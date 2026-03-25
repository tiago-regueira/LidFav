# --------------------------------------------------------------------------- #

# Importando os pacotes

library( tidyverse )
library( lidR )
library( terra )
library( sf )

# --------------------------------------------------------------------------- #

# 0 Discussão

# --------------------------------------------------------------------------- #

## 1 Insumos

calcula_indicadores_favelas_nucleos <- function( nomes ) {
  
  library(sf)
  library(dplyr)
  library(terra)
  library(lidR)
  
  #### Dados base
  favelas <- st_read("./1_Insumos/Favelas_Habita_20260307.gpkg") |>
    mutate(tipo = "favela")
  
  nucleos <- st_read("./1_Insumos/Nucleos_Habita_20260325.gpkg") |>
    mutate(tipo = "nucleo")
  
  assentamentos <- bind_rows(favelas, nucleos)
  
  nomes_validos <- assentamentos$nome
  
  nomes_invalidos <- setdiff(nomes, nomes_validos)
  
  if (length(nomes_invalidos) > 0) {
    stop(
      paste0(
        "Os seguintes nomes não foram encontrados nas bases de favelas/núcleos:\n",
        paste(nomes_invalidos, collapse = ", ")
      )
    )
  }
  
  malha_mdt_mds <- st_read("./1_Insumos/Quadricula_MDT_MDS_2020.geojson") |>
    st_set_crs('EPSG:31983')
  
  pasta_base <- "./2_Processamento"
  
  processa_um <- function(nome_area) {
    
    message("Processando: ", nome_area)
    
    area_sel <- assentamentos |>
      filter(nome == nome_area)
    
    if (nrow(area_sel) == 0) {
      warning("Nome não encontrado: ", nome_area)
      return(NULL)
    }
    
    tipo_area <- area_sel$tipo[1]
    
    pasta_area <- file.path(pasta_base, nome_area)
    pasta_result <- file.path(pasta_area, "Resultados")
    rds_path <- file.path(pasta_result, "indicadores.rds")
    
    #### CACHE: se já existe RDS, retorna direto
    if (file.exists(rds_path)) {
      return(readRDS(rds_path))
    }
    
    favela_buffer <- st_buffer(area_sel, 120)
    
    dir.create(pasta_area, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(pasta_area, "MDT"), showWarnings = FALSE)
    dir.create(file.path(pasta_area, "MDS"), showWarnings = FALSE)
    dir.create(pasta_result, showWarnings = FALSE)
    
    #### Seleção de quadrículas
    malha_selecao <- st_intersection(malha_mdt_mds, favela_buffer)
    
    malha_download <- malha_selecao |>
      st_drop_geometry() |>
      pull(cd_quadricula)
    
    extensoes <- c("MDS", "MDT")
    
    #### Download (com verificação)
    for (extensao in extensoes) {
      
      pasta_destino <- file.path(pasta_area, extensao)
      
      for (codigo in malha_download) {
        
        zip_destino <- file.path(pasta_destino, paste0(codigo, ".zip"))
        
        if (!file.exists(zip_destino)) {
          
          url <- paste0(
            "http://download.geosampa.prefeitura.sp.gov.br/PaginasPublicas/",
            "downloadArquivo.aspx?orig=DownloadMapaArticulacao",
            "&arq=", extensao, "_2020%5C", codigo, ".zip",
            "&arqTipo=MAPA_ARTICULACAO"
          )
          
          download.file(url, zip_destino, mode = "wb")
          unzip(zip_destino, exdir = pasta_destino)
          
          Sys.sleep(1)
        }
      }
    }
    
    #### Processamento
    for (extensao in extensoes) {
      
      pasta_destino <- file.path(pasta_area, extensao)
      
      ctg <- readLAScatalog(pasta_destino) |>
        clip_roi(favela_buffer)
      
      if (extensao == "MDS") {
        r <- rasterize_canopy(ctg, res = 1, p2r())
      }
      
      if (extensao == "MDT") {
        ctg$Classification[ctg$Classification == 8] <- 2L
        r <- rasterize_terrain(ctg, res = 1, tin())
      }
      
      r <- mask(r, vect(area_sel))
      
      writeRaster(
        r,
        file.path(pasta_destino, paste0(extensao, "_mosaico.tif")),
        overwrite = TRUE
      )
    }
    
    #### Limpeza LAS/LAZ
    for (extensao in extensoes) {
      
      pasta_destino <- file.path(pasta_area, extensao)
      
      arquivos_las <- list.files(
        pasta_destino,
        pattern = "\\.(las|laz)$",
        full.names = TRUE
      )
      
      file.remove(arquivos_las)
    }
    
    #### CHM
    mds <- rast(file.path(pasta_area, "MDS", "MDS_mosaico.tif"))
    mdt <- rast(file.path(pasta_area, "MDT", "MDT_mosaico.tif"))
    
    mds_area <- crop(mds, area_sel)
    mdt_area <- crop(mdt, area_sel)
    
    chm <- mds_area - mdt_area
    
    writeRaster(
      chm,
      file.path(pasta_result, "CHM.tif"),
      overwrite = TRUE
    )
    
    #### Indicadores
    edificacoes <- chm > 2
    
    area_total <- as.numeric(st_area(area_sel))
    area_pixel <- prod(res(chm))
    
    n_pixels <- global(edificacoes, "sum", na.rm = TRUE)[1]
    area_construida <- n_pixels * area_pixel
    
    TO <- area_construida / area_total
    
    altura_edif <- mask(chm, edificacoes, maskvalues = FALSE)
    volume_pixel <- altura_edif * area_pixel
    volume_total <- global(volume_pixel, "sum", na.rm = TRUE)[1,1]
    
    IV <- volume_total / area_total
    
    resultado <- data.frame(
      nome = nome_area,
      tipo = tipo_area,
      area_construida = area_construida,
      TO = TO,
      volume_total = volume_total,
      IV = IV
    )
    
    #### Salva cache
    saveRDS(resultado, rds_path)
    
    return(resultado)
  }
  
  #### Iteração
  resultados <- lapply(nomes, processa_um) |>
    bind_rows()
  
  return(resultados)
}