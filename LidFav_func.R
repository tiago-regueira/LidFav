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

calcula_indicadores_favelas <- function(nomes_favelas) {
  
  library(sf)
  library(dplyr)
  library(terra)
  library(lidR)
  
  #### Dados base
  favelas <- st_read("./1_Insumos/Favelas_Habita_20260307.gpkg")
  
  malha_mdt_mds <- st_read("./1_Insumos/Quadricula_MDT_MDS_2020.geojson") |>
    st_set_crs('EPSG:31983')
  
  pasta_base <- "./2_Processamento"
  
  processa_uma_favela <- function(nome_favela) {
    
    message("Processando: ", nome_favela)
    
    favela_selecionada <- favelas |>
      filter(nome == nome_favela)
    
    favela_buffer <- st_buffer(favela_selecionada, 120)
    
    pasta_favela <- file.path(pasta_base, nome_favela)
    
    dir.create(pasta_favela, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(pasta_favela, "MDT"), showWarnings = FALSE)
    dir.create(file.path(pasta_favela, "MDS"), showWarnings = FALSE)
    dir.create(file.path(pasta_favela, "Resultados"), showWarnings = FALSE)
    
    #### Seleção de quadrículas
    malha_selecao <- st_intersection(malha_mdt_mds, favela_buffer)
    
    malha_download <- malha_selecao |>
      st_drop_geometry() |>
      pull(cd_quadricula)
    
    extensoes <- c('MDS', 'MDT')
    
    #### Download
    for (extensao in extensoes) {
      
      pasta_destino <- file.path(pasta_favela, extensao)
      
      for (codigo in malha_download) {
        
        url <- paste0(
          "http://download.geosampa.prefeitura.sp.gov.br/PaginasPublicas/",
          "downloadArquivo.aspx?orig=DownloadMapaArticulacao",
          "&arq=", extensao, "_2020%5C", codigo, ".zip",
          "&arqTipo=MAPA_ARTICULACAO"
        )
        
        zip_destino <- file.path(pasta_destino, paste0(codigo, ".zip"))
        
        download.file(url, zip_destino, mode = "wb")
        unzip(zip_destino, exdir = pasta_destino)
        
        Sys.sleep(1)
      }
    }
    
    #### Processamento
    for (extensao in extensoes) {
      
      pasta_destino <- file.path(pasta_favela, extensao)
      
      ctg <- readLAScatalog(pasta_destino) |>
        clip_roi(favela_buffer)
      
      if (extensao == "MDS") {
        r <- rasterize_canopy(ctg, res = 1, p2r())
      }
      
      if (extensao == "MDT") {
        ctg$Classification[ctg$Classification == 8] <- 2L
        r <- rasterize_terrain(ctg, res = 1, tin())
      }
      
      r <- mask(r, vect(favela_selecionada))
      
      writeRaster(
        r,
        file.path(pasta_destino, paste0(extensao, "_mosaico.tif")),
        overwrite = TRUE
      )
    }
    
    #### Remover LAS/LAZ
    for (extensao in extensoes) {
      pasta_destino <- file.path(pasta_favela, extensao)
      
      arquivos_las <- list.files(
        pasta_destino,
        pattern = "\\.(las|laz)$",
        full.names = TRUE
      )
      
      file.remove(arquivos_las)
    }
    
    #### CHM
    mds <- rast(file.path(pasta_favela, "MDS", "MDS_mosaico.tif"))
    mdt <- rast(file.path(pasta_favela, "MDT", "MDT_mosaico.tif"))
    
    mds_favela <- crop(mds, favela_selecionada)
    mdt_favela <- crop(mdt, favela_selecionada)
    
    chm <- mds_favela - mdt_favela
    
    writeRaster(
      chm,
      file.path(pasta_favela, "Resultados", "CHM.tif"),
      overwrite = TRUE
    )
    
    #### Indicadores
    edificacoes <- chm > 2
    
    area_favela <- as.numeric(st_area(favela_selecionada))
    area_pixel <- prod(res(chm))
    
    n_pixels_edif <- global(edificacoes, "sum", na.rm = TRUE)[1]
    area_edificada <- n_pixels_edif * area_pixel
    
    TI <- area_edificada / area_favela
    
    altura_edif <- mask(chm, edificacoes, maskvalues = FALSE)
    volume_pixel <- altura_edif * area_pixel
    volume_total <- global(volume_pixel, "sum", na.rm = TRUE)[1,1]
    
    IV <- volume_total / area_favela
    
    data.frame(
      favela = nome_favela,
      area_edificada = area_edificada,
      TI = TI,
      volume_total = volume_total,
      IV = IV
    )
  }
  
  #### Iteração
  resultados <- lapply(nomes_favelas, processa_uma_favela) |>
    bind_rows()
  
  return(resultados)
}
