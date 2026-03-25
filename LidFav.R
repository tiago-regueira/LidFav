# --------------------------------------------------------------------------- #

# 0 Discussão

# --------------------------------------------------------------------------- #

## 1 Insumos

setwd( this.path::this.dir() )

source( "./LidFav_func.R" )

assentamentos_pesquisa <- c( 'JARDIM NOVA PANTANAL',
                             'JARDIM NAKAMURA I',
                             'VILA RUBI',
                             'JARDIM IMBUIAS I',
                             'CHÁCARA DO CONDE II',
                             'JARDIM ARACATI I' )

indicadores <- calcula_indicadores_favelas_nucleos( assentamentos_pesquisa )

indicadores
