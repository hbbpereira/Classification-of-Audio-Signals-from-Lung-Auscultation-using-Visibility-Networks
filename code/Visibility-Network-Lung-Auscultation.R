# ---------------------------------------------------------------------------------
# Redes de visibilidade para ausculta pulmonar              -----------------------
# Visibility networks from lung auscultation                -----------------------
# Cálculo de propriedades de redes selecionadas             -----------------------
# Calculation of properties of selected networks            -----------------------
# Autores(authors): Dirceu Melo                             -----------------------
#                   Hernane Borges de Barros Pereira        -----------------------
#            Criação(creation): 08/08/2023                  -----------------------
# Última atualização(last update): 07/10/2024               -----------------------
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# Instalação de pacotes necessários                         -----------------------
# Installing required packages.                             -----------------------
# ---------------------------------------------------------------------------------

# install.packages(c("igraph","readr","tuneR","cowplot","poweRlaw","brainGraph"))
# library(igraph)
# library('somebm')
# library(cowplot)
# library(poweRlaw)
# library(readr)
# library(tuneR)
# library(brainGraph)

# ---------------------------------------------------------------------------------
# Preparação de diretório de trabalho                       -----------------------
# Working directory preparation                             -----------------------
# ---------------------------------------------------------------------------------

# Definição de diretório de trabalho
# Working directory definition
setwd("/Users/hbbpereira/Downloads/Ausculta")

# Carregar os nomes dos arquivos com extensão WAV
# Load file names with WAV extension
input <- dir(pattern='.wav') 

# ---------------------------------------------------------------------------------
# Transformação de wav em séries de desvios padrões         -----------------------
# Transformation of wav into series of standard deviations  -----------------------
# ---------------------------------------------------------------------------------

# Definição de variáveis
# Definition of variables
L <- length(input)
train_audio <- NULL
train_audio_norm <- NULL
variancia <- NULL
variancia_acumulada <- NULL
serie_variancias <- NULL

for (i in 1:L){
  
  # Quantidade de segundos do arquivo WAV; arquivos de ausculta pulmonar t = 19s
  # Number of seconds in the WAV file; lung auscultation files t = 19s
  N <- 837900
  
  # Janela para sample rate de 44,1 khz
  # # Window for 44.1 kHz sample rate
  m <- 441 
  r <- 1
  k <- 1
  train_audio[[i]] = readWave(input[i], from=1, to=N, units=c("samples"))
  train_audio_norm[[i]] <- normalize(train_audio[[i]], unit=c("8"), center=TRUE,level=1, rescale=TRUE)
  
  while (r < N){
    variancia[[k]] <- var(train_audio_norm[[i]][c(r:(r+m-1))]@left)
    r <- r+m
    k <- k+1
  }
  
  serie_variancias[[i]] <- variancia
  
}

# ---------------------------------------------------------------------------------
# Transformação de uma série de desvios padrões em rede de visibilidade -----------
# transformation of a series of standard deviations into visibility network  ------
# ---------------------------------------------------------------------------------

# Definição de variáveis
# Definition of variables
y <- serie_variancias
x <- 1:length(variancia)
tab.pontos <- list()
tab.pontos.geral <- list()
tab.geral <- list()
tg.max <- list()
tx.var <- list()
tam <- length(variancia)
g <- list()
vert <- list()
arest <- list()
dens <- list()
gr_m <- list()
mod <- list()
diam <- list()
cmm <- list()
c.aglom <- list()
efici_glob <- list()
fit1 <- list()
gama1 <- list()
assort <- list()
central_grau <- list()
central_prox <- list()
centr_autovetor <- list()
centr_autovalor <- list()

# Aplicação do critério de visibilidade
# Application of visibility criteria
for(j in 1:L){
  
  tab.pontos <- NULL
  
  for(i in 1:(tam-1)){
    tg.max [[j]] <- (as.numeric(y[[j]][i+1])-as.numeric(y[[j]][i]))/(x[i+1]-x[i]) 
    tab.pontos<- rbind(tab.pontos, c(x[i],x[i+1]))
    
    for(k in (i+1):tam){
      tx.var[[j]] <- (as.numeric(y[[j]][k])-as.numeric(y[[j]][i]))/(x[k]-x[i])
      if(as.numeric(tx.var[[j]]) > as.numeric(tg.max[[j]])){
        tg.max[[j]] <- tx.var[[j]]
        tab.pontos <- rbind(tab.pontos, c(x[i],x[k]))
      } 
    }
  }
  
  # Estrutura de dados geral para as redes de visibilidade e suas propriedades
  # General data structure for visibility networks and their properties 
  tab.pontos.geral[[j]] <- tab.pontos
  g[[j]] <- graph_from_edgelist(tab.pontos.geral[[j]], directed = FALSE)
  
  # Calcular as propriedades das redes de visibilidade
  # Calculate the properties of visibility networks
  vert <- vcount(g[[j]])
  arest <- ecount(g[[j]])
  gr_m <-mean(degree(g[[j]],mode="all"))
  dens <- edge_density(g[[j]], loops = FALSE)
  mod <-modularity(cluster_louvain(g[[j]]))
  diam <- diameter(g[[j]]) 
  cmm <- mean_distance(g[[j]],directed = FALSE)
  c.aglom <- transitivity(g[[j]], type = "average")
  efici_glob <- efficiency(g[[j]], type=c("global")) 
  fit1 <- fit_power_law(degree(g[[j]])+1, 10)
  gama1 <- fit1$alpha
  assort <- assortativity_degree(g[[j]], directed=F) 
  central_grau <- centr_degree(g[[j]], mode="in", normalized=T)$centralization 
  central_prox <- centr_clo(g[[j]], mode="all", normalized=T)$centralization 
  centr_autovetor <- centr_eigen(g[[j]], directed=T, normalized=T)$centralization
  centr_autovalor <- centr_eigen(g[[j]], directed=T, normalized=T)$value

  # Preencher a estrutura de dados geral 
  # Populate the general data structure 
  tab.prop.redes.todas<- data.frame(vert, arest,gr_m,dens, mod, diam, cmm,c.aglom, efici_glob,gama1,assort,central_grau,central_prox,centr_autovetor, centr_autovalor)
  tab.geral[[j]]<-rbind(tab.prop.redes.todas)
  dados <- do.call(rbind.data.frame, tab.geral)
  
}

# Salvar em CSV as propriedades das redes de visibilidade
# Save visibility network properties to CSV
write.csv2(dados,"visibility-network-properties.csv")
