#----- Pacotes

if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(plotly)) {install.packages("plotly"); library(plotly)}
if(!require(knitr)) {install.packages("knitr"); library(knitr)}
if(!require(likert)){install.packages("likert"); require(likert)}
if(!require(factoextra)){install.packages("factoextra"); require(factoextra)}
if(!require(dendextend)){install.packages("dendextend"); require(dendextend)}
if(!require(heatmaply)){install.packages("heatmaply"); require(heatmaply)}
if(!require(psych)){install.packages("psych"); require(psych)}
if(!require(psy)){install.packages("psy"); require(psy)}
if(!require(nFactors)){install.packages("nFactors"); require(nFactors)}

#----- Funções

tab <- function(x){cbind(table(x), prop.table(table(x))*100)}

Val1 <- function(MDF1){
  faPR1 <- principal(MDF1, nfactors=1,  rotate="varimax")
  itens<- dim(MDF1)[2]
  AVE<- mean(faPR1$communality)
  KMO<- kmo(MDF1)$KMO
  AC<-cronbach(MDF1)$alpha
  DG<- (sum(principal(MDF1,rotate="varimax")$loadings)^2)/ ((sum(principal(MDF1,rotate="varimax")$loadings)^2) + (sum(1-(principal(MDF1,rotate="varimax")$communality)^2)))
  tab<-cbind(itens, AVE, AC, DG, KMO, 1,1,1,1)
  colnames(tab)<- c("itens","AVE","AC","DG","KMO","noc","naf","nparallel","nkaiser")
  tab
}

Val <- function(MDF1){
  nf<- function(boi){
    boi1<-eigen(cor(boi, method="spearman"), TRUE)
    eig <-  boi1$values
    N <- dim(boi)[1]
    res <- nSeScree(eig)
    eigenvalues<- eigen(cor(boi,method="spearman"), TRUE)$values
    nsubjects <- dim(boi)[1]
    variables <- length(eigenvalues) # Computes the number of variables
    rep <- 1000 # Number of replications for PA analysis
    cent <- 0.95 # Centile value of PA analysis
    aparallel <- parallel(var = variables,
                          subject = nsubjects,
                          rep = rep,
                          cent = cent)$eigen$qevpea # The 95 centile
    results <- nScree(x=eigenvalues, aparallel=aparallel)
    results$Components}
  faPR1 <- principal(MDF1, nfactors=1,  rotate="varimax")
  itens<- dim(MDF1)[2]
  AVE<- mean(faPR1$communality)
  KMO<- kmo(MDF1)$KMO
  AC<-cronbach(MDF1)$alpha
  DG<- (sum(principal(MDF1,rotate="varimax")$loadings)^2)/ ((sum(principal(MDF1,rotate="varimax")$loadings)^2) + (sum(1-(principal(MDF1,rotate="varimax")$communality)^2)))
  ND<- nf(MDF1)
  cbind(itens, AVE, AC, DG, KMO, ND)
}

kmo <- function(x)
{
  x <- subset(x, complete.cases(x)) # Omit missing values
  r <- cor(x) # Correlation matrix
  r2 <- r^2 # Squared correlation coefficients
  i <- solve(r) # Inverse matrix of correlation matrix
  d <- diag(i) # Diagonal elements of inverse matrix
  p2 <- (-i/sqrt(outer(d, d)))^2 # Squared partial correlation coefficients
  diag(r2) <- diag(p2) <- 0 # Delete diagonal elements
  KMO <- sum(r2)/(sum(r2)+sum(p2))
  MSA <- colSums(r2)/(colSums(r2)+colSums(p2))
  return(list(KMO=KMO, MSA=MSA))
}

fa <- function(x){
  Mat<- cbind(principal(x,1, rotate="varimax")$loadings
              ,principal(x,1, rotate="varimax")$communality
              ,principal(x,1, rotate="varimax")$weights)
  Mat1<- Mat[order(abs(Mat[,1]), decreasing = TRUE),]
  colnames(Mat1)<- c("C.F.","Com.","Peso")
  Mat1
}


table.aux <- function(x,y,  more=F){
  tab<- table(factor(x), factor(y))
  tab1<- round(prop.table(tab,1)*100,2)
  tab2<- cbind(tab[,1], tab1[,1], tab[,2], tab1[,2])
  p_valor<- round(chisq.test(tab)$p.value,3)
  tab5<- cbind(tab2,p_valor)
  colnames(tab5) <- c(rep(c("N", "%"), 2), "Valor-p")
  tab5}

basic.stats <- function(x,more=F) {
  stats <- list()
  
  clean.x <- x[!is.na(x)]
  
  #stats$N <- length(x)
  #stats$NAs <- stats$N-length(clean.x)
  stats$N_validos <- round(length(clean.x),3)
  stats$Média <- round(mean(clean.x),3)
  stats$E.P <- round(sd(clean.x)/sqrt(length(clean.x)),3)
  stats$Q1 <- round(fivenum(clean.x)[2],3)
  stats$Q2<- round(fivenum(clean.x)[3],3)
  stats$Q3 <- round(fivenum(clean.x)[4],3)
  t1<- unlist(stats)
  names(t1)<- c("N válidos", "Média", "E.P", "1Q", "2Q", "3Q")
  t1
}


basic <- function(x,more=F) {
  stats <- list()
  
  clean.x <- x[!is.na(x)]
  
  #stats$N <- length(x)
  #stats$NAs <- stats$N-length(clean.x)
  stats$N_validos <- round(length(clean.x),3)
  stats$Média <- round(mean(clean.x),3)
  stats$D.P <- round(sd(clean.x),3)
  stats$Mín. <- round(min(clean.x),3)
  stats$Q1 <- round(fivenum(clean.x)[2],3)
  stats$Q2<- round(fivenum(clean.x)[3],3)
  stats$Q3 <- round(fivenum(clean.x)[4],3)
  stats$Máx. <- round(max(clean.x),3)
  t1<- unlist(stats)
  names(t1)<- c("N válidos", "Média", "D.P.","Mín.", "1Q", "2Q", "3Q","Máx.")
  t1
}

m.b <- function(x){
  media<-c()
  for (i in 1:1000) {
    boot<-sample(x, replace=TRUE)
    media[i] <- mean(boot, na.rm =TRUE)
  }
  LI<- quantile(media, probs=c(0.025), na.rm =TRUE)
  LS<- quantile(media, probs=c(0.975), na.rm =TRUE)
  valores<- c(LI,LS)
  valores
}

basic.np <- function(x,more=F) {
  stats <- list()
  
  clean.x <- x[!is.na(x)]
  
  stats$N_validos <- round(length(clean.x),3)
  stats$Média <- round(mean(clean.x),3)
  stats$D.P. <- round(sd(clean.x),3)
  stats$L.I.<- round(as.numeric(m.b(clean.x)[1]),3)
  stats$L.S.<- round(as.numeric(m.b(clean.x)[2]),3)
  stats$Mín. <- round(min(clean.x),3)
  stats$Q1 <- round(fivenum(clean.x)[2],3)
  stats$Q2<- round(fivenum(clean.x)[3],3)
  stats$Q3 <- round(fivenum(clean.x)[4],3)
  stats$Máx. <- round(max(clean.x),3)
  t1<- unlist(stats)
  names(t1)<- c("N válidos", "Média", "D.P.","L.I", "L.S","Mín.", "1Q", "2Q", "3Q","Máx.")
  t1
}

plotar_likert <- function(data, titulo) {
  plot(likert(data), centered = T, plot.percents = F) + 
    ggtitle(titulo)+
    ylab("Porcentagem") +
    guides(fill=guide_legend("Escala Likert"))
  
}

kruskal_b <- function(y, z, more=F){
  tab <-matrix(NA, length(levels(factor(z))), 6)
  for(i in 1:length(levels(factor(z)))){ 
    desc<- tapply(y, factor(z),  basic.stats)[i]
    desc1<- unlist(desc)
    for(j in 1:6){ 
      tab[i,j] <-desc1[j]
    }
  }
  p_valor<- rep(kruskal.test(y~factor(z))$p.value, length(levels(factor(z))))
  tab<- cbind(tab, p_valor)
  
  tab_new <- c(tab[1, c(2,3)], tab[2, c(2,3)], tab[3, c(2,3, 7)]) 
  
}


#----- Leitura base de dados

dados <- read.csv2("wiki4HE.csv", na.strings = "?")
dados %>% head() 
dados %>% dim()

#----- Análise de dados faltantes

plot_na <- dados %>% 
  dplyr::select(-c("OTHER_POSITION", "OTHERSTATUS")) %>% 
  summarise_all(function(x) sum(is.na(x)) / n()) %>% 
  gather("variavel", "prop_na") %>% 
  filter(prop_na > 0) %>% 
  arrange(prop_na) %>% 
  ggplot(aes(x = variavel, y = prop_na*100)) +
  ylab("Percentual de dados ausentes") +
  xlab("Variáveis") +
  geom_col(fill="#2EB4C4") +
  coord_flip()

ggplotly(plot_na) %>%
  config(displayModeBar = F)

ind_aior_10 <- table(apply(dados, 1, function(x) sum(is.na(x)) / length(x)) > 10)

#----- Imputação dos dados

dados[ ,c(1, 5, 11:53)] <- apply(dados[ ,c(1, 5, 11:53)], 2, function(x) ifelse(is.na(x), round(mean(x, na.rm=T)), x))

#----- Análise Exploratória

tab_desc1 <- dados %>% 
  dplyr::select(c("GENDER", "DOMAIN", "PhD", "UNIVERSITY", 
                  "UOC_POSITION", "OTHER_POSITION",
                  "OTHERSTATUS", "USERWIKI")) %>% 
  apply(., 2, tab) %>% 
  do.call(rbind,.) %>% 
  round(.,2) %>% 
  as.data.frame() 

colnames(tab_desc1) <- c("N", "%")

rownames(tab_desc1) <- c(
  "Masculino",
  "Feminino",
  
  "Artes e Humanidades",
  "Ciências",
  "Ciências da Saúde",
  "Engenharia e Arquitetura",
  "Direito e Política",
  "6",
  
  "PhD-Não",
  "PhD-Sim",
  
  "Uni-UOC",
  "Uni-UPF",
  
  "Pos_Acad: Professor",
  "Pos_Acad: Associado",
  "Pos_Acad: Assistente",
  "Pos_Acad: Palestrante",
  "Pos_Acad: instrutor",
  "Pos_Acad: Adjunto",
  
  
  "Outra_pos: sim",
  "Outra_pos: não",
  
  "Outra_Pos: Professor",
  "Outra_Pos: Associado",
  "Outra_Pos: Assistente",
  "Outra_Pos: Palestrante",
  "Outra_Pos: instrutor",
  "Outra_Pos: Adjunto",
  "7",
  
  "Usa_wiki: não",
  "Usa_wiki: sim"
)

tab_desc2 <- dados %>% 
  dplyr::select(c("AGE", "YEARSEXP")) %>% 
  apply(., 2, basic) %>% 
  t() %>% 
  as.data.frame()  

rownames(tab_desc2) <- c("Idade", "Anos de Experiência")

tab_desc1 
tab_desc2 

#----- Gráficos análise exploratória

# Organizando o nome das variáveis para o plot

dados$GENDER <- factor(dados$GENDER, labels = c("Masculino", "Feminino"))

dados$DOMAIN <- factor(dados$DOMAIN, labels = c("Artes e Humanidades",
                                                "Ciências",
                                                "Ciências da Saúde",
                                                "Engenharia e Arquitetura",
                                                "Direito e Política",
                                                "6"))

dados$PhD <- factor(dados$PhD, labels = c("Não", "Sim"))

dados$UNIVERSITY <- factor(dados$UNIVERSITY, labels = c("UOC", "UPF"))

dados$UOC_POSITION <- factor(dados$UOC_POSITION, labels = c("Professor",
                                                            "Associado",
                                                            "Assistente",
                                                            "Palestrante",
                                                            "Instrutor",
                                                            "Adjunto"))

dados$OTHER_POSITION <- factor(dados$OTHER_POSITION, labels = c("Sim", "Não"))

dados$OTHERSTATUS <- factor(dados$OTHERSTATUS, labels = c("Professor",
                                                          "Associado",
                                                          "Assistente",
                                                          "Palestrante",
                                                          "Instrutor",
                                                          "Adjunto",
                                                          "7"))

dados$USERWIKI <- factor(dados$USERWIKI, labels = c("Não", "Sim"))

# Plot 

plot1 <- dados %>% 
  drop_na(GENDER) %>% 
  ggplot(aes(as.factor(GENDER))) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill="#2EB4C4") + 
  xlab("Sexo") + 
  ylab("Frequência")+ 
  ggtitle("") + 
  scale_y_continuous(labels=scales::percent)


plot2 <- dados %>% 
  drop_na(DOMAIN) %>% 
  ggplot(aes(as.factor(DOMAIN))) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill="#2EB4C4") + 
  xlab("Dimíno") + 
  ylab("Frequência")+ 
  ggtitle("Dimíno") + 
  scale_y_continuous(labels=scales::percent) +
  coord_flip()

plot3 <- dados %>% 
  drop_na(PhD) %>% 
  ggplot(aes(as.factor(PhD))) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill="#2EB4C4") + 
  xlab("PhD") + 
  ylab("Frequência")+ 
  ggtitle("") + 
  scale_y_continuous(labels=scales::percent)

plot4 <- dados %>% 
  drop_na(UNIVERSITY) %>% 
  ggplot(aes(as.factor(UNIVERSITY))) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill="#2EB4C4") + 
  xlab("Universidade") + 
  ylab("Frequência")+ 
  ggtitle("Universidade") + 
  scale_y_continuous(labels=scales::percent)

plot5 <- dados %>% 
  drop_na(UOC_POSITION) %>% 
  ggplot(aes(as.factor(UOC_POSITION))) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill="#2EB4C4") + 
  xlab("UOC_POSITION") + 
  ylab("Frequência")+ 
  ggtitle("UOC_POSITION") + 
  scale_y_continuous(labels=scales::percent) +
  coord_flip()

plot6 <- dados %>% 
  drop_na(OTHER_POSITION) %>% 
  ggplot(aes(as.factor(OTHER_POSITION))) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill="#2EB4C4") + 
  xlab("OTHER_POSITION") + 
  ylab("Frequência")+ 
  ggtitle("OTHER_POSITION") + 
  scale_y_continuous(labels=scales::percent)

plot7 <- dados %>% 
  drop_na(OTHERSTATUS) %>% 
  ggplot(aes(as.factor(OTHERSTATUS))) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill="#2EB4C4") + 
  xlab("OTHERSTATUS") + 
  ylab("Frequência")+ 
  ggtitle("OTHERSTATUS") + 
  scale_y_continuous(labels=scales::percent) +
  coord_flip()

plot8 <- dados %>% 
  drop_na(USERWIKI) %>% 
  ggplot(aes(as.factor(USERWIKI))) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill="#2EB4C4") + 
  xlab("Usa Wikipédia") + 
  ylab("Frequência")+ 
  ggtitle("Usa Wikipédia") + 
  scale_y_continuous(labels=scales::percent)

subplot(ggplotly(plot1), ggplotly(plot3)) %>%
  config(displayModeBar = F)

ggplotly(plot2) %>%
  config(displayModeBar = F)

subplot(ggplotly(plot4), ggplotly(plot6)) %>%
  config(displayModeBar = F)

ggplotly(plot5) %>%
  config(displayModeBar = F)

ggplotly(plot7) %>%
  config(displayModeBar = F)

ggplotly(plot8) %>%
  config(displayModeBar = F)

#----- Análise Exploratória dos Itens 

set.seed(220120)
tab_desc3 <- dados %>% 
  dplyr::select(-c("GENDER", "DOMAIN", "PhD", "UNIVERSITY", 
                   "UOC_POSITION", "OTHER_POSITION",
                   "OTHERSTATUS", "USERWIKI", "AGE", "YEARSEXP")) %>% 
  apply(., 2, basic.np) %>% 
  t() %>% 
  round(.,2) %>% 
  as.data.frame() 

tab_desc3 

# Gráfico Likert dos Itens

data_aux_plot <- dados %>% 
  dplyr::select(-c("GENDER", "DOMAIN", "PhD", "UNIVERSITY", 
                   "UOC_POSITION", "OTHER_POSITION",
                   "OTHERSTATUS", "USERWIKI", "AGE", "YEARSEXP")) 

data_aux_plot <- lapply(data_aux_plot, factor, levels = c("1", "2", "3", "4", "5")) %>% as.data.frame()

plot_likert <- plotar_likert(data_aux_plot, "")

ggplotly(plot_likert) %>%
  config(displayModeBar = F)

#----- Machine Learning - Modelo Não Supervisionado

# Análise de Agrupamento

data_cluster <- dados %>% 
  dplyr::select(-c("GENDER", "DOMAIN", "PhD", "UNIVERSITY", 
                   "UOC_POSITION", "OTHER_POSITION",
                   "OTHERSTATUS", "USERWIKI"))

data_cluster <- apply(data_cluster, 2, function(x) scale(x))

# Definindo o número de grupos

set.seed(1654651)
wss <- (nrow(data_cluster)-1)*sum(apply(data_cluster,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(data_cluster,
                                     centers=i)$withinss)

plot(1:15, wss, type="b", xlab="Nº de Clusters",
     ylab="Soma dos quadrados dos grupos", pch = 16, col="#2EB4C4")

dados$Grupo <- kmeans(data_cluster, 3, iter.max = 10, nstart = 1)$cluster 

tab(dados$Grupo) %>% as.data.frame()

plot_grupo <- dados %>% 
  ggplot(aes(as.factor(Grupo))) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill="#2EB4C4") + 
  xlab("Grupo") + 
  ylab("Frequência")+ 
  ggtitle("Grupo") + 
  scale_y_continuous(labels=scales::percent)

ggplotly(plot_grupo)

### Comparação entre os grupos 

tab_comp <- rbind(
  kruskal_b(dados$ENJ1, dados$Grupo),
  kruskal_b(dados$ENJ2, dados$Grupo)) %>% 
  as.data.frame() %>% 
  round(.,2)

colnames(tab_comp) <- c(rep(c("Média", "E.P."), 3), "Valor-p")
rownames(tab_comp) <- c("ENJ1", "ENJ2") 

tab_comp %>% as.data.frame() 

dados$Grupo <- factor(dados$Grupo, labels = c("1", "2", "3"))

box1 <- ggplot(dados, aes(x=Grupo, y=ENJ1)) +
  geom_boxplot(fill='#2EB4C4', color="black")+
  theme_classic()

box2 <- ggplot(dados, aes(x=Grupo, y=ENJ2)) +
  geom_boxplot(fill='#2EB4C4', color="black")+
  theme_classic()

tab_wiki_grupo <- table.aux(dados$Grupo, dados$USERWIKI)
tab_wiki_grupo %>% as.data.frame() 

#----- Análise Fatorial

# Definindo os constructos

UP <- dados %>% dplyr::select(c("PU1", "PU2", "PU3"))
FUP <- dados %>% dplyr::select(c("PEU1", "PEU2", "PEU3"))
PP <- dados %>% dplyr::select(c("ENJ1", "ENJ2"))
QU <- dados %>% dplyr::select(c("Qu1", "Qu2", "Qu3", "Qu4", "Qu5"))
VS <- dados %>% dplyr::select(c("Vis1", "Vis2", "Vis3"))
IS <- dados %>% dplyr::select(c("Im1", "Im2", "Im3"))
AC <- dados %>% dplyr::select(c("SA1", "SA2", "SA3"))
UC <- dados %>% dplyr::select(c("Use1", "Use2", "Use3", "Use4", "Use5"))
PF <- dados %>% dplyr::select(c("Pf1", "Pf2", "Pf3"))
RT <- dados %>% dplyr::select(c("JR1", "JR2"))
IC <- dados %>% dplyr::select(c("Inc1", "Inc2", "Inc3", "Inc4"))
EXP <- dados %>% dplyr::select(c("Exp1", "Exp2", "Exp3", "Exp4", "Exp5"))

# Modelo Inicial #

tab_fat1 <- rbind(
  fa(UP),
  fa(FUP),
  fa(PP),
  fa(QU), #
  fa(VS),
  fa(IS),
  fa(AC),
  fa(UC),
  fa(PF),
  fa(RT),
  fa(IC),
  fa(EXP)) %>%  as.data.frame()

# Modelo Final #

tab_fat2 <- rbind(
  fa(UP),
  fa(FUP),
  fa(PP),
  fa(QU[,-c(4,5)]), #
  fa(VS),
  fa(IS),
  fa(AC),
  fa(UC),
  fa(PF),
  fa(RT),
  fa(IC),
  fa(EXP)) %>% as.data.frame()

tab_fat2 

tab_Valt2 <- rbind(
  Val(UP),
  Val(FUP),
  Val1(PP),
  Val(QU[,-c(4,5)]), #
  Val(VS),
  Val(IS),
  Val(AC),
  Val(UC),
  Val(PF),
  Val1(RT),
  Val(IC),
  Val(EXP)) %>% as.data.frame()

tab_Val <- tab_Valt2[,-c(6,8,9)] %>% round(.,2)
rownames(tab_Val) <- c("UP","FUP","PP","QU", "VS", "IS", "AC", "UC", 
                       "PF", "RT", "IC", "EXP") 
tab_Val 

dados$UPi <- apply(UP, 1, mean)
dados$FUPi <- apply(FUP, 1, mean)
dados$PPi <- apply(PP, 1, mean)
dados$QUi <- apply(QU[,-c(4,5)], 1, mean)
dados$VSi <- apply(VS, 1, mean)
dados$ISi <- apply(IS, 1, mean)
dados$ACi <- apply(AC, 1, mean)
dados$UCi <- apply(UC, 1, mean)
dados$PFi <- apply(PF, 1, mean)
dados$RTi <- apply(RT, 1, mean)
dados$ICi <- apply(IC, 1, mean)
dados$EPXi <- apply(EXP, 1, mean)

tab_desc_ind <- dados %>% 
  dplyr::select(c("UPi", "FUPi", "PPi", "QUi", "VSi", "ISi", 
                  "ACi", "UCi", "PFi", "RTi", "ICi", "EPXi")) %>% 
  apply(., 2, basic) %>% 
  t() %>% 
  as.data.frame() %>% 
  round(.,2)

tab_desc_ind 
