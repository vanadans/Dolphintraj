##########################analisis manuscrito dolphin movements######################

# Start ----
#setwd("D:/Vana/Documentos/manuscritos/Random walk/PlosONE")
#setwd("C:/Users/54280/Documents/manuscritos/random walk/PLOSONE")


####Adehabitat####

######con adehabitat se extraen las metricas a un archivo dataframe################

# load('.RData')

packagesOfInterest <- c("rgdal", "spatial",  "adehabitatLT", "ggspatial", "magrittr", "ggpubr")
idx <- which(is.na(match(packagesOfInterest, installed.packages()[,"Package"] )))
if(length(idx)) install.packages(packagesOfInterest[idx])

idx <- which(is.na(match( paste("package:",packagesOfInterest,sep=""), search())))
if(length(idx)) sapply(idx, function(i) library(packagesOfInterest[i], character.only = TRUE))
rm(idx)


#####el archivo original es "delfinestodosgolfos" pero se limpio, y eso esta en el script analisistrayectorias########

ifile <- "dolphinstrajectories.csv"
x <- read.table(file = ifile, header = TRUE, as.is = TRUE, dec = ",", sep = ";", quote = "", fileEncoding = "latin1")
head(x)
dim(x)
str(x)

xwithloc<-x

# Converting time data ----
# For adehabitat date needs to be transformed into an object of the class POSIXct.
da <- paste(paste(xwithloc$Año, 
                  ifelse(nchar(xwithloc$Mes) ==1, paste("0", xwithloc$Mes, sep = ""), xwithloc$Mes), 
                  ifelse(nchar(xwithloc$Día) ==1, paste("0", xwithloc$Día, sep = ""),xwithloc$Día), sep = "-"), xwithloc$hora)
head(da)
# este paso agrega una colummna con fecha hora y diferencia en uso horario, todo en la misma columna, el objeto x pasa a 6948 x 17.
xwithloc$date <- as.POSIXct(strptime(da,"%Y-%m-%d %H:%M:%S"), tz = "America/Buenos_Aires")


# para cambiar el nombre de la columna IdPrimerAvistaje
idx <- match("IdPrimerAvistaje", names(xwithloc))
if(length(idx)) names(xwithloc)[idx] <- "ID"

#####################################

# Checking lon-lat ranges
format(range(xwithloc$long), digits = 7, scientific = FALSE)
format(range(xwithloc$lat), digits = 7, scientific = FALSE)

####################

# Converting into spatial objects ----
ID <- unique(xwithloc$ID)
coords <- SpatialPoints(xwithloc[c("long", "lat")], proj4string = CRS("+init=epsg:4326"))
xspdf  <-  SpatialPointsDataFrame(coords, xwithloc)

# Project data into planar units ----
# For selecting the projection
hist(xwithloc$long)
median(xwithloc$long)

xspdf <- spTransform(xspdf, CRS("+init=epsg:5346"))  #GKF4



# Converting into trajectories ----
# xtraj <- as.ltraj(xy = x[,c("long", "lat")], date = x$date, id = x$ID)
xtraj <- as.ltraj(xy = coordinates(xspdf), id = xspdf$ID, date = xspdf@data$date, infolocs = xspdf)
options(max.print = 2500)
xtraj

is.regular(xtraj) ####TRUE




#########################################################################
# Convert trajectory to dataframe ---
xtrajdf <- ld(xtraj)
head(xtrajdf)


#################
dat<-xtrajdf

###agrupo en 3 clases de comportamiento

var<-ifelse(dat$Estado=="Alimentación","Feeding",ifelse(dat$Estado=="Traslado","Traveling","Other"))
var
dat$Estado2<-var

###agrupando los meses####
var<-ifelse(dat$Mes== 12,"Summer",ifelse(dat$Mes<5, "Summer","Winter"))
var
dat$Season<-var





####################################################################################################
######modelo mixto con glmmTMB#####

library(glmmTMB)

library(lme4)
library(bbmle) ###for AICtab
library(ggplot2)###Diagnosis
library(DHARMa)##Diagnosis


####segun pdf sobre todos los metodos hay una recomendacion de usar tt-1 para el termino ar1###

tt<-as.factor(dat$Intervalo)
class(tt)

#######se corre el modelo completo con todas las interacciones dobles, termino aleatorio, y autocorrelacion temporal, se estima con REML, y funcion gamma

REMLgammafull<-glmmTMB(dist~ cos(rel.angle) * Estado2 + Estado2 * Especie + Especie * Season + 
                         Estado2 * Season + cos(rel.angle)* Especie + cos(rel.angle)* Season + Tamaño.grupo * Estado2 + Tamaño.grupo * Especie + 
                         Tamaño.grupo * Season + Tamaño.grupo * cos(rel.angle) + (1|ID)+ ar1(tt-1|ID), REML = TRUE, family = Gamma(log), data = dat)
summary(REMLgammafull)
#
#
# ###sacando el termino aleatorio para evaluar su efecto
#
REMLgamma<-glmmTMB(dist~ cos(rel.angle) * Estado2 + Estado2 * Especie + Especie * Season + 
                     Estado2 * Season + cos(rel.angle)* Especie + cos(rel.angle)* Season + 
                     Tamaño.grupo * Estado2 + Tamaño.grupo * Especie + Tamaño.grupo * Season + 
                     Tamaño.grupo * cos(rel.angle) + ar1(tt-1|ID), REML = TRUE, 
                   family = Gamma(log), data = dat)
summary(REMLgamma)
#
#
#
# ####mismos pasos considerando funcion gauss
#
REMLgaussfull<-glmmTMB(dist~ cos(rel.angle) * Estado2 + Estado2 * Especie + 
                         Especie * Season + Estado2 * Season + 
                         cos(rel.angle)* Especie + cos(rel.angle)* Season + 
                         Tamaño.grupo * Estado2 + Tamaño.grupo * Especie + 
                         Tamaño.grupo * Season + Tamaño.grupo * cos(rel.angle) + 
                         (1|ID)+ ar1(tt-1|ID), REML = TRUE, family = gaussian(), data = dat)
summary(REMLgaussfull)
#
# ###para evaluar el efecto del termino aleatorio
REMLgauss<-glmmTMB(dist~ cos(rel.angle) * Estado2 + Estado2 * Especie + 
                     Especie * Season + Estado2 * Season + cos(rel.angle)* Especie + 
                     cos(rel.angle)* Season + Tamaño.grupo * Estado2 + 
                     Tamaño.grupo * Especie + Tamaño.grupo * Season + 
                     Tamaño.grupo * cos(rel.angle) + ar1(tt-1|ID), REML = TRUE, 
                   family = gaussian(), data = dat)
summary(REMLgauss)
#
#
#

####para comparar

AICtab(REMLgaussfull, REMLgauss, REMLgammafull, REMLgamma)



#
#

#################poner a prueba interacciones de variables
####parto del beyond optimal model con distribucion gamma termino aleatorio y autocorrelacion temporal a 1 time lag, y las interacciones dobles de 5 variable: Especie, Estado, turning angle, season y tamaño de grupo, se estima por LME (no por REML)

#
beyondoptimal<-glmmTMB(dist~ cos(rel.angle) * Estado2
                        + Estado2 * Especie
                        + Especie * Season
                        + Estado2 * Season
                        + cos(rel.angle)* Especie
                        + cos(rel.angle)* Season
                        + Tamaño.grupo * Estado2
                        + Tamaño.grupo * Especie
                        + Tamaño.grupo * Season
                        + Tamaño.grupo * cos(rel.angle)
                        + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(beyondoptimal)
#
# #####diagnosis

beyondoptimal_simres <- simulateResiduals(beyondoptimal)
plot(beyondoptimal_simres)
#
ggplot(data = NULL, aes(y = resid(beyondoptimal, type = "pearson"), x = predict(beyondoptimal))) + geom_point()
#
# ###para explorar correlacion entre variables explicatorias
#
cov2cor(vcov(beyondoptimal)$cond)
#
# ##para explorar outliers
#
outliers(beyondoptimal_simres, lowerQuantile = 0, upperQuantile = 1,
          return = c("index", "logical"))
#
#
# #####################a partir del beyond optimal model saco interacciones dobles de a una
#
#
# ####Turning angle*species
#
glmmA<-glmmTMB(dist~ cos(rel.angle) * Estado2
           + Estado2 * Especie
           + Especie * Season
           + Estado2 * Season
           + cos(rel.angle)* Season
           + Tamaño.grupo * Estado2
           + Tamaño.grupo * Especie
           + Tamaño.grupo * Season
           + Tamaño.grupo * cos(rel.angle)
           + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(glmmA)
#
# ###Species*season
glmmB<-glmmTMB(dist~ cos(rel.angle) * Estado2
                        + Estado2 * Especie
                        + Estado2 * Season
                        + cos(rel.angle)* Especie
                        + cos(rel.angle)* Season
                        + Tamaño.grupo * Estado2
                        + Tamaño.grupo * Especie
                        + Tamaño.grupo * Season
                        + Tamaño.grupo * cos(rel.angle)
                        + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(glmmB)
#
# ###Estado2 * Season
glmmC<-glmmTMB(dist~ cos(rel.angle) * Estado2
           + Estado2 * Especie
           + Especie * Season
           + cos(rel.angle)* Especie
           + cos(rel.angle)* Season
           + Tamaño.grupo * Estado2
           + Tamaño.grupo * Especie
           + Tamaño.grupo * Season
           + Tamaño.grupo * cos(rel.angle)
           + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(glmmC)
#
# ####Turning angle * Estado2
#
glmmD<-glmmTMB(dist~ Estado2 * Especie
           + Especie * Season
           + Estado2 * Season
           + cos(rel.angle)* Especie
           + cos(rel.angle)* Season
           + Tamaño.grupo * Estado2
           + Tamaño.grupo * Especie
           + Tamaño.grupo * Season
           + Tamaño.grupo * cos(rel.angle)
           + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(glmmD)
#
# ####Estado2 * Especie
#
glmmE<-glmmTMB(dist~ cos(rel.angle) * Estado2
           + Especie * Season
           + Estado2 * Season
           + cos(rel.angle)* Especie
           + cos(rel.angle)* Season
           + Tamaño.grupo * Estado2
           + Tamaño.grupo * Especie
           + Tamaño.grupo * Season
           + Tamaño.grupo * cos(rel.angle)
           + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(glmmE)
#
#
# #####turning angl*season
glmmF<-glmmTMB(dist~ cos(rel.angle) * Estado2
           + Estado2 * Especie
           + Especie * Season
           + Estado2 * Season
           + cos(rel.angle)* Especie
           + Tamaño.grupo * Estado2
           + Tamaño.grupo * Especie
           + Tamaño.grupo * Season
           + Tamaño.grupo * cos(rel.angle)
           + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(glmmF)
#
# ########Turning angle*Tamaño.grupo
#
glmmG<-glmmTMB(dist~ cos(rel.angle) * Estado2
                + Estado2 * Especie
                + Especie * Season
                + Estado2 * Season
                + cos(rel.angle)* Especie
                + cos(rel.angle)* Season
                + Tamaño.grupo * Estado2
                + Tamaño.grupo * Especie
                + Tamaño.grupo * Season
                + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(glmmG)
#
# #####Especie* Tamaño.grupo
#
glmmH<-glmmTMB(dist~ cos(rel.angle) * Estado2
           + Estado2 * Especie
           + Especie * Season
           + Estado2 * Season
           + cos(rel.angle)* Especie
           + cos(rel.angle)* Season
           + Tamaño.grupo * Estado2
           + Tamaño.grupo * Season
           + Tamaño.grupo * cos(rel.angle)
           + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(glmmH)
#
# #######Season*Tamaño.grupo
glmmI<-glmmTMB(dist~ cos(rel.angle) * Estado2
           + Estado2 * Especie
           + Especie * Season
           + Estado2 * Season
           + cos(rel.angle)* Especie
           + cos(rel.angle)* Season
           + Tamaño.grupo * Estado2
           + Tamaño.grupo * Especie
           + Tamaño.grupo * cos(rel.angle)
           + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(glmmI)
#
#
# ####Estado2*Tamaño.grupo
glmmJ<-glmmTMB(dist~ cos(rel.angle) * Estado2
           + Estado2 * Especie
           + Especie * Season
           + Estado2 * Season
           + cos(rel.angle)* Especie
           + cos(rel.angle)* Season
           + Tamaño.grupo * Especie
           + Tamaño.grupo * Season
           + Tamaño.grupo * cos(rel.angle)
           + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(glmmJ)
#
# #
#
#para armar tabla 3 del ms
AIC(beyondoptimal,glmmA, glmmB, glmmC, glmmD, glmmE, glmmF, glmmG, glmmH, glmmI, glmmJ)


#
#
# ####modelo nuevo sacando ineracciones no significativas #new beyond optimal model

beyondoptimal2<-glmmTMB(dist~  Estado2 * Especie
                         + cos(rel.angle)* Especie
                         + cos(rel.angle)* Season
                         + Tamaño.grupo * Estado2
                         + Tamaño.grupo * Especie
                         + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
 summary(beyondoptimal2)
#
#
# #####Tamaño.grupo * Especie
beyondoptimal2A<-glmmTMB(dist~  Estado2 * Especie
                         + cos(rel.angle)* Especie
                         + cos(rel.angle)* Season
                         + Tamaño.grupo * Estado2
                         + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(beyondoptimal2A)
#
# ######Tamaño.grupo * Estado2
#
beyondoptimal2B<-glmmTMB(dist~  Estado2 * Especie
                         + cos(rel.angle)* Especie
                         + cos(rel.angle)* Season
                         + Tamaño.grupo * Especie
                         + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(beyondoptimal2B)
#
# ####+ cos(rel.angle)* Season
beyondoptimal2C<-glmmTMB(dist~  Estado2 * Especie
                         + cos(rel.angle)* Especie
                         + Tamaño.grupo * Estado2
                         + Tamaño.grupo * Especie
                         + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(beyondoptimal2C)

# ##### + cos(rel.angle)* Especie
beyondoptimal2D<-glmmTMB(dist~  Estado2 * Especie
                         + cos(rel.angle)* Season
                         + Tamaño.grupo * Estado2
                         + Tamaño.grupo * Especie
                         + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(beyondoptimal2D)
#
# ######Estado2 * Especie
beyondoptimal2E<-glmmTMB(dist~  cos(rel.angle)* Especie
                         + cos(rel.angle)* Season
                         + Tamaño.grupo * Estado2
                         + Tamaño.grupo * Especie
                         + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(beyondoptimal2E)
#
AIC(beyondoptimal2,beyondoptimal2A, beyondoptimal2B, beyondoptimal2C, beyondoptimal2D, beyondoptimal2E)

####Para tabla 3######
AICtab(beyondoptimal2,beyondoptimal2A, beyondoptimal2B, beyondoptimal2C, beyondoptimal2D, beyondoptimal2E)
anova(beyondoptimal2,beyondoptimal2D)

########la intraccion turning angle*Species es marginalmente significativa

###FINAL#######
####modelo final####
###Este es el modelo que queda despues de sacar las interacciones dobles no significativas

Finalselectedmodel<-glmmTMB(dist~  Estado2 * Especie
                         + cos(rel.angle)* Season
                         + Tamaño.grupo * Estado2
                         + Tamaño.grupo * Especie
                         + (1|ID)+ ar1(tt-1|ID), family = Gamma(log), data = dat)
summary(Finalselectedmodel)

#####diagnosis del modelo final

Finalselectedmodel_simres <- simulateResiduals(Finalselectedmodel)
plot(Finalselectedmodel_simres)
ggplot(data = NULL, aes(y = resid(Finalselectedmodel, type = "pearson"), x = predict(Finalselectedmodel))) + geom_point()


residuals(Finalselectedmodel)

toto <- recalculateResiduals(Finalselectedmodel_simres, group = tt, aggregateBy = mean)
#####Error in aggregate.data.frame(as.data.frame(x), ...) : 
  ###arguments must have same length


acf(toto$scaledResiduals)
pacf(toto$scaledResiduals)


plot(residuals(Finalselectedmodel), Finalselectedmodel_simres$fittedResiduals)

# testTemporalAutocorrelation(toto$scaledResiduals, levels(toto$group))

# Step length
par(mfrow=c(1,1))
# ylimit <- qnorm((1 + 0.95)/2)/sqrt(length(pr$stepRes[!is.na(pr$stepRes)])) + 1

acf(Finalselectedmodel_simres$fittedResiduals,lag.max = 50)
pacf(Finalselectedmodel_simres$fittedResiduals,lag.max = 50)

# name.plot <- paste0("./Figures/", paste0("F_mod_", i, "_acf_step.png"))
# dev.copy(png, name.plot,  width = 800, height = 500)
# dev.off()


ggplot(data = NULL, aes(y = resid(Finalselectedmodel, type = "pearson"), x = predict(Finalselectedmodel))) + geom_point()

                           #could not find function "%>%"
                                                                              # > 
######el mismo modelo sin correlacion temporal####

finalselectedmodel0<-glmmTMB(dist~  Estado2 * Especie
                         + cos(rel.angle)* Season
                         + Tamaño.grupo * Estado2
                         + Tamaño.grupo * Especie
                         + (1|ID), family = Gamma(log), data = dat)
summary(finalselectedmodel0)

AIC(Finalselectedmodel, finalselectedmodel0)

#####diagnosis del Modelo sin correlacion temporal

finalselectedmodel0_simres <- simulateResiduals(finalselectedmodel0)
plot(finalselectedmodel0_simres)

ggplot(data = NULL, aes(y = resid(finalselectedmodel0, type = "pearson"), x = predict(finalselectedmodel0))) + geom_point()

toto0 <- recalculateResiduals(finalselectedmodel0_simres, group = tt, aggregateBy = mean)

acf(toto0$scaledResiduals)
pacf(toto0$scaledResiduals)

plot(levels(tt), toto0$scaledResiduals, type = "b")


acf(finalselectedmodel0_simres$fittedResiduals,lag.max = 50)
pacf(finalselectedmodel0_simres$fittedResiduals,lag.max = 50)

#########profile likelihood ratio confidence intervals of the parameters of the final select model to be added to S2 Table 
####intervalos de confianza en S2 Table Supp Inf
packageVersion("glmmTMB")
packageVersion("TMB")

toto <- glmmTMB:::confint.glmmTMB(Finalselectedmodel, method = "profile", parallel = "snow", ncpus = 4)
# write.table(toto, file="Lprofile")
# toto %>% tibble::as_tibble() %>% mutate(param = rownames(toto)) %>% write.csv(.,file = 'toto.csv', row.names = F)
#Error in toto %>% as_tibble() %>% mutate(param = rownames(toto)) %>% write.csv(.,  : 

##################Figures in the ms


####Figure 2

##################
library(ggplot2)
library(magrittr)
library(ggspatial)
library(dplyr)

delfin1 <- dat[which(dat$ID == 10013),]

# asignar estado de comportamiento a cada vector 
## establecer los segmentos
delfin1_secuencia <- as.data.frame(with(delfin1,cbind(embed(y,2),embed(x,2))))
colnames(delfin1_secuencia) <- c('yend','y','xend','x')

## extraer los comportamientos que seran asignados a cada vector
# delfin1_estados <- delfil
#   delfin1 %>% 
#   dplyr::select(x, Estado2)  
  # remover la ultima observacion porque no se corresponde con ningun vector
  # el ultimo punto en el grafico se plotea de todas maneras como un triangulo

delfin1_estados <- delfin1[1:nrow(delfin1)-1,c("x", "Estado2")]

## unir la secuencia de vectores con el estado comportamental. Se unen con la primera posicion x, como hace adehabitat
delfin1_secuencia <- merge(x = delfin1_estados, y = delfin1_secuencia, by = 'x')
delfin1_secuencia <- delfin1_secuencia %>% 
  transform(Estado2 = factor(Estado2,
                             levels = c(
                               'Feeding',  'Traveling',  'Other'
                             )))

commonplot <- ggplot(data = delfin1_secuencia) +
  geom_segment(aes(x = x,y = y,xend = xend,yend = yend,colour = Estado2)) +
  geom_point( aes(x = x , y = y), col = 'gray70', cex = 0.85) +
  # scale_color_viridis_d(option = "C", end = 0.7) +
  geom_point(x = delfin1$x[1], y = delfin1$y[1], cex= 2, shape = 22, fill = "black") +
  geom_point(x = delfin1$x[length(delfin1$x)], y = delfin1$y[length(delfin1$x)], cex= 2, shape = 24, fill = "black") +
  labs(color = "Behavior") +
  theme_bw() +
  annotation_scale()  +
  theme(legend.title = element_blank())


delfin2 <- dat[which(dat$ID == 271),]

# asignar estado de comportamiento a cada vector 
## establecer los segmentos
delfin2_secuencia <- as.data.frame(with(delfin2,cbind(embed(y,2),embed(x,2))))
colnames(delfin2_secuencia) <- c('yend','y','xend','x')

## extraer los comportamientos que seran asignados a cada vector

delfin2_estados <- delfin2[1:nrow(delfin2)-1,c("x", "Estado2")]
# delfin2_estados <- 
#   delfin2 %>% 
#   select(x, Estado2) %>% 
#   # remover la ultima observacion porque no se corresponde con ningun vector
#   # el ultimo punto en el grafico se plotea de todas maneras como un triangulo
#   slice(., 1:(n() - 1)) 

## unir la secuencia de vectores con el estado comportamental. Se unen con la primera posicion x, como hace adehabitat
delfin2_secuencia <- merge(x = delfin2_estados, y = delfin2_secuencia, by = 'x')
delfin2_secuencia <- delfin2_secuencia %>% 
  transform(Estado2 = factor(Estado2,
                             levels = c(
                               'Feeding',  'Traveling',  'Other'
                             )))

duskyplot <- ggplot(data = delfin2_secuencia) +
  geom_segment(aes(x = x,y = y,xend = xend,yend = yend,colour = Estado2)) +
  geom_point( aes(x = x , y = y), col = 'gray70', cex = 0.85) +
  # scale_color_viridis_d(option = "C", end = 0.7) +
  geom_point(x = delfin2$x[1], y = delfin2$y[1], cex = 2, shape = 22, fill = "black") +
  geom_point(x = delfin2$x[length(delfin2$x)], y = delfin2$y[length(delfin2$x)], cex = 2, shape = 24, fill = "black") +
  labs(color = "Behavior") +
  theme_bw() +
  annotation_scale()  +
  theme(legend.title = element_blank())

# library(patchwork)
# 
# patchwork::wrap_plots(duskyplot, commonplot, ncol = 2)
# cowplot::plot_grid(duskyplot, commonplot, ncol = 2)

p.dolphin.tracks <- ggpubr::ggarrange(duskyplot,
          commonplot,
          labels = NULL, ncol = 2, common.legend = T,legend = "bottom" )

# ggsave(plot = p.dolphin.tracks, filename = 'Fig_dolphin_tracks.png', height = 8, width = 13)

x <- bind_rows(
  delfin1[c(1, nrow(delfin1)), c('date', 'Season', 'Especie', 'ID')],
  delfin2[c(1, nrow(delfin2)), c('date', 'Season', 'Especie', 'ID')])

# write.csv(x,
#   file = 'Fig_dolphin_tracks_Start_End.csv', row.names = FALSE)

##################
###################Figure 3 and 4
dat$Especie <- factor(dat$Especie, levels = c("Lagenorhynchus obscurus", "Delphinus delphis"),
                      labels = c("Dusky dolphin", "Common dolphin"))

############Figure 3

ggplot(data = dat, aes(x = dist)) +
  geom_histogram(aes(y = stat(density)),
                 colour = 1, fill = "white") +
  geom_density(lwd = 1, colour = 4,
               alpha = 0.25) +
  xlab("Step length") +
  facet_wrap(Especie ~Estado2) +
  theme_bw()
# ggsave("StepLength_Behavior_Especie.png")

###Figure 4

ggplot(data = dat, aes(x = rel.angle)) +
  geom_histogram(aes(y = stat(density)),
                 colour = 1, fill = "white") +
  geom_density(lwd = 1, colour = 4,
               alpha = 0.25) +
  xlab("Turning angle") +
  facet_wrap(Especie ~Estado2) +
  theme_bw()
# ggsave("TurnAngle_Behavior_Especie.png")

##################


##########figure 5 a 8 separando por especie######
dat<-xtrajdf
###agrupando comportamientos pero sin agrupar meses

var<-ifelse(dat$Estado=="Alimentación","Feeding",ifelse(dat$Estado=="Traslado","Traveling","Other"))
var
dat$Estado2<-var

####subset para cada especie
Lo <- subset(dat, dat$Especie=="Lagenorhynchus obscurus")
dim(Lo)
idx <- which(is.na(Lo$dist))
Lodistdata<- Lo[-idx,]
dim(Lodistdata)
max((Lodistdata$dist))
mean((Lodistdata$dist))
sd((Lodistdata$dist))
median(Lodistdata$dist)

max(Lodistdata$Tamaño.grupo)

idx<- which(is.na(Lo$rel.angle))
Lorelangledata<-Lo[-idx,]
dim(Lorelangledata)
Lotraj<-dl(Lo)
Lodisttraj<-dl(Lodistdata)
Loangtraj<-dl(Lorelangledata)

Dd <- subset (dat, dat$Especie=="Delphinus delphis")
dim(Dd)
idx <- which(is.na(Dd$dist))
Dddistdata<- Dd[-idx,]
dim(Dddistdata)
max(Dddistdata$dist)
mean(Dddistdata$dist)
sd(Dddistdata$dist)
median(Dddistdata$dist)

max(Dddistdata$Tamaño.grupo)
min(Dddistdata$Tamaño.grupo)

idx<- which(is.na(Dd$rel.angle))
Ddrelangledata<-Dd[-idx,]
dim(Ddrelangledata)
Ddtraj<-dl(Dd)
Dddisttraj<-dl(Dddistdata)
Ddangtraj<-dl(Ddrelangledata)


#################Figure 5
par(mfrow=c(1,2))

Lox<- Lo
Lox$Tamaño.grupo<-factor(Lox$Tamaño.grupo, levels=c("<10", "10 a 20", "20 a 50", "50 a 70", "70 a 100", ">100"))
boxplot(dist~Lox$Tamaño.grupo, data=Lox, main="Dusky dolphin", xlab="Group size", ylab="Step length (m)")

Ddx<- Dd
Ddx$Tamaño.grupo<-factor(Ddx$Tamaño.grupo, levels=c("<10", "10 a 20", "20 a 50", "50 a 70", "70 a 100", ">100"))
boxplot(dist~Ddx$Tamaño.grupo, data=Ddx, main="Common dolphin", xlab="Group size", ylab="Step length (m)")

#dist~behavior en cada especie######Figure 6
boxplot(dist~Lo$Estado2, data=Lo, main="Dusky dolphin", xlab="Behavioral state", ylab="Step length (m)")
boxplot(dist~Dd$Estado2, data=Dd, main="Common dolphin", xlab="Behavioral state", ylab="Step length (m)")



###Turning angle en funcion de behavior 
###para discriminar por especie cada grafico#####Figure 7

travelingLo<- subset(Lo, Lo$Estado2=="Traveling")
feedingLo<-subset(Lo, Lo$Estado2=="Feeding")
otherLo<-subset(Lo, Lo$Estado2=="Other")

travelingDd<- subset(Dd, Dd$Estado2=="Traveling")
feedingDd<-subset(Dd, Dd$Estado2=="Feeding")
otherDd<-subset(Dd, Dd$Estado2=="Other")

par(mfrow=c(1,3))
plot(dist~travelingLo$rel.angle, data=travelingLo, main="Traveling", xlab="Turning angle (rad)", ylab="Step length (m)", col="grey")
points(dist~travelingDd$rel.angle, data=travelingDd, main="Traveling", xlab="Turning angle (rad)", ylab="Step length (m)", col="black")

plot(dist~feedingLo$rel.angle, data=feedingLo, main="Feeding", xlab="Turning angle (rad)", ylab="Step length (m)", col="grey")
points(dist~feedingDd$rel.angle, data=feedingDd, main="Traveling", xlab="Turning angle (rad)", ylab="Step length (m)", col="black")

plot(dist~otherLo$rel.angle, data=otherLo, main="Other", xlab="Turning angle (rad)", ylab="Step length (m)", col="grey")
points(dist~otherDd$rel.angle, data=otherDd, main="Traveling", xlab="Turning angle (rad)", ylab="Step length (m)", col="black")




#Dist~Mes en cada especie############Figure 8
par(mfrow=c(1,2))
boxplot(dist~ Lo$Mes, data=Lo, main="Dusky dolphin", xlab= "Month",ylab="Step length (m)")
boxplot(dist~Dd$Mes, data=Dd, main="Common dolphin", xlab= "Month",ylab="Step length (m)")

##################################################################
#####para graficar predichos en funcion d e interacciones####
####Figure 9
#Creo el set de datos para mostrar el modelo ajustado 
####los observados son los datos con lo que se corrio el modelo agrupando comportamientos, meses y redefiniendo el intervalo como tt
dat<-xtrajdf

###agrupo comportamientos
var<-ifelse(dat$Estado=="Alimentación","Feeding",ifelse(dat$Estado=="Traslado","Traveling","Other"))
var
dat$Estado2<-var

##agrupo meses
var<-ifelse(dat$Mes== 12,"Summer",ifelse(dat$Mes<5, "Summer","Winter"))
var
dat$Season<-var

#redefino tt
tt<-as.factor(dat$Intervalo)
class(tt)


####primero creo un grid de todas las combinaciones entre variables que quiero mostrar

####para armar grilla de prediccion###


newData1<-expand.grid(rel.angle = rep(seq(from=-3.14, to= 3.14, 0.19625),1),
                      Especie = rep (c( "Lagenorhynchus obscurus","Delphinus delphis"), 1), 
                      Season = rep (c( "Summer", "Winter"), 1),
                      Estado2 = rep (c("Feeding", "Traveling", "Other"), 1),
                      Tamaño.grupo = rep (c("<10", "10 a 20", "20 a 50", "50 a 70", "70 a 100", ">100"), 1),
                      ID= rep((NA), 1),
                      tt = rep((NA),1 ))

newData1
dim(newData1)
head(newData1)



predict<-predict(Finalselectedmodel, newdata= newData1, allow.new.levels = TRUE, type="response", re.from=NULL)

########Grafico final con predicciones del modelo#######



p<- ggplot( newData1, aes(x = rel.angle, y =predict, colour= Tamaño.grupo))+ 
  scale_color_manual(values=c('grey', 'orange', 'red', 'violet','blue', 'black'))+
  geom_line()+
  scale_x_continuous(breaks = c(-3,0, 3))+ #Si querés que el valor sea -2, 0, 2 entonces borrá esta lína
  theme(legend.position = "bottom")
p




#######################Figura 9 en el ms#######

p<-p+ facet_grid(Especie~Season+Estado2 , scales = "fixed", margins=FALSE)

p

######etiquetas y ejes

p<- p+labs(x = "Turning angle (rad)", y = "Step length (m)") 
p + ggtitle("" ) +
  scale_color_manual(values=c('grey', 'orange', 'red', 'violet','blue', 'black'), 
                     name = "Group size", labels = c('<10','10-20','20-50','50-70','70-100','>100'))

delfines <- c(`Lagenorhynchus obscurus` = 'Dusky dolphin',
              `Delphinus delphis` = 'Common dolphin'
)

p<-p+ facet_grid(Especie~Season+Estado2 , scales = "fixed", margins=FALSE, labeller = labeller(Especie = delfines)) +
  labs(color= 'Group Size')

p

# ggsave("Figure 9.png")


