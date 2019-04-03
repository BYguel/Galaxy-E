#!/usr/bin/env Rscript



##################################################################################################################################
############## FUNCTION TO CALCULATE AND PLOT CSI CTI AND CTRI OF communities  function:csi_cti_ctri    ##############################
##################################################################################################################################

#### Based on Romain Lorrillière R script
#### Modified by Alan Amosse and Benjamin Yguel for integrating within Galaxy-E

### made with R version 3.5.1

suppressMessages(library(RODBC))  ##Version: 1.3-15
suppressMessages(library(reshape))  ##Version: 0.8.8
suppressMessages(library(data.table))  ##Version: 1.12.0
suppressMessages(library(rgdal))  ##Version: 1.3-4
suppressMessages(library(lubridate))  ##Version: 1.7.4
suppressMessages(library(RPostgreSQL))  ##Version: 0.6-2
suppressMessages(library(doBy))  ##Version: 4.6-2
suppressMessages(library(arm))  ##Version: 1.10-1
suppressMessages(library(ggplot2))  ##Version: 3.1.0
suppressMessages(library(scales))  ##Version: 1.0.0
suppressMessages(library(mgcv))  ##Version: 1.8-24
suppressMessages(library(visreg))  ##Version: 2.5-0
suppressMessages(library(plyr))  ##Version: 1.8.4
suppressMessages(library(lme4))  ##Version: 1.1-18-1
suppressMessages(library(lmerTest))  ##Version: 3.1-0




###########
#delcaration des arguments et variables/ declaring some variables and load arguments

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
    stop("At least 5 arguments must be supplied :\n- An input dataset filtered (.tabular). May come from the filter rare species tool.\n- A species detail table (.tabular)\n- A species ssi/sti table.\n- A table with plots coordinates.\n- table with csi calculated before 2001.\n\n", call.=FALSE) #si pas d'arguments -> affiche erreur et quitte / if no args -> error and exit1
} else {
    Datafiltered<-args[1] ###### Nom du fichier avec extension ".typedefichier", peut provenir de la fonction "FiltreEspeceRare" / file name without the file type ".filetype", may result from the function "FiltreEspeceRare"    
    tabSpecies<-args[2] ###### Nom du fichier avec extension ".typedefichier", fichier mis à disposition dans Galaxy-E avec specialisation à l'habitat des especes et si espece considérée comme indicatrice / file name without the file type ".filetype", file available in Galaxy-E containing habitat specialization for each species and whether or not they are considered as indicator  
    tabssi<-args[3] ##### Nom du fichier avec extension ".typedefichier", fichier mis à disposition dans Galaxy-E avec degre de specialisation de l espece et affinite thermique /file name without the file type ".filetype", file available in Galaxy-E containing specilalization degree as well as thermic preferences
    coordCarre<-[4] #### Nom du fichier avec extension ".typedefichier", fichier mis à disposition dans Galaxy-E avec les coordonnees gps des carres /file name without the file type ".filetype", file available in Galaxy-E containing gps coordinates of the plots
	csibefore2001<-[5] #### Nom du fichier avec extension ".typedefichier", fichier mis à disposition dans Galaxy-E avec les calcul de csi avant 2001 sur base de données différentes / file name without the file type ".filetype", file available in Galaxy-E containing csi calculated for years before 2001
}


#Import des données / Import data 
tabCLEAN <- read.table(Datafiltered,sep="\t",dec=".",header=TRUE) #### charge le fichier de données d abondance / load abundance of species
tabsp <- read.table(tabSpecies,sep="\t",dec=".",header=TRUE)   #### charge le fichier de donnees sur nom latin, vernaculaire et abbreviation, espece indicatrice ou non / load the file with information on species specialization and if species are indicators

spTrait=read.csv2("tabssi") ############# species_indicateur_fonctionnel.csv  doit posseder au moins une colonne ssi, une colonne thermal_niche_mean  /
coordCarre=read.csv2("coordCarre") ######## carre.csv  charge les coordonnées des carrés qui sont utilisés comme covariable  / 
csibefore2001=read.csv2("csibefore2001")##### csi_init.csv  

######### Calcul du csi (indice de specialisation à l'habitat) / Calculation of the csi (index of habitat specialization)
tabCLEAN$ssi=spTrait$ssi[match(tabCLEAN$espece,spTrait$pk_species)] ### recupere donnee de ssi par espece calcule a partir des cv entre habitat obtenu a partir des occurences par carre et descriptif des habitats de chaque carre / retrieve ssi data for each species based on their variation coeffficent between habitat described for each plots
tabCLEAN$longitude_grid_wgs84=coordCarre$longitude_grid_wgs84[match(tabCLEAN$carre,coordCarre$pk_carre)] #### recupere les coordonnees carres
tabCLEAN$latitude_grid_wgs84=coordCarre$latitude_grid_wgs84[match(tabCLEAN$carre,coordCarre$pk_carre)]#### recupere les coordonnees carres

ssicarre=aggregate(ssi*abond~annee+carre,tabCLEAN,sum) ### ssi totale par annee et carre pondere par les abondances / ssi total per year and per plots weighted by abundances
abcarre=aggregate(abond~annee+carre,tabCLEAN,sum) ### somme des abondances totales par annee et carre / sum of total abundance per year and plots
csi=ssicarre[,3]/abcarre[,3]  #### le ssi moyen par carre = csi par carre et annee / mean ssi per plots = csi per year and plots

dd1=data.frame(csi,ssicarre$carre,ssicarre$annee) 
names(dd1)[2]="carre"
names(dd1)[3]="year"
######### Calcul du cti (indice de preference thermique) / Calculation of the cti (index of thermal preference)
tabCLEAN$sti=spTrait$sti[match(tabCLEAN$espece,spTrait$pk_species)] ### recupere donnee de sti par espece calcule a partir temperature moyenne obtenu a partir des occurences par carre et descriptif des habitats de chaque carre / retrieve sti data for each species based on the temperature mean between habitat described for each plots

sticarre=aggregate(sti*abond~annee+carre,tabCLEAN,sum) ### sti totale par annee et carre pondere par les abondances / ssi total per year and per plots weighted by abundances
abcarre=aggregate(abond~annee+carre,tabCLEAN,sum) ### somme des abondances totales par annee et carre / sum of total abundance per year and plots
cti=sticarre[,3]/abcarre[,3]  #### le sti moyen par carre = cti par carre et annee / mean sti per plots = cti per year and plots

dd2=data.frame(cti,sticarre$carre,sticarre$annee)  #############
names(dd2)[2]="carre"
names(dd2)[3]="year"
dd3=merge(dd1,dd2,by=c("carre","year")) ################ 

######### Calcul du ctri (indice de niveau trophique) avec utilisation de l'exponentiel du niveau trophique / Calculation of the ctri (index of trophic level) using the exponential of the trophic level (1=herbivorous birds, 2=insectivorous/nematophagous birds, 3=birds predator)
tabCLEAN$stri=spTrait$exp_stri[match(tabCLEAN$espece,spTrait$pk_species)] ### recupere donnee de stri par espece calcule a partir temperature moyenne obtenu a partir des occurences par carre et descriptif des habitats de chaque carre / retrieve sti data for each species based on the temperature mean between habitat described for each plots

stricarre=aggregate(stri*abond~annee+carre,tabCLEAN,sum) ### sti totale par annee et carre pondere par les abondances / ssi total per year and per plots weighted by abundances
abcarre=aggregate(abond~annee+carre,tabCLEAN,sum) ### somme des abondances totales par annee et carre / sum of total abundance per year and plots
ctri=stricarre[,3]/abcarre[,3]  #### le stri moyen par carre = ctri par carre et annee / mean stri per plots = ctri per year and plots

dd4=data.frame(ctri,stricarre$carre,stricarre$annee)  #############  
names(dd4)[2]="carre"
names(dd4)[3]="year"
dd=merge(dd3,dd4,by=c("carre","year")) ################ 


dd$longitude_grid_wgs84=coordCarre$longitude_grid_wgs84[match(dd$carre,coordCarre$pk_carre)] #### recupere coordonnées gps / retrieve gps coordinates
dd$latitude_grid_wgs84=coordCarre$latitude_grid_wgs84[match(dd$carre,coordCarre$pk_carre)]  #### recupere coordonnées gps / retrieve gps coordinates
dd$id_plot=dd$carre  ### id_plot nom données aux carrés dans le script /id_plot is use as the name of the plot in the following script



############################# The function


csi_cti_ctri <- function(dd=dd,ic=TRUE, indic="csi",methode="gam", ####### indic: choisir indicateur "csi" "cti" ou "ctri" ; methode: choisir modele "gam" ou "lmer" ; ic pour calcul de l'interval de confiance plus rapide sans mais moins fiable / indic is the choice of an indicator: "csi", "cti" or "ctri" ; methode is the statistical model use for the analysis "gam" or "lmer" ; ic is for the calculation of confidence interval faster without but less reliable
                          firstYear=NULL,lastYear=NULL,altitude=800,departement=NULL,onf=TRUE,distance_contact=NULL, #### altitude, departement onf, distance de contact = Argument non utilise, se trouvait dans requete sql / altitude, departement onf, distance de contact = not use anymore was in a postgres request
                         spExcluPassage1=c("MOTFLA","SAXRUB","ANTPRA","OENOEN","PHYTRO"),# (Prince et al. 2013 Env. Sc. and Pol.) + "OENOEN","PHYTRO" avis d'expert F. Jiguet, #### Argument non utilise, se trouvait dans requete sql / not use anymore was in a postgres request
                         seuilAbondance=.99,init_1989 = FALSE,plot_smooth=TRUE, ###### init_1989 si TRUE, option que pour csi et besoin du fichier des csi calculés sur les données avant 2001 (pas forcement fiable car protocole un peu different) / init_1989 if TRUE, only working for csi, and use calculation of csi based on data before 2001 (protocol was bit different, not totally reliable)
                          champSp = "code_sp", sp=NULL,champsHabitat=FALSE, #### Argument non utilise, se trouvait dans requete sql / not use anymore was in a postgres request
                          anglais=FALSE,seuilSignif=0.05,##### #### anglais=FALSE Argument non utilise, se trouvait dans requete sql / not use anymore was in a postgres request
                          couleur="#4444c3",
                          titreY=indic,titreX="Années",titre=indic,
                          savePostgres=FALSE,output=FALSE,   ##### OPTION "output" pour afficher le resultat dans R  / OPTION "output" is only to show the result in the R window 
                          operateur=c("Lorrilliere Romain","lorrilliere@mnhn.fr"), encodingSave="ISO-8859-1",fileName="dataCSI",id="France"){ ####### nom des fichiers de sorties et de l'operateur / name of the output files and of the operator

start <- Sys.time()
            dateExport <- format(start,"%Y-%m-%d")

if(is.null(firstYear)) firstYear <- 2001
if(is.null(lastYear)) lastYear <- 9999
if(is.null(altitude)) altitude <- 10000




 annee <- sort(unique(dd$year))
            nban <- length(annee)
            pasdetemps <-nban-1
			
			
if(methode == "gam") {

            cat("Methode: gam\n")
			
			
####browser()
            ## Utilisation des modèles GAMM pour obtenir les tendances d evolution par an du csi cti ou ctri !!!! Marche pas si peu de données !!!!  / Use of GAMM model for the estimation of the annual variations of the csi cti or ctri  !!! does not work with few data !!! 
            cat("\nEstimation de la variation annuelle ",indic,"~ factor(year)+s(longitude_grid_wgs84,latitude_grid_wgs84,bs='sos')\n",sep="")
if (indic=="csi"){
        gammf <- gamm(csi~ factor(year)+s(longitude_grid_wgs84,latitude_grid_wgs84,bs="sos"), data=dd,random=reStruct(object = ~ 1| id_plot, pdClass="pdDiag"),correlation=corAR1(form=~year)) #### spline sur les coordonnées, effet aleatoire sur les carres, methode autoregressive sur l'année N-1  / spline on the gps coordinates, random effect on the plots, autoregressive method on the year-1 
}
if(indic=="cti"){
        gammf <- gamm(cti~ factor(year)+s(longitude_grid_wgs84,latitude_grid_wgs84,bs="sos"), data=dd,random=reStruct(object = ~ 1| id_plot, pdClass="pdDiag"),correlation=corAR1(form=~year)) #### 

}
if(indic=="ctri"){
        gammf <- gamm(ctri~ factor(year)+s(longitude_grid_wgs84,latitude_grid_wgs84,bs="sos"), data=dd,random=reStruct(object = ~ 1| id_plot, pdClass="pdDiag"),correlation=corAR1(form=~year)) #### 

}

            sgammf<-summary(gammf$gam)
            coefdata=coefficients(gammf$gam) ### recupere les coefficient de regression de la variable "annee" / retrieve the regression coefficient of the variable "year"
            coefannee <- c(0,sgammf$p.coeff[2:nban])  ### meme chose que au dessus / same as before
            erreuran <- c(0,sgammf$se[2:nban])### recupere les erreurs standard des coefficient de regression de la variable "annee" / retrieve the standard errors of the regression coefficient of the variable "year"
            pval <-  c(1,sgammf$p.pv[2:nban])### recupere les p value de la variable "annee" / retrieve the p value of the variable year

            
			## calcul des intervalles de confiance  / confidence interval calculation
			
        if(ic) {
            # gammf.sim <- sim(gammf)  ######################  VERSION ROMAIN mais fct sim() ne marche pas avec GAMM / old version using function sim() but did not work with Gamm models
            # ic_inf_sim <- c(0,tail(apply(coef(gammf.sim), 2, quantile,.025),pasdetemps))
            # ic_sup_sim <- c(0,tail(apply(coef(gammf.sim), 2, quantile,.975),pasdetemps))
			icalpha05 <- as.data.frame(confint(gammf$gam))[2:nban,1:2]  ########## VERSION BENJ 
			ic_inf_sim <- icalpha05[,1]
			ic_inf_sim <- c("NA",ic_inf_sim[1:nban-1])
			ic_sup_sim <- icalpha05[,2]
			ic_sup_sim <- c("NA",ic_sup_sim[1:nban-1])
        } else
        {
            ic_inf <- "not assessed"
            ic_sup <- "not assessed"
        }

       

            tabfgamm <- data.frame(model = "gamm factor(year) plot",annee,coef=coefannee,se = erreuran,pval,signif=pval<seuilSignif,Lower_ci=ic_inf_sim,upper_ci=ic_sup_sim) #### recupère les resultats des modèles avec interval de confiance / retrieve results of the models used with confidence interval 






if(plot_smooth) {   #### Representation graphique de l'evolution annuelle des indicateurs  / Graphical representation of the annual evolution of the indicators
            cat("\nGam pour la figure ",indic,"~s(year)\n",sep="")
            ## create a sequence of temperature that spans your temperature  #####not use anymore 
            ## http://zevross.com/blog/2014/09/15/recreate-the-gam-partial-regression-smooth-plots-from-r-package-mgcv-with-a-little-style/  #### method for the plot

if (indic=="csi"){
            gammgg <- gamm(csi~s(year), data=dd,random=reStruct(object = ~ 1| id_plot, pdClass="pdDiag"),correlation=corAR1(form=~year))  #### spline sur l'année, effet aleatoire des carres sur ordonnée à l'origine, methode autoregressive sur l'année N-1  / spline on the year, random effect of the plots on the intercept, autoregressive method on the year-1 
}
if (indic=="cti"){
            gammgg <- gamm(cti~s(year), data=dd,random=reStruct(object = ~ 1| id_plot, pdClass="pdDiag"),correlation=corAR1(form=~year))  #### 
}
if (indic=="ctri"){
            gammgg <- gamm(ctri~s(year), data=dd,random=reStruct(object = ~ 1| id_plot, pdClass="pdDiag"),correlation=corAR1(form=~year))  ####  
}

            maxyear<-max(dd$year)
            minyear<-min(dd$year)
            year.seq<-sort(unique(c(minyear:maxyear,(seq(minyear, maxyear,length=1000)))))
            year.seq<-data.frame(year=year.seq)

                                        # predict only the temperature term (the sum of the   ########### ???? not use anymore
                                        # term predictions and the intercept gives you the overall########### ???? not use anymore
                                        # prediction)########### ???? not use anymore

            preds<-predict(gammgg$gam, newdata=year.seq, type="terms", se.fit=TRUE)


                                        # set up the temperature, the fit and the upper and lower########### ???? not use anymore
                                        # confidence interval########### ???? not use anymore

            year <-year.seq$year
            realYear <- sort(unique(dd$year))
            fit<-as.vector(preds$fit)
            init <- fit[1]

            # fit.up95<-fit-1.96*as.vector(preds$se.fit)    ###########  not use anymore

            # fit.low95<-fit+1.96*as.vector(preds$se.fit)

           # ggGamData <- data.frame(year=year, csi=fit,ic_low95 = fit.low95, ic_up95 = fit.up95)

        fit <- fit - init
        fit.up95 <- fit.up95 - init
        fit.low95 <- fit.low95 - init
		
		
if (indic=="csi"){
        ggGamData <- data.frame(year=year, csi=fit,ic_low95 = fit.low95, ic_up95 = fit.up95)  ####### Recupère les resultats des modèles / retrieve the results of the models 
}
if (indic=="cti"){
        ggGamData <- data.frame(year=year, cti=fit,ic_low95 = fit.low95, ic_up95 = fit.up95)
}       
if (indic=="ctri"){
        ggGamData <- data.frame(year=year, ctri=fit,ic_low95 = fit.low95, ic_up95 = fit.up95)
} 
 ## The ggplot:
if (indic=="csi"){
            gg <- ggplot(data=ggGamData,aes(x=year,y=csi))
           }
if (indic=="cti"){
            gg <- ggplot(data=ggGamData,aes(x=year,y=cti))
           }
if (indic=="ctri"){
            gg <- ggplot(data=ggGamData,aes(x=year,y=ctri))
           }		   
		   gg <- gg + geom_ribbon(aes(ymin=ic_low95, ymax=ic_up95),fill = couleur,alpha=.2)+ geom_line(size=1,colour=couleur)
            gg <- gg +  geom_point(data = subset(ggGamData,year %in% realYear),size=3,colour=couleur) + geom_point(data = subset(ggGamData,year %in% realYear),size=1.5,colour="white")
            gg <- gg + labs(y=titreY,x=titreX,title=titre)+scale_x_continuous(breaks=pretty_breaks())

            ggsave(paste("Output/fig",indic,"_plot",id,".png",sep=""),gg)
            cat("\n--> Output/fig",indic,"_plot",id,".png\n",sep="")

            tabPredict <- subset(ggGamData,year %in% realYear)########### Tableau des resultats pour ne prendre que les valeurs d'IC pour l'année pas entre les années (spline sur annee) !!!plus utilisé!! / Table of the results not taking confidence interval between year but at each year (because of the spline of year)
            colnames(tabPredict)[1:2] <- c("annee",paste(indic,"_predict",sep=""))
            ##tabgamm <- merge(tabfgamm,tabPredict,by="annee") #### Desactivation car merge sortie de modèles différents (le modèle dont on tire les coef de regression pour année, avec spline sur les coordonnées geo vs celui pour faire la figure avec splin sur année uniquement)  / not use anymore (as before) because use the results of the restricted model with the spline on the year while the better analysis is on full model with the spline on gps coordinates
			tabgamm <-  tabfgamm  #### remplace la ligne au dessus  / replace the line above

} else {
    if(init_1989) { ########### Pour utiliser les données csi avant 2001 / in order to use csi data from before 2001
        histo <- read.csv2("csi_init.csv")  ###############   Csi calculé avant 2001 / Data of the csi calculated before 2001 
        init <- subset(histo,annee == 2001)$csi

        tabgamm <-  tabfgamm
        tabgamm$coef <-  tabgamm$coef+init
        histo <- subset(histo,annee<=2001)
        tabgamm <- subset(tabgamm,annee>2001)
        tabHisto <- data.frame(model= tabgamm$model[1],annee = histo$annee,
                               coef = histo$csi,se=histo$se,pval=NA,signif=NA)
        tabgamm <- rbind(tabgamm,tabHisto)
        tabgamm <- tabgamm[order(tabgamm$annee),]

        tabgamm$ic_up95<-tabgamm$coef-1.96*as.vector(tabgamm$se)

           tabgamm$ic_low95<-tabgamm$coef+1.96*as.vector(tabgamm$se)

    }
#    browser()
         ## The ggplot:
            gg <- ggplot(data=tabgamm,aes(x=annee,y=coef))
           gg <- gg + geom_ribbon(aes(ymin=ic_low95, ymax=ic_up95),fill = couleur,alpha=.2)+ geom_line(size=1,colour=couleur)
            gg <- gg +  geom_point(size=3,colour=couleur)
            gg <- gg + labs(y=titreY,x=titreX,title=titre)+scale_x_continuous(breaks=pretty_breaks())

            ggsave(paste("Output/fig",indic,"_plot",id,".png",sep=""),gg)
        cat("\n--> Output/fig",indic,"_plot",id,".png\n",sep="")


    write.csv(tabgamm,paste("Output/",indic,"_gammPlot_",id,".csv",sep=""),row.names=FALSE) ####
        cat("\n  --> Output/",indic,"_gammPlot_",id,".csv\n",sep="")





}

            cat("\nEstimation de la tendence  ",indic,"~ year+s(longitude_grid_wgs84,latitude_grid_wgs84,bs='sos')\n")
if (indic=="csi"){
            gammc <- gamm(csi~year+s(longitude_grid_wgs84,latitude_grid_wgs84,bs="sos"), data=dd,random=reStruct(object = ~ 1| id_plot, pdClass="pdDiag"),correlation=corAR1(form=~year))### spline sur les coordonnées, effet aleatoire des carres sur ordonnée à l'origine, methode autoregressive sur l'année N-1  / spline on the gps coordinates, random effect of the plots on intercept, autoregressive method on the year-1
}
if (indic=="cti"){
            gammc <- gamm(cti~year+s(longitude_grid_wgs84,latitude_grid_wgs84,bs="sos"), data=dd,random=reStruct(object = ~ 1| id_plot, pdClass="pdDiag"),correlation=corAR1(form=~year))
}
if (indic=="ctri"){
            gammc <- gamm(ctri~year+s(longitude_grid_wgs84,latitude_grid_wgs84,bs="sos"), data=dd,random=reStruct(object = ~ 1| id_plot, pdClass="pdDiag"),correlation=corAR1(form=~year))
}
            sgammc <-summary(gammc$gam)



            coefannee <- sgammc$p.coeff[2]  #### coefficient de regression de la variable année / regression coefficient of the variable "year"
            ## erreur standard / standard error
            erreuran <- sgammc$se[2]
            ## p value
            pval <-  sgammc$p.pv[2]
			
			#### Calcul des intervalles de confiances / calculation of the confidence intervals

        if(ic) {
            # gammc.sim <- sim(gammc)######################  VERSION ROMAIN mais fct sim() marche pas avec Gamm / old version using function sim() but did not work with Gamm models
            # ic_inf_sim <- c(0,tail(apply(coef(gammc.sim), 2, quantile,.025),pasdetemps))
            # ic_sup_sim <- c(0,tail(apply(coef(gammc.sim), 2, quantile,.975),pasdetemps))
			icalpha05 <- as.data.frame(confint(gammc$gam))[2,1:2]  ########## VERSION BENJ
			ic_inf_sim <- icalpha05[,1]
			ic_sup_sim <- icalpha05[,2]
        } else
        {
            ic_inf_sim <- "not assessed"
            ic_sup_sim <- "not assessed"
        }
if (indic=="csi"){
            tabcgamm <- data.frame(model = "gamm numeric(year) plot",annee = NA,coef = coefannee,se = erreuran,pval,signif = pval<seuilSignif, csi_predict= NA , ic_low95 = ic_inf_sim, ic_up95 = ic_sup_sim)#### recupère les resultats des modèles avec interval de confiance / retrieve results of the models used with confidence interval
}
if (indic=="cti"){
            tabcgamm <- data.frame(model = "gamm numeric(year) plot",annee = NA,coef = coefannee,se = erreuran,pval,signif = pval<seuilSignif, cti_predict= NA , ic_low95 = ic_inf_sim, ic_up95 = ic_sup_sim)#### 
}
if (indic=="ctri"){
            tabcgamm <- data.frame(model = "gamm numeric(year) plot",annee = NA,coef = coefannee,se = erreuran,pval,signif = pval<seuilSignif, ctri_predict= NA , ic_low95 = ic_inf_sim, ic_up95 = ic_sup_sim)#### 
}

            tabgamm <- tabgamm[,colnames(tabcgamm)]


            tabgamm <- rbind(tabgamm,tabcgamm)

            write.csv(tabgamm,paste("Output/",indic,"_gammPlot_",id,".csv",sep=""),row.names=FALSE)
        cat("\n  --> Output/",indic,"_gammPlot_",id,".csv\n",sep="")
		
		}
		
		
if (methode == "lmer") {
            cat("Method : lmer \n")


###################
# browser()
   ### Utilisation des modèles mixtes pour obtenir les tendances d evolution par an du csi cti ou ctri / Use of mixte model for the estimation of the annual variations of the csi cti or ctri 
            cat("\nEstimation de la variation annuelle lmer(",indic,"~ factor(year)+(1|id_plot)\n",sep="")
if (indic=="csi"){            
			md.f <- lmer(csi~ factor(year)+(1|id_plot),data=dd)  ##### effet aleatoire liés aux carrés sur l'ordonnée à l'origine / random effects of plots on intercept 
}
if (indic=="cti"){            
			md.f <- lmer(cti~ factor(year)+(1|id_plot),data=dd)
}
if (indic=="ctri"){            
			md.f <- lmer(ctri~ factor(year)+(1|id_plot),data=dd)
}
            smd.f<-summary(md.f)    
            coefdata.f <-  as.data.frame(smd.f$coefficients)
            coefdata.f <- data.frame(model="Annual fluctuation", variable = rownames(coefdata.f),coefdata.f)

            ggdata <<- data.frame(year=c(1989,as.numeric(substr(coefdata.f$variable[-1],13,16))),
                                 estimate=c(0,coefdata.f$Estimate[-1]),
                                 se=c(0,coefdata.f$Std..Error[-1]))   #####################  resultat du modèle / results of the models
            #ggdata$estimate <-  ggdata$estimate
             #ggdata$se.supR <- ggdata$estimate +  ggdata$se ############################################################################## METHODE ROMAIN 
             #ggdata$se.infR <- ggdata$estimate -  ggdata$se
            # ggdata$estimate2 <- c(coefdata.f$Estimate[1],coefdata.f$Estimate[1] + coefdata.f$Estimate[-1])
			# ggdata$se.sup2 <- ggdata$estimate2 +  ggdata$se
            # ggdata$se.inf2 <- ggdata$estimate2 -  ggdata$se
		
		prof <- profile(md.f)  #### Nouvel interval de confiance avec utilisation du logarithme des ecarts types / logarithms of standard deviations are used, while varianceProf converts from the standard-deviation to the variance scale
		MODconfint <- confint(prof) #### plus rapide de passer par la fonction profile mais pas indispensable fonctionne aussi directement sur modele mixte md.f  / more rapid using both function profile and confint but works also directly on output of the model 
		se.sup <- MODconfint[2:nban+2,2]
		ggdata$se.sup <- c(0,se.sup) 
		se.inf <- MODconfint[2:nban+2,1]
        ggdata$se.inf <- c(0,se.inf)   			
coefdata.f$se.inf <- ggdata$se.inf
coefdata.f$se.sup <- ggdata$se.sup
browser()
            #gg <<- ggplot(ggdata,aes(x=year,y=estimate))+ geom_ribbon(ymin=ggdata$se.infR,ymax=ggdata$se.supR,alpha=.25)+geom_errorbar(ymin=ggdata$se.infR,ymax=ggdata$se.supR,width=0,alpha=.25)+ geom_point() + geom_line() + ylim(min(ggdata$se.infR),max(ggdata$se.supR)) + labs(x="Years",y=paste(indic," variation",sep="")) #####  AVEC INTERVAL ROMAIN
			gg <<- ggplot(ggdata,aes(x=year,y=estimate))+ geom_ribbon(ymin=ggdata$se.inf,ymax=ggdata$se.sup,alpha=.25)+geom_errorbar(ymin=ggdata$se.inf,ymax=ggdata$se.sup,width=0,alpha=.25)+ geom_point() + geom_line() + ylim(min(ggdata$se.inf),max(ggdata$se.sup)) + labs(x="Years",y=paste(indic," variation",sep="")) #####

            ggfile <- paste("output/",indic,id,".png",sep="")
            ggsave(ggfile,gg)


############ Estimation de la tendance sur la periode étudiée  / Trends estimation on the time period studied
			cat("\nEstimation de la tendance lmer(",indic,"~ factor(year)+(1|id_plot)\n",sep="")
if (indic=="csi"){            
			md.c <- lmer(csi~ year+(1|id_plot),data=dd)##### effet aleatoire liés aux carrés sur l'ordonnée à l'origine / random effects of plots on intercept 
}
if (indic=="cti"){            
			md.c <- lmer(cti~ year+(1|id_plot),data=dd)
}
if (indic=="ctri"){            
			md.c <- lmer(ctri~ year+(1|id_plot),data=dd)
}        
            smd.c<-summary(md.c)

            coefdata.c <-  as.data.frame(smd.c$coefficients)
profc=profile(md.c) ######### Ajout des intervalles de confiances / addition of the confidence intervals
MODconfint=confint(profc)
se.inf=MODconfint[4,1]
se.sup=MODconfint[4,2]
			coefdata.c <- data.frame(model = "Linear trend", variable = rownames(coefdata.c),coefdata.c,se.inf,se.sup)

            coefdata <- rbind(coefdata.c,coefdata.f)


write.csv(coefdata,paste("Output/lmer_coefficient_",indic,id,".csv",sep=""),row.names=FALSE)
write.csv(ggdata,paste("Output/ggdata_",indic,id,".csv",sep=""),row.names=FALSE)


            smd.file <- paste("output/summary_lmer_",indic,"_",id,".txt",sep="")



#####################






        }
		
		
		
		
		
		
		}
		
		
################## 
###  Do your analysis

csi_cti_ctri(dd=dd,indic="csi",ic=TRUE,plot_smooth = TRUE,methode="lmer")  ##### exemple pour l'indicateur csi avec interval de confiance et utilisant modele mixte / example for the csi index with confidence interval using the mixte model 