#!/usr/bin/env Rscript

##################################################################################################################################################
############## CALCULATE AND PLOT EVOLUTION OF SPECIES POPULATION BY SPECIALIZATION GROUP  function:analyse.Groupe  ##############################
##################################################################################################################################################

#### Based on Romain Lorrillière R script
#### Modified by Alan Amosse and Benjamin Yguel for integrating within Galaxy-E

suppressMessages(library(lme4))
suppressMessages(library(ggplot2))
suppressMessages(library(speedglm))
suppressMessages(library(arm))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))

check_file<-function(dataset,err_msg,vars,nb_vars){
    if(ncol(dataset)!=nb_vars){ #Verifiction de la présence du bon nb de colonnes, si c'est pas le cas= message d'erreur / checking for right number of columns in the file if not = error message
        cat("\nerr nb var\n") 
        stop(err_msg, call.=FALSE)
    }

    for(i in vars){
        if(!(i %in% names(dataset))){
            stop(err_msg,call.=FALSE)
        }
    }
}

###########
#delcaration des arguments et variables/ declaring some variables and load arguments

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=6) {
    stop("The tool need the following inputs :\n\n- A yearly species variations data set (.tabular). It may come from the main glm tool.\n- A species global tendencies dataset (.tabular). It may come from the main glm tool.\n- A species table filtered (.tabular). It may come from the Filter rare species tool.\n- An id to fix output repository name.\n- A list of species to exclude, can be empty.\n- A bias file.\n\n", call.=FALSE) #si pas d'arguments -> affiche erreur et quitte / if no args -> error and exit1
} else {
    donnees<-args[1] ###### Nom du fichier avec extension "***variationsAnnuellesEspece***.tabular", peut provenir de la fonction "mainglm" / file name without the file type "***variationsAnnuellesEspece***.tabular", may result from the function "mainglm"    
    donneesTrend <- args[2] ####### Nom du fichier avec extension "***tendanceGlobalEspece***.tabular", peut provenir de la fonction "mainglm" / / file name without the file type "***tendanceGlobalEspece***.tabular", may result from the function "mainglm"    
    tabSpecies<-args[3] ###### Nom du fichier avec extension ".typedefichier", peut provenir de la fonction "FiltreEspeceRare" / file name without the file type ".filetype", may result from the function "FiltreEspeceRare"  
    id<-args[4]  ##### nom du dossier de sortie des resultats / name of the output folder
    spExclude <- args [5] ##### liste d'espece qu on veut exclure de l analyse  / list of species that will be excluded
    tBiais <-args [6] ##########   fichier contenant le biais de détéction en fonction des occurances, obtenu à partir d'un modéle théorique de dynamique de pop et de survey / the file containing the detection bias depending on occurance data obtained with theoretical model of population dynamic and survey
}




#Import des données / Import data 
tBiais=read.table(tBiais,sep="\t",dec=".",header=TRUE) ###### charge le fichier contenant le biais de détéction en fonction des occurances, obtenu à partir d'un modéle théorique de dynamique de pop et de survey / load the file containing the detection bias obtained with theoretical model of population dynamic and survey
donnees <-  read.table(donnees,sep="\t",dec=".",header=TRUE) #### charge le fichier de resultat sur les tendances annuelles par espèce / load annual population evolution trend for each species obtained with the function mainglm
donneesTrend <- read.table(donneesTrend,sep="\t",dec=".",header=TRUE)#### charge le fichier de resultat sur les tendances sur la periode etudiée par espèce / load population evolution trend on the studied period for each species obtained with the function mainglm
tabsp <- read.table(tabSpecies,sep="\t",dec=".",header=TRUE)   #### charge le fichier de donnees sur nom latin, vernaculaire et abbreviation, espece indicatrice ou non / load the file with information on species specialization and if species are indicators




vars_donnees<-c("id","code_espece","nom_espece","indicateur","annee","abondance_relative","IC_inferieur","IC_superieur","erreur_standard","p_value","significatif","nb_carre","nb_carre_presence","abondance")
err_msg_donnees<-"\nThe yearly species variation dataset doesn't have the right format. It need to have following 14 variables :\n- id\n- code_espece\n- nom_espece\n- indicateur\n- annee\n- abondance_relative\n- IC_inferieur\n- IC_superieur\n- erreur_standard\n- p_value\n- significatif\n- nb_carre\n- nb_carre_presence\n- abondance\n"

vars_donneesTrend<-c("id","code_espece","nom_espece","indicateur","nombre_annees","premiere_annee","derniere_annee","tendance","IC_inferieur","IC_superieur","pourcentage_variation","erreur_standard","p_value","significatif","categorie_tendance_EBCC","mediane_occurrence","valide","raison_incertitude")
err_msg_donneesTrend<-"\nThe species global tendances dataset doesn't have the right format. It need to have following 18 variables :\n- id\n- code_espece\n- nom_espece\n- indicateur\n- nombre_annees\n- premiere_annee\n- derniere_annee\n- tendance\n- IC_inferieur\n- IC_superieur\n- pourcentage_variation\n- erreur_standard\n- p_value\n- significatif\n- categorie_tendance_EBCC\n mediane_occurrence\n valide\n raison_incertitude\n"

vars_tabsp<-c("espece","nom","nomscientific","indicateur","specialisation")
err_msg_tabsp<-"\nThe species dataset filtered doesn't have the right format. It need to have the following 4 variables :\n- espece\n- nom\n- nomscientific\n- indicateur\n- specialisation\n"


check_file(donnees,err_msg_donnees,vars_donnees,14)
check_file(donneesTrend,err_msg_donneesTrend,vars_donneesTrend,18)
check_file(tabsp,err_msg_tabsp,vars_tabsp,5)


spsFiltre=unique(levels(donnees$code_espece)) #### Recupère la liste des especes du tabCLEAN qui ont été sélectionnée et qui ont passé le filtre / retrieve species name that were selected and then filtered before

spExclude=subset (tabsp, !(espece %in% spsFiltre)) #### liste des espèces exclu par le filtre ou manuellement / List of species excluded manually or by the filter from the analyses 
tabsp=subset (tabsp, (espece %in% spsFiltre)) #### Enlève les espèces qui n'ont pas passé le filtre ou exclu manuellement pour les analyses / keep only selected species and species with enough data
sp=as.character(tabsp$espece)  ##### liste des espece en code ou abbreviation gardées pour les analyses ### arg de la fonction  DECLARE AUSSI APRES DS FONCTION  / list of the code or abbreviation of the species kept for the analyses
tabsp=data.frame(tabsp,sp)### rajoute une colonne identique appelé sp / add new column called sp


## creation d'un dossier pour y mettre les resultats / create folder for the output of the analyses   ###### NORMALEMENT DOIT ËTRE DEJ2 CREER POUR LES SORTIES TENDANCES PAR SPS DONC PAS SUR QU IL FAUT REFAIRE CETTE ETAPE

dir.create(paste("Output/",id,sep=""),recursive=TRUE,showWarnings=FALSE)
cat(paste("Create Output/",id,"\n",sep=""))
dir.create(paste("Output/",id,"/Incertain/",sep=""),recursive=TRUE,showWarnings=FALSE)
cat(paste("Create Output/",id,"Incertain/\n",sep=""))




## Analyse par groupe de specialisation Ã  partir des resulats de variation d'abondance par especes
## id identifiant de la session
## ICfigureGroupeSp affichage des intervalles de confiances sur la figure
## correctionAbondanceNull correction des abondance NULL
analyseGroupe <- function(id="france",tabsp=tabsp,donnees=donnees,donneesTrend=donneesTrend,ICfigureGroupeSp=TRUE,powerWeight=2,
                          correctionAbondanceNull = 0.000001,
                          groupeNom = c("generaliste","milieux batis","milieux forestiers","milieux agricoles"),
                          groupeCouleur = c("black","firebrick3","chartreuse4","orange")) {
    
    ## donnees tendances globales
    donneesTrend <- subset(donneesTrend, select = c(code_espece,valide,mediane_occurrence))
	
    ## table de reference espece
    tabsp <- subset(tabsp, select= c(sp,nom,indicateur, specialisation))
    donnees <- merge(donnees,donneesTrend,by="code_espece")
    donnees <- merge(donnees,tabsp,by.x="code_espece",by.y="sp")
    ## table de correspondance de biais en fonction des medianes des occuerences
	
    
    nameFileSpe <-  paste("Output/",id,"/variationsAnnuellesGroupes_",id, ############# Declare le fichier de sortie des variations annuelles par groupe / declare the name of the outputfile for annual population evolution trend by group 
                          ".tabular",sep="" )
    nameFileSpepng <-  paste("Output/",id,"/variationsAnnuellesGroupes_",id, ############# Declare le fichier de sortie graphique des variations annuelles par groupe / declare the name of the graphical output file for annual population evolution trend by group
                             ".png",sep="" )
    
    grpe <- donnees$specialisation
    
    
    ff <- function(x,y) max(which(y<=x)) ## fonction pour recherche valeur plus petite que seuil d'occurence / function to look for values under the occurence threshold
     
    IncertW <- ifelse(donnees$valide=="Incertain",tBiais$biais[sapply(as.vector(donnees$mediane_occurrence),ff,y=tBiais$occurrenceMed)],1) ## pr verifier poids de l'espèce dans analyse, récupére seuil occurence minimum pour lequel tendance pas bonne, et compare avec mediane occurence des données  / to check the weight of species in the analysis, this retrieve occurence threshold with wich real occurence measured on data are compared in order to verify the accuracy of the trend measurment
    ## poids du Ã  la qualitÃ© de l'estimation
                                        #   erreur_stW <- 1/((donnees$erreur_st+1)^powerWeight)
                                        #	erreur_stW <- ifelse( is.na(donnees$IC_superieur),0,erreur_stW)
    erreur_stW <- ifelse(is.na(donnees$IC_superieur),0,1)##### 
    ## calcul du poids total de chaque espèce / calcul of the weight of each species 
    W <- IncertW * erreur_stW
    
    ## variable de regroupement pour les calculs par groupe de specialisation et par an / variables gathered to identify group for the calculation (per specialization and per year)
    grAn <- paste(donnees$specialisation,donnees$annee,sep="_")
    ## data frame pour le calcul / dataframe made for the calcul
    dd <- data.frame(grAn,annee = donnees$annee, grpe,W,ab=donnees$abondance_relative,ICinf= donnees$IC_inferieur, ICsup= ifelse(is.na(donnees$IC_superieur),10000,donnees$IC_superieur)) 
    ## table resumer de tous les poids / table to sum up the weights of each species depending on the incertainty in the calcul of the poulation evolution trends
    ddd <- data.frame(code_espece = donnees$code_espece,nom_espece = donnees$nom_espece,annee = donnees$annee, 
                      groupe_indicateur = grpe,
                      poids_erreur_standard = round(erreur_stW,3), poids_incertitude = round(IncertW,3),poids_final = round(W,3),
                      abondance_relative=donnees$abondance_relative,
                      IC_inferieur= donnees$IC_inferieur, 
                      IC_superieur= ifelse(is.na(donnees$IC_superieur),10000,donnees$IC_superieur),
                      valide = donnees$valide, mediane_occurrence = donnees$mediane_occurrence) 

    nomFileResum <- paste("Resultats/",id,"/donneesGroupes_",id, ###### declaration du nom du repertoire et des fichiers de sortie / declaring the name of the output folder and files  
                          ".tabular",sep="" )
    write.table(ddd,nomFileResum,row.names=FALSE,sep="\t",dec=".")
    cat(" <--",nomFileResum,"\n")
    
    ## calcul des moyennes pondÃ©rÃ© par groupe par an et pour les abondance et les IC	/ calcul of weighted means per specialization group and per year for the abundance and confidence interval
    for(j in 5:7) dd[,j] <- ifelse(dd[,j]==0,correctionAbondanceNull,dd[,j])	
    ag <- apply(dd[,5:7], 2,  ######## sur les abondances relatives, les ICinf et ICsup
                function(x) {
                    sapply(split(data.frame(dd[,1:4], x), dd$grAn),  ###### fait les moyennes pondérés par groupe grAn / calculate the weighted mean by group grAn
                           function(y) round(geometriqueWeighted(y[,5], w = y$W),3))
                })
    ##	gg <- subset(dd,as.character(dd$grAn)=="milieux forestier_2014")  #############################################################

    ag <- ifelse(is.na(ag),1,ag)
    ag <- as.data.frame(ag)
    ag$grAn <-  rownames(ag)
    dbon <- subset(donnees,valide=="bon")
    dIncert <- subset(donnees,valide=="Incertain")
    ## calcul nombre d'espece "bonne" pour le calcul / calculating the number of species with low level of incertainty, "good" species 
    bon <- tapply(dbon$nom,dbon$specialisation,FUN=function(X)length(unique(X)) )
    bon <- ifelse(is.na(bon),0,bon)
    tbon <- data.frame(groupe=names(bon),bon)
    ## calcul nombre d'especes "incertaines" pour le calcul / calculating the number of species with high level of incertainty, "bad" species
    Incert <- tapply(dIncert$nom,dIncert$specialisation,FUN=function(X)length(unique(X)) )
    Incert <- ifelse(is.na(Incert),0,Incert)
    tIncert <- data.frame(groupe=names(Incert),Incertain=Incert)

    tIncert <- merge(tIncert,tbon,by="groupe")
    
    ## table de données avec les moyennes ponderees par groupe / table of the data with the weighted mean by group 
    da <- merge(unique(dd[,1:3]),ag,by="grAn")[,-1]
    colnames(da) <- c("annee","groupe","abondance_relative","IC_inferieur","IC_superieur")

    da$annee <- as.numeric(da$annee)
    da <-  merge(da,tIncert,by="groupe") #### ajoute le nombre d'espece "incertaines" et "bonne" aux resultats  / add the number of "good" and "bad" species to the overall resutls
    da <- subset(da, groupe != "non")
    colnames(da)[6:7] <-  c("nombre_especes_incertaines","nombre_espece_bonnes")
    a <- data.frame(id,da)
    write.table(da,file=nameFileSpe,row.names=FALSE,quote=FALSE,sep="\t",dec=".")

    cat(" <--",nameFileSpe,"\n")
    yearsrange <- c(min(da$annee),max(da$annee))
    
    ## figure par ggplot2  / plots with ggplot2
    titre <- paste("Variation de l'indicateur groupe de spÃ©cialisation",sep="")

    vecCouleur <- setNames(groupeCouleur,groupeNom)
                                        #browser()
    p <- ggplot(data = da, mapping = aes(x = annee, y = abondance_relative, colour=groupe,fill=groupe))
    p <- p + geom_hline(aes(yintercept = 1), colour="white", alpha=1,size=1.2) 
    if(ICfigureGroupeSp)
        p <- p + geom_ribbon(mapping=aes(ymin=IC_inferieur,ymax=IC_superieur),linetype=2,alpha=.1,size=0.1) 
    p <- p + geom_line(size=1.5)
    p <- p +  ylab("") + xlab("AnnÃ©e")+ ggtitle(titre) 
    if(!is.null(groupeNom)) p <- p + scale_colour_manual(values=vecCouleur, name = "" )+
                                scale_x_continuous(breaks=unique(da$annee))
    if(!is.null(groupeNom)) p <- p +  scale_fill_manual(values=vecCouleur, name="")
    p <- p +  theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank()) 
    ggsave(nameFileSpepng, p,width=17,height=10,units="cm")

                                        #   cat(" <==",nameFileSpepng,"\n")
    
    ## calul pour chaque groupe une pente de regression d'evolution des abondances sur la periode étudiée / calculating for each group the regression slope for the abundance evolution on the studied period
    vecSpe <- unique(da$groupe)
    datasum <- data.frame(groupe=NULL,tendance=NULL,pourcentage_variation=NULL)
    for(spe in 1:4){
                                        # print(spe)
        subtab <- subset(da,groupe==vecSpe[spe])
        if(nrow(subtab)>1) {
            sumlm <- summary(lm(abondance_relative~annee,data=subtab)) ##### recupère les resultats du modèle linéaire / retrieve the results of the linear model
            subdatasum <- data.frame(groupe=vecSpe[spe],
                                     tendance=round(sumlm$coefficients[2,1],3),
                                     pourcentage_variation=round(sumlm$coefficients[2,1]*(nrow(subtab)-1)*100,3)) #### assemble les resultats pour en faire une sortie  /  bring together the results for an output file
            datasum <- rbind(datasum,subdatasum)
            
        }
        
    }
    datasum <- merge(datasum,tIncert,by="groupe") #### 
    datasum <- data.frame(id,datasum)
                                        #datasum$cat_tendance_EBCC <- affectCatEBCC(trend,pVal,ICinf,ICsup
    namefilesum <- paste("Resultats/",id,"/tendancesGlobalesGroupes_",id,
                         ".tabular",sep="" )
    write.table(datasum,file=namefilesum,row.names=FALSE,quote=FALSE,sep="\t",dec=".")
    cat(" <--",namefilesum,"\n")
}




################## 
###  Do your analysis
analyseGroupe(id,tabsp=tabsp,donnees=donnees,donneesTrend=donneesTrend,ICfigureGroupeSp=TRUE,groupeNom = groupeNom,groupeCouleur=groupeCouleur)
