#!/usr/bin/env Rscript

######################################################################################################################
############## CALCULATE AND PLOT EVOLUTION OF SPECIES POPULATION  function:main.glm    ##############################
######################################################################################################################

#### Based on Romain Lorrillière R script
#### Modified by Alan Amosse and Benjamin Yguel for integrating within Galaxy-E

suppressMessages(library(lme4))
suppressMessages(library(ggplot2))
suppressMessages(library(speedglm))
suppressMessages(library(arm))
suppressMessages(library(reshape))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))

source("main_glm_functions.r")


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) { #TODO !=5 ? 6? to fix later
    stop("At least 5 arguments must be supplied :\n- Input data filtered (csv file)\n- Tab species (csv file)\n- True/False to perform a standard glm with IC\n- List of species to perform the analys on\n- First year\n- Last year\n\n  ", call.=FALSE) #TODO a mettre à jour
} else {
    Datafilteredfortrendanalysis<-args[1] #datafile may result from the function "FiltreEspeceRare"    
    tabSpecies<-args[2] # tab of species with information on species specialization and if species are indicators
#    spExclude <- args [4] # list of species that will be excluded #Pas pris en compte ce truc ? #Plutot à faire dans les étapes de filter avant
#    AssessIC <-args[3] # TRUE or FALSE perform a "standard" glm with confidance interval or speedglm without CI #A mettre en relation avec l'envie de figure.If figure=True -> AssessIC=True
#    Listsp<-args[4] #To analyse a few species ? faster
#    firstYear<-args[5] #Set a Ystart-Yend analyse 
#    lastYear<-args[6]      
#    figure<-args[7] #True/False create figure
}




# Import data 
tabCLEAN <- read.csv(Datafilteredfortrendanalysis,sep=";",dec=".") 
tabsp <- read.csv(tabSpecies,sep=";",dec=".")  

if(ncol(tabCLEAN)<4){ #Check file integrity
    stop("The file don't have at least 4 variables", call.=FALSE)
}



firstYear<-2007
lastYear<-2009
id="France" #fix id to have same output dir at every execution
assessIC=TRUE #TRUE Si on veux figure sinon erreur #Error: Discrete value supplied to continuous scale
listSp=NULL
annees<- firstYear:lastYear
figure=TRUE
description=TRUE
tendanceSurFigure=TRUE
tendanceGroupeSpe=FALSE
seuilOccu=14
seuilAbond=NA
ecritureStepByStep=FALSE

dir.create(paste("Output/",id,sep=""),recursive=TRUE)
cat(paste("Create Output/",id,"\n",sep=""))
dir.create(paste("Output/",id,"/Incertain/",sep=""),recursive=TRUE)
cat(paste("Create Output/",id,"Incertain/\n",sep=""))

main.glm(id,tabCLEAN,assessIC,listSp,tabsp,annees=NULL,figure,description,tendanceSurFigure,tendanceGroupSpe,seuilOccu,seuilAbond,ecritureStepByStep)
