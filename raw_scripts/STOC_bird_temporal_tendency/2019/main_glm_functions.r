#!/usr/bin/env Rscript


## renvoie la categorie EBCC de la tendance en fonction
## trend l'estimateur de la tendance
## pVal la p value
## ICinf ICsup l intervalle de confiance a 95 pourcent
affectCatEBCC <- function(trend,pVal,ICinf,ICsup){
    catEBCC <- ifelse(pVal>0.05,
               ifelse(ICinf < 0.95 | ICsup > 1.05,"Incertain","Stable"),
               ifelse(trend<1,
               ifelse(ICsup<0.95,"Fort déclin","Déclin modéré"),
               ifelse(ICinf>1.05,"Forte augmentation","Augmentation modérée")))
    return(catEBCC)
}


## fonction general de calcul de la variation temporelle et de la tedence generale
## la fonction genere aussi les graphiques
main.glm <- function(id,donneesAll,assessIC= TRUE,listSp=NULL,tabsp,annees=NULL,figure=TRUE,description=TRUE,tendanceSurFigure=TRUE,tendanceGroupSpe = FALSE,
                     seuilOccu=14,seuilAbond=NA,ecritureStepByStep=FALSE) {

    ##  donneesAll=data;listSp=sp;annees=firstYear:lastYear;figure=TRUE;description=TRUE;tendanceSurFigure=TRUE;tendanceGroupSpe = FALSE;
    ##                   seuilOccu=14;seuilAbond=NA;ecritureStepByStep=TRUE
    
                                        #donneAll = data;listSp=NULL;annees=NULL;echantillon=1;methodeEchantillon=NULL;
                                        #figure=TRUE;description=TRUE;tendanceSurFigure=TRUE;tendanceGroupSpe = FALSE;
                                        #seuilOccu=14;seuilAbond=NA;ecritureStepByStep=FALSE
                                        # browser()
    require(arm)
    require(ggplot2)

    filesaveAn <-  paste("Output/",id,"/variationsAnnuellesEspece_",id,".csv",
                         sep = "")
    filesaveTrend <-  paste("Output/",id,"/tendanceGlobalEspece_",id,".csv",
                            sep = "")

    fileSaveGLMs <-  paste("Output/",id,"/listGLM_",id,sep = "")

    
    ## seuil de significativite
    seuilSignif <- 0.05
    
    ## tabsp table de reference des especes
    rownames(tabsp) <- tabsp$espece
    
    
    ##vpan vecteur des panels de la figure
    vpan <- c("Variation abondance")
    if(description) vpan <- c(vpan,"Occurrences","Abondances brutes")
                                        # nomfile1 <- paste("Donnees/carre2001-2014REGION_",id,".csv",sep="")

    ## des variable annees
    annee <- sort(unique(donneesAll$annee))
    nbans <- length(annee)
    pasdetemps <- nbans-1
    firstY <- min(annee)
    lastY <- max(annee)

    ## Ordre de traitement des especes
    spOrdre <- aggregate(abond~espece,data=donneesAll,sum)
    spOrdre <- merge(spOrdre,tabsp,by="espece")
    
    spOrdre <- spOrdre[order(as.numeric(spOrdre$indicateur),spOrdre$abond,decreasing = TRUE),]
    
    
    listSp <- spOrdre$espece
    i <- 0
    nbSp <- length(listSp)
                                        #	browser()
    ## analyse par espece
    ### browser()
    ## affichage des especes conservées pour l'analyse
    cat("\n",nbSp," Espèces conservées pour l'analyse\n\n",sep="")
    rownames(tabsp) <- tabsp$espece
    tabCons <- data.frame(Code_espece = listSp, nom_espece = tabsp[as.character(listSp),"nom"])
    print(tabCons)  
    cat("\n\n",sep="")
    flush.console()




    ## initialisation de la liste de sauvegarde

    


    ##browser()
    
    for (sp in listSp) {
    ###        if(sp=="PHYCOL")            browser()

        i <- i + 1
        ## d data pour l'espece en court    
        d <- subset(donneesAll,espece==sp)
        ## info sp
        nomSp <- as.character(tabsp[sp,"nom"])
        cat("\n(",i,"/",nbSp,") ",sp," | ", nomSp,"\n",sep="")
        flush.console()
        #### shortlist espece fait partie des especes indiactrice reconnue à l'echelle national
        ### shortlist <- tabsp[sp,"shortlist"]
        ## indic espèce utilisé pour le calcul des indicateurs par groupe de specialisation 
        indic <- tabsp[sp,"indicateur"]

        ## Occurrence
        ## nb_carre nombre de carre suivie par annee
        nb_carre = tapply(rep(1,nrow(d)),d$annee,sum)
        ## nb_carre_presence nombre de carre de presence par annee
        nb_carre_presence = tapply(ifelse(d$abond>0,1,0),d$annee,sum)
        ## tab2 table de resultat d'analyse
        tab2 <- data.frame(annee=rep(annee,2),val=c(nb_carre,nb_carre_presence),LL = NA,UL=NA,
                           catPoint=NA,pval=NA,
                           courbe=rep(c("carre","presence"),each=length(annee)),panel=vpan[2])
        tab2$catPoint <- ifelse(tab2$val == 0,"0",ifelse(tab2$val < seuilOccu,"infSeuil",NA))
        
        ## abondance brut
        ## abond abondance par annee
        abond <- tapply(d$abond,d$annee,sum)
        ## tab3 tab3 pour la figure
        tab3 <- data.frame(annee=annee,val=abond,LL = NA,UL=NA,catPoint=NA,pval=NA,courbe=vpan[3],panel=vpan[3])
        tab3$catPoint <- ifelse(tab3$val == 0,"0",ifelse(tab3$val < seuilAbond,"infSeuil",NA))

        ## GLM variation d abondance
       formule <- as.formula("abond~as.factor(carre)+as.factor(annee)")
       if(assessIC) {
           glm1 <- glm(formule,data=d,family=quasipoisson)
       } else {
           glm1 <- try(speedglm(formule,data=d,family=quasipoisson()))
           if(class(glm1)[1]=="try-error")
               glm1 <- glm(formule,data=d,family=quasipoisson) 
       }
       sglm1 <- summary(glm1)
       sglm1 <- coefficients(sglm1)
       sglm1 <- tail(sglm1,pasdetemps)
       coefan <- as.numeric(as.character(sglm1[,1]))
        ## coefannee vecteur des variation d'abondance par annee back transformee
        coefannee <- c(1,exp(coefan))
        erreuran <- as.numeric(as.character(sglm1[,2]))
        ## erreur standard back transformee
        erreurannee1 <- c(0,erreuran*exp(coefan))
        pval <- c(1,as.numeric(as.character(sglm1[,4])))
        
        ## calcul des intervalle de confiance
        if(assessIC) {
        glm1.sim <- sim(glm1)
        ic_inf_sim <- c(1,exp(tail(apply(coef(glm1.sim), 2, quantile,.025),pasdetemps)))
        ic_sup_sim <- c(1,exp(tail(apply(coef(glm1.sim), 2, quantile,.975),pasdetemps)))
        } else {
            ic_inf_sim <- NA
            ic_sup_sim <- NA
 
        }
        
        
        ## tab1 table pour la realisation des figures
        tab1 <- data.frame(annee,val=coefannee,
                           LL=ic_inf_sim,UL=ic_sup_sim,
                           catPoint=ifelse(pval<seuilSignif,"significatif",NA),pval,
                           courbe=vpan[1],
                           panel=vpan[1])
        ## netoyage des intervalle de confiance superieur très très grande
        if(assessIC) {
        tab1$UL <- ifelse( nb_carre_presence==0,NA,tab1$UL)
        tab1$UL <-  ifelse(tab1$UL == Inf, NA,tab1$UL)
        tab1$UL <-  ifelse(tab1$UL > 1.000000e+20, NA,tab1$UL)
        tab1$UL[1] <- 1
        tab1$val <-  ifelse(tab1$val > 1.000000e+20,1.000000e+20,tab1$val)
        }
        ## indice de surdispersion
        ## browser()
        if(assessIC) dispAn <- glm1$deviance/glm1$null.deviance else dispAn <- glm1$deviance/glm1$nulldev


        ## tabAn table de sauvegarde des resultats      
        tabAn <- data.frame(id,code_espece=sp, nom_espece = nomSp,indicateur = indic,annee = tab1$annee,
                            abondance_relative=round(tab1$val,3),
                            IC_inferieur = round(tab1$LL,3), IC_superieur = round(tab1$UL,3),
                            erreur_standard = round(erreurannee1,4),
                            p_value = round(tab1$pval,3),significatif = !is.na(tab1$catPoint),
                            nb_carre,nb_carre_presence,abondance=abond)
        
        ## GLM tendance generale sur la periode
        formule <- as.formula(paste("abond~ as.factor(carre) + annee",sep=""))
          #  browser()
    
       
         if(assessIC) {
             md2 <- glm(formule,data=d,family=quasipoisson) }
        else {
                md2 <- try(speedglm(formule,data=d,family=quasipoisson()),silent=TRUE)

                if(class(md2)[1]=="try-error")
                    md2 <- glm(formule,data=d,family=quasipoisson)
            }

        
        smd2 <- summary(md2)
        smd2 <- coefficients(smd2)
        smd2 <- tail(smd2,1)
       
        ## tendences sur la periode
        coefan <- as.numeric(as.character(smd2[,1]))
        trend <- round(exp(coefan),3)
        ## pourcentage de variation sur la periode
        pourcentage <- round((exp(coefan*pasdetemps)-1)*100,2)
        pval <- as.numeric(as.character(smd2[,4]))
        
        erreuran <- as.numeric(as.character(smd2[,2])) 
        ## erreur standard 
        erreurannee2 <- erreuran*exp(coefan)
        
        
        ## calcul des intervalle de confiance
        LL <- NA
        UL <- NA
        if(assessIC) {
            md2.sim <- sim(md2)
            LL <- round(exp(tail(apply(coef(md2.sim), 2, quantile,.025),1)),3)
            UL <- round(exp(tail(apply(coef(md2.sim), 2, quantile,.975),1)),3)
        } else {
            LL <- NA
            UL <- NA
        }
        
        ## tab1t table utile pour la realisation des figures 
        tab1t <- data.frame(Est=trend,
                            LL , UL,
                            pourcent=pourcentage,signif=pval<seuilSignif,pval)
        
        
        trendsignif <- tab1t$signif
        pourcent <- round((exp(coefan*pasdetemps)-1)*100,3)
        ## surdispersion

          if(assessIC) dispTrend <- md2$deviance/md2$null.deviance else dispTrend <- md2$deviance/md2$nulldev


        
        ## classement en categorie incertain
        # browser()
        if(assessIC) {
        if(dispTrend > 2 | dispAn > 2 | median( nb_carre_presence)<seuilOccu) catIncert <- "Incertain" else catIncert <-"bon"
        vecLib <-  NULL
        if(dispTrend > 2 | dispAn > 2 | median( nb_carre_presence)<seuilOccu) {
            if(median( nb_carre_presence)<seuilOccu) {
                vecLib <- c(vecLib,"espece trop rare")
            }
            if(dispTrend > 2 | dispAn > 2) {
                vecLib <- c(vecLib,"deviance")
            }
        }
        raisonIncert <-  paste(vecLib,collapse=" et ")
        } else {
            catIncert <- NA
            raisonIncert <- NA
        }
        
        
        
        ## affectation des tendence EBCC
        catEBCC <- NA
        if(assessIC)  catEBCC <- affectCatEBCC(trend = as.vector(trend),pVal = pval,ICinf=as.vector(LL),ICsup=as.vector(UL)) else catEBCC <- NA
        ## table complete de resultats
        # browser()
        tabTrend <- data.frame(
            id,code_espece=sp,nom_espece = nomSp,indicateur = indic,
            nombre_annees = pasdetemps,premiere_annee = firstY,derniere_annee = lastY,
            tendance = as.vector(trend) ,  IC_inferieur=as.vector(LL) , IC_superieur = as.vector(UL),pourcentage_variation=as.vector(pourcent),
            erreur_standard = as.vector(round(erreurannee2,4)), p_value = round(pval,3),
            significatif = trendsignif,categorie_tendance_EBCC=catEBCC,mediane_occurrence=median( nb_carre_presence) ,
            valide = catIncert,raison_incertitude = raisonIncert)


        if(assessIC)  listGLMsp <- list(list(glm1,glm1.sim,md2,md2.sim)) else  listGLMsp <- list(list(glm1,md2))
        names(listGLMsp)[[1]] <-sp 
        fileSaveGLMsp <- paste(fileSaveGLMs,"_",sp,".Rdata",sep="")
        
        save(listGLMsp,file=fileSaveGLMsp)
        cat("--->",fileSaveGLMsp,"\n")
        flush.console()

        if(sp==listSp[1]) {
            glmAn <- tabAn
            glmTrend <- tabTrend
        } else  {
            glmAn <- rbind(glmAn,tabAn)
            glmTrend <- rbind(glmTrend,tabTrend)
        }
	## les figures     
        if(figure) {
            ## table complete pour la figure en panel par ggplot2
            ## table pour graphe en panel par ggplot2
            if(description)	dgg <- rbind(tab1,tab2,tab3) else dgg <- tab1
            ## les figures     
            ggplot.espece(dgg,tab1t,id,serie=NULL,sp,valide=catIncert,nomSp,description,tendanceSurFigure,seuilOccu=14,vpan = vpan)
            
        }
        
        if(ecritureStepByStep) {

            write.csv2(glmAn,filesaveAn,row.names=FALSE,quote=FALSE)
            cat("--->",filesaveAn,"\n")
            write.csv2(glmTrend,filesaveTrend,row.names=FALSE,quote=FALSE)
            cat("--->",filesaveTrend,"\n")
            
            flush.console()

        }
        
        
    }
    
    write.csv2(glmAn,filesaveAn,row.names=FALSE,quote=FALSE)
    cat("--->",filesaveAn,"\n")
    write.csv2(glmTrend,filesaveTrend,row.names=FALSE,quote=FALSE)
    cat("--->",filesaveTrend,"\n")
    
    flush.console()
  
}




#PLot function called during main.glm execution
ggplot.espece <- function(dgg,tab1t,id,serie=NULL,sp,valide,nomSp=NULL,description=TRUE,
                          tendanceSurFigure=TRUE,seuilOccu=14, vpan) {

    require(ggplot2)
    figname<- paste("Output/",id,"/",ifelse(valide=="Incertain","Incertain/",""),
                    sp,"_",id,serie, ".png",
                    sep = "")
    ## coordonnée des ligne horizontal de seuil pour les abondances et les occurences
    hline.data1 <- data.frame(z = c(1), panel = c(vpan[1]),couleur = "variation abondance",type="variation abondance")
    hline.data2 <- data.frame(z = c(0,seuilOccu), panel = c(vpan[2],vpan[2]),couleur = "seuil",type="seuil")
    hline.data3 <- data.frame(z = 0, panel = vpan[3] ,couleur = "seuil",type="seuil")  
    hline.data <- rbind(hline.data1,hline.data2,hline.data3)
    titre <- paste(nomSp)#,"\n",min(annee)," - ",max(annee),sep="")

    ## texte de la tendence
    tab1 <- subset(dgg,panel =="Variation abondance")
    pasdetemps <- max(dgg$annee) - min(dgg$annee) + 1
    txtPente1 <- paste(tab1t$Est,
                       ifelse(tab1t$signif," *",""),"  [",tab1t$LL," , ",tab1t$UL,"]",
                       ifelse(tab1t$signif,paste("\n",ifelse(tab1t$pourcent>0,"+ ","- "),
                                                 abs(tab1t$pourcent)," % en ",pasdetemps," ans",sep=""),""),sep="")
    ## table du texte de la tendence
    tabTextPent <- data.frame(y=c(max(c(tab1$val,tab1$UL),na.rm=TRUE)*.9),
                              x=median(tab1$annee),
                              txt=ifelse(tendanceSurFigure,c(txtPente1),""),
                              courbe=c(vpan[1]),panel=c(vpan[1]))
    ## les couleurs
    vecColPoint <- c("#ffffff","#eeb40f","#ee0f59")
    names(vecColPoint) <- c("significatif","infSeuil","0")
    vecColCourbe <- c("#3c47e0","#5b754d","#55bb1d","#973ce0")
    names(vecColCourbe) <- c(vpan[1],"carre","presence",vpan[3])
    vecColHline <- c("#ffffff","#e76060")
    names(vecColHline) <- c("variation abondance","seuil")

    col <- c(vecColPoint,vecColCourbe,vecColHline)
    names(col) <- c(names(vecColPoint),names(vecColCourbe),names(vecColHline))

    ## si description graphique en 3 panels
    if(description) {
       
        p <- ggplot(data = dgg, mapping = aes(x = annee, y = val))
        ## Titre, axes ...
        p <- p + facet_grid(panel ~ ., scale = "free") +
            theme(legend.position="none",
                  panel.grid.minor=element_blank(),
                  panel.grid.major.y=element_blank())  +
            ylab("") + xlab("Année")+ ggtitle(titre) +
            scale_colour_manual(values=col, name = "" ,
                                breaks = names(col))+
            scale_x_continuous(breaks=min(dgg$annee):max(dgg$annee))
        p <- p + geom_hline(data =hline.data,mapping = aes(yintercept=z, colour = couleur,linetype=type ),
                            alpha=1,size=1.2)

        p <- p + geom_ribbon(mapping=aes(ymin=LL,ymax=UL),fill=col[vpan[1]],alpha=.2) 
        p <- p + geom_pointrange(mapping= aes(y=val,ymin=LL,ymax=UL),fill=col[vpan[1]],alpha=.2)
        p <- p + geom_line(mapping=aes(colour=courbe),size = 1.5)
        p <- p + geom_point(mapping=aes(colour=courbe),size = 3)
        p <- p + geom_point(mapping=aes(colour=catPoint,alpha=ifelse(!is.na(catPoint),1,0)),size = 2)
        p <-  p + geom_text(data=tabTextPent, mapping=aes(x,y,label=txt),parse=FALSE,color=col[vpan[1]],fontface=2, size=4)
        ggsave(figname, p,width=16,height=21, units="cm")
    } else {

        p <- ggplot(data = subset(dgg,panel=="Variation abondance"), mapping = aes(x = annee, y = val))
        ## Titre, axes ...
        p <- p + facet_grid(panel ~ ., scale = "free") +
            theme(legend.position="none",
                  panel.grid.minor=element_blank(),
                  panel.grid.major.y=element_blank())  +
            ylab("") + xlab("Année")+ ggtitle(titre) +
            scale_colour_manual(values=col, name = "" ,
                                breaks = names(col))+
            scale_x_continuous(breaks=min(dgg$annee):max(dgg$annee))
        p <- p + geom_hline(data =subset(hline.data,panel=="Variation abondance"),mapping = aes(yintercept=z, colour = couleur,linetype=type ),
                            alpha=1,size=1.2)
        
        p <- p + geom_ribbon(mapping=aes(ymin=LL,ymax=UL),fill=col[vpan[1]],alpha=.2) 
        p <- p + geom_pointrange(mapping= aes(y=val,ymin=LL,ymax=UL),fill=col[vpan[1]],alpha=.2)
        p <- p + geom_line(mapping=aes(colour=courbe),size = 1.5)
        p <- p + geom_point(mapping=aes(colour=courbe),size = 3)
        p <- p + geom_point(mapping=aes(colour=catPoint,alpha=ifelse(!is.na(catPoint),1,0)),size = 2)
        p <-  p + geom_text(data=tabTextPent, mapping=aes(x,y,label=txt),parse=FALSE,color=col[vpan[1]],fontface=2, size=4)
        ggsave(figname, p,width=15,height=9,units="cm")
    }
}
