########################required packages########################
#################################################################


  craigs_packages <- c("dplyr", "ggplot2", "vegan", "tidyr", "stringr", "ggrepel", 
                       "gridExtra", "tibble","plyr","ggpubr","kableExtra",
                       "performance", "sjPlot", "sjmisc", "labdsv", "cluster", 
                       "scatterpie", "dendextend", "wesanderson", "psych")
  
  lapply(craigs_packages, require, character.only = TRUE)

############################palettes#############################
#################################################################

  wood_palette <- c("#73bbd2", "#ea8932","#eb3324", "#469d8a")  

  collection_set_colors<-  c("set1" = "grey70",
                  "set2" = "grey10")
      collection_set_pallette.c <- scale_color_manual(values = collection_set_colors)
      collection_set_pallette.f <- scale_fill_manual(values = collection_set_colors)

  cluster_pal_colors <- c("1"="#73bbd2",
                          "2"="#469d8a",
                          "3"="#ea8932",
                           "4"="#eb3324",
                          "X1"="#73bbd2",
                          "X2"="#469d8a",
                          "X3"="#ea8932",
                          "X4"="#eb3324")
      
      cluster_pallette.c <- scale_color_manual(values = cluster_pal_colors)
      cluster_pallette.f <- scale_fill_manual(values = cluster_pal_colors)
      
  second_stage_colors <- c("1"="#73bbd2",
                           "2"="#ea8932",
                           "3"="#eb3324",
                           "4"="#469d8a",
                           "X1"="#73bbd2",
                           "X2"="#ea8932",
                           "X3"="#eb3324",
                           "X4"= "#469d8a")

  second_stage_pal.c <- scale_color_manual(values = second_stage_colors)
  second_stage_pal.f <- scale_fill_manual(values = second_stage_colors)
  

#############################themes##############################
#################################################################
    theme_craig <- function () { 
      theme_bw(base_size=12, base_family="Helvetica") %+replace% 
        theme(
          # change stuff here
          axis.line = element_line(colour = "darkgrey"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          strip.background = element_blank(),
          axis.title = element_text(face="bold"),
          plot.title = element_text(lineheight=.8, face="bold", hjust = 0),
          legend.position="none")
                                  }

    theme_craig_legend <- function () { 
      theme_bw(base_size=12, base_family="Helvetica") %+replace% 
        theme(
          # change stuff here
          axis.line = element_line(colour = "darkgrey"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          strip.background = element_blank(),
          axis.title = element_text(face="bold"),
          plot.title = element_text(lineheight=.8, face="bold", hjust = 0))
                                  }

############################load data############################
#################################################################                                                 
    setwd("~/Dropbox/Return of the Woodfall/Wood Falls/WF Data")
    #setwd("C:/Users/rdixo/Dropbox/Return of the Woodfall/Wood Falls/WF Data")
    
    logbytaxa<- data.frame(read.csv("woodfall_logbytaxa.csv", 
                                    header=TRUE, stringsAsFactors=FALSE))
    logsummary <-  data.frame(read.csv("woodfall_logs.csv", 
                                       header=TRUE, stringsAsFactors=FALSE))
    
########################Data manipulation########################
#################################################################     
    NCOL <- ncol(logbytaxa)
    
    #removing stations and species with zero total abundance
    #removing unwanted taxa
        logbytaxa <- logbytaxa %>%
          replace(is.na(.), 0) %>% #replace all NAs (missing data cells) with zeros
          mutate(rsum = rowSums(.[2:NCOL])) %>% #remove species (rows) with sums of zeros
          filter(rsum>0) %>%
          select(-rsum) %>%
          filter(!ID %in% c("xylo-unid", "gas-juv", "gas-unid", "idas-juv","poly-unid", "poly-juv",
                        "sip-unid","iso-unid","holo-unid","holo-juv", "deca-unid", "gal-UNID", 
                        "poly-larva", "amph-unid", "fish-1", "wtf-2", "wtf-5", "wtf-6", "wtf-7", "wtf-10",
                        "capit-unid", "mys-1", "biv-unid",  "idas-unid"
                        )) 


    #transpose for vegan
        #first remember the names
           n <- logbytaxa$ID
      
        #transpose all but the first column (name)
            logbytaxa_t <- as.data.frame(t(logbytaxa[,-1]))
            colnames(logbytaxa_t) <- n
            logbytaxa_t$log_id <- factor(row.names(logbytaxa_t))
          
            str(logbytaxa_t) # Check the column types
          
            logbytaxa_t<- logbytaxa_t[grep("L.Y", rownames(logbytaxa_t)), ]
            
            logbytaxa_t  <- logbytaxa_t[, colSums(logbytaxa_t != 0) > 0] 
    
      #diversity time
          NCOL <- ncol(logbytaxa_t)-1
        
          logdiversity_wt <-  logbytaxa_t %>%
            mutate(Abundance = rowSums(logbytaxa_t[,1:NCOL]),                 
                   S = specnumber( logbytaxa_t[,1:NCOL]),
                   H = diversity(logbytaxa_t[,1:NCOL],index="shannon"),
                   Simp = diversity(logbytaxa_t[,1:NCOL],index="simpson"),
                   log10Abundance=log10(Abundance)) %>%
            select(Abundance, log10Abundance, S, H, Simp, log_id)
         
        #joining diversity sheets together
            logdiversity_wt$log_id<-str_replace_all(logdiversity_wt$log_id, "[.]", "-")
            logdiversity_wt <- left_join(logdiversity_wt, logsummary, by=c("log_id"="Log.ID"))  
          
        #creating new variables
            logdiversity_wt <- logdiversity_wt %>%
              mutate (log10Mass = log10(wood.mass.final),
                      size = as.factor(notes),
                      hardness <- as.numeric(revalue(wood.type, 
                                                           c("Celtis laevigata"="3910", 
                                                             "Pinus echinata"="3070", 
                                                             "Pinus elliottii"="3380",
                                                             "Quercus rubra"="5430",
                                                             "Quercus virginiana"="12920",
                                                             "Carya illinoiensis"="8100",
                                                             "Magnolia grandiflora"="4540",     
                                                             "Quercus alba"="5990",
                                                             "Salix nigra"="1920",
                                                             "Taxodium distichum"="2270",
                                                             "Ulmus americana"="2830"))),
                      logHard = log10(hardness),
                      shortnames = revalue(wood.type, 
                                                  c("Celtis laevigata"="", 
                                                    "Pinus echinata"="ShP", 
                                                    "Pinus elliottii"="SlP",
                                                    "Quercus rubra"="RO",
                                                    "Quercus virginiana"="LO",
                                                    "Carya illinoiensis"="Pec",
                                                    "Magnolia grandiflora"="SM",     
                                                    "Quercus alba"="WO",
                                                    "Salix nigra"="BW",
                                                    "Taxodium distichum"="BC",
                                                    "Ulmus americana"="AE")),
                      set.numeric= as.numeric(revalue(collection_set, 
                                                             c("set1"="21", 
                                                               "set2"="31"))))%>%
            select(log_id,wood.type, shortnames,
                    Abundance, log10Abundance, 
                    S, H, Simp,
                    deployment.row, deployment.column, 
                   size, logHard, log10Mass, collection_set, set.numeric, wood.mass.final)
      
         
  
#######################First set of plots########################
#################################################################            
    #plotting with wood size#
    
    

    sizeA <- ggplot(data=logdiversity_wt, aes(y=log10Abundance, x=log10Mass))+
      geom_text_repel(aes(label=shortnames), size=4, max.overlaps = Inf)+
      #geom_point(pch=21, size=4, alpha=1, color="black")+
      xlab("log10 Woodfall Size (kg)")+
      ylab("log10 Abundance")+
      theme_craig()+
              ggtitle("A")
                
    sizeS <- ggplot(data=logdiversity_wt, aes(y=S, x=log10Mass))+
      geom_text_repel(aes(label=shortnames), size=4, max.overlaps = Inf)+
      #geom_point(pch=21, size=4, alpha=1, color="black")+
      xlab("log10 Woodfall Size (kg)")+
      theme_craig()+
      ggtitle("B")
    
    sizeJ <- ggplot(data=logdiversity_wt, aes(y=Simp, x=log10Mass))+
      geom_text_repel(aes(label=shortnames), size=4, max.overlaps = Inf)+
      #geom_point(pch=21, size=4, alpha=1, color="black")+
      xlab("log10 Woodfall Size (kg)")+
      theme_craig()+
      ggtitle("C")
    
  #plotting with hardness#
  
    

    hardA <- ggplot(data=logdiversity_wt, aes(y=log10Abundance, x=logHard))+
      geom_text_repel(aes(label=shortnames), size=4, max.overlaps = Inf)+
      #geom_point(pch=21, size=4, alpha=1, color="black")+
      xlab("log10 Jenka Hardness (N)")+
      ylab("log10 Abundance")+
      theme_craig()+
      ggtitle("D")
    
    hardS <- ggplot(data=logdiversity_wt, aes(y=S, x=logHard))+
      geom_text_repel(aes(label=shortnames), size=4, max.overlaps = Inf)+
      #geom_point(pch=21, size=4, alpha=1, color="black")+
      xlab("log10 Jenka Hardness (N)")+
      theme_craig()+
      ggtitle("E")
    
    hardJ <- ggplot(data=logdiversity_wt, aes(y=Simp, x=logHard))+
      geom_text_repel(aes(label=shortnames), size=4, max.overlaps = Inf)+
      #geom_point(pch=21, size=4, alpha=1, color="black")+
      xlab("log10 Jenka Hardness (N)")+
      theme_craig()+
      ggtitle("F")
    
    
    grid.arrange(sizeA, sizeS, sizeJ,
                 hardA, hardS, hardJ,
                 nrow=2)
  
    #S vs I#
      ggplot(data=logdiversity_wt, aes(y=S, x=log10Abundance))+
        geom_text(aes(label=shortnames,color=logHard))+
        theme_craig()+
        xlab("log10 Abundance")+
        scale_fill_gradientn(colours = wood_palette)+
        scale_colour_gradientn(colours = wood_palette)+
        stat_smooth(method = "lm", formula = y ~ x,se=F, color="grey40")
      
      p3 <- ggplot(data=logdiversity_wt, aes(y=S, x=log10Abundance, 
                                       color=collection_set, group=collection_set, fill=collection_set))+
        geom_text(aes(label=shortnames,color=collection_set))+
        theme_craig()+
        xlab("log10 Abundance")+
        collection_set_pallette.f+
        collection_set_pallette.c+
        stat_smooth(method = "lm", formula = y ~ x,se=F)+
        ggtitle("C")
  
    #collection sets#
    p_set_S<-ggplot(data=logdiversity_wt, aes(x=S, group=collection_set, color=collection_set, fill=collection_set))+
      geom_density(alpha=0.3)+
      theme_craig()+
      collection_set_pallette.c+
      collection_set_pallette.f+
      ggtitle("B")
    
    p_set_A<-ggplot(data=logdiversity_wt, aes(x=log10Abundance, group=collection_set, color=collection_set, fill=collection_set))+
      geom_density(alpha=0.3)+
      theme_craig()+
      xlab("log10 Abundance")+
      collection_set_pallette.c+
      collection_set_pallette.f+
      ggtitle("A")
    
    
    grid.arrange(p_set_A, p_set_S,p3, nrow=1)
######################Univariate Models##########################
#################################################################
  ##S models
    S_model_int3 <- glm(S~(logHard+collection_set+log10Mass+log10Abundance)^3, data=logdiversity_wt)
    S_model_int2 <- glm(S~(logHard+collection_set+log10Mass+log10Abundance)^2, data=logdiversity_wt)
    S_model <- glm(S~logHard+collection_set+log10Mass+log10Abundance, data=logdiversity_wt)
    S_model_sz <- glm(S~log10Mass, data=logdiversity_wt)
    S_model_hardness <- glm(S~logHard, data=logdiversity_wt)
    S_model_set <- glm(S~collection_set, data=logdiversity_wt)
    S_model_abundance <- glm(S~log10Abundance, data=logdiversity_wt)
    S_model_abundance_sz <- glm(S~log10Abundance+log10Mass, data=logdiversity_wt)

    S_model_table<-AIC(S_model_int3, S_model_int2, S_model, S_model_sz, S_model_hardness, S_model_set, S_model_abundance,S_model_abundance_sz )


      row.names(S_model_table) <- c("Three-way interactins", "Two-way interactions", "No interactions", 
                                  "Wood Mass Only", 
                                  "Wood Hardness Only",
                                  "Time Set Only",
                                  "Abundance Only",
                                  "Abundance and Wood Mass")
      
      S_model_table  %>%
        kable(caption = "Table: Summary statistics for models", digits=2) %>%
        kable_classic_2(full_width = F)
      
      summary(S_model)
      r2_nagelkerke(S_model)
      S_model2 <- glm(S~logHard+collection_set+log10Mass, data=logdiversity_wt)
      summary(S_model2)
      r2_nagelkerke(S_model2)
      
      #plot models
      r2 <- plot_model(S_model, type = "pred", terms = c("log10Mass", "log10Abundance", "collection_set"))+
        theme_craig_legend()+
        scale_color_manual(values = c("grey70", "grey50", "grey10"))+
        scale_fill_manual(values = c("grey70", "grey50", "grey10"))

      
##Evenness models
      Simp_model_int3 <- glm(Simp~(logHard+collection_set+log10Mass+log10Abundance)^3, data=logdiversity_wt)
      Simp_model_int2 <- glm(Simp~(logHard+collection_set+log10Mass+log10Abundance)^2, data=logdiversity_wt)
      Simp_model <- glm(Simp~logHard+collection_set+log10Mass+log10Abundance, data=logdiversity_wt)
      Simp_model_sz <- glm(Simp~log10Mass, data=logdiversity_wt)
      Simp_model_hardness <- glm(Simp~logHard, data=logdiversity_wt)
      Simp_model_set <- glm(Simp~collection_set, data=logdiversity_wt)
      Simp_model_abundance <- glm(Simp~log10Abundance, data=logdiversity_wt)
      Simp_model_abundance_sz <- glm(Simp~log10Abundance+log10Mass, data=logdiversity_wt)
        
      Simp_model_table<-AIC(Simp_model_int3, Simp_model_int2, Simp_model, 
                            Simp_model_sz, Simp_model_hardness, Simp_model_set, 
                            Simp_model_abundance,Simp_model_abundance_sz )
        
        
        row.names(Simp_model_table) <- c("Three-way interactins", "Two-way interactions", "No interactions", 
                                    "Wood Mass Only", 
                                    "Wood Hardness Only",
                                    "Time Set Only",
                                    "Abundance Only",
                                    "Abundance and Wood Mass")
        
        Simp_model_table  %>%
        kable(caption = "Table: Summary statistics for models", digits=2) %>%
        kable_classic_2(full_width = F)
        
        summary(Simp_model)
        r2_nagelkerke(Simp_model)
        
        #plot models
        r3 <- plot_model(Simp_model, type = "pred", terms = c("log10Mass", "log10Abundance", "collection_set"))+
        theme_craig_legend()+
        scale_color_manual(values = c("grey70", "grey50", "grey10"))+
        scale_fill_manual(values = c("grey70", "grey50", "grey10"))
  
      
  ##A models
    A_model_int3 <- glm(log10Abundance~(logHard+collection_set+log10Mass)^3, data=logdiversity_wt)
    A_model_int2 <- glm(log10Abundance~(logHard+collection_set+log10Mass)^2, data=logdiversity_wt)
    A_model <- glm(log10Abundance~logHard+collection_set+log10Mass, data=logdiversity_wt)
    A_model_sz <- glm(log10Abundance~log10Mass, data=logdiversity_wt)
    A_model_hardness <- glm(log10Abundance~logHard, data=logdiversity_wt)
    A_model_set <- glm(log10Abundance~collection_set, data=logdiversity_wt)

    
    A_model_table<-AIC(A_model_int3, A_model_int2, A_model, A_model_sz, A_model_hardness, A_model_set)
    
    
    row.names(A_model_table) <- c("Three-way interactions", "Two-way interactions", "No interactions", 
                                  "Wood Mass Only", 
                                  "Wood Hardness Only",
                                  "Time Set Only")
    
    A_model_table  %>%
      kable(caption = "Table: Summary statistics for models", digits=2) %>%
      kable_classic_2(full_width = F)
    
    
  
    summary(A_model)
    r2_nagelkerke(A_model)
  
    #plot models
    r1 <- plot_model(A_model, type = "pred", terms = c("log10Mass", "logHard", "collection_set"))+
      theme_craig_legend()+
      scale_color_manual(values = c("grey70", "grey50", "grey10"))+
      scale_fill_manual(values = c("grey70", "grey50", "grey10"))
    
    grid.arrange(r1, r2, r3, nrow=3)
    
########################Multivariate#############################
#################################################################
    row.names(logbytaxa_t) <- logbytaxa_t$log_id
    logbytaxa_t <-logbytaxa_t %>% select(-log_id)
    logbytaxa_t <-  logbytaxa_t[, colSums(logbytaxa_t != 0) > 0] 
    
    
    dist <- logdiversity_wt %>%
      select(deployment.row, deployment.column)
    
    #rename for convenience
    spe <- logbytaxa_t 
    
    #environmental
    env<-logdiversity_wt %>%
      select(log10Mass, logHard, set.numeric)
    
    
    
    # Species indicator values (Dufrene and Legendre)
    # ***********************************************
    # Divide the sites into 4 groups depending 
    #hardness(2), wood.mass.final (1), collection set (3)
    factor = env[,1]
    
    
    das.D1 <- dist(data.frame(das=factor, row.names=logdiversity_wt$log_id))
    dasD1.kmeans <- kmeans(das.D1, centers=2, nstart=100)
    dasD1.kmeans$cluster
    # Indicator species for this typology of the sites
    (iva <- indval(spe, dasD1.kmeans$cluster)) 
    
    boxplot(factor~dasD1.kmeans$cluster)
    
    # Table of the significant indicator species
    gr <- iva$maxcls[iva$pval <= 0.05]
    iv <- iva$indcls[iva$pval <= 0.05]
    pv <- iva$pval[iva$pval <= 0.05]
    
    rf1 <- iva$relfrq[iva$pval <= 0.05,1]
    rf2 <- iva$relfrq[iva$pval <= 0.05,2]
    
    ra1 <- iva$relabu[iva$pval <= 0.05,1]
    ra2 <- iva$relabu[iva$pval <= 0.05,2]
    
    fr <- apply(spe > 0, 2, sum)[iva$pval <= 0.05]
    fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr, relfrq1=rf1, relfrq2=rf2, relabu1=ra1, relabu2=ra2)
    fidg <- fidg[order(fidg$group, -fidg$indval),]
    fidg
  
    
    
    # Fuzzy c-means clustering of the species data
    # *************************************************

    spe.norm <- decostand(spe, "normalize") 
    spe.ch <- vegdist (spe.norm, "euc")
    
    k <-4	# Choose the number of clusters
    spe.fuz <- fanny(spe.ch, k=k, memb.exp=1.5)
    summary(spe.fuz)
    
    # Site fuzzy membership
    spe.fuz$membership
    # Nearest crisp clustering
    spe.fuz$clustering
    spefuz.g <- spe.fuz$clustering
    
    # Silhouette plot
    plot(silhouette(spe.fuz), main="Silhouette plot - Fuzzy clustering", 
         cex.names=0.8,col = c("#73bbd2", "#469d8a", "#ea8932", "#eb3324"))
    

    
    # Ordination of fuzzy clusters (PCoA)
    dc.pcoa <- cmdscale(spe.ch, eig=TRUE)
    dc.scores <- scores(dc.pcoa, choices=c(1,2))

    #get species scores
    species.scores <- data.frame(wascores(dc.pcoa$points[,1:2], spe))
      #limit to upper quantiles
    species.scores2 <- species.scores %>%
      filter(
          X1<= quantile (species.scores$X1, probs=.10) | 
          X1>= quantile (species.scores$X1, probs=.90) | 
          X2<= quantile (species.scores$X2, probs=.10) | 
          X2>= quantile (species.scores$X2, probs=.90) )  %>%
      rownames_to_column('name')
    
    species.scores3 <- species.scores %>%
      filter(row.names(species.scores) %in% c("xylo-5","barn-1","xylo-1","nem-1","platy-1", 
                             "xylo-7", "xylo-12", "xylo-2", "poly-57", "poly-11", 
                             "gal-1", "poly-51", "gas-1", "gas-6 ", "ast-2", "poly-102", 
                             "poly-16", "amph-1", "holo-2", "poly-6", "poly-101", "poly-26",  
                             "xylo-15", "deca-1", "sib-1"))   %>%
      rownames_to_column('name')

    
    
    #new fuzzy cluster plot

    logdiversity_wt$DIM1 <-scores(dc.pcoa)[,1]
    logdiversity_wt$DIM2 <-scores(dc.pcoa)[,2]
    
    
    #create memberships dataframe for pie chats
    membership <- data.frame(spe.fuz$membership)
    membership$DIM1 <-scores(dc.pcoa)[,1]
    membership$DIM2 <-scores(dc.pcoa)[,2]  
    membership$size <- logdiversity_wt$log10Mass
    membership$hard <- logdiversity_wt$logHard
    membership$set  <- logdiversity_wt$collection_set
    membership$cluster <- spe.fuz$clustering
    
    #create centroids for wood types
    wt_centroids <- logdiversity_wt %>%
      dplyr::group_by(shortnames) %>%
      dplyr::summarise(xDIM1=mean(DIM1),
                xDIM2=mean(DIM2))
    
    
    #create hulls for clusters
    hull_data <- 
      membership %>%
      group_by(cluster) %>% 
      slice(chull(DIM1,DIM2)) %>%
      dplyr::select(cluster, DIM1, DIM2)
  
    

    #environmental fits
    vpcoa.env <- envfit(dc.pcoa, env)
    vec.sp.df<-as.data.frame(vpcoa.env$vectors$arrows*sqrt(vpcoa.env$vectors$r))
    vec.sp.df$env<-c("Mass", "Hardness", "Time")
    
    #And Finally the Plot Together
    ggplot(data=logdiversity_wt, aes(DIM1, DIM2))+
      geom_scatterpie(aes(x=DIM1, y=DIM2, r=.03), data=membership, cols=c("X1", "X2", "X3", "X4"))+
      coord_equal()+
      geom_text_repel(aes(label=shortnames), point.size = 10, force=6)+
      geom_text_repel(data=species.scores2, 
                      aes(x=X1, y=X2, label=name), point.size = 10, force=6, color="red4")+
      geom_polygon(data = hull_data,
                   aes(fill = as.factor(cluster),
                       color = as.factor(cluster)),
                   alpha = 0.3,
                   show.legend = FALSE)+
      #geom_text_repel(data=wt_centroids, aes(x=xDIM1, y=xDIM2, label=wood.type), color="blue", alpha=0.7, size=7)+
      geom_segment(data=vec.sp.df,aes(x=0,xend=Dim1,y=0,yend=Dim2),
                   arrow = arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE)+
      geom_text(data=vec.sp.df,aes(x=Dim1,y=Dim2,label=env),size=5)+
      cluster_pallette.c+
      cluster_pallette.f+
      theme_craig()
    
    
    #Broken Up
    PCOA1 <- ggplot(data=logdiversity_wt, aes(DIM1, DIM2))+
      geom_scatterpie(aes(x=DIM1, y=DIM2, r=.06), data=membership, cols=c("X1", "X2", "X3", "X4"))+
      coord_equal()+
      #geom_text_repel(aes(label=shortnames), point.size = 10, force=6,  max.overlaps = Inf)+
      geom_polygon(data = hull_data,
                   aes(fill = as.factor(cluster),
                       color = as.factor(cluster)),
                   alpha = 0.3,
                   show.legend = FALSE)+
      cluster_pallette.c+
      cluster_pallette.f+
      ggtitle("A. Fuzzy Custers")+
      theme_craig()
    
    
    PCOA2 <-   ggplot(data=membership, aes(x=DIM1, y=DIM2))+
      geom_point(pch=21, alpha=.8, size=3,aes(fill = as.factor(cluster)))+
      coord_equal()+
      geom_segment(data=vec.sp.df,aes(x=0,xend=Dim1,y=0,yend=Dim2),arrow = arrow(length = unit(0.5, "cm")),colour="grey20")+
      geom_text(data=vec.sp.df,aes(x=Dim1,y=Dim2,label=env),size=5)+
      cluster_pallette.c+
      cluster_pallette.f+
      ggtitle("B. Factor Loadings")+
      theme_craig()
    
    PCOA3 <-  ggplot(data=membership, aes(x=DIM1, y=DIM2))+
      geom_point(pch=21, alpha=.3, size=3,aes(fill = as.factor(cluster)))+
      coord_equal()+
      geom_text_repel(data=species.scores3, 
                      aes(x=X1, y=X2, label=name), color="red4", max.overlaps = Inf)+
      cluster_pallette.c+
      cluster_pallette.f+
      ggtitle("C. Species Loadings")+
      theme_craig()
    
    PCOA4 <-  ggplot(data=membership, aes(x=DIM1, y=DIM2))+
      geom_point(pch=21, alpha=.3, size=3,aes(fill = as.factor(cluster)))+
      geom_text_repel(data=wt_centroids, aes(x=xDIM1, y=xDIM2, label=shortnames), color="grey10", size=4)+
      coord_equal()+
      cluster_pallette.c+
      cluster_pallette.f+
      ggtitle("D. Wood Type Centroids")+
      theme_craig()
    
    PCOA5 <- ggplot(data=membership, aes(x=DIM1, y=DIM2))+
      geom_point(pch=21, alpha=.8, aes(fill = as.factor(cluster), size=logdiversity_wt$size))+
      coord_equal()+
      ggtitle("log10 Mass")+
      cluster_pallette.c+
      cluster_pallette.f+
      theme_craig()
    
    PCOA6 <-   ggplot(data=membership, aes(x=DIM1, y=DIM2))+
      geom_point(pch=21, alpha=.8, aes(fill = as.factor(cluster), size=hard))+
      coord_equal()+
      cluster_pallette.c+
      cluster_pallette.f+
      ggtitle("log10 Hardness")+
      theme_craig()
    
    PCOA7 <- ggplot(data=membership, aes(x=DIM1, y=DIM2))+
      geom_point(pch=21, alpha=.8, aes(fill = as.factor(cluster), size=set))+
      coord_equal()+
      cluster_pallette.c+
      cluster_pallette.f+
      ggtitle("Collection Set")+
      theme_craig()
    
    grid.arrange(PCOA1, PCOA3, PCOA4, nrow=1)
    grid.arrange(PCOA2, PCOA5, PCOA6, PCOA7, nrow=2)
    
    #plot for paper
    grid.arrange(PCOA1, PCOA2, PCOA3, PCOA4, nrow=2)
    
    
    
    ggplot(data=logdiversity_wt, aes(DIM1, DIM2))+
      geom_scatterpie(aes(x=DIM1, y=DIM2, r=.06), data=membership, cols=c("X1", "X2", "X3", "X4"))+
      coord_equal()+
      geom_polygon(data = hull_data,
                   aes(fill = as.factor(cluster),
                       color = as.factor(cluster)),
                   alpha = 0.3,
                   show.legend = FALSE)+
      cluster_pallette.c+
      cluster_pallette.f+
      ggtitle("Fuzzy Custers")+
      theme_craig()
    
    
#######################MANOVA for Clusters#######################
#################################################################
    clusters.man <- manova(cbind(X1, X2, X3, X4)
                                 ~size+hard+set, data = membership)
                      
    summary(clusters.man, tol=0)
    summary.aov(clusters.man)
    
    x1_hard<-ggplot(data=membership, aes(y=X1, x=hard))+
      geom_point(pch=21, fill="#73bbd2", color="black", alpha=.6, size=3)+
      theme_craig()+
      geom_smooth(method=lm, se=FALSE, color="#73bbd2")+
      ylab("Representation in Cluster 1")+
      xlab("Log10 Jenka Hardness")+
      ggtitle("C")
    
    x2_hard<-ggplot(data=membership, aes(y=X2, x=hard))+
      geom_point(pch=21, fill="#469d8a", color="black", alpha=.6, size=3)+
      theme_craig()+
      geom_smooth(method=lm, se=FALSE, color="#469d8a")+
      ylab("Representation in Cluster 2")+
      xlab("Log10 Jenka Hardness")+
      ggtitle("A")
    
    x3_hard<-ggplot(data=membership, aes(y=X3, x=hard))+
      geom_point(pch=21, fill="#ea8932", color="black", alpha=.6, size=3)+
      theme_craig()+
      geom_smooth(method=lm, se=FALSE, color="#ea8932")+
      ylab("Representation in Cluster 3")+
      xlab("Log10 Jenka Hardness")+
      ggtitle("B")
    
    x3_size<-ggplot(data=membership, aes(y=X3, x=size))+
      geom_point(pch=21, fill="#ea8932", color="black", alpha=.6, size=3)+
      theme_craig()+
      geom_smooth(method=lm, se=FALSE, color="#ea8932")+
      ylab("Representation in Cluster 3")+
      xlab("Log10 Woodfall Size (kg)")+
      ggtitle("D")
    
    x2_set<-ggplot(data=membership, aes(y=X2, x=set))+
      geom_boxplot(pch=21, fill="#469d8a", color="black", alpha=.6)+
      theme_craig()+
      ylab("Representation in Cluster 2")+
      xlab("Set")+
      ggtitle("E")
    
    x4_set<-ggplot(data=membership, aes(y=X4, x=as.factor(set)))+
      geom_boxplot(pch=21, fill="#eb3324", color="black", alpha=.6)+
      theme_craig()+
      ylab("Representation in Cluster 4")+
      xlab("Set")+
      ggtitle("F")
    
    grid.arrange(x2_hard,x3_hard,x1_hard,
                 x3_size, x2_set, x4_set, nrow=2)
    
    
####2nd Stage Analyses####
    spe.poly <- spe %>%
      select(contains('poly'))
    
    spe.xylo <- spe %>%
      select(contains('xylo'))
    
    spe.mol <- spe %>%
      select(contains(c('biv', 'gas' ,'idas' ,'limp' , 'polyplac' , 'aplac')))
    
    spe.arth <- spe %>%
      select(contains(c('tan', 'amph' ,'deca' ,'iso' , 'gal' , 'cum', 'barn')))
    
    spe.cnid <- spe %>%
      select(contains(c('hyd', 'anem')))
    
    spe.ech <- spe %>%
      select(contains(c('ophi', 'ast' ,'holo')))
    
    spe.oth <- spe %>%
      select(contains(c('oli', 'sip' ,'sib', 'nem', 'platy')))
    
    spe.bin <- spe
      spe.bin[spe.bin > 0] = 1
      
    
  
    hard = env[,2]
    mass = env[,1]
    time = env[,3]
    
    
    emp_matrices <- c("env",
                      "hard", "mass", "time", 
                        "spe.poly", "spe.xylo", "spe.mol",
                      "spe.arth", "spe.cnid", "spe.ech", 
                      "spe.oth", "spe.bin")
    
      
    
    emp_distances <- c("dist.env",
                          "dist.hard", "dist.mass", "dist.time", 
                          "dist.spe.poly", "dist.spe.xylo", "dist.spe.mol",
                          "dist.spe.arth", "dist.spe.cnid",  "dist.spe.ech",
                          "dist.spe.oth", "dist.spe.bin")   
    
    emp_labels <- c("env",
                    "hard", "mass", "time", 
                    "poly", "xylo", "mol",
                    "arth", "cnid", "ech", 
                    "oth", "bin")
    

    
    for (i in emp_matrices) {
        assign(paste0("dist.", i), dist(get(i), method="euclidean"))
             }    
    
    
    for(k in emp_distances){
      
      for(j in emp_distances){
        
        assign(paste0(k, ".", j), mantel(get(k), get(j)))
        
      }
      
    }


  second.stage.cor = matrix(nrow = length(emp_distances), ncol = length(emp_distances))
  colnames(second.stage.cor) = emp_distances
  rownames(second.stage.cor) = emp_distances
  
  
  for(m in emp_distances){
    for(n in emp_distances){
          pair = paste0(m, ".", n)
          attach(get(pair))
          second.stage.cor[m,n] = statistic
          detach(get(pair))
      
    }
      }
  

  # Fuzzy c-means clustering of 2nd Stage
  # *************************************************
  rownames(second.stage.cor) <- emp_labels
  
  colnames(second.stage.cor) <- emp_labels
  
  second.stage.dist <- cor2dist(second.stage.cor)
 

  
  k <-4	# Choose the number of clusters
  second.stage.fuz <- fanny(second.stage.dist, k=k, memb.exp=1.5)
  summary(second.stage.fuz)
  
  
  # Silhouette plot
  plot(silhouette(second.stage.fuz), main="Silhouette plot - Fuzzy clustering", 
       cex.names=0.8, col=c("#73bbd2", "#ea8932","#eb3324", "#469d8a"))  
  
  
  
  # Ordination of fuzzy clusters (PCoA)
  dc.pcoa <- cmdscale(second.stage.dist, eig=TRUE)
  dc.scores <- scores(dc.pcoa, choices=c(1,2))
  
  
  #create memberships dataframe for pie chats
  membership <- data.frame(second.stage.fuz$membership)
  membership$DIM1 <-scores(dc.pcoa)[,1]
  membership$DIM2 <-scores(dc.pcoa)[,2]  
  membership$cluster <- second.stage.fuz$clustering
  membership$group <- emp_labels
  
  
  
  #create hulls for clusters
  hull_data <- 
    membership %>%
    group_by(cluster) %>% 
    slice(chull(DIM1,DIM2)) %>%
    dplyr::select(cluster, DIM1, DIM2)
  
  #And Finally the Plot Together
  ggplot(data=membership, aes(DIM1, DIM2))+
    geom_scatterpie(aes(x=DIM1, y=DIM2, r=.03), data=membership, cols=c("X1", "X2", "X3", "X4"))+
    coord_equal()+
    geom_text_repel(aes(label=group), point.size = 10, force=6)+
    geom_polygon(data = hull_data,
                 aes(fill = as.factor(cluster),
                     color = as.factor(cluster)),
                 alpha = 0.3,
                 show.legend = FALSE)+
    second_stage_pal.c+
    second_stage_pal.f+
    theme_craig()
  

  

  
  
    
    
    

#######################Nestedness########################
################################################################# 
  
  require(betapart)
  

  
  #analyses
  logbytaxa_bin <-  logbytaxa_t %>%
    mutate_if(is.numeric, ~1 * (. != 0)) 
  
  wt.betapart<-betapart.core(logbytaxa_bin)
  wt.betapart.dist1<-beta.pair(wt.betapart, index.family="sor")
  
  wt.beta.SIM1<-wt.betapart.dist1$beta.sim
  wt.beta.SNE1<-wt.betapart.dist1$beta.sne
  wt.beta.SOR1<-wt.betapart.dist1$beta.sor
  

  
  #get everything in the same data.frame
  
  wt.betaS1<-NULL
  wt.betaS1$SIM<-as.vector(wt.beta.SIM1)
  wt.betaS1$SNE<-as.vector(wt.beta.SNE1)
  wt.betaS1$SOR<-as.vector(wt.beta.SOR1)
  wt.betaS1$HARD<-as.vector(dist.hard)
  wt.betaS1$MASS<-as.vector(dist.mass)
  wt.betaS1<-data.frame(wt.betaS1)

  
  #plots
  p1<-ggplot(wt.betaS1, aes(HARD, SIM)) +
    geom_point(size=2, alpha=0.7, pch=21, color="black", fill="grey")+
    geom_smooth(method=lm,se=FALSE, color="black")+
    ylab("Replacement")+
    xlab("Difference in Wood Fall Hardness")+
    theme_craig() +
    ggtitle("A")
  
  
  p2<-ggplot(wt.betaS1, aes(HARD, SNE)) +
    geom_point(size=2, alpha=0.7, pch=21, color="black", fill="grey")+
    geom_smooth(method=lm,se=FALSE, color="black")+
    ylab("Nestedness")+
    xlab("Difference in Wood Fall Hardness")+
    theme_craig() +
    ggtitle("C")
  
  
  p3<-ggplot(wt.betaS1, aes(HARD, SOR)) +
    geom_point(size=2, alpha=0.7, pch=21, color="black", fill="grey")+
    geom_smooth(method=lm,se=FALSE, color="black")+
    ylab("Overall")+
    xlab("Difference in Wood Fall Hardness")+
    theme_craig() +
    ggtitle("E")
  
  
  p4<-ggplot(wt.betaS1, aes(MASS, SIM)) +
    geom_point(size=2, alpha=0.7, pch=21, color="black", fill="grey")+
    geom_smooth(method=lm,se=FALSE, color="black")+
    ylab("Replacement")+
    xlab("Difference in Wood Fall Size")+
    theme_craig()+
    ggtitle("B")
  
  
  p5<-ggplot(wt.betaS1, aes(MASS, SNE)) +
    geom_point(size=2, alpha=0.7, pch=21, color="black", fill="grey")+
    geom_smooth(method=lm,se=FALSE, color="black")+
    ylab("Nestedness")+
    xlab("Difference in Wood Fall Size")+
    theme_craig()+
    ggtitle("D")
  
  
  p6<-ggplot(wt.betaS1, aes(MASS, SOR)) +
    geom_point(size=2, alpha=0.7, pch=21, color="black", fill="grey")+
    geom_smooth(method=lm,se=FALSE, color="black")+
    ylab("Overall")+
    xlab("Difference in Wood Fall Size")+
    theme_craig()+
    ggtitle("F")
  
  grid.arrange(p1,p4,
               p2,p5,
               p3,p6,
               nrow=3)
  #hardness
    #replacment
    mantel(dist.hard,wt.beta.SIM1, method="pearson", permutations=9999, na.rm = FALSE)
    
    #nestedness
    mantel(dist.hard, wt.beta.SNE1, method="pearson", permutations=9999, na.rm = FALSE)
    
    #overall
    mantel(dist.hard,wt.beta.SOR1, method="pearson", permutations=9999, na.rm = FALSE)
    
  #mass
    #replacment
    mantel(dist.mass,wt.beta.SIM1, method="pearson", permutations=9999, na.rm = FALSE)
    
    #nestedness
    mantel(dist.mass, wt.beta.SNE1, method="pearson", permutations=9999, na.rm = FALSE)
    
    #overall
    mantel(dist.mass,wt.beta.SOR1, method="pearson", permutations=9999, na.rm = FALSE)
    

    
      
      
    
  