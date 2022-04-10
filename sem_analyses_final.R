################################packages################################
    require(dplyr)
    require(ggplot2)
    require(vegan)
    require(tidyr)
    require(stringr)
    library(ggrepel)
    require(gridExtra)
    library(tibble)
    library(kableExtra)
    library(broom)
    library(scales)


################################palettes################################
    study_color_pallette <- scale_color_manual(values = c("Quercus rubra" = "indianred1",
                                                          "Quercus virginiana" = "navajowhite2",
                                                          "Pinus echinata" = "darkseagreen4" , 
                                                          "Pinus elliottii" = "darkseagreen4",
                                                          "Celtis laevigata" = "lemonchiffon2", 
                                                          "1" = "grey30", 
                                                          "2" = "grey60"))
    study_color_pallette2 <- scale_fill_manual(values = c("Quercus rubra" = "indianred1",
                                                          "Quercus virginiana" = "navajowhite2",
                                                          "Pinus echinata" = "darkseagreen4" , 
                                                          "Pinus elliottii" = "darkseagreen4",
                                                          "Celtis laevigata" = "lemonchiffon2", 
                                                          "1" = "grey30", 
                                                          "2" = "grey60"))

################################create my custom theme################################
    theme_woodfall <- function () { 
      theme_bw(base_size=12) %+replace% 
        theme(
          # change stuff here
          axis.line = element_line(colour = "darkgrey"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.position="none",
          strip.background = element_blank(),
          plot.title = element_text(lineheight=.8, face="bold", hjust = 0))
    }



################################load data################################

    logbytaxa<- data.frame(read.csv("woodfall_logbytaxa.csv", header=TRUE, stringsAsFactors=FALSE))
    
    logsummary <-  data.frame(read.csv("woodfall_logs.csv", header=TRUE, stringsAsFactors=FALSE))
    
    logsummary_WF1 <- logsummary%>%
      filter(site=="WF1")
    

    
    NCOL <- ncol(logbytaxa)
    
    #filter stations and species with zero total abundance
    #filter unidentified species and juveniles
    logbytaxa_clean <- logbytaxa %>%
      replace(is.na(.), 0) %>% #replace all NAs (missing data cells) with zeros
      mutate(rsum = rowSums(.[2:NCOL])) %>% #remove species (rows) with sums of zeros
      filter(rsum>0) %>%
      select(-rsum) %>%
      filter(!ID %in% c("xylo-unid", "gas-juv", "gas-unid", "idas-juv","poly-unid", "poly-juv",
                        "sip-unid","iso-unid","holo-unid","holo-juv", "deca-unid", "gal-UNID", 
                        "poly-larva", "amph-unid"))
    
    
    logbytaxa_clean <- logbytaxa_clean[, colSums(logbytaxa_clean != 0) > 0] #remove stations (cols) with sums of zero
    
    logbytaxa_WF1 <- logbytaxa_clean %>%
      select( -"X.1") %>% 
      select(-contains("L.Y."))
    
    
    #time to transpose for vegan
    # first remember the names
    n <- logbytaxa_WF1$ID
    
    # transpose all but the first column (name)
    logbytaxa_WF1_vegan <- as.data.frame(t(logbytaxa_WF1[,-1]))
    colnames(logbytaxa_WF1_vegan) <- n
    logbytaxa_WF1_vegan$log_id <- factor(row.names(logbytaxa_WF1_vegan))



################################ diversity time################################
    logdiversity_WF1 <- logbytaxa_WF1_vegan %>%
      mutate(Abundance = rowSums(logbytaxa_WF1_vegan[,1:161]),                 
             S = specnumber(logbytaxa_WF1_vegan[,1:161]),
             H = diversity(logbytaxa_WF1_vegan[,1:161],index="shannon"),
             Simp = diversity(logbytaxa_WF1_vegan[,1:161],index="simpson"),
             log10Abundance=log10(Abundance)) %>%
      select(Abundance, log10Abundance, S, H, Simp, log_id)
    
    
    logdiversity_WF1$log_id<-str_replace_all(logdiversity_WF1$log_id, "[.]", "-")
    
    logdiversity_WF1 <- left_join(logdiversity_WF1, logsummary_WF1, by=c("log_id"="Log.ID"))  %>%
      select(Abundance, log10Abundance, S, H, Simp, log_id, wood.type, wood.mass.final, wood.mass.final) %>%
      mutate(log10Mass = log10(wood.mass.final))

################################MULTIVARIATE################################
    logbytaxa_WF1_vegan <-logbytaxa_WF1_vegan  %>% select(-log_id)
    
    #Hellinger pre-transformation of the species matrix
    log.h <- decostand (logbytaxa_WF1_vegan, "hellinger")
    log.h.pca <- rda(log.h )
    log.h.pca.summ <- summary(log.h.pca)
    
    #pull coordinates and loadings
    df1  <- data.frame(log.h.pca.summ$sites[,1:2])
    logdiversity_WF1$PC1 <- df1[,1]
    logdiversity_WF1$PC2 <- df1[,2]
    
    #pull loadings
    log.loadings  <- data.frame(log.h.pca.summ$species[,1:2])
    log.loadings2 <- log.loadings %>%
      rownames_to_column('Species') %>%
      filter(PC1 > quantile(log.loadings$PC1, 0.98) | PC1 < quantile(log.loadings$PC1, 0.02) |
               PC2 > quantile(log.loadings$PC2, 0.98) | PC2 < quantile(log.loadings$PC2, 0.02) )
    
    
    #figure
    pc_plot <- ggplot(data = logdiversity_WF1, aes(x=PC1, y=PC2, label=wood.type)) + 
      geom_point(pch=21, color="black",aes(size=log10(wood.mass.final), fill = wood.type))+
      geom_text_repel(data=log.loadings2, aes(x=PC1, y=PC2, label=Species), color="grey50")+
      #geom_text_repel(aes(label=log_id, color = wood.type),size=4, vjust=-2) +
      geom_hline(yintercept=0, linetype="dotted") +
      geom_vline(xintercept=0, linetype="dotted") +
      study_color_pallette+
      study_color_pallette2+
      theme_woodfall()+
      theme(legend.position="none")+
      ggtitle("A")




    
    ################################body size################################
    
 
    
    speciesmasses <-  data.frame(read.csv("woodfall_taxa_masses.csv", header=TRUE, stringsAsFactors=FALSE))  
    speciesmasses$averagesize=(speciesmasses$Biomass..mg./speciesmasses$Individuals)
    
    mass_summary <-speciesmasses %>%
      group_by(Species.ID) %>%
      summarise(size=mean(averagesize)) 
    
    wood_triplet <- logbytaxa_WF1_vegan %>%
      rownames_to_column(var = "logID") %>%
      gather(species, abundance, -logID) %>%
      filter(abundance > 0)   %>%
      mutate(species = sub("(.)", "\\U\\1", species, perl=TRUE), #capitalize first letter 
             logID = str_replace_all(logID, "[.]", "-"))
    
    wood_triplet <- left_join(wood_triplet, mass_summary, by=c("species"="Species.ID"))
    
    wood_triplet <- left_join(wood_triplet, logdiversity_WF1, by=c("logID"="log_id"))
    #%>%
      #filter(logID!="L-G-34")
    
    
    require(moments)
    
    comm_size_summary <- wood_triplet %>%
      group_by(logID) %>%
      summarise(
        kurtS=kurtosis(log10(size), na.rm = TRUE),
        skewS=skewness(log10(size), na.rm = TRUE),
        meanS=mean(log10(size), na.rm = TRUE),
        medS=median(log10(size), na.rm = TRUE),
        devS=sd(log10(size), na.rm = TRUE)
      )
    
    logdiversity_WF1 <- left_join(comm_size_summary, logdiversity_WF1, by=c("logID"="log_id"))



################################SEM################################


    #individual piece wise models
    logdiversity_WF1$wood.type <- as.factor(logdiversity_WF1$wood.type)
    
    logdiversity_WF1$wood.type.bin <- logdiversity_WF1$wood.type
    levels(logdiversity_WF1$wood.type.bin) <- c("0", "1")
    logdiversity_WF1$wood.type.bin <- as.numeric(logdiversity_WF1$wood.type.bin)
    
    #full model
    
        #richness
        S_mod <- glm(S ~ wood.type.bin + log10Mass + PC1 + PC2 + skewS + log10Abundance +
                       Simp, data=logdiversity_WF1)
        #abundance
        Ab_mod <- glm(log10Abundance ~ wood.type.bin + log10Mass + PC1 + PC2 + skewS,
                      data=logdiversity_WF1)
        #simpson
        Simp_mod <- glm(Simp ~ wood.type.bin + log10Mass + PC1 + PC2 + skewS + log10Abundance,
                        data=logdiversity_WF1)
        #PC1
        PC1_mod <- glm(PC1 ~ wood.type.bin + log10Mass + skewS ,
                       data=logdiversity_WF1)
        #PC2
        PC2_mod <- glm(PC2 ~ wood.type.bin + log10Mass + skewS,
                       data=logdiversity_WF1)
        #Skewed S
        skewS_mod <- glm(skewS ~ wood.type.bin + log10Mass,
                         data=logdiversity_WF1)
        #SEM Model
        library(piecewiseSEM)
        
        #full model
        model_whole <- psem (S_mod, Ab_mod, Simp_mod, PC1_mod, PC2_mod, skewS_mod)
        model_whole
        summary(model_whole)
        plot(model_whole, node_attrs = list(shape = "circle",  color = "grey50", fontsize="5",fillcolor = "grey"))  
    
      #ecological model
        
        #richness
        S_mod_eco <- glm(S ~ PC1 + PC2 + skewS + log10Abundance + Simp, 
                         data=logdiversity_WF1)
        #abundance
        Ab_mod_eco <- glm(log10Abundance ~ wood.type.bin + log10Mass + PC1 + PC2 + skewS,
                      data=logdiversity_WF1)
        #simpson
        Simp_mod_eco <- glm(Simp ~ wood.type.bin + log10Abundance,
                        data=logdiversity_WF1)
        #PC1
        PC1_mod_eco <- glm(PC1 ~ wood.type.bin + log10Mass,
                       data=logdiversity_WF1)
        #PC2
        PC2_mod_eco <- glm(PC2 ~ wood.type.bin + log10Mass,
                       data=logdiversity_WF1)
        #Skewed S
        skewS_mod_eco <- glm(skewS ~ wood.type.bin + log10Mass,
                         data=logdiversity_WF1)
        model_eco <- psem (S_mod_eco, Ab_mod_eco, Simp_mod_eco, PC1_mod_eco, PC2_mod_eco, skewS_mod_eco)
        model_eco
        summary(model_eco)
        
      #reduced model
        
        #richness
        S_mod_red <- glm(S ~ log10Abundance + Simp, data=logdiversity_WF1)
        #abundance
        Ab_mod_red <- glm(log10Abundance ~ log10Mass + PC2,
                          data=logdiversity_WF1)
        #simpson
        Simp_mod_red <- glm(Simp ~ wood.type.bin + log10Abundance + PC1 + PC2 + skewS,
                            data=logdiversity_WF1)
        #PC1
        PC1_mod_red <- glm(PC1 ~ wood.type.bin,
                           data=logdiversity_WF1)
        #PC2
        PC2_mod_red <- glm(PC2 ~ log10Mass,
                           data=logdiversity_WF1)
        #Skewed S
        skewS_mod_red <- glm(skewS ~ wood.type.bin,
                             data=logdiversity_WF1)
        
        #full model
        model_red <- psem (S_mod_red, Ab_mod_red, Simp_mod_red, PC1_mod_red, PC2_mod_red, skewS_mod_red)
        model_red
        summary(model_red)
 
    
    
    
################################New figures################################
require(ggeffects)

    #richness
      ggeS1 <- ggemmeans(S_mod_red2, "Simp [all]" )
      
        s1 <- ggplot(ggeS1, aes(x, predicted)) +
          geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                      color = NA, alpha = 0.2) +
          geom_line() +
          geom_point(data = get_residual_point_data(ggeS1,S_mod_red2),cex=4, pch=21 ,color="black", alpha=0.5) +
          theme_woodfall()+
          xlab("Evenness")+
          ylab("Richness Partial Residials")+
          ggtitle("B")
      
      ggeS2 <- ggemmeans(S_mod_red2, "log10Abundance [all]" )
      
        s2 <- ggplot(ggeS2, aes(x, predicted)) +
          geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                      color = NA, alpha = 0.2) +
          geom_line() +
          geom_point(data = get_residual_point_data(ggeS2,S_mod_red2),cex=4, pch=21,color="black", alpha=0.5) +
          theme_woodfall()+
          xlab("Log10 Abundance")+
          ylab("Richness Partial Residials")+
          ggtitle("A")

      grid.arrange(s2, s1, nrow=1)
      
      pdf(file="figure3.pdf",width=8, height=4)
      par(mar=c(5,3,2,2)+0.1) #removes space from around edges of pdf
      grid.arrange(s2, s1, nrow=1)
      dev.off()
    
      #abundance 
      ggeA1 <- ggemmeans(Ab_mod_red3, "log10Mass [all]" )

        A1 <- ggplot(ggeA1, aes(x, predicted)) +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                    color = NA, alpha = 0.2) +
        geom_line() +
        geom_point(data = get_residual_point_data(ggeA1,Ab_mod_red2),cex=4, pch=21,fill="grey50",color="black", alpha=0.5) +
        theme_woodfall()+
        xlab("Log10 Wood Fall Mass")+
        ylab("Log10 Abundance Partial Residials")+
        ggtitle("A")
        
        
      ggeA3 <- ggemmeans(Ab_mod_red2, "PC2 [all]" )
        
        A3 <- ggplot(ggeA3, aes(x, predicted)) +
          geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                      color = NA, alpha = 0.2) +
          geom_line() +
          geom_point(data = get_residual_point_data(ggeA3,Ab_mod_red2),cex=4, pch=21,fill="grey50",color="black", alpha=0.5) +
          theme_woodfall()+
          xlab("PC2")+
          ylab("Log10 Abundance Partial Residials")+
          ggtitle("B")
        
    
        pdf(file="figure4.pdf",width=8, height=4)
        par(mar=c(5,3,2,2)+0.1) #removes space from around edges of pdf
        grid.arrange(A1,A3, nrow=1)
        dev.off()    
        
      
        
      #evenness
      
      ggeE1 <- ggemmeans(Simp_mod, "log10Abundance [all]" )
      
      E1 <- ggplot(ggeE1, aes(x, predicted)) +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                    color = NA, alpha = 0.2) +
        geom_line() +
        geom_point(data = get_residual_point_data(ggeE1,Simp_mod),cex=4, pch=21,fill="grey50",color="black", alpha=0.5) +
        theme_woodfall()+
        xlab("Log10 Abundance")+
        ylab("Evenness Partial Residials")+
        ggtitle("A")
      
      ggeE2 <- ggemmeans(Simp_mod, "skewS [all]" )
      
      E2 <- ggplot(ggeE2, aes(x, predicted)) +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                    color = NA, alpha = 0.2) +
        geom_line() +
        geom_point(data = get_residual_point_data(ggeE2,Simp_mod),cex=4, pch=21,fill="grey50",color="black", alpha=0.5) +
        theme_woodfall()+
        xlab("Skewness of Body Size")+
        ylab("Evenness Partial Residials")+
        ggtitle("B")
      
      ggeE3 <- ggemmeans(Simp_mod, "PC1 [all]" )
      
      E3 <- ggplot(ggeE3, aes(x, predicted)) +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                    color = NA, alpha = 0.2) +
        geom_line() +
        geom_point(data = get_residual_point_data(ggeE3,Simp_mod),cex=4, pch=21,fill="grey50",color="black", alpha=0.5) +
        theme_woodfall()+
        xlab("PC1")+
        ylab("Evenness Partial Residials")+
        ggtitle("C")
      
      ggeE4 <- ggemmeans(Simp_mod, "PC2 [all]" )
      
      E4 <- ggplot(ggeE4, aes(x, predicted)) +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                    color = NA, alpha = 0.2) +
        geom_line() +
        geom_point(data = get_residual_point_data(ggeE4,Simp_mod),cex=4, pch=21,fill="grey50",color="black", alpha=0.5) +
        theme_woodfall()+
        xlab("PC2")+
        ylab("Evenness Partial Residials")+
        ggtitle("D")
      
      ggeE5 <- ggemmeans(Simp_mod, "wood.type.bin [all]" )
      
      E5 <- ggplot(ggeE5, aes(x, predicted)) +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                    color = NA, alpha = 0.2) +
        geom_line() +
        geom_point(data = get_residual_point_data(ggeE5,Simp_mod),cex=4, pch=21,fill="grey50",color="black", alpha=0.5) +
        theme_woodfall()+
        xlab("Wood Type")+
        ylab("Evenness Partial Residials")+
        ggtitle("E")
      
      pdf(file="figure5.pdf",width=11.5, height=8)
      par(mar=c(5,3,2,2)+0.1) #removes space from around edges of pdf
      grid.arrange(E1,E2, E3, E4, E5, nrow=2)
      dev.off()    
      
  
      
      ##PC1 and PC2
      
      ggePC1.1 <- ggemmeans(PC1_mod_red, "wood.type.bin [all]" )
      
      PC1.1 <- ggplot(ggePC1.1, aes(x, predicted)) +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                    color = NA, alpha = 0.2) +
        geom_line() +
        geom_point(data = get_residual_point_data(ggePC1.1,PC1_mod_red),cex=4, pch=21,fill="grey50",color="black", alpha=0.5) +
        theme_woodfall()+
        xlab("Wood type")+
        ylab("PC1 Partial Residials")+
        ggtitle("B")
      
      
      ggePC2 <- ggemmeans(PC2_mod_red, "log10Mass [all]")
      PC2.1 <- ggplot(ggePC2, aes(x, predicted)) +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                    color = NA, alpha = 0.2) +
        geom_line() +
        geom_point(data = get_residual_point_data(ggePC2,PC2_mod_red),cex=4, pch=21,fill="grey50",color="black", alpha=0.5) +
        theme_woodfall()+
        xlab("Log10 Wood Fall Mass")+
        ylab("PC2 Partial Residials")+
        ggtitle("C")
      
      lay <- rbind(c(1,1),
                   c(2,3))
     
      pdf(file="figure6.pdf",width=8, height=8)
      par(mar=c(5,3,2,2)+0.1) #removes space from around edges of pdf
      grid.arrange(pc_plot, PC1.1, PC2.1, layout_matrix = lay)
      dev.off()    
      
      
      #body size
      
      
      ggeSkew <- ggemmeans(skewS_mod_red, "wood.type.bin [all]" )
      
      Skew1 <- ggplot(ggeSkew, aes(x, predicted)) +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                    color = NA, alpha = 0.2) +
        geom_line() +
        geom_point(data = get_residual_point_data(ggeSkew,skewS_mod_red),cex=4, pch=21,fill="grey50",color="black", alpha=0.5) +
        theme_woodfall()+
        xlab("Wood Type")+
        ylab("Skewness of Body Size Partial Residials")
      
      pdf(file="figure7.pdf",width=4, height=4)
      par(mar=c(5,3,2,2)+0.1) #removes space from around edges of pdf
      Skew1
      dev.off()    
      
#######plot#######
      mp_S <- ggplot(data=logdiversity_WF1,aes(x=log10Mass, y=S))+
        geom_point(cex=4, alpha=0.8, pch=21, color="black", aes(fill = wood.type))+
        study_color_pallette+
        study_color_pallette2+
        theme_woodfall()+
        geom_smooth(method=lm,se=FALSE, color="grey50")+
        xlab("Log10 Wood Fall Mass")
      
      mp_A <- ggplot(data=logdiversity_WF1,aes(x=log10Mass, y=log10Abundance))+
        geom_point(cex=4, alpha=0.8, pch=21, color="black", aes(fill = wood.type))+
        study_color_pallette+
        study_color_pallette2+
        theme_woodfall()+
        geom_smooth(method=lm,se=FALSE, color="grey50")+
        xlab("Log10 Wood Fall Mass")+
        ylab("Log10 Abundance")
      
      mp_SA <- ggplot(data=logdiversity_WF1,aes(y=S, x=log10Abundance))+
        geom_point(cex=4, alpha=0.8, pch=21, color="black", aes(fill = wood.type))+
        study_color_pallette+
        study_color_pallette2+
        theme_woodfall()+
        geom_smooth(method=lm,se=FALSE, color="grey50")+
        ylab("S")+
        xlab("Log10 Abundance")

      
      pdf(file="figure1.pdf",width=8, height=4)
      par(mar=c(5,3,2,2)+0.1) #removes space from around edges of pdf
      grid.arrange(mp_S, mp_A, mp_SA, nrow=1, respect=TRUE)
      dev.off()
      
      mpS_lm <- lm(S~log10Mass, data=logdiversity_WF1)
      mpS_lm <- lm(S~log10Mass+wood.type, data=logdiversity_WF1)
      summary(mpS_lm)
      
      mpA_lm <- lm(log10Abundance~log10Mass, data=logdiversity_WF1)
      mpA_lm <- lm(log10Abundance~log10Mass+wood.type, data=logdiversity_WF1)
      summary(mpA_lm)
     
      mpSA_lm <- lm(S~log10Abundance, data=logdiversity_WF1)
      summary(mpSA_lm)
      