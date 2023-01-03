library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)

### read in study and filemap ###

# Note that original transfection groups C and D are labeled E anf F (and vice versa) to group homo- and heterodimers together
study = read.csv("./Figure_2/study.csv",
                 stringsAsFactors = F, skipNul = T)

study = as.data.table(study)


filemap = read.csv("./Figure_2/filemap.csv",
                   stringsAsFactors = F, skipNul = T)

filemap = as.data.table(filemap)

### Note all 7 TRnos have AR2 in cols 1-12 and AR6 in cols 13-24; only AR2 data used for manuscript

### read in data by TR and map to study ###

dt0 = NULL

for (i in 1:nrow(filemap)){
  
  sub.name = paste("./Figure_2/data/",
                   "TRno", filemap[i, TR], ".csv", sep = "")
  
  sub.data = read.csv(sub.name, stringsAsFactors = F, skipNul = T)
  sub.data$Destination.Well = paste(substr(sub.data$Well, 1, 1), as.numeric(substr(sub.data$Well, 2, 3)), sep = "")
  
  temp = study
  temp$Destination.Well = temp$wells
  temp$RLU = sub.data[,3][match(temp$Destination.Well, sub.data$Destination.Well)]
  temp$apid = filemap[i, TR]
  temp$replicate = filemap[i, replicate]
  temp$medium = filemap[i, medium]
  
  dt0 = rbind(dt0, temp)
}

dt0 <- as.data.table(dt0)
setkey(dt0, chnm)

# Agonist Panel
dt1 = dt0[coli %in% 1:6 & medium == "ESM",]
dt1[, bval := median(RLU[chnm == "DMSO"]), by = apid]
dt1[, nval := log(RLU/bval, 2)]

#plot dose-response data

source("./tcpl_Lite/tcplFit.R")
source("./tcpl_Lite/tcpl_Fit_Lite_Sample_Data.R")
source("./tcpl_Lite/tcplObjCnst.R")
source("./tcpl_Lite/tcplObjHill.R")
source("./tcpl_Lite/tcplObjGnls.R")

hill_curve <- function(hill_tp, hill_ga, hill_gw, lconc){
  return(hill_tp/(1+10^((hill_ga - lconc)*hill_gw)))
}
##########
#Sample data is a data frame with 5 chemicals. Change to data table, then loop over chemical name through tcplFit
#tcplFit needs something called bmad, but it doesn't impact the results, 
#but can limit the number of curves attempted to be fit, so set smaller than you think it really is
bmad <- .001
##########

dt1$spid = dt1$chnm
dt1$logc = log10(dt1$final.uM)
dt1$resp = dt1$nval

dt1 = dt1[chnm != "DMSO",] #remove vehicle wells


  sub.dt <- dt1

  dat_steve <- sub.dt
  setkey(dat_steve, spid)
  dat_steve <- dat_steve[,spid := as.factor(spid)]
  
  #Get hill parameters for each spid
  dat_hill <- data.table(spid = unique(dat_steve[,spid]))
  
  for(chem in dat_hill[,spid]){
    pipefit <- tcplFit_Lite(dat_steve[spid == chem,logc], dat_steve[spid == chem,resp], bmad)
    
    aic.vals <- c(pipefit$hill_aic, pipefit$gnls_aic, pipefit$cnst_aic)
    
    hill_wtp <- ifelse(min(aic.vals) == pipefit$hill_aic, pipefit$hill_tp, 
                       ifelse(min(aic.vals) == pipefit$gnls_aic, pipefit$gnls_tp, 0))
    hill_wga <- ifelse(min(aic.vals) == pipefit$hill_aic, pipefit$hill_ga, 
                       ifelse(min(aic.vals) == pipefit$gnls_aic, pipefit$gnls_ga, 0))
    hill_wgw <- ifelse(min(aic.vals) == pipefit$hill_aic, pipefit$hill_gw, 
                       ifelse(min(aic.vals) == pipefit$gnls_aic, pipefit$gnls_gw, 0))
    
    dat_hill[spid == chem, hill_tp := hill_wtp]
    dat_hill[spid == chem, hill_ga := hill_wga]
    dat_hill[spid == chem, hill_gw := hill_wgw]
  }
  
  # The palette with black:
  cbbPalette <- c("#000000", "#ff0000", "#008800", "#0000ff") 
  
  curves_plot <- ggplot(dat_steve, aes(x=logc, y=resp, color = spid, shape = spid))
  
  #plot the fitted curves, color by spid. This code may need to be modified if you have more colors than in cbbPalette
  k <- 1L
  for(chem in dat_hill[,spid]){
    chemical_fit_hill_tp <- dat_hill[spid == chem, hill_tp]
    chemical_fit_hill_ga <- dat_hill[spid == chem, hill_ga]
    chemical_fit_hill_gw <- dat_hill[spid == chem, hill_gw]
    curves_plot <- curves_plot +
      stat_function(fun = hill_curve, args=list(hill_tp = chemical_fit_hill_tp, 
                                                hill_ga = chemical_fit_hill_ga, 
                                                hill_gw = chemical_fit_hill_gw),
                    color = cbbPalette[k], 
                    size=1) 
    k <- k + 1L
  }
  

  curves_plot <- curves_plot +
    theme_bw() +
    geom_point(size = 2.5) + 
    scale_colour_manual(values=cbbPalette, name = "") +
    scale_shape_manual(values= c(0,1), name = "") +
    ylab("Fold Change (log-2)") +
    xlab("Concentration (log-10 uM)") +
    ggtitle("") + theme (axis.title=element_text(size=14,face="bold"),
                                plot.title=element_text(size=18,face="bold", hjust = 0.5)) +
    geom_segment(aes(x = dat_hill[spid == "R1881", hill_ga], xend = dat_hill[spid == "R1881", hill_ga], 
                     y = 0, yend = dat_hill[spid == "R1881", hill_tp/2]), color = "#000000", linetype = 5) +
    geom_segment(aes(x = -4, xend = dat_hill[spid == "R1881", hill_ga], 
                     y = dat_hill[spid == "R1881", hill_tp/2], yend = dat_hill[spid == "R1881", hill_tp/2]), 
                 color = "#000000", linetype = 5) +
    geom_segment(aes(x = dat_hill[spid == "Testosterone", hill_ga], xend = dat_hill[spid == "Testosterone", hill_ga], 
                     y = 0, yend = dat_hill[spid == "Testosterone", hill_tp/2]), color = "#ff0000", linetype = 5) +
    geom_segment(aes(x = -4, xend = dat_hill[spid == "Testosterone", hill_ga], 
                     y = dat_hill[spid == "Testosterone", hill_tp/2], yend = dat_hill[spid == "Testosterone", hill_tp/2]), 
                 color = "#ff0000", linetype = 5) + 
    geom_segment(aes(x = -2, y = 6.3, xend = -2, yend = 5.8),arrow = arrow(length = unit(0.25, "cm")), size = 1, color = "#0000ff")
  


  Figure_2A = curves_plot

  write.csv(dat_hill, "./Figure_2/agonist_fits.csv",
            row.names = F)


# Antagonist Panel
dt2 = dt0[coli %in% 7:12 & medium != "DMEM",] #this will include cols 7-12 of Expt 13 (TRnos 1787,1793,1801)
dt2[, bval := median(RLU[chnm == "DMSO"]), by = apid]
dt2[, pval := min(RLU), by = .(apid)]
dt2[, nval := ((RLU-bval)/(pval-bval))*100]

dt2$spid = dt2$chnm
dt2$logc = log10(dt2$final.uM)
dt2$resp = dt2$nval

dt2 = dt2[chnm != "DMSO",] #remove vehicle wells

sub.dt <- dt2

dat_steve <- sub.dt
setkey(dat_steve, spid)
dat_steve <- dat_steve[,spid := as.factor(spid)]

#Get hill parameters for each spid
dat_hill <- data.table(spid = unique(dat_steve[,spid]))

for(chem in dat_hill[,spid]){
  pipefit <- tcplFit_Lite(dat_steve[spid == chem,logc], dat_steve[spid == chem,resp], bmad)
  
  aic.vals <- c(pipefit$hill_aic, pipefit$gnls_aic, pipefit$cnst_aic)
  
  hill_wtp <- ifelse(min(aic.vals) == pipefit$hill_aic, pipefit$hill_tp, 
                     ifelse(min(aic.vals) == pipefit$gnls_aic, pipefit$gnls_tp, 0))
  hill_wga <- ifelse(min(aic.vals) == pipefit$hill_aic, pipefit$hill_ga, 
                     ifelse(min(aic.vals) == pipefit$gnls_aic, pipefit$gnls_ga, 0))
  hill_wgw <- ifelse(min(aic.vals) == pipefit$hill_aic, pipefit$hill_gw, 
                     ifelse(min(aic.vals) == pipefit$gnls_aic, pipefit$gnls_gw, 0))
  
  dat_hill[spid == chem, hill_tp := hill_wtp]
  dat_hill[spid == chem, hill_ga := hill_wga]
  dat_hill[spid == chem, hill_gw := hill_wgw]
}

# The palette:
cbbPalette <- c("#ff0000", "#ff9999", "#0000ff", "#55aaff", "#4c9a2a", "#acdf87") 

curves_plot <- ggplot(dat_steve, aes(x=logc, y=resp, color = spid, shape = spid))

#plot the fitted curves, color by spid. This code may need to be modified if you have more colors than in cbbPalette
k <- 1L
for(chem in dat_hill[,spid]){
  chemical_fit_hill_tp <- dat_hill[spid == chem, hill_tp]
  chemical_fit_hill_ga <- dat_hill[spid == chem, hill_ga]
  chemical_fit_hill_gw <- dat_hill[spid == chem, hill_gw]
  curves_plot <- curves_plot +
    stat_function(fun = hill_curve, args=list(hill_tp = chemical_fit_hill_tp, 
                                              hill_ga = chemical_fit_hill_ga, 
                                              hill_gw = chemical_fit_hill_gw),
                  color = cbbPalette[k], 
                  size=1) 
  k <- k + 1L
}


curves_plot <- curves_plot +
  theme_bw() +
  geom_point(size = 2.5) + 
  scale_colour_manual(values=cbbPalette, name = "") +
  scale_shape_manual(values= c(0,1,2,5,6,15), name = "") +
  ylab("% Inhibition") +
  xlab("Concentration (log-10 uM)") +
  ggtitle("") + theme (axis.title=element_text(size=14,face="bold"),
                       plot.title=element_text(size=18,face="bold", hjust = 0.5)) +
  scale_y_reverse()

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 1)),
  colhead = list(fg_params=list(cex = 1)),
  rowhead = list(fg_params=list(cex = 0)))

Figure_2B = curves_plot

write.csv(dat_hill, "./Figure_2/antagonist_fits.csv",
          row.names = F)

Figure_2 = ggarrange(plotlist = list(Figure_2A + theme(legend.position='bottom'), 
                                     Figure_2B + theme(legend.position='bottom')), labels = LETTERS[1:2], ncol = 2, nrow = 1)

ggsave(plot = Figure_2, width = 300, height = 150, device = "png", dpi = 600, units = "mm", filename = "Figure_2.png", 
       path = "./Figure_2/Draft_Figures")

