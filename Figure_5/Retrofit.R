library(data.table)
library(ggplot2)

### read in study and filemap ###

# Note that original transfection groups C and D are labeled E anf F (and vice versa) to group homo- and heterodimers together
study = read.csv("./Figure_5/study.csv",
                 stringsAsFactors = F, skipNul = T)

study = as.data.table(study)


filemap = read.csv("./Figure_5/filemap.csv",
                   stringsAsFactors = F, skipNul = T)

filemap = as.data.table(filemap)

### Note all 7 TRnos have AR2 in cols 1-12 and AR6 in cols 13-24; only AR2 data used for manuscript

### read in data by TR and map to study ###

dt0 = NULL

for (i in 1:nrow(filemap)){
  
  sub.name = paste("./Figure_5/data/",
                   "TRno", filemap[i, TR], ".csv", sep = "")
  
  sub.data = read.csv(sub.name, stringsAsFactors = F, skipNul = T)
  sub.data$Destination.Well = paste(substr(sub.data$Well, 1, 1), as.numeric(substr(sub.data$Well, 2, 3)), sep = "")
  
  temp = study
  temp$Destination.Well = temp$wells
  temp$RLU = sub.data[,3][match(temp$Destination.Well, sub.data$Destination.Well)]
  temp$apid = filemap[i, TR]
  temp$replicate = filemap[i, replicate]

  dt0 = rbind(dt0, temp)
}

dt0 <- as.data.table(dt0)
setkey(dt0, chnm)

# Agonist Panel
dt0[, bval := median(RLU[final.uM < 0.03]), by = .(apid,chnm,biogroup)]
dt0[, pval := min(RLU), by = .(apid)]
dt0[, nval := ((RLU-bval)/(pval-bval))*100]

dt0 = dt0[biogroup != "CYP2D6"] #removed to unclutter graph since it is identical to CYP3A4

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

# FLUT with Bgal and CYPs
dt0$spid = paste(dt0$chnm, dt0$biogroup, sep = " + ")
dt0$spid = factor(dt0$spid, levels = c("2-Hydroxyflutamide + NoRNA", "Flutamide + NoRNA", "Flutamide + Bgal", "Flutamide + CYP1A2",
                                       "Flutamide + CYP3A4"))
dt0$logc = log10(dt0$final.uM)
dt0$resp = dt0$nval

dt0 = dt0[chnm != "DMSO",] #remove vehicle wells
#dt0 = dt0[apid != 1801,] #N3 was noisy replicate out of range

sub.dt <- dt0

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
cbbPalette <- c("#000000","#999999","#00aaff","#ff0000","#3d8c40") 

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
  scale_colour_manual(values=cbbPalette, name = "", labels = c("2-Hydroxyflutamide", "Flutamide",
                                                               "Flutamide + Bgal", "Flutamide + CYP1A2", "Flutamide + CYP3A4")) +
  scale_shape_manual(values= c(1,1,1,16,16), name = "", labels = c("2-Hydroxyflutamide", "Flutamide",
                                                               "Flutamide + Bgal", "Flutamide + CYP1A2", "Flutamide + CYP3A4")) +
  coord_cartesian(ylim = c(105,-5)) +
  ylab("% Inhibition") +
  xlab("Concentration (log-10 uM)") +
  ggtitle("") + theme (axis.title=element_text(size=14,face="bold"),
                       plot.title=element_text(size=18,face="bold", hjust = 0.5)) +
  scale_y_reverse()


file.path <- "./Figure_5/Draft_Figures"
ggsave(plot = curves_plot, width = 200, height = 150, device = "png", dpi = 600, units = "mm", filename = "Figure_5.png",
       path = file.path)

write.csv(dat_hill, "./Figure_5/antagonist_fits.csv",
          row.names = F)
