library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)
library(reshape2)

assay.vol.ul = 40

### read in study and filemap ###
study = read.csv("./Figure_3/study.csv",
                 stringsAsFactors = F, skipNul = T)

study = as.data.table(study)


filemap = read.csv("./Figure_3/filemap.csv",
                   stringsAsFactors = F, skipNul = T)

filemap = as.data.table(filemap)


### read in data by TR and map to study ###

dt0 = NULL

for (i in 1:nrow(filemap)){
  
  sub.name = paste("./Figure_3/data/",
                   "TRno", filemap[i, TR], ".csv", sep = "")
  
  sub.data = read.csv(sub.name, stringsAsFactors = F, skipNul = T)
  sub.data$Destination.Well = paste(substr(sub.data$Well, 1, 1), as.numeric(substr(sub.data$Well, 2, 3)), sep = "")
  
  temp = study[Destination.Plate.Barcode == paste("apid_", filemap[i, apid], sep = ""),]
  temp$rowi = substr(temp$Destination.Well, 1, 1)
  temp$coli = substr(temp$Destination.Well, 2, nchar(temp$Destination.Well))
  temp$RLU = sub.data[,3][match(temp$Destination.Well, sub.data$Destination.Well)]
  temp$apid = filemap[i, TR]
  temp$assay = filemap[i, assay]
  temp$replicate = filemap[i, replicate]

  dt0 = rbind(dt0, temp)
}

dt0[, logc := log10(Transfer.Volume*stock.mM/(assay.vol.ul+Transfer.Volume/1000))]
setkey(dt0, chnm)

#remove affected wells from Echo exceptions reports
exceptions1 = read.csv("./Figure_3/validation4.n1.exceptionsreport.csv",
                       stringsAsFactors = F, skipNul = T)
exceptions1$replicate = 1

exceptions2 = read.csv("./Figure_3/validation4.n2.exceptionsreport.csv",
                       stringsAsFactors = F, skipNul = T)
exceptions2$replicate = 2

exceptions3 = read.csv("./Figure_3/validation4.n3.exceptionsreport.csv",
                       stringsAsFactors = F, skipNul = T)
exceptions3$replicate = 3

exceptions = as.data.table(rbind(exceptions1, exceptions2, exceptions3))
exceptions.summ = exceptions[, unique(Source.Well), by = .(Source.Well, Destination.Plate.Barcode, replicate)]

delete.rows = NULL
for (n in 1:nrow(exceptions.summ)){
  
  temp.row = exceptions.summ[n,]
  temp.rows = which(dt0$Source.Plate.Barcode == "DPID_128" & dt0$Source.Well == temp.row[,Source.Well] & 
                      dt0$Destination.Plate.Barcode == temp.row[,Destination.Plate.Barcode] & dt0$replicate == temp.row[,replicate])
  delete.rows = c(delete.rows, temp.rows)
}

dt0 = dt0[-delete.rows,]



#Load tcplFit_Lite

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


# Agonist Mode
dt1 = dt0[assay == "AR_agonist",]
dt1[, bval := median(RLU[coli %in% c(12,13,22,23)]), by = apid] #bval using two lowest test concentrations by assay plate
dt1[, pval := median(RLU[rowi == "B" & coli %in% c(9:23)]), by = apid]
dt1[, nval := log(RLU/bval, 2)]
dt1[, resp := nval]
dt1[, spid := chnm]

bmad1 = dt1[coli %in% c(12,13,22,23), mad(nval)]

dt1 = dt1[!chnm %in% c("DMSO", ""),] # remove vehicle and empty wells

agonist.fits = NULL

for (i in 1:length(unique(dt1$chnm))){

chem.title = unique(dt1$chnm)[i]
sub.dt = dt1[chnm == chem.title,]
sub.med = sub.dt[, .(max_med = median(resp)), by = logc]

dat_steve <- sub.dt
setkey(dat_steve, spid)
dat_steve <- dat_steve[,spid := as.factor(spid)]

#Get hill parameters for each spid
dat_hill <- data.table(spid = unique(dat_steve[,spid]))

for(chem in dat_hill[,spid]){
  pipefit = tcplFit_Lite(dat_steve[spid == chem,logc], dat_steve[spid == chem,resp], bmad)
  
  aic.vals = c(pipefit$hill_aic, pipefit$gnls_aic, pipefit$cnst_aic)
  
  win.model = ifelse(min(aic.vals) == pipefit$hill_aic, "Hill", 
                     ifelse(min(aic.vals) == pipefit$gnls_aic, "Gain-Loss", "Constant"))
  
  hill_wtp <- ifelse(win.model == "Hill", pipefit$hill_tp, ifelse(win.model =="Gain-Loss", pipefit$gnls_tp, 0))
  hill_wga <- ifelse(win.model == "Hill", pipefit$hill_ga, ifelse(win.model =="Gain-Loss", pipefit$gnls_ga, 0))
  hill_wgw <- ifelse(win.model == "Hill", pipefit$hill_gw, ifelse(win.model =="Gain-Loss", pipefit$gnls_gw, 0))
  
  dat_hill[spid == chem, win_model := win.model]
  dat_hill[spid == chem, hill_tp := hill_wtp]
  dat_hill[spid == chem, hill_ga := hill_wga]
  dat_hill[spid == chem, hill_gw := hill_wgw]
  dat_hill[spid == chem, max_med := max(sub.med[, max_med], na.rm = T)]
}

# The palette:
cbbPalette <- c("#0000ff") 

curves_plot <- ggplot(dat_steve, aes(x=logc, y=resp, color = spid))

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
  coord_cartesian(ylim = c(0,max(dt1$resp))) +
  scale_colour_manual(values=cbbPalette, guide = "none") +
  ylab("Fold Change (log-2)") +
  xlab("Concentration (log-10 uM)") +
  geom_hline(yintercept = 5*bmad1, color = "#ff0000", linetype = 5) +
  ggtitle(chem.title) + theme (axis.title=element_text(size=14,face="bold"),
                       plot.title=element_text(size=18,face="bold", hjust = 0.5))


file.name <- paste("chem", i, "png", sep = ".")
file.path <- "./Figure_3/plots/AR_agonist"
ggsave(plot = curves_plot, filename = file.name, 
       path = file.path)

agonist.fits = rbind(agonist.fits, dat_hill)

}

write.csv(agonist.fits, "./Figure_3/agonist_fits.csv",
          row.names = F)

agonist.hits = agonist.fits[win_model != "Constant" & max_med > 5*bmad1 & spid != "R1881",][order(max_med, decreasing = T),spid]

agonist.hitc = agonist.fits
agonist.hitc$agonist_ga = ifelse(agonist.hitc$hill_ga == 0, 2.5, agonist.hitc$hill_ga)
agonist.hitc$agonist_ga = replace(agonist.hitc$agonist_ga, is.na(agonist.hitc$agonist_ga), 2.5)
agonist.hitc[, hitc := ifelse(win_model != "Constant" & max_med > 5*bmad1 & spid != "R1881", 1, 0)]
agonist.hitc$agonist_ga = replace(agonist.hitc$agonist_ga, agonist.hitc$hitc == 0, 2.5)


write.csv(agonist.hitc, "./Figure_3/agonist.hitc.csv",
          row.names = F)

# Figure 3A
target.chems = c("R1881", "17beta-Trenbolone", "Cyproterone acetate", "Bicalutamide", "5alpha-Dihydrotestosterone", "Hydroxyflutamide", 
                 "Spironolactone", "Bisphenol A")
fig3a = dt1[chnm %in% target.chems,]

plot_list = list()

for (i in 1:length(target.chems)){
  
  chem.title = target.chems[i]
  dat_steve <- fig3a[chnm == chem.title,]
  
  y.label = ifelse(i %in% c(1,5), "Fold Change (log-2)", "")
  x.label = ""
  
  setkey(dat_steve, spid)
  dat_steve <- dat_steve[,spid := as.factor(spid)]
  
  #Get hill parameters for each spid
  dat_hill <- data.table(spid = unique(dat_steve[,spid]))
  
  for(chem in dat_hill[,spid]){
    pipefit = tcplFit_Lite(dat_steve[spid == chem,logc], dat_steve[spid == chem,resp], bmad)
    
    aic.vals = c(pipefit$hill_aic, pipefit$gnls_aic, pipefit$cnst_aic)
    
    win.model = ifelse(min(aic.vals) == pipefit$hill_aic, "Hill", 
                       ifelse(min(aic.vals) == pipefit$gnls_aic, "Gain-Loss", "Constant"))
    
    hill_wtp <- ifelse(win.model == "Hill", pipefit$hill_tp, ifelse(win.model =="Gain-Loss", pipefit$gnls_tp, 0))
    hill_wga <- ifelse(win.model == "Hill", pipefit$hill_ga, ifelse(win.model =="Gain-Loss", pipefit$gnls_ga, 0))
    hill_wgw <- ifelse(win.model == "Hill", pipefit$hill_gw, ifelse(win.model =="Gain-Loss", pipefit$gnls_gw, 0))
    
    dat_hill[spid == chem, hill_tp := hill_wtp]
    dat_hill[spid == chem, hill_ga := hill_wga]
    dat_hill[spid == chem, hill_gw := hill_wgw]

  }
  
  # The palette:
  cbbPalette <- c("#0000ff") 
  
  curves_plot <- ggplot(dat_steve, aes(x=logc, y=resp, color = spid))
  
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
    coord_cartesian(ylim = c(0,max(dt1$resp))) +
    scale_colour_manual(values=cbbPalette, guide = "none") +
    ylab(y.label) +
    xlab(x.label) +
    geom_hline(yintercept = 5*bmad1, color = "#ff0000", linetype = 5) +
    ggtitle(chem.title) + theme (axis.title=element_text(size=14,face="bold"),
                                 plot.title=element_text(size=18,face="bold", hjust = 0.5))

    
    plot_list[[i]] = curves_plot

}

Figure_3A = ggarrange(plotlist = plot_list, labels = c(LETTERS[c(1,3,5,7)], LETTERS[c(2,4,6,8)]), ncol = 4, nrow = 2)



# Antagonist Mode
dt2 = dt0[assay == "AR_antagonist",]
dt2[, bval := median(RLU[coli %in% c(12,13,22,23)]), by = apid] #bval using two lowest test concentrations by assay plate
dt2[, pval := median(RLU[rowi == "B" & coli %in% c(9:12)]), by = apid]
dt2[, nval := 100*(RLU - bval)/(pval - bval)]
dt2[, resp := nval]
dt2[, spid := assay]

bmad2 = dt2[coli %in% c(12,13,22,23), mad(nval)]

dt2 = dt2[!chnm %in% c("DMSO", ""),] # remove vehicle and empty wells

# Viability
dt3 = dt0[apid %in% c(2609:2613,2639:2643,2669:2673),] #paired only viability paired with antagonist plates (not agonist)
dt3[, bval := median(RLU[coli %in% c(12,13,22,23)]), by = apid] #bval using two lowest test concentrations by assay plate
dt3[, pval := median(RLU[rowi == "A" & coli %in% c(9:12)], na.rm = T), by = apid]
dt3[, nval := 100*(RLU - bval)/(pval - bval)]
dt3[, resp := nval]
dt3[, spid := assay]

bmad3 = dt3[coli %in% c(12,13,22,23), mad(nval)]

dt3 = dt3[!chnm %in% c("DMSO", ""),] # remove vehicle and empty wells

antagonist.fits = NULL

dt23 = rbind(dt2,dt3)

for (i in 1:length(unique(dt23$chnm))){
  
  chem.title = unique(dt2$chnm)[i]
  sub.dt = dt23[chnm == chem.title,]
  sub.med = sub.dt[, .(max_med = median(resp)), by = .(logc, assay)]
  
  dat_steve <- sub.dt
  setkey(dat_steve, spid)
  dat_steve <- dat_steve[,spid := as.factor(spid)]
  
  #Get hill parameters for each spid
  dat_hill <- data.table(spid = unique(dat_steve[,spid]))
  
  for(chem in dat_hill[,spid]){
    pipefit = tcplFit_Lite(dat_steve[spid == chem,logc], dat_steve[spid == chem,resp], bmad)
    
    aic.vals = c(pipefit$hill_aic, pipefit$gnls_aic, pipefit$cnst_aic)
    
    win.model = ifelse(min(aic.vals) == pipefit$hill_aic, "Hill", 
                       ifelse(min(aic.vals) == pipefit$gnls_aic, "Gain-Loss", "Constant"))
    
    hill_wtp <- ifelse(win.model == "Hill", pipefit$hill_tp, ifelse(win.model =="Gain-Loss", pipefit$gnls_tp, 0))
    hill_wga <- ifelse(win.model == "Hill", pipefit$hill_ga, ifelse(win.model =="Gain-Loss", pipefit$gnls_ga, 0))
    hill_wgw <- ifelse(win.model == "Hill", pipefit$hill_gw, ifelse(win.model =="Gain-Loss", pipefit$gnls_gw, 0))
    
    dat_hill[spid == chem, win_model := win.model,]
    dat_hill[spid == chem, hill_tp := hill_wtp]
    dat_hill[spid == chem, hill_ga := hill_wga]
    dat_hill[spid == chem, hill_gw := hill_wgw]
    dat_hill[spid == chem, max_med := max(sub.med[assay == chem, max_med], na.rm = T)]
  }
  
  # The palette with red:
  cbbPalette <- c("#ff0000", "#666666") 
  
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
    coord_cartesian(ylim = c(max(dt2$resp),0)) +
    scale_colour_manual(values=cbbPalette, name = "Assay", labels = c("AR antagonist", "Cytotoxicity")) +
    scale_shape_manual(values= c(16,1), name = "Assay", labels = c("AR antagonist", "Cytotoxicity")) +
    ylab("% Response") +
    xlab("Concentration (log-10 uM)") +
    geom_hline(yintercept = 3*bmad2, color = "#ff0000", linetype = 5) +
    ggtitle(chem.title) + theme (axis.title=element_text(size=14,face="bold"),
                                 plot.title=element_text(size=18,face="bold", hjust = 0.5)) + 
    scale_y_reverse()

  
  file.name <- paste("chem", i, "png", sep = ".")
  file.path <- "./Figure_3/plots/AR_antagonist"
  ggsave(plot = curves_plot, filename = file.name, 
         path = file.path)
  
  temp = data.table(chnm = chem.title, dat_hill)
  antagonist.fits = rbind(antagonist.fits, temp)
  
}

write.csv(antagonist.fits, "./Figure_3/antagonist_fits.csv",
          row.names = F)

antagonist.hits = antagonist.fits[!chnm %in%  c("BICAL", "DCLN") &  spid == "AR_antagonist" & win_model != "Constant" & max_med > 3*bmad2, chnm]

# AR antagonist vs. Cytotoxicity Selectivity
antagonist.fits[, hitc := ifelse(spid == "AR_antagonist" & max_med > 3*bmad2, 1, ifelse(spid == "Viability" & max_med > 3*bmad3, 1, 0))]
antagonist.fits[, adj_ga := ifelse(spid == "AR_antagonist" & hitc == 1, hill_ga, ifelse(spid == "Viability" & hitc == 1, hill_ga, 2.5))]

delta_ga = antagonist.fits[, .(antag_hitc = hitc[spid == "AR_antagonist"], antag_ga = hill_ga[spid == "AR_antagonist"], delta_ga = adj_ga[spid == "Viability"] - adj_ga[spid == "AR_antagonist"]), by = chnm]

selective.antagonists = delta_ga[chnm != "BICAL" & antag_hitc == 1 & delta_ga > 0.5, chnm]
delta_ga[, select.hitc := ifelse(chnm %in% selective.antagonists, 1, 0)]
delta_ga[, final_ga := ifelse(select.hitc == 1, antag_ga, 2.5)]
delta_ga = delta_ga[!chnm %in% c("BICAL", "DCLN"),]

write.csv(delta_ga, "./Figure_3/delta_ga.csv",
          row.names = F)

antagonist.hitc = data.table(chnm = delta_ga[, chnm], hitc = delta_ga[, antag_hitc], 
                             antagonist_ga = ifelse(delta_ga[, antag_hitc] == 1, delta_ga[, antag_ga], 2.5))

write.csv(antagonist.hitc, "./Figure_3/antagonist.hitc.csv",
          row.names = F)

# Figure 3B
target.chems = c("BICAL", "Flutamide", "Cyproterone acetate", "17-Methyltestosterone", "DCLN", "Spironolactone", 
                 "Hydroxyflutamide", "Simazine")
fig3b = dt23[chnm %in% target.chems,]

plot_list = list()

for (i in 1:length(target.chems)){
  
  chem.title = target.chems[i]
  sub.dt = dt23[chnm == chem.title,]

  y.label = ifelse(i %in% c(1,5), "% Response", "")
  x.label = ifelse(i %in% 5:8, "Concentration (log-10 uM)", "")
  
  dat_steve <- sub.dt
  setkey(dat_steve, spid)
  dat_steve <- dat_steve[,spid := as.factor(spid)]
  
  #Get hill parameters for each spid
  dat_hill <- data.table(spid = unique(dat_steve[,spid]))
  
  for(chem in dat_hill[,spid]){
    pipefit = tcplFit_Lite(dat_steve[spid == chem,logc], dat_steve[spid == chem,resp], bmad)
    
    aic.vals = c(pipefit$hill_aic, pipefit$gnls_aic, pipefit$cnst_aic)
    
    win.model = ifelse(min(aic.vals) == pipefit$hill_aic, "Hill", 
                       ifelse(min(aic.vals) == pipefit$gnls_aic, "Gain-Loss", "Constant"))
    
    hill_wtp <- ifelse(win.model == "Hill", pipefit$hill_tp, ifelse(win.model =="Gain-Loss", pipefit$gnls_tp, 0))
    hill_wga <- ifelse(win.model == "Hill", pipefit$hill_ga, ifelse(win.model =="Gain-Loss", pipefit$gnls_ga, 0))
    hill_wgw <- ifelse(win.model == "Hill", pipefit$hill_gw, ifelse(win.model =="Gain-Loss", pipefit$gnls_gw, 0))
    
    dat_hill[spid == chem, hill_tp := hill_wtp]
    dat_hill[spid == chem, hill_ga := hill_wga]
    dat_hill[spid == chem, hill_gw := hill_wgw]
  }
  
  # The palette with red:
  cbbPalette <- c("#ff0000", "#666666") 
  
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
    coord_cartesian(ylim = c(max(dt2$resp),0)) +
    scale_colour_manual(values=cbbPalette, guide = "none") +
    scale_shape_manual(values= c(16,1), guide = "none") +
    ylab(y.label) +
    xlab(x.label) +
    geom_hline(yintercept = 3*bmad2, color = "#ff0000", linetype = 5) +
    ggtitle(chem.title) + theme (axis.title=element_text(size=14,face="bold"),
                                 plot.title=element_text(size=18,face="bold", hjust = 0.5)) + 
    scale_y_reverse()
  
  plot_list[[i]] = curves_plot
  
}

Figure_3B = ggarrange(plotlist = plot_list, labels = c(LETTERS[c(9,11,13,15)], LETTERS[c(10,12,14,16)]), ncol = 4, nrow = 2)

Figure_3 = ggarrange(plotlist = list(Figure_3A, Figure_3B), ncol = 1, nrow = 2)

ggsave(plot = Figure_3, width = 400, height = 400, device = "png", dpi = 600, units = "mm", filename = "Figure_3.png", 
       path = "./Figure_3/Draft_Figures")


