library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)

assay.vol.ul = 45

### read in studies and filemap ###
ag.study = read.csv("./Figure_6_and_7/agonist_study.csv",
                 stringsAsFactors = F, skipNul = T)

ag.study = as.data.table(ag.study)

antag.study = read.csv("./Figure_6_and_7/antagonist_study.csv",
                    stringsAsFactors = F, skipNul = T)

antag.study = as.data.table(antag.study)


filemap = read.csv("./Figure_6_and_7/filemap.csv",
                   stringsAsFactors = F, skipNul = T)

filemap = as.data.table(filemap)


### read in data by TR and map to study ###

dt0 = NULL

for (i in 1:nrow(filemap)){
  
  sub.name = paste("./Figure_6_and_7/data/",
                   "TRno", filemap[i, TR], ".csv", sep = "")
  
  sub.data = read.csv(sub.name, stringsAsFactors = F, skipNul = T)
  sub.data$Destination.Well = paste(sub.data$Well.Row, sub.data$Well.Col, sep = "")
  
  temp.study = if(filemap[i,mode] == "agonist"){ag.study}else{antag.study}
  temp = temp.study[cpid == filemap[i, cpid],]
  temp$rowi = substr(temp$Destination.Well, 1, 1)
  temp$coli = substr(temp$Destination.Well, 2, nchar(temp$Destination.Well))
  temp$RLU = sub.data[,4][match(temp$Destination.Well, sub.data$Destination.Well)]
  temp$apid = filemap[i, TR]
  temp$assay = filemap[i, assay]
  temp$mode = filemap[i, mode]
  temp$biogroup = filemap[i, biogroup]
  temp$replicate = filemap[i, replicate]
  
  dt0 = rbind(dt0, temp)
}

dt0[, logc := log10(Transfer.Volume*stock.mM/(assay.vol.ul+Transfer.Volume/1000))]
setkey(dt0, chnm)

dt0$biogroup = factor(dt0$biogroup, levels = c("CYP1A2", "CYP2A6",  "CYP2B6",  "CYP2C8",  "CYP2C9",  "CYP2C19",
                                               "CYP2D6",  "CYP2E1",  "CYP2J2",  "CYP3A4", "Bgal", "No_RNA"))

dt0 = dt0[biogroup != "CYP2E1-WT",] #only using original CYP2E1 mRNA sequence data
dt0[, Destination.Plate.Barcode := paste(biogroup, cpid, sep = "_")]


#remove affected wells from Echo exceptions reports
exceptions = read.csv("./Figure_6_and_7/exceptions.csv",
                       stringsAsFactors = F, skipNul = T)
exceptions = as.data.table(exceptions)
exceptions = exceptions[!grep(x = Destination.Plate.Barcode, pattern = "CYP2E1.wt"),] #only using original CYP2E1 mRNA sequence data

exceptions.summ = exceptions[, .(Source.Well = unique(Source.Well)), by = .(Source.Well, Destination.Plate.Barcode, replicate, mode)]

delete.rows = NULL
for (n in 1:nrow(exceptions.summ)){
  
  temp.row = exceptions.summ[n,]
  temp.rows = which(dt0$Source.Plate.Barcode == "DPID_128" & dt0$Source.Well == temp.row[,Source.Well] & 
                      dt0$Destination.Plate.Barcode == temp.row[,Destination.Plate.Barcode] & dt0$replicate == temp.row[,replicate] &
                      dt0$mode == temp.row[,mode])
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
dt1 = dt0[mode == "agonist" & assay == "AR",]
dt1[, bval := median(RLU[coli %in% c(12,13,22,23)]), by = apid] #bval using two lowest test concentrations by assay plate
dt1[, pval := median(RLU[rowi == "B" & coli %in% c(9:12)]), by = apid]
dt1[, nval := log(RLU/bval, 2)]
dt1[, resp := nval]
dt1[, spid := biogroup]

bmad1 = dt1[coli %in% c(12,13,22,23), mad(nval)]

dt1 = dt1[!chnm %in% c("DMSO", ""),] # remove vehicle and empty wells

agonist.fits = NULL

plot_list = list()

for (i in 1:length(unique(dt1$chnm))){
  
  for(j in 3:12){
  
  chem.title = paste(unique(dt1$chnm)[i], unique(dt1$spid)[j], sep = ": ")
  sub.dt = dt1[chnm == unique(dt1$chnm)[i] & spid %in%  unique(dt1$spid)[c(1,2,j)],]
  sub.med = sub.dt[, .(max_med = median(resp)), by = .(logc, spid)]
  
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
    
    dat_hill[spid == chem, chnm := unique(dt1$chnm)[i]]
    dat_hill[spid == chem, win_model := win.model]
    dat_hill[spid == chem, hill_tp := hill_wtp]
    dat_hill[spid == chem, hill_ga := hill_wga]
    dat_hill[spid == chem, hill_gw := hill_wgw]
    dat_hill[spid == chem, max_med := max(sub.med[spid == chem, max_med], na.rm = T)]
  }
  
  # The palette with black:
  cbbPalette <- c("#0000ff", "#000000", "#999999") 
  
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
    scale_colour_manual(values=cbbPalette, name = "") +
    ylab("Fold Change (log-2)") +
    xlab("Concentration (log-10 uM)") +
    geom_hline(yintercept = 5*bmad1, color = "#ff0000", linetype = 5) +
    ggtitle(chem.title) + theme (axis.title=element_text(size=8,face="bold"),
                                 plot.title=element_text(size=12,face="bold", hjust = 0.5))

  plot_list[[(i-1)*10+(j-2)]] = curves_plot
  

    
  file.name <- paste("chem", i, unique(dt1$spid)[j], "png", sep = ".")
  file.path <- "./Figure_6_and_7/plots/agonist"
  ggsave(plot = curves_plot, filename = file.name, 
         path = file.path)
  
  agonist.fits = rbind(agonist.fits, dat_hill)
  }
  
}

ggsave(
  filename = "agonist_cyp_plots.pdf", 
  path = "./Figure_6_and_7/plots",
  plot = marrangeGrob(plot_list, nrow=5, ncol=2), 
  width = 9, height = 9
)

#remove repeats for No_RNA and Bgal and add hitc
agonist.fits = agonist.fits[, .(win_model = unique(win_model), hill_tp = median(hill_tp), hill_ga = median(hill_ga), 
                                hill_gw = median(hill_gw), max_med = median(max_med)), by = .(chnm, spid)]
agonist.fits[, hitc := ifelse(max_med > 5*bmad1, 1, 0)]

write.csv(agonist.fits, "./Figure_6_and_7/agonist_fits.csv",
          row.names = F)



# Figure 6 (Agonist with CYP Metabolism)
target.chems = c("17alpha-Estradiol", "Mestranol", "Flutamide", "Prochloraz", "Danazol")
target.biogroups = c("CYP1A2", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2J2", "Bgal", "No_RNA")
target.biogroups = factor(target.biogroups, levels = c("CYP1A2", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2J2", "Bgal", "No_RNA"))
fig6 = dt1[chnm %in% target.chems & biogroup %in% target.biogroups,]

plot_list = list()

for (i in 1:length(target.chems)){
  for (j in 1:(length(target.biogroups)-2)){
  
    sub.dt = fig6[chnm == target.chems[i] & biogroup %in% target.biogroups[c(j,6,7)],]
    biotitle = ifelse(i == 1, as.character(target.biogroups[j]), "")
    y.label = ifelse(j == 1, as.character(target.chems[i]), "")
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
    
    # The palette with black:
    cbbPalette <- c("#0000ff", "#000000", "#999999") 
    
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
      coord_cartesian(ylim = c(0,max(dt1$resp)), xlim = c(min(fig6$logc), max(fig6$logc))) +
      scale_colour_manual(values=cbbPalette, labels = c("CYP", "Bgal", "No RNA"), name = "") +
      scale_x_continuous(breaks = c(-3,-2,-1,0,1,2)) +
      ylab(y.label) +
      xlab("") +
      geom_hline(yintercept = 5*bmad1, color = "#ff0000", linetype = 5) +
      ggtitle(biotitle) + theme (axis.title=element_text(size=18,face="bold"),
                                   plot.title=element_text(size=18,face="bold", hjust = 0.5)) +
      theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
            legend.text=element_text(size=20), legend.key.size = unit(3,"line")) +
      guides(color = guide_legend(override.aes = list(size = 4)))
  
  plot_list[[5*(i-1)+j]] = curves_plot
  
  } 
}

Figure_6 = ggarrange(plotlist = plot_list, labels = LETTERS[1:25], ncol = 5, nrow = 5, common.legend = TRUE, legend = "bottom")

ggsave(plot = Figure_6, width = 500, height = 500, device = "png", dpi = 600, units = "mm", filename = "Figure_6.png", 
       path = "./Figure_6_and_7/Draft_Figures")




### Antagonist Mode ###
dt2 = dt0[mode == "antagonist" & assay == "AR",]
dt2[, bval := median(RLU[coli %in% c(12,13,22,23)]), by = apid] #bval using two lowest test concentrations by assay plate
dt2[, pval := median(RLU[rowi == "B" & coli %in% c(9:12)]), by = apid]
dt2[, nval := 100*(RLU - bval)/(pval - bval)]
dt2[, resp := nval]
dt2[, spid := biogroup]

bmad2 = dt2[coli %in% c(12,13,22,23), mad(nval)]

dt2 = dt2[!chnm %in% c("DMSO", ""),] # remove vehicle and empty wells


antagonist.fits = NULL

plot_list = list()

for (i in 1:length(unique(dt2$chnm))){
  for(j in 3:12){
    
    chem.title = paste(unique(dt2$chnm)[i], unique(dt2$biogroup)[j], sep = ": ")
    sub.dt = dt2[chnm == unique(dt2$chnm)[i] & biogroup %in%  unique(dt2$biogroup)[c(1,2,j)],]
    sub.med = sub.dt[, .(max_med = median(resp)), by = .(logc, spid)]
    
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
    
    dat_hill[spid == chem, chnm := unique(dt2$chnm)[i]]
    dat_hill[spid == chem, win_model := win.model,]
    dat_hill[spid == chem, hill_tp := hill_wtp]
    dat_hill[spid == chem, hill_ga := hill_wga]
    dat_hill[spid == chem, hill_gw := hill_wgw]
    dat_hill[spid == chem, max_med := max(sub.med[spid == chem, max_med], na.rm = T)]
  }
  
  # The palette with black:
  cbbPalette <- c("#ff0000", "#000000", "#999999") 
  
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
    scale_colour_manual(values=cbbPalette, name = "") +
    scale_shape_manual(values=c(16,16,1), name = "") +
    ylab("% Inhibition") +
    xlab("Concentration (log-10 uM)") +
    geom_hline(yintercept = 3*bmad2, color = "#ff0000", linetype = 5) +
    ggtitle(chem.title) + theme (axis.title=element_text(size=8,face="bold"),
                                 plot.title=element_text(size=12,face="bold", hjust = 0.5)) + 
    scale_y_reverse() 
  
    plot_list[[(i-1)*10+(j-2)]] = curves_plot
  
  file.name <- paste("chem", i, unique(dt2$biogroup)[j], "png", sep = ".")
  file.path <- "./Figure_6_and_7/plots/antagonist"
  ggsave(plot = curves_plot, filename = file.name, 
         path = file.path)
  
  antagonist.fits = rbind(antagonist.fits, dat_hill)
  }
}

ggsave(
  filename = "antgonist_cyp_plots.pdf", 
  path = "./Figure_6_and_7/plots",
  plot = marrangeGrob(plot_list, nrow=5, ncol=2), 
  width = 9, height = 9
)


#remove repeats for No_RNA and Bgal and add hitc
antagonist.fits = antagonist.fits[, .(win_model = unique(win_model), hill_tp = median(hill_tp), hill_ga = median(hill_ga), 
                                hill_gw = median(hill_gw), max_med = median(max_med)), by = .(chnm, spid)]
antagonist.fits[, hitc := ifelse(max_med > 3*bmad2, 1, 0)]

write.csv(antagonist.fits, "./Figure_6_and_7/antagonist_fits.csv",
          row.names = F)


# Figure 7 (Antagonist with CYP Metabolism)
target.chems = c("Flutamide", "Hydroxyprogesterone caproate", "Testosterone propionate", "17alpha-Ethinylestradiol", "Equilin")
target.biogroups = c("CYP1A2", "CYP2C8", "CYP2C9", "CYP2C19", "CYP3A4", "Bgal", "No_RNA")
target.biogroups = factor(target.biogroups, levels = c("CYP1A2", "CYP2C8", "CYP2C9", "CYP2C19", "CYP3A4", "Bgal", "No_RNA"))
fig7 = dt2[chnm %in% target.chems & biogroup %in% target.biogroups,]

plot_list = list()

for (i in 1:length(target.chems)){
  for (j in 1:(length(target.biogroups)-2)){
    
    sub.dt = fig7[chnm == target.chems[i] & biogroup %in% target.biogroups[c(j,6,7)],]
    biotitle = ifelse(i == 1, as.character(target.biogroups[j]), "")
    y.label = ifelse(j == 1, as.character(target.chems[i]), "")
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
    
    # The palette with black:
    cbbPalette <- c("#ff0000", "#000000", "#999999") 
    
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
      coord_cartesian(ylim = c(max(dt2$resp),0), xlim = c(min(fig7$logc), max(fig7$logc))) +
      scale_colour_manual(values=cbbPalette, labels = c("CYP", "Bgal", "No RNA"), name = "") +
      scale_shape_manual(values=c(16,16,1), labels = c("CYP", "Bgal", "No RNA"), name = "") +
      scale_x_continuous(breaks = c(-3,-2,-1,0,1,2)) +
      ylab(y.label) +
      xlab("") +
      geom_hline(yintercept = 3*bmad2, color = "#ff0000", linetype = 5) +
      ggtitle(biotitle) + theme (axis.title=element_text(size=16,face="bold"),
                                 plot.title=element_text(size=18,face="bold", hjust = 0.5)) +
      scale_y_reverse(breaks = c(0,20,40,60,80,100)) +
      theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
            legend.text=element_text(size=20), legend.key.size = unit(3,"line")) +
      guides(color = guide_legend(override.aes = list(size = 4)))
    
    plot_list[[5*(i-1)+j]] = curves_plot
    
  } 
}

Figure_7 = ggarrange(plotlist = plot_list, labels = LETTERS[1:25], ncol = 5, nrow = 5, common.legend = TRUE, legend = "bottom")

ggsave(plot = Figure_7, width = 500, height = 500, device = "png", dpi = 600, units = "mm", filename = "Figure_7.png", 
       path = "./Figure_6_and_7/Draft_Figures")

