library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)

### data retrieved for all spids with chmn == 128 reference chnm from InvitroDB on 05/02/2022 
      #for aeids 740, 741 (OT_AR_ARSRC1 8 and 16hr) 
      #and 2386, 2387 (UPITT_HCI_U2OS_AR_TIF2_Nucleoli antagonist and agonist)
      #and 2390, 2391 (UPITT_HCI_U2OS_AR_TIF2_Cytoplasm_Ratio antagonist and agonist)

AR.mc3 = as.data.table(read.csv("./Figure_4/data/AR.mc3.csv", 
                                stringsAsFactors = F, skipNul = T))

AR.mc5 = as.data.table(read.csv("./Figure_4/data/AR.mc5.csv", 
                                stringsAsFactors = F, skipNul = T))

AR2.agonist = as.data.table(read.csv("./Figure_3/agonist.hitc.csv", 
                                     stringsAsFactors = F, skipNul = T))
AR2.agonist = AR2.agonist[!spid %in% c("R1881", "BICAL", "DCLN"),]

AR2.antagonist = as.data.table(read.csv("./Figure_3/delta_ga.csv", 
                                     stringsAsFactors = F, skipNul = T))
AR2.antagonist = AR2.antagonist[!chnm %in% c("R1881", "BICAL", "DCLN"),]

chem.use.class = as.data.table(read.csv("./Figure_4/Ref_chem_use_classes.csv", 
                                        stringsAsFactors = F, skipNul = T))

palette12col = c("#ff0000", "#ff9999", "#0000ff", "#55aaff", "#4c9a2a", "#acdf87", 
                 "#a45ee5", "#bf9b30", "#000000", "#555555", "#999999", "#bbbbbb")

AR2.combined = data.table(chnm = AR2.agonist[,spid], ag.hitc = AR2.agonist[,hitc], ag.modl_ga = AR2.agonist[,agonist_ga])
AR2.combined[, antag.hitc := AR2.antagonist$select.hitc[match(AR2.combined$chnm, AR2.antagonist$chnm)]]
AR2.combined[, antag.modl_ga := AR2.antagonist$final_ga[match(AR2.combined$chnm, AR2.antagonist$chnm)]]
AR2.combined[, hitc := ifelse(ag.hitc == 0 & antag.hitc == 0, 0, 1)]
AR2.combined[, modl_ga := min(ag.modl_ga, antag.modl_ga), by = chnm]
AR2.combined[, use.class := chem.use.class$class[match(AR2.combined$chnm, chem.use.class$chnm)]]

## Odyssey Thera 8-hour

OT8.mc5 = AR.mc5[aeid == 740 & logc_max > 1.477,] #high-dose testing only (21 samples retested at logc_max == 0 generating duplicates stats)
OT8.AR2 = OT8.mc5[, .(OT8.hitc = median(hitc), OT8.ga = median(modl_ga[hitc == 1])), by = chmn]
OT8.AR2$OT8.ga = replace(OT8.AR2$OT8.ga, is.na(OT8.AR2$OT8.ga), 2.5)
OT8.AR2$AR2.agonist.hitc = AR2.agonist$hitc[match(OT8.AR2$chmn, AR2.agonist$spid)]
OT8.AR2$AR2.agonist.ga = AR2.agonist$agonist_ga[match(OT8.AR2$chmn, AR2.agonist$spid)]
OT8.AR2$AR2.antagonist.hitc = AR2.antagonist$select.hitc[match(OT8.AR2$chmn, AR2.antagonist$chnm)]
OT8.AR2$AR2.antagonist.ga = AR2.antagonist$final_ga[match(OT8.AR2$chmn, AR2.antagonist$chnm)]
OT8.AR2[, use.class := chem.use.class$class[match(OT8.AR2$chmn, chem.use.class$chnm)]]
OT8.AR2$use.class = factor(OT8.AR2$use.class, levels = c("drug (androgen)", "drug (anti-androgen)", "drug (estrogen)", "drug (anti-estrogen)", "drug (progestin)", "drug (anti-progestin)", 
                                                         "drug (anti-aromatase)", "drug (cortisol)", "cosmetic",  "industrial", "mycotoxin","pesticide"))

OT8.AR2.ag.plot = ggplot(OT8.AR2, aes(x = AR2.agonist.ga, y = OT8.ga, color = use.class)) + geom_point(size = 3) + theme_bw() +
                    scale_color_manual(values = palette12col, name = "") +
                    coord_cartesian(xlim = c( -3.5,2.5), ylim = c( -3.5,2.5)) +
  scale_x_continuous(breaks = seq(from = -3, to = 3, by =1)) +
  scale_y_continuous(breaks = seq(from = -3, to = 3, by =1)) +
                    xlab("AR2 AC50 (log-10 uM)") + ylab("OT AC50 (log-10 uM)") +
                    geom_abline(slope = 1, intercept = 0, linetype = 5) + 
                    ggtitle("AR2 Agonist vs. OT (8-hour)")

OT8.AR2.antag.plot = ggplot(OT8.AR2, aes(x = AR2.antagonist.ga, y = OT8.ga, color = use.class)) + geom_point(size = 3) + theme_bw() +
                      scale_color_manual(values = palette12col, name = "") +
                      coord_cartesian(xlim = c( -3.5,2.5), ylim = c( -3.5,2.5)) +   
  scale_x_continuous(breaks = seq(from = -3, to = 3, by =1)) +   
  scale_y_continuous(breaks = seq(from = -3, to = 3, by =1)) +
                      xlab("AR2 AC50 (log-10 uM)") + ylab("OT AC50 (log-10 uM)") +
                      geom_abline(slope = 1, intercept = 0, linetype = 5) + 
                      ggtitle("AR2 Antagonist vs. OT (8-hour)")

OT8.combo = OT8.mc5[, .(OT8.hitc = median(hitc), OT8.ga = median(modl_ga[hitc == 1])), by = chmn]
OT8.combo$OT8.ga = replace(OT8.AR2$OT8.ga, is.na(OT8.AR2$OT8.ga), 2.5)
OT8.combo$AR2.hitc = AR2.combined$hitc[match(OT8.AR2$chmn, AR2.combined$chnm)]
OT8.combo$AR2.ga = AR2.combined$modl_ga[match(OT8.AR2$chmn, AR2.combined$chnm)]
OT8.combo$use.class = chem.use.class$class[match(OT8.AR2$chmn, chem.use.class$chnm)]
OT8.combo$use.class = factor(OT8.combo$use.class, levels = c("drug (androgen)", "drug (anti-androgen)", "drug (estrogen)", "drug (anti-estrogen)", "drug (progestin)", "drug (anti-progestin)", 
                                                         "drug (anti-aromatase)", "drug (cortisol)", "cosmetic",  "industrial", "mycotoxin","pesticide"))


write.csv(OT8.combo, "./Figure_4/OT8_AR_combined.csv",
          row.names = F)

OT8.AR2.combo.plot = ggplot(OT8.combo, aes(x = AR2.ga, y = OT8.ga, color = use.class)) + geom_point(size = 3) + theme_bw() +
  scale_color_manual(values = palette12col, name = "") +
  coord_cartesian(xlim = c( -3.5,2.5), ylim = c( -3.5,2.5)) +   
  scale_x_continuous(breaks = seq(from = -3, to = 3, by =1)) +   
  scale_y_continuous(breaks = seq(from = -3, to = 3, by =1)) +
  xlab("AR2 AC50 (log-10 uM)") + ylab("OT AC50 (log-10 uM)") +
  geom_abline(slope = 1, intercept = 0, linetype = 5) + 
  ggtitle("AR2 Combined vs. OT (8-hour)")


ggsave(plot = OT8.AR2.combo.plot, width = 125, height = 125, device = "png", dpi = 600, units = "mm", filename = "OT8.AR.combo.png", 
       path = "./Figure_4/plots")

## Odyssey Thera 16-hour
OT16.mc5 = AR.mc5[aeid == 741 & logc_max > 1.477,] #high-dose testing only (21 samples retested at logc_max == 0 generating duplicates stats)
OT16.AR2 = OT16.mc5[, .(OT16.hitc = median(hitc), OT16.ga = median(modl_ga[hitc == 1])), by = chmn]
OT16.AR2$OT16.ga = replace(OT16.AR2$OT16.ga, is.na(OT16.AR2$OT16.ga), 2.5)
OT16.AR2$AR2.agonist.hitc = AR2.agonist$hitc[match(OT16.AR2$chmn, AR2.agonist$spid)]
OT16.AR2$AR2.agonist.ga = AR2.agonist$agonist_ga[match(OT16.AR2$chmn, AR2.agonist$spid)]
OT16.AR2$AR2.antagonist.hitc = AR2.antagonist$select.hitc[match(OT16.AR2$chmn, AR2.antagonist$chnm)]
OT16.AR2$AR2.antagonist.ga = AR2.antagonist$final_ga[match(OT16.AR2$chmn, AR2.antagonist$chnm)]
OT16.AR2[, use.class := chem.use.class$class[match(OT16.AR2$chmn, chem.use.class$chnm)]]
OT16.AR2$use.class = factor(OT16.AR2$use.class, levels = c("drug (androgen)", "drug (anti-androgen)", "drug (estrogen)", "drug (anti-estrogen)", "drug (progestin)", "drug (anti-progestin)", 
                                                         "drug (anti-aromatase)", "drug (cortisol)", "cosmetic",  "industrial", "mycotoxin","pesticide"))

OT16.AR2.ag.plot = ggplot(OT16.AR2, aes(x = AR2.agonist.ga, y = OT16.ga, color = use.class)) + geom_point(shape = 1, size = 3) + theme_bw() +
  scale_color_manual(values = palette12col, name = "") +
                    coord_cartesian(xlim = c( -3.5,2.5), ylim = c( -3.5,2.5)) +   
  scale_x_continuous(breaks = seq(from = -3, to = 3, by =1)) +   
  scale_y_continuous(breaks = seq(from = -3, to = 3, by =1)) +
                    xlab("AR2 AC50 (log-10 uM)") + ylab("OT AC50 (log-10 uM)") +
                    geom_abline(slope = 1, intercept = 0, linetype = 5) + 
                    ggtitle("AR2 Agonist vs. OT (16-hour)")

OT16.AR2.antag.plot = ggplot(OT16.AR2, aes(x = AR2.antagonist.ga, y = OT16.ga, color = use.class)) + geom_point(shape = 1, size = 3) + theme_bw() +
  scale_color_manual(values = palette12col, name = "") +
                        coord_cartesian(xlim = c( -3.5,2.5), ylim = c( -3.5,2.5)) +   
  scale_x_continuous(breaks = seq(from = -3, to = 3, by =1)) +   
  scale_y_continuous(breaks = seq(from = -3, to = 3, by =1)) +
                        xlab("AR2 AC50 (log-10 uM)") + ylab("OT AC50 (log-10 uM)") +
                        geom_abline(slope = 1, intercept = 0, linetype = 5) + 
                        ggtitle("AR2 Antagonist vs. OT (16-hour)")


OT16.combo = OT16.mc5[, .(OT16.hitc = median(hitc), OT16.ga = median(modl_ga[hitc == 1])), by = chmn]
OT16.combo$OT16.ga = replace(OT16.AR2$OT16.ga, is.na(OT16.AR2$OT16.ga), 2.5)
OT16.combo$AR2.hitc = AR2.combined$hitc[match(OT16.AR2$chmn, AR2.combined$chnm)]
OT16.combo$AR2.ga = AR2.combined$modl_ga[match(OT16.AR2$chmn, AR2.combined$chnm)]
OT16.combo$use.class = chem.use.class$class[match(OT16.AR2$chmn, chem.use.class$chnm)]
OT16.combo$use.class = factor(OT16.combo$use.class, levels = c("drug (androgen)", "drug (anti-androgen)", "drug (estrogen)", "drug (anti-estrogen)", "drug (progestin)", "drug (anti-progestin)", 
                                                         "drug (anti-aromatase)", "drug (cortisol)", "cosmetic",  "industrial", "mycotoxin","pesticide"))


write.csv(OT16.combo, "./Figure_4/OT16_AR_combined.csv",
          row.names = F)

OT16.AR2.combo.plot = ggplot(OT16.combo, aes(x = AR2.ga, y = OT16.ga, color = use.class)) + geom_point(size = 3) + theme_bw() +
  scale_color_manual(values = palette12col, name = "") +
  coord_cartesian(xlim = c( -3.5,2.5), ylim = c( -3.5,2.5)) +   
  scale_x_continuous(breaks = seq(from = -3, to = 3, by =1)) +   
  scale_y_continuous(breaks = seq(from = -3, to = 3, by =1)) +
  xlab("AR2 AC50 (log-10 uM)") + ylab("OT AC50 (log-10 uM)") +
  geom_abline(slope = 1, intercept = 0, linetype = 5) + 
  ggtitle("AR2 Combined vs. OT (16-hour)")

ggsave(plot = OT16.AR2.combo.plot, width = 125, height = 125, device = "png", dpi = 600, units = "mm", filename = "OT16.AR.combo.png", 
       path = "./Figure_4/plots")



## UPitt Nucleoli antagonist
UPitt.Nuc.antag.mc5 = AR.mc5[aeid == 2386 & logc_max > 1.477,] #high-dose testing only (2 samples retested at logc_max == 0 generating duplicates stats)
UPitt.Nuc.antag.summ = UPitt.Nuc.antag.mc5[, .(med.hitc = median(hitc), med.modl_ga = median(modl_ga)), by = chmn]
UPitt.Nuc.antag.summ$med.modl_ga = replace(UPitt.Nuc.antag.summ$med.modl_ga, UPitt.Nuc.antag.summ$med.hitc == 0, 2.5)
UPitt.Nuc.antag.summ$AR2.antagonist.hitc = AR2.antagonist$select.hitc[match(UPitt.Nuc.antag.summ$chmn, AR2.antagonist$chnm)]
UPitt.Nuc.antag.summ$AR2.antagonist.ga = AR2.antagonist$final_ga[match(UPitt.Nuc.antag.summ$chmn, AR2.antagonist$chnm)]
UPitt.Nuc.antag.summ[, use.class := chem.use.class$class[match(UPitt.Nuc.antag.summ$chmn, chem.use.class$chnm)]]
UPitt.Nuc.antag.summ$use.class = factor(UPitt.Nuc.antag.summ$use.class, levels = c("drug (androgen)", "drug (anti-androgen)", "drug (estrogen)", "drug (anti-estrogen)", "drug (progestin)", "drug (anti-progestin)", 
                                                               "drug (anti-aromatase)", "drug (cortisol)", "cosmetic",  "industrial", "mycotoxin","pesticide"))


write.csv(UPitt.Nuc.antag.summ, "./Figure_4/UPitt.Nuc.AR2.antagonist.csv",
          row.names = F)

UPitt.Nuc.antag.plot = ggplot(UPitt.Nuc.antag.summ, aes(x = AR2.antagonist.ga, y = med.modl_ga, color = use.class)) + geom_point(size = 3) + theme_bw() +
  scale_color_manual(values = palette12col, name = "") +
                        coord_cartesian(xlim = c( -3.5,2.5), ylim = c( -3.5,2.5)) +   
  scale_x_continuous(breaks = seq(from = -3, to = 3, by =1)) +   
  scale_y_continuous(breaks = seq(from = -3, to = 3, by =1)) +
                        xlab("AR2 AC50 (log-10 uM)") + ylab("UPitt AC50 (log-10 uM)") +
                        geom_abline(slope = 1, intercept = 0, linetype = 5) + 
                        ggtitle("AR2 Antagonist vs. UPitt Nucleoli Antagonist")

ggsave(plot = UPitt.Nuc.antag.plot, width = 125, height = 125, device = "png", dpi = 600, units = "mm", filename = "UPitt.Nuc.antag.png", 
       path = "./Figure_4/plots")


## UPitt Nucleoli agonist
UPitt.Nuc.ag.mc5 = AR.mc5[aeid == 2387 & logc_max > 1.477,] #high-dose testing only (2 samples retested at logc_max == 0 generating duplicates stats)
UPitt.Nuc.ag.summ = UPitt.Nuc.antag.mc5[, .(med.hitc = median(hitc), med.modl_ga = median(modl_ga)), by = chmn]
UPitt.Nuc.ag.summ$med.modl_ga = replace(UPitt.Nuc.ag.summ$med.modl_ga, UPitt.Nuc.ag.summ$med.hitc == 0, 2.5)
UPitt.Nuc.ag.summ$AR2.agonist.hitc = AR2.agonist$hitc[match(UPitt.Nuc.ag.summ$chmn, AR2.agonist$spid)]
UPitt.Nuc.ag.summ$AR2.agonist.ga = AR2.agonist$agonist_ga[match(UPitt.Nuc.ag.summ$chmn, AR2.agonist$spid)]
UPitt.Nuc.ag.summ[, use.class := chem.use.class$class[match(UPitt.Nuc.ag.summ$chmn, chem.use.class$chnm)]]
UPitt.Nuc.ag.summ$use.class = factor(UPitt.Nuc.ag.summ$use.class, levels = c("drug (androgen)", "drug (anti-androgen)", "drug (estrogen)", "drug (anti-estrogen)", "drug (progestin)", "drug (anti-progestin)", 
                                                               "drug (anti-aromatase)", "drug (cortisol)", "cosmetic",  "industrial", "mycotoxin","pesticide"))

write.csv(UPitt.Nuc.ag.summ, "./Figure_4/UPitt.Nuc.AR2.agonist.csv",
          row.names = F)

UPitt.Nuc.ag.plot = ggplot(UPitt.Nuc.ag.summ, aes(x = AR2.agonist.ga, y = med.modl_ga, color = use.class)) + geom_point(size = 3) + theme_bw() +
  scale_color_manual(values = palette12col, name = "") +
                      coord_cartesian(xlim = c( -3.5,2.5), ylim = c( -3.5,2.5)) +   
  scale_x_continuous(breaks = seq(from = -3, to = 3, by =1)) +   
  scale_y_continuous(breaks = seq(from = -3, to = 3, by =1)) +
                      xlab("AR2 AC50 (log-10 uM)") + ylab("UPitt AC50 (log-10 uM)") +
                      geom_abline(slope = 1, intercept = 0, linetype = 5) + 
                      ggtitle("AR2 Agonist vs. UPitt Nucleoli Agonist")

ggsave(plot = UPitt.Nuc.ag.plot, width = 125, height = 125, device = "png", dpi = 600, units = "mm", filename = "UPitt.Nuc.ag.png", 
       path = "./Figure_4/plots")

# Figure 4
Figure_4 = ggarrange(OT8.AR2.combo.plot + theme(legend.position = "none") + xlab(""), OT16.AR2.combo.plot + xlab("")  + ylab("") + theme(legend.position = "none"), 
                     UPitt.Nuc.ag.plot+ theme(legend.position = "none"), UPitt.Nuc.antag.plot + ylab("") + theme(legend.position = "none"), 
                     labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")

ggsave(plot = Figure_4, width = 300, height = 300, device = "png", dpi = 600, units = "mm", filename = "Figure_4.png", 
       path = "./Figure_4/Draft_Figures")