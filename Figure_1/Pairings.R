library(data.table)
library(ggplot2)
library(grid)

### read in study and filemap ###

# Note that original transfection groups C and D are labeled E anf F (and vice versa) to group homo- and heterodimers together
study = read.csv("./Figure_1/study.csv", stringsAsFactors = F, skipNul = T)

study = as.data.table(study)


filemap = read.csv("./Figure_1/filemap.csv",
                   stringsAsFactors = F, skipNul = T)

filemap = as.data.table(filemap)

### read in data by TR and map to study ###

dt0 = NULL

for (i in 1:nrow(filemap)){
  
  sub.name = paste("./Figure_1/data/",
                   "TRno", filemap[i, TR], ".csv", sep = "")
  
  sub.data = read.csv(sub.name, stringsAsFactors = F, skipNul = T)
  sub.data$Destination.Well = paste(substr(sub.data$Well, 1, 1), as.numeric(substr(sub.data$Well, 2, 3)), sep = "")
  
  temp = study
  temp$Destination.Well = temp$wells
  temp$RLU = sub.data[,3][match(temp$Destination.Well, sub.data$Destination.Well)]
  temp$chnm = ifelse(temp$rowi %in% LETTERS[1:4], "DMSO", filemap[i, Agonist])
  temp$apid = filemap[i, TR]
  temp$agonist = filemap[i, Agonist]
  temp$rep = filemap[i, rep]
  
  dt0 = rbind(dt0, temp)
}

dt0 <- as.data.table(dt0)
setkey(dt0, chnm)
dt0[, bval := median(RLU[chnm == "DMSO"]), by = .(apid, TG)]
dt0[, nval := RLU/bval]
dt1 <- dt0[!is.na(Final.nM),]

# Plot TG group median response by R1881 concentration
summ = dt1[!is.na(Final.nM), .(group.med = median(nval), group.mad = mad(nval)), by = .(TG, Final.nM)]
summ$Final.nM = factor(summ$Final.nM, levels = c("10", "1", "0.1"))

Fig1A.plot = ggplot(summ, aes(x = group.med - 1, y = TG, fill = as.factor(Final.nM))) + 
              geom_errorbar(aes(xmin = group.med-group.mad-1, xmax = group.med+group.mad-1), width = 0.5, position=position_dodge(0.8)) +
              geom_bar(stat = "identity", position=position_dodge(0.8), width = 0.75) + 
              scale_fill_manual(values = c("#00aaff", "#ffb855", "#ff0000"), name = "R1881 (nM)", guide = guide_legend(reverse = TRUE)) + 
              xlab("Median Fold Change Above Vehicle") + 
              ylab("") + 
              theme_bw() + 
              scale_y_discrete(limits=rev) +
              geom_text(aes(x = -0.35, y = 10.5, label = "AR Homodimers"), size = 4, angle = 90, color = "#555555") + 
              geom_segment(aes(x = -0.2, xend = -0.2, y = 8.7, yend = 12.3), size = 0.8, color = "#555555") +
              geom_text(aes(x = 2.85, y = 4.5, label = "AR + SRC1 Heterodimers"), size = 4, angle = 270, color = "#555555") + 
              geom_segment(aes(x = 2.7, xend = 2.7, y = 0.7, yend = 8.3), size = 0.8, color = "#555555") +
              theme(legend.position = c(0.85, 0.3))

ggsave(filename = "Draft_Figure_1A.png", plot = Fig1A.plot, device = "png", 
       path = "./Figure_1/Draft_Figures", 
       width = 150, height = 150, units = "mm", dpi = 600)

# Plot TG group cummulative rZ across R1881 concentrations

bmads = dt0[is.na(Final.nM), mad(nval), by = TG]
summ$bmad = bmads$V1[match(summ$TG, bmads$TG)]
summ[, rZ := 1-(3*(group.mad + bmad)/abs(group.med - 1))]

rZprime = summ[, .(cumm.rZ = sum(rZ)), by = TG]

Fig1B.plot = ggplot(rZprime, aes(x = cumm.rZ, y = TG)) + 
  geom_bar(stat = "identity", fill = "#999999") + 
  coord_cartesian(xlim = c(0,0.15)) +
  xlab("Cummulative rZ Score") + 
  ylab("") + 
  theme_bw() + 
  scale_y_discrete(limits=rev)

ggsave(filename = "Draft_Figure_1B.png", plot = Fig1B.plot, device = "png", 
       path = "./Figure_1/Draft_Figures", 
       width = 150, height = 150, units = "mm", dpi = 600)

grid.newpage()
grid.draw(cbind(ggplotGrob(Fig1A.plot), ggplotGrob(Fig1B.plot), size = "max"))

# manually save
