#code martin
install.packages("BiodiversityR")

install.packages(pkgs=c("BiodiversityR", "vegan", 
                        "Rcmdr", "MASS", "mgcv", 
                        "cluster", "RODBC", "rpart", "effects", "multcomp", 
                        "ellipse", "maptree", "sp", "splancs", "spatial", 
                        "akima", "nnet", "dismo", "raster", "rgdal", 
                        "bootstrap", "PresenceAbsence", 
                        "maxlike", "gbm", "randomForest", "gam", "earth", "mda", 
                        "kernlab", "e1071", "glmnet", "sem", "rgl", "relimp", 
                        "lmtest", "leaps", "Hmisc", "colorspace", "aplpack", 
                        "abind", "XLConnect", "car", "markdown", "knitr",
                        "geosphere", "maptools", "rgeos", "ENMeval", "red"), 
                 dependencies=c("Depends", "Imports"))
setwd("C:/Users/psanchez/OneDrive - Leitat/Desktop/MICROSERVICES/16S")

library("vegan")
library(BiodiversityR)
library(microeco)
library("tibble")
library("dplyr")
library(phyloseq); library(tidyverse)
library("mirlyn")
library("BiocManager")
library(data.table)
library(file2meco)


#Import data at OTU level----
otus <- read.table('OTU_16S.txt', header = T, row.names = 1, sep = '\t')#dataset at OTU level provided by Martin
#otus <- otus[rowSums(otus)>0,]
names(otus) <- gsub('X', '', names(otus))
dim(otus)

#the taxonomy
library(data.table)
#to import a table with different number of columns
tax <- read.table("taxonomy_16S.txt", sep = ',', header=,col.names = paste0("V", seq_len(8)), fill = TRUE)
otus$samples <- row.names(otus)
tax<-as.data.frame(tax)
tax <- tax[match (otus$samples, tax$V1),]
otus$samples <- NULL
row.names(tax)<- tax$V1
row.names(tax) == row.names(otus)
tax$V1<-NULL

names(tax) <- c('Kingdom', 'Phylum', 'Class',  'Order', 'Family','Genus')
tax$Kingdom <- gsub('[0-9]', '', tax$Kingdom)
tax$Kingdom <- gsub('d:', '', tax$Kingdom)
tax$Phylum <- gsub('[0-9]', '', tax$Phylum)
tax$Phylum <- gsub('p:', '', tax$Phylum)
tax$Class <- gsub('[0-9]', '', tax$Class)
tax$Class <- gsub('c:', '', tax$Class)
tax$Order <- gsub('[0-9]', '', tax$Order)
tax$Order <- gsub('o:', '', tax$Order)
tax$Family <- gsub('[0-9]', '', tax$Family)
tax$Family <- gsub('f:', '', tax$Family)
tax$Genus <- gsub('[0-9]', '', tax$Genus)
tax$Genus <- gsub('g:', '', tax$Genus)


#to confirm both datasets have the same observation = ASV
otus.filt <- otus[row.names(otus)%in%row.names(tax),]
taxfilt<-tax[row.names(tax)%in%row.names (otus.filt),]
row.names(taxfilt)==row.names(otus.filt)
#todo ok?? seguimos

#Introducing the LUCAS datasets 
LUCASmtdt <- read.table ('LUCAS-SOIL-2018_mtdt.txt', header = T, sep ='\t')#metadata of all sites from LUCAS
barcodes <- read.table('LUCAS_ID_fastaq_number.txt', header = T, sep='\t')#all barcodes corresponding to LUCAS POINTID
LUCASfull <- read.table('LUCAS-SOIL-2018_old.txt', header = T, sep = '\t')#dataset with physicochemical parameters
#select samples present in both datasets
LUCASmtdt2 <- LUCASmtdt[LUCASmtdt$POINTID%in%barcodes$LUCAS_ID,]
barcodes2 <- barcodes[barcodes$LUCAS_ID%in%LUCASmtdt2$POINTID,]
barcodes2<-barcodes2[match(LUCASmtdt2$POINTID, barcodes2$LUCAS_ID),]

#we loss 40 samples, which barcodes are these?
barcodes3 <- barcodes[!(barcodes$LUCAS_ID%in%LUCASmtdt2$POINTID),]
barcodes3<-barcodes3[match(LUCASmtdt2$POINTID, barcodes3$LUCAS_ID),]
write.table(barcodes3, "barcodes_non-identified.txt",sep = "\t", row.names= FALSE)
LUCASmtdt2<- cbind(LUCASmtdt2, barcodes2)
#select only LC0 related to croplands
mtdtcropland <-LUCASmtdt2[LUCASmtdt2$LC0_Desc == 'Cropland',] 

#adding some more parameters to cropland metadata
LUCASfull <- LUCASfull[LUCASfull$POINTID%in%mtdtcropland$POINTID,]
LUCASfull<-LUCASfull[match(mtdtcropland$POINTID, LUCASfull$POINTID),]
mtdtcropland2 <- cbind(mtdtcropland, LUCASfull$LC1_Desc)
colnames(mtdtcropland2)[14]<- "LC1_Desc"


#select the OTUS that correspond only to cropland sites
otus.filt2 <- otus.filt[names(otus.filt)%in%mtdtcropland$Barcode_ID]
#to extract variable order
paste(unique(names(otus.filt2)), collapse = "','")
target <- c('103','105','106','108','11','110','112','117','12','124','129','130','133','145','148','149','151','155','158','163','164','165','166','167','168','17','171','172','173','176','177','178','18','181','184','188','191','198','201','203','204','205','207','208','21','212','215','216','219','221','224','226','227','228','229','23','232','233','235','238','24','240','244','247','256','259','26','261','262','265','268','278','284','288','290','293','294','296','299','30','300','303','306','310','312','319','320','324','326','327','330','338','339','34','340','341','342','344','345','346','347','348','350','358','361','362','363','366','367','370','378','38','380','382','384','385','391','392','393','395','398','399','402','403','404','406','410','413','418','423','427','428','434','436','440','441','447','455','459','46','460','462','466','468','47','471','472','473','476','478','479','480','482','484','486','487','488','489','49','494','495','5','501','502','504','505','508','510','512','513','514','517','519','521','523','525','526','527','529','53','532','534','535','537','539','541','543','545','547','55','550','551','553','557','558','559','56','562','566','57','577','579','584','585','589','591','592','593','597','599','6','60','603','604','608','61','610','612','615','616','617','618','619','622','626','629','632','638','64','640','641','642','643','644','645','646','65','654','656','66','663','673','675','679','680','682','683','684','689','69','691','693','694','695','697','698','699','7','70','700','704','707','708','710','713','715','719','721','73','730','731','733','735','736','740','741','744','750','752','755','76','762','764','769','77','770','776','780','781','782','786','791','793','800','801','802','805','806','807','809','812','813','816','819','825','830','831','835','837','839','840','843','844','847','848','853','858','86','862','864','874','875','878','881','882','883','884','885','93','95','97','99')
mtdtcropland2<-mtdtcropland2[match(target, mtdtcropland2$Barcode_ID),]
mtdtcropland2$Barcode_ID
row.names(mtdtcropland2) <- mtdtcropland2$POINTID 
names(otus.filt2) == mtdtcropland2$Barcode_ID #TRUE.both dataset has the same order
names(otus.filt2) <- mtdtcropland2$POINTID
#to eliminate OTUS that are not present in cropland sites
otus.filt2 <- otus.filt2[rowSums(otus.filt2)>0,]
dim(otus.filt2)
taxfilt2<-tax[row.names(taxfilt)%in%row.names (otus.filt2),]
row.names(taxfilt2)==row.names(otus.filt2)

#generate a PHYLOSEQ OBJECT: TO SEE THGE RAREFACTUON (myrlin package) ----
tax.m <- as.matrix(taxfilt2); rownames(tax.m) <- row.names(taxfilt2); colnames(tax.m) <- names(taxfilt2)#OJO, row.names de la metadata tienen que ser los mismos que los names de la otu table
physeq <- phyloseq(otu_table(otus.filt2, taxa_are_rows = T), tax_table(tax.m), sample_data(mtdtcropland2))
dim(physeq@tax_table)
rarefaction_df <- phyloseq_to_df(physeq)
#to extract only otutable and tax
rarefaction_df_asv <- get_asv_table(rarefaction_df, OTU="Row.names",Sample="sample", Species = "Genus")
#Creation of rarefaction curve data frame
Rarefy_whole_rep_example <- rarefy_whole_rep(physeq, rep = 1)
library(vegan)
Rarecurve_ex <- mirlyn::rarecurve(Rarefy_whole_rep_example)
pdf("rarecurve.pdf", height = 6, width = 25)
Rarecurve_ex 
dev.off()
#1000 reads seems to be the most optimal number to rarefy.

#calculate abundance / alpha div / beta div normalization at 1.000 (OTUS)----
#Normalization with rarefy_samples
microeco <- phyloseq2meco(physeq)
#microeco <- microtable$new(otu_table = otus.filt2, sample_table = mtdtcropland2, tax_table = taxfilt2)
dim(microeco$sample_table)
microeco$rarefy_samples(sample.size = 50000)#30 out of 332 samples and 584 out of 3997 OTUs were removed 
microeco$tidy_dataset()
microeco$sample_sums()
spdf <- read.table('new_spdf.txt', header = T, row.names = 1, sep = '\t')#dataset at OTU level provided by Martin
microeco$sample_table$ac_zone <- spdf$ac.zone

#To analyse only cereals ----
microeco$sample_sums()
microeco$cal_alphadiv(PD = F)
microeco$alpha_diversity
alphadiv_table <- cbind(microeco$alpha_diversity, microeco$sample_table)
otus_rarefy <- microeco$otu_table
otus_rarefy <- as.data.frame(t(otus_rarefy))
row.names(alphadiv_table) == row.names(otus_rarefy)
mtdtcropland_rarefy <- mtdtcropland2[row.names(mtdtcropland2)%in%alphadiv_table$POINTID,]
mtdtcropland_rarefy$POINTID == alphadiv_table$POINTID 
alphadiv_table$LC1 <- mtdtcropland_rarefy$LC1
alphadiv_table$NUTS_0
alphas__cereals <- subset (alphadiv_table, LC1 == "Common wheat" | LC1 == "Durum wheat" |LC1 == "Barley" |LC1 == "Mix of cereals" |LC1 == "Maize" |LC1 == "Oats" |LC1 == "Rye"| LC1 == "Triticale"  )
alphadiv_table$ac_zone

alphas__cereals %>% 
  ggplot(aes(ac_zone, Shannon, color = ac_zone, fill = ac_zone, group = ac_zone)) +
  geom_point(aes(color = ac_zone, group = ac_zone), position = position_jitterdodge(), size = 4, alpha = 0.2) +
  #geom_violin(width=1) +
  #geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  #facet_wrap(~LC1, scales = "free_x", drop = TRUE)+
  geom_boxplot(aes(fill = ac_zone), position = 'dodge', outlier.colour = NA, alpha = .4) + stat_compare_means(size = 3, label.y = 2)+
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA), 
        panel.border = element_rect(fill=NA), 
        strip.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) + guides(fill = F, shape = F) +
  labs(x = NULL, y = 'Shannon index value', shape = "Month", color = 'Agroclimatic zone (only cereals)')

alphas__cereals$ac_zone <- factor(alphas__cereals$ac_zone, levels = c("BOS", "NEM", "CON", "PAN", "NMA", "SMA", "MED"))

alphas__cereals %>% 
  ggplot(aes(NUTS_0, Observed, color = NUTS_0, fill = NUTS_0, group = NUTS_0)) +
  geom_point(aes(color = NUTS_0, group = NUTS_0), position = position_jitterdodge(), size = 4, alpha = 0.2) +
  #geom_violin(width=1) +
  #geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  facet_wrap(~ac_zone, scales = "free_x", drop = TRUE)+
  geom_boxplot(aes(fill = NUTS_0), position = 'dodge', outlier.colour = NA, alpha = .4) +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA), 
        panel.border = element_rect(fill=NA), 
        strip.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) + guides(fill = F, shape = F) +
  labs(x = NULL, y = 'AD Observed', shape = "Month", color = 'Region (NUTS_0')
dev.off()



microeco__cereals<- clone(microeco)
microeco__cereals$sample_table <- subset (microeco__cereals$sample_table, LC1_Desc == "Common wheat" | LC1_Desc == "Durum wheat" |LC1_Desc == "Barley" |LC1_Desc == "Mix of cereals" |LC1_Desc == "Maize" |LC1_Desc == "Oats" |LC1_Desc == "Rye"| LC1_Desc == "Triticale"  )
microeco__cereals$tidy_dataset()
#most abundant phyla and genera
microeco__cereals$sample_table$full<- interaction(microeco__cereals$sample_table$ac_zone, microeco__cereals$sample_table$NUTS_0)
microeco__cereals$cal_abund()
microeco__cereals$sample_sums()

#TRY TO SHOW RAREFACTION
my_rarefy <- function (x,sa)
  t1 <- mecodev::trans_rarefy$new(microeco__cereals, alphadiv = "Shannon", depth = c(0, 100, 1000, 10000,50000))
t1$plot_rarefy(color = "NUTS_3", show_point = FALSE, add_fitting = FALSE)

#beta diversity(Bray)/Cannonical PCoA ----
otus_rarefy <- microeco__cereals$otu_table
byc.otus <- vegan::vegdist(t(otus_rarefy), method = 'bray', binary = "FALSE")
braymatrix<- as.matrix(byc.otus)
byc.otus.pcoa <- ape::pcoa(byc.otus)
byc.otus.scores <- as.data.frame(byc.otus.pcoa$vectors)
byc.otus.variance <- as.data.frame(byc.otus.pcoa$values$Relative_eig*100)
byc.otus.variance_corr <- as.data.frame(byc.otus.pcoa$values$Rel_corr_eig*100)
byc.otus.variance_cum.corr <- byc.otus.pcoa$values$Cum_corr_eig*100
Axis1<-byc.otus.scores$Axis.1
Axis2<-byc.otus.scores$Axis.2
Axis3<-byc.otus.scores$Axis.3
byc.otus.scores.top3 <- cbind(Axis1, Axis2, Axis3, microeco__cereals$sample_table)


# run diagnostics(cannonical PCoA con BiodiversityR)----
axes <- 50 # number of axis to inspect, increase if you don't see a plateau
nc <- nrow(as.matrix(byc.otus))
success <- data.frame(m=numeric(nc), class.success=numeric(nc))
#doubts: data = design --> metadata,any specific format?
#loop toprocess the first 50 axis of beta diversity matrix
for (i in 1:axes){
  cap <- CAPdiscrim(byc.otus ~ constraining_factor, data = microeco__cereals$sample_table, m = i, add = TRUE)
  success[i,1] <- cap$m
  success[i,2] <- 100/length(cap$group)*length(which(cap$group == cap$CV))
}
# plot success rate for each axis
plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels=success$m, pos=1, cex=0.6)
# take number of axis where the curves starts to plateau, to run the final CAP (let's say 15 for this example)

# run final CAP
CAPdiscrim(OTU.bray ~ constraining_factor,, data=design, m = 15, permutations = 999, add = TRUE)
success <- cbind(data.frame(ISS.cap$group), data.frame(ISS.cap$CV))
colnames(success) <- c("source","classified")
rownames(success) <- rownames(ISS.cap$PCoA)
success <- success[order(success$source),]
