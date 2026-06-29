# gene reference - select key genes reported in RRMM
setwd("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/")
gene_ref2<- read.delim("./2025_manuscript/nat_gen_script/zold/all_myeloma_and_cancer_driver_genes.txt", sep="\t")
driver_myeloma<- c("TNFRSF17","CD38", "IKZF3",
                   "CUL4A","DDB1","RBX1","PAX5",
                   "IKZF1","CRBN","POU2AF1","TNFRSF13B")
gene_ref<- gene_ref2[gene_ref2$region %in%driver_myeloma, ]
gene_ref$gene_class<-"LOSS"
colnames(gene_ref)<- c("chrom","start","end","region","gene_class")

all_ref<- rbind.data.frame(sv_hotspot2, gene_ref)

write.table(all_ref,"~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/all_myeloma_and_cancer_driver_genes.txt", sep="\t", quote=F)


setwd("~/Library/CloudStorage/Box-Box/RRMMproject/")
all_cnv_37<- read.delim("all_cnv_rrmm_hg37_def.txt")

moffitt_cnv<- read.delim("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2024_final_calls/all_cnv_moffitt_CART_hg37.txt")
moffitt_cnv2<-moffitt_cnv[,c("sample","chr","start","end","total","minor")]
colnames(moffitt_cnv2)<- c("sample","chr","start","end","tot","min")

setwd("~/OneDrive - Memorial Sloan Kettering Cancer Center/File di Ziccheddu, Bachisio - BiSpecific_Calgary_MSKCC/")
setwd("./analysis/liftover/")
list.files()

CNV_calgary<- read.delim("all_cnv_biS_hg37.txt", sep="\t")
CNV_calgary<-CNV_calgary[,-ncol(CNV_calgary)]
CNV_calgary<- CNV_calgary[! CNV_calgary$sample %in% moffitt_cnv2$sample,]
colnames(CNV_calgary)<- colnames(moffitt_cnv2)

### MSK
setwd("~/OneDrive - Memorial Sloan Kettering Cancer Center/MSKCC_project/MSK_TCE/fra_files_isabl_TMP-selected/")
CNV<- read.delim("MSK_TCE_CNV.txt", sep="\t")
CNV$chrom<- gsub("chr","", CNV$chrom)
colnames(CNV)<- colnames(moffitt_cnv2)
CNV<- CNV[! CNV$sample %in% moffitt_cnv2$sample,]


ALL_CNV<- rbind.data.frame(CNV, CNV_calgary, moffitt_cnv2)
ALL_CNV$chr<- gsub("chr","", ALL_CNV$chr)
# ALL_CNV<- ALL_CNV[ALL_CNV$sample %in%clin_all$sample, ]

write.table(ALL_CNV, "~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/ALL_CNV_FINAL.txt",
            sep="\t", quote=F)



########################################################################
##
## upload genomic data
##
#######################################################################

dnds_all2<- readRDS("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2024_final_calls/output_dnds_all.RDS")
# dnds_all_SIG<-dnds_all2$sel_cv
# dnds_all_SIG[dnds_all_SIG$qglobal_cv<0.1,1:10]
dnds_all_ANN<- dnds_all2$annotmuts
all_jco$sample[!all_jco$sample%in% dnds_all_ANN$sampleID]





# setwd("~/OneDrive - Memorial Sloan Kettering Cancer Center/File di Ziccheddu, Bachisio - BiSpecific_Calgary_MSKCC/")
# setwd("./analysis/liftover/")
# ns_calgary<- read.delim("all_snv_nonsyn_biS_hg37.txt", sep="\t")
# colnames(ns_calgary)[1]<-"sample"
# ns_calgary<-ns_calgary[,-c(2:3)]


dnds_all2<- readRDS("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2024_final_calls/output_dnds_all.RDS")
# dnds_all_SIG<-dnds_all2$sel_cv
# dnds_all_SIG[dnds_all_SIG$qglobal_cv<0.1,1:10]
dnds_all_ANN<- dnds_all2$annotmuts
dnds_all_ANN<-dnds_all_ANN[dnds_all_ANN$impact !="Synonymous",]
colnames(dnds_all_ANN)[1]<-"sample"

setwd("~/OneDrive - Memorial Sloan Kettering Cancer Center/MSKCC_project/MSK_TCE/fra_files_isabl_TMP-selected/")
driver_dnds<- read.delim("MSK_TCE_annot_dnds.txt")
colnames(driver_dnds)[1]<-"sample"
driver_dnds<-driver_dnds[driver_dnds$impact !="Synonymous",]
driver_dnds<-driver_dnds[! driver_dnds$sample %in% dnds_all_ANN$sample,]

all_dnds<-rbind.data.frame(driver_dnds, dnds_all_ANN)
clin_all$sample[!clin_all$sample %in%all_dnds$sample]
write.table(all_dnds, "~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/ALL_NS_SNV_FINAL.txt", sep="\t", quote=F)




moffitt_sv<- read.delim("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2024_final_calls/all_sv_moffitt_CART_hg37.txt")
moffitt_sv2<-moffitt_sv[,c(1:8)]
colnames(moffitt_sv2)<-c("sample" , "chrom1" , "pos1"   , "chrom2",  "pos2"   , "SVTYPE" , "strand1", "strand2")

### Calgary
setwd("~/OneDrive - Memorial Sloan Kettering Cancer Center/File di Ziccheddu, Bachisio - BiSpecific_Calgary_MSKCC/")
setwd("./analysis/liftover/")
list.files()
SV_calgary<- read.delim("all_sv_biS_hg37.txt", sep="\t")
SV_calgary2<- as.data.frame(str_split_fixed(SV_calgary$code, "_", 8))
colnames(SV_calgary2)<-c("sample","chrom1","pos1","chrom2","pos2","SVTYPE","strand1","strand2")
SV_calgary22<-SV_calgary2[!SV_calgary2$sample %in% moffitt_sv2$sample,]

setwd("~/OneDrive - Memorial Sloan Kettering Cancer Center/MSKCC_project/MSK_TCE/fra_files_isabl_TMP-selected/")
SV<-  read.delim("MSK_TCE_SV.txt", sep="\t")
SV$sample<- gsub("_merged_svs.pass.annotated.tsv", "", SV$sample)
SV$chrom1<- gsub("chr","", SV$chrom1)
SV$chrom2<- gsub("chr","", SV$chrom2)

ALL_SV_TCE<- rbind.data.frame(SV, SV_calgary22,moffitt_sv2 )
ALL_SV_TCE$chrom1<- gsub("chr","", ALL_SV_TCE$chrom1)
ALL_SV_TCE$chrom2<- gsub("chr","", ALL_SV_TCE$chrom2)
ALL_SV_TCE<-ALL_SV_TCE[,c(1:8)]
ALL_SV_TCE

opposite<- ALL_SV_TCE[ALL_SV_TCE$pos2<ALL_SV_TCE$pos1 & ALL_SV_TCE$chrom1==ALL_SV_TCE$chrom2,]
opposite<- opposite[,c(1,4,5,2,3,6:8)]
colnames(opposite)<- colnames(ALL_SV_TCE)
ALL_SV_TCE_no_opp<- ALL_SV_TCE[!rownames(ALL_SV_TCE) %in% rownames(opposite),]
ALL_SV_TCE2<- rbind.data.frame(ALL_SV_TCE_no_opp, opposite)
ALL_SV_TCE2<-ALL_SV_TCE2[!is.na(ALL_SV_TCE2$SVTYPE),]
ALL_SV_TCE2[ALL_SV_TCE2$SVTYPE=="INS",]
clin_all$sample[!clin_all$sample %in%ALL_SV_TCE$sample]
ALL_SV_TCE2$SVTYPE[ALL_SV_TCE2$SVTYPE =="INS"]<-"DEL"
write.table(ALL_SV_TCE2, "~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/ALL_SV_FINAL.txt", sep="\t", quote=F)



all_tce_CART<- read.delim("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/final_matrix_genomic_drivers.txt")
all_tce_CART$POU2AF1<-0
all_tce_CART$POU2AF1[all_tce_CART$sample %in% c("64451_CG22_2338","S98_EXP1828DNA10_CAR384","SL521899")]<-1 

all_tce_CART$TNFRSF13B<-0
all_tce_CART$TNFRSF13B[all_tce_CART$sample %in% c("IID_H212047_T01_01_WG01","P1814","SL521899")]<-1 
write.table(all_tce_CART, "~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/final_matrix_genomic_drivers.txt", sep="\t", quote=F)
all_tce_CART[all_tce_CART$sample %in% c("64451_CG22_2338","S98_EXP1828DNA10_CAR384","SL521899"),]
  