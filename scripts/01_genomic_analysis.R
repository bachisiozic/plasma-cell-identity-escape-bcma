
library(stringr) # done
library(deconstructSigs)
# library(dndscv)
library(mmsig)
library(quantsmooth)
library(dplyr) # done
require(plyr) # done
library(GenomicRanges)
library(IRanges)
library(pheatmap) # done
library(survminer) # done
library(RColorBrewer)  # done
library(survival) # done
library(tidyr) # done
library(data.table) # done
library(ggplot2) # done
library(survminer)  # done
library('GenomicFeatures') # done
library(IRanges) # done
library(maftools) # done
library(ComplexHeatmap) # done
library(colorspace) # done
library(reshape) # done
library(tidyverse)  # done
library(readxl) # done
library(data.table) # done
library(biomaRt) # done
library(forestmodel)
set.seed(1)


####################################################################################
##
## functions
##
####################################################################################
run_km_screen <- function(data,
                          start_col = 17,
                          time_col = "pfs_time",
                          event_col = "pfs_code",
                          min_cases = 3,
                          p_cutoff = 0.05) {
  
  # Make a copy
  df <- data
  
  # Convert all values >1 to 1 in tested columns
  df[, start_col:ncol(df)][df[, start_col:ncol(df)] > 1] <- 1
  
  # Results container
  results <- data.frame(
    p_value = numeric(),
    variable = character(),
    n_positive = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through variables
  for (i in start_col:ncol(df)) {
    
    current_var <- df[, i]
    
    if (length(unique(current_var)) > 1 &&
        sum(current_var == 1, na.rm = TRUE) >= min_cases) {
      
      df$test <- current_var
      
      fit <- surv_fit(
        Surv(df[[time_col]], df[[event_col]]) ~ test,
        data = df
      )
      
      km_p_value <- as.numeric(
        gsub("p = ", "", surv_pvalue(fit)[, 2])
      )
      
      results <- rbind(
        results,
        data.frame(
          p_value = km_p_value,
          variable = colnames(df)[i],
          n_positive = sum(current_var == 1, na.rm = TRUE)
        )
      )
      
      if (km_p_value < p_cutoff) {
        print(
          c(
            p_value = km_p_value,
            variable = colnames(df)[i],
            n_positive = sum(current_var == 1, na.rm = TRUE)
          )
        )
      }
    }
  }
  
  results$p_adj <- p.adjust(results$p_value, method = "BH")
  results <- results[order(results$p_value), ]
  
  return(results)
}


plot_km <- function(data,
                    group_var,
                    time_var = "pfs_time",
                    event_var = "pfs_code",
                    break_time = 180,
                    palette = c("dodgerblue",
                                "firebrick",
                                "forestgreen",
                                "orange",
                                "black")) {
  
  formula <- as.formula(
    paste0("Surv(", time_var, ", ", event_var, ") ~ ", group_var)
  )
  
  fit <- surv_fit(formula, data = data)
  
  ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    palette = palette,
    break.time.by = break_time,
    conf.int = FALSE,
    risk.table = TRUE,
    xlab = "Days",
    ylab = "Proportion"
  )
}

### run fisher's exact test to assess refractoriness using sliding window
run_fisher_screen <- function(jco_clin_genomic,
                              refractory_time,
                              feature_cols = NULL,
                              fish_p_value_cart = NA) {
  
  # define columns safely inside function
  if (is.null(feature_cols)) {
    feature_cols <- 17:ncol(jco_clin_genomic)
  }
  
  ref_jco_int <- list()
  
  for (i in feature_cols) {
    
    if (length(unique(jco_clin_genomic[, i])) > 1 &&
        sum(jco_clin_genomic[, i] == 1, na.rm = TRUE) > 2) {
      
      fish_p_value <- fisher.test(
        table(
          jco_clin_genomic[, i],
          jco_clin_genomic$refractory
        )
      )$p.value
      
      ref_jco_int[[length(ref_jco_int) + 1]] <- data.frame(
        refra_time = refractory_time,
        driver = colnames(jco_clin_genomic)[i],
        p_all = fish_p_value,
        p_CART = fish_p_value_cart
      )
    }
  }
  
  do.call(rbind, ref_jco_int)
}

####################################################################################
##
## upload clinical files
##
####################################################################################

setwd("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/")

# upload gene driver reference list that are recurrently involved by SV (i.e. Rustad et al BCD 2020) 
# or reported as enriched in RRMM or impact clinical outcomes with immunotherapies
# (e.g. Diamond et al. Blood 2025, Maura et al. Nat Cancer 2023)

gene_ref<- read.delim("all_myeloma_and_cancer_driver_genes.txt")


### upload clinical data TCE and cart

clin_all<- read.delim("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/all_tce_cart_clin.txt")
clin_all$note<-""
clin_all$note[clin_all$sample %in%c("P2123","P2289")] <-"post-TCE_ref" # these are the Two patients treated with TCE underwent WGS at progression after one month; they were included in the subsequent clinical analyses, as clonal evolution is generally not observed within such a short timeframe


# figure 1a
clin_all$code<-1
plot_km(  data = clin_all[clin_all$cohort=="CART",],group_var = "code")
# figure 1b
KK<-plot_km(  data = clin_all[clin_all$cohort=="TCE",],group_var = "code")
pdf(paste0("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/new_figure/tce_pfs.pdf"), width = 5, height = 6)
print(KK, newpage = FALSE)
dev.off()


####################################################################################
##
## upload Maura et al. JCO 2024 classification 
##
####################################################################################
#  
all_jco<- read.delim("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/jco_matrix_cart_tce.txt", sep="\t")
head(all_jco)
length(unique(all_jco$sample))
# [1] 99 --> P2336 is the patient that was treated with CART and TCE with a 10-month interval between the two treatments, and the same sample was used for both groups 

# Ccomplex annotation based on  Maura et al. JCO 2024 - https://github.com/UM-Myeloma-Genomics/GCP_MM
all_jco$complex_new<- 1
all_jco$complex_new[all_jco$clusters.new %in% c("CCND1_Simple","HRD_RAS","Simple","HRD_gain")]<-0

###########################################################################################
##
## upload genomic data
##
###########################################################################################

ALL_CNV<- read.delim( "~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/ALL_CNV_FINAL.txt")
all_dnds<- read.delim("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/ALL_NS_SNV_FINAL.txt")
ALL_SV<- read.delim("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/ALL_SV_FINAL.txt")

all_jco$sample[!all_jco$sample%in% all_dnds$sample]
all_jco$sample[!all_jco$sample%in% ALL_SV$sample]
all_jco$sample[!all_jco$sample%in% ALL_CNV$sample]

#######################################################################
##
## JCO driver impact on clinic
##
#######################################################################

## in this analysis each genomic events is counted as one indipendently if it is monoallelic or biallelic
jco_clin_genomic<- merge(clin_all,all_jco, by="sample" )
jco_clin_genomic$Gain_Amp1q[jco_clin_genomic$Gain_Amp1q==2]<-1 
jco_clin_genomic[,16:(ncol(jco_clin_genomic)-1)][jco_clin_genomic[,16:(ncol(jco_clin_genomic)-1)]>1]<-1

jco_clin_sig <- run_km_screen(jco_clin_genomic)
jco_clin_sig_CART <- run_km_screen(jco_clin_genomic[jco_clin_genomic$cohort=="CART",])
jco_clin_sig_tce <- run_km_screen(jco_clin_genomic[jco_clin_genomic$cohort=="TCE",])

# print some results
jco_clin_sig[jco_clin_sig$p_adj<0.12,]
# p_value    variable n_positive      p_adj
# 45 0.0001905497         ATM         10 0.01734002
# 58 0.0022804866        CYLD         40 0.07734618
# 9  0.0025498741  Gain_Amp1q         61 0.07734618
# 28 0.0056271701        TET2         17 0.10266242
# 1  0.0056407926   chr3.gain         25 0.10266242
# 48 0.0089992664       ARID2         14 0.11229838
# 25 0.0091207153       SP140         24 0.11229838
# 90 0.0111064333 complex_new         91 0.11229838

## impact of complex genomic profile on PFS - all TCE but one were complex
plot_km(  data = jco_clin_genomic,group_var = "complex_new")

## figure 2a - impact of complex genomic profile on PFS
plot_km(  data = jco_clin_genomic[jco_clin_genomic$cohort=="CART",],group_var = "complex_new")

###################################################
##
## JCO classification BIALLELIC EVENTS
##
###################################################

# in this analysis you only counts the biallelic or high level CNV
jco_clin_genomic_1_2<- merge(clin_all,all_jco, by="sample" )
# jco_clin_genomic_1_2$Gain_Amp1q[jco_clin_genomic_1_2$Gain_Amp1q==2]<-1
jco_clin_genomic_1_2[,16:(ncol(jco_clin_genomic_1_2)-1)][jco_clin_genomic_1_2[,16:(ncol(jco_clin_genomic_1_2)-1)]==1]<-0
sig_jco_1_2 <- run_km_screen(jco_clin_genomic_1_2)
sig_jco_1_2[sig_jco_1_2$p_adj<0.1 ,]

# p_value    variable n_positive      p_adj
# 12 0.01110643 complex_new         91 0.09029226
# 10 0.01997301     MAF.IGH          7 0.09029226
# 7  0.02562095        TP53         14 0.09029226
# 1  0.03009742  Gain_Amp1q         23 0.09029226

## Figure 5A
plot_km(  data = jco_clin_genomic,group_var = "Gain_Amp1q")

# Figure 
plot_km(  data = jco_clin_genomic,group_var = " MAF.IGH")


###################################################
##
## JCO classification for refractoriness
##
###################################################
jco_clin_genomic<- merge(clin_all,all_jco, by="sample" )
jco_clin_genomic$Gain_Amp1q[jco_clin_genomic$Gain_Amp1q==2]<-1 
jco_clin_genomic[,16:(ncol(jco_clin_genomic)-1)][jco_clin_genomic[,16:(ncol(jco_clin_genomic)-1)]>1]<-1
ref_time<- c(30,60,90, 120,150,180, 270, 360) # adjust definition of refractory based on different time brackets. 

ref_jco=list()
for(j in (1:5)){
ref_jco[[j]] <- run_fisher_screen(
  jco_clin_genomic = jco_clin_genomic,
  refractory_time = ref_time[j]
)
ref_jco2<- do.call(rbind.data.frame, ref_jco)
ref_jco2$fdr<- p.adjust(ref_jco2$p_all, method = "fdr")
ref_jco2_sign<- ref_jco2[ref_jco2$p_all<0.05,]
table(ref_jco2_sign$refra_time, ref_jco2_sign$driver)

#           Amp_8q24.21 ARID2 ATM chr15.gain chr19.gain chr3.gain chr9.gain CNV.Sig complex_new CYLD Del_10p15.3 Del_11q22.1 Gain_Amp1q KRAS
# 180           0     0   1          0          0         1         0       0           1    0           0           1          1    0
# 270           1     1   1          0          0         1         0       1           1    1           0           0          1    1
# 30            0     1   0          0          0         0         0       0           0    1           1           0          1    0
# 360           0     0   1          1          1         1         0       1           1    1           0           0          1    0
# 90            0     0   1          0          0         0         1       0           0    1           0           0          1    0
# 
# MAF.IGH MAX NCOR1 TP53
# 180       1   1     0    0
# 270       1   0     0    1
# 30        0   0     0    0
# 360       0   0     0    0
# 90        1   1     1    0


################################################
##
## biallelic loss of BCMA
##
###############################################

# idenitfy and manually inspect biallelic loss of BCMA

BCMA_ANNOT<-all_dnds[all_dnds$gene =="TNFRSF17",]
BCMA_CNV<- ALL_CNV[ALL_CNV$chr==17 & 
          ALL_CNV$start<12059181 &
          ALL_CNV$end> 12059067 &
          ALL_CNV$min<0.5,]
BCMA_SV<- ALL_SV[ALL_SV$chrom1==17 & 
         ALL_SV$chrom2==17 &
         ALL_SV$pos1<12059181 &
         ALL_SV$pos2> 12059067 ,]

BCMA_ANNOT[BCMA_ANNOT$sample =="55289_CG22_1561",] ## neutral mutation
# sample   chr      pos ref mut     gene strand ref_cod mut_cod ref3_cod mut3_cod aachange
# 1041746 55289_CG22_1561 chr16 11965421   C   T TNFRSF17      1       C       T      TCC      TTC     P33S
# ntchange codonsub   impact             pid
# 1041746     C97T  CCT>TCT Missense ENSP00000053243

# IID_H212040_T01_01_WG01 --> loss BCMA after Bela and before Tec - tec was given without knowing and in fact it was refractory
# IID_H135317_T02_01_WG01 --> loss BCMA after Bela and before Cilta was given without knowing and in fact it was refractory 632

## final annotation for each sample combining SV, CNV, Nonsynonimouse mutations and manual inspection

del_mut_bcma<- c("P2226", "P2798", "P2510" ,"P3298", "P1883" ,"P2085", "P3101")# del mut

mut_only<- c("P1924","P2969") #mut only

non_extra_mut_no_del<-("S105_EXP1828DNA17_CAR_NA")

bi_del<- c("P1759","P2825","P2918","S120_EXP1828DNA32_CAR182", "S143_EXP1828DNA54_CAR511",
           "S157_EXP1828DNA68_CAR182","S163_EXP1828DNA74_CAR182",
           "IID_H212040_T01_01_WG01","IID_H135317_T02_01_WG01",
           "IID_H212061_T01_01_WG01")

all_pts<- c(del_mut_bcma, mut_only, non_extra_mut_no_del, bi_del)
# 
# 
jco_clin_genomic2<-jco_clin_genomic
jco_clin_genomic2$BCMA<- NA
jco_clin_genomic2$BCMA[jco_clin_genomic2$sample %in% bi_del]<- "bi_del"

jco_clin_genomic2$BCMA[jco_clin_genomic2$sample=="P1759"]<-"subclonal_bi_del" ## based on single cell WGS (Lee et al. Nat med 2023)
jco_clin_genomic2$BCMA[jco_clin_genomic2$sample=="P2825"]<-"subclonal_bi_del" ## based on single cell WGS (Lee et al. Nat med 2023)
jco_clin_genomic2$BCMA[paste(jco_clin_genomic2$BCMA , jco_clin_genomic2$cohort) =="NA CART" ]<-"CART WT"
jco_clin_genomic2$BCMA[paste(jco_clin_genomic2$BCMA , jco_clin_genomic2$cohort) =="NA TCE" ]<-"TCE WT"

# table(jco_clin_genomic2$BCMA, jco_clin_genomic2$CYLD)

fit.null3 <- surv_fit(Surv(pfs_time, pfs_code) ~ BCMA, 
                      data = jco_clin_genomic2)

ggsurvplot(fit.null3, data = jco_clin_genomic2, 
                   pval = TRUE,palette=c("dodgerblue",
                                         "firebrick","forestgreen",
                                         "orange","black"),
                   break.time.by = 180,
                   conf.int = FALSE,risk.table = T,
                   xlab = "Days",ylab="Proportion")


# Figure 1F
loss_bcma<-read.delim("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2024_manuscript/2025_biallelic/biallelic_heatmap_2025.txt")
loss_bcma$CART[loss_bcma$CART %in% c( "idecel","carvykti" )]<-"CART"
loss_bcma$pre_bcma[loss_bcma$pre_bcma %in% c("no", "NO PRIOR BCMA")]<-"none"
loss_bcma$pre_bcma[loss_bcma$pre_bcma %in% c( "CAR-T","yes CART","CART"  )]<-"CART"
color_annot_file_first<- list("CART" = c("CART" = "cornflowerblue", "TCE"="orange3"),
                              "pre_bcma" = c("none" =  "grey90",
                                             "CART" = "violet",
                                             "Blenrep"="purple"))

loss_bcma2<-loss_bcma[order(loss_bcma$BCMA),]
rownames(loss_bcma2)<- loss_bcma2$pts
pheatmap(t(loss_bcma2[,c("TRAF3","TRAF2","CYLD")]), show_colnames = T,cluster_cols = F,cluster_rows = F,
         annotation_col = loss_bcma2[,c("CART","pre_bcma")],
         annotation_colors = color_annot_file_first,
         col=c("grey90","coral", "brown4","dodgerblue4","dodgerblue","grey50"))



################################################
##
## nonsynonymous mutations 
##
#################################################

all_dnds_driver<-all_dnds[all_dnds$gene %in% gene_ref$region,]
all_ns_annot_interval_mat<- as.data.frame.matrix(table(all_dnds_driver$sample, all_dnds_driver$gene))
all_ns_annot_interval_mat$sample<- rownames(all_ns_annot_interval_mat)
clin_ns_annot_interval<- join(clin_all,all_ns_annot_interval_mat,  by="sample")
clin_ns_annot_interval[is.na(clin_ns_annot_interval)]<-0
dim(clin_ns_annot_interval)
clin_ns_annot_interval[,16:ncol(clin_ns_annot_interval)][clin_ns_annot_interval[,16:ncol(clin_ns_annot_interval)]>1]<-1
ns_annot_driver <- run_km_screen(clin_ns_annot_interval)
ns_annot_driver_CART <- run_km_screen(clin_ns_annot_interval[clin_ns_annot_interval$cohort=="CART",])
ns_annot_driver_tce <- run_km_screen(clin_ns_annot_interval[clin_ns_annot_interval$cohort=="CART",])

## refractory
ref_ns_focal=list()
for(j in (1:5)){
  ref_ns_focal[[j]] <- run_fisher_screen(
    jco_clin_genomic = clin_ns_annot_interval,
    refractory_time = ref_time[j]
  )
}
ref_ns_focal2<- do.call(rbind.data.frame, ref_ns_focal)
ref_ns_focal2[ref_ns_focal2$p_all<0.05,]

## TP53 is the non synonymous mutation consistently impacting clinical outcomes both PFS and refractory usigng different temporal cut off

################################################
##
## loss of function events
##
#################################################

## Reference
gene_ref$start<- as.numeric(as.character(gene_ref$start))
gene_ref$end<- as.numeric(as.character(gene_ref$end))
gr0 = with(gene_ref, GRanges(chrom, IRanges(start=(start)-10000, end=(end)+10000)))
values(gr0) <- gene_ref[,c("region","gene_class")]

## focal SV
ALL_SV$pos2<- as.numeric(as.character(ALL_SV$pos2))
ALL_SV$pos1<- as.numeric(as.character(ALL_SV$pos1))
SV_focal<- ALL_SV[ALL_SV$chrom1 == ALL_SV$chrom2 & (ALL_SV$pos2-ALL_SV$pos1)<10000000,]
SV_focal$chrom1<- gsub("chr","", SV_focal$chrom1)
SV_focal$chrom2<- gsub("chr","", SV_focal$chrom2)
gr1 = with(SV_focal, GRanges(chrom1, IRanges(start=pos1-10000, end=pos2+10000)))
values(gr1) <- SV_focal[,c("sample","SVTYPE")]
ranges_SV <- merge(as.data.frame(gr0),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
ranges_SV_1 <- ranges_SV[with(ranges_SV, startB <= startA & endB >= endA),]
ranges_SV_2 <- ranges_SV[with(ranges_SV, startB <= startA & endB < endA & startA < endB),]
ranges_SV_3 <- ranges_SV[with(ranges_SV, startB > startA & endB >= endA & endA > startB),]
ranges_SV_4 <- ranges_SV[with(ranges_SV, startB > startA & endB < endA),]
ranges_del_dup<- rbind.data.frame(ranges_SV_1, ranges_SV_2, ranges_SV_3, ranges_SV_4)
ranges_del_sv<-ranges_del_dup[ranges_del_dup$SVTYPE!="DUP",]
ranges_del_sv$code_sv<- paste(ranges_del_sv$sample,ranges_del_sv$seqnames, ranges_del_sv$startA, ranges_del_sv$endA, ranges_del_sv$SVTYPE)
## clin impact
all_sv_focal_interval_mat<- as.data.frame.matrix(table(ranges_del_sv$sample, ranges_del_sv$region))
all_sv_focal_interval_mat$sample<- rownames(all_sv_focal_interval_mat)
clin_sv_focal_interval<- join(clin_all,all_sv_focal_interval_mat,  by="sample")
clin_sv_focal_interval[is.na(clin_sv_focal_interval)]<-0
dim(clin_sv_focal_interval)
clin_sv_focal_interval[,16:ncol(clin_sv_focal_interval)][clin_sv_focal_interval[,16:ncol(clin_sv_focal_interval)]>1]<-1
sv_focal_driver <- run_km_screen(clin_sv_focal_interval)
sv_focal_driver_CART <- run_km_screen(clin_sv_focal_interval[clin_sv_focal_interval$cohort=="CART",])
sv_focal_driver_tce <- run_km_screen(clin_sv_focal_interval[clin_sv_focal_interval$cohort!="CART",])

# "CDKN2C|FAF1" , "GLCCI1" ,  "TNFRSF17" --> signigicant genes

## refractory
ref_sv_focal=list()  
for(j in (1:5)){
ref_sv_focal[[j]] <- run_fisher_screen(
  jco_clin_genomic = clin_sv_focal_interval,
  refractory_time = ref_time[j]
)
}
ref_sv_focal2<- do.call(rbind.data.frame, ref_sv_focal)
ref_sv_focal2[ref_sv_focal2$p_all<0.05,]

# For example this shows that GLCCI1 focal deletions associated with both shorter PFS and refractoriness.

ref_sv_focal2$fdr<- p.adjust(ref_sv_focal2$p_all, method = "fdr")

#deletions using CNV <5 mb - focal deletions and non synonymous mutations are key element to distinguish a large aspecific event (e.g 1q gain) from a specific gene (Maura et al. JCO 2024). 
CNV_any_del<- ALL_CNV[ALL_CNV$min<0.5 & ALL_CNV$tot<5,]
CNV_any_del<- CNV_any_del[! is.na(CNV_any_del$start),]
CNV_any_del<-CNV_any_del[(CNV_any_del$end-CNV_any_del$start)<5000000,]
gr1 = with(CNV_any_del, GRanges(chr, IRanges(start=start-10000, end=end+10000)))
values(gr1) <- CNV_any_del[,c("sample","tot","min")]
ranges_cnv_del <- merge(as.data.frame(gr0),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
ranges_cnv_del_1 <- ranges_cnv_del[with(ranges_cnv_del, startB <= startA & endB >= endA),]
ranges_cnv_del_2 <- ranges_cnv_del[with(ranges_cnv_del, startB <= startA & endB < endA & startA < endB),]
ranges_cnv_del_3 <- ranges_cnv_del[with(ranges_cnv_del, startB > startA & endB >= endA & endA > startB),]
ranges_cnv_del_4 <- ranges_cnv_del[with(ranges_cnv_del, startB > startA & endB < endA),]
ranges_cnv_del<- rbind.data.frame(ranges_cnv_del_1, ranges_cnv_del_2, ranges_cnv_del_3, ranges_cnv_del_4)
ranges_cnv_del$code_cnv<- paste(ranges_cnv_del$sample, ranges_cnv_del$seqnames, ranges_cnv_del$startA, ranges_cnv_del$endA, ranges_cnv_del$tot, ranges_cnv_del$min)
ranges_cnv_del$cnv_type<- "DEL"
## clin impact del  CNV <5 mb
all_cnv_focal_interval_mat<- as.data.frame.matrix(table(ranges_cnv_del$sample, ranges_cnv_del$region))
all_cnv_focal_interval_mat$sample<- rownames(all_cnv_focal_interval_mat)
clin_cnv_focal_interval<- join(clin_all,all_cnv_focal_interval_mat,  by="sample")
clin_cnv_focal_interval[is.na(clin_cnv_focal_interval)]<-0
dim(clin_cnv_focal_interval)
clin_cnv_focal_interval[,16:ncol(clin_cnv_focal_interval)][clin_cnv_focal_interval[,16:ncol(clin_cnv_focal_interval)]>1]<-1
cnv_focal_driver <- run_km_screen(clin_cnv_focal_interval)
cnv_focal_driver_CART <- run_km_screen(clin_cnv_focal_interval[clin_cnv_focal_interval$cohort=="CART",])
cnv_focal_driver_tce <- run_km_screen(clin_cnv_focal_interval[clin_cnv_focal_interval$cohort!="CART",])

# this showed  "POU2AF1","CD38" and CCSER1"        

## refractory
ref_cnv_focal_focal=list()
for(j in (1:5)){
  ref_cnv_focal_focal[[j]] <- run_fisher_screen(
    jco_clin_genomic = clin_cnv_focal_interval,
    refractory_time = ref_time[j]
  )
}
ref_cnv_focal_focal2<- do.call(rbind.data.frame, ref_cnv_focal_focal)
ref_cnv_focal_focal2[ref_cnv_focal_focal2$p_all<0.05,]
         

#deletions using CNV any size
CNV_any_del<- ALL_CNV[ALL_CNV$min<0.5 & ALL_CNV$tot<5,]
CNV_any_del<- CNV_any_del[! is.na(CNV_any_del$start),]
# CNV_any_del<-CNV_any_del[(CNV_any_del$end-CNV_any_del$start)<5000000,]
gr1 = with(CNV_any_del, GRanges(chr, IRanges(start=start-10000, end=end+10000)))
values(gr1) <- CNV_any_del[,c("sample","tot","min")]
ranges_cnv_del <- merge(as.data.frame(gr0),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
ranges_cnv_del_1 <- ranges_cnv_del[with(ranges_cnv_del, startB <= startA & endB >= endA),]
ranges_cnv_del_2 <- ranges_cnv_del[with(ranges_cnv_del, startB <= startA & endB < endA & startA < endB),]
ranges_cnv_del_3 <- ranges_cnv_del[with(ranges_cnv_del, startB > startA & endB >= endA & endA > startB),]
ranges_cnv_del_4 <- ranges_cnv_del[with(ranges_cnv_del, startB > startA & endB < endA),]
ranges_cnv_del<- rbind.data.frame(ranges_cnv_del_1, ranges_cnv_del_2, ranges_cnv_del_3, ranges_cnv_del_4)
ranges_cnv_del$code_cnv<- paste(ranges_cnv_del$sample, ranges_cnv_del$seqnames, ranges_cnv_del$startA, ranges_cnv_del$endA, ranges_cnv_del$tot, ranges_cnv_del$min)
ranges_cnv_del$cnv_type<- "DEL"
## clin impact 
all_cnv_any_interval_mat<- as.data.frame.matrix(table(ranges_cnv_del$sample, ranges_cnv_del$region))
all_cnv_any_interval_mat$sample<- rownames(all_cnv_any_interval_mat)
clin_cnv_any_interval<- join(clin_all,all_cnv_any_interval_mat,  by="sample")
clin_cnv_any_interval[is.na(clin_cnv_any_interval)]<-0
dim(clin_cnv_any_interval)
clin_cnv_any_interval[,16:ncol(clin_cnv_any_interval)][clin_cnv_any_interval[,16:ncol(clin_cnv_any_interval)]>1]<-1
cnv_any_driver <- run_km_screen(clin_cnv_any_interval)
cnv_any_driver_CART <- run_km_screen(clin_cnv_any_interval[clin_cnv_any_interval$cohort=="CART",])
cnv_any_driver_tce <- run_km_screen(clin_cnv_any_interval[clin_cnv_any_interval$cohort!="CART",])
## refractory
ref_cnv_any<-list()
for(j in (1:5)){
  ref_cnv_any[[j]] <- run_fisher_screen(
    jco_clin_genomic = clin_cnv_any_interval,
    refractory_time = ref_time[j]
  )
}
ref_cnv_any2<- do.call(rbind.data.frame, ref_cnv_any)
ref_cnv_any2[ref_cnv_any2$p_all<0.05,]

# SV brk within the gene
ALL_SV$codesv<- paste(ALL_SV$sample, ALL_SV$chrom1, ALL_SV$pos1, ALL_SV$chrom1, ALL_SV$pos2)
SV_1<- ALL_SV[,c("codesv","sample","chrom1","pos1","SVTYPE")]
SV_2<- ALL_SV[,c("codesv","sample","chrom2","pos2","SVTYPE")]
colnames(SV_2)<- colnames(SV_1)
all_ig_conogenic<- rbind.data.frame(SV_1, SV_2)
all_ig_conogenic$pos1<- as.numeric(as.character(all_ig_conogenic$pos1))
gr1 = with(all_ig_conogenic, GRanges(chrom1, IRanges(start=pos1-10000, end=pos1+10000)))
values(gr1) <- all_ig_conogenic[,c("sample","SVTYPE","codesv")]
ranges_sv_brk_gist <- merge(as.data.frame(gr0),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
ranges_sv_brk_gist_1 <- ranges_sv_brk_gist[with(ranges_sv_brk_gist, startB <= startA & endB >= endA),]
ranges_sv_brk_gist_2 <- ranges_sv_brk_gist[with(ranges_sv_brk_gist, startB <= startA & endB < endA & startA < endB),]
ranges_sv_brk_gist_3 <- ranges_sv_brk_gist[with(ranges_sv_brk_gist, startB > startA & endB >= endA & endA > startB),]
ranges_sv_brk_gist_4 <- ranges_sv_brk_gist[with(ranges_sv_brk_gist, startB > startA & endB < endA),]
ranges_sv_brk_gist<- rbind.data.frame(ranges_sv_brk_gist_1, ranges_sv_brk_gist_2, ranges_sv_brk_gist_3, ranges_sv_brk_gist_4)
## clin impact 
all_sv_brk_interval_mat<- as.data.frame.matrix(table(ranges_sv_brk_gist$sample, ranges_sv_brk_gist$region))
all_sv_brk_interval_mat$sample<- rownames(all_sv_brk_interval_mat)
clin_sv_brk_interval<- join(clin_all,all_sv_brk_interval_mat,  by="sample")
clin_sv_brk_interval[is.na(clin_sv_brk_interval)]<-0
dim(clin_sv_brk_interval)
clin_sv_brk_interval[,16:ncol(clin_sv_brk_interval)][clin_sv_brk_interval[,16:ncol(clin_sv_brk_interval)]>1]<-1
sv_brk_driver <- run_km_screen(clin_sv_brk_interval)
sv_brk_driver_CART <- run_km_screen(clin_sv_brk_interval[clin_sv_brk_interval$cohort=="CART",])
sv_brk_driver_tce <- run_km_screen(clin_sv_brk_interval[clin_sv_brk_interval$cohort!="CART",])
## refractory
ref_brk<-list()
for(j in (1:5)){
  ref_brk[[j]] <- run_fisher_screen(
    jco_clin_genomic = clin_sv_brk_interval,
    refractory_time = ref_time[j]
  )
}
ref_brk2<- do.call(rbind.data.frame, ref_brk)
ref_brk2[ref_brk2$p_all<0.05,]



################################################
##
## gain of function events
##
#################################################

## focal SV DUP
ALL_SV$pos2<- as.numeric(as.character(ALL_SV$pos2))
ALL_SV$pos1<- as.numeric(as.character(ALL_SV$pos1))
SV_focal_gain<- ALL_SV[ALL_SV$chrom1 == ALL_SV$chrom2 & (ALL_SV$pos2-ALL_SV$pos1)<10000000,]
SV_focal_gain$chrom1<- gsub("chr","", SV_focal_gain$chrom1)
SV_focal_gain$chrom2<- gsub("chr","", SV_focal_gain$chrom2)
gr1 = with(SV_focal_gain, GRanges(chrom1, IRanges(start=pos1-10000, end=pos2+10000)))
values(gr1) <- SV_focal_gain[,c("sample","SVTYPE")]
ranges_SV <- merge(as.data.frame(gr0),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
ranges_SV_1 <- ranges_SV[with(ranges_SV, startB <= startA & endB >= endA),]
ranges_SV_2 <- ranges_SV[with(ranges_SV, startB <= startA & endB < endA & startA < endB),]
ranges_SV_3 <- ranges_SV[with(ranges_SV, startB > startA & endB >= endA & endA > startB),]
ranges_SV_4 <- ranges_SV[with(ranges_SV, startB > startA & endB < endA),]
ranges_del_dup<- rbind.data.frame(ranges_SV_1, ranges_SV_2, ranges_SV_3, ranges_SV_4)
ranges_del_sv<-ranges_del_dup[ranges_del_dup$SVTYPE=="DUP",]
ranges_del_sv$code_sv<- paste(ranges_del_sv$sample,ranges_del_sv$seqnames, ranges_del_sv$startA, ranges_del_sv$endA, ranges_del_sv$SVTYPE)
## clin impact
all_sv_focal_gain_interval_mat<- as.data.frame.matrix(table(ranges_del_sv$sample, ranges_del_sv$region))
all_sv_focal_gain_interval_mat$sample<- rownames(all_sv_focal_gain_interval_mat)
clin_sv_focal_gain_interval<- join(clin_all,all_sv_focal_gain_interval_mat,  by="sample")
clin_sv_focal_gain_interval[is.na(clin_sv_focal_gain_interval)]<-0
dim(clin_sv_focal_gain_interval)
clin_sv_focal_gain_interval[,16:ncol(clin_sv_focal_gain_interval)][clin_sv_focal_gain_interval[,16:ncol(clin_sv_focal_gain_interval)]>1]<-1
sv_focal_gain_driver <- run_km_screen(clin_sv_focal_gain_interval)
sv_focal_gain_driver_CART <- run_km_screen(clin_sv_focal_gain_interval[clin_sv_focal_gain_interval$cohort=="CART",])
sv_focal_gain_driver_tce <- run_km_screen(clin_sv_focal_gain_interval[clin_sv_focal_gain_interval$cohort!="CART",])

# example of false positive
#
# BTG2 is coming out from this analysis, but only three patients overlapping with other drivers and with 1q gain
# S98_EXP1828DNA10_CAR384 --> chromothripsis on chr 1q
# P2123 focal gain on prior gain
# P2123 focal gain on prior gain 

## refractory
ref_sv_focal_gain<-list()
for(j in (1:5)){
  ref_sv_focal_gain[[j]] <- run_fisher_screen(
    jco_clin_genomic = clin_sv_focal_gain_interval,
    refractory_time = ref_time[j]
  )
}
ref_sv_focal_gain2<- do.call(rbind.data.frame, ref_sv_focal_gain)
ref_sv_focal_gain2[ref_sv_focal_gain2$p_all<0.05,]

#gainetions using CNV <5 mb - focal gain and non synonymous mutations are key element to distinguish a large non-specific event (e.g 1q gain) from a specific gene (Maura et al. JCO 2024). 
CNV_any_gain<- ALL_CNV[ALL_CNV$tot>2.5,]
CNV_any_gain<- CNV_any_gain[! is.na(CNV_any_gain$start),]
CNV_any_gain<-CNV_any_gain[(CNV_any_gain$end-CNV_any_gain$start)<5000000,]
gr1 = with(CNV_any_gain, GRanges(chr, IRanges(start=start-10000, end=end+10000)))
values(gr1) <- CNV_any_gain[,c("sample","tot","min")]
ranges_cnv_gain <- merge(as.data.frame(gr0),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
ranges_cnv_gain_1 <- ranges_cnv_gain[with(ranges_cnv_gain, startB <= startA & endB >= endA),]
ranges_cnv_gain_2 <- ranges_cnv_gain[with(ranges_cnv_gain, startB <= startA & endB < endA & startA < endB),]
ranges_cnv_gain_3 <- ranges_cnv_gain[with(ranges_cnv_gain, startB > startA & endB >= endA & endA > startB),]
ranges_cnv_gain_4 <- ranges_cnv_gain[with(ranges_cnv_gain, startB > startA & endB < endA),]
ranges_cnv_gain<- rbind.data.frame(ranges_cnv_gain_1, ranges_cnv_gain_2, ranges_cnv_gain_3, ranges_cnv_gain_4)
ranges_cnv_gain$code_cnv<- paste(ranges_cnv_gain$sample, ranges_cnv_gain$seqnames, ranges_cnv_gain$startA, ranges_cnv_gain$endA, ranges_cnv_gain$tot, ranges_cnv_gain$min)
ranges_cnv_gain$cnv_type<- "gain"
## clin impact del  CNV <5 mb
all_cnv_focal_interval_mat<- as.data.frame.matrix(table(ranges_cnv_gain$sample, ranges_cnv_gain$region))
all_cnv_focal_interval_mat$sample<- rownames(all_cnv_focal_interval_mat)
clin_cnv_gain<- join(clin_all,all_cnv_focal_interval_mat,  by="sample")
clin_cnv_gain[is.na(clin_cnv_gain)]<-0
dim(clin_cnv_gain)
clin_cnv_gain[,16:ncol(clin_cnv_gain)][clin_cnv_gain[,16:ncol(clin_cnv_gain)]>1]<-1
cnv_focal_driver <- run_km_screen(clin_cnv_gain)
cnv_focal_driver_CART <- run_km_screen(clin_cnv_gain[clin_cnv_gain$cohort=="CART",])
cnv_focal_driver_tce <- run_km_screen(clin_cnv_gain[clin_cnv_gain$cohort!="CART",])

# TNFRSF13B is coming out from this analysis

## refractory
clin_cnv_gain_focal<-list()
for(j in (1:5)){
  ref_cnv_focal_focal[[j]] <- run_fisher_screen(
    jco_clin_genomic = clin_cnv_gain,
    refractory_time = ref_time[j]
  )
}
clin_cnv_gain_focal2<- do.call(rbind.data.frame, clin_cnv_gain_focal)
clin_cnv_gain_focal2[clin_cnv_gain_focal2$p_all<0.05,]


#gainetions using CNV any size
CNV_any_gain<- ALL_CNV[ALL_CNV$tot>2.5,]
CNV_any_gain<- CNV_any_gain[! is.na(CNV_any_gain$start),]
# CNV_any_gain<-CNV_any_gain[(CNV_any_gain$end-CNV_any_gain$start)<5000000,]
gr1 = with(CNV_any_gain, GRanges(chr, IRanges(start=start-10000, end=end+10000)))
values(gr1) <- CNV_any_gain[,c("sample","tot","min")]
ranges_cnv_gain <- merge(as.data.frame(gr0),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
ranges_cnv_gain_1 <- ranges_cnv_gain[with(ranges_cnv_gain, startB <= startA & endB >= endA),]
ranges_cnv_gain_2 <- ranges_cnv_gain[with(ranges_cnv_gain, startB <= startA & endB < endA & startA < endB),]
ranges_cnv_gain_3 <- ranges_cnv_gain[with(ranges_cnv_gain, startB > startA & endB >= endA & endA > startB),]
ranges_cnv_gain_4 <- ranges_cnv_gain[with(ranges_cnv_gain, startB > startA & endB < endA),]
ranges_cnv_gain<- rbind.data.frame(ranges_cnv_gain_1, ranges_cnv_gain_2, ranges_cnv_gain_3, ranges_cnv_gain_4)
ranges_cnv_gain$code_cnv<- paste(ranges_cnv_gain$sample, ranges_cnv_gain$seqnames, ranges_cnv_gain$startA, ranges_cnv_gain$endA, ranges_cnv_gain$tot, ranges_cnv_gain$min)
ranges_cnv_gain$cnv_type<- "gain"
## clin impact 
all_cnv_any_interval_mat<- as.data.frame.matrix(table(ranges_cnv_gain$sample, ranges_cnv_gain$region))
all_cnv_any_interval_mat$sample<- rownames(all_cnv_any_interval_mat)
clin_cnv_any_interval<- join(clin_all,all_cnv_any_interval_mat,  by="sample")
clin_cnv_any_interval[is.na(clin_cnv_any_interval)]<-0
dim(clin_cnv_any_interval)
clin_cnv_any_interval[,16:ncol(clin_cnv_any_interval)][clin_cnv_any_interval[,16:ncol(clin_cnv_any_interval)]>1]<-1
cnv_any_driver <- run_km_screen(clin_cnv_any_interval)
cnv_any_driver_CART <- run_km_screen(clin_cnv_any_interval[clin_cnv_any_interval$cohort=="CART",])
cnv_any_driver_tce <- run_km_screen(clin_cnv_any_interval[clin_cnv_any_interval$cohort!="CART",])
## refractory
ref_cnv_any<-list()
for(j in (1:5)){
  ref_cnv_any[[j]] <- run_fisher_screen(
    jco_clin_genomic = clin_cnv_any_interval,
    refractory_time = ref_time[j]
  )
}
ref_cnv_any2<- do.call(rbind.data.frame, ref_cnv_any)
ref_cnv_any2[ref_cnv_any2$p_all<0.05,]

## this shows motsly 1q gain events

# SV brk_gain within the gene - gain of function
ALL_SV$codesv<- paste(ALL_SV$sample, ALL_SV$chrom1, ALL_SV$pos1, ALL_SV$chrom1, ALL_SV$pos2)
SV_1<- ALL_SV[,c("codesv","sample","chrom1","pos1","SVTYPE")]
SV_2<- ALL_SV[,c("codesv","sample","chrom2","pos2","SVTYPE")]
colnames(SV_2)<- colnames(SV_1)
all_ig_conogenic<- rbind.data.frame(SV_1, SV_2)
all_ig_conogenic$pos1<- as.numeric(as.character(all_ig_conogenic$pos1))
gr1 = with(all_ig_conogenic, GRanges(chrom1, IRanges(start=pos1-500000, end=pos1+500000)))
values(gr1) <- all_ig_conogenic[,c("sample","SVTYPE","codesv")]
ranges_sv_brk_gain_gist <- merge(as.data.frame(gr0),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
ranges_sv_brk_gain_gist_1 <- ranges_sv_brk_gain_gist[with(ranges_sv_brk_gain_gist, startB <= startA & endB >= endA),]
ranges_sv_brk_gain_gist_2 <- ranges_sv_brk_gain_gist[with(ranges_sv_brk_gain_gist, startB <= startA & endB < endA & startA < endB),]
ranges_sv_brk_gain_gist_3 <- ranges_sv_brk_gain_gist[with(ranges_sv_brk_gain_gist, startB > startA & endB >= endA & endA > startB),]
ranges_sv_brk_gain_gist_4 <- ranges_sv_brk_gain_gist[with(ranges_sv_brk_gain_gist, startB > startA & endB < endA),]
ranges_sv_brk_gain_gist<- rbind.data.frame(ranges_sv_brk_gain_gist_1, ranges_sv_brk_gain_gist_2, ranges_sv_brk_gain_gist_3, ranges_sv_brk_gain_gist_4)
## clin impact 
all_sv_brk_gain_interval_mat<- as.data.frame.matrix(table(ranges_sv_brk_gain_gist$sample, ranges_sv_brk_gain_gist$region))
all_sv_brk_gain_interval_mat$sample<- rownames(all_sv_brk_gain_interval_mat)
clin_sv_brk_gain_interval<- join(clin_all,all_sv_brk_gain_interval_mat,  by="sample")
clin_sv_brk_gain_interval[is.na(clin_sv_brk_gain_interval)]<-0
dim(clin_sv_brk_gain_interval)
clin_sv_brk_gain_interval[,16:ncol(clin_sv_brk_gain_interval)][clin_sv_brk_gain_interval[,16:ncol(clin_sv_brk_gain_interval)]>1]<-1
sv_brk_gain_driver <- run_km_screen(clin_sv_brk_gain_interval)
sv_brk_gain_driver_CART <- run_km_screen(clin_sv_brk_gain_interval[clin_sv_brk_gain_interval$cohort=="CART",])
sv_brk_gain_driver_tce <- run_km_screen(clin_sv_brk_gain_interval[clin_sv_brk_gain_interval$cohort!="CART",])
## refractory
ref_brk_gain<-list()
for(j in (1:5)){
  ref_brk_gain[[j]] <- run_fisher_screen(
    jco_clin_genomic = clin_sv_brk_gain_interval,
    refractory_time = ref_time[j]
  )
}
ref_brk_gain2<- do.call(rbind.data.frame, ref_brk_gain)
ref_brk_gain2[ref_brk_gain2$p_all<0.05,]


###########################################################
###
### Bulk RNA seq analysis
###
############################################################

## RNA seq was not available for these WGS, but large RNAseq can still be helpful to 
## investigate certain associations

moffitt_rna<- readxl::read_xlsx("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2024_final_calls/CAR-T_RNAseq_Data.xlsx")
moffitt_rna<- as.data.frame(moffitt_rna)

moffitt_rna_commpass<-moffitt_rna[grep("MMRF", moffitt_rna$Sample_ID),] 
moffitt_rna_moffitt<-moffitt_rna[-grep("MMRF", moffitt_rna$Sample_ID),] 

moffitt_rna_commpass$TNFRSF17 <- scale(moffitt_rna_commpass$`TNFRSF17_Log2_TPM+1`)
moffitt_rna_commpass$XBP1 <- scale(moffitt_rna_commpass$`XBP1_Log2_TPM+1`)
moffitt_rna_commpass$CD38 <- scale(moffitt_rna_commpass$`CD38_Log2_TPM+1`)
moffitt_rna_commpass$TP53 <- scale(moffitt_rna_commpass$`TP53_Log2_TPM+1`)
moffitt_rna_commpass$SP140 <- scale(moffitt_rna_commpass$`SP140_Log2_TPM+1`)

moffitt_rna_moffitt$TNFRSF17 <- scale(moffitt_rna_moffitt$`TNFRSF17_Log2_TPM+1`)
moffitt_rna_moffitt$XBP1 <- scale(moffitt_rna_moffitt$`XBP1_Log2_TPM+1`)
moffitt_rna_moffitt$CD38 <- scale(moffitt_rna_moffitt$`CD38_Log2_TPM+1`)
moffitt_rna_moffitt$TP53 <- scale(moffitt_rna_moffitt$`TP53_Log2_TPM+1`)
moffitt_rna_moffitt$SP140 <- scale(moffitt_rna_moffitt$`SP140_Log2_TPM+1`)
rna_zscore<- rbind.data.frame(moffitt_rna_moffitt, moffitt_rna_commpass)

rna_zscore<- read.delim("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/rna_seq_bcma.txt")
rna_zscore_rr<- rna_zscore[rna_zscore$Disease_State %in% c("ERRMM","LRRMM"),]
rna_zscore_dg<- rna_zscore[!rna_zscore$Disease_State %in% c("ERRMM","LRRMM"),]

## Figure Suppl. 8A
par(mfrow=c(1,1))
summary(lm(rna_zscore_dg$TNFRSF17~rna_zscore_dg$TP53))
plot(rna_zscore_dg$TP53,rna_zscore_dg$TNFRSF17,
     , pch=21, bg="olivedrab2", col="grey30", cex=2,
     ylab=c("TNFRSF17"), xlab=c("TP53"))

## Figure 5C
summary(lm(rna_zscore_rr$TNFRSF17~rna_zscore_rr$TP53))
plot(rna_zscore_rr$TP53,rna_zscore_rr$TNFRSF17,
     , pch=21, bg="goldenrod1", col="grey30", cex=2,
     ylab=c("TNFRSF17"), xlab=c("TP53"))
abline(lm(rna_zscore_rr$TNFRSF17~rna_zscore_rr$TP53))

## Suppl. Figure 7G
summary(lm(rna_zscore_dg$TNFRSF17~rna_zscore_dg$XBP1))
plot(rna_zscore_dg$XBP1,rna_zscore_dg$TNFRSF17,
     , pch=21, bg="olivedrab2", col="grey30", cex=2,
     ylab=c("TNFRSF17"), xlab=c("XBP1"))
abline(lm(rna_zscore_dg$TNFRSF17~rna_zscore_dg$XBP1))

## Suppl. Figure 7F
summary(lm(rna_zscore_rr$TNFRSF17~rna_zscore_rr$XBP1))
plot(rna_zscore_rr$XBP1,rna_zscore_rr$TNFRSF17,
     , pch=21, bg="goldenrod1", col="grey30", cex=2,
     ylab=c("TNFRSF17"), xlab=c("XBP1"))

## Suppl. Figure 7D
summary(lm(rna_zscore_dg$TNFRSF17~rna_zscore_dg$CD38))
plot(rna_zscore_dg$CD38,rna_zscore_dg$TNFRSF17,
     , pch=21, bg="olivedrab2", col="grey30", cex=2,
     ylab=c("TNFRSF17"), xlab=c("TP53"))
abline(lm(rna_zscore_dg$TNFRSF17~rna_zscore_dg$CD38))

text("R2=0.2463", x=1, y=2, adj=0)
text("p<0.000001", x=1, y=1.5, adj=0)
dev.off()

## Figure 4D
summary(lm(rna_zscore_rr$TNFRSF17~rna_zscore_rr$CD38))
plot(rna_zscore_rr$CD38,rna_zscore_rr$TNFRSF17,
     , pch=21, bg="goldenrod1", col="grey30", cex=2,
     ylab=c("TNFRSF17"), xlab=c("CD38"))
abline(lm(rna_zscore_rr$TNFRSF17~rna_zscore_rr$CD38))

# Example of how RNAseq can help to guide selection of distinct genomic drivers
#
# SP140 is a gene whose role in myeloma is poorly characterized. A mixture of large and focal deletions seems to associate with poor outcomes, but 
# non synonymous mutations are neutral.For this type of genes whose role in anti-BCMA therapy is well known, we check the expression. 

summary(lm(rna_zscore_rr$SP140~rna_zscore_rr$TNFRSF17))
plot(rna_zscore_rr$SP140,rna_zscore_rr$TNFRSF17,
     , pch=21, bg="goldenrod1", col="grey30", cex=2,
     ylab=c("TNFRSF17"), xlab=c("CD38"))

# SP140-BCMA Adjusted R-squared:  0.03497 

#############################################################
##
## example of how to check and manually curate driver genes
##
##############################################################
# 
# we check the size of the deletion (focal vs large) and the position of the SV breakpoints
# proportion of focal events, events with biallelic loss and/or nonsynonymous mutations are key to define a specific gene as driver of resistance 

# POU2AF1 is particularly important to check as its deletions have never been reported in MM.
all_cnv_any_interval_mat[all_cnv_any_interval_mat$POU2AF1>0,]$sample
all_sv_focal_interval_mat[all_sv_focal_interval_mat$POU2AF1>0,]$sample
all_cnv_focal_interval_mat[all_cnv_focal_interval_mat$POU2AF1>0,]$sample

# "64451_CG22_2338"     --> focal gain without SV supproting it     
# "S98_EXP1828DNA10_CAR384" --> deletions with brk on the gene
# "SL521899" --> deletion after gain  

# "FT-SA243202"     --> whole arm event
# "IID_H212040_T01_01_WG01"    --> whole arm event with FAT3 mutation (biallelic)
# "IID_H212047_T01_01_WG01"  --> whole arm
# "IID_H212048_T01_01_WG01"  --> whole arm with ATM biallelic mut
# "P2337" --> whole arm event with FAT3 mutation (biallelic)
# "P3234"        --> gain of function event on POU2AF1 on whole arm CN-LOH
# "P3277"     --> gain of function event on POU2AF1 on whole arm CN-LOH
# "S120_EXP1828DNA32_CAR182"  whole arm and possible gain of function

#############################################################
##
## final matrix
##
##############################################################
col_sel_features<- c("CDKN2C","CCSER1","GLCCI1", "Gain_Amp1q",   "RPL5" ,   "ATM"  , "TP53","MAF",
                     "XBP1" ,    "CD38"    , "IKZF3" ,  "TNFRSF17","POU2AF1", "TNFRSF13B")

all_tce_CART<- read.delim("~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/final_matrix_genomic_drivers2.txt")
# write.table(all_tce_CART, "~/OneDrive - Memorial Sloan Kettering Cancer Center/CART_Moffitt/2025_manuscript/nat_gen_script/final_matrix_genomic_drivers2.txt",
#             sep="\t", quote=F)

## define complex genomics
all_tce_CART$complex_genomoics<- 0
all_tce_CART$complex_genomoics[rowSums(all_tce_CART[,c("CDKN2C","CCSER1","GLCCI1",
                                                       "Gain_Amp1q",   "RPL5" ,
                                                       "ATM"  , "MAF","TP53")])>0]<- "complex"
## define loss of plasma cell genes
all_tce_CART$plasmacell<- 0
all_tce_CART$plasmacell[rowSums(all_tce_CART[,c( "XBP1" ,"CD38", "IKZF3" ,
                                                 "TNFRSF17","TNFRSF13B",
                                                 "POU2AF1")])>0]<- "plasmacell"
# figure 4A
color_annot_final<- list("refractory" = c("0" = "cornflowerblue", "1"="firebrick2"),
                              "pfs_code" = c("0" =  "lightskyblue1",
                                             "1" ="mediumpurple4"),
                              "os_code" = c( "0" ="darkolivegreen1", "1" = "forestgreen"),
                              "emd" = c("No"="goldenrod","Yes"="darkolivegreen4"),
                              "pre_bcma" = c("No"="cornflowerblue", "Yes"="darkblue"),
                              "cohort"=c("CART" = "grey15", "TCE"="grey90"))


pheatmap(t(all_tce_CART[,col_sel_features]), show_colnames = T,cluster_cols = T,cluster_rows = F,
         annotation_col = all_tce_CART[,c("emd","pfs_code","os_code",
                                          "refractory" ,"cohort")],
         annotation_colors = color_annot_final,
         col=c("grey90", "red", "brown4","dodgerblue"))



# figure 4G
fit.null3 <- surv_fit(Surv(pfs_time, pfs_code) ~ plasmacell ,
                      data =all_tce_CART)

fig<-ggsurvplot(fit.null3, data =all_tce_CART,
                pval = TRUE,risk.table = TRUE, 
                palette=c("dodgerblue",
                          "firebrick"),
                break.time.by = 180)

print(fig)


# figure 5D
fit.null3 <- surv_fit(Surv(pfs_time, pfs_code) ~ complex_genomoics ,
                      data =all_tce_CART)

fig<-ggsurvplot(fit.null3, data =all_tce_CART,
                pval = TRUE,risk.table = TRUE, 
                palette=c("dodgerblue",
                          "firebrick"),
                break.time.by = 180)
print(fig)

all_tce_CART$final_code<-1
all_tce_CART$final_code[paste(all_tce_CART$complex_genomoics, all_tce_CART$plasmacell)=="complex plasmacell"]<-0
cox_model<- coxph(Surv(pfs_time, pfs_code) ~ complex_genomoics+ plasmacell + emd + age + pre_bcma + num_lines + cohort, 
                  data = all_tce_CART)
# Figure 5E
summary(cox_model)
forest_model(cox_model)


## refractory
all_tce_CART_ref<- all_tce_CART
all_tce_CART_ref$plasmacell[all_tce_CART_ref$plasmacell!=0]<-1
all_tce_CART_ref$complex_genomoics[all_tce_CART_ref$complex_genomoics!=0]<-1
final_referactory<-list()
for(j in (1:length(ref_time))){
  final_referactory[[j]] <- run_fisher_screen(
    jco_clin_genomic = all_tce_CART_ref,
    feature_cols = 14:38,
    refractory_time = ref_time[j]
  )
}
final_referactory2<- do.call(rbind.data.frame, final_referactory)
final_referactory2<- final_referactory2[final_referactory2$driver %in% c(col_sel_features,
                                                                         "complex_genomoics",
                                                                         "plasmacell"),]
# Suppl. Table 
final_referactory2$driver<- factor(final_referactory2$driver, levels = unique(final_referactory2$driver))
refr_sig_final<- final_referactory2[final_referactory2$p_all<0.05,]
genomic_events<- c(col_sel_features,"complex_genomoics",  "plasmacell")
genomic_events[! genomic_events%in% refr_sig_final$driver]
genomic_events[genomic_events%in% refr_sig_final$driver]


# sensitivity statement excluding the two TCE patients whose WGS was obtained only at early progression.
all_tce_CART_only_pre<- all_tce_CART[!all_tce_CART$sample %in% c("P2123","P2289"),]
sensitivity_pre <- run_km_screen(all_tce_CART_only_pre, start_col = 14)
sensitivity_pre[sensitivity_pre$p_value>0.05,]

# only RPL5 is >0.05 after removing the two patients whose samples were collected at progression after rapid progression from the treatment start (i.e. primary refractory).
#
# 18 0.13047795     RPL5         32 0.15420121


####################################################################
###
### plot minially deleted regions - for example CD38
###
#####################################################################

ALL_CNV_cd38<- ALL_CNV[ALL_CNV$sample %in% all_tce_CART$sample[all_tce_CART$CD38>0],] ## select deleted
ALL_CNV_cd382<-ALL_CNV_cd38[ALL_CNV_cd38$sample %in% all_tce_CART$sample,]# pre treatment sample only
ALL_CNV_cd382<-ALL_CNV_cd382[ALL_CNV_cd382$chr==4 & ALL_CNV_cd382$min==0 &
                               ALL_CNV_cd382$start<gene_ref$end[gene_ref$region=="CD38"] &
                               ALL_CNV_cd382$end>gene_ref$start[gene_ref$region=="CD38"],]


sv_cd38<- ALL_SV[ALL_SV$chrom1==4 & ALL_SV$chrom2==4 & ALL_SV$pos1<gene_ref$end[gene_ref$region=="CD38"]& 
                 ALL_SV$pos2>gene_ref$start[gene_ref$region=="CD38"],]
sv_cd38<-sv_cd38[sv_cd38$sample %in% all_tce_CART$sample,]
sv_cd38<-sv_cd38[sv_cd38$chrom1==4 & sv_cd38$chrom2==4 &
                   sv_cd38$pos1<gene_ref$end[gene_ref$region=="CD38"] &
                   sv_cd38$pos2>gene_ref$start[gene_ref$region=="CD38"],]

### two patients have focal deletions not capture by CNV callers. 
sv_cd38[! sv_cd38$sample %in%ALL_CNV_cd382$sample, ]

ALL_CNV_cd382$color<- "coral"
ALL_CNV_cd382$color[ALL_CNV_cd382$tot<0.6]<- "black"
ALL_CNV_cd382$order<-0
ALL_CNV_cd382$order[ALL_CNV_cd382$sample %in% all_tce_CART$sample[all_tce_CART$cohort=="CART"]]<-1



## add focal SV coordinates
ALL_CNV_cd38_cnv_sv<-rbind.data.frame(ALL_CNV_cd382, 
                                        c("S186_CAR_344_CD138_pos_CAR344",4,15741285,15804966, 1, 0,"coral", 1 ),
                                        c("IID_H201460_T03_01_WG01",4,12393088,15311603, 1, 0,"coral", 0))
ALL_CNV_cd38_cnv_sv$start<- as.numeric(as.character(ALL_CNV_cd38_cnv_sv$start))
ALL_CNV_cd38_cnv_sv$end<- as.numeric(as.character(ALL_CNV_cd38_cnv_sv$end))
ALL_CNV_cd38_cnv_sv<- ALL_CNV_cd38_cnv_sv[order(ALL_CNV_cd38_cnv_sv$order, (ALL_CNV_cd38_cnv_sv$end-ALL_CNV_cd38_cnv_sv$start)),]
length(unique(ALL_CNV_cd38_cnv_sv$sample))
length(unique(ALL_CNV_cd38$sample))

### Figure 4B
plot(NULL, 
     xlim = c(1, max(ALL_CNV_cd38_cnv_sv$end)),  # Set x-axis limits (adjust as needed)
     ylim = c(1, nrow(ALL_CNV_cd38_cnv_sv)),  # Set y-axis limits (adjust as needed)
     xlab = "",       # Remove x-axis label
     ylab = "",       # Remove y-axis label
     xaxt = "n",      # Remove x-axis ticks and numbers
     yaxt = "n",      # Remove y-axis ticks and numbers
     bty = "n",       # Remove the box around the plot
     main = "",      # Remove the title
     axes = FALSE)   # Remove both axes entirely

for(i in (1:nrow(ALL_CNV_cd38_cnv_sv))){
  segments(x0 =ALL_CNV_cd38_cnv_sv$start[i], x1 = ALL_CNV_cd38_cnv_sv$end[i], 
           y0=i, y1 = i, col=ALL_CNV_cd38_cnv_sv$color[i], lwd = 5)
}


par(new=TRUE, xpd=T)
paintCytobands(gsub("chr","",4), pos = c(0, -0.4), units = "bases", width = 0.3, cex.leg = 0.1,
               bands = "major", legend=F, length.out =max(ALL_CNV_cd38_cnv_sv$end))

