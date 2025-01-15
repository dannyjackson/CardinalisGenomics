# C1_genelists.sh
# analyzing gene lists

cd /xdisk/mcnew/dannyjackson/cardinals/relevantgenelists

cp /xdisk/mcnew/dannyjackson/cardinals/fst/relevantgenenames_noca_windowed_fst.txt .
cp /xdisk/mcnew/dannyjackson/cardinals/fst/relevantgenenames_pyrr_windowed_fst.txt .
cp /xdisk/mcnew/dannyjackson/cardinals/raisd/noca/relevantgenenames_noca_urban_raisd.txt .
cp /xdisk/mcnew/dannyjackson/cardinals/raisd/noca/relevantgenenames_noca_rural_raisd.txt .
cp /xdisk/mcnew/dannyjackson/cardinals/raisd/pyrr/relevantgenenames_pyrr_urban_raisd.txt .
cp /xdisk/mcnew/dannyjackson/cardinals/raisd/pyrr/relevantgenenames_pyrr_rural_raisd.txt .

cp /xdisk/mcnew/dannyjackson/cardinals/tajima/nocaurban.taj.relevantgenenames_50kb.txt .
cp /xdisk/mcnew/dannyjackson/cardinals/tajima/nocarural.taj.relevantgenenames_50kb.txt .
cp /xdisk/mcnew/dannyjackson/cardinals/tajima/pyrrurban.taj.relevantgenenames_50kb.txt .
cp /xdisk/mcnew/dannyjackson/cardinals/tajima/pyrrrural.taj.relevantgenenames_50kb.txt .





# fst
import pandas as pd
import numpy as np

pyrr = pd.read_csv('relevantgenenames_pyrr_windowed_fst.txt')
noca = pd.read_csv('relevantgenenames_noca_windowed_fst.txt')
np.intersect1d(pyrr, noca)
# array(['CARCAR_R14856', 'Dgkq', 'Nudt12', 'Pip5k1c'], dtype=object)

pyrr_df = np.array(pyrr)
noca_df = np.array(noca)

dif1 = np.setdiff1d(pyrr_df, noca_df)
dif2 = np.setdiff1d(noca_df, pyrr_df)

print(list(dif1))
print(list(dif2))
# dif 1
# Aimp1,Alcam,Anxa7_0,Anxa7_1,Avd_1,Avd_3,Baz1b,CARCAR_R00463,CARCAR_R00752,CARCAR_R04365,CARCAR_R04519,CARCAR_R07473,CARCAR_R12379,CARCAR_R12385,CARCAR_R15577,CARCAR_R15737,Capn5_0,Cck,Ccl4_2,Ccl5_2,Cdh22,Cdkal1,Cfl2,Chst11_1,Cmbl_0,Col20a1,Creb3l3_1,Ctnna3_0,Dap,Dazl,Dmbx1,Dnai1,Dtnbp1_0,Eefsec,Eif3f,Exoc2,Fgf18,Fkbp6,Fzd9,Gimd1,Gpatch2l,Gpr83_2,Ift122,Ift27,Jarid2,Kif24,Klhl17,March6,Mre11,Mtmr4,Ndrg1,Otud7a,Panx1,Pctp,Pla2g4b,Pla2g4e_1,Prvm,Rftn1,Rho,Ruvbl1,Scfd1,Sdk1,Sec61g,Six2,Slc12a2,Snx6,Syt1_1,Tln1,Tpk1,Trak1_0,Trim39,Ubap1,Wisp1,Zfhx4
# dif 2
# Anks1b,Arl8a,CARCAR_R04343,CARCAR_R08758,CARCAR_R10590,CARCAR_R14533,CARCAR_R14987,Cadn,Cdc42,Cdh18,Cdon,Cl004,Cldn16_1,Cldn19_1,Cplx3,Crim1,Ctdp1,Cxcr5,Dchs2,Ddx6,Dpagt1,Dpp10,Dyrk4,Fam219b,Foxp2,Gli3,Gucy1a1,H2afx_0,Hmbs,Hoxa2,Hoxa3,Hoxa5,Hoxa6,Il1r2,Immp2l_0,Lrp8,Magohb,Megf11,Mpi,Mtco2_4,Myo10_1,Myo10_2,Ncald,Nck2,Nrx3a,P3h1,Ptpn7,Rad51ap1,Scamp2,Slc6a5,Tab2,Tmcc3,Tnik_2,Ulk3,Vps11


# local analysis
install.packages('ggvenn')
# library(venn)
# venn(5, ilab=TRUE, zcolor = "style")

# fst
library(ggvenn)

df_pyrr<-read.csv("relevantgenenames_pyrr_windowed_fst.txt", header=FALSE)
df_noca<-read.csv("relevantgenenames_noca_windowed_fst.txt", header=FALSE)
colnames(df_pyrr) <- "Gene"
colnames(df_noca) <- "Gene"

x <- list(
pyrr = df_pyrr$Gene,
noca = df_noca$Gene
)


pdf(file = "species_fst.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("pyrr", "noca"),
    fill_color = c("#4EAFAF", "#FF817E"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()


both_df = intersect(df_pyrr$Gene, df_noca$Gene)


write.csv(both_df, file="intersection_fst.csv", row.names=FALSE)

sed -i 's/\"//g' intersection_fst.csv
tail -n +2 intersection_fst.csv > intersection_fst.csv.tmp
mv intersection_fst.csv.tmp intersection_fst.csv

# raisd

library(ggvenn)

df_pyrr_urban<-read.csv("relevantgenenames_pyrr_urban_raisd.txt", header=FALSE)
df_pyrr_rural<-read.csv("relevantgenenames_pyrr_rural_raisd.txt", header=FALSE)
df_noca_urban<-read.csv("relevantgenenames_noca_urban_raisd.txt", header=FALSE)
df_noca_rural<-read.csv("relevantgenenames_noca_rural_raisd.txt", header=FALSE)

colnames(df_pyrr_urban) <- "Gene"
colnames(df_pyrr_rural) <- "Gene"
colnames(df_noca_urban) <- "Gene"
colnames(df_noca_rural) <- "Gene"

x <- list(
pyrr_urban = df_pyrr_urban$Gene,
pyrr_rural = df_pyrr_rural$Gene,
noca_urban = df_noca_urban$Gene,
noca_rural = df_noca_rural$Gene
)


pdf(file = "species_raisd.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("pyrr_urban", "pyrr_rural", "noca_rural", "noca_urban"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()


both_df = intersect(df_pyrr$Gene, df_noca$Gene)


write.csv(both_df, file="intersection_fst.csv", row.names=FALSE)

# fst and raisd
library(ggvenn)

df_pyrr_urban<-read.csv("relevantgenenames_pyrr_urban_raisd.txt", header=FALSE)
df_pyrr_rural<-read.csv("relevantgenenames_pyrr_rural_raisd.txt", header=FALSE)
df_noca_urban<-read.csv("relevantgenenames_noca_urban_raisd.txt", header=FALSE)
df_noca_rural<-read.csv("relevantgenenames_noca_rural_raisd.txt", header=FALSE)

df_pyrr<-read.csv("relevantgenenames_pyrr_windowed_fst.txt", header=FALSE)
df_noca<-read.csv("relevantgenenames_noca_windowed_fst.txt", header=FALSE)
colnames(df_pyrr) <- "Gene"
colnames(df_noca) <- "Gene"

colnames(df_pyrr_urban) <- "Gene"
colnames(df_pyrr_rural) <- "Gene"
colnames(df_noca_urban) <- "Gene"
colnames(df_noca_rural) <- "Gene"



x <- list(
pyrr_urban = df_pyrr_urban$Gene,
noca_urban = df_noca_urban$Gene,
pyrr = df_pyrr$Gene,
noca = df_noca$Gene
)


pdf(file = "species_raisd_fst.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("noca", "noca_urban", "pyrr_urban", "pyrr"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()


intersect(df_noca$Gene, df_pyrr_urban$Gene)
# "Hoxa2" "Hoxa3" "Hoxa5" "Hoxa6"
intersect(df_noca$Gene, df_noca_urban$Gene)
# "Hoxa2" "Hoxa3" "Tab2" 
intersect(df_noca$Gene, df_pyrr$Gene)
# "CARCAR_R14856" "Dgkq"          "Nudt12"        "Pip5k1c"
intersect(df_noca_urban$Gene, df_pyrr_urban$Gene)

# unique to urban (noca)
raisd_noca_urban_unique <- c("Adgrg4","Ahcyb","Aifm3","Ap5s1","Apmap","Arhgef6","Arid1b","Arx_1","Asip","Atp9b","Barx1b","Bhlhe23","CARCAR_R00615","CARCAR_R04350","CARCAR_R05153","CARCAR_R06321","CARCAR_R08420","CARCAR_R09689","CARCAR_R14549","CARCAR_R15815","Cd247","Cd2ap","Cd40lg","Cdk17","Chid1","Cntnap5","Col17a1","Csgalnact1","Cst7","Ctnnd2_1","Dicer1","Dlg2","Ehd4","Fcgbp_0","Fcgbp_1","Gabra5","Gtf2f2","Ha1f_2","Htr3a_0","Ikzf1","Itgb3","Klhl5","Mboat1","Mcc","Mid1","Mrps28","Nkx28","Nol4_0","Olig3","Parp12_1","Pim1_0","Pole","Pou2f1","Ptprr","Rbmx","Rprml","Sh2d4a","Slc10a2","Slc25a30","Slc30a9","Slc7a11_0","Slk","Snx10","Sult1b1_0","Sult1d1","Syndig1","Sytl5","Tab2","Tbx20","Thsd7a","Tmem33","Tpd52","Tpt1","Ubl3","Ugt2c1","Usp28","Vgll3","Wdr19","Wnk2","Zbtb16a","Zc3h12d","Zscan20_1","Zscan32_1")


# unique to urban (pyrr)
raisd_pyrr_urban_unique<-c("A2ml1_2","Acr_7","Arhgap27_1","CARCAR_R04345","CARCAR_R04347","CARCAR_R06009","CARCAR_R06010","CARCAR_R06691","CARCAR_R15075","Ccnd2","Chst1","Cnksr2_1","Col6a3","Cst7","Ctnnbl1","Derl2_0","Dhx33","Drd1","Fgf23","Fgf6","Fuk","Gabpa","Ha1f_2","Hoxa1","Hoxa2","Hoxa3","Hoxa5","Hoxa6","Hoxa7","Itpripl1_0","Itpripl1_1","Kcnc1","Kcnk16_0","Mia3_0","Mia3_1","Mis12","Nr0b1","Pcdh10_1","Pdpr","Pds5a","Phf20","Pim1_1","Plekhm1","Prl_0","Rfx4","Rimklb","Rusc2","Satb2","Slc12a7","Slc30a10","Slc8a3","Smoc1","St3gal2","Tmem263_1","Trmt61b","Ube2k","Uvrag","Vstm2l")

intersect(raisd_noca_urban_unique, raisd_pyrr_urban_unique)
# "Cst7"(innate immune)  "Ha1f_2"




# tajima

library(ggvenn)

df_pyrr_urban<-read.csv("pyrrurban.taj.relevantgenenames_50kb.txt", header=FALSE)
df_pyrr_rural<-read.csv("pyrrrural.taj.relevantgenenames_50kb.txt", header=FALSE)
df_noca_urban<-read.csv("nocaurban.taj.relevantgenenames_50kb.txt", header=FALSE)
df_noca_rural<-read.csv("nocarural.taj.relevantgenenames_50kb.txt", header=FALSE)

colnames(df_pyrr_urban) <- "Gene"
colnames(df_pyrr_rural) <- "Gene"
colnames(df_noca_urban) <- "Gene"
colnames(df_noca_rural) <- "Gene"

x <- list(
pyrr_urban = df_pyrr_urban$Gene,
pyrr_rural = df_pyrr_rural$Gene,
noca_urban = df_noca_urban$Gene,
noca_rural = df_noca_rural$Gene
)


pdf(file = "species_tajima.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("pyrr_urban", "pyrr_rural", "noca_rural", "noca_urban"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()


both_df = intersect(df_pyrr_urban$Gene, df_noca_urban$Gene)


write.csv(both_df, file="intersection_fst.csv", row.names=FALSE)































# local analyses
# make plots

# BiocManager::install("org.Gg.eg.db")
library(org.Gg.eg.db)
library(limma)
library(rrvgo)

# noca_fst_GOterms_pantherslim_FDR
go_list <- c("GO:0009952","GO:0003002","GO:0007389")

simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/noca_fst_GOterms_pantherslim_FDR.scatter.png", width = 7, height =7, units = "in", res = 300)
scatterPlot(simMatrix, reducedTerms)
dev.off()

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/noca_fst_GOterms_pantherslim_FDR.tree.png", width = 7, height =7, units = "in", res = 300)

treemapPlot(reducedTerms)
dev.off()

# noca_fst_GOterms_pantherslim_nocorrection
go_list <- c("GO:0033962","GO:0061709","GO:0036297","GO:0006739","GO:0017148","GO:0009225","GO:0006783","GO:0016237","GO:0007224","GO:0034727","GO:0042168","GO:0089718","GO:0009952","GO:0006778","GO:0046148","GO:0043009","GO:0140029","GO:0048568","GO:0003002","GO:0033013","GO:0042440","GO:0007157","GO:0007032","GO:0046928","GO:1901379","GO:0007389","GO:0019722","GO:0010506","GO:0042594","GO:0000422","GO:0009991","GO:0048598","GO:0031667","GO:0001501","GO:0007033","GO:0000045","GO:0006417","GO:0034248","GO:1903008","GO:0000725","GO:0000724","GO:0016236","GO:0098742","GO:0018105","GO:0006914","GO:0061919","GO:0018209","GO:0009790","GO:0006310","GO:0006887","GO:0007417","GO:0051248","GO:0022411","GO:0006302","GO:0009887","GO:0031329","GO:0032940","GO:0010608","GO:0046903","GO:0140694","GO:0098609","GO:0007005","GO:0141124","GO:0140352","GO:0009605","GO:0070925","GO:0007155","GO:0007275","GO:0048731","GO:0032501","GO:0022607","GO:0048856","GO:0006357","GO:0032502","GO:0031323","GO:0019222","GO:0051171","GO:0006355","GO:2001141","GO:0080090","GO:0010468","GO:0010556","GO:0031326","GO:0009889","GO:0060255","GO:0050896","GO:0050794","GO:0050789","GO:0065007","GO:0008150","UNCLASSIFIED","GO:0050794","GO:0050789","GO:0065007","GO:0008150","UNCLASSIFIED")

simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/noca_fst_GOterms_pantherslim_nocorrection.scatter.png", width = 7, height =7, units = "in", res = 300)
scatterPlot(simMatrix, reducedTerms)
dev.off()

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/noca_fst_GOterms_pantherslim_nocorrection.tree.png", width = 7, height =7, units = "in", res = 300)

treemapPlot(reducedTerms)
dev.off()


# pyrr_fst_GOterms_pantherslim_FDR


simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/pyrr_fst_GOterms_pantherslim_FDR.scatter.png", width = 7, height =7, units = "in", res = 300)
scatterPlot(simMatrix, reducedTerms)
dev.off()

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/pyrr_fst_GOterms_pantherslim_FDR.tree.png", width = 7, height =7, units = "in", res = 300)

treemapPlot(reducedTerms)
dev.off()


# pyrr_fst_GOterms_pantherslim_nocorrection
go_list <- c("GO:0042138","GO:0007586","GO:0007095","GO:0031204","GO:0071108","GO:0036158","GO:0044818","GO:0044773","GO:0006739","GO:0030042","GO:0006620","GO:0010972","GO:0044774","GO:0006303","GO:0035721","GO:0070536","GO:1902750","GO:0010507","GO:0010389","GO:0000077","GO:0042770","GO:0071482","GO:0034198","GO:0043162","GO:0007416","GO:0007019","GO:0000723","GO:0032200","GO:1902749","GO:0031570","GO:0055078","GO:0070286","GO:0035567","GO:0051261","GO:0007156","GO:0007093","GO:1901991","GO:0031330","GO:0061982","GO:0010506","GO:0009895","GO:0045930","GO:0000075","GO:0009416","GO:0051606","GO:1901988","GO:0098742","GO:0034329","GO:0009314","GO:0010948","GO:0032984","GO:0006417","GO:0034248","GO:1903046","GO:1901990","GO:0060070","GO:0016055","GO:0198738","GO:0000725","GO:0000724","GO:0045786","GO:0046434","GO:0051321","GO:0006338","GO:0006325","GO:0007155","GO:0050808","GO:0031329","GO:0098609","GO:0071824","GO:0010628","GO:1901987","GO:0034330","GO:0007346","GO:0006310","GO:0006650","GO:0009628","GO:0009894","GO:0046488","GO:0022411","GO:0006302","GO:0046486","GO:0043933","GO:0006644","GO:0009605","GO:0006974","GO:0033554","GO:0019637","GO:0022607","GO:0007010","GO:0006950","GO:0016043","GO:0044085","GO:0071840","GO:0007165","GO:0023052","GO:0007154","GO:0051716","GO:0050896","GO:0009987","GO:0050794","GO:0008150","UNCLASSIFIED","GO:0050896","GO:0009987","GO:0050794","GO:0008150","UNCLASSIFIED")

simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/pyrr_fst_GOterms_pantherslim_nocorrection.scatter.png", width = 7, height =7, units = "in", res = 300)
scatterPlot(simMatrix, reducedTerms)
dev.off()

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/pyrr_fst_GOterms_pantherslim_nocorrection.tree.png", width = 7, height =7, units = "in", res = 300)

treemapPlot(reducedTerms)
dev.off()




# raisd 
import pandas as pd
import numpy as np

urban_noca = pd.read_csv('relevantgenenames_noca_urban_raisd.txt')
rural_noca = pd.read_csv('relevantgenenames_noca_rural_raisd.txt')
np.intersect1d(urban_noca, rural_noca)

urban_noca_df = np.array(urban_noca)
rural_noca_df = np.array(rural_noca)

dif1_noca = np.setdiff1d(urban_noca_df, rural_noca_df)
dif2_noca = np.setdiff1d(rural_noca_df, urban_noca_df)

print(list(dif1_noca))
print(list(dif2_noca))

urban_pyrr = pd.read_csv('relevantgenenames_pyrr_urban_raisd.txt')
rural_pyrr = pd.read_csv('relevantgenenames_pyrr_rural_raisd.txt')
np.intersect1d(urban_pyrr, rural_pyrr)

urban_pyrr_df = np.array(urban_pyrr)
rural_pyrr_df = np.array(rural_pyrr)

dif1_pyrr = np.setdiff1d(urban_pyrr_df, rural_pyrr_df)
dif2_pyrr = np.setdiff1d(rural_pyrr_df, urban_pyrr_df)
list(set.intersection(dif1_pyrr, dif1_noca))
print(np.intersect1d(dif1_pyrr, dif1_noca))

print(list(dif1_pyrr))
print(list(dif2_pyrr))



# unique to urban (noca)
Adgrg4,Ahcyb,Aifm3,Ap5s1,Apmap,Arhgef6,Arid1b,Arx_1,Asip,Atp9b,Barx1b,Bhlhe23,CARCAR_R00615,CARCAR_R04350,CARCAR_R05153,CARCAR_R06321,CARCAR_R08420,CARCAR_R09689,CARCAR_R14549,CARCAR_R15815,Cd247,Cd2ap,Cd40lg,Cdk17,Chid1,Cntnap5,Col17a1,Csgalnact1,Cst7,Ctnnd2_1,Dicer1,Dlg2,Ehd4,Fcgbp_0,Fcgbp_1,Gabra5,Gtf2f2,Ha1f_2,Htr3a_0,Ikzf1,Itgb3,Klhl5,Mboat1,Mcc,Mid1,Mrps28,Nkx28,Nol4_0,Olig3,Parp12_1,Pim1_0,Pole,Pou2f1,Ptprr,Rbmx,Rprml,Sh2d4a,Slc10a2,Slc25a30,Slc30a9,Slc7a11_0,Slk,Snx10,Sult1b1_0,Sult1d1,Syndig1,Sytl5,Tab2,Tbx20,Thsd7a,Tmem33,Tpd52,Tpt1,Ubl3,Ugt2c1,Usp28,Vgll3,Wdr19,Wnk2,Zbtb16a,Zc3h12d,Zscan20_1,Zscan32_1

# unique to rural (noca)
A2ml1_2, Acyp1, Ak9_0, Ambra1, Antxr2_0, Antxr2_1, Anxa7_0, Anxa7_1, Atp6v1d, CARCAR_R02215, CARCAR_R09899, CARCAR_R12248, CARCAR_R13196, CARCAR_R14918, Cab39l, Chp1, Chrm4, Clic3, Crim1, Csrp2, Ctnna1, Dst, Eed, Egr1, Eif2b2, Eif3f, Entpd2_0, Epcam, Fa2h, Fgf18, Fig4, Fndc3a, Fzd1, Gnptg, Gtf2a1, Hdac11, Hdac9, Hecw1, Hikeshi, Hnrnpa2b1, Hoxa10_1, Hoxa5, Hoxa6, Hoxa7, Hoxa9, Iars2_1, Il21r, Il4r, Il9r, Ino80, Kdm3b, Lcor, Lpp, Mbnl3, Mcat, Mlh3_0, Mlh3_1, Mlkl, Mpp5, Neb_7, Neb_8, Nfe2l3, Npm1, Oip5, Paxx, Pgf, Plcl2_0, Plekhg1, Ptprm_0, Rab19, Rab3gap2, Reep2, Rfwd3, Rgs7, Rimklb, Rpl18a, Rps6kl1, Sapcd2, Setdb2, Sipa1l2, Slc37a3, Slc8a3, Slit1, Smoc1, Ston2, Tmem141, Tnfrsf11a, Tsnax, Tspo_1, Tsr3, Ttll1, Ube2i, Vwf_1, Wdr59, Yd021_0, Ythdc2_2, Zbtb14, Zc2hc1c, Zdhhc17

# unique to urban (pyrr)
A2ml1_2, Acr_7, Arhgap27_1, CARCAR_R04345, CARCAR_R04347, CARCAR_R06009, CARCAR_R06010, CARCAR_R06691, CARCAR_R15075, Ccnd2, Chst1, Cnksr2_1, Col6a3, Cst7, Ctnnbl1, Derl2_0, Dhx33, Drd1, Fgf23, Fgf6, Fuk, Gabpa, Ha1f_2, Hoxa1, Hoxa2, Hoxa3, Hoxa5, Hoxa6, Hoxa7, Itpripl1_0, Itpripl1_1, Kcnc1, Kcnk16_0, Mia3_0, Mia3_1, Mis12, Nr0b1, Pcdh10_1, Pdpr, Pds5a, Phf20, Pim1_1, Plekhm1, Prl_0, Rfx4, Rimklb, Rusc2, Satb2, Slc12a7, Slc30a10, Slc8a3, Smoc1, St3gal2, Tmem263_1, Trmt61b, Ube2k, Uvrag, Vstm2l

# unique to rural (pyrr)
Abhd11, Ahcyb, Antxr1_1, Ap1g1, Asip, Atxn1l, Babam2, Bcl11b, Btbd11, CARCAR_R00003, CARCAR_R00027, CARCAR_R00109, CARCAR_R02090, CARCAR_R02917, CARCAR_R02918, CARCAR_R02919, CARCAR_R03448, CARCAR_R04578, CARCAR_R08420, CARCAR_R09899, CARCAR_R10282, CARCAR_R14885, CARCAR_R14906, CARCAR_R15737, CARCAR_R15738, Chd6, Clec18a, Cmtm3, Cmtm4, Dhodh, Dhx38, Dlg1_0, Dlg2, Dll1_1, Dync1li2, Eif1b, Elovl4_1, Ephb2, F13a1, Fbln7_1, Fgf18, Glg1, Guca1c, H101_0, H2eb1, Ift27, Ist1, Limch1, Lmx1b_1, Lnpk, Man2b2, Mcat, Mei1, Mtmr4, Nuak1, Nwd2, Phox2b, Prlh, Prvm, Rap1gap2, Rbp3, Rd3, Ric8b, Snu13, Snx10, Srl, Stxbp5l, Syt14, Tmem19, Traf5, Tspo_1, Ttll1, Wnk2, Xrcc6, Zdhhc1, Zfc3h1, Zkscan7_0, Znf271, Znf774, Znf821

# local analysis
install.packages('ggvenn')
# library(venn)
# venn(5, ilab=TRUE, zcolor = "style")

# post-introduction selection
library(ggvenn)

df_pyrr<-read.csv("relevantgenenames_pyrr_fst.txt", header=FALSE)
df_noca<-read.csv("relevantgenenames_noca_fst.txt", header=FALSE)
colnames(df_pyrr) <- "Gene"
colnames(df_noca) <- "Gene"

x <- list(
pyrr = df_pyrr$Gene,
noca = df_noca$Gene
)


pdf(file = "species.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("pyrr", "noca"),
    fill_color = c("#4EAFAF", "#FF817E"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()


both_df = intersect(df_pyrr$Gene, df_noca$Gene)


write.csv(both_df, file="intersection_fst.csv", row.names=FALSE)


# NOCA_RAISD
# make plots

# BiocManager::install("org.Gg.eg.db")
library(org.Gg.eg.db)
library(limma)
library(rrvgo)

# noca_raisd_uniqueurban_GOterms_pantherslim_nocorrection
go_list <- c("GO:0030183","GO:0070527","GO:0042113","GO:2000117","GO:0030168","GO:0015721","GO:0032438","GO:0097194","GO:0043113","GO:0035721","GO:0034109","GO:0030098","GO:0051965","GO:0071786","GO:0000077","GO:0042770","GO:0002695","GO:0050866","GO:1903131","GO:0006271","GO:0007599","GO:0007596","GO:0050650","GO:0031570","GO:0050817","GO:0060078","GO:0006024","GO:0006284","GO:0050728","GO:0001775","GO:0000075","GO:1901988","GO:0010948","GO:0007423","GO:0051960","GO:0031503","GO:0016579","GO:0045786","GO:0070646","GO:0120031","GO:0030031","GO:1901987","GO:0060271","GO:0051240","GO:0006869","GO:0044782","GO:0120036","GO:0030030","GO:0006468","GO:0065008","GO:0006357","GO:0016043","GO:0071840","GO:0006468","GO:0065008","GO:0006357","GO:0016043","GO:0071840")

simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/noca_raisd_uniqueurban_GOterms_pantherslim_nocorrection.scatter.png", width = 7, height =7, units = "in", res = 300)
scatterPlot(simMatrix, reducedTerms)
dev.off()

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/noca_raisd_uniqueurban_GOterms_pantherslim_nocorrection.tree.png", width = 7, height =7, units = "in", res = 300)

treemapPlot(reducedTerms)
dev.off()


# noca_raisd_allurban_GOterms_pantherslim_nocorrection
go_list <- c("GO:0043030","GO:0030183","GO:0007076","GO:0051402","GO:0006120","GO:0002695","GO:0050866","GO:0070527","GO:0042130","GO:0042113","GO:0006312","GO:2000117","GO:0030888","GO:0040001","GO:0030168","GO:0032981","GO:0015721","GO:0032438","GO:0097194","GO:0000281","GO:0043113","GO:0046330","GO:0035721","GO:0030261","GO:0051653","GO:0051293","GO:0002683","GO:0000077","GO:0042770","GO:0002694","GO:0031570","GO:0050865","GO:0043009","GO:0048568","GO:0001775","GO:0061640","GO:0000910","GO:0051960","GO:0000075","GO:0051241","GO:0044057","GO:0048598","GO:0000070","GO:0001501","GO:0051301","GO:1902850","GO:0000819","GO:0140014","GO:0009952","GO:1901988","GO:0010948","GO:0009790","GO:0007423","GO:0003002","GO:0031503","GO:0006915","GO:0016579","GO:0045786","GO:0007059","GO:2000026","GO:0120031","GO:0030031","GO:0030198","GO:0043062","GO:0045229","GO:0051640","GO:0060271","GO:0120036","GO:0030030","GO:0006468","GO:0006974","GO:0016310","GO:0006357","GO:0006355","GO:2001141","GO:0051252","GO:0019219","GO:0010468","GO:0010556","GO:0051171","GO:0031326","GO:0009889","GO:0080090","GO:0032502","GO:0031323","GO:0060255","GO:0016043","GO:0019222","GO:0071840","GO:0050794","GO:0050789","GO:0065007","GO:0008150","UNCLASSIFIED","GO:0044249","GO:1901576","GO:0009058","GO:0008150","UNCLASSIFIED","GO:0044249","GO:1901576","GO:0009058")


simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/noca_raisd_allurban_GOterms_pantherslim_nocorrection.scatter.png", width = 7, height =7, units = "in", res = 300)
scatterPlot(simMatrix, reducedTerms)
dev.off()

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/noca_raisd_allurban_GOterms_pantherslim_nocorrection.tree.png", width = 7, height =7, units = "in", res = 300)

treemapPlot(reducedTerms)
dev.off()


# noca_raisd_allrural_GOterms_pantherslim_nocorrection
go_list <- c("GO:0043030","GO:0030183","GO:0007076","GO:0051402","GO:0006120","GO:0002695","GO:0050866","GO:0070527","GO:0042130","GO:0042113","GO:0006312","GO:2000117","GO:0030888","GO:0040001","GO:0030168","GO:0032981","GO:0015721","GO:0032438","GO:0097194","GO:0000281","GO:0043113","GO:0046330","GO:0035721","GO:0030261","GO:0051653","GO:0051293","GO:0002683","GO:0000077","GO:0042770","GO:0002694","GO:0031570","GO:0050865","GO:0043009","GO:0048568","GO:0001775","GO:0061640","GO:0000910","GO:0051960","GO:0000075","GO:0051241","GO:0044057","GO:0048598","GO:0000070","GO:0001501","GO:0051301","GO:1902850","GO:0000819","GO:0140014","GO:0009952","GO:1901988","GO:0010948","GO:0009790","GO:0007423","GO:0003002","GO:0031503","GO:0006915","GO:0016579","GO:0045786","GO:0007059","GO:2000026","GO:0120031","GO:0030031","GO:0030198","GO:0043062","GO:0045229","GO:0051640","GO:0060271","GO:0120036","GO:0030030","GO:0006468","GO:0006974","GO:0016310","GO:0006357","GO:0006355","GO:2001141","GO:0051252","GO:0019219","GO:0010468","GO:0010556","GO:0051171","GO:0031326","GO:0009889","GO:0080090","GO:0032502","GO:0031323","GO:0060255","GO:0016043","GO:0019222","GO:0071840","GO:0050794","GO:0050789","GO:0065007","GO:0008150","UNCLASSIFIED","GO:0044249","GO:1901576","GO:0009058","GO:0008150","UNCLASSIFIED","GO:0044249","GO:1901576","GO:0009058")


simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/noca_raisd_allrural_GOterms_pantherslim_nocorrection.scatter.png", width = 7, height =7, units = "in", res = 300)
scatterPlot(simMatrix, reducedTerms)
dev.off()

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/noca_raisd_allrural_GOterms_pantherslim_nocorrection.tree.png", width = 7, height =7, units = "in", res = 300)

treemapPlot(reducedTerms)
dev.off()



# PYRR RAISD

# unique to urban (pyrr) genes
A2ml1_2,Acr_7,Arhgap27_1,CARCAR_R04345,CARCAR_R04347,CARCAR_R06009,CARCAR_R06010,CARCAR_R06691,CARCAR_R15075,Ccnd2,Chst1,Cnksr2_1,Col6a3,Cst7,Ctnnbl1,Derl2_0,Dhx33,Drd1,Fgf23,Fgf6,Fuk,Gabpa,Ha1f_2,Hoxa1,Hoxa2,Hoxa3,Hoxa5,Hoxa6,Hoxa7,Itpripl1_0,Itpripl1_1,Kcnc1,Kcnk16_0,Mia3_0,Mia3_1,Mis12,Nr0b1,Pcdh10_1,Pdpr,Pds5a,Phf20,Pim1_1,Plekhm1,Prl_0,Rfx4,Rimklb,Rusc2,Satb2,Slc12a7,Slc30a10,Slc8a3,Smoc1,St3gal2,Tmem263_1,Trmt61b,Ube2k,Uvrag,Vstm2l

# pyrr_raisd_uniqueurban_GOterms_pantherslim_nocorrection
go_list <- c('GO:0043030','GO:0051382','GO:0051383','GO:0045943','GO:0007064','GO:2000117','GO:0006356','GO:0009225','GO:1902808','GO:0098754','GO:0007062','GO:0009952','GO:0002695','GO:0009636','GO:0050866','GO:0003002','GO:0008543','GO:0043009','GO:0044344','GO:0071774','GO:0045931','GO:0006688','GO:0050728','GO:0007389','GO:0048568','GO:0008361','GO:0006040','GO:0044843','GO:0000082','GO:0031348','GO:0050727','GO:0065004','GO:2000045','GO:0030488','GO:1902806','GO:0002683','GO:0048598','GO:0000070','GO:0001501','GO:0000819','GO:0140014','GO:0009887','GO:0098813','GO:0051276','GO:0008284','GO:0010628','GO:0001934','GO:0043410','GO:0009790','GO:0031401','GO:0000280','GO:0001932','GO:0070848','GO:0071363','GO:0007059','GO:1903047','GO:0098771','GO:0048285','GO:0071824','GO:0071805','GO:0006813','GO:0031399','GO:0030334','GO:0022402','GO:2000145','GO:0040012','GO:0140694','GO:0042127','GO:0042325','GO:0007049','GO:0042327','GO:0043408','GO:0010562','GO:0045937','GO:0051174','GO:0019220','GO:0055080','GO:0000278','GO:0050801','GO:0065003','GO:0051246','GO:0019725','GO:1902533','GO:0098662','GO:0098655','GO:0048513','GO:0030001','GO:0098660','GO:0043933','GO:0006812','GO:0006357','GO:0009653','GO:0034220','GO:0006355','GO:2001141','GO:0010557','GO:0010468','GO:0031328','GO:0051171','GO:0009891','GO:0010556','GO:0080090','GO:0031326','GO:0009889','GO:0051252','GO:0031323','GO:0060255','GO:0019219','GO:0019222','GO:0032502','GO:0048856','GO:0007275','GO:0065007','GO:0050794','GO:0050789','GO:0008150','UNCLASSIFIED','GO:0065007','GO:0050794','GO:0050789','GO:0008150','UNCLASSIFIED')

simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/pyrr_raisd_uniqueurban_GOterms_pantherslim_nocorrection.scatter.png", width = 7, height =7, units = "in", res = 300)
scatterPlot(simMatrix, reducedTerms)
dev.off()

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/pyrr_raisd_uniqueurban_GOterms_pantherslim_nocorrection.tree.png", width = 7, height =7, units = "in", res = 300)

treemapPlot(reducedTerms)
dev.off()

# pyrr_raisd_uniqueurban_GOterms_pantherslim_FDR
go_list <- c('GO:0051382','GO:0051383','GO:0007064','GO:0009952','GO:0003002','GO:0007389','GO:0009887','GO:0051276','GO:0022402','GO:0006357','GO:0006355','GO:2001141','GO:0010468','GO:0051171','GO:0010556','GO:0080090','GO:0031326','GO:0009889','GO:0051252','GO:0031323','GO:0060255','GO:0019219','GO:0019222')


simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/pyrr_raisd_uniqueurban_GOterms_pantherslim_FDR.scatter.png", width = 7, height =7, units = "in", res = 300)
scatterPlot(simMatrix, reducedTerms)
dev.off()

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/pyrr_raisd_uniqueurban_GOterms_pantherslim_FDR.tree.png", width = 7, height =7, units = "in", res = 300)

treemapPlot(reducedTerms)
dev.off()


# pyrr_raisd_allurban_GOterms_pantherslim_FDR
go_list <- c('GO:0051382','GO:0051383','GO:0009952','GO:0003002','GO:0007389','GO:0006357','GO:0006355','GO:2001141','GO:0010468','GO:0051171','GO:0010556','GO:0031326','GO:0009889','GO:0080090','GO:0031323','GO:0060255','GO:0050794','GO:0065007','GO:0050789')

simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/pyrr_raisd_allurban_GOterms_pantherslim_FDR.scatter.png", width = 7, height =7, units = "in", res = 300)
scatterPlot(simMatrix, reducedTerms)
dev.off()

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/pyrr_raisd_allurban_GOterms_pantherslim_FDR.tree.png", width = 7, height =7, units = "in", res = 300)

treemapPlot(reducedTerms)
dev.off()

# pyrr_raisd_allurban_GOterms_pantherslim_nocorrection

go_list <- c('GO:0043030','GO:0051382','GO:0051383','GO:0051402','GO:0099565','GO:0045943','GO:0007064','GO:0002695','GO:0050866','GO:0070527','GO:0042130','GO:0006312','GO:2000117','GO:0030888','GO:0030168','GO:0006356','GO:0009225','GO:0046330','GO:1903038','GO:1902808','GO:0098754','GO:0034109','GO:0050868','GO:0007062','GO:0009952','GO:0042129','GO:0002683','GO:0051965','GO:0051250','GO:0006298','GO:0009636','GO:0003002','GO:0002694','GO:0046328','GO:0008543','GO:0050865','GO:0043009','GO:0044344','GO:0071774','GO:0007389','GO:0048568','GO:0065004','GO:0051960','GO:0051241','GO:0044057','GO:0048598','GO:0000070','GO:0001501','GO:0000819','GO:0140014','GO:0009887','GO:0043410','GO:0042127','GO:0070848','GO:0071363','GO:0007059','GO:0098813','GO:0051276','GO:2000026','GO:0071805','GO:0006813','GO:0008284','GO:0042391','GO:0010628','GO:0001934','GO:0009790','GO:0007189','GO:0031401','GO:0043408','GO:0022402','GO:0001932','GO:1903047','GO:0071824','GO:1902533','GO:0007049','GO:0031399','GO:0010647','GO:0023056','GO:0006281','GO:0098662','GO:0042325','GO:0098655','GO:0051174','GO:0019220','GO:0030001','GO:0098660','GO:0006812','GO:0051246','GO:0034220','GO:0006357','GO:0010557','GO:0031328','GO:0009891','GO:0065008','GO:0006355','GO:2001141','GO:0007166','GO:0010468','GO:0051171','GO:0010556','GO:0031326','GO:0009889','GO:0080090','GO:0051252','GO:0019219','GO:0031323','GO:0060255','GO:0048522','GO:0048856','GO:0032502','GO:0019222','GO:0048518','GO:0050794','GO:0050896','GO:0051716','GO:0065007','GO:0050789','GO:0007165','GO:0023052','GO:0007154','GO:0008150','UNCLASSIFIED','GO:0007165','GO:0023052','GO:0007154','GO:0008150','UNCLASSIFIED')

simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/pyrr_raisd_allurban_GOterms_pantherslim_nocorrection.png", width = 7, height =7, units = "in", res = 300)
scatterPlot(simMatrix, reducedTerms)
dev.off()

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/pyrr_raisd_allurban_GOterms_pantherslim_nocorrection.tree.png", width = 7, height =7, units = "in", res = 300)

treemapPlot(reducedTerms)
dev.off()

# pyrr_raisd_allrural_GOterms_pantherslim_nocorrection
go_list <- c('GO:0051402','GO:0099565','GO:0072376','GO:0006896','GO:0070527','GO:0042130','GO:0007599','GO:0007596','GO:0006312','GO:0050817','GO:0030888','GO:0030168','GO:0032438','GO:0043113','GO:0006303','GO:0046330','GO:1903038','GO:0034109','GO:0050878','GO:0050868','GO:0042129','GO:0051965','GO:0071786','GO:0051250','GO:0006298','GO:0006206','GO:0042060','GO:0000723','GO:0032200','GO:0009611','GO:0006892','GO:0051960','GO:0006302','GO:0042127','GO:0006468','GO:0016310','GO:0006281','GO:0051668','GO:0006974','GO:0007167','GO:0006259','GO:0007166','GO:0065008','GO:0006796','GO:0006793','GO:0051716','GO:0007165','GO:0050896','GO:0023052','GO:0007154','GO:0050794','GO:0065007','GO:0050789','GO:0023052','GO:0007154','GO:0050794','GO:0065007','GO:0050789')

simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/pyrr_raisd_allrural_GOterms_pantherslim_nocorrection.png", width = 7, height =7, units = "in", res = 300)
scatterPlot(simMatrix, reducedTerms)
dev.off()

png(file = "/Users/danjack/Documents/CardinalGenomics/Paper/GeneListAnalyses/pyrr_raisd_allrural_GOterms_pantherslim_nocorrection.tree.png", width = 7, height =7, units = "in", res = 300)

treemapPlot(reducedTerms)
dev.off()