#Draw violinplot by scCustomize package
library(scCustomize)
library(qs)

####################### FIBROBLAST MARKERS ####
#annotation platform
#Fascia | Disease only
markers.to.plot <- c("ITGA10", "CCN3", "DPP4", "CDH13", "PRG4", "CRTAC1", "PCOLCE2", "LGR5")
Annot_plat_Fascia <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Fascia") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/062525_FB_Fascia.pdf", width = 15, height = 15)
Annot_plat_Fascia
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

p <- Stacked_VlnPlot(
  seurat_object = seurat_merged,
  features = c("IBSP","ITGAV","ITGB5"),
  x_lab_rotate = TRUE,
  #colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6")
)

# Force collection of all legends and show on the right
p + plot_layout(guides = "collect") & theme(legend.position = "right")

#F1 Superficial | Healthy same as Disease but disease has no COL13A1
markers.to.plot <- c("APCDD1", "COL18A1", "COL23A1", "COL13A1", "COMP", "NKD2", "RSPO1", "AXIN2")
Annot_plat_Superficial <- DotPlot(TN.combined_DF, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Superficial Fibroblasts") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/061825_All_Superficial.pdf", width = 15, height = 15)
Annot_plat_Superficial
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))


#F1 CRABP1 | Disease only
markers.to.plot <- c("CRABP1", "CYP26B1", "TNFRSF21", "CXCL1", "WNT5A", "COL18A1",
                     "COL23A1", "COL13A1", "NKD2", "AXIN2", "RSPO1")
Annot_plat_CRABP <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Superficial CRABP1+") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/062525_FB_Sup_CRABP1.pdf", width = 15, height = 15)
Annot_plat_CRABP
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))


#F2 Universal
markers.to.plot <- c("CD34", "PI16", "MFAP5", "DPP4", "PCOLCE2", "LGR5", "SLPI", "CD70", "CTHRC1")
Annot_plat_Universal <- DotPlot(TN.combined_DF, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Universal Fibroblasts") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/061825_All_UniversalFB.pdf", width = 15, height = 15)
Annot_plat_Universal
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

markers.to.plot <- c("PI16", "DPP4", "PCOLCE2", "MFAP5", "CD70", "LGR5")
Annot_plat_Universal_Dis <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Universal Fibroblasts") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/071125_FB_UniversalFB_Dis.pdf", width = 15, height = 15)
Annot_plat_Universal_Dis
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))


markers.to.plot <- c("FAP", "SPP1", "IL6", "COL1A1", "COL2A1", "SOX9", "RUNX2", "ACAN", "TNMD")

#F2/3 Perivascular | Healthy same as Disease but disease has no PPARG
markers.to.plot <- c("CXCL12", "APOE", "EFEMP1", "APOC1", "C7", "PLA2G2A", "PPARG", "MYOC", "GDF10")
Annot_plat_Perivascular <- DotPlot(TN.combined_DF, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Perivascular (Universal) Fibroblasts") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/061825_All_Peri.pdf", width = 15, height = 15)
Annot_plat_Perivascular
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

markers.to.plot <- c("CXCL12", "APOE", "C7", "PLA2G2A", "EFEMP1", "GDF10", "MYOC")
Annot_plat_Perivascular_Dis <- DotPlot(TN.combined_DF, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Perivascular (Universal) Fibroblasts") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/Markers/070125_All_Dis_Peri.pdf", width = 15, height = 15)
Annot_plat_Perivascular_Dis
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))


#F3 FRC-like | CCL19+
markers.to.plot <- c("CCL19", "CD74", "CH25H", "TNFSF13B", "IL33",
                     "IRF8", "IL15", "VCAM1", "HLA-DRB1", "HLA-DRA")
Annot_plat_Ret_CCL19 <- DotPlot(TN.combined_DF, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="CCL19+ Fibroblast reticular-like cell") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/061825_All_Ret_CCL19.pdf", width = 15, height = 15)
Annot_plat_Ret_CCL19
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

markers.to.plot <- c("CCL19", "CD74", "CH25H", "TNFSF13B", "IL33",
                     "HLA-DRA", "IRF8", "COX4I2", "RBP5", "ADAMDEC1",
                     "CXCL9", "CXCL10", "APOE", "CXCL12")
Annot_plat_Ret_CCL19_Dis <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="CCL19+ Fibroblast reticular-like cell") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/Markers/070125_FB_Dis_Ret_CCL19.pdf", width = 15, height = 15)
Annot_plat_Ret_CCL19_Dis
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

#F4 HF DPEP
markers.to.plot <- c("DPEP1", "MYL4", "MEF2C", "COL11A1")
Annot_plat_HF_DPEP <- DotPlot(TN.combined_DF, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Dermal sheath_DPEP1+ (hair-follicle associated)") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/061825_All_HF_DPEP.pdf", width = 15, height = 15)
Annot_plat_HF_DPEP
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

markers.to.plot <- c("MEF2C", "MYL4", "COL11A1", "POSTN", "DPEP1")
Annot_plat_HF_DPEP_Dis <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Dermal sheath_DPEP1+ (hair-follicle associated)") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/Markers/070125_FB_Dis_HF_DPEP.pdf", width = 15, height = 15)
Annot_plat_HF_DPEP_Dis
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

#F4  HF TNN
markers.to.plot <- c("TNN", "COCH", "TNMD", "MKX", "NRG3", "SLITRK6")
Annot_plat_HF_TNN <- DotPlot(TN.combined_DF, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Hair-follicle associated (TNN+COCH+)") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/061825_All_HF_TNN.pdf", width = 15, height = 15)
Annot_plat_HF_TNN
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

markers.to.plot <- c("COCH", "CRABP1", "COL24A1", "RSPO4", "SLITRK6", "NRG3", "MKX", "TNMD")
Annot_plat_HF_TNN_Dis <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Hair-follicle associated (TNN+COCH+)") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/Markers/070125_FB_Dis_HF_TNN.pdf", width = 15, height = 15)
Annot_plat_HF_TNN_Dis
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

#F4 HF DP | HHIP+
markers.to.plot <- c("CORIN", "HHIP", "BMP7", "WNT5A", "LEF1")
Annot_plat_HF_DP <- DotPlot(TN.combined_DF, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Hair-follicle associated (dermal papilla)") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/061825_All_HF_DP.pdf", width = 15, height = 15)
Annot_plat_HF_DP
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))


markers.to.plot <- c("CRABP1", "COL24A1", "RSPO4", "RSPO3", "BMP7", "WNT5A",
                     "LEF1", "SOX18", "HHIP")
Annot_plat_HF_DP_Dis <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Hair-follicle associated (dermal papilla)") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/Markers/070125_FB_Dis_HF_DP.pdf", width = 15, height = 15)
Annot_plat_HF_DP_Dis
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))


#F5 Schwann-like Fibroblast | NGFR+
markers.to.plot <- c("NGFR", "ITGA6", "SCN7A", "CDH19", "CLDN1", "SFRP4", "TENM2")
Annot_plat_SC_NGFR <- DotPlot(TN.combined_DF, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Schwann-like fibroblast (NGFR+)") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/061825_All_SC_NGFR.pdf", width = 15, height = 15)
Annot_plat_SC_NGFR
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))


markers.to.plot <- c("NGFR", "TM4SF1", "SFRP4", "ANGPTL7", "ITGA6", "CDH19", "CLDN1", "EBF2", "OLFM2", "SCN7A")
Annot_plat_SC_NGFR_Dis <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Schwann-like fibroblast (NGFR+)") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/Markers/070125_FB_Dis_SC_NGFR.pdf", width = 15, height = 15)
Annot_plat_SC_NGFR_Dis
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))


#F5 Schwann-like Fibroblast | RAMP1+ | Almost same byt SCN7A (H) and COL26A1&FMO2 (D)
markers.to.plot <- c("RAMP1", "RELN", "PLEKHA6", "IGFBP2", "FGFBP2", "SCN7A")
Annot_plat_SC_RAMP1 <- DotPlot(TN.combined_DF, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Schwann-like fibroblast (RAMP1+)") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/061825_All_SC_RAMP1.pdf", width = 15, height = 15)
Annot_plat_SC_RAMP1
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

markers.to.plot <- c("RAMP1", "IGFBP2", "RELN", "COL26A1", "PLEKHA6", "FMO2", "FGFBP2")
Annot_plat_SC_RAMP1_Dis <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Schwann-like fibroblast (RAMP1+)") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/Markers/070125_FB_Dis_SC_RAMP1.pdf", width = 15, height = 15)
Annot_plat_SC_RAMP1_Dis
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

#F6 | Inflammatory Myofibroblasts
markers.to.plot <- c("IL11", "IL24", "CXCL5", "CXCL6", "CXCL8", "MMP9", "WNT2", "COL10A1", "MMP1", "MMP3")
Annot_plat_In_Myo <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Inflammatory Myofibroblast") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/062525_FB_Myo.pdf", width = 15, height = 15)
Annot_plat_In_Myo
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))


#F7 | Universal Myofibroblasts
markers.to.plot <- c("ACTA2", "TAGLN", "CTHRC1", "RUNX2", "KIF26B", "SULF1", "ADAM12", "COL8A1", "LRRC15",
                     "CCN4", "ASPN", "POSTN", "TNC", "COL3A1", "WNT2", "COL10A1")
Annot_plat_Myo <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Myofibroblast") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/062525_FB_MyoFB.pdf", width = 15, height = 15)
Annot_plat_Myo
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

#F8 | Fascia-like Myofibroblasts
markers.to.plot <- c("ACAN", "ITGA10", "CDH2", "DPP4", "CCN3")
Annot_plat_FL_Myo <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Fascia-like Myofibroblast") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/062525_FB_Fascia_Myo.pdf", width = 15, height = 15)
Annot_plat_FL_Myo
dev.off()

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))


markers.to.plot <- c("COL23A1", "MFAP5", "POSTN", "CCL19")
Annot_PPRX1_Douetal_2025 <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Fibroblast Markers") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/062625_FB_markers.pdf", width = 15, height = 15)
Annot_PPRX1_Douetal_2025
dev.off()

markers.to.plot <- c("IL4R", "IL13RA1", "COL1A1", "FN1", "CCL19", "APOE", "CTHRC1", "ANGPTL1", "TAGLN", "ACTA2", "APCDD1", "COL18A1", "ITGA6", "ITGB4")
Annot_FB_KGD_markers_2025 <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Fibroblast Markers") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/062625_FB_KGD_markers.pdf", width = 15, height = 15)
Annot_FB_KGD_markers_2025
dev.off()


########### Other Genes ####

markers.to.plot <- c("COL10A1", "COL11A1", "COL12A1", "COL14A1", "COL16A1", "ASPN")
Chondro_Osteo_M <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Chondro/Osteogenic Markers") + theme(plot.title = element_text(hjust = 0.5, size=24))

pdf("D:/Jojie/Analysis_Outputs/scRNA-20250617/061825_Chondro_Osteo_M.pdf", width = 15, height = 15)
Chondro_Osteo_M
dev.off()

markers.to.plot <- c("IGF1R")
IGFR <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="IGF1R in Cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
IGFR

FeaturePlot(FB.subgroup, features = "IGF1R",
            combine = T,cols =  {c("gray", "darkred")},  pt.size = 0.6, split.by = "orig.ident1")


markers.to.plot <- c("DCC", "UNC5C", "FIBIN")
Bio_Rebio <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Mechanosensitive Genes") + theme(plot.title = element_text(hjust = 0.5, size=24))
Bio_Rebio


Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

Stacked_VlnPlot(seurat_object = FB.subgroup, features = c("CD74", "PIEZO1", "PIEZO2"), x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

DotPlot(FB.subgroup, features = c("CD74", "PIEZO1", "PIEZO2"), cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="CD74 in Fibroblasts") + theme(plot.title = element_text(hjust = 0.5, size=24))

FeaturePlot(FB.subgroup, features = c("TXNDC5"),
            cols = c("lightgrey", "brown3"), split.by = "orig.ident1")

FeaturePlot(TN.combined_DF, features = c("CCN3"),
            cols = c("lightgrey", "brown3"), split.by = "orig.ident1")

Stacked_VlnPlot(seurat_object = FB.subgroup, features = "TXNDC5", x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

DotPlot(FB.subgroup, features = "TXNDC5", cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="TXNDC5 in Fibroblasts") + theme(plot.title = element_text(hjust = 0.5, size=24))


########### Vol Plot #####
fib_healthy <- subset(FB.subgroup, subset = orig.ident1 == "1_Healthy_skin")
fib_adj <- subset(FB.subgroup, subset = orig.ident1 == "2_Adjacent_skin")
fib_keloid <- subset(FB.subgroup, subset = orig.ident1 == "3_Keloid_skin")
fib_reg <- subset(FB.subgroup, subset = orig.ident1 == "4_Regression_skin")

fib_H_v_K <- merge(fib_healthy, fib_keloid)
fib_A_v_K <- merge(fib_adj, fib_keloid)
Idents(fib_H_v_K) <- fib_H_v_K$orig.ident1
Idents(fib_A_v_K) <- fib_A_v_K$orig.ident1

deg <- FindMarkers(FB.subgroup, ident.1 = 1, ident.2 = 3)
deg_up <- deg[deg$avg_log2FC > 0.25 & deg$p_val_adj < 0.05 & deg$pct.1 > 0.1, ]
deg_down <- deg[deg$avg_log2FC < -0.25 & deg$p_val_adj < 0.05 & deg$pct.2 > 0.1, ]
write.csv(deg_up, "D:/Jojie/Analysis_Outputs/scRNA-20250617/DEGup_1v3.csv")
write.csv(deg_down, "D:/Jojie/Analysis_Outputs/scRNA-20250617/DEGdown_1v3.csv")

deg_clean <- deg
deg_clean$p_val_adj[deg_clean$p_val_adj == 0] <- min(deg_clean$p_val_adj[deg_clean$p_val_adj > 0], na.rm = TRUE) * 1e-1
deg_clean <- deg_clean[complete.cases(deg_clean) & is.finite(deg_clean$avg_log2FC), ]


library(EnhancedVolcano)
EnhancedVolcano(deg_clean,
                lab = rownames(deg),
                x = "avg_log2FC",
                y = "p_val_adj",
                xlim = c(-1.5,1.5),
                pCutoff = 0.05,
                FCcutoff = 0.25,
                title = "Cluster 1 vs Cluster 3",
                pointSize = 2.5,
                labSize = 3.5)

EnhancedVolcano(deg_clean,
                lab = rownames(deg),
                x = "avg_log2FC",
                y = "p_val_adj",
                xlim = c(-2,2),
                pCutoff = 0.05,
                FCcutoff = 0.25,
                title = "Adjacent vs Keloid",
                subtitle = NULL,
                pointSize = 2.5,
                labSize = 3.5,
                selectLab = c("COL11A1", "COL12A1", "COL14A1", "COL16A1", "RUNX2", "ASPN"))

EnhancedVolcano(deg_clean,
                lab = rownames(deg),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.25,
                col = c('blue',     # downregulated
                        'gray',     # non-significant
                        'gray',     # non-significant
                        'red'),     # upregulated
                colAlpha = 0.7,
                pointSize = 2.0,
                labSize = 3.5,
                title = "Healthy vs Keloid",
                subtitle = "Volcano Plot",
                caption = "log2FC cutoff = 0.25; adj p < 0.05")





####### Module Scores #####

exprMatrix <- GetAssayData(FB.subgroup, assay = "RNA", slot = "data")

migration_genes <- c("MMP2", "MMP9", "VIM", "SNAI1", "FN1") #change
geneSets <- list(migration_genes = migration_genes)
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 1, plotStats = FALSE)

# Calculate AUC scores
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

# Add AUC score to Seurat object metadata
FB.subgroup$Migration_AUCell <- as.numeric(getAUC(cells_AUC)["migration_genes", ])

FeaturePlot(FB.subgroup, features = "Migration_AUCell")
VlnPlot(FB.subgroup, features = "Migration_AUCell")
VlnPlot(FB.subgroup, features = "Migration_AUCell", split.by = "orig.ident1",
        pt.size = FALSE, cols = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

####### Also Module Scores ####

library(org.Hs.eg.db)

go_term_angiogenesis <- "GO:0001525"  # Angiogenesis
go_term_migration <- "GO:0016477"     # Cell migration
go_term_apoptosis <- "GO:0008626"     # Granzyme
go_term_migration2 <- "GO:0030334"    # Positive reg of migration
go_term_migration3 <- "GO:0006929"    # Substrate-dependent cell migration

angiogenesis_genes <- bitr(go_term_angiogenesis, fromType = "GO", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
migration1_genes <- bitr(go_term_migration, fromType = "GO", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
migration2_genes <- bitr(go_term_migration2, fromType = "GO", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
migration2_genes <- bitr(go_term_migration3, fromType = "GO", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
apoptosis_genes <- bitr(go_term_apoptosis, fromType = "GO", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL

angiogenesis_genes <- intersect(angiogenesis_genes, rownames(FB.subgroup))
migration1_genes <- intersect(migration1_genes, rownames(FB.subgroup))
migration2_genes <- intersect(migration2_genes, rownames(FB.subgroup))
migration3_genes <- intersect(migration2_genes, rownames(FB.subgroup))
apoptosis_genes <- intersect(apoptosis_genes, rownames(FB_Clusters_2))

FB.subgroup <- AddModuleScore(
  object = FB.subgroup,
  features = list(migration1_genes, migration2_genes, migration3_genes),
  name = c("Migration_Score", "Regulation_of_Migration_Score", "Substrate_Dependent_Migration_Score")
)

VlnPlot(FB.subgroup, features = c("Substrate_Dependent_Migration_Score"), pt.size = 0)

#AUCell
#prepare data
FB.subgroup.norm <- NormalizeData(FB.subgroup, normalization.method = "RC", scale.factor = 1e6)
expression_data <- as.matrix(FB.subgroup.norm@assays$RNA@data)

#prepare genes
geneSets <- list(Angiogenesis = angiogenesis_genes)
geneSets <- list(
  Angiogenesis = c("VEGFA", "ANGPT1", "PDGFB")
)

# Create AUCell rankings
cells_rankings <- AUCell_buildRankings(expression_data, nCores = 1, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
auc_scores <- as.data.frame(getAUC(cells_AUC))
colnames(auc_scores) <- make.names(colnames(auc_scores))
FB_Clusters_2 <- AddMetaData(FB_Clusters_2, metadata = auc_scores)
auc_scores <- FB_Clusters_2@meta.data$auc_scores
VlnPlot(seurat_obj, features = c("Angiogenesis_Score", "Migration_Score"), 
        pt.size = 0, geom = "boxplot")
VlnPlot(FB_Clusters_2, features = c("Angiogenesis_Score"))



########## PSEUDOTIME ####

FeaturePlot(FB.subgroup, features = c("MKI67", "PCNA", "MCM2"),
                         cols = c("lightgray", "darkred"), split.by = "orig.ident1")

FeaturePlot(FB.subgroup, features = c("S100B", "CCN3", "NGFR", "RAMP1"),
            cols = c("lightgray", "darkred"), split.by = "orig.ident1")

FeaturePlot(FB.subgroup, features = c("PLP1"),
            cols = c("lightgray", "darkred"), split.by = "orig.ident1")

markers.to.plot <- c("MKI67", "PCNA", "MCM2") #proliferating markers

library(monocle)

data <- as(as.matrix(FB.subgroup@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = FB.subgroup@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
FB.subgroup_mono2 <- newCellDataSet(data,
                                    phenoData = pd,
                                    featureData = fd,
                                    lowerDetectionLimit = 0.25,
                                    expressionFamily = negbinomial.size())
FB.subgroup_mono2
pData(FB.subgroup_mono2)
fData(FB.subgroup_mono2)

#Estimate size factors and dispersions
FB.subgroup_mono2 <- estimateSizeFactors(FB.subgroup_mono2)
FB.subgroup_mono2 <- estimateDispersions(FB.subgroup_mono2)

#Filtering low-quality cells 
FB.subgroup_mono2 <- detectGenes(FB.subgroup_mono2, min_expr = 0.1)
print(head(pData(FB.subgroup_mono2)))
expressed_genes <- row.names(subset(fData(FB.subgroup_mono2),num_cells_expressed >= 5))

L <- log(exprs(FB.subgroup_mono2[expressed_genes,]))
mL <- apply(L,1,function(x){mean(x[is.finite(x)])})
sdL <- apply(L,1,function(x){sd(x[is.finite(x)])})
Lstd <- (L-mL)/sdL

library(reshape)
#melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
#                                 ^^^^^

melted_dens_df <- melt(as.matrix(Lstd))
#Warning messages:
#  1: In type.convert.default(X[[i]], ...) :
#  'as.is' should be specified by the caller; using TRUE
#2: In type.convert.default(X[[i]], ...) :
#  'as.is' should be specified by the caller; using TRUE

qplot(value, geom = "density", data = melted_dens_df) + stat_function(fun = dnorm, size = 0.5, color = 'red') +xlab("Standardized log(FPKM)") + ylab("Density")

#Get more RAM
gc()
memory.limit(9999999999)
gc()

# Step 3: Select ordering genes using HVG
ordering_genes <- VariableFeatures(FB.subgroup)

# Confirm overlap with monocle_cds genes
ordering_genes <- ordering_genes[ordering_genes %in% row.names(fData(FB.subgroup_mono2))]

# Set ordering filter
FB.subgroup_mono2 <- setOrderingFilter(FB.subgroup_mono2, ordering_genes)

# Step 4: Dimensionality reduction and trajectory building

FB.subgroup_mono2 <- reduceDimension(FB.subgroup_mono2, max_components = 2, method = 'DDRTree')
FB.subgroup_mono2 <- orderCells(FB.subgroup_mono2)

# Find root using
table(pData(FB.subgroup_mono2)$State)

# Select own root cluster
FB.subgroup_mono2 <- orderCells(FB.subgroup_mono2, root_state = 3)


# Step 5: Visualize
plot_cell_trajectory(FB.subgroup_mono2, color_by = "Pseudotime")
plot_cell_trajectory(FB.subgroup_mono2, color_by = "State")
plot_cell_trajectory(FB.subgroup_mono2, color_by = "seurat_clusters")
plot_cell_trajectory(FB.subgroup_mono2, color_by = "orig.ident1")

plot_cell_trajectory(FB.subgroup_mono2, color_by = "seurat_clusters") +
  facet_wrap(~seurat_clusters, nrow = 3)

diff_test_res <- differentialGeneTest(FB.subgroup_mono2[expressed_genes, ],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_genes <- row.names(subset(diff_test_res, qval < 0.05))

plot_pseudotime_heatmap(FB.subgroup_mono2[sig_genes[1:50], ], #changed from 100 to 50
                        num_clusters = 5,
                        show_rownames = TRUE)



##### Gene changes by pseudotime branch or cluster #####

FB.subgroup_mono2 <- orderCells(FB.subgroup_mono2)
branch_1_cds <- FB.subgroup_mono2[, pData(FB.subgroup_mono2)$State == 1]

branch_1_cds <- setOrderingFilter(branch_1_cds, ordering_genes)  # use same genes as full CDS

diff_test_branch <- differentialGeneTest(branch_1_cds,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)")

# Top ranked genes by q-value
top_genes <- row.names(subset(diff_test_branch, qval < 0.05))
head(diff_test_branch[order(diff_test_branch$qval), 1:4])

write.csv(diff_test_branch, file = "D:/Jojie/Analysis_Outputs/scRNA-20250617/branch_1_genes.csv")

# Optional: plot top ones
plot_genes_in_pseudotime(branch_1_cds[top_genes[1:6], ])
plot_genes_in_pseudotime(branch_1_cds[c("APOD", "SELENOP", "CFD", "PTGDS", "SDC1", "APOE1"), ])
plot_genes_in_pseudotime(branch_1_cds[c("CCN3", "SFRP4", "NGFR", "RAMP1"), ])

plot_genes_in_pseudotime(branch_1_cds[c("CALB2"), ])


branch_2_cds <- FB.subgroup_mono2[, pData(FB.subgroup_mono2)$State == 2]

branch_2_cds <- setOrderingFilter(branch_2_cds, ordering_genes)  # use same genes as full CDS

diff_test_branch <- differentialGeneTest(branch_2_cds,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)")

# Top ranked genes by q-value
top_genes <- row.names(subset(diff_test_branch, qval < 0.05))
head(diff_test_branch[order(diff_test_branch$qval), 1:4])

write.csv(diff_test_branch, file = "D:/Jojie/Analysis_Outputs/scRNA-20250617/branch_2_genes.csv")

# Optional: plot top ones
plot_genes_in_pseudotime(branch_2_cds[top_genes[1:6], ])
plot_genes_in_pseudotime(branch_2_cds[c("APOD", "SELENOP", "CFD", "PTGDS", "SDC1", "APOE1"), ])
plot_genes_in_pseudotime(branch_2_cds[c("CCN3", "SFRP4", "NGFR", "RAMP1"), ])

plot_genes_in_pseudotime(branch_2_cds[c("CALB2"), ])


plot_genes_in_pseudotime(FB.subgroup_mono2[top_genes[1:6], ])

#### Slingshot Pseudotime ####
BiocManager::install("slingshot")

# Assuming FB.subgroup is already clustered and has UMAP
DefaultAssay(FB.subgroup) <- "RNA"

library(SingleCellExperiment)
library(slingshot)

FB.sce <- as.SingleCellExperiment(FB.subgroup)
# Check the UMAP name
reducedDimNames(FB.sce)  # should include "UMAP"

# Run slingshot using clusters and UMAP coordinates
FB.sce <- slingshot(FB.sce, clusterLabels = "seurat_clusters", reducedDim = "UMAP")

# Plot colored by clusters
plot(reducedDims(FB.sce)$UMAP, 
     col = brewer.pal(8, "Set2")[FB.sce$seurat_clusters],
     pch = 16, asp = 1,
     main = "Slingshot Trajectory on UMAP")

# Add trajectory lines and arrows
lines(SlingshotDataSet(FB.sce), lwd = 2, col = 'black', curveType = "lineage")


library(RColorBrewer)

# Basic UMAP scatter plot
plot(reducedDims(FB.sce)$UMAP,
     col = brewer.pal(8, "Set2")[FB.sce$seurat_clusters],
     pch = 16, asp = 1,
     main = "Slingshot Trajectory with Cluster Labels")

# Add trajectory lines/arrows
lines(SlingshotDataSet(FB.sce), lwd = 2, col = 'black')

# Get cluster centers
umap_coords <- reducedDims(FB.sce)$UMAP
cluster_labels <- colData(FB.sce)$seurat_clusters
centroids <- aggregate(umap_coords, by = list(cluster_labels), FUN = mean)

# Add text labels at cluster centroids
text(centroids[,2], centroids[,3], labels = centroids[,1], cex = 1, font = 2)

# Step 4: Extract pseudotime (optional)
pseudotime <- slingPseudotime(FB.sce)

# Add to Seurat object
FB.subgroup$pseudotime <- pseudotime[, 1]  # first lineage
FB.subgroup$pseudotime <- pseudotime[, 2]  # second lineage
FB.subgroup$pseudotime <- pseudotime[, 3]  # third lineage
FB.subgroup$pseudotime <- pseudotime[, 4]  # fourth lineage

FeaturePlot(FB.subgroup, features = "pseudotime", reduction = "umap")

FB.sce.2 <- slingshot(FB.sce, clusterLabels = "seurat_clusters", reducedDim = "UMAP", start.clus ="0")
pseudotime.2 <- slingPseudotime(FB.sce.2)
FB.subgroup$pseudotime.2 <- pseudotime.2[, 1]  # first lineage
FeaturePlot(FB.subgroup, features = "pseudotime.2", reduction = "umap")

# Basic UMAP scatter plot
plot(reducedDims(FB.sce.2)$UMAP,
     col = brewer.pal(8, "Set2")[FB.sce$seurat_clusters],
     pch = 16, asp = 1,
     main = "Slingshot Trajectory with Cluster Labels")

# Add trajectory lines/arrows
lines(SlingshotDataSet(FB.sce.2), lwd = 2, col = 'black')

# Get cluster centers
umap_coords <- reducedDims(FB.sce.2)$UMAP
cluster_labels <- colData(FB.sce.2)$seurat_clusters
centroids <- aggregate(umap_coords, by = list(cluster_labels), FUN = mean)

# Add text labels at cluster centroids
text(centroids[,2], centroids[,3], labels = centroids[,1], cex = 1, font = 2)


#### Create Loupe ####
library(loupeR)

create_loupe_from_seurat(
  FB.subgroup,
  output_dir = "D:/Jojie/Analysis_Outputs/scRNA-20250617/",
  output_name = "FB_subgroup_loupe",
  dedup_clusters = FALSE,
  feature_ids = NULL,
  executable_path = NULL,
  force = FALSE
)


########

saveRDS(FB.subgroup, file = "D:/Jojie/Analysis_Outputs/scRNA-20250617/072125_FBClusters.rds")


markers.to.plot <- c("TGFB1", "ITGB1", "COL1A1", "TRPC1", "TRPC3", "TRPC4", "TRPC5", "TRPC6", "PIEZO1", "PIEZO2")

markers.to.plot <- c("ITGAV", "ITGB5")
markers.to.plot <- c("POSTN", "ASPN", "MMP11", "MMP14", "MT-ND4", "MT-CO2", "MT-CYB")

Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))

markers.to.plot <- c("POSTN", "ASPN", "APCDD1", "WIF1", "WISP2", "SLPI", "CCL19", "APOE")
FB_general <- DotPlot(FB.subgroup, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="Fibroblast Classifications") + theme(plot.title = element_text(hjust = 0.5, size=24))
FB_general

markers.to.plot <- c("SFRP2", "DPP4", "FMO1", "LSP")
DotPlot(FB.subgroup, features = c("SFRP1","SYVN1"), cols = c("white", "darkred"), dot.scale = 8) + 
  RotatedAxis() + labs(title="SYVN and SFRP1") + theme(plot.title = element_text(hjust = 0.5, size=24))

markers.to.plot <- c("SFRP2", "DPP4", "FMO1", "LSP")
Stacked_VlnPlot(seurat_object = FB.subgroup, features = markers.to.plot, x_lab_rotate = FALSE,
                split.by = "orig.ident1", colors_use = c("#658354", "#D8B863", "#CC183c", "#92b1b6"))