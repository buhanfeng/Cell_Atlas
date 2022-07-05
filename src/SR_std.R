

# R exe: /usr/local/bin/R
# R exe: /usr/local/bin/Rscript

# input from cmd:
# 1: data dir as data_dir, default '/root'
# 2: output dir as output_dir, default '/root'
# 3: resolution as relu, 0.2,0.05
# 4: logfc threshold as logfc_threshold, 1
# 5: min pct as min_pct, 0.2
data_dir = '/root'
output_dir = '/root'
relu = '0.2,0.05'

# output
# 1: vlnplot.pdf
# 2: var_features.csv
# 3: var_features.pdf
# 4: pca_st.pdf
# 5: JackStrawPlot.pdf
# 6: ElbowPlot.pdf
# 6: 
# 6: 
# 6: 
# 6: 

library(generics)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(magrittr)

args=commandArgs(T)
print(args[1])
print(args[2])
print(args[3])
print(args[4])

data_dir = args[1]
output_dir = args[2]
resolution = print(args[3])

vlnp_name = 'vlnplot.pdf'
vlnp_path = paste(output_dir, vlnp_name, sep = '/')
var_ft_name = 'var_features.csv'
var_ft_path = paste(output_dir, var_ft_name, sep = '/')
var_ft_p_name = 'var_features.pdf'
var_ft_p_path = paste(output_dir, var_ft_p_name, sep = '/')
pca_st_name = 'pca_st.pdf'
pca_st_path = paste(output_dir, pca_st_name, sep = '/')
JackStrawPlot_name = 'JackStrawPlot.pdf'
JackStrawPlot_path = paste(output_dir, JackStrawPlot_name, sep = '/')
ElbowPlot_name = 'ElbowPlot.pdf'
ElbowPlot_path = paste(output_dir, ElbowPlot_name, sep = '/')

#read 10X data
SRO.data = Read10X(data.dir = data_dir)

#create seurat object
SRO = CreateSeuratObject(counts = SRO.data, min.cells = 3, min.features = 200)

#QC
#list top ten gene
#get mt gene
SRO[["percent.mt"]] = PercentageFeatureSet(object = SRO, pattern = "^mt-")
#draw vlnplot
vp = VlnPlot(object = SRO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(vp ,filename = 'vp.pdf')

#obtain quanlity data
SRO = subset(x = SRO, subset = nFeature_RNA>200 & nFeature_RNA<12000 & percent.mt<25 & nCount_RNA<60000)
#normalization
#$SRO = NormalizeData(object = SRO, normalization.method = "LogNormalize", scale.factor = 10000)
SRO = NormalizeData(object = SRO)

#find high variable gene
SRO = FindVariableFeatures(object = SRO, selection.method = "vst", nfeatures = 2000)

# get variable features
variable_features = VariableFeatures(object = SRO)
vout = data.frame(variable_features)
write.csv(vout, var_ft_path, row.names = F, quote = F)

#Identity the 10 most highly variable genes
top10 = head(x = VariableFeatures(object = SRO), 10)

#plot variable features with and without labels
plot1 = VariableFeaturePlot(object = SRO)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
varp = CombinePlots(plots = list(plot1, plot2))
ggsave(varp ,filename = var_ft_p_path)

#PCA
all.gene = rownames(x = SRO)
#SRO = ScaleData(object = SRO, features = all.gene)
SRO = ScaleData(object = SRO)
#dimension reduction
SRO = RunPCA(object = SRO, features = VariableFeatures(object = SRO))
pcap = DimPlot(object = SRO, reduction = "pca")
ggsave(pcap ,filename = pca_st_path)

SRO = Run

#identify avarible dimension
SRO = JackStraw(object = SRO, num.replicate = 100)
SRO = ScoreJackStraw(object = SRO, dims = 1:20)
jackp = JackStrawPlot(object = SRO, dims = 1:20)
ggsave(jackp ,filename = JackStrawPlot_path)
elp = ElbowPlot(SRO)
ggsave(elp ,filename = ElbowPlot_path)

#cluster
#construct KNN map
SRO = FindNeighbors(object = SRO, dims = 1:20)

relu_list = strsplit(relu)
for(res in relu_list){
  res = as.numeric(res)
}

#process cluster, resolution determine by cell count
SRO = FindClusters(object = SRO, resolution = 0.2)
head(Idents(SRO), 5)
SRO = RunUMAP(object = SRO, dims = 1:20)
umapp = DimPlot(object = SRO, reduction = "umap")

SRO = RunTSNE(object = SRO, dims = 1:20)
tsnep = DimPlot(object = SRO, reduction = "tsne")


# get cluster data
cluster_data = SRO@active.ident
write.csv(cluster_data, 'D:\\desktop\\advance\\projects\\zebrafish_B_cell_and_Ig_gene_analysis\\3.5df_zibrafish_CHT_run1\\compare\\3.5df_zibrafish_CHT_2000_1-20_0d2_cluster_data.csv')

SRO.markers = FindAllMarkers(SRO, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 1)
SRO.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.csv(SRO.markers, 
          'D:\\desktop\\advance\\projects\\zebrafish_B_cell_and_Ig_gene_analysis\\3.5df_zibrafish_CHT_run1\\compare\\3.5df_zibrafish_CHT_2000_1-20_0d2_0d25_1_all_marker.csv')

FeaturePlot(object = SRO, features = c("lyz","ighz2"))

# do heat map
top10 = SRO.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
DoHeatmap(SRO, features = top10$gene)

all_barcode = SRO@assays$RNA@data@Dimnames[[2]]
all_indent = SRO@active.ident
if(is.null(all_barcode) | length(all_barcode)==0){
  return(message('barcode list is none'))
}else{
  exp_matrix = as.matrix(SRO@assays$RNA@counts)
  # exp_matrix = t(exp_matrix)
}
if(is.null(barcode_list)){
  if(is.null(indent)){
    barcode_list = all_barcode
  }else{
    barcode_list = c()
    for(ind in indent){
      barcode_list = c(barcode_list, all_barcode[all_indent == ind])
    }
  }
}
exp_matrix = exp_matrix[ , barcode_list]
return(exp_matrix)
