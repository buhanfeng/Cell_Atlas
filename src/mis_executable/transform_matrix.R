

# R exe: /usr/local/bin/R
# R exe: /usr/local/bin/Rscript

# input from cmd:
# 1: transform type as tt, default NA
# 2: from data as fd, default NA
# 3: destination as dd, default NA

args=commandArgs(T)
print('[transform_matrix.R] is running...')
print(args[1])
print(args[2])
print(args[3])

tt = args[1]
fd = args[2]
dd = args[3]

 
library(generics)
library(Seurat)
library(dplyr)
library(patchwork)
library(reshape2)


crfmt2mat_wth = function(fd, dd){
  SRO.data = Read10X(data.dir = fd)
  SRO = CreateSeuratObject(counts = SRO.data, min.cells = 1, min.features = 1)
  all_barcode = SRO@assays$RNA@data@Dimnames[[2]]
  all_indent = SRO@active.ident
  if(is.null(all_barcode) | length(all_barcode)==0){
    return(message('barcode list is none'))
  }else{
    exp_matrix = as.matrix(SRO@assays$RNA@counts)
    exp_matrix = as.data.frame(exp_matrix)
    write.table(exp_matrix, file = dd, quote = F, sep = '\t')
  }
}


mat_wth2mat_lg = function(fd, dd){
  temp =  read.csv(fd, sep = '\t')
  temp$gene = rownames(temp)
  out = melt(temp, id="gene", variable.name="cell",value.name="expression")
  write.table(out, file = dd, quote = F, sep = '\t', row.names = F)
}


if(tt == 'crfmt2mat_wth'){
  crfmt2mat_wth(fd, dd)
}else if(tt == 'mat_wth2mat_lg'){
  mat_wth2mat_lg(fd, dd)
}



















