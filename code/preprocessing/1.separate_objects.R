# functions to check UMI count distribution per dataset
mcell_plot_umis_per_cell(mat_id_XXX)
plot_cum_umi_distr(mat_id_XXX)


############### SCHN_02 ##################
setwd("/DATA/project")

source("Functions_separate_objects.R")
tbk_reload()    # change dir names in functions_separate_objects file if needed

# load data
data = load("HN_full_dataset_scdb/ori_objects_online_only/GEX_VDJ_5821_GEX_S1.Rda")
mat_seurat_S1 = counts
data = load("HN_full_dataset_scdb/ori_objects_online_only/GEX_VDJ_5821_GEX_S2.Rda")
mat_seurat_S2 = counts

# export seurat object to a metacell object
sce_S1 = as.SingleCellExperiment(mat_seurat_S1)
sce_S2 = as.SingleCellExperiment(mat_seurat_S2)
mat_S1 = scm_import_sce_to_mat(sce_S1)
mat_S2 = scm_import_sce_to_mat(sce_S2)

mat_S1@cells = paste0('HN_02_S1_', mat_S1@cells)
rownames(mat_S1@cell_metadata) = mat_S1@cells
colnames(mat_S1@mat) = mat_S1@cells

mat_S2@cells = paste0('HN_02_S2_', mat_S2@cells)
rownames(mat_S2@cell_metadata) = mat_S2@cells
colnames(mat_S2@mat) = mat_S2@cells

mat_S1@cell_metadata[, "orig.ident"] = "HN_02_S1"
mat_S2@cell_metadata[, "orig.ident"] = "HN_02_S2"
mat_S1@cell_metadata$seq_run = "GCF-5821"
mat_S2@cell_metadata$seq_run = "GCF-5821"

scdb_add_mat(mat = mat_S1, id = mat_id_S1)
scdb_add_mat(mat = mat_S2, id = mat_id_S2)

mat_names = c(mat_id_S1, mat_id_S2)
for (i in mat_names){
  mat = scdb_mat(i)
  print(dim(mat@mat))
  
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO1"   |
                      mat@cell_metadata$hash.ID == "HTO2", "patient"] = "Pat10"
  
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO3"   |
                      mat@cell_metadata$hash.ID == "HTO4", "patient"] = "Pat15"
  
  
  
  mat@cell_metadata$HT = substr(mat@cell_metadata$hash.ID, 1, 5)
  mat@cell_metadata$new_hash_ID = mat@cell_metadata$hash.ID
  
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO1"   |
                      mat@cell_metadata$hash.ID == "HTO3", "condition"] = "pre"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO2"   |
                      mat@cell_metadata$hash.ID == "HTO4", "condition"] = "post"
  
  mat@cell_metadata[mat@cell_metadata$hash.ID == "Negative", "condition"] = "neg"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "Doublet", "condition"] = "doub"
  #mat@cell_metadata[mat@cell_metadata$condition == "-8h-aPD1", "condition"] = "8h-aPD1"
  
  scdb_add_mat(id = i,  mat = mat)
}

mat_S1 = scdb_mat(mat_id_S1)  
mat_S2 = scdb_mat(mat_id_S2)  






############### SCHN_04 ####################
setwd("/DATA/project")

# load data
data = load("HN_full_dataset_scdb/ori_objects_online_only/GEX_VDJ_5954.Rda")
mat_seurat_S1 = counts
#data = load("HN_02_scdb/5821_S2.Rda")
#mat_seurat_S2 = counts

# export seurat object to a metacell object
sce_S1 = as.SingleCellExperiment(mat_seurat_S1)
#sce_S2 = as.SingleCellExperiment(mat_seurat_S2)
mat_S1 = scm_import_sce_to_mat(sce_S1)
#mat_S2 = scm_import_sce_to_mat(sce_S2)

mat_S1@cells = paste0('HN_04_S1_', mat_S1@cells)
rownames(mat_S1@cell_metadata) = mat_S1@cells
colnames(mat_S1@mat) = mat_S1@cells

#mat_S2@cells = paste0('HN_02_S2_', mat_S2@cells)
#rownames(mat_S2@cell_metadata) = mat_S2@cells
#colnames(mat_S2@mat) = mat_S2@cells

mat_S1@cell_metadata[, "orig.ident"] = "HN_04_S1"
#mat_S2@cell_metadata[, "orig.ident"] = "HN_02_S2"
mat_S1@cell_metadata$seq_run = "GCF-5954"
#mat_S2@cell_metadata$seq_run = "GCF-5821"

scdb_add_mat(mat = mat_S1, id = mat_id_S1)
#scdb_add_mat(mat = mat_S2, id = mat_id_S2)

mat_names = c(mat_id_S1)
for (i in mat_names){
  mat = scdb_mat(i)
  print(dim(mat@mat))
  
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO1-NR-pre"   |
                      mat@cell_metadata$hash.ID == "HTO2-NR-post", "patient"] = "Pat30"
  
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO4-R-pre", "patient"] = "Pat22"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO8-R-post", "patient"] = "Pat17"
  
  
  mat@cell_metadata$HT = substr(mat@cell_metadata$hash.ID, 1, 5)
  mat@cell_metadata$new_hash_ID = mat@cell_metadata$hash.ID
  
  mat@cell_metadata$condition = substr(mat@cell_metadata$new_hash_ID, 8,  50)
  mat@cell_metadata[mat@cell_metadata$hash.ID == "Negative", "condition"] = "neg"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "Doublet", "condition"] = "doub"
  mat@cell_metadata[mat@cell_metadata$condition == "-pre", "condition"] = "pre"
  mat@cell_metadata[mat@cell_metadata$condition == "-post", "condition"] = "post"
  
  scdb_add_mat(id = i,  mat = mat)
}

mat_S1 = scdb_mat(mat_id_S1)  
#mat_S2 = scdb_mat(mat_id_S2)  




############## SCHN_05 ############

# load data
data = load("HN_full_dataset_scdb/ori_objects_online_only/GEX_VDJ_6009_S1.Rda")
mat_seurat_S1 = counts
data = load("HN_full_dataset_scdb/ori_objects_online_only/GEX_VDJ_6009_S2.Rda")
mat_seurat_S2 = counts

# export seurat object to a metacell object
sce_S1 = as.SingleCellExperiment(mat_seurat_S1)
sce_S2 = as.SingleCellExperiment(mat_seurat_S2)
mat_S1 = scm_import_sce_to_mat(sce_S1)
mat_S2 = scm_import_sce_to_mat(sce_S2)

mat_S1@cells = paste0('HN_05_S1_', mat_S1@cells)
rownames(mat_S1@cell_metadata) = mat_S1@cells
colnames(mat_S1@mat) = mat_S1@cells

mat_S2@cells = paste0('HN_05_S2_', mat_S2@cells)
rownames(mat_S2@cell_metadata) = mat_S2@cells
colnames(mat_S2@mat) = mat_S2@cells

mat_S1@cell_metadata[, "orig.ident"] = "HN_05_S1"
mat_S2@cell_metadata[, "orig.ident"] = "HN_05_S2"
mat_S1@cell_metadata$seq_run = "GCF-6009"
mat_S2@cell_metadata$seq_run = "GCF-6009"

scdb_add_mat(mat = mat_S1, id = mat_id_S1)
scdb_add_mat(mat = mat_S2, id = mat_id_S2)

mat_names = c(mat_id_S1, mat_id_S2)
for (i in mat_names){
  mat = scdb_mat(i)
  print(dim(mat@mat))
  
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO9-pre" |
                      mat@cell_metadata$hash.ID == "HTO10-post", "patient"] = "Pat31"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO1-pre"   |
                      mat@cell_metadata$hash.ID == "HTO2-post", "patient"] = "Pat34"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO7-pre"   |
                      mat@cell_metadata$hash.ID == "HTO8-post", "patient"] = "Pat29"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO4-pre"   |
                      mat@cell_metadata$hash.ID == "HTO6-post", "patient"] = "Pat37"
  
  
  
  mat@cell_metadata$HT = substr(mat@cell_metadata$hash.ID, 1, 5)
  mat@cell_metadata$new_hash_ID = mat@cell_metadata$hash.ID
  
  mat@cell_metadata$condition = substr(mat@cell_metadata$new_hash_ID, 6,  50)
  mat@cell_metadata[mat@cell_metadata$condition == "ive", "condition"] = "neg"
  mat@cell_metadata[mat@cell_metadata$condition == "et", "condition"] = "doub"
  mat@cell_metadata[mat@cell_metadata$condition == "-post", "condition"] = "post"
  
  scdb_add_mat(id = i,  mat = mat)
}

mat_S1 = scdb_mat(mat_id_S1)  
mat_S2 = scdb_mat(mat_id_S2)  




################# SCHN_06 ##############

# load data
data = load("HN_full_dataset_scdb/ori_objects_online_only/GEX_VDJ_6013.Rda")
mat_seurat_S1 = counts
#data = load("HN_05_scdb/GEX_VDJ_6009_S2.Rda")
#mat_seurat_S2 = counts

# export seurat object to a metacell object
sce_S1 = as.SingleCellExperiment(mat_seurat_S1)
#sce_S2 = as.SingleCellExperiment(mat_seurat_S2)
mat_S1 = scm_import_sce_to_mat(sce_S1)
#mat_S2 = scm_import_sce_to_mat(sce_S2)

mat_S1@cells = paste0('HN_06_S1_', mat_S1@cells)
rownames(mat_S1@cell_metadata) = mat_S1@cells
colnames(mat_S1@mat) = mat_S1@cells

#mat_S2@cells = paste0('HN_05_S2_', mat_S2@cells)
#rownames(mat_S2@cell_metadata) = mat_S2@cells
#colnames(mat_S2@mat) = mat_S2@cells

mat_S1@cell_metadata[, "orig.ident"] = "HN_06_S1"
#mat_S2@cell_metadata[, "orig.ident"] = "HN_05_S2"
mat_S1@cell_metadata$seq_run = "GCF-6013"
#mat_S2@cell_metadata$seq_run = "GCF-6009"

scdb_add_mat(mat = mat_S1, id = mat_id_S1)
#scdb_add_mat(mat = mat_S2, id = mat_id_S2)

mat_names = c(mat_id_S1)
for (i in mat_names){
  mat = scdb_mat(i)
  print(dim(mat@mat))
  
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO1-pre" |
                      mat@cell_metadata$hash.ID == "HTO2-post", "patient"] = "Pat38"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO4-pre"   |
                      mat@cell_metadata$hash.ID == "HTO6-post", "patient"] = "Pat12"
  
  mat@cell_metadata$HT = substr(mat@cell_metadata$hash.ID, 1, 5)
  mat@cell_metadata$new_hash_ID = mat@cell_metadata$hash.ID
  
  mat@cell_metadata$condition = substr(mat@cell_metadata$new_hash_ID, 6,  50)
  mat@cell_metadata[mat@cell_metadata$condition == "ive", "condition"] = "neg"
  mat@cell_metadata[mat@cell_metadata$condition == "et", "condition"] = "doub"
  #mat@cell_metadata[mat@cell_metadata$condition == "-8h-aPD1", "condition"] = "8h-aPD1"
  
  scdb_add_mat(id = i,  mat = mat)
}

mat_S1 = scdb_mat(mat_id_S1)  
#mat_S2 = scdb_mat(mat_id_S2)  





################# SCHN_07 ###############

# load data
data = load("HN_full_dataset_scdb/ori_objects_online_only/GEX_VDJ_6027.Rda")
mat_seurat_S1 = counts
#data = load("HN_05_scdb/GEX_VDJ_6009_S2.Rda")
#mat_seurat_S2 = counts

# export seurat object to a metacell object
sce_S1 = as.SingleCellExperiment(mat_seurat_S1)
#sce_S2 = as.SingleCellExperiment(mat_seurat_S2)
mat_S1 = scm_import_sce_to_mat(sce_S1)
#mat_S2 = scm_import_sce_to_mat(sce_S2)

mat_S1@cells = paste0('HN_07_S1_', mat_S1@cells)
rownames(mat_S1@cell_metadata) = mat_S1@cells
colnames(mat_S1@mat) = mat_S1@cells

#mat_S2@cells = paste0('HN_05_S2_', mat_S2@cells)
#rownames(mat_S2@cell_metadata) = mat_S2@cells
#colnames(mat_S2@mat) = mat_S2@cells

mat_S1@cell_metadata[, "orig.ident"] = "HN_07_S1"
#mat_S2@cell_metadata[, "orig.ident"] = "HN_05_S2"
mat_S1@cell_metadata$seq_run = "GCF-6027"
#mat_S2@cell_metadata$seq_run = "GCF-6009"

scdb_add_mat(mat = mat_S1, id = mat_id_S1)
#scdb_add_mat(mat = mat_S2, id = mat_id_S2)

mat_names = c(mat_id_S1)
for (i in mat_names){
  mat = scdb_mat(i)
  print(dim(mat@mat))
  
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO1-pre" |
                      mat@cell_metadata$hash.ID == "HTO2-post", "patient"] = "Pat26"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO4-pre"   |
                      mat@cell_metadata$hash.ID == "HTO6-post", "patient"] = "Pat33"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO7-pre" |
                      mat@cell_metadata$hash.ID == "HTO8-post", "patient"] = "Pat04"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO9-pre"   |
                      mat@cell_metadata$hash.ID == "HTO10-post", "patient"] = "Pat21"
  
  mat@cell_metadata$HT = substr(mat@cell_metadata$hash.ID, 1, 5)
  mat@cell_metadata$new_hash_ID = mat@cell_metadata$hash.ID
  
  mat@cell_metadata$condition = substr(mat@cell_metadata$new_hash_ID, 6,  50)
  mat@cell_metadata[mat@cell_metadata$condition == "ive", "condition"] = "neg"
  mat@cell_metadata[mat@cell_metadata$condition == "et", "condition"] = "doub"
  mat@cell_metadata[mat@cell_metadata$condition == "-post", "condition"] = "post"
  
  scdb_add_mat(id = i,  mat = mat)
}

mat_S1 = scdb_mat(mat_id_S1)  



############### SCHN_08 ############

# load data
data = load("HN_full_dataset_scdb/ori_objects_online_only/GEX_VDJ_6031_S1.Rda")
mat_seurat_S1 = counts
data = load("HN_full_dataset_scdb/ori_objects_online_only/GEX_VDJ_6031_S2.Rda")
mat_seurat_S2 = counts

# export seurat object to a metacell object
sce_S1 = as.SingleCellExperiment(mat_seurat_S1)
sce_S2 = as.SingleCellExperiment(mat_seurat_S2)
mat_S1 = scm_import_sce_to_mat(sce_S1)
mat_S2 = scm_import_sce_to_mat(sce_S2)

mat_S1@cells = paste0('HN_08_S1_', mat_S1@cells)
rownames(mat_S1@cell_metadata) = mat_S1@cells
colnames(mat_S1@mat) = mat_S1@cells

mat_S2@cells = paste0('HN_08_S2_', mat_S2@cells)
rownames(mat_S2@cell_metadata) = mat_S2@cells
colnames(mat_S2@mat) = mat_S2@cells

mat_S1@cell_metadata[, "orig.ident"] = "HN_08_S1"
mat_S2@cell_metadata[, "orig.ident"] = "HN_08_S2"
mat_S1@cell_metadata$seq_run = "GCF-6031"
mat_S2@cell_metadata$seq_run = "GCF-6031"

scdb_add_mat(mat = mat_S1, id = mat_id_S1)
scdb_add_mat(mat = mat_S2, id = mat_id_S2)

mat_names = c(mat_id_S1, mat_id_S2)
for (i in mat_names){
  mat = scdb_mat(i)
  print(dim(mat@mat))
  
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO9-pre" |
                      mat@cell_metadata$hash.ID == "HTO10-post", "patient"] = "Pat39"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO1-pre"   |
                      mat@cell_metadata$hash.ID == "HTO2-post", "patient"] = "Pat27"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO7-pre", "patient"] = "Pat17"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO8-post", "patient"] = "Pat22"
  mat@cell_metadata[mat@cell_metadata$hash.ID == "HTO4-pre"   |
                      mat@cell_metadata$hash.ID == "HTO6-post", "patient"] = "Pat36"
  
  
  
  mat@cell_metadata$HT = substr(mat@cell_metadata$hash.ID, 1, 5)
  mat@cell_metadata$new_hash_ID = mat@cell_metadata$hash.ID
  
  mat@cell_metadata$condition = substr(mat@cell_metadata$new_hash_ID, 6,  50)
  mat@cell_metadata[mat@cell_metadata$condition == "ive", "condition"] = "neg"
  mat@cell_metadata[mat@cell_metadata$condition == "et", "condition"] = "doub"
  mat@cell_metadata[mat@cell_metadata$condition == "-post", "condition"] = "post"
  
  scdb_add_mat(id = i,  mat = mat)
}

mat_S1 = scdb_mat(mat_id_S1)  
mat_S2 = scdb_mat(mat_id_S2)  


