########## #########
### METACELL GENERATION PIPELINE
########## #########


###### load data and source code ######

# load from laptop
setwd("/DATA/project")

# load from MC0123
#setwd("/Users/a.vd.leun/surfdrive/Metacell/IMCISION/HN_full_dataset/")

source("Functions_hn_full.R")
tbk_reload()

# load previously generate metacell objects #
mat_HN_02_S1 = scdb_mat("HN_02_S1")
mat_HN_02_S2 = scdb_mat("HN_02_S2")
mat_HN_04_S1 = scdb_mat("HN_04_S1")
mat_HN_05_S1 = scdb_mat("HN_05_S1")
mat_HN_05_S2 = scdb_mat("HN_05_S2")
mat_HN_06_S1 = scdb_mat("HN_06_S1")
mat_HN_07_S1 = scdb_mat("HN_07_S1")
mat_HN_08_S1 = scdb_mat("HN_08_S1")
mat_HN_08_S2 = scdb_mat("HN_08_S2")

full_mat = scm_merge_mats(scmat1 = mat_HN_02_S1, scmat2 = mat_HN_02_S2)
full_mat = scm_merge_mats(scmat1 = full_mat, scmat2 = mat_HN_04_S1)
full_mat = scm_merge_mats(scmat1 = full_mat, scmat2 = mat_HN_05_S1)
full_mat = scm_merge_mats(scmat1 = full_mat, scmat2 = mat_HN_05_S2)
full_mat = scm_merge_mats(scmat1 = full_mat, scmat2 = mat_HN_06_S1)
full_mat = scm_merge_mats(scmat1 = full_mat, scmat2 = mat_HN_07_S1)
full_mat = scm_merge_mats(scmat1 = full_mat, scmat2 = mat_HN_08_S1)
full_mat = scm_merge_mats(scmat1 = full_mat, scmat2 = mat_HN_08_S2)

scdb_add_mat(mat_id, full_mat)

mat_list = list(rownames(mat_HN_02_S1@mat),
                rownames(mat_HN_02_S2@mat),
                rownames(mat_HN_04_S1@mat),
                rownames(mat_HN_05_S1@mat), 
                rownames(mat_HN_05_S2@mat), 
                rownames(mat_HN_06_S1@mat), 
                rownames(mat_HN_07_S1@mat),
                rownames(mat_HN_08_S1@mat),
                rownames(mat_HN_08_S1@mat))

non_overlapping_genes = lapply(1:length(mat_list), 
                               function(n) setdiff(mat_list[[n]], unlist(mat_list[-n])))
exclude_genes = unlist(non_overlapping_genes)
write.csv(x = exclude_genes, "HN_full_dataset_scdb/exclude_genes.csv")

# collapse clonotype IDs - create df with unique cdr3s AA SEQS!!! and add new clonotypeid NO NTs IN THIS DATASET
nt <- data.frame(cdr3s_nt_check=unique(full_mat@cell_metadata$chain_cdr3), new_clonotype_id=paste("clonotype",1:length(unique(full_mat@cell_metadata$chain_cdr3))-1, sep=""), new_clononum=1:length(unique(full_mat@cell_metadata$chain_cdr3))-1 )

# m atch the original cdr3s within the unique set -> create dataframe of results
newctdf <- nt[(match(full_mat@cell_metadata$chain_cdr3, nt$cdr3s_nt_check)),]
# add rownames from original metadata
rownames(newctdf) <- rownames(full_mat@cell_metadata)

# add dataframe to metadata
full_mat@cell_metadata[rownames(newctdf), "unique_clone_ID"] =  newctdf[rownames(newctdf), "new_clononum"]

# add TCR frequencies to metadata
tcr_counts = plyr::count(full_mat@cell_metadata$unique_clone_ID)
for (i in tcr_counts$x){
  full_mat@cell_metadata[full_mat@cell_metadata$unique_clone_ID == i, "TCR_count"] = tcr_counts[tcr_counts$x == i ,"freq"]
}
full_mat@cell_metadata[full_mat@cell_metadata$unique_clone_ID == 0, c("unique_clone_ID", "TCR_count")] =  NA
full_mat@cell_metadata$pat_clone_ID = paste0(full_mat@cell_metadata$patient, "_", full_mat@cell_metadata$unique_clone_ID)

scdb_add_mat(mat_id, full_mat)


# add TCRs to lateral_gset
lateral_gset_mel = scdb_gset(lateral_gset_id)

tcr_genes = grep("^TRB[V|D|J|]|^TRA[V|D|J|]|^TRD[V|D|J|]|^TRG[V|D|J|]", full_mat@genes, v=T, perl=T)

lateral_gset_tcrs_id = paste0(lateral_gset_id, "_tcrs")
scdb_add_gset(lateral_gset_tcrs_id, lateral_gset_mel)
mcell_gset_add_gene(lateral_gset_tcrs_id, tcr_genes, 4)

##### TAKE TO THE SERVER ######

