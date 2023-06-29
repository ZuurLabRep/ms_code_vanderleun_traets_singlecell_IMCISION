### metacell formation on large datasets ###

# load code and input form server
setwd("/DATA/project")

source("./Functions_hn_coley.R")
tbk_reload()

# load merged dataset
mat = scdb_mat(mat_id)
print(dim(mat@mat))

# plot UMI stats of unfiltered data
mcell_plot_umis_per_cell(mat_id)  

# plot ECDF umi counts
plot_cum_umi_distr(mat_id)

# clean up data by blacklisting gene lists and excluding cells with low UMI count or high mito gene expression
mel_build_blacklist_gsets_by_gene_nms_ann(mat_id)  # check mitocarte file in right directory - hardcoded
mito_gset_id = "mito_gset_updated"
plot_mito_genes(mat_id, mito_gset_id) 

clean_up_data_set(mat_id = mat_id, 
                  max_mito_f = 0.5, 
                  min_umis_post_gene_ignore = 800) ### CHECK MIN UMIS ####
print(dim(mat@mat))

# derive gene stats for all genes
mcell_add_gene_stat(gstat_id, clean_mat_id, force = T)
gstat = scdb_gstat(gstat_id)
dim(gstat)

# define most differential genes and plot their statistics
generate_feats_gset(gstat_id, feats_gset_id) #### CHECK TVM ####
length(names(feats_gset@gene_set))

# cross check feature genes with genes that are not overlapping between datasets
exclude_genes = read.csv(file = "HN_full_dataset_scdb/exclude_genes.csv")
exclude_genes = as.vector(exclude_genes$x)

feats_gset = scdb_gset(feats_gset_id)
remove_genes = feats_gset@gene_set[names(feats_gset@gene_set) %in% exclude_genes]
#write.csv(x = remove_genes, file = "HN_full_dataset_scdb/excluded_genes_from_feats.csv")

# check for additional lateral genes
check_lateral_genes(gset_id = feats_gset_id, mat_id = clean_mat_id, K = 70, 
                    filt_clusts = c(2,13,16,23,26,29,35),
                    new_lat_gset_id = lat_mod_id)

#XX/01/2021 c(2,15,19,56,50,35,38,39,45,26,27,60)
lat_mod_gset = scdb_gset(lat_mod_id)

feats_gset = scdb_gset(feats_gset_id)

# check expression of modules between exps
for(i in unique(feats_gset@gene_set)){
 gset = names(feats_gset@gene_set[feats_gset@gene_set == i])

  cc_mat = mat@mat[intersect(gset, rownames(mat@mat)),]

  cc_umis = colSums(cc_mat)
  cc_fraction = cc_umis/colSums(mat@mat)

  md = mat@cell_metadata[names(cc_fraction), ]


  metacell:::.plot_start(scfigs_fn(paste0("HN_full_dataset_module_comp/", mat_id), sprintf("%s_by_%s", i, "experiment")), 600, 400)
  mar.default <- c(8,4,2,2) + 0.1
  par(mar = mar.default)

  boxplot(cc_fraction[rownames(md)] ~  md[, "orig.ident"], las = 2, ylab = "fraction of UMIs", ylim = c(0,0.1), main = sprintf("%s_%s_by_%s", i, "MEL57", "experiment"))
  stripchart(cc_fraction[rownames(md)] ~ md[, "orig.ident"], ylim = c(0,0.1), vertical = TRUE,  method = "jitter", add = TRUE, pch = 20, cex=0.1, col = 'blue')

  dev.off()

 }

# add non-overlapping genes (remove_genes) to lateral geneset
add_lat = c(names(remove_genes), 
            rownames(mat@mat[grep("IFI", rownames(mat@mat)),]),
            rownames(mat@mat[grep("HIST", rownames(mat@mat)),]),
            rownames(mat@mat[grep("HSP", rownames(mat@mat)),]),
            c("RPLP1", "RPLP0", "EEF1A1", "TPT1"))

lat_gset = gset_add_genes(lat_mod_gset, add_lat, 0)
scdb_add_gset(id = lat_gset_hn_id, lat_gset)

# filter feature genes for lateral genes
filter_feats_gset(feats_gset_id, lat_gset_hn_id)
length(names(filt_feats_gset@gene_set))

# plot mean of top variable genes in feats
metacell:::.plot_start(scfigs_fn(mat_id, "mean_top_feature_genes"), 1200,1200)
genes = gstat[names(filt_feats_gset@gene_set),]
plot(tail(genes[order(genes[,"ds_mean"]),"ds_mean"], 50), pch = 19, cex = 0.5, main = "mean_top_feats")
text(x = tail(genes[order(genes[,"ds_mean"]),"ds_mean"], 50), labels = rownames(tail(genes[order(genes[,"ds_mean"]),], 50)), cex = 1, pos = 2)
dev.off()
write.csv(genes[,c("name", "ds_mean")], file = "HN_full_dataset_figs/filt_feats_expression.csv")

high_expr = rownames(genes[genes$ds_mean > 4.8,])
lat_gset = gset_add_genes(lat_gset, high_expr, 100)
scdb_add_gset(id = lat_gset_hn_id, lat_gset)

# filter feature genes for lateral genes
filter_feats_gset(feats_gset_id, lat_gset_hn_id)
length(names(filt_feats_gset@gene_set))


# build a cell graph using balanced knn graph on given gene features
mcell_add_cgraph_from_mat_bknn(clean_mat_id, filt_feats_gset_id, graph_id, K=100, dsamp=T)

# compute metacell using resampling iterations of graph cover k-means-like approach
mcell_coclust_from_graph_resamp(coc_id <<- paste0(exp_nm, "_coc500"), graph_id, min_mc_size=20, p_resamp=0.75, n_resamp=500)

# build a metacell cover from co-clust data through filtering un-balanced edges and running graph cover
mcell_mc_from_coclust_balanced(mc_id = mc_id, coc_id=paste0(exp_nm, "_coc500"), mat_id = clean_mat_id, K=30, min_mc_size=30, alpha=2)

mc = scdb_mc(mc_id)
length(names(mc@mc))
length(names(mc@annots))

# split metacells and exclude outlier cells (see parameters in function)
#mcell_plot_outlier_heatmap(mc_id, clean_mat_id, T_lfc=100)

mcell_mc_split_filt(filt_mc_id, 
                    mc_id, 
                    clean_mat_id,
                    T_lfc=3000, plot_mats=F)

mc <<- scdb_mc(filt_mc_id)


###### generate 2D and confusion maps and plot gmod expression per metacell
mc_colorize_default(filt_mc_id)
generate_metacell_plots(mc_id = filt_mc_id,
                        graph_id = graph_id,
                        mc2d_id = mc2d_id,
                        supmc_file = NULL,
                        marks_file = NULL)

generate_metacell_plots(mc_id = filt_mc_id,
                        graph_id = graph_id,
                        mc2d_id = mc2d_id,
                        supmc_file = "config/sup_mc_colors.txt",
                        marks_file = "config/colorize.txt")

col2group = get_mc_col2group(mc)
group2col = get_mc_group2col(mc)

# generate marker geneset & plot marker genes heatmap
generate_marker_gset(gset_id = markers_gset_id, mc_id = filt_mc_id, blacklist_gset_id = lat_gset_hn_id, filt_gset_id = lat_gset_hn_id, mat_id = clean_mat_id)
marker_gset = scdb_gset(markers_gset_id)

# look at marker expression per gene module (filt feats genes)
filt_feats_gset = scdb_gset(filt_feats_gset_id)
View(filt_feats_gset)
for (i in unique(unname(filt_feats_gset@gene_set))) {
  mcell_gset_remove_clusts(gset_id = filt_feats_gset_id, filt_clusts = i, new_id = paste0(filt_feats_gset_id, "_", i), reverse = T)
  mcell_mc_plot_marks(mc_id = filt_mc_id, gset_id = paste0(filt_feats_gset_id, "_", i), mat_id = clean_mat_id,  plot_cells = F, 
                      fig_fn = paste0("HN_full_dataset_figs/HN_full_dataset_heatmark_mods/HN_full_dataset_heatmark_mod_", i, ".png"))
}

# generate lfp table
lfp = log2(mc@mc_fp)
generate_lfp_table(filt_mc_id, gstat_id, clean_mat_id)


# add metacells and metacell_groups to metadata
mat = scdb_mat(clean_mat_id)
mat@cell_metadata[names(mc@mc), 'mc_group'] = col2group[mc@colors[mc@mc]]
mat@cell_metadata[names(mc@mc), 'mc'] = as.data.frame(mc@mc)
mat@cell_metadata$pat_condition = paste0(mat@cell_metadata$patient, "_", mat@cell_metadata$condition)
mat@cell_metadata$pat_mc_group_condition = paste0(mat@cell_metadata$mc_group, "_", mat@cell_metadata$pat_condition)

mat@cell_metadata[is.na(mat@cell_metadata$frequency1), "frequency1"] = 0
mat@cell_metadata[mat@cell_metadata$frequency1 < 1 , "clonal_exp"] = "no"
mat@cell_metadata[mat@cell_metadata$frequency1 == 1 , "clonal_exp"] = "low"
mat@cell_metadata[mat@cell_metadata$frequency1 > 1 & mat@cell_metadata$frequency1 < 5 , "clonal_exp"] = "int"
mat@cell_metadata[mat@cell_metadata$frequency1 >= 5 , "clonal_exp"] = "hi"
mat@cell_metadata$factor = paste0(mat@cell_metadata$patient, "_", mat@cell_metadata$clonal_exp)
mat@cell_metadata[mat@cell_metadata$frequency1 > 5, "expanded"] = paste0(mat@cell_metadata[mat@cell_metadata$frequency1 >5, "unique_clone_ID"], "_", mat@cell_metadata[mat@cell_metadata$frequency1 >5, "patient"])

mat@cell_metadata[!is.na(mat@cell_metadata$patient), 'response'] = "NR"
mat@cell_metadata[(mat@cell_metadata$patient == "Pat15" |
                     mat@cell_metadata$patient == "Pat17" |
                     mat@cell_metadata$patient == "Pat22" |
                     mat@cell_metadata$patient == "Pat29" |
                     mat@cell_metadata$patient == "Pat31" |
                     mat@cell_metadata$patient == "Pat12" |
                     mat@cell_metadata$patient == "Pat04" |
                     mat@cell_metadata$patient == "Pat21" |
                     mat@cell_metadata$patient == "Pat39") & 
                    !is.na(mat@cell_metadata$patient) , 'response'] = "RE"

mat@cell_metadata$exp_pat_condition = paste0(mat@cell_metadata$orig.ident, "_", mat@cell_metadata$pat_condition)
mat@cell_metadata$condition_response = paste0(mat@cell_metadata$condition, "_", mat@cell_metadata$response)

scdb_add_mat(id = clean_mat_id, mat = mat)



#plot UMI distribution cleaned mat
mcell_plot_umis_per_cell(clean_mat_id)  

# plot UMI counts per mc_group
umis_counts(clean_mat_id, filt_mc_id, "mc", width =8000)

mcell_mc2d_plot_by_factor(mc2d_id = mc2d_id, mat_id = clean_mat_id, meta_field = "pat_condition", single_plot = F, ncols = 4)
