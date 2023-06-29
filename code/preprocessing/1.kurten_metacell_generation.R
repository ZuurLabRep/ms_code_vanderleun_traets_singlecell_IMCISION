################
#### Kurten ####
################

# preprocessing of data from Kürten CHL, Kulkarni A, Cillo AR, Santos PM, Roble AK, Onkar S, Reeder C, Lang S, Chen X, Duvvuri U, Kim S, Liu A, Tabib T, Lafyatis R, Feng J, Gao SJ, Bruno TC, Vignali DAA, Lu X, Bao R, Vujanovic L, Ferris RL. Investigating immune and non-immune cell interactions in head and neck tumors by single-cell RNA sequencing. Nat Commun. 2021 Dec 17;12(1):7338. doi: 10.1038/s41467-021-27619-4. PMID: 34921143; PMCID: PMC8683505.
# using feature genes generated in the single cell IMCISION processing with Metacell
setwd("/DATA/project")
source("./code/Figures_general_functions.R")
source("./code/Functions_reprocessing_metacell.R")
tbk_reload()

if(!dir.exists("figures/Figures_Kurten")) dir.create("Figures_Kurten/") 
scfigs_init("figures/Figures_Kurten/")
fig_path = .scfigs_base

if(!dir.exists("data/Data_Kurten")) dir.create("data/Data_Kurten/") 
scdb_init("data/Data_Kurten/", force_reinit=T)

exp_nm <<- "Kurten_full_dataset"
mat_id <<- exp_nm
p_mat_id <<- paste0("p_", exp_nm)
clean_mat_id <<- paste0(exp_nm, "_clean")
gstat_id <<- clean_mat_id
feats_gset_id <<- paste0(clean_mat_id, "_feats")
filt_feats_gset_id <<- paste0(feats_gset_id, "_filtered")
graph_id <<- paste0(exp_nm, "_graph")
markers_gset_id <<- paste0(clean_mat_id, "_markers")
mc_id<<- paste0(exp_nm, "_mc")
filt_mc_id<<- paste0(mc_id,"_f")
mc2d_id <<- paste0(mc_id, "_2dproj")
lateral_gset_mel_id <<- "mel_lateral"
lat_mod_id <<- "lateral_modules"
lat_gset_hn_id <<- "HN_lateral"

## create seurat object from GEO Ferris data, only CD45p
library(Seurat)
library(scCustomize)

setwd("~/Documents/Zuur_lab/Single_cell_data/IMCISION_final_vanderleun_rev2/Final_code Github/data/Data_Kurten")
sample_IDs<-system("ls", intern=T)
sample_IDs<-sample_IDs[grep("*CD45p_*", sample_IDs)]
sample_names_fix <- as.data.frame(do.call(rbind,strsplit(sample_IDs,"_"))[,1:3])
sample_names_fix$names <- apply( sample_names_fix[ , 1:3 ] , 1 , paste , collapse = "_" )

sample_names_fix = unique(sample_names_fix$names)

sample_name_x = sample_names_fix[1]
data <- Read10X_GEO(data.dir = "./", sample.names =paste0(sample_name_x, "_"))
seurat_object_merged <- CreateSeuratObject(counts = data[[paste0(sample_name_x, "_")]]$`Gene Expression`,project=sample_name_x)

for (sample_name_x in sample_names_fix[2:18]){
  data <- Read10X_GEO(data.dir = "./", sample.names =paste0(sample_name_x, "_"))
  seurat_object_tmp <- CreateSeuratObject(counts = data[[paste0(sample_name_x, "_")]]$`Gene Expression`,project=sample_name_x)
  seurat_object_merged = merge(seurat_object_merged,seurat_object_tmp,project="All_sets_CD45_Ferrris")
}

seurat_object_merged = subset(x = seurat_object_merged, subset = orig.ident %in% sample_names_fix[c(1:11,15)])
seurat_object_merged = subset(x = seurat_object_merged, subset = orig.ident != "GSM5017027_HN03_CD45p") # outlier in data

# from seurat
setwd("~/Documents/Zuur_lab/Single_cell_data/IMCISION_final_vanderleun_rev2/Final_code Github/")
sce_ferris = as.SingleCellExperiment(seurat_object_merged)
mat_ferris = scm_import_sce_to_mat(sce_ferris)

scdb_add_mat(mat = mat_ferris, id = mat_id)

mcell_plot_umis_per_cell(mat_id)
plot_cum_umi_distr(mat_id)

## load dataset check
mat = scdb_mat(mat_id)
print(dim(mat@mat))

library(tgstat) 
options(tgs_max.processes=5)

## 2. process metacell
mat = scdb_mat(mat_id)
print(dim(mat@mat))

mcell_plot_umis_per_cell(mat_id,min_umis_cutoff=800)  #800
plot_cum_umi_distr(mat_id)

mel_build_blacklist_gsets_by_gene_nms_ann(mat_id)  # check mitocarte file in right directory - hardcoded
mito_gset_id = "mito_gset_updated"
plot_mito_genes(mat_id, mito_gset_id) 

# also removes cells of bad quality < 800 reads
clean_up_data_set(mat_id = mat_id, 
                  max_mito_f = 0.4, ## check mito
                  min_umis_post_gene_ignore = 800) ### CHECK MIN UMIS #### 
test = scdb_mat(clean_mat_id)
print(dim(test@mat))

# derive gene stats for all genes
mcell_add_gene_stat(gstat_id, clean_mat_id, force = T)
gstat = scdb_gstat(gstat_id)
dim(gstat)

# define most differential genes and plot their statistics
generate_feats_gset(gstat_id, feats_gset_id,T_vm = 0.95) #### CHECK TVM ####   0.95
length(names(feats_gset@gene_set))
# replace feats
feats_gset_id = "HN_full_dataset_clean_feats_filtered"
feats_gset = scdb_gset(feats_gset_id) # replace with head and neck feats!!

# cross check feature genes with genes that are not overlapping between datasets
feats_gset@gene_set = feats_gset@gene_set[names(feats_gset@gene_set) %in% rownames(test@mat)]
# not properly in data
feats_gset@gene_set = feats_gset@gene_set[names(feats_gset@gene_set) %in% gstat$name]
scdb_add_gset(feats_gset_id, feats_gset)

# check for additional lateral genes
check_lateral_genes(gset_id = feats_gset_id, mat_id = clean_mat_id, K = 50,
                    filt_clusts = c(),
                    new_lat_gset_id = lat_mod_id)
# save lateral genes
lat_mod_gset = scdb_gset(lat_mod_id)

feats_gset = scdb_gset(feats_gset_id)
#remove_genes = feats_gset@gene_set[names(feats_gset@gene_set) %in% exclude_genes]

# add non-overlapping genes (remove_genes) to lateral geneset
load("./data/gset.mel_lateral.Rda")
object@gene_set = object@gene_set[names(object@gene_set) %in% rownames(test@mat)]
add_lat = c(names(object@gene_set),
            rownames(mat@mat[grep("IFI", rownames(mat@mat)),]),
            rownames(mat@mat[grep("HIST", rownames(mat@mat)),]),
            rownames(mat@mat[grep("HSP", rownames(mat@mat)),]),
            c("RPLP1", "RPLP0", "EEF1A1", "TPT1"))

lat_gset = gset_add_genes(lat_mod_gset, add_lat, 0)
scdb_add_gset(id = lat_gset_hn_id, lat_gset)

# filter feature genes for lateral genes
filter_feats_gset(feats_gset_id, lat_gset_hn_id)
length(names(filt_feats_gset@gene_set))

metacell:::.plot_start(scfigs_fn(mat_id, "mean_top_feature_genes"), 1200,1200)
genes = gstat[names(filt_feats_gset@gene_set),]
plot(tail(genes[order(genes[,"ds_mean"]),"ds_mean"], 50), pch = 19, cex = 0.5, main = "mean_top_feats")
text(x = tail(genes[order(genes[,"ds_mean"]),"ds_mean"], 50), labels = rownames(tail(genes[order(genes[,"ds_mean"]),], 50)), cex = 1, pos = 2)
dev.off()
write.csv(genes[,c("name", "ds_mean")], file = "./data/Data_Ferris/filt_feats_expression.csv")

high_expr = rownames(genes[genes$ds_mean > 3,])
lat_gset = gset_add_genes(lat_gset, high_expr, 100)
scdb_add_gset(id = lat_gset_hn_id, lat_gset)

# filter feature genes for lateral genes
filter_feats_gset(feats_gset_id, lat_gset_hn_id)
length(names(filt_feats_gset@gene_set))

library(tgstat) 
options(tgs_max.processes=10)

## construct
mcell_add_cgraph_from_mat_bknn(clean_mat_id, filt_feats_gset_id, graph_id, K=100, dsamp=T) # K=100

# compute metacell using resampling iterations of graph cover k-means-like approach
mcell_coclust_from_graph_resamp(coc_id <<- paste0(exp_nm, "_coc500"), graph_id, min_mc_size=20, p_resamp=0.75, n_resamp=500) 

# build a metacell cover from co-clust data through filtering un-balanced edges and running graph cover
mcell_mc_from_coclust_balanced(mc_id = mc_id, coc_id=paste0(exp_nm, "_coc500"), mat_id = clean_mat_id, K=35, min_mc_size=35, alpha=2)  #K=30, min_mc_size=30,

mc = scdb_mc(mc_id)
length(names(mc@mc))
length(names(mc@annots))

# remove outliers
mcell_plot_outlier_heatmap(mc_id, clean_mat_id, T_lfc=3000)
mcell_mc_split_filt(filt_mc_id,  mc_id,  clean_mat_id, T_lfc=3000, plot_mats=F)

mc <<- scdb_mc(filt_mc_id)

# use filt mc
mc_hc = mcell_mc_hclust_confu(mc_id = filt_mc_id, graph_id = graph_id)
mc_sup = mcell_mc_hierarchy(mc_id = filt_mc_id, mc_hc = mc_hc, T_gap = 0.14)

mcell_mc_plot_hierarchy(mc_id = filt_mc_id, graph_id = graph_id,
                        mc_order = mc_hc$order,
                        sup_mc = mc_sup,
                        width = 1200, height = 2400,
                        min_nmc=2, show_mc_ids = T)

lfp = log2(mc@mc_fp)
generate_lfp_table(filt_mc_id, gstat_id, clean_mat_id)

mc_colorize_default(filt_mc_id)
mc <<- scdb_mc(filt_mc_id)

# check expression of genes of interest for each metacell
plt = function(gene1, gene2, lfp, colors)
{
  plot(lfp[gene1, ], lfp[gene2, ], pch=21, cex=3, bg=colors, xlab=gene1, ylab=gene2)
  text(lfp[gene1, ], lfp[gene2, ], colnames(lfp))
}

# plot genes (rerun mc filter)
gene1 <- "GZMB"
gene2 <- "ZNF683"
pdf(file=paste0("./Arnon_figs/",gene1,"_",gene2,".pdf"), width=5, height=5) #
plt(gene1 = gene1, gene2 = gene2, lfp = lfp, colors = mc@colors)
dev.off()


supid <- c(4,11,16,20,23,28,35,41,46,55,63,72,93,111,78,77,109,82)
color <- c("#E2ADF2","#73D9D5","#780585","brown","#DDAF08","#D97373","#0042b0","#59b56f",
           "#F56105","#C20088","#4C005C","red","#A1A1A1","grey","#e6e600","#8b0000","#0089e0","#ffa500") #,"#0E30DD","#D97373","#1A5606","#73D993",
# "#F56105","#DDB1B1","#DDAF08","#73D9D5","#780585"            DDAF08  0E30DD
name <- c("cell1","cell2","cell3","cell4","cell5","cell6","Treg_GIMAP",
          "Treg_TNF","cell9","Dysf_ZNF","Dysf_GZMB","cell11","cell12","cell13","cDC_CLEC9A","cDC_LAMP3","pDC","cDC_CD1C")
supmc_tab <- as.data.frame(cbind(supid,color,name))
write.table(supmc_tab,"./Data_Ferris/CD45_tumor_supmc.txt",sep="\t")

# plot confusion matrix
mc_colorize_sup_hierarchy(mc_id = filt_mc_id,
                          supmc = mc_sup,
                          supmc_key = "./Data_Ferris/CD45_tumor_supmc.txt")
#gene_key = "./S1_pilot_scdb/S1_marks.txt")

mcell_mc_plot_hierarchy(mc_id = filt_mc_id,
                        graph_id = graph_id,
                        mc_order = mc_hc$order,
                        sup_mc = mc_sup,
                        width=1800, height=2400, min_nmc=2, show_mc_ids= T)

# 2d plot
mcell_mc2d_force_knn(mc2d_id, filt_mc_id, graph_id, ignore_mismatch = T)
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id=mc2d_id, legend_pos = "topleft",sc_cex=2)

# reload!
mc <<- scdb_mc(filt_mc_id)

# replace with new mat
col2group <- get_mc_col2group(mc)
group2col <- get_mc_group2col(mc)
mat = scdb_mat(clean_mat_id)
mat@cell_metadata[names(mc@mc), 'mc_group'] <- col2group[mc@colors[mc@mc]]
mat@cell_metadata[names(mc@mc), 'mc'] <- as.data.frame(mc@mc)

scdb_add_mat(id = clean_mat_id, mat = mat)

