tbk_reload = function() {
  library("metacell")
  library("dplyr")
  #library("flowCore")
  library("Matrix")
  library("tidyr")
  library("ggplot2")
  library("data.table")
  library("Seurat")
  library("pheatmap")
  
  if(!dir.exists("HN_04_scdb")) dir.create("HN_04_scdb/") 
  scdb_init("HN_04_scdb/", force_reinit=T)
  
  if(!dir.exists("HN_02_figs")) dir.create("HN_02_figs/") 
  scfigs_init("HN_02_figs/")
  
  exp_nm <<- "HN_04"
  mat_id <<- exp_nm
  mat_id_S1 <<- paste0(exp_nm, "_S1")
  mat_id_S2 <<- paste0(exp_nm, "_S2")
  p_mat_id <<- paste0("p_", exp_nm)
  clean_mat_id <<- paste0(exp_nm, "_clean")
  gstat_id <<- clean_mat_id
  feats_gset_id <<- paste0(clean_mat_id, "_feats")
  graph_id <<- paste0(exp_nm, "_graph")
  markers_gset_id <<- paste0(clean_mat_id, "_markers")
  mc_id<<- paste0(exp_nm, "_mc")
  filt_mc_id<<- paste0(mc_id,"_f")
  mc2d_id <<- paste0(mc_id, "_2dproj")
  lateral_gset_mel_id <<- "mel_lateral"
}

## plot cumulative distribution of total cell UMI counts
plot_cum_umi_distr = function(mat_id) {
  
  mat = scdb_mat(mat_id)
  uc <<- Matrix::colSums(mat@mat)
  metacell:::.plot_start(scfigs_fn(mat_id, "cumulative_UMI_distr"), w = 450, h = 450)
  
  plot(ecdf(log2(uc)), xlim= c(0, 15), xlab = "total cell UMIs (log2)", ylab = "fraction of data", main = "ECDF of UMI counts")
  dev.off()
}