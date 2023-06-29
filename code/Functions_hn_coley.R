# install/update metacell
#if (!require("BiocManager")) install.packages('BiocManager')
#BiocManager::install('metacell',  site_repository = 'tanaylab.github.io/repo', update = FALSE)


# load packages and object names
tbk_reload = function() {
  library("metacell")
  library("dplyr")
  #library("flowCore")
  library("Matrix")
  library("tidyr")
  library("ggplot2")
  library("data.table")
  #library("Seurat")
  library("pheatmap")
  
  #tgconfig::override_params("../config/TBK.yaml", "metacell")
  
  if(!dir.exists("HN_full_dataset_scdb")) dir.create("HN_full_dataset_scdb/") 
  scdb_init("HN_full_dataset_scdb/", force_reinit=T)
  
  if(!dir.exists("HN_full_dataset_figs")) dir.create("HN_full_dataset_figs/") 
  scfigs_init("HN_full_dataset_figs/")
  
  exp_nm <<- "HN_full_dataset"
  mat_id <<- exp_nm
  p_mat_id <<- paste0("p_", exp_nm)
  clean_mat_id <<- "HN_full_dataset_clean_v3"
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
}

# plot cumulative distribution of total cell UMI counts
plot_cum_umi_distr = function(mat_id) {
  
  mat = scdb_mat(mat_id)
  uc <<- Matrix::colSums(mat@mat)
  metacell:::.plot_start(scfigs_fn(mat_id, "cumulative_UMI_distr"), w = 450, h = 450)
  
  plot(ecdf(log2(uc)), xlim= c(0, 15), xlab = "total cell UMIs (log2)", ylab = "fraction of data", main = "ECDF of UMI counts")
  dev.off()
}

# generate blacklist gsets
mel_build_blacklist_gsets_by_gene_nms_ann = function(all_id) 
{
  full_m = scdb_mat(all_id)
  mito_gset_id = "mito_gset_updated"
  ig_gset_id = "ig_gset_updated"
  ncrna_gset_id = "ncrna_gset_updated"
  ncrp_gset_id = "ncrp_gset_updated"
  
  # mitochondrial gene set
  mt_cands = grep("^MT-|^MTRN|^MTAT|^MTND|^MRP", full_m@genes, v=T, perl=T)
  mito = data.table::fread("./config_files/MitoCarta2_human.txt", header=T, sep="\t", stringsAsFactors=F)
  mt_both = intersect(mt_cands, mito$Symbol)
  mt_cands = setdiff(mt_cands, mt_both)
  mitocarta = setdiff(mito$Symbol, mt_both)
  mt_genes = c(rep('regexp MT', length(mt_cands)), rep('MitoCarta2', length(mitocarta)), rep('MitoCarta2 and regexp MT', length(mt_both)))
  names(mt_genes) = c(mt_cands, mitocarta, mt_both)
  scdb_add_gset(mito_gset_id, gset_new_gset(mt_genes, 'mitochondrial genes'))
  
  # IG gene set
  ig_nms = grep("^IGK|^IGL|^IGJ|^IGH|^IGBP|^IGSF", full_m@genes, v=T, perl=T)
  ig_genes = rep("IG", length(ig_nms))
  names(ig_genes) = ig_nms
  scdb_add_gset(ig_gset_id, gset_new_gset(ig_genes, 'IG genes'))
  
  # ncRNA gene set
  ncrna_nms = c('MALAT1', 'XIST', 'NEAT1', 'hsa-mir-6723')
  ncrna_genes = rep("ncRNA", length(ncrna_nms))
  names(ncrna_genes) = ncrna_nms
  scdb_add_gset(ncrna_gset_id, gset_new_gset(ncrna_genes, "ncRNA genes"))
  
  # RP pseudo genes set
  ncrp_nms = grep("^RP.[0-9]+", full_m@genes, v=T, perl=T)
  ncrp_genes = rep("ncRP", length(ncrp_nms))
  names(ncrp_genes) = ncrp_nms
  scdb_add_gset(ncrp_gset_id, gset_new_gset(ncrp_genes, "RP##- genes"))
  
}

# plot mito gene fraction 
plot_mito_genes = function(mat_id, gset_id){
  mat = scdb_mat(mat_id)
  mito_genes = names(scdb_gset(gset_id)@gene_set)
  
  metacell:::.plot_start(scfigs_fn(mat_id, "fraction_mito_vs_total_UMIs"), 400,400)
  
  mito_genes_inmat = intersect(mito_genes, rownames(mat@mat))
  mito_f = Matrix::colSums(mat@mat[mito_genes_inmat,])/colSums(mat@mat)
  
  plot(log2(uc), mito_f)
  dev.off()
  
}

# filter dataset 
clean_up_data_set = function(mat_id, min_umis_post_gene_ignore = 500, max_mito_f = 0.6, force_new=T){
  mat = scdb_mat(mat_id)
  
  ig_genes <<- names(scdb_gset("ig_gset_updated")@gene_set)
  mito_genes <<- names(scdb_gset("mito_gset_updated")@gene_set)
  ncRNA_genes <<- names(scdb_gset("ncrna_gset_updated")@gene_set)
  ncRP_genes <<- names(scdb_gset("ncrp_gset_updated")@gene_set)
  bad_genes <<- unique(c(ig_genes, mito_genes, ncRNA_genes, ncRP_genes))
  
  uc = Matrix::colSums(mat@mat)
  
  mt_genes = intersect(mito_genes, mat@genes)
  mito_f = Matrix::colSums(mat@mat[mt_genes, ]) / uc
  
  mcell_mat_ignore_genes(clean_mat_id, mat_id, bad_genes, reverse=F)
  
  clean_mat = scdb_mat(clean_mat_id)
  
  mcell_mat_ignore_cells(clean_mat_id, clean_mat_id, union(clean_mat@ignore_cells, names(uc)[mito_f >= max_mito_f | Matrix::colSums(clean_mat@mat) <= min_umis_post_gene_ignore | mat@cell_metadata$HTO_classification.global != "Singlet"]))
  
  mat <<- scdb_mat(clean_mat_id)
}

# generate and filter feats genesets
generate_feats_gset = function(gstat_id, gset_id, T_vm = 0.08) {
  
  mcell_gset_filter_varmean(gstat_id, gset_id, T_vm, force_new=T)
  mcell_gset_filter_cov(gstat_id, gset_id, T_tot=100, T_top3=2)
  
  feats_gset <<- scdb_gset(gset_id)
  gstat <<- scdb_gstat(gstat_id)  
  
  mcell_plot_gstats(gstat_id, gset_id)
}


# look for additional lateral genes besides the mel_gset:
check_lateral_genes = function(gset_id, mat_id, K, filt_clusts = NULL, new_lat_gset_id){
  
  gset = scdb_gset(gset_id)
  
  lat_mat_id = paste0(mat_id, "_lateral")
  mcell_mat_ignore_genes(lat_mat_id, mat_id, names(gset@gene_set), reverse = T )
  mcell_gset_split_by_dsmat(gset_id = gset_id, mat_id = lat_mat_id, K = K, force = F)
  
  mcell_plot_gset_cor_mats(gset_id = gset_id, scmat_id = lat_mat_id, downsamp = T)
  mcell_gset_remove_clusts(gset_id = gset_id, filt_clusts = filt_clusts, new_id = new_lat_gset_id, reverse = T)
  
}


filter_feats_gset = function(gset_id, filt_id){
  feats_gset = scdb_gset(gset_id)
  filt_gset = scdb_gset(filt_id)
  
  feats_gset = gset_new_restrict_gset(feats_gset, filt_gset, inverse = T, desc = "cgraph gene markers w/o lateral genes")
  scdb_add_gset(paste0(feats_gset_id, "_filtered"), gset = feats_gset)
  
  filt_feats_gset_id <<- paste0(feats_gset_id, "_filtered")
  filt_feats_gset <<- scdb_gset(filt_feats_gset_id)
}


# split outliers cells

split_outliers = function (mc_id, mat_id, filt_mc_id) {
  mcell_plot_outlier_heatmap(mc_id, clean_mat_id, T_lfc=100)
  
  mcell_mc_split_filt(filt_mc_id, 
                      mc_id, 
                      clean_mat_id,
                      T_lfc=3000, plot_mats=T)
  
  mc <<- scdb_mc(filt_mc_id) 
  
}


# generate confusion map and 2D projection
colorize_by_confusion_mat = function(mc_id, graph_id, supmc_file=NULL, marks_file=NULL, res=NULL, show_mc_ids=F) 
{
  # Cluster metacells by confusion matrix
  if (is.null(res)) {	
    mc_hc = mcell_mc_hclust_confu(mc_id=mc_id,
                                  graph_id=graph_id)
  }
  else {
    mc_hc = res$mc_hc
  }
  
  # Annotate metacell clusters (globally and locally enriched genes)
  if (is.null(res)) {
    mc_sup <<- mcell_mc_hierarchy(mc_id=mc_id,
                                  mc_hc=mc_hc, T_gap=0.04)
  }
  else {
    mc_sup = res$mc_sup
  }
  
  # Colorize metacells based on the manually created supmc_file (and optionally by the marks_file)
  if (!is.null(supmc_file)) {
    mc_colorize_sup_hierarchy(mc_id = mc_id,
                              supmc = mc_sup,
                              supmc_key = supmc_file,
                              gene_key = marks_file)
  }
  
  # generate metacell clusters heatmap
  mcell_mc_plot_hierarchy(mc_id=mc_id,
                          graph_id=graph_id,
                          mc_order=mc_hc$order,
                          sup_mc = mc_sup,
                          width=4800, heigh=8000, min_nmc=2, show_mc_ids=show_mc_ids)
  
  list(mc_hc=mc_hc, mc_sup=mc_sup)
}


generate_metacell_plots = function(mc_id, graph_id, supmc_file = NULL, marks_file = NULL, mc2d_id) {
  
  conf_res <<- colorize_by_confusion_mat(mc_id, graph_id)
  conf_res <<- colorize_by_confusion_mat(mc_id, graph_id, 
                                         supmc_file = supmc_file,
                                         marks_file = marks_file,
                                         res = conf_res)
  
  mc <<- scdb_mc(mc_id)
  
  mcell_mc2d_force_knn(mc2d_id, mc_id, graph_id, ignore_mismatch = T)
  tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
  tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
  tgconfig::set_param("mcell_mc2d_T_edge",0.03, "metacell")
  mcell_mc2d_plot(mc2d_id, plot_edges = F)
}


# get groups by colors
get_mc_col2group = function(mc, white_is_undet=T) {
  col2group = as.character(mc@color_key$group)
  names(col2group) = as.character(mc@color_key$color)
  if (white_is_undet) {
    if ('white' %in% mc@colors) {
      col2group = c(col2group, c('white'='Undet'))
    }
  }
  col2group
}

# get colors by groups
get_mc_group2col = function(mc, white_is_undet=T) {
  group2col = as.character(mc@color_key$color)
  names(group2col) = as.character(mc@color_key$group)
  if (white_is_undet) {
    if ('white' %in% mc@colors) {
      group2col = c(group2col, c('Undet'='white'))
    }
  }
  
  group2col
}


# generate list of marker genes
generate_marker_gset = function(gset_id, mc_id, blacklist_gset_id, filt_gset_id, mat_id){
  
  mcell_gset_from_mc_markers(gset_id = gset_id, mc_id = mc_id, blacklist_gset_id = blacklist_gset_id)
  mcell_gset_from_mc_markers(gset_id = paste0(gset_id, "_lateral"), mc_id = mc_id, filt_gset_id = filt_gset_id)
  
  mcell_mc_plot_marks(mc_id = mc_id, gset_id = gset_id, mat_id = mat_id, lateral_gset_id = paste0(gset_id, "_lateral"), plot_cells = T)
}

plt = function(nm1, nm2, lfp, cols, ofn=NULL, x=NULL, y=NULL, show_mc_ids=T, cex=3, cex.lab=1, add_grid=F, main="", xlim=NULL, ylim=NULL) {
  if (is.null(x)) {
    x = lfp[nm1, ]
  }
  if (is.null(y)) {
    y = lfp[nm2, ]
  }
  if (!is.null(ofn)) {
    .plot_start(ofn, 450, 450)
  }
  if (is.null(xlim)) {
    xlim = range(x)
  }
  if (is.null(ylim)) {
    ylim = range(y)
  }
  plot(x, y, pch=21, cex=cex, bg=cols, xlab=nm1, ylab=nm2, cex.lab=cex.lab, main=main, xlim=xlim, ylim=ylim)
  if (show_mc_ids) {
    text(x, y, colnames(lfp), cex=cex/4)
  }
  if (add_grid) {
    grid(col='black', lwd=0.5)
  }
  if (!is.null(ofn)) {
    dev.off()
  }
}


# export lfp table
generate_lfp_table = function(mc_id, gstat_id, mat_id){
  mc = scdb_mc(mc_id)
  lfp <<- log2(mc@mc_fp)
  mcell_mc_export_tab(mc_id, gstat_id, mat_id, metadata_fields = c('amp_batch_id'))
  
}


# plot UMI counts per cell for each mc_group
umis_counts = function(mat_id, mc_id, condition, width = 800, height = 400) {
  
  mat = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)
  
  metacell:::.plot_start(scfigs_fn(mc_id, sprintf("UMI_counts_per_%s", condition)), width, height)
  
  umi_counts = Matrix::colSums(mat@mat[,names(mc@mc)])
  
  plot_umis = ggplot(data = as.data.frame(umi_counts), 
                     mapping = aes(x = as.character(mat@cell_metadata[names(mc@mc), condition]),
                                   y = log2(umi_counts)))  + geom_violin()
  plot(plot_umis)
  
  dev.off()
}


celltype_proportions = function(mat_id, mc_id, condition, response, lead_state, order) {
  
  mat = scdb_mat(clean_mat_id)
  mc = scdb_mc(mc_id)
  md = mat@cell_metadata[names(mc@mc),]
  md = md[!is.na(md$mc_group) & md$condition == condition,]
  
  md = mat@cell_metadata[names(mc@mc),]
  
  n_cells_mc = as.data.frame(t(table(md[md$response == response, "patient"], 
                                     md[md$response == response, "cell_type"])))
  
  n_cells_mc_wide = spread(n_cells_mc, Var1, Freq)
  rownames(n_cells_mc_wide) = n_cells_mc_wide$Var2
  n_cells_mc_wide = n_cells_mc_wide[,-1]
  
  mc_prop <<- n_cells_mc_wide / Matrix::rowSums(n_cells_mc_wide)
  mc_prop = mc_prop[order(mc_prop[,lead_state]),]
  mc_prop = mc_prop[,order]

  tot_cells = rowSums(n_cells_mc_wide)
  
  mc_prop_r <- reshape2::melt(t(as.matrix(mc_prop)))
  
  mc_prop_r$patient <- factor(mc_prop_r$Var2,levels=unique(mc_prop_r[["Var2"]]))
  
  mc_prop_r = merge(mc_prop_r,n_cells_mc)
  
  width_fig= 3.2
  if (response == "RE"){
    width_fig = 4.3
  }
  
  pdf(file=paste0("./HN_manuscript_figures/Ordered_patients_barplots_",response,"_pre.pdf"), width=width_fig, height=2.8)
  p = ggplot(mc_prop_r,aes(x=patient,y=value,fill=Var1)) + geom_bar(position="fill", stat="identity") +
    ggtitle(paste0(response,"_pre CD8 ordered >")) + #scale_fill_manual(values=c(wesanderson::wes_palette("Moonrise3")[1:5])) +# +coord_flip() +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_text(aes(label=Freq), size = 2.3, position = position_stack(vjust = 0.5))
  print(p)
  dev.off()
}


DE_wilcoxon = function(on_t_df,pre_t_df){
  on_t_cells = rownames(on_t_df)
  pre_t_cells = rownames(pre_t_df)
  
  umis = as.matrix(mat_ds[,c(on_t_cells, pre_t_cells)])
  umis_n = sweep(umis,2,colSums(umis), "/") * 1000 # normalized
  
  # Fold-change
  logFC = as.vector(apply(umis, 1, function(x) log2((mean(x[on_t_cells]))/(mean(x[pre_t_cells])))))
  logFCx = as.vector(apply(umis, 1, function(x) log2(mean(x[on_t_cells]))))
  logFCy = as.vector(apply(umis, 1, function(x) log2(mean(x[pre_t_cells]))))
  
  # Wilcoxon-Mann-Whitney test
  wilc_adj = p.adjust(as.vector(apply(umis_n, 1, function(x) wilcox.test(x[on_t_cells], x[pre_t_cells])$p.value)), method = "fdr")
  wilc = apply(umis_n, 1, function(x) wilcox.test(x[on_t_cells], x[pre_t_cells])$p.value)
  
  pb <- data.frame(rownames(umis), logFC, logFCx, logFCy, wilc, wilc_adj)
  colnames(pb) <- c("gene", "fc", "x", "y", "pval", "pval_adj")
  pb$pval[pb$pval=="NaN"] <- NA
  pb <- na.omit(pb)
  pb <- pb[-which(abs(pb$fc) == "Inf"),]
  rownames(pb) <- pb$gene
  
  pb$pval[pb$pval == 0] <- 10e-280
  return(pb)
}


make_line_boxplot = function(cell_x, ylim_ax,prepost_frac_c,prepost_final_RE,prepost_final_NR){
  cell_x_counts <- prepost_frac_c[prepost_frac_c$mc_group == cell_x,]
  cell_x_counts$pat_cond <- paste0(cell_x_counts$patient,cell_x_counts$condition)
  prepost_final_RE$pat_cond <- paste0(prepost_final_RE$patient,prepost_final_RE$condition)
  prepost_final_RE_plot <- merge(prepost_final_RE,cell_x_counts,by= "pat_cond")
  pdf(file=paste0("./HN_manuscript_figures/",cell_x,"_post_pre_RE.pdf"), width=2.5, height=2.3) #
  p = ggplot(prepost_final_RE_plot[prepost_final_RE_plot$patient.x != "Pat27",], aes(x=condition.x, y=.data[[cell_x]], fill = condition.x)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(data=prepost_final_RE_plot[prepost_final_RE_plot$patient.x == "Pat27",],
               aes(x=condition.x, y=.data[[cell_x]],group = condition.x, size = count),
               pch = 21, alpha = 0.5, fill="white") +
    geom_line(data=prepost_final_RE_plot[prepost_final_RE_plot$patient.x == "Pat27",],aes(group = patient.x), color = "lightgrey", alpha = 0.5)+
    geom_point(aes(group = condition.x, size = count, fill=condition.x), pch = 21, alpha = 0.5,
               position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0))+ 
    geom_line(aes(group = patient.x), color = "lightgrey", alpha = 0.5)+
    ggtitle(paste0(cell_x, " ","RE")) +
    scale_fill_manual(values=c("#89CC9B","#116936")) +
    scale_size(name   = "# cells",
               breaks = c(10,20,40,80,160,320,640),
               labels = c(10,20,40,80,160,320,640),range=c(1,4),
               limits = c(0,500))+
    ylim(ylim_ax)+ ylab("fraction")+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(p)
  dev.off()
  print(prepost_final_RE_plot[prepost_final_RE_plot$condition.y == "pre",][[cell_x]])
  print(prepost_final_RE_plot[prepost_final_RE_plot$condition.y == "post",][[cell_x]])
  pval = wilcox.test(prepost_final_RE_plot[prepost_final_RE_plot$condition.y == "pre",][[cell_x]],
                     prepost_final_RE_plot[prepost_final_RE_plot$condition.y == "post",][[cell_x]])
  print("RE:")
  print(pval)
  
  prepost_final_NR$pat_cond <- paste0(prepost_final_NR$patient,prepost_final_NR$condition)
  prepost_final_NR_plot <- merge(prepost_final_NR,cell_x_counts,by= "pat_cond")
  pdf(file=paste0("./HN_manuscript_figures/",cell_x,"_post_pre_NR.pdf"), width=2.5, height=2.3) #
  p=ggplot(prepost_final_NR_plot, aes(x=condition.x, y=.data[[cell_x]], fill = condition.x)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(data=prepost_final_NR_plot[prepost_final_NR_plot$patient.x == "Pat27",],
               aes(x=condition.x, y=.data[[cell_x]],group = condition.x, size = count),
               pch = 21, alpha = 0.5, fill="white") +
    geom_line(data=prepost_final_NR_plot[prepost_final_NR_plot$patient.x == "Pat27",],aes(group = patient.x), color = "lightgrey", alpha = 0.5)+
    geom_point(aes(group = condition.x, size = count, fill=condition.x), pch = 21, alpha = 0.5,
               position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0))+ 
    geom_line(aes(group = patient.x), color = "lightgrey", alpha = 0.5)+
    ggtitle(paste0(cell_x, " ","NR")) +
    scale_fill_manual(values=c("#DC7985","#9B1F23")) +
    scale_size(name   = "# cells",
               breaks = c(10,20,40,80,160,320,640),
               labels = c(10,20,40,80,160,320,640),range=c(1,4),
               limits = c(0,800))+
    ylim(ylim_ax)+ ylab("fraction")+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(p)
  dev.off()
  pval = wilcox.test(prepost_final_NR_plot[prepost_final_NR_plot$condition.y == "pre",][[cell_x]],
                     prepost_final_NR_plot[prepost_final_NR_plot$condition.y == "post",][[cell_x]])
  print("NR:")
  print(pval)
}


diff_expr = function(mc, mat_ds, mcs1=NULL, mcs2=NULL, reg=5, min_max_umi=50, nms1=NULL, nms2=NULL, filter_outlier_genes=F, compare_top_mc_to_n_highest=3, max_top_to_n_highest_ratio=3, verbose=T, geo_mean=F, geo_mean_per_cell=F)
{
  if (is.null(nms1)) {
    nms1 = names(mc@mc)[mc@mc %in% mcs1]
  }
  if (is.null(nms2)) {
    nms2 = names(mc@mc)[mc@mc %in% mcs2]
  }
  nms1 = intersect(colnames(mat_ds), nms1)
  nms2 = intersect(colnames(mat_ds), nms2)
  if (verbose) {
    message(sprintf("comparing %d vs %d cells", length(nms1), length(nms2)))
  }
  
  if (geo_mean) {
    df = data.frame(row.names=rownames(mat_ds), gene=rownames(mat_ds), mu1=apply(mat_ds[, nms1], 1, function(y) {exp(mean(log(1+y)))-1}), mu2=apply(mat_ds[, nms2], 1, function(y) {exp(mean(log(1+y)))-1}), stringsAsFactors = F)
    df$tot1 = df$mu1 * length(nms1)
    df$tot2 = df$mu2 * length(nms2)
  } else {
    df = data.frame(row.names=rownames(mat_ds), gene=rownames(mat_ds), tot1=Matrix::rowSums(mat_ds[, nms1]), tot2=Matrix::rowSums(mat_ds[, nms2]), stringsAsFactors = F)
  }
  
  norm_by = min(sum(df$tot1), sum(df$tot2))
  df$tot1 = df$tot1 / sum(df$tot1) * norm_by
  df$tot2 = df$tot2 / sum(df$tot2) * norm_by
  df = df[pmax(df$tot1, df$tot2) >= min_max_umi, ]
  
  if (geo_mean && geo_mean_per_cell) {
    df$enr = log2( (df$mu1 + reg) / (df$mu2 + reg))
  } else {
    df$enr = log2( (df$tot1 + reg) / (df$tot2 + reg))
  }
  df = df[order(df$enr, decreasing=T), ]
  
  if (filter_outlier_genes) {
    fp = mc@mc_fp[intersect(rownames(mc@mc_fp), df$gene), ]
    if (!is.null(mcs1)) {
      fp = fp[, mcs1]
    }
    gmax = apply(fp, 1, max)
    gnext = apply(fp, 1, function(v) head(tail(sort(v), n=compare_top_mc_to_n_highest), n=1) )
    df[rownames(fp), 'out_r'] =  gmax/gnext
    to_filt = !is.na(df$out_r) & df$out_r > max_top_to_n_highest_ratio & df$enr > 0
    if (sum(to_filt) > 0) {
      if (verbose) {
        message(sprintf("filtering %d outlier genes:", sum(to_filt)))
      }
      print(df[to_filt, ])
      df = df[!to_filt, ]
    }
  }
  df	        
}



