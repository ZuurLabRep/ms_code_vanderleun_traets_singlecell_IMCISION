### general functions for analysis and plotting ###

# load packages and metacell object names
tbk_reload = function() {
  library("metacell")
  library("dplyr")
  library("Matrix")
  library("tidyr")
  library("ggplot2")
  library("data.table")
  library("pheatmap")
  
  if(!dir.exists("data")) dir.create("data/") 
  scdb_init("data/", force_reinit=T)
  
  if(!dir.exists("figures")) dir.create("figures/") 
  scfigs_init("figures/")
  
  exp_nm <<- "HN_full_dataset"
  mat_id <<- exp_nm
  clean_mat_id <<- "HN_full_dataset_clean_v3"
  mc_id<<- paste0(exp_nm, "_mc")
  filt_mc_id<<- paste0(mc_id,"_f")
  mc2d_id <<- paste0(mc_id, "_2dproj")
  
  exp_nm_kurten <<- "Kurten_full_dataset"
  mat_id_kurten <<- exp_nm_kurten
  clean_mat_id_kurten <<- "Kurten_full_dataset_clean"
  mc_id_kurten<<- paste0(exp_nm_kurten, "_mc")
  filt_mc_id_kurten<<- paste0(mc_id_kurten,"_f")
  mc2d_id_kurten <<- paste0(mc_id_kurten, "_2dproj")
}

# function to plot 2D projection of cells with their corresponding cell state color
# input: mc2d projection metacell object, clean matrix metacell object, metacell info object, colors assigned to cell states
# output: 2D projection plot of cells colored by cell state
plot_coord_mc_colors = function(mc2d_name,mat_name,mc_name,mc_colors){
  mat = scdb_mat(mat_name)
  mc = scdb_mc(mc_name)
  mc2d = scdb_mc2d(mc2d_name)
  
  # coordinates cells
  xco <- data.frame(mc2d@sc_x, cellnms = names(mc2d@sc_x))
  yco <- data.frame(mc2d@sc_y, cellnms = names(mc2d@sc_y))
  coo <- merge(xco, yco, by="cellnms")
  coo <- na.omit(coo)
  
  cells <- rownames(mat@cell_metadata)
  cells <- cells[cells %in% names(mc@mc)]
  cell_codes <- intersect(cells, coo$cellnms)
  
  # annotate cells with colors of mc groups
  # replace " " with "_" to match color annotation file
  coo_plot <- as.data.frame(coo)
  rownames(coo_plot) = coo_plot$cellnms
  coo_plot <- coo_plot[,-1]
  coo_plot <- coo_plot[intersect(rownames(coo_plot), cell_codes),]
  coo_plot[cell_codes, "mc_group"] = mat@cell_metadata[cell_codes, "mc_group"]
  coo_plot$mc_group <- gsub(" ", "_", coo_plot$mc_group)
  colnames(mc_colors)[3] <- "mc_group"
  coo_plot$cell.id <- rownames(coo_plot)
  coo_colors <- plyr::join(coo_plot,mc_colors,by="mc_group")
  
  coo_colors$mc_group <- factor(coo_colors$mc_group,levels=mc_colors$mc_group)
  
  # plot mc_groups
  pdf(file=paste0(fig_path,"/2dprojection_mc_groups.pdf"), width=6.7, height=4)
  p = ggplot()+
    ggtitle(label = sprintf("mc2d_mc_groups")) + 
    geom_point(data = coo_colors, aes(x=mc2d.sc_x, y = mc2d.sc_y, fill = mc_group),colour = "black",pch=21, size=1.3,stroke = 0.1)+
    scale_fill_manual(values=mc_colors$color) +guides(colour = guide_legend(override.aes = list(size=3)))+
    theme_bw() + theme(legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    guides(fill = guide_legend(override.aes = list(size=2.5,stroke=0.3)))
  print(p)
  dev.off()
}

# function to plot gene expression of one gene on 2D projection of cells from metacell objects
# input: mc2d projection metacell object, clean matrix metacell object, prepared clean matrix, metacell info object, gene of interest
# output: gene expression projected on 2D projection plot of cells
plot_coord_gene_expression = function(mc2d_name,mat_name,umis_n,mc_name,gene_oi){
  mc = scdb_mc(mc_name)
  mat = scdb_mat(mat_name)
  submat = mat@cell_metadata[names(mc@mc),]
  mc2d = scdb_mc2d(mc2d_name)

  # gene expression
  genes_cells <- as.data.frame(umis_n[which(rownames(umis_n) %in% gene_oi),])
  genes_cells$gene_oi <- rowSums(genes_cells)
  genes_cells$lognorm <- log2(genes_cells$gene_oi + 1)
  
  md <- mat@cell_metadata
  md <- md[!is.na(md$mc_group),]
  genes_cells$cell.id <- rownames(genes_cells)
  md$cell.id <- rownames(md)
  
  all_merged <- merge(genes_cells,md,by="cell.id")

  # coordinates
  xco <- data.frame(mc2d@sc_x, cellnms = names(mc2d@sc_x))
  yco <- data.frame(mc2d@sc_y, cellnms = names(mc2d@sc_y))
  coo <- merge(xco, yco, by="cellnms")
  coo <- na.omit(coo)
  
  # align with meta data
  cells <- rownames(mat@cell_metadata)
  cells <- cells[cells %in% names(mc@mc)]
  cell_codes <- intersect(cells, coo$cellnms)
  
  # prepare coordinates for plotting
  coo_c <- as.data.frame(coo)
  rownames(coo_c) <- coo_c$cellnms
  coo_c <- coo_c[,-1]
  coo_c <- coo_c[intersect(rownames(coo_c), cell_codes),]
  
  coo_c$cell.id <- rownames(coo_c)
  to_plot <- merge(coo_c,all_merged,by="cell.id")
  
  png(file=paste0(fig_path,"/2dprojection_gene_",gene_oi,".png"), width=1000, height=1000)
  p = ggplot()+geom_point(data = to_plot[to_plot$lognorm == 0,], aes(x=mc2d.sc_x, y = mc2d.sc_y),fill = "lightgrey",pch=21, size=6,stroke = 0) +
    geom_point(data = to_plot[to_plot$lognorm !=0,], aes(x=mc2d.sc_x, y = mc2d.sc_y, fill = lognorm),colour = "black",pch=21, size=6,stroke = 0) +
    scale_fill_gradient2(low="#ffe7ca", mid="#f48622", high="#3f2800",midpoint = max(to_plot$lognorm)/2, 
                         limits = c(0, max(to_plot$lognorm)))+
    theme(legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank())
  print(p)
  dev.off()
}


# function to plot 2D projection of cells with their corresponding cell state color, baseline RE vs NR
# input: mc2d projection metacell object, clean matrix metacell object, metacell info object, colors assigned to cell states
# output: two 2D projection plots of cells colored by cell state, baseline in RE and NR
plot_coord_mc_colors_baseline_resp = function(mc2d_name,mat_name,mc_name,mc_colors){
  mat = scdb_mat(mat_name)
  mc = scdb_mc(mc_name)
  
  mc2d = scdb_mc2d(mc2d_name)
  xco = data.frame(mc2d@sc_x, cellnms = names(mc2d@sc_x))
  yco = data.frame(mc2d@sc_y, cellnms = names(mc2d@sc_y))
  coo = merge(xco, yco, by="cellnms")
  coo = na.omit(coo)
  
  cells = rownames(mat@cell_metadata)
  
  cells = cells[cells %in% names(mc@mc)]
  cell_codes = intersect(cells, coo$cellnms)
  
  coo_df = as.data.frame(coo)
  rownames(coo_df) = coo_df$cellnms
  coo_df = coo_df[,-1]
  coo_df = coo_df[intersect(rownames(coo_df), cell_codes),]
  coo_df[cell_codes, "mc_group"] = mat@cell_metadata[cell_codes, "mc_group"]
  coo_df[cell_codes, "condition"] = mat@cell_metadata[cell_codes, "condition"]
  coo_df[cell_codes, "response"] = mat@cell_metadata[cell_codes, "response"]
  coo_df[cell_codes, "patient"] = mat@cell_metadata[cell_codes, "patient"]
  coo_df$mc_group <- gsub(" ", "_", coo_df$mc_group)
  colnames(mc_colors)[3] <- "mc_group"
  coo_df$cell.id <- rownames(coo_df)
  coo_colors <- plyr::join(coo_df,mc_colors,by="mc_group")
  
  coo_colors$condition_response <- paste0(coo_colors$condition,"_",coo_colors$response)
  
  yes_colors_pre_NR <- coo_colors[coo_colors$condition_response %in% c("pre_NR"),]
  yes_colors_pre_RE <- coo_colors[coo_colors$condition_response %in% c("pre_RE"),]
  
  yes_colors_pre_RE$mc_group <- factor(yes_colors_pre_RE$mc_group,levels=mc_colors$mc_group)
  
  # plot condition/response
  pdf(file=paste0(fig_path,"/2dprojection_mc_groups_pre_RE_fix_v2.pdf"), width=6.7, height=4)
  p = ggplot()+
    ggtitle(label = sprintf("mc2d_mc_groups_RE")) + 
    geom_point(data = yes_colors_pre_NR, aes(x=mc2d.sc_x, y = mc2d.sc_y),size=1,color="#e3e3e3")+
    geom_point(data = yes_colors_pre_RE, aes(x=mc2d.sc_x, y = mc2d.sc_y, fill = mc_group),colour = "black",pch=21, size=1.2,stroke = 0.1)+
    scale_fill_manual(values=mc_colors$color) +guides(colour = guide_legend(override.aes = list(size=3)))+
    theme_bw() + theme(legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  print(p)
  dev.off()
  #write.table(yes_colors_pre_RE,"./HN_manuscript_figures/2dprojection_mc_groups_pre_RE_fix_v2.txt",quote=F,sep="\t")
  
  yes_colors_pre_NR$mc_group <- factor(yes_colors_pre_NR$mc_group,levels=mc_colors$mc_group)
  
  pdf(file=paste0(fig_path,"/2dprojection_mc_groups_pre_NR_fix_v2.pdf"), width=6.7, height=4)
  p = ggplot()+
    ggtitle(label = sprintf("mc2d_mc_groups_NR")) + 
    geom_point(data = yes_colors_pre_RE, aes(x=mc2d.sc_x, y = mc2d.sc_y),size=1,color="#e3e3e3")+
    geom_point(data = yes_colors_pre_NR, aes(x=mc2d.sc_x, y = mc2d.sc_y,  fill = mc_group),colour = "black",pch=21, size=1.2,stroke = 0.1)+
    scale_fill_manual(values=mc_colors$color) +guides(colour = guide_legend(override.aes = list(size=3)))+
    theme_bw() + theme(legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  print(p)
  dev.off()
  #write.table(yes_colors_pre_NR,"./HN_manuscript_figures/2dprojection_mc_groups_pre_NR_fix_v2.txt",quote=F,sep="\t")
}


# function to plot gene expression per cell state, on baseline
# input: clean matrix metacell object, metacell info object, prepared clean matrix, gene of interest, colors assigned to cell states
# output: gene expression in violin plots per cell state
plot_violin_gene_expression = function(mat_name=clean_mat_id,umis_n=umis_n,mc_name=filt_mc_id,gene_oi=gene_oi,mc_colors=mc_colors){
  genes_cells <- as.data.frame(umis_n[which(rownames(umis_n) %in% gene_oi),])
  
  genes_cells$gene_oi <- rowSums(genes_cells)
  genes_cells$lognorm = log2(genes_cells$gene_oi + 1) # log2 of gene exression
  
  md = mat@cell_metadata
  genes_cells$cell.id <- rownames(genes_cells)
  md$cell.id <- rownames(md)
  
  # only baseline
  md = md[md$condition == "pre",]
  
  all_merged <- merge(genes_cells,md,by="cell.id")
  
  mc_colors$name <- gsub("_", " ", mc_colors$name)
  all_merged$mc_group = factor(all_merged$mc_group,levels=mc_colors$name)
  
  # violin plot
  pdf(file=paste0(fig_path,"/Expression_",gene_oi,"_baseline_response.pdf"), width=10, height=3.5)
  p=ggplot(all_merged, aes(x = mc_group, y = lognorm, fill = mc_group)) +
    geom_violin( adjust =1,trim=TRUE, scale = "width")    +   
    geom_point( aes(x = mc_group, y = lognorm, fill = mc_group),position = position_jitterdodge(dodge.width=1,jitter.width=3.3),size=0.6)+
    scale_fill_manual(values=mc_colors$color)+
    theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(p)
  dev.off()
}


# function to plot cell type proportions, based in response and timepoint
# input: clean matrix metacell object, metacell info object, condition (timepoint), response, order
# output: boxplot with celltype proportions
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
  mc_nr_r<- reshape2::melt(t(as.matrix(n_cells_mc_wide)))
  colnames(mc_nr_r)[3] = "tot_nr"
  mc_prop_r$patient <- factor(mc_prop_r$Var2,levels=unique(mc_prop_r[["Var2"]]))
  mc_prop_nr_r = plyr::join(mc_prop_r,mc_nr_r)
  
  pdf(file= paste0(fig_path,"/Ordered_patients_barplots_",response,"_",condition,".pdf"), width=4.5, height=3.2)
  p=ggplot(mc_prop_nr_r,aes(x=patient,y=value,fill=Var1,label=tot_nr)) + geom_bar(position="fill", stat="identity") +
    ggtitle(paste0(response,"_",condition," CD8 ordered >")) +  
    geom_text(position = position_stack(vjust = 0.5),
              color="black", size=3.5) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  print(p)
  dev.off()
}


# function differentially expressed genes between group 1 and  2, based on wilcoxon
# input: group 1, group 2
# output: volcano plot
DE_wilcoxon = function(gr1_df,gr2_df){
  gr1_cells = rownames(gr1_df)
  gr2_cells = rownames(gr2_df)
  
  umis = as.matrix(mat_ds[,c(gr1_cells, gr2_cells)])
  umis_n = sweep(umis,2,colSums(umis), "/") * 1000 # scaled
  
  # Fold-change
  logFC = as.vector(apply(umis, 1, function(x) log2((mean(x[gr1_cells]))/(mean(x[gr2_cells])))))
  logFCx = as.vector(apply(umis, 1, function(x) log2(mean(x[gr1_cells]))))
  logFCy = as.vector(apply(umis, 1, function(x) log2(mean(x[gr2_cells]))))
  
  # Wilcoxon-Mann-Whitney test
  wilc_adj = p.adjust(as.vector(apply(umis_n, 1, function(x) wilcox.test(x[gr1_cells], x[gr2_cells])$p.value)), method = "fdr")
  wilc = apply(umis_n, 1, function(x) wilcox.test(x[gr1_cells], x[gr2_cells])$p.value)
  
  pb <- data.frame(rownames(umis), logFC, logFCx, logFCy, wilc, wilc_adj)
  colnames(pb) <- c("gene", "fc", "x", "y", "pval", "pval_adj")
  pb$pval[pb$pval=="NaN"] <- NA
  pb <- na.omit(pb)
  pb <- pb[-which(abs(pb$fc) == "Inf"),]
  rownames(pb) <- pb$gene
  
  pb$pval[pb$pval == 0] <- 10e-280
  return(pb) # return data frame with wilcoxon output and fold change
}

# function to plot line boxplot, RE and NR
# input: , ylim, combined data frame with fractions, data frame RE, data frame NR
# output: line boxplot
make_line_boxplot = function(cell_x, ylim_ax,prepost_frac_c,prepost_final_RE,prepost_final_NR){
  cell_x_counts <- prepost_frac_c[prepost_frac_c$mc_group == cell_x,]
  cell_x_counts$pat_cond <- paste0(cell_x_counts$patient,cell_x_counts$condition)
  prepost_final_RE$pat_cond <- paste0(prepost_final_RE$patient,prepost_final_RE$condition)
  prepost_final_RE_plot <- merge(prepost_final_RE,cell_x_counts,by= "pat_cond")
  pdf(file=paste0(fig_path,"/",cell_x,"_post_pre_RE.pdf"), width=2.5, height=2.3) #
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
  pdf(file=paste0(fig_path,"/",cell_x,"_post_pre_NR.pdf"), width=2.5, height=2.3) #
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


# function calculate diff expressed genes between two groups, using geom mean
# input: mc, downsampled mat, groups to compare
# output: data frame with diff expressed genes
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

# function to plot UMAP with the expression of indv genes in log2 scale
# input: monocle subset, scaled matrix, gene of interest
# output: UMAP with gene expression indv gene
plot_UMAP_expression = function(cds_subset,umis_n,gene_oi,cell_type){
  genes_cells <- as.data.frame(umis_n[which(rownames(umis_n) %in% gene_oi),])
  
  genes_cells$gene_exp <- rowSums(genes_cells)
  genes_cells$lognorm = log2(genes_cells$gene_exp + 1)
  md = cds_subset@colData
  
  genes_cells$cell.id <- rownames(genes_cells)
  md$cell.id <- rownames(md)
  all_merged <- merge(genes_cells,md,by="cell.id")
  
  xco = umap_1
  yco = umap_2
  coo = as.data.frame(cbind(xco, yco))
  coo$cell.id = umap_cell
  coo$mc = umap_mc
  
  to_plot = as.data.frame(merge(coo, all_merged, by = "cell.id"))
  
  pdf(file=paste0(fig_path,"/2dprojection_gene_",gene_oi,"_Tregs.pdf"),width=4.6, height=4.5)
  p = ggplot()+geom_point(data = to_plot[to_plot$lognorm == 0,], aes(x=xco, y = yco),fill = "#440154FF",pch=21, size=2.5,stroke = 0) +
    geom_point(data = to_plot[to_plot$lognorm !=0,], aes(x=xco, y = yco, fill = lognorm),colour = "black",pch=21, size=2.5,stroke = 0) +
    scale_fill_gradient2(low="#440154FF", mid="#21908CFF", high="#ffe000",midpoint = max(to_plot$lognorm)/2, 
                         limits = c(0, max(to_plot$lognorm)))+
    theme(legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank())
  print(p)
  dev.off()
  
  print(cell_type)
  if (cell_type == "treg"){
    col_cells = c("#b3d4f0", "#a14f00","#107538","#0042b0","#59b56f")
  }
  if (cell_type == "dysf"){
    col_cells = c("#E5D3B3", "#4C005C","#C20088","#ff7f7f","#c8a2c8")
  }
  
  pdf(file=paste0(fig_path,"/2dprojection_gene_mc_group_for_legend_Tregs.pdf"), width=6, height=5)
  p = ggplot()+geom_point(data = to_plot, aes(x=xco, y = yco,fill=mc_group,color=mc_group),pch=21, size=2.5,stroke = 0) +
    theme(legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank()) +
          scale_fill_manual(values=col_cells) +
          scale_color_manual(values=col_cells)
  print(p)
  dev.off()
}


# function to generate NK cell signatures based on:
# Yang, C., et al. Heterogeneity of human bone marrow and blood natural killer cells defined by single-cell transcriptome. Nat Commun 10, 3931 (2019)
# output: two gene signatures based on a CD56bright, mature & active signatures of NK cells
generate_NKcell_signatures  = function(){
  blood_bright = read.table("./data/signatures/CD56bright_blood_NK.txt")
  bm_bright = read.table("./data/signatures/CD56bright_BM_NK.txt")
  CD56_bright = unique(c(blood_bright$V1,bm_bright$V1))
  CD56_bright = CD56_bright[-grep("RPL",CD56_bright)]
  CD56_bright = CD56_bright[-grep("RPS",CD56_bright)]
  
  bm_mature = read.table("./data/signatures/Mature_BM_NK.txt")
  blood_mature = read.table("./data/signatures/Mature_blood_NK.txt")
  bm_active = read.table("./data/signatures/Active_blood_NK.txt")
  blood_active = read.table("./data/signatures/Active_BM_NK.txt")
  active_mature = unique(c(bm_mature$V1,blood_mature$V1,bm_active$V1,blood_active$V1))
  
  return(list("CD56_bright" = CD56_bright,"active_mature" = active_mature))
}



# function to load Adaptive data
# input: working directory with the adaptive data
# output: large df with all TCRs (including template counts)
load_adpative_data_bulkTCR  = function(working_dir){
  couple_samples <- read.csv("./data/CFMPB559 koppeling CFnrs & AKPnrs.csv")
  
  working_dir_old <- getwd()
  setwd(working_dir)
  
  sample_IDs <- system("ls", intern=T)
  
  couple_samples <- couple_samples[c(1:18),]
  couple_samples$Total_blood_filtered <- NA
  couple_samples$Productive_simpson_clonality <- NA
  
  all_templates_samples <- as.data.frame(matrix(ncol = 3,nrow=0))
  colnames(all_templates_samples) <- c("cdr3_rearrangement","Patient","condition")
  
  for (x in 1:18){
    print(sample_IDs[x])
    
    sample_TCR = read.table(sample_IDs[x],header=T,sep="\t")
    
    sample_id <- do.call(rbind,strsplit(sample_IDs[x],"_"))[,1]
    sample_inf <- couple_samples[couple_samples$CFSample == sample_id,]
    sample_inf$Patient <- paste0("Pat",sample_inf$Patient)
    sample_TCR$Patient <- sample_inf$Patient
    sample_TCR$condition <- sample_inf$Timepoint
    
    sample_TCR <- sample_TCR[sample_TCR$cdr3_rearrangement != "",]
    sample_TCR <- sample_TCR[sample_TCR$v_gene != "",]
    sample_TCR <- sample_TCR[sample_TCR$j_gene != "",]
    sample_TCR <- sample_TCR[,c("cdr3_rearrangement","templates","Patient","condition")]
    all_templates_sample <- sample_TCR |> uncount(templates)
    
    all_templates_samples <- rbind(all_templates_samples,all_templates_sample)
  }
  all_templates_samples$test <- paste0(all_templates_samples$Patient,"_",all_templates_samples$condition)
  setwd(working_dir_old)
  
  return (all_templates_samples)
}


# function to plot TCR coverage within metacells
# input: mat_id, path for the figure
# output: figure showing clone counts per metacell
TCR_coverage = function(mat_id, width = 1600, heigth = 800) {
  mat = scdb_mat(mat_id)
  y = vector(length = length(unique(mat@cell_metadata$mc))-1)
  names(y) = 1:(length(unique(mat@cell_metadata$mc))-1)
  names_y = vector()
  tot_cells = vector()
  for (i in unique(mat@cell_metadata$mc)){
    df = mat@cell_metadata[which(mat@cell_metadata$mc == i),]
    y[i] = nrow(df[which(df$clonotype_id != "NA"),])/nrow(df)
    #names_y = na.omit(c(names_y, i))
    tot_cells[i] = nrow(mat@cell_metadata[which(mat@cell_metadata$mc == i),])
  }
  plot = y[order(mc@colors)]
  
  metacell:::.plot_start(scfigs_fn(dir="figures/Figure_S1/",mat_id, sprintf("clonotype_counts_per_metacell")), width, heigth)
  barplot(plot, col = mc@colors[as.numeric(names(plot))], cex.names = 0.5, las = 2, main = "clonotype_counts_per_metacell", ylab = "fraction_clonotype/total")
  mtext(tot_cells, 3, cex = 0.7, at=seq(1, by=1.2, length=length(tot_cells)), las=2, line=0.5)
  dev.off()
}



