#!/usr/bin/env Rscript
####################################################
#Goal:		CPM and proportion bar plot for sample-specific isoforms
#Author:	Yang Xu
#E-mail:	yangax@pennmedicine.upenn.edu
####################################################

args = commandArgs(trailingOnly=TRUE)
library(ComplexHeatmap)
library(circlize)
#library(ggplot2)
#library(RColorBrewer)
#library(dplyr)
#library(viridis)

target_trans = args[1]
out_dir = args[2]
window_size = args[3]
start_pos = as.numeric(args[4])
end_pos = as.numeric(args[5])
inf_name = paste(out_dir, '/7_3_calculate_score_matrix_WindowSize_', window_size,'_', target_trans, '_with_topology_z_score.txt',sep='')

fig_data <- read.table(inf_name, header = T,sep='\t', row.names=1, check.names = FALSE)
outf_name = paste(out_dir, '/7_4_heatmap_', target_trans,'_score_line_z_score.png', sep='')
outf_pdf_name = paste(out_dir, '/7_4_heatmap_', target_trans, '_score_line_z_score.pdf', sep = '')
#dim(fig_data)

###### subset specific range #######
fig_data_trans <- t(fig_data)
fig_data_trans <- fig_data_trans[as.numeric(rownames(fig_data_trans)) >= start_pos & as.numeric(rownames(fig_data_trans)) <= end_pos, ]
fig_data <- t(fig_data_trans)
#dim(fig_data)


topology_row <- as.character(fig_data[(row.names(fig_data) %in% c('Location')),])
#mean_row <- as.numeric(fig_data[(row.names(fig_data) %in% c('Log2FC_mean')),])
median_row <- as.numeric(fig_data[(row.names(fig_data) %in% c('Log2FC_median')),])
fig_data <- fig_data[!(row.names(fig_data) %in% c('Location','Log2FC_median','Log2FC_mean')),]
#dim(fig_data)


df <- data.frame(matrix(as.numeric(unlist(fig_data)), nrow=dim(fig_data)[1], byrow=FALSE), stringsAsFactors=FALSE)
rownames(df) <- rownames(fig_data)
colnames(df) <- colnames(fig_data)

max_value = ceiling(max(max(df),-min(df)))
max_value_ha = ceiling(max(as.numeric(as.character(median_row))))
min_value_ha = floor(min(as.numeric(as.character(median_row))))
font_size = 6

#length(median_row)

#col_fun_ha = colorRamp2(c(-max_value_ha, 0, max_value_ha), c('#1b998b', "white", '#f46036'))

#col_fun = colorRamp2(c(0, max_value/2, max_value), c("#235789", "white", "#C1292E"))  #blue red
#col_fun = colorRamp2(c(0, max_value/2, max_value), c("white", "yellow", "red"))    #isoCirc
#col_fun = colorRamp2(c(0, max_value/4, max_value/2, max_value*3/4, max_value), viridis(5, option="magma"))  #viridis
#col_fun = colorRamp2(c(0, max_value/4, max_value/2, max_value*3/4, max_value), heat.colors(5))  #heat
col_fun = colorRamp2(c(-max_value, 0, max_value), c("#235789", "white", "#C1292E"))
col_fun_ha = colorRamp2(c(-max_value_ha, 0, max_value_ha), c("#235789", "white", "#C1292E"))

if (sum(is.na(topology_row)) == length(topology_row)){
ha = HeatmapAnnotation(
    Specificity = anno_lines(median_row, gp=gpar(col = '#8E44AD', lwd=2), add_points=FALSE, pt_gp=gpar(col = '#8E44AD'), pch=c(1), height = unit(0.8, "in"), ylim = c(min_value_ha, max_value_ha), extend = 0.05),
	annotation_name_gp = gpar(fontsize = font_size+1),
	simple_anno_size_adjust = TRUE,
	annotation_name_side = "left"
)
}else{
ha = HeatmapAnnotation(
	Specificity = anno_lines(median_row, gp=gpar(col = '#8E44AD', lwd=2), add_points=FALSE, pt_gp=gpar(col = '#8E44AD'), pch=c(1), height = unit(0.8, "in"), ylim = c(min_value_ha, max_value_ha), extend = 0.05),
	Topology = topology_row,
	col = list(
		Topology = c("Extracellular" = "#F6D7A7", "TM" = "#C8E3D4", "Cytoplasmic"="#87AAAA")
	),
	annotation_legend_param = list(
		Topology = list(
			title = "Topology",
			labels = c('Extracellular','TM','Cytoplasmic'),
			at = c('Extracellular','TM','Cytoplasmic'),
			labels_gp = gpar(fontsize = font_size),
			title_gp = gpar(fontsize = font_size),
			rect_gp = gpar(col = "white", lwd = 0.5),
			legend_direction = "horizontal",
			color_bar = "discrete",
			nrow = 1,
			title_position = "topleft",
			legend_width = unit(1.5, "in")
		)
	),
	annotation_name_gp = gpar(fontsize = font_size+1),
	height = unit(0.8, "in"),
	simple_anno_size_adjust = TRUE,
	annotation_name_side = "left"
	#gp = gpar(col = "white", lwd = 0.01)
)
}

#lgd = Legend(col_fun = col_fun, title = "Tumor-specificity score", direction = "horizontal", title_position = "lefttop")
#rect_gp = gpar(col = "white", lwd = 0.5)

f_fig = Heatmap(df,
	heatmap_legend_param = list(
		title = "Z-score of expression", title_gp = gpar(fontsize = font_size), labels_gp = gpar(fontsize = font_size), color_bar = "continuous", legend_direction="horizontal", title_position = "topleft", legend_width = unit(1.5, "in")
		), 
		rect_gp = gpar(col = "white", lwd = 0.5), col = col_fun, column_title = "Amino acid position along the protein", column_title_side = "bottom", row_title="Samples", column_title_gp = gpar(fontsize = font_size+2, fontface = "bold"), row_title_gp = gpar(fontsize = font_size+2, fontface = "bold"), cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", column_names_side = "bottom", show_column_names=FALSE, row_names_gp = gpar(fontsize = font_size), column_names_gp = gpar(fontsize = font_size+1), column_names_rot = 90, top_annotation = ha)

fig_2 <- draw(f_fig, heatmap_legend_side = "top", annotation_legend_side = "top", merge_legend =TRUE)
png(file=outf_name, width=4, height=7, units='in',res=1200)
#draw(f_fig, heatmap_legend_side = "top", annotation_legend_side = "top", merge_legend =TRUE)
fig_2
dev.off()

pdf(file=outf_pdf_name, width=4, height=7)
#draw(f_fig, heatmap_legend_side = "top", annotation_legend_side = "top", merge_legend =TRUE)
fig_2
dev.off()

#ggsave(outf_pdf_name, plot=f_fig, width=8, height=6, dpi=300)
