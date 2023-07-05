#!/usr/bin/env Rscript
####################################################
#Goal:		CPM and proportion bar plot for sample-specific isoforms
#Author:	Yang Xu
#E-mail:	yangax@pennmedicine.upenn.edu
####################################################

args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(viridis)

target_trans = args[1]
out_dir = args[2]
window_size = args[3]
request_pos = args[4]
inf_name = paste(out_dir, '/7_5_barplot_WindowSize_', window_size,'_', target_trans, '_reshaped.txt',sep='')

fig_data <- read.table(inf_name, header = T,sep='\t')
outf_name = paste(out_dir, '/7_6_Bar_', target_trans,'_exp.png', sep='')
outf_pdf_name = paste(out_dir, '/7_6_Bar_', target_trans, '_exp.pdf', sep = '')
fig_data$Sample <- factor(fig_data$Sample, levels=c(rev(unique(as.character(fig_data$Sample)))))
fig_data$Group <- factor(fig_data$Group, levels=c('Tumor','Tissue'))
dim(fig_data)
fig_data <- fig_data[fig_data$AA_index == as.numeric(request_pos), ]
dim(fig_data)


font_size = 7
f_fig <- ggplot(fig_data, aes(x=Sample, y=CPM))
f_fig <- f_fig + geom_bar(aes(fill=Group), color='black', size=0.1, stat='identity')
f_fig <- f_fig + scale_fill_manual(name='Group', values=c('#ffae57','#b877bd'), breaks=c('Tumor','Tissue'),labels=c('Tumor','Tissue'))
#f_fig <- f_fig + scale_colour_manual(name="Subtype",values=c("#0000F9","#02a620","#e68e19","#FF0018"),labels=c("Luminal","HER2 enriched","Basal A","Basal B"),breaks=c("Luminal","HER2_Amp","Basal_A","Basal_B"))
f_fig <- f_fig + theme(panel.grid.major = element_line(color = "white",size=0.01),panel.grid.minor = element_line(color = "white",size=0.01), axis.title.y=element_text(size = font_size+1), axis.text.y=element_text(size = font_size+2), axis.title.x=element_text(size = font_size), axis.text.x=element_text(size = font_size+1), plot.title = element_blank(), legend.position = 'top', legend.title = element_text(size = font_size), legend.text = element_text(size = font_size), panel.border = element_rect(fill=NA, size=0.3, colour = "black"), panel.background = element_blank(), legend.box="vertical", legend.direction="horizontal", plot.margin = unit(c(5,5,5,5), "pt"))
f_fig <- f_fig + labs(x="Sample", y="Sum of the expression of transcripts encoding peptide\nwith the highest specificity score (CPM)", title="")
#f_fig <- f_fig + scale_y_continuous(limits=c(0,8), breaks=c(0,2,4,6,8))
f_fig <- f_fig + coord_flip()
f_fig <- f_fig + guides(fill=guide_legend(title.position="top", title.hjust =0))


png(file=outf_name, width=3, height=8, units='in',res=300)
f_fig
dev.off()
ggsave(outf_pdf_name, plot=f_fig, width=3, height=8, dpi=300)

