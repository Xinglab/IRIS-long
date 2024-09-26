args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(123)

input_gene = args[1]
outf_dir = args[2]
inf_name = paste(outf_dir, "/", input_gene, "_TCGA_GTEx_gene_exp.txt", sep="")
file_keyword <- input_gene

fig_data <- read.table(inf_name, header = T, sep='\t')
outf_name = paste(outf_dir, "/Gene_exp_box_",file_keyword,".png", sep="")
outf_pdf_name = paste(outf_dir, "/Gene_exp_box_",file_keyword,".pdf", sep="")


fig_data$Type <- factor(fig_data$Type, level=c("TCGA","GTEx"), labels=c('TCGA', 'GTEx'))
fig_data_2 <- fig_data[order(fig_data$Type, fig_data$Group),]
fig_data$Group <- factor(fig_data$Group, level=as.character(unique(fig_data_2$Group)))
#fig_data <- fig_data[fig_data$SE_PSI>0,]

font_size = 8
f_A <- ggplot(fig_data, aes(x=Group, y=Normalized_count, fill=Type))
#f_A <- f_A + geom_violin(trim=TRUE)
f_A <- f_A + geom_boxplot(outlier.size = 0.1)
f_A <- f_A + theme(panel.grid.major = element_line(color = "white",size=0.01),panel.grid.minor = element_line(color = "white",size=0.01), plot.title=element_text(size=font_size, hjust=0.5), axis.text.x=element_text(size = font_size-2, angle=90, vjust=0.5, hjust=1),axis.title.x = element_blank(), axis.title.y = element_text(size = font_size), axis.text.y=element_text(size=font_size, colour = 'black'), legend.title = element_text(size = font_size-1), legend.text = element_text(size = font_size-2), panel.border = element_rect(fill=NA, size=0.3, colour ="black"), panel.background = element_blank(),legend.key = element_rect(fill = "white", color = NA),legend.position='right')
f_A <- f_A + labs(x="Group", y="Expression of gene\n(DESeq2 normalized count)", title=paste("Expression of ",file_keyword, " in TCGA and GTEx samples", sep=''))
#f_A

def_width = 8.5
if(length(as.character(unique(fig_data_2$Group))) < 40){
	def_width = 6
}
png(file=outf_name, width=def_width, height=3, units='in',res=1500)
f_A
dev.off()
ggsave(outf_pdf_name, plot=f_A, width=def_width, height=3, dpi=300)
