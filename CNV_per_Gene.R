
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(org.Hs.eg.db)
library(rhdf5)


###### Step1 : Loading & Filtering ##### 
###### load node_unmerged_cnv_calls.bed & cnv_data.h5 from Cellrnager-dna (10X Genomics) output ######
### Load bed_file ###
setwd("/path_to_working_directory/")

bed = read.table("node_unmerged_cnv_calls.bed", header=FALSE, sep="\t")
colnames(bed) = c("chrom", "start", "end", "id", "copy_number", "event_confidence")
bed=bed[which(bed$chrom=="chrXX"), ]  ##### chrXX correspond to the gene of interest chr locus, ie for TNFRSF17 this would be chr16 ######

### Create data.frame with col = groups and rows = cells ###
h5 = H5Fopen("cnv_data.h5")

group = as.data.frame(h5$"/tree/is_cell_in_group", header=FALSE)
colnames(group) = c(as.numeric(h5$"/constants/num_cells"):(length(unique(bed$id))-1))
colnames(group) = c(as.numeric(h5$"/constants/num_cells"):(as.numeric(h5$"/constants/num_cells")*2-2))
rownames(group) = c(0:(as.numeric(h5$"/constants/num_cells")-1))

#### Exclude very large nodes ### (arbitrarily set up at over 80% of all cells)

group = group[,!colSums(group)>=(as.numeric(h5$"/constants/num_cells")*0.8)]

##### Step 2: Convert bins into genes + translation to SYMBOL #####

cnv_gr <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
g <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only=FALSE)

olaps <- findOverlaps(g, cnv_gr)

df_gene <- tibble(entrezgene = names(g)[queryHits(olaps)], 
                  copy_number = mcols(cnv_gr)$copy_number[subjectHits(olaps)],
                  id = mcols(cnv_gr)$id[subjectHits(olaps)],
                  event_confidence = mcols(cnv_gr)$event_confidence[subjectHits(olaps)])

entrezgene_symbol_map <- as.list(org.Hs.egSYMBOL)
entrezgene_symbol_map <- lapply(entrezgene_symbol_map, `[`, 1)

df_gene <- dplyr::filter(df_gene, entrezgene %in% names(entrezgene_symbol_map)) %>% 
  dplyr::mutate(gene_id = unlist(entrezgene_symbol_map[entrezgene])) %>% 
  dplyr::select(gene_id, entrezgene, copy_number, id, event_confidence) %>% 
  drop_na()
df_gene = as.data.frame(df_gene)

###### Step3 : CNV analysis #####

#Recover nodes with each cell
cells_nodes = apply(group==1,1,function(a) as.integer(colnames(group))[a])

#studied_gene, XXXXX corresponf to the gene studied (for example TNFRSF17)
studied_gene = df_gene[which(df_gene$gene_id=="XXXXX"), c("id","copy_number", "event_confidence")]

#keep, for each cell, the event with the highest confidence.
output = data.frame("id"= 0, 
                    "copy_number" = studied_gene[which(studied_gene$id %in% c(0, cells_nodes[["0"]])),][which.max(studied_gene[which(studied_gene$id %in% c(0,cells_nodes[["0"]])),]$event_confidence), "copy_number"], 
                    "event_confidence" = studied_gene[which(studied_gene$id %in% c(0, cells_nodes[["0"]])),][which.max(studied_gene[which(studied_gene$id %in% c(0,cells_nodes[["0"]])),]$event_confidence), "event_confidence"])

for(i in 1:(as.numeric(h5$"/constants/num_cells")-1)) {
  copy_number = studied_gene[which(studied_gene$id %in% c(i ,cells_nodes[[i+1]])),]
  output[i+1,"id"]=i
  output[i+1,"copy_number"]= copy_number[which.max(copy_number$event_confidence), "copy_number"]
  output[i+1,"event_confidence"]= copy_number[which.max(copy_number$event_confidence), "event_confidence"]
}

##### Step 3 : Filtering #####

### Take out noisy cells ###

noisy_df = as.data.frame(h5$"/per_cell_summary_metrics/is_noisy", header=FALSE)
colnames(noisy_df) = "noisy"
noisy_df$id = c(0:(as.numeric(h5$"/constants/num_cells")-1))

noisy_cells = noisy_df[which(noisy_df$noisy==1),"id"]

output_f = output[!(rownames(output) %in% noisy_cells),]

### Filter out normal cells based on ploidy

normal_df = as.data.frame(h5$"/per_cell_summary_metrics/mean_ploidy", header=FALSE)
colnames(normal_df) = "ploidy"
normal_df$id = c(0:(as.numeric(h5$"/constants/num_cells")-1))

normal_df = normal_df[which(normal_df$ploidy>=1.95 & normal_df$ploidy <=2.05 ),"id"]
output_f_nn = output_f[!(rownames(output_f) %in% normal_df),]


############## Last step : visualisation (example below for TNFRSF17 gene)
MM_TNFRSF17 = output_f_nn
MM_TNFRSF17$sample = "TNFRSF17"
MM_TNFRSF17[which(MM_TNFRSF17$copy_number>=4),"copy_number"] = ">=4"

library(ggplot2)
library(tidyr)
library(patchwork)


dloupe.palette <- c("0"="#2570B1", "1"="#5fb1ff", "2"="#ffffff", "3"="#F0823C", ">=4"="#C31718")
level_order_cnv <- c(">=4", "3", "2", "1", "0")

table(MM_TNFRSF17$copy_number)
Prop_table <- table(MM_TNFRSF17$copy_number)
prop.table(Prop_table)
prop.table(Prop_table) %>% '*' (100) %>% round(1)
Copy_Number <- as.vector(table(MM_TNFRSF17$copy_number))
TNFRSF17_prop_table <- prop.table(Prop_table) %>% '*' (100) %>% round(1)
write.csv(TNFRSF17_prop_table, "/path_to_directory_to_store_csv/TNFRSF17_prop_table.csv")

p1 = MM_TNFRSF17 %>%
  
  ggplot(aes (x = factor(sample), fill = factor(copy_number, level_order_cnv), y = (..count..)/sum(..count..))) +
  geom_bar(position="fill", colour = "black", size=1) +
  # coord_flip() +
  scale_fill_manual(values=dloupe.palette) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  theme_classic(base_size = 16) +
  guides(fill = guide_legend(title = "Copy Number")) +
  theme(legend.position="left") +
  scale_x_discrete(limits=c("TNFRSF17")) +
  xlab("") +
  ylab("Cells (%)")


p1 