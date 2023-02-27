library(stringr)
library(tidyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(rhdf5)
library(stringr)
library(karyoploteR)
library(ggplot2)
library(ggpubr)
library(cowplot)
setwd("path to working directory")

source("path to /utils.R")

# load barcodes from dloupe manual selection
heatmap.csv.dir <- c("path to file heatmap csv files from 10X Genomic dloupe")
heatmap.csv.files <- list.files(path = heatmap.csv.dir, pattern = "*heatmap_copy number.csv")
heatmap.csv.data <- lapply(paste0(heatmap.csv.dir, heatmap.csv.files), read.csv)
heatmap.csv.IDs <- str_extract(heatmap.csv.files, "P[[:digit:]]{4}") # (samples are coded as P followed by 4 digits, ie P2510)
names(heatmap.csv.data) <- heatmap.csv.IDs
heatmap.barcodes <- lapply(heatmap.csv.data, function(x){ unlist( strsplit(as.character(pull(x, barcodes)), ";") ) })

dir.create("cnv-barcode-filtered")
mapply(function(id, csv){ 
  write.table(csv, file = paste0("cnv-barcode-filtered/", id, "_CNV_barcodes.txt"), quote = F, row.names = F, col.names = F)
}, id=heatmap.csv.IDs, csv=heatmap.barcodes)


# load cnv data from hdf5r
cnv.data.root <- "path to directory with 10xCNV cnv_data.h5 files"
cnv.h5.files <- paste0(cnv.data.root, heatmap.csv.IDs, "_CNV/outs/cnv_data.h5")
all( file.exists(cnv.h5.files) )
# TRUE


cnv.data <- mapply(function(h5.cnv, heatmap.barcodes){ get.cnv.data(h5.cnv, 100, heatmap.barcodes) }, cnv.h5.files, heatmap.barcodes)


# rows: genomic coordinates; columns: cell_id
cnv.matrix <- do.call(cbind, cnv.data)
length(unlist(heatmap.barcodes)) == ncol(cnv.matrix)

# define metadata
samples.list.expanded <- rep( names( sapply(heatmap.barcodes, length) ), sapply(heatmap.barcodes, length))
cnv.short.metadata <- data.frame(Sample=heatmap.csv.IDs, 
                                 ID=c("Patient","Patient"), 
                                 Cond=c("Pre", "Post"))
cnv.short.metadata$ID.Cond <- paste(cnv.short.metadata$ID, cnv.short.metadata$Cond, sep="-")
cnv.cell.metadata <- data.frame(row.names = paste(colnames(cnv.matrix), samples.list.expanded, sep="_"), Sample=samples.list.expanded)
cnv.cell.metadata <- cnv.cell.metadata %>% tibble::rownames_to_column("cell_id") %>% left_join(cnv.short.metadata, by = c("Sample"="Sample"))
colnames(cnv.matrix) <- cnv.cell.metadata$cell_id

cnv.pos.metadata <- data.frame(row.names = rownames(cnv.matrix), 
                               position = rownames(cnv.matrix),
                               chr=str_split_fixed(rownames(cnv.matrix), pattern = ":|-", n = 3)[,1],
                               start=str_split_fixed(rownames(cnv.matrix), pattern = ":|-", n = 3)[,2],
                               end=str_split_fixed(rownames(cnv.matrix), pattern = ":|-", n = 3)[,3])

#### now we have cnv.matrix with the data, cnv.pos.metadata with genomic loc. and cnv.cell.metadata with sample/condition info

saveRDS(cnv.matrix, "cnv.matrix.rds")
saveRDS(cnv.pos.metadata, "cnv.pos.metadata.rds")
saveRDS(cnv.cell.metadata, "cnv.cell.metadata.rds")
# 
# cnv.matrix <- readRDS("cnv.matrix.rds")
# cnv.pos.metadata <- readRDS("cnv.pos.metadata.rds")
# cnv.cell.metadata <- readRDS("cnv.cell.metadata.rds")

# Step2 -------------------------------------------------------------------

library(tidyr)
library(dplyr)
library(stringr)
library(karyoploteR)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(cowplot)

Chr.palette <- rep(brewer.pal(8, "Greys")[c(5, 2)], length.out=24)
Condition.palette <- brewer.pal(3, "Set2")[c(1, 2)]

Patient.ID <- "Patient"
Patient.ID.Pre <- paste0(Patient.ID, "-Pre")
Patient.ID.Post <- paste0(Patient.ID, "-Post")

sectors <- paste0("chr", 1:22)

cnv.matrix <- readRDS("cnv.matrix.rds")
cnv.pos.metadata <- readRDS("cnv.pos.metadata.rds")
cnv.cell.metadata <- readRDS("cnv.cell.metadata.rds")

cnv.pos.metadata$position <- as.character(cnv.pos.metadata$position)
cnv.pos.metadata$start <- as.numeric(as.character( cnv.pos.metadata$start ))
cnv.pos.metadata$end <- as.numeric(as.character(cnv.pos.metadata$end))

dloupe.palette <- c("0"="#2570B1", "1"="#8BAFD2", "2"="#F0F0F0", "3"="#FBD5B9", "4"="#F9BC8B", 
                    "5"="#F0823C", "6"="#F0823C", "7"="#EB562B", "8"="#E52A19", "9"="#D42119",
                    "10"="#C31718", "11"="#A21814", "12"="#801810", "13" = "#61170F", "14"="#41150E", "NaN"="#CDCCCC")

cell_ids_pre = subset(cnv.cell.metadata, ID.Cond == Patient.ID.Pre) %>% pull(cell_id)
cell_ids_post = subset(cnv.cell.metadata, ID.Cond == Patient.ID.Post) %>% pull(cell_id)
cell_ids <- c(cell_ids_pre, cell_ids_post)
positions = subset( cnv.pos.metadata) %>% pull(position) # chr %in% c("chr18", "chr19")
reduced.matrix <- cnv.matrix[ positions, cell_ids]
reduced.matrix <- round(reduced.matrix)

CNV.mat.reduced.df <- as.data.frame( t( reduced.matrix ) ) %>% 
  tibble::rownames_to_column("cell_id") %>% 
  gather(Pos, Copies, -cell_id) %>%
  separate(Pos, into = c("Chr", "Start", "End"), sep="[^[:alnum:]]+", remove = FALSE) %>%
  separate(cell_id, into = c("CELL_ID", "LIBRARY"), sep="_", remove = FALSE)
CNV.mat.reduced.df$Chr <- factor(CNV.mat.reduced.df$Chr, levels=sectors)
CNV.mat.reduced.df$Start <- as.numeric(CNV.mat.reduced.df$Start)
CNV.mat.plot <- ggplot(CNV.mat.reduced.df, aes(Start, factor(cell_id, levels = rev(cell_ids)), fill=as.factor(Copies))) + 
  geom_tile() +
  scale_fill_manual(values=dloupe.palette) +
  scale_x_continuous(expand = c(0,0))+
  facet_grid(LIBRARY~Chr, scales = "free", space = "free_x") + #facet by group
  # theme_pubclean() + 
  theme(legend.position = "none",
        axis.title = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text = element_blank(), 
        panel.grid = element_blank(),
        panel.border=element_rect(fill = NA, size=0.1, colour="lightgrey"),
        strip.background = element_blank(), strip.text = element_blank(),
        panel.spacing = unit(0.1, "lines"), plot.margin=unit(c(0, 0, 0, 0), "cm"))

CNV.chr.df <- CNV.mat.reduced.df %>% dplyr::select(Start, Chr) %>% group_by(Start, Chr) %>% distinct()
CNV.chr.plot <- ggplot(CNV.chr.df, aes(Start, 1, fill=Chr)) +
  geom_tile() +
  scale_fill_manual(values=Chr.palette) +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  facet_grid(.~ factor(Chr, labels=""), scales = "free", space = "free") +  #facet by group
  # theme_pubclean() +
  theme(legend.position = "none", 
        axis.title = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.text = element_blank(),
        panel.border=element_rect(fill = NA,size=0.1, colour="lightgrey"),
        panel.grid = element_blank(), panel.background = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 5),
        panel.spacing = unit(0.1, "lines"), plot.margin=unit(c(0, 0, 0, 0), "cm"))

CNV.sample.df <- CNV.mat.reduced.df %>% dplyr::select(cell_id, LIBRARY) %>% distinct()
CNV.sample.plot <- ggplot(CNV.sample.df, aes(1, factor(cell_id, levels = rev(cell_ids)), fill=LIBRARY)) +
  geom_tile() +
  scale_fill_manual(values=Condition.palette) +
  facet_grid(factor(LIBRARY, labels=c(Patient.ID.Pre, Patient.ID.Post))~., scales = "free_y", switch = "y") +  #facet by group
  # theme_pubclean()+  
  scale_x_continuous(expand = c(0,0))+
  theme(legend.position = "none", axis.title = element_blank(),panel.grid = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.text = element_blank(),
        panel.border=element_rect(fill = NA, size=0.1, colour="lightgrey"),
        strip.background = element_blank(), strip.text = element_text(size = 5, angle = 180),
        panel.spacing = unit(0.1, "lines"), plot.margin=unit(c(0, 0, 0, 0), "cm"))

Plot.Put.Together <- ggdraw(xlim = c(0, 116), ylim=c(0,58), clip = "on") +   
  draw_plot(CNV.mat.plot, x=7, y=0, width = 108, height = 51) +
  draw_plot(CNV.sample.plot, x=0, y=0, width = 7, height = 51) + 
  draw_plot(CNV.chr.plot, x=7, y=51, width = 108, height = 7)
# Plot.Put.Together
ggsave(Plot.Put.Together, filename = paste0("Plot.Put.Together.", Patient.ID, ".png"), width = 116, height = 58, units = "mm", dpi = 600)
