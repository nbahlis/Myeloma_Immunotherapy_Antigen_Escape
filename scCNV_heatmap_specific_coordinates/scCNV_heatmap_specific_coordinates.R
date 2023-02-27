# this script takes heatmap csv exports from loupe cnv and saves figures
library(tidyverse)
library(ggpubr)
library(cowplot)

dloupe.palette <- c("0"="#2570B1", "1"="#8BAFD2", "2"="#F0F0F0", "3"="#FBD5B9", "4"="#F9BC8B", 
                    "5"="#F0823C", "6"="#F0823C", "7"="#EB562B", "8"="#E52A19", "9"="#D42119",
                    "10"="#C31718", "11"="#A21814", "12"="#801810", "13" = "#61170F", "14"="#41150E", "NaN"="#CDCCCC")



# FILES
# cnv files are extracted from 10X Genomiocs Loupe scDNA Browser after uploading cellranger dna processed .dloupe file 
#and selecting genomic coordinate of interest (ie TNFRSSF17 as gene or genomic coordinates chr16:11860001-chr16:12260000)
# and exporting heatmap copy numver data as CSV

cnv.files.path <- list.files(path="~/Path_to_CNV_plots/", pattern = ".*.csv", full.names = TRUE)

# SAMPLE LABELS

cnv.files.names <- str_extract(string = basename(cnv.files.path), pattern = "\\w+")


# load heatmap files and annotate
if(all(file.exists(cnv.files.path))) {
  heatmap.files <- lapply(cnv.files.path, read_csv)
} else {
  stop("Not all heatmap files were found")
}

# rename list objects
if(length(heatmap.files)==length(cnv.files.names)) {
  names(heatmap.files) <- cnv.files.names
} else {
  stop("there should be as many names as heatmap files")
}

# pivot to longer format before binding together only if they all match the expected format
# also annotate position for graphical representations
if(all(unlist(lapply(heatmap.files, function(x) c("node_id", "barcodes", "num_cells", "num_noisy") %in% colnames(x))))) {
  long.cnv.data <- heatmap.files %>%
    map(~ ( mutate(.x,
                   loupe_index=row_number(), 
                   loupe_cell_number=cumsum(num_cells),
                   loupe_cell_before=lag(loupe_cell_number, default = 0),
                   total_cells=sum(num_cells)))) %>%
    map(~ ( pivot_longer(.x, 
                         cols = -c(node_id, barcodes, num_cells, num_noisy, loupe_index, loupe_cell_number, loupe_cell_before, total_cells), 
                         names_to = "segment",
                         values_to = "copy_number"))) %>%
    bind_rows(.id = "sample_name")
} else {
  stop("check heatmap file format / must match expected from loupe output")
}

# translate segment name to genomic coordinates, number of cells to sample, and window size
long.cnv <- long.cnv.data %>%
  separate(segment, into = c("chr", "start", "end"), sep = "[:-]", remove = FALSE) %>%
  mutate_at(.vars = vars(start, end), str_remove_all, pattern=",") %>% 
  mutate_at(.vars = vars(start, end), as.numeric) %>%
  group_by(chr) %>%
  mutate(chrn=str_remove(chr, pattern = "chr"),
         sample_name_number_cells=paste0(sample_name, " (n=", format(total_cells, big.mark = ","), " cells)"),
         window_coord=paste0(chr, ":", format(min(start), big.mark=","), "-", format(max(end), big.mark=","))) %>%
  ungroup()

# ggplot object
# key is that tile height should correspond to cluster size (=N cells)
genomic.window.plot <- ggplot(long.cnv, aes(xmin=start, xmax=end, ymin=loupe_cell_before, ymax=loupe_cell_number, fill=as.factor(copy_number))) +
  geom_rect() +
  scale_fill_manual(values=dloupe.palette) +
  scale_x_continuous(labels = function(x) format(x, big.mark=","), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(sample_name_number_cells ~ window_coord, space = "free_y", scales = "free_y", switch = "y") +
  theme(legend.position = "none", panel.grid = element_blank(), panel.background = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(fill = NA), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank())



# I make my own legend
gglegend <- data.frame(CN=names(dloupe.palette), HEX=dloupe.palette) %>%
  dplyr::filter(CN %in% c(0:5)) %>%
  ggplot(aes(1, CN, fill=CN)) +
  scale_fill_manual(values=dloupe.palette)+
  geom_tile() +
  theme_pubclean() +
  labs(x="Copy\nnumber")+
  theme(aspect.ratio = 3, axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_blank(), legend.position = "none", axis.title.x = element_text(size=9))


library(ggbio)
library(EnsDb.Hsapiens.v86)
ensdb <- EnsDb.Hsapiens.v86
# below an example how to plot heatmap for chr16 (11860001, 12260000)
gr <- GRanges(seqnames = 16, IRanges(11860001, 12260000), strand = "+")
p2 <- autoplot(ensdb, GRangesFilter(gr), names.expr = "gene_name") + theme_null()
p1 <- genomic.plot.with.legend <- plot_grid( genomic.window.plot, gglegend, p2@ggplot, rel_widths = c(6, 1))
p1
