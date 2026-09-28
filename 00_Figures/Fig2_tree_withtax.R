.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(data.table)
library(treeio)

# load tree generated with GTDB-Tk and IQTREE as described in the methods section of the manuscript
tree <- read.tree('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/classification_tree/tree/250429_zoogloea_gtdb_midas_HQ_R226/gtdbtk.bac120.user_msa.fasta.treefile')
tr <- phytools::midpoint.root(tree)

# load metadata with proposed names for the genomes
propnames <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/metadata/proposed_R226_new_names_updated.tsv', col_names = T)  
tr$tip.label <- gsub("_SRR.*|_ASM.*|_TT.*|.fa", "", tr$tip.label)
propnames[] <- lapply(propnames, function(x) if(is.factor(x) || is.character(x)) gsub("Zoogloea", "Z.", as.character(x)) else x)

pp <- ggtree(tr ,size = 0.6 ) 
p.genome <-pp %<+% propnames 

# Put bootstrap values into a new dataframe and assign them to a status category for coloring the nodes in the tree
df_boot <- p.genome %>% as.treedata %>% as_tibble
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 90), "90-100 %",
                         ifelse((df_boot$label >= 70),"70-90 %","0-70 %"))
df_boot <- df_boot[, c("node","status")]


w <- p.genome %<+% df_boot + geom_nodepoint(aes(color=status), size=1.75) +
  scale_color_manual(values = c("90-100 %" = "grey0", "70-90 %"="gray50", "0-70 %"="gray75"),
                     na.translate = FALSE, drop = F) +
  labs( color="Bootstrap") + xlim(-0.15, NA)+
  theme(legend.position = c(0.275, 0.50))


ps <- w + new_scale_color() +
  geom_treescale(x=0, y=20, offset=1) 
ps <- ps %>% 
     rotate(99)



ps$data$GTDBr226_species <- gsub("s__","", ps$data$GTDBr226_species)
ps$data$midas_species <- gsub("_"," ", ps$data$midas_species)
ps$data$label <- gsub("GCA_009026175.1","GCF_009026175.1", ps$data$label)
ps$data$label <- gsub('barcode02_Z_caeni','GCA_986281895*', ps$data$label)
ps$data$label <- gsub('barcode03_Z_oleivorans','GCA_986281795*', ps$data$label)
ps$data$label <- gsub('barcode04_Z_resiniphila','GCA_986340685*', ps$data$label)


text.size=2.8

rem <- c("s:|;")
xq <- 0.237

ps_new <- ps +
  geom_segment(aes(x=0.244, y=y, xend=0.52, #yend=ylim, 
                   color = if_else(is.na(GTDBr226_species), 'g__Thaura',proposed_name_new)),
               data=ps$data %>% filter(isTip),
               alpha = 1, linewidth=4.1)+
  geom_tippoint(mapping=aes(x= x+0.0028,color = if_else(is.na(GTDBr226_species), 'g__Thaura',proposed_name_new)),
                size=3, stroke=0, shape=15) +
  geom_tiplab(mapping = aes(label = label),
              geom='text',size = text.size, linesize = 0, offset = 0.009,lineheight=0, linetype='blank',
              alpha = 1, align = T)+
  geom_tiplab(mapping = aes(label = GTDBr226_species,
                            fontface = if_else(str_detect(GTDBr226_species, 'sp'), 'plain', 'italic')), 
              geom='text',
              size = text.size, linesize = 0, offset = 0.097,lineheight=0,linetype='blank',
              
              alpha = 1, align = T) +
  geom_tiplab(mapping = aes(label = str_remove_all(midas_species,rem),
                            fontface = if_else(str_detect(midas_species, 'Z.'), 'italic', 'plain'))
              , geom='text',
              size = text.size, linesize = 0, offset = 0.17,lineheight=0,linetype='blank',
              alpha = 1, align = T) +
  geom_tiplab(mapping = aes(label = proposed_name_new), geom='text',
              fontface = 'italic',
              size = text.size, linesize = 0, offset = 0.23,lineheight=0,linetype='blank',
              alpha = 1, align = T) +
  geom_text(data = data.frame(
    x = c(xq+0.008, xq+0.096, xq+0.169, xq+0.229),  # the same as your offsets
    y = max(ps$data$y) + 1.5,       # place headers *above* the top tip
    label = c("Accession", "GTDB R226", "MIDAS 5.3", "Proposed Name")
  ),
  aes(x = x, y = y, label = label),
  inherit.aes = FALSE,
  size = 3, fontface = "bold", hjust = 0)



p.genomec <- collapse(ps_new, node = 65, 'min', fill='grey75')

col <- c('transparent','#e8d9c8','#e8d9c8','grey80','#e8d9c8','grey80','#e8d9c8','grey80','#e8d9c8','grey80','#e8d9c8','grey80','#e8d9c8','#e8d9c8','grey80','grey80','#e8d9c8','grey80','#e8d9c8','grey80','grey80')

s <- p.genomec  +
  scale_color_manual('', values = col
                     , guide = "none"
                     )

ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/Fig2_total_tree_wtax_r226.tiff',
       s, dpi=1200,
       compression = "lzw",
       height = 8,
       width=15)

ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/Fig2_total_tree_wtax_r226.jpeg',
       s, dpi=400,
       height = 8,
       width=15)


ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/Fig2_total_tree_wtax_r226.pdf',
       s, device = cairo_pdf,
       height = 8,
       width=13)

