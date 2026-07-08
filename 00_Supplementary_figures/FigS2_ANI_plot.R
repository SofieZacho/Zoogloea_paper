.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

# Load the data
data <- read_tsv("/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/classification_tree/ANI/250122_fastani/fastANI.tsv", col_names = c("query", "reference", "ANI", "alignment_coverage_nr", "total_query_fragments"))

# remove '/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/symlinks/Zoogloea_thau_HQ_midas_gtdb/' from column 1 and 2 - also remove '.fa'
data$query <- gsub("/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/symlinks/Zoogloea_thau_HQ_midas_gtdb/", "", data$query)
data$reference <- gsub("/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/symlinks/Zoogloea_thau_HQ_midas_gtdb/", "", data$reference)
data$query <- gsub("_SRR.*|_ASM.*|_TT.*|.fa", "", data$query)
data$reference <- gsub("_SRR.*|_ASM.*|_TT.*|.fa", "", data$reference)

propnames <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/metadata/proposed_R226_new_names_updated.tsv', col_names = T)  
propnames[] <- lapply(propnames, function(x) if(is.factor(x) || is.character(x)) gsub("Zoogloea", "Z.", as.character(x)) else x)
propnames <- propnames %>% select(genome, proposed_name_new)

# load metadata
# load meta data and combine and manipulate
metagtdb <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/metadata/meta_thau_zoo_manipulated.tsv')
metamidas <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/metadata/meta_zoo_midas.tsv')

mg <- metagtdb %>% select(accession, Completeness = checkm2_completeness, Contamination = checkm2_contamination,
                          Genome_Size=genome_size, contigs = contig_count, MAG_status = mimag_high_quality, 
                          tax_gtdb = gtdb_taxonomy) %>% 
  mutate(MAG_status = if_else(MAG_status == T, 'HQ','MQ'),
         tax_midas = NA, source = 'gtdb')
mm <- metamidas %>% select(accession = bin, Completeness, Contamination,Genome_Size,
                           contigs, MAG_status, tax_gtdb, tax_midas) %>% mutate(source = 'MiDAS_catalog')

meta <- rbind(mg, mm) %>%
  separate(tax_gtdb, into = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ';') %>% 
  rename(tip_label = accession) %>%
  # if 'barcode' is part of tip label, make source = 'This study'
  mutate(source = if_else(str_detect(tip_label, 'barcode'), 'This study', source),
         genome = tip_label) %>% 
  mutate(species = if_else(genome %in% 'barcode02_Z_caeni',
                           's__', # as it was manipulated to have that taxonomy
                           species))



# merge data and metadata



# plot
colspp <- c( "grey90",'black','black','black', "#EAC5AA", "#84541E","#D18E35", "#C6A07E","#f2dfcb", "#E0C37E", "#80CDC1", "#34978F",  "#4E8077", "#B7A27C", "#95825A", "#9c7605")
cols <- c("#F5F5F5",'#f7ecdf',"#f2dfcb", "#EAC5AA", "#80CDC1", "#34978F",  "#4E8077")
ggplot(data) +
  geom_tile(aes(x=query, y=reference, fill = ANI)) +
  scale_fill_gradientn(colors = cols)


# put it on ggtree
t_genome <- '/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/classification_tree/tree/250120_zoogloea_gtdb_midas_HQ/MSA_zoogloea_HQ_IQtree.faa.treefile'
treeg <- read.tree(t_genome)
trg <- phytools::midpoint.root(treeg)
############################################################
############ plot full genome tree ########################
#remove everything in tip labels after '_SRR', '_ASM', and also remove 'RS_' in front of accession
trg$tip.label <- gsub('_SRR.*|_ASM.*|_TT.*|RS_', '', trg$tip.label)
p <- ggtree(trg ,size = 0.5 ) 

p.genome <- p %<+% meta 
p.genomec <-collapse(p.genome, node = 65, 'max', fill='grey80')

ps <- p.genomec+ geom_treescale() +
  geom_tiplab(mapping = aes(label = paste0(species,'::',genome),
                            fill = species),
              geom='label',
              size = 2, linesize = 0.2, offset = 0.01,
              alpha = 0.8, align = T)+theme(legend.position = 'none') +
  scale_fill_manual(values = colspp)



datap <- data %>% rename(tip_label = query)
datap$reference <- factor(datap$reference, levels = ps$data %>% filter(isTip) %>% select(label,y) %>% arrange(y)%>% pull(label)) 

thau = c('GCF_002245655.1','GCF_003030465.1','GCF_000443165.1')
datap <- datap %>% filter(!reference %in% thau,
                          !tip_label %in% thau)

ps + new_scale_fill() + geom_fruit(data= datap,
                                   geom = geom_tile,
                                   mapping = aes(
                                     y=tip_label, 
                                     x=reference,
                                     group = label,
                                     fill=ANI
                                   ),
                                   offset= 0.2, pwidth = 1.2,
                                   axis.params=list(
                                     axis       = "x",
                                     text.size  = 1.8,
                                     hjust      = 0,
                                     vjust      = 0.5, text.angle = 90
                                   )#width=0.8 #before -0.03
) + 
  scale_fill_gradientn(colors = cols) +
  theme(legend.position = 'right') 

ps$data$label <- gsub('barcode02_Z_caeni','GCA_986281895*', ps$data$label)
ps$data$label <- gsub('barcode03_Z_oleivorans','GCA_986281795*', ps$data$label)
ps$data$label <- gsub('barcode04_Z_resiniphila','GCA_986340685*', ps$data$label)



####################################
############ only heatmap ####
dataheat <- datap
dataheat <- dataheat %>% left_join(., propnames, by = c('tip_label'='genome')) %>% 
  left_join(., propnames, by = c('reference'='genome'))



dataheat$tip_label <- gsub('barcode02_Z_caeni','GCA_986281895*', dataheat$tip_label)
dataheat$tip_label <- gsub('barcode03_Z_oleivorans','GCA_986281795*', dataheat$tip_label)
dataheat$tip_label <- gsub('barcode04_Z_resiniphila','GCA_986340685*',dataheat$tip_label)

dataheat$reference <- gsub('barcode02_Z_caeni','GCA_986281895*', dataheat$reference)
dataheat$reference <- gsub('barcode03_Z_oleivorans','GCA_986281795*', dataheat$reference)
dataheat$reference <- gsub('barcode04_Z_resiniphila','GCA_986340685*',dataheat$reference)


# 1. Get the correct tree-based order
correct_order <- ps$data %>%
  filter(isTip) %>%
  select(label, y) %>%
  left_join(., propnames, by = c('label'='genome')) %>% 
  arrange(y) %>%
  pull(label)

correct_order_species <- ps$data %>%
  filter(isTip) %>%
  select(label, y) %>%
  left_join(., propnames, by = c('label'='genome')) %>% 
  arrange(y) %>%
  pull(proposed_name_new) %>% unique()

# 2. Set factors for tip_label and reference based on correct order
dataheat <- dataheat %>%
  mutate(
    tip_label = factor(tip_label, levels = correct_order),
    reference = factor(reference, levels = correct_order),
    proposed_name_new.x = factor(proposed_name_new.x, levels = correct_order_species),
    proposed_name_new.y = factor(proposed_name_new.y, levels = rev(correct_order_species))
  )

dataheat$tip_label <- gsub("GCA_009026175.1","GCF_009026175.1",dataheat$tip_label)
dataheat$reference <- gsub("GCA_009026175.1","GCF_009026175.1",dataheat$reference)

heat_clean <- ggplot(dataheat) +
  geom_tile(aes(x = tip_label, y = reference, fill = ANI), color = 'grey40', size = 0.021) +
  scale_fill_gradientn(colors = cols) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=7, color = 'black'),
    axis.text.y =element_text(size=7, color = 'black'),
    strip.text.x = element_text(angle = 90, hjust = 0, face = 'italic', color = 'black'),
    strip.text.y = element_text(angle = 0, hjust=0, face = 'italic', color = 'black'),
    panel.spacing = unit(0.05, "lines"),
    panel.grid = element_blank()
  ) +
  facet_grid(rows =vars(proposed_name_new.y),
             cols = vars(proposed_name_new.x),
             scales = 'free', space = 'free')
heat_clean

ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS2_ANI_propnames_facet.jpeg',
       heat_clean, dpi=600,
       height = 9,
       width=10)







