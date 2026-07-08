.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")
library(tidyverse)
library(ggtree)



# data paths
zoogloea <- read_tsv('/home/bio.aau.dk/zs85xk/projects/epsSMASH/epsSMASH_analyses/Zooglan_Sofie/region_counts.tsv', skip = 1)
zoogloea$record <- str_remove_all(zoogloea$record, "_ASM.*|_SRR.*|_TT.*|\\.fa")


# #hits$fasta <- str_remove_all(hits$fasta, "_TT.*")
# rethink <- read_tsv('/home/bio.aau.dk/zs85xk/projects/epsSMASH/epsSMASH_analyses/REThiNk_catalogue/data/epsSMASH_results/region_counts.tsv', skip = 1)
# rethink$record <- gsub('.gbff','',rethink$record)

hit.genomes <- read_tsv('~/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/AS_MAGs/250224_An/hit_genomes.txt',
                        col_names = 'genome')


# read trees
propnames <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/metadata/proposed_R226_new_names_updated.tsv') %>% 
  select(genome, proposed_name = proposed_name_new)
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
  mutate(source = if_else(str_detect(tip_label, 'barcode'), 'This study', source),
         genome = tip_label) %>% 
  mutate(species = if_else(genome %in% 'barcode02_Z_caeni',
                           's__', # as it was manipulated to have that taxonomy
                           species))
meta.zoogloea <- meta %>% left_join( ., propnames, by = 'genome')






# filter to relevsnt genomes
rethink.hit.genomes <- rethink %>% 
  filter(record %in% hit.genomes$genome)

meta <- read_tsv('/projects/PHN/MiDAS/04-Analysis/Status-analysis/103-samples/100-add_checkm1/03-MAGs.final-version/000-MAGs-metadata/MiDAS-22277.hqMAGs.metadata.tsv')%>%
  separate(tax_gtdb.r220, into = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ';')
# 
# # remove columns with no hits (meaning column sum is 0) - only for numeric columns
# rethink.hit.genomes.filt <- rethink.hit.genomes %>% 
#   select(where(~ !is.numeric(.) || (is.numeric(.) && sum(.) > 0))) %>% 
#   select(#-contains("putative"),
#          -description)
# 
# 
# lv <-names(colSums(rethink.hit.genomes.filt %>% select(where (~ is.numeric(.) ))) %>% sort())
# 
# rethink.hit.genomes.filt.long <- rethink.hit.genomes.filt %>% 
#   pivot_longer(., cols = -record) %>% 
#   mutate(name = factor(name, levels = rev(lv)))
# 
# rethink.hit.genomes.filt.long.meta <- rethink.hit.genomes.filt.long %>% 
#   left_join(., meta, by = c('record'='bin'))
# 
# pp <-rethink.hit.genomes.filt.long.meta %>% filter(!name %in% 'total_count') %>% 
# ggplot() +
#   geom_tile(aes(x= paste0( genus, species,record), y = name, fill = factor(value)),
#             size= 0.02, color = 'grey20')+
#   scale_fill_manual('Number',values = c('grey95','#EAC5AA', '#d98c6b')) +
#   theme_minimal()+
#   theme(
#     axis.text.x = element_text(size=3.5, angle = 90, hjust=1, vjust= 0.5),
#     axis.title = element_blank(),
#     strip.text.x = element_text(angle=90, hjust=0, face = 'bold', size = 8),
#     strip.clip = 'off'
#   )+
#   facet_grid(cols=vars( family), scales = 'free', space= 'free' )
# 
# 
# ggsave('~/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/AS_MAGs/epsSMASH/000_plots/hit_genomes_epshits.jpeg',
#        pp,
#        width = 13,
#        height = 5.5)
# 


# only zoogloea
# remove columns with no hits (meaning column sum is 0) - only for numeric columns
zoogloea.filt <- zoogloea %>% 
  select(#-contains("putative"),
    -description, -hybrid) 



zoogloea.filt.long <- zoogloea.filt %>% 
  pivot_longer(., cols = -record) 

putative_sum <- zoogloea.filt.long %>%
  filter(grepl("putative", name, ignore.case = TRUE)) %>%
  group_by(record) %>%
  summarise(name = "putative", value = sum(value), .groups = "drop")



zoogloea.filt.long.meta <- zoogloea.filt.long %>%
  filter(!grepl("putative", name, ignore.case = TRUE)) %>%
  bind_rows(putative_sum) %>%
  arrange(record) %>% 
  left_join(., meta.zoogloea, by = c('record'='tip_label')) %>% 
  filter(!genus %in% 'g__Thauera')%>% filter(!name %in% 'total_count') %>% 
  mutate(name = factor(name, levels = c('putative','pel','zooglan-like')))

zoogloea.filt.long.meta$record <- gsub('barcode02_Z_caeni','GCA_986281895*', zoogloea.filt.long.meta$record)
zoogloea.filt.long.meta$record <- gsub('barcode03_Z_oleivorans','GCA_986281795*', zoogloea.filt.long.meta$record)
zoogloea.filt.long.meta$record <- gsub('barcode04_Z_resiniphila','GCA_986340685*',zoogloea.filt.long.meta$record)


p_tile <- ggplot(zoogloea.filt.long.meta) +
  geom_tile(aes(x= record, y = name, fill = factor(value)),
            size= 0.02, color = 'grey20')+
  geom_text(aes(x= record, y = name, label= value),
            size=2.5)+
  scale_fill_manual(values = c('grey98','#EAC5AA', '#d98c6b','#cf6c42')) +
  facet_grid(cols = vars(proposed_name),
             space = 'free', scales = 'free')+
  theme_minimal()+
  theme(
    legend.background = element_rect(fill = 'transparent'),
    axis.title = element_blank(),
    axis.text.y = element_text(size=8, angle =0, hjust=1, vjust = 0.5, face='bold', color = 'black'),
    axis.text.x = element_text(size=6, angle =90, hjust=1, vjust = 0.5), 
    strip.text.x = element_text(angle=90, hjust=0, face = 'bold', size = 6),
    strip.clip = 'off',
    #strip.text.x = element_text(size=10, face = 'bold'),
    # panel.spacing.y=unit(0.045, "lines"),
    # panel.spacing.x=unit(0.7, "lines"),
    panel.grid = element_blank(),
    legend.position = 'none'
    #strip.background = element_rect(fill = 'grey90')
  ) 


p_tile



library(gggenes)


bigscape <- read_tsv('~/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/Zoogloea/250402_cluster_gene_clusters/bigscape/250403_bigscape/output_files/2025-04-03_16-17-52_full.network')


gene_clusters <- read_tsv('/home/bio.aau.dk/zs85xk/projects/epsSMASH/epsSMASH_analyses/Zooglan_Sofie/gene_info_zoogloea.tsv', col_names = T)
gene_functions <- read_tsv('/home/bio.aau.dk/zs85xk/projects/epsSMASH/epsSMASH/epssmash/outputs/html/gene_functions.tsv', 
                           col_names = c('Query','gene_function'))


# read trees
propnames <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/metadata/proposed_R226_new_names_updated.tsv', col_names = T)  
propnames[] <- lapply(propnames, function(x) if(is.factor(x) || is.character(x)) gsub("Zoogloea", "Z.", as.character(x)) else x)
propnames <- propnames %>% select(genome, proposed_name_new)


### only putative
gene_clusters.put <- gene_clusters %>% 
  filter(str_detect(product, "putative")) %>%
  mutate(
    start = str_extract(location, "(?<=\\[)\\d+"),
    stop = str_extract(location, "(?<=:)\\d+"),
    strand = str_extract(location, "(?<=\\()[-+]") # Extracts "+" or "-"
  ) %>%
  mutate(across(c(start, stop), as.integer)) 

gene_clusters.put <- gene_clusters.put %>% 
  mutate(GBK=paste0(genome_id, "__",contig_id,".",region_id))
#putative.GBK <- test$GBK %>% unique()



gc.put.filt <- gene_clusters.put %>% 
  group_by(genome_id, contig_id, region_id, gene_id) %>%
  slice_min(order_by = `e-value`, with_ties = TRUE) %>%  # Get rows with the min e-value
  slice(1)


# remove stuff from genome names
gc.put.filt$genome_id <- str_remove_all(gc.put.filt$genome_id , "_(ASM|SRR|TT).*")
ord <- read_tsv('~/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/Zoogloea/250402_cluster_gene_clusters/bigscape/250403_bigscape/output_files/2025-04-03_16-17-52_c0.75/mix/mix_clustering_c0.75.tsv')
gc.put.filt.meta <- gc.put.filt %>% 
  left_join(., gene_functions, by = "Query") %>% 
  left_join(., propnames, by = c('genome_id'='genome')) %>% 
  left_join(., ord, by = 'GBK')
  




gc.put.filt.meta$genome_id <- gsub("GCA_009026175.1","GCF_009026175.1",gc.put.filt.meta$genome_id)
gc.put.filt.meta$genome_id <- gsub('barcode02_Z_caeni','GCA_986281895*', gc.put.filt.meta$genome_id)
gc.put.filt.meta$genome_id <- gsub('barcode03_Z_oleivorans','GCA_986281795*', gc.put.filt.meta$genome_id)
gc.put.filt.meta$genome_id <- gsub('barcode04_Z_resiniphila','GCA_986340685*',gc.put.filt.meta$genome_id)

library(glue)
gc.put.filt.meta <- gc.put.filt.meta %>%
  mutate(y_label = glue("italic('{proposed_name_new}')*' ({genome_id}::{contig_id}::{str_remove_all(region_id, 'region')})'"))

p_geme <- gc.put.filt.meta %>% mutate(stranded = if_else(strand == '+', 1, 0)) %>% 
  ggplot(., aes(xmin = start, xmax = stop, y = y_label ,
                fill = gene_function, forward = stranded)) +
  geom_gene_arrow(alpha = 1,
                  arrowhead_height = grid::unit(2.5, "mm"),
                  arrow_body_height = grid::unit(2, "mm"),
                  arrowhead_width = grid::unit(1, "mm")) + 
  geom_text(
    aes(x=if_else(stop > start,
                  start+100,
                  stop +100), 
        label = Query),
    #nudge_y = 0.03, 
    #angle = 5,
    hjust= 0,
    size= 1)+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=5),
        strip.text.y = element_text(angle=0, vjust=0.5, face = 'bold', size = 8)) +
    facet_grid(rows=vars(Family), scales = 'free', space= 'free' ) +
 # scale_fill_viridis_d(option ="F", na.value = 'grey90')+
  scale_fill_manual( "Gene function",
    values = rev(c( "#34978F","#2F5D5B",
    "#9FB8AD","#A47559","#D9A066","#EAC5AA","#84541E")), na.value = 'grey95')+
  scale_y_discrete(labels = function(x) parse(text = x))

complot <- cowplot::plot_grid(p_tile,
                              p_geme,
                              # labs(tag = 'B'),
                              labels = 'AUTO',
                              label_fontface = "plain",
                              ncol=1,
                              rel_heights = c(0.4,1))
ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS14_epsSMASH_zoogloea_GCfamily_putative.jpeg',
       complot, dpi = 600,
       height = 10,
       width = 11)



  