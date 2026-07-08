.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")
library(tidyverse)
library(ggtree)



# data paths
rethink <- read_tsv('/home/bio.aau.dk/zs85xk/projects/epsSMASH/epsSMASH_analyses/REThiNk_catalogue/data/epsSMASH_results/region_counts.tsv', skip = 1) # data from doi: https://doi.org/10.64898/2025.12.21.693542 (Daugberg et al., 2025)
rethink$record <- gsub('.gbff','',rethink$record)

# genomes with Zoogloea polysaccharide gene cluster
hit.genomes <- read_tsv('~/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/AS_MAGs/250224_An/hit_genomes.txt',
                        col_names = 'genome')


# filter to relevant genomes
rethink.hit.genomes <- rethink %>% 
  filter(genome_id %in% hit.genomes$genome)
n_distinct(rethink.hit.genomes$genome_id)

meta <- read_tsv('/projects/PHN/MiDAS/04-Analysis/Status-analysis/103-samples/100-add_checkm1/03-MAGs.final-version/000-MAGs-metadata/MiDAS-22277.hqMAGs.metadata.tsv')%>%
  separate(tax_gtdb.r220, into = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ';')

# remove columns with no hits (meaning column sum is 0) - only for numeric columns
rethink.hit.genomes.filt <- rethink.hit.genomes %>% 
  select(where(~ !is.numeric(.) || (is.numeric(.) && sum(.) > 0))) %>% 
  select(#-contains("putative"),
         -description)


lv <-names(colSums(rethink.hit.genomes.filt %>% select(where (~ is.numeric(.) ))) %>% sort())

rethink.hit.genomes.filt.long <- rethink.hit.genomes.filt %>% 
  pivot_longer(., cols = -record) %>% 
  mutate(name = factor(name, levels = rev(lv)))

rethink.hit.genomes.filt.long.meta <- rethink.hit.genomes.filt.long %>% 
  left_join(., meta, by = c('record'='bin'))

pp <-rethink.hit.genomes.filt.long.meta %>% filter(!name %in% 'total_count') %>% 
ggplot() +
  geom_tile(aes(x= paste0( genus, species,record), y = name, fill = factor(value)),
            size= 0.02, color = 'grey20')+
  scale_fill_manual('Number',values = c('grey95','#EAC5AA', '#d98c6b')) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(size=3.5, angle = 90, hjust=1, vjust= 0.5),
    axis.title = element_blank(),
    strip.text.x = element_text(angle=90, hjust=0, face = 'bold', size = 8),
    strip.clip = 'off'
  )+
  facet_grid(cols=vars( family), scales = 'free', space= 'free' )


ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS15_hit_genomes_epshits.jpeg',
       pp,
       width = 13,
       height = 5.5)
