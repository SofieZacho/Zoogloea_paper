.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")
library(tidyverse)
library(ggtree)



# data paths
rethink <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/epsSMAH/260724_AS_MAGs_epsSMASH/region_counts.tsv', skip = 1) # data from epsSMASH run on MiDAS MAGs with Zoogloea polysaccharide gene cluster detected
rethink$record <- gsub('.fa','',rethink$record)


meta <- read_tsv('/projects/PHN/MiDAS/04-Analysis/Status-analysis/103-samples/100-add_checkm1/03-MAGs.final-version/000-MAGs-metadata/MiDAS-22277.hqMAGs.metadata.tsv')%>%
  separate(tax_gtdb.r220, into = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ';')

# # remove columns with no hits (meaning column sum is 0) - only for numeric columns
rethink.hit.genomes.filt <- rethink %>%
  select(where(~ !is.numeric(.) || (is.numeric(.) && sum(.) > 0))) %>%
  select(#-contains("putative"),
         -description)


lv <-names(colSums(rethink.hit.genomes.filt %>% select(where (~ is.numeric(.) ))) %>% sort())

rethink.hit.genomes.filt.long <- rethink.hit.genomes.filt %>%
  pivot_longer(., cols = -record) %>%
  mutate(name = factor(name, levels = rev(lv)))



rethink.hit.genomes.filt.long.meta <- rethink.hit.genomes.filt.long  %>%
  left_join(., meta, by = c('record'='bin'))





pp <-rethink.hit.genomes.filt.long.meta %>% filter(!name %in% 'total_count') %>% 
  ggplot() +
  geom_tile(aes(x= paste0( genus, species,record), y = name, fill = factor(value)),
            size= 0.02, color = 'grey20')+
  scale_fill_manual('Number',values = c('grey95','#EAC5AA', '#d98c6b','hotpink4')) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(size=3.5, angle = 90, hjust=1, vjust= 0.5),
    axis.text.y=element_text(size=9, color = 'black'),
    axis.title = element_blank(),
    strip.text.x = element_text(angle=90, hjust=0, face = 'bold', size = 8),
    strip.clip = 'off'
  )+
  facet_grid(cols=vars( family), scales = 'free', space= 'free' )


ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS15_hit_genomes_epshits.jpeg',
       pp, dpi = 600,
       width = 13,
       height = 5.5)
