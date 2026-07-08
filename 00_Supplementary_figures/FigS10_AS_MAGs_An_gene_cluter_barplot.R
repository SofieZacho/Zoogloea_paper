.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(data.table)
library(treeio)
library(gggenes)


# load hits from 

polymer <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/AS_MAGs/250224_An/hmmsearch_output/250224_hmmsearch_results.tsv',
                    col_names = T)

pol <- polymer %>% mutate(contig = sub(".*(contig_.*)_.*$", "\\1", Gene_ID))

approved <- read_tsv('/projects/PHN/MiDAS/04-Analysis/Status-analysis/103-samples/100-add_checkm1/03-MAGs.final-version/000-MAGs-metadata/22277.hqMAGs.IDs', col_names = F) %>% 
  rename(bin=X1)
approved$bin <- gsub('.fa','',approved$bin )

metadata <- read_tsv('/projects/PHN/MiDAS/04-Analysis/Status-analysis/103-samples/100-add_checkm1/01-MAGs/000-MAGs-metadata/MiDAS-MAGs-metadata.table.updated-2024-07/MiDAS-27160.hqMAGs.metadata.tsv')%>%
  separate(tax_gtdb.r220, into = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ';')

pattern <- c('wzi','wza','etk','wzc',
             'zooP',
             'zooSA',
             'zooM3',
             'zooM2','zooGT1','epsH',
             #'zooGT2',
             'asnH','wzy') # found in many of the hits - really conserved

# In the same genome
com_pattern <- pol %>%
  filter(Genome %in% approved$bin) %>% 
  group_by(Genome) %>%
  filter(all(pattern %in% HMM)) %>%
  ungroup()
c.pat.met <- com_pattern %>% left_join(., metadata, by = c('Genome'='bin'))
c.pat.met.tax <- c.pat.met %>% select(Genome, genus) %>% unique() %>%  count(genus) %>% arrange(-n) %>% rename(MAG.with.core.genes =n)




#### sort stuff
com_sorted <- com_pattern %>% 
  group_by(Genome, contig) %>% 
  arrange(Start_Position, .by_group = T) %>%
  group_by(Gene_ID) %>%
  filter(`E-value` == min(`E-value`))


n <- 15000
com_filtered <- com_sorted  %>% 
  group_by(Genome, contig) %>% 
  mutate(prior =lag(End_Position), after = lead(Start_Position)) %>% 
  mutate(dif_bef = Start_Position - lag(End_Position),
         dif_aft = lead(Start_Position) - End_Position) %>% 
  filter(!if_all(c(dif_bef, dif_aft), is.na)) %>% #removes hits from contigs with only 1 hit
  filter(dif_bef < n | dif_aft < n | is.na(dif_bef) & dif_aft < n | is.na(dif_aft) & dif_bef < n)


com_plot<- com_filtered %>% 
  ungroup() %>% 
  group_by(Genome) %>% 
  mutate(
    operon_start = ifelse(dif_bef > n | is.na(dif_bef), "start", NA),
    operon = cumsum(!is.na(operon_start)),  # Count occurrences of "start" within each group
  ) %>%
  group_by(Genome, operon) %>%
  filter(n_distinct(HMM) >1) %>%
  rename(tip_label = Genome)

# In the same gene cluster only
com_pattern_operon <- com_plot %>%
  group_by(tip_label, contig,operon) %>%
  filter(all(pattern %in% HMM)) %>%
  ungroup()

# join with metadata
com_pattern_operon_meta <- com_pattern_operon %>% left_join(., metadata, by = c('tip_label' = 'bin')) 

hit <- com_pattern_operon_meta%>% select(tip_label, 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species',tax_midas_v5.3) %>% unique()
hit.tax <- hit %>% select(tip_label, genus) %>% unique() %>%  count(genus) %>% arrange(-n) %>% rename(MAG.with.cluster=n)
tot.tax <- c.pat.met.tax %>% left_join(., hit.tax, by = 'genus') %>% 
  left_join(., met.tax, by = 'genus') 


# other plot
dat <- tot.tax %>% 
  filter(!genus %in% 'g__',
         MAG.nr >2
  ) %>% 
  mutate(genus = fct_reorder(genus, MAG.nr),
         p.MAG.nr =100* MAG.nr/MAG.nr,
         p.MAG.with.core.genes = 100*MAG.with.core.genes/MAG.nr,
         p.MAG.with.cluster = 100*MAG.with.cluster/MAG.nr,
         p.MAG.with.cluster= replace_na(p.MAG.with.cluster, 0),
         p.MAG.with.core.genes= replace_na(p.MAG.with.core.genes, 0)) %>% 
  mutate(p.MAG.nr = p.MAG.nr-p.MAG.with.core.genes,
         p.MAG.with.core.genes = p.MAG.with.core.genes-p.MAG.with.cluster) %>% # Order by decreasing percentage
  pivot_longer(cols = c(p.MAG.nr,p.MAG.with.core.genes,p.MAG.with.cluster), names_to = 'group', values_to = 'number' ) %>% 
  mutate(group = factor(group, levels= c('p.MAG.nr','p.MAG.with.core.genes','p.MAG.with.cluster'))) 



# Order groups to ensure proper stacking
# Compute sorting order:
genus_order <- dat %>%
  filter(group %in% c("p.MAG.with.cluster", "p.MAG.with.core.genes")) %>%
  group_by(genus) %>%
  summarize(
    cluster_value = sum(number[group == "p.MAG.with.cluster"]),  # Value for cluster
    core_value = sum(number[group == "p.MAG.with.core.genes"])   # Value for core genes
  ) %>%
  arrange(desc(cluster_value), desc(core_value)) %>%
  pull(genus)

# Apply ordering to genus
dat$genus <- factor(dat$genus, levels = genus_order)


bar <- dat %>% filter(group %in% c('p.MAG.with.core.genes','p.MAG.with.cluster')) %>% 
  mutate(group = if_else(group == 'p.MAG.with.core.genes', 'core.genes.detected','core.genes.in.cluster')) %>% 
  ggplot(., aes(x = genus)) +
  geom_col(aes(y = number, fill = group),position = 'stack') +  # Stacks bars and scales to proportions
  #scale_y_continuous(labels = scales::percent_format()) +  # Show proportions as %
  scale_fill_manual(values = c( "#F5C9B0", "#EC9E75")) +  # Custom colors
  theme_minimal(base_size = 10) +
  geom_label(aes(label = MAG.nr, y = 50), size=2)+
  labs(x = "Genus", y = "Percentage of genomes in genus [%]", fill = "", tag = "A") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), axis.title.x= element_blank(),
        legend.position = 'top') +
  scale_y_continuous(expand = c(0,0))



####################################################################################################
####### without wzi as core gene ###################################################################
####################################################################################################
pattern <- c('wza','etk','wzc',
             'zooP',
             'zooSA',
             'zooM3',
             'zooM2','zooGT1','epsH',
             'asnH','wzy') # found in many of the hits - really conserved

# In the same genome
com_pattern <- pol %>%
  filter(Genome %in% approved$bin) %>% 
  group_by(Genome) %>%
  filter(all(pattern %in% HMM)) %>%
  ungroup()
c.pat.met <- com_pattern %>% left_join(., metadata, by = c('Genome'='bin'))
c.pat.met.tax <- c.pat.met %>% select(Genome, genus) %>% unique() %>%  count(genus) %>% arrange(-n) %>% rename(MAG.with.core.genes =n)




#### sort stuff
com_sorted <- com_pattern %>% 
  group_by(Genome, contig) %>% 
  arrange(Start_Position, .by_group = T) %>%
  group_by(Gene_ID) %>%
  filter(`E-value` == min(`E-value`))


n <- 15000
com_filtered <- com_sorted  %>% 
  group_by(Genome, contig) %>% 
  mutate(prior =lag(End_Position), after = lead(Start_Position)) %>% 
  mutate(dif_bef = Start_Position - lag(End_Position),
         dif_aft = lead(Start_Position) - End_Position) %>% 
  filter(!if_all(c(dif_bef, dif_aft), is.na)) %>% #removes hits from contigs with only 1 hit
  filter(dif_bef < n | dif_aft < n | is.na(dif_bef) & dif_aft < n | is.na(dif_aft) & dif_bef < n)


com_plot<- com_filtered %>% 
  ungroup() %>% 
  group_by(Genome) %>% 
  mutate(
    operon_start = ifelse(dif_bef > n | is.na(dif_bef), "start", NA),
    operon = cumsum(!is.na(operon_start)),  # Count occurrences of "start" within each group
  ) %>%
  group_by(Genome, operon) %>%
  filter(n_distinct(HMM) >1) %>%
  rename(tip_label = Genome)

# In the same gene cluster only
com_pattern_operon <- com_plot %>%
  group_by(tip_label, contig,operon) %>%
  filter(all(pattern %in% HMM)) %>%
  ungroup()

# join with metadata
com_pattern_operon_meta <- com_pattern_operon %>% left_join(., metadata, by = c('tip_label' = 'bin')) 

hit <- com_pattern_operon_meta%>% select(tip_label, 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species',tax_midas_v5.3) %>% unique()
hit.tax <- hit %>% select(tip_label, genus) %>% unique() %>%  count(genus) %>% arrange(-n) %>% rename(MAG.with.cluster=n)
tot.tax <- c.pat.met.tax %>% left_join(., hit.tax, by = 'genus') %>% 
  left_join(., met.tax, by = 'genus') 


# other plot
dat <- tot.tax %>% 
  filter(!genus %in% 'g__',
         MAG.nr >2
  ) %>% 
  mutate(genus = fct_reorder(genus, MAG.nr),
         p.MAG.nr =100* MAG.nr/MAG.nr,
         p.MAG.with.core.genes = 100*MAG.with.core.genes/MAG.nr,
         p.MAG.with.cluster = 100*MAG.with.cluster/MAG.nr,
         p.MAG.with.cluster= replace_na(p.MAG.with.cluster, 0),
         p.MAG.with.core.genes= replace_na(p.MAG.with.core.genes, 0)) %>% 
  mutate(p.MAG.nr = p.MAG.nr-p.MAG.with.core.genes,
         p.MAG.with.core.genes = p.MAG.with.core.genes-p.MAG.with.cluster) %>% # Order by decreasing percentage
  pivot_longer(cols = c(p.MAG.nr,p.MAG.with.core.genes,p.MAG.with.cluster), names_to = 'group', values_to = 'number' ) %>% 
  mutate(group = factor(group, levels= c('p.MAG.nr','p.MAG.with.core.genes','p.MAG.with.cluster'))) 



# Order groups to ensure proper stacking
# Compute sorting order:
genus_order <- dat %>%
  filter(group %in% c("p.MAG.with.cluster", "p.MAG.with.core.genes")) %>%
  group_by(genus) %>%
  summarize(
    cluster_value = sum(number[group == "p.MAG.with.cluster"]),  # Value for cluster
    core_value = sum(number[group == "p.MAG.with.core.genes"])   # Value for core genes
  ) %>%
  arrange(desc(cluster_value), desc(core_value)) %>%
  pull(genus)

# Apply ordering to genus
dat$genus <- factor(dat$genus, levels = genus_order)


bar.without.wzi <- dat %>% filter(group %in% c('p.MAG.with.core.genes','p.MAG.with.cluster')) %>% 
  mutate(group = if_else(group == 'p.MAG.with.core.genes', 'core.genes.detected','core.genes.in.cluster')) %>% 
  ggplot(., aes(x = genus)) +
  geom_col(aes(y = number, fill = group),position = 'stack') +  # Stacks bars and scales to proportions
  #scale_y_continuous(labels = scales::percent_format()) +  # Show proportions as %
  scale_fill_manual(values = c( "#F5C9B0", "#EC9E75")) +  # Custom colors
  theme_minimal(base_size = 10) +
  geom_label(aes(label = MAG.nr, y = 50), size=2)+
  labs(x = "Genus", y = "Percentage of genomes in genus [%]", fill = "", tag = "B") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), axis.title.x= element_blank(),
        legend.position = 'none') +
  scale_y_continuous(expand = c(0,0))

psave <- cowplot::plot_grid(bar,bar.without.wzi,
                            ncol=1, rel_heights = c(1,0.9))

psave

ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS10_an_gene_cluster_bar_withoutwzi_approved_15000bp.jpeg',
       psave, dpi=600,
       height = 8,
       width = 15)





