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

pattern <- c('wza','etk','wzc',
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
met.tax <- metadata %>% select(bin, genus) %>% filter(bin %in% approved$bin) %>% unique() %>%  count(genus) %>% arrange(-n) %>% rename(MAG.nr = n)

n_distinct(com_pattern$Genome)


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

n_distinct(com_filtered$Genome)


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
n_distinct(com_plot$tip_label)



# order of genes
prot_vector <- c(
  "asnB", "zooTRP", "wzxC", "zooGT8", "zooGT7", "zooGH", "zooM1", 
  "zooGT6", "zooGT5", "epsH", "wzy", "capK", "zooGT4", "zooGT3", 
  "asnH", "zooGT2", "zooGT1", "zooM2", "zooM3", "zooSA", "zooP", 
  "wzc", "etk", "wza", "wzi", "lolD"
)

# factorise
com_plot$HMM <- factor(com_plot$HMM, 
                       levels = rev(prot_vector))
com_plot %>% ggplot(aes(x=HMM,y= paste0(tip_label),
                        fill=HMM)) +
  geom_tile() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=2))


# In the same gene cluster only
com_pattern_operon <- com_plot %>%
  group_by(tip_label, contig,operon) %>%
  filter(all(pattern %in% HMM)) %>%
  ungroup()

n_distinct(com_pattern_operon$tip_label)




com_pattern_operon_meta <- com_pattern_operon %>% left_join(., metadata, by = c('tip_label' = 'bin')) 

hit <- com_pattern_operon_meta%>% select(tip_label, 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species',tax_midas_v5.3) %>% unique()
hit.tax <- hit %>% select(tip_label, genus) %>% unique() %>%  count(genus) %>% arrange(-n) %>% rename(MAG.with.cluster=n)
tot.tax <- c.pat.met.tax %>% left_join(., hit.tax, by = 'genus') %>% 
  left_join(., met.tax, by = 'genus') 


## plot of overview of nr. MAGs, nr. MAGs with core genes, nr. MAGs with core genes in cluster
tot.tax %>% 
  filter(!genus %in% 'g__',
         MAG.nr >2) %>% 
  mutate(genus = fct_reorder(genus, MAG.nr, .fun = desc)) %>%  # Order by decreasing percentage
  pivot_longer(cols = -genus, names_to = 'group', values_to = 'number') %>% 
  mutate(group = factor(group, levels= c('MAG.nr','MAG.with.core.genes','MAG.with.cluster'))) %>% 
  ggplot(aes(x = genus, fill = group, y = number)) +
  geom_col(position = position_dodge(), width = 0.9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank())

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



###### try here
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


bar <- ggplot(dat, aes(x = genus)) +
  geom_col(aes(y = number, fill = group),position = 'stack') +  # Stacks bars and scales to proportions
  #scale_y_continuous(labels = scales::percent_format()) +  # Show proportions as %
  scale_fill_manual(values = c("#FFF5EB",  "#F5C9B0", "#EC9E75")) +  # Custom colors
  theme_minimal(base_size = 10) +
  geom_label(aes(label = MAG.nr, y = 50), size=2)+
  labs(x = "Genus", y = "Proportion", fill = "Group") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))





com_pattern_operon_meta %>% 
  group_by(tip_label,operon) %>%
  filter(!any(HMM == "wzi")) %>% 
  ggplot(aes(x=HMM,y= paste0(family, genus, species,tip_label),
                                       fill=HMM)) +
  geom_tile() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=3))




com_pattern_operon_meta %>% count(HMM)


###################################### start her #####################################
######################################################################################
# for gene cluster plot
darrow <- com_pattern_operon_meta%>% 
  group_by(tip_label,operon) %>%
  filter(!any(HMM == "wzi")) %>% 
  filter(HMM != "wza" | `E-value` == min(`E-value`[HMM == "wza"])) %>% 
  mutate(wza_strand = ifelse(any(HMM == "wza" & Strandedness == 1), 1, 0)) %>% 
  mutate(start_relative = if_else(Strandedness == 1, 
                                  Start_Position - Start_Position[HMM == "wza"],
                                  End_Position - End_Position[HMM == "wza"]),
         end_relative = if_else(Strandedness == 1, 
                                End_Position- Start_Position[HMM == "wza"],
                                Start_Position - End_Position[HMM == "wza"])) %>% 
  mutate(s = if_else(wza_strand == 1, start_relative, -start_relative),
         e = if_else(wza_strand == 1, end_relative, -end_relative))


# Reorder factor so genomes appear in desired order (e.g., by genus)
df <- darrow %>% ungroup() %>% 
  arrange(order, family,genus, species, tip_label,operon) %>%
  mutate(tip_label = factor(tip_label, levels = unique(tip_label)),
         tip_y = as.numeric(tip_label))

xpos <- -45000
# Compute segment positions for each genus
genus_segments <- df %>%
  group_by(genus) %>%
  summarise(
    y_min = min(tip_y) -0.3,
    y_max = max(tip_y) +0.3, # line position
    label_x = xpos-500, # label position
    label_y = (y_min + y_max) / 2,
    .groups = "drop"
  ) %>% filter(!genus %in% c('g__')) %>% 
  mutate(genus = str_remove_all(genus, 'g__'))

genus_segments$genus <- sapply(genus_segments$genus, function(x) {
  if (grepl("^Ca_", x)) {
    rest <- sub("^Ca_", "", x)
    paste0("italic('Ca.')~'", rest, "'")
  } else if (!grepl("[0-9]", x)) {
    paste0("italic('", x, "')")
  } else {
    paste0("'", x, "'")
  }
})


colors <- c(
  "#4BA9A3", "#6CB5AB", "#8AC4B5", "#A5D1C0", "#BEDFCE",  # Soft teals to seafoam
  "#C9D9C6", "#D9D1B8", "#E6D8B8", "#F2E0C3", "#F9EACF",  # Sand and cream
  "#FCEFCF", "#F9E6B0", "#F2D58F", "#EBC27A", "#E1AE67",  # Pale gold to warm amber
  "#D29B61", "#C78B5B", "#B97A53", "#A96C4A", "#986042",  # Warm tans to brown
  "#87613F", "#765139", "#654430", "#553A2A", "#483024", "#3A261D" # Driftwood fade
)

p_gene <- darrow %>% 
  ggplot(., aes(xmin = s, xmax = e, y = paste0(order,'::',family,'::', genus, '::',species,'::', tip_label) ,
                fill = HMM,label = HMM)) +
  geom_gene_arrow(alpha = 0.8,
                  arrowhead_height = grid::unit(1, "mm"),
                  arrow_body_height = grid::unit(0.8, "mm"),
                  arrowhead_width = grid::unit(0.5, "mm"),
                  linetype = 'blank') + 
  theme_minimal() +
  # Genus grouping lines
  geom_segment(data = genus_segments ,
               aes(x = xpos, xend = xpos, y = y_min, yend = y_max),
               inherit.aes = FALSE, color = "grey30", size = 0.6) +
  # Genus labels
  geom_text(data = genus_segments ,
            aes(x = label_x, y = label_y, label = genus),
            inherit.aes = FALSE, hjust = 1, size = 2.1, 
            parse=T)+
  scale_fill_manual('Gene', values = colors,
                    guide = guide_legend(ncol = 1))+
  coord_cartesian(xlim = c(-40000, 70000), clip = "off") +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        plot.margin = margin(5.5, 5, 5.5, 100) ,
        panel.grid.major = element_line(linewidth = 0.2)
  ) 


# 
ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS11_an_gene_cluster_approved_withoutwzi_15000bp.jpeg',
       p_gene, dpi=800,
       height = 14,
       width = 10)



