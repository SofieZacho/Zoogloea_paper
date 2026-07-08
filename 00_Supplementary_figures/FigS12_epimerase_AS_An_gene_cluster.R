

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
## epimerase
wecB <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/AS_MAGs/250224_An/epimerase/250304_hmmsearch_wecB/250304_hmmsearch_wecB_results.tsv') %>% 
  mutate(contig = sub(".*(contig_.*)_.*$", "\\1", Gene_ID))

polall <- rbind(pol, wecB)

approved <- read_tsv('/projects/PHN/MiDAS/04-Analysis/Status-analysis/103-samples/100-add_checkm1/03-MAGs.final-version/000-MAGs-metadata/22277.hqMAGs.IDs', col_names = F) 

approved$X1 <- gsub('.fa','',approved$X1 )
approved <- approved%>% 
  rename(bin=`X1`)

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
com_pattern <- polall %>%
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
  "asnH", "zooGT2", "zooGT1", "zooM2", "zooM3","wecB_TIGR00236", "zooSA", "zooP", 
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



com_pattern_operon_meta %>% count(HMM)


###################################### start her #####################################
######################################################################################
# for gene cluster plot
darrow <- com_pattern_operon_meta%>% 
  group_by(tip_label,operon) %>% 
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

xpos <- -25000
# Compute segment positions for each genus
genus_segments <- df %>%
  group_by(genus) %>%
  summarise(
    y_min = min(tip_y) -0.3,
    y_max = max(tip_y) +0.3, # line position
    label_x = xpos-500, # label position
    label_y = (y_min + y_max) / 2,
    .groups = "drop"
  )%>% filter(!genus %in% c('g__')) %>% 
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

p_gene <- darrow%>% filter(s > -20000, s < 60000) %>%
  ggplot(., aes(xmin = s, xmax = e, y = paste0(order,'::',family,'::', genus, '::',species,'::', tip_label) ,
                fill = if_else(HMM %in% 'wecB_TIGR00236', 'wecB', 'other'))) +
  geom_gene_arrow(alpha = 0.8,
                  arrowhead_height = grid::unit(1, "mm"),
                  arrow_body_height = grid::unit(0.8, "mm"),
                  arrowhead_width = grid::unit(0.5, "mm"),
                  linetype = 'blank') +
  theme_minimal() +
  geom_segment(data = genus_segments,
               aes(x = xpos, xend = xpos, y = y_min, yend = y_max),
               inherit.aes = FALSE, color = "grey30", size = 0.6) +
  # Genus labels
  geom_text(data = genus_segments,
            aes(x = label_x, y = label_y, label = genus),
            inherit.aes = FALSE, hjust = 1, size = 2.2,
            parse=T) +
  coord_cartesian(xlim = c(-20000, 60000), clip = "off") +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        plot.margin = margin(5.5, 5, 5.5, 100) ,
        panel.grid.major = element_line(linewidth = 0.2)
  ) +
  scale_fill_manual('HMM',values= c('grey80','hotpink4'))


ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS12_an_gene_cluster_wecB_color_approved_15000bp.jpeg',
       p_gene, dpi = 600,
       height = 10,
       width = 9)




