.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(data.table)
library(treeio)
library(gggenes)


# load hits from HMMsearch of the Zoogloea-polysaccharide gene cluster in AS MAGs as described in the methods section of the manuscript.
polymer <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/AS_MAGs/250224_An/hmmsearch_output/250224_hmmsearch_results.tsv',
                    col_names = T)

pol <- polymer %>% mutate(contig = sub(".*(contig_.*)_.*$", "\\1", Gene_ID))

# data from MiDAS genome paper
approved <- read_tsv('/projects/PHN/MiDAS/04-Analysis/Status-analysis/103-samples/100-add_checkm1/03-MAGs.final-version/000-MAGs-metadata/22277.hqMAGs.IDs', col_names = c('bin')) 
approved$bin <- gsub('.fa','',approved$bin )
metadata <- read_tsv('/projects/PHN/MiDAS/04-Analysis/Status-analysis/103-samples/100-add_checkm1/01-MAGs/000-MAGs-metadata/MiDAS-MAGs-metadata.table.updated-2024-07/MiDAS-27160.hqMAGs.metadata.tsv')%>%
  separate(tax_gtdb.r220, into = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ';')

pattern <- c('wzi','wza','etk','wzc','zooP','zooSA','zooM3', 'zooM2','zooGT1','epsH',
             'asnH','wzy') 

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

p_gene <- df %>% filter(s > -20000, s < 60000) %>% 
  ggplot(., aes(xmin = s, xmax = e, y = paste0(order,'::',family,'::', genus, '::',species,'::', tip_label) ,
                fill = HMM,label = HMM)) +
  geom_gene_arrow(alpha = 1,
                  arrowhead_height = grid::unit(1, "mm"),
                  arrow_body_height = grid::unit(0.8, "mm"),
                  arrowhead_width = grid::unit(0.5, "mm"),
                  linetype = 'blank') +
  #geom_text( aes(x= (s+e)/2, label = HMM ), size=0.6) +
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
  scale_fill_manual('Gene',values=colors,
                    guide = guide_legend(ncol = 1))+
  # scale_fill_viridis_d(option ="F") +
  coord_cartesian(xlim = c(-20000, 60000), clip = "off") +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        plot.margin = margin(5.5, 5, 5.5, 100) ,
        panel.grid.major = element_line(linewidth = 0.2)
        ) 


ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/Fig6_an_gene_cluster_approved_15000bp_genuslab_with_labels.jpeg',
       p_gene, dpi=600,
       height = 10,
       width = 9)









# Genome tree for all MAGs in MiDAS genome paper
t_genome <- '/projects/PHN/MiDAS/200-Rethink/Go-Nagoya/10-others/Sofie/midas.all.hqMAGs-gtdbtk.bac120.decorated.tree'

# Tree file
tree <- read.tree(t_genome)
tree$tip.label <- gsub("'","",tree$tip.label)

############################################################
############ plot full genome tree ########################
p <- ggtree(tree ,size = 0.01, layout = 'circular' ) 

n_distinct(p$data %>% filter(isTip) %>% pull(label))


meta.approved <- metadata %>% filter(bin %in% approved$bin)
meta.approved %>% count(phylum) %>% arrange(-n)

# color top phyla
chip <- meta.approved %>% mutate(
  phyl_col = case_when(phylum == 'p__Pseudomonadota' ~ 'Pseudomonadota', 
                       phylum == 'p__Bacteroidota' ~ 'Bacteroidota',
                       phylum == 'p__Myxococcota' ~'Myxococcota',
                       phylum == 'p__Actinomycetota' ~ 'Actinomycetota',
                       phylum == 'p__Patescibacteria' ~ 'Patescibacteria',
                       phylum == 'p__Chloroflexota' ~ 'Chloroflexota',
                       phylum == 'p__Planctomycetota' ~ 'Planctomycetota',
                       phylum == 'p__Acidobacteriota' ~ 'Acidobacteriota',
                       phylum == 'p__Verrucomicrobiota' ~ 'Verrucomicrobiota',
                       .default = 'Other'
  )) %>% mutate(phyl_col = factor(phyl_col, 
                                  levels = c('Pseudomonadota','Bacteroidota', 'Actinomycetota', 'Patescibacteria','Chloroflexota','Planctomycetota', 'Myxococcota', 
                                             'Verrucomicrobiota',
                                             'Acidobacteriota','Other'))) %>% 
  select(tip_label = bin, phyl_col)

# join with metadata from MiDAS genome paper to get the geographical information for the MAGs
geo_info <- read_tsv('/projects/PHN/MiDAS/04-Analysis/Status-analysis/103-samples/100-add_checkm1/03-MAGs.final-version/000-MAGs-metadata/geo.info/geo.info.tsv', 
                     col_names = c('plant.ID', 'continent', 'country', 'city'))
hits.geo <- hit %>%
  mutate(plant.ID = str_extract(tip_label, "[A-Z]+\\d*-\\d+"),
         hit = 'gene cluster') %>% 
  left_join(., geo_info)
p.genome <- p %<+% hits.geo
#p.genome <- p %<+% rep_hit_names


species_count <- com_pattern_operon_meta %>% 
  group_by(family) %>%
  summarise(unique_genomes = n_distinct(tip_label)) %>%
  mutate(label = paste0(family, " (", unique_genomes, ")"))

# Create a named vector for renaming legend labels
family_labels <- setNames(species_count$label, species_count$family)


fam_colors <- c( "#A96C4A", "#d08c7e", "#f6b26b", "#F9EACF",  "#e5c06f","#ccd5ae","#a3b18a", "#A5D1C0", '#4BA9A3','transparent')
# Plot with updated legend labels
ps <- p.genome + geom_treescale() +
  geom_tippoint(mapping = aes(color = family), size = 2, na.rm = TRUE) +
  scale_color_manual(name = "Family (Genome count)", labels = family_labels, 
                       values = fam_colors, na.value = 'transparent')


ps <- p.genome + geom_treescale() +
  geom_tippoint(mapping = aes(color = hit), size = 2, na.rm = TRUE) +
  scale_color_manual(name = "Zoogloea gene cluster", 
                       values = '#a03a3a', na.value = 'transparent')


# phyla_colors = c('#F5F5F5','#452929','#84541E','#9c7605', '#6a170e','#bf4e4e','#D18E35', '#E0C37E' ,'#80CDC1','#34978F')
# phyla_colors <- rev(phyla_colors)


# Plot with updated legend labels
ps <- p.genome + geom_treescale() +
  geom_tippoint(mapping = aes(fill = order, color = order), size = 2.2, na.rm = TRUE, shape=21) +
  scale_fill_manual(name = "Order", 
                     values = c("#F9EACF","#c8a96a", "#c65e03"),na.value = "transparent",
                    na.translate = FALSE) +
  scale_color_manual(name = "", 
                    values = c("black","black", "black"),na.value = "transparent",
                    na.translate = FALSE,
                    guide= 'none')


phyla_colors <- c("#5B8FA2", "#A0AFC0", "#C1C9D6", "#E1DDE6", "#9BA8A1", "#B6B8AC", "#D9D4CF", "#A3A1B0", "#768793",'transparent')

#phyla_colors <- c("#5B8FB8", "#A174A0", "#D4A373", "#9DBFAF", "#728F81", "#B38B6D", "#8493A7", "#A6B6C9", "#779CAB","#F5F5F5")

psave <- ps + new_scale_fill() + geom_fruit(
  data = chip,
  geom = geom_tile,
  mapping = aes(
    y=tip_label,
    group=label,
    fill=phyl_col
  ),
  alpha = 0.8,
  width=0.2,
  offset= -0.08 #before -0.03
) +
  scale_fill_manual('Phylum',
                    values = phyla_colors)+
  labs(tag='A')+
  theme(plot.tag.position = c(0.06, 0.88),
        plot.tag = element_text(size = 20),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA))



tree_data <- ps$data %>% filter(isTip) %>% 
  select(label, x, y)

phyla_labels <- chip %>%
  group_by(phyl_col) %>%
  slice(n()) %>%      # or slice(1)
  ungroup() %>% 
  left_join(tree_data, by = c("tip_label" = "label")) %>%
  filter(!phyl_col %in% 'Other') %>% 
  mutate(
    angle = atan2(y, x),
    angle_deg = angle * 180 / pi,
    hjust = ifelse(angle_deg > 90 | angle_deg < -90, 1, 0),
    angle_deg = ifelse(angle_deg > 90 | angle_deg < -90,
                       angle_deg + 180,
                       angle_deg),
  ) %>% 
  select(-x,-y) %>% 
  mutate(hjust = case_when(phyl_col == 'Pseudomonadota' ~ -0.03, 
                           phyl_col == 'Bacteroidota' ~ 0.95,
                           phyl_col == 'Myxococcota' ~ 0.85,
                           phyl_col == 'Actinomycetota' ~ -0.04,
                           phyl_col == 'Patescibacteria' ~ 0.3,
                           phyl_col == 'Chloroflexota' ~ 0.1,
                           phyl_col == 'Planctomycetota' ~ 1,
                           phyl_col == 'Acidobacteriota' ~ 0.95,
                           phyl_col == 'Verrucomicrobiota' ~ 1))



psave <- ps +
  new_scale_fill() +
  geom_fruit(
    data = chip,
    geom = geom_tile,
    mapping = aes(
      y = tip_label,
      group = label,
      fill = phyl_col
    ),
    width = 0.2,
    offset = -0.08,
    alpha = 0.8
  ) +
  scale_fill_manual(values = phyla_colors, guide = 'none')+
  geom_fruit(
    data = phyla_labels,
    geom = geom_text,
    mapping = aes(y = tip_label, label = phyl_col, hjust=hjust),
    offset = 0.1, angle=0,
     size=3
  )+
 # labs(tag='A')+
  theme(plot.tag.position = c(0.06, 0.88),
        plot.tag = element_text(size = 20),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA))



psave

####################### sankey plot

library(ggsankey)


# Count occurrences for each column
taxonomic_cols <- c('class', 'order', 'family', 'genus', 'species')
df_with_counts <- hits.geo %>%
  mutate(across(all_of(taxonomic_cols), ~ paste0(.x, " (", ave(seq_along(.x), .x, FUN = length))))
df_with_counts_meta <-meta.approved.hits.geo %>%
  mutate(across(all_of(taxonomic_cols), 
                list(count = ~ ave(seq_along(.x), .x, FUN = length)))) %>% 
  select(bin,class_count,order_count,family_count,genus_count,species_count)

df_with_counts_tot <- df_with_counts %>% left_join(., df_with_counts_meta, by = c('tip_label'='bin')) %>% 
  mutate(across(all_of(taxonomic_cols), ~ paste0(.x, "/", get(paste0(cur_column(), "_count")), ")")))

# Ensure df.san.hits is sorted so that higher-frequency nodes appear first
level_nodes <- df_with_counts_tot %>%
  make_long("class", "order", "family", "genus") %>%
  group_by(x, node) %>%
  summarise(count = n()) %>%
  arrange(x, desc(count)) 


levl <- df_with_counts_tot %>%
  arrange(class, order, family, genus) %>%
  distinct(class, order, family, genus) %>%
  unlist(use.names = FALSE) %>% unique()



df.san.hits <- df_with_counts_tot %>%
  filter(!genus %in% 'g__ (9/2957)') %>% 
   make_long(  "class"    ,      "order"    ,     "family"     ,    "genus") %>% 
  mutate(node = factor(node, levels = levl))
  #make_long(  "class"    ,      "order"    ,     "family"     ,   'continent')

san_colors=c('#dad7cd','#E0C37E',
             '#e9edc9','#ccd5ae','#adc178','#a3b18a' )


san_colors=c('#f0efeb',#'#f0efeb','#f0efeb','#f0efeb',
             "#b23a48","#e62989", "#ff8e72",
             
             
             "#A96C4A", "#d08c7e", "#f6b26b", "#F9EACF",  "#e5c06f","#ccd5ae","#a3b18a", "#A5D1C0", '#4BA9A3',
             
             
             "#A96C4A", "#A96C4A", "#A96C4A", "#A96C4A", "#A96C4A",   
             "#d08c7e", 
             "#f6b26b",  
             "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF",
             "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF",
             
             "#e5c06f",
             "#ccd5ae",
             "#a3b18a","#a3b18a","#a3b18a",
             '#A5D1C0',
           '#4BA9A3' )


san_colors=c('#f0efeb',#'#f0efeb','#f0efeb','#f0efeb',
             "#F9EACF","#c8a96a", "#c65e03",
             
             
             "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF",  "#F9EACF","#F9EACF","#c8a96a", "#c65e03", '#c65e03',
             
             
             "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF",   
             "#F9EACF", 
             "#F9EACF",  
             "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF",
             "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF", "#F9EACF",
             
             "#c8a96a",
             "#c8a96a",
             "#c8a96a","#c8a96a","#c8a96a",
             '#c65e03',
             '#c65e03' )




psan <- ggplot(df.san.hits, aes(x = x, 
                                next_x = next_x, 
                                node = node, 
                                next_node = next_node,
                                fill = node,
                                label = node)) +
  geom_sankey(type = 'alluvial',
              smooth = 6,
              width = 0, 
              space= 7) +
  geom_sankey_text( type = 'alluvial',hjust=0, vjust=0.25,
                     space= 7, color = 'black',
                    size=3)+
  theme_sankey(base_size = 18)+
  #labs(tag='B')+
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.tag.position = c(0.1, 0.65),
        plot.tag = element_text(size = 15),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = margin(t = 100, r = 10, b = 10, l = -50)) +
  scale_fill_manual(values =san_colors) +
  scale_color_manual(values =san_colors)

psan

plot_5 <- cowplot::ggdraw() +
  cowplot::draw_plot(psan) +                       # main plot
  cowplot::draw_plot(
    psave,
    x = 0.02,      # left position (0 = left edge)
    y = 0.42,      # bottom position
    width = 0.62,
    height = 0.62
  )+ # A label (top-left of main panel)
  cowplot::draw_label(
    "A",
    x = 0.01,
    y = 0.98,
    size = 18,
    fontface = "plain"
  ) +
  
  # B label (top-left of inset OR second panel position)
  cowplot::draw_label(
    "B",
    x = 0.01,
    y = 0.43,
    size = 18,
    fontface = "plain"
  )



ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/Fig5_hits_AS_MAGS_approved_15000bp.jpeg',
       plot_5, dpi = 800,
       height = 8,
       width = 13)



