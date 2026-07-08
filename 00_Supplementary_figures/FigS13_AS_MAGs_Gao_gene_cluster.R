.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(data.table)
library(treeio)
library(gggenes)


# load hits from 

polymer <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/AS_MAGs/250305_Gao/hmmsearch_output/250305_hmmsearch_results_gao.tsv',
                    col_names = T)

pol <- polymer %>% mutate(contig = sub(".*(contig_.*)_.*$", "\\1", Gene_ID))

approved <- read_tsv('/projects/PHN/MiDAS/04-Analysis/Status-analysis/103-samples/100-add_checkm1/03-MAGs.final-version/000-MAGs-metadata/22277.hqMAGs.IDs', col_names = F) %>% 
  rename(bin=X1)
approved$bin <- gsub('.fa','',approved$bin )

metadata <- read_tsv('/projects/PHN/MiDAS/04-Analysis/Status-analysis/103-samples/100-add_checkm1/01-MAGs/000-MAGs-metadata/MiDAS-MAGs-metadata.table.updated-2024-07/MiDAS-27160.hqMAGs.metadata.tsv')%>%
  separate(tax_gtdb.r220, into = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ';')

pattern <- c("epsB2", "prsK", "prsR", "prsT", "ugd") # found in many of the hits - really conserved

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
prot_vector <- c("epsB2", "prsK", "prsR", "prsT", "ugd","mltE", "degQ2"
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
  scale_fill_manual(values = c("#FFF5EB",  "#F5C9B0", "#EC9E75")) +  # Custom colors
  theme_minimal(base_size = 10) +
  geom_label(aes(label = MAG.nr, y = 50), size=1)+
  labs(x = "Genus", y = "Proportion", fill = "Group") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))


ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/AS_MAGs/250305_Gao/000_plots/gao_gene_cluster_bar_approved_15000bp.jpeg',
       bar,
       height = 6,
       width = 20)




com_pattern_operon_meta %>% ggplot(aes(x=HMM,y= paste0(family, genus, species,tip_label),
                                       fill=HMM)) +
  geom_tile() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=3))




com_pattern_operon_meta %>% count(HMM)


###################################### start her #####################################
######################################################################################
# for gene cluster plot

darrow <- com_pattern_operon_meta %>% 
  group_by(tip_label, operon) %>% 
  
  # Determine if any epsB2 is on the forward strand
  mutate(wza_strand = ifelse(any(HMM == "epsB2" & Strandedness == 1), 1, 0)) %>%
  
  # Identify the first epsB2 based on strand direction
  mutate(first_epsB2_start = ifelse(wza_strand == 1, 
                                    min(Start_Position[HMM == "epsB2"], na.rm = TRUE), 
                                    max(Start_Position[HMM == "epsB2"], na.rm = TRUE)),
         first_epsB2_end = ifelse(wza_strand == 1, 
                                  min(End_Position[HMM == "epsB2"], na.rm = TRUE), 
                                  max(End_Position[HMM == "epsB2"], na.rm = TRUE))) %>%
  
  # Compute relative positions based on first epsB2
  mutate(start_relative = if_else(Strandedness == 1, 
                                  Start_Position - first_epsB2_start,
                                  End_Position - first_epsB2_end),
         end_relative = if_else(Strandedness == 1, 
                                End_Position - first_epsB2_start,
                                Start_Position - first_epsB2_end)) %>% 
  
  # Flip coordinates if the strand is reversed
  mutate(s = if_else(wza_strand == 1, start_relative, -start_relative),
         e = if_else(wza_strand == 1, end_relative, -end_relative))


darrow<- darrow %>% rename(Gene = HMM)

p_gene <- darrow %>% 
  ggplot(., aes(xmin = s, xmax = e, y = tip_label,
                fill = Gene,label = Gene)) +
  geom_gene_arrow(alpha = 0.9,
                  arrowhead_height = grid::unit(0.5, "mm"),
                  arrow_body_height = grid::unit(1, "mm"),
                  arrowhead_width = grid::unit(0.2, "mm")) + 
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 8))  +
  # geom_text(data = darrow %>% mutate(pos = (s+e)/2 -300),
  #           aes(x=pos, label = Gene),
  #           #nudge_y = 0.03, angle = 75,
  #           hjust= 0,
  #           size= 1.5)+
  scale_fill_viridis_d(option ="D") # +
#scale_x_continuous(limits = c(-500,12000))



ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/AS_MAGs/250305_Gao/000_plots/gao_gene_cluster_approved_15000bp.jpeg',
       p_gene,
       height = 20,
       width = 10)

poly_zoo <- read_tsv( '/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/AS_MAGs/250224_An/hit_genomes.txt',
                      col_names = F) %>% pull(1)
darrowplot <- darrow %>% filter(tip_label %in% poly_zoo)  #only hit genomes

# Reorder factor so genomes appear in desired order (e.g., by genus)
df <- darrowplot %>% ungroup() %>% 
  arrange(order, family,genus, species, tip_label,operon) %>%
  mutate(tip_label = factor(tip_label, levels = unique(tip_label)),
         tip_y = as.numeric(tip_label))

xpos <- -30000
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

colors_gene <- c(
  "#CCD5AE",  # sage
  "#A3B18A",  # olive green
  "#6B8E6E",  # forest sage
  "#E9EDC9",  # pale cream
  "#D4A373",  # muted ochre
  "#B07D62",  # warm brown
  "#7C9EB2"   # muted blue
)

colors_gene <- c(
  "#E9EDC9",  # light cream
  "#B8C480",  # light olive
  "#7AA974",  # medium green
  "#3F7D5C",  # dark green
  "#7C9EB2",  # muted teal
  "#D4A373",  # muted blue
  "#B07D62"   # muted ochre
)
p_gene_poly <- darrow %>% filter(tip_label %in% poly_zoo) %>% 
  ggplot(., aes(xmin = s, xmax = e, y = paste0(order,'::', family, '::', genus, '::', species, '::', tip_label),
                fill = Gene,label = Gene)) +
  geom_gene_arrow(alpha = 0.9,
                  arrowhead_height = grid::unit(2.5, "mm"),
                  arrow_body_height = grid::unit(2, "mm"),
                  arrowhead_width = grid::unit(1.2, "mm")) + 
  theme_minimal() +
  geom_text(data = darrowplot %>% mutate(pos = (s+e)/2 -500),
            aes(x=pos, label = Gene),
            #nudge_y = 0.03, angle = 75,
            hjust= 0,
            size= 1.5)+
  geom_segment(data = genus_segments ,
               aes(x = xpos, xend = xpos, y = y_min, yend = y_max),
               inherit.aes = FALSE, color = "grey30", size = 0.6) +
  # Genus labels
  geom_text(data = genus_segments ,
            aes(x = label_x, y = label_y, label = genus),
            inherit.aes = FALSE, hjust = 1, size = 3, 
            parse=T)+
 # scale_fill_viridis_d(option ="E")+
  scale_fill_manual(values = colors_gene) +
  coord_cartesian(xlim = c(-25000, 50000), clip = "off") +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        plot.margin = margin(5.5, 5, 5.5, 100) ,
        panel.grid.major = element_line(linewidth = 0.2)
  ) 

p_gene_poly

















# with polysaccharide
poly_zoo <- read_tsv( '/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/AS_MAGs/250224_An/hit_genomes.txt',
                      col_names = F)
poly_zoo <- poly_zoo %>% 
  mutate(hit = 'polysaccharide') %>% 
  rename(tip_label=1)

hit.7 <- hit %>% select(tip_label) %>% 
  mutate(hit7 = '7-gene cluster') 
hit.5 <- com_pattern%>% select(tip_label=Genome) %>% unique() %>% 
  mutate(hit7.essential = '5 essential genes')

# pepcterm
pepA <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/AS_MAGs/250305_Gao/pepA/250307_hmmsearch_pepA/hmmsearch_pepA_results.tsv') %>% 
  mutate(contig = sub(".*(contig_.*)_.*$", "\\1", Gene_ID))

pepA.count <-pepA %>% count(Genome)%>% rename('nr pepA' = n) %>% 
  mutate(pepA = 'yes')

poly_zoo.com <- poly_zoo %>% left_join(., hit.7, by = 'tip_label') %>% 
  left_join(., hit.5, by = 'tip_label') %>% 
  left_join(., pepA.count, by = c('tip_label' ='Genome')) %>% 
  left_join(., metadata, by = c('tip_label' ='bin'))






####################### sankey plot

library(ggsankey)


# Count occurrences for each column
taxonomic_cols <- c('class', 'order', 'family', 'genus', 'species')
df_with_counts <- hit %>%
  mutate(across(all_of(taxonomic_cols), ~ paste0(.x, " (", ave(seq_along(.x), .x, FUN = length))))
df_with_counts_meta <-metadata  %>%
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

df.san.hits <- df_with_counts_tot %>%
  make_long(  "class"    ,      "order"    ,     "family"     ,    "genus") %>% 
  mutate(node = factor(node, levels = rev(level_nodes$node)))
#make_long(  "class"    ,      "order"    ,     "family"     ,   'continent')



san_colors=c('#dad7cd','#E0C37E',
             '#e9edc9','#ccd5ae','#adc178','#a3b18a' ,'#588157')

psan <- ggplot(df.san.hits, aes(x = x, 
                                next_x = next_x, 
                                node = node, 
                                next_node = next_node,
                                fill = factor(node),
                                label = node)) +
  geom_sankey(type = 'alluvial') +
  geom_sankey_label(alpha = 0.9, type = 'alluvial')+
  theme_sankey(base_size = 12)+
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank()) +
  scale_fill_manual(values = scales::seq_gradient_pal(san_colors)(seq(0,1,length.out=n_distinct(df.san.hits$node)+4)))




df.san <- poly_zoo.com %>%
  mutate(hit7 = if_else(is.na(hit7), 'no', 'yes'),
         pepA = if_else(is.na(pepA), 'no', pepA)
         )%>%
  mutate(across(c(hit, hit7.essential), 
                ~ paste0(.x, " (", n(), ")"))) %>%
  group_by(hit7) %>% 
  mutate(hit7 = paste0(hit7, " (", n(), ")")) %>% 
  ungroup() %>%
  group_by(pepA) %>% 
  mutate(pepA = paste0(pepA, " (", n(), ")")) %>% 
  ungroup() %>% 
  rename('Genes from\n7-gene cluster' =hit7.essential,
         'Zoogloea\npolysaccharide\ngene cluster' =hit,
         'Genes in cluster' =hit7) %>% 
  make_long('Zoogloea\npolysaccharide\ngene cluster', 'Genes from\n7-gene cluster', 'Genes in cluster','pepA')

pp <-ggplot(df.san, aes(x = x, 
                   next_x = next_x, 
                   node = node, 
                   next_node = next_node,
                   fill = node,
                   label = node),
            linewidth=0) +
  geom_sankey() +
  geom_sankey_label()+
  theme_sankey(base_size = 16) +
  theme(legend.position = 'none',
        axis.title = element_blank())  +
  scale_fill_manual(values = #c( '#e9edc9','#ccd5ae','#adc178','#a3b18a' ,'#588157','#a3b18a' )
                      c('#ccd5ae','#8b9d6b','#a3b18a','#e9edc9','#a3b18a','#8b9d6b','#a3b18a'  )
                    ) 


complot <- cowplot::plot_grid(p_gene_poly,
                              pp,
                                # labs(tag = 'B'),
                                labels = 'AUTO',
                              label_fontface = "plain",
                              ncol=1,
                              rel_heights = c(1,0.175))
ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS13_gao_hitpolysaccharide_gene_cluster_approved_15000bp_sankey.jpeg',
       complot, dpi = 600,
       height = 15,
       width = 14)

