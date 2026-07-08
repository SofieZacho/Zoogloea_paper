.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")

library(gggenes)
library(data.table)
library(Biostrings)
library(tidyverse)

############### taxonomy +length queries  ######################################################################################################################
cwd <- "/home/bio.aau.dk/kl42gg/projects/rethink/epsSMASH/epsProtocol"
tax_file <- file.path(cwd, "gtdb/gtdb_taxonomy.tsv")
tax <- read_tsv(tax_file, col_names = c("genome", "taxonomy"), col_types = "cc")
tax_sep <- tax %>%
  separate(taxonomy, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") 

queries <- readAAStringSet('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/gene_cluster_literature/zooglan.faa')
query_data <- data.frame(
  query = names(queries),
  query_length = nchar(queries)) %>%
  group_by(query) %>%
  summarise(query_length = mean(query_length))
###############################################################################################################################################

filtered_genomes <- read_csv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/GTDBr220/250227_An/core_genes_csv_without_wzi')


conserved <- c('wzi',
               'wza','etk','wzc','zooP','zooSA','zooM3','zooM2','zooGT1','epsH', 'asnH','wzy')


data_core_tax <- filtered_genomes %>% left_join(., tax_sep, by ='genome')
# investigate taxonomy
st <- data_core_tax %>% select("domain", "phylum", "class", "order", "family", "genus", "species") %>% 
  unique()
st %>% ungroup() %>%  count(genus) %>% arrange(-n)

n_distinct(data_core_tax$genome)

data_core_tax %>% group_by(genome) %>% count(name) %>% group_by(name) %>%  summarise(me = mean(n))
data_core_tax <- data_core_tax %>% left_join(., query_data, by = c('name'='query'))

########## plot
# order of genes
prot_vector <- c(
  "asnB", "zooTRP", "wzxC", "zooGT8", "zooGT7", "zooGH", "zooM1", 
  "zooGT6", "zooGT5", "epsH", "wzy", "capK", "zooGT4", "zooGT3", 
  "asnH", "zooGT2", "zooGT1", "zooM2", "zooM3", "zooSA", "zooP", 
  "wzc", "etk", "wza", "wzi", "lolD"
)

# factorise
data_core_tax$name <- factor(data_core_tax$name, 
                             levels = rev(prot_vector))
data_core_tax %>% ggplot(aes(x=name,y=species,fill=name)) +
  geom_tile() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=3))

n_distinct(data_core_tax$genome)


########################################################################
################## try filter regions ############################
########################################################################
com_sorted <- data_core_tax %>% 
  mutate(start = as.numeric(start),
         end=as.numeric(end)) %>% 
  group_by(genome, contig) %>% 
  arrange(start, .by_group = T) %>%
  group_by(protID) %>%
  filter(evalue == min(evalue))

n <- 15000
com_filtered <- com_sorted  %>% 
  group_by(genome, contig) %>% 
  mutate(prior =lag(end), after = lead(start)) %>% 
  mutate(dif_bef = start - lag(end),
         dif_aft = lead(start) - end) %>% 
  filter(!if_all(c(dif_bef, dif_aft), is.na)) %>% #removes hits from contigs with only 1 hit
  filter(dif_bef < n | dif_aft < n | is.na(dif_bef) & dif_aft < n | is.na(dif_aft) & dif_bef < n)


com_plot<- com_filtered %>% 
  ungroup() %>% 
  group_by(genome) %>% 
  mutate(
    operon_start = ifelse(dif_bef > n | is.na(dif_bef), "start", NA),
    operon = cumsum(!is.na(operon_start)),  # Count occurrences of "start" within each group
  ) %>%
  group_by(genome, operon) %>%
  filter(n_distinct(name) >2) # filter out hits from operons with 3 or less genes 
n_distinct(com_plot$genome)
n_distinct(com_plot$contig)

com_plot %>% ggplot(aes(x=name,y=species,fill=name)) +
  geom_tile() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=3))

## filter to one contig
com_plot_region <- com_plot %>%
  group_by(genome, contig,operon) %>%
  filter(all(conserved %in% name)) %>%
  #filter(length(unique(query)) >= min_genes) %>%
  ungroup()
n_distinct(com_plot_region$genome)


com_plot_region %>% ggplot(aes(x=name,y=species,fill=name)) +
  geom_tile() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=5),
        axis.text.x = element_text(angle=90))








### count per genome how many of each gene


com_plot%>% ungroup() %>% select("domain", "phylum", "class", "order", "family", "genus", "species") %>% 
  unique() %>% ungroup() %>%  count(family) %>% arrange(-n)
com_plot_region %>% ungroup() %>% select("domain", "phylum", "class", "order", "family", "genus", "species") %>% 
  unique() %>% ungroup() %>%  count(family) %>% arrange(-n)




#################################################################################################################################################################
########################################################
################ re-iterate filtering as part of pattern ################
########################################################
pattern <- c('wzi','wza','etk','wzc','zooP','zooSA','zooM3','zooM2','zooGT1','epsH','asnH','wzy') # found in many of the hits - really conserved

data_refilt <- com_plot %>%
  group_by(genome) %>%
  filter(all(pattern %in% name)) %>%
  #filter(length(unique(query)) >= min_genes) %>%
  ungroup()

n_distinct(data_refilt$genome)
n_distinct(com_plot_region$genome)


################## try filter regions ############################

# In the same gene cluster
com_plot_refilt <- com_plot %>%
  group_by(genome, contig,operon) %>%
  filter(all(pattern %in% name)) %>% 
  ungroup()
n_distinct(com_plot_refilt$genome)

com_plot_refilt$name <- factor(com_plot_refilt$name, 
                               levels = rev(prot_vector))

com_plot_refilt %>% ggplot(aes(x=name,y=species,fill=name)) +
  geom_tile() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=2))




### count per genome how many of each gene


com_plot_refilt %>% select("domain", "phylum", "class", "order", "family", "genus", "species") %>% 
  unique() %>% ungroup() %>%  count(family) %>% arrange(-n)

data_refilt %>% select("domain", "phylum", "class", "order", "family", "genus", "species") %>% 
  unique() %>% ungroup() %>%  count(family) %>% arrange(-n)

taxi<- com_plot_refilt %>% select("domain", "phylum", "class", "order", "family", "genus", "species") %>% 
  unique() 


######## all genes found on same contig


# for gene cluster plot
darrow <- com_plot_refilt%>% 
  group_by(genome,operon) %>% 
  filter(name != "etk" | evalue == min(evalue[name == "etk"])) %>% 
  mutate(wza_strand = ifelse(any(name == "etk" & strand == 1), 1, 0)) %>% 
  mutate(start_relative = if_else(strand == 1, 
                                  start - start[name == "etk"],
                                  end - end[name == "etk"]),
         end_relative = if_else(strand == 1, 
                                end- start[name == "etk"],
                                start - end[name == "etk"])) %>% 
  mutate(s = if_else(wza_strand == 1, start_relative, -start_relative),
         e = if_else(wza_strand == 1, end_relative, -end_relative))




###### plot re-filtered hits on gtdb tree
# make tree
library(ggtree)
library(treeio)
library(ggtreeExtra)
library(ggnewscale)
tree <- read.tree('/databases/GTDB/gtdbtk_packages/GTDB-TK_release220_2024_04_18/pplacer/gtdb_r220_bac120.refpkg/gtdb_r220_bac120_decorated_unrooted.tree')



p <- ggtree(tree, layout = "fan",
            size = 0.0025)


tax_GTDB <- read_tsv("/databases/GTDB/gtdbtk_packages/GTDB-TK_release220_2024_04_18/taxonomy/bac120_taxonomy_r220_reps.tsv",
                     col_names = F) %>%
  rename(genome = X1) %>% 
  separate(X2, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") 

dbroad <- com_plot_region %>% select(genome) %>% unique() %>% mutate(Hit = 'Zoogloea gene cluster') ## change here
con <- tax_GTDB %>% left_join(., dbroad, by = c("genome")) %>% 
  rename(tip_label = genome)


family_count <- con %>% filter(Hit %in% 'Zoogloea gene cluster') %>%  
  group_by(family) %>%
  mutate(unique_genomes = n_distinct(tip_label)) %>%
  mutate(fam.label = paste0(family, " (", unique_genomes, ")")) %>% ungroup() %>% select(tip_label, fam.label)

# Create a named vector for renaming legend labels
con <- left_join(con, family_count, by = 'tip_label')


chip <- con %>% mutate(
  phyl_col = case_when(phylum == 'p__Pseudomonadota' ~ 'Pseudomonadota', 
                       phylum == 'p__Bacillota_A' ~ 'Bacillota_A',
                       phylum == 'p__Bacillota' ~'Bacillota',
                       phylum == 'p__Bacillota_I' ~'Bacillota_I',
                       phylum == 'p__Desulfobacterota' ~'Desulfobacterota',
                       phylum == 'p__Bacteroidota' ~ 'Bacteroidota',
                       phylum == 'p__Actinomycetota' ~ 'Actinomycetota',
                       phylum == 'p__Patescibacteria' ~ 'Patescibacteria',
                       phylum == 'p__Chloroflexota' ~ 'Chloroflexota',
                       #phylum == 'p__Cyanobacteriota' ~ 'Cyanobacteriota',
                       #phylum == 'p__Planctomycetota' ~ 'Planctomycetota',
                       #phylum == 'p__Acidobacteriota' ~ 'Acidobacteriota',
                       #phylum == 'p__Verrucomicrobiota' ~ 'Verrucomicrobiota',
                       .default = 'Other'
  )) %>% mutate(phyl_col = factor(phyl_col, 
                                  levels = c('Pseudomonadota','Bacillota_A','Bacteroidota', 'Actinomycetota', 'Patescibacteria', 'Bacillota', 'Bacillota_I',
                                             'Chloroflexota','Desulfobacterota','Other'))) 





phyla_colors = c('#F5F5F5','#452929','#84541E','#9c7605', '#6a170e','#bf4e4e','#D18E35', '#E0C37E' ,'#80CDC1','#34978F')
phyla_colors <- rev(phyla_colors)


tip_data <- chip %>% select(tip_label, phylum,class, order, family,genus, species, Hit,phyl_col,fam.label) 

p_save <- p %<+% tip_data 

# combine
p_savephyl <- p_save +
  geom_tippoint(
                data = ~subset(.x, !is.na(fam.label)),
                aes(fill = fam.label),
                shape=21,
                size=2.2)+
  scale_fill_discrete(
    na.value = NA,
    'Family (Genome count)') 


p_savephyl_zooglan <- p_savephyl + new_scale_fill() + geom_fruit(
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
                    values = phyla_colors)

ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS9_an_gtdbr220_15000bp.jpeg',
       p_savephyl_zooglan, dpi=500)


# ggsave('~/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/GTDBr220/250227_An/000_plots/an_gtdbr220_tree_gene_cluster_15000bp.jpeg',
#        p_savephyl_zooglan)







