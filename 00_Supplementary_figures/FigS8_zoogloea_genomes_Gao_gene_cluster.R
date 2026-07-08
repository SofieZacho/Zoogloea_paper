.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(data.table)
library(treeio)
library(gggenes)


# load hits from 

polymer <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/Zoogloea/250305_Gao/00_Gao_r220_total_blast_hits.tsv',
                    col_names = F)
names <- c('genome','qseqid','protein', 'pident', 'length','evalue',
           'bitscore',
           'mismatch', 'gaps','qstart', 'qend', 'qlen',  'sstart', 'send','slen', 'full_qseq')
# put col names on
colnames(polymer) <- names

polymer$genome <- str_remove_all(polymer$genome, "_ASM.*")
polymer$genome <- str_remove_all(polymer$genome, "_SRR.*")
polymer$genome <- str_remove_all(polymer$genome, "_TT.*")

# load hits from annotation
hits <- read_tsv("/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/annotation/KEGG_zoogloea_25_01_24_HQ/00_condensed_HQ_dram.tsv") %>%
  filter(fasta != 'fasta')
hits$fasta <- str_remove_all(hits$fasta, "_ASM.*")
hits$fasta <- str_remove_all(hits$fasta, "_SRR.*")
hits$fasta <- str_remove_all(hits$fasta, "_TT.*")

hits <- hits %>% rename(qseqid=1)
com <- polymer %>% left_join(hits, by = c('genome'='fasta','qseqid'))


#### sort stuff
com_sorted <- com %>% 
  group_by(genome, scaffold) %>% 
  arrange(start_position, .by_group = T) %>%
  group_by(qseqid) %>%
  filter(evalue == min(evalue))

n <- 15000
com_filtered <- com_sorted  %>% 
  group_by(genome, scaffold) %>% 
  mutate(prior =lag(end_position), after = lead(start_position)) %>% 
  mutate(dif_bef = start_position - lag(end_position),
         dif_aft = lead(start_position) - end_position) %>% 
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
  filter(n_distinct(protein) >1) %>%
  rename(tip_label = genome)

n_distinct( com_plot$tip_label)
ss <- com_plot %>% select(tip_label, operon) %>% unique()

# order of genes
prot_vector <- c(
  "epsB2", "prsK", "prsR", "prsT", "ugd", "mltE", "degQ2"
)

# factorise
com_plot$protein <- factor(com_plot$protein, 
                           levels = prot_vector)
com_plot %>% ggplot(aes(x=protein,y= paste0(tip_label),
                        fill=protein)) +
  geom_tile() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=2))



# read trees
t_genome <- '/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/classification_tree/tree/250429_zoogloea_gtdb_midas_HQ_R226/gtdbtk.bac120.user_msa.fasta.treefile'
# Taxonomy file
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



# Tree file
tree <- read.tree(t_genome)
trg <- phytools::midpoint.root(tree)
############################################################
############ plot full genome tree ########################
#remove everything in tip labels after '_SRR', '_ASM', and also remove 'RS_' in front of accession
trg$tip.label <- gsub('_SRR.*|_ASM.*|_TT.*|RS_', '', trg$tip.label)
p <- ggtree(trg ,size = 0.5 ) 


df_boot <- trg %>% as.treedata %>% as_tibble
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 90), "90-100 %",
                         ifelse((df_boot$label >= 70),"70-90 %","0-70 %"))
df_boot <- df_boot[, c("node","status")]

w <- ggtree(trg) 
w <- w %<+% df_boot + geom_nodepoint(aes(color=status), size=1.5) +
  scale_color_manual(values = c("90-100 %" = "grey0", "70-90 %"="gray50", "0-70 %"="gray75"),
                     na.translate = FALSE) +
  labs( color="Bootstrap")


propnames <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/metadata/proposed_R226_new_names_updated.tsv', col_names = T)  
propnames[] <- lapply(propnames, function(x) if(is.factor(x) || is.character(x)) gsub("Zoogloea", "Z.", as.character(x)) else x)
propnames <- propnames %>% rename(midas_species_label = midas_species)


p.genome <- w %<+% meta %<+% propnames
p.genomec <-collapse(p.genome, node = 65, 'max', fill='grey80')

p.genomec <- p.genomec %>% rotate(99) 
p.genomec$data$genome <- gsub("GCA_009026175.1","GCF_009026175.1", p.genomec$data$genome)
#colspp <- c( "grey90",'black','black','black', "#EAC5AA", "#84541E","#D18E35", "#C6A07E","#f2dfcb", "#E0C37E", "#80CDC1", "#34978F",  "#4E8077", "#B7A27C", "#95825A", "#9c7605")
colspp <- c(    '#9dB8bD', "#84541E","#D18E35",'#ddcea1',"#f2dfa2", "#D9A066","#E0C37E", "#C6A07E","#8dB8AD",  "#4E8077", 
                  "#B7A27C", "#95825A", "#9c7605",
                  "#EAC5a1", '#c6b46f',"#EAC5AA","#BFA5A0", 
                  "#A3B18A","#9FB8AD")


p.genomec$data$genome  <- gsub('barcode02_Z_caeni','GCA_986281895*', p.genomec$data$genome )
p.genomec$data$genome  <- gsub('barcode03_Z_oleivorans','GCA_986281795*', p.genomec$data$genome )
p.genomec$data$genome  <- gsub('barcode04_Z_resiniphila','GCA_986340685*', p.genomec$data$genome )

p.genomec$data$label <- gsub('barcode02_Z_caeni','GCA_986281895*', p.genomec$data$label)
p.genomec$data$label <- gsub('barcode03_Z_oleivorans','GCA_986281795*', p.genomec$data$label)
p.genomec$data$label <- gsub('barcode04_Z_resiniphila','GCA_986340685*', p.genomec$data$label)


ps <- p.genomec+ geom_treescale(x = 0.02, y = 12) + new_scale_color()+
  geom_tippoint(mapping=aes(color=proposed_name_new),
                size=2, show.legend = F,
                shape = 15)+
  geom_tiplab(mapping = aes(label = paste0(proposed_name_new,' (',genome,')')),
              geom='text',
              size = 1.5, linesize = 0, offset = 0.12, show.legend = F, colour = 'black',
              alpha = 0.8, align = T, hjust=1)+ #theme(legend.position = 'none') +
  scale_color_manual(values = colspp) 




############################################################################################################
######################################################
######### filter based on conserved pattern #########
############################################################################################################
pattern <- c("epsB2", "prsK", "prsR", "prsT", "ugd", "mltE", "degQ2") # found in many of the hits - really conserved

com_plot$tip_label <- gsub('barcode02_Z_caeni','GCA_986281895*', com_plot$tip_label)
com_plot$tip_label <- gsub('barcode03_Z_oleivorans','GCA_986281795*', com_plot$tip_label)
com_plot$tip_label <- gsub('barcode04_Z_resiniphila','GCA_986340685*', com_plot$tip_label)


# In the same gene cluster
com_plot_refilt <- com_plot %>%
  group_by(tip_label, scaffold,operon) %>%
  filter(all(pattern %in% protein)) %>%
  ungroup()
n_distinct(com_plot_refilt$tip_label)



########################################################################################################### 

darrow <- com_plot_refilt%>% 
  group_by(tip_label,operon) %>% 
  mutate(wza_strand = ifelse(any(protein == "epsB2" & strandedness == 1), 1, 0)) %>% 
  mutate(start_relative = if_else(strandedness == 1, 
                                  start_position - start_position[protein == "epsB2"],
                                  end_position - end_position[protein == "epsB2"]),
         end_relative = if_else(strandedness == 1, 
                                end_position- start_position[protein == "epsB2"],
                                start_position - end_position[protein == "epsB2"])) %>% 
  mutate(s = if_else(wza_strand == 1, start_relative, -start_relative),
         e = if_else(wza_strand == 1, end_relative, -end_relative))


darrow <- com_plot_refilt %>% 
  group_by(tip_label, operon) %>% 
  
  # Determine if any epsB2 is on the forward strand
  mutate(wza_strand = ifelse(any(protein == "epsB2" & strandedness == 1), 1, 0)) %>%
  
  # Identify the first epsB2 based on strand direction
  mutate(first_epsB2_start = ifelse(wza_strand == 1, 
                                    min(start_position[protein == "epsB2"], na.rm = TRUE), 
                                    max(start_position[protein == "epsB2"], na.rm = TRUE)),
         first_epsB2_end = ifelse(wza_strand == 1, 
                                    min(end_position[protein == "epsB2"], na.rm = TRUE), 
                                    max(end_position[protein == "epsB2"], na.rm = TRUE))) %>%
  
  # Compute relative positions based on first epsB2
  mutate(start_relative = if_else(strandedness == 1, 
                                  start_position - first_epsB2_start,
                                  end_position - first_epsB2_end),
         end_relative = if_else(strandedness == 1, 
                                end_position - first_epsB2_start,
                                start_position - first_epsB2_end)) %>% 
  
  # Flip coordinates if the strand is reversed
  mutate(s = if_else(wza_strand == 1, start_relative, -start_relative),
         e = if_else(wza_strand == 1, end_relative, -end_relative))


darrow<- darrow %>% rename(Gene = protein,
                           Percent_identity = pident)

labs <- ps$data %>% filter(isTip, genus %in% 'g__Zoogloea') %>% arrange(y)%>% pull(label)


darrow$tip_label <- gsub("GCA_009026175.1","GCF_009026175.1", darrow$tip_label)


pssa <- ps+ new_scale_fill()+ geom_fruit(pwidth = 4, offset = 0.001,
                                         data = darrow %>% filter(!s<0),
                                         geom = geom_gene_arrow,
                                         mapping = aes(xmin = (s+14^4)*0.00001, xmax = (e+14^4)*0.00001, y = tip_label, 
                                                       fill = Gene, 
                                                     #  alpha = Percent_identity
                                                       ), 
                                         arrowhead_height = grid::unit(2.2, "mm"),
                                         arrow_body_height = grid::unit(1.8, "mm"),arrowhead_width = grid::unit(1.8, "mm"),
                                         size=0.2
) +
  scale_fill_viridis_d(option ="E") +
  theme(
    legend.position = c(0.09, 0.62),
    legend.key.size = unit(0.4, 'cm'), #change legend key size
    legend.key.height = unit(0.4, 'cm'), #change legend key height
    legend.key.width = unit(0.4, 'cm'), #change legend key width
    legend.spacing.y = unit(0.02, 'cm'),
    legend.title = element_text(size=8), #change legend title font size
    legend.text = element_text(size=6),
    legend.background = element_rect(fill = 'transparent'))





ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS8_tree_genes_gao_withoutalpha.jpeg',
       pssa, dpi=600,
       height = 5.8,
       width = 5)















