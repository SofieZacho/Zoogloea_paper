.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")
library(tidyverse)
library(ggtree)
library(ggnewscale)
library(data.table)
library(treeio)
library(gggenes)
library(ggtreeExtra)


# load hits from DIAMOND blastp search of the Zoogloea-polysaccharide gene cluster in Zoogloea genomes as described in the methods section of the manuscript.

polymer <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/EPS_gene_cluster_search/HQ_genomes/Zoogloea/250203_An/00_An_r220_total_blast_hits.tsv',
                    col_names = F)
names <- c('genome','qseqid','protein', 'pident', 'length','evalue',
           'bitscore',
           'mismatch', 'gaps','qstart', 'qend', 'qlen',  'sstart', 'send','slen', 'full_qseq')
# put col names on
colnames(polymer) <- names

polymer$genome <- str_remove_all(polymer$genome, "_(ASM|SRR|TT).*")

# load hits from annotation using DRAM as described in the methods section of the manuscript
hits <- read_tsv("/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/annotation/KEGG_zoogloea_25_01_24_HQ/00_condensed_HQ_dram.tsv") %>%
  filter(fasta != 'fasta')
hits$fasta <- str_remove_all(hits$fasta, "_(ASM|SRR|TT).*")

# load metadata with proposed names for the genomes
propnames <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/metadata/proposed_R226_new_names_updated.tsv', col_names = T)  
propnames[] <- lapply(propnames, function(x) if(is.factor(x) || is.character(x)) gsub("Zoogloea", "Z.", as.character(x)) else x)
propnames <- propnames %>% rename(midas_species_label = midas_species)

hits <- hits %>% rename(qseqid=1)
com <- polymer %>% left_join(hits, by = c('genome'='fasta','qseqid')) %>% 
  left_join(., propnames, by = 'genome')


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
  "asnB", "zooTRP", "wzxC", "zooGT8", "zooGT7", "zooGH", "zooM1", 
  "zooGT6", "zooGT5", "epsH", "wzy", "capK", "zooGT4", "zooGT3", 
  "asnH", "zooGT2", "zooGT1", "zooM2", "zooM3", "zooSA", "zooP", 
  "wzc", "etk", "wza", "wzi", "lolD"
)

# factorise
com_plot$protein <- factor(com_plot$protein, 
                           levels = rev(prot_vector))
com_plot %>% ggplot(aes(x=protein,y= paste0(tip_label),
                        fill=protein)) +
  geom_tile() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=2))



# load tree generated with GTDB-Tk and IQTREE as described in the methods section of the manuscript
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
                           species)) %>% 
  left_join(., propnames, by = c('tip_label'='genome'))



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

p.genome <- w %<+% meta 
p.genomec <-collapse(p.genome, node = 65, 'max', fill='grey80')

p.genomec <- p.genomec %>% rotate(99) 
p.genomec$data$genome <- gsub("GCA_009026175.1","GCF_009026175.1", p.genomec$data$genome)

p.genomec$data$genome  <- gsub('barcode02_Z_caeni','GCA_986281895*', p.genomec$data$genome )
p.genomec$data$genome  <- gsub('barcode03_Z_oleivorans','GCA_986281795*', p.genomec$data$genome )
p.genomec$data$genome  <- gsub('barcode04_Z_resiniphila','GCA_986340685*', p.genomec$data$genome )

p.genomec$data$label <- gsub('barcode02_Z_caeni','GCA_986281895*', p.genomec$data$label)
p.genomec$data$label <- gsub('barcode03_Z_oleivorans','GCA_986281795*', p.genomec$data$label)
p.genomec$data$label <- gsub('barcode04_Z_resiniphila','GCA_986340685*', p.genomec$data$label)

colspp <- c( "grey90",'black','black','black', "#EAC5AA", "#84541E","#D18E35", "#C6A07E","#f2dfcb", "#E0C37E", "#80CDC1", "#34978F",  "#4E8077", "#B7A27C", "#95825A", "#9c7605")
col <- c('#4E8077', "#D18E35", "#84593A","#E0C37E", 
         "#C6A07E","#D9A066",  
         "#f2dfa2","#95825A","#B7A27C",  "#9c7605",
         "#9FB8AD", '#e8c262',"#34978F","#A3B18A","#BFA5A0",  "#8dB8AD",'#A47559',
         "#ffe599","#9dB8bD")



ps <- p.genomec+ geom_treescale(x = 0.02, y = 12) + new_scale_color()+
  geom_tippoint(mapping=aes(color=proposed_name_new), 
                size=2, show.legend = F, alpha = 0.9, 
                shape = 15) +
  geom_tiplab(mapping = aes(label = paste0("italic('", proposed_name_new, "')~' (", genome, ")'")),
              geom='text', parse=T,
              size = 2, linesize = 0, offset = 0.13, show.legend = F, colour = 'black',
              alpha = 0.8, align = T, hjust=1) + #theme(legend.position = 'none') +
 # scale_color_viridis_d(option = 'F')
  scale_color_manual(values = col) 





############################################################################################################
######################################################
######### filter based on conserved pattern #########
############################################################################################################
pattern <- c('wza','etk','wzc',
             'zooSA','zooM2','zooGT1','epsH',
             'asnH','wzy') # found in many of the hits - really conserved
com_plot$tip_label <- gsub('barcode02_Z_caeni','GCA_986281895*', com_plot$tip_label)
com_plot$tip_label <- gsub('barcode03_Z_oleivorans','GCA_986281795*', com_plot$tip_label)
com_plot$tip_label <- gsub('barcode04_Z_resiniphila','GCA_986340685*', com_plot$tip_label)
# In the same gene cluster
com_plot_refilt <- com_plot %>%
  group_by(tip_label, scaffold,operon) %>%
  filter(all(pattern %in% protein)) %>%
  ungroup()
n_distinct(com_plot_refilt$tip_label)

com$genome <- gsub('barcode02_Z_caeni','GCA_986281895*', com$genome)
com$genome <- gsub('barcode03_Z_oleivorans','GCA_986281795*', com$genome)
com$genome <- gsub('barcode04_Z_resiniphila','GCA_986340685*', com$genome)


com$protein <- factor(com$protein, 
                      levels = rev(prot_vector))


colors <- c(
  "#4BA9A3", "#6CB5AB", "#8AC4B5", "#A5D1C0", "#BEDFCE",  # Soft teals to seafoam
  "#C9D9C6", "#D9D1B8", "#E6D8B8", "#F2E0C3", "#F9EACF",  # Sand and cream
  "#FCEFCF", "#F9E6B0", "#F2D58F", "#EBC27A", "#E1AE67",  # Pale gold to warm amber
  "#D29B61", "#C78B5B", "#B97A53", "#A96C4A", "#986042",  # Warm tans to brown
  "#87613F", "#765139", "#654430", "#553A2A", "#483024", "#3A261D" # Driftwood fade
)

heat <- ps + new_scale_fill()+ 
  geom_fruit(pwidth = 1, offset = 0.3,
                                          data = com_plot_refilt,
                                          geom = geom_tile,
                                          mapping = aes(
                                            y=tip_label,x=protein, fill = protein), 
                                          size = 0.03, color = 'grey50', axis.params=list(
                                            axis       = "x", 
                                            text.size  = 2, text.angle =90,
                                            hjust      = -0)) + 
  geom_fruit(pwidth = 1, offset = 0.07,
               data = com,
               geom = geom_tile,
               mapping = aes(
                 y=genome,x=protein, fill = protein), 
               size = 0.03, color = 'grey50', axis.params=list(
                 axis       = "x", #title = 'Metabolic reconstruction',
                 text.size  = 2, text.angle =90, 
                 hjust      = -0)
) +
  scale_fill_manual('Gene',values=colors)+
 # scale_fill_viridis_d(option ="F") +
  theme(
    legend.position = c(0.1, 0.62),
    legend.key.size = unit(0.4, 'cm'), #change legend key size
    legend.key.height = unit(0.4, 'cm'), #change legend key height
    legend.key.width = unit(0.4, 'cm'), #change legend key width
    legend.spacing.y = unit(0.02, 'cm'),
    legend.title = element_text(size=8), #change legend title font size
    legend.text = element_text(size=6),
    legend.background = element_rect(fill = 'transparent'))

ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS6_heat_inBGC_genes_an.jpeg',
       heat, dpi=600,
       height = 7,
       width = 10.5)










# remove one wzc for plotting
com_plot_refilt_un <- com_plot_refilt %>% filter(!qseqid %in% 'GCA_016714095.1_ASM1671409v1_genomic_JADJNV010000013.1_141') #splittedwzc gene

darrow <- com_plot_refilt%>% 
  group_by(tip_label,operon) %>% 
  mutate(wza_strand = ifelse(any(protein == "wzc" & strandedness == 1), 1, 0)) %>% 
  mutate(start_relative = if_else(strandedness == 1, 
                                  start_position - start_position[protein == "wzc"],
                                  end_position - end_position[protein == "wzc"]),
         end_relative = if_else(strandedness == 1, 
                                end_position- start_position[protein == "wzc"],
                                start_position - end_position[protein == "wzc"])) %>% 
  mutate(s = if_else(wza_strand == 1, start_relative, -start_relative),
         e = if_else(wza_strand == 1, end_relative, -end_relative))


darrow<- darrow %>% rename(Gene = protein,
                           Percent_identity = pident)

labs <- ps$data %>% filter(isTip, genus %in% 'g__Zoogloea') %>% arrange(y)%>% pull(label)




############################################################################################################
############################################################################################################
######### extract full clusters ########################################################################
############################################################################################################

# start is always lolD
lolD_pos <- darrow %>% filter(Gene %in% 'lolD') %>% 
  group_by(tip_label) %>% filter(s==min(s))

leng=38000
lolD_pos_label <- lolD_pos %>% select(tip_label, start_BGC=qseqid, scaffold, start_gene_cluster=start_position, start_gene_cluster_minus=end_position,
                                      strandedness) %>% 
  mutate(end_gene_cluster = if_else(strandedness == 1,
                                    start_gene_cluster+ leng,
                                    start_gene_cluster_minus- leng))

hits$fasta <- gsub('barcode02_Z_caeni','GCA_986281895*', hits$fasta)
hits$fasta <- gsub('barcode03_Z_oleivorans','GCA_986281795*', hits$fasta)
hits$fasta<- gsub('barcode04_Z_resiniphila','GCA_986340685*', hits$fasta)

combi <- lolD_pos_label %>% left_join(., hits, by =c('tip_label' = 'fasta', 'scaffold'='scaffold'))


combifilt <- combi %>% group_by(tip_label) %>% 
  filter(if_else(strandedness.x == 1,
                 start_position >= start_gene_cluster -5000,
                 start_position <= start_gene_cluster_minus +5000),
         if_else(strandedness.x == 1,
                 end_position < end_gene_cluster,
                 end_position > end_gene_cluster))  %>% 
  mutate(st = if_else(strandedness.x == 1, start_position-start_gene_cluster,end_position-start_gene_cluster_minus),
         en = if_else(strandedness.x == 1, end_position-start_gene_cluster,start_position-start_gene_cluster_minus)) %>% 
  mutate(s = if_else(strandedness.x == 1, st, -st),
         e = if_else(strandedness.x == 1, en, -en))%>% 
  mutate(s. = if_else(strandedness.x != strandedness.y, e, s),
         e. = if_else(strandedness.x != strandedness.y, s, e))


ss <-combifilt %>% ungroup() %>%  count(pfam_hits) %>% filter(n >1) %>% pull(pfam_hits)

combifilt.hits <- combifilt %>% left_join(.,darrow %>% select(tip_label,qseqid, Gene, Percent_identity), by = c('qseqid','tip_label'))


ds <- combifilt %>% select(tip_label, start_gene_cluster,start_gene_cluster_minus, end_gene_cluster,start_position, end_position,st,en,s,e,s.,e.,gene_position, 
                           strandedness.x, strandedness.y)







################ ################ ################ ################ ################ ################ 
################ tree with full clusters - homologs + other genes ################ ################ 
################ ################ ################ ################ ################ ################ ################ 
################ ################ ################ ################ ################ ################ ################ 


pnew <- ps+ new_scale_fill()+geom_fruit( pwidth = 10, offset =-100,
                                         data=combifilt.hits,
                                         geom = geom_gene_arrow,
                                         mapping= aes(xmin=(s.+14.5^4)*0.00001, xmax = (e.+14.5^4)*0.00001, y = tip_label,
                                                      fill = Gene),
                                         arrowhead_height = grid::unit(2.3, "mm"),
                                         arrow_body_height = grid::unit(1.8, "mm"),arrowhead_width = grid::unit(1, "mm"),
                                         size=0.2
) +
  
  scale_fill_manual(values = colors, na.value= 'grey95') +
  theme(
    legend.position = c(0.09, 0.72),
    legend.key.size = unit(0.4, 'cm'), #change legend key size
    legend.key.height = unit(0.4, 'cm'), #change legend key height
    legend.key.width = unit(0.4, 'cm'), #change legend key width
    legend.spacing.y = unit(0.02, 'cm'),
    legend.title = element_text(size=8), #change legend title font size
    legend.text = element_text(size=6),
    legend.background = element_rect(fill = 'transparent')) 




# Gene cluster as in query (Zoogloea resiniphila)
qplot <- combifilt.hits %>% filter(tip_label %in% 'GCA_986340685*',
                                   !is.na(Gene)) %>%  
  ggplot()+
  geom_gene_arrow( aes(xmin=s., xmax = e., y = tip_label,
                       fill = Gene),
                   arrowhead_height = grid::unit(5, "mm"),
                   arrow_body_height = grid::unit(3.5, "mm"),
                   arrowhead_width = grid::unit(2, "mm"), size=0.2
  ) +
  geom_text(data = combifilt.hits %>% filter(tip_label %in% 'GCA_986340685*',
                                             !is.na(Gene)) %>%  
              mutate(
                #pos = (s.+e.)/2 -400
                pos = if_else(s.<e., (s.+e.)/2 -450, (s.+e.)/2 -300)     
              ),
            aes(x=pos, label = Gene, y = tip_label),
            angle = 20,fontface='bold',
            hjust= -0.5, vjust = -0.8,
            size= 3, color='black')+
  theme_void()+
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_manual(values = colors, na.value= 'grey95') 



qspprep <- qplot+  labs(tag = 'A') + pnew+theme(legend.position = 'none')+ labs(tag = 'B') + 
  patchwork::plot_layout(heights = c(0.8,7.2)) 

ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS7_spprep_non_tree_genes_others_an.jpeg',
       qspprep, dpi = 600,
       height = 8,
       width = 10)



### spp rep
meta$tip_label <- gsub('barcode02_Z_caeni','GCA_986281895*', meta$tip_label)
meta$tip_label<- gsub('barcode03_Z_oleivorans','GCA_986281795*', meta$tip_label)
meta$tip_label<- gsub('barcode04_Z_resiniphila','GCA_986340685*', meta$tip_label)

metaHQ <- meta %>% filter(tip_label %in% hits$fasta)
metaHQ <- metaHQ %>% group_by(species) %>% arrange(-Completeness, .by_group = T)
spprep <- c( 'GCA_986281895*',
             '095-US2-13.bin.1.298',
             '072-HK-05.bin.1.452',
             '042-US1-105.bin.1.238',
             '001-SG-02.bin.2.148',
             '045-US1-87.bin.2.288', # only one
             
             
             'GCF_012927165.1',
             'GCA_986281795*',
             'GCF_009469605.1', #oryzae
             '048-UY-15.bin.1.146', # ramigera
             
             'GCF_002028455.1', #only one
             '098-UY-13.bin.1.166', #sp008011
             'GCA_009026175.1', #only one
             '028-AU-34.bin.1.573', #sp01565
             
             '034-IT-09.bin.1.324', # sp01672
             '033-IN-11.bin.1.356', #two midas
             '090-US1-124.bin.1.734', #only one
             'GCA_986340685*', # sp017309145
             '032-ES-07.bin.1.336'
)

combifilt.hits.spprep <- combifilt.hits %>% filter(tip_label %in% spprep)

# order
# Assuming combifilt.hits.spprep is your dataframe
combifilt.hits.spprep.meta <- combifilt.hits.spprep%>% left_join(., metaHQ) 
# %>% 
#   left_join(., propnames, by = c('tip_label'='genome'))
combifilt.hits.spprep.meta$tip_label <- gsub("GCA_009026175.1","GCF_009026175.1", combifilt.hits.spprep.meta$tip_label)

combifilt.hits.spprep.meta$tip_label %>% unique()
combifilt.hits.spprep.meta$proposed_name_new %>% unique()

# Compute mean perc_id per species
levl <- combifilt.hits.spprep.meta %>%
  group_by(tip_label) %>%
  summarise(mean_perc_id = mean(Percent_identity, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_perc_id)) %>% pull(tip_label)



# Create a factor column based on sorted species order
combifilt.hits.spprep.meta <- combifilt.hits.spprep.meta %>%
  mutate(tip_label = factor(tip_label, levels = rev(levl))) %>%
  mutate(y_label = paste0("italic('", proposed_name_new, "')~' (", tip_label, ")'")) %>% 
  mutate(y_label = factor(y_label, levels = unique(y_label)))



qplot <-combifilt.hits %>% filter(tip_label %in% 'GCA_986340685*',
                                  !is.na(Gene)) %>%  
  ggplot()+
  geom_gene_arrow( aes(xmin=s., xmax = e., y = tip_label,
                       fill = Gene),
                   arrowhead_height = grid::unit(5, "mm"),
                   arrow_body_height = grid::unit(3.5, "mm"),
                   arrowhead_width = grid::unit(2, "mm"), size=0.2
  ) +
  geom_text(data = combifilt.hits %>% filter(tip_label %in% 'GCA_986340685*',
                                             !is.na(Gene)) %>%  
              mutate(
                #pos = (s.+e.)/2 -400
                pos = if_else(s.<e., (s.+e.)/2 -450, (s.+e.)/2 -300)     
              ),
            aes(x=pos, label = Gene, y = tip_label),
            angle = 20,fontface='bold',
            hjust= -0.42, vjust = -0.7,
            size= 2.5, color='black')+
  theme_void()+
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_manual(values = colors, na.value= 'grey95')  +
  coord_cartesian(xlim = c(-3300, 39000), clip = "off") 



pgene.spprep <-combifilt.hits.spprep.meta  %>% 
  ggplot() +
  geom_gene_arrow(aes(xmin=s., xmax = e., y = y_label,
                      fill = Gene), 
                  arrowhead_height = grid::unit(4.6, "mm"),
                  arrow_body_height = grid::unit(4, "mm"),
                  arrowhead_width = grid::unit(1.8, "mm"), 
                  size=0.2 )+
  
  geom_text(data = combifilt.hits.spprep.meta %>% 
              filter(Gene %in% c('zooTRP')) %>% 
              mutate(pos = if_else(s.<e., s+30, s +200)),
            aes(x=pos, label = Gene, y = y_label),
            color = 'black',
            fontface ='bold', hjust= 0,size= 1.45)+
  
  geom_text(data = combifilt.hits.spprep.meta %>% mutate(
    pos = if_else(s.<e., s+30, 
                  if_else(Gene %in% c('asnB'), s+600,s +200) )    ),
  aes(x=pos, label = Gene, y = y_label,
      color = if_else(Gene %in% c('zooGT7','zooGT8','zooTRP','wzxC','asnB'),
                      'white','black'
                      )),
  fontface ='bold', hjust= 0,size= 1.4)+

  scale_color_manual(values = c('black','grey95')) +
  scale_fill_manual(values=colors, na.value= 'grey95')+
  scale_y_discrete(labels = function(x) parse(text = x))+
  theme_minimal() +
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_text(color = 'black')) +
  coord_cartesian(xlim = c(-3300, 39000), clip = "off") 



# put plots together
qspprep<-cowplot::plot_grid(qplot, pgene.spprep, ncol = 1, labels= c('A','B'),
                            label_fontface = "plain", label_size = 12,
                            rel_heights = c(1,7.2),
                            label_y = c(1.05, 1.04),
                            align = 'v')
qspprep

ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/Fig4_spprep_tree_genes_others_an.jpeg',
       qspprep, dpi=600,
       height = 6,
       width = 12.5)
















