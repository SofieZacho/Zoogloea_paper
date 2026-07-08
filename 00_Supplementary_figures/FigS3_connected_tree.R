
.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(data.table)
library(treeio)

# read trees
t_path16S <- '/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/classification_tree/tree_16s/250121_zoogloea_gtdb_midas_HQ/output/msa_trim_10pct/trim_msa_16s_gtdb_midas_combined_1000bp_HQ.fna_10pct.fa.treefile'
t_genome <- '/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/classification_tree/tree/250429_zoogloea_gtdb_midas_HQ_R226/gtdbtk.bac120.user_msa.fasta.treefile'


tree16s <- read.tree(t_path16S)
tr16s <- phytools::midpoint.root(tree16s)
treeg <- read.tree(t_genome)
trg <- phytools::midpoint.root(treeg)


propnames <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/metadata/proposed_R226_new_names_updated.tsv', col_names = T)  
trg$tip.label <- gsub("_SRR.*|_ASM.*|_TT.*|.fa", "", trg$tip.label)
propnames[] <- lapply(propnames, function(x) if(is.factor(x) || is.character(x)) gsub("Zoogloea", "Z.", as.character(x)) else x)
propnames <- propnames %>% rename(midas_species_label = midas_species)


pp <- ggtree(trg ,size = 0.6 ) 
p.genome <-pp %<+% propnames 



meta16s <- data.frame(tip_label =tr16s$tip.label)
sintax <- fread(input = '/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/classification_tree/tree_16s/241218_zoogloea_gtdb_midas/output/MiDAS_classification/16S_all_MiDAS5.3.tsv',
                sep = "\t", fill = TRUE, header = FALSE, data.table = TRUE)[, c(1, 2, 3)]
colnames(sintax) <- c('tip_label','midas_16S_tax','ID')
sintax <- sintax %>% mutate(tip_label = sub(" .*", "", tip_label))
# remove everything after the not last but second to last '_' in tip labels (including both '_') 
meta16s$genome <- gsub('_[^_]*$', '', gsub('_[^_]*$', '', meta16s$tip_label))
meta16s$genome <- gsub('_SRR.*|_ASM.*|_TT.*|RS_', '', meta16s$genome)
meta16s <- meta16s %>% left_join(., propnames, by = 'genome') %>% 
  left_join(., sintax, by = 'tip_label') %>% 
  separate(midas_16S_tax, into = c('midas_domain', 'midas_phylum', 'midas_class', 'midas_order', 'midas_family', 'midas_genus', 'midas_species'), sep = ',') %>% 
  mutate(
    midas_species = if_else(ID < 98.7,
                            paste0('*',midas_species),
                            midas_species)) 



# check outliers and remove
p16 <- ggtree(tr16s ,size = 0.3 ) 
pp16 <- p16 %<+% meta16s 

ps16midas <- pp16+ geom_treescale() +
  geom_tiplab(mapping = aes(label = paste0(midas_genus,'::',midas_species,'::',genome),
                            fill = midas_species),geom='label',
              size = 2, linesize = 0.2) +theme(legend.position = 'none') +
  geom_tippoint(mapping=aes(color=midas_species),size=1, stroke=0) 

psav16 <- ps16midas+ #geom_nodelab(aes(label = node), size=2, geom= 'label') +
  geom_tiplab(aes(label = node), size=2, geom= 'label')


# remove 16s from other genera
tips.other.genera <- c(1,2,3,4,5,6,7,8,9) 
rem16 <- ape::drop.tip(tr16s, tips.other.genera)
# list which genomes these are from:
setdiff( tr16s$tip.label, rem16$tip.label) %>% sort()
# out on metadata
prem16 <- ggtree(rem16 ,size = 0.2 ) %<+% meta16s 








############################################################
############ plot ##########################################
############################################################
p.genome$data$label <- gsub("GCA_009026175.1","GCF_009026175.1", p.genome$data$label)

p.genome$data$label <- gsub('barcode02_Z_caeni','GCA_986281895*', p.genome$data$label)
p.genome$data$label <- gsub('barcode03_Z_oleivorans','GCA_986281795*', p.genome$data$label)
p.genome$data$label <- gsub('barcode04_Z_resiniphila','GCA_986340685*', p.genome$data$label)


df_boot <- p.genome %>% as.treedata %>% as_tibble
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 90), "90-100 %",
                         ifelse((df_boot$label >= 70),"70-90 %","0-70 %"))
df_boot <- df_boot[, c("node","status")]


w <- p.genome %<+% df_boot + geom_nodepoint(aes(color=status), size=1.5) +
  scale_color_manual(values = c("90-100 %" = "grey0", "70-90 %"="gray50", "0-70 %"="gray75"),
                     na.translate = FALSE) +
  labs( color="Bootstrap")



ps <- w + new_scale_color() +
  geom_treescale(x=0, y=20, offset=1)  
  # geom_nodelab(aes(label = node), size =3, color='green3')

ps <- ps %>% 
  rotate(99) %>% rotate(70)
   # rotate(68) %>%
   # rotate(69)%>%
   # #rotate(93)%>%
   # rotate(95)

library(glue)

ps_new <- ps +
  geom_segment(aes(x=0.234, y=y, xend=0.32, #yend=ylim, 
                   color = if_else(is.na(proposed_name_new), 'g__Thauera',proposed_name_new)),
               data=ps$data %>% filter(isTip),
               alpha = 0.8, linewidth=5.7)+
  geom_tiplab(mapping = aes(label =
                              glue("italic('{proposed_name_new}')*' ({label})'")),
    geom='text', parse = TRUE,linetype='blank',
    size = 2.4, linesize = 0, offset = 0.003,
   alpha = 1, align = T) 

# collapse outgroup
p.genomec <-collapse(ps_new, node = 65, 'min', fill='grey80')

prem16r <- prem16 %>% 
  rotate(333) %>%
  rotate(369)

prem16rc <-collapse(prem16r, node = 253, 'min', fill='grey80')

#  rotate(379)
d.genome <- ps$data
d.ribo <- prem16r$data
d.ribo$genome <- gsub("GCA_009026175.1","GCF_009026175.1", d.ribo$genome)


## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d.ribo$x <- max(d.ribo$x) - d.ribo$x + max(d.genome$x) + 0.2
d.ribo$y <- d.ribo$y * (max(d.genome$y) / max(d.ribo$y))

varx <- 0.48

ppp <- p.genomec +  geom_segment(aes(x=varx, y=y, xend=0.61, #yend=ylim, 
                                     color = if_else(midas_genus %in% 'g:Thauera',
                                                     'g:Thauera',
                                                     if_else(!str_detect(midas_species,'^\\*s'),
                                                             proposed_name_new, 'No_spp'))),
                                 data=d.ribo %>% filter(isTip),
                                 alpha = 0.8, linewidth=1.45) +
  geom_tree(data=d.ribo) 



dd <- bind_rows(d.genome %>% mutate(genome =label), d.ribo) %>% # draw lines between genomes and 16S pairs
  filter(isTip) %>%
  mutate(
    xlim = if_else(x>max(d.genome$x), x, NA),
    ylim = if_else(x>max(d.genome$x), y, NA),
    x = if_else(x<=max(d.genome$x), x, NA),
    y = if_else(x<=max(d.genome$x), y, NA)) #%>%
#select(x, y, xlim, ylim, genome)

dd <- full_join((dd %>%
                   select(genome, x, y) %>%
                   filter(!is.na(x), !is.na(y)) %>%
                   distinct()),
                (dd %>%
                   select(genome, xlim, ylim) %>%
                   filter(!is.na(xlim), !is.na(ylim)) %>%
                   distinct())) %>% 
  left_join(., propnames, by = 'genome')


# try to add midas label on the right
# Compute segment positions for each genus
species_segments <- d.ribo %>%
  filter(isTip) %>% 
  group_by(midas_species) %>%
  filter(!midas_genus %in% 'g:Thauera') %>% 
  summarise(
    y_min = min(y) -0.03,
    y_max = max(y) +0.03, # line position
    label_x = 0.616, # label position
    label_y = (y_min + y_max) / 2,
    .groups = "drop"
  ) %>%
  filter(!str_detect(midas_species, "^\\*s"),
         !str_detect(midas_species, "^s:Thauera"),
         !str_detect(midas_species, "s:midas_s_256;"))

col <- c(
  'transparent','transparent','transparent',
  
  "#8db7d2","#5e62a9","#434279",
  "#A3B18A","#9FB8AD",         "#c45161","#ff6361","#e094a0", "#b9b5c9", "#f2b6c0",
  
  "#0B3C5D", "#2A9D8F", "#5E60CE", "#144552", "#3AAFB9", "#1B2A41", "#72B5A4", "#B5E2DE","#8dB8AD",  "#4E8077", 
  
         "#8B1E3F", "#C1440E", "#A63F03", "#D95D39", "#B23A48", "#802420", "#E07A5F", "#F4A261", 
         
         "#95825A","#D18E35",'#ddcea1', "#84593A","#D9A066","#E0C37E", "#C6A07E",
  "#f2dfa2","#B7A27C",  "#9c7605",
         "#EAC5a1", '#c6b46f',"#EAC5AA","#BFA5A0" ) 

# 16 midas species

s <- ppp + geom_segment(aes(x=0.321, y=y, xend=if_else(xlim> varx, varx-0.0009, xlim), yend=ylim, 
                            color = if_else(is.na(proposed_name_new), 'g__Thauera',proposed_name_new)),
  data=dd,
  alpha = 0.75, linewidth=1.3,
  lineend = 'round')+ 
  geom_segment(data = species_segments,
               aes(x = label_x, xend = label_x, y = y_min, yend = y_max, color= midas_species),
               inherit.aes = FALSE, size = 0.7)+
  geom_text(data = species_segments,
            aes(x = label_x+0.003, y = label_y, label = str_remove_all(midas_species, "s:|;")),
            inherit.aes = FALSE, hjust = 0, size = 3) +
  scale_color_manual('', values = col) + 
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(0, 0.7), clip = "off") 

 
ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS3compare_tree_hq_propnames.jpeg',
       s,
       height = 11,
       width=16.5
       )

ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS3compare_tree_hq_propnames.pdf',
       s,
       height = 11,
       width=16.5
)







