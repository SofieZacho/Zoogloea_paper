.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(data.table)
library(treeio)
library(gggenes)

# load metadata with proposed names for the genomes
propnames <- read_tsv('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/metadata/proposed_R226_new_names_updated.tsv', col_names = T)  
propnames[] <- lapply(propnames, function(x) if(is.factor(x) || is.character(x)) gsub("Zoogloea", "Z.", as.character(x)) else x)
propnames <- propnames %>% rename(midas_species_label = midas_species)



# load hits from annotation using DRAM as described in the methods section of the manuscript
hits <- read_tsv("/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/generated/annotation/KEGG_zoogloea_25_01_24_HQ/00_condensed_HQ_dram.tsv") %>%
    filter(fasta != 'fasta')
hits$fasta <- str_remove_all(hits$fasta, "_ASM.*")
hits$fasta <- str_remove_all(hits$fasta, "_SRR.*")
hits$fasta <- str_remove_all(hits$fasta, "_TT.*")


# load KO IDs of interest can be found in the Supplementary Data File 3 of the manuscript
KOs <- readxl::read_xlsx('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/supplementary/KOs_SZV.xlsx',
                         skip = 2) %>% 
  filter(!Gene %in% c('nxrA', 'nxrB')) %>% 
  filter(!Gene %in% 'sat' | !Metabolism %in% 'DSR') %>% 
  filter(!Metabolism %in% c('ASR','DSR', 'Fructose'), # as none Zoogloea has this
      !Gene %in% c('ugpA', 'ugpB','ugpE',
                   'hycE','hyfB','hyfC','hyfE', 'hyfF', 'mbhJ',
                   'glgB other',
                  'amoA', 'amoB','amoC','hao',
                  'sor','sorA','sorB','SUOX',
                   #'pta', #acetate 
                   'porA','porB' # pyruvate something - not detected in any
                   ))




# subset to wanted KOs
hits_sub <- hits %>% filter(ko_id %in% unique(KOs$KO))
combinations <- expand.grid(unique(hits_sub$fasta), unique(KOs$KO)) %>%
  rename('fasta'=1,'ko_id'=2)


# number of found genes
n_distinct(hits_sub$ko_id)

# count number of hits for each KO
hits_sub_count <- hits_sub %>% 
  group_by(fasta) %>% 
  count(ko_id) 


com <- combinations %>% left_join(hits_sub_count) %>%
  left_join(KOs, by = c('ko_id'= 'KO')) %>%
  #left_join(meta, by = c('fasta' = 'tip_label')) %>%
  rename(tip_label=fasta)

# geom_tile
tile <- com %>% select(tip_label,ko_id,Gene, n,Overall_metabolism, Metabolism) %>% 
  arrange(Overall_metabolism) %>% 
  mutate(Gene = factor(Gene, levels = unique(Gene)))


tile <- tile %>% left_join(., propnames, by = c('tip_label'='genome'))


tile <- tile %>% 
  filter(!is.na(proposed_name_new)) %>% 
    mutate(tip_label = factor(x= tip_label, levels = propnames$genome),
           proposed_name_new = factor(x= proposed_name_new, levels = unique(propnames$proposed_name_new)))
tile$tip_label <- gsub("GCA_009026175.1","GCF_009026175.1", tile$tip_label)

tile$tip_label <- gsub('barcode02_Z_caeni','GCA_986281895*', tile$tip_label)
tile$tip_label <- gsub('barcode03_Z_oleivorans','GCA_986281795*', tile$tip_label)
tile$tip_label <- gsub('barcode04_Z_resiniphila','GCA_986340685*', tile$tip_label)
  

metcol <- c("Carbon sources and processing"     =    "#EAC5AA",
                 
                 "Nitrogen cycle" ="#34978F"  ,       
                  "PHA, glycogen & amino acid storage"     ='#B7A27C', 'Storage polymers'=    '#B7A27C',
               "AA storage"       ='#95825A', 
                 "Aromatic aa"  ='#9c7605',    "Cyanophycin"   ='#84541E',    "Polyphosphate"   ='#E0C37E',  
                 "ASR" ='#cf6c42','DSR'='#dead9f', 
                 'SOX'= '#d98c6b', 'Sulfur metabolism' = '#d98c6b'
  )
  
## subset of genes
tile.subset <- tile %>% filter(Metabolism %in% c('N fixation',
                                                 'Denitrification','ANR','DNRA',
                                                 'PHA','Glycogen','Cyanophycin',
                                                 'Polyphosphate','Sulfur oxidation'),
                               !Gene %in% c('nasA','nasB', 'narB','nirA','anfG',
                                            'glgX')) %>% 
  mutate(Overall_metabolism = if_else(Overall_metabolism %in% 'PHA, glycogen & amino acid storage',
                                      'Storage polymers',
                                      Overall_metabolism))


# metabolism
p_tile <- ggplot(data=tile.subset)+ 
  geom_tile(
    mapping = aes(
      y=tip_label, x=Gene, fill = ifelse(is.na(n), NA, Overall_metabolism)), 
    size = 0.03, color = 'grey10'
  )+ 
  scale_fill_manual('Metabolism',values = metcol, na.value = "grey98", breaks = unique(tile.subset$Overall_metabolism))+
  geom_text( mapping = aes(y=tip_label, x=Gene, label = ifelse(is.na(n), '', n)),
             size= 2.1
  )  +
  facet_grid(cols = vars(Overall_metabolism), rows = vars(proposed_name_new),
             space = 'free', scales = 'free')+
  theme_minimal()+
  theme(
    legend.background = element_rect(fill = 'transparent'),
    axis.title = element_blank(),
    axis.text.x = element_text(size=7.5, angle =90, hjust=1, vjust = 0.5, face = 'italic', color = 'black'),
    axis.text.y= element_text(size=7, color = 'black'), 
    strip.text.y = element_text(angle=0, hjust=0, face = 'bold', size = 7, color = 'black'),
    strip.clip = 'off',
    strip.text.x = element_text(size=10, face = 'bold', color = 'black'),
    panel.spacing.y=unit(0.045, "lines"),
    panel.spacing.x=unit(0.7, "lines"),
    panel.grid = element_blank(),
    legend.position = 'none'
  ) 


ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/Fig3_heat_simple_subset_propnames.tiff',
       p_tile,dpi=1200,
       compression = "lzw",
       height = 7.5,
       width = 9.5)

ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/Fig3_heat_simple_subset_propnames.jpeg',
       p_tile,dpi=600,
       height = 7.5,
       width = 9.5)



### supplementary Figures S4 and S5


# column with continent. Metadata from MiDAS genome paper
geo_info <- read_tsv('/projects/PHN/MiDAS/04-Analysis/Status-analysis/103-samples/100-add_checkm1/03-MAGs.final-version/000-MAGs-metadata/geo.info/geo.info.tsv', 
                     col_names = c('plant.ID', 'continent', 'country', 'city'))


tile.subset.geo <- tile.subset %>% 
  mutate(plant.ID = str_extract(tip_label, "[A-Z]+\\d*-\\d+")) %>% 
  left_join(., geo_info) %>% 
  filter(!is.na(continent))


# metabolism
p_tile_continent <- ggplot(data=tile.subset.geo)+ 
  geom_tile(
    mapping = aes(
      y=paste0("italic('", proposed_name_new, "')~' (", tip_label, ")'"),
      x=Gene, fill = ifelse(is.na(n), NA, Overall_metabolism)), 
    size = 0.03, color = 'grey10'
  )+ 
  scale_fill_manual('Metabolism',values = metcol, na.value = "grey95", breaks = unique(tile.subset$Overall_metabolism))+
  geom_text( mapping = aes(y=paste0("italic('", proposed_name_new, "')~' (", tip_label, ")'"),
                           x=Gene, label = ifelse(is.na(n), '', n)),
             size= 2.1
  )  +
  facet_grid(cols = vars(Overall_metabolism), rows = vars(continent),
             space = 'free', scales = 'free')+
  theme_minimal()+
  theme(
    legend.background = element_rect(fill = 'transparent'),
    axis.title = element_blank(),
    axis.text.x = element_text(size=7.5, angle =90, hjust=1, vjust = 0.5, face = 'italic', color = 'black'),
    axis.text.y= element_text(size=7, color = 'black'), 
    strip.text.y = element_text(angle=0, hjust=0, face = 'bold', size = 8, color = 'black'),
    strip.clip = 'off',
    strip.text.x = element_text(size=10, face = 'bold', color = 'black'),
    panel.spacing.y=unit(0.045, "lines"),
    panel.spacing.x=unit(0.7, "lines"),
    panel.grid = element_blank(),
    legend.position = 'none'
    #strip.background = element_rect(fill = 'grey90')
  ) +
  scale_y_discrete(labels = function(x) parse(text = x))





tile.subset.geo.m2 <- tile.subset.geo %>% group_by(plant.ID) %>% 
  filter(n_distinct(tip_label) > 1)

p_tile_country_continent <-ggplot(data=tile.subset.geo.m2)+ 
  geom_tile(
    mapping = aes(
      y=paste0("italic('", proposed_name_new, "')~' (", tip_label, ")'"),
      x=Gene, fill = ifelse(is.na(n), NA, Overall_metabolism)), 
    size = 0.03, color = 'grey10'
  )+ 
  scale_fill_manual('Metabolism',values = metcol, na.value = "grey95", breaks = unique(tile.subset$Overall_metabolism))+
  geom_text( mapping = aes(y=paste0("italic('", proposed_name_new, "')~' (", tip_label, ")'"),
                           x=Gene, label = ifelse(is.na(n), '', n)),
             size= 2.1
  )  +
  facet_grid(cols = vars(Overall_metabolism), rows = vars(continent,country, plant.ID),
             space = 'free', scales = 'free')+
  theme_minimal()+
  theme(
    legend.background = element_rect(fill = 'transparent'),
    axis.title = element_blank(),
    axis.text.x = element_text(size=7.5, angle =90, hjust=1, vjust = 0.5, face = 'italic', color = 'black'),
    axis.text.y= element_text(size=7, color = 'black'), 
    strip.text.y = element_text(angle=0, hjust=0, face = 'bold', size = 8, color = 'black'),
    strip.clip = 'off',
    strip.text.x = element_text(size=10, face = 'bold', color = 'black'),
    panel.spacing.y=unit(0.045, "lines"),
    panel.spacing.x=unit(0.7, "lines"),
    panel.grid = element_blank(),
    legend.position = 'none'
    #strip.background = element_rect(fill = 'grey90')
  ) +
  scale_y_discrete(labels = function(x) parse(text = x))

ptotcon <- cowplot::plot_grid(p_tile_continent + labs(tag= 'A') + theme(plot.tag = element_text(face='plain') , axis.text.x = element_blank()),
                              p_tile_country_continent+ labs(tag= 'B') + theme(plot.tag = element_text(face='plain'), strip.text.x =  element_blank()),
                              ncol=1,
                              rel_heights = c(1,0.6),
                              align = 'v')
ptotcon

ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS5_heat_simple_subset_geo_combined.jpeg',
       ptotcon, dpi = 600,
       height = 10,
       width = 10)





#### arrow plot of nitrogen genes in Z. caeni. Genes of interest are described in the Supplementary Data File 3 of the manuscript. The genes are annotated using DRAM as described in the methods section of the manuscript.
KOs_nitrogen <- readxl::read_xlsx('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/supplementary/KOs_SZV.xlsx', skip=2)  %>% 
  filter(Overall_metabolism %in% 'Nitrogen cycle',!KO %in% c('K00371','K00370') 
         )%>% 
  filter(!Description %in% 'nrfA; nitrite reductase (cytochrome c-552) [EC:1.7.2.2]')


# subset to wanted KOs
hits_sub <- hits %>% filter(ko_id %in% unique(KOs_nitrogen$KO))%>%
  left_join(KOs_nitrogen, by = c('ko_id'= 'KO')) 


# 1️⃣ Prepare genes
genes <- hits_sub %>%
  filter(fasta == "barcode02_Z_caeni") %>%
  arrange(scaffold, start_position) %>%
  mutate(strandedness = strandedness == 1,
         scaffold_num = as.numeric(as.factor(scaffold)))

# 2️⃣ Identify gene blocks (>10 kb gaps)
genes <- genes %>%
  group_by(scaffold_num) %>%
  mutate(
    gap = start_position - lag(end_position, default = first(start_position)),
    new_block = if_else(gap > 10000, 1, 0),
    block_id = cumsum(new_block) + 1
  ) %>%
  ungroup()

# 3️⃣ Compute compressed positions
block_offsets <- genes %>%
  group_by(block_id) %>%
  summarise(
    block_start = min(start_position),
    block_end = max(end_position),
    block_len = block_end - block_start,
    .groups = "drop"
  ) %>%
  mutate(offset = lag(cumsum(block_len + 1000), default = 0))

genes <- genes %>%
  left_join(block_offsets, by = "block_id") %>%
  mutate(
    start_compressed = start_position - block_start + offset,
    end_compressed = end_position - block_start + offset,
    pos_compressed = (start_compressed + end_compressed)/2
  )

# 4️⃣ Prepare "slash" labels at block boundaries
slash_labels <- block_offsets %>%
  filter(row_number() != 1) %>%
  mutate(
    x = offset -500,                  
    y = max(genes$scaffold_num) + 0, 
    label = "//"
  )

# 5️⃣ Prepare x-axis ticks with real genomic positions
axis_labels <- block_offsets %>%
  rowwise() %>%
  mutate(
    ticks = list(seq(block_start, block_end, length.out = 3))  # 3 ticks per block
  ) %>%
  unnest(ticks) %>%
  mutate(
    pos_compressed = ticks - block_start + offset,
    label = ticks
  )

cols <-c( "#34978F" ,"#45A59E" ,"#57B3AD", "#69C1BC", "#7ACCCA" ,"#2E8780" ,"#29726D", "#235F5A" ,"#41a39a" ,"#63B0AA", "#28907D" ,"#287A72")

# 6️⃣ Plot
supplot <- ggplot(genes, aes(xmin = start_compressed, xmax = end_compressed, y = scaffold_num, fill = Gene)) +
  geom_gene_arrow(aes(forward = strandedness),
                  alpha = 1,
                  arrowhead_height = grid::unit(6,"mm"),
                  arrow_body_height = grid::unit(6,"mm"),
                  arrowhead_width = grid::unit(1.5,"mm"),
                  linetype = "blank") +
  geom_text(aes(x = pos_compressed, label = Gene),
            size = 2.5, angle = 30, hjust = 0.5, fontface = "bold") +
  geom_text(data = slash_labels,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 5.5, fontface = "bold") +
  scale_x_continuous(
    breaks = axis_labels$pos_compressed,
    labels = axis_labels$label
  ) +
  scale_y_continuous(
    breaks = genes$scaffold_num,
    labels = genes$scaffold
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 8, angle=45, hjust=1, vjust=1, color = 'black'),
    axis.text.y = element_blank(),
    legend.position = "none", 
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = 'black', linewidth = 0.2)
  ) +
  scale_fill_manual(values = cols)


ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS4_sup_nitrogen_genes.jpeg',
       supplot, dpi=600,
       height = 1.3,
       width = 10)




