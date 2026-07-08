.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")
library(ampvis2)
library(data.table)
library(tidyverse)
library(patchwork)

#load data (from MiDAS 4 paper https://doi.org/10.1038/s41467-022-29438-7)
seq_metadata <- read.csv("/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/Dueholm2021/2020-07-03_Sequencing_metadata.txt", sep="\t")
wwtp_m <- read.csv("/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/Dueholm2021/DataS1_210413.txt", sep="\t")
wwtp_m$Mean <- "mean"
wwtp_m$abs_Latitude <- abs(wwtp_m$Latitude)

#V13
seq_V13 <- seq_metadata[ , which(names(seq_metadata) %in% c("V13_seq_id","WWTP_id"))]

V13metadata <- merge.data.frame(seq_V13, wwtp_m, by ="WWTP_id", all.x = TRUE)
V13metadata$WWTP_ID <- V13metadata$WWTP_id
V13metadata <- V13metadata[,-1]

#load data (from MiDAS 4 paper https://doi.org/10.1038/s41467-022-29438-7)
d13 <- amp_load("/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/Dueholm2021/V13_ASVtab.txt", 
                taxonomy = "/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/raw/Dueholm2021/V13ASV_vs_MiDAS_4.8.sintax",
                metadata = V13metadata)

# Remove samples with low read count
d13n <- amp_subset_samples(d13, minreads = 10000, normalise = TRUE)

#Subset AS only in 4 basic process_type(s)
d13nAS <- amp_subset_samples(d13n, Plant_type == "Activated sludge")
d13nAS4PT <-  amp_subset_samples(d13nAS, Process_type %in% c("C", "C,N", "C,N,DN", "C,N,DN,P"))  

colnames(d13nAS4PT$metadata)
d13nAS4PT$metadata %>% pull(WWTP_ID) %>% n_distinct()
d13nAS4PT$metadata %>% pull(Country) %>% n_distinct()
d13nAS4PT$metadata %>% pull(V13_seq_id) %>% n_distinct()

meta <- d13nAS4PT$metadata

# heatmap of top genera
col <- c("#f4e7d9","#80CDC1","#4E8077","#c24e44", "#990000")

heat_genus <- amp_heatmap(d13nAS4PT,
            group_by = c("Mean"),
            tax_aggregate = "Genus",
            measure = "mean",
            tax_show = 25,
            normalise = FALSE,
            textmap = T,
            tax_empty = 'best',
            )

rownames(heat_genus) <- sapply(rownames(heat_genus), function(x) {
  if (grepl("^Ca_", x)) {
    rest <- sub("^Ca_", "", x)
    paste0("italic('Ca.')~'", rest, "'")
  } else if (!grepl("[0-9]", x)) {
    paste0("italic('", x, "')")
  } else {
    paste0("'", x, "'")
  }
})

heat_genus$Genus <- rownames(heat_genus)
heat_genus$Genus <- factor(heat_genus$Genus, levels = rownames(heat_genus))
melted_genus_df <- reshape2::melt(heat_genus, id.vars = "Genus")

# Create the bins and corresponding colors
bins <- c(0, 0.01, 0.1, 1, 10, Inf)
melted_genus_df$bin <- cut(melted_genus_df$value, breaks = bins, labels = FALSE, right = FALSE)

plotg <- ggplot(melted_genus_df, aes(x='Global mean', y = Genus, fill = as.factor(bin))) +
  geom_tile() +
  geom_text(aes(label = round(value, digits = 2)),
            size = 2.5) +
 scale_y_discrete(limits = rev, drop = T, labels = scales::parse_format()) +
   scale_fill_manual(values = col, 
                    breaks = c(1,2,3,4,5), 
                    labels = c("0-0.01", "0.01-0.1", "0.1-1", "1-10", ">10")) +
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks.x = element_blank(),
    legend.position = 'none'
   ) 


# Species distribution of Zoogloea pr continent
d13_zoogloea <- amp_subset_taxa(d13nAS4PT, tax_vector = 'g__Zoogloea')

# textmap
gV4 <- amp_heatmap(d13_zoogloea,
            group_by = c("Continent"),
            tax_aggregate = "Species",
            measure = "mean",
            tax_show = 7000,
            normalise = FALSE,
            textmap = T,
            tax_empty = 'best'
            )


# --- 1. Extract ASV rows ---
unclassified <- gV4[grep("ASV", rownames(gV4)), ]

# --- 2. Remove ASV rows from gV4 ---
gV4_noASV <- gV4[-grep("ASV", rownames(gV4)), ]

# --- 3. Filter rows with any value > 0.01 ---
gV4_filtered <- gV4_noASV %>% filter_all(any_vars(. > 0.01))
gV4_filtered_continent <- gV4_noASV %>% filter_all(any_vars(. > 0.01))

# --- Identify rows that were filtered out ---
gV4_not_filtered <- gV4_noASV[!rownames(gV4_noASV) %in% rownames(gV4_filtered), ]

# --- Compute sum of filtered-out species for "Other" ---
other_row <- colSums(gV4_not_filtered)

# --- 5. Compute sum of all original rows for "Zoogloea (sum)" ---
zoogloea_sum <- colSums(gV4)

# --- 6. Compute sum of ASV rows for "Unclassified" ---
unclassified_sum <- colSums(unclassified)

# --- 7. Bind rows together ---
new <- rbind(
  'Zoogloea (sum)' = zoogloea_sum,
  gV4_filtered,
  'Other' = other_row,
  'Unclassified' = unclassified_sum
)

#new <- rbind('Zoogloea (sum)' =colSums(gV4),gV4,'Unclassified' =colSums(unclassified))
rownames(new) <- sapply(rownames(new), function(x) {
  if (x == "Zoogloea (sum)") {
    "italic('Zoogloea')*' (sum)'"  # preserves the space and parentheses
  } else if (grepl("^Zoogloea_", x)) {
    species <- gsub("^Zoogloea_", "", x)
    paste0("italic('Zoogloea')*' '*italic('", species, "')")
  } else {
    paste0("'", x, "'")  # keep other names plain
  }
})



# number of WWTPs in every continent
d13nAS4PT$metadata %>% select(WWTP_ID, Continent) %>% unique() %>% count(Continent)
 
# Create the bins and corresponding colors
bins <- c(0, 0.01, 0.1, 1, 10, Inf)
#colors <- c("gray90","#FFF5EB", "#F5C9B0", "#EC9E75", "#E2733B")

new$Genus <- rownames(new)
new$Genus <- factor(new$Genus, levels = rownames(new))

melted_df <- reshape2::melt(new, id.vars = "Genus")
melted_df$bin <- cut(melted_df$value, breaks = bins, labels = FALSE, right = FALSE)




# Reorder the levels of Genus to introduce a gap
melted_df$Genus <- factor(melted_df$Genus, 
                          levels = c(levels(melted_df$Genus)[1],  # Top row
                                     " ",  # Spacer
                                     levels(melted_df$Genus)[-1]))  # Rest

# Plot with modified Genus
plot <- ggplot(melted_df, aes(x = variable, y = Genus, fill = as.factor(bin))) +
  geom_tile() + theme_minimal()+
  geom_text(aes(label = round(value, digits = 2)), size = 2.5) +
  labs(y = "Species", fill = "% relative\n abundance", x = "") +
  scale_y_discrete(limits = rev, drop = FALSE, # Keep the blank space
                   labels = scales::parse_format()) + 
  scale_fill_manual(values = col, 
                    breaks = c(1,2,3,4,5), 
                    labels = c("0-0.01", "0.01-0.1", "0.1-1", "1-10", ">10")) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank(), panel.grid = element_blank()
  ) +
  guides(fill = guide_legend(reverse = TRUE,
                             override.aes = list(color = "black", size = 0.5)))





# process design overall
gV4 <- amp_heatmap(d13_zoogloea,
            group_by = c("Process_type"),
            tax_aggregate = "Species",
            measure = "mean",
            tax_show = 7000,
            normalise = FALSE,
            textmap = T,
            tax_empty = 'best'
            )


# --- 1. Extract ASV rows ---
unclassified <- gV4[grep("ASV", rownames(gV4)), ]

# --- 2. Remove ASV rows from gV4 ---
gV4_noASV <- gV4[-grep("ASV", rownames(gV4)), ]

# --- 3. Filter rows with any value > 0.01 ---
gV4_filtered <- gV4_noASV[!rownames(gV4_noASV) %in% rownames(gV4_filtered_continent), ]
# --- Identify rows that were filtered out ---
gV4_not_filtered <- gV4_noASV[rownames(gV4_filtered_continent), ]

# --- Compute sum of filtered-out species for "Other" ---
other_row <- colSums(gV4_filtered)

# --- 5. Compute sum of all original rows for "Zoogloea (sum)" ---
zoogloea_sum <- colSums(gV4)

# --- 6. Compute sum of ASV rows for "Unclassified" ---
unclassified_sum <- colSums(unclassified)

# --- 7. Bind rows together ---
new <- rbind(
  'Zoogloea (sum)' = zoogloea_sum,
  gV4_not_filtered,
  'Other' = other_row,
  'Unclassified' = unclassified_sum
)



#new <- rbind('Zoogloea (sum)' =colSums(gV4),gV4,'Unclassified' =colSums(unclassified))
rownames(new) <- sapply(rownames(new), function(x) {
  if (x == "Zoogloea (sum)") {
    "italic('Zoogloea')*' (sum)'"  # preserves the space and parentheses
  } else if (grepl("^Zoogloea_", x)) {
    species <- gsub("^Zoogloea_", "", x)
    paste0("italic('Zoogloea')*' '*italic('", species, "')")
  } else {
    paste0("'", x, "'")  # keep other names plain
  }
})



tokeep <- plot$data$Genus %>% unique()



# Create the bins and corresponding colors
bins <- c(0, 0.01, 0.1, 1, 10, Inf)

new$Genus <- rownames(new)
new$Genus <- factor(new$Genus, levels = rownames(new))
melted_df <- reshape2::melt(new, id.vars = "Genus") 

melted_df$bin <- cut(melted_df$value, breaks = bins, labels = FALSE, right = FALSE)

# Reorder the levels of Genus to introduce a gap

melted_df <- melted_df %>% filter(Genus %in% tokeep)
melted_df$Genus <- factor(melted_df$Genus, 
                          levels = c(levels(tokeep)))  # Rest tokeep

# Plot with modified Genus
plotproces <- ggplot(melted_df, aes(x = variable, y = Genus, fill = as.factor(bin))) +
  geom_tile() + theme_minimal()+
  geom_text(aes(label = round(value, digits = 2)), size = 2.5) +
  labs(y = "Species", fill = "% relative\n abundance", x = "") +
  scale_y_discrete(limits = rev, drop = FALSE, labels = scales::parse_format()) +  # Keep the blank space
  scale_fill_manual(values = col, 
                    breaks = c(1,2,3,4,5), 
                    labels = c("0-0.01", "0.01-0.1", "0.1-1", "1-10", ">10")) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank(), panel.grid = element_blank()
  ) +
  guides(fill = guide_legend(reverse = TRUE,
                             override.aes = list(color = "black", size = 0.5)))










d_exp_zoogloea <- amp_export_long(d13_zoogloea)

d_exp_caeni <- d_exp_zoogloea %>% filter(Species == 's__Zoogloea_caeni') %>% 
  group_by(V13_seq_id) %>% mutate(zoo_caeni = sum(count))

d_exp_caeni_new <- d_exp_caeni %>% ungroup()  %>% select(V13_seq_id,Continent, Country, Latitude, Longitude, City, WWTP_ID, zoo_caeni) %>% unique() %>% 
  group_by(Country) %>% mutate(zoo_caeni_country = mean(zoo_caeni))

d_exp_caeni_country <- d_exp_caeni_new %>% select(Continent, Country, Latitude, Longitude, zoo_caeni_country) %>% unique()

library(maps)

mp <- NULL
mapWorld <- borders("world", colour = "#333333", fill="#edeadd", size=0.1)
mp <- ggplot() + mapWorld



d_exp_caeni <- d_exp_zoogloea %>% #filter(Species == 's__Zoogloea_caeni') %>% 
  group_by(V13_seq_id) %>% mutate(zoo_sum = sum(count))

d_exp_caeni_new <- d_exp_caeni %>% ungroup()  %>% select(V13_seq_id,Continent, Country, Latitude, Longitude, City, WWTP_ID, zoo_sum) %>% unique() %>% 
  group_by(Country) %>% mutate(zoo_caeni_country = mean(zoo_sum))

dfplot <- d_exp_caeni_new %>% group_by(WWTP_ID) %>% mutate(zoo_wwtp = mean(zoo_sum)) %>% select(-zoo_sum,-V13_seq_id) %>% unique()


bins <- c(0, 0.01, 0.1, 1, 10, Inf)
colors <- c("gray90","#FFF5EB", "#F5C9B0", "#EC9E75", "#E2733B")

dfplot$bin <- cut(dfplot$zoo_wwtp, breaks = bins, labels = FALSE, right = FALSE)

#col <- c("#4E8077","#80CDC1","#f4e7d9","#c24e44", "#990000")
col <- c("#f4e7d9","#80CDC1","#4E8077","#c24e44", "#990000")
mp1 <- mp + geom_point(data = dfplot,aes(x=Longitude, y=Latitude, fill = as.factor(bin)),
                       shape =21, alpha=0.95,
                       size=1.6) +
  scale_fill_manual('Relative\nabundance (%)',
    values = col, 
                    breaks = c(1,2,3,4,5), 
                    labels = c("0-0.01", "0.01-0.1", "0.1-1", "1-10", ">10")) +
  guides(fill = guide_legend( nrow=1,
                              byrow = TRUE,
    # keywidth  = unit(0.35, "cm"),
    # keyheight = unit(0.35, "cm"),
    default.unit = "cm",
    override.aes = list(shape = 22, size = 6)  # shape 22 is a square
  )) +
  theme_minimal()+
  theme(axis.title = element_blank(),
        legend.position = 'bottom',
        legend.title = element_text(face='bold'),
        legend.spacing.x = unit(0.1, 'cm'),   # horizontal spacing between keys
    legend.spacing.y = unit(0.001, 'cm'),   # vertical spacing between rows
    legend.key.size = unit(0.2, 'cm'),    # size of the squares
    legend.margin = margin(t = 1, r = 1, b = 1, l = 1),
    legend.box.spacing = unit(0.5, "cm"),
    legend.text = element_text(margin = margin(l = -0.5)))+ 
  coord_quickmap(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE)




#play 
# create ordered WWTP_ID levels grouped by continent + country
library(grid)
# poster
gV4 <- amp_heatmap(d13_zoogloea,
            group_by = c("WWTP_ID"),
            tax_aggregate = "Species",
            measure = "mean",
            tax_show = 7000,
            normalise = FALSE,
            textmap = T,
            tax_empty = 'best'
            )

# new <- gV4  %>% # textmap from ampvis
#    filter_all(.,any_vars(. >0.5))

# --- 1. Extract ASV rows ---
unclassified <- gV4[grep("ASV", rownames(gV4)), ]

# --- 2. Remove ASV rows from gV4 ---
gV4_noASV <- gV4[-grep("ASV", rownames(gV4)), ]

# --- 3. Filter rows with any value > 0.01 ---
gV4_filtered <- gV4_noASV %>% filter_all(any_vars(. > 0.1))

# --- Identify rows that were filtered out ---
gV4_not_filtered <- gV4_noASV[!rownames(gV4_noASV) %in% rownames(gV4_filtered), ]

# --- Compute sum of filtered-out species for "Other" ---
other_row <- colSums(gV4_not_filtered)

# --- 5. Compute sum of all original rows for "Zoogloea (sum)" ---
zoogloea_sum <- colSums(gV4)

# --- 6. Compute sum of ASV rows for "Unclassified" ---
unclassified_sum <- colSums(unclassified)

# --- 7. Bind rows together ---
new <- rbind(
  gV4_filtered,
  'Other' = other_row,
  'Unclassified' = unclassified_sum
)

#new <- rbind('Zoogloea (sum)' =colSums(gV4),gV4,'Unclassified' =colSums(unclassified))
rownames(new) <- sapply(rownames(new), function(x) {
  if (x == "Zoogloea (sum)") {
    "italic('Zoogloea')*' (sum)'"  # preserves the space and parentheses
  } else if (grepl("^Zoogloea_", x)) {
    species <- gsub("^Zoogloea_", "", x)
    paste0("italic('Zoogloea')*' '*italic('", species, "')")
  } else {
    paste0("'", x, "'")  # keep other names plain
  }
})




# Create the bins and corresponding colors
bins <- c(0, 0.01, 0.1, 1, 10, Inf)

new$Genus <- rownames(new)
new$Genus <- factor(new$Genus, levels = rownames(new))
melted_df <- reshape2::melt(new, id.vars = "Genus")

melted_df$bin <- cut(melted_df$value, breaks = bins, labels = FALSE, right = FALSE)

meta_wwtpid <- meta %>% select(-V13_seq_id) %>% unique()
melted_df_meta <- melted_df %>% left_join(., meta_wwtpid, by = c('variable'='WWTP_ID')) %>% 
  mutate(country_tl = substr(variable, 1, 2)) 

melted_df_meta$Continent = factor(melted_df_meta$Continent, levels=c("North America" ,"South America","Europe" , "Africa"     ,         "Asia"   , "Oceania"        ))

# --- Order country facets by continent ---
melted_df_meta <- melted_df_meta %>%
  mutate(country_tl = factor(country_tl,
                             levels = melted_df_meta %>%
                               distinct(country_tl, Continent) %>%
                               arrange(Continent, country_tl) %>%
                               pull(country_tl)))


library(ggh4x)
# --- Base plot ---
plotbiogeo <- ggplot(melted_df_meta,
       aes(x = variable, y = Genus, fill = as.factor(bin))) +
  geom_tile(color = 'white',
    linewidth = 0.1
    ) +
  
  ggh4x::facet_nested( ~ Continent+country_tl, scales = "free", space = "free", 
                       solo_line = T,
               nest_line = element_line(color="grey20", linewidth = 0.4), resect=unit(2, "pt"),
               , strip = ggh4x::strip_nested(clip = "off", 
                                             text_x = elem_list_text(angle=c(rep(0,6), rep(90,30))))) +
  
  labs(y = "Species", fill = "% relative\n abundance", x = "") +
  scale_y_discrete(limits = rev, drop = FALSE, labels = scales::parse_format(), expand = c(0,0)) +
  # scale_y_discrete(limits = rev, drop = TRUE, expand = c(0,0)) +
  scale_fill_manual(values = col,
                    breaks = c(1,2,3,4,5),
                    labels = c("0-0.01", "0.01-0.1", "0.1-1", "1-10", ">10")) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = TRUE,
                             override.aes = list(color = "black", size = 0.5))) +
  coord_cartesian(clip = "off") +
  
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    strip.clip = 'off',
    strip.text = element_text(face = "bold", size = 8,),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(t = 3, r = 5, b = 5, l = 5),
    panel.spacing = unit(0.1, "lines", data = NULL)
  )



complot <- cowplot::plot_grid(plotg +theme(text = element_text(color = 'black'),
                                           axis.text.y = element_text(size=8)),
                              
                              mp1+theme(text = element_text(color = 'black')),
                              nrow = 1, rel_widths = c(1.65,7),
                              labels = c("A", "B"), label_size = 14, label_fontface = "plain")
heats <- cowplot::plot_grid(plot+theme(legend.position = 'none',
                                       text = element_text(color = 'black'))+ labs(tag = 'D'),
                              plotproces+theme(legend.position = 'none', axis.text.y=element_blank(),
                                               text = element_text(color = 'black'))+ labs(tag = 'E') ,
                            align = 'h', nrow=1, rel_widths = c(1,0.45))
complot1 <- cowplot::plot_grid(complot,
                               heats,
                            
                              nrow = 1, rel_widths = c(2,1.15))

compf <- cowplot::plot_grid(complot1+theme(text = element_text(color = 'black')),
                           plotbiogeo + labs(tag = 'C')+theme(legend.position = 'none'),
                           nrow = 2, rel_heights = c(1,0.6))





ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/Fig1_global_abundance_all_high_switch.jpeg',
       compf, width = 12.5, height = 8, dpi=800,
       bg = 'white')
ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/Fig1_global_abundance_all_high_switch.svg',
       compf, width = 12.5, height = 8,
       bg = 'white')





## other systems than AS
#Subset AS only in 4 basic process_type(s)
d13nOther <- amp_subset_samples(d13n, !Plant_type %in% c("Activated sludge",""))

# Count samples per Plant_type (unique WWTP_IDs)
plant_counts <- d13nOther$metadata %>%
  ungroup() %>%
  select(WWTP_ID, Plant_type) %>%
  distinct() %>%
  count(Plant_type)

# Build labels with n
new_labels <- paste0(plant_counts$Plant_type, " (n=", plant_counts$n, ")")
names(new_labels) <- plant_counts$Plant_type

# Heatmap data
heat_genus <- amp_heatmap(
  d13nOther,
  group_by = c("Plant_type"),
  tax_aggregate = "Genus",
  measure = "mean",
  tax_show = 25,
  normalise = FALSE,
  textmap = TRUE,
  tax_empty = "best"
)

heat_genus$Genus <- rownames(heat_genus)
heat_genus$Genus <- factor(heat_genus$Genus, levels = rownames(heat_genus))

melted_df <- reshape2::melt(heat_genus, id.vars = "Genus")

# Apply new x-axis labels
melted_df$variable <- factor(
  melted_df$variable,
  levels = names(new_labels),
  labels = new_labels
)

# Plot
plot11 <- ggplot(melted_df, aes(x = variable, y = Genus, fill = value)) +
  geom_tile() +
  theme_minimal() +
  geom_text(aes(label = round(value, 2)), size = 2.25) +
  labs(y = "Species", fill = "% relative\n abundance", x = "") +
  scale_y_discrete(limits = rev, drop = FALSE) +
  scale_fill_gradient(low = "#FFF5EE", high = "#E2733B") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = 'black'),
    axis.text.y= element_text(color = 'black'),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )


ggsave('/home/bio.aau.dk/kl42gg/projects/rethink/zoogloea_paper/data/plots/supplementary/FigS1_other_systems_global_abundance.jpeg',
       dpi=600,
       plot11, width = 5, height = 5)

