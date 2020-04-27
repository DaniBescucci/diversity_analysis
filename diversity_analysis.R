
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(PMCMRplus)
library(here)

observed_otus <- read_tsv('tourmaline/04-diversity/dada2-pe/observed_otus/alpha-diversity.tsv')

observed_otus <- dplyr::rename(observed_otus, `sample-id` = X1)

shannon <- read_tsv('tourmaline/04-diversity/dada2-pe/shannon_div/alpha-diversity.tsv')

metadata <- read_tsv('tourmaline/00-data/luminal_metadata.csv')

merged <- left_join(shannon, metadata) %>% left_join(., observed_otus)

# Alpha Div Visualizations ----------------------------------------------------------

merged$Day <- factor(merged$Day, cecum
                     levels=c("Day2",
                              "Day6",
                              "Day10")
                     )

merged$Tissue <- factor(merged$Tissue, 
                     levels=c("Ileum",
                              "Cecum",
                              "Sp-Colon")
)

ggplot(merged, aes(x=Treatment, y=shannon)) + 
  #geom_boxplot() +
  geom_boxplot() +
  facet_grid(Day ~ Tissue)



# Alpha diversity statistical analysis ------------------------------------

# Choose only day 2

merged_day2_cecum <- merged %>%
  filter(Day=="Day2") %>%
  filter(Tissue=="Cecum")

merged_day2_cecum$Treatment <- as.factor(merged_day2_cecum$Treatment)

#merge_by_day <- split(merged, merged$Day)

merged_kruskal <- kruskal.test(shannon ~ Treatment,data=merged_day2_cecum)

posthoc.kruskal.nemenyi.test(Treatment~, data=merged_day2_cecum, 
                             dist = "Chisquare")

# Wilcox test / Mann-Whitney test

wilcox.test(shannon ~ Treatment, 
            alternative="greater",
            data = merged_day2_cecum)


# Relative Abundance Visualizations ---------------------------------------

taxa_counts <- Sys.glob(here("tourmaline",
                             "04-diversity",
                             "dada2-pe",
                             "taxa_barplot_viz",
                             "level*.csv")
                        )

taxa_level7 <- read_csv(here("tourmaline",
                             "04-diversity",
                             "dada2-pe",
                             "taxa_barplot_viz",
                             "level-7.csv"))

taxa_level3 <- read_csv(here("tourmaline",
                             "04-diversity",
                             "dada2-pe",
                             "taxa_barplot_viz",
                             "level-3.csv"))

taxa_level2 <- read_csv(here("tourmaline",
                             "04-diversity",
                             "dada2-pe",
                             "taxa_barplot_viz",
                             "level-2.csv"))


taxa_count_df <- taxa_counts %>%
  map_df(
    ~ read_csv(.x)
  )

taxa_ileum_day2 <- taxa_count_df %>%
  filter(Day=="2") %>%
  filter(Tissue=="Ileum")

taxa_ileum_day2_lev7 <- taxa_level7 %>%
  filter(Day=="2") %>%
  filter(Tissue=="Ileum")

taxa_ileum_day2_lev3 <- taxa_level3 %>%
  filter(Day==2) %>%
  filter(Tissue=="Ileum")

taxa_ileum_day2_lev2 <- taxa_level2 %>%
  filter(Day==2) %>%
  filter(Tissue=="Ileum")

taxa_cecum_day2_lev2 <- taxa_level2 %>%
  filter(Day==2) %>%
  filter(Tissue=="Cecum")

taxa_cecum_lev2 <- taxa_level2 %>%
  filter(Tissue=="Cecum")

taxa_ileum_day2_lev7 <- taxa_level7 %>%
  tidyr::gather(., key="taxon",value="counts", 2:398)

taxa_ileum_day2_lev2 <- taxa_ileum_day2_lev2 %>%
  tidyr::gather(., key="taxon",value="counts", 2:20)

taxa_cecum_day2_lev2 <- taxa_cecum_day2_lev2 %>%
  tidyr::gather(., key="taxon",value="counts", 2:20)

taxa_cecum_lev2 <- taxa_cecum_lev2 %>%
  tidyr::gather(., key="taxon",value="counts", 2:20)

taxa_ileum_day2_lev2 <- taxa_ileum_day2_lev2 %>%
  filter(taxon != "D_0__Bacteria;__")

taxa_cecum_day2_lev2 <- taxa_cecum_day2_lev2 %>%
  filter(taxon != "D_0__Bacteria;__")

taxa_cecum_lev2 <- taxa_cecum_lev2 %>%
  filter(taxon != "D_0__Bacteria;__")


newCrayolaPalette6 <- c(
  "#CB7119", 
  "#404E5A",
  "#5F4F3A",
  "#D6AEDD",
  "#FEBAAD",
  "#6F7285",
  "#803790",
  "#0095B6",
  "#FFCD48",
  "#F653A6"
)

newCrayolaPalette7 <- c(
  "#FE6F5E",
  "#346114",
  "#CB7119",
  "#D6AEDD",
  "#0095B6",
  "#F653A6",
  "#803790",
  "#FFCD48",
  "#631F41",
  "#404E5A"
) 

dichromat_palette <- c(
  sample(dichromat::colorschemes$Categorical.12, 5)
)

new_palette <- c(newCrayolaPalette6, newCrayolaPalette7, dichromat_palette)

ggplot(taxa_ileum_day2_lev2, 
       aes(x=Treatment, y=counts, fill=taxon))+
  geom_bar(stat = "identity", position="fill") +
  scale_fill_manual(values = new_palette)


write_csv(taxa_ileum_day2_lev2, 'taxa_ileum_day2_lev2.csv')

day2 <- merged %>%
  filter(Day == 2)

day10 <- merged %>%
  filter(Day == 10)

day6 <- merged %>%
  filter(Day == 6)

ggplot(taxa_cecum_day2_lev2, 
       aes(x=Treatment, y=counts, fill=taxon))+
  geom_bar(stat = "identity", position="fill") +
  scale_fill_manual(values = new_palette)

ggplot(taxa_cecum_lev2, 
       aes(x=Treatment, y=counts, fill=taxon))+
  geom_bar(stat = "identity", position="fill") +
  scale_fill_manual(values = new_palette) +
  facet_grid(Day ~.)



# Other plots -------------------------------------------------------------


ggplot(day10, aes(x=Tissue, y=shannon)) + 
  geom_boxplot() + 
  facet_grid(. ~ Treatment)


# Ordination --------------------------------------------------------------

# Import data from QIIME2 into R

library(qiime2R)
library(phyloseq)
library(ggplot2)

# Example to import the QZA file of the Unweighted Unifrac

unweighted_unifrac <- read_qza("~/inglis/luminal_content/tourmaline/04-diversity/dada2-pe/unweighted_unifrac_pcoa_results.qza")
weighted_unifrac <- read_qza("~/inglis/luminal_content/tourmaline/04-diversity/dada2-pe/weighted_unifrac_pcoa_results.qza")

# Read QIIME objects as phyloseq
# Attempt ordination manually with PhyloSeq


luminal_qiime <- qiime2R::qza_to_phyloseq(
  features = "~/inglis/luminal_content/tourmaline/02-denoised/dada2-pe/table.qza",
  tree = "~/inglis/luminal_content/tourmaline/03-repseqs/dada2-pe/rooted_tree.qza",
  taxonomy = "~/inglis/luminal_content/tourmaline/03-repseqs/dada2-pe/taxonomy.qza",
  metadata = "~/inglis/luminal_content/tourmaline/00-data/luminal_metadata.csv"
)

luminal_qiime@sam_data$Day <- factor(
  luminal_qiime@sam_data$Day, 
  levels=c("Day2","Day6","Day10")
)

luminal_unifrac <- UniFrac(luminal_qiime, weighted = T, parallel = T)

luminal_ordinate <- ordinate(
  luminal_qiime, 
  method= "PCoA", 
  distance = "luminal_unifrac"
  )

luminal_ord <- plot_ordination(
  luminal_qiime, luminal_ordinate, 
  type = "samples", color="Tissue", shape="Treatment") +
  geom_point(size=5) +
  facet_grid(. ~ Day)
  #stat_ellipse(geom="polygon", alpha=0.3, aes(fill=Treatment))

luminal_ord



