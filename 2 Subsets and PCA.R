library(tidyverse, warn.conflicts=FALSE)
library(dplyr)
library(stringr)
library(readxl)
library(ggplot2)
library(patchwork)
library(jpeg)


merged_PFF <- read.csv("Datasets/merged_PFF.csv")

UniProt_OlinkID <- read.csv("Datasets/UniProt_OlinkIDs.csv")



# ------------------------------------------
# Method 1: Impute LOD/2
# ------------------------------------------

LODcols <- c("UniProt", "ProteinName", "GeneName", "DilutionFactor", "Panel", "LOD")

LODs <- read_excel("olink-LOD-data.xlsx", 
                   range = cell_cols("A:F"), col_names = LODcols)
LODs <- LODs[-c(1:2),]            # remove header rows
LODs$LOD <- as.numeric(LODs$LOD)  # convert LOD to numeric for calculations


# Imputation loop
# ----------------
cols_to_remove <- c()
merged_PFF_impute <- merged_PFF #needed to not impute on the original dataset

for (j in 146:ncol(merged_PFF_impute)) {
  
  # Get UniProt ID from column name
  uniprot_id <- colnames(merged_PFF_impute)[j]
  
  # Look up LOD value
  lod_value <- LODs$LOD[LODs$UniProt == uniprot_id]
  
  if (length(lod_value) == 1 && !is.na(lod_value)) {
    # Impute NA values with LOD/2
    na_indices <- which(is.na(merged_PFF_impute[[j]]))
    merged_PFF_impute[na_indices, j] <- lod_value / 2
  } 
  
  else {
    # No valid LOD found, mark column for removal
    cols_to_remove <- c(cols_to_remove, j)
  }
}

# Remove columns with no LOD match
if (length(cols_to_remove) > 0) {
  merged_PFF_impute <- merged_PFF_impute[, -cols_to_remove]
}
  # this removed 535 columns


any(is.na(merged_PFF_impute[, 146:ncol(merged_PFF_impute)]))
  # FALSE


write.csv(merged_PFF_impute, "Datasets/merged_PFF_imputed.csv", 
          row.names=FALSE)


# PCA
# -------------------
# n=1036 (subjects), j=2404 (cytokines)
pcaobject <- prcomp(merged_PFF_impute[,-c(1:145)], center = T, scale. = T, retx = T)
variance_explained <- pcaobject$sdev^2 / sum(pcaobject$sdev^2)
percent_var <- round(variance_explained * 100, 1)

pca_df <- data.frame(PC1 = pcaobject$x[, 1],
                     PC2 = pcaobject$x[, 2],
                     Diagnosis = merged_PFF_impute$Diag_Short)

# PC1 explains 11.6%, PC2 explains 8.06%

pcs <- pcaobject$x
dim(pcs) # 1036 by 1036

PCA_impute <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Diagnosis)) +
  geom_point(size = 1.5, alpha = 0.7) +
  labs(x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)")) +
  theme_minimal() + 
  scale_color_brewer(palette = "Set1")


ggsave("Plots/PCA_impute.jpg", plot = PCA_impute)


# ------------------------------------------
# Method 2: Remove subjects with at least 20% NAs, then all cytokines with missing values
# ------------------------------------------

# Missingness fraction for each subject
# --------------------------------

# Count number of missing
missing <- c()
  
for (i in 1:nrow(merged_PFF)) {
  missing <- c(missing, sum(is.na(merged_PFF[i, 146:ncol(merged_PFF)])))
}

# Calculate fraction of missingness per subject
missing_frac <- missing/(ncol(merged_PFF)-145) #2939-145
df_missing_frac <- data.frame(cbind(missing_frac, merged_PFF[,"SSID"]))
class(df_missing_frac$missing_frac) <- "numeric"
colnames(df_missing_frac)[2] <- "SSID"

# Sort by decreasing fraction
df_missing_frac <- df_missing_frac[order(-df_missing_frac$missing_frac), ]
df_missing_frac$SSID <- factor(df_missing_frac$SSID, levels = df_missing_frac$SSID)

# How many subjects have at least 20% missing
sum(df_missing_frac$missing_frac >= 0.20)

# Bar plot
bar_plot <- ggplot(df_missing_frac, aes(x = SSID, y = missing_frac)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0.2) +
  labs(title = "Missingness Fraction by Subject",
       x = "Subject",
       y = "Missing Fraction") + 
  theme_minimal() +
  theme(axis.text.x = element_blank())


ggsave("Plots/missingness_fraction.jpg", bar_plot)


# Create dataframe
# --------------------

# Remove subjects with at least 20% missing
merged_PFF_rm_subjects <- merged_PFF[rowMeans(is.na(merged_PFF[, 146:ncol(merged_PFF)])) < 0.20, ]

# Remove cytokines with any missing values
df_first_145 <- merged_PFF_rm_subjects[, 1:145]
df_rest <- merged_PFF_rm_subjects[, 146:ncol(merged_PFF_rm_subjects)]
df_rest_noNA <- df_rest[, colSums(is.na(df_rest)) == 0]
merged_PFF_rm_subjects <- cbind(df_first_145, df_rest_noNA)


any(is.na(merged_PFF_rm_subjects[, 146:ncol(merged_PFF_rm_subjects)]))
# FALSE

write.csv(merged_PFF_rm_subjects, "Datasets/merged_PFF_remove_subjects.csv", 
          row.names=FALSE)


# PCA
# -------------------
# n=1009 (subjects), j=94 (cytokines)
pcaobject <- prcomp(merged_PFF_rm_subjects[,-c(1:145)], center = T, scale. = T, retx = T)

variance_explained <- pcaobject$sdev^2 / sum(pcaobject$sdev^2)
percent_var <- round(variance_explained * 100, 1)

pca_df <- data.frame(PC1 = pcaobject$x[, 1],
                     PC2 = pcaobject$x[, 2],
                     Diagnosis = merged_PFF_rm_subjects$Diag_Short)

# PC1 explains 22.1%, PC2 explains 20.1%

pcs <- pcaobject$x
dim(pcs) # 1009 by 94

PCA_remove_subjects <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Diagnosis)) +
  geom_point(size = 1.5, alpha = 0.7) +
  labs(x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)")) +
  theme_minimal() + 
  scale_color_brewer(palette = "Set1")

ggsave("Plots/PCA_remove_subjects.jpg", plot = PCA_remove_subjects)



# ------------------------------------------
# Method 3: Remove all cytokines with missing values
# ------------------------------------------

# Remove cytokines with any missing values
merged_PFF_rm_na <- merged_PFF[, colSums(is.na(merged_PFF)) == 0]

any(is.na(merged_PFF_rm_na[, 146:ncol(merged_PFF_rm_na)]))
# FALSE

write.csv(merged_PFF_rm_na, "Datasets/merged_PFF_remove_na.csv", 
          row.names=FALSE)



# PCA
# --------------

# n=1036 (subjects), j=62 (cytokines)
pcaobject <- prcomp(merged_PFF_rm_na[,-c(1:145)], center = T, scale. = T, retx = T)

variance_explained <- pcaobject$sdev^2 / sum(pcaobject$sdev^2)
percent_var <- round(variance_explained * 100, 1)

pca_df <- data.frame(PC1 = pcaobject$x[, 1],
                     PC2 = pcaobject$x[, 2],
                     Diagnosis = merged_PFF_rm_na$Diag_Short)

# PC1 explains 22.0%, PC2 explains 21.1%

PCA_remove_NAs <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Diagnosis)) +
  geom_point(size = 1.5, alpha = 0.7) +
  labs(x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)")) +
  theme_minimal() + 
  scale_color_brewer(palette = "Set1")


ggsave("Plots/PCA_remove_NAs.jpg", plot = PCA_remove_NAs)



# ------------------------------------------
# Check whether the same cytokines are left for methods 2 and 3
# ------------------------------------------

# Get column names
colnames_df1 <- colnames(merged_PFF_rm_subjects[ , 146:ncol(merged_PFF_rm_subjects)])
colnames_df2 <- colnames(merged_PFF_rm_na[ , 146:ncol(merged_PFF_rm_na)])

# Find common column names
common_cols <- intersect(colnames_df1, colnames_df2)

# Number of common column names
length(common_cols)
  # 62 cytokines are the same (all of method 3)





# ------------------------------------------
# Method 4: Look for IPF subtypes
# ------------------------------------------

# Filter imputed dataset for just IPF
# --------------
merged_PFF_IPF <- filter(merged_PFF_impute, Diag_Short == "IPF")
# 914 X 2549

write.csv(merged_PFF_IPF, "Datasets/merged_PFF_IPF.csv", 
          row.names=FALSE)



# PCA
# --------------

# n=1036 (subjects), j=2404 (cytokines)
pcaobject <- prcomp(merged_PFF_IPF[,-c(1:145)], center = T, scale. = T, retx = T)

variance_explained <- pcaobject$sdev^2 / sum(pcaobject$sdev^2)
percent_var <- round(variance_explained * 100, 1)

pca_df <- data.frame(PC1 = pcaobject$x[, 1],
                     PC2 = pcaobject$x[, 2],
                     Diagnosis = merged_PFF_IPF$Diag_Short)

# PC1 explains 12%, PC2 explains 8.4%

PCA_IPF <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Diagnosis)) +
  geom_point(size = 1.5, alpha = 0.7) +
  labs(x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)")) +
  theme_minimal() + 
  scale_color_brewer(palette = "Set1")


ggsave("Plots/PCA_IPF.jpg", PCA_IPF)




# ------------------------------------------
# Method 5: Randomly select 61 IPF patients
# ------------------------------------------

# Loop to randomly sample 61 patients with IPF 10 times

PCA_plots <- list()

for (i in 1:10) {

  # Subset IPF rows and randomly sample 61
  IPF_subset <- merged_PFF_impute %>% 
    filter(Diag_Short == "IPF") %>% 
    sample_n(61)
  
  # Subset RA and SSc rows
  non_IPF_all <- merged_PFF_impute %>% 
    filter(Diag_Short == "RA" | Diag_Short == "SSc")
  
  # Combine sampled IPF with non-IPF rows
  merged_PFF_IPF_sample <- bind_rows(IPF_subset, non_IPF_all)
  
  
  # PCA
  # --------------
  pcaobject <- prcomp(merged_PFF_IPF_sample[,-c(1:145)], center = T, scale. = T, retx = T)
  variance_explained <- pcaobject$sdev^2 / sum(pcaobject$sdev^2)
  percent_var <- round(variance_explained * 100, 1)
  
  pca_df <- data.frame(PC1 = pcaobject$x[, 1],
                       PC2 = pcaobject$x[, 2],
                       Diagnosis = merged_PFF_IPF_sample$Diag_Short)

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Diagnosis)) +
    geom_point(size = 1.5, alpha = 0.7) +
    labs(x = paste0("PC1 (", percent_var[1], "%)"),
         y = paste0("PC2 (", percent_var[2], "%)")) +
    theme_minimal() + 
    scale_color_brewer(palette = "Set1") +
    theme(plot.margin = margin(t=5, r=5, b=5, l=5, unit = "mm"))

  PCA_plots[[i]] <- p
}

PCA_IPF_random_sample <- wrap_plots(PCA_plots, nrow = 5, ncol = 2)

ggsave("Plots/PCA_IPF_random_sample.jpg", plot = PCA_IPF_random_sample, 
       device = "jpeg", width = 10, height = 18, dpi = 300)



