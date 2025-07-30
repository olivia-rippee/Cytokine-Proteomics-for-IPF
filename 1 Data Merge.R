library(tidyverse, warn.conflicts=FALSE)
library(dplyr)
library(stringr)
library(readxl)


# ---------------------------
# MERGE 1
# ---------------------------

# Read in data
assay_data <- read.csv("assay_data.csv", 
                       sep=";")

clinical_data <- read_excel("clinical_data.xlsx", 
                       range = cell_cols("A:J"),
                       col_types = c("text", "text", "text", "text", "numeric", "text", 
                                     "text", "text", "numeric", "text")) # otherwise ID is numeric

colnames(clinical_data) <- c("Box_ID", "Well_ID","Unique_Sample_ID", "Sample_type",
                             "Sample_volume", "Diagnosis", "Site", 
                             "Sex", "Age_draw", "Race")


# Set QC fail NPX to NA
assay_data$NPX[assay_data$QC_Warning == "WARN"] <- NA


# order assay data by increasing Sample ID
assay_data <- assay_data[order(assay_data$SampleID),]

# Remove index variable
assay_data <- assay_data[ , !(names(assay_data) %in% c("Index"))]

# Remove controls
assay_data <- assay_data[!grepl("^CONTROL", assay_data$SampleID), ]

# Dataframe for assay warning
assay_warning_wider <- pivot_wider(assay_data,
                                   names_from = OlinkID,
                                   id_cols = SampleID, 
                                   values_from = Assay_Warning)

write.csv(assay_warning_wider, "Datasets/assay_warning.csv", 
          row.names=FALSE)


# Pivot wider using Sample ID as row names, OlinkID as cytokine (column) name, and NPX as cell values 
assay_data_wider <- pivot_wider(assay_data,
                                names_from = OlinkID,
                                id_cols = SampleID, 
                                values_from = NPX)


# Merge clinical and assay data
merged_data <- merge(clinical_data, assay_data_wider,
                     by.x = "Unique_Sample_ID", by.y = "SampleID", 
                     all.x = TRUE)

write.csv(merged_data, "Datasets/merged_data.csv", 
          row.names = FALSE)


# ---------------------------
# MERGE 2
# ---------------------------

# Load in merged_data (from above)
# ---------------------------------
merged_data <- read.csv("Datasets/merged_data.csv")
merged_data$Race <- str_to_title(merged_data$Race)
merged_data$Sex <- str_to_title(merged_data$Sex)
merged_data$Race <- case_when(
  merged_data$Race == "White" ~ "White",
  merged_data$Race == "Black" ~ "Black",
  merged_data$Race == "Other/Uk" ~ "Other/Unknown",
  merged_data$Race == "Asian" ~ "Asian",
  T ~ "Other/Unknown")



# Load in new clinical data
# ---------------------------------
clinical_data_2 <- read.csv("clinical_data_2.csv")

clinical_data_2$Race <- case_when(
  clinical_data_2$Race == "White" ~ "White",
  clinical_data_2$Race == "Black" ~ "Black",
  clinical_data_2$Race == "Other" ~ "Other/Unknown",
  clinical_data_2$Race == "Asian" ~ "Asian",
  T ~ "Other/Unknown")

clinical_data_2$Ethnicity <- case_when(
  clinical_data_2$Ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino",
  clinical_data_2$Ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
  T ~ "Unknown")


# Add shorthand for diagnosis to clinical_data_2
# ---------------------------------
clinical_data_2 <- clinical_data_2 %>% 
  mutate(Diag_Short = case_when(str_detect(tolower(Primary_Diag_Enroll), "idiopathic pulmonary fibrosis") ~ "IPF",
                                str_detect(tolower(Primary_Diag_Enroll), "scleroderma") ~ "SSc",
                                str_detect(tolower(Primary_Diag_Enroll), "rheumatoid arthritis") ~ "RA",
                                TRUE ~ "Other")) %>%
  relocate(Diag_Short, .before="Primary_Diag_Enroll")


# Load in IDs
# ---------------------------------

IDs <- read_excel("SSID to Plasma ID.xlsx",
                  range = cell_cols("A:K"), col_type = c("text", "text", "text", "text", 
                                        "text", "numeric", "text", "text", 
                                        "text", "numeric", "text"))
colnames(IDs)[1:2] <- c("SSID", "Unique_Sample_ID")
IDs$Sex <- str_to_title(IDs$Sex)
IDs$Race <- str_to_title(IDs$Race)
IDs$Race <- case_when(
  IDs$Race == "White" ~ "White",
  IDs$Race == "Black" ~ "Black",
  IDs$Race == "Hisp" ~ "Hisp",
  IDs$Race == "Other/Uk" ~ "Other/Unknown",
  IDs$Race == "Asian" ~ "Asian",
  T ~ "Other/Unknown")



### SKIP ###

# Isolate ID lists
#----------------------
#merged_dat <- merged_data[,c("Unique_Sample_ID", "Site", "Diagnosis")]
#write.csv(merged_dat, "Datasets/assay_IDs.csv", 
#          row.names = FALSE)

#clin_dat <- clinical_data_2[,c("SSID", "Site", "Diag_Short")] %>%
#  filter(Diag_Short %in% c("IPF", "SSc", "RA"))
#write.csv(clin_dat, "Datasets/clinical_IDs.csv", 
#          row.names = FALSE)

### END SKIP ###




# Drop Site variable (numeric here, char in merged_data) for merge. Would need Site decoder to merge by Site.
clinical_data_2$Site <- NULL


# Create merged dataset with just PFF Site
# ------------------------------------------

# Merge clinical_data_2 with IDs 
# NOTE/FIX: (HISP is in IDs but not clinical_2)
clinical_data_3 <- merge(clinical_data_2, IDs, by = c("SSID", "Sex")) %>%
  relocate(Race.y, .after="Race.x")
  #[unique(clinical_data_3$Site) --> [1] "PFF"]


# Merge clinical with assay data 
merged_PFF <- merge(clinical_data_3, merged_data, by = c("Unique_Sample_ID", "Sex")) %>% 
  relocate(Race, .before="Race.x")


# Remove diagnoses other than IPF, RA, and SSc
merged_PFF <- filter(merged_PFF, Diag_Short != "Other")



# Rename OlinkID as UniProt
# ---------------------------

UniProt_OlinkID <- assay_data[, c("OlinkID", "UniProt")]
UniProt_OlinkID <- UniProt_OlinkID[!duplicated(UniProt_OlinkID), ]

id_map <- setNames(UniProt_OlinkID$UniProt, UniProt_OlinkID$OlinkID)

colnames(merged_PFF) <- ifelse(colnames(merged_PFF) %in% names(id_map),
                               id_map[colnames(merged_PFF)],
                               colnames(merged_PFF))

write.csv(UniProt_OlinkID, "Datasets/UniProt_OlinkIDs.csv", 
          row.names=FALSE)

write.csv(merged_PFF, "Datasets/merged_PFF.csv",
          row.names=FALSE)

