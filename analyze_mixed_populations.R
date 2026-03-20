library(purrr)
library(ggplot2)
library(flowCore)
library(stringr)
library(dplyr)
library(uwot)
library(FlowSOM)
library(class)

# ---------------------------------------------------------------------
# 1. Get all the .fcs files in the working directory. In the example, 
# there are 21 single color controls and one mixed file to be analyzed.
# ---------------------------------------------------------------------
## Required info to be filled
setwd("E:/Brown/Lab/Thesis/Data/double color/data")
sample_to_analysis = "16mix_ad"

fcs_files <- list.files(pattern = "\\.fcs$", full.names = FALSE)

read_and_process <- function(filename) {
  fcs <- read.FCS(filename, transformation = FALSE, truncate_max_range = FALSE)
  exprs_data <- exprs(fcs)
  # Arcsinh transform
  cofactor <- 150
  exprs_data <- asinh(exprs_data / cofactor)
  df <- as.data.frame(exprs_data)
  try({
    markers <- markernames(fcs)
    names(df)[match(names(markers), names(df))] <- markers
  }, silent = TRUE)
  df$Sample_ID <- str_remove(filename, "\\.fcs")
  return(df)
}

combined_df <- map_dfr(fcs_files, read_and_process)

print(dim(combined_df))
head(combined_df)


# ---------------------------------------------------------------------
# 2. Applying downsampling strategy and run umap to separate different 
# groups
# ---------------------------------------------------------------------
# Downsampling
set.seed(1)
df_subset <- combined_df %>%
  dplyr::group_by(Sample_ID) %>%
  dplyr::mutate(n_cells = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(keep_n = ifelse(Sample_ID == sample_to_analysis, 20000, 2000)) %>%
  dplyr::group_by(Sample_ID) %>%
  dplyr::sample_n(size = min(keep_n, n()))

# choosing channels for UMAP
markers_to_use <- colnames(df_subset)[!colnames(df_subset) %in% 
                                        c("Time", "FSC-A", "FSC-H", "SSC-A",
                                          "Sample_ID", "n_cells", "keep_n",
                                          "SSC-H", "SSC-B-H", "SSC-B-A")]
print("Channels to be analyzed:")
print(markers_to_use)

# run UMAP
set.seed(1)
umap_out <- umap(df_subset[, markers_to_use], 
                 n_neighbors = 30, min_dist = 0.2)
df_subset$UMAP_1 <- umap_out[, 1]
df_subset$UMAP_2 <- umap_out[, 2]
df_controls <- df_subset %>% dplyr::filter(Sample_ID != sample_to_analysis)

# ---------------------------------------------------------------------
# 3. Plotting for all the controls
# ---------------------------------------------------------------------
df_labels <- df_controls %>%
  group_by(Sample_ID) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
unique_samples <- unique(df_controls$Sample_ID) 
sample_map <- data.frame(
  Sample_ID = unique_samples
) %>%
  mutate(Group_Num = 1:n()) %>% 
  mutate(Legend_Label = paste0(Group_Num, ". ", Sample_ID))
sample_map$Legend_Label <- factor(
  sample_map$Legend_Label, 
  levels = sample_map$Legend_Label[order(sample_map$Group_Num)]
)
df_controls <- df_controls %>% left_join(sample_map, by = "Sample_ID")
df_labels <- df_labels %>% left_join(sample_map, by = "Sample_ID")

# Raw figure
ggplot(df_controls, aes(x = UMAP_1, y = UMAP_2, color = Sample_ID)) +
  geom_point(size = 0.5, alpha = 0.6) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 3)))
# Figure with detailed labels
ggplot(df_controls, aes(x = UMAP_1, y = UMAP_2, color = Legend_Label)) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_text(
    data = df_labels, 
    aes(label = Group_Num), 
    color = "black",
    fontface = "bold",
    size = 4
  ) +
  theme_bw() +
  guides(color = guide_legend(
    title = "Groups",
    override.aes = list(size = 4, alpha = 1),
    ncol = 1 
  )) +
  scale_color_discrete()


# ---------------------------------------------------------------------
# 4. KNN for mix separation
# ---------------------------------------------------------------------
train_df <- df_subset %>% dplyr::filter(Sample_ID != sample_to_analysis)
test_df <- df_subset %>% dplyr::filter(Sample_ID == sample_to_analysis)
train_matrix <- train_df[, markers_to_use]
test_matrix  <- test_df[, markers_to_use]
train_labels <- train_df$Sample_ID

print(paste("Training size:", nrow(train_matrix)))
print(paste("Test size:", nrow(test_matrix)))

set.seed(123)
# Run kNN
predicted_labels <- knn(train = train_matrix, 
                        test  = test_matrix, 
                        cl    = train_labels, 
                        k     = 15)

test_df$Predicted_Identity <- predicted_labels
head(test_df[, c("Sample_ID", "Predicted_Identity")])

composition_stats <- test_df %>%
  count(Predicted_Identity) %>%
  mutate(Percentage = n / sum(n) * 100) %>%
  arrange(desc(Percentage))
print(composition_stats)

# Lollipop Chart for prediction
ggplot(composition_stats, aes(x = reorder(Predicted_Identity, Percentage), y = Percentage)) +
  geom_segment(aes(xend = Predicted_Identity, yend = 0), color = "grey") +
  geom_point(size = 4, color = "steelblue") +
  geom_text(aes(label = round(Percentage, 1)), vjust = -1, size = 3.5) +
  scale_y_continuous(limits = c(0, max(composition_stats$Percentage) + 5)) +
  coord_flip() +
  theme_bw() +
  labs(title = "Composition of the 14-mix Sample",
       subtitle = "Based on kNN classification (k=15)",
       x = "Cell Type",
       y = "Percentage (%)")

# ---------------------------------------------------------------------
# 5. Plotting for the mixed sample
# ---------------------------------------------------------------------
test_df <- test_df %>%
  left_join(sample_map, by = c("Predicted_Identity" = "Sample_ID"))
test_df$Legend_Label <- factor(
  test_df$Legend_Label,
  levels = sample_map$Legend_Label[order(sample_map$Group_Num)]
)
test_labels <- test_df %>%
  group_by(Group_Num, Legend_Label) %>%
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

# Raw figure
ggplot(test_df, aes(x = UMAP_1, y = UMAP_2, color = Predicted_Identity)) +
  geom_point(size = 0.5, alpha = 0.6) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Predicted Identities projected on sample UMAP",
       color = "Predicted Type")
# Figure with detailed labels
ggplot(test_df, aes(x = UMAP_1, y = UMAP_2, color = Legend_Label)) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_text(
    data = test_labels,
    aes(label = Group_Num),
    color = "black",
    fontface = "bold",
    size = 4,
    check_overlap = TRUE
  ) +
  theme_bw() +
  guides(color = guide_legend(
    title = "Predicted Type",
    override.aes = list(size = 4, alpha = 1),
    ncol = 1
  )) +
  labs(
    title = "Predicted Identities projected on sample UMAP",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_discrete()



