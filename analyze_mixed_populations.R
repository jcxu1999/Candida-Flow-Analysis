library(purrr)
library(ggplot2)
library(flowCore)
library(stringr)
library(dplyr)
library(uwot)
library(FlowSOM)
library(class)
library(caret)

# ---------------------------------------------------------------------
# 1. Get all the .fcs files in the working directory. In the example, 
# there are 21 single color controls and one mixed file to be analyzed.
# ---------------------------------------------------------------------
## Required info to be filled
setwd("//path_to_data")
sample_to_analysis = "21mix_equal"

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
set.seed(10)
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
set.seed(10)
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
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
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
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  ) +
  guides(color = guide_legend(
    title = "Groups",
    override.aes = list(size = 4, alpha = 1),
    ncol = 2
  )) +
  scale_color_discrete()

# ---------------------------------------------------------------------
# 4. Model Validation, Parameter Optimization, and Cross-Validation
# ---------------------------------------------------------------------
set.seed(123)
df_controls_split <- df_subset %>%
  dplyr::filter(Sample_ID != sample_to_analysis) %>%
  dplyr::group_by(Sample_ID) %>%
  dplyr::mutate(split = sample(c("train", "test"), 
                               size = n(), replace = TRUE, 
                               prob = c(0.7, 0.3))) %>%
  dplyr::ungroup()

train_df <- df_controls_split %>% dplyr::filter(split == "train")
test_df_internal  <- df_controls_split %>% dplyr::filter(split == "test")

train_matrix <- train_df[, markers_to_use]
test_matrix_internal  <- test_df_internal[, markers_to_use]

train_labels <- factor(train_df$Sample_ID)
test_labels_internal <- factor(test_df_internal$Sample_ID, levels = levels(train_labels))

# cross-validation
k_values <- c(5, 10, 15, 20, 25, 30)
cv_results <- data.frame(k = k_values, accuracy = NA)

set.seed(123)
folds <- createFolds(train_labels, k = 5, list = TRUE)

for (i in seq_along(k_values)) {
  k_val <- k_values[i]
  fold_acc <- c()
  for (fold in folds) {
    cv_train_x <- train_matrix[-fold, ]
    cv_train_y <- train_labels[-fold]
    cv_test_x  <- train_matrix[fold, ]
    cv_test_y  <- train_labels[fold]
    
    pred <- knn(train = cv_train_x, test = cv_test_x, cl = cv_train_y, k = k_val)
    fold_acc <- c(fold_acc, mean(pred == cv_test_y))
  }
  cv_results$accuracy[i] <- mean(fold_acc)
}

print(cv_results)
best_k <- cv_results$k[which.max(cv_results$accuracy)]
print(paste("Optimal k:", best_k))

ggplot(cv_results, aes(x = k, y = accuracy)) +
  geom_line() + geom_point(size = 3) +
  theme_bw() +
  labs(title = "5-fold CV Accuracy vs. k",
       x = "k (number of neighbors)", y = "Mean CV Accuracy") +
  scale_y_continuous(limits = c(0.99, 1))

# Confusion matrix
set.seed(123)
pred_internal <- knn(train = train_matrix, test = test_matrix_internal, 
                     cl = train_labels, k = best_k)

cm <- confusionMatrix(pred_internal, test_labels_internal)

# Overall accuracy
print(cm$overall["Accuracy"])

# Accuracy for each color population
print(cm$byClass[, c("Sensitivity", "Specificity", "Precision", "F1", "Balanced Accuracy")])

# Confusion matrix
cm_df <- cm$table
cm_prop <- as.data.frame(prop.table(cm_df, margin = 2))
cm_counts <- as.data.frame(cm_df)
cm_data <- merge(cm_counts, cm_prop, by = c("Prediction", "Reference"))
cm_data$Percentage <- cm_data$Freq.y * 100  # Covert to 0-100%
ggplot(cm_data, aes(x = Reference, y = Prediction, fill = Percentage)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.1f", Percentage)), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)) +
  labs(title = paste0("Confusion Matrix (k = ", best_k, ")"),
       subtitle = paste("Overall Accuracy:", round(cm$overall["Accuracy"] * 100, 2), "%"),
       x = "True Population", y = "Predicted Population")

# ---------------------------------------------------------------------
# 5. KNN for mix separation
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
                        k     = best_k)

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
       subtitle = paste0("Based on kNN classification (k=", best_k, ")"),
       x = "Cell Type",
       y = "Percentage (%)")

# ---------------------------------------------------------------------
# 6. Plotting for the mixed sample
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
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  ) +
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
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  ) +
  guides(color = guide_legend(
    title = "Predicted Type",
    override.aes = list(size = 4, alpha = 1),
    ncol = 2
  )) +
  labs(
    title = "Predicted Identities projected on sample UMAP",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_discrete()



