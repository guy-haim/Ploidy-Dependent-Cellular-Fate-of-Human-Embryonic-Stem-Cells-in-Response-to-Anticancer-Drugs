library(stringr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggsignif)
library(ggpubr)

setwd(choose.dir())
foci_number <- read.csv("Summary 1st-1_250 2nd-1_500 UT and 2.5uM bleo (Red threshold-1000) (ImageJ analysis).csv")

foci_number[foci_number$Count<=5,4] <- "0-5"
# foci_number[foci_number$Count>=1 & foci_number$Count<=5,4] <- "1-5"
foci_number[foci_number$Count>=6 & foci_number$Count<=10,4] <- "6-10"
foci_number[foci_number$Count>=11 & foci_number$Count<=15,4] <- "11-15"
foci_number[foci_number$Count>=16 & foci_number$Count<=20,4] <- "16-20"
foci_number[foci_number$Count>20 & foci_number$Count<=30,4] <- "21-30"
foci_number[foci_number$Count>30,4] <- ">30"

# Treatment-Ploidy Foci Distribution Analysis
# This script creates pie charts for foci distribution across treatment-ploidy combinations

# Create group variable for foci counts
foci_number$group <- cut(foci_number$Count, 
                         breaks = c(0, 2, 5, 10, 15, 20, 30, Inf),
                         labels = c("0-1","2-5", "6-10", "11-15", "16-20", "21-30", ">30"),
                         include.lowest = TRUE,
                         right = TRUE)

# Calculate percentages for each treatment-ploidy combination
pie_data <- foci_number %>%
  group_by(treatment, ploidy, group) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(treatment, ploidy) %>%
  mutate(percentage = Count / sum(Count) * 100,
         label = paste0(group, "\n", round(percentage, 1), "%"))

# Get unique treatment-ploidy combinations
treatment_ploidy_combos <- pie_data %>%
  select(treatment, ploidy) %>%
  distinct() %>%
  arrange(treatment, ploidy)

cat("\nCreating pie charts for", nrow(treatment_ploidy_combos), "treatment-ploidy combinations:\n")
print(treatment_ploidy_combos)

# Create pie chart for each treatment-ploidy combination
for (i in 1:nrow(treatment_ploidy_combos)) {
  current_treatment <- treatment_ploidy_combos$treatment[i]
  current_ploidy <- treatment_ploidy_combos$ploidy[i]
  
  # Filter data for current combination
  combo_data <- pie_data %>% 
    filter(treatment == current_treatment & ploidy == current_ploidy)
  
  # Create pie chart
  p <- ggplot(combo_data, aes(x = "", y = percentage, fill = group)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(title = paste("Foci Distribution"),
         subtitle = paste(current_treatment, "-", current_ploidy, "Ploidy"),
         fill = "Foci Group") +
    geom_text(aes(label = label), 
              position = position_stack(vjust = 0.5),
              size = 3.5) +
    scale_fill_viridis_d(option = "magma", begin = 0.5, end = 0.975, direction = -1) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 12))
  
  # Display plot
  print(p)
  
  # Save plot with descriptive filename
  filename <- paste0("Pie_chart_", 
                     gsub(" ", "_", current_treatment), 
                     "_", 
                     gsub(" ", "_", current_ploidy), 
                     "_ploidy.pdf")
  
  ggsave(filename = filename,
         plot = p,
         width = 8,
         height = 6,
         dpi = 300)
  
  cat("Saved:", filename, "\n")
}

# Create summary statistics table
cat("\n=== Summary statistics by treatment and ploidy ===\n")
summary_stats <- foci_number %>%
  group_by(treatment, ploidy) %>%
  summarise(
    n_cells = n(),
    mean_foci = mean(Count, na.rm = TRUE),
    median_foci = median(Count, na.rm = TRUE),
    sd_foci = sd(Count, na.rm = TRUE),
    min_foci = min(Count, na.rm = TRUE),
    max_foci = max(Count, na.rm = TRUE),
    .groups = 'drop'
  )
print(summary_stats)

# Optional: Create a comparison table showing distribution percentages
cat("\n=== Distribution percentages by group ===\n")
distribution_table <- pie_data %>%
  select(treatment, ploidy, group, percentage) %>%
  arrange(treatment, ploidy, group)
print(distribution_table)

# Alternative: Create all pie charts in one figure using facet_wrap
# Create combined label for faceting
pie_data <- pie_data %>%
  mutate(treatment_ploidy = paste0(treatment, "\n", ploidy, " Ploidy"))

all_pie <- ggplot(pie_data, aes(x = "", y = percentage, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  facet_wrap(~treatment_ploidy, ncol = 2) +
  theme_void() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  labs(title = "Foci Distribution by Treatment and Ploidy",
       fill = "Foci Group") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5),
            size = 2.5) +
  scale_fill_viridis_d(option = "magma", begin = 0.5, end = 0.975, direction = -1)

print(all_pie)

# Save the combined figure
ggsave(filename = "Pie_chart_all_treatment_ploidy_combinations.pdf",
       plot = all_pie,bg = "white",
       width = 12,
       height = 8,
       dpi = 300)

cat("\nSaved combined figure: Pie_chart_all_treatment_ploidy_combinations.pdf\n")

# BAR PLOTS
# Calculate percentages for each treatment-ploidy combination
bar_data <- foci_number %>%
  group_by(treatment, ploidy, group) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(treatment, ploidy) %>%
  mutate(percentage = Count / sum(Count) * 100)

# Get unique treatment-ploidy combinations
treatment_ploidy_combos <- bar_data %>%
  select(treatment, ploidy) %>%
  distinct() %>%
  arrange(treatment, ploidy)

cat("\nCreating bar plots for", nrow(treatment_ploidy_combos), "treatment-ploidy combinations:\n")

# Create bar plot for each treatment-ploidy combination
for (i in 1:nrow(treatment_ploidy_combos)) {
  current_treatment <- treatment_ploidy_combos$treatment[i]
  current_ploidy <- treatment_ploidy_combos$ploidy[i]
  
  # Filter data for current combination
  combo_data <- bar_data %>% 
    filter(treatment == current_treatment & ploidy == current_ploidy)
  
  # Create bar plot
  p <- ggplot(combo_data, aes(x = group, y = percentage, fill = group)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_text(aes(label = paste0(round(percentage, 1), "%")), 
              vjust = -0.5, size = 4) +
    theme_minimal() +
    labs(title = "Foci Distribution",
         subtitle = paste(current_treatment, "-", current_ploidy, "Ploidy"),
         x = "Foci Group",
         y = "Percentage (%)") +
    scale_fill_viridis_d(option = "magma", begin = 0.5, end = 0.975, direction = -1) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          legend.position = "none") +
    ylim(0, max(combo_data$percentage) * 1.15)
  
  # Display plot
  print(p)
  
  # Save plot
  filename <- paste0("Barplot_", 
                     gsub(" ", "_", current_treatment), 
                     "_", 
                     gsub(" ", "_", current_ploidy), 
                     "_ploidy.pdf")
  
  ggsave(filename = filename,
         plot = p,
         width = 8,
         height = 6,
         dpi = 300)
  
  cat("Saved:", filename, "\n")
}

# Alternative: Create all bar plots in one figure using facet_wrap
# Create combined label for faceting
bar_data <- bar_data %>%
  mutate(treatment_ploidy = paste0(treatment, "\n", ploidy, " Ploidy"))

all_bars <- ggplot(bar_data, aes(x = group, y = percentage, fill = group)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            vjust = -0.5, size = 2.5) +
  facet_wrap(~treatment_ploidy, ncol = 2, scales = "free_y") +
  theme_minimal() +
  labs(title = "Foci Distribution by Treatment and Ploidy",
       x = "Foci Group",
       y = "Percentage (%)") +
  scale_fill_viridis_d(option = "magma", begin = 0.5, end = 0.975, direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 9, face = "bold"),
        legend.position = "none")

print(all_bars)

# Save the combined bar plot figure
ggsave(filename = "Barplot_all_treatment_ploidy_combinations.pdf",
       plot = all_bars,
       width = 14,
       height = 10,
       dpi = 300)

cat("\nSaved combined bar plot: Barplot_all_treatment_ploidy_combinations.pdf\n")

# FOCI DISTRIBUTION BY PLOIDY PER TREATMENT WITH STATISTICAL TESTS
cat("\n=== Creating foci distribution plots by ploidy per treatment with statistics ===\n")

# Load required packages
library(ggpubr)
library(rstatix)

# Get unique treatments
treatments <- unique(foci_number$treatment)

# Create plots for each treatment
for (treat in treatments) {
  cat("\nProcessing treatment:", treat, "\n")
  
  # Filter data for current treatment
  treatment_data <- foci_number %>% filter(treatment == treat)
  
  # Get unique ploidy groups for this treatment
  ploidy_groups <- unique(treatment_data$ploidy)
  n_groups <- length(ploidy_groups)
  
  # Perform appropriate statistical test
  if (n_groups == 2) {
    # Wilcoxon test for 2 groups
    stat_test <- treatment_data %>%
      wilcox_test(Count ~ ploidy) %>%
      add_significance() %>%
      add_xy_position(x = "ploidy") %>%
      mutate(p.adj.signif = p.signif)  # Add p.adj.signif column for consistency
    test_name <- "Wilcoxon test"
  } else if (n_groups > 2) {
    # Kruskal-Wallis test for >2 groups
    kw_test <- treatment_data %>%
      kruskal_test(Count ~ ploidy)
    
    # Dunn's test for pairwise comparisons
    stat_test <- treatment_data %>%
      dunn_test(Count ~ ploidy, p.adjust.method = "bonferroni") %>%
      add_significance() %>%
      add_xy_position(x = "ploidy")
    
    test_name <- "Kruskal-Wallis with Dunn's post-hoc"
    cat("Kruskal-Wallis p-value:", kw_test$p, "\n")
  } else {
    cat("Only one ploidy group, skipping statistical test\n")
    stat_test <- NULL
    test_name <- "No test (single group)"
  }
  
  # Print statistical results
  if (!is.null(stat_test)) {
    cat("\nPairwise comparisons:\n")
    print(stat_test)
  }
  
  # Violin plot with boxplot overlay and statistics
  ploidy_violin <- ggplot(treatment_data, aes(x = ploidy, y = Count, fill = ploidy)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    theme_minimal() +
    labs(title = paste("Foci Count Distribution by Ploidy -", treat),
         subtitle = test_name,
         x = "Ploidy",
         y = "Number of Foci") +
    scale_fill_manual(values = c("royalblue","red3")) +
    theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
          legend.position = "none")
  
  # Add statistical annotations if available
  if (!is.null(stat_test)) {
    ploidy_violin <- ploidy_violin + 
      stat_pvalue_manual(stat_test, label = "p.adj.signif", 
                         tip.length = 0.01, hide.ns = FALSE, size = 4,
                         inherit.aes = FALSE)
  }
  
  print(ploidy_violin)
  ggsave(filename = paste0("Foci_distribution_by_ploidy_violin_", gsub(" ", "_", treat), ".pdf"),
         plot = ploidy_violin,
         width = 10,
         height = 7,
         dpi = 300)
  
  # Density plot
  ploidy_density <- ggplot(treatment_data, aes(x = Count, fill = ploidy)) +
    geom_density(alpha = 0.6) +
    theme_minimal() +
    labs(title = paste("Foci Count Density Distribution by Ploidy -", treat),
         subtitle = test_name,
         x = "Number of Foci",
         y = "Density",
         fill = "Ploidy") +
    scale_fill_manual(values = c("royalblue","red3")) +
    theme(axis.text = element_text(size = 11),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
          legend.position = "right",
          legend.title = element_text(face = "bold"))
  
  print(ploidy_density)
  ggsave(filename = paste0("Foci_distribution_by_ploidy_density_", gsub(" ", "_", treat), ".pdf"),
         plot = ploidy_density,
         width = 10,
         height = 6,
         dpi = 300)
  
  # Histogram with facets
  ploidy_histogram <- ggplot(treatment_data, aes(x = Count, fill = ploidy)) +
    geom_histogram(bins = 30, color = "black", alpha = 0.7) +
    facet_wrap(~ploidy, ncol = 1, scales = "free_y") +
    theme_minimal() +
    labs(title = paste("Foci Count Distribution by Ploidy -", treat),
         subtitle = test_name,
         x = "Number of Foci",
         y = "Frequency") +
    scale_fill_manual(values = c("royalblue","red3")) +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
          strip.text = element_text(size = 11, face = "bold"),
          legend.position = "none")
  
  print(ploidy_histogram)
  ggsave(filename = paste0("Foci_distribution_by_ploidy_histogram_", gsub(" ", "_", treat), ".pdf"),
         plot = ploidy_histogram,
         width = 10,
         height = 8,
         dpi = 300)
  
  # Statistical comparison for this treatment
  cat("\n=== Statistical summary by ploidy for", treat, "===\n")
  ploidy_stats <- treatment_data %>%
    group_by(ploidy) %>%
    summarise(
      n_cells = n(),
      mean = mean(Count, na.rm = TRUE),
      median = median(Count, na.rm = TRUE),
      sd = sd(Count, na.rm = TRUE),
      min = min(Count, na.rm = TRUE),
      max = max(Count, na.rm = TRUE),
      q25 = quantile(Count, 0.25, na.rm = TRUE),
      q75 = quantile(Count, 0.75, na.rm = TRUE),
      .groups = 'drop'
    )
  print(ploidy_stats)
}

# Combined plot: All treatments with ploidy comparison
cat("\n=== Creating combined plot for all treatments with statistics ===\n")

# Perform statistical tests for each treatment
stat_test_combined <- foci_number %>%
  group_by(treatment) %>%
  dunn_test(Count ~ ploidy, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "ploidy")

combined_violin <- ggplot(foci_number, aes(x = ploidy, y = Count, fill = ploidy)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  facet_wrap(~treatment, ncol = 2, scales = "free") +
  stat_pvalue_manual(stat_test_combined, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = FALSE, size = 3,
                     inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = "Foci Count Distribution by Ploidy Across Treatments",
       subtitle = "Dunn's test with Bonferroni correction",
       x = "Ploidy",
       y = "Number of Foci") +
  scale_fill_manual(values = c("royalblue","red3")) +
  theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "none")

print(combined_violin)
ggsave(filename = "Foci_distribution_by_ploidy_all_treatments.pdf",
       plot = combined_violin,
       width = 12,
       height = 10,
       dpi = 300)

cat("\nSaved all foci distribution plots by ploidy per treatment with statistical tests\n")

# KOLMOGOROV-SMIRNOV TEST FOR 2n vs 3n COMPARISON
cat("\n=== Kolmogorov-Smirnov test for 2n vs 3n comparison ===\n")

# Filter for only 2n and 3n ploidy
ks_data <- foci_number %>%
  filter(ploidy %in% c("2n", "3n"))

# Perform KS test for each treatment
ks_results <- list()

for (treat in unique(ks_data$treatment)) {
  cat("\nProcessing treatment:", treat, "\n")
  
  # Filter data for current treatment
  treatment_ks_data <- ks_data %>% filter(treatment == treat)
  
  # Check if both 2n and 3n exist for this treatment
  available_ploidy <- unique(treatment_ks_data$ploidy)
  
  if (length(available_ploidy) < 2) {
    cat("Skipping - only", available_ploidy, "available\n")
    next
  }
  
  # Get 2n and 3n data
  data_2n <- treatment_ks_data %>% filter(ploidy == "2n") %>% pull(Count)
  data_3n <- treatment_ks_data %>% filter(ploidy == "3n") %>% pull(Count)
  
  # Perform KS test
  ks_test <- ks.test(data_2n, data_3n)
  
  # Store results
  ks_results[[treat]] <- data.frame(
    treatment = treat,
    n_2n = length(data_2n),
    n_3n = length(data_3n),
    D_statistic = ks_test$statistic,
    p_value = ks_test$p.value,
    significance = case_when(
      ks_test$p.value <= 0.0001 ~ "****",
      ks_test$p.value <= 0.001 ~ "***",
      ks_test$p.value <= 0.01 ~ "**",
      ks_test$p.value <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
  
  cat("KS test D =", round(ks_test$statistic, 4), ", p-value =", 
      format.pval(ks_test$p.value, digits = 3), "\n")
  
  # Create comparison plot with KS test result
  ks_plot <- ggplot(treatment_ks_data, aes(x = Count, color = ploidy, fill = ploidy)) +
    geom_density(alpha = 0.3, linewidth = 1) +
    annotate("text", x = Inf, y = Inf, 
             label = paste0("KS test\nD = ", round(ks_test$statistic, 4),
                            "\np = ", format.pval(ks_test$p.value, digits = 3),
                            "\n", ks_results[[treat]]$significance),
             hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
    theme_minimal() +
    labs(title = paste("Distribution Comparison: 2n vs 3n -", treat),
         subtitle = "Kolmogorov-Smirnov test",
         x = "Number of Foci",
         y = "Density",
         color = "Ploidy",
         fill = "Ploidy") +
    scale_color_manual(values = c("2n" = "royalblue", "3n" = "red3")) +
    scale_fill_manual(values = c("2n" = "royalblue", "3n" = "red3")) +
    theme(axis.text = element_text(size = 11),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
          legend.position = "right",
          legend.title = element_text(face = "bold"))
  
  print(ks_plot)
  ggsave(filename = paste0("KS_test_2n_vs_3n_", gsub(" ", "_", treat), ".pdf"),
         plot = ks_plot,
         width = 10,
         height = 6,
         dpi = 300)
  
  # Create cumulative distribution plot
  ks_cdf_plot <- ggplot(treatment_ks_data, aes(x = Count, color = ploidy)) +
    stat_ecdf(linewidth = 1.5) +
    annotate("text", x = Inf, y = 0.5, 
             label = paste0("KS test\nD = ", round(ks_test$statistic, 4),
                            "\np = ", format.pval(ks_test$p.value, digits = 3),
                            "\n", ks_results[[treat]]$significance),
             hjust = 1.1, vjust = 0.5, size = 5, fontface = "bold") +
    theme_minimal() +
    labs(title = paste("Cumulative Distribution: 2n vs 3n -", treat),
         subtitle = "Kolmogorov-Smirnov test",
         x = "Number of Foci",
         y = "Cumulative Probability",
         color = "Ploidy") +
    scale_color_manual(values = c("2n" = "#440154", "3n" = "#FDE724")) +
    theme(axis.text = element_text(size = 11),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
          legend.position = "right",
          legend.title = element_text(face = "bold"))
  
  print(ks_cdf_plot)
  ggsave(filename = paste0("KS_test_CDF_2n_vs_3n_", gsub(" ", "_", treat), ".pdf"),
         plot = ks_cdf_plot,
         width = 10,
         height = 6,
         dpi = 300)
}

# Combine all KS test results
ks_summary <- do.call(rbind, ks_results)
rownames(ks_summary) <- NULL

cat("\n=== Summary of Kolmogorov-Smirnov tests ===\n")
print(ks_summary)

# Create summary barplot of KS statistics
if (nrow(ks_summary) > 0) {
  ks_summary_plot <- ggplot(ks_summary, aes(x = treatment, y = D_statistic, fill = significance)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_text(aes(label = significance), vjust = -0.5, size = 5, fontface = "bold") +
    theme_minimal() +
    labs(title = "Kolmogorov-Smirnov Test Summary: 2n vs 3n",
         subtitle = "D statistic measures maximum difference between cumulative distributions",
         x = "Treatment",
         y = "D Statistic",
         fill = "Significance") +
    scale_fill_manual(values = c("****" = "#d62728", "***" = "#ff7f0e", 
                                 "**" = "#ffbb78", "*" = "#fdd835", "ns" = "#cccccc")) +
    theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray30"),
          legend.position = "right")
  
  print(ks_summary_plot)
  ggsave(filename = "KS_test_summary_2n_vs_3n.pdf",
         plot = ks_summary_plot,
         width = 10,
         height = 6,
         dpi = 300)
}

cat("\nSaved Kolmogorov-Smirnov test results and plots\n")

stat=compare_means(Count~ploidy, data = foci_number, method="t.test", p.adjust.method = "BH",group.by = "treatment")
ggbarplot(foci_number, x= "treatment", y= "Count",label = T,lab.nb.digits = 3,lab.vjust = -2,fill="ploidy",palette = c("royalblue","red3"),width = 0.7, position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, size = 4.5,x = "treatment", label = "p.format",tip.length = 0.0025, y.position = c(1.5,6),hide.ns = F)+
  theme_classic()+ylab("Mean foci number")+xlab("Treatment")+
  theme(text = element_text(size = 13.5,face = "bold"), legend.position = "top")

