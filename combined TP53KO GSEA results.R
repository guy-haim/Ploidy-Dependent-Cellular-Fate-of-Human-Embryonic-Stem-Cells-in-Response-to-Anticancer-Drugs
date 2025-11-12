library(stringr)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggfortify)
library(pheatmap)
library(data.table)
library(zoo)
library(gplots)
library(PCAtools)
library(reshape2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(FSA)
library(ggsignif)
library(ggpubr)

setwd(choose.dir())


### heatmap of GSEA terms in treated diploid and triploid Vs their untreated controls:
GSEA_TP53KO_2n_vs_TP53KO_3n_5FU_down <- read.csv("All downregulated GSEA results of TP53KO 5FU-treated triploid vs diploid hESCs.csv")
GSEA_TP53KO_2n_vs_TP53KO_3n_5FU_up <- read.csv("All upregulated GSEA results of TP53KO 5FU-treated triploid vs diploid hESCs.csv")
GSEA_WT_3n_vs_TP53KO_3n_5FU_down <- read.csv("All downregulated GSEA results of both WT and TP53-KO triploid hESCs treated with 5FU.csv")
GSEA_WT_3n_vs_TP53KO_3n_5FU_up <- read.csv("All upregulated GSEA results of both WT and TP53-KO triploid hESCs treated with 5FU.csv")
GSEA_WT_2n_vs_TP53KO_2n_5FU_down <- read.csv("All downregulated GSEA results of both WT and TP53-KO diploid hESCs treated with 5FU.csv")
GSEA_WT_2n_vs_TP53KO_2n_5FU_up <- read.csv("All upregulated GSEA results of both WT and TP53-KO diploid hESCs treated with 5FU.csv")
GSEA_TP53KO_2n_vs_TP53KO_3n_Bleo_down <- read.csv("All downregulated GSEA results of TP53KO Bleo-treated triploid vs diploid hESCs.csv")
GSEA_TP53KO_2n_vs_TP53KO_3n_Bleo_up <- read.csv("All upregulated GSEA results of TP53KO Bleo-treated triploid vs diploid hESCs.csv")
GSEA_WT_3n_vs_TP53KO_3n_Bleo_down <- read.csv("All downregulated GSEA results of both WT and TP53-KO triploid hESCs treated with Bleo.csv")
GSEA_WT_3n_vs_TP53KO_3n_Bleo_up <- read.csv("All upregulated GSEA results of both WT and TP53-KO triploid hESCs treated with Bleo.csv")
GSEA_WT_2n_vs_TP53KO_2n_Bleo_down <- read.csv("All downregulated GSEA results of both WT and TP53-KO diploid hESCs treated with Bleo.csv")
GSEA_WT_2n_vs_TP53KO_2n_Bleo_up <- read.csv("All upregulated GSEA results of both WT and TP53-KO diploid hESCs treated with Bleo.csv")
GSEA_TP53KO_2n_vs_TP53KO_3n_Cis_down <- read.csv("All downregulated GSEA results of TP53KO Cis-treated triploid vs diploid hESCs.csv")
GSEA_TP53KO_2n_vs_TP53KO_3n_Cis_up <- read.csv("All upregulated GSEA results of TP53KO Cis-treated triploid vs diploid hESCs.csv")
GSEA_WT_3n_vs_TP53KO_3n_Cis_down <- read.csv("All downregulated GSEA results of both WT and TP53-KO triploid hESCs treated with Cis.csv")
GSEA_WT_3n_vs_TP53KO_3n_Cis_up <- read.csv("All upregulated GSEA results of both WT and TP53-KO triploid hESCs treated with Cis.csv")
GSEA_WT_2n_vs_TP53KO_2n_Cis_down <- read.csv("All downregulated GSEA results of both WT and TP53-KO diploid hESCs treated with Cis.csv")
GSEA_WT_2n_vs_TP53KO_2n_Cis_up <- read.csv("All upregulated GSEA results of both WT and TP53-KO diploid hESCs treated with Cis.csv")
GSEA_TP53KO_2n_vs_TP53KO_3n_Pacli_down <- read.csv("All downregulated GSEA results of TP53KO Pacli-treated triploid vs diploid hESCs.csv")
GSEA_TP53KO_2n_vs_TP53KO_3n_Pacli_up <- read.csv("All upregulated GSEA results of TP53KO Pacli-treated triploid vs diploid hESCs.csv")
GSEA_WT_3n_vs_TP53KO_3n_Pacli_down <- read.csv("All downregulated GSEA results of both WT and TP53-KO triploid hESCs treated with Pacli.csv")
GSEA_WT_3n_vs_TP53KO_3n_Pacli_up <- read.csv("All upregulated GSEA results of both WT and TP53-KO triploid hESCs treated with Pacli.csv")
GSEA_WT_2n_vs_TP53KO_2n_Pacli_down <- read.csv("All downregulated GSEA results of both WT and TP53-KO diploid hESCs treated with Pacli.csv")
GSEA_WT_2n_vs_TP53KO_2n_Pacli_up <- read.csv("All upregulated GSEA results of both WT and TP53-KO diploid hESCs treated with Pacli.csv")



# Function to combine GSEA enrichment and depletion results
# Selects the best result for each pathway based on p-value and FDR
combine_gsea_results <- function(enrichment_df, depletion_df, 
                                 primary_metric = "pval", 
                                 secondary_metric = "padj",
                                 add_source_column = TRUE) {
  
  # Check if both dataframes have the same pathways
  if (!identical(sort(enrichment_df$pathway), sort(depletion_df$pathway))) {
    stop("Error: Enrichment and depletion dataframes must contain the same pathways")
  }
  
  # Ensure both dataframes are ordered by pathway for proper comparison
  enrichment_df <- enrichment_df[order(enrichment_df$pathway), ]
  depletion_df <- depletion_df[order(depletion_df$pathway), ]
  
  # Initialize combined dataframe
  combined_df <- enrichment_df[0, ]  # Empty dataframe with same structure
  
  if (add_source_column) {
    combined_df$source <- character(0)
  }
  
  # Compare each pathway and select the best result
  for (i in 1:nrow(enrichment_df)) {
    enrich_row <- enrichment_df[i, ]
    deplete_row <- depletion_df[i, ]
    
    # Handle NA values - if one has NA and other doesn't, choose the non-NA
    enrich_primary_na <- is.na(enrich_row[[primary_metric]])
    deplete_primary_na <- is.na(deplete_row[[primary_metric]])
    
    if (enrich_primary_na && !deplete_primary_na) {
      selected_row <- deplete_row
      source_label <- "depletion"
    } else if (!enrich_primary_na && deplete_primary_na) {
      selected_row <- enrich_row
      source_label <- "enrichment"
    } else if (enrich_primary_na && deplete_primary_na) {
      # Both are NA, use enrichment as default
      selected_row <- enrich_row
      source_label <- "enrichment"
    } else {
      # Both have valid values, compare primary metric
      if (enrich_row[[primary_metric]] < deplete_row[[primary_metric]]) {
        selected_row <- enrich_row
        source_label <- "enrichment"
      } else if (enrich_row[[primary_metric]] > deplete_row[[primary_metric]]) {
        selected_row <- deplete_row
        source_label <- "depletion"
      } else {
        # Primary metrics are equal, use secondary metric
        enrich_secondary_na <- is.na(enrich_row[[secondary_metric]])
        deplete_secondary_na <- is.na(deplete_row[[secondary_metric]])
        
        if (enrich_secondary_na && !deplete_secondary_na) {
          selected_row <- deplete_row
          source_label <- "depletion"
        } else if (!enrich_secondary_na && deplete_secondary_na) {
          selected_row <- enrich_row
          source_label <- "enrichment"
        } else if (!enrich_secondary_na && !deplete_secondary_na) {
          if (enrich_row[[secondary_metric]] <= deplete_row[[secondary_metric]]) {
            selected_row <- enrich_row
            source_label <- "enrichment"
          } else {
            selected_row <- deplete_row
            source_label <- "depletion"
          }
        } else {
          # Both secondary metrics are NA, default to enrichment
          selected_row <- enrich_row
          source_label <- "enrichment"
        }
      }
    }
    
    # Add source information if requested
    if (add_source_column) {
      selected_row$source <- source_label
    }
    
    # Add to combined dataframe
    combined_df <- rbind(combined_df, selected_row)
  }
  
  return(combined_df)
} 

GSEA_WT_2n_vs_TP53KO_2n_5FU <- combine_gsea_results(GSEA_WT_2n_vs_TP53KO_2n_5FU_up,GSEA_WT_2n_vs_TP53KO_2n_5FU_down)
GSEA_WT_3n_vs_TP53KO_3n_5FU <- combine_gsea_results(GSEA_WT_3n_vs_TP53KO_3n_5FU_up,GSEA_WT_3n_vs_TP53KO_3n_5FU_down)
GSEA_TP53KO_2n_vs_TP53KO_3n_5FU <- combine_gsea_results(GSEA_TP53KO_2n_vs_TP53KO_3n_5FU_up,GSEA_TP53KO_2n_vs_TP53KO_3n_5FU_down)
GSEA_WT_2n_vs_TP53KO_2n_Bleo <- combine_gsea_results(GSEA_WT_2n_vs_TP53KO_2n_Bleo_up,GSEA_WT_2n_vs_TP53KO_2n_Bleo_down)
GSEA_WT_3n_vs_TP53KO_3n_Bleo <- combine_gsea_results(GSEA_WT_3n_vs_TP53KO_3n_Bleo_up,GSEA_WT_3n_vs_TP53KO_3n_Bleo_down)
GSEA_TP53KO_2n_vs_TP53KO_3n_Bleo <- combine_gsea_results(GSEA_TP53KO_2n_vs_TP53KO_3n_Bleo_up,GSEA_TP53KO_2n_vs_TP53KO_3n_Bleo_down)
GSEA_WT_2n_vs_TP53KO_2n_Cis <- combine_gsea_results(GSEA_WT_2n_vs_TP53KO_2n_Cis_up,GSEA_WT_2n_vs_TP53KO_2n_Cis_down)
GSEA_WT_3n_vs_TP53KO_3n_Cis <- combine_gsea_results(GSEA_WT_3n_vs_TP53KO_3n_Cis_up,GSEA_WT_3n_vs_TP53KO_3n_Cis_down)
GSEA_TP53KO_2n_vs_TP53KO_3n_Cis <- combine_gsea_results(GSEA_TP53KO_2n_vs_TP53KO_3n_Cis_up,GSEA_TP53KO_2n_vs_TP53KO_3n_Cis_down)
GSEA_WT_2n_vs_TP53KO_2n_Pacli <- combine_gsea_results(GSEA_WT_2n_vs_TP53KO_2n_Pacli_up,GSEA_WT_2n_vs_TP53KO_2n_Pacli_down)
GSEA_WT_3n_vs_TP53KO_3n_Pacli <- combine_gsea_results(GSEA_WT_3n_vs_TP53KO_3n_Pacli_up,GSEA_WT_3n_vs_TP53KO_3n_Pacli_down)
GSEA_TP53KO_2n_vs_TP53KO_3n_Pacli <- combine_gsea_results(GSEA_TP53KO_2n_vs_TP53KO_3n_Pacli_up,GSEA_TP53KO_2n_vs_TP53KO_3n_Pacli_down)


##############################################################################
combined_GSEA <- Reduce(function(x, y) merge(x, y, all=T, by = "pathway"),
                        list(GSEA_TP53KO_2n_vs_TP53KO_3n_5FU[,c(2,4,7)],
                             GSEA_TP53KO_2n_vs_TP53KO_3n_Bleo[,c(2,4,7)],
                             GSEA_TP53KO_2n_vs_TP53KO_3n_Cis[,c(2,4,7)],
                             GSEA_TP53KO_2n_vs_TP53KO_3n_Pacli[,c(2,4,7)]))
colnames(combined_GSEA) <- c("pathway","5FU_p.adj","5FU_NES","Bleo_p.adj","Bleo_NES","Cis_p.adj","Cis_NES","Pacli_p.adj","Pacli_NES")

relevant_terms <- c("P53","TP53","G2","G1","_M_","_S_","PHASE","CHECKPOINT","GTSE","DNA","MYC","E2F","MTOR","OXIDATIVE_PHOSPHORYLATION","SENESCENCE",
                    "PATTERN_SPECIFICATION","SKELETAL_SYSTEM_DEVELOPMENT","UV","APOP","DAMAGE","MITOTIC","MITOSIS","SPINDLE","PI3","AKT","P16","CDKN1A",
                    "CDKN2A","CHROMOSOME_SEGREGATION","REPAIR","PCNA","P21","P27","TNFA","NECRO","DOUBLE_STRAND_BREAK","MICROTUBULE","PHASE_TRANSITION")



##### generate a unified GSEA terms list and display the relevant pathways in 3n vs 2n TP53-KO cells:
# Function to create summary of GSEA pathways by treatment
create_gsea_pathway_summary <- function(combined_gsea_df, pathway_groups) {
  
  # Initialize results dataframe
  summary_results <- data.frame()
  
  # Get unique treatments from column names (assuming format like "Treatment_5FU", "Treatment_Bleo", etc.)
  treatment_cols <- colnames(combined_gsea_df)
  treatments <- c("5FU", "Bleo", "Cis", "Pacli")  # Define your treatments
  
  # Process each pathway group
  for (group_name in names(pathway_groups)) {
    
    cat("Processing pathway group:", group_name, "\n")
    
    # Get pathways for this group
    group_pathways <- pathway_groups[[group_name]]
    
    # Find rows that match any of the pathways in this group
    matching_rows <- combined_gsea_df[grepl(paste(group_pathways, collapse = "|"), 
                                            combined_gsea_df$pathway, ignore.case = TRUE), ]
    
    if (nrow(matching_rows) == 0) {
      cat("  No pathways found for group:", group_name, "\n")
      next
    }
    
    cat("  Found", nrow(matching_rows), "pathways for group:", group_name, "\n")
    
    # Process each treatment
    for (treatment in treatments) {
      
      # Find columns for this treatment (padj and NES)
      padj_col <- paste0(treatment, "_p.adj")  # Adjust column naming as needed
      nes_col <- paste0(treatment, "_NES")    # Adjust column naming as needed
      
      # Alternative column naming patterns (adjust based on your actual column names)
      if (!padj_col %in% colnames(matching_rows)) {
        # Try alternative naming patterns
        potential_padj_cols <- grep(paste0(treatment, ".*padj|padj.*", treatment), 
                                    colnames(matching_rows), ignore.case = TRUE, value = TRUE)
        potential_nes_cols <- grep(paste0(treatment, ".*NES|NES.*", treatment), 
                                   colnames(matching_rows), ignore.case = TRUE, value = TRUE)
        
        if (length(potential_padj_cols) > 0) padj_col <- potential_padj_cols[1]
        if (length(potential_nes_cols) > 0) nes_col <- potential_nes_cols[1]
      }
      
      # Check if columns exist
      if (!padj_col %in% colnames(matching_rows) || !nes_col %in% colnames(matching_rows)) {
        cat("  Warning: Columns not found for treatment", treatment, "\n")
        cat("  Looking for:", padj_col, "and", nes_col, "\n")
        cat("  Available columns:", paste(colnames(matching_rows), collapse = ", "), "\n")
        next
      }
      
      # Get valid (non-NA) padj values for this treatment
      valid_rows <- !is.na(matching_rows[[padj_col]])
      
      if (sum(valid_rows) == 0) {
        cat("  No valid padj values for treatment", treatment, "in group", group_name, "\n")
        next
      }
      
      # Find row with minimum padj value
      min_padj_idx <- which.min(matching_rows[[padj_col]][valid_rows])
      # Get the actual row index in the original dataframe
      actual_idx <- which(valid_rows)[min_padj_idx]
      best_row <- matching_rows[actual_idx, ]
      
      # Extract the best values
      best_padj <- best_row[[padj_col]]
      best_nes <- best_row[[nes_col]]
      best_pathway <- best_row$pathway
      
      # Create summary row
      summary_row <- data.frame(
        General_Pathway = group_name,
        Treatment = treatment,
        Best_Pathway = best_pathway,
        Min_padj = best_padj,
        Corresponding_NES = best_nes,
        stringsAsFactors = FALSE
      )
      
      # Add to results
      summary_results <- rbind(summary_results, summary_row)
      
      cat("  ", treatment, "- Best pathway:", substr(best_pathway, 1, 50), "..., padj:", 
          round(best_padj, 6), ", NES:", round(best_nes, 3), "\n")
    }
  }
  
  return(summary_results)
}

# Define pathway groups with common themes
pathway_groups <- list(
  "Cell_Cycle_G1/S_G2/M_Phase_Transition" = c("REACTOME_MITOTIC_G1_PHASE_AND_G1_S_TRANSITION", "REACTOME_M_PHASE","GOBP_CELL_CYCLE_G1_S_PHASE_TRANSITION", 
                                              "REACTOME_G1_S_SPECIFIC_TRANSCRIPTION", "HALLMARK_G2M_CHECKPOINT", "GOBP_REGULATION_OF_MITOTIC_CELL_CYCLE",
                                              "GOBP_REGULATION_OF_MITOTIC_CELL_CYCLE_PHASE_TRANSITION","GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION",
                                              "GOBP_REGULATION_OF_CELL_CYCLE_PHASE_TRANSITION","GOBP_REGULATION_OF_CELL_CYCLE_G2_M_PHASE_TRANSITION",
                                              "GOBP_REGULATION_OF_CELL_CYCLE_G1_S_PHASE_TRANSITION","GOBP_POSITIVE_REGULATION_OF_MITOTIC_CELL_CYCLE_PHASE_TRANSITION",
                                              "GOBP_POSITIVE_REGULATION_OF_MITOTIC_CELL_CYCLE","GOBP_POSITIVE_REGULATION_OF_CELL_CYCLE_PHASE_TRANSITION",
                                              "GOBP_NEGATIVE_REGULATION_OF_MITOTIC_CELL_CYCLE","GOBP_MITOTIC_CELL_CYCLE_PHASE_TRANSITION",
                                              "REACTOME_MITOTIC_PROPHASE","REACTOME_MITOTIC_PROMETAPHASE","REACTOME_CONDENSATION_OF_PROPHASE_CHROMOSOMES"),
  "DNA_Repair" = c("WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK", "WP_DNA_IRDAMAGE_AND_CELLULAR_RESPONSE_VIA_ATR", "REACTOME_HOMOLOGY_DIRECTED_REPAIR", 
                   "REACTOME_DNA_DOUBLE_STRAND_BREAK_RESPONSE", "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR","KEGG_MISMATCH_REPAIR","GOBP_DNA_INTEGRITY_CHECKPOINT_SIGNALING",
                   "KEGG_MEDICUS_REFERENCE_MISMATCH_REPAIR","HALLMARK_DNA_REPAIR","GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
                   "GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR","GOBP_REGULATION_OF_DNA_REPAIR","GOBP_POSITIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR",
                   "GOBP_POSITIVE_REGULATION_OF_DNA_REPAIR","GOBP_NUCLEOTIDE_EXCISION_REPAIR","GOBP_INTERSTRAND_CROSS_LINK_REPAIR","GOBP_DOUBLE_STRAND_BREAK_REPAIR"),
  "DNA_Replication" = c("REACTOME_DNA_REPLICATION", "GOBP_DNA_TEMPLATED_DNA_REPLICATION", "GOBP_DNA_REPLICATION"),
  "Senescence" = c("WP_SENESCENCEASSOCIATED_SECRETORY_PHENOTYPE_SASP", "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP", "REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE", 
                   "REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE","REACTOME_CELLULAR_SENESCENCE"),
  "Pattern_Specification_Process" = c("GOBP_PATTERN_SPECIFICATION_PROCESS"),
  "Microtubule_Organization/Spindle_Function" = c("GOBP_REGULATION_OF_MICROTUBULE_CYTOSKELETON_ORGANIZATION", "GOBP_MICROTUBULE_POLYMERIZATION_OR_DEPOLYMERIZATION", "GOCC_MICROTUBULE",
                                 "GOBP_MICROTUBULE_DEPOLYMERIZATION","GOCC_SPINDLE","WP_MICROTUBULE_CYTOSKELETON_REGULATION")
)

# Example usage:
# Assuming your combined_gsea_df has columns like:
# pathway, 5FU_padj, 5FU_NES, Bleo_padj, Bleo_NES, Cis_padj, Cis_NES, Pacli_padj, Pacli_NES

# Create the summary
gsea_summary <- create_gsea_pathway_summary(combined_GSEA, pathway_groups)

# View results
ggplot(data = gsea_summary, aes(General_Pathway, -log2(Min_padj)*sign(Corresponding_NES), group = Treatment)) +
  geom_col(aes(fill = Treatment), width = 0.5, position = position_dodge(0.6)) +
  coord_flip() + 
  scale_fill_manual(values = c("goldenrod1", "orange", "tan3", "chocolate3")) +
  labs(x = "Pathway", y = expression(-log[2](FDR))) + 
  theme_classic() + 
  ggtitle("TP53-KO 3n/2n Unified GSEA") +
  scale_x_discrete(limits = c("Microtubule_Organization/Spindle_Function",
                              "DNA_Replication",
                              "Cell_Cycle_G1/S_G2/M_Phase_Transition",
                              "DNA_Repair",
                              "Pattern_Specification_Process")) +
  geom_hline(yintercept = -log2(0.05), linetype = 2, color = "grey", size = 0.25) +
  geom_hline(yintercept = log2(0.05), linetype = 2, color = "grey", size = 0.25) +
  geom_hline(yintercept = -log2(1), linetype = 1, color = "grey", size = 0.25) +
  theme(axis.text.y = element_text(size = 10)) +
  scale_y_continuous(limits = c(-30, 30),n.breaks = 6)
# ggsave(filename="GSEA_relevant_general_terms_in_3n_TP53KO_bars.pdf")
