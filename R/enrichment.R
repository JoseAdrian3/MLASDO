#' Enrichment analysis
#'
#' This function performs enrichment analysis based on the most important genes
#' in the dataset using the elbow rule method. Omic may be transcriptomic,
#' proteomic or metabolomic.
#'
#' @param importance Data frame containing variables and their importance scores in $importance
#' @param saving_name String for the saving directory and file prefix.
#' @param text String to append to file names.
#' @param enrichment String with the enrichment type: "transcriptomic" or "metabolomic".
#' @param omic_data_changed Data frame of omics data with the changed condition.
#' @param omic_data_maintained Data frame of omics data with the original condition.
#' @param bool String indicating condition: "changed" or other.
#'
#' @return Saves enrichment analysis results and elbow data.
#' @export
enrichment <- function(importance, saving_name, text, enrichment, omic_data_changed, omic_data_maintained, bool) {

  # Create directory path for saving results
  dir_path <- paste(saving_name, "enrichment", sep = "/")

  if (!is.null(importance)) {

    # Sort importance by descending order
    sorted_importance <- importance[order(-importance$importance), ]

    # Create x and y values for elbow point calculation
    x_values <- 1:nrow(sorted_importance)
    y_values <- sorted_importance$importance

    importance_elbow <- data.frame(
      x = x_values,
      y = y_values
    )

    # Find elbow point (default to 50 if not found)
    elbow <- find_curve_elbow(importance_elbow)
    if (is.null(elbow) || length(elbow) == 0) {
      elbow <- 50
    }

    # Select top genes based on elbow
    all_genes <- sorted_importance$variable[1:elbow]

    # Calculate mean values for changed and maintained conditions
    mean_changed <- colMeans(omic_data_changed[, -((ncol(omic_data_changed) - 1):ncol(omic_data_changed))])
    mean_maintained <- colMeans(omic_data_maintained[, -((ncol(omic_data_maintained) - 1):ncol(omic_data_maintained))])

    # Subset mean values for selected genes
    mean_changed_genes <- mean_changed[all_genes]
    mean_maintained_genes <- mean_maintained[all_genes]

    # Identify genes based on the condition ("changed" or "maintained")
    if (bool == "changed") {
      genes <- all_genes[mean_changed_genes < mean_maintained_genes]
    } else {
      genes <- all_genes[mean_changed_genes > mean_maintained_genes]
    }

    if (enrichment == "transcriptomic") {

      # Perform transcriptomic enrichment analysis using gprofiler2
      enrich <- gprofiler2::gost(list(genes),
                                 correction_method = "g_SCS",
                                 sources = c("GO", "KEGG", "REAC", "HP"),
                                 organism = "hsapiens",
                                 exclude_iea = TRUE)$result

      if (!is.null(enrich)) {
        enrich <- enrich %>%
          mutate(gene_ratio = intersection_size / term_size,
                 log_p_value = -log10(p_value))
      }

      # Save enrichment results to a file
      saveRDS(enrich, file.path(dir_path, paste0(saving_name, "_enrichment_", text, "_", bool, ".rds")))

    } else if (enrichment == "metabolomic") {

      # Perform metabolomic enrichment analysis
      m_set <- InitDataObjects("conc", "msetora", FALSE)

      m_set <- Setup.MapData(m_set, genes)
      m_set <- CrossReferencing(m_set, "hmdb")
      m_set <- CreateMappingResultTable(m_set)
      m_set <- SetMetabolomeFilter(m_set, FALSE)
      m_set <- SetCurrentMsetLib(m_set, "smpdb_pathway", 0)
      m_set <- CalculateHyperScore(m_set)

      # Plot results and save to a file
      m_set <- PlotORA(m_set, file.path(dir_path, paste0(saving_name, "_enrichment_", text, "_", bool)), "bar", "png", 72, width = NA)

    } else {
      enrich <- NULL
    }

    # Save elbow data to a file
    saveRDS(elbow, file.path(dir_path, paste0(saving_name, "_elbow_", text, ".rds")))
  }
}
