#' Clustering function
#'
#' This function performs hierarchical clustering, performs feature importance
#' based on random forest model training and enrichment analysis most importance features.
#' It saves all clustering, features importance and enrichmnent results and outputs for downstream analysis.
#'
#' @param omic_data Transcriptomic data frame.
#' @param class_variable Stringname representing the class in the data.
#' @param changed_clinic_data Clinic dataset with MLASDO changed class variable.
#' @param id_column String name representing the sample ids column name.
#' @param saving_name String name used for naming saved output files.
#' @param changed_class Name of the change to be specifically analyzed 'Control2Case' or 'Case2Control'.
#' @param class Secondary change class: 'Case' for 'Case2Control' or 'Control' for 'Control2Case'.
#' @param mtry Number of variables to split at each node in random forest.
#' @param n_cores Number of cores to use in random forest training.
#' @param seed Seed for reproducibility.
#' @param enrichment Logical flag indicating if enrichment analysis should be performed.
#'
#' @return Clustering, feature importance and enrichment RDS files in a specified directory.
#'
#' @export
clustering <- function(omic_data, class_variable, changed_clinic_data, id_column, saving_name, changed_class, class, mtry, n_cores, seed, enrichment) {

  ### FIRST PART: CLUSTERING ###

  # Create a working copy of the data and subset by relevant classes
  omic_data_clustering <- omic_data
  omic_data_clustering[[class_variable]] <- changed_clinic_data[[class_variable]]
  omic_data_clustering <- omic_data_clustering[omic_data_clustering[[class_variable]] %in% c(changed_class, class), ]

  # Remove the ID column
  omic_data_clustering <- omic_data_clustering[, -which(names(omic_data_clustering) %in% c(id_column))]

  # Separate predictors and response variables
  x <- omic_data_clustering[, -ncol(omic_data_clustering)]
  labels <- as.factor(omic_data_clustering[, ncol(omic_data_clustering)])

  # Scale data and compute distance matrix
  x_scaled <- scale(x)
  distance <- as.dist(x_scaled)

  # Perform hierarchical clustering
  hc <- hclust(distance)
  clusters <- cutree_dynamic(dendro = hc, dist_m = as.matrix(distance), min_cluster_size = 40)
  dend <- hc %>% as.dendrogram

  # Save clustering results
  dir_path <- paste(saving_name, "clustering", sep = "/")
  saveRDS(x_scaled, paste(dir_path, paste0(saving_name, "_x_", changed_class, ".rds", sep = ""), sep = "/"))
  saveRDS(labels, paste(dir_path, paste0(saving_name, "_labels_", changed_class, ".rds", sep = ""), sep = "/"))
  saveRDS(hc, paste(dir_path, paste0(saving_name, "_hc_", changed_class, ".rds", sep = ""), sep = "/"))
  saveRDS(clusters, paste(dir_path, paste0(saving_name, "_clusters_", changed_class, ".rds", sep = ""), sep = "/"))
  saveRDS(dend, paste(dir_path, paste0(saving_name, "_dend_", changed_class, ".rds", sep = ""), sep = "/"))

  ### SECOND PART: RANDOM FOREST MODELS ###

  rf_importances <- list()
  rf_mc <- list()

  # Reorder data according to the dendrogram
  ordered_indices <- order.dendrogram(dend)
  x_ordered <- x[ordered_indices, ]
  clusters_ordered <- clusters[ordered_indices]
  labels_ordered <- labels[ordered_indices]
  omic_data_clustering_ordered <- omic_data_clustering[ordered_indices, ]

  # Train Random Forest models for valid clusters
  valid_clusters <- clusters_ordered[clusters_ordered != 0]
  num_clusters <- length(unique(valid_clusters))

  for (cluster_num in 1:num_clusters) {

    # Create binary response variable for the cluster
    cluster_indices <- which(clusters_ordered == cluster_num)
    response <- rep(0, length(labels_ordered))
    response[cluster_indices] <- 1
    response[-cluster_indices] <- 0

    omic_data_cluster <- data.frame(response = factor(response),
                                    omic_data_clustering_ordered)
    omic_data_cluster <- omic_data_cluster[!is.na(omic_data_cluster$response), ]

    # Define predictors (X) and response (Y)
    x <- omic_data_cluster[, !names(omic_data_cluster) %in% c("response")]
    y <- as.factor(omic_data_cluster[, "response"])
    levels(y) <- c("out_cluster", "in_cluster")

    # Define a custom metric for balanced accuracy
    metricas <- function(data, lev = levels(as.factor(data$obs)), model = NULL) {
      c(
        ba = (sensitivity(data[, "pred"], data[, "obs"], positive = "in_cluster") +
                specificity(data[, "pred"], data[, "obs"], negative = "out_cluster")) / 2
      )
    }

    # Train control settings for Random Forest
    train_control <- train_control(method = "cv",
                                   number = 10,
                                   save_predictions = "all",
                                   class_probs = TRUE,
                                   sampling = "down",
                                   return_resamp = "all",
                                   allow_parallel = TRUE,
                                   summary_function = metricas)

    # Model grid and trees
    mygrid <- expand.grid(mtry = c(mtry), splitrule = c("gini"), min.node.size = c(1))
    trees <- c(500)

    # Train Random Forest models
    cl <- make_cluster(n_cores)
    register_do_parallel(cl)

    rf_cv_10 <- list()
    for (tree in trees) {
      rf_model <- train(
        x,
        y,
        method = "ranger",
        tr_control = train_control,
        tune_grid = mygrid,
        num.trees = tree,
        classification = TRUE,
        seed = seed,
        metric = "ba",
        importance = "permutation"
      )
      rf_cv_10[[paste(tree, "trees")]] <- rf_model
    }
    stop_cluster(cl)

    # Save model results and variable importance
    best_model_trees <- aggregate(metric_value ~ trees, data = datos_combinados, fun = max)
    best_config <- best_model_trees[which.max(best_model_trees$metric_value), "trees"]
    best_rf_model <- rf_cv_10[[best_config]]

    pred <- predict(best_rf_model, x)
    rf_mc[[paste0("cluster_", cluster_num)]] <- confusion_matrix(pred, y, positive = "in_cluster")

    importance_df <- var_imp(best_rf_model, scale = TRUE)$importance
    importance_df <- as.data.frame(importance_df)
    importance_df$variable <- rownames(importance_df)
    names(importance_df) <- c("importance", "variable")
    sorted_importance_df <- importance_df[order(-importance_df$importance), ]
    rf_importances[[paste0("cluster_", cluster_num)]] <- sorted_importance_df
  }

  # Save Random Forest importance and confusion matrices
  saveRDS(rf_importances, file = paste(dir_path, paste0(saving_name, "_importances_", changed_class, ".rds", sep = ""), sep = "/"))
  saveRDS(rf_mc, file = paste(dir_path, paste0(saving_name, "_cm_", changed_class, ".rds", sep = ""), sep = "/"))

  ### THIRD PART: ENRICHMENT ANALYSIS ###

  rf_enrichments <- list()
  rf_elbows <- list()
  rf_genes <- list()

  for (cluster_name in names(rf_importances)) {

    importance <- rf_importances[[cluster_name]]
    sorted_importance <- importance[order(-importance$importance), ]

    # Identify the elbow point for variable importance
    x_values <- 1:nrow(sorted_importance)
    y_values <- sorted_importance$importance
    importance_elbow <- data.frame(x = x_values, y = y_values)
    elbow <- find_curve_elbow(importance_elbow)
    all_genes <- sorted_importance$variable[1:elbow]

    # Perform enrichment analysis
    enrich <- gprofiler2::gost(list(all_genes),
                               correction_method = "g_scs",
                               sources = c("go", "kegg", "reac", "hp"),
                               organism = "hsapiens",
                               exclude_iea = TRUE)$result

    if (!is.null(enrich)) {
      enrich <- enrich %>%
        mutate(gene_ratio = intersection_size / term_size,
               log_p_value = -log10(p_value))
    }

    cluster_num <- gsub("cluster_", "", cluster_name)
    rf_enrichments[[paste0("cluster_", cluster_num)]] <- enrich
    rf_elbows[[paste0("cluster_", cluster_num)]] <- elbow
    rf_genes[[paste0("cluster_", cluster_num)]] <- all_genes
  }

  # Save enrichment analysis results
  saveRDS(rf_enrichments, file = paste(dir_path, paste0(saving_name, "_enrichments_", changed_class, ".rds", sep = ""), sep = "/"))
  saveRDS(rf_elbows, file = paste(dir_path, paste0(saving_name, "_elbows_", changed_class, ".rds", sep = ""), sep = "/"))
  saveRDS(rf_genes, file = paste(dir_path, paste0(saving_name, "_genes_", changed_class, ".rds", sep = ""), sep = "/"))
}
