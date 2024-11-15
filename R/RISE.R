#' Find Neighbors Function
#'
#' `FindNeighbors` identifies neighboring spots based on their 2D coordinates,
#' which is common in spatial transcriptomics (ST) data.
#' Returns a boolean matrix indicating neighbors for each spot.
#'
#' @param coordinates A nx2 matrix with x and y coordinates for each spot, where `n` is the number of spots.
#' @param thresh Side length (0 to 1) of the square used to define neighborhood range; larger values include more neighbors.
#' @return A nxn boolean matrix indicating neighbor relationships across spots.
#' @export
FindNeighbors <- function(coordinates, thresh = 0.05) {

  # Input checks
  if (!is.matrix(coordinates) && !is.data.frame(coordinates)) {
    stop("Error: 'coordinates' must be a matrix or data frame.")
  }

  if (ncol(coordinates) != 2) {
    stop("Error: 'coordinates' must have exactly two columns (x and y coordinates).")
  }

  if (any(is.na(coordinates))) {
    stop("Error: 'coordinates' contains missing values. Please remove or impute missing data.")
  }

  if (!is.numeric(thresh) || thresh <= 0 || thresh > 1) {
    stop("Error: 'thresh' must be a numeric value between 0 and 1.")
  }

  # Standardize coordinates to range [0, 1] for both x and y dimensions
  xmin <- min(coordinates[,1])
  xmax <- max(coordinates[,1])
  ymin <- min(coordinates[,2])
  ymax <- max(coordinates[,2])
  coordinates[,'x'] <- (coordinates[,1] - xmin)/(xmax - xmin)
  coordinates[,'y'] <- (coordinates[,2] - ymin)/(ymax - ymin)

  # Pre-allocate logical matrix for neighbors
  n <- nrow(coordinates)
  neighbors <- matrix(FALSE, n, n, dimnames = list(rownames(coordinates), rownames(coordinates)))

  # Calculate distance between all points using vectorized operation
  for (i in seq_len(n)) {
    neighbors[, i] <- abs(coordinates[, 'x'] - coordinates[i, 'x']) < thresh &
      abs(coordinates[, 'y'] - coordinates[i, 'y']) < thresh
  }

  # Remove self-neighbor entries
  diag(neighbors) <- FALSE

  return(neighbors)
}



#' RISE
#'
#' RISE uses predictions from a neighbor-based regression model to impute zeros in spatial transcriptomic (ST) data.
#' @param spot_express A pxn matrix of spot expression, with rows as genes and columns as spots.
#' @param coordinates An nx2 matrix with x and y coordinates of each spot, where rows in `coordinates` align with columns in `spot_express`.
#' @param thresh Side length (0 to 1) of the square used to identify neighboring spots; larger values include more neighbors.
#' @param norm A boolean value indicating if `spot_express` needs normalization. The default is TRUE, which normalizes `spot_express` by the scaling factor method.
#' @param cores Number of CPU cores to use for parallel processing. The default is 2.
#' @return A pxn matrix of spot expression with zeros imputed by RISE.
#' @export
RISE <- function(spot_express, coordinates, thresh = 0.05, norm = TRUE, cores = 2) {

  library(parallel)
  if (ncol(coordinates) != 2) {
    stop("Error: 'coordinates' must be an nx2 matrix with x and y coordinates.")
  }
  if (ncol(spot_express) != nrow(coordinates)) {
    stop("Error: Number of spots in 'spot_express' and 'coordinates' must match.")
  }
  if (!is.numeric(thresh) || thresh <= 0 || thresh > 1) {
    stop("Error: 'thresh' must be a numeric value between 0 and 1.")
  }

  # Normalization (if not already normalized)
  if (norm) {
    scaling_factor <- median(colSums(spot_express))
    spot_express <- log2(t(t(spot_express) / colSums(spot_express) * scaling_factor) + 1)
  }

  # Neighbors identification and expression prediction
  neigh_mat <- FindNeighbors(coordinates, thresh)
  pred_express <- matrix(NA, nrow = nrow(spot_express), ncol = ncol(spot_express),
                         dimnames = dimnames(spot_express))
  cl <- parallel::makeCluster(cores)  # Create a cluster with the specified number of cores
  clusterExport(cl, varlist = c("spot_express", "neigh_mat", "impute_spot"), envir = environment())
  pred_list <- parLapply(cl, seq_len(ncol(spot_express)), impute_spot)
  stopCluster(cl)

  # Impute zeros in original expression matrix with predicted values
  pred_express <- do.call(cbind, pred_list)
  spot_express[spot_express == 0] <- pred_express[spot_express == 0]

  return(spot_express)
}


# Define function for each spot's imputation
impute_spot <- function(spot) {
  neighbor_spots <- which(neigh_mat[, spot])

  if (length(neighbor_spots) > 0) {
    X <- spot_express[, neighbor_spots, drop = FALSE]
    y <- spot_express[, spot]

    # Fit model and predict values
    model <- lm(y ~ ., data = as.data.frame(X))
    predict(model, newdata = as.data.frame(X))
  } else {
    rep(NA, nrow(spot_express))  # Return NA if no neighbors
  }
}
