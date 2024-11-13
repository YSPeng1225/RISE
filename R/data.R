#' Example Spot Expression Data and Coordinates
#'
#' This example dataset contains a high-sparsity gene expression matrix (`spot_express`)
#' and corresponding 2D spatial coordinates (`coordinates`). The data is generated
#' using the following R code:
#'
#' ```
#' p <- 100  # number of genes
#' n <- 50   # number of cells
#'
#' # Generate gene expression matrix with high sparsity
#' set.seed(123)
#' spot_express <- matrix(sample(c(0, rpois(p * n, lambda = 2)), p * n, replace = TRUE, prob = c(0.95, 0.05)),
#'                        nrow = p, ncol = n)
#' rownames(spot_express) <- paste0("gene", 1:p)
#' colnames(spot_express) <- paste0("spot", 1:n)
#'
#' # Generate corresponding coordinates of each spot
#' coordinates <- data.frame(x = round(runif(n, min = 0, max = 1), 2),
#'                           y = round(runif(n, min = 0, max = 1), 2))
#' rownames(coordinates) <- paste0("spot", 1:n)
#' ```
#'
#' @format A list with two elements:
#' \describe{
#'   \item{spot_express}{A `p x n` matrix of gene expression values with genes as rows and spots as columns.}
#'   \item{coordinates}{A `n x 2` matrix of x and y coordinates for each spot, corresponding to columns in `spot_express`.}
#' }
#' @source Generated using the code above.
"spot_example"
