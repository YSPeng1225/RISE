
# RISE Package

`RISE` is an R package designed for tackling zero inflation in
sequence-based spatial transcriptomic (ST) data by leveraging
neighborhood-based regression model.

## Installation

To install RISE from Github:

```r
devtools::install_github("YSPeng1225/RISE")
library(RISE)
```

To install RISE manually:

``` r
install.packages("path/to/RISE", repos = NULL, type = "source")
library(RISE)
```

## Example Data

To help users understand the data structure, RISE includes an example
dataset, `spot_example`, with two elements:

- `spot_express`: A sparse matrix of gene expression values with genes
  as rows and spots as columns.

- `coordinates`: A matrix containing x and y coordinates for each spot.

Load and print the example data:

``` r
data(spot_example)
print(spot_example$spot_express[1:10,1:10])
print(spot_example$coordinates[1:10,])
```

## RISE Usage

RISE mainly relies on two functions, `FindNeighbors` and `RISE`, with
the first focusing on neighboring spot identification and the second
imputing the zeros using neighbor-based regression predictions.

The role of the `FindNeighbors` function is to let users adjust the number
of spot’s neighbors through the parameter `thresh` (range from 0 to 1), with
larger `thresh` including more neighbors. `FindNeighbors` returns a
logistic matrix containing the neighbors' indicator, which facilitates the
direct application in `RISE`. However, there are simple examples to
further explore the neighboring spots' information.

``` r
# neighbors identification
neighbors <- FindNeighbors(spot_example$coordinates, thresh = 0.3)

# To summarize the number of neighbors
summary(rowSums(neighbors))
hist(rowSums(neighbors), breaks = 200)

# to obtain the neighbors' names
spots <- colnames(neighbors)
neighbor_names <- apply(neighbors, 1, function(x) spots[x])
neighbor_names[["spot29"]]
```

Once the value of `thresh` to identify neighbors is determined by the
user, the `RISE` function can be used to gain the imputed gene expression
matrix.

``` r
spot_imputed <- RISE(spot_example$spot_express, spot_example$coordinates, thresh = 0.3)
# After imputation
spot_imputed[1:5,1:5]
# Before imputation
spot_example$spot_express[1:5,1:5]
```
