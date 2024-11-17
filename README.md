
# RISE Package

`RISE` is an R package designed for tackling the zero inflation problem in
sequence-based spatial transcriptomic (ST) data by leveraging
neighborhood-based regression models.

## Installation

You can install RISE directly from Github

```r
devtools::install_github("YSPeng1225/RISE")
library(RISE)
```

or you can first download RISE and install it manually in your local computer

``` r
install.packages("path/to/RISE", repos = NULL, type = "source")
library(RISE)
```

## Example Data

To help users understand the data structure, RISE includes an example
dataset, `spot_example`, with two elements:

- `spot_express`: a sparse matrix of gene expression values with genes
  as rows and spots as columns.

- `coordinates`: a matrix containing spatial x and y coordinates for spots.

Load and print the example data:

``` r
data(spot_example)
print(spot_example$spot_express[1:10,1:10])
print(spot_example$coordinates[1:10,])
# normalization
norm_data <- log2(t(t(spot_example$spot_express)/colSums(spot_example$spot_express) * median(colSums(spot_example$spot_express))) + 1)

```

## RISE Usage

RISE mainly relies on two functions, `FindNeighbors` and `RISE`, with
the previous focusing on the neighborhood selection and the latter
aiming to imputing zeros using regression models.

In `FindNeighbors` function, users can adjust the neighbor number
of spots through the parameter `thresh` (range from 0 to 1), with
larger `thresh` including more neighbors. `FindNeighbors` returns a
Boolean matrix containing the neighbors' indicators, which facilitates the
direct application in `RISE`. 

``` r
# neighbors identification
neighbors <- FindNeighbors(spot_example$coordinates, thresh = 0.3)

# to summarize the number of neighbors
summary(rowSums(neighbors))
hist(rowSums(neighbors), breaks = 200)

# to obtain the neighbors' names
spots <- colnames(neighbors)
neighbor_names <- apply(neighbors, 1, function(x) spots[x])
neighbor_names[["spot29"]]
```

Once the value of `thresh` is determined by the
user, the `RISE` function can be used to obtain the imputed gene expression
matrix.

``` r
spot_imputed <- RISE(spot_example$spot_express, spot_example$coordinates, thresh = 0.3)
# after imputation
spot_imputed[1:5,1:5]
# before imputation
norm_data[1:5,1:5]
```
