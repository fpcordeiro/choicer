# Halton draws for mixed logit

Create halton normal draws in appropriate format for mixed logit
estimation

## Usage

``` r
get_halton_normals(S, N, K_w)
```

## Arguments

- S:

  Number of draws for each choice situation

- N:

  number of choice situations

- K_w:

  dimension of random coefficients (number of columns in W matrix)

## Value

K_w x S x N array with halton standard normal draws

## Examples

``` r
draws <- get_halton_normals(S = 50, N = 10, K_w = 2)
dim(draws)  # 2 x 50 x 10
#> [1]  2 50 10
```
