
<!-- README.md is generated from README.Rmd. Please edit that file -->

# linearalgebra

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/JamesFrech/linearalgebra.svg?branch=main)](https://travis-ci.com/JamesFrech/linearalgebra)
<!-- badges: end -->

The goal of linearalgebra is to provide various functions to perform
tasks dealing with linear algebra. Included are functions that can
return the rref of a matrix, find the determinant of a matrix, and more.

## Installation

You can install the released version of linearalgebra from github with:

``` r
install_github("JamesFrech/linearalgebra")
```

## Usage

Here are some examples of the functions included in this package:

The first function returns the transpose of a matrix.

``` r
transpose <- function(matrix){
  t_matrix <- matrix(nrow = ncol(matrix), ncol = nrow(matrix))
  # make the rows of the matrix the columns of the transpose
  for(i in 1:nrow(matrix)){
    t_matrix[,i] <- matrix[i,]
  }
  return(t_matrix)
}

transpose(matrix(1:6, nrow = 2))
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    3    4
#> [3,]    5    6
```

This function multiplies matrices together.

``` r
m1 <- matrix(1:6, ncol = 2)
m2 <-  matrix(8:13, nrow = 2)
print(m1)
#>      [,1] [,2]
#> [1,]    1    4
#> [2,]    2    5
#> [3,]    3    6
print(m2)
#>      [,1] [,2] [,3]
#> [1,]    8   10   12
#> [2,]    9   11   13

multiply_matrices <- function(A, B){
  if(ncol(A) != nrow(B)){
    stop("Number of columns in matrix A must match number of rows in matrix B.")
  }
  result <- matrix(nrow = dim(A)[1], ncol = dim(B)[2])
  # Each i,j spot of the result is the sum of the i of the first matrix multiplied by each of the j's of the second
  for(i in 1:dim(A)[1]){
    for(j in 1:dim(B)[2]){
      result[i,j] <- sum(A[i,]*B[,j])
    }
  }
  return(result)
}

multiply_matrices(m1, m2)
#>      [,1] [,2] [,3]
#> [1,]   44   54   64
#> [2,]   61   75   89
#> [3,]   78   96  114
```

Functions row\_swap, row\_replace, and arrange\_rows are usable as well.
However they are just meant to be used when called in rref. Refer to
their documentation for examples on these functions.

``` r
row_swap <- function(matrix, row1, row2){
  # copy the first desired row as a temporary row to the end of the matrix
  matrix <- rbind(matrix, matrix[row1,])
  # make the first desired row a copy of the second desired row
  matrix[row1,] <- matrix[row2,]
  # make the second desired row a copy of the temporary row created
  matrix[row2,] <- matrix[nrow(matrix),]
  # select all rows of the matrix other than the temporary (last) row
  matrix <- matrix[1:(nrow(matrix)-1),]
  return(matrix)
}

row_replace <- function(matrix, row1, row2){
  position <- 1
  # find the first position where neither row is equal to 0
  while((matrix[row1, position] == 0) | (matrix[row2, position] == 0)){
    position <- position + 1
  }
  # find the constant to multiply the first row by in order to cancel out the first term in the second row
  x <- -matrix[row2, position]/matrix[row1, position]
  # add the corresponding spot in the first row times the constant to each element of the second row
  for(j in seq_along(matrix[row2,])){
    matrix[row2, j] <- matrix[row2, j] + x*matrix[row1, j]
  }
  return(matrix)
}

arrange_rows <- function(matrix){
  column <- 1
  # first row is the first row that you will start with in each column to rearrange the order.
  first_row <- 1
  new_matrix <- matrix
  old_matrix <- new_matrix
  while(column <= ncol(matrix)){
    i <- first_row
    while(i < nrow(matrix)){
      if(new_matrix[i, column] < new_matrix[i + 1, column]){
        old_matrix <- new_matrix
        new_matrix <- row_swap(old_matrix, i, i + 1)
        if((i > 1) & (i > first_row)){
          i <- i - 1
        }
      }
      else{
        i <- i + 1
      }
    }
    for(i in first_row:nrow(new_matrix)){
      if(new_matrix[i, column] == 0){
        first_row <- i
        break
      }
    }
    column <- column + 1
  }
  return(new_matrix)
}

rref <- function(matrix){
  # arrange rows so that 0's are at bottom
  new_matrix <- arrange_rows(matrix)
  column <- 1
  first_row <- 0
  previous <- -1
  while(column <= ncol(new_matrix)){
    for(i in 1:nrow(new_matrix)){
      if(abs(new_matrix[i, column]) < 0.001){
        new_matrix[i, column] <- 0
      }
    }
    # find the pivot row for each column if there is a pivot
    if(first_row != nrow(new_matrix)){
      for(i in (first_row+1):nrow(new_matrix)){
        if(new_matrix[i, column] != 0){
          first_row <- i
          break
        }
      }
    }
    # Can only have one pivot per row.
    if((first_row == previous)){
      column <- column + 1
      next
    }
    # make the pivot = 1
    if(new_matrix[first_row, column] != 1){
      new_matrix[first_row,] <- new_matrix[first_row,]/new_matrix[first_row, column]
    }
    # replace numbers below pivot with 0
    if(first_row != nrow(new_matrix)){
      for(i in (first_row+1):nrow(new_matrix)){
        if(new_matrix[i, column] != 0){
          new_matrix <- row_replace(new_matrix, first_row, i)
        }
      }
    }
    # replace numbers above pivot with 0
    if(first_row != 1){
      for(i in 1:(first_row-1)){
        if(new_matrix[i, column] != 0){
          new_matrix <- row_replace(new_matrix, first_row, i)
        }
      }
    }
    # keep track of the last pivot row
    previous <- first_row
    column <- column + 1
  }
  return(new_matrix)
}

rref(rbind(c(-3, 6, -1, 1, -7), c(1, -2, 2, 3, -1), c(2, -4, 5, 8, -4)))
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1   -2    0   -1    3
#> [2,]    0    0    1    2   -2
#> [3,]    0    0    0    0    0
```

This function creates the identity matrix of desired size.

``` r
i_matrix <- function(size){
  matrix <- matrix(0, nrow= size, ncol = size)
  # make diagonals 1
  for(i in 1:nrow(matrix)){
    matrix[i, i] <- 1
  }
  return(matrix)
}

i_matrix(5)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    0    0    0    0
#> [2,]    0    1    0    0    0
#> [3,]    0    0    1    0    0
#> [4,]    0    0    0    1    0
#> [5,]    0    0    0    0    1
```

This function computes the determinant of a matrix.

``` r
det <- function(matrix){
  if(ncol(matrix) != nrow(matrix)){
    stop("Must be an n by n matrix")
  }
  # det = ad-bc for 2x2
  if((ncol(matrix) == 2) & (nrow(matrix) == 2)){
    return(matrix[1,1]*matrix[2,2]-matrix[1,2]*matrix[2,1])
  }
  # need to separate into smaller parts if more than 2x2
  part <- rep(0, ncol(matrix))
  for(i in 1:ncol(matrix)){
    rows <- 2:nrow(matrix)
    columns <- setdiff(1:ncol(matrix), i)
    part[i] <- matrix[1, i]*det(matrix[rows, columns])
    # even number parts are subtracted
    if(i %% 2 == 0){
      part[i] <- -1*part[i]
    }
  }
  return(sum(part))
}

det(matrix(1:49, ncol = 7))
#> [1] 0
```

More functions may come in the future.

## License

This package is licensed under the GPL-3 license.
