############################
# Linear Algebra Functions #
############################

#' transpose
#'
#' @description transposes a given matrix
#'
#' @param matrix the matrix you would like transposed
#'
#' @examples \dontrun{transpose(matrix(1:6, nrow = 2))}
#'
#' @return the transpose of the matrix you inputed
#'
#' @export
transpose <- function(matrix){
  t_matrix <- matrix(nrow = ncol(matrix), ncol = nrow(matrix))
  # make the rows of the matrix the columns of the transpose
  for(i in 1:nrow(matrix)){
    t_matrix[,i] <- matrix[i,]
  }
  return(t_matrix)
}

#' multiply_matrices
#'
#' @description Allows for the multiplication of matrices
#'
#' @param A The first matrix
#' @param B The second matrix
#'
#' @note Number of columns in matrix A must match the number of rows in matrix B
#'
#' @examples \dontrun{multiply_matrices(matrix(1:6, ncol = 2), matrix(8:13, nrow = 2))}
#'
#' @return the resulting matrix
#'
#' @export
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

#' row_swap
#'
#' @description swaps two rows in a matrix
#'
#' @param matrix The matrix you would like the rows to be swapped in
#' @param row1 The first row you would like swapped
#' @param row2 The second row you would like swapped
#'
#' @examples \dontrun{row_swap(matrix(1:9, nrow = 3), 1, 3)}
#'
#' @return The matrix that results after the rows are swapped
#'
#' @export
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

#' row_replace
#'
#' @description replaces a row by multiplying by another row by a constant and adding it to the row you want to replace.
#'
#' @param matrix The matrix where you would like a row replaced
#' @param row1 The row you use as a reference to replace the other row
#' @param row2 The row to be replaced
#'
#' @examples \dontrun{row_replace(matrix(1:20, nrow = 4), 2, 4)}
#'
#' @return The resulting matrix with the replaced row
#'
#' @export
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

#' arrange_rows
#'
#' @description arranges rows of matrix in descending order starting with the leftmost position of each row.
#'
#' @param matrix The matrix in which you would like to rearrange the rows
#'
#' @examples \dontrun{arrange_rows(rbind(c(0, 0, 1), c(1, 0, 0), c(0, 1, 1)))}
#'
#' @return The matrix with its rows rearranged
#'
#' @export
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

#' rref
#'
#' @description algorithm for getting the row reduced echelon form of a matrix
#'
#' @param matrix The matrix you would like the rref of
#'
#' @examples \dontrun{rref(rbind(c(1, 3, 4, 6, 7), c(0, 0, 4, 5, 3), c(10, 4, 5, 8, 9)))}
#'
#' @return The row reduced echelon form of the given matrix
#'
#' @export
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

#' i_matrix
#'
#' @description creates an identity matrix
#'
#' @param size The number of rows/columns for the identity matrix
#'
#' @examples \dontrun{i_matrix(5)}
#'
#' @return The identity matrix of the given size
#'
#' @export
i_matrix <- function(size){
  matrix <- matrix(0, nrow= size, ncol = size)
  # make diagonals 1
  for(i in 1:nrow(matrix)){
    matrix[i, i] <- 1
  }
  return(matrix)
}

#' is_invertible
#'
#' @description tells whether a given matrix is invertible
#'
#' @param matrix The matrix you want to know if it is invertible
#'
#' @examples \dontrun{is_invertible(matrix(1:25, ncol = 5))}
#'
#' @return TRUE if it is invertible, FALSE if not invertible
#'
#' @export
is_invertible <- function(matrix){
  if(ncol(matrix) != nrow(matrix)){
    return(FALSE)
  }
  else if(det(matrix) == 0){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

#' inverse_matrix
#'
#' @description creates the inverse of a matrix
#'
#' @param matrix The matrix you want the inverse of
#'
#' @examples \dontrun{inverse_matrix(matrix(1:36, ncol = 6))}
#'
#' @return The inverse of the given matrix
#'
#' @export
inverse_matrix <- function(matrix){
  new_matrix <- cbind(matrix, i_matrix(nrow(matrix)))
  new_matrix <- rref(new_matrix)
  first_column <- (ncol(new_matrix)/2)+1
  inverse <- new_matrix[, first_column:ncol(new_matrix)]
  return(inverse)
}

#' det
#'
#' @description computes the determinant of a matrix
#'
#' @param matrix The matrix you want the determinant of
#'
#' @note Number of rows must equal number of columns
#'
#' @examples \dontrun{det(matrix(1:49, ncol = 7))}
#' @examples \dontrun{det(i_matrix(5))}
#'
#' @return The determinant of the given matrix
#'
#' @export
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
