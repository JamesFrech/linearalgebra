% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Linear Algebra Functions.R
\name{multiply_matrices}
\alias{multiply_matrices}
\title{multiply_matrices}
\usage{
multiply_matrices(A, B)
}
\arguments{
\item{A}{The first matrix}

\item{B}{The second matrix}
}
\value{
the resulting matrix
}
\description{
Allows for the multiplication of matrices
}
\note{
Number of columns in matrix A must match the number of rows in matrix B
}
\examples{
\dontrun{multiply_matrices(matrix(1:6, ncol = 2), matrix(8:13, nrow = 2))}

}
