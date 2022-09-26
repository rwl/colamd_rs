/// Copyright (c) 1998-2016, Timothy A. Davis
/// Copyright (c) 2022, Richard W. Lincoln
/// All Rights Reserved.
///
/// Approximate minimum degree column ordering algorithms.
///
/// `colamd`: an approximate minimum degree column ordering algorithm,
/// for LU factorization of symmetric or unsymmetric matrices,
/// QR factorization, least squares, interior point methods for
/// linear programming problems, and other related problems.
///
/// `symamd`: an approximate minimum degree ordering algorithm for Cholesky
/// factorization of symmetric matrices.
///
/// Colamd computes a permutation `Q` such that the Cholesky factorization of
/// `(AQ)'(AQ)` has less fill-in and requires fewer floating point operations
/// than `A'A`. This also provides a good ordering for sparse partial
/// pivoting methods, `P(AQ) = LU`, where `Q` is computed prior to numerical
/// factorization, and `P` is computed during numerical factorization via
/// conventional partial pivoting with row interchanges. Colamd is the
/// column ordering method used in SuperLU, part of the ScaLAPACK library.
/// It is also available as built-in function in MATLAB Version 6,
/// available from MathWorks, Inc. (http://www.mathworks.com). This
/// routine can be used in place of colmmd in MATLAB.
///
/// Symamd computes a permutation `P` of a symmetric matrix A such that the
/// Cholesky factorization of `PAP'` has less fill-in and requires fewer
/// floating point operations than `A`. Symamd constructs a matrix `M` such
/// that `M'M` has the same nonzero pattern of `A`, and then orders the columns
/// of `M` using colmmd. The column ordering of `M` is then returned as the
/// row and column ordering `P` of `A`.
///
/// The authors of the code itself are Stefan I. Larimore and Timothy A.
/// Davis (davis at cise.ufl.edu), University of Florida. The algorithm was
/// developed in collaboration with John Gilbert, Xerox PARC, and Esmond
/// Ng, Oak Ridge National Laboratory.
///
/// This work was supported by the National Science Foundation, under
/// grants DMS-9504974 and DMS-9803599.
mod codes;
mod col;
mod colamd;
mod colamd2;
mod debug;
mod internal;
mod report;
mod row;
mod stats;
mod symamd;

pub use crate::colamd::colamd;
