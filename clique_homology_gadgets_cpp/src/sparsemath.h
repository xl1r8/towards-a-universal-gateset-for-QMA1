#pragma once

#include "eigen.h"
#include <Eigen/SPQRSupport>
#include <optional>

using namespace Eigen;
constexpr const double EPS = 1e-6;

Index rank(const SparseMatrix<double> &A);

inline Index corank(const SparseMatrix<double> &A) { return A.cols() - rank(A); };

std::optional<VectorXd> smallest_eigenvalues(const SparseMatrix<double> &A, int k);

bool has_nullspace_dim(const SparseMatrix<double> &A, int k);

bool vec_in_image(const SparseMatrix<double> &A, const VectorXd &v);

bool vec_in_image(const SPQR<SparseMatrix<double>> &spqr, const SparseMatrix<double> &A, const VectorXd &v);

int count_zeros(const VectorXd &v, double tol = 1e-9);