#include "sparsemath.h"
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <iostream>

using namespace Spectra;

Index rank(const SparseMatrix<double> &A) {
  SPQR spqr(A);
  if (spqr.info() != Success) {
    std::cerr << "SPQR failed\n";
    return -1;
  }
  return spqr.rank();
}

std::optional<VectorXd> smallest_eigenvalues(const SparseMatrix<double> &A, int k) {
  SparseSymMatProd<double> op(A);
  SymEigsSolver eigs(op, k, 2 * k);
  eigs.init();
  auto nconv = eigs.compute(SortRule::SmallestAlge, 2000, 1e-14, SortRule::SmallestAlge);
  if (eigs.info() == CompInfo::Successful) {
    return eigs.eigenvalues();
  }
  std::cerr << "eigensolver unsuccessful; nconv = " << nconv << '\n';
  return {};
}

bool has_nullspace_dim(const SparseMatrix<double> &A, int d) {
  auto eigs = smallest_eigenvalues(A, std::min(2 * d, (int)A.rows() - 1));
  return eigs && count_zeros(*eigs) == d;
}

int count_zeros(const VectorXd &v, double tol) {
  int count = 0;
  for (Index i = 0; i < v.rows(); ++i) {
    if (std::abs(v[i]) < tol) {
      count++;
    }
  }
  return count;
}

bool vec_in_image(const SparseMatrix<double> &A, const VectorXd &v) {
  SPQR spqr(A);
  return vec_in_image(spqr, A, v);
}

bool vec_in_image(const SPQR<SparseMatrix<double>> &spqr, const SparseMatrix<double> &A, const VectorXd &v) {
  if (spqr.info() != Success) {
    std::cerr << "SPQR error\n";
    return false;
  }
  auto x = spqr.solve(v);
  if (spqr.info() != Success) {
    return false;
  }
  return (A * x - v).isZero(EPS);
}