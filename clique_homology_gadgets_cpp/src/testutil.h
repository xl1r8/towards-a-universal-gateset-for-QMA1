#pragma once
#include "gadgets.h"
#include "itertools.h"
#include "sparsemath.h"
#include "gtest/gtest.h"

#define V(s...) g.simplices_to_vec({{s}})
#define B(b...) g.embed_bits({b})

inline void verify_cycle(const QubitGraph &g, const Graph &g_cycle, const StrStrMap &f, const BoolVecVec &states) {
  auto mc = g_cycle.maximal_cliques(f);
  auto mc2 = g.cliques_for_states(states);
  EXPECT_EQ(mc, mc2);
}

inline void verify_gadget(QubitGraph &g, const VectorXd &v, const std::vector<VectorXd> &basis,
                          bool compute_d1_rank = true) {
  int dim = (1 << g.n) - 1;
  int k = 2 * g.n - 1;

  EXPECT_EQ(g.h->euler_characteristic(), dim);

  auto &d1 = g.h->diff(k);
  auto &d2 = g.h->diff(k + 1);

  auto &L = g.h->laplacian(k);
  EXPECT_TRUE(has_nullspace_dim(L, dim));
  std::cout << "checked Nullspace of L" << std::endl;

  SPQR spqr_d2(d2);
  auto r2 = spqr_d2.rank();
  std::cout << "rank(d2) = " << r2 << std::endl;

  if (compute_d1_rank) {
    auto cr1 = corank(d1);
    EXPECT_EQ(cr1 - r2, dim);
    std::cout << "corank(d1) = " << cr1 << std::endl;
  }

  EXPECT_TRUE(vec_in_image(spqr_d2, d2, v));
  EXPECT_TRUE((d1 * v).isZero(EPS));
  EXPECT_EQ(basis.size(), dim);

  for (const auto &w : basis) {
    EXPECT_TRUE((d1 * w).isZero(EPS));
  }

  SparseMatrix<double> M(d2);
  M.conservativeResize(d2.rows(), d2.cols() + basis.size());
  for (int j = 0; j < basis.size(); ++j) {
    for (int i = 0; i < M.rows(); ++i) {
      auto x = basis[j][i];
      if (x != 0.) {
        M.insert(i, j + d2.cols()) = x;
      }
    }
  }

  // verify that the basis vectors are independent modulo Im(d2) by checking $Span(basis) \cap Im(d2) = 0$
  EXPECT_EQ(rank(M), r2 + dim);
  std::cout << "checked rank(M)" << std::endl;
}
