#include "testutil.h"

TEST(CliqueHomologyGadgets, Diff) {
  QubitGraph g(1);
  g.make_homology();

  // for (auto &x: g.g.label_to_id) {
  //     std::cout << x.first << ": " << x.second << '\n';
  // }
  // std::cout << '\n';

  auto v = V(g.x[0], g.a4[0]);
  auto &d = g.h->diff(1);
  auto v2 = d * v;
  VectorXd v2_ = V(g.a4[0]) - V(g.x[0]);
  EXPECT_TRUE((v2 - v2_).isZero(EPS));
}

TEST(CliqueHomologyGadgets, Diff2) {
  QubitGraph g(2);
  g.make_homology();

  // for (auto &x: g.g.label_to_id) {
  //     std::cout << x.first << ": " << x.second << '\n';
  // }
  // std::cout << '\n';

  auto v = V(g.x[0], g.a4[0]);
  auto &d = g.h->diff(1);
  VectorXd v2 = d * v;
  VectorXd v2_ = V(g.a4[0]) - V(g.x[0]);
  EXPECT_TRUE((v2 - v2_).isZero(EPS));

  v = V(g.x[0], g.a4[0], g.x[1], g.a3[1]);
  auto &d3 = g.h->diff(3);
  v2 = d3 * v;
  v2_ = V(g.a4[0], g.x[1], g.a3[1]) - V(g.x[0], g.x[1], g.a3[1]) + V(g.x[0], g.a4[0], g.a3[1]) -
        V(g.x[0], g.a4[0], g.x[1]);
  EXPECT_TRUE((v2 - v2_).isZero(EPS));
}

TEST(CliqueHomologyGadgets, Laplacian) {
  QubitGraph g(1);
  g.make_homology();
  auto L = g.h->laplacian(1);
  EXPECT_EQ(L.rows(), L.cols());
  EXPECT_EQ(L.rows(), g.h->simplex_to_id[2].size());
  EXPECT_EQ(corank(L), 2);
}

TEST(CliqueHomologyGadgets, Laplacian2) {
  int max_n = 3;
#ifdef NDEBUG
  max_n = 5;
#endif
  for (int n = 2; n <= max_n; ++n) {
    QubitGraph g(n);
    g.make_homology();
    // std::cout << "made homology\n";
    auto L = g.h->laplacian(2 * n - 1);
    EXPECT_TRUE(L.isApprox(L.transpose()));
    // std::cout << L.rows() << '\n';
    // std::cout << "made laplacian\n";
    // EXPECT_EQ(corank(L), 1 << n);
    std::cout << "Laplacian2 n = " << n << '\n';
    EXPECT_TRUE(has_nullspace_dim(L, 1 << n));
    if (n <= 3)
      EXPECT_EQ(corank(L), 1 << n);
    if (n <= 4) {
      auto &d = g.h->diff(2 * n - 1);
      EXPECT_EQ(corank(d), 1 << n);
    }
  }
}

TEST(itertools, product) {
  std::vector<std::vector<int>> input = {{1, 2}, {3, 4}, {5, 6}};
  std::vector<std::vector<int>> expected_output = {{1, 3, 5}, {1, 3, 6}, {1, 4, 5}, {1, 4, 6},
                                                   {2, 3, 5}, {2, 3, 6}, {2, 4, 5}, {2, 4, 6}};

  auto result = cartesian_product(input);

  std::sort(result.begin(), result.end());
  std::sort(expected_output.begin(), expected_output.end());

  EXPECT_EQ(result, expected_output);
}

TEST(itertools, flat_product) {
  IntVecVecVec input = {{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}, {{9}}};
  std::vector<std::vector<int>> expected_output = {{1, 2, 5, 6, 9}, {1, 2, 7, 8, 9}, {3, 4, 5, 6, 9}, {3, 4, 7, 8, 9}};

  auto result = flat_cartesian_product(input);

  std::sort(result.begin(), result.end());
  std::sort(expected_output.begin(), expected_output.end());

  EXPECT_EQ(result, expected_output);
}

TEST(CliqueHomologyGadgets, FillCycle_ket0) {
  QubitGraph g(1);
  Graph g_cycle;
  g_cycle.add_edges(g.bit_simplices(false, 0));
  auto g_fill = g.g.fill_cycle(g_cycle);
  g.g = g_fill;
  g.make_homology();
  auto &L = g.h->laplacian(1);
  EXPECT_EQ(corank(L), 1);

  auto ket0 = B(false);
  auto ket1 = B(true);

  verify_gadget(g, ket0, {ket1});
}

TEST(CliqueHomologyGadgets, FillCycle_ketminus) {
  QubitGraph g(1);
  Graph g_cycle;
  std::string x1 = "0.x1";
  // g_cycle.add_edges(
  //         {{g.a2[0], g.a3[0]},
  //          {g.a3[0], g.x[0]},
  //          {g.x[0],  g.b4[0]},
  //          {g.b4[0], g.b2[0]},
  //          {g.b2[0], g.b3[0]},
  //          {g.b3[0], x1},
  //          {x1,      g.a4[0]},
  //          {g.a4[0], g.a2[0]}}
  // );
  g_cycle.add_edges({{g.a2[0], g.a3[0]},
                     {g.a3[0], g.x[0]},
                     {g.x[0], g.b3[0]},
                     {g.b3[0], g.b2[0]},
                     {g.b2[0], g.b4[0]},
                     {g.b4[0], x1},
                     {x1, g.a4[0]},
                     {g.a4[0], g.a2[0]}});
  auto g_fill = g.g.fill_cycle(g_cycle, {{x1, g.x[0]}});
  g.g = g_fill;

  g.make_homology();
  auto &L = g.h->laplacian(1);
  EXPECT_EQ(corank(L), 1);

  auto ket0 = B(false);
  auto ket1 = B(true);

  VectorXd ket_plus = ket0 + ket1, ket_minus = ket0 - ket1;

  verify_gadget(g, ket_minus, {ket_plus});
}

TEST(CliqueHomologyGadgets, FillCycle_ket0_ket1) {
  QubitGraph g(2);

  Graph g_cycle_0, g_cycle_1;
  std::string x1 = "1.x1";
  g_cycle_0.add_edges(g.bit_simplices(false, 0));
  g_cycle_1.add_edges(g.bit_simplices(true, 1));

  auto g_cycle = g_cycle_0.make_join(g_cycle_1);
  auto g_fill = g.g.fill_cycle(g_cycle, {{x1, g.x[1]}});
  g.g = g_fill;
  g.make_homology();
  auto &L = g.h->laplacian(3);
  EXPECT_EQ(corank(L), 3);

  auto v00 = B(0, 0);
  auto v01 = B(0, 1);
  auto v10 = B(1, 0);
  auto v11 = B(1, 1);
  verify_gadget(g, v01, {v00, v10, v11});
}

TEST(CliqueHomologyGadgets, FillCycle_ket0_ketminus) {
  QubitGraph g(2);

  Graph g_cycle_0, g_cycle_minus;
  std::string x1 = "1.x1";
  g_cycle_0.add_edges(g.bit_simplices(false, 0));
  g_cycle_minus.add_edges({{g.a2[1], g.a3[1]},
                           {g.a3[1], g.x[1]},
                           {g.x[1], g.b3[1]},
                           {g.b3[1], g.b2[1]},
                           {g.b2[1], g.b4[1]},
                           {g.b4[1], x1},
                           {x1, g.a4[1]},
                           {g.a4[1], g.a2[1]}});

  auto g_cycle = g_cycle_0.make_join(g_cycle_minus);
  auto g_fill = g.g.fill_cycle(g_cycle, {{x1, g.x[1]}});
  g.g = g_fill;

  // std::cout << '\n';
  // for (const auto &e: g.g.sorted_edges()) {
  //     std::cout << e.first << ' ' << e.second << '\n';
  // }
  // std::cout << '\n';

  g.make_homology();
  auto &L = g.h->laplacian(3);
  EXPECT_EQ(corank(L), 3);

  VectorXd v = B(0, 0) - B(0, 1);

  VectorXd w1 = B(0, 0) + B(0, 1);
  VectorXd w2 = B(1, 0);
  VectorXd w3 = B(1, 1);
  verify_gadget(g, v, {w1, w2, w3});
}

TEST(CliqueHomologyGadgets, Gadget_00) {
  QubitGraph g(2);
  auto cycle = Graph::from_edges(g.bit_simplices(false, 0)).make_join(Graph::from_edges(g.bit_simplices(false, 1)));
  g.fill_cycle(cycle);
  auto &L = g.h->laplacian(3);
  EXPECT_EQ(corank(L), 3);
  EXPECT_EQ(g.h->euler_characteristic(), 3);
}

TEST(CliqueHomologyGadgets, Gadget_00_11) {
  QubitGraph g(2);
  auto [g_cycle, f] = g.cycle_flip2(0, 0);

  g.fill_cycle(g_cycle, f);
  auto &L = g.h->laplacian(3);
  EXPECT_EQ(corank(L), 3);

  VectorXd v = B(0, 0) - B(1, 1);
  VectorXd w1 = B(0, 0) + B(1, 1);
  VectorXd w2 = B(0, 1) + B(1, 0);
  VectorXd w3 = B(0, 1) - B(1, 0);
  verify_gadget(g, v, {w1, w2, w3});
}

TEST(CliqueHomologyGadgets, Gadget_00_11_plus) {
  QubitGraph g(2);
  auto [g_cycle, f] = g.cycle_flip2(0, 0, false);
  g.fill_cycle(g_cycle, f);
  auto &L = g.h->laplacian(3);
  EXPECT_EQ(corank(L), 3);

  VectorXd v = B(0, 0) + B(1, 1);
  VectorXd w1 = B(0, 0) - B(1, 1);
  VectorXd w2 = B(0, 1) + B(1, 0);
  VectorXd w3 = B(0, 1) - B(1, 0);
  verify_gadget(g, v, {w1, w2, w3});
}

TEST(CliqueHomologyGadgets, Gadget_01_10) {
  QubitGraph g(2);
  auto [g_cycle, f] = g.cycle_flip2(0, 1);
  g.fill_cycle(g_cycle, f);
  auto &L = g.h->laplacian(3);
  EXPECT_EQ(corank(L), 3);

  VectorXd v = B(0, 1) - B(1, 0);
  VectorXd w1 = B(0, 0) + B(1, 1);
  VectorXd w2 = B(0, 1) + B(1, 0);
  VectorXd w3 = B(0, 0) - B(1, 1);
  verify_gadget(g, v, {w1, w2, w3});
}

TEST(CliqueHomologyGadgets, Gadget_0000_1100) {
  QubitGraph g(4);
  auto [g_cycle, f] = g.cycle_flip2(0, 0);
  g_cycle = g_cycle.make_join(Graph::from_edges(g.bit_simplices(false, 2)))
                .make_join(Graph::from_edges(g.bit_simplices(false, 3)));
  g.fill_cycle(g_cycle, f);

  VectorXd v = B(0, 0, 0, 0) - B(1, 1, 0, 0);

  std::vector<VectorXd> basis;
  for (bool a : {false, true}) {
    for (bool b : {false, true}) {
      if (a || b) {
        basis.emplace_back(B(0, 0, a, b) - B(1, 1, a, b));
      }
      basis.emplace_back(B(0, 0, a, b) + B(1, 1, a, b));
      basis.emplace_back(B(0, 1, a, b) + B(1, 0, a, b));
      basis.emplace_back(B(0, 1, a, b) - B(1, 0, a, b));
    }
  }

  verify_gadget(g, v, basis, false);
}

TEST(CliqueHomologyGadgets, Gadget_00p10p11) {
  QubitGraph g(2);
  auto [g_cycle, f] = g.cycle_00p10p11();
  g.fill_cycle(g_cycle, f);

  verify_cycle(g, g_cycle, f, {{0, 0}, {1, 0}, {1, 1}});

  VectorXd v = B(0, 0) + B(1, 0) + B(1, 1);
  VectorXd w1 = B(0, 0) - B(1, 1);
  VectorXd w2 = B(1, 0) - B(1, 1);
  VectorXd w3 = B(0, 1);
  verify_gadget(g, v, {w1, w2, w3});
}

TEST(CliqueHomologyGadgets, Gadget_00p10p11_0) {
  QubitGraph g(3);
  auto [g_cycle, f] = g.cycle_00p10p11();

  g_cycle = g_cycle.make_join(Graph::from_edges(g.bit_simplices(false, 2)));
  g.fill_cycle(g_cycle, f);

  VectorXd v = B(0, 0, 0) + B(1, 0, 0) + B(1, 1, 0);
  std::vector<VectorXd> basis;
  for (bool a : {false, true}) {
    if (a)
      basis.emplace_back(B(0, 0, a) + B(1, 0, a) + B(1, 1, a));
    basis.emplace_back(B(0, 0, a) - B(1, 1, a));
    basis.emplace_back(B(1, 0, a) - B(1, 1, a));
    basis.emplace_back(B(0, 1, a));
  }
  EXPECT_EQ(g.h->euler_characteristic(), 7);
  verify_gadget(g, v, basis);
}

TEST(CliqueHomologyGadgets, Gadget_H_trans_00m10m11) {
  QubitGraph g(2);
  auto [g_cycle, f] = g.cycle_Htrans_00m10m11();
  g.fill_cycle(g_cycle, f);
  auto &L = g.h->laplacian(3);
  EXPECT_EQ(corank(L), 3);

  VectorXd v = B(0, 0) - B(1, 0) - B(1, 1);
  VectorXd w1 = B(0, 0) + B(1, 1);
  VectorXd w2 = B(1, 0) - B(1, 1);
  VectorXd w3 = B(0, 1);
  verify_gadget(g, v, {w1, w2, w3});
}

TEST(CliqueHomologyGadgets, Gadget_00m11m11) {
  QubitGraph g(2);
  auto [g_cycle, f] = g.cycle_00m11m11();
  g.fill_cycle(g_cycle, f);
  auto &L = g.h->laplacian(3);
  EXPECT_EQ(corank(L), 3);

  VectorXd v = B(0, 0) - 2 * B(1, 1);
  VectorXd w1 = 2 * B(0, 0) + B(1, 1);
  VectorXd w2 = B(1, 0);
  VectorXd w3 = B(0, 1);
  verify_gadget(g, v, {w1, w2, w3});
}

TEST(CliqueHomologyGadgets, Gadget_H_trans_01m10p11) {
  QubitGraph g(2);
  auto [g_cycle, f] = g.cycle_Htrans_01m10p11();
  g.fill_cycle(g_cycle, f);
  auto &L = g.h->laplacian(3);
  EXPECT_EQ(corank(L), 3);

  VectorXd v = B(0, 1) - B(1, 0) + B(1, 1);
  VectorXd w1 = B(0, 1) + B(1, 0);
  VectorXd w2 = B(1, 0) + B(1, 1);
  VectorXd w3 = B(0, 0);
  verify_gadget(g, v, {w1, w2, w3});
}

TEST(CliqueHomologyGadgets, Gadget_00p11p11p11) {
  QubitGraph g(2);
  auto [g_cycle, f] = g.cycle_00p11p11p11();
  g.fill_cycle(g_cycle, f);
  auto &L = g.h->laplacian(3);

  VectorXd v = B(0, 0) + 3 * B(1, 1);
  VectorXd w1 = 3 * B(0, 0) - B(1, 1);
  VectorXd w2 = B(0, 1);
  VectorXd w3 = B(1, 0);
  verify_gadget(g, v, {w1, w2, w3});
}
