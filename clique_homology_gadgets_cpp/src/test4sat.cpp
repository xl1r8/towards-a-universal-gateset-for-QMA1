#include "testutil.h"

// it suffices to check 0000 since the gadgets for other 4-qubit computational basis states are the same, up to renaming
// vertices
TEST(fourSAT, state_0000) {
  QubitGraph g(4);
  auto cycle = Graph::from_edges(g.bit_simplices(false, 0))
                   .make_join(Graph::from_edges(g.bit_simplices(false, 1)))
                   .make_join(Graph::from_edges(g.bit_simplices(false, 2)))
                   .make_join(Graph::from_edges(g.bit_simplices(false, 3)));
  verify_cycle(g, cycle, {}, {{0, 0, 0, 0}});

  g.fill_cycle(cycle);

  VectorXd v = B(0, 0, 0, 0);

  std::vector<VectorXd> basis;
  for (int i = 1; i < 16; ++i) {
    basis.push_back(B(bool(i & 1), bool((i >> 1) & 1), bool((i >> 2) & 1), bool((i >> 3) & 1)));
  }
  verify_gadget(g, v, basis);
}

// |11> - |01>
TEST(fourSAT, state_11m01) {
  QubitGraph g(2);

  auto [g_cycle_minus, f] = g.cycle_0m1();
  auto g_cycle = g_cycle_minus.make_join(Graph::from_edges(g.bit_simplices(true, 1)));

  verify_cycle(g, g_cycle, f, {{1, 1}, {0, 1}});

  g.fill_cycle(g_cycle, f);

  VectorXd v = B(1, 1) - B(0, 1);
  VectorXd w1 = B(1, 1) + B(0, 1);
  VectorXd w2 = B(0, 0);
  VectorXd w3 = B(1, 0);
  verify_gadget(g, v, {w1, w2, w3});
}

// |0000> - |1100>
TEST(fourSAT, state_0000m1100) {
  QubitGraph g(4);

  auto [g_cycle, f] = g.cycle_flip2(0, 0);
  g_cycle = g_cycle.make_join(Graph::from_edges(g.bit_simplices(false, 2)))
                .make_join(Graph::from_edges(g.bit_simplices(false, 3)));
  verify_cycle(g, g_cycle, f, {{0, 0, 0, 0}, {1, 1, 0, 0}});

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

// |1011> - |0111>
TEST(fourSAT, state_1011m0111) {
  QubitGraph g(4);

  auto [g_cycle, f] = g.cycle_flip2(1, 0);
  g_cycle = g_cycle.make_join(Graph::from_edges(g.bit_simplices(true, 2)))
                .make_join(Graph::from_edges(g.bit_simplices(true, 3)));
  verify_cycle(g, g_cycle, f, {{1, 0, 1, 1}, {0, 1, 1, 1}});

  g.fill_cycle(g_cycle, f);

  VectorXd v = B(1, 0, 1, 1) - B(0, 1, 1, 1);

  std::vector<VectorXd> basis;
  for (bool a : {false, true}) {
    for (bool b : {false, true}) {
      if (!a || !b) {
        basis.emplace_back(B(1, 0, a, b) - B(0, 1, a, b));
      }
      basis.emplace_back(B(0, 0, a, b) + B(1, 1, a, b));
      basis.emplace_back(B(1, 0, a, b) + B(0, 1, a, b));
      basis.emplace_back(B(0, 0, a, b) - B(1, 1, a, b));
    }
  }

  verify_gadget(g, v, basis, false);
}

// (|01> - |10> + |11> = -(|10> - |11> - |01>))|0>
TEST(fourSAT, state_010m100p110) {
  QubitGraph g(3);
  auto [g_cycle, f] = g.cycle_Htrans_01m10p11();
  g_cycle = g_cycle.make_join(Graph::from_edges(g.bit_simplices(false, 2)));
  verify_cycle(g, g_cycle, f, {{0, 1, 0}, {1, 0, 0}, {1, 1, 0}});

  g.fill_cycle(g_cycle, f);

  VectorXd v = B(0, 1, 0) - B(1, 0, 0) + B(1, 1, 0);
  std::vector<VectorXd> basis{B(0, 1, 1) - B(1, 0, 1) + B(1, 1, 1)};
  for (bool a : {false, true}) {
    basis.emplace_back(B(0, 1, a) + B(1, 0, a));
    basis.emplace_back(B(1, 0, a) + B(1, 1, a));
    basis.emplace_back(B(0, 0, a));
  }
  verify_gadget(g, v, basis);
}

// (|00> - |10> - |11> = -(|10> + |11> - |00>))|1>
TEST(fourSAT, state_001m101m111) {
  QubitGraph g(3);
  auto [g_cycle, f] = g.cycle_Htrans_00m10m11();
  g_cycle = g_cycle.make_join(Graph::from_edges(g.bit_simplices(true, 2)));
  verify_cycle(g, g_cycle, f, {{0, 0, 1}, {1, 0, 1}, {1, 1, 1}});

  g.fill_cycle(g_cycle, f);

  VectorXd v = B(0, 0, 1) - B(1, 0, 1) - B(1, 1, 1);
  std::vector<VectorXd> basis{B(0, 0, 0) - B(1, 0, 0) - B(1, 1, 0)};
  for (bool a : {false, true}) {
    basis.emplace_back(B(0, 0, a) + B(1, 1, a));
    basis.emplace_back(B(1, 0, a) - B(1, 1, a));
    basis.emplace_back(B(0, 1, a));
  }
  verify_gadget(g, v, basis);
}
