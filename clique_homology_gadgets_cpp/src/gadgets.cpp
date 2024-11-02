#include "gadgets.h"
#include "itertools.h"

#include <format>
#include <iostream>

VertexNames::VertexNames(int n, const std::string &suffix) : n{n} {
  for (int i = 0; i < n; ++i) {
    x.push_back(std::format("{}.xx", i)); // don't apply suffix so that the copies share the x vertices
    a2.push_back(std::format("{}.a2{}", i, suffix));
    a3.push_back(std::format("{}.a3{}", i, suffix));
    a4.push_back(std::format("{}.a4{}", i, suffix));
    b2.push_back(std::format("{}.b2{}", i, suffix));
    b3.push_back(std::format("{}.b3{}", i, suffix));
    b4.push_back(std::format("{}.b4{}", i, suffix));
    qubit_vertices.push_back({x[i], a2[i], a3[i], a4[i], b2[i], b3[i], b4[i]});
  }
}

QubitGraph::QubitGraph(int n) : VertexNames(n) {
  for (int i = 0; i < n; ++i) {
    g.add_edge(x[i], a3[i]);
    g.add_edge(a3[i], a2[i]);
    g.add_edge(a2[i], a4[i]);
    g.add_edge(a4[i], x[i]);
    g.add_edge(x[i], b3[i]);
    g.add_edge(b3[i], b2[i]);
    g.add_edge(b2[i], b4[i]);
    g.add_edge(b4[i], x[i]);

    for (const auto &u : qubit_vertices[i]) {
      for (int j = 0; j < i; ++j) {
        for (const auto &v : qubit_vertices[j]) {
          g.add_edge(u, v);
        }
      }
    }
  }
}

void QubitGraph::make_homology() {
  h = std::make_unique<CliqueHomology>(g);
  h->process();
}

int QubitGraph::gadget_dim() const {
  int d = static_cast<int>(h->simplex_to_id[2 * n].size());
  return d;
};

VectorXd QubitGraph::simplices_to_vec(const StrVecVec &simplices) const {
  assert(!simplices.empty());
  int k = static_cast<int>(simplices[0].size());
  int d = static_cast<int>(h->simplex_to_id[k].size());

  VectorXd v(d);
  v.setZero();
  // std::cout << "simplices_to_vec\n" << v << "\n\n";
  for (const auto &labels : simplices) {
    IntVec indices;
    indices.reserve(k);
    assert(k == labels.size());
    for (const auto &l : labels) {
      indices.push_back(g.label_to_id.at(l));
    }
    auto [simplex, sign] = Simplex::make_with_orientation(std::move(indices));
    auto i = h->simplex_to_id[k][simplex];
    v[i] = sign;
  }

  return v;
}

VectorXd QubitGraph::embed_bits(const BoolVec &bits) const {
  assert(bits.size() == n);

  StrVecVecVec local_gadgets;
  for (int i = 0; i < n; ++i) {
    local_gadgets.push_back(bit_simplices(bits[i], i));
  }

  auto prod = flat_cartesian_product(local_gadgets);

  // for (const auto &s: prod) {
  //     for (const auto &l: s) {
  //         std::cout << l << ' ';
  //     }
  //     std::cout << '\n';
  // }
  // std::cout << '\n';

  return simplices_to_vec(prod);
}

StrVecVec QubitGraph::cliques_for_states(const BoolVecVec &states) const {
  StrVecVec cliques;

  for (const auto &bits : states) {
    assert(bits.size() == n);
    StrVecVecVec local_gadgets;
    for (int i = 0; i < n; ++i) {
      local_gadgets.push_back(bit_simplices(bits[i], i));
    }
    auto prod = flat_cartesian_product(local_gadgets);
    for (auto &c : prod) {
      std::sort(c.begin(), c.end());
      cliques.push_back(c);
    }
  }

  std::sort(cliques.begin(), cliques.end());

  return cliques;
}

VectorXd QubitGraph::embed_qubits(const std::vector<std::array<double, 2>> &bits) const {
  BoolVecVec b;
  for (int i = 0; i < n; ++i) {
    b.push_back({false, true});
  }
  const auto bit_strings = cartesian_product(b);
  VectorXd v(gadget_dim());
  v.setZero();
  for (const auto &s : bit_strings) {
    double a = 1;
    for (int i = 0; i < n; ++i) {
      a *= bits[i][s[i]];
    }
    if (a != 0) {
      v += a * embed_bits(s);
    }
  }
  return v;
}

StrVecVec VertexNames::bit_simplices(bool bit, int i) const {
  if (bit)
    return {{x[i], b3[i]}, {b3[i], b2[i]}, {b2[i], b4[i]}, {b4[i], x[i]}};
  else
    return {{x[i], a3[i]}, {a3[i], a2[i]}, {a2[i], a4[i]}, {a4[i], x[i]}};
}

void QubitGraph::fill_cycle(const Graph &cycle, const StrStrMap &f, const std::string &v0) {
  auto g_fill = g.fill_cycle(cycle, f);
  g = g_fill;
  make_homology();
}

[[nodiscard]] std::pair<Graph, StrStrMap> QubitGraph::cycle_flip2(bool b0, bool b1, bool minus) const {
  StrVec dx; // dummy x
  for (int i = 1; i < 5; ++i) {
    dx.push_back(std::format("x{}", i));
  }

  auto dx1 = dx;
  if (!minus) // [Footnote 14, KK23]
    std::swap(dx1[1], dx1[3]);

  Graph g_cycle;
  // Figure 12 in [KK23]
  for (const auto &b : dx) {
    g_cycle.add_edge(x[0], b); // connect top vertex
    g_cycle.add_edge(x[1], b); // connect bottom vertex
  }
  for (int i = 0; i < dx.size(); ++i) {
    g_cycle.add_edge(dx[i], dx[(i + 1) % dx.size()]); // connect cycle x1 x2 x3 x4 x1
  }

  // Loops from Figure 11 for |00> and |11>. Swap for the other cases
  StrVec loop00{a3[0], a4[1], a4[0], a3[1]};
  StrVec loop11{b3[0], b4[1], b4[0], b3[1]};
  if (b0) {
    std::swap(loop00[0], loop11[0]);
    std::swap(loop00[2], loop11[2]);
  }
  if (b1) {
    std::swap(loop00[1], loop11[1]);
    std::swap(loop00[3], loop11[3]);
  }

  // Figure 13
  auto g_cycle_00 = make_part2(bit_simplices(b0, 0), bit_simplices(b1, 1), dx, loop00);
  auto g_cycle_11 = make_part2(bit_simplices(!b0, 0), bit_simplices(!b1, 1), dx1, loop11);

  g_cycle = g_cycle.make_union(g_cycle_00).make_union(g_cycle_11);
  g_cycle.delete_edge(x[0], x[1]);

  StrStrMap f; /* [Eq. 14, KK23] */
  for (const auto &xi : dx) {
    f[xi] = x[0];
  }

  return {g_cycle, f};
}

// build one of the parts in (e.g.) Figure 16 of [KK23]
Graph QubitGraph::make_part2(const StrVecVec &simpl0, const StrVecVec &simpl1, const StrVec &dx,
                             const StrVec &loop) const {
  assert(dx.size() == 4);
  assert(loop.size() == 4);
  auto gc = Graph::from_edges(simpl0).make_join(Graph::from_edges(simpl1));

  for (int i = 0; i < 4; ++i) {
    gc.add_edge(x[0], dx[i]);                    // connect top vertex
    gc.add_edge(x[1], dx[i]);                    // connect bottom vertex
    gc.add_edge(dx[i], dx[(i + 1) % dx.size()]); // connect cycle
    for (int j = 0; j < 2; ++j) {
      gc.add_edge(loop[(i - j + 4) % 4], dx[i]);
    }
  }

  return gc;
};

#define BOILERPLATE_3_2qubits                                                                                          \
  VertexNames n_(2, "_");                                                                                              \
  StrVec dx; /* dummy x */                                                                                             \
  for (int i = 1; i < 7; ++i) {                                                                                        \
    dx.push_back(std::format("x{}", i));                                                                               \
  }                                                                                                                    \
                                                                                                                       \
  StrStrMap f; /* [Eq. 15, KK23] */                                                                                    \
  for (const auto &xi : dx) {                                                                                          \
    f[xi] = x[0];                                                                                                      \
  }                                                                                                                    \
  for (int i = 0; i < 2; ++i) {                                                                                        \
    f[n_.a2[i]] = a2[i];                                                                                               \
    f[n_.a3[i]] = a3[i];                                                                                               \
    f[n_.a4[i]] = a4[i];                                                                                               \
    f[n_.b2[i]] = b2[i];                                                                                               \
    f[n_.b3[i]] = b3[i];                                                                                               \
    f[n_.b4[i]] = b4[i];                                                                                               \
  }

#define BOILERPLATE_4_2qubits                                                                                          \
  VertexNames n_[] = {{2, "_1"}, {2, "_2"}};                                                                           \
  StrVec dx; /* dummy x */                                                                                             \
  for (int i = 1; i < 9; ++i) {                                                                                        \
    dx.push_back(std::format("x{}", i));                                                                               \
  }                                                                                                                    \
                                                                                                                       \
  StrStrMap f; /* [Eq. 16, KK23] */                                                                                    \
  for (const auto &xi : dx) {                                                                                          \
    f[xi] = x[0];                                                                                                      \
  }                                                                                                                    \
  for (int i = 0; i < 2; ++i) {                                                                                        \
    for (int j = 0; j < 2; ++j) {                                                                                      \
      f[n_[j].a2[i]] = a2[i];                                                                                          \
      f[n_[j].a3[i]] = a3[i];                                                                                          \
      f[n_[j].a4[i]] = a4[i];                                                                                          \
      f[n_[j].b2[i]] = b2[i];                                                                                          \
      f[n_[j].b3[i]] = b3[i];                                                                                          \
      f[n_[j].b4[i]] = b4[i];                                                                                          \
    }                                                                                                                  \
  }                                                                                                                    \
  /* [Figure 19, KK23]                                                                                                 \
   * The extra vertices x9 and x10 add spurious holes, but we can fix the gadget by just connecting x2-x7 and x4-x6    \
   */                                                                                                                  \
  Graph g_cycle;                                                                                                       \
  g_cycle.add_edge(dx[1], dx[6]);                                                                                      \
  g_cycle.add_edge(dx[3], dx[5]);

// |00> + |10> + |11>
[[nodiscard]] std::pair<Graph, StrStrMap> QubitGraph::cycle_00p10p11() const {
  BOILERPLATE_3_2qubits;

  // 00
  auto g_cycle = make_part2(bit_simplices(false, 0), bit_simplices(false, 1), {dx[0], dx[1], dx[2], dx[3]},
                            {a3[0], a4[1], a4[0], a3[1]});
  // 11
  g_cycle = g_cycle.make_union(make_part2(bit_simplices(true, 0), bit_simplices(true, 1), {dx[1], dx[4], dx[5], dx[2]},
                                          {b3[0], b4[1], b4[0], b3[1]}));

  // 10
  g_cycle = g_cycle.make_union(make_part2(n_.bit_simplices(true, 0), n_.bit_simplices(false, 1),
                                          {dx[4], dx[0], dx[3], dx[5]}, {n_.b3[0], n_.a4[1], n_.b4[0], n_.a3[1]}));

  g_cycle.delete_edge(x[0], x[1]);
  return {g_cycle, f};
}

// |00> - |10> - |11>
[[nodiscard]] std::pair<Graph, StrStrMap> QubitGraph::cycle_Htrans_00m10m11() const {
  BOILERPLATE_3_2qubits;

  // 00
  auto g_cycle = make_part2(bit_simplices(false, 0), bit_simplices(false, 1), {dx[0], dx[1], dx[2], dx[3]},
                            {a3[0], a4[1], a4[0], a3[1]});
  // 11
  g_cycle = g_cycle.make_union(make_part2(bit_simplices(true, 0), bit_simplices(true, 1), {dx[1], dx[2], dx[5], dx[4]},
                                          {b3[0], b4[1], b4[0], b3[1]}));

  // 10
  g_cycle = g_cycle.make_union(make_part2(n_.bit_simplices(true, 0), n_.bit_simplices(false, 1),
                                          {dx[4], dx[5], dx[3], dx[0]}, {n_.b3[0], n_.a4[1], n_.b4[0], n_.a3[1]}));

  g_cycle.delete_edge(x[0], x[1]);

  return {g_cycle, f};
}

// |00> - |11> - |11>
[[nodiscard]] std::pair<Graph, StrStrMap> QubitGraph::cycle_00m11m11() const {
  BOILERPLATE_3_2qubits;

  // 00
  auto g_cycle = make_part2(bit_simplices(false, 0), bit_simplices(false, 1), {dx[0], dx[1], dx[2], dx[3]},
                            {a3[0], a4[1], a4[0], a3[1]});
  // 11
  g_cycle = g_cycle.make_union(make_part2(bit_simplices(true, 0), bit_simplices(true, 1), {dx[1], dx[2], dx[5], dx[4]},
                                          {b3[0], b4[1], b4[0], b3[1]}));

  // 11
  g_cycle = g_cycle.make_union(make_part2(n_.bit_simplices(true, 0), n_.bit_simplices(true, 1),
                                          {dx[4], dx[5], dx[3], dx[0]}, {n_.b3[0], n_.b4[1], n_.b4[0], n_.b3[1]}));

  g_cycle.delete_edge(x[0], x[1]);

  return {g_cycle, f};
}

// |01> - |10> + |11>
[[nodiscard]] std::pair<Graph, StrStrMap> QubitGraph::cycle_Htrans_01m10p11() const {
  BOILERPLATE_3_2qubits;

  // 01
  auto g_cycle = make_part2(bit_simplices(false, 0), bit_simplices(true, 1), {dx[0], dx[1], dx[2], dx[3]},
                            {a3[0], b4[1], a4[0], b3[1]});
  // 10
  g_cycle = g_cycle.make_union(make_part2(bit_simplices(true, 0), bit_simplices(false, 1), {dx[1], dx[2], dx[5], dx[4]},
                                          {b3[0], a4[1], b4[0], a3[1]}));

  // 11
  g_cycle = g_cycle.make_union(make_part2(n_.bit_simplices(true, 0), n_.bit_simplices(true, 1),
                                          {dx[4], dx[0], dx[3], dx[5]}, {n_.b3[0], n_.b4[1], n_.b4[0], n_.b3[1]}));

  g_cycle.delete_edge(x[0], x[1]);

  return {g_cycle, f};
}

[[nodiscard]] std::pair<Graph, StrStrMap> QubitGraph::cycle_00p11p11p11() const {
  BOILERPLATE_4_2qubits;

  // for (const auto &kv : f) {
  //   std::cout << "f: " << kv.first << " -> " << kv.second << "\n";
  // }

  // [Figure 18, KK23]:
  // (a) 00
  g_cycle = g_cycle.make_union(make_part2(bit_simplices(false, 0), bit_simplices(false, 1),
                                          {dx[0], dx[1], dx[2], dx[3]}, {a3[0], a4[1], a4[0], a3[1]}));
  // (b) 11
  g_cycle = g_cycle.make_union(make_part2(bit_simplices(true, 0), bit_simplices(true, 1), {dx[1], dx[4], dx[5], dx[2]},
                                          {b3[0], b4[1], b4[0], b3[1]}));
  // (c) 11
  g_cycle = g_cycle.make_union(make_part2(n_[0].bit_simplices(true, 0), n_[0].bit_simplices(true, 1),
                                          {dx[4], dx[6], dx[7], dx[5]},
                                          {n_[0].b3[0], n_[0].b4[1], n_[0].b4[0], n_[0].b3[1]}));
  // (d) 11
  g_cycle = g_cycle.make_union(make_part2(n_[1].bit_simplices(true, 0), n_[1].bit_simplices(true, 1),
                                          {dx[6], dx[0], dx[3], dx[7]},
                                          {n_[1].b3[0], n_[1].b4[1], n_[1].b4[0], n_[1].b3[1]}));
  g_cycle.delete_edge(x[0], x[1]);

  // for (const auto &e : g_cycle.sorted_edges()) {
  //   std::cout << "e = " << e.first << ", " << e.second << "\n";
  // }

  return {g_cycle, f};
}

[[nodiscard]] std::pair<Graph, StrStrMap> QubitGraph::cycle_0m1() const {
  Graph g_cycle;
  std::string x1 = "0.x1";
  g_cycle.add_edges({{a2[0], a3[0]},
                     {a3[0], x[0]},
                     {x[0], b3[0]},
                     {b3[0], b2[0]},
                     {b2[0], b4[0]},
                     {b4[0], x1},
                     {x1, a4[0]},
                     {a4[0], a2[0]}});
  return {g_cycle, {{x1, x[0]}}};
}
