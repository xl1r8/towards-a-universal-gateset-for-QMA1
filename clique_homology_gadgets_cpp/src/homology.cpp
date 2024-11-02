#include "homology.h"
#include <algorithm>
#include <iostream>
#include <utility>

extern "C" {
#include "thirdparty/cliquer-1.21/cliquer.h"
}

CliqueHomology::CliqueHomology(Graph g) : graph{std::move(g)} {}

void CliqueHomology::process() {
  simplex_to_id.emplace_back();
  Simplex empty(IntVec{});
  simplex_to_id[0][empty] = 0;

  auto n = graph.size();
  auto cg = graph_new(n);
  for (int i = 0; i < n; ++i) {
    for (int j : graph.adj[i]) {
      SET_ADD_ELEMENT(cg->edges[i], j);
    }
  }
  auto user_function = [](set_t s, graph_t *g, clique_options *opt) {
    auto *ch = static_cast<CliqueHomology *>(opt->user_data);
    auto &gg = ch->graph;
    IntVec clique;
    for (int i = 0; i < SET_MAX_SIZE(s); i++) {
      if (SET_CONTAINS(s, i)) {
        clique.push_back(i);
      }
    }
    std::sort(clique.begin(), clique.end());
    //        for (const auto &v: clique) {
    //            std::cout << gg.id_to_label[v] << ", ";
    //        }
    //        std::cout << '\n';

    Simplex simplex(clique, false);
    int n = simplex.n;
    while (ch->simplex_to_id.size() <= n + 1)
      ch->simplex_to_id.emplace_back();
    auto id = static_cast<int>(ch->simplex_to_id[n].size());
    ch->simplex_to_id[n][simplex] = id;

    return 1;
  };
  clique_options options = {reorder_by_default, nullptr, nullptr, nullptr, user_function, this, nullptr, 0};
  clique_unweighted_find_all(cg, 1, 1000, false, &options);
  graph_free(cg);
}

const SparseMatrix<double> &CliqueHomology::diff(int k, bool print) {
  if (!diff_cache[k]) {
    auto n1 = static_cast<Index>(simplex_to_id[k + 1].size());
    auto n0 = static_cast<Index>(simplex_to_id[k].size());
    // printf("generate diff(%d): %lld x %lld\n", k, n1, n0);
    SparseMatrix<double> M(n0, n1);
    // M.reserve(n1 * (k + 1));
    std::vector<Triplet<double>> entries;
    entries.reserve(n1 * (k + 1));
    for (const auto &sd : simplex_to_id[k + 1]) {
      auto j = sd.second;
      if (print) {
        print_simplex(std::cout, sd.first);
        std::cout << " -> ";
      }
      for (int l = 0; l < k + 1; ++l) {
        auto s = sd.first.remove(l);
        auto i = simplex_to_id[k][s];
        // M.insert(i, j) = (l % 2) == 0 ? 1 : -1;
        entries.emplace_back(i, j, (l % 2) == 0 ? 1 : -1);
        if (print) {
          if (l > 0) {
            std::cout << ((l % 2 == 0) ? " + " : " - ");
          }
          print_simplex(std::cout, s);
        }
      }
      if (print)
        std::cout << '\n';
    }
    if (print)
      std::cout << '\n';
    M.setFromTriplets(entries.begin(), entries.end());
    diff_cache[k] = M;
    // std::cout << "made diff(" << k << ")\n";
  }
  return diff_cache[k].value();
}

SparseMatrix<double> &CliqueHomology::laplacian(int k) {
  if (!laplacian_cache[k]) {
    auto L = diff(k).transpose() * diff(k) + diff(k + 1) * diff(k + 1).transpose();
    laplacian_cache[k] = L;
  }
  return laplacian_cache[k].value();
}

void CliqueHomology::print_simplex(std::ostream &stream, const Simplex &s) const {
  stream << '[';
  for (int i = 0; i < s.n; ++i) {
    stream << graph.id_to_label[s.vertices[i]];
    if (i < s.n - 1) {
      stream << ", ";
    }
  }
  stream << ']';
}

int CliqueHomology::euler_characteristic() const {
  int sum = 0;
  for (int i = 0, j = 1; i < simplex_to_id.size(); ++i, j *= -1) {
    auto x = j * static_cast<int>(simplex_to_id[i].size());
    // std::cout << "Euler " << x << "\n";
    sum += x;
  }
  return sum;
}
