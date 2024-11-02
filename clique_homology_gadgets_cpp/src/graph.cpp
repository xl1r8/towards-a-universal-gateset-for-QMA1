#include "graph.h"

#include <algorithm>
#include <cassert>
#include <format>
#include <iostream>
#include <set>
#include <unordered_set>

extern "C" {
#include "thirdparty/cliquer-1.21/cliquer.h"
}

int Graph::add_vertex(const std::string &v) {
  auto it = label_to_id.find(v);
  if (it != label_to_id.end()) {
    return it->second;
  }
  int i = static_cast<int>(label_to_id.size());
  label_to_id[v] = i;
  id_to_label.push_back(v);
  adj.emplace_back();
  return i;
}

void Graph::add_edge(const std::string &u, const std::string &v) {
  int i = add_vertex(u);
  int j = add_vertex(v);
  add_edge(i, j);
}

void Graph::add_edges(const StrVecVec &ee) {
  for (const auto &e : ee) {
    assert(e.size() == 2);
    add_edge(e[0], e[1]);
  }
}

IntPair make_sorted_pair(int i, int j) { return i <= j ? std::make_pair(i, j) : std::make_pair(j, i); }

void Graph::add_edge(int i, int j) {
  if (!has_edge(i, j)) {
    adj[i].push_back(j);
    adj[j].push_back(i);
    edges.insert(make_sorted_pair(i, j));
  }
}

bool Graph::has_edge(const std::string &u, const std::string &v) const {
  return has_edge(label_to_id.at(u), label_to_id.at(v));
}

bool Graph::has_edge(int i, int j) const { return edges.contains(make_sorted_pair(i, j)); }

Graph Graph::make_union(const Graph &o) const {
  Graph ug{*this};
  for (int i = 0; i < o.adj.size(); ++i) {
    for (int j : o.adj[i]) {
      ug.add_edge(o.id_to_label[i], o.id_to_label[j]);
    }
  }
  return ug;
}

Graph Graph::make_join(const Graph &o) const {
#ifndef NDEBUG
  for (const auto &l : o.id_to_label) {
    assert(!label_to_id.contains(l));
  }
#endif
  Graph jg = make_union(o);
  int n1 = size();
  int n2 = o.size();
  for (int i = 0; i < n1; ++i) {
    for (int j = n1; j < n1 + n2; ++j) {
      jg.add_edge(i, j);
    }
  }
  return jg;
}

Graph Graph::thicken() const {
  Graph g;
  const auto n = size();
  for (int i = 0; i < n; ++i) {
    const auto &u = id_to_label[i];
    const auto u0 = std::format("{}.0", u);
    const auto u1 = std::format("{}.1", u);
    for (auto j : adj[i]) {
      const auto &v = id_to_label[j];
      const auto v0 = std::format("{}.0", v);
      const auto v1 = std::format("{}.1", v);
      if (u < v) {
        g.add_edge(u0, v0);
        g.add_edge(u1, v1);
        g.add_edge(u0, v1);
      }
    }
  }
  return g;
}

int Graph::size() const { return static_cast<int>(adj.size()); }

constexpr auto &get_or_key(const auto &map, const auto &key) {
  const auto it = map.find(key);
  return it == map.cend() ? key : it->second;
}

Graph Graph::fill_cycle(const Graph &cycle, const StrStrMap &f, const std::string &v0) const {
  StrStrMap merge;

  auto g_thicken = cycle.thicken();
  for (const auto &v : cycle.id_to_label) {
    // std::cout << "v = " << v << "\n";
    const auto v_0 = std::format("{}.0", v);
    const auto v_1 = std::format("{}.1", v);
    g_thicken.add_edge(v_0, v_1);
    g_thicken.add_edge(v_1, v0);
    merge[v_0] = get_or_key(f, v);
  }

  auto g_fill = make_union(g_thicken);
  return g_fill.merge_vertices(merge);
}

Graph Graph::merge_vertices(const StrStrMap &map) const {
#ifndef NDEBUG
  // for (const auto &label: label_to_id) {
  //     std::cout << "MERGE " << label.first << "\n";
  // }
  // for (const auto &label: id_to_label) {
  //     std::cout << "MERGE2 " << label << "\n";
  // }
  for (const auto &kv : map) {
    // std::cout << "MAP: " << kv.first << " -> " << kv.second << "\n";
    assert(label_to_id.contains(kv.first));
    assert(label_to_id.contains(kv.second));
  }
#endif

  Graph g_merge;
  for (int i = 0; i < size(); ++i) {
    const auto &u = get_or_key(map, id_to_label[i]);
    for (int j : adj[i]) {
      const auto &v = get_or_key(map, id_to_label[j]);
      if (u != v) {
        g_merge.add_edge(u, v);
      }
    }
  }
  return g_merge;
}

std::vector<std::pair<std::string, std::string>> Graph::sorted_edges() const {
  std::vector<std::pair<std::string, std::string>> ee;
  ee.reserve(edges.size());
  for (const auto &e : edges) {
    auto &b = ee.emplace_back(id_to_label[e.first], id_to_label[e.second]);
    if (b.second < b.first) {
      std::swap(b.first, b.second);
    }
  }
  std::sort(ee.begin(), ee.end());
  return ee;
}

Graph Graph::from_edges(const StrVecVec &edges) {
  Graph g;
  for (const auto &e : edges) {
    assert(e.size() == 2);
    g.add_edge(e[0], e[1]);
  }
  return g;
}

void Graph::delete_edge(const std::string &u, const std::string &v) {
  delete_edge(label_to_id.at(u), label_to_id.at(v));
}

void Graph::delete_edge(int i, int j) {
  assert(has_edge(i, j));
  edges.erase(make_sorted_pair(i, j));
  std::erase(adj[i], j);
  std::erase(adj[j], i);
}

StrVecVec Graph::maximal_cliques(const StrStrMap &f) const {
  auto n = size();
  auto cg = graph_new(n);
  for (int i = 0; i < n; ++i) {
    for (int j : adj[i]) {
      SET_ADD_ELEMENT(cg->edges[i], j);
    }
  }
  struct Context {
    const Graph *g{};
    StrVecVec cliques;
    const StrStrMap *f{};
  } c;
  c.g = this;
  c.f = &f;

  auto user_function = [](set_t s, graph_t *g, clique_options *opt) {
    auto *c = static_cast<Context *>(opt->user_data);
    std::set<std::string> clique;
    bool dup_free = true;
    for (int i = 0; i < SET_MAX_SIZE(s); i++) {
      if (SET_CONTAINS(s, i)) {
        auto [it, ins] = clique.insert(get_or_key(*c->f, c->g->id_to_label[i]));
        dup_free &= ins;
      }
    }
    if (dup_free) {
      c->cliques.emplace_back(clique.begin(), clique.end());
    }
    return 1;
  };
  clique_options options = {reorder_by_default, nullptr, nullptr, nullptr, user_function, &c, nullptr, 0};
  clique_unweighted_find_all(cg, 0, 0, true, &options);
  graph_free(cg);

  std::sort(c.cliques.begin(), c.cliques.end());

  // for (const auto &clique : c.cliques) {
  //   std::cout << "Clique: ";
  //   for (const auto &s : clique) {
  //     std::cout << s << ", ";
  //   }
  //   std::cout << '\n';
  // }

  return c.cliques;
}
