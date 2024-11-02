#pragma once

#include "types.h"
#include <unordered_set>

template <> struct std::hash<IntPair> {
  size_t operator()(const IntPair &p) const noexcept { return ((size_t)p.first << 32) | (unsigned int)p.second; }
};

struct Graph {
  IntVecVec adj;
  StrVec id_to_label;
  StrIntMap label_to_id;
  std::unordered_set<IntPair> edges;

  int add_vertex(const std::string &v);

  void add_edge(const std::string &u, const std::string &v);

  void add_edges(const StrVecVec &ee);

  void add_edge(int i, int j);

  static Graph from_edges(const StrVecVec &edges);

  [[nodiscard]] bool has_edge(const std::string &u, const std::string &v) const;

  [[nodiscard]] bool has_edge(int i, int j) const;

  void delete_edge(const std::string &u, const std::string &v);

  void delete_edge(int i, int j);

  [[nodiscard]] Graph make_union(const Graph &o) const;

  [[nodiscard]] Graph make_join(const Graph &o) const;

  [[nodiscard]] Graph thicken() const;

  [[nodiscard]] Graph fill_cycle(const Graph &cycle, const StrStrMap &f = {}, const std::string &v0 = "v0") const;

  [[nodiscard]] Graph merge_vertices(const StrStrMap &map) const;

  [[nodiscard]] int size() const;

  [[nodiscard]] std::vector<std::pair<std::string, std::string>> sorted_edges() const;

  [[nodiscard]] StrVecVec maximal_cliques(const StrStrMap &f = {}) const;
};