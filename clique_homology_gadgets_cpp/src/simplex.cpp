#include "simplex.h"
#include <algorithm>
#include <cassert>

Simplex::Simplex(const Simplex &s) : n{s.n}, vertices{s.vertices} {}

Simplex::Simplex(Array &&v, int n, bool sort) : n{n}, vertices{v} {
  if (sort)
    sort_vertices();
}

Simplex::Simplex(const IntVec &v, bool sort) : n{static_cast<int>(v.size())}, vertices{} {
  assert(v.size() <= SIMPLEX_MAX_SIZE);
  std::ranges::copy(v, vertices.begin());
  if (sort)
    sort_vertices();
}

void Simplex::sort_vertices() {
  auto a = vertices.data();
  std::sort(a, a + n);
}

std::pair<Simplex, int> Simplex::make_with_orientation(IntVec v) {
  int sign = 1;
  for (int i = 0; i < v.size(); ++i) {
    for (int j = 0; j < v.size() - i - 1; ++j) {
      if (v[j] > v[j + 1]) {
        sign *= -1;
        std::swap(v[j], v[j + 1]);
      }
    }
  }
  auto simplex = Simplex(v, false);
  return {simplex, sign};
}

Simplex Simplex::remove(int i) const {
  Simplex s(*this);
  for (int j = i + 1; j < s.n; ++j) {
    s.vertices[j - 1] = s.vertices[j];
  }
  s.vertices[s.n - 1] = 0;
  s.n--;
  return s;
}

std::strong_ordering Simplex::operator<=>(const Simplex &o) const {
  if (auto c = n <=> o.n; c != 0)
    return c;
  return vertices <=> o.vertices;
}

bool Simplex::operator==(const Simplex &o) const { return n == o.n && vertices == o.vertices; }