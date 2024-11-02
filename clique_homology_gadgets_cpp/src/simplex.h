#pragma once

#include "types.h"
#include <array>
#include <limits>

constexpr size_t SIMPLEX_MAX_SIZE = 10;

struct Simplex {
  using Array = std::array<int, SIMPLEX_MAX_SIZE>;
  int n;
  Array vertices{};

  Simplex(const Simplex &s);
  Simplex(Array &&v, int n, bool sort = true);
  explicit Simplex(const IntVec &v, bool sort = true);

  void sort_vertices();

  static std::pair<Simplex, int> make_with_orientation(IntVec v);

  [[nodiscard]] Simplex remove(int i) const;

  std::strong_ordering operator<=>(const Simplex &o) const;
  bool operator==(const Simplex &o) const;
};

template <typename T> void hash_combine(std::size_t &seed, const T &v);

template <typename T> void hash_combine(std::size_t &seed, const T &v) {
  static constexpr auto digits = std::numeric_limits<std::size_t>::digits;
  static_assert(digits == 64 || digits == 32);

  if constexpr (digits == 64) {
    std::size_t x = seed + 0x9e3779b9 + std::hash<T>()(v);
    const std::size_t m = 0xe9846af9b1a615d;
    x ^= x >> 32;
    x *= m;
    x ^= x >> 32;
    x *= m;
    x ^= x >> 28;
    seed = x;
  } else {
    std::size_t x = seed + 0x9e3779b9 + std::hash<T>()(v);
    const std::size_t m1 = 0x21f0aaad;
    const std::size_t m2 = 0x735a2d97;
    x ^= x >> 16;
    x *= m1;
    x ^= x >> 15;
    x *= m2;
    x ^= x >> 15;
    seed = x;
  }
}

template <> struct std::hash<Simplex> {
  std::size_t operator()(const Simplex &s) const noexcept {
    size_t h = std::hash<int>()(s.n);
    for (int i = 0; i < s.n; ++i) {
      hash_combine(h, s.vertices[i]);
    }
    return h;
  }
};
