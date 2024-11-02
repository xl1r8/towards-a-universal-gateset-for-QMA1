#pragma once

#include <vector>

template <typename T> std::vector<std::vector<T>> cartesian_product(const std::vector<std::vector<T>> &input) {
  if (input.empty())
    return {{}};

  std::vector<std::vector<T>> result = {{}};
  for (const auto &pool : input) {
    std::vector<std::vector<T>> new_result;
    for (const auto &x : result) {
      for (const auto &y : pool) {
        auto new_combination = x;
        new_combination.push_back(y);
        new_result.push_back(new_combination);
      }
    }
    result = std::move(new_result);
  }
  return result;
}

template <typename T>
std::vector<std::vector<T>> flat_cartesian_product(const std::vector<std::vector<std::vector<T>>> &input) {
  if (input.empty())
    return {{}};

  std::vector<std::vector<T>> result = {{}};
  for (const auto &pool : input) {
    std::vector<std::vector<T>> new_result;
    for (const auto &x : result) {
      for (const auto &y : pool) {
        auto new_combination = x;
        new_combination.reserve(x.size() + y.size());
        new_combination.insert(new_combination.end(), y.begin(), y.end());
        new_result.push_back(new_combination);
      }
    }
    result = std::move(new_result);
  }
  return result;
}
