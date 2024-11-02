#pragma once

// fixes a weird compilation bug where std::ptrdiff_t is not the same as int64_t
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int64_t
#include <Eigen/Sparse>