#pragma once

using num_t = double;

constexpr int N_PARTICLES{100000}; // Change to adjust particle count.

constexpr int WIDTH = 3840 * 1.5;
constexpr int HEIGHT = 2160 * 1.5;
constexpr num_t DEPTH{3840};

constexpr num_t G{6.67408e-11};
constexpr num_t PI{3.141592653589793238462643383279502884197169399};
constexpr num_t softening{1e-3};

constexpr num_t dt{0.01};
