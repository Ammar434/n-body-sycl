#pragma once
#include "../include/constants.hpp"
#include <array>
#include <cmath>
#include <random>
#include <vector>

class Derivative {
public:
  std::array<num_t, 3> dp; // change in position
  std::array<num_t, 3> dv; // change in velocity
};

class Particle {
public:
  num_t pos[3];
  num_t vel[3];
  num_t acc[3];

  num_t mass;
  num_t energy;

  num_t local_density;

  // virtual void Particle() {};

  void update_position() {

    for (int i = 0; i < 3; ++i) {
      vel[i] += acc[i] * dt;
      pos[i] += vel[i] * dt;
      acc[i] = 0;
    }
  }

  virtual ~Particle() {}
};