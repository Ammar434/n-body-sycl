#pragma once

#include "particle.hpp"
#include <array>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sycl/sycl.hpp>

class Simulation {
public:
  Simulation();
  virtual ~Simulation();

  void moveSyclRK4(sycl::queue &q);
  void moveSyclEuler(sycl::queue &q);

  void computeKineticEnergy(sycl::queue &q);
  void computePotentialEnergy(sycl::queue &q);

  void outputData(const std::string &filename, int step);

  const Particle *get_particles() const;
  const float *get_density() const;

private:
  std::unique_ptr<Particle[]> particle_array;
  std::unique_ptr<float[]> density;
  float kineticEnergy;
  float potentialEnergy;

  void initializeParticles();
  void initializeGalaxyParticles(std::mt19937 &gen,
                                 std::uniform_real_distribution<> &radius_dist,
                                 std::uniform_real_distribution<> &theta_dist,
                                 std::uniform_real_distribution<> &height_dist,
                                 std::exponential_distribution<> &mass_dist,
                                 num_t disk_radius, num_t disk_height,
                                 num_t total_mass);

  void initializeCentralParticle(num_t total_mass);
};
