#include "../include/simulation.hpp"
#include <algorithm>
#include <cmath>
#include <execution>
#include <random>

Simulation::Simulation() { initializeParticles(); }

Simulation::~Simulation() {
  particle_array.release();
  density.release();
}

void Simulation::initializeParticles() {
  particle_array = std::make_unique<Particle[]>(N_PARTICLES);
  density = std::make_unique<float[]>(N_PARTICLES);

  const num_t total_mass = 1e22;
  const num_t disk_radius = WIDTH;
  const num_t disk_height = disk_radius / 10.0; // Flatter disk

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> radius_dist(0, 1);
  std::uniform_real_distribution<> theta_dist(0, 2 * PI);
  std::uniform_real_distribution<> height_dist(-1, 1);
  std::exponential_distribution<> mass_dist(1.0 / (total_mass / N_PARTICLES));

  initializeGalaxyParticles(gen, radius_dist, theta_dist, height_dist,
                            mass_dist, disk_radius, disk_height, total_mass);
  initializeCentralParticle(total_mass);
}

void Simulation::initializeGalaxyParticles(
    std::mt19937 &gen, std::uniform_real_distribution<> &radius_dist,
    std::uniform_real_distribution<> &theta_dist,
    std::uniform_real_distribution<> &height_dist,
    std::exponential_distribution<> &mass_dist, num_t disk_radius,
    num_t disk_height, num_t total_mass) {
  for (int i = 1; i < N_PARTICLES; ++i) {
    num_t r = disk_radius * std::sqrt(radius_dist(gen)); // Disk distribution
    num_t theta = theta_dist(gen);
    num_t height = disk_height * height_dist(gen);

    particle_array[i].pos[0] = WIDTH / 2 + r * std::cos(theta);
    particle_array[i].pos[1] = HEIGHT / 2 + r * std::sin(theta);
    particle_array[i].pos[2] = height;

    // Calculate circular orbit velocity
    num_t distance_from_center = std::sqrt(r * r + height * height);
    num_t v_circular = std::sqrt(G * total_mass / distance_from_center);

    // Add some random perturbation to the velocity
    std::normal_distribution<> velocity_perturbation(0, v_circular * 0.1);

    particle_array[i].vel[0] =
        v_circular * std::sin(theta) + velocity_perturbation(gen);
    particle_array[i].vel[1] =
        -v_circular * std::cos(theta) + velocity_perturbation(gen);
    particle_array[i].vel[2] = velocity_perturbation(gen);

    for (int j = 0; j < 3; ++j) {
      particle_array[i].acc[j] = 0;
      // particle_array[i].vel[j] = 0;
    }

    particle_array[i].mass = mass_dist(gen);
  }
}

void Simulation::initializeCentralParticle(num_t total_mass) {
  particle_array[0].mass = total_mass * 0.4;

  particle_array[0].pos[0] = WIDTH / 2;
  particle_array[0].pos[1] = HEIGHT / 2;
  particle_array[0].pos[2] = 0;

  for (int j = 0; j < 3; ++j) {
    particle_array[0].vel[j] = 0;
    particle_array[0].acc[j] = 0;
  }
}

SYCL_EXTERNAL Derivative
evaluate(const sycl::accessor<Particle, 1, sycl::access::mode::read_write>
             &particles,
         const Particle &p, num_t dt, const Derivative &d) {
  Particle state;
  for (int i = 0; i < 3; ++i) {
    state.pos[i] = p.pos[i] + d.dp[i] * dt;
    state.vel[i] = p.vel[i] + d.dv[i] * dt;
  }
  state.mass = p.mass;

  Derivative output;
  for (int i = 0; i < 3; ++i) {
    output.dp[i] = state.vel[i];
  }

  std::array<num_t, 3> acc = {0, 0, 0};

  for (int j = 0; j < N_PARTICLES; ++j) {
    const Particle &other = particles[j];
    if (&p != &other) {
      std::array<num_t, 3> diff;
      num_t distanceSqr = softening * softening;

      for (int k = 0; k < 3; ++k) {
        diff[k] = other.pos[k] - state.pos[k];
        distanceSqr += diff[k] * diff[k];
      }

      num_t distanceInv = sycl::rsqrt(distanceSqr);
      num_t c = sycl::min(distanceSqr / 100000000, num_t(1));
      num_t forceGravity = other.mass * (1.0 - sycl::sin(c * (M_PI / 2)));
      num_t factor = G * forceGravity * distanceInv * distanceInv * distanceInv;

      for (int k = 0; k < 3; k++) {
        acc[k] += diff[k] * factor;
      }
    }
  }

  for (int i = 0; i < 3; ++i) {
    output.dv[i] = acc[i];
  }

  return output;
}

void Simulation::moveSyclRK4(sycl::queue &q) {
  sycl::buffer<Particle, 1> particleBuf(particle_array.get(),
                                        sycl::range<1>(N_PARTICLES));
  sycl::buffer<float, 1> densityBuf(density.get(), sycl::range<1>(N_PARTICLES));

  try {
    q.submit([&](sycl::handler &h) {
      sycl::accessor particles(particleBuf, h, sycl::read_write);
      sycl::accessor density_acc(densityBuf, h, sycl::write_only);

      h.parallel_for(sycl::range<1>(N_PARTICLES), [=](sycl::id<1> idx) {
        int i = idx[0];
        if (i > 0) {

          Particle &p = particles[i];
          /*
            // Calculate accelerations
            std::array<num_t, 3> acc = {0, 0, 0};
            for (int j = 0; j < N_PARTICLES; ++j) {
              if (i != j) {
                const Particle &other = particles[j];

                // Euler method
                std::array<num_t, 3> diff;
                num_t distanceSqr = softening * softening;

                for (int k = 0; k < 3; ++k) {
                  diff[k] = other.pos[k] - p.pos[k];
                  distanceSqr += diff[k] * diff[k];
                }

                num_t distanceInv = sycl::rsqrt(distanceSqr);
                num_t c = sycl::min(distanceSqr / 100000, num_t(1));
                num_t forceGravity = other.mass * (1.0 - sycl::sin(c * (PI /
            2))); num_t factor = G * forceGravity * distanceInv * distanceInv *
            distanceInv;

                for (int k = 0; k < 3; k++) {
                  acc[k] += diff[k] * factor;
                }
              }
            }

            // Update position
            if (i > 0) {
              for (int k = 0; k < 3; k++) {
                p.vel[k] += acc[k] * dt;
                p.pos[k] += p.vel[k] * dt;
              }
            }
            */

          // RK4 integration
          Derivative k1 = evaluate(particles, p, 0.0f, Derivative());
          Derivative k2 = evaluate(particles, p, dt * 0.5f, k1);
          Derivative k3 = evaluate(particles, p, dt * 0.5f, k2);
          Derivative k4 = evaluate(particles, p, dt, k3);

          // Update position and velocity
          for (int j = 0; j < 3; ++j) {
            p.pos[j] +=
                (k1.dp[j] + 2.0f * k2.dp[j] + 2.0f * k3.dp[j] + k4.dp[j]) * dt /
                6.0f;
            p.vel[j] +=
                (k1.dv[j] + 2.0f * k2.dv[j] + 2.0f * k3.dv[j] + k4.dv[j]) * dt /
                6.0f;
          }

          // Calculate density
          const float radius_square = 10000.0f; // 100.0f * 100.0f
          int count = 0;

          for (int j = 0; j < N_PARTICLES; ++j) {
            if (i == j)
              continue;

            const Particle &other = particles[j];
            float dx = p.pos[0] - other.pos[0];
            float dy = p.pos[1] - other.pos[1];
            float dz = p.pos[2] - other.pos[2];
            float distance_sq = dx * dx + dy * dy + dz * dz;

            if (distance_sq < radius_square) {
              count++;
            }
          }

          density_acc[i] =
              static_cast<float>(count) / (4.0f / 3.0f * M_PI * 1000000.0f);
        }
      });
    });

    q.wait_and_throw();

  } catch (sycl::exception const &e) {
    std::cout << "SYCL exception caught: " << e.what() << std::endl;
  }
}

void Simulation::moveSyclEuler(sycl::queue &q) {
  sycl::buffer<Particle, 1> particleBuf(particle_array.get(),
                                        sycl::range<1>(N_PARTICLES));
  sycl::buffer<float, 1> densityBuf(density.get(), sycl::range<1>(N_PARTICLES));

  try {
    q.submit([&](sycl::handler &h) {
      sycl::accessor particles(particleBuf, h, sycl::read_write);
      sycl::accessor density_acc(densityBuf, h, sycl::write_only);

      h.parallel_for(sycl::range<1>(N_PARTICLES), [=](sycl::id<1> idx) {
        int i = idx[0];
        if (i > 0) {

          Particle &p = particles[i];
          // Calculate accelerations
          std::array<num_t, 3> acc = {0, 0, 0};
          for (int j = 0; j < N_PARTICLES; ++j) {
            if (i != j) {
              const Particle &other = particles[j];

              // Euler method
              std::array<num_t, 3> diff;
              num_t distanceSqr = softening * softening;

              for (int k = 0; k < 3; ++k) {
                diff[k] = other.pos[k] - p.pos[k];
                distanceSqr += diff[k] * diff[k];
              }

              num_t distanceInv = sycl::rsqrt(distanceSqr);
              num_t c = sycl::min(distanceSqr / 100000, num_t(1));
              num_t forceGravity = other.mass * (1.0 - sycl::sin(c * (PI / 2)));
              num_t factor =
                  G * forceGravity * distanceInv * distanceInv * distanceInv;

              for (int k = 0; k < 3; k++) {
                acc[k] += diff[k] * factor;
              }
            }
          }

          // Update position
          if (i > 0) {
            for (int k = 0; k < 3; k++) {
              p.vel[k] += acc[k] * dt;
              p.pos[k] += p.vel[k] * dt;
            }
          }

          // Calculate density
          const float radius_square = 10000.0f; // 100.0f * 100.0f
          int count = 0;

          for (int j = 0; j < N_PARTICLES; ++j) {
            if (i == j)
              continue;

            const Particle &other = particles[j];
            float dx = p.pos[0] - other.pos[0];
            float dy = p.pos[1] - other.pos[1];
            float dz = p.pos[2] - other.pos[2];
            float distance_sq = dx * dx + dy * dy + dz * dz;

            if (distance_sq < radius_square) {
              count++;
            }
          }

          density_acc[i] =
              static_cast<float>(count) / (4.0f / 3.0f * M_PI * 1000000.0f);
        }
      });
    });

    q.wait_and_throw();

  } catch (sycl::exception const &e) {
    std::cout << "SYCL exception caught: " << e.what() << std::endl;
  }
}

const Particle *Simulation::get_particles() const {
  return particle_array.get();
}

const float *Simulation::get_density() const { return density.get(); }

void Simulation::outputData(const std::string &filename, int step) {
  std::ofstream outFile(filename, std::ios::app);
  outFile << std::setprecision(6) << std::fixed;

  outFile << "Step: " << step << " " << kineticEnergy << " " << potentialEnergy
          << "\n";
  outFile.close();
}

void Simulation::computeKineticEnergy(sycl::queue &q) {
  sycl::buffer<Particle, 1> particleBuf(particle_array.get(),
                                        sycl::range<1>(N_PARTICLES));
  sycl::buffer<float, 1> energyBuf(&kineticEnergy, sycl::range<1>(1));

  q.submit([&](sycl::handler &h) {
     sycl::accessor particles(particleBuf, h, sycl::read_only);
     sycl::accessor energy(energyBuf, h, sycl::write_only, sycl::no_init);

     h.parallel_for(sycl::range<1>(N_PARTICLES), [=](sycl::id<1> idx) {
       int i = idx[0];
       float velocity_squared = particles[i].vel[0] * particles[i].vel[0] +
                                particles[i].vel[1] * particles[i].vel[1] +
                                particles[i].vel[2] * particles[i].vel[2];
       float particle_energy = 0.5f * particles[i].mass * velocity_squared;

       sycl::atomic_ref<float, sycl::memory_order::relaxed,
                        sycl::memory_scope::device,
                        sycl::access::address_space::global_space>
           atomic_energy(energy[0]);

       atomic_energy.fetch_add(particle_energy);
     });
   }).wait();
}

void Simulation::computePotentialEnergy(sycl::queue &q) {
  sycl::buffer<Particle, 1> particleBuf(particle_array.get(),
                                        sycl::range<1>(N_PARTICLES));
  sycl::buffer<float, 1> energyBuf(&potentialEnergy, sycl::range<1>(1));

  q.submit([&](sycl::handler &h) {
     sycl::accessor particles(particleBuf, h, sycl::read_only);
     sycl::accessor energy(energyBuf, h, sycl::write_only, sycl::no_init);

     h.parallel_for(sycl::range<1>(N_PARTICLES), [=](sycl::id<1> idx) {
       int i = idx[0];
       float local_energy = 0.0f;

       for (int j = 0; j < N_PARTICLES; j++) {
         if (i != j) {
           float dx = particles[j].pos[0] - particles[i].pos[0];
           float dy = particles[j].pos[1] - particles[i].pos[1];
           float dz = particles[j].pos[2] - particles[i].pos[2];
           float distance_squared = dx * dx + dy * dy + dz * dz;
           float distance = sycl::rsqrt(distance_squared);

           if (distance > 0) {
             local_energy +=
                 G * particles[i].mass * particles[j].mass * distance;
           }
         }
       }

       sycl::atomic_ref<float, sycl::memory_order::relaxed,
                        sycl::memory_scope::device,
                        sycl::access::address_space::global_space>
           atomic_energy(energy[0]);
       atomic_energy.fetch_add(local_energy);
     });
   }).wait();
}