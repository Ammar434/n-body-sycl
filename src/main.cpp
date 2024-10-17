#include "../include/constants.hpp"
#include "../include/screen.hpp"
#include <sycl/sycl.hpp>

int main() {

  // Initialize SYCL queue
  sycl::queue q{sycl::gpu_selector_v};
  sycl::device dev = q.get_device();

  std::cout << "device name : "
            << q.get_device().get_info<sycl::info::device::name>() << "\n";
  std::cout << "global_mem_size : "
            << q.get_device().get_info<sycl::info::device::global_mem_size>()
            << "\n";
  std::cout << "local_mem_size : "
            << q.get_device().get_info<sycl::info::device::local_mem_size>()
            << "\n";
  std::cout
      << "max_work_group_size: "
      << q.get_device().get_info<sycl::info::device::max_work_group_size>()
      << "\n";

  Screen screen;
  Simulation simulation;

  sycl::buffer<Particle, 1> display_buf(simulation.get_particles(),
                                        sycl::range<1>(N_PARTICLES));

  while (!screen.quit_program() && screen.frameCount < 100) {
    // Update SDL window with new particle positions and colors.
    screen.update(simulation);

    // Manipulate particle positions for next iteration.
    simulation.moveSyclRK4(q);

    simulation.computeKineticEnergy(q);
    simulation.computePotentialEnergy(q);

    std::cout << "Generated " << screen.frameCount << std::endl;
    simulation.outputData("simulation_data_rk4.txt", screen.frameCount);
  }
  return 0;
}
