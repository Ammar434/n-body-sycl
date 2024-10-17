
#include "../include/screen.hpp"

Screen::Screen() : m_window{nullptr}, m_renderer{nullptr}, m_texture{nullptr} {
  // Initialize all required SDL functionality and create SDL objects.
  init_SDL();
  init_window_renderer();
  init_texture();

  // Check if requirement are meet to save image
  if (IMG_Init(IMG_INIT_JPG) < 0) {
    std::cerr << "IMG_Init Error: " << SDL_GetError() << std::endl;
  }
}

void Screen::init_SDL() {
  // Initialize SDL library. Returns 0 if successful.
  // Must be called before using any SDL functionality.
  if (!SDL_Init(SDL_INIT_VIDEO)) {
    std::cout << "SDL_Init Error: \n" << SDL_GetError() << std::endl;
    exit(1);
  }
}
void Screen::init_window_renderer() {

  if (!SDL_CreateWindowAndRenderer("Particle Simulation ", WIDTH, HEIGHT,
                                   SDL_WINDOW_HIDDEN, &m_window, &m_renderer)) {
    std::cout << "SDL_CreateWindow or Renderer Error: \n"
              << SDL_GetError() << std::endl;
    exit(2);
  }
}

void Screen::init_texture() {
  m_texture = SDL_CreateTexture(m_renderer, SDL_PIXELFORMAT_RGBA8888,
                                SDL_TEXTUREACCESS_STREAMING, WIDTH, HEIGHT);
  if (!m_texture) {
    std::cout << "SDL_CreateTexture Error: \n" << SDL_GetError() << std::endl;
    SDL_DestroyRenderer(m_renderer);
    SDL_DestroyWindow(m_window);
    exit(3);
  }
}

void Screen::update(Simulation &simulation) {

  SDL_SetRenderDrawColor(m_renderer, 0, 0, 0, 255);
  SDL_RenderClear(m_renderer);

  const Particle *const particles = simulation.get_particles();
  const float *const density = simulation.get_density();
  const float *max_density = std::max_element(density, density + N_PARTICLES);

  for (size_t i = 0; i < N_PARTICLES; ++i) {
    const Particle &particle = particles[i];
    SDL_FPoint screenPos = world_to_screen(particle.pos[0], particle.pos[1]);

    // Get color based on density
    SDL_Color color = calculate_color(density[i], *max_density);
    if (i == 0) { // Central particle (supermassive black hole)
      fill_circle(m_renderer, screenPos.x, screenPos.y, 10, 255, 255, 255, 255);
    } else {
      // Use particle mass to determine size (larger mass = larger particle)
      // float size = 2 + std::exp(particle.mass / N_PARTICLES);
      fill_circle(m_renderer, screenPos.x, screenPos.y, 5, color.r, color.g,
                  color.b, color.a);
    }
  }

  SDL_RenderPresent(m_renderer);
  saveRendererToImage();

  const float initial_zoom_speed = 0.01f;
  const float decay_rate = 0.005f;
  float zoom_factor = initial_zoom_speed * std::exp(-decay_rate * frameCount);
  m_zoom /= (1 + zoom_factor);

  frameCount++;
}

// Function to save SDL renderer content as a JPEG file
void Screen::saveRendererToImage() {
  std::ostringstream oss;
  oss << "images/frame_" << std::setw(6) << std::setfill('0') << frameCount
      << ".jpg";

  SDL_Surface *surface = SDL_RenderReadPixels(m_renderer, NULL);
  if (!surface) {
    std::cerr << "SDL_RenderReadPixels Error: " << SDL_GetError() << std::endl;
    return;
  }

  // Save the surface as a JPEG file
  if (!IMG_SaveJPG(surface, oss.str().c_str(), 70)) {
    std::cerr << "IMG_SaveJPG Error: " << SDL_GetError() << std::endl;
  }

  // Free the surface
  SDL_DestroySurface(surface);
}

bool Screen::quit_program() {
  // Check for SDL events. If the window is closed, quit the program.
  while (SDL_PollEvent(&m_event)) {
    switch (m_event.type) {
    case SDL_EVENT_QUIT:
      return true;
    case SDL_EVENT_MOUSE_WHEEL:
      handle_zoom(m_event);
      break;
    case SDL_EVENT_MOUSE_BUTTON_DOWN:
    case SDL_EVENT_MOUSE_BUTTON_UP:
    case SDL_EVENT_MOUSE_MOTION:
      handle_drag(m_event);
      break;
    }
  }
  return false;
}

void Screen::handle_zoom(SDL_Event &event) {
  float zoom_speed = 0.1f;
  SDL_FPoint mouseWorldPos =
      screen_to_world(event.wheel.mouse_x, event.wheel.mouse_y);

  if (event.wheel.y > 0) {
    m_zoom *= (1 + zoom_speed);
  } else if (event.wheel.y < 0) {
    m_zoom /= (1 + zoom_speed);
  }

  // m_zoom = std::max(m_zoom, 0.1f);

  SDL_FPoint newMouseWorldPos =
      screen_to_world(event.wheel.mouse_x, event.wheel.mouse_y);
  m_offsetX += mouseWorldPos.x - newMouseWorldPos.x;
  m_offsetY += mouseWorldPos.y - newMouseWorldPos.y;
}

void Screen::handle_drag(SDL_Event &event) {
  switch (event.type) {
  case SDL_EVENT_MOUSE_BUTTON_DOWN:
    if (event.button.button == SDL_BUTTON_LEFT) {
      m_isDragging = true;
      m_lastMouseX = event.button.x;
      m_lastMouseY = event.button.y;
    }
    break;
  case SDL_EVENT_MOUSE_BUTTON_UP:
    if (event.button.button == SDL_BUTTON_LEFT) {
      m_isDragging = false;
    }
    break;
  case SDL_EVENT_MOUSE_MOTION:
    if (m_isDragging) {
      int dx = event.motion.x - m_lastMouseX;
      int dy = event.motion.y - m_lastMouseY;
      m_offsetX -= dx / m_zoom;
      m_offsetY -= dy / m_zoom;
      m_lastMouseX = event.motion.x;
      m_lastMouseY = event.motion.y;
    }
    break;
  }
}
