#pragma once
#include "constants.hpp"
#include "helper.hpp"
#include "simulation.hpp"
#include <SDL3/SDL.h>
#include <SDL3_image/SDL_image.h>
#include <iomanip>

#include <iostream>
class Screen {
public:
  Screen();
  virtual ~Screen() {
    // Destroy all SDL related objects.
    SDL_DestroyRenderer(m_renderer);
    SDL_DestroyTexture(m_texture);
    SDL_DestroyWindow(m_window);
    SDL_Quit();
  };
  void update(Simulation &);

  bool quit_program();
  void handle_zoom(SDL_Event &event);
  void handle_drag(SDL_Event &event);
  void saveRendererToImage();

  SDL_FPoint world_to_screen(float x, float y) {
    return {(x - m_offsetX) * m_zoom + WIDTH / 2,
            (y - m_offsetY) * m_zoom + HEIGHT / 2};
  }

  SDL_FPoint screen_to_world(int x, int y) {
    return {(x - WIDTH / 2) / m_zoom + m_offsetX,
            (y - HEIGHT / 2) / m_zoom + m_offsetY};
  }
  int frameCount = 0;

private:
  void init_SDL();
  void init_window_renderer();
  void init_texture();

  void set_pixel(SDL_Renderer *rend, int x, int y, Uint8 r, Uint8 g, Uint8 b,
                 Uint8 a) {
    SDL_SetRenderDrawColor(rend, r, g, b, a);
    SDL_RenderPoint(rend, x, y);
  }

  void fill_circle(SDL_Renderer *gRenderer, int cx, int cy, int radius, Uint8 r,
                   Uint8 g, Uint8 b, Uint8 a) {

    for (double dy = 1; dy <= radius; dy += 1.0) {

      double dx = floor(sqrt((2.0 * radius * dy) - (dy * dy)));
      SDL_SetRenderDrawColor(gRenderer, r, g, b, a);
      SDL_RenderLine(gRenderer, cx - dx, cy + dy - radius, cx + dx,
                     cy + dy - radius);
      SDL_RenderLine(gRenderer, cx - dx, cy - dy + radius, cx + dx,
                     cy - dy + radius);
    }
  }

  // Helper function to convert HSV to RGB
  SDL_Color hsv_to_rgb(float h, float s, float v) {
    float c = v * s;
    float x = c * (1 - std::abs(std::fmod(h / 60.0f, 2) - 1));
    float m = v - c;
    float r, g, b;

    if (h >= 0 && h < 60) {
      r = c, g = x, b = 0;
    } else if (h >= 60 && h < 120) {
      r = x, g = c, b = 0;
    } else if (h >= 120 && h < 180) {
      r = 0, g = c, b = x;
    } else if (h >= 180 && h < 240) {
      r = 0, g = x, b = c;
    } else if (h >= 240 && h < 300) {
      r = x, g = 0, b = c;
    } else {
      r = c, g = 0, b = x;
    }

    return {static_cast<Uint8>((r + m) * 255),
            static_cast<Uint8>((g + m) * 255),
            static_cast<Uint8>((b + m) * 255), 255};
  }

  // Enhanced function to calculate color based on density
  SDL_Color calculate_color(float density, float max_density) {
    float normalized_density = std::min(density / max_density, 1.0f);

    float hue = 240 - normalized_density * 240; // Blue (240) to Red (0)
    float saturation = 0.8f + normalized_density * 0.2f;
    float value = 0.7f + normalized_density * 0.3f;

    return hsv_to_rgb(hue, saturation, value);
  }

private:
  SDL_Window *m_window{nullptr};
  SDL_Renderer *m_renderer{nullptr};
  SDL_Texture *m_texture{nullptr};
  SDL_Event m_event;

  float m_zoom = 0.11;
  float m_offsetX = WIDTH / 2;
  float m_offsetY = HEIGHT / 2;
  bool m_isDragging = false;
  int m_lastMouseX = 0;
  int m_lastMouseY = 0;
};