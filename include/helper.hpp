#pragma once
#include "constants.hpp"
#include <SDL3/SDL.h>
#include <iostream>

// void set_pixel(SDL_Renderer *rend, int x, int y, Uint8 r, Uint8 g, Uint8 b,
//                Uint8 a) {
//   SDL_SetRenderDrawColor(rend, r, g, b, a);
//   SDL_RenderPoint(rend, x, y);
// }

// void fill_circle(SDL_Renderer *gRenderer, int cx, int cy, int radius, Uint8
// r,
//                  Uint8 g, Uint8 b, Uint8 a) {

//   static const int BPP = 4;

//   for (double dy = 1; dy <= radius; dy += 1.0) {

//     double dx = floor(sqrt((2.0 * radius * dy) - (dy * dy)));
//     int x = cx - dx;
//     SDL_SetRenderDrawColor(gRenderer, r, g, b, a);
//     SDL_RenderLine(gRenderer, cx - dx, cy + dy - radius, cx + dx,
//                    cy + dy - radius);
//     SDL_RenderLine(gRenderer, cx - dx, cy - dy + radius, cx + dx,
//                    cy - dy + radius);
//   }
// }

// // Helper function to calculate color based on distance
// SDL_Color calculate_color(float density, float max_density) {
//   float normalized_density = std::min(density / max_density, 1.0f);

//   // Create a gradient from blue (low density) to red (high density)
//   uint8_t r = static_cast<uint8_t>(255 * normalized_density);
//   uint8_t b = static_cast<uint8_t>(255 * (1.0f - normalized_density));
//   uint8_t g =
//       static_cast<uint8_t>(255 * (1.0f - std::abs(2 * normalized_density -
//       1)));

//   return {r, g, b, 255};
// }
