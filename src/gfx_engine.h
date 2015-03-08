/*
Hand coded 3d software renderer project.
Graphics engine code.
Seoras Macdonald
seoras1@gmail.com
2015
*/
#pragma once

#ifdef __linux__
    #include "SDL2/SDL.h"
#elif _WIN32
    #include "SDL.h"
#endif 

#include "engine_types.h"

//Projection Constants
static const int VIEW_WIDTH    = 1366;
static const int VIEW_HEIGHT   = 768;
static const int Z_FAR         = 500;
static const int Z_NEAR        = 10;
static const float FOV_X       = 10000;//1280;//1.5f;
static const float FOV_Y       = 10000;//960;//1.5f;

//IDEA have pixelBuffer as static variable in gfx_engine, use a function to set it?
void drawRect (SDL_Rect rect, uint32_t color, PixelBuffer pixelBuffer);
void drawVector (Vector3 vector, uint32_t color, PixelBuffer pixelBuffer);
void drawLine (Vector3 start, Vector3 end, uint32_t color, PixelBuffer pixelBuffer);
void rasterizePolygon (Triangle poly, uint32_t color, PixelBuffer pixelBuffer);
void draw (PixelBuffer pixelBuffer, Entity camera, Entity* entityList, int entityCount,
           bool shouldDrawWireframe, bool shouldDrawSurfaces);
Matrix4 mulMatrix4 (Matrix4 mat1, Matrix4 mat2);
Vector3 transform (Matrix4 matrix, Vector3 vector);
Mesh loadMeshFromFile (char* fileName);
