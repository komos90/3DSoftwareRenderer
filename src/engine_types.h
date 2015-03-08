#pragma once

#include <stdint.h>
#include <stdbool.h>

typedef struct 
{
    uint32_t* pixels;
    float* zBuffer;
    int width;
    int height;
} PixelBuffer;

typedef struct
{
    int x;
    int y;
    int z;
} Vector3Int;

typedef struct 
{
    float x;
    float y;
    float z;
} Vector3;

typedef struct
{
    float x;
    float y;
    float z;
    float w;
} Vector4;

typedef struct
{
    float x;
    float y;
    float z;
    float w;
    float h;
    float d;
} Box;

typedef struct
{
    Vector3 vectors[3];
} Triangle;

typedef struct
{
    float values[16];
} Matrix4;

typedef struct
{
    int polyCount;
    Triangle* polygons;
} Mesh;

typedef struct
{
    Vector3 position;
    Vector3 rotation; //TODO Quaternion?
    Vector3 scale;
    Mesh mesh;
    Box collisionBox;
    uint32_t color; //TEMP
} Entity;

typedef struct
{
    int length;
    Entity* data;
} EntityArray;

Vector3 vector3Add(Vector3 vec1, Vector3 vec2);
Vector3 vector3ScalarAdd(Vector3 vec1, float scalar);
Vector3 vector3Floor(Vector3 vec);
bool doBoxesCollide(Box box1, Box box2);
