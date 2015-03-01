#pragma once

typedef struct 
{
	uint32_t* pixels;
	int32_t* zBuffer;
	int width;
	int height;
} PixelBuffer;

typedef struct
{
    int x;
    int y;
} Vector2Int;

typedef struct 
{
	float x;
	float y;
	float z;
} Vector3;

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
	//Vector3 origin;
	int polyCount;
	Triangle* polygons;
} Mesh;

typedef struct
{
	Vector3 position;
	Vector3 rotation; //TODO Quaternion?
	Vector3 scale;
	Mesh mesh;
} Entity;

typedef struct
{
    int length;
    Entity* data;
} EntityArray;