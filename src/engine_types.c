#include <math.h>

#ifdef __linux__
    #include "SDL2/SDL.h"
#elif _WIN32
    #include "SDL.h"
#endif

#include "engine_types.h"


Vector3 vector3Add(Vector3 vec1, Vector3 vec2)
{
    Vector3 tmp;
    tmp.x = vec1.x + vec2.x;
    tmp.y = vec1.y + vec2.y;
    tmp.z = vec1.z + vec2.z;
    return tmp;
}

Vector3 vector3Sub(Vector3 vec1, Vector3 vec2)
{
    Vector3 tmp;
    tmp.x = vec1.x - vec2.x;
    tmp.y = vec1.y - vec2.y;
    tmp.z = vec1.z - vec2.z;
    return tmp;
}

Vector3 vector3ScalarAdd(Vector3 vec1, float scalar)
{
    Vector3 tmp;
    tmp.x = vec1.x + scalar;
    tmp.y = vec1.y + scalar;
    tmp.z = vec1.z + scalar;
    return tmp;
}

Vector3 vector3Floor(Vector3 vec)
{
    Vector3 tmp = {floor(vec.x), floor(vec.y), floor(vec.z)};
    return tmp;
}

float vector3Dot(Vector3 vec1, Vector3 vec2)
{
    return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
}

Vector3 vector3Cross(Vector3 vec1, Vector3 vec2)
{
    Vector3 returnVec = {
        .x=vec1.y * vec2.z - vec1.z * vec2.y,
        .y=vec1.x * vec2.z - vec1.z * vec2.x,
        .z=vec1.x * vec2.y - vec1.y * vec2.x
    };
    return returnVec;
}

float vector3Abs(Vector3 vec)
{
    return sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

Vector3 vector3Normalize(Vector3 vec)
{
    float magnitude = sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
    vec.x /= magnitude;
    vec.y /= magnitude;
    vec.z /= magnitude;
    return vec;
}

Vector3 getTriangleNormal(Triangle tri) {
    Vector3 vec1 = vector3Sub(tri.vectors[0], tri.vectors[1]);
    Vector3 vec2 = vector3Sub(tri.vectors[0], tri.vectors[2]);
    return vector3Normalize(vector3Cross(vec1, vec2));
}

Vector3 getTriangleCenter(Triangle tri)
{
    Vector3 sumVec = vector3Add(vector3Add(tri.vectors[0], tri.vectors[1]), tri.vectors[2]);
    sumVec.x /= 3;
    sumVec.y /= 3;
    sumVec.z /= 3;
    return sumVec;
}

bool doBoxesCollide(Box box1, Box box2)
{

    //SDL_Log("%f, %f, %f, %f, %f, %f", box1.x, box1.y, box1.z, box2.x, box2.y, box2.z);
    //SDL_Log("%d", box1.x + box1.w > box2.x && box1.x < box2.x + box2.w);
    return box1.x + box1.w > box2.x && box1.x < box2.x + box2.w &&
           box1.y + box1.h > box2.y && box1.y < box2.y + box2.h &&
           box1.z + box1.d > box2.z && box1.z < box2.z + box2.d;
}
