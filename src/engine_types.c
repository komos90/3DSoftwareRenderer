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

bool doBoxesCollide(Box box1, Box box2)
{

    //SDL_Log("%f, %f, %f, %f, %f, %f", box1.x, box1.y, box1.z, box2.x, box2.y, box2.z);
    //SDL_Log("%d", box1.x + box1.w > box2.x && box1.x < box2.x + box2.w);
    return box1.x + box1.w > box2.x && box1.x < box2.x + box2.w &&
           box1.y + box1.h > box2.y && box1.y < box2.y + box2.h &&
           box1.z + box1.d > box2.z && box1.z < box2.z + box2.d;
}
