/*
Hand coded 3d software renderer project.
Seoras Macdonald
seoras1@gmail.com
2015
*/
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#ifdef __linux__
    #include <SDL2/SDL.h>
    #define M_PI 3.14159265358979323846
#elif _WIN32
    #include <SDL.h>
#endif

#include "engine_types.h"
#include "load_level.h"
#include "meshes.h"
#include "gfx_engine.h"

static const int SCREEN_WIDTH  = 854;//300;//640;
static const int SCREEN_HEIGHT = 480;//300;//480;

//Temp Globals

bool collidedWithSomething(Entity entity, EntityArray otherEntities)
{
    for (int i = 0; i < otherEntities.length; i++)
    {
        Entity otherEntity = otherEntities.data[i];

        Box box1;
        box1.x = entity.collisionBox.x + entity.position.x;
        box1.y = entity.collisionBox.y + entity.position.y;
        box1.z = entity.collisionBox.z + entity.position.z;
        box1.w = entity.collisionBox.w;
        box1.h = entity.collisionBox.h;
        box1.d = entity.collisionBox.d;

        Box box2;
        box2.x = otherEntity.collisionBox.x + otherEntity.position.x;
        box2.y = otherEntity.collisionBox.y + otherEntity.position.y;
        box2.z = otherEntity.collisionBox.z + otherEntity.position.z;
        box2.w = otherEntity.collisionBox.w;
        box2.h = otherEntity.collisionBox.h;
        box2.d = otherEntity.collisionBox.d;

        //SDL_Log("%f, %f, %f, %f, %f, %f", box1.w, box1.h, box1.d, box2.w, box2.h, box2.d);

        if (doBoxesCollide(box1, box2))
        {
            //SDL_Log("Collide.");
            return true;
        }
    }
    return false;
}

int main( int argc, char* args[] )
{
    //Initialise SDL ====
    if( SDL_Init( SDL_INIT_EVERYTHING) < 0 )
    {
        printf( "SDL could not initialize! SDL_Error: %s\n", SDL_GetError() );
        return 0;
    }
    SDL_Window* window = SDL_CreateWindow(
        "Pixel buffer Playground :P", SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, 
        SDL_WINDOW_FULLSCREEN_DESKTOP);
    if (window == NULL)
    {
        printf( "Window could not be created! SDL_Error: %s\n", SDL_GetError() );
        return 0;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);
    SDL_Texture* screenTexture = SDL_CreateTexture(
        renderer, SDL_PIXELFORMAT_ARGB8888,
        SDL_TEXTUREACCESS_STREAMING, SCREEN_WIDTH,
        SCREEN_HEIGHT);
    //grab cursor
    SDL_SetRelativeMouseMode(true);

    //Allocate pixel and z-buffer
    uint32_t* pixels = (uint32_t*) malloc(SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(uint32_t));
    float*  zBuffer = (float*) malloc(SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(float));
    PixelBuffer pixelBuffer = {pixels, zBuffer, SCREEN_WIDTH, SCREEN_HEIGHT};

    initGfxEngine();

    //Initialise Meshes and Entities ====
    //Load meshes
    cubeMesh = loadMeshFromFile("../res/meshes/cube.raw");
    plane = loadMeshFromFile("../res/meshes/plane.raw");
    monkey = loadMeshFromFile("../res/meshes/monkey.raw");
    monkeyHd = loadMeshFromFile("../res/meshes/monkeyhd.raw");
    monkeySuperHd = loadMeshFromFile("../res/meshes/monkeysuperhd.raw");

    //Load level from file and add level entities to entity list
    //Level level = loadLevel("../res/levels/level2.lvl");
    //EntityArray entities = createLevelEntities(level);
    EntityArray entities;
    entities.data = (Entity*)malloc(sizeof(Entity));
    Entity temp = {.position = {400, 0, 200}, .mesh=monkeyHd, .scale={100, 100, 100}, .rotation={M_PI/2,0,0}};
    entities.data[0] = temp;
    entities.length = 1;

    //Initialise Entities
    Entity camera = {.position={100, 0, 100}, .rotation={0, M_PI/2},
                     .collisionBox={-10,-1000,-10,20,2000,20}};

    //Get input devices' states
    SDL_Joystick* gamePad = SDL_JoystickOpen(0);
    const uint8_t* keyState = SDL_GetKeyboardState(NULL);    
    
    bool running = true;
    bool paused = false;
    bool shouldDrawWireframe = false;
    bool shouldDrawSurfaces = true;

    //tmp frame rate counter
    int framerateDisplayDelayCounter = 0;

    //Main Loop ====
    while(running)    
    {
        int frameStartTime = SDL_GetTicks();
        
        //SDL Event Loop
        SDL_Event e;
        while (SDL_PollEvent(&e))
        {
            switch(e.type)
            {
            case SDL_QUIT:
                running = false;
                break;
            case SDL_KEYDOWN:
                switch(e.key.keysym.sym)
                {
                    case SDLK_RETURN:
                        if (e.key.keysym.mod & KMOD_ALT)
                        {
                            uint32_t flags = SDL_GetWindowFlags(window);
                            if (flags & SDL_WINDOW_FULLSCREEN ||
                                flags & SDL_WINDOW_FULLSCREEN_DESKTOP)
                            {
                                SDL_SetWindowFullscreen(window, 0);
                            }
                            else
                            {
                                SDL_SetWindowFullscreen(
                                    window, SDL_WINDOW_FULLSCREEN_DESKTOP);
                            }
                        }
                        break;
                    case SDLK_SPACE:
                        camera.position.y += 50;
                        break;
                    case SDLK_1:
                        shouldDrawWireframe = !shouldDrawWireframe;
                        break;
                    case SDLK_2:
                        shouldDrawSurfaces = !shouldDrawSurfaces;
                        break;
                }
                break;
            case SDL_MOUSEMOTION:
                camera.rotation.y -= e.motion.xrel * 0.001;
               
                //Keep look up and down bounded 
                float newCamRotX = camera.rotation.x - e.motion.yrel * 0.001;
                if (newCamRotX > -M_PI/2 && newCamRotX < M_PI/2)
                    camera.rotation.x = newCamRotX;
                break;
            }
        }
        //Joystick input
        //HACK Should decouple Input and Game Logic
        {
            const int JOYSTICK_DEAD_ZONE = 8000;
            int moveVel = 3;
            if( SDL_JoystickGetAxis(gamePad, 0) < -JOYSTICK_DEAD_ZONE )
            {
                camera.position.z += moveVel * cosf(camera.rotation.y + M_PI/2);
                camera.position.x += moveVel * sinf(camera.rotation.y + M_PI/2);
            }
            //Right of dead zone
            else if( SDL_JoystickGetAxis(gamePad, 0) > JOYSTICK_DEAD_ZONE )
            {
                camera.position.z += moveVel * cosf(camera.rotation.y - M_PI/2);
                camera.position.x += moveVel * sinf(camera.rotation.y - M_PI/2);
            }
            //Left of dead zone
            if( SDL_JoystickGetAxis(gamePad, 1) < -JOYSTICK_DEAD_ZONE )
            {
                camera.position.z += moveVel * cosf(camera.rotation.y);
                camera.position.x += moveVel * sinf(camera.rotation.y);
            }
            //Right of dead zone
            else if( SDL_JoystickGetAxis(gamePad, 1) > JOYSTICK_DEAD_ZONE )
            {
                camera.position.z += moveVel * cosf(camera.rotation.y + M_PI);
                camera.position.x += moveVel * sinf(camera.rotation.y + M_PI);
            }
            //Left of dead zone
            if( SDL_JoystickGetAxis(gamePad, 2) < -JOYSTICK_DEAD_ZONE )
            {
                camera.rotation.y += 0.04 * -SDL_JoystickGetAxis(gamePad, 2) / 32767.f;
            }
            //Right of dead zone
            else if( SDL_JoystickGetAxis(gamePad, 2) > JOYSTICK_DEAD_ZONE )
            {
                camera.rotation.y -= 0.04 * SDL_JoystickGetAxis(gamePad, 2) / 32767.f;
            }
        }
        //Keyboard Input
        { 
            Vector3 oldCameraPos = camera.position;
            int moveVel = 3; 
            if (keyState[SDL_SCANCODE_A])
            {
                camera.position.z += moveVel * cosf(camera.rotation.y + M_PI/2);
                camera.position.x += moveVel * sinf(camera.rotation.y + M_PI/2);
            }
            if (keyState[SDL_SCANCODE_D])
            {
                camera.position.z += moveVel * cosf(camera.rotation.y - M_PI/2);
                camera.position.x += moveVel * sinf(camera.rotation.y - M_PI/2);
            }
            if (keyState[SDL_SCANCODE_S])
            {
                camera.position.z += moveVel * cosf(camera.rotation.y + M_PI);
                camera.position.x += moveVel * sinf(camera.rotation.y + M_PI);
            }
            if (keyState[SDL_SCANCODE_W])
            {
                camera.position.z += moveVel * cosf(camera.rotation.y);
                camera.position.x += moveVel * sinf(camera.rotation.y);
            }
            if (keyState[SDL_SCANCODE_LEFT])
            {
                camera.rotation.y += 0.02;
            }
            if (keyState[SDL_SCANCODE_RIGHT])
            {
                camera.rotation.y -= 0.02;
            }
            //For each other entity do collision detection. Only update newPos if
            //player does not collide with anything.
            if (collidedWithSomething(camera, entities))
            {
                camera.position = oldCameraPos;
            }
        }
        
        if(!paused)
        {
            //entities.data[0].rotation.x += 0.01;
            entities.data[0].rotation.y += 0.01;
        }    
       
        //Send game entities to gfx engine to be rendered 
        draw(pixelBuffer, camera, entities.data, entities.length, shouldDrawWireframe, shouldDrawSurfaces);

        //Render the pixel buffer to the screen
        SDL_UpdateTexture(screenTexture, NULL, pixelBuffer.pixels, SCREEN_WIDTH * sizeof(uint32_t));        
        SDL_RenderCopy(renderer, screenTexture, NULL, NULL);
        SDL_RenderPresent(renderer);

        //Clear the pixel buffer
        memset((void*)pixelBuffer.pixels, 200, SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(uint32_t));
        //Clear the z-buffer
        for (int i = 0; i < pixelBuffer.width * pixelBuffer.height; i++)
        {
            pixelBuffer.zBuffer[i] = INT_MIN;
        }

        //Lock to 60 fps
        int delta = SDL_GetTicks() - frameStartTime;
        if (delta < 1000/60)
        {
            SDL_Delay(1000/60 - delta);
        }
        framerateDisplayDelayCounter = (framerateDisplayDelayCounter + 1) % 30;
        //if (framerateDisplayDelayCounter == 0)
            //SDL_Log("FPS: %f", 1000.f/(SDL_GetTicks() - frameStartTime));
    }
    return 0;
}
