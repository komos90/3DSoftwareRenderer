/*
Hand coded 3d software renderer project.
Graphics engine code.
Seoras Macdonald
seoras1@gmail.com
2015
*/
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#ifdef __linux__
    #include <SDL2/SDL.h>
#elif __Win32
    #include <SDL.h>
#endif

#include "gfx_engine.h"
#include "engine_types.h"

void drawRect(SDL_Rect rect, uint32_t color, PixelBuffer pixelBuffer)
{
    for (int yy = rect.y; yy < rect.y + rect.h; ++yy)
    {
        for (int xx = rect.x; xx < rect.x + rect.w; ++xx)
        {
            pixelBuffer.pixels[yy * pixelBuffer.width + xx] = color;
        }
    }
}

//TODO use floats for depth buffer
//TODO Get proper per pixel z values for each polygon
void drawVector(Vector3 vector, uint32_t color, PixelBuffer pixelBuffer)
{
    if (vector.z < pixelBuffer.zBuffer[(int)vector.y * pixelBuffer.width + (int)vector.x])
    {
        pixelBuffer.pixels[(int)vector.y  * pixelBuffer.width + (int)vector.x] = color;
        pixelBuffer.zBuffer[(int)vector.y  * pixelBuffer.width + (int)vector.x] = vector.z;
    }
}

void drawLine(Vector3 start, Vector3 end, uint32_t color,
              PixelBuffer pixelBuffer)
{
    int x0 = (int)start.x;
    int y0 = (int)start.y;
    int x1 = (int)end.x;
    int y1 = (int)end.y;
    
    //Distance between x0 and x1, y0 and y1
    int dx = abs(x1 - x0);
    int dy = abs(y1 - y0);
    //Determine whether to step backwards or forwards through the line
    int sx = (x0 < x1) ? 1 : -1;
    int sy = (y0 < y1) ? 1 : -1;

    int err = dx - dy;

    for(;;)
    {
        //Draw the current pixel
        Vector3 tmp = {x0, y0};
        drawVector(tmp, color, pixelBuffer);

        //Break if we reach the end of the line.
        if ((x0 == x1) && (y0 == y1)) break;

        int e2 = 2 * err;
        //Step x
        if (e2 > -dy)
        { 
            err -= dy;
            x0 += sx;
        }
        //Step y
        if (e2 < dx)
        {
            err += dx;
            y0 += sy;
        }
    }
}

void rasterizePolygon(Triangle poly, uint32_t color,
                      PixelBuffer pixelBuffer)
{
    int topIndex = 0;
    int leftIndex = 0;
    int rightIndex = 0;
    //Find top vertex
    for (int i = 1; i < 3; i++)
    {
        if (poly.vectors[i].y < poly.vectors[topIndex].y)
        {
            topIndex = i;
        }
        else if (poly.vectors[i].y == poly.vectors[topIndex].y)
        {
            if (poly.vectors[i].x < poly.vectors[(i+1)%3].x ||
                poly.vectors[i].x < poly.vectors[(i+2)%3].x)
            {
                topIndex = i;
            }
        }
    }
    //Find left and right vertices
    leftIndex  = (topIndex + 2) % 3;
    rightIndex  = (topIndex + 1) % 3;
    
    //Initilise vertices for triangle drawing
    Vector3Int topLO = {(int)poly.vectors[topIndex].x,
                        (int)poly.vectors[topIndex].y,
                        (int)poly.vectors[topIndex].z};
    Vector3Int left  = {(int)poly.vectors[leftIndex].x,
                        (int)poly.vectors[leftIndex].y,
                        (int)poly.vectors[leftIndex].z};
    Vector3Int right = {(int)poly.vectors[rightIndex].x,
                        (int)poly.vectors[rightIndex].y,
                        (int)poly.vectors[rightIndex].z};
    Vector3Int topRO  = topLO;
    Vector3Int topR  = topLO;
    Vector3Int topL  = topLO;

    //Line drawing variables for left line
    //NOTE Removed syL and syR as the lines are always increasing
    // i.e. syL = syR = 1 in all cases 
    int dxL = abs(left.x-topL.x);
    int dyL = abs(left.y-topL.y);
    int sxL = topL.x<left.x ? 1 : -1;
    int errL = (dxL>dyL ? dxL : -dyL)/2;
    int e2L;

    //Line drawing variables for right line
    int dxR = abs(right.x-topR.x);
    int dyR = abs(right.y-topR.y);
    int sxR = topR.x<right.x ? 1 : -1;
    int errR = (dxR>dyR ? dxR : -dyR)/2;
    int e2R;

    //z-buffer
    float zL = topL.z;
    float zR = topL.z;

    for(;;)
    {
        //Draw current scanline
        for(;;)
        {
            //Handle breakpoint on right line
            if (topR.x == right.x && topR.y == right.y)
            {
                topRO = right;
                right = left;
                dxR = abs(right.x-topR.x);
                dyR = abs(right.y-topR.y);
                sxR = topR.x<right.x ? 1 : -1;
                errR = (dxR>dyR ? dxR : -dyR)/2;
            }
            if (topL.y == topR.y) break;
            e2R = errR;
            if (e2R >-dxR) { errR -= dyR; topR.x += sxR; }
            if (e2R < dyR) { errR += dxR; topR.y += 1; }
            
            //calculate zR value by interpolating between top.z and right.z
            if(topRO.y == right.y) 
                zR = topRO.z;
            else
                zR = topRO.z + (right.z - topRO.z) * ((float)(topR.y - right.y) / (float)(topRO.y - right.y));
        }
        //Fill scanline
        for (int i = topL.x; i < topR.x; i++)
        {
            //calculate z value by interpolating between zL and zR
            float curZ = zL + (zR -zL) * ((float)(i - topL.x) / (float)(topR.x - topL.x));
            Vector3 tmp = {(int)i, (int)topL.y, curZ};
            //uint32_t outCol = ((uint8_t)(curZ)) + ((uint8_t)(curZ) >> 8) + ((uint8_t)(curZ) >> 16); 
            drawVector(tmp, color, pixelBuffer);
        } 
        if (topL.x == left.x && topL.y == left.y)
        {
            if (right.y <= topL.y)
                break;
            else
            {
                topLO = left;
                left = right;
                dxL = abs(left.x-topL.x);
                dyL = abs(left.y-topL.y);
                sxL = topL.x<left.x ? 1 : -1;
                errL = (dxL>dyL ? dxL : -dyL)/2;
            }
        }
        e2L = errL;
        if (e2L >-dxL) { errL -= dyL; topL.x += sxL; }
        if (e2L < dyL) { errL += dxL; topL.y += 1; }

        //calculate zL value by interpolating between top.z and left.z
        if (topLO.y == left.y)
            zL = topLO.z;
        else
            zL = topLO.z + (left.z - topLO.z) * ((float)(topL.y - left.y) / (float)(topLO.y - left.y));
    } 
}

Matrix4 mulMatrix4(Matrix4 mat1, Matrix4 mat2)
{
    Matrix4 result = {{0}};
    for (int i = 0; i < 16; i++)
        for (int j = 0; j < 4; j++)
            result.values[i] += mat1.values[(i/4)*4 + j] * mat2.values[i%4 + j*4];
    return result;
}

Vector3 transform(Matrix4 matrix, Vector3 vector) 
{
    Vector3 result;
    float w = 1.f;
    result.x = matrix.values[0] * vector.x + matrix.values[1] * vector.y +
               matrix.values[2] * vector.z + matrix.values[3] * 1.f;
    result.y = matrix.values[4] * vector.x + matrix.values[5] * vector.y +
               matrix.values[6] * vector.z + matrix.values[7] * 1.f;
    result.z = matrix.values[8] * vector.x + matrix.values[9] * vector.y +
               matrix.values[10] * vector.z + matrix.values[11] * 1.f;
           w = matrix.values[12] * vector.x + matrix.values[13] * vector.y +
               matrix.values[14] * vector.z + matrix.values[15] * 1.f;

    if (w != 0.f && w != 1.f) {
        result.x /= w;
        result.y /= w;
        result.z /= w;
    }
    //SDL_Log("%f, %f, %f, %f", result.x, result.y, result.z, result.w);
    return result;
}

bool isInsideRect(Vector3 vector, SDL_Rect rect)
{
    return vector.x > rect.x && vector.x < rect.x + rect.w &&
           vector.y > rect.y && vector.y < rect.y + rect.h &&
           vector.z > 100.f;
}

bool isOnRectBorder(Vector3 vector, SDL_Rect rect)
{
    return ((vector.x == rect.x || vector.x == rect.x + rect.w) &&
            vector.y >= rect.y && vector.y <= rect.y + rect.h) ||
           ((vector.y == rect.y || vector.y == rect.y + rect.h) &&
            vector.x >= rect.x && vector.x <= rect.x + rect.w);
}

float get2DMagnitude(Vector3 vector)
{
    return sqrt(vector.x * vector.x + vector.y * vector.y);
}

Vector3 getIntersect(Vector3 start, Vector3 end, SDL_Rect rect)
{
    int x0 = (int)start.x;
    int y0 = (int)start.y;
    int x1 = (int)end.x;
    int y1 = (int)end.y;
    
    //Distance between x0 and x1, y0 and y1
    int dx = abs(x1 - x0);
    int dy = abs(y1 - y0);
    //Determine whether to step backwards or forwards through the line
    int sx = (x0 < x1) ? 1 : -1;
    int sy = (y0 < y1) ? 1 : -1;

    int err = dx - dy;

    for(;;)
    {
        //Check for rect boundary
        Vector3 tmp = {x0, y0, 0};
        if (isOnRectBorder(tmp, rect))
        {
            if (end.z == start.z)
                tmp.z = start.z;
            else
                tmp.z = start.z + (end.z - start.z) * ((get2DMagnitude(tmp) - (get2DMagnitude(start))) / (get2DMagnitude(end) - get2DMagnitude(start)));
                //SDL_Log ("%f, %f, %f", start.z, tmp.z, end.z);
            return tmp;
        }

        //Break if we reach the end of the line.
        if ((x0 == x1) && (y0 == y1)) break;

        int e2 = 2 * err;
        //Step x
        if (e2 > -dy)
        { 
            err -= dy;
            x0 += sx;
        }
        //Step y
        if (e2 < dx)
        {
            err += dx;
            y0 += sy;
        }
    }
    Vector3 error = {0, 0, 0};
    return error;
}

//Has some temp debug parameters
void draw(
    PixelBuffer pixelBuffer, Entity* camera,
    Entity* entityList, int entityCount,
    bool shouldDrawWireframe,
    bool shouldDrawSurfaces)
{
    for (int k = 0; k < entityCount; k++) 
    {
        Entity* entity = &entityList[k];
        uint32_t lineColor = 0xffffffff;
        uint32_t fillColor = 0x00555555;

        //Scale Matrix
        Vector3 s = entity->scale;
        Matrix4 scaleMat = {{
            s.x, 0  , 0  , 0,
            0  , s.y, 0  , 0,
            0  , 0  , s.z, 0,
            0  , 0  , 0  , 1
        }};
        //ZAxis Rotation Matrix
        float a = entity->rotation.z;
        Matrix4 zRotMat = {{
             cosf(a), sinf(a), 0, 0,
            -sinf(a), cosf(a), 0, 0,
             0      , 0      , 1, 0,
             0      , 0      , 0, 1
        }};
        //YAxis Rotation Matrix
        a = entity->rotation.y;
        Matrix4 yRotMat ={{
             cosf(a), 0, sinf(a), 0,
             0      , 1, 0      , 0,
            -sinf(a), 0, cosf(a), 0,
             0      , 0, 0      , 1
        }};
        //XAxis Rotation Matrix
        a = entity->rotation.x;
        Matrix4 xRotMat = {{
            1,  0      , 0      , 0,
            0,  cosf(a), sinf(a), 0,
            0, -sinf(a), cosf(a), 0,
             0,  0      , 0      , 1
        }};
        //World translation matrix
        Vector3 p = entity->position;
        Matrix4 worldTranslate = {{
            1, 0, 0, p.x,
            0, 1, 0, p.y,
            0, 0, 1, p.z,
             0, 0, 0, 1
        }};
        //Camera Z Rotation Matrix
        a = camera->rotation.z;
        Matrix4 cameraZRotation = {{
             cosf(-a), sinf(-a), 0, 0,
            -sinf(-a), cosf(-a), 0, 0,
             0       , 0       , 1, 0,
             0       , 0       , 0, 1
        }};
        //Camera Y Rotation Matrix
        a = camera->rotation.y;
        Matrix4 cameraYRotation = {{
             cosf(-a), 0, sinf(-a), 0,
             0       , 1, 0       , 0,
            -sinf(-a), 0, cosf(-a), 0,
             0       , 0, 0       , 1
        }};
        //Camera X Rotation
        a = camera->rotation.x;
        Matrix4 cameraXRotation = {{
            1,  0       , 0       , 0,
            0,  cosf(-a), sinf(-a), 0,
            0, -sinf(-a), cosf(-a), 0,
             0,  0       , 0       , 1
        }};
        //Camera translate Matrix    
        p = camera->position;
        Matrix4 cameraTranslate = {{
            1, 0, 0, p.x,
            0, 1, 0, p.y,
            0, 0, 1, p.z,
             0, 0, 0, 1
        }};
        
        //perspective projection matrix
        float e11 = atanf((FOV_X/VIEW_WIDTH)/2);
        float e22 = atanf((FOV_Y/VIEW_HEIGHT)/2);
        float e33 = -(Z_FAR + Z_NEAR) / (Z_FAR - Z_NEAR);
        float e34 = (-2 * (Z_FAR*Z_NEAR))/(Z_FAR - Z_NEAR);
        Matrix4 perspectiveProjection = {{
            e11, 0  , 0  , 0  ,
            0  , e22, 0  , 0  ,
            0  , 0  , e33, e34,
            0  , 0  , -1 , 0
        }};
        
        int w = pixelBuffer.width;
        int h = pixelBuffer.height;
        float d = 10000000.f;
        Matrix4 correctForScreen = {{
            w/2, 0  ,0  , w/2,
            0  , h/2,0  , h/2,
            0  , 0  ,d/2, d/2,
             0  , 0  ,0  , 1
        }};
        
        //Combine matrices into one transformation matrix
        //Model Space -> World Space
        Matrix4 finalTransform = mulMatrix4(xRotMat, scaleMat);
        finalTransform = mulMatrix4(yRotMat, finalTransform);    
        finalTransform = mulMatrix4(zRotMat, finalTransform);    
        finalTransform = mulMatrix4(worldTranslate, finalTransform);
        //World Space -> View Space    
        finalTransform = mulMatrix4(cameraTranslate, finalTransform);    
        finalTransform = mulMatrix4(cameraYRotation, finalTransform);    
        finalTransform = mulMatrix4(cameraXRotation, finalTransform);    
        //finalTransform = mulMatrix4(cameraZRotation, finalTransform);    
        //View Space -> Projection Space
        finalTransform = mulMatrix4(perspectiveProjection, finalTransform);    

        //For each polygon
        for (int i = 0; i < entity->mesh.polyCount; i++)
        {    
            Triangle displayPoly;

            //For each vertex
            for (int j = 0; j < 3; j++)
            {
                //Convert triangular Vector3 polygon to up to five vector4s
                displayPoly.vectors[j].x = entityList[k].mesh.polygons[i].vectors[j].x;
                displayPoly.vectors[j].y = entityList[k].mesh.polygons[i].vectors[j].y;
                displayPoly.vectors[j].z = entityList[k].mesh.polygons[i].vectors[j].z;

                //Apply all transformations =====
                displayPoly.vectors[j] = transform(finalTransform, displayPoly.vectors[j]);
                
                //Projection Space -> Screen Friendly
                displayPoly.vectors[j] = transform(correctForScreen, displayPoly.vectors[j]);

            }

            //1.clip by screen rectangle
            //Parameters Vector3 vector, SDL_Rect screenRect
            //Returns Vector3[5] pentagon, int polygonSize
            Vector3 clippedPolygon[5];
            int clippedPolyTriangleCount = 0;
            {
                //for each vertex in Triangle
                //1.Check vertex[i] and vertex[(i+1) % polyLen] is inside
                //the screenRect and set is1Inside and is2Inside bits occordingly
                //(use a bitfield)
                //2.go through the Sutheran-Hodgman cases for the 2 vertices
                //For finding the intersection use Bresenham line algo.
                //until line is on the screenRect boundary, then add that point
                //to outPoly
                
                int clippedPolyLen = 0;
                SDL_Rect screenRect = {0, 0, pixelBuffer.width, pixelBuffer.height};
                for (int j = 0; j < 3; j++)
                {
                    uint8_t insideBitField = 0;
                    if (isInsideRect(displayPoly.vectors[j], screenRect))
                        insideBitField |= 0x1;
                    if (isInsideRect(displayPoly.vectors[(j + 1) % 3], screenRect))
                        insideBitField |= 0x2;

                    //Switch over the Sutherland-Hodgman cases
                    switch (insideBitField)
                    {
                        case 0:
                        //Both outside
                        break;
                        case 1:
                        //First vector inside
                        clippedPolygon[clippedPolyLen++] = displayPoly.vectors[j];
                        clippedPolygon[clippedPolyLen++] = getIntersect(displayPoly.vectors[j], displayPoly.vectors[(j + 1) % 3], screenRect);
                        break;
                        case 2:
                        //Second vector inside
                        clippedPolygon[clippedPolyLen++] = getIntersect(displayPoly.vectors[j], displayPoly.vectors[(j + 1) % 3], screenRect);
                        break;
                        case 3:
                        //Both inside
                        clippedPolygon[clippedPolyLen++] = displayPoly.vectors[j];
                        break;
                    }
                }
                clippedPolyTriangleCount = clippedPolyLen - 2;
            }

            Triangle clippedTriangles[3];
            for (int j = 0; j < clippedPolyTriangleCount; j++)
            {
                //2.convert pentagon(at most) to 3 triangles
                //Returns Triangle[3] clippedPolygons
                clippedTriangles[j].vectors[0] = clippedPolygon[0];
                clippedTriangles[j].vectors[1] = clippedPolygon[1 + j];
                clippedTriangles[j].vectors[2] = clippedPolygon[2 + j];

                if(shouldDrawSurfaces)
                {
                    rasterizePolygon(clippedTriangles[j], fillColor, pixelBuffer);
                }
                fillColor = ~fillColor;
                //Only draw lines between vectors that haven't been culled
                if(shouldDrawWireframe)
                {
                    for (int k = 0; k < 3; k++)
                        drawLine(clippedTriangles[j].vectors[k], clippedTriangles[j].vectors[(k + 1) % 3], lineColor, pixelBuffer);        
                }
            }
        }
    }
}

Mesh loadMeshFromFile(char* fileName)
{
    FILE* file = fopen(fileName, "r");
    Mesh mesh = {0};
    int lineCount = 0;

    if (file == NULL)
    {
        SDL_Log("Could not open mesh file.");
        SDL_Log(fileName);
    }

    //Save the pos beginning of file
    fpos_t filePos;
    fgetpos(file, &filePos);
    //Count the number of lines in the file
    {
        int ch;
        while (EOF != (ch=getc(file)))
               if (ch=='\n')
                ++lineCount;
    }

    mesh.polyCount = lineCount;
    mesh.polygons = (Triangle*)malloc(lineCount * sizeof(Triangle));

    //Go back to beginning of file
    fsetpos(file, &filePos);
    char line[256] = {0};
    int i = 0;
    int k = 0;
    if (fgets(line, 255, file) == NULL) SDL_Log("fgets error in loadMeshFromFile()!");
    while (!feof(file))
    {
        float vertices[9] = {0};
        //Split line by spaces and store floats
        k = 0;
        while (k < 9)
        {
            if (k == 0)
                vertices[k] = atof(strtok(line, " "));
            else
                vertices[k] = atof(strtok(NULL, " "));
            k++;
        }
        //Store the loaded vertices into the return polygon
        mesh.polygons[i].vectors[0].x = vertices[0];
        mesh.polygons[i].vectors[0].y = vertices[1];
        mesh.polygons[i].vectors[0].z = vertices[2];
        mesh.polygons[i].vectors[1].x = vertices[3];
        mesh.polygons[i].vectors[1].y = vertices[4];
        mesh.polygons[i].vectors[1].z = vertices[5];
        mesh.polygons[i].vectors[2].x = vertices[6];
        mesh.polygons[i].vectors[2].y = vertices[7];
        mesh.polygons[i].vectors[2].z = vertices[8];
        if (fgets(line, 255, file) == NULL) SDL_Log("fgets error in loadMeshFromFile()!");
        i++;
    }
    return mesh;
}

