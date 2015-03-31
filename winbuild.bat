@echo off
gcc -std=c99 -O3 src/*.c -IC:/MinGW/SDL/include/SDL2 -LC:/MinGW/SDL/lib -Wall -lmingw32 -lSDL2main -lSDL2 -o bin/wingame
pause
