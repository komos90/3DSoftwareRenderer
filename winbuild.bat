@echo off
@rem gcc -std=c99 -O3 src/*.c -IC:/MinGW/SDL/include/SDL2 -LC:/MinGW/SDL/lib -Wall -lmingw32 -lSDL2main -lSDL2 -o bin/wingame
clang -m32 -Ilibs/SDL2-2.0.7/include -Llibs/SDL2-2.0.7/lib/x86 -lSDL2main -lSDL2 -Xlinker /subsystem:windows src/*.c -o bin/wingame.exe

pause
