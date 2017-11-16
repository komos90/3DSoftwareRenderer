#/bin/bash
clang -std=c99 -Wall -O3 -g src/*.c -lSDL2 -lm -o bin/game.elf
