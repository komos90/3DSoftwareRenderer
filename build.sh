#/bin/bash

gcc -std=c99 -D_BSD_SOURCE -Wall src/*.c -lSDL2 -lm -o bin/game.elf

