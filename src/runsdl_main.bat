@echo off
mkdir "build/win" 2>NUL
gcc main.c graphics.c -IC:\sdl2\include\SDL2 -LC:\sdl2\lib -g -ggdb -lmingw32 -lSDL2main -lSDL2 -lSDL2_image -lSDL2_mixer -lSDL2_ttf -o build/win/main
cd "build/win"
.\main

