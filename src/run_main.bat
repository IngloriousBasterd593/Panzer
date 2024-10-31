@echo off
mkdir "build" 2>NUL
gcc main.c glad.c -I..\sdl2\include\SDL2 -I..\glad\include -I..\glew-2.1.0\include -I..\glfw-3.4.bin.WIN64\include -L..\sdl2\lib -L..\glew-2.1.0\lib\Release\x64 -L..\glfw-3.4.bin.WIN64\lib-mingw-w64 -g -ggdb -lmingw32 -lSDL2main -lSDL2 -lSDL2_image -lSDL2_mixer -lSDL2_ttf -lglew32 -lglfw3 -lopengl32 -lgdi32 -o build/main
cd build
.\main


