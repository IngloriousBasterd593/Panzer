@echo off
mkdir "build/win" 2>NUL
gcc main.c glad.c graphics.c -IC:\Panzer\sdl2\include\SDL2 -IC:\Panzer\glad\include -IC:\Panzer\glew-2.1.0\include -IC:\Panzer\glfw-3.4.bin.WIN64\include -LC:\Panzer\sdl2\lib -LC:\Panzer\glew-2.1.0\lib\Release\x64 -LC:\Panzer\glfw-3.4.bin.WIN64\lib-mingw-w64 -g -ggdb -lmingw32 -lSDL2main -lSDL2 -lSDL2_image -lSDL2_mixer -lSDL2_ttf -lglew32 -lglfw3 -lopengl32 -lgdi32 -o build/win/main
cd "build/win"
.\main


