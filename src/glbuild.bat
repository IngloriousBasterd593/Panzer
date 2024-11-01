@echo off
mkdir "glbuild" 2>NUL
gcc glsusy.c glad.c -I..\glad\include -I..\glew-2.1.0\include -I..\glfw-3.4.bin.WIN64\include -L..\glew-2.1.0\lib\Release\x64 -L..\glfw-3.4.bin.WIN64\lib-mingw-w64 -g -ggdb -lmingw32  -lglew32 -lglfw3 -lopengl32 -lgdi32 -o build/main 
cd glbuild
.\glsusy