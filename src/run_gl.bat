@echo off
mkdir "glbuild/glwin" 2>NUL
gcc glmain.c -I C:\glfw-3.3.8.bin.WIN64\include -I C:\glew-2.1.0\include -L C:\glfw-3.3.8.bin.WIN64\lib-mingw-w64 -L C:\glew-2.1.0\lib\Release\x64 -lglfw3dll -lglew32 -lopengl32 -o  glbuild/glwin/glmain  
cd "glbuild/glwin"
.\glmain