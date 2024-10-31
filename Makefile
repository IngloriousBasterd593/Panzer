CC = gcc

INCLUDE_DIRS = -I..\sdl2\include\SDL2 -I..\glad\include -I..\glew-2.1.0\include -I..\glfw-3.4.bin.WIN64\include
LIB_DIRS = -L..\sdl2\lib -L..\glew-2.1.0\lib\Release\x64 -L..\glfw-3.4.bin.WIN64\lib-mingw-w64

LIBS = -lmingw32 -lSDL2main -lSDL2 -lSDL2_image -lSDL2_mixer -lSDL2_ttf -lglew32 -lglfw3 -lopengl32 -lgdi32

CFLAGS = -g -ggdb

SOURCES = src/main.c src/glad.c

OUTPUT = build\main.exe

all: $(OUTPUT)

$(OUTPUT): $(SOURCES)
	@if not exist build ( mkdir build ) 
	$(CC) $(SOURCES) $(INCLUDE_DIRS) $(LIB_DIRS) $(CFLAGS) -o $(OUTPUT) $(LIBS)

clean:
	del build\main.exe  
	rmdir /s /q build   
