CC = g++
FLAGS = -Wall -O2
LIBS = -lSDL2 -lGLEW 
FRAMEWORKS = -framework OpenGL
TARGET = SPH
OBJ = OpenGL.o Solver.o Renderer.o Timer.o

.PHONY:
$(TARGET): main.cpp $(OBJ)
	$(CC) -o $(TARGET) main.cpp $(OBJ) $(FLAGS) $(LIBS) $(FRAMEWORKS) 

%.o: %.cpp
	$(CC) -c $< $(FLAGS)
