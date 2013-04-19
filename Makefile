CC = g++
FLAGS = -Wall -O2
LIBS = -lSDL2 -lGLEW -lopencv_core -lopencv_highgui
FRAMEWORKS = -framework OpenGL
TARGET = SPH
OBJ = OpenGL.o Solver.o Renderer.o Timer.o VideoWriter.o

$(TARGET): main.cpp $(OBJ)
	$(CC) -o $(TARGET) main.cpp $(OBJ) $(FLAGS) $(LIBS) $(FRAMEWORKS) 
%.o: %.cpp
	$(CC) -c $< $(FLAGS)

.PHONY: clean
clean:
	rm $(OBJ) $(TARGET)
