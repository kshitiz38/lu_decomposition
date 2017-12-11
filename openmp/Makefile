CC = g++
FLAGS = -fopenmp -Wall -c -fPIC
COMMON = -I./
TARGET  = luDecompose.so

all: $(TARGET)
$(TARGET): LUDecomp.cpp
	$(CC) $(COMMON) $(FLAGS) -o $(TARGET) -fPIC -shared LUDecomp.cpp
