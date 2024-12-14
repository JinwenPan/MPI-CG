CC = mpic++
CFLAGS  = -O3 -std=c++17 -Wall -Wshadow -DNDEBUG -Werror
TARGET = cgsolver

all: main.cpp 
	$(CC) $(CFLAGS) main.cpp Grid.cpp -o $(TARGET)