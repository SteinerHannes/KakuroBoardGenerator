CC = g++
CFLAGS = -std=c++17 -Wall -O3
TARGET = kakuro_generator
SRC = main.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)