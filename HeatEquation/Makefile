CC = cc
SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)
BIN = ./bin
TARGET = calc.out

DFLAGS = -lm 
CFLAGS = -I./include -g -std=c99 -Wall -O3

all: dir $(BIN)/$(TARGET)

dir: ${BIN}

${BIN}:
	mkdir -p $(BIN)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(BIN)/$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^

.PHONY: clean
clean:
	rm -f $(OBJ) $(BIN)/$(TARGET)
