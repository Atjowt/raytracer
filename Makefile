.PHONY: all run clean

NAME := raytracer
CC := gcc
CFLAGS := -std=c99 -pedantic -Wall -O2 -lm
SOURCES := raytracer.c

all: compile_commands.json

compile_commands.json: $(NAME)
	bear --output $@ -- make $<

$(NAME): $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $^

run: $(NAME)
	./$<
	feh ./output.ppm

clean:
	rm -f ./compile_commands.json ./$(NAME)
