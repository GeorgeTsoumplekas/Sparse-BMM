SHELL := /bin/bash


# ===============================
# COMMANDS

CFLAGS = -O3
CC = gcc
RM = rm -f

# ===============================
# TARGETS

EXECUTABLES = serial blocked parallel parallel_blocked filtered_serial filtered_parallel

default:all

all: $(EXECUTABLES)

serial: test_serial.c
	$(CC) $(CFLAGS) -o test_serial test_serial.c -lm

blocked: test_blocked.c
	$(CC) $(CFLAGS) -o test_blocked test_blocked.c -lm

parallel: test_parallel.c
	$(CC) $(CFLAGS) -o test_parallel test_parallel.c -lm -fopenmp

parallel_blocked: test_blocked_parallel.c
	$(CC) $(CFLAGS) -o test_blocked_parallel test_blocked_parallel.c -lm -fopenmp

filtered_serial: test_serial_filtered.c
	$(CC) $(CFLAGS) -o test_serial_filtered test_serial_filtered.c -lm

filtered_parallel: test_parallel_filtered.c
	$(CC) $(CFLAGS) -o test_parallel_filtered test_parallel_filtered.c -lm -fopenmp

.PHONY: clean

# =================================
# CLEAN

clean:
	$(RM) *.o *~ test_serial test_blocked test_parallel test_blocked_parallel test_serial_filtered test_parallel_filtered
