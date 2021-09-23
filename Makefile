SHELL := /bin/bash


# ===============================
# COMMANDS

CFLAGS = -O3
CC = gcc
RM = rm -f

# ===============================
# TARGETS

EXECUTABLES = serial blocked

default:all

all: $(EXECUTABLES)

serial: test_serial.c
	$(CC) $(CFLAGS) -o test_serial test_serial.c -lm

blocked: test_blocked.c
	$(CC) $(CFLAGS) -o test_blocked test_blocked.c -lm

.PHONY: clean

# =================================
# CLEAN

clean:
	$(RM) *.o *~ $(EXECUTABLES)
