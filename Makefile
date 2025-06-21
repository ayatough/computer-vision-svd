CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -g
LDFLAGS = -lm

SRCDIR = .
SOURCES = matrix.c test_svd.c test_8point.c test_homography.c test_pca.c
OBJECTS = $(SOURCES:.c=.o)
TARGET = test_svd
TARGET_8POINT = test_8point
TARGET_HOMOGRAPHY = test_homography
TARGET_PCA = test_pca

.PHONY: all clean run run-8point run-homography run-pca test-all

all: $(TARGET) $(TARGET_8POINT) $(TARGET_HOMOGRAPHY) $(TARGET_PCA)

$(TARGET): matrix.o test_svd.o
	$(CC) matrix.o test_svd.o -o $(TARGET) $(LDFLAGS)

$(TARGET_8POINT): matrix.o test_8point.o
	$(CC) matrix.o test_8point.o -o $(TARGET_8POINT) $(LDFLAGS)

$(TARGET_HOMOGRAPHY): matrix.o test_homography.o
	$(CC) matrix.o test_homography.o -o $(TARGET_HOMOGRAPHY) $(LDFLAGS)

$(TARGET_PCA): matrix.o test_pca.o
	$(CC) matrix.o test_pca.o -o $(TARGET_PCA) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

matrix.o: matrix.c matrix.h

test_svd.o: test_svd.c matrix.h

test_8point.o: test_8point.c matrix.h

test_homography.o: test_homography.c matrix.h

test_pca.o: test_pca.c matrix.h

run: $(TARGET)
	./$(TARGET)

run-8point: $(TARGET_8POINT)
	./$(TARGET_8POINT)

run-homography: $(TARGET_HOMOGRAPHY)
	./$(TARGET_HOMOGRAPHY)

run-pca: $(TARGET_PCA)
	./$(TARGET_PCA)

test-all: $(TARGET) $(TARGET_8POINT) $(TARGET_HOMOGRAPHY) $(TARGET_PCA)
	@echo "Running SVD tests..."
	./$(TARGET)
	@echo "Running 8-point algorithm tests..."
	./$(TARGET_8POINT)
	@echo "Running homography algorithm tests..."
	./$(TARGET_HOMOGRAPHY)
	@echo "Running PCA tests..."
	./$(TARGET_PCA)

clean:
	rm -f $(OBJECTS) $(TARGET) $(TARGET_8POINT) $(TARGET_HOMOGRAPHY) $(TARGET_PCA)

debug: CFLAGS += -DDEBUG -g3
debug: $(TARGET)

install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

.SUFFIXES: .c .o