CC = gcc
NVCC =
TARGET = test

CFLAGS = -I../include -g
LIBS = -L/../lib -lcheck -lsubunit -lgsl -lgslcblas -lm

# Source code from Python, which is what we'll be testing
PYTHON_SOURCES = \
	../source/xlog.c \
	../source/rdpar.c \
	../source/synonyms.c \
	../source/matrix_cpu.c

# This is sources to be compiled with NVCC
CUDA_SOURCES = ../source/matrix_gpu.cu

# Source code for the testing framework and unit test files
TARGET_SOURCES = source/main.c
TEST_SOURCES = source/tests/test_matrix.c

# This is all sources, to be compiled with CC
SOURCES = $(TARGET_SOURCES) $(TEST_SOURCES) $(PYTHON_SOURCES)

# Object files, again keeping C and CUDA sources separate
OBJS = $(patsubst %.c, %.o, $(SOURCES))

ifneq ($(NVCC),)
	CFLAGS += -DCUDA_ON
	LIBS += -lcudart -lcusolver
	CUDA_OBJS = $(patsubst %.cu, %.o, $(CUDA_SOURCES))
else
	CUDA_OBJS =
endif

all: startup clean $(TARGET)

startup:
	@echo 'ALL SOURCES =' $(SOURCES) $(CUDA_SOURCES)

cuda: startup $(CUDA_SOURCES)
ifneq ($(NVCC),)
	$(NVCC) -c $(CUDA_SOURCES) -DCUDA_ON -o $(CUDA_OBJS)
endif

# TODO: we should store objects in a directory local to test
$(TARGET): startup cuda $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(CUDA_OBJS) -o $(TARGET) $(LIBS)

clean:
	rm -f $(TARGET) $(OBJS) $(CUDA_OBJS)

.PHONY: clean
