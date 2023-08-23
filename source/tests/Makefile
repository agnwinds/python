CC = mpicc
NVCC =
TARGET = py_unit_test

CFLAGS = -I../../include -DMPI_ON
NVCCFLAGS = -DCUDA_ON
LIBS = -L../../lib -lcunit -lgsl -lgslcblas -lm

# Source code from Python, which is what we'll be testing
PYTHON_SOURCES = \
	../xlog.c \
	../rdpar.c \
	../synonyms.c \
	../matrix_cpu.c

# This is sources to be compiled with NVCC
CUDA_SOURCES = ../matrix_gpu.cu

# Source code for the testing framework and unit test files
TARGET_SOURCES = source/main.c
TEST_SOURCES = source/tests/test_matrix.c

# This is all sources, to be compiled with CC
SOURCES = $(TARGET_SOURCES) $(TEST_SOURCES) $(PYTHON_SOURCES)

# Object files, again keeping C and CUDA sources separate
OBJS = $(patsubst %.c, %.o, $(SOURCES))

ifneq ($(NVCC),)
	CFLAGS += -DCUDA_ON
	CUDA_LIBS = -lcudart -lcusolver
	CUDA_OBJS = $(patsubst %.cu, %.o, $(CUDA_SOURCES))
else
	CUDA_OBJS =
	CUDA_LIBS =
endif

all: startup clean $(TARGET)

startup:
	@echo 'ALL SOURCES =' $(SOURCES) $(CUDA_SOURCES)

cuda: $(CUDA_SOURCES)
ifneq ($(NVCC),)
	$(NVCC) -c $(NVCCFLAGS) $(CUDA_SOURCES) -o $(CUDA_OBJS)
endif

$(TARGET): cuda $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(CUDA_OBJS) -o $(TARGET) $(CUDA_LIBS) $(LIBS)
	mv $(TARGET) ../../bin

check: $(TARGET)
	../../bin/$(TARGET)

clean:
	rm -f ../../bin/$(TARGET) $(OBJS) $(CUDA_OBJS)

.PHONY: clean
