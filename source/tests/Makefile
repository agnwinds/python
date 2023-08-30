CC = mpicc
NVCC =
TARGET = py_unit_test

C_FLAGS = -I../../include
NVCC_FLAGS = -DCUDA_ON

LIBS = -L../../lib -lcunit -lgsl -lgslcblas -lm -lstdc++

# Source code from Python, which is what we'll be testing
#PYTHON_SOURCES = \
#	../xlog.c \
#	../rdpar.c \
#	../synonyms.c \
#	../matrix_cpu.c \

# This is less than ideal, but it's far too complicated and time consuming to
# figure out what tests needs which source files. So as we've done in the past,
# we'll just include ALL of Python's source
PYTHON_SOURCES = \
	agn.c anisowind.c atomic_extern_init.c atomicdata.c atomicdata_init.c  \
	atomicdata_sub.c bands.c bb.c bilinear.c brem.c cdf.c charge_exchange.c  \
	compton.c continuum.c cooling.c corona.c cv.c cylind_var.c cylindrical.c  \
	density.c diag.c dielectronic.c direct_ion.c disk.c disk_init.c  \
	disk_photon_gen.c emission.c estimators_macro.c estimators_simple.c  \
	extract.c frame.c gradv.c gridwind.c homologous.c hydro_import.c import.c  \
	import_calloc.c import_cylindrical.c import_rtheta.c import_spherical.c  \
	ionization.c knigge.c levels.c lines.c macro_accelerate.c macro_gen_f.c  \
	macro_gov.c matom.c matom_diag.c matrix_cpu.c matrix_ion.c  \
	models_extern_init.c para_update.c parse.c partition.c paths.c  \
	phot_util.c photon2d.c photon_gen.c photon_gen_matom.c pi_rates.c  \
	python_extern_init.c radiation.c random.c rdpar.c recipes.c recomb.c  \
	resonate.c reverb.c roche.c rtheta.c run.c saha.c setup.c setup_disk.c  \
	setup_domains.c setup_files.c setup_line_transfer.c setup_reverb.c  \
	setup_star_bh.c shell_wind.c signal.c spectra.c spectral_estimators.c  \
	spherical.c stellar_wind.c sv.c synonyms.c time.c trans_phot.c  \
	vvector.c walls.c wind.c wind2d.c wind_sum.c wind_updates2d.c wind_util.c  \
	windsave.c windsave2table_sub.c xlog.c xtest.c zeta.c get_models.c rdpar_init.c

PYTHON_SOURCES := $(patsubst %,../%, $(PYTHON_SOURCES))

# This is sources to be compiled with NVCC
CUDA_SOURCES = ../matrix_gpu.cu

# Source code for the testing framework
TARGET_SOURCES = main.c

# Source code for unit tests
TEST_SOURCES = \
	tests/test_matrix.c \
	tests/test_compton.c

# This is all sources, to be compiled with CC
SOURCES = $(TARGET_SOURCES) $(TEST_SOURCES) $(PYTHON_SOURCES)

# Object files, again keeping C and CUDA sources separate
OBJS = $(patsubst %.c, %.o, $(SOURCES))

ifeq ($(CC), mpicc)
	C_FLAGS += -DMPI_ON
endif

ifneq ($(NVCC),)
	CUDA_OBJS = $(patsubst %.cu, %.o, $(CUDA_SOURCES))
	LIBS += -lcudart -lcusolver
	C_FLAGS += -DCUDA_ON
	OBJS += $(CUDA_OBJS)
endif

all: startup clean $(TARGET)

startup:
	@echo "C_FLAGS = " $(C_FLAGS)
	@echo "NVCC_FLAGS =" $(NVCC_FLAGS)
	@echo "LIBS =" $(LIBS)

cuda: $(CUDA_SOURCES)
ifneq ($(NVCC),)
	$(NVCC) -c $(NVCC_FLAGS) $(CUDA_SOURCES) -o $(CUDA_OBJS)
endif

$(TARGET): cuda $(OBJS)
	$(CC) $(C_FLAGS) $(OBJS) -o $(TARGET) $(LIBS)
	mv $(TARGET) ../../bin

check: startup $(TARGET)
	../../bin/$(TARGET)

clean:
	rm -f *.o
	rm -f tests/*.o
	rm -f ../*.o

.PHONY: clean
