# Makefile
TARGET    = d2q9-bgk.llvm
BUILD_DIR = ./build
SRC_DIRS  = ./src

SRCS := $(shell find $(SRC_DIRS) -name *.cpp -or -name *.c -or -name *.s)
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS  := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS  = $(INC_FLAGS) -MMD -MP
CPPFLAGS += -DNDEBUG
CPPFLAGS += -std=c++20 -Wall -g

# Clang setup
# CXX       = clang
# CPPFLAGS += -O3 -march=native -mavx512f -mavx512f -mavx512dq -mavx512bw
# CPPFLAGS += --gcc-toolchain=/apps/2021/gcc/10.2/

# Intel setup
CXX       = icpx
CPPFLAGS += -Ofast -march=native -fno-strict-aliasing
CPPFLAGS += -mavx512f -mavx512dq -mavx512bw
CPPFLAGS += -mprefer-vector-width=512
#CPPFLAGS += -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize

LDFLAGS   = -lm

# MPI vars
MPI_CXX_FLAGS = -I/opt/mpi/openmpi/4.1.4.2/include -I/opt/mpi/openmpi/4.1.4.2/fortran-intel-ilp64/include -g -m64
MPI_LD_FLAGS  = -g -m64 -L/usr/lib64 -L/opt/mpi/openmpi/4.1.4.2/lib -L/opt/mpi/openmpi/4.1.4.2/fortran-intel-ilp64/lib -lmpi

FINAL_STATE_FILE      = ./final_state.dat
AV_VELS_FILE          = ./av_vels.dat
REF_FINAL_STATE_FILE  = check/128x128.final_state.dat
REF_AV_VELS_FILE      = check/128x128.av_vels.dat

all: $(TARGET)

run: $(TARGET)
	@./run

$(TARGET): $(OBJS)
	$(CXX) $(CPPFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(MPI_LD_FLAGS)

# c++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@ $(MPI_CXX_FLAGS)

check:
	python3 check/check.py --ref-av-vels-file=$(REF_AV_VELS_FILE) --ref-final-state-file=$(REF_FINAL_STATE_FILE) --av-vels-file=$(AV_VELS_FILE) --final-state-file=$(FINAL_STATE_FILE)

.PHONY: all check clean

clean:
	$(RM) -f $(TARGET)
	$(RM) -r $(BUILD_DIR)

-include $(DEPS)

MKDIR_P ?= mkdir -p