# Makefile
TARGET    = d2q9-bgk.llvm
BUILD_DIR = ./build
SRC_DIRS  = ./src

SRCS := $(shell find $(SRC_DIRS) -name *.cpp -or -name *.c -or -name *.s)
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS  := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CXX       = clang
CPPFLAGS  = $(INC_FLAGS) -MMD -MP
CPPFLAGS += -DNDEBUG
CPPFLAGS += -std=c++2b -Wall -g
CPPFLAGS += -O3 -march=native
CPPFLAGS += --gcc-toolchain=/apps/2021/gcc/10.2/
LDFLAGS   = -lm

FINAL_STATE_FILE      = ./final_state.dat
AV_VELS_FILE          = ./av_vels.dat
REF_FINAL_STATE_FILE  = check/128x128.final_state.dat
REF_AV_VELS_FILE      = check/128x128.av_vels.dat

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CPPFLAGS) $(OBJS) -o $@ $(LDFLAGS)

# c++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

check: $(TARGET)
	python check/check.py --ref-av-vels-file=$(REF_AV_VELS_FILE) --ref-final-state-file=$(REF_FINAL_STATE_FILE) --av-vels-file=$(AV_VELS_FILE) --final-state-file=$(FINAL_STATE_FILE)

.PHONY: all check clean

clean:
	$(RM) -f $(TARGET)
	$(RM) -r $(BUILD_DIR)

-include $(DEPS)

MKDIR_P ?= mkdir -p