TARGETS = sysjitter
include rules.mk


# Only supports x86 and x86_64 processors at the moment.
processor := $(shell uname -m || uname -p)
ifneq ($(processor),x86_64)
# Need to tell gcc we have a reasonably recent cpu to get the atomics.
CFLAGS += -march=i686
endif


sysjitter: LIBS := -lpthread
