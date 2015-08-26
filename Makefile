TARGETS = sysjitter
include rules.mk


is_ppc		:= $(shell (uname -m || uname -p) | grep ppc)
is_x86		:= $(shell (uname -m || uname -p) | grep i.86)
is_x86_64	:= $(shell (uname -m || uname -p) | grep x86_64)

ifneq ($(is_x86),)
# Need to tell gcc we have a reasonably recent cpu to get the atomics.
CFLAGS += -march=i686
endif


sysjitter: LIBS := -lpthread
