
ifndef CFLAGS
CFLAGS		= -O2 -Wall
endif

ifdef NDEBUG
CPPFLAGS	+= -DNDEBUG
endif

INCLUDE		+= -I.


%: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(INCLUDE) $< $(CLINK) $(LIBS) -o $@


all: $(TARGETS)

clean:
	rm -f *.o $(TARGETS)
