# Build the C implementation of Hornet.
#
#   make                                   # htslib from pkg-config (or -lhts)
#   make HTSLIB_PREFIX=$CONDA_PREFIX       # htslib installed in a prefix (e.g. conda)
#   make HTSLIB_DIR=/path/to/htslib-1.x    # statically link a built htslib source tree
#
# Then copy ./hornet anywhere on your PATH (or `make install PREFIX=...`).

CC      ?= cc
CFLAGS  ?= -O3 -g -Wall
PREFIX  ?= /usr/local

ifdef HTSLIB_DIR
  HTS_CPPFLAGS = -I$(HTSLIB_DIR)
  HTS_LIBS     = $(HTSLIB_DIR)/libhts.a $(shell sed -n 's/^static_libs=//p' $(HTSLIB_DIR)/htslib.pc.tmp 2>/dev/null || echo -lz -lm -lbz2 -llzma -lcurl)
else ifdef HTSLIB_PREFIX
  HTS_CPPFLAGS = -I$(HTSLIB_PREFIX)/include
  HTS_LIBS     = -L$(HTSLIB_PREFIX)/lib -Wl,-rpath,$(HTSLIB_PREFIX)/lib -lhts
else
  HTS_CPPFLAGS = $(shell pkg-config --cflags htslib 2>/dev/null)
  HTS_LIBS     = $(shell pkg-config --libs htslib 2>/dev/null || echo -lhts)
endif

LIBS = $(HTS_LIBS) -lz -lm -lpthread

hornet: hornet.c
	$(CC) $(CFLAGS) $(HTS_CPPFLAGS) -o $@ $< $(LIBS)

install: hornet
	install -d $(PREFIX)/bin
	install -m 755 hornet $(PREFIX)/bin/hornet

clean:
	rm -f hornet

.PHONY: install clean
