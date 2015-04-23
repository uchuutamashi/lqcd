IDIR =/usr/local/include
CC=gcc
CFLAGS=-I$(IDIR) -std=c99

ODIR=obj
LDIR =/usr/local/lib

LIBS=-lm -lgsl -lgslcblas

_DEPS =
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = su2_gauge_hb.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

su2_gauge_hb: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
