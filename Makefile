
FILES =	Makefile  .cvsignore *-plot.sh jobs.sh

default: main

examples.zip: ${FILES}
	@zip -9 -u -y $@ $^

clean:
	-/bin/rm -f *.o *.tmp plot.eps PI[0-9]* \
		*.G.* *.A.* *.Aalpha.* *.Abeta.* \
		*.Gx.* *.Gy.* *.Gz.* *.x.* *.y.* *.z.* *.b.* *.x0.* *~

distclean: clean
	-/bin/rm -f team21 *.o *~

lib:
	@(cd ../src; $(MAKE))

PLASMA_OFILES = utils2D.o medit2D.o vtk2D.o  equil2D.o



main: $(PLASMA_OFILES) 
	$(LINKER) $(LDFLAGS) -o $@ $(PLASMA_OFILES) $(LIBS)


include ${PHG_MAKEFILE_INC}

CFLAGS +=-Wno-unused -Wno-implicit
#CFLAGS +=-D _HWENO


.PHONY: default all clean distclean lib
