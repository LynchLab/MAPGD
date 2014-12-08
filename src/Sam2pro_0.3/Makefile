CC=gcc
CFLAGS= -O3 -Wall -Wshadow -pedantic -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -DVER32 \
	-I/opt/local/include/ -L/opt/local/lib/ # -m64 -g

# The source files, object files, libraries and executable name.
SRCFILES= sam2pro.c eprintf.c interface.c tab.c stringUtil.c
OBJFILES= sam2pro.o eprintf.o interface.o tab.o stringUtil.o
LIBS= -lm

EXECFILE= ../../bin/sam2pro
DIRECTORY= Sam2pro

VERSION= 0.3

# The make rule for the executable
.PHONY : all
all : $(EXECFILE)

$(EXECFILE) : $(OBJFILES)
	$(CC) $(CFLAGS) -o $(EXECFILE) $(OBJFILES) $(LIBS)	

doc:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)

# Other Standard make rules
lint : 
	lint $(SRCFILES) | more

clean:
	rm -f *.o ../bin/sam2pro

tarfile:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)
	mkdir $(DIRECTORY)_$(VERSION)
	cp -rf $(SRCFILES) *.h  filter.sh Makefile README COPYING \
	../Doc/sam2proDoc.pdf $(DIRECTORY)_$(VERSION)
	cp ../Data/test.sam $(DIRECTORY)_$(VERSION)
	tar cvzfh $(EXECFILE)_$(VERSION).tgz $(DIRECTORY)_$(VERSION)
	mv $(EXECFILE)_$(VERSION).tgz ../
	/bin/rm -r $(DIRECTORY)_$(VERSION)
