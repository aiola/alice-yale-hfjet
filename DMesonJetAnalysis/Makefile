CC=`root-config --cxx`
LD=`root-config --ld`
CFLAGS=-c -g -Wall `root-config --cflags` -fPIC
LDFLAGS=-shared `root-config --libs` -lEG 
SOURCES=MassFitter.cxx
OBJECTS=$(SOURCES:.cxx=.o)
LIBRARY=MassFitter.so

$(LIBRARY): $(OBJECTS)
	$(LD) $(OBJECTS) $(OBJECTS:.o=_Dict.o) -o $@ $(LDFLAGS)

%.o: %.cxx
	rootcint -f $(@:.o=_Dict.cxx) -c -I`root-config --incdir` $(@:.o=.h)
	$(CC) $(CFLAGS) $(@:.o=_Dict.cxx) $(@:.o=.cxx)
	
clean:
	rm ./*~ ./*.o ./MassFitter.so
