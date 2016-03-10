CC=clang++
LD=clang++
CFLAGS=-c -g -Wall `root-config --cflags`
LDFLAGS=`root-config --glibs` -leg -shared
SOURCES=DJetCorrBase.cxx DJetCorrAnalysis.cxx DJetCorrAnalysisParams.cxx MassFitter.cxx HistoStyler.cxx
OBJECTS=$(SOURCES:.cxx=.o)
LIBRARY=DJetCorr.so

$(LIBRARY): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(OBJECTS:.o=_Dict.o) -o $@

%.o:
	rootcint -f $(@:.o=_Dict.cxx) -c -I`root-config --incdir` $(@:.o=.h)
	$(CC) $(CFLAGS) $(@:.o=_Dict.cxx) $(@:.o=.cxx)
	
clean:
	rm ./*~ ./*.o ./DJetCorr.so
   