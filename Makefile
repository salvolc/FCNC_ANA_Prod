#CopyPaste Sachen damit hoffentlich root einigermaßen zufrieden läuft
#--------------------------------------------------------------------
# Root variables (SYSTEM)
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := -lMinuit $(shell root-config --libs)
ROOTTEST	 := -I$(shell root-config --prefix)/test
EXROOT = -L/home/salv/Dokumente/Masterarbeit/MadGraph/MG5_aMC_v2_6_0/ExRootAnalysis/ -lExRootAnalysis

DELPHES = -L/home/salv/Dokumente/Masterarbeit/MadGraph/Delphes-3.4.1/ -lDelphes

# Root variables (.local) 
# ROOTCFLAGS :=
# ROOTLIBS := 

CXXFLAGS     =  -g -Wall -fPIC -Wno-deprecated -ggdb
CXXFLAGS    += $(ROOTCFLAGS) $(ROOTTEST) 
CXXFLAGS    += -I. -I./include -I/home/salv/Dokumente/Masterarbeit/MadGraph/MG5_aMC_v2_6_0/ExRootAnalysis/ExRootAnalysis -I/home/salv/Dokumente/Masterarbeit/MadGraph/Delphes-3.4.1/
LIBS 		= $(ROOTLIBS) -L./ -lTreePlayer -lMinuit $(EXROOT) $(DELPHES)
#--------------------------------------------------------------------
plots = plots/photon/tua_RH_Photon_PT.pdf 
pres = Presentation.pdf
SRC = $(wildcard *.cxx)
#$(info $$plots is [${plots}])
OOBJ = $(SRC:%.cxx=%.o)
#$(info $$OOBJ is [${OOBJ}])

EX = $(wildcard *.cpp)
PRO = $(EX:%.cpp=%)
#$(info $$PRO is [${PRO}])

all: $(PRO) $(OOBJ) $(plots) $(pres)

%.o: %.cxx %.h
	g++ -c $< -o $@ $(CXXFLAGS)

$(PRO): *.cpp $(OOBJ)
	g++ $^ -o $@ -fopenmp $(CXXFLAGS) $(LIBS)

$(plots): $(PRO) aufgabe.py
	./$(PRO)
	python3.5 aufgabe.py

py: aufgabe.py
	python3.5 aufgabe.py

light: $(PRO)
	./$(PRO)

pres: $(plots)
	cd Vortrag/ && make
	cp Vortrag/build/Presentation.pdf ./

#aufgabe: aufgabe.cpp

#$(data): aufgabe
#	./aufgabe

#$(plots): $(data) aufgabe.py
#	python3.5 aufgabe.py


clean:
	@rm -rf $(OOBJ) $(PRO) 

