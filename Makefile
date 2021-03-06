TARGET = vetopulseMain 
CC = g++
ROOTCINT = $(ROOTSYS)/bin/rootcint
DICTNAME = vetopulse_dict
SRCS = $(addsuffix .C, $(TARGET))
DIR = .
SDIR = /home/hqian/montecarlo/g4ds10/Linux-g++/SLADLoadClass
SOURCES = $(DIR)/vetopulsemc.cc $(DIR)/vetopulsereal.cc $(DIR)/vetopulsebasic.cc $(DIR)/vetopulsefitfunc.cc $(DIR)/vetopulsetotalfit.cc $(DIR)/odselector/odselector.C $(DIR)/DSTtreeSelector/DSTtreeSelector.C $(DIR)/SLADDSTSelector/SLADDSTSelector.C $(DIR)/vetopulseslad.cc $(SDIR)/SLADLoadClass.cc 
HEADERS = $(DIR)/vetopulsemc.hh $(DIR)/vetopulsereal.hh $(DIR)/vetopulsebasic.hh $(DIR)/vetopulsefitfunc.hh $(DIR)/vetopulsetotalfit.hh $(DIR)/odselector/odselector.h $(DIR)/DSTtreeSelector/DSTtreeSelector.h $(DIR)/SLADDSTSelector/SLADDSTSelector.h $(DIR)/vetopulseslad.hh $(SDIR)/SLADLoadClass.hh
DEPS = $(DIR)/reconmcvar.hh $(DIR)/vetoreadcfg.hh $(SDIR)/sladtpcstruct.hh $(SDIR)/sladodstruct.hh 
OBJS = $(addsuffix .o, $(notdir $(basename $(SRCS))))
SOBJS = $(addsuffix .o, $(notdir $(basename $(SOURCES))))
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags) #-Wall -fPIC -g 
ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs) -lProof -lProofPlayer #-lRooFit -lRooFitCore
ROOTLDFLAGS = $(shell $(ROOTSYS)/bin/root-config --ldflags)

all: $(DICTNAME).C
	$(CC) $(ROOTCFLAGS) $(ROOTLIBS) $(SRCS) $(SOURCES) $^ -o $(TARGET) 

$(DICTNAME).C: $(HEADERS) $(DEPS)
	$(ROOTCINT) -f $@ -c  $^ LinkDef.h 

.PHONY: clean

clean:
	rm -rf $(DIR)/*_dict.h $(DIR)/*_dict.C $(TARGET)