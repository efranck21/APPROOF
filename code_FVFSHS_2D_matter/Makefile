APPLI1 = run

CSRC1 = $(shell echo Mesh/*.cpp)
CSRC2 =  $(shell echo AnaSol/*.cpp)
CSRC3 = $(shell echo BC/*.cpp)
CSRC4 =  $(shell echo Diagnostics/*.cpp)
CSRC5 = $(shell echo Init/*.cpp)
CSRC6 =  $(shell echo PlotDatas/*.cpp)
CSRC7 = $(shell echo Flux/*.cpp)
CSRC8 =  $(shell echo Source/*.cpp)
CSRC9 =  $(shell echo Functions/*.cpp)
CSRC10 =  $(shell echo Timeschemes/*.cpp)
CSRC11 =  $(shell echo ParamPhysic/*.cpp)
CSRC12 =  $(shell echo HighOrder/*.cpp)
CSRC0 =  $(shell echo *.cpp)
CC = g++ -g -I./hppfiles -I./hppParam -I.


LIB = -lm 

INCLUDES = 
WARNINGS = -Wall -Wsign-compare -mmacosx-version-min=10.7

 OPT =  -O2 -DNDEBUG $(INCLUDES) $(WARNINGS)

all:	$(APPLI1)
	@echo OK.

.cpp.o:
	$(CC) $(OPT) $(UI) -c $*.cpp -Wall -o $*.o
COBJ0 = $(CSRC0:.cpp=.o)
COBJ1 = $(CSRC1:.cpp=.o)
COBJ2 = $(CSRC2:.cpp=.o)
COBJ3 = $(CSRC3:.cpp=.o)
COBJ4 = $(CSRC4:.cpp=.o)
COBJ5 = $(CSRC5:.cpp=.o)
COBJ6 = $(CSRC6:.cpp=.o)
COBJ7 = $(CSRC7:.cpp=.o)
COBJ8 = $(CSRC8:.cpp=.o)
COBJ9 = $(CSRC9:.cpp=.o)
COBJ10 = $(CSRC10:.cpp=.o)
COBJ11 = $(CSRC11:.cpp=.o)
COBJ12 = $(CSRC12:.cpp=.o)
COBJ0 = $(CSRC0:.cpp=.o)

$(APPLI1):$(COBJ1) $(COBJ2) $(COBJ3) $(COBJ4) $(COBJ5) $(COBJ6) $(COBJ7) $(COBJ8) $(COBJ9) $(COBJ10) $(COBJ11) $(COBJ12) $(COBJ0)
	$(CC) $(OPT) $(UI) -o $(APPLI1) $^ $(LIB)

clean:
	@echo Cleaning up...
	@rm -f _tout.c nohup.out core gmon.out
	@rm -f Mesh/*.o BC/*.o Init/*.o AnaSol/*.o Diagnostics/*.o ParamPhysic/*.o
	@rm -f PlotDatas/*.o Flux/*.o Source/*.o
	@rm -f Functions/*.o Timeschemes/*.o HighOrder/*.o
	@rm -f *.o *~ $(APPLI1)
	@echo Done.

cleandata:
	@rm -rf DATA 
	@mkdir DATA
