#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cassert>
using namespace std;

#include "settings.h"  //gloval variables etc
#include "DimacsBinary.h" //reading DIMACS format
#include "classic.h"

int main(int argc, char *argv[]){
  if (argc < 2) {
    fprintf(stderr,
	    "ERR: Usage is %s Dimacs *Binary* Graph [optional: seed]\n",
	    argv[0]);
    exit(10);
  }

#ifdef DB
#warning "Debug Option On"
  BB = true;
#endif

#ifdef V
#warning "VERBOSE ON"
#endif
  
  //initialize input graphs and other settings
  auto nVertices = 0;
  auto nEdges = 0;
  auto nAnts = 0;
  auto nCycles=0;
  auto breakCycles=0;  
  auto rSizeLimit=0;
  auto moveLimit=0;
  auto nRLFSetLimit=0;
  initialize(argc, argv, nVertices, nEdges, nAnts,
	     bestResult,
	     nCycles, breakCycles, moveLimit,
	     rSizeLimit, nRLFSetLimit);
  
  //Recursive Largest First with Xconstraint
  int XRLFColor[nVertices];
  int nColorsXRLF = XRLF(XRLFColor, nVertices, nEdges, nRLFSetLimit);

  int tempColorAssigned[nVertices];
  for(int i=0;i <nVertices ;++i){tempColorAssigned[i]=-1;}

  setUpColorClasses(tempColorAssigned, XRLFColor, nColorsXRLF, nVertices);

  //Ant Algorithm
  auto bestCycle = antsOps(nColorsXRLF, tempColorAssigned,
			   nVertices, nAnts,
			   nCycles, breakCycles, moveLimit, rSizeLimit);

  printSol(BB==true, bestResult, bestCycle, nVertices, nEdges);
  
  cleanUp(nVertices);
  return 0; 
}
