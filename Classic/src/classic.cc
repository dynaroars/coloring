#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cassert>
using namespace std;

#include "DimacsBinary.h" //reading DIMACS format
#include "settings.h"  //gloval variables etc

/************* GLOBAL ***************/
time_t seed_t = 0;
vertex **pVertices;
vector<ant *> vAnts;

auto BB = false ;//debug option

//results stores here
vector<int> bestColorResult;

#include "classic.h"

int main(int argc, char *argv[]){
  if (argc < 2) {
    fprintf(stderr,
	    "ERR: Usage is %s input_file.b [optional: seed]\n",
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

  //seed
  if(argc == 3){
    seed_t = atoi(argv[2]);
  }else {
    seed_t = time(0);}
  srand(seed_t);
  
  //initialize input graphs and other settings
  auto nVertices = 0;
  auto nEdges = 0;

  char const * inputFile = argv[1]; 
  auto nAnts = 0;
  auto nCycles=0;
  auto breakCycles=0;  
  auto rSizeLimit=0;
  auto moveLimit=0;
  auto nRLFSetLimit=0;
  int bestResult = initialize(inputFile,
			      nVertices, nEdges, nAnts,
			      nCycles, breakCycles, moveLimit,
			      rSizeLimit, nRLFSetLimit);
  
  //Recursive Largest First with Xconstraint
  int XRLFColor[nVertices];
  XRLF(XRLFColor, bestResult, nVertices, nEdges, nRLFSetLimit);

  int tempColorAssigned[nVertices];
  for(int i=0; i <nVertices ;++i){tempColorAssigned[i]=-1;}

  setUpColorClasses(tempColorAssigned, XRLFColor, bestResult, nVertices);

  //Ant Algorithm
  auto bestCycle = antsOps(bestResult, tempColorAssigned,
			   nVertices, nAnts,
			   nCycles, breakCycles, moveLimit, rSizeLimit);

  printSol(BB==true, bestResult, bestCycle, nVertices, nEdges);
  
  cleanUp(nVertices);
  return 0; 
}
