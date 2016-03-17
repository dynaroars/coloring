#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <algorithm>
using namespace std;

#include "DimacsBinary.h" //reading DIMACS format
#include "settings.h"  //gloval variables etc

/************* GLOBAL ***************/
Vertex **pVertices;
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

  time_t seed_t = (argc == 3) ? atoi(argv[2]) : time(0);
  srand(seed_t);
  
  //initialize input graphs and other settings
  const char *inputFile = argv[1];  
  auto nVertices = 0, nEdges = 0;
  readDIMACSBinaryFormat(inputFile, nVertices, nEdges);
  if(BB){
    printf("Graph %s, nVertices %d, nEdges %d\n",inputFile,nVertices,nEdges);
  }
  pVertices = initVerticesAndEdges(nVertices, nEdges);  
  for(auto i = 0; i < nVertices; ++i) bestColorResult.push_back(i);

  //Recursive Largest First with Xconstraint
  int XRLFColor[nVertices];
  auto bestResult = XRLF(XRLFColor, nVertices, nEdges);

  int tempColorAssigned[nVertices];
  for(int i=0; i < nVertices ;++i){tempColorAssigned[i]=-1;}

  setUpColorClasses(tempColorAssigned, XRLFColor, bestResult, nVertices);

  //Ant Algorithm
  auto bestCycle = antsOps(bestResult, tempColorAssigned, nVertices);
  printSol(BB==true, seed_t, bestResult, bestCycle, nVertices, nEdges);
  
  //cleanup
  for (auto i = 0; i < nVertices; ++i){delete pVertices[i];} delete pVertices;
  return 0; 
}