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
auto BB = false ;//debug option
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
    printf("Graph %s, nVertices %d, nEdges %d\n",
	   inputFile, nVertices, nEdges);
  }
  Vertex **pVertices = initVerticesAndEdges(nVertices, nEdges);
  Sol sol;
  XRLF(sol, pVertices, nVertices); //Recursive Largest First with Xconstraint
  vector<int> ttmp_colors = setUpColorClasses(sol.colors, sol.nColors);
  antsOps(sol, ttmp_colors, pVertices, nVertices);
  printSol(sol, seed_t, nEdges, BB==true);  
  
  //cleanup
  for(auto i = 0; i < nVertices; ++i){delete pVertices[i];} delete pVertices;
  return 0; 
}
