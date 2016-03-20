#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cassert>
using namespace std;

#include "config.h"  //gloval variables etc
#include "DimacsBinary.h" //reading DIMACS format
#include "alg.h"  //main algorithm

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
  initialize(argc, argv);
  
  //Recursive Largest First with Xconstraint
  int XRLFColor[nVertices];
  int nColorsXRLF = XRLF(XRLFColor, nVertices, nEdges);

  int tempColorAssigned[nVertices];
  for(int i=0;i <nVertices ;++i){tempColorAssigned[i]=-1;}

  setUpColorClasses(nColorsXRLF,tempColorAssigned,XRLFColor);

  //Ant Algorithm
  AntsOps(nColorsXRLF, tempColorAssigned);

  printSol(BB==true);
  
  cleanUp();
  return 0; 
}
