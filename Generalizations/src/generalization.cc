#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <cmath>
using namespace std;
#include "MT_rand.h"
#include "config.h"
#include "commonUtils.h"
#include "DimacsAscii.h"
#include "BMCP2MCPtrans.h"
#include "generalization.h"
#include "Pheromone.h"
#include "JOLT.h"
#include "LocalOpt.h"
#include "AntsOps.h"


void RW_Test();
void Confs_Test();

int main(int argc, char *argv[]){

  if (argc < 3) {fprintf(stderr,"ERR: Usage is %s Dimacs *Ascii* Graph type [optional: seed]\n",argv[0]); exit(10);}

#ifdef DB
#warning "DEBUG OPTION ON"
  BB=true;
#endif

#ifdef V
#warning "VERBOSE ON"
#endif


  inputFile=argv[1];
  PTYPE=atoi(argv[2]);//note GCP (choice 0) uses different program, not this one ... avail choices 1,2,3  
  
  initGraph(argc,argv);//initialize graphs

  int bestConfigA[nVertices]; //results saves here

  int nDistinctColorsIG=IG(bestConfigA,nVertices,nEdges);


  int setupColorConfigA[nVertices];copyArray(bestConfigA,setupColorConfigA);
  setUpColorClasses(nDistinctColorsIG,setupColorConfigA);


  AntsOps(setupColorConfigA,bestConfigA);//run the ants


  printSol(bestConfigA,(argc==4)?TOCHECK:NOCHECK);
  
  cleanUp();//clean up mem 
  return 0; 
}


////////////////////////////////////////////////////////////////////////////////////////////////



// DEBUG AND TESTING CODE ,  NOT IMPORTANT AND NOT USED, ONLY USED WHEN I DEBUG THINGS 
// I JUST LEAVE THEM HERE JUST IN CASE OF REUSING THEM


//   Code to test to see if roulette wheel gives answers proportionally to score
void RW_Test(){
  vector<double>a; 
  for(int i = 0 ; i < 10 ; ++i){
	a.push_back((genrand_int32()%100)/100.00);
  }
  printVector("a ",a);
  
  double ScoreSum=0.0;
  for(int i = 0 ; i < a.size();++i){ScoreSum+= a.at(i) ;}
  
  
  int freqq[a.size()];  for(int i = 0 ; i < a.size() ; ++i){freqq[i]=0;}
  for(int i = 0 ; i < 100 ; ++i){
	freqq[chooseColor_RouletteWheel(a,ScoreSum)]++;	
  }
  
  for(int i = 0 ; i < a.size() ;++i){
	printf("color %d, score %f, totalScore %f (percentile %d), chosen %d times\n",
		   i,a.at(i),ScoreSum,(int)((a.at(i)/ScoreSum)*100.00),freqq[i]);
  }
  
  assert(1==0);
  /////////////////////////////
}




void Confs_Test(){
  vector<vector<int> >aa; 
  for(int i = 0 ; i < 10 ; ++i){
	vector<int>a;
	for(int j = 0 ; j < 5 ; ++j){
	  a.push_back((genrand_int32()%100));
	}
	aa.push_back(a);
  }
  print2DVector("aa",aa);


  const int INDEX = 5;
  const int TSIZE=aa.at(INDEX).size();
  int b[TSIZE];
  for(int i = 0 ; i < TSIZE ; ++i){
	b[i]=aa.at(INDEX).at(i)+1;
  }

  appendToAllConfs(b,TSIZE,aa);

  print2DVector("aa_prime",aa);
  
}
