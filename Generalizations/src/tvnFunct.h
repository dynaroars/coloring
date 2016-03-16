//***** Function Utitilities *****//
void cleanUp(){ for (int i=0;i<nVertices;++i){delete pVertices[i];} delete pVertices;}


int getDistinctColors(const int v[]){
  if(!distinctColorsV.empty())distinctColorsV.clear();

  vector<int>ColorCount;

  for(int i = 0 ; i < nVertices;++i){
	while(ColorCount.size()<=v[i]){
	  ColorCount.push_back(0);
	}
	ColorCount.at(v[i])++;

	if(ColorCount.at(v[i])==1){//if it's the first time added
	  distinctColorsV.push_back(v[i]);
	}
  }

  return distinctColorsV.size();
}

int getConflictOfVertexW(const int &v, const int color[]){
  int nConflicts=0;
  vector<int>adjV = pVertices[v]->adj;
  for (int i = 0 ; i < adjV.size(); ++i){

	if(BB)assert(edgeMatrix.at(v).at(adjV.at(i)));


	if(abs(color[v]-color[adjV.at(i)])<weightMatrix[v][adjV.at(i)]){
	  /*printf("conflict, v %d, colored %d , v %d colored %d, weight[%d][%d]=%d\n",
			 pVertex->v_id,color[pVertex->v_id],adjV.at(i),color[adjV.at(i)],
			 pVertex->v_id,adjV.at(i),weightMatrix[pVertex->v_id][adjV.at(i)]);
	  */
	  nConflicts++;
	}

  }
  //printf("vertex %d has %d conflict\n",pVertex->v_id,nConflicts);
  return nConflicts;
}


int getConflictOfVertexAHEAD(const int &v, const int color[],const int &AHEAD_DIST){
  
  if(AHEAD_DIST==1){
	return getConflictOfVertexW(v,color);
  }
  else{
	int nConflicts=0;
	vector<int>adjV = pVertices[v]->adj;
	for (int i = 0 ; i < adjV.size(); ++i){
	  nConflicts+=getConflictOfVertexAHEAD(adjV.at(i),color,AHEAD_DIST-1);
	}

	return nConflicts;
  }
}




int updateConflictTableAHEAD(const int color[],int conflictsTableAHEAD[], const int &AHEAD_DIST){
  int totalConflicts=0;
  for(int i = 0 ; i < nVertices ; ++i){
	conflictsTableAHEAD[i]=getConflictOfVertexAHEAD(i,color,AHEAD_DIST);
	totalConflicts+=conflictsTableAHEAD[i];
  }
  return totalConflicts;
}


int updateConflictTableW(const int color[], int conflictsTable[]){
  int totalConflicts=0;
  for(int i = 0 ; i < nVertices ; ++i){
	conflictsTable[i]=getConflictOfVertexW(i,color);
	totalConflicts+=conflictsTable[i];
  }
  return totalConflicts;
}

//return false if there's no conflicts within its adj vertices
bool isConflictAmongAdj(const int &v, const int color[]){
  
  vector<int>adjV=pVertices[v]->adj;
  for(int i = 0 ;i < adjV.size(); ++i){
	if(getConflictOfVertexW(adjV.at(i),color))return true;
  }

  return false;
}






void printSol(const int bestConfig[],const int &MODE){

  int nthreads=1;//temp value, will be change when doing multi-threads

  if(BB){int tempConflicts[nVertices];assert(updateConflictTableW(bestConfig,tempConflicts)==0);  }

  maxK=getMaxK(bestConfig);

  printf("nthreads: %d maxK: %d vertices: %d edges: %d %s: () bestIndex: %d seed: %u\n",
		 nthreads,maxK+1,nVertices,nEdges,PTYPE_name[PTYPE].c_str(),bestCycle,seed_t);

  if (MODE==TOCHECK){//write to file
	printf("Outputing sol to \" %s \" \n",outputFile);
	fstream os(outputFile,ios::out);
	os << "c FILENAME: " << inputFile << "\n";
	os << "c bestcolor: " << maxK +1<< "\n";
    for(int i = 0; i < nVertices ; ++i){
      os << "s " << bestConfig[i]+1 << "\n";  //+1  means the index of color starts from 1 instead of 0     
    }

  }

	
#ifdef V
  printSeed("End",1);  //print out the seed and some rand for debugging purpose
#endif

}




//***** Initialize *****//
void initVerticesAndEdges(const int &nVertices, int &nEdges){


  pVertices = new vertex* [nVertices];  

  int edgeCount=0; 
  
  for (int i = 0 ;  i < nVertices ; ++i){
	vertex *vPtr = new vertex;	vPtr->v_id = i ; 
	pVertices[i]=vPtr;
	  
	
	for (int j = 0 ; j < i ;++j){
	  if(BB)assert(j!=i);
	  if (edgeMatrix.at(i).at(j)){
		pVertices[i]->adj.push_back(j);
		pVertices[j]->adj.push_back(i);
		edgeCount++;
		if(BB)assert(weightMatrix.at(i).at(j)!=0);
		if(weightMatrix.at(i).at(j)>maxWeight)maxWeight=weightMatrix.at(i).at(j);
		//printf("e %d %d %d\n",i+1,j+1,weightMatrix.at(i).at(j));
	  }
	  
	}
  }
  
  
  //there are cases when the Dimacs graphs incorrectly have duplicate edges (i.e., edge i,j then edge j,i)
  if(edgeCount!=nEdges){
	if(BB)printf("WARNING: edgecount %d != nEdges %d , setting nEdges=edgeCount!!!!\n",edgeCount,nEdges);
	nEdges=edgeCount;
  }
}



void initGraph(const int &argc,char *argv[]){

  randInit(argc,argv);

  read_graph_DIMACS_ascii(inputFile,nVertices,nEdges,edgeMatrix,weightMatrix,ku);
  


  if(PTYPE==GCP||PTYPE==MCP){//ignore weight matrix, weight = 1
	for(int i =  0 ;  i < weightMatrix.size() ; ++i){
	  for(int j = 0 ; j < weightMatrix.at(i).size();++j){
		weightMatrix.at(i).at(j)=1;
	  }
	}
  }

  if(PTYPE==GCP||PTYPE==BCP){//ignore ku, ku = 1 (each set size = 1)
	for(int i = 0 ; i < ku.size() ; ++i){
	  ku.at(i)=1;
	}
  }



  if(PTYPE==BMCP||PTYPE==MCP){//clique translation
	initSuperSet(nVertices,nEdges);
  }
  
  



  initVerticesAndEdges(nVertices,nEdges);
  maxK=nVertices*2;


  nCycles=nVertices*nCyclesFactor;if(nCycles>4000)nCycles=4000;//if(nCycles<1000)nCycles=1000;
  nAnts=(int)(nVertices*nAntsPercent);if(nAnts>100)nAnts=100;//if(nAnts<60)nAnts=60;
  joltCycles=(int)(nCycles*nJoltPercent); if(joltCycles>300)joltCycles=300;
  breakCycles=int(nCycles/nBreakCycleFactor); if(breakCycles>1500)breakCycles=1500;
  nMoves=(nVertices>100)?20+(int)(nVertices/nAnts):(int)(nMovesPercent*nVertices);
  rSizeLimit = (int)(nMoves/rSizeLimitFactor);if(rSizeLimit<1)rSizeLimit=1;



#ifdef V
  printf("\n*******************\n");
  printf("Graph %s, nVertices %d, nEdges %d, maxWeight %d, maxK %d\n",inputFile,nVertices,nEdges,maxWeight,maxK);
  printf("nCycles %d, ",nCycles);
  printf("nAnts %d\n",nAnts);
  printf("joltCycles %d, ",joltCycles);
  printf("breakCycles %d, ",breakCycles);
  printf("qChangeCycles %d, ",Q_CHANGE_CYCLE);
  printf("nMoves %d, ",nMoves);
  printf("rSizeLimit %d, ",rSizeLimit);
  printSeed("\nTesting",0);
  printf("*******************\n");

#endif

}


void setUpColorClasses(const int &nDistinctColors, int setupColorConfigA[]){


  //printArray("old colorA: ",setupColorConfigA);

  int betaLimit=(int)ceil(nDistinctColors*BETA);

  //  printVector("distinctColorsV: ", distinctColorsV);

  vector<int>availV,selectedV;
  for(int i=0;i<distinctColorsV.size();++i)availV.push_back(distinctColorsV.at(i));

  int smallestAmongSelected=nVertices;//some temp big number
  int tNum;
  while(selectedV.size()<betaLimit){
	int tIndex=genrand_int32()%availV.size();

	tNum=availV.at(tIndex);
	selectedV.push_back(tNum);
	availV.erase(availV.begin()+tIndex);


	//avoid (but still select it) choosing the first color, 0 , 
	//since that limit the choice of deleting colors less than it
	if (tNum&&tNum<smallestAmongSelected){smallestAmongSelected=tNum;}
  }

  if(BB){
	printVector("selectedV: ", selectedV);
	printVector("availV: ", availV);
  }



  if(BB)printf("Smallest %d, maxK %d\n",smallestAmongSelected,maxK+1);


  /*Recoloring all vertices with colors <= smallestAmongSelected*/

  //  printArray("old colorA: ",setupColorConfigA,nVertices);

  for(int i=0; i<nVertices ;++i){
	if(setupColorConfigA[i]<smallestAmongSelected){
	  setupColorConfigA[i]=0;
	}
	else{
	  setupColorConfigA[i] -= smallestAmongSelected;

	  //also limit the "maxK"
	  //	  setupColorConfigA[i] = ((double)setupColorConfigA[i]>(maxK*2)/3)?(int)(maxK*2/3):setupColorConfigA[i];

	  if(setupColorConfigA[i]>(int)maxK*2/3){
		int oldC=setupColorConfigA[i];
		
		//IMPORTANT:  todo :  change 2/3 to a variable called VFactor and make sure the 
		//in AntsOp , the init avail color (q) will be maxK * this same VFactor
		
		setupColorConfigA[i]=(int)maxK*2/3;
		//printf("maxK %d, color for %d is %d (old) %d (new)\n",maxK,i,oldC,setupColorConfigA[i]);
	  }
	}



  }//recoloring loop


  //printArray("new colorA: ",setupColorConfigA);



  if(BB){
	int conflictsTable[nVertices];
	printf("SetupColorClasses results in %d conflicts\n",updateConflictTableW(setupColorConfigA,conflictsTable));

	printSeed("End of SetupColorClasses",0);
  }

}//end of setUpColorClasses(const int &nDistinctColors, const int origColorA[])






void reColorMoreThanQ(const int &nColors, int colorTable[], 
					  int conflictsTable[],int conflictsTableAHEAD[]){
  if(BB)assert(nColors>0);//if it's 0 then only 1 color
  
  //orig , KA's
  //  for(int i = 0; i < nVertices; ++i){
  //	colorTable[i]--;
  //	if(colorTable[i]<0){colorTable[i]=0;}
  //  }

  for(int i = 0; i < nVertices; ++i){
  	if(colorTable[i]>nColors){
	  colorTable[i]=genrand_int32()%nColors;
	}
  }


  updateConflictTableW(colorTable,conflictsTable);  
  updateConflictTableAHEAD(colorTable,conflictsTableAHEAD,AHEAD_DIST);  

  if(BB && AHEAD_DIST==1)for(int i=0;i<nVertices;++i)assert(conflictsTable[i]==conflictsTableAHEAD[i]);
}
