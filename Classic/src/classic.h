//***** Utilities *****//
void printSeed(const char *c){
  printf("%s : seed %u, %d %d %d\n",c,seed_t,rand(),rand(),rand());
}

int getDistinctColors(const vector<int> &v){
  auto distinctColors = 0;
  vector<int>ColorCount;
  for(auto &x : v){
    while(ColorCount.size() <= x){
      ColorCount.push_back(0);
    }
    ColorCount.at(x)++;

    if(ColorCount.at(x) == 1){//if it's the first time added
      distinctColors++;
    }
  }
  return distinctColors;
}

int getConflictOfVertex(const vertex *pVertex, const int color[]){
  int nConflicts=0;
  vector<int>adjV = pVertex->adj;
  for (auto &x : adjV){
    if(BB) assert(pVertices[x]->id == x);
    if(BB) assert(getDIMACSBinaryEdgeSwap(pVertex->id, x));

    if(color[pVertex->id]==color[x]){
      nConflicts++;
    }
  }
  //printf("vertex %d has %d conflict\n",pVertex->id,nConflicts);
  return nConflicts;
}

int updateConflictTable(const int color[], int conflict[],
			const int &nVertices){
  int totalConflicts = 0;
  for(auto i = 0 ; i < nVertices ; ++i){
    conflict[i] = getConflictOfVertex(pVertices[i],color);
    totalConflicts += conflict[i];
  }
  return totalConflicts;
}

void cleanUp(const int &nVertices){
  for (auto i = 0; i < nVertices; ++i){delete pVertices[i];} delete pVertices;
  for (auto i = 0; i < vAnts.size(); ++i){delete vAnts.at(i);}
}

void printSol(const bool &write,
	      const int &bestResult,
	      const int &bestCycle,	      
	      const int &nVertices, const int &nEdges){
  int nthreads=1; //temp value, will be change when doing multi-threads
  if(BB) assert(getDistinctColors(bestColorResult) == bestResult);
  printf("nthreads: %d colors: %d vertices: %d edges %d "
	 "colors_table: () bestIndex: %d seed: %u\n",
	 nthreads, bestResult, nVertices, nEdges, bestCycle, seed_t);
  
  if (write){//write to file
    const char * solFile = "sol.txt";
    printf("sol written to '%s'\n", solFile);
    fstream os(solFile, ios::out);

#ifndef KA  
    for(auto i = 0; i < nVertices ; ++i){
      //+1  b/c color index starts from 1 instead of 0
      os << bestColorResult.at(i)+1 << "\n";  
    }
#else //write to file using KA's format, just so able to run KA's verifier
#warning "Solution output in KA style"
    os << "c FILENAME: " << inputFile << "\n";
    os << "c bestcolor: " << bestResult << "\n";
    os << "c maxVertices: " << nVertices << "\n";
    os << "c actual Edges: "<< nEdges << "\n";
    os << "p edges " <<nVertices << " " << nEdges << " " << bestResult <<"\n";

    for(auto i = 0; i < nVertices ; ++i){
      os << "v " << i+1 << " " << bestColorResult.at(i)+1 << "\n";
    }

#endif
  }

#ifdef V
  printSeed("End");  //print out the seed and some rand for debugging purpose
#endif

}

int compare_degree (const void *v1, const void *v2 ){
  vertex* a1 = *(vertex**)v1;  vertex* a2 = *(vertex**)v2;
  return a2->adj.size() - a1->adj.size();
}

void setGreedyClique(vector<vertex *>&theGreedyClique,const int &nVertices){
  vertex *pVertex ;
  vertex **sv = new vertex *[nVertices];
  int status[nVertices];
  for (int i=0;  i < nVertices; i++) {
    pVertex = pVertices[i];
    status[i] = 0;
    sv[i] = pVertex;
  }
  qsort((void*) sv, (size_t) nVertices, sizeof(vertex*), compare_degree);

  for (auto i = 0; i < nVertices; i++){
    pVertex = sv[i];

    if (status[sv[i]->id] == 0){
      status[sv[i]->id] = 1;
      theGreedyClique.push_back(pVertex);

      for (int j = i+1 ; j < nVertices;++j){
	if (!getDIMACSBinaryEdgeSwap(pVertex->id,sv[j]->id)){
	  status[sv[j]->id]=-1;
	}
      }//end for j
    }//end     if (pVertex->status == POSSIBLY_IN_OUTPUT){
  }

#ifdef V
  printf("Clique Size %d\n",theGreedyClique.size());
#endif

  delete sv;
}

//***** Initialize *****//
void initAnts(const int &nAnts){
  for(auto i = 0 ; i < nAnts ; ++i){
    ant *a = new ant();
    a->id = i;
    a->current = nullptr;
    a->old = nullptr;
    vAnts.push_back(a);
  }
}

void initVerticesAndEdges(const int &nVertices, int &nEdges){
  pVertices = new vertex* [nVertices];  
  auto edgeCount = 0; 
  for(auto i = 0;  i < nVertices; ++i){
    vertex *vPtr = new vertex;
    vPtr->id = i ;
    vPtr->numAnts = 0;

    pVertices[i]= vPtr;

    //todo: when done, change EdgeSwap to just Edge for performance
    for (auto j = 0 ; j < i ;++j){
	  
      if (getDIMACSBinaryEdgeSwap(i,j)){
	//printf("e %d %d\n",i,j);
	
	pVertices[i]->adj.push_back(j);
	pVertices[j]->adj.push_back(i);
	edgeCount++;
      }//end if (getDIMACSBinaryEdgeSwap(i,j))
    }//end for j 
  }//end for i


  /*check for cases when Dimacs graphs incorrectly 
    have duplicate edges (i.e., edge i,j then edge j,i)*/
  
  if(edgeCount != nEdges){
    if(BB)
      printf("W: edgecount %d != nEdges %d, setting nEdges=edgeCount!!!!\n",
	     edgeCount,nEdges);
    nEdges = edgeCount;
  }
}

void copyBestResult(const int color[], const int &nVertices){
  if(BB)assert(bestColorResult.size()==nVertices);
  for(auto i = 0 ; i < nVertices ; ++i){
    bestColorResult.at(i) = color[i];
  }
}

int initialize(char const * inputFile,
		int &nVertices, int &nEdges, int &nAnts,
		int &nCycles, int &breakCycles, 
		int &moveLimit,
		int &rSizeLimit,
		int & nRLFSetLimit){
  
  //read_graph_DIMACS_ascii(inputFile,nVertices,nEdges);
  readDIMACSBinaryFormat(inputFile, nVertices, nEdges);
  for(auto i = 0; i < nVertices; ++i){
    bestColorResult.push_back(i);
  }
  
  //initialize data structures
  initVerticesAndEdges(nVertices, nEdges);
  
  //other settings
  nCycles = nVertices * nCyclesFactor;
  if(nCycles > 4000) nCycles = 4000;
  
  nAnts = (int)(nVertices * nAntsPercent);
  if(nAnts > 100) nAnts = 100;
  
  breakCycles = int(nCycles / 2);
  if(breakCycles > 1500) breakCycles = 1500;
    
  moveLimit = (nVertices > 100) ?
    20 + (int)(nVertices / nAnts) : (int) (moveLimitPercent*nVertices);
  
  rSizeLimit = (int)(moveLimit / rSizeLimitFactor);
  if(rSizeLimit < 1) rSizeLimit = 1;
  
  nRLFSetLimit = (int)(nVertices * RLFSetPercent);

  if(BB){
    printf("*******************\n");
    printf("Graph %s, nVertices %d, nEdges %d\n",inputFile,nVertices,nEdges);
    printf("nCycles %d, ",nCycles);
    printf("nAnts %d\n",nAnts);
    printf("breakCycles %d, ",breakCycles);
    printf("moveLimit %d, ",moveLimit);
    printf("rSizeLimit %d, ",rSizeLimit);
    printf("nRLFSetLimit %d",nRLFSetLimit);
    printf("\n*******************\n");
  }

  return bestColorResult.size();
}


//***** XLRF *****//
void updateDegreeB(const vector<int> &adj, const bool W[], int degreeB[]){
  vertex *tV; //temp vertex
  for(auto &x : adj){
    tV=pVertices[x];
    for(auto &y : tV->adj){
      if (W[y]){//if it's uncolored and safe, update its blacklist degree
	degreeB[y]++;
      }
    }
  }
}

const int getLargestDegreeB(const bool W[], const int degreeB[], const int &nVertices){
  int chosenV=-1, maxDeg=0;
  for(auto i = 0; i < nVertices ;++i){
    if(W[i] && degreeB[i]>maxDeg){
	maxDeg = degreeB[i]; 
	chosenV = i;
    }
  }

  if(chosenV==-1){
    int adjCount=0;
    for(auto i = 0 ; i < nVertices ;++i){
      if(W[i]){
	if(chosenV==-1){
	  chosenV=i; //desperately seeking for a vertex
	} 
	adjCount=pVertices[i]->adj.size();
	if(adjCount>maxDeg){
	  maxDeg=adjCount;
	  chosenV=i;
	}
      }
    }
  }

  return chosenV;
}

void markW(const vector<int> &adj, bool W[], int &vSize){
  for(auto &x : adj){
    if (W[x]){ //if it's NOT colored and is considered
      W[x] = false;//no longer consider it
      vSize--;
    }
  }
}

void markBlackList(const vector<int> &adj, const bool isColored[],
		   bool B[], int &vSize){
  for(auto &x : adj){
    if (!isColored[x] && !B[x]){//if it's UNcolored and not yet in blacklist
      B[x]=true;
      vSize++;
    }
  }//and for 
}


void XRLF(int colorAssigned[], int &bestResult,
	 const int &nVertices, const int &nEdges,
	 const int &nRLFSetLimit){

  bool isColored[nVertices];//stable set
  bool W[nVertices];//uncolored but can be included in stable set
  bool B[nVertices];//uncolored but can NO longer belong to stable set
  int degreeB[nVertices]; //degree of B that is inside W 

  auto wSize = nVertices;
  int bSize ;

  auto currentColor = -1;
  auto nVerticesWithThisColored = 0;
  
  for(auto i = 0 ; i < nVertices ;++i){
    colorAssigned[i] = -1; 
    isColored[i] = false;
    W[i] = true; //initially all can be included
  }
  
  int maxTest, chosenV, vAdjSize;   //temp vals
  while((nVerticesWithThisColored < nVertices) && (wSize > 0)){
    currentColor++; //current color 
    nVerticesWithThisColored=0; //no vertex has this color yet 

    //reset Blacklist
    for(auto i = 0 ; i < nVertices ; ++i){
      B[i] = false;
      degreeB[i] = 0;
    }
    bSize = 0;
    chosenV = -1;   //reset chosen

    //adj density based, vertex w/ largest density chosen for MXRLF    
    maxTest = 0;
    for(auto i = 0; i < nVertices; ++i){
      if (W[i]){//only check UNcolored and safe vertices
	vAdjSize = pVertices[i]->adj.size();
	if (maxTest < vAdjSize){
	  maxTest = vAdjSize;
	  chosenV = i;
	}
      }
    }
	
    if(BB) assert(chosenV!=-1);

    colorAssigned[chosenV] = currentColor;//incre then assign color

    isColored[chosenV] = true;
    nVerticesWithThisColored++;

    W[chosenV]=false;
    wSize--; //no longer considered

    markW(pVertices[chosenV]->adj, W, wSize);//marked W[adjV] to false 

    while(wSize > 0){
      //mark "UNcolored" adj ones to black list
      markBlackList(pVertices[chosenV]->adj,isColored,B,bSize);


      if(nVerticesWithThisColored > nRLFSetLimit){
	bSize=0;
	chosenV=-1;//reset
	for(auto i = 0 ; i < nVertices ; ++i){
	  if(!isColored[i]){//if it's not colored,add it to blacklist
	    B[i]=true;
	    bSize++;
	  }
	}

#ifdef V
	printf("break out of while loop\n");
#endif
	break; //break out of the loop
		
      }

      //update degree of B in respect to W
      updateDegreeB(pVertices[chosenV]->adj,W,degreeB);
      chosenV = getLargestDegreeB(W,degreeB, nVertices);

      colorAssigned[chosenV]=currentColor;//incre then assign color

      isColored[chosenV]=true;
      nVerticesWithThisColored++;

      W[chosenV]=false;
      wSize--; //no longer considered
	  
      markW(pVertices[chosenV]->adj,W,wSize);//marked W[adjV] to false 
	  
      //printf("coloring %d with color %d\n", chosenV,currentColor);
	  
    }//end while(wSize>0)
	
    markBlackList(pVertices[chosenV]->adj,isColored,B,bSize);
	
    wSize = 0;
    for(auto i=0; i < nVertices; ++i){
      W[i] = B[i];
      //reset the count of 1's in W                                                               
      if(B[i]) wSize++;
    }
  }//while loop

  //save solution
  for(auto i = 0; i < nVertices; ++i){
    if (colorAssigned[i]<0){//if vertex is uncolored
      currentColor++;
      colorAssigned[i]= currentColor;
      printf("current color %d\n",currentColor);
    }
  }
  
  int XRLF_Colors = currentColor + 1; //since color index starts from 0 ;

  //save the best (lowest) colors
  //this won't even need to be check if this function is called 
  //in the first time (since bestNumColor == nVertices init)
  if(BB) assert(XRLF_Colors <= bestResult);
  if (XRLF_Colors < bestResult){
    bestResult = XRLF_Colors;
    copyBestResult(colorAssigned, nVertices);
  }
}

void setUpColorClasses(int colorAssigned[],
		       const int colorOrig[],
		       const int &colorsUsedXRLF,
		       const int &nVertices){

  int betaLimit=(int)(colorsUsedXRLF*BETA);
  int vertexRecoloredCounter=0;
  int betaReColorIndex=-1;
  bool toBeRecolor[nVertices];
  vector<int>keepColorsV,unkeepColorsV;
  

  for(auto i = 0 ; i < nVertices ;++i) toBeRecolor[i]=true;
  
  for(auto j = 0 ; j< colorsUsedXRLF ;++j){
    if(rand()%1000<=BETA*1000){//if choose to reIndex this color
      keepColorsV.push_back(j);
    }
    else{
      unkeepColorsV.push_back(j);
    }
	
  }

  //  printf("before, beta limit %d, keep %d, unkeep %d\n",betaLimit,keepColorsV.size(),unkeepColorsV.size());  
  if(BB)assert(keepColorsV.size()+unkeepColorsV.size()==colorsUsedXRLF);

  while(keepColorsV.size() > betaLimit){
    int tIndex = rand()%keepColorsV.size();
    unkeepColorsV.push_back(keepColorsV.at(tIndex));
    keepColorsV.erase(keepColorsV.begin()+tIndex);
  }
  while(keepColorsV.size() + 1 <= betaLimit){
    int tIndex = rand()%unkeepColorsV.size();
    keepColorsV.push_back(unkeepColorsV.at(tIndex));
    unkeepColorsV.erase(unkeepColorsV.begin()+tIndex);
  }

  //  printf("after, beta limit %d, keep %d, unkeep %d\n",betaLimit,keepColorsV.size(),unkeepColorsV.size());  
  if(BB)assert(keepColorsV.size()+unkeepColorsV.size()==colorsUsedXRLF);


  //recolor & reindex
  for(auto &x: keepColorsV){
    betaReColorIndex++;
    for(auto i = 0 ; i < nVertices ;++i){
      if(colorOrig[i] == x){
	//printf("choosing color %d==%d, vertex %d\n",j,colorAssigned[i],i);
	if(BB)assert(toBeRecolor[i]);
	vertexRecoloredCounter++;
	colorAssigned[i] = betaReColorIndex;
	toBeRecolor[i] = false;
      }
    }
  }
	   
  int lg=0;
  if(BB){
    lg=0;
    for(int i=0;i <nVertices ;++i){
      if(colorAssigned[i]>lg)lg=colorAssigned[i];
    }

#ifdef V
    printf("Q BETA is %d, percent %f, XG %d, "
	   "recolored vertices %d, largest Color Index %d\n",
	   keepColorsV.size(),
	   (double)keepColorsV.size()/(double)colorsUsedXRLF,
	   colorsUsedXRLF,vertexRecoloredCounter,lg);
#endif

    if(lg==0){
      printf("lg = %d\n",lg);
      printf("(QColor)largest color index %d\n",lg);
      assert(false);
    }
  }

  vertexRecoloredCounter=0; //reset  , can be removed , only for debug purpose


  //random color distribution
  int deltaNumColors=(int)(colorsUsedXRLF*DELTA); //limit the numbers to 0:availColor  

  //probably I can just pick the colors RANDOMLY
  for(auto i = 0 ; i < nVertices ;++i){
    if(toBeRecolor[i]){//if can be colored
      if(BB)assert(colorAssigned[i]==-1);
      colorAssigned[i]=rand()%(deltaNumColors);
      vertexRecoloredCounter++;
    }
  }

  lg=0;
  if(BB){
    lg=0;
    for(int i=0;i <nVertices ;++i){
      if(colorAssigned[i]>lg)lg=colorAssigned[i];
    }

#ifdef V
    printf("Q DELTA is %d, percent %f, XG %d, "
	   "recolored vertices %d, largest Color Index %d\n",
	   deltaNumColors,
	   (double)deltaNumColors/(double)colorsUsedXRLF,
	   colorsUsedXRLF,vertexRecoloredCounter,lg);
#endif

    if(lg==0){
      printf("lg = %d\n",lg);
      printf("(QColor)largest color index %d\n",lg);
      assert(false);
    }
  }
}

int chooseInitialMoveKA(const MoveOpts &MOVE_METHOD,
			const int conflictTable[],
			const int &nVertices){
  int movePos=-1;

  if(MOVE_METHOD == MoveOpts::Random){movePos=rand()%nVertices; }

  if(MOVE_METHOD == MoveOpts::MaxConflict){
    movePos=0;

    for(int i=1;i<nVertices;++i){
      if(BB)assert(conflictTable[i]!=-1);

      if(conflictTable[i] > conflictTable[movePos]){
	movePos=i;//tvn todo: what if more than 1 has max conflicts ?
      }
    }//end for

    if(BB){
      for(auto i = 0 ;  i< nVertices;++i)
	assert(conflictTable[i]<=conflictTable[movePos]);
    }
  }

  return movePos;
}

void color(const ant* theAnt, 
	   int colorsTable[], int conflictsTable[],
	   const int &maxNColor,
	   const int &nVertices){
  
  //  printf("Max Color %d\n",maxNColor);

  auto currPos = -1, oldPos = -1; 

  if(BB)assert(theAnt->current != nullptr);

  if(theAnt->old != nullptr) oldPos=theAnt->old->id;
  if(theAnt->current != nullptr) currPos=theAnt->current->id;

  vector<int> adjV = theAnt->current->adj;
  
  //if(BB)assert(currPos!=oldPos);
  
  int oldColor=colorsTable[currPos];
  int colorInAdj[maxNColor];
  for(auto i = 0 ; i < maxNColor; ++i){colorInAdj[i]=0;}

  for(auto &x : adjV){
    if(BB){
      if(colorsTable[x] >= maxNColor){
	printf("color in vertex %d is %d, maxNColor %d\n",
	       x, colorsTable[x], maxNColor);	  
	assert(false);
      }
    }
    colorInAdj[colorsTable[x]]++;
  }

  //avoid choosing color from old POS
  if(getDIMACSBinaryEdgeSwap(oldPos,currPos)){
    //printf("Will not choose color %d\n", colorsTable[oldPos]);
    colorInAdj[colorsTable[oldPos]] = 2 * nVertices;
  }


  vector<int>leftOver;
  for(auto i = 0 ; i < maxNColor; ++i){
     //if no adj has this color, then push it to leftOver
    if(!colorInAdj[i]) leftOver.push_back(i);
  }

  if(leftOver.size()){
    colorsTable[currPos] = leftOver.at(rand()%leftOver.size());
    conflictsTable[currPos] = 0;
    
    if(BB)assert(getConflictOfVertex(pVertices[currPos],colorsTable)==0);

    //assert(oldColor!=colorsTable[currPos]);
  }
  else{
    int minCount=colorInAdj[0];
    for(auto i = 1 ; i < maxNColor; ++i){
      if(colorInAdj[i] < minCount){
	minCount = colorInAdj[i];
      }
    }
	
    if(BB){ for(auto i = 0 ; i < maxNColor ; ++i){ assert(colorInAdj[i]>=minCount);}}
    vector<int> MinIndexVector;
    for(auto i = 0 ; i < maxNColor ; ++i){
      if(colorInAdj[i]==minCount)MinIndexVector.push_back(i);
    }
    colorsTable[currPos]=MinIndexVector.at(rand()%MinIndexVector.size());
    conflictsTable[currPos]=getConflictOfVertex(pVertices[currPos],colorsTable);
  }

  //  printf("after %d, %d, color %d\n",getConflictOfVertex(pVertices[currPos],colorsTable),getTotalConflict(colorsTable),colorsTable[currPos]);

  //update only the conflicts LOCALLY
  if(colorsTable[currPos]!=oldColor){//if changed the Color
	
    //decrement since the colors has changed
    for(auto &x : adjV){
      if(colorsTable[x] == oldColor){
	conflictsTable[x]--;
      }
    }
	
    //if there's still conflict, and if they are the same color , incre
    if(conflictsTable[currPos]){//if still conflict  , TVN: not sure if this should be done , doesn't make too much sense
      for(auto &x : adjV){
	if(colorsTable[x] == colorsTable[currPos]){
	  conflictsTable[x]++;
	}
      }
    }//if still conflict

	
  }
}

//move the ants, very frequently used operation
int move(const int &curPos, const MoveOpts &MOVE_METHOD, 
	 const vector<int> &recentlyVisited,bool considerVertex[],
	 const int conflictTable[],
	 const int &nVertices){
  
  int movePos=-1;
  vector<int> adjV = pVertices[curPos]->adj;
  vector<int>chosenOne;

  if(!adjV.empty()){//if it has some adj vertices
    for(auto i=0 ; i < nVertices ;++i) considerVertex[i] = false;//reset
    //consider these adj's
    for(auto &x : adjV) considerVertex[x] = true; 
    //don't considered ones recently visited
    for(auto &x : recentlyVisited) considerVertex[x] = false;
	
    if(MOVE_METHOD == MoveOpts::Random){	
      for(auto &x : adjV){
	if(considerVertex[x]){
	  if(BB) assert(getDIMACSBinaryEdgeSwap(curPos, x));
	  chosenOne.push_back(x);
	}
      }
      //if(BB){if(chosenOne.empty())printf("chosenOne is empty");}
    }
    if(MOVE_METHOD == MoveOpts::MaxConflict){
      int maxConflict = 0;
      for(auto &x : adjV){
	if(considerVertex[x] && conflictTable[x] > maxConflict){
	    maxConflict = conflictTable[x];
	}
      }//end for
	  
      if(maxConflict){//if there are conflicts
	if(BB)assert(chosenOne.empty());
	
	for (auto &x : adjV){
	  if(considerVertex[x] && conflictTable[x] == maxConflict){
	      chosenOne.push_back(x);
	  }
	}
	
	if(BB)assert(!chosenOne.empty()); 
      }
	  
    }
  }
  else{//if this vertex has no adj vertices
    for(auto i = 0; i < nVertices ;++i) considerVertex[i] = true;//reset
    //don't considered ones recently visited
    for(auto &x : recentlyVisited) considerVertex[x] = false;
    for(auto i = 0; i < nVertices ;++i){
      if(considerVertex[i]) chosenOne.push_back(i);
    }
	
  }

  movePos=chosenOne.empty() ?
    adjV.at(rand() % adjV.size()) :
    chosenOne.at(rand() % chosenOne.size());  

  if(BB)assert(movePos!=-1);
  return movePos;
}

void reColorMoreThanQ(int colorTable[], int conflictsTable[],
		      const int &nColors,
		      const int &nVertices){
  if(BB)assert(nColors>0);//if it's 0 then only 1 color
  
  for(auto i = 0; i < nVertices; ++i){
    if(colorTable[i] >= nColors){ 
      colorTable[i] = rand() % nColors;
    }
  }
  updateConflictTable(colorTable,conflictsTable, nVertices);  
}

const int antsOps(int &bestResult, int currentColor[],
		  const int &nVertices, const int &nAnts,
		  const int &nCycles, const int &breakCycles,
		  const int &moveLimit,
		  const int &rSizeLimit){
  
  initAnts(nAnts);//add the ants
  auto alphaNumColors = (int)(bestResult * ALPHA);
  auto bestCycle = 0;
  
  int lg=0;
  if(BB){
    lg=0;
    for(auto i = 0;i < nVertices ;++i){
      if(currentColor[i]>lg)lg=currentColor[i];
    }

#ifdef V
    printf("Q ALPHA is %d, percent %f, XG %d, largest Color Index %d\n",
	   alphaNumColors,(double)alphaNumColors/(double)bestResult,bestResult,lg);
#endif

    if(lg==0){
      printf("lg = %d\n",lg);
      printf("(QColor)largest color index %d\n",lg);
      assert(false);
    }
  }

  int conflictsTable[nVertices];
  int totalConflicts = updateConflictTable(currentColor,conflictsTable, nVertices);  
  //  printf("Total Conflict %d, %d\n",totalConflicts, getTotalConflict(conflictsTable));

  vector<int> recentlyVisited ;
  bool considerVertices[nVertices];

  int changedCycle=0;
  int moveSoFar;

  for (auto iCycles = 0; iCycles < nCycles; ++iCycles) {
    for(auto iAnts = 0; iAnts<nAnts; ++iAnts){
      auto *theAnt = vAnts.at(iAnts);
      moveSoFar = 0;

      if(!recentlyVisited.empty())recentlyVisited.clear();//resest
      //assert(theAnt->current==NULL);
      theAnt->current=pVertices[chooseInitialMoveKA(MOVE_METHOD,
						    conflictsTable,
						    nVertices)];
      moveSoFar++; if(BB)assert(moveSoFar==1);
      color(theAnt, currentColor,conflictsTable, alphaNumColors, nVertices);

      int movePos;
      while(moveSoFar < moveLimit){ 
	if(BB)assert(theAnt->current != nullptr); 
	theAnt->old = theAnt->current;
	for(int distI=0; distI < HOW_FAR; ++distI){
	  ///move to Max Conflict adj lastly
	  movePos = move(theAnt->current->id,
			 (distI==HOW_FAR-1)?MoveOpts::MaxConflict:MoveOpts::Random,
			 recentlyVisited,considerVertices,conflictsTable, nVertices);
	  theAnt->current=pVertices[movePos];

	  //printf("-> [%d] old %d, current %d\n",iAnts,theAnt->old->id,theAnt->current->id);

	  //if reach maximum then remove from the beginning, todo, use a diff structure
	  if(recentlyVisited.size()== rSizeLimit)
	    recentlyVisited.erase(recentlyVisited.begin());
	  recentlyVisited.push_back(theAnt->current->id);  

	}//HOW_FAR

	color(theAnt, currentColor, conflictsTable, alphaNumColors, nVertices);
	moveSoFar++;

      }//while(moveSoFar<moveLimit)
    }//for ant loop


    totalConflicts = updateConflictTable(currentColor, conflictsTable, nVertices);  

    if(!totalConflicts && alphaNumColors<bestResult){
      bestResult=alphaNumColors;
      copyBestResult(currentColor, nVertices);
      bestCycle = iCycles;
      
#ifdef V
      printf("*********** Achieve %d colors at cycle %d ************\n",
	     bestResult, iCycles);
#endif
      
      if(alphaNumColors > 1){
	alphaNumColors--;
#ifdef V
	printf("Cycle [last changed %d, current %d], "
	       "Q_CHANGE_CYCLE %d, decrease nColors from %d to %d\n",
	       changedCycle,iCycles,Q_CHANGE_CYCLE,alphaNumColors+1,alphaNumColors);
#endif
	changedCycle = iCycles;
	reColorMoreThanQ(currentColor, conflictsTable, alphaNumColors, nVertices);
	}
      else{
	fprintf(stderr,"E: attempting to have only 1 color, "
		"only possible if there's no edge .. haha\n");
	assert(false);
      }
    }//if no conflict and better than best result


    if (iCycles - changedCycle == Q_CHANGE_CYCLE){
      int tIncr = (int)(bestResult- alphaNumColors/4);

      if(tIncr==0 || alphaNumColors+tIncr >= bestResult){
	tIncr = 1;
      }
	  
      if(alphaNumColors + tIncr < bestResult){
	alphaNumColors += tIncr;
	changedCycle = iCycles;
#ifdef V
	printf("Cycle [last changed %d, current %d], "
	       "Q_CHANGE_CYCLE %d, increasing nColors from %d to %d\n",
	       changedCycle, iCycles, Q_CHANGE_CYCLE,
	       alphaNumColors-tIncr, alphaNumColors);
#endif
      }
    }

    if(iCycles - changedCycle == breakCycles){
#ifdef V
      printf("Break at cycle %d, last changed cycle %d, breakCycle %d\n",
	     iCycles, changedCycle, breakCycles);
#endif
      break;
    }
  }//cycle loop

  return bestCycle;
}
