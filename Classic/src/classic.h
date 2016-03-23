void printSeed(const char *c, const time_t &seed_t){
  printf("%s : seed %u, %d %d %d\n",c,seed_t,rand(),rand(),rand());
}

const int getDistinctColors(const vector<int> &v){
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

void printSol(const Sol &sol,
	      const time_t &seed_t,	      
	      const int &nEdges,
	      const bool &write_to_file){

  int nthreads=1; //temp value, will be change when doing multi-threads
  if(BB) assert(getDistinctColors(sol.colors) == sol.nColors);
  printf("nthreads: %d colors: %d vertices: %d edges: %d "
	 "bestCycle: %d seed: %u\n",
	 nthreads, sol.nColors, sol.colors.size(), nEdges, sol.iCycle, seed_t);
  
  if (write_to_file){
    const char * solFile = "sol.txt";
    printf("sol written to '%s'\n", solFile);
    fstream os(solFile, ios::out);
    for(auto &color : sol.colors){      
      os << color + 1 << "\n";  //+1  b/c color index starts from 1 instead of 0
    }
  }

#ifdef V
  printSeed("End", seed_t);  //print out the seed and some rand for debugging purpose
#endif
}

const int getConflictOfVertex(Vertex *pVertex, const vector<int> &colors){
  int nConflicts=0;
  for (auto &x : pVertex->adj){
    if(BB) assert(getDIMACSBinaryEdgeSwap(pVertex->id, x));
    if(colors.at(pVertex->id) == colors.at(x)){
      nConflicts++;
    }
  }
  return nConflicts;
}

const int updateConflictTable(const vector<int> &colors,
			      int conflict[],
			      Vertex **pVertices){
  int totalConflicts = 0;
  for(auto i = 0; i < colors.size(); ++i){
    conflict[i] = getConflictOfVertex(pVertices[i], colors);
    totalConflicts += conflict[i];
  }
  return totalConflicts;
}

Vertex **initVerticesAndEdges(const int &nVertices, int &nEdges){
  Vertex **vertices = new Vertex *[nVertices];  
  auto edgeCount = 0; 
  for(auto i = 0;  i < nVertices; ++i){
    Vertex *vPtr = new Vertex;
    vPtr->id = i;
    vertices[i]= vPtr;

    for (auto j = 0; j < i; ++j){
      if (getDIMACSBinaryEdgeSwap(i,j)){
	vertices[i]->adj.push_back(j);
	vertices[j]->adj.push_back(i);
	edgeCount++;
      }
    }
  }

  /*check for cases when Dimacs graphs incorrectly 
    have duplicate edges (i.e., edge i,j then edge j,i)*/
  if(edgeCount != nEdges){
    if(BB)
      printf("W: edgecount %d != nEdges %d, setting nEdges=edgeCount!!!!\n",
	     edgeCount,nEdges);
    nEdges = edgeCount;
  }
  return vertices;
}

void updateSol(Sol &sol, const vector<int> &colors, const int &nColors,
	       const int &iCycle){
  sol.nColors = nColors;
  sol.colors = colors;
  sol.iCycle = iCycle;
}

//***** XLRF *****//
void updateDegreeB(const vector<int> &adj, const bool W[], int degreeB[],
		   Vertex **pVertices){
  for(auto &x : adj){
    for(auto &y : pVertices[x]->adj){
      if (W[y]){//if it's uncolored and safe, update its blacklist degree
	degreeB[y]++;
      }
    }
  }
}

const int getLargestDegreeB(const bool W[], const int degreeB[],
			    Vertex **pVertices,
			    const int &nVertices){
  int chosenV=-1, maxDeg=0;
  for(auto i = 0; i < nVertices; ++i){
    if(W[i] && degreeB[i] > maxDeg){
	maxDeg = degreeB[i]; 
	chosenV = i;
    }
  }

  if(chosenV==-1){
    int adjCount=0;
    for(auto i = 0; i < nVertices; ++i){
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

void XRLF(Sol &sol, Vertex **pVertices, const int &nVertices){
  vector<int> colors (nVertices);
  bool isColored[nVertices];//stable set
  bool W[nVertices];//uncolored but can be included in stable set
  bool B[nVertices];//uncolored but can NO longer belong to stable set
  int degreeB[nVertices]; //degree of B that is inside W 

  auto wSize = nVertices;
  int bSize;

  auto currentColor = -1;
  auto nVerticesWithThisColored = 0;
  auto nRLFSetLimit = (int)(nVertices * RLFSetPercent);
  
  for(auto i = 0; i < nVertices; ++i){
    colors.at(i) = -1; 
    isColored[i] = false;
    W[i] = true; //initially all can be included
  }

  int maxTest, chosenV, vAdjSize;   //temp vals
  while((nVerticesWithThisColored < nVertices) && (wSize > 0)){
    currentColor++; //current color 
    nVerticesWithThisColored = 0; //no vertex has this color yet 

    //reset Blacklist
    for(auto i = 0; i < nVertices; ++i){
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
	
    if(BB) assert(chosenV != -1);
    colors.at(chosenV) = currentColor;//incre then assign color
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
	for(auto i = 0; i < nVertices; ++i){
	  if(!isColored[i]){//if it's not colored,add it to blacklist
	    B[i]=true;
	    bSize++;
	  }
	}

#ifdef V
	printf("break out of while loop\n");
#endif
	break; //break out of loop
      }

      //update degree of B in respect to W
      updateDegreeB(pVertices[chosenV]->adj, W, degreeB, pVertices);
      chosenV = getLargestDegreeB(W, degreeB, pVertices, nVertices);

      colors.at(chosenV) = currentColor;//incre then assign color

      isColored[chosenV]=true;
      nVerticesWithThisColored++;

      W[chosenV]=false;
      wSize--; //no longer considered
	  
      markW(pVertices[chosenV]->adj,W,wSize);//marked W[adjV] to false 
	  
      //printf("coloring %d with color %d\n", chosenV,currentColor);
	  
    }//end while(wSize>0)
	
    markBlackList(pVertices[chosenV]->adj, isColored, B, bSize);
	
    wSize = 0;
    for(auto i=0; i < nVertices; ++i){
      W[i] = B[i];
      //reset the count of 1's in W                                                               
      if(B[i]) wSize++;
    }
  }//while loop

  //save solution
  for(auto i = 0; i < nVertices; ++i){
    if (colors.at(i) < 0){//if vertex is uncolored
      currentColor++;
      colors.at(i) = currentColor;
      printf("current color %d\n",currentColor);
    }
  }
  
  auto nColors = currentColor + 1; //since color index starts from 0;
  auto iCycle = -1; //not in any cycle yet
  updateSol(sol, colors, nColors, -1);
}

vector<int> setUpColorClasses(const vector<int> &colors, const int &nColors){

  const auto nVertices = colors.size();
  int betaLimit=(int)(nColors*BETA);
  int vertexRecoloredCounter=0;
  int betaReColorIndex=-1;
  bool toBeRecolor[nVertices];
  vector<int>keepColorsV, unkeepColorsV;

  vector<int>tmp_colors(nVertices);
  for(auto i = 0; i < nVertices; ++i){
    toBeRecolor[i]=true;
    tmp_colors.at(i) = -1;
  }
  
  for(auto i = 0; i< nColors; ++i){
    if(rand() % 1000 <= BETA*1000){
      keepColorsV.push_back(i);
    }
    else{
      unkeepColorsV.push_back(i);
    }
  }

  //  printf("before, beta limit %d, keep %d, unkeep %d\n",betaLimit,keepColorsV.size(),unkeepColorsV.size());  
  if(BB)assert(keepColorsV.size()+unkeepColorsV.size()==nColors);

  while(keepColorsV.size() > betaLimit){
    int tIndex = rand() % keepColorsV.size();
    unkeepColorsV.push_back(keepColorsV.at(tIndex));
    keepColorsV.erase(keepColorsV.begin() + tIndex);
  }
  while(keepColorsV.size() + 1 <= betaLimit){
    int tIndex = rand() % unkeepColorsV.size();
    keepColorsV.push_back(unkeepColorsV.at(tIndex));
    unkeepColorsV.erase(unkeepColorsV.begin() + tIndex);
  }

  if(BB)assert(keepColorsV.size() + unkeepColorsV.size() == nColors);

  //recolor & reindex
  for(auto &x: keepColorsV){
    betaReColorIndex++;
    for(auto i = 0; i < nVertices; ++i){
      if(colors.at(i) == x){
	if(BB)assert(toBeRecolor[i]);
	vertexRecoloredCounter++;
	tmp_colors.at(i) = betaReColorIndex;
	toBeRecolor[i] = false;
      }
    }
  }
	   
  if(BB){
    auto lg = *max_element(tmp_colors.begin(), tmp_colors.end());
#ifdef V
    printf("Q BETA is %d, percent %f, XG %d, "
	   "recolored vertices %d, largest Color Index %d\n",
	   keepColorsV.size(),
	   (double)keepColorsV.size()/(double)nColors,
	   nColors,vertexRecoloredCounter,lg);
#endif
    assert (lg != 0);
  }

  vertexRecoloredCounter=0; //reset  , can be removed , only for debug purpose
  //random color distribution
  auto deltaNumColors=(int)(nColors * DELTA); //limit the numbers to 0:availColor  

  //probably I can just pick the colors RANDOMLY
  for(auto i = 0; i < nVertices; ++i){
    if(toBeRecolor[i]){//if can be colored
      if(BB)assert(tmp_colors.at(i) == -1);
      tmp_colors.at(i) = rand() % deltaNumColors;
      vertexRecoloredCounter++;
    }
  }

  
  if(BB){
    auto lg = *max_element(tmp_colors.begin(), tmp_colors.end());

#ifdef V
    printf("Q DELTA is %d, percent %f, XG %d, "
	   "recolored vertices %d, largest Color Index %d\n",
	   deltaNumColors,
	   (double)deltaNumColors/(double)nColors,
	   nColors,vertexRecoloredCounter,lg);
#endif

    assert (lg != 0);
  }
  
  return tmp_colors;
}

int selectInitMove(const MoveOpts &move_method,
		   const int conflictTable[],
		   const int &nVertices){
  int movePos=-1;

  if(move_method == MoveOpts::Random){
    movePos = rand()%nVertices;
  }
  if(move_method == MoveOpts::MaxConflict){
    movePos = 0;

    for(auto i = 1;i<nVertices; ++i){
      if(BB)assert(conflictTable[i]!=-1);

      if(conflictTable[i] > conflictTable[movePos]){
	movePos=i;//tvn todo: what if more than 1 has max conflicts ?
      }
    }//end for

    if(BB){
      for(auto i = 0;  i< nVertices; ++i)
	assert(conflictTable[i] <= conflictTable[movePos]);
    }
  }

  return movePos;
}

void color(const Ant &ant, 
	   vector<int> &colors,
	   int conflictsTable[],
	   const int &maxNColor,
	   Vertex **pVertices,
	   const int &nVertices){

  if(BB)assert(ant.current != nullptr);  

  auto currPos = -1, oldPos = -1; 
  if(ant.old != nullptr) oldPos = ant.old->id;
  if(ant.current != nullptr) currPos = ant.current->id;

  vector<int> adjV = ant.current->adj;
  int oldColor = colors.at(currPos);
  int colorInAdj[maxNColor];
  for(auto i = 0; i < maxNColor; ++i){colorInAdj[i]=0;}

  for(auto &x : adjV){
    if(BB){
      if(colors.at(x) >= maxNColor){
	printf("color in vertex %d is %d, maxNColor %d\n",
	       x, colors.at(x), maxNColor);	  
	assert(false);
      }
    }
    colorInAdj[colors.at(x)]++;
  }

  //avoid choosing color from old POS
  if(getDIMACSBinaryEdgeSwap(oldPos,currPos)){
    colorInAdj[colors.at(oldPos)] = 2 * nVertices;
  }

  vector<int> leftOver;
  for(auto i = 0; i < maxNColor; ++i){
     //if no adj has this color, then push it to leftOver
    if(!colorInAdj[i]) leftOver.push_back(i);
  }

  if(leftOver.size()){
    auto color = leftOver.at(rand() % leftOver.size());
    colors.at(currPos) = color;     
    conflictsTable[currPos] = 0;
    
    if(BB)assert(getConflictOfVertex(pVertices[currPos], colors)==0);
  }
  else{
    int minCount=colorInAdj[0];
    for(auto i = 1; i < maxNColor; ++i){
      if(colorInAdj[i] < minCount){
	minCount = colorInAdj[i];
      }
    }
	
    if(BB){ for(auto i = 0; i < maxNColor; ++i){ assert(colorInAdj[i]>=minCount);}}
    vector<int> MinIndexVector;
    for(auto i = 0; i < maxNColor; ++i){
      if(colorInAdj[i]==minCount)MinIndexVector.push_back(i);
    }
    auto color = MinIndexVector.at(rand() % MinIndexVector.size());
    colors.at(currPos) = color;    
    conflictsTable[currPos] = getConflictOfVertex(pVertices[currPos], colors);
  }

  //update only the conflicts LOCALLY
  if(colors.at(currPos) != oldColor){//if changed the Color
	
    //decrement since the colors has changed
    for(auto &x : adjV){
      if(colors.at(x) == oldColor){
	conflictsTable[x]--;
      }
    }
	
    //if there's still conflict, and if they are the same color , incre
    if(conflictsTable[currPos]){//if still conflict  , TVN: not sure if this should be done , doesn't make too much sense
      for(auto &x : adjV){
	if(colors.at(x) == colors.at(currPos)){
	  conflictsTable[x]++;
	}
      }
    }//if still conflict
  }
}

//move the ants, very frequently used operation
int move(const int &curPos, const MoveOpts &move_method, 
	 const vector<int> &recentlyVisited,
	 const int conflictTable[],
	 Vertex **pVertices,
	 const int &nVertices){
  
  int movePos=-1;
  vector<int> adjV = pVertices[curPos]->adj;
  vector<int>chosenOne;
  bool considerVertex[nVertices];

  if(!adjV.empty()){//if it has some adj vertices

    //these 3 loops are fast (so don't try to optimize)
    for(auto i=0; i < nVertices; ++i) considerVertex[i] = false;//reset
    for(auto &x : adjV) considerVertex[x] = true; 
    for(auto &x : recentlyVisited) considerVertex[x] = false;
	
    if(move_method == MoveOpts::Random){	
      for(auto &x : adjV){
	if(considerVertex[x]){
	  if(BB) assert(getDIMACSBinaryEdgeSwap(curPos, x));
	  chosenOne.push_back(x);
	}
      }
    }
    if(move_method == MoveOpts::MaxConflict){
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
  else{
    //if has no adj vertices, chose all ones not recently visited
    for(auto i = 0; i < nVertices; ++i) considerVertex[i] = true;
    for(auto &x : recentlyVisited) considerVertex[x] = false;
    for(auto i = 0; i < nVertices; ++i){
      if(considerVertex[i]) chosenOne.push_back(i);
    }
  }

  movePos=chosenOne.empty() ?
    adjV.at(rand() % adjV.size()) :
    chosenOne.at(rand() % chosenOne.size());  

  if(BB)assert(movePos!=-1);
  return movePos;
}

void rand_color(vector<int> &colors, const int &nColors){
  if(BB)assert(nColors>0);//if it's 0 then only 1 color
  for(auto i = 0; i < colors.size(); ++i){
    if(colors.at(i) >= nColors){
      colors.at(i) = rand() % nColors;
    }
  }
}

void antsOps(Sol &sol,
	     vector <int> &cur_colors,
	     Vertex **pVertices,
	     const int &nVertices){

  auto nAnts = (int)(nVertices * nAntsPercent);
  if(nAnts > 100) nAnts = 100;

  auto moveLimit = (nVertices > 100) ?
    20 + (int)(nVertices / nAnts) : (int) (moveLimitPercent * nVertices);
  
  auto rSizeLimit = (int)(moveLimit / rSizeLimitFactor);
  if(rSizeLimit < 1) rSizeLimit = 1;


  auto nCycles = nVertices * nCyclesFactor;
  if(nCycles > 4000) nCycles = 4000;
  
  auto breakCycles = int(nCycles / 2);
  if(breakCycles > 1500) breakCycles = 1500;      

  auto alphaNumColors = (int)(sol.nColors * ALPHA);
  auto bestCycle = 0;

  //add ants
  vector<Ant> ants(nAnts);
  for(auto i = 0; i < nAnts; ++i){
    auto ant = Ant();
    ant.id = i;
    ant.current = nullptr;
    ant.old = nullptr;
    ants.at(i) = ant;
  }
    
  if(BB){
    auto lg = *max_element(cur_colors.begin(), cur_colors.end());
    
#ifdef V
    printf("Q ALPHA is %d, percent %f, XG %d, largest Color Index %d\n",
	   alphaNumColors,(double)alphaNumColors/(double)sol.nColors,sol.nColors,lg);
#endif
    assert (lg != 0);
  }

  int conflictsTable[nVertices];
  updateConflictTable(cur_colors, conflictsTable, pVertices);  
  int totalConflicts = -1;

  vector<int> recentlyVisited;
  int changedCycle=0;
  int moveSoFar;

  for (auto iCycle = 0; iCycle < nCycles; ++iCycle) {
    for(auto &ant : ants){
      moveSoFar = 0;
      if(!recentlyVisited.empty())
	recentlyVisited.clear();

      auto pos = selectInitMove(MOVE_METHOD, conflictsTable, nVertices);
      ant.current = pVertices[pos];
      moveSoFar++;
      if(BB)assert(moveSoFar==1);
      color(ant, cur_colors, conflictsTable, alphaNumColors,
	    pVertices, nVertices);

      int movePos;
      while(moveSoFar < moveLimit){ 
	if(BB)assert(ant.current != nullptr); 
	ant.old = ant.current;
	for(int distI=0; distI < HOW_FAR; ++distI){
	  ///move to Max Conflict adj lastly
	  movePos = move(ant.current->id,
			 (distI==HOW_FAR-1)?MoveOpts::MaxConflict:MoveOpts::Random,
			 recentlyVisited, conflictsTable,
			 pVertices, nVertices);
	  ant.current=pVertices[movePos];

	  //printf("-> [%d] old %d, current %d\n",iAnts,ant.old->id,ant.current->id);

	  //if reach maximum then remove from the beginning, todo, use a diff structure
	  if(recentlyVisited.size()== rSizeLimit)
	    recentlyVisited.erase(recentlyVisited.begin());
	  recentlyVisited.push_back(ant.current->id);  

	}//HOW_FAR

	color(ant, cur_colors, conflictsTable, alphaNumColors,
	      pVertices, nVertices);
	moveSoFar++;

      }//while(moveSoFar<moveLimit)
    }//for ant loop

    totalConflicts = updateConflictTable(cur_colors, conflictsTable, pVertices);  
    if(!totalConflicts && alphaNumColors < sol.nColors){
      updateSol(sol, cur_colors, alphaNumColors, iCycle);
      
#ifdef V
      printf("*********** Achieve %d colors at cycle %d ************\n",
	     sol.nColors, iCycle);
#endif
      
      if(alphaNumColors > 1){
	alphaNumColors--;
#ifdef V
	printf("Cycle [last changed %d, current %d], "
	       "Q_CHANGE_CYCLE %d, decrease nColors from %d to %d\n",
	       changedCycle, iCycle, Q_CHANGE_CYCLE,
	       alphaNumColors+1, alphaNumColors);
#endif
	changedCycle = iCycle;
	rand_color(cur_colors, alphaNumColors);
	updateConflictTable(cur_colors, conflictsTable, pVertices);
	}
      else{
	fprintf(stderr,"E: attempting to have only 1 color, "
		"only possible if there's no edge .. haha\n");
	assert(false);
      }
    }//if no conflict and better than best result


    if (iCycle - changedCycle == Q_CHANGE_CYCLE){
      int tIncr = (int)(sol.nColors - alphaNumColors/4);

      if(tIncr==0 || alphaNumColors+tIncr >= sol.nColors){
	tIncr = 1;
      }
	  
      if(alphaNumColors + tIncr < sol.nColors){
	alphaNumColors += tIncr;
	changedCycle = iCycle;
#ifdef V
	printf("Cycle [last changed %d, current %d], "
	       "Q_CHANGE_CYCLE %d, increasing nColors from %d to %d\n",
	       changedCycle, iCycle, Q_CHANGE_CYCLE,
	       alphaNumColors-tIncr, alphaNumColors);
#endif
      }
    }

    if(iCycle - changedCycle == breakCycles){
#ifdef V
      printf("Break at cycle %d, last changed cycle %d, breakCycle %d\n",
	     iCycle, changedCycle, breakCycles);
#endif
      break;
    }
  }//cycle loop

  //clean up
  //for (auto &ant : ants) delete ant;
}
