/*movePos in both case can result in -1 if there's no conflict*/
int placeAnt(const int &MOVE_METHOD, const int conflictTable[],const int conflictTableAHEAD[],
	     const vector<int> &recentlyVisited,bool considerVertex[]){

  int movePos=-1;
  int maxConflictVal=0;
  vector<int>maxConflictV;

  
  for (int i=0 ; i < nVertices; ++i)considerVertex[i]=true;//reset
  for(int i=0 ; i < recentlyVisited.size(); ++i) {
    considerVertex[recentlyVisited.at(i)]=false;//don't considered ones recently visited
  }
  
  for(int i=0;i<nVertices;++i){
    if(considerVertex[i]){
      if(conflictTable[i]>maxConflictVal){
	maxConflictVal=conflictTable[i];
	maxConflictV.push_back(i);   //01302006  is this a bug here ? basically I just push back everything that is > currentMax
      }
    }
  }

  if(maxConflictV.size()>1){//several , need to break the tie
    maxConflictVal=0;
    vector<int>maxConflictVAHEAD;
    int aV;
    for(int i = 0 ; i < maxConflictV.size();++i){
      aV=maxConflictV.at(i);
      if(conflictTableAHEAD[aV]>maxConflictVal){
	maxConflictVal=conflictTableAHEAD[aV];
	maxConflictVAHEAD.push_back(aV);
      }
    }

    if(maxConflictVAHEAD.size()){movePos=maxConflictVAHEAD.at(genrand_int32()%maxConflictVAHEAD.size());}
    else{if(BB)assert(1==0);}

  }
  if(maxConflictV.size()==1){
    movePos=maxConflictV.at(0);
  }

  //if movePos==-1 then that means there's no fonclit
  if(BB){if(movePos==-1) assert(maxConflictV.empty());}

  return movePos;
}




//move the ants, very frequently used operation
int  moveToNeighbor(const int &curPos,const int &MOVE_METHOD, 
		    const vector<int> &recentlyVisited,bool considerVertex[],
		    const int conflictTable[]){
  
  int movePos=-1;
  vector<int> adjV=pVertices[curPos]->adj;
  vector<int>chosenOne;

  int tempV;


  if(!adjV.empty()){//if it has some adj vertices
  
    for(int i=0;i < nVertices ;++i)considerVertex[i]=false;//reset	
    for (int i=0 ; i < adjV.size(); ++i)considerVertex[adjV.at(i)]=true; //consider these adj's
    for(int i=0 ; i < recentlyVisited.size(); ++i) {
      if(considerVertex[recentlyVisited.at(i)]){
	considerVertex[recentlyVisited.at(i)]=false;//don't considered ones recently visited
      }
    }
	
	
	
    if(MOVE_METHOD==RANDOM_CHOOSE){	
      for(int i=0 ; i < adjV.size();++i){
	tempV=adjV.at(i);
	if(considerVertex[tempV]&&conflictTable[tempV]){
	  if(BB) assert(edgeMatrix.at(curPos).at(tempV));

	  chosenOne.push_back(tempV);
	}
      }
      //if(BB){if(chosenOne.empty())printf("chosenOne is empty");}
    }
	
    if(MOVE_METHOD==MAX_DEG){
	  
      int maxConflict=0;
      int vid;
	  
      for(int i = 0 ; i < adjV.size(); ++i){
	vid=adjV.at(i);
	if(considerVertex[vid]){
	  if(conflictTable[vid]>maxConflict){
	    maxConflict=conflictTable[vid];
	  }
	}
      }//end for
	  
	  
      if(maxConflict){//if there are conflicts
		
	if(BB)assert(chosenOne.empty());

	for(int i = 0 ; i < adjV.size(); ++i){
	  vid=adjV.at(i);
	  if(considerVertex[vid]){
	    if(conflictTable[vid]==maxConflict)chosenOne.push_back(vid);
	  }
	}
	if(BB)assert(!chosenOne.empty()); 
      }//if maxConflict != 0 
	  
    }

  }
  else{//if this vertex has no adj vertices

    for(int i=0 ; i < nVertices ;++i)considerVertex[i]=true;//reset	
    for(int i=0 ; i < recentlyVisited.size(); ++i) {
      considerVertex[recentlyVisited.at(i)]=false;//don't considered ones recently visited
    }
	
    for(int i=0 ; i < nVertices ;++i){
      if(considerVertex[i])chosenOne.push_back(i);
    }
	
  }
	


  if(chosenOne.empty()){

    for(int i=0 ; i < nVertices ;++i)considerVertex[i]=true;//reset	
    for(int i=0 ; i < recentlyVisited.size(); ++i) {
      considerVertex[recentlyVisited.at(i)]=false;//don't considered ones recently visited
    }

    for(int i = 0 ; i < nVertices ; ++i){
      if(considerVertex[i]&&conflictTable[i]){
	chosenOne.push_back(i);
      }
    }
  }
  



  movePos=chosenOne.empty()?-1:chosenOne.at(genrand_int32()%chosenOne.size());  




  if(BB){if(movePos!=-1)assert(conflictTable[movePos]>0);}


  return movePos;
}



//color the vertex
void colorTheVertex(const int &currPos,const int &maxNColor,
		    int tempConfigA[], int conflictsTable[], int conflictsTableAHEAD[]){
  

  vector<int> adjV=pVertices[currPos]->adj;

  
  int oldColor=tempConfigA[currPos];
  
  int freq[maxNColor+1];for(int i = 0 ; i <= maxNColor; ++i){freq[i]=0;}
  

  int lowerC,upperC;
	
  int aV;
  for(int i = 0 ; i < adjV.size() ; ++i){
    aV = adjV.at(i);  //adjacent Vertex
    lowerC = tempConfigA[aV]-weightMatrix[aV][currPos]+1;if (lowerC<0) lowerC=0;
    upperC = tempConfigA[aV]+weightMatrix[aV][currPos]-1;if(upperC>maxNColor)upperC=maxNColor;

    for(int j = lowerC ; j <= upperC ; ++j)freq[j]++;

  }//end for all neighbor


  vector<int>forbiddenSet;
  for(int g = 0 ; g <= maxNColor ; ++ g ){
    if(freq[g]!=0)forbiddenSet.push_back(g);
  }
  






  if(forbiddenSet.empty()){
    tempConfigA[currPos]=oldColor==0?1:0;//no forbidden, assign the smallest, actually can be anything since no forbidden
  }
  else{

    int maxInforbiddenSet =  forbiddenSet.at(forbiddenSet.size()-1);  if(BB)assert(maxInforbiddenSet<=maxNColor);
    int aSize =maxInforbiddenSet+1;
	  
    bool tempA[aSize];
    for(int g=0;g<aSize;++g) tempA[g]=true;//consider these

    for(int g=0;g<forbiddenSet.size();++g){
      tempA[forbiddenSet.at(g)]=false ; // don't consider these
    }
	
    //assert(tempA[maxInforbiddenSet]==false);

	
    vector<int>diffSet;
    for(int g=0;g<maxInforbiddenSet;++g){
      if(tempA[g] && g!=oldColor){diffSet.push_back(g);}
    }


    if(diffSet.empty()){//if all colors are used, the next valid color is maxInForbidden + 1
      if(maxInforbiddenSet < maxNColor){
	tempConfigA[currPos] = maxInforbiddenSet+1;
      }
      else{		//todo 2 ,  get the one w/ least freq (i.g., the # overlappsed conflict range is the least);

		
	int leastFreqVal = -1;
	for(int g=0;g<=maxNColor;++g){
	  if(g!=oldColor){
	    if(leastFreqVal==-1){//first time
	      leastFreqVal=freq[g];
	    }
	    else{
	      if(freq[g]<leastFreqVal){
		leastFreqVal=freq[g];
	      }
	    }
	  }
	}

	vector<int>leastFreqIndices;
	for(int g=0;g<=maxNColor;++g){
	  if(g!=oldColor){
	    if(freq[g]==leastFreqVal)leastFreqIndices.push_back(g);
	  }
		  
	}
		
	if(BB)assert(leastFreqIndices.size());
	tempConfigA[currPos]=leastFreqIndices.at(genrand_int32()%leastFreqIndices.size());
	if(BB)assert(tempConfigA[currPos]!=oldColor);

      }
	  
    }
    else{//if diffSet is not empty, i.e., some intervals,  then choose the middle of the largest avail color interval.

      /*Ex diffset   2 3 4 5     9 10 11       12    15 16   then will choose 3 or 4 most likely*/


      vector<vector<int> >intervalSet;
      for(int i =0 ; i < diffSet.size();++i){
	int dsVal = diffSet.at(i);
	if(i==0){
	  vector<int>tempSet;
	  tempSet.push_back(dsVal);
	  intervalSet.push_back(tempSet);
	}
	else{
	  if(BB)assert(dsVal-diffSet.at(i-1)!=0);

	  if(dsVal-diffSet.at(i-1)>1) {//not continous
	    vector<int>tempSet;
	    tempSet.push_back(dsVal);
	    intervalSet.push_back(tempSet);
	  }
	  else{//=1 , continuous
	    intervalSet.at(intervalSet.size()-1).push_back(dsVal);
	  }
	}
      }
		  
      //	  printVector("diffSet",diffSet);

      //	  for(int i=0; i < intervalSet.size();++i)printVector("intervalSet",intervalSet.at(i));


      int largestISSize=0; //largest Interval Set Size
      for(int i=0; i < intervalSet.size();++i){
	if(intervalSet.at(i).size()>largestISSize){
	  largestISSize=intervalSet.at(i).size();
	}
      }
	  
	  

      vector<int>largestISIndices;
      for(int i=0;i<intervalSet.size();++i){
	if(intervalSet.at(i).size()==largestISSize){
	  largestISIndices.push_back(i);
	}
      }

      if(BB)assert(largestISIndices.size());
	  

      /*e.g., the set is 4 5 7 6, size 4 then the index is 2-1=2 */
      /*e.g., the set is 3 4 5 7 6, size 5 then the index is ceil(2.5)-1=2 */
      int midIndex=ceil(largestISSize/2.0)-1;  


      tempConfigA[currPos]=oldColor;//temp assign
      while(tempConfigA[currPos]==oldColor){
	int randIndex=largestISIndices.at(genrand_int32()%largestISIndices.size());
	tempConfigA[currPos]=intervalSet.at(randIndex).at(midIndex);
      }
	  

    }//if diffSet NOT empty
	
  }//if forbiddenSet NOT empty
  
  
  if(BB)assert(tempConfigA[currPos]!=oldColor);



  //update conflicts
  conflictsTable[currPos]=getConflictOfVertexW(currPos,tempConfigA);
  conflictsTableAHEAD[currPos]=getConflictOfVertexAHEAD(currPos,tempConfigA,AHEAD_DIST);

  
  int newColor=tempConfigA[currPos];
  
  int neighborColor, weightD;
  for(int i= 0 ; i<adjV.size();++i){
    neighborColor = tempConfigA[adjV.at(i)];
    weightD=weightMatrix.at(currPos).at(adjV.at(i));
	
    //if previous there were conflicts 
    if (abs(oldColor-neighborColor)<weightD){
      //and now there's no 
      if(abs(newColor-neighborColor)>=weightD){
	conflictsTable[adjV.at(i)]--;
	conflictsTableAHEAD[adjV.at(i)]--;
	//assert(conflictsTable[adjV.at(i)]>=0);
      }
    }
    else{//if prevous there where NO conflicts and now there is 
      if(abs(newColor-neighborColor)<weightD){
	conflictsTable[adjV.at(i)]++;
	conflictsTableAHEAD[adjV.at(i)]++;
      }
    }
	
    //assert(getConflictOfVertexW(adjV.at(i),tempConfigA)==conflictsTable[adjV.at(i)]);
	
  }//adjV loop

}//enb colorTheVertex






void AntsOps(int tempConfigA[], int bestConfigA[]){


  //printVector("bes colorA: ", bestColorResult);
  //  printArray("tempConfigA:", tempConfigA,nVertices);

  
  int alphaMaxK= (int)(maxK * ALPHA); //- qSmallest;

  if(alphaMaxK <= maxWeight){
    alphaMaxK = maxWeight + 1; 
    //could happen w/ geom20b
    //	printf("alphaMax %d is <= than maxWeight %d\n",alphaMaxK,maxWeight);
    //	assert(1==0);
  }


  if(BB)printf("alphaMaxK(q) is %d\n",alphaMaxK+1);


  


  int conflictsTable[nVertices],conflictsTableAHEAD[nVertices];
  int totalConflicts = updateConflictTableW(tempConfigA,conflictsTable);  
  int totalConflictsAHEAD = updateConflictTableAHEAD(tempConfigA,conflictsTableAHEAD,AHEAD_DIST); 

  if(BB && AHEAD_DIST==1)for(int i=0;i<nVertices;++i)assert(conflictsTable[i]==conflictsTableAHEAD[i]);

  if(BB)printf("Total Conflict %d\n",totalConflicts);

  int changedCycle=0;

  vector<int> recentlyVisited ;
  bool considerVertices[nVertices];
  
  //store the maxK achieved by Ants operations
  int AntsMaxKAchieved=0,LOptMaxKAchieved=0; 
  int ltempConfigA[nVertices];//temp stores config returns by local Opt  



  if(BB){
    printf("NANTS %d, NCYCLES %d, MOVE_LIMIT %d, nDists %d, R_SIZE_LIMIT %d\n",
	   nAnts, nCycles, nMoves, nDists, rSizeLimit);
  }


  int currPos=-1; 

  bool JUSTBREAK = false;

  for (int iCycles=0 ; iCycles < nCycles; ++iCycles) {

    for(int iAnts=0; iAnts<nAnts; ++iAnts){
	  
      if(!recentlyVisited.empty())recentlyVisited.clear();//resest
	  
      currPos=placeAnt(MAX_DEG,conflictsTable,conflictsTableAHEAD,recentlyVisited,considerVertices);
	  
      if(currPos==-1) {break;}//no conflicts among vertices 
	  

      //otherwise, (there is at least a conflict)
	  
      if(BB)assert(getConflictOfVertexW(currPos,tempConfigA)!=0);

      colorTheVertex(currPos,alphaMaxK,tempConfigA,conflictsTable,conflictsTableAHEAD);
	  

      for (int iMoves=0;iMoves<nMoves;++iMoves){ 
		
	for(int iDists=0;iDists<nDists;++iDists){

		  
	  if(iDists==nDists-1){
	    if(isConflictAmongAdj(currPos,tempConfigA)){
			  
	      currPos==moveToNeighbor(currPos,MAX_DEG,
				      recentlyVisited,considerVertices,conflictsTable);
			  
	    }
	    else{
	      currPos=placeAnt(MAX_DEG,conflictsTable,conflictsTableAHEAD,recentlyVisited,considerVertices);
			  
	    }
			
	  }
	  else{//if iDist is not the last
	    currPos=moveToNeighbor(currPos, /*seems like trouble is here*/
				   RANDOM_CHOOSE,
				   recentlyVisited,considerVertices,conflictsTable);
	  }
		  

	  if(currPos==-1){JUSTBREAK=true;break;}

	  if(recentlyVisited.size()==rSizeLimit)recentlyVisited.erase(recentlyVisited.begin());
	  recentlyVisited.push_back(currPos);  

	}//for iDists loop

	if(JUSTBREAK){
	  if(BB)assert(currPos==-1);
	  JUSTBREAK=false;//reset
	  break;
	}
	colorTheVertex(currPos,alphaMaxK,tempConfigA,conflictsTable,conflictsTableAHEAD);
	//		if(iCycles-changedCycle>joltCycles){
	//		  printf("attempting to recolor %d, alphaMaxK %d\n",currPos,alphaMaxK);
	//		}

      }//for iMoves loop
    }//for ant loop



    totalConflicts=	updateConflictTableW(tempConfigA,conflictsTable);  
    totalConflictsAHEAD=updateConflictTableAHEAD(tempConfigA,conflictsTableAHEAD,AHEAD_DIST);  

    if(BB && AHEAD_DIST==1){for(int i=0;i<nVertices;++i)assert(conflictsTable[i]==conflictsTableAHEAD[i]);}

    if(!totalConflicts){ 

      if(BB) printf("\n\nCycle %d, Total conflict is 0\n",iCycles);

      AntsMaxKAchieved=getMaxK(tempConfigA);

      if(BB)assert(AntsMaxKAchieved<maxK);

      //localOpt
      LOptMaxKAchieved=localOpt(tempConfigA,ltempConfigA,nVertices,nEdges);
	  
      if(LOptMaxKAchieved<AntsMaxKAchieved){
	maxK=LOptMaxKAchieved;
	copyArray(ltempConfigA,tempConfigA);
	//		printf("LOptMaxKAchieved %d < AntsMaxKAchieved %d at cycle %d\n",LOptMaxKAchieved,AntsMaxKAchieved,iCycles);
      }
      else{
	maxK=AntsMaxKAchieved;
      }

      //assert(alphaMaxK>=maxK);
	  
      alphaMaxK=maxK;
      copyArray(tempConfigA,bestConfigA);
      bestCycle=iCycles;

#ifdef V
      printf("*********** Achieve maxK %d at cycle %d ************\n",maxK+1,iCycles);
	  
      if(maxK!=getMaxK(bestConfigA)){
	printf("maxK %d, getMaxK %d\n",maxK,getMaxK(bestConfigA));
	if(BB)assert(maxK==getMaxK(bestConfigA));
      }
	  
      if(BB)assert (updateConflictTableW(bestConfigA,conflictsTable)==0);
#endif
	  
      if(alphaMaxK>1){
	alphaMaxK--;
#ifdef V
	printf("Cycle [last changed %d, current %d], Q_CHANGE_CYCLE %d, decrease nColors from %d to %d\n",
	       changedCycle,iCycles,Q_CHANGE_CYCLE,alphaMaxK+1+1,alphaMaxK+1);
#endif
	changedCycle = iCycles;
	reColorMoreThanQ(alphaMaxK,tempConfigA,conflictsTable,conflictsTableAHEAD);
      }
      else{
	fprintf(stderr,"ERROR: attempting to have only 1 color, only possible if there's no edge .. haha\n");
	if(BB)assert(1==0);
      }
	  
    }//if no conflict
	
	
    //if no improvement for a long time
    if (iCycles-changedCycle == Q_CHANGE_CYCLE){
      int tIncr = (int)(maxK-alphaMaxK/4);

      if(tIncr==0||alphaMaxK+tIncr>maxK){tIncr=1;}
	  
      if(alphaMaxK+tIncr<maxK){
	alphaMaxK+=tIncr;
	changedCycle=iCycles;
	if(BB)printf("Cycle [last changed %d, current %d], Q_CHANGE_CYCLE %d, increasing nColors from %d to %d\n",
		     changedCycle,iCycles,Q_CHANGE_CYCLE,(alphaMaxK-tIncr)+1,alphaMaxK+1);
      }
    }

    //if no improvement for a BREAK_CYCLES of time
    if(iCycles-changedCycle==breakCycles){
      if(BB)printf("Break at cycle %d, last changed cycle %d, breakCycle %d\n",iCycles,changedCycle,breakCycles);
      break;
    }
	
	
  }//cycle loop

}
