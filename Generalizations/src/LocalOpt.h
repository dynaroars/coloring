
/*
  This functions "GREEDILY" assign colors to vertices.
  It will return a VALID colors set results, but most certainly not OPTIMAL
  TODO: do something like XRLF, only get 70% or so percent of the vertices colored
  using greedy method. the 30 % left = > random
*/



int chooseV(const int color[], const int inColor[],const int &OPTION){
  
  int chosenV=-1;

  if(OPTION == MAX_DEG){
    int maxTest=0,vAdjSize;
    for(int i = 0;i<nVertices;++i){
      if (color[i]==-1){//only choose UNcolored
	vAdjSize=pVertices[i]->adj.size();
	if (maxTest<vAdjSize){
	  maxTest=vAdjSize;  chosenV=i;
	}
      }
    }//for (i=0:nVertices)
  }

  if(OPTION==RANDOM_CHOOSE){
    vector<int>availV;
    for(int i = 0 ; i < nVertices ; ++i){
      if(color[i]==-1) availV.push_back(i);
    }
    chosenV=availV.empty()?-1:availV[genrand_int32()%availV.size()];
  }


  if(OPTION==MAX_COLOR){
    int maxC=-1;
    for(int i = 0;i<nVertices;++i){
      if (color[i]==-1){//only choose UNcolored
	if (maxC<inColor[i]){
	  maxC=inColor[i]; chosenV=i;
	}
      }
    }//for (i=0:nVertices)
  }
	
  return chosenV;
}

void getIG(int colorV[],const int &nVertices, const int &nEdges){


  int nUnColored=nVertices;//nothing is colored yet
  int currentColor=-1; //color Index starts with 0 

  
  //init
  for(int i = 0 ; i < nVertices ;++i){colorV[i]=-1;}


  vector<int>forbiddenSet;
  int maxInForbidden, chosenV, vAdjSize, lowerC, upperC;



  while(nUnColored> 0){

    chosenV=chooseV(colorV,NULL,IG_OPT);

    if(BB)assert(chosenV<nVertices);
    if(chosenV==-1){
      if(BB)printf("chosenV is -1 , break\n");
      break; //break from while loop, nothing else to choose 
    }



    forbiddenSet.clear();maxInForbidden=-1;	//reset
	
    vAdjSize=pVertices[chosenV]->adj.size();
    for(int i = 0 ; i < vAdjSize ; ++i){
      int aV = pVertices[chosenV]->adj.at(i);  //adjacent Vertex
	  
      if(colorV[aV]!=-1){//if this neighbored is colored


	if(BB)assert(chosenV!=aV);
	if(BB)assert(weightMatrix[aV][chosenV]==weightMatrix[aV][chosenV]);
		
		
	lowerC=colorV[aV]-weightMatrix[aV][chosenV] + 1;
	if(lowerC<0)lowerC=0;
		
	upperC=colorV[aV]+weightMatrix[aV][chosenV] - 1;
		
	/*		printf("lowerC %d , upperC %d, vertex %d , color %d, weight(%d,%d) is %d\n",
			lowerC+1,upperC+1,aV+1,colorV[aV]+1,aV+1,chosenV+1,weightMatrix[aV][chosenV]);
	*/
		  
		
	for(int j = lowerC ; j <= upperC ; ++j){
		  
	  bool pushBack = true;
	  for(int g=0;g<forbiddenSet.size();++g){
	    if(forbiddenSet.at(g)==j){pushBack=false;break;}
	  }
		  
	  if(pushBack){forbiddenSet.push_back(j);if(j>maxInForbidden)maxInForbidden=j;}
		  
	  if(BB){//making sure no duplicates in forbiddenSet
	    for(int g=0;g<forbiddenSet.size();++g){
	      for(int g1=0;g1<forbiddenSet.size();++g1){
		if(g!=g1){assert(forbiddenSet.at(g)!=forbiddenSet.at(g1));}
	      }
	    }
	  }
		  
	}//lowerC to upperC loop
		
      }//if neighboris colored
	  
    }//end for all neighbor
	



	

    //The following code block determines and assign the color to the chosen vertex (chosenV)
    if(forbiddenSet.empty()){
      currentColor=0;//no forbidden, assign the smallest, actually can be anything since no forbidden
    }
    else{
      /*
	Given a SORTED forbidden set, determines if the smallest number that is NOT in that forbidden set
	e.g,  forbiddenSet = "0 1 2 3 5 7 8 9 10 13 15" 
	the chosen color will be "4", since it's the smallest that's NOT in the forbidden set.
		
	If the forbiddenSet = " 0 1 2 3 4 5 6 7 8 9 10 " , then the next color is 10+1 , or 11.
      */
	  
      //printVector("forbidden",forbiddenSet);	
	  

      //since index start from 0,  so even if the largest in forbidden is 0, 
      //that still means cannot consider that one lement
      //this damn BUG costs me 4 HOURS !!!! 
      //the BUG itself is easy to fix but realizing the BUG was caused by this is not
      //the damn bool tempA[] was exceeded and modified values from bool isColored[] !!!
	  
      int aSize = maxInForbidden+1;
	  
      bool tempA[aSize];
      for(int g=0;g<aSize;++g){
	tempA[g]=true;//consider these
      }
	  
	  
      for(int g=0;g<forbiddenSet.size();++g){
	tempA[forbiddenSet.at(g)]=false ; // don't consider these
      }
	  

	
	
      int TColor = -1;
      for(int g=0;g<maxInForbidden;++g){
	if(tempA[g]){TColor=g ; break;}
      }
	  
      //determine current color
      currentColor=TColor>-1?TColor:maxInForbidden+1;
      //printf("currentColor %d\n",currentColor);
	  
    }
	
    //assigning color

    //	if(colorV[chosenV]==-1 || isColored[chosenV]==false){
    //	  printf("vertex %d already COLORED(%d), %d\n",chosenV+1,isColored[chosenV],colorV[chosenV]+1);
    //	};




    if(BB)assert(currentColor!=-1);
    colorV[chosenV]=currentColor;
    nUnColored--;

    //printf("%d, chosen %d colored %d\n",nUnColored,chosenV+1,colorV[chosenV]+1);
	


  }//while





  //if there's any uncolored, color it with the smallest color index
  for(int i = 0 ; i < nVertices ;++i){if(colorV[i]==-1) colorV[i]=0; }




  //get the "DISTINCT" colors, e.g., if the colors used are 1, 100, 
  //then the nColors is just 2 .
  
}

int IG(int bestConfigA[],const int &nVertices, const int &nEdges){
  
  int colorV[nVertices], t_maxK=-1;
  if(IG_OPT==MAX_DEG){
    getIG(colorV,nVertices,nEdges);
    t_maxK=getMaxK(colorV);    
    if(t_maxK<maxK){
      maxK=t_maxK;
      copyArray(colorV,bestConfigA);//copy to result vector
    }
	
  }
  if(IG_OPT==RANDOM_CHOOSE){
    for(int i = 0 ; i < 20; ++i){
      getIG(colorV,nVertices,nEdges);
      t_maxK=getMaxK(colorV);    
      if(BB)printf("iteration %d , maxK %d\n",i,t_maxK);
      if(t_maxK<maxK){
	maxK=t_maxK;
	copyArray(colorV,bestConfigA);//copy to result vector
	if(BB)printf("--- new maxK %d\n",maxK);
      }
    }//for i
  }

  

  baseScore_MaxK=maxK;

  addPheromoneScore(maxK,0,bestConfigA);//0 conflict

  //  for(int i = 0 ; i < nVertices ;++i){
  //	printVector("",pVertices[i]->phColor);
  //  }
  


  if(BB){int conflictsTable[nVertices];	assert(updateConflictTableW(colorV,conflictsTable)==0);}  


  int bestResult=getDistinctColors(bestConfigA);
  if(BB)printf("*** Ig distinct colors %d, maxK %d\n",bestResult,maxK);

  return bestResult;

}//end IG()



int localOpt(const int inColor[], int lColor[], const int &nVertices, const int &nEdges){

  for(int i = 0 ; i < nVertices ;++i){lColor[i]=-1;}
  int currentColor=-1; //color Index starts with 0 

  

  vector<int>forbiddenSet;
  int chosenV,lowerC,upperC,nUnColored=nVertices,maxInForbidden,vAdjSize;
  
  while(nUnColored> 0){
    chosenV=chooseV(lColor,inColor,MAX_COLOR);
    if(BB)assert(chosenV!=-1);

    forbiddenSet.clear();maxInForbidden=-1;	//reset
    vAdjSize=pVertices[chosenV]->adj.size();

    for(int i = 0 ; i < vAdjSize ; ++i){
      int aV = pVertices[chosenV]->adj.at(i);  //adjacent Vertex
	  
      if(lColor[aV]!=-1){//if this neighbored is colored
	lowerC = 0 ; //smallest avail

	int tmpD = lColor[aV]-weightMatrix[aV][chosenV] + 1; ;
	if (lowerC<tmpD) lowerC=tmpD;
		
	upperC=lColor[aV]+weightMatrix[aV][chosenV] - 1;
		
	/*		printf("lowerC %d , upperC %d, vertex %d , color %d, weight(%d,%d) is %d\n",
			lowerC+1,upperC+1,aV+1,colorV[aV]+1,aV+1,chosenV+1,weightMatrix[aV][chosenV]);
	*/
		  
		
	for(int j = lowerC ; j <= upperC ; ++j){
		  
	  bool pushBack = true;
	  for(int g=0;g<forbiddenSet.size();++g){
	    if(forbiddenSet.at(g)==j){pushBack=false;break;}
	  }
		  
	  if(pushBack){forbiddenSet.push_back(j);if(j>maxInForbidden)maxInForbidden=j;}
		  
	  if(BB){//making sure no duplicates in forbiddenSet
	    for(int g=0;g<forbiddenSet.size();++g){
	      for(int g1=0;g1<forbiddenSet.size();++g1){
		if(g!=g1){assert(forbiddenSet.at(g)!=forbiddenSet.at(g1));}
	      }
	    }
	  }
		  
	}//lowerC to upperC loop
		
      }//if neighboris colored
	  
    }//end for all neighbor



    if(forbiddenSet.empty()){
      currentColor=0;//no forbidden, assign the smallest, actually can be anything since no forbidden
    }
    else{
	  
      int aSize = maxInForbidden+1;
	  
      bool tempA[aSize];
      for(int g=0;g<aSize;++g){tempA[g]=true;}
      for(int g=0;g<forbiddenSet.size();++g){tempA[forbiddenSet.at(g)]=false ; }
	  
	  
      int TColor = -1;
      for(int g=0;g<maxInForbidden;++g){
	if(tempA[g]){TColor=g ; break;}
      }
	  
      currentColor=TColor>-1?TColor:maxInForbidden+1;
    }


    lColor[chosenV]=currentColor;
    nUnColored--;


  }//while
  

  //if there's any uncolored, color it with the smallest color index
  for(int i = 0 ; i < nVertices ;++i){
    if(lColor[i]==-1)  lColor[i]=0;

  }

  if(BB){int tempConflicts[nVertices];assert(updateConflictTableW(lColor,tempConflicts)==0); }
  
  return   getMaxK(lColor);
  
}

