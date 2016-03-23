//add Pheromone to vertices

/*CHECK TO SEE WHEN I REACH the point of no better improvement, then 
  how much 'difference' in the config will each cycle give ?  does it 
  even give a diff config at all  ?*/
void addScore(const int &confMaxK,const double &score,const int config[]){
  if(BB){assert(score!=0.0);  int testest = getMaxK(config);  assert(testest==confMaxK); 
    int tempConflicts[nVertices];assert(updateConflictTableW(config,tempConflicts)==0);  }

  int bCounter[nVertices];for(int i = 0 ; i < nVertices ;++i)bCounter[i]=0; for(int i = 0 ; i < nVertices ; ++i){ vector<double>&phC = pVertices[i]->phColor;	for(int j = 0 ; j < phC.size();++j){ if(phC.at(j))bCounter[i]++;	}  }
  

  //printf("maxK %d, score %f\n",confMaxK,score);


  //example if phC's color range is 0 ... c, then the size is c+1
  //therefore if input confMaxK =k then the size of the is exactly k+1,
  int vSize=confMaxK+1;


  for(int i = 0 ; i < nVertices ; ++i){ 
    vector<double>&phC = pVertices[i]->phColor;

    assert(config[i]<=confMaxK);

    while(phC.size()<config[i]+1){phC.push_back(0.0);}	
    while(phC.size()>vSize){ phC.pop_back();}

    phC[config[i]]+=score;	//assign the score

  }
  //  printf("testing, vertex[0]phColor size %d\n",pVertices[0]->phColor.size());
  //  printf("testing, vertex[0] has at color %d has score %f",45,pVertices[0]->phColor[45]);
  //printf("done\n");


  int aCounter[nVertices];for(int i = 0 ; i < nVertices ;++i)aCounter[i]=0;
  for(int i = 0 ; i < nVertices ; ++i){ 
    vector<double>&phC = pVertices[i]->phColor;
    for(int j = 0 ; j < phC.size();++j){
      if(phC.at(j))aCounter[i]++;
    }
  }

  if(BB){
    for(int i = 0 ; i < nVertices ; ++i){
      //WHY THE HELL DOESN't THIS ABBORT ?  
      //it should as in the case beforeCounter[i] has a confMaxK+1 color
      //it will be deleted
      //UPDATE: well it can but it rarely happens
      //assert(bCounter[i]<=aCounter[i]);
      //printf("v %d, colors avail before %d, after %d\n",i,bCounter[i],aCounter[i]);
    }
  }

  //update the Sum
  

  for(int i = 0 ; i < nVertices ; ++i){

    double &sumScore=pVertices[i]->phColorSumScore;
    //printf("vertex %d, before score %f, ",i,pVertices[i]->phColorSumScore);
    sumScore=0;
    vector<double>&phC = pVertices[i]->phColor;
    for(int j = 0 ; j < phC.size();++j){
      if(phC.at(j))sumScore+=phC.at(j);
    }
    //	printf("after score %f\n",pVertices[i]->phColorSumScore);
  }


}

//NOTE IMPORTANT:  have to check to make sure not ADDING similar 
//configs !!!
void addPheromoneScore(const int &confMaxK, const int &confTotalConflicts,
		       const int config[]){

  //  if(!configExist(config)){
  double scoreK=(double)baseScore_MaxK/confMaxK * scoreKFactor;
  double scoreC=(double)(1.0/(1.0+confTotalConflicts)) * scoreConflictFactor;
  double configScore=scoreK+scoreC ;
  addScore(confMaxK,configScore,config);//remember this config
  //  }
}





