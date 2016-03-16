int getEdgesInClique(const int &cliqueSize){ return (cliqueSize*(cliqueSize-1))/2;}


int setupClique(const vector<int> &subSet,const int &theWeight){
  //  printVector("subset member",subSet);
  int tempCount=0;
  int nV1,nV2;
  //a clique among the subset
  for(int i = 0 ; i < subSet.size(); ++i){
	nV1=subSet.at(i);
	for(int j = 0 ; j <i ;++j){
	  nV2=subSet.at(j);
	  assert(nV1!=nV2);
	  edgeMatrix.at(nV1).at(nV2)=true;tempCount++;
	  edgeMatrix.at(nV2).at(nV1)=true;
	  weightMatrix.at(nV1).at(nV2)=theWeight;
	  weightMatrix.at(nV2).at(nV1)=theWeight;
	}
  }
  int shouldbe=getEdgesInClique(subSet.size());
  if(tempCount!=shouldbe)printf("size %d, tempCount %d, shouldbe %d\n",subSet.size(),tempCount,shouldbe);
  assert(tempCount==shouldbe);

  return shouldbe;

}

int merge2Cliques(const vector<int> &subSet1, const int &subSet1Weight, 
				  const vector<int> &subSet2,const int &subSet2Weight, 
				  const int &connectingWeight){

  int tempCount1=setupClique(subSet1,subSet1Weight);
  int tempCount2=setupClique(subSet2,subSet1Weight);
  int tempCount3=0;

  //connect them together
  int nV1, nV2;
  for(int i = 0 ; i < subSet1.size(); ++i){
	nV1=subSet1.at(i);
	for(int j = 0 ; j < subSet2.size(); ++j){
	  nV2=subSet2.at(j); assert(nV1 != nV2);
	  tempCount3++;
	  if(i==0&&j==0){
		assert(edgeMatrix.at(nV1).at(nV2)==true);
		assert(weightMatrix.at(nV1).at(nV2)==connectingWeight);
	  }
	  else{
		assert(edgeMatrix.at(nV1).at(nV2)==false);
		edgeMatrix.at(nV1).at(nV2)=true;
		edgeMatrix.at(nV2).at(nV1)=true;
		weightMatrix.at(nV1).at(nV2)=connectingWeight;
		weightMatrix.at(nV2).at(nV1)=connectingWeight;
	  }

	}
  }
  
  //
  int shouldbe=getEdgesInClique(subSet1.size()+subSet2.size());
  
  //  assert(tempCount==shouldbe);
  return shouldbe;
}

void initSuperSet(int &nVertices, int &nEdges){

  vector<vector<int> >superSet;

  int currIndex=nVertices-1;  //since index starts at 0
  for(int i = 0 ; i < ku.size() ; ++i){
	vector<int>tempSet;
	tempSet.push_back(i);//itself
	for(int j = 0 ; j < ku.at(i)-1 ; ++j){//its member (excluding itself hence -1)
	  currIndex+=1;
	  tempSet.push_back(currIndex);
	}

	superSet.push_back(tempSet);
	if(BB){
	  printf("set %d ",i);
	  printVector("members ",superSet.at(i));
	}
  }


  //assertion
  assert(superSet.size()==nVertices&&superSet.size()==ku.size());

  int new_nVertices=0;
  for(int i = 0 ; i < superSet.size() ;++i){
	assert(superSet.at(i).size()==ku.at(i));
	new_nVertices+=ku.at(i);
  }

  assert(new_nVertices==currIndex+1);

  
  //printf("old nVertices %d, new nVertices %d\n",nVertices,new_nVertices);
  int extraVertices=new_nVertices-nVertices;

  //
  //  for(int i = 0 ; i < superSet.size() ;++i){
  //	printVector("subset",superSet.at(i));
  //  }
  
  
  //resetup edgeMatrix and weightMatrix, basically recreated a new_nVertices x new_nVertices size matrix
  //while reserving old information

  //for the old, add more extraVertices (to the row)
  for(int i = 0 ; i < edgeMatrix.size(); ++i){
	for(int j = 0 ; j < extraVertices ; ++j){
	  edgeMatrix.at(i).push_back(false);//default val
	  weightMatrix.at(i).push_back(-10);//default
	}
  }

  //for the new , add extraVertices (to the col)
  for(int i = 0 ; i < extraVertices ; ++i){
	vector<bool>tempEMatrix(new_nVertices,false);
	edgeMatrix.push_back(tempEMatrix);
	vector<bool>tempMatrix(new_nVertices,false);

	vector<int>tempWMatrix(new_nVertices,-10);
	weightMatrix.push_back(tempWMatrix);
  }

  //assertion check
  for(int i = 0 ; i < edgeMatrix.size(); ++i){
	assert(edgeMatrix.at(i).size()==new_nVertices);
  }
  assert(edgeMatrix.size()==new_nVertices);
  //end assertion check

  
  int tempCount=0;
  for(int i = 0 ; i < superSet.size();++i){
	setupClique(superSet.at(i),weightMatrix.at(i).at(i));
	int ttt=getEdgesInClique(superSet.at(i).size());
	tempCount+=ttt;
	//	printf("n %d %d cliqueSize %d\n",i+1,superSet.at(i).size(),ttt);
  }

  //printf("new edges (NOT YET MERGE) %d\n",tempCount);


  
  int tempCount2=0;

  for (int i = 0 ;  i < nVertices ; ++i){
	for (int j = 0 ; j < i ;++j){
	  if (edgeMatrix.at(i).at(j)){
		int newE = merge2Cliques(superSet.at(i),weightMatrix.at(i).at(i),
								 superSet.at(j),weightMatrix.at(j).at(j),
								 weightMatrix.at(i).at(j));
		if(BB)printf("attempting to merge set %d (size %d, w %d )and set %d (size %d, w %d), connectW %d, #of new edges %d\n",
			   i,superSet.at(i).size(),weightMatrix.at(i).at(i),
			   j,superSet.at(j).size(), weightMatrix.at(j).at(j),
			   weightMatrix.at(i).at(j),
			   newE);
		tempCount2+=newE;
	  }
	}
  }
  

  if(BB)printf("new edges (AFTER MERGE) %d\n",tempCount2);

  
  nVertices=new_nVertices;



}

