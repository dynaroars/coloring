//***** Common Utilities *****//
void randInit(const int &argc,char *argv[]){if(argc==4){seed_t=atoi(argv[3]);}else{seed_t=time(NULL);}init_genrand(seed_t);}

void printSeed(const char *c, const int &OPTION){ 
  if(OPTION==TOCHECK)printf("%s : seed %u, %d %d %d\n",c,seed_t,genrand_int32(),genrand_int32(),genrand_int32());
  else	printf("%s : seed %u\n",c,seed_t);
}

template<class T>
void printArray(const char *c, const T a[],const int &aSize){
  printf("%s: ",c);for(int i = 0 ; i < aSize ; ++i)	cout << a[i] << " "; printf("\n");
}

template<class T>
void printVector(const char *c, const vector<T>&v){
  if(c!=""){
    printf("%s: ",c);
  }
  for(int i = 0 ; i < v.size() ; ++i)cout << v.at(i) << " "; printf("\n");
}


template<class T>
void print2DVector(const char *c, const vector<vector<T> >&v){
  printf("%s:\n",c);
  for(int i = 0 ; i < v.size() ; ++i){
    printVector("",v.at(i));
  }
}

void copyArray(const int A[],int B[]){for(int i = 0 ; i < nVertices ; ++i)B[i]=A[i];}


int getMaxK(const int a[]){ 
  int mK_index=0; for(int i = 1 ; i < nVertices ;++i)if(a[i]>a[mK_index])mK_index=i; return a[mK_index];
}


//return the color (int) proportionally based on its avail colors scores

int chooseColor_RouletteWheel(const vector<double>&score, const double &tSum){
  
  int theColor = -1;
  
  int tneedle = genrand_int32()%int(tSum*100);
  double needle = (float)(tneedle/100.00);

  //  printf("tSum=%f, need= %f\n",tSum,needle);

  for (int i = 0 ; i < score.size() ;++i){
    if (needle-score.at(i) > 0.0){
      needle-=score.at(i);
    }
    else{
      theColor=i;
      break;
    }
  }

  //  printf("chosen color is %d\n",theColor);

  assert(theColor!=-1);
  return theColor;

}


bool isConfsSimilar(const vector<int>&confA, const vector<int>&confB){
  if(confA.size()!=confB.size())return false;

  for(int i = 0 ;i < confB.size();++i){
    if(confA.at(i)!=confB.at(i))return false;
  }

  return true;
}

bool isConfsSimilar(const int confA[], const vector<int>&confB){
  for(int i = 0 ;i < confB.size();++i){
    if(confA[i]!=confB.at(i))return false;
  }
  return true;
}


//vector version
bool isConfExist(const vector<int>&conf ,  const vector<vector<int> >&allConfs){
  if(allConfs.empty()){
    return false;
  }
  else{
    for(int i = 0 ; i < allConfs.size();++i){
      if(isConfsSimilar(conf,allConfs.at(i)))return true;
    }
    return false;
  }
}


//array version
bool isConfExist(const int conf[] ,  const vector<vector<int> >&allConfs){
  if(allConfs.empty()){
    return false;
  }
  else{
    for(int i = 0 ; i < allConfs.size();++i){
      if(isConfsSimilar(conf,allConfs.at(i)))return true;
    }
    return false;
  }
}


void appendToAllConfs(const int conf[], const int &confSize ,vector<vector<int> >&allConfs){
  if(!isConfExist(conf,allConfs)){
    vector<int>newConf;
    for(int i = 0 ; i < confSize ; ++i){
      newConf.push_back(conf[i]);
    }
    allConfs.push_back(newConf);
  }
}

template<class T>
void getNonZeroIndices(const vector<T> &inputV, vector<int> &outputV){
  for(int i = 0 ; i < inputV.size(); ++i){
    if(inputV.at(i))outputV.push_back(i);
  }
}
