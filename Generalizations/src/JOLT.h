/*shuffling*/
void jolt(int tempConfigA[],const int &JOLT_PERCENT, const int &alphaMaxK){
  for(int i =0; i < nVertices ; ++i){
    if(genrand_int32()%100 <= JOLT_PERCENT){
      tempConfigA[i]=genrand_int32()%alphaMaxK;
    }
  }
}
