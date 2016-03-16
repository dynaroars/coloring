well have to take care of the weight thing 
weight neighbor  -> all cliques of subset (i) to all cliques of subjset (j) have the weights of weight(i,j)
  weight self  -> all cliques of subset(i) have the weight of weight(i,i) among themselves ... have to check this 
is there weight(i,i) ???


const double nAntsPercent = 1 ;   
const int nCyclesFactor = 6 ; //nCycles = this * nVertices,
const double nBreakCycleFactor = 1;
const double nMovesPercent=.33;
const double nJoltPercent=.25;
const double rSizeLimitFactor=2.0;

nCycles=nVertices*nCyclesFactor;if(nCycles>4000)nCycles=4000;//if(nCycles<1000)nCycles=1000;
nAnts=(int)(nVertices*nAntsPercent);if(nAnts>100)nAnts=100;//if(nAnts<60)nAnts=60;
joltCycles=(int)(nCycles*nJoltPercent); if(joltCycles>300)joltCycles=300;
breakCycles=int(nCycles/nBreakCycleFactor); if(breakCycles>1500)breakCycles=1500;
nMoves=(int)nVertices*nMovesPercent;
rSizeLimit = (int)(nMoves/rSizeLimitFactor);if(rSizeLimit<1)rSizeLimit=1;


for (int iCycles=0 ; iCycles < nCycles; ++iCycles) {

  for(int iAnts=0; iAnts<nAnts; ++iAnts){
	  
	move to vertex where there is a MAX _# of conflicts(_if cannot find then _break, no conflict)
	color that vertex
	
	
	for (int iMoves=0;iMoves<nMoves;++iMoves){ 
		
	  for(int iDists=0;iDists<nDists;++iDists){
		
		  
		if(iDists==nDists-1){
		  move to vertex among neighbors that has the MAX conflicts,  if_cannot_do that 
		  then find any vertex that has MAX Conflicts
		  
		}
		else{
		  move to a RANDOM vertex among neighbors
		}
		
		
		if_cannot move to any vertex (all clear) then_break;
		
		
		push the vertex to tabuList, remove oldest vertex from TabuList in the_case TabuList is full
		
	  }//for iDists loop
	  

	  color_the_vertex
	  
	}//for iMoves loop
  }//for ant loop
  
  
  
  totalConflicts=updateConflictTableW(tempConfigA,conflictsTable);  
  totalConflictsAHEAD=updateConflictTableAHEAD(tempConfigA,conflictsTableAHEAD,AHEAD_DIST);  


  if(!totalConflicts){ 

	run Local Opt 
	save the best
	reduce attemptMaxK
  }//if no conflict
	
	

  if (iCycles-changedCycle == Q_CHANGE_CYCLE){
	increase attempMaxK
  }


  if(iCycles-changedCycle==joltCycles){
	jolt 10 percent randomly 
  }
	



  if(iCycles-changedCycle==breakCycles){
	break;
  }
	
	
}//cycle loop
