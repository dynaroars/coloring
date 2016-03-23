
/*
  This file contains Parameters, Data structure and Global Variables used in program
*/

/************* PARAMS ***************/
enum {GCP,BCP,MCP,BMCP};int PTYPE=-1;
string PTYPE_name[4]={"GCP","BCP","MCP","BMCP"};
enum {NOCHECK,TOCHECK};
enum {RANDOM_CHOOSE,MAX_DEG,MAX_COLOR};

const int IG_OPT=RANDOM_CHOOSE;


const double nAntsPercent = 1 ;   //change to 1.0 later
const int nCyclesFactor = 6 ; //nCycles = this * nVertices,  
const double nBreakCycleFactor = 1; //change to 1.0 later
const double nMovesPercent=.33;
const double nJoltPercent=.25;
const double rSizeLimitFactor=2.0;
const double scoreKFactor = 1.0;
const double scoreConflictFactor = 10.0;



const double ALPHA = .8 ; 
const double BETA = .2 ; //choosing BETA*IG colors as candidate to select the smallest one to detel



const int nDists=2;

const int MOVE_METHOD = RANDOM_CHOOSE;
const int Q_CHANGE_CYCLE=5;
const int AHEAD_DIST = 2 ;

/************* DATA TYPE ***************/
struct vertex{
  int v_id;
  vector<int>adj;
  vector<int>edgeList;
  vector<double>phColor;
  double phColorSumScore; // stores sum of all values in phColor

};



/************* GLOBAL ***************/
time_t seed_t = 0;
int nVertices=0;
int nEdges=0;
int nAnts=0;
int nCycles=0;
int joltCycles=0;
int breakCycles=0;
int nMoves=0;
int rSizeLimit=0;
int nRLFSetLimit=0;
char* inputFile="";
char* outputFile="SOLUTION.S";


int bestCycle=0;


vertex **pVertices;



bool BB=false ;//debug option


//results stores here

int maxWeight=0,maxK=0;

vector<int>distinctColorsV;


int bestResult=-1;

vector<vector<int> >weightMatrix;
vector<vector<bool> >edgeMatrix;
vector<int>ku;

//store the maxK obtained from intial greedy, this will be the based 
//for scoring, i.e., all later configs with maxK_prime will have score 
//initScore_MaxK / maxK_prime  ....  that is , the lower the maxK_prime
//the higher the score.
int baseScore_MaxK=0; 
vector<vector<int> >allConfs;
