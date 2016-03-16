/************* PARAMS ***************/
/*
  This file contains Parameters, Data structure and Global Variables used in program
*/

/************* PARAMS ***************/
enum class MoveOpts {Random, MaxConflict};

const auto nAntsPercent = .2 ; 
const auto nCyclesFactor = 6 ; //nCycles = this * nVertices, TVN todo: this is too much cycles
const auto nJoltPercent = .5 ; 
const auto nBreakCycleFactor = 2.5;
const auto moveLimitPercent=.25;
const auto rSizeLimitFactor=3;
const auto RLFSetPercent = .7 ; 
const auto XLRF_METHOD = 0 ; //0 = adj density based, vertex w/ largest density chosen for MXRLF

const auto ALPHA = .7 ; 
const auto BETA = .5 ; 
const auto DELTA = .7 ; 

const auto HOW_FAR = 2;

const auto MOVE_METHOD = MoveOpts::Random;
const auto Q_CHANGE_CYCLE=20;

/************* DATA TYPE ***************/
struct vertex{
  int id;
  vector<int> adj;
  vector<int> edgeList;
  int numAnts; 
};

struct ant{
  int id;
  vertex *current;
  vertex *old;
};

/************* GLOBAL ***************/
time_t seed_t = 0;
auto nVertices=0;
auto nEdges=0;
auto nAnts=0;
auto nCycles=0;
auto nJolts=0;
auto breakCycles=0;
auto moveLimit=0;
auto rSizeLimit=0;
auto nRLFSetLimit=0;
auto bestCycle=0;

vertex **pVertices;
vector<ant *> vAnts;

auto BB = false ;//debug option

//results stores here
vector<int> bestColorResult;
auto bestResult = -1;
