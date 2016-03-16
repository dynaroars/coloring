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
