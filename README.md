# AntColor: An Ant-based Algorithm for Graph Coloring problems

**AntColor** implements a very efficient, heuristic ant-based algorithm for the (classical) Graph Coloring problem. **AntColor** also supports several popular generalizations, namely the Bandwidth Coloring, Multi Coloring, and Bandwidth Multi Coloring problems.

Code is written in C/C++ and released under the MIT license.

## Classic Graph Coloring

```
$ cd Classic/src
$ g++ classic.cc -std=c++11  -o classic
$ ./classic ../examples/DSJC125.1.col.b 1  #optional -D DB : turn on debugging and output sol to sol.txt
nthreads: 1 colors: 5 vertices: 125 edges: 736 bestCycle: 160 seed: 1
```

In the above example, `../examples/DSJC125.1.col.b` is the input graph in DIMACS *binary* format, and the [optional] integer `1` is the seed. The result `colors: 5` indicates this graph can be colored using 5 colors.


## MISCS
You might also be interested in

* **Benchmark graphs**: [https://github.com/dynaroars/npbench/](https://github.com/dynaroars/npbench/)
* **[Converter](https://github.com/dynaroars/npbench/tree/master/instances/converter)**: convert from DIMACS binary to ASCII and vice versa

## Publications
1. Bui, T., T. Nguyen, C. Patel, and K. Phan, "An Ant-Based Algorithm for Coloring Graphs," Journal of Discrete Applied Mathematics, Vol. 156(2), 2008, pp. 190 --- 200.
1. T. Nguyen., M.S. Thesis: "On the Graph Coloring Problem and Its Generalizations," Penn State University, 2006
1. Bui, T and T. Nguyen, "An Agent-Based Algorithm for Generalized Graph Colorings," Proc. 8th Annual Conference on Genetic and Evolutionary Computation Conference (GECCO), 2006, pp. 19 --- 26

