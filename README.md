# RQMC-simulation-hub
This repository contains tools and examples for the simulation of stochastic models with Monte Carlo and randomized quasi-Monte Carlo. Most of the provided packages depend on the library [Stochastic Simulation in Java (SSJ)](http://simul.iro.umontreal.ca/ssj/), get it [here](https://github.com/umontreal-simul/ssj), and have been used to simulate benchmark problems in our research papers.

The packages included in this repository are 

* [tauleap](#Tauleap): The simulation of chemical reaction networks with tau-leaping and Array-RQMC (see our [paper](https://arxiv.org/abs/2009.00337)).
* [cde](#CDE): Density estimation for stochastic models via conditioning (see our [paper](https://arxiv.org/abs/1906.04607)). **To be added soon.** 

The repository and this list will be updated continuously to provide the code for other simulation problems that arise from my research.

## Tauleap
The package offers *abstract classes* to easily implement your own chemical reaction network, as well as your own sorting algorithms for Array-RQMC. We also provide instances of each of these abstract classes to demonstrate how they can be used in practice.  Finally, we also include a class containing a main procedure that was used to produce the results from the [reference paper](https://arxiv.org/abs/2009.00337).

### Chemical reaction networks
They are implemented in the abstract class [ChemicalReactionNetwork](https://github.com/FlorianPuchhammer/rqmc-simulation-hub/blob/main/src/main/java/tauleap/ChemicalReactionNetwork.java). This means that all the important procedures for this type of object, like computing the next tau-leap iteration, are already implemented and to create a specific instance of a chemical reaction network, only system specific parameters need to be given. In particular, this concerns the system's propensity functions, a performance measure and similar things, as well as initial values, reaction rates, the stoichiometric matrix, tau, the simulation time, etc.

To see examples for implementations of a [ChemicalReactionNetwork](https://github.com/FlorianPuchhammer/rqmc-simulation-hub/blob/main/src/main/java/tauleap/ChemicalReactionNetwork.java), see the classes [PKA](https://github.com/FlorianPuchhammer/rqmc-simulation-hub/blob/main/src/main/java/tauleap/PKA.java), [ReversibleIsomerization](https://github.com/FlorianPuchhammer/rqmc-simulation-hub/blob/main/src/main/java/tauleap/ReversibleIsomerization.java), [SchloeglSystem1D](https://github.com/FlorianPuchhammer/rqmc-simulation-hub/blob/main/src/main/java/tauleap/SchloeglSystem1D.java), and [SchloeglSystem2D](https://github.com/FlorianPuchhammer/rqmc-simulation-hub/blob/main/src/main/java/tauleap/SchloeglSystem2D.java). For specialized instances, in which the states are transformed into the unit cube, see [PKA01](https://github.com/FlorianPuchhammer/rqmc-simulation-hub/blob/main/src/main/java/tauleap/PKA01.java) and [SchloeglSystem01](https://github.com/FlorianPuchhammer/rqmc-simulation-hub/blob/main/src/main/java/tauleap/SchloeglSystem01.java), for instance. Note that the latter two classes require the presence of certain csv-files to construct the transformation. For more information refer to the comments in the code.

In most of these classes, the observed performance measure is the copy number of one species of the system. More sophisticated performance measures can easily be implemented by adapting the function *getPerformance().*

### Sorting algorithms
Most of the multivariate sorting algorithms are already implemented in the dependency [SSJ](http://simul.iro.umontreal.ca/ssj/). This package contains the specialized sorting method, in which the potentially multi-dimensional state of a Markov chain is mapped  via a *score function* or *importance function* to the real numbers, where sorting can be trivially done by value. This kind of sort is implemented in the abstract class [MultDimToOneDimSort](https://github.com/FlorianPuchhammer/rqmc-simulation-hub/blob/main/src/main/java/tauleap/MultDimToOneDimSort.java). The only abstract method is the score function itself. Examples for this type of sort are given in [PKAOSLASort](https://github.com/FlorianPuchhammer/rqmc-simulation-hub/blob/main/src/main/java/tauleap/PKAOSLASort.java) and in [SchloeglOSLASort](https://github.com/FlorianPuchhammer/rqmc-simulation-hub/blob/main/src/main/java/tauleap/SchloeglOSLASort.java), which gives the so-called *one-step look ahead importance function (OSLAIF)* for the [PKA](https://github.com/FlorianPuchhammer/rqmc-simulation-hub/blob/main/src/main/java/tauleap/PKA.java) and the [SchloeglSystem2D](https://github.com/FlorianPuchhammer/rqmc-simulation-hub/blob/main/src/main/java/tauleap/SchloeglSystem2D.java), respectively. 

### Runnable files
The class [TestChemicalReactionNetwork](https://github.com/FlorianPuchhammer/rqmc-simulation-hub/blob/main/src/main/java/tauleap/TestChemicalReactionNetwork.java) contains the *main()* method that was used for all experiments in the [reference paper](https://arxiv.org/abs/2009.00337). It simulates a chemical reaction network with tau-leaping and Array-RQMC for each specified sort, each specified RQMC-point set, for various point set sizes, and performs an indicated number of replications to estimate several performance metrics. To use this class, carefully follow the comments provided in the code.

## CDE
To be added soon.
