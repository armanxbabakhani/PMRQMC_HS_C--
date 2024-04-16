# Permutation Matrix Representation Quantum Monte Carlo for High Spins

## Table of Contents

- [Introduction](#introduction)
- [The Code](#thecode)
- [Contributions](#contributions)
- [License](#license)

## Introduction

This program PMRQMC.cpp simulates the equilibrium quantities for a specified Hamiltonian (using an input file). The only simulation parameter is the inverse temperature of the system that can be altered directly by editing PMRQMC.cpp file (it is set to $\beta = 1$ by default).
This code takes in an input file specifying the Hamiltonian in sum of pauli strings.
The code consists of two parts, descibed as follows.

## The Code

### A. Data to PMR (DatatoPMR.h)

This header file provides all the necessary functions to convert the pauli string input data provided by the user into PMR representation. The data structure describing the PMR formalism is the following 

```
typedef vector<vector<complex<double>>> Coeffs;
typedef vector<vector<int>> ParticlePerm;
typedef vector<vector<pair<int,int>>> TotalPerms;
typedef pair<int,int> ParticleDiag;

struct TotalDiag {
    int ztype , k , particle;
};

typedef vector<vector<vector<ParticleDiag>>> ParticleDVecs;
typedef vector<vector<vector<TotalDiag>>> TotalDVecs;
typedef pair<vector<int> , pair<vector<ParticleDiag> , vector<complex<double>>>> PauliCDPs;

struct PDdata {
    int twoSplusOne , NumOfParticles;
    TotalPerms Permutations;
    TotalDVecs Diagonals;
    Coeffs Coefficients;
    vector<vector<TotalDiag>> D0;
    vector<complex<double>> D0Coeffs;
};
```
In the DatatoPMR.hpp ParticlePerm data type is used to extract the permutation for a particular particle. For example, $X_1 ^2$ would produce a $P^0$, $P^1$, and a $P^2$ term. So, the powers of the permutation would be stored in the form of a ParticlePerm data structure as [0,1,2] .

## Contributions

The main contributors to the code are Lev Barash and Arman Babakhani. The main contributions to designing the algorithm has been done by Prof. Itay Hen.

## License

Specify the project's license.
