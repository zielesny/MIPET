# Mesoscopic Interaction Parameter Estimation with Tinker (MIPET)

## Content of this document
* [Description](#Description)
* [Contents of this repository](#Contents-of-this-repository)
  * [Sources](#Sources)
  * [Molecules](#Molecules)
  * [Visualization Tools](#Visualization-Tools)
  * [Tutorial/Usage](#tutorialusage)
  * [Installation](#Installation)
    * [Application](#application)
    * [Source Code](#source-code)
    * [Tinker - Software Tools for Molecular Design](#tinker---software-tools-for-molecular-design)

## Description
MIPET ('Mesoscopic Interaction Parameter Estimation with Tinker') is an open software project for calculating non-bonding intermolecular interaction energies and determining coordination numbers. With the results, MIPET can generate a "particle set" for Dissipative Particle Dynamics (DPD) simulation. Tools from the <a href="https://dasher.wustl.edu/tinker/">Tinker</a> are used for both the calculation of interaction energies and coordination numbers.

## Content of this repository 
### Sources 
The <a href="https://github.com/zielsny/MIPET/tree/main/resources/de/whs/ibci/mipet">"src"</a> subfolder contains all source code packages. The tests can be found in the <a href="https://github.com/zielsny/MIPET/tree/main/test/de/whs/ibci">"test"</a> subfolder.

### Molecules
The <a href="https://github.com/zielsny/MIPET/tree/main/Molecules"> Molecules</a> subfolder contains Tinker .xyz files for different molecules for the force fields AMOEBA09, MM3, MMFF, OPLSAA and SMOOTHAA as well as their SMILES codes in the <a href="https://github.com/zielsny/MIPET/tree/main/Molecules/SMILES">SMILES</a> subfolder.
The Tinker .xyz file is the basic Tinker file type for Cartesian coordinates. After a tile line it contains one line for each atom in the structure. Each atom respectively its line contains the sequential number within the structure, an atomic symbol or name, X-, Y-, and Z-coordinates, the force field atom type number of the atom, and a list of connected atoms [<a href="https://dasher.wustl.edu/tinker/downloads/tinker-guide.pdf">Tinker guide</a>]. 
To create a new Tinker .xyz file a mol file of the molecule can be transformed into a Tinker .xyz file using the Tinker tool *mol2xyz*. The corresponding force field atom type numbers of the atoms must be added manually.

### Visualization Tools
The subfolder <a href="https://github.com/zielsny/MIPET/tree/main/Visualization_Tools">Visualization Tools</a> includes multiple Mathematica notebooks designed to aid in evaluating and visualizing the results of the pipeline. Each notebook includes a short introductory tutorial on its usage. It is important to note that a licensed version of Mathematica 13 (or higher) is required to use these notebooks.

### Tutorial/Usage
Not yet implemented

## Installation
### Application
Pre-compiled and executable MIPET distributions can be found attached to the <a href="">marked releases</a>.

**Windows:**  
A convenient Windows OS installer executable for MIPET is available (click <a href="">here</a> to automatically download the installer .exe of the latest version). Download the installer .exe, start, and follow the instructions to install MIPET. Note that the installation includes a full Java Runtime Environment (JRE).  

**Linux:**  
...

### Source Code
Not yet implemented

### Tinker - Software Tools for Molecular Design
In order to run MIPET, the Tinker molecular modeling tool package has to be added. For this, the Tinker executables and the files for the atomic parameters of version 8.10.2 have to be downloaded (click <a href="https://dasher.wustl.edu/tinker/">here</a> to go to the Tinker website). Choose the correct version for your OS. Please note the current <a href="https://dasher.wustl.edu/tinker/downloads/license.pdf">licence</a>. The executables and the parameter files have to be placed in a subfolder named "tinker", which you have to create in the root directory of MIPET.
Exampled directory structure:  

    ├── MIPET                   # MIPET root directory  
    │   ├── ...              
    │   ├── tinker              # Frequently asked questions  
    │   │  ├── params           # Directory to put parameter files  
    │   │  ├── alchemy.exe      # Tinker executable  
    │   │  ├── analyze.exe      # Tinker executable  
    │   │  └── ...              # more Tinker executables  
    │   └── ...                 

## Dependencies
* Java Development Kit (JDK) version 21
    * [AdoptOpenJDK](https://adoptopenjdk.net) (as one possible source of the JDK)
* Chemistry Development Kit (CDK) version 2.9
    * [Chemistry Development Kit on GitHub](https://cdk.github.io/)
* Tinker - Software Tools for Molecular Design version 8.10.2
    * [Tinker](https://dasher.wustl.edu/tinker/)
* JUnit version 4.13.2
    * [JUnit 4](https://junit.org/junit4/)
* Apache commons-lang3 version 3.12.0
    * [commons-lang3](https://commons.apache.org/proper/commons-lang/)

## Acknowledgments 
### Developers
* Felix Baensch
* Mirco Daniel
* Veit Hucke

### Initialization, conceptualization, and supervision:
* Achim Zielesny
* Christoph Steinbeck

### MIPET was developed at:
Zielesny Research Group  
Westphalian University of Applied Sciences  
August-Schmidt-Ring 10  
D-45665 Recklinghausen Germany
