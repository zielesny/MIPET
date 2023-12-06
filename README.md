[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/zielesny/MIPET/graphs/commit-activity)
[![GitHub issues](https://img.shields.io/github/issues/zielesny/MIPET.svg)](https://GitHub.com/zielesny/MIPET/issues/)
[![GitHub contributors](https://img.shields.io/github/contributors/zielesny/MIPET.svg)](https://github.com/zielesny/MIPET/graphs/contributors/)
[![GitHub release](https://img.shields.io/github/release/zielesny/MIPET.svg)](https://github.com/zielesny/MIPET/releases/)

# Mesoscopic Interaction Parameter Estimation with Tinker (MIPET)

## Content of this document
* [Description](#Description)
* [Contents of this repository](#Contents-of-this-repository)
  * [Sources](#Sources)
  * [Molecules](#Molecules)
  * [Mathematica notebooks for visualization](#mathematica-notebooks-for-visualization)
  * [Tutorial](#tutorial)
  * [Publication 2024](#publication-2024)
    * [C10E4-water bilayer formation study](#c10e4-water-bilayer-formation-study)
    * [Gaussian job files](#gaussian-job-files)
  * [Legacy](#legacy)
* [Installation](#Installation)
    * [Application](#application)
    * [Tinker - Software Tools for Molecular Design](#tinker---software-tools-for-molecular-design)
* [Dependencies](#dependencies)
* [Acknowledgments](#acknowledgments)

## Description
MIPET (Mesoscopic Interaction Parameter Estimation with Tinker) is an open-source software project that calculates 
non-bonding intermolecular interaction energies and determines coordination numbers. The results of these calculations 
are used by MIPET to generate a "particle set" for the execution of a Dissipative Particle Dynamics (DPD) simulation.
<a href="https://dasher.wustl.edu/tinker/">Tinker</a> tools are utilized for computing both interaction energies and 
coordination numbers.

## Content of this repository 
### Sources 
The <a href="https://github.com/zielsny/MIPET/tree/main/resources/de/whs/ibci/mipet">"src"</a> subfolder contains all
source code packages. The tests can be found in the 
<a href="https://github.com/zielsny/MIPET/tree/main/test/de/whs/ibci">"test"</a> subfolder.

### Molecules
The <a href="https://github.com/zielsny/MIPET/tree/main/Molecules"> Molecules</a> subfolder contains Tinker .xyz files 
for different molecules for the force fields AMOEBA09, MM3, MMFF, OPLSAA and SMOOTHAA as well as their SMILES codes in 
the <a href="https://github.com/zielsny/MIPET/tree/main/Molecules/SMILES">SMILES</a> subfolder. The Tinker .xyz file is 
the basic Tinker file type for Cartesian coordinates. After a tile line it contains one line for each atom in the structure. 
Each atom respectively its line contains the sequential number within the structure, an atomic symbol or 
name, X-, Y-, and Z-coordinates, the force field atom type number of the atom, and a list of connected atoms 
[<a href="https://dasher.wustl.edu/tinker/downloads/tinker-guide.pdf">Tinker guide</a>]. To create a new Tinker .xyz file 
a mol file of the molecule can be transformed into a Tinker .xyz file using the Tinker tool *mol2xyz*. The corresponding 
force field atom type numbers of the atoms must be added manually. A more detailed instruction can be found
in the <a href="https://github.com/zielesny/MIPET/tree/main/Tutorial">Tutorial</a>.

### Mathematica notebooks for visualization
The subfolder <a href="https://github.com/zielesny/MIPET/tree/main/Mathematica%20notebooks%20for%20visualization">Mathematica notebooks for visualization</a> includes
multiple Mathematica notebooks designed to aid in evaluating and visualizing the results of the pipeline. Each notebook
includes a short introductory tutorial on its usage. It is important to note that a licensed version of 
Mathematica 13 (or higher) is required to use these notebooks.

### Tutorial
The <a href="https://github.com/zielesny/MIPET/tree/main/Tutorial">"Tutorial" folder</a> contains a detailed tutorial 
on how to install and use MIPET, together with tutorial job files.

### Publication 2024
The subfolder <a href="https://github.com/zielesny/MIPET/tree/main/Publication%202024">Publication 2024<a>  contains supplementary and supporting materials for our forthcoming article on the calculation pipeline, set to be published in 2024.

#### C10E4-water bilayer formation study
The subfolder <a href="https://github.com/zielesny/MIPET/tree/main/Publication%202024/C10E4-water%20bilayer%20formation%20study">C10E4-water bilayer formation study</a>
contains a bunch of job files for DPD simulations performed with MFsim. Detailed information can be found in our 
scientific articles, that will be published 2024.

#### Gaussian job files
The subfolder <a hre="https://github.com/zielesny/MIPET/tree/main/Publication%202024/Gaussian%20job%20files"> Gaussian job files</a>
contains a bunch of job files for DFT calculations that were carried out as part of the development of MIPET using Gaussian.
Detailed information can be found in our scientific articles, that will be published 2024.

### Legacy
The <a href="https://github.com/zielesny/MIPET/tree/main/Legacy"> Legacy </a> subfolder contains legacy source code and 
files which were developed as part of the development of MIPET.

## Installation
### Application
Pre-compiled and executable MIPET distributions can be found attached to the 
<a href="https://github.com/zielesny/MIPET/releases">marked releases</a>.

**Windows:**  
A convenient Windows OS installer executable for MIPET is available 
(click <a href="https://github.com/zielesny/MIPET/releases/download/MIPET/MIPET_win_x86_64.exe">here</a> 
to automatically download the installer .exe of the latest version). Download the installer .exe, start, and follow the
instructions to install MIPET. Note that the installation includes a full Java Runtime Environment (JRE). A more 
detailed instruction can be found in the <a href="https://github.com/zielesny/MIPET/tree/main/Tutorial">Tutorial</a>.

**Linux:**  
Every release has a zip file attached which contains the MIPET Java archive with all dependencies (except Tinker) and a
full Java runtime environment 
(click <a href="https://github.com/zielesny/MIPET/releases/download/MIPET/MIPET_linux_x86_64.tar.gz">here</a> 
to automatically download the zip file of the latest version). A more detailed instruction can be found
in the <a href="https://github.com/zielesny/MIPET/tree/main/Tutorial">Tutorial</a>.

### Tinker - Software Tools for Molecular Design
In order to run MIPET, the Tinker molecular modeling tool package has to be added. For this, the Tinker executables and 
the files for the atomic parameters of version 8.10.2 must be downloaded 
(click <a href="https://dasher.wustl.edu/tinker/downloads/">here</a> to go to the Tinker website). Choose the correct version for
your OS. The atomic parameters can be found in the "Tinker Package Distribution" archive. 
Inside is a directory named "params" in which all atomic parameter sets can be found.
Please note the current <a href="https://dasher.wustl.edu/tinker/downloads/license.pdf">licence</a>. 
The executables and the parameter files have to be placed in a subfolder named "tinker", which you have to create in the
root directory of MIPET.\
A more detailed instruction can be found in the <a href="https://github.com/zielesny/MIPET/tree/main/Tutorial">Tutorial</a>.


Example directory structure:  

    ├── MIPET                   # MIPET root directory  
    │   ├── ...              
    │   ├── tinker              # Tinker directory
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
