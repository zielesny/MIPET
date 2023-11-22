# Tutorial

Last update: 20th November 2023\
Reference version of MIPET 1.0.0.0\

## Contents
* [Introduction](#introduction)
* [Installation and start-up](#installation-and-start-up)
  * [Windows](#windows)
  * [Linux](#linux)
  * [Tinker](#tinker)

## Introduction
MIPET (Mesoscopic Interaction Parameter Estimation with Tinker) is an automated pipeline for comprehensive calculation of intermolecular interaction energies and coordination numbers based on molecular force-fields using the <a href="https://dasher.wustl.edu/tinker/">Tinker</a> molecular modeling package.
A scientific article describing MIPET can be found <a href="TODO: add link to paper"> here </a>. The MIPET GitHub repository were you can download the software, browse the source code, and report issues can be found <a href="https://github.com/zielesny/MIPET">here</a>.

## Installation and start-up
Pre-compiled and executable MIPET distributions can be found attached to the <a href="TODO: add link to releases">marked releases</a>. Please note that MIPET only supports x64 architecture in general.

### Windows
A convenient dedicated Windows OS installer executable for MIPET is available <a href="TODO: add link to release"> here </a>, it comprises a full Java runtime environment. Download the installer executable, start, and follow the instructions to install the calculation pipeline. Follow the instructions for the [Tinker](#tinker) installation. To run the pipeline, double-click the batch file "MIPET.bat". 

Should this installation or the execution of the batch files not work for you, try the guidelines for Linux described below. As an alternative way, they should also work on Windows.

### Linux
Every release has a zip file attached which contains the MIPET Java archive with all dependencies (except Tinker) and a full Java runtime environment (click <a href="TODO add link to zip file for linux">here</a> to automatically download the zip file of the latest version). After extracting the zip file, follow the instructions for the [Tinker](#tinker) installation. To run MIPET, execute the JAR from the command line. 

### Tinker
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
