# Tutorial

Last update: 20th November 2023\
Reference version of MIPET 1.0.0.0\

## Contents
* [Introduction](#introduction)
* [Installation and start-up](#installation-and-start-up)
  * [Windows](#windows)
  * [Linux](#linux)
  * [Tinker](#tinker)
* [Change calculation properties, molecules and force fields](#change-calculation-properties-molecules-and-force-fields) 
  * [Configure water model](#configure-water-model)
* [Create and configure molecules](#create-and-configure-molecules)

## Introduction
MIPET (Mesoscopic Interaction Parameter Estimation with Tinker) is an automated pipeline for comprehensive calculation
of intermolecular interaction energies and coordination numbers based on molecular force-fields using the 
<a href="https://dasher.wustl.edu/tinker/">Tinker</a> molecular modeling package.
A scientific article describing MIPET can be found <a href="TODO: add link to paper"> here </a>. 
The MIPET GitHub repository were you can download the software, browse the source code, and report issues can be found 
<a href="https://github.com/zielesny/MIPET">here</a>.

## Installation and start-up
Pre-compiled and executable MIPET distributions can be found attached to the 
<a href="https://github.com/zielesny/MIPET/releases">marked releases</a>. Please note that MIPET only supports x64 
architecture in general.\
With this tutorial and on the first run after a fresh installation a short "proof-of-work" calculation is configured. 
It should take about 10 to 15 min.
The tutorial files can be found in the <a href="https://github.com/zielesny/MIPET/tree/main/Tutorial">Tutorial</a> folder 
of MIPET and must be copied to the appropriate 
places to repeat the tutorial calculation. 

### Windows
A convenient dedicated Windows OS installer executable for MIPET is available
<a href="https://github.com/zielesny/MIPET/releases/download/MIPET/MIPET_win_x86_64.exe"> here </a>, 
it comprises a full Java runtime environment. Download the installer executable, start, and follow the instructions to
install the calculation pipeline. Follow the instructions for the [Tinker](#tinker) installation. To run the pipeline, 
a scratch directory must be configured in the MIPET.properties file, which can be found in the "resources" directory. 
Then simply double-click on the batch file "MIPET.bat".

Should this installation or the execution of the batch files not work for you, try the guidelines for Linux described 
below. As an alternative way, they should also work on Windows.

### Linux
Every release has a zip file attached which contains the MIPET Java archive with all dependencies (except Tinker) and a
full Java runtime environment 
(click <a href="https://github.com/zielesny/MIPET/releases/download/MIPET/MIPET_linux_x86_64.tar.gz">here</a> to 
automatically download the zip file of the latest version). 
After extracting the zip file, follow the instructions for the [Tinker](#tinker) installation. To run the pipeline, 
a scratch directory must be configured in the "MIPET.properties" file, which can be found in the "resources" directory. 
Then execute the JAR from the command line.

### Tinker
In order to run MIPET, the Tinker molecular modeling tool package has to be added. For this, the Tinker executables and 
the files for the atomic parameters of version 8.10.2 have to be downloaded 
(click <a href="https://dasher.wustl.edu/tinker/">here</a> to go to the Tinker website). Choose the correct version for 
your OS. Please note the current <a href="https://dasher.wustl.edu/tinker/downloads/license.pdf">licence</a>. 
The executables and the parameter files have to be placed in a subfolder named "tinker", which you have to create in the
root directory of MIPET.
Exampled directory structure:

    ├── MIPET                   # MIPET root directory  
    │   ├── ...              
    │   ├── tinker              # Frequently asked questions  
    │   │  ├── params           # Directory to put parameter files  
    │   │  ├── alchemy.exe      # Tinker executable  
    │   │  ├── analyze.exe      # Tinker executable  
    │   │  └── ...              # more Tinker executables  
    │   └── ...                 

## Change calculation properties, molecules and force fields
The "MIPET.properties" file,  which can be found in the "resources" directory, allows to change properties for a 
calculation. For example, to adjust the number of evenly spaced points on the spheres for the calculation of the 
intermolecular energy, the "MIPETSphereNodeNumber" setting can be changed. Each property in this file has a short comment
describing what it does. To change the used molecules the file "MIPET.job" holds a section named "Particles". The molecules 
to use must be listed in this section by using their abbreviations. To not consider a molecule during the pipeline a
preceded hashtag is needed. Configured and ready-to-use molecules can be found in the "Molecules" folder, seperated by 
different force fields. To create and configure new molecules for the calculation pipeline follow the instruction of the
section [Create and configure molecules](#create-and-configure-molecules).\
Different force fields can be used to calculate the intermolecular energy and determine the coordination number. The 
force field to be used for calculating the intermolecular energy must be preceded by an asterisk and the force field 
to be used for calculating the coordination number must be preceded by a dollar sign.

### Configure water model
To configure the water model for the force fields OPLSAA and SMOOTHAA just copy and paste the file of the chosen water 
model inside the force field directory and rename the new file to "H2O.xyz". By default, after a fresh installation, and 
if nothing was changed, the TIP5P water model is chosen.

## Create and configure molecules

