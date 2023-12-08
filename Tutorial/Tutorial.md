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
MIPET (Mesoscopic Interaction Parameter Estimation with Tinker) is an automated pipeline that comprehensively calculates
intermolecular interaction energy and coordination numbers using molecular force-fields, all through the <a href="https://dasher.wustl.edu/tinker/">Tinker</a>
molecular modeling package. A scientific article about MIPET is forthcoming. To download the software, browse its source
code, and report issues, access the <a href="https://github.com/zielesny/MIPET">MIPET GitHub repository</a>.

## Installation and start-up
Pre-built MIPET distributions that are ready to be executed can be accessed through the 
<a href="https://github.com/zielesny/MIPET/releases">marked releases</a>. Please note that MIPET only supports x64 architecture.\
This tutorial configures a short "proof-of-work" calculation on the first run after a fresh installation.  The calculation
typically lasts around 10 to 15 minutes. The tutorial files are located in the Tutorial folder of MIPET and should be 
duplicated to the relevant locations to replicate the tutorial calculation.

### Windows
A convenient dedicated Windows OS installer executable for MIPET is readily available for
<a href="https://github.com/zielesny/MIPET/releases/download/MIPET/MIPET_win_x86_64.exe"> download </a>, 
which includes a full Java runtime environment. Download and run the installation program and follow the installation 
prompts. Once MIPET is installed successfully, proceed with the installation of [Tinker](#tinker) by following the 
provided instructions. To execute the pipeline, it is imperative to configure a scratch directory (use "/" instead of "\" as delimiter) within the MIPET.properties file, which can be located in the MIPET root directory. A directory on a RAM
disk would be ideal. If the default install directory is utilized, modifying the MIPET.properties file necessitates 
admin privileges. As a workaround, copy the file to a writable directory, personalize it, and copy it back to its 
original location. Subsequently, double-click on the "MIPET.bat" batch file.

If this installation or the execution of the batch files does not function properly, refer to the Linux instructions 
provided below. They can also be applied as an alternative method on Windows.

### Linux
Every version includes a zip file that comprises the MIPET Java archive containing all dependencies except Tinker, 
along with a comprehensive Java runtime environment. By clicking 
<a href="https://github.com/zielesny/MIPET/releases/download/MIPET/MIPET_linux_x86_64.tar.gz">here</a>, 
you can automatically download the zip file of the latest version. Follow the instructions for installing [Tinker](#tinker)
after extracting the zip file.  Before running the pipeline, configure a scratch directory (use "/" instead of "\" as delimiter) in the "MIPET.properties" file, located in the root directory of MIPET. A directory on a RAM disk would be 
ideal. Then, execute the JAR from the command line.

### Tinker
Please note the current <a href="https://dasher.wustl.edu/tinker/downloads/license.pdf">licence</a> of Tinker.
To use MIPET, add the Tinker molecular modeling tool package. For this task, you must download Tinker executables 
version 8.10.2 and the atomic parameter files. Visit the Tinker downloads website by clicking <a href="https://dasher.wustl.edu/tinker/downloads/">here</a>
to access the downloads. Make sure to choose the appropriate version for your operating system:

  * <a href="https://dasher.wustl.edu/tinker/downloads/bin-windows-8.10.2.zip"> Windows </a>
  * <a href="https://dasher.wustl.edu/tinker/downloads/bin-linux-8.10.2.tar.gz"> Linux</a>.

The atomic parameter sets are located in the "Tinker Package Distribution" archive:

  * <a href="https://dasher.wustl.edu/tinker/downloads/tinker-8.10.2.zip"> Windows </a>
  * <a href="https://dasher.wustl.edu/tinker/downloads/tinker-8.10.2.tar.gz"> Linux</a>.

Inside the "params" directory, all atomic parameter sets can be found. To set up MIPET, create a subfolder named 
"tinker" (if it does not already exists) in the root directory and place the executables and parameter files inside. 
Example directory structure:

    ├── MIPET                   # MIPET root directory  
    │   ├── ...              
    │   ├── tinker              # Tinker directory
    │   │  ├── params           # Directory to put parameter files  
    │   │  ├── alchemy.exe      # Tinker executable  
    │   │  ├── analyze.exe      # Tinker executable  
    │   │  └── ...              # more Tinker executables  
    │   └── ...                 

### Results
The results of a MIPET run are available in the directory specified in the "MIPET.properties" file under "MIPET.Directory.result". In addition, the subdirectory "IE" contains the particle sets created for Dissipative Particle Dynamics (DPD) simulations with <a href="https://github.com/zielesny/MFsim">MFsim</a> as well as the results of the intermolecular energy calculation.  
The numerical results of intermolecular energy calculation are in the .dat files located in the subfolders of the dimers within the "IE" folder. The "CN" folder has the coordination number determination results divided into subfolders named after the dimers. A .dat file with the numerical results is also available here.
Mathematica notebooks for visual inspection of results are available in the <a href="https://github.com/zielesny/MIPET/tree/main/Mathematica%20notebooks%20for%20visualization">Mathematica notebooks for visualization</a> subfolder. Each notebook includes a brief description of its functionality and instructions for usage. A license for Mathematica 13 or newer is mandatory for accessing the notebooks.

## Change calculation properties, molecules and force fields
The "MIPET.properties" file, located in the "resources" directory, enables users to modify calculation properties. For
instance, by adjusting the "MIPETSphereNodeNumber" setting, the number of evenly spaced points on the spheres involved 
in the intermolecular energy calculation can be customized. Each property in this file is briefly commented on to describe
its function. To alter the utilized molecules, take a look at the "Particles" section of the "MIPET.job" file. List the 
molecules to be used in this section using their abbreviations. To exclude a molecule from the pipeline, add a hashtag 
before it. Ready-to-use molecules are available in the "Molecules" folder, separated by different force fields. To create
and set up new molecules in the calculation pipeline, refer to the guidelines in the 
[Create and configure molecules](#create-and-configure-molecules) section.\
Various force fields can be utilized to compute intermolecular energy and determine coordination numbers. Precede the 
chosen force field for calculating intermolecular energy with an asterisk and the force field for calculating the 
coordination number with a dollar sign.

### Configure water model
To set up the water model for the OPLSAA and SMOOTHAA force fields, duplicate the chosen water model file inside the 
force field directory, and rename the new copy to "H2O.xyz". The TIP5P water model is automatically selected by default
after a clean installation, and if no changes were made.

## Create and configure molecules
To create a Tinker xyz file for a molecule to be used in the calculation pipeline, a pdb file with at least 2D coordinates
of the desired molecule should be available. Utilizing the Tinker tool "pdbxyz" will transform the pdb file to a Tinker 
xyz file. It is also necessary to have a key file that refers to the force field's desired parameter file such as OPLSAA:
    
    PARAMETERS "<Path to Tinker directory>\Tinker\params\oplsaa.prm"  

The generated Tinker xyz file for methanol using this method appears as follows:

    6    MeOH
    1  C      0.680000   -0.025000    0.000000     0     2     3     4     5
    2  H      1.047000    1.011000   -0.000000     0     1
    3  H      1.009000   -0.547000   -0.904000     0     1
    4  H      1.009000   -0.547000    0.904000     0     1
    5  O     -0.746000    0.124000    0.000000     0     1     6
    6  H     -1.172000   -0.756000    0.000000     0     5

Replace the first row after the coordinates with the corresponding numbers for the atom types from the force field:

    6    MeOH
    1  C      0.680000   -0.025000    0.000000     99     2     3     4     5
    2  H      1.047000    1.011000   -0.000000     98     1
    3  H      1.009000   -0.547000   -0.904000     98     1
    4  H      1.009000   -0.547000    0.904000     98     1
    5  O     -0.746000    0.124000    0.000000     96     1     6
    6  H     -1.172000   -0.756000    0.000000     97     5

This file should be copied to the force field directory within the "Molecules" directory.
