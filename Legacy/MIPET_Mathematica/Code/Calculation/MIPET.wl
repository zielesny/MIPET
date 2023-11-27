(* ::Package:: *)

(*
-----------------------------------------------------------------------
Mesoscopic Interaction Parameter Estimation with Tinker (MIPET)
Version 0.9.6 for Mathematica 12.0
-----------------------------------------------------------------------

Authors: Mirco Daniel, Felix Baensch
Institute for Bioinformatics and Chemoinformatics
Westphalian University of Applied Sciences, Recklinghausen, Germany

Copyright 2023 Daniel, Baensch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL) as 
published by the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
Lesser General Public License (LGPL) for more details.

You should have received a copy of the GNU Lesser General Public 
License along with this program. If not, see 
<http://www.gnu.org/licenses/>. 
-----------------------------------------------------------------------
*)
(* 
Summary:
This package calculates the intermolecular (non-bonded) energy of two molecules 
and determines the coordination number of a molecule A surrounded by solvent 
molecules B.

For the calculation of intermolecular energy of two same or of two different 
molecules some steps are required. First of all, tinker xyz-files for both 
molecules with correct atom types are necessary. To do this the 6th column of
tinker xyz-file has to fill up with correct atom types. You can find them in 
the .prm file in tinker/params directory.

Next step is the optimization of the molecule structure. For this MIPET calls 
tinker's tool 'optimize' and 'scan'. The tool 'scan' does the conformational 
search and the optimized structure are saved in code/calculation/optXYZ directory.
MIPET uses the structure with lowest energy.

To find out the intermolecular energy of two molecules all kind of configurations 
have to be anlysed. MIPET tries select finite configurations of the two molecules. 
The higher the better the estimation of intermolecular energy but there is more 
computing time required. Firstly, the unweighted center of the two molecules is 
calculated. Then fixed number of rotamers by rotating around the center are calculated. 
For that rotamers equaly distributed sphere nodes on the sphere are used. Two of 
the sphere nodes and both centers of two molecules are broght to same line and the 
intermolecular (non-bonded) energy with tinker tool 'analyze' is calculated. To cover 
all configurations the second molecule is roted around the center to center axis.
For example the analyse of 2*81 nodes and 60 degree there have to be 81*81*6 
configurations calculated. These intermolecular energies are weighted with boltzmann 
distribution. After that the center to center distance is changed.

To find the minimum of the intermolecular energy the center to center distance is roughly
and then more finely changed. Lastly, the structure with the lowest intermolecular energy
is optimized with tinker's 'optimize' tool.    

Coordination number determination uses the Tinker tool 'dynamic'.
The coordination number in a dynamic simulation is the number of neighbor molecules
around of the solved molecule. The distance must be smaller than the sum of the vdW-radii
of adjacent atoms from solved molecule and its neighbor molecule as well as capture radius.
Because of the dynamic nature of the simulation the coordination number is fluctuating 
over the time. MIPET observes the coordination numbers from simulation and take the 
mean value.    
To determine the box length of the simulation box, the density of solvent molecule 
has to be known. Many interested molecules are gas at room temperature. The coordination 
number of gas molecules would be zero. So, it is not helpful to take the density of gas.
Moreover the observed fragments are part of a big molecule and its density is in most cases 
unknown. MIPET assumes that all fragments have the same ratio of molecule volume (mole volume
of one molecule) and vdW volume (calculated by VabcVolume methode in CDK package) from water. 
With this assumption the density of simulation box is calculated. *)
(*
Dependencies:
- JLink
- Chemistry Development Kit (CDK) Version 2.8
- Tinker Verison 8.10.3
*)



(* ::Section:: *)
(*Package*)


BeginPackage["MIPET`"]
<< JLink`;
ParallelNeeds["JLink`"];


(* ::Section:: *)
(*Default Options*)


Options[MIPETOptions] = 
	{
		(* Temperatures *)
		MIPETOptionTemperature -> 300,
    
		(* CPU core number *)
		MIPETOptionCPUCoreNumber -> 4,

		(* Boolean value for whether to use a Fibonacci sphere algorithm to determine the *)
		(* sphere node coordinates. If false the coordinates of Fliege et.al. will be used *)
		MIPETOptionsUseFibonacciSphereAlgorithm -> True,
  	
		(* Configuration number of one particle *)
		MIPETOptionSphereNodeNumber -> 4,
    
		(* Rotation configuration number *)
		MIPETOptionRotationNumber -> 5,

		(* Minimum atom distance to prevent tinker error output *)
		MIPETOptionMinAtomDistance -> 1.0,

		(* Force field *)
		MIPETOptionForceFieldName -> "MMFF",
    
		(* Jar Filename *)
		MIPETOptionJarFileName -> "cdk-2.8.jar",
  	
		(* Tinker's scan program *)
		(* 0: automatic selection *)
		(* 1: manual selection of angles to rotate *)
		(* 2: manual selection of angles to freeze *)
		MIPETOptionScanProgram -> 0,
  	
		(* Number of search directions *)
		MIPETOptionNumberSearchDirection -> 5,
  	
		(* Energy threshold for local minima *)
		MIPETOptionEnergyThreshold -> 100,
  	
		(* RMS Gradient for atom criterion *)
		MIPETOptionRmsGradient -> 0.0001,
    
		(* RMS Gradient for BFGS minimization of solvent box *)
		MIPETOptionBfgsRmsGradient -> 0.01,
    
		(* Maximum iteration steps for BFGS minimization of solvent box *)
		MIPETOptionsBfgsMaxIteration -> 100,
    
		(*RMS Gradient for Tinker Optimize - to optimize particles xyz file *)
		MIPETOptionOptimizeRmsGradient -> 0.01,
    
		(* Dielectric constant *)
		MIPETOptionDielectricConstant -> 1,
	  	
		(* Solvent molecule number *)
		MIPETOptionSolventMoleculeNumber -> 400,
    
		(* Step number *)
		MIPETOptionStepNumber -> 1000,
		
		(* Steps per Tinker dynamic round for Cij calculation  *)
		MIPETOptionsStepsPerRound -> 2000,
    
		(* Time step length in fs. If remains unchanged, number of steps equals the simulation length in femtoseconds *)
		MIPETOptionTimeStep -> 1.0,
    
		(* Printinterval in ps (time between dumps) *)
		MIPETOptionPrintinterval -> 0.25,
    
		(* WarmUp step number *)
		MIPETOptionWarmUpStepNumber -> 1000,
    
		(* Tinker Dynamic Simulation Type *)
		(* Available Statistical Mechanical Ensembles: *)
		(* (1) Microcanonical (NVE) *)
		(* (2) Canonical (NVT) *)
		(* (3) Isoenthalpic-Isobaric (NPH) *)
		(* (4) Isothermal-Isobaric (NPT) *)
		MIPETOptionSimulationType -> 2,
    
		(* Catch radius in Angstrom *)
		MIPETOptionCatchRadius -> {.5},
    
		(* Lower boundary for search for intermolecular energy minimum *) 
		MIPETOptionLowerBoundary -> 3,
    
		(* Upper boundary for search for intermolecular energy minimum *) 
		MIPETOptionUpperBoundary -> 10,
    
		(* Prescan step size for search for intermolecular energy minimum *)
		MIPETOptionPrescanStepSize -> 0.5,
    
		(* Step size for search for intermolecular energy minimum *)
		MIPETOptionStepSize -> 0.1,
    
		(* Step size for a more precise search for intermolecular energy minimum *)
		MIPETOptionPreciseStepSize -> 0.01,
    
		(* Boolean value for whether or not to create an overview notebook. *)
		MIPETOptionCreateOverview -> True,
		
		(* Boolean value for whether or not to export charts *)
		MIPETOptionExportCharts -> True,			
		
		(* Title for parameter set *)
		MIPETOptionParameterSetTitle -> "Parameter",
		
		(* Title abbreviation for parameter set *)
		MIPETOptionParameterSetTitleAbr -> "Prm",
		
		(* Fraction of energy values used for the Boltzmann distribution *)
		(* - fraction = 0.0: Only E(nonbonded, min) is used (no "Boltzmann average" is calculated) *)
  		(* - 0.0 < fraction < 1.0: All configurational E(nonbonded) values are sorted ascending and the lower "numberOfValues*fractionForAverage" E(nonbonded) values are used for "Boltzmann average" calculation *)
    	(* 	- Example: For 144x144x16 = 331776 E(nonbonded) values for a specific molecule distance r and a fractionForAverage of 0.25 the lowest Round(331776x0.25) = 82944 E(nonbonded) values are used for "Boltzmann average" calculation only *)
  		(* - fraction = 1.0: all configurational E(nonbonded) values are used for "Boltzmann average" calculation *)
		MIPETOptionsBoltzmannFraction -> 0.0,
		
		(* Boolean value for whether notebooks should be created that show steps of the Zij determination or not - increases computing time *)
		MIPETOptionsExportZijShow -> False,
		
		(* Booelan value for whether to use cell list mehtod by Tildesley to determine Zij or use brute force methoed *)
		MIPETOptionsUseCellList -> True
	}


(* ::Section:: *)
(*Declarations*)


RunMIPET::usage = 
	"RunMIPET[oldParticles, newParticles]";
Off[Decrement::rvalue];
Off[General::munfl];
SetSystemOptions["CheckMachineUnderflow"->False];



(* ::Section:: *)
(*Functions*)


(* TODO_AZ: Complex data structure, i.e. all matrices, must be properly commented otherwise their structure can not be comprehended *)
(* TODO_FB: Suppress popups caused by calling external commands *)

Begin["`Private`"]

(* --------------------------------------------------- *)
(* Private method compilation overview: *)
(* -------------------------------- *)
(* Compiled methods: *)
(* --------------- *)
(* - GetMoleculeIndexByAtomIndex *)
(* - GetNeighborParticleNumbers *)
(* - IsTooClose *)
(* ------------------------ *)
(* Uncompiled methods (compilation not necessary): *)
(* ------------------------ *)
(* - GetNumberExtension *)
(* - ConvertTinkerXYZToXYZ *)
(* - ConvertTinkerXYZToXYZInMemory *)
(* - FormatXYZData *)
(* - GetInterMolecularEnergy *)
(* - GetVdwVolume *)
(* --------------------------------------------------- *)

(* "[... ]all results returned to the master kernel by the slaves are checked for aborted evaluations using MemberQ[res, $Aborted]. *)
(* Here, res is a large matrix in the form of a packed array [some] 160MB in size, and the unpacking of this by MemberQ accounts for the poor performance *)
(* and considerable memory consumption of this example. The peak memory consumption does not persist, however, *)
(* since after the absence of $Aborted has been verified, the intermediate (unpacked) results are discarded." *)
(* - Oleksandr R. on mathematica stackexchange *)
(* https://mathematica.stackexchange.com/questions/2886/transferring-a-large-amount-of-data-in-parallel-calculations, Accessed 10. Juni 2022  *)
withModifiedMemberQ[expr_] :=
 Module[{doneQ, unmatchable},
  Internal`InheritedBlock[{MemberQ},
   Unprotect[MemberQ];
   (* Can uncomment this if we want to print out the MemberQ calls:
   mq:MemberQ[args___]/;(Print@HoldForm[mq];True):=mq;
   *)
   MemberQ[list_, patt_Symbol, args___] /; !TrueQ[doneQ] :=
    Block[
    	{
    		doneQ = True
    	},
     	MemberQ[
      		Unevaluated[list] /. _List?Developer`PackedArrayQ -> {unmatchable},
      		Unevaluated[patt],
      		args
     	]
    ];
   Protect[MemberQ];
   expr
  ]
 ];
SetAttributes[withModifiedMemberQ, HoldAllComplete];

(* Get Van der Waals Volume for given SMILES from CDK *)
(* Note: Compilation is not necessary due to external library call  *)
GetVdwVolume[
	smilesString_,
	cdkJarFullPath_
	] := 
		Block[
			{
				smilesParser,
				particle,
				vabcVolume
			},
			ReinstallJava[ClassPath -> cdkJarFullPath];
			LoadJavaClass["org.openscience.cdk.silent.SilentChemObjectBuilder"];
			LoadJavaClass["org.openscience.cdk.tools.manipulator.AtomContainerManipulator"];
			LoadJavaClass["org.openscience.cdk.geometry.volume.VABCVolume"];
			smilesParser = JavaNew["org.openscience.cdk.smiles.SmilesParser", SilentChemObjectBuilder`getInstance[]];
			particle = smilesParser@parseSmiles[smilesString];
			AtomContainerManipulator`percieveAtomTypesAndConfigureAtoms[particle];
			vabcVolume = VABCVolume`calculate[particle];
			Return[vabcVolume];
		]	
	
(* Get atomic mass for given SMILES from CDK *)
(* Note: Compilation is not necessary due to external library call  *)
GetAtomicMass[
	smilesString_,
	cdkJarFullPath_
	] := 
		Block[
			{
				smilesParser,
				particle,
				atomicMass
			},
			ReinstallJava[ClassPath -> cdkJarFullPath];
			LoadJavaClass["org.openscience.cdk.silent.SilentChemObjectBuilder"];
			LoadJavaClass["org.openscience.cdk.tools.manipulator.AtomContainerManipulator"];						
			smilesParser = JavaNew["org.openscience.cdk.smiles.SmilesParser", SilentChemObjectBuilder`getInstance[]];
			particle = smilesParser@parseSmiles[smilesString];				
			AtomContainerManipulator`percieveAtomTypesAndConfigureAtoms[particle];
			atomicMass = AtomContainerManipulator`getMass[particle, 1];
			Return[atomicMass];
		]	

(* Returns a list with the indices of the neighbouring cells of the specified solubleCellIndex, including solubleCellIndex itself. *)        		
GetNeighborCells[
   solubleCellIndex_,
   cellsInARow_,
   indexCube_
   ] :=
  		Block[
		   {
		    neighborCells,
		    tmpCellIndex
		    },
		   neighborCells = Last[
		      Reap[
		       tmpCellIndex = 1;
		       Do[
		        If[
		         solubleCellIndex == tmpCellIndex,
		         Sow[indexCube[[i, j, k]]];
		         Sow[indexCube[[i, j, k + 1]]];
		         Sow[indexCube[[i, j, k - 1]]];
		         
		         Sow[indexCube[[i, j + 1, k]]];
		         Sow[indexCube[[i, j + 1, k + 1]]];
		         Sow[indexCube[[i, j + 1, k - 1]]];
		         
		         Sow[indexCube[[i, j - 1, k]]];
		         Sow[indexCube[[i, j - 1, k + 1]]];
		         Sow[indexCube[[i, j - 1, k - 1]]];
		         
		         Sow[indexCube[[i + 1, j, k]]];
		         Sow[indexCube[[i + 1, j, k + 1]]];
		         Sow[indexCube[[i + 1, j, k - 1]]];
		         
		         Sow[indexCube[[i + 1, j + 1, k]]];
		         Sow[indexCube[[i + 1, j + 1, k + 1]]];
		         Sow[indexCube[[i + 1, j + 1, k - 1]]];
		         
		         Sow[indexCube[[i + 1, j - 1, k]]];
		         Sow[indexCube[[i + 1, j - 1, k + 1]]];
		         Sow[indexCube[[i + 1, j - 1, k - 1]]];
		         
		         Sow[indexCube[[i - 1, j, k]]];
		         Sow[indexCube[[i - 1, j, k + 1]]];
		         Sow[indexCube[[i - 1, j, k - 1]]];
		         
		         Sow[indexCube[[i - 1, j + 1, k]]];
		         Sow[indexCube[[i - 1, j + 1, k + 1]]];
		         Sow[indexCube[[i - 1, j + 1, k - 1]]];
		         
		         Sow[indexCube[[i - 1, j - 1, k]]];
		         Sow[indexCube[[i - 1, j - 1, k + 1]]];
		         Sow[indexCube[[i - 1, j - 1, k - 1]]];
		         ];
		        tmpCellIndex++,
		        {i, 2, cellsInARow + 1}, {j, 2, cellsInARow + 1}, {k, 2, cellsInARow + 1}
		        ]
		       ]
		      ] // Flatten;
		   Return[neighborCells]
		   ];	   

(* Returns the items of the specified cell *)   
GetCellItems[
   cellListHead_,
   cellListList_,
   cellIndex_
   ] :=
	  Block[
	  	{
	    	cellItems,
	    	pos
	    },
	   	cellItems = 
	    	Flatten[
	     		Last[
	      			Reap[
	       				Sow[pos = cellListHead[[cellIndex]]];
	       				While[
					        pos != 0,
	         				pos = cellListList[[pos]]; 
					        Sow[pos]
	         			]
	       			]
	      		]
	     	];
	   	Return[cellItems]
	];
	
(* Returns the index of the molecule to which the atom with the specified index belongs. *)
(* Deprecated *)
GetMoleculeIndexByAtomIndex = 
	Compile[
		{
			{index, _Integer},
			{nSolventAtomsPerMolecule, _Integer}
		},
		Block[
			{
				atomIndex
			},
			If[
				Mod[index, nSolventAtomsPerMolecule] === 0,
				atomIndex = IntegerPart[index/nSolventAtomsPerMolecule],
				atomIndex = IntegerPart[index/nSolventAtomsPerMolecule] + 1
			];
			atomIndex
		],
		CompilationTarget -> "C",
		RuntimeAttributes -> {Listable}
	];

(* Creates and returns an index cube with dimensions "cellsInARow" + 2  *)
CreateIndexCube[
		cellsInARow_
	] :=
		Block[
			{
				cellIndex,
				indexCube
			},
			cellIndex = 1.;
			indexCube = Table[0., {i, cellsInARow + 2.}, {j, cellsInARow + 2.}, {k, cellsInARow + 2.}];
			Do[
				indexCube[[i, j, k]] = cellIndex++,
  				{i, 2., cellsInARow + 1.}, {j, 2., cellsInARow + 1.}, {k, 2., cellsInARow + 1.}
  			];
  			(* bottom & top face *)
			Do[
 				indexCube[[1, j, k]] = indexCube[[cellsInARow + 1, j, k]];
 				indexCube[[cellsInARow + 2, j, k]] = indexCube[[2, j, k]],
 				{j, 2, cellsInARow + 1}, {k, 2, cellsInARow + 1}
 			];
			(* left & right face *)
			Do[
 				indexCube[[i, j, 1]] = indexCube[[i, j, cellsInARow + 1]];
 				indexCube[[i, j, cellsInARow + 2]] = indexCube[[i, j, 2]],
 				{i, cellsInARow + 2}, {j, 2, cellsInARow + 1}
 			];
			(* front & back face *)
			Do[
				indexCube[[i, 1, k]] = indexCube[[i, cellsInARow + 1, k]];
			 	indexCube[[i, cellsInARow + 2, k]] = indexCube[[i, 2, k]],
			 	{i, cellsInARow + 2}, {k, cellsInARow + 2}
 			];
 			indexCube 					
		];	     
		
(* Creates and returns cell list and a list of all cells that contains elements of the soluble molecule; consists of HEAD and LIST *)
(* Returns {HEAD, LIST, solubleCells} *)
CreateCellList[
		nCellsPerRow_,
		cellI_,
		coordinates_,
		nAllAtoms_,
		nSolubleAtoms_
	]:=
		Block[
			{
				cellListHead,
				cellListList,
				solubleCells,
				icell
			},
			cellListHead = Developer`ToPackedArray[ConstantArray[0, {nCellsPerRow * nCellsPerRow * nCellsPerRow}]];
			cellListList = Developer`ToPackedArray[ConstantArray[0, {nAllAtoms}]];
			solubleCells = Developer`ToPackedArray[ConstantArray[0, {nSolubleAtoms}]];
			Do[
				icell = 1 + IntegerPart[coordinates[[i, 1]] * cellI] + 
					IntegerPart[coordinates[[i, 2]] * cellI] * nCellsPerRow +
					IntegerPart[coordinates[[i, 3]] * cellI] * nCellsPerRow * nCellsPerRow;
				If[
					icell > Length[cellListHead],
					Continue[];
				];
				If[
					i <= nSolubleAtoms,
					solubleCells[[i]] = icell
				];
				cellListList[[i]] = cellListHead[[icell]];
				cellListHead[[icell]] = i,
				{i, Length[coordinates]}
			];
			solubleCells = Select[solubleCells, UnequalTo[0]];
			Return[{cellListHead, cellListList, solubleCells}]		
		];			

(* Determine file extension for xyz-files *)
GetNumberExtension[
	n_?IntegerQ
	] :=
		Block[
			{
				number = n
			},
			If[
				number >= 100,
				Return[ToString[number]],
				Return[StringRepeat["0", 3 - StringLength[ToString[number]]] <> ToString[number]]
			]
		];
		

(* Convert tinker-xyz file to xyz file *)
(* Note: Compilation is not necessary due to file operations *)
ConvertTinkerXYZToXYZ[
	tinkerxyzFileName_, 
	xyzFileName_
	] := 
		Block[
			{
				inputFileName = tinkerxyzFileName,
				outputFileName = xyzFileName,
				atomNumber,
				first3Lines,
				tinkerxyzData,
				xyzDataHeader,
				xyzDataBody, (* [[i, j]] i: atoms, j: xyz-coordinate *)
				xyzData, (* [[i, j]] i: header + atoms, j: xyz-coordinate *)
				stream,
				startRow
			},
			If[
				!FileExistsQ[inputFileName],
				Print["Error: tinker xyz file not found."];
				Return[];
			];
			stream = OpenRead[inputFileName];
			first3Lines = ReadList[stream, Record, 3, RecordLists -> True];
			startRow = 0;
			If[
				First[StringSplit[first3Lines[[2, 1]] ]] == "1",
				startRow = 2,
				startRow = 3
			];
			SetStreamPosition[stream, 0];
			atomNumber = Read[stream, Number];
			SetStreamPosition[stream, 0];
			xyzDataHeader = ReadList[stream, {Number, Word}, startRow - 1, RecordLists -> True];
			xyzDataBody = 
				Table[
					ReadList[stream, {Number, Word, Real, Real, Real, Record}, atomNumber], 
					1
				][[1, All, 2;;5]];
			Close[stream];
			xyzData = Join[xyzDataHeader[[1]], xyzDataBody];
			xyzData = FormatXYZData[xyzData];
			Export[outputFileName, xyzData, "Table"];
		];

(* In-memroy conversion from Tinker-XYZ to XYZ file *)
(* Note: Compilation is not necessary *)
ConvertTinkerXYZToXYZInMemory[
	tinkerXYZ_ (* [[i, j]] i: rows, j: xyz-coordinate *)
	] := 
	Block[
		{
			tinkerXYZData = tinkerXYZ, (* [[i, j]] i: rows, j: xyz-coordinate *)
			xyzHead,
			xyzBody,
			xyzData (* [[i, j]] i: head + atoms, j: xyz-coordinate *)
		},
		xyzHead = tinkerXYZData[[1]];
		xyzBody = tinkerXYZData[[2;;All, 2;;5]];
		xyzData = Insert[xyzBody, xyzHead, 1];
		Return[xyzData];
	];

ConvertXYZDataToXYZ[
	tinkerXYZ_, (* [[i, j]] i: rows, j: xyz-coordinate *)
	solubleAtoms_,
	solventAtoms_,
	atomCount_,
	nIterations_
	] :=
	Block[
		{
			xyzData, (* [[i, j]] i: atoms, j: xyz-coordinate *)			
			i
		},
		xyzData = tinkerXYZ;
			Do[			
				If[
					i - 1 < Length[solubleAtoms],
					xyzData[[i]] = Insert[xyzData[[i]], solubleAtoms[[Mod[(i - 1), Length[solubleAtoms]] + 1]], 1],					
					xyzData[[i]] = Insert[xyzData[[i]], solventAtoms[[Mod[(i - 1 - Length[solubleAtoms]), Length[solventAtoms]] + 1]], 1]
				],
				{i, atomCount}
			];
		xyzData = TableForm[xyzData, TableSpacing -> {0, 4}]//ToString;
		xyzData = ImportString[xyzData, "XYZ"];		
		Return[xyzData];		
	];	

(* Format xyz data *)
(* Note: Compilation is not necessary *)
FormatXYZData[
	xyzData_ (* [[i, j]] i: atoms, j: xyz-coordinate *)
	] := 
		Block[
			{
				xyz = xyzData, (* [[i, j]] i: atoms, j: xyz-coordinate *)
				startRow,
				i
			},
			If[
				StringMatchQ[ToString[xyz[[2, 1]] ], RegularExpression["[A-Z,a-z]{1,2}"]],
				startRow = 2,
				startRow = 3;
			];
			Do[
				xyz[[i, 2]] = ToString[DecimalForm[xyz[[i, 2]], {8, 6}]];
				xyz[[i, 3]] = ToString[DecimalForm[xyz[[i, 3]], {8, 6}]];
				xyz[[i, 4]] = ToString[DecimalForm[xyz[[i, 4]], {8, 6}]],
				{i, startRow, Length[xyz]}
			];
			Return[xyz];
		];

(* Find out whether two atoms are too close *)
(* EuclideanDistance seems to be not compilable. *)
IsTooClose = 
	Compile[
		{
			{atomCoordinates1, _Real, 2}, (* [[i, j]] i: atoms1, j: xyz-coordinates *) 
			{atomCoordinates2, _Real, 2}, (* [[i, j]] i: atoms2, j: xyz-coordinates *)
			{minimumAtomDistance, _Real}
		},
		Block[
			{
				coordinates1 = atomCoordinates1,
				coordinates2 = atomCoordinates2,
				minAtomDistanceSquared = minimumAtomDistance*minimumAtomDistance,
				distance,
				deltaX, deltaY, deltaZ,
				isTooClose = False
			},	
			Do[
				{deltaX, deltaY, deltaZ} = coordinates1[[i]] - coordinates2[[j]];
				distance = {deltaX, deltaY, deltaZ}.{deltaX, deltaY, deltaZ};
				If[
					distance < minAtomDistanceSquared,
					isTooClose = True;
					Break[];					
				],
				{i, Length[coordinates1]}, {j, Length[coordinates2]}
			];
			isTooClose
		],
		CompilationTarget -> "C"
	];

(* Determine the intermolecular energies for the given list of distances and boltzmann weighted the energies with the given boltzmann factor *)
(* Returns a list with all minimum energies for each distances, the overall boltzmann weighted minimum and its corresponding distance *)
(* Note: Compilation is not necessary *)
GetInterMolecularEnergy[
	tinkerAnalyze_,
	distances_,
	nCpuCores_,
	minAtomDistance_,
	keyName_,
	xyzData1_,
	xyzData2_,
	nAtoms1_,
	nAtoms2_,
	rotatedCoordinates1_,
	rotatedCoordinates2_
	] :=
		Block[
			{
				searchPatternEnergy = "Intermolecular Energy :",
				energyData,
				energyDataAtDistance,
				nConfigurations1,
				nConfigurations2,
				nDistances,
				exportData,
				exportData1,
				exportData2,
				rotatedData,
				energyMinimum,
				energyMinimumIndex,
				stream
			},
			nConfigurations1 = Length[rotatedCoordinates1];
			nConfigurations2 = Length[rotatedCoordinates2];
			nDistances = Length[distances];

			exportData1 = xyzData1;
			exportData1[[1, 1]] = nAtoms1 + nAtoms2; (*1st row:number of atoms+name*)
			exportData2 = Drop[xyzData2, 1];
			(*1. column:index of atoms*)
			(*2. column:atom name*)
			(*3-5. column:xyz-coordinates*)
			(*6. column:atom type*)
			(*7+column:connectivity to other atoms*)
			exportData2[[All, 1]] = Table[i, {i, nAtoms1 + 1, nAtoms1 + nAtoms2}];
			exportData2[[All, 7 ;;]] = exportData2[[All, 7 ;;]] + nAtoms1;
			
			DistributeDefinitions[exportData, exportData1, exportData2, nConfigurations1, nConfigurations2, distances, rotatedData, rotatedCoordinates1, rotatedCoordinates2, nAtoms1, nAtoms2, IsTooClose];
			LaunchKernels[nCpuCores];
			ParallelDo[					
				exportData = Table[{}, nConfigurations1, nConfigurations2];
				rotatedData = ((## + {distances[[j]], 0., 0.})&/@#)&/@rotatedCoordinates2;
				Do[
					exportData1[[2;;nAtoms1 + 1, 3;;5]] = rotatedCoordinates1[[k]]; 
					Do[
						If[
							!IsTooClose[rotatedCoordinates1[[k]], rotatedData[[l]], minAtomDistance],
							exportData2[[All, 3;;5]] = rotatedData[[l]];
							exportData[[k,l]] = Join[exportData1, exportData2];
						],
						{l, nConfigurations2}
					],
					{k, nConfigurations1}
				];
				exportData = StringRiffle[Flatten[Join[exportData], 2], "\n", "\t"];
				Export["output." <> ToString[j], exportData, "Text"],
				{j, nDistances}
			];
			CloseKernels[];
			(*calculate intermolecular energy using TINKER analyze*)
			DistributeDefinitions[keyName, tinkerAnalyze];
			LaunchKernels[nCpuCores];
			ParallelTable[
   					ReadList["!" <> tinkerAnalyze <> "output." <> ToString[j] <> " -k " <> keyName <> " " <> "e >output" <> ToString[j] <> ".out"],
   					{j, nDistances}
   			];
			CloseKernels[];
			(*read intermolecular energy*)
			energyData = Table[{}, nDistances];
			energyDataAtDistance = Table[{}, nDistances];
			Do[
			  stream = OpenRead["output" <> ToString[i] <> ".out"];
			  energyData[[i]] = 
			  	Flatten[
			    	StringCases[
			     		FindList[stream, searchPatternEnergy],
			     		x : NumberString :> ToExpression[x]
			     	]
			    ];
			  energyDataAtDistance[[i]] = {distances[[i]], energyData[[i]]};
			  Close[stream],
			  {i, nDistances}
			];
			energyMinimum = 100.;
			energyMinimumIndex = 0;
			Do[
				If[
					Min[energyData[[i]] ] <= energyMinimum,
					energyMinimum = Min[energyData[[i]] ];
					energyMinimumIndex = i
				],
				{i, Length[energyData]}
			];
			Return[{energyMinimum, distances[[energyMinimumIndex]], energyDataAtDistance, energyMinimumIndex}];
		];

(* Returns the indices of the neighbor particles *)
(* Compiled *)
GetNeighborParticleNumbers = 
	Compile[
		{
			{coordinates, _Real, 2},
			{allAtomTypes, _Integer, 1},
			{nSolubleAtoms, _Integer},
			{solventAtomIndices, _Integer, 1},
			{catchRadius, _Real},
			{boxLength, _Real},
			{vdWRadii, _Real, 1}
		},
		Block[
			{
				solubleCoords,
				solventCoords,
				nFirstSolventAtom = nSolubleAtoms + 1,
				distances,
				indices,
				nSolventAtoms ,
				deltaX, deltaY, deltaZ, result, oneHalf = 0.5
			},
			solubleCoords = coordinates[[;; nSolubleAtoms]];
			solventCoords = coordinates[[#]] & /@ solventAtomIndices;
			nSolventAtoms = Length[solventAtomIndices];
			distances = Table[0., {nSolubleAtoms}, {nSolventAtoms}];
			Do[
				{deltaX, deltaY, deltaZ} = solubleCoords[[i]] - solventCoords[[j]];
				If[deltaX > boxLength*oneHalf, deltaX -= boxLength];
				If[deltaX <= -boxLength*oneHalf, deltaX += boxLength];
				If[deltaY > boxLength*oneHalf, deltaY -= boxLength];
				If[deltaY <= -boxLength*oneHalf, deltaY += boxLength];
				If[deltaZ > boxLength*oneHalf, deltaZ -= boxLength];
				If[deltaZ <= -boxLength*oneHalf, deltaZ += boxLength];
			 	distances[[i, j]] = {deltaX, deltaY, deltaZ} . {deltaX, deltaY, deltaZ},
			 	{i, nSolubleAtoms}, {j, nSolventAtoms}
			];
			indices = Table[0, {nSolubleAtoms}, {nSolventAtoms}];
			Do[
				If[
					distances[[i, j]] <= 
					((vdWRadii[[allAtomTypes[[i]] ]] + vdWRadii[[allAtomTypes[[solventAtomIndices[[j]] ]] ]] + catchRadius) *
						(vdWRadii[[allAtomTypes[[i]] ]] + vdWRadii[[allAtomTypes[[solventAtomIndices[[j]] ]] ]] +  catchRadius)
			    	),
			  		indices[[i, j]] = solventAtomIndices[[j]]
			  	],
			 	{i, nSolubleAtoms}, {j, nSolventAtoms}
			];
			indices
		],
		CompilationTarget -> "C"
	];

GetSolubleCoordinates[	
	dataXYZ_,
	shotCount_,
	solubleAtoms_
	] :=
	Block[
		{
			solubleCoords,
			solubleAtomCount,
			i
		},
		solubleAtomCount = Length[solubleAtoms];
		solubleCoords = dataXYZ[[;;solubleAtomCount]];
			
		Do[
			solubleCoords[[i]] = Insert[solubleCoords[[i]], solubleAtoms[[Mod[(i - 1), Length[solubleAtoms]] + 1]], 1],
			{i, solubleAtomCount}
		];
		Return[solubleCoords];	
	];

GetSolventCoordinates[
	dataXYZ_,
	indices_,
	solubleAtoms_,
	solventAtoms_	
	] :=
	Block[
		{
			solventCoords,
			solventAtomCount,
			solubleAtomCount,
			i,
			j
		},
		solventAtomCount = Length[solventAtoms];
		solubleAtomCount = Length[solubleAtoms];
		Do[
			If[
				(1 + solubleAtomCount + (indices[[i]] - 1) * solventAtomCount) > Length[dataXYZ] || (solubleAtomCount + indices[[i]] * solventAtomCount) > Length[dataXYZ],
				x = 5		
			],
			{i, Length[indices]}
		];
		solventCoords = 
			Table[
				dataXYZ[[(1 + solubleAtomCount + (indices[[i]] - 1) * solventAtomCount);;(solubleAtomCount + indices[[i]] * solventAtomCount) ]],
				{i, Length[indices]}
			];
		Do[
			solventCoords[[i, j]] = Insert[solventCoords[[i, j]], solventAtoms[[Mod[(j - 1), solventAtomCount] + 1]], 1],
			{i, Length [solventCoords]}, {j, Length[solventCoords[[i]] ]}
		];
		Return[solventCoords];
	];
	
GetCoordinationNumberSphereApproximation = 
	Compile[
		{
			{vdWVolumei, _Real},
			{vdWVolumej, _Real},
			{factor, _Real}
		},
		Block[
			{
				coordinationNumber,
				numerator
			},
			numerator = (Power[(CubeRoot[vdWVolumei] + 2 * CubeRoot[vdWVolumej]), 3] - vdWVolumei);
			coordinationNumber = factor * numerator / vdWVolumej;
			coordinationNumber
		],
		CompilationTarget -> "C"	
	];	

ExportParameterSet[
	title_,
	titleAbr_,
	forceFieldName_,
	resultsDir_,
	smiles_,
	cdkJarFullPath_,
	waterVolumeRatio_,
	boltzmannFraction_,
	catchRadius_,
	energySearchPattern_,
	parameters_
	]:=
		Block[
			{
				aij,
				densestPackingFator = 0.74,
				randomPackingFactor = 0.64,
				boltzmannConst = 3.297*10^-27, (* kcal/mole *)
				coordinationNumbers,
				coordinationNumber11,
				coordinationNumber12,
				coordinationNumber21,
				coordinationNumber22,
				calculationData,
				chiNumerator,
				defaultCharge,
				denseZij,
				descriptionStringLength,
				file,
				firstLine,
				fragments,
				fragmentIndex,
				graphicsRadius,
				header,
				inputDataStream,
				inputPath,
				inputPathNames,
				intermolecularEnergy11,
				intermolecularEnergy12,
				intermolecularEnergy22,
				massDPD,
				outputStream,
				pairs,
				particle1Row,
				particle2Row,
				particleDescription,
				particleInteractions,
				particleInteractionsWithZijEqualsOne,
				particleInteractionsWithDenseZijSphereApprox,
				particleInteractionsWithRandomZijSphereApprox,
				particleName1,
				particleName2,
				particleSmiles,
				randomZij,
				standardColor,
				temperature,
				volumeParticleI,
				volumeParticleJ,
				versionNr = "1.0.0.0"
			},
			SetDirectory[resultsDir];
			inputPathNames = Select[FileNames["*",""], DirectoryQ];
			fragments = Table["", Length[inputPathNames]*2];
			If[
				Length[inputPathNames] == 0,
				Return[];
			];
			fragmentIndex = 1;
			Do[
				particleName1 = StringSplit[inputPathNames[[i]], "_"][[1]];
				fragments[[fragmentIndex++]] = particleName1;
				particleName2 = StringSplit[inputPathNames[[i]], "_"][[2]];
				fragments[[fragmentIndex++]] = particleName2,
				{i, Length[inputPathNames]}
			];
			fragments = Sort[DeleteDuplicates[fragments]];
			pairs = Sort[Tuples[fragments, 2]];
			(* Get temperatures *)
			inputPath = FileNameJoin[{resultsDir, inputPathNames[[1]]}];
			SetDirectory[inputPath];
			inputDataStream = OpenRead[inputPathNames[[1]] <> ".dat"];
			temperature =
				StringCases[
					Find[
						inputDataStream, "temperature [K] ="
					],
					x : NumberString :> ToExpression[x]
				]; (* in K *)
			Close[inputDataStream];
			(* /Get temperatures *)
			(* particle desrciptions *)
			particleDescription = Table[Null, {Length[fragments] + 1}, 8];
			massDPD = 1; (* set to 1 for all *)
			graphicsRadius = 0.5; (* set to 0.5 for all *)
			defaultCharge = 0; (* check TODO *)
			standardColor = "black";
			particleDescription[[1, 1]] = "# Particle"; (* Particle *)
			particleDescription[[1, 2]] = "Name"; (* Name TODO: get full name*)
			particleDescription[[1, 3]] = "Mass[DPD]"; (* Mass[DPD] *) 
			particleDescription[[1, 4]] = "Charge"; (* Charge *) 
			particleDescription[[1, 5]] = "Mass[g/mol]"; (* Mass[g/mol] *) 
			particleDescription[[1, 6]] = "Volume[A^3]"; (* Volume[A^3] *) 
			particleDescription[[1, 7]] = "Graphics-Radius"; (* Graphics-Radius *) 
			particleDescription[[1, 8]] = "Standard-Color";(* Standard-Color *) 
			Do[
				particleDescription[[i + 1, 1]] = ToString[fragments[[i]]]; (* Particle *)
				particleDescription[[i + 1, 2]] = ToString[fragments[[i]]]; (* Name TODO: get full name*)
				particleDescription[[i + 1, 3]] = ToString[massDPD]; (* Mass[DPD] *) 
				particleDescription[[i + 1, 4]] = ToString[defaultCharge]; (* Charge *) 
				particleDescription[[i + 1, 5]] = ToString[GetAtomicMass[smiles[fragments[[i]] ], cdkJarFullPath ]]; (* Mass[g/mol] *) 
				particleDescription[[i + 1, 6]] = ToString[GetVdwVolume[smiles[fragments[[i]] ], cdkJarFullPath ] * waterVolumeRatio]; (* Volume[A^3] *) 
				particleDescription[[i + 1, 7]] = ToString[graphicsRadius]; (* Graphics-Radius *) 
				particleDescription[[i + 1, 8]] = ToString[standardColor], (* Standard-Color *)
				{i, Length[fragments]}
			];
			(* /particle descriptions *)
			(* particle interactions *)
			particleInteractions = Table[{}, Length[inputPathNames]];
			particleInteractionsWithZijEqualsOne = Table[{}, Length[inputPathNames]];
			particleInteractionsWithDenseZijSphereApprox = Table[{}, Length[inputPathNames]];
			particleInteractionsWithRandomZijSphereApprox = Table[{}, Length[inputPathNames]];
			(* collect data to calculate aij parameters *)
			calculationData = Table[Null, Length[inputPathNames], 4];
			(* 
				{
					{
						particlePair {particle1, particle2}, 
						Eij {{boltzmannFraction1, bf2, bf3, ... for temperature 1}, {boltzmannFraction1, bf2, bf3, ... for temperature 2},,...}, 
						Zij {T1,T2,...}, 
						Zji {T1,T2,...}
					},
				 	...
			   } 
			*)
			Do[
				inputPath = FileNameJoin[{resultsDir, inputPathNames[[i ]] }];
				particleName1 = StringSplit[inputPathNames[[i]], "_"][[1]];
			  	particleName2 = StringSplit[inputPathNames[[i]], "_"][[2]];
			  	calculationData[[i, 1]] = {particleName1, particleName2};
			  	SetDirectory[inputPath];
			  	inputDataStream = OpenRead[inputPathNames[[i]] <> ".dat"];
			  	(* energy values *)
			  	calculationData[[i, 2]] = 
			  		ToExpression[
			  			Flatten /@ Map[
			  				StringSplit[#, ","]&,
			  				StringSplit[
			  					Find[inputDataStream, energySearchPattern], "="
			  				][[2]]
			  			]
			  		];
			  	(* Zij *)	
			  	SetStreamPosition[inputDataStream, 1];
			  	calculationData[[i,3]] = 
			  		ToExpression[
			  			Flatten /@ Map[
			  				StringSplit[#, ","]&,
			  				StringSplit[
			  					Find[
			  						inputDataStream,
			  						StringJoin[{"Mean(",particleName1,"/",particleName2,") ="}]			  						
			  					],
			  					"="
			  				][[2]]
			  			]
			  		];
			  	(* Zji *)
			  	SetStreamPosition[inputDataStream, 1];	
			  	calculationData[[i,4]] = 
			  		ToExpression[
			  			Flatten /@ Map[
			  				StringSplit[#, ","]&,
			  				StringSplit[
			  					Find[
			  						inputDataStream,
			  						StringJoin[{"Mean(",particleName2,"/",particleName1,") ="}]			  						
			  					],
			  					"="
			  				][[2]]
			  			]
			  		];
			  	Close[inputDataStream],			  
				{i, Length[inputPathNames]}
			];
			(* calculate aij - see DPD theory pdf *)
			Do[
				particleName1 = calculationData[[i, 1, 1 ]];
				particleName2 = calculationData[[i, 1, 2 ]];
				particle1Row = Position[calculationData, {particleName1, particleName1}][[1, 1 ]];
				particle2Row = Position[calculationData, {particleName2, particleName2}][[1, 1 ]];
				intermolecularEnergy12 = calculationData[[i, 2 ]];
				intermolecularEnergy11 = calculationData[[particle1Row, 2 ]];
				intermolecularEnergy22 = calculationData[[particle2Row, 2 ]];
				coordinationNumber12 = calculationData[[i, 3]];
				coordinationNumber21 = calculationData[[i, 4]];
				coordinationNumber11 = calculationData[[particle1Row, 3 ]];
				coordinationNumber22 = calculationData[[particle2Row, 3 ]];
				chiNumerator = (
					coordinationNumber12 * intermolecularEnergy12 +
					coordinationNumber21 * intermolecularEnergy12 -
					coordinationNumber11 * intermolecularEnergy11 - 
					coordinationNumber22 * intermolecularEnergy22
				  );
				aij = temperature / 12 + 1.7483 * chiNumerator;			
				particleInteractions[[i]] = {particleName1 <> "_" <> particleName2, ToString[aij[[1]] ]},			
				{i, Length[calculationData]}
			];
			(* calculate aij with Zij = 1 - see DPD theory pdf *)
			Do[
				particleName1 = calculationData[[i, 1, 1 ]];
				particleName2 = calculationData[[i, 1, 2 ]];
				particle1Row = Position[calculationData, {particleName1, particleName1}][[1, 1 ]];
				particle2Row = Position[calculationData, {particleName2, particleName2}][[1, 1 ]];
				intermolecularEnergy12 = calculationData[[i, 2 ]];
				intermolecularEnergy11 = calculationData[[particle1Row, 2 ]];
				intermolecularEnergy22 = calculationData[[particle2Row, 2 ]];
				coordinationNumber12 = 1.;
				coordinationNumber21 = 1.;
				coordinationNumber11 = 1.;
				coordinationNumber22 = 1.;
				chiNumerator = (
					coordinationNumber12 * intermolecularEnergy12 +
					coordinationNumber21 * intermolecularEnergy12 -
					coordinationNumber11 * intermolecularEnergy11 - 
					coordinationNumber22 * intermolecularEnergy22
				  );
				aij = temperature / 12 + 1.7483 * chiNumerator;			
				particleInteractionsWithZijEqualsOne[[i]] = {particleName1 <> "_" <> particleName2, ToString[aij[[1]] ]},			
				{i, Length[calculationData]}
			];
			(* calculate aij with Zij of sphere approx - see DPD theory pdf *) (* TODO sphere approximation *)
			coordinationNumbers = CreateDataStructure["HashTable"];
			Do[
			  	volumeParticleI = 			  		
			  		GetVdwVolume[smiles[pairs[[i, 1]] ], cdkJarFullPath];
			  	volumeParticleJ = 
			  		GetVdwVolume[smiles[pairs[[i, 2]] ], cdkJarFullPath];				  	
			  	denseZij = GetCoordinationNumberSphereApproximation[volumeParticleI, volumeParticleJ, densestPackingFator];	
			  	randomZij = GetCoordinationNumberSphereApproximation[volumeParticleI, volumeParticleJ, randomPackingFactor];	
			  	coordinationNumbers["Insert", {pairs[[i, 1]], pairs[[i, 2]] } -> {denseZij, randomZij}],
				{i, Length[pairs]}
			];
			Do[
				particleName1 = calculationData[[i, 1, 1 ]];
				particleName2 = calculationData[[i, 1, 2 ]];
				particle1Row = Position[calculationData, {particleName1, particleName1}][[1, 1 ]];
				particle2Row = Position[calculationData, {particleName2, particleName2}][[1, 1 ]];
				intermolecularEnergy12 = calculationData[[i, 2 ]];
				intermolecularEnergy11 = calculationData[[particle1Row, 2 ]];
				intermolecularEnergy22 = calculationData[[particle2Row, 2 ]];
				aij = Table[
					coordinationNumber11 = coordinationNumbers["Lookup", {particleName1, particleName1}][[j]];
					coordinationNumber22 = coordinationNumbers["Lookup", {particleName2, particleName2}][[j]];
					coordinationNumber12 = coordinationNumbers["Lookup", {particleName1, particleName2}][[j]];
					coordinationNumber21 = coordinationNumbers["Lookup", {particleName2, particleName1}][[j]];
					chiNumerator = (
						coordinationNumber12 * intermolecularEnergy12 +
						coordinationNumber21 * intermolecularEnergy12 -
						coordinationNumber11 * intermolecularEnergy11 - 
						coordinationNumber22 * intermolecularEnergy22
					);
					temperature / 12 + 1.7483 * chiNumerator,
					{j, 2} (* 2 cause of denseZij and random Zij *)
				];					
				particleInteractionsWithDenseZijSphereApprox[[i]] = {particleName1 <> "_" <> particleName2, ToString[aij[[1, 1]] ]};
				particleInteractionsWithRandomZijSphereApprox[[i]] = {particleName1 <> "_" <> particleName2, ToString[aij[[2, 1]] ]},			
				{i, Length[calculationData]}
			];			
			(* / particle interactions*)		
			(* particle SMILES *)
  			particleSmiles = Table[Null, Length[fragments] + 1, 2];  
  			particleSmiles[[1, 1]] = "# Particle";
  			particleSmiles[[1, 2]] = "SMILES"; 			
  			Do[
  				particleSmiles[[i+1, 1]] = ToString[fragments[[i]] ];
  				particleSmiles[[i+1, 2]] = ToString[smiles[fragments[[i]] ]],
  				{i, Length[fragments]}
  			];
			(* Create and write to file *)
			SetDirectory[resultsDir];				
			file = FileNameJoin[StringJoin[titleAbr, "_EijFraction_", ToString[boltzmannFraction], "_catchRadius_", ToString[catchRadius[[1]] ], ".txt"]];
			WriteParameterSet[
				file,
				title,
				titleAbr,
				temperature,
				parameters,
				particleDescription,
				particleInteractions,
				particleSmiles,
				forceFieldName
			];
			file = FileNameJoin[StringJoin[titleAbr, "_Zij=1_EijFraction_", ToString[boltzmannFraction], "_catchRadius_", ToString[catchRadius[[1]] ], ".txt"]];
			WriteParameterSet[
				file,
				title,
				titleAbr,
				temperature,
				parameters,
				particleDescription,
				particleInteractionsWithZijEqualsOne,
				particleSmiles,
				forceFieldName
			];
			file = FileNameJoin[StringJoin[titleAbr, "_ZijDensestSphere_EijFraction_", ToString[boltzmannFraction], "_catchRadius_", ToString[catchRadius[[1]] ], ".txt"]];
			WriteParameterSet[
				file,
				title,
				titleAbr,
				temperature,
				parameters,
				particleDescription,
				particleInteractionsWithDenseZijSphereApprox,
				particleSmiles,
				forceFieldName
			];
			file = FileNameJoin[StringJoin[titleAbr, "_ZijRandomSphere_EijFraction_", ToString[boltzmannFraction], "_catchRadius_", ToString[catchRadius[[1]] ], ".txt"]];
			WriteParameterSet[
				file,
				title,
				titleAbr,
				temperature,
				parameters,
				particleDescription,
				particleInteractionsWithRandomZijSphereApprox,
				particleSmiles,
				forceFieldName
			];
		];

WriteParameterSet[
	fileName_,
	title_,
	titleAbr_,
	temperature_,
	parameters_,
	particleDescription_,
	particleInteractions_,
	particleSmiles_,
	forceFieldName_
] :=
	Block[
		{
			descriptionStringLength,
			firstLine,
			outputStream
		},
		outputStream = OpenWrite[fileName];
		descriptionStringLength = 53;
		WriteLine[outputStream, "# Particle set for MFSim created by MIPET4Mathematica"];
		WriteLine[outputStream, StringPadRight["# Force Field:", descriptionStringLength] <> ToString[forceFieldName] ];		
		WriteLine[outputStream, StringPadRight["# CPU cores:", descriptionStringLength] <> ToString[parameters[[1]] ]];
		WriteLine[outputStream, StringPadRight["# Temperature:",descriptionStringLength] <> ToString[temperature[[1]] ]];
		WriteLine[outputStream, StringPadRight["# Sphere nodes calculated with Fibonacci algorithm:", descriptionStringLength] <> ToString[parameters[[2]] ]];
		WriteLine[outputStream, StringPadRight["# Sphere node number:", descriptionStringLength] <> ToString[parameters[[3]] ]];
		WriteLine[outputStream, StringPadRight["# Sphere rotation number:", descriptionStringLength] <> ToString[parameters[[4]]] ];
		WriteLine[outputStream, StringPadRight["# Fraction for Boltzmann averaging:", descriptionStringLength] <> ToString[parameters[[5]] ]];
		WriteLine[outputStream, StringPadRight["# Used cell index method for MD simulation analysis:", descriptionStringLength] <> ToString[parameters[[6]] ]];
		WriteLine[outputStream, StringPadRight["# Number of solvent molecules:", descriptionStringLength] <> ToString[parameters[[7]] ]];
		WriteLine[outputStream, StringPadRight["# MD step number for warm up:", descriptionStringLength] <> ToString[parameters[[8]] ]];
		WriteLine[outputStream, StringPadRight["# MD step number:", descriptionStringLength] <> ToString[parameters[[9]] ]];
		WriteLine[outputStream, StringPadRight["# Catch radius for MD simulation analysis:", descriptionStringLength] <> ToString[parameters[[10, 1]] ]];
		WriteLine[outputStream, "\n"];
		(* Write particle Description *)
		WriteLine[outputStream, "[Particle description]"];
		Do[
			WriteLine[outputStream, StringRiffle[StringPadRight[#, 17] & /@ particleDescription[[k]] ]],				
			{k, Length[particleDescription]}
		];
		WriteLine[outputStream, "[/Particle description]\n"];
		(*Write particle interactions*)
		WriteLine[outputStream, "[Particle interactions]\n" <> "# Repulsion parameters a(ij) for different temperatures (in K)"];
		firstLine = {"Pair", ToString[temperature[[1]] ]};
		WriteLine[outputStream, StringRiffle[StringPadRight[#, 20]&/@firstLine]];
		Do[
			WriteLine[outputStream, StringRiffle[StringPadRight[#, 20]&/@particleInteractions[[k]] ]],
			{k, Length[particleInteractions]}
		];
		WriteLine[outputStream, "[/Particle interactions]\n"];
		(* Write particle SMILES *)		
		WriteLine[outputStream, "[Particle SMILES]"];
		Do[
			WriteLine[outputStream, StringRiffle[StringPadRight[#, 20]&/@particleSmiles[[k]], "\t" ]],
			{k, Length[particleSmiles]}
		];
		WriteLine[outputStream, "[/Particle SMILES]"];
		Close[outputStream];		
	];

ExportParticleSetSphereApproximantion[
	title_,
	titleAbr_,
	resultsDir_,
	smiles_,
	waterVolumeRatio_,
	cdkJarFullPath_,
	boltzmannFraction_,
	energySearchPattern_
	]:=
		Block[
			{
				densestPackingFator = 0.74,
				randomPackingFactor = 0.64,
				energyDataTable,
				energyDataTableHeader,
				chiNumerator,
				coordinationNumbers,
				coordinationNumbersTable,
				coordination11,
				coordination12,
				coordination21,
				coordination22,
				particleInteractions,
				parameterSet,
				inputDataStream,
				inputPathNames,
				inputPath,
				fragments,
				fragmentIndex,
				pairs,
				keys,
				particleName1,
				particleName2,
				volumeParticleI,
				volumeParticleJ,
				denseZij,
				randomZij,
				energyData,
				energies,
				energy11,
				energy12,
				energy22,
				aij,
				i,
				j,
				k
			},
			SetDirectory[resultsDir];
			inputPathNames = Select[FileNames["*",""], DirectoryQ];
			fragments = Table["", Length[inputPathNames]*2];
			If[
				Length[inputPathNames] == 0,
				Return[];
			];
			fragmentIndex = 1;
			Do[
				particleName1 = StringSplit[inputPathNames[[i]], "_"][[1]];
				fragments[[fragmentIndex++]] = particleName1;
				particleName2 = StringSplit[inputPathNames[[i]], "_"][[2]];
				fragments[[fragmentIndex++]] = particleName2,
				{i, Length[inputPathNames]}
			];
			fragments = Sort[DeleteDuplicates[fragments]];
			pairs = Sort[Tuples[fragments, 2]];
			coordinationNumbers = CreateDataStructure["HashTable"];
			Do[
			  	volumeParticleI = 			  		
			  		GetVdwVolume[smiles[pairs[[i, 1]] ], cdkJarFullPath];
			  	volumeParticleJ = 
			  		GetVdwVolume[smiles[pairs[[i, 2]] ], cdkJarFullPath];				  	
			  	denseZij = GetCoordinationNumberSphereApproximation[volumeParticleI, volumeParticleJ, densestPackingFator];	
			  	randomZij = GetCoordinationNumberSphereApproximation[volumeParticleI, volumeParticleJ, randomPackingFactor];	
			  	coordinationNumbers["Insert", {pairs[[i, 1]], pairs[[i, 2]] } -> {denseZij, randomZij}],
				{i, Length[pairs]}
			];
			energyData = CreateDataStructure["HashTable"];
			Do[
				inputPath = FileNameJoin[{resultsDir, inputPathNames[[i]] }];
				particleName1 = StringSplit[inputPathNames[[i]], "_"][[1]];
			  	particleName2 = StringSplit[inputPathNames[[i]], "_"][[2]];		  
			  	SetDirectory[inputPath];
			  	inputDataStream = OpenRead[inputPathNames[[i]] <> ".dat"];
			  	(* energy values *)
			  	energies =
  			  		ToExpression[
   			  			Map[
     			  				StringSplit[#, ","] &,
     			  				StringSplit[
     			  					Find[
     			  						inputDataStream, 
        								energySearchPattern
        							],
        							 "="
       			  				][[2]]
     			  		]
   			  		];
   			  	energyData["Insert", {particleName1, particleName2} -> energies];
   			  	If[
   			  		particleName1 != particleName2,
   			  		energyData["Insert", {particleName2, particleName1} -> energies];
   			  	];	
			  	Close[inputDataStream],
				{i, Length[inputPathNames]}
			];
			particleInteractions = Table[Null, Length[inputPathNames], 2];
			Do[
				particleName1 = StringSplit[inputPathNames[[i]], "_"][[1]];
			  	particleName2 = StringSplit[inputPathNames[[i]], "_"][[2]];	
				Do[
					energy11 = energyData["Lookup", {particleName1, particleName1}];
					energy22 = energyData["Lookup", {particleName2, particleName2}];
					energy12 = energyData["Lookup", {particleName1, particleName2}]; 
					coordination11 = coordinationNumbers["Lookup", {particleName1, particleName1}][[j]];
					coordination22 = coordinationNumbers["Lookup", {particleName2, particleName2}][[j]];
					coordination12 = coordinationNumbers["Lookup", {particleName1, particleName2}][[j]];
					coordination21 = coordinationNumbers["Lookup", {particleName2, particleName1}][[j]];
					chiNumerator = (coordination12*energy12 + coordination21*energy12 - coordination11*energy11 - coordination22*energy22);
					aij = 300/12 + 1.7483*chiNumerator;
					particleInteractions[[i,j]] = aij,
					{j, 2}
				];
				particleInteractions[[i]] = Insert[particleInteractions[[i]], particleName2, 1];
				particleInteractions[[i]] = Insert[particleInteractions[[i]], particleName1, 1];
				particleInteractions[[i]] = Flatten[particleInteractions[[i]] ],
				{i, Length[inputPathNames]}
			];
			particleInteractions = Sort[particleInteractions];
			keys = Sort[coordinationNumbers["Keys"]];			
			coordinationNumbersTable = 
				Table[					
					{
						keys[[i, 1]],
						keys[[i, 2]],
						coordinationNumbers["Lookup", keys[[i]] ]
					},
					{i, Length[keys]}
				];		
			coordinationNumbersTable = Insert[coordinationNumbersTable, {"particle i", "particle j", "densest packing", "random packing"}, 1];
			Do[
  				coordinationNumbersTable[[i]] = Flatten[coordinationNumbersTable[[i]] ],
  				{i, Length[coordinationNumbersTable]}
  			];
  			energyDataTableHeader = {boltzmannFraction};
			energyDataTableHeader = Insert[energyDataTableHeader, "particle j", 1];
			energyDataTableHeader = Insert[energyDataTableHeader, "particle i", 1];	
  			keys = energyData["Keys"];
			energyDataTable = 
				Table[					
					{
						keys[[i, 1]],
						keys[[i, 2]],
						energyData["Lookup", keys[[i]] ]
					},
					{i, Length[keys]}
				];		
			energyDataTable = Join[{energyDataTableHeader}, Sort[energyDataTable]];
			Do[
  				energyDataTable[[i]] = Flatten[energyDataTable[[i]] ],
  				{i, Length[energyDataTable]}
  			];			
			parameterSet = coordinationNumbersTable;
			parameterSet = Insert[parameterSet, "coordination numbers", 1];
			parameterSet = Insert[parameterSet, "energyData", -1];			
			Do[
				parameterSet = Insert[parameterSet, energyDataTable[[j]], -1],
				{j, Length[energyDataTable]}
			];	
			parameterSet = Insert[parameterSet, "particle interactions", -1];
			parameterSet = 
					Insert[parameterSet, StringJoin["interaction boltzmann fraction ", ToString[boltzmannFraction]], -1];
				Do[
					parameterSet =
					Insert[
						parameterSet,
						Flatten[{
							particleInteractions[[j, 1 ;; 2]],
							particleInteractions[[j, 2 ;; 3]]
						}],
						-1
					],
					{j, Length[particleInteractions]}
				];
			Export[FileNameJoin[{resultsDir,"parameterSet.txt"}], parameterSet, "Table"];					
		];

RunMIPET[
	(* Calculated particles *)
	oParticles_,	
	(* New particles *)
	nParticles_,
	(* Scratch directory *)
	scratchDir_,	
	(* Options *)
	opts___
	] := 
		Block[
			{
				oldParticles = oParticles,
				newParticles = nParticles,
				scratchDirectory = scratchDir,
				adjustedCoordinate1,
				adjustedCoordinate2,
				nAllAtoms,
				allAtomsNames,
				allDistances,
				allParticles, (* list of all particles, union of newParticles and oldParticles *)
				allParticlePairs, (* //tuples of length 2 of allParticles *)
				allSolventAtomsIndices,
				allSolventAtomsNames,
				nSolubleAtoms,
				nSolventAtoms,
				averageCharts,
				avogadroConst,
				boltzmannFraction,
				boxLength,
				boxLengthHalf,
				boxMinimizationTime,
				calculationDirectory,
				catchRadius,
				cdkJarFullPath,
				cdkJarName,
				cellI,
				cellLength,
				cellListCellsPerRow,
				cellListHead,
				cellListIndexCube,
				cellListList,
				cellListMaxDistances,
				cellListNeighborCells,
				cellListNeighborCellsItems,
				cellListNeighborCellsSolventItems,
				cellListResult,
				cellListSolubleCellsList,
				centerCoordinate1,
				centerCoordinate2,
				combinedNb,
				configurationNumber1,
				configurationNumber2,
				conformationSearchTime,
				coordinateData1,
				coordinateData2,
				coordinationCalculationTime,
				cpuCoreNumber,
				createOverview,
				dataName,
				dielectricConstant,
				directoryName,
				distances,
				dynamicAnalysisTime,
				dynamicWarmUpTime,
				elementList,
				elementSymbolToAtomicNumberTable,
				energyCalculationTime,
				energyDataFractioned,	
				energyDataHistogramTable,
				energyDataSorted,
				energyDataTable,
				energyListLinePlots,
				energyMinimumIndex,
				energyMinimumValueLines,
				energyMinimumXyzData,
				energyPreciseListLinePlots,		
				energyThreshold,
				exportCharts,
				exportZijShow,
				forceField,
				fractionToMax,
				gasConstant,
				histograms,
				importData1,
				importData2,
				index,
				indexShift,
				inputDirectory,
				inputPath,
				inputdirectoryName,
				inputdirectoryNamePosition, 
				inputdirectoryNameString, (* String representation of particle pair for the inputdirectory, whereas e.g Et_Me and Me_Et lead to the same dir *) 
				interMolecularEnergyResult,
				isSameParticle,
				isScratchDirectoryCreated,
				keyFileString,
				keyFileStringOrigin,				
				keyName,
				librariesDirectory,
				lineIndex,
				linesToSkip,
				lowerBoundary,
				minAtomDistance,
				minDistance,
				minDistanceGlobal,
				minEnergy,
				minEnergyGlobal,
				minimizeMaxIteration,
				minimumList,
				minimumList1,
				minimumList2,
				preciseStepSize,
				morePreciseXyData,
				moleculeLookupTable,
				nDec,
				nDynamicIteration,
				nDynamicSteps,
				nDynamicWarmupIteration,
				nEntryLines,
				nIterations,
				nRun,
				nRunIteration,
				nSearchDirection,
				nAllSolventAtoms,
				nSolventMolecules,
				nStepsPerRound,
				nTotalSolventatoms,
				nbFile,
				nbFilesList,
				neighborIndices,
				neighborParticleMean,
				neighborParticleNumbers,
				newXYZData,
				notebook,
				numberAtoms1,
				numberAtoms2,
				oldParticlePairs,
				optimizedData,
				optimizedMinimum,
				optXYZDirectory,
				optXYZdirectoryNames,
				optimizeRmsGradient,
				outputData,
				outputDirectory,
				outputDirectoryNames,
				outputPath,
				parameterFileTitle,
				parameterFileTitleAbr,
				particleLogFile,
				particleName,
				particlePair,
				particlePairString, (* String representation of particle pair (e.g. "Et_Et", whereas "Et_Me" is not equal to "Me_Et") *)
				particlePairTime,
				particles,
				phi,
				preciseXyData,
				preciseXyLinePlots,
				precisionGoal,
				prescanStepSize,
				scanStepSize,
				preciseScanStepSize,
				printInterval,
				radiusSphereCoordinates,
				randomseed,
				rawDynamicTime,
				resultDirectory,
				rmsGradient,
				rmsMinimizeGradient,
				rotatedCoordinate1,
				rotatedCoordinate2,
				rotationCalculationTime,
				rotationMatrices1,
				rotationMatrices2,
				rotationNumber,
				scanProgram,
				scratchDirectoryParticlePair,
				simulationType,
				smiles,
				smilesData,
				smilesDirectory,
				solubleAtomNames,
				solubleAtomNamesString,
				solubleCoordsFirst,
				solubleCoordsFirstDetail,
				solubleCoordsLast,
				solubleCoordsLastDetail,
				solubleData,
				solubleName,
				solventAtomNames,
				solventAtomNamesString,
				solventAtomsStartIndex,
				solventCoordsFirst,
				solventCoordsFirstDetail,
				solventCoordsLast,
				solventCoordsLastDetail,
				solventData,
				solventMoleculeNumber,
				solventName,
				sourceDirectory,
				sphereNodeCoordinates,
				sphereNodeNumber,
				startRow,
				stepSize,
				stream,
				tempNb,
				temperature,
				theta,
				timeStep,
				tinkerAnalyze,
				tinkerDynamic,
				tinkerMinimize,
				tinkerOptimize,
				tinkerScan,
				tinkerXYZEdit,
				tinkerXYZpdb,
				totalHistogram,
				totalTime,
				upperBoundary,
				useCellList,
				useFibonacciSphereAlgorithm,
				vdwSolventVolume,
				vdwVolume1,
				vdwVolume2,
				vdWRadii,
				warmUpNeighborIndices,
				warmUpPrintInterval,
				warmUpStepNumber,
				warmUpTimeStep,
				waterVolumeRatio,
				weight,
				weightedAverageList,
				weightedMinEnergyTextString = "(Weighted) Minimum Intermolecular Energy [kcal/mole] = ",
				xyData,
				xyzData1,
				xyzData2,
				xyzDataList,
				xyzDataListInput,				
				xyzDataWarmUpCoords,
				xyzDataWarmUpCoordsInput,
				xyDataWeighted,
				xyzFirst,
				xyzLast,
				xyLinePlots,
				xyzNames,
				xyzTemp,
				xSphereCoordinate,
				ySphereCoordinate,
				zSphereCoordinate
			},			
			(* Check parameters *)
			(* Create scratch directory if it does not already exist *)
			If[
				MatchQ[scratchDirectory, ""],
				Print["No valid scratch directory was specified."];
				Return[]
			];
			If[
				DirectoryQ[scratchDirectory],
				DeleteDirectory[scratchDirectory, DeleteContents -> True]
			];			
			If[
				!DirectoryQ[scratchDirectory], 
				isScratchDirectoryCreated = CreateDirectory[scratchDirectory];
				If[
					MatchQ[ToString[isScratchDirectoryCreated], "$Failed"],
					Print["\"" <> scratchDirectory <> "\"" <> " is not a valid path for the scratch directory."];
					Return[];
				]
			];
			(* Check if newParticles is not empty *)
			If[
				Length[newParticles] <= 0,
				Print["No new particle was specified"];
				Return[]
			];				
			(* Start measurement of total calculation time *)
			totalTime = Now;	
			(* Initialization - Options that can be adjusted in the Start Notebook are set here. In addition, paths and constants are set here. *)
			cpuCoreNumber = MIPETOptionCPUCoreNumber/.{opts}/.Options[MIPETOptions];
			temperature = MIPETOptionTemperature/.{opts}/.Options[MIPETOptions];			
			sphereNodeNumber = MIPETOptionSphereNodeNumber/.{opts}/.Options[MIPETOptions];
			rotationNumber = MIPETOptionRotationNumber/.{opts}/.Options[MIPETOptions];
			minAtomDistance = MIPETOptionMinAtomDistance/.{opts}/.Options[MIPETOptions];
			scanProgram = MIPETOptionScanProgram/.{opts}/.Options[MIPETOptions]; 
			nSearchDirection = MIPETOptionNumberSearchDirection/.{opts}/.Options[MIPETOptions];
			energyThreshold = MIPETOptionEnergyThreshold/.{opts}/.Options[MIPETOptions];
			rmsGradient = MIPETOptionRmsGradient/.{opts}/.Options[MIPETOptions]; 
			rmsMinimizeGradient = MIPETOptionBfgsRmsGradient/.{opts}/.Options[MIPETOptions];
			minimizeMaxIteration = MIPETOptionsBfgsMaxIteration/.{opts}/.Options[MIPETOptions];
			forceField = MIPETOptionForceFieldName/.{opts}/.Options[MIPETOptions];		
			keyFileStringOrigin = 
				StringJoin[{"# Force Field Selection\nPARAMETERS\t\"", NotebookDirectory[], "Tinker\\params\\", ToLowerCase[forceField], ".prm\"\n"}];
			dielectricConstant = MIPETOptionDielectricConstant/.{opts}/.Options[MIPETOptions];
			smilesDirectory = FileNameJoin[{NotebookDirectory[], "Molecules", "SMILES" }];
			calculationDirectory = FileNameJoin[{NotebookDirectory[], "Calculation"}];
			cdkJarName = MIPETOptionJarFileName/.{opts}/.Options[MIPETOptions];
			librariesDirectory = FileNameJoin[{NotebookDirectory[], "Libraries"}];			
			sourceDirectory = FileNameJoin[{NotebookDirectory[], "Molecules", forceField}];
			inputDirectory = FileNameJoin[{calculationDirectory, "Input"}];
			optXYZDirectory = FileNameJoin[{calculationDirectory, "OptXYZ", forceField}];
			outputDirectory = FileNameJoin[{ParentDirectory[NotebookDirectory[]], "Results"}];
			resultDirectory = FileNameJoin[{outputDirectory, forceField}];
			lowerBoundary = MIPETOptionLowerBoundary/.{opts}/.Options[MIPETOptions];
			upperBoundary = MIPETOptionUpperBoundary/.{opts}/.Options[MIPETOptions];
			prescanStepSize = MIPETOptionPrescanStepSize/.{opts}/.Options[MIPETOptions];
			stepSize = MIPETOptionStepSize/.{opts}/.Options[MIPETOptions];
			preciseStepSize = MIPETOptionPreciseStepSize/.{opts}/.Options[MIPETOptions];
			precisionGoal = 3;
			isSameParticle = False;
			avogadroConst = 6.022*10^23; (* 1/mol *)
			gasConstant = 1.9872*10^-3;(*gas constant kB=[kcal/(mol*K)]*)
			solventMoleculeNumber = MIPETOptionSolventMoleculeNumber/.{opts}/.Options[MIPETOptions];			
			nDynamicSteps =  MIPETOptionStepNumber/.{opts}/.Options[MIPETOptions];
			timeStep = MIPETOptionTimeStep/.{opts}/.Options[MIPETOptions];
			printInterval = timeStep*10^(-3);			
			simulationType = MIPETOptionSimulationType/.{opts}/.Options[MIPETOptions];
			catchRadius = MIPETOptionCatchRadius/.{opts}/.Options[MIPETOptions];
			optimizeRmsGradient = MIPETOptionOptimizeRmsGradient/.{opts}/.Options[MIPETOptions];			
			warmUpStepNumber = MIPETOptionWarmUpStepNumber/.{opts}/.Options[MIPETOptions];
			warmUpTimeStep = 1.0;
			warmUpPrintInterval = N[warmUpStepNumber * 10^(-3)];			
			nDynamicWarmupIteration = Ceiling[(warmUpStepNumber * warmUpTimeStep * 10. ^(-15))/(warmUpPrintInterval * 10.^(-12))]; (* timeStep in fs and printIntervall in ps*)
			nDynamicIteration = Ceiling[(nDynamicSteps * timeStep * 10. ^(-15))/(printInterval * 10.^(-12))]; (* timeStep in fs and printIntervall in ps*)		
			tinkerAnalyze = FileNameJoin[{NotebookDirectory[], "Tinker", "analyze.exe "}];
			tinkerScan = FileNameJoin[{NotebookDirectory[], "Tinker", "scan.exe "}];
			tinkerXYZEdit = FileNameJoin[{NotebookDirectory[], "Tinker", "xyzedit.exe "}];
			tinkerMinimize = FileNameJoin[{NotebookDirectory[], "Tinker", "minimize.exe " }];
			tinkerDynamic = FileNameJoin[{NotebookDirectory[], "Tinker", "dynamic.exe "}];
			tinkerOptimize = FileNameJoin[{NotebookDirectory[], "Tinker", "optimize.exe "}];
			tinkerXYZpdb = FileNameJoin[{NotebookDirectory[], "Tinker", "xyzpdb.exe "}];
			nbFilesList = {};
			createOverview = MIPETOptionCreateOverview/.{opts}/.Options[MIPETOptions];
			exportCharts = MIPETOptionExportCharts/.{opts}/.Options[MIPETOptions];
			parameterFileTitle = MIPETOptionParameterSetTitle/.{opts}/.Options[MIPETOptions];
			parameterFileTitleAbr = MIPETOptionParameterSetTitleAbr/.{opts}/.Options[MIPETOptions];
			boltzmannFraction = MIPETOptionsBoltzmannFraction/.{opts}/.Options[MIPETOptions];
			exportZijShow = MIPETOptionsExportZijShow/.{opts}/.Options[MIPETOptions];
			useCellList = MIPETOptionsUseCellList/.{opts}/.Options[MIPETOptions];
			nStepsPerRound = MIPETOptionsStepsPerRound/.{opts}/.Options[MIPETOptions];
			randomseed = 123456789;
			useFibonacciSphereAlgorithm = MIPETOptionsUseFibonacciSphereAlgorithm/.{opts}/.Options[MIPETOptions];
			
			(* Load smiles data - Loads the SMILES code of the "known" molecules*)
			SetDirectory[smilesDirectory];
			smilesData = Import["Smiles.dat", "Table"];	
			(* Puts SMILES into a "table" *)
			Do[
				smiles[smilesData[[i, 1]] ] = smilesData[[i, 2]],
				{i, Length[smilesData]}
			];

			(* Load CDK and get vdw radii of H, C, N, O, P, S *)
			(* "The class path is the set of directories in which the Java runtime looks for classes." -Wolfram Language & System Documentation Center *)
			cdkJarFullPath = FileNameJoin[{ librariesDirectory, cdkJarName}];
			ReinstallJava[ClassPath -> cdkJarFullPath];
			LoadJavaClass["org.openscience.cdk.tools.periodictable.PeriodicTable"];
			(* 
			- Use integer codes so that vdWRadii becomes a pure integer (!) array \[Rule] Number of protons/atomic number
			- Make conversion array for atomic number to element symbol string for data interplay
			Specifiy for all atoms of periodic table, e.g. (please check syntax)
			*)
			vdWRadii = Table[0.0, {i, 92}];
			Do[
				elementSymbolToAtomicNumberTable[PeriodicTable`getSymbol[i]] = i;
				vdWRadii[[i]] = PeriodicTable`getVdwRadius[PeriodicTable`getSymbol[i]],
				{i, Length[vdWRadii]}
			];			
			vdWRadii = vdWRadii/.{Null -> 0.0};
			(* Calculate water volume ratio - ratio of Vparticle and Vvdw of water*)
            waterVolumeRatio = (30.0036 / GetVdwVolume["O", cdkJarFullPath]); 
			(* Prepair jobs *)
			(* Create and set result directory and select contained results*)
			If[
				!DirectoryQ[resultDirectory], 
				CreateDirectory[resultDirectory]
			];
			SetDirectory[resultDirectory];
			outputDirectoryNames = Select[FileNames["*", "", Infinity], DirectoryQ];

			(* Create log file - Logfile for each particle pair*)
			gblLogFileName = FileNameJoin[{resultDirectory, "log.txt"}];
			If[
				!FileExistsQ[gblLogFileName],
				CreateFile[gblLogFileName];
				gblLogFile = OpenAppend[gblLogFileName],
				gblLogFile = OpenAppend[gblLogFileName];
				WriteString[gblLogFile, "\n"];
			];
			WriteString[gblLogFile, StringForm["CPU cores: `` \n", ToString[cpuCoreNumber]]];	
			WriteString[gblLogFile, StringForm["Sphere node number: `` \n", ToString[sphereNodeNumber]]];
			WriteString[gblLogFile, StringForm["Sphere rotation number: `` \n", ToString[rotationNumber]]];
			WriteString[gblLogFile, StringForm["Warm up steps for dynamic simulation: `` \n", ToString[warmUpStepNumber]]]; 
			WriteString[gblLogFile, StringForm["Steps for dynamic simulation: `` \n", ToString[nDynamicSteps]]];
			
			(* Sort particles lexicographically *)
			newParticles = LexicographicSort[newParticles];
			oldParticles = LexicographicSort[oldParticles];
			(* Create input folders for each particle pair if it does not yet exist and copies xyz file to folder *)
			If[
				DirectoryQ[inputDirectory],
				DeleteDirectory[inputDirectory, DeleteContents -> True]
			];	
			If[
				!DirectoryQ[inputDirectory],
				CreateDirectory[inputDirectory]
			];
			SetDirectory[inputDirectory];
			Do[
				directoryName = newParticles[[i]] <> "_" <> oldParticles[[j]];
				If[
					!DirectoryQ[directoryName] && !MemberQ[outputDirectoryNames, directoryName],
					CreateDirectory[directoryName]
				];
				CopyFile[
					FileNameJoin[{sourceDirectory, newParticles[[i]] <> ".xyz"}],
					FileNameJoin[{inputDirectory, directoryName, newParticles[[i]] <> ".xyz"}], 
					OverwriteTarget -> True
				];
				CopyFile[
					FileNameJoin[{sourceDirectory, oldParticles[[j]] <> ".xyz"}],
					FileNameJoin[{inputDirectory, directoryName, oldParticles[[j]] <> ".xyz"}], 
					OverwriteTarget -> True
				],
				{i, Length[newParticles]}, {j, Length[oldParticles]}
			];
			Do[
				directoryName = newParticles[[i]] <>"_" <> newParticles[[i]];
				If[
					!DirectoryQ[directoryName] && !MemberQ[outputDirectoryNames, directoryName],
					CreateDirectory[directoryName]
				];
				CopyFile[
					FileNameJoin[{sourceDirectory, newParticles[[i]] <> ".xyz"}],
					FileNameJoin[{inputDirectory, directoryName, newParticles[[i]] <> ".xyz"}],
					OverwriteTarget ->True],
					{i, Length[newParticles]}
			];
			Do[
				Do[
					directoryName = newParticles[[i]] <> "_" <> newParticles[[j]];
					If[
						!DirectoryQ[directoryName] && !MemberQ[outputDirectoryNames, directoryName],
						CreateDirectory[directoryName]
					];
					CopyFile[
						FileNameJoin[{sourceDirectory, newParticles[[i]] <> ".xyz"}],
						FileNameJoin[{inputDirectory, directoryName, newParticles[[i]] <> ".xyz"}],
						OverwriteTarget -> True
					];
					CopyFile[
						FileNameJoin[{sourceDirectory, newParticles[[j]] <> ".xyz"}],
						FileNameJoin[{inputDirectory, directoryName, newParticles[[j]] <> ".xyz"}],
						OverwriteTarget -> True
					],
					{j, i - 1}
				],
				{i, 2, Length[newParticles]}
			];
			(* Load surface coordinates *)
			(* The coordinates for equidistantly distributed points on a sphere from Technical University of Dortmund are used, thanks to J. Fliege and U. Maier *)
			(* http://www.mathematik.uni-dortmund.de/lsx/research/projects/fliege/nodes/nodes.html *)
			rotationCalculationTime = Now;
			If[
				useFibonacciSphereAlgorithm,
				(* TRUE *)
				(* Fibonacci sphere algorithm is used for determination of the *) 
				(* golden angle in radians *)
				phi = Pi * (3. - Sqrt[5.]);
				nDec = sphereNodeNumber - 1.;
				index = Range[0., nDec];
				ySphereCoordinate = N[1. - 2.  index/nDec];
				radiusSphereCoordinates = Sqrt[1. - ySphereCoordinate * ySphereCoordinate];
				theta = phi * index;
				xSphereCoordinate = Cos[theta] * radiusSphereCoordinates;
				zSphereCoordinate = Sin[theta] * radiusSphereCoordinates;
				sphereNodeCoordinates = Transpose[{xSphereCoordinate, ySphereCoordinate, zSphereCoordinate}],
				(* FALSE *)
				(* sphere node coordinates from Fliege et. al. will be used *)					
				SetDirectory[FileNameJoin[{librariesDirectory, "SphereNodes"}]];
				sphereNodeCoordinates = ToExpression[Import["SphereNodes" <> ToString[sphereNodeNumber] <> ".txt", "List"]];
			];			
			configurationNumber1 = sphereNodeNumber;
			configurationNumber2 = sphereNodeNumber * rotationNumber;
			rotationMatrices1 = Table[0., configurationNumber1];
			rotationMatrices2 = Table[0., configurationNumber2];
			(* Calculate rotation matrices used to rotate the particle/atom coordinates *)
			Do[
				rotationMatrices1[[i]] = RotationMatrix[{sphereNodeCoordinates[[i]], {1., 0., 0.}}],
				{i, configurationNumber1}
			];
			iIndex = 0;
			Do[
				Do[
					iIndex++;
					rotationMatrices2[[iIndex]] = 
						RotationMatrix[2. * Pi * j / rotationNumber, {-1., 0., 0.}] . RotationMatrix[{sphereNodeCoordinates[[i]], {-1., 0., 0.}}],
					{j, rotationNumber}
				],
				{i, configurationNumber1}
			];
			WriteString[gblLogFile, StringForm["Time for calculation of rotation matrices: `` \n", ToString[Now - rotationCalculationTime]]];
			rotationCalculationTime =.;
			(* Start of the calculation *)
			SetDirectory[inputDirectory];
			inputdirectoryName = Select[FileNames["*", "", Infinity], DirectoryQ];		
			gblEnergyCalculationTime=Now;
			outputData = Table[Null, Length[inputdirectoryName], 24]; (* 24 := magic number, do not touch *)
			energyListLinePlots = Table[Null, Length[inputdirectoryName]];
			energyPreciseListLinePlots = Table[Null, Length[inputdirectoryName]];
			Do[		
				particlePairTime = Now;
				(* Splits inputdirectoryName to retrieve the names of both particles of the pair *)
				Do[
					particles = StringSplit[inputdirectoryName[[1]], "_"];
					If[
						particles[[1]] != particles[[2]],
						inputdirectoryName = RotateLeft[inputdirectoryName, 1],
						Break[];
					],
					{i, Length[inputdirectoryName]}
				];
				(* Create key file for Tinker analyze TODO: add water model specification here *)
				particlePair = inputdirectoryName[[f]];
				inputPath = FileNameJoin[{inputDirectory, particlePair}];
				SetDirectory[inputPath];
				xyzNames = FileNames["*.xyz"];
				keyFileString = 
					keyFileStringOrigin <>
					"DIELECTRIC\t" <> ToString[dielectricConstant];			
				Export[particlePair <> ".key", keyFileString, "Text"];
				keyName = FileNameJoin[{inputPath, particlePair}] <> ".key";
				(* Create particle pair log file *)
				SetDirectory[resultDirectory];
				If[
					!DirectoryQ[particlePair],
					CreateDirectory[particlePair],
					ReadList["!attrib " <> FileNameJoin[{inputDirectory, particlePair}] <> " -r"];
				];
				outputPath = FileNameJoin[{resultDirectory, particlePair}];
				particleLogFile = FileNameJoin[{outputPath, particlePair <> "_log.txt"}];
				If[
					!FileExistsQ[particleLogFile],
					CreateFile[particleLogFile];
					particleLogFile = OpenAppend[particleLogFile],
					particleLogFile = OpenAppend[particleLogFile];
					WriteString[particleLogFile, "\n"];
				];
				(* Splits xyzName if length > 1 and particles are not equal and sets nRun (number fo runs for determination of coordination number) to 2 *)
				If[
					Length[xyzNames]== 1,
					isSameParticle = True;
					nRun = 1,
					xyzNames[[1]] = StringSplit[particlePair, "_"][[1]] <> ".xyz";
					xyzNames[[2]] = StringSplit[particlePair, "_"][[2]] <> ".xyz";
					isSameParticle = False;
					nRun = 2;
				];			
				(* Conformational search via Tinker Scan *)
				(* Create optimized xyz files for the particles via Tinker optimize, if not already present.*)
				Do[
					conformationSearchTime = Now;
					If[
						!DirectoryQ[optXYZDirectory], 
						CreateDirectory[optXYZDirectory]
					];
					SetDirectory[optXYZDirectory];
					(*optXYZdirectoryNames = Select[FileNames["*", "", Infinity], DirectoryQ];*)
					particleName = FileBaseName[xyzNames[[i]]];
					
					If[
						FileExistsQ[FileNameJoin[{optXYZDirectory, particleName, particleName <> ".xyz"}]],
						SetDirectory[particleName];
						CopyFile[
						    xyzNames[[i]], 
							FileNameJoin[{inputDirectory, particlePair, xyzNames[[i]] }],
							OverwriteTarget -> True
						], (* else *)
						CreateDirectory[particleName];
						SetDirectory[particleName];
						CopyFile[
							FileNameJoin[{sourceDirectory, xyzNames[[i]] }],
							xyzNames[[i]]
						];
						elementList = Import[particleName <> ".xyz", "Table"][[2;;All, 2]];
						Export[particleName <> ".key", keyFileString, "Text"];
						
						(* Optimization of particle via Tinker Optimize *)
						ReadList[
							"!" <>
							tinkerOptimize <>
							xyzNames[[i]] <>
							" " <>
							ToString[optimizeRmsGradient] <>
							" > " <> FileNameJoin[{outputPath, particleName <> "_optimize.log"}]
						];
						RenameFile[
							FileBaseName[xyzNames[[i]] ] <> ".xyz_2", xyzNames[[i]],
							OverwriteTarget -> True
						];
						
						(* Run Tinker Scan with optimized particles *)
						ReadList[
							"!" <>
							tinkerScan <>
							xyzNames[[i]] <>
							" " <>
							ToString[scanProgram] <> 
							" " <>
							ToString[nSearchDirection] <> 
							" " <>
							ToString[energyThreshold] <> 
							" " <>
							ToString[rmsGradient] <> 
							" > " <>
							particleName <> 
							".out"
						];						
						minimumList = FindList[particleName <> ".out", "Potential Surface Map       Minimum"];
						If[
							Length[minimumList] > 0,
							minimumList1 = ToExpression[StringTrim[StringTake[minimumList, {42, 46}]]];
							minimumList2 = ToExpression[StringTrim[StringTake[minimumList, {58, 68}]]];
							minimumList = SortBy[Transpose[{minimumList1, minimumList2}], Last];
							xyzTemp = Import[particleName <> "." <> "arc", "Table"]; (*arc -> xyzExtension*)
							tmpStartIndex = 2;
							tmpEndIndex = Length[elementList] + 1;
							tmpHead = xyzTemp[[1]];
							
							Do[
								xyzTemp[[tmpStartIndex;;tmpEndIndex, 2]] = elementList;
								xyzOutput = ConvertTinkerXYZToXYZInMemory[xyzTemp[[tmpStartIndex - 1;;tmpEndIndex]]];
								Export[
									FileBaseName[xyzNames[[i]]] <> "_o" <> ToString[j] <> ".xyz", 
									xyzOutput, 
									"Table", 
									OverwriteTarget -> True
								];
								tmpStartIndex = tmpEndIndex + 2;
								tmpEndIndex += Length[elementList] + 1,						
								{j, Length[minimumList]}
							];
							
							Export[xyzNames[[i]],
							  xyzTemp[[1;;Length[elementList] + 1]],
							  "Table",
							  OverwriteTarget -> True
							];
							CopyFile[
								FileBaseName[xyzNames[[i]]] <> ".xyz",
								FileNameJoin[{inputDirectory, particlePair, xyzNames[[i]]}],
								OverwriteTarget -> True
							];
							
							(* Creates notebook and adds particle configuration of lowest intermolecular energy *)
							notebook = 
								CreateDocument[
									{
										Table[
											{
												MoleculePlot3D[Import[
													particleName <> 
													"_o" <> ToString[j] <> ".xyz",
													"XYZ"
												]],
												minimumList[[j, 1]],
												minimumList[[j, 2]]
											},
											{j, Length[minimumList]}
										]
									},
									Visible->False
								];
							NotebookEvaluate[notebook, InsertResults -> True];
							NotebookSave[notebook, FileNameJoin[{optXYZDirectory, particleName, particleName}]];
							NotebookClose[notebook],
							
							(* Length[minimumList] == 0 *)
							ConvertTinkerXYZToXYZ[
								particleName <> ".xyz",
								particleName <> "_o1." <> "xyz"
							];
							newXYZData = Import[particleName <> "_o1." <> "xyz", "Table"];
							newXYZData[[2;;All, 1]] = elementList;
							(*DeleteFile[particleName <> "_o1." <> "xyz"];*)
							Export[
								particleName <> "_o1." <> "xyz", 
								newXYZData, 
								"Table",
								OverwriteTarget -> True
							];
							newXYZData = Import[particleName <> ".xyz", "Table"];
							newXYZData[[2;;All, 2]] = elementList;
							Export[
								particleName <> ".xyz", 
								newXYZData, 
								"Table",
								OverwriteTarget -> True
							];
							CopyFile[
								xyzNames[[i]],
								FileNameJoin[{inputDirectory, particlePair, xyzNames[[i]]}],
								OverwriteTarget -> True
							];
							SetDirectory[FileNameJoin[{optXYZDirectory, particleName}]];
							notebook = 
								CreateDocument[
									{
										MoleculePlot3D[Import[particleName <> "_o1." <> "xyz", "XYZ"]]
									},
									Visible->False
								];
							NotebookEvaluate[notebook, InsertResults -> True];
							NotebookSave[notebook, FileNameJoin[{optXYZDirectory, particleName, particleName}]];
							NotebookClose[notebook];
						];
						WriteString[particleLogFile, StringForm["Time for Conformational Search via Tinker Scan: ``\n", ToString[Now - conformationSearchTime]]];
						conformationSearchTime =.;
					],
					{i, Length[xyzNames]}
				];

				(* Import tinker - xyz-data *)
				SetDirectory[inputPath];
				xyzData1 = Import[xyzNames[[1]], "Table"];
				If[
					isSameParticle == True,
					xyzData2 = xyzData1,
					xyzData2 = Import[xyzNames[[2]], "Table"];
				];
				(* Import coordinate data and adjust to the center (0,0,0) *)
				importData1 = Drop[xyzData1, 1];
				numberAtoms1 = Dimensions[importData1][[1]];
				coordinateData1 = importData1[[All, 3;;5]];
				centerCoordinate1 = 
					Sum[
						coordinateData1[[i]], 
						{i, numberAtoms1}
					] / numberAtoms1;
				adjustedCoordinate1 = 
					Table[
						coordinateData1[[i]] - centerCoordinate1,
						{i, numberAtoms1}
					];
				If[
					isSameParticle == True,
					(* true *)
					importData2 = importData1;
					numberAtoms2 = numberAtoms1;
					adjustedCoordinate2 = adjustedCoordinate1,
					(* false *)
					importData2 = Drop[xyzData2, 1];
					numberAtoms2 = Dimensions[importData2][[1]];
					coordinateData2 = importData2[[All, 3;;5]];
					centerCoordinate2 = 
						Sum[
							coordinateData2[[i]],
							{i, numberAtoms2}
						] / numberAtoms2;
					adjustedCoordinate2 = 
						Table[
							coordinateData2[[i]] - centerCoordinate2,
							{i, numberAtoms2}
						];
				];				
				(* Calculates the rotated atom coordinates using the rotation matrices *)
				SetDirectory[scratchDirectory];
				rotatedCoordinate1 = Table[0., configurationNumber1];
				rotatedCoordinate2 = Table[0., configurationNumber2];
				Do[
					rotatedCoordinate1[[i]] =
						Table[
							rotationMatrices1[[i]] . adjustedCoordinate1[[j]],
							{j, numberAtoms1}
						],
					{i, configurationNumber1}
				];
				Do[
					rotatedCoordinate2[[i]] =
						Table[
							rotationMatrices2[[i]] . adjustedCoordinate2[[j]],
							{j, numberAtoms2}
						],
					{i, configurationNumber2}
				];
				(* Scans for the lowest intermolecular energy between particles *)				
                DeleteFile[keyName];
                keyFileString = 
					keyFileStringOrigin <>
					"DIELECTRIC\t" <> ToString[dielectricConstant] <>
					"\nNONBONDTERM ONLY";			
				Export[particlePair <> ".key", keyFileString, "Text"];
				keyName = particlePair <> ".key";
				(* Prescan *)			
				energyCalculationTime = Now;				
				distances = Range[lowerBoundary, upperBoundary, prescanStepSize];				
				minEnergyGlobal = 100.;
				minDistanceGlobal = 100.;				
				interMolecularEnergyResult = 
					GetInterMolecularEnergy[
						tinkerAnalyze,
						distances,
						cpuCoreNumber,
						minAtomDistance,
						keyName,
						xyzData1,
						xyzData2,
						numberAtoms1,
						numberAtoms2,
						rotatedCoordinate1,
						rotatedCoordinate2
					];	
				If[
					minEnergyGlobal >= interMolecularEnergyResult[[1]],
					minEnergyGlobal = interMolecularEnergyResult[[1]];
					minDistanceGlobal = interMolecularEnergyResult[[2]]
				];	
				energyDataTable = interMolecularEnergyResult[[3]];
				allDistances = distances;
				(*Precise scan*)
				distances = 
					Range[
						minDistanceGlobal - prescanStepSize + stepSize,
						minDistanceGlobal + prescanStepSize - stepSize,
						stepSize
					];
				interMolecularEnergyResult = 
					GetInterMolecularEnergy[
						tinkerAnalyze,
						distances,
						cpuCoreNumber,
						minAtomDistance,
						keyName,
						xyzData1,
						xyzData2,
						numberAtoms1,
						numberAtoms2,
						rotatedCoordinate1,
						rotatedCoordinate2
					];	
				If[
					minEnergyGlobal >= interMolecularEnergyResult[[1]],
					minEnergyGlobal = interMolecularEnergyResult[[1]];
					minDistanceGlobal = interMolecularEnergyResult[[2]]
				];	
				energyDataTable = Join[energyDataTable, interMolecularEnergyResult[[3]] ];
				allDistances = Join[allDistances, distances];
				(* More precise scan*)
				distances = 
					Range[
						minDistanceGlobal - stepSize + preciseStepSize,
						minDistanceGlobal + stepSize - preciseStepSize,
						preciseStepSize
					];
				interMolecularEnergyResult = 
					GetInterMolecularEnergy[
						tinkerAnalyze,
						distances,
						cpuCoreNumber,
						minAtomDistance,
						keyName,
						xyzData1,
						xyzData2,
						numberAtoms1,
						numberAtoms2,
						rotatedCoordinate1,
						rotatedCoordinate2
					];	
				If[
					minEnergyGlobal >= interMolecularEnergyResult[[1]],
					minEnergyGlobal = interMolecularEnergyResult[[1]];
					minDistanceGlobal = interMolecularEnergyResult[[2]]
				];	
				energyDataTable = Join[energyDataTable, interMolecularEnergyResult[[3]] ];	
				allDistances = Join[allDistances, distances];
				allDistances = DeleteDuplicates[allDistances];
				(*weightedAverageList = interMolecularEnergyResult[[1]];
				morePreciseXyData = 
					Table[
						{distances[[i]], weightedAverageList[[i]]},
						{i, Length[distances]}
					];
				If[
					minWeightedEnergy >= interMolecularEnergyResult[[2]],
					minWeightedEnergy = interMolecularEnergyResult[[2]];
					minDistance = interMolecularEnergyResult[[3]]
				];
				energyDataTable = Join[energyDataTable, interMolecularEnergyResult[[4]] ];*)
				(* energyDataTable is a list of lists:  *)
				(* {
					{distances1, {n*n*m energy values for distances1}},
					{distances2, {n*n*m energy values for distances2}},
					{distances3, {n*n*m energy values for distances3}},
					....
				} *)
				(* So energyDataTable[[i, 1]] returns a distance and energyDataTable[[i, 2]] returns a list of energy vales *)
				energyDataTable = DeleteDuplicates[energyDataTable];
				energyDataTable = SortBy[energyDataTable, 1];
				xyData =
					Table[
						{
							energyDataTable[[i,1]],
							Min[energyDataTable[[i,2]]]
						},
						{i, Length[energyDataTable]}
					];
				(* Export global minimum energy configuration *)
				energyMinimumIndex = interMolecularEnergyResult[[4]];
				stream = OpenRead["output" <> ToString[energyMinimumIndex] <> ".out"];
				energyMinimumValueLines =
	  				Flatten[
	   					StringCases[
	    					FindList[stream, "Intermolecular Energy :"],
	    					x : NumberString :> ToExpression[x]
	    				]
	   				];
				Close[stream];
				lineIndex = FirstPosition[energyMinimumValueLines, minEnergyGlobal][[1]];						
				nEntryLines = numberAtoms1 + numberAtoms2 + 1;
				stream = OpenRead["output." <> ToString[energyMinimumIndex]];
				Skip[stream, Record, nEntryLines*(lineIndex - 1)];
				energyMinimumXyzData = ReadList[stream, Record, numberAtoms1 + numberAtoms2 + 1];
				Close[stream];
	  			Export["output.0", energyMinimumXyzData, "Text", OverwriteTarget -> True]; 
	  			Export["output0.out", minEnergyGlobal, "Table", OverwriteTarget -> True];
	  			(* Minimize lowest energy conformation *)
				keyFileString = 
					keyFileStringOrigin <> 
					"DIELECTRIC\t" <> ToString[dielectricConstant];				
				Export[particlePair <> ".key", keyFileString, "Text"];
				keyName = FileNameJoin[{scratchDirectory, particlePair}] <> ".key";			
				ReadList[
					"!" <>
					tinkerOptimize <> 
					"output.0" <> 
					" -k " <> 
					keyName <> 
					" " <> 
					ToString[optimizeRmsGradient] <>
					" > optimizeEnergy.log",
					Record
				];			
				ReadList[
					"!" <>
					tinkerAnalyze <> 
					"output.xyz" <> 
					" -k " <> 
					keyName <> 
					" e" <> 
					" >output0_opt.txt"
				];
				If[
					FileExistsQ["output0_opt.txt"],
					optimizedData = OpenRead["output0_opt.txt"];
					optimizedMinimum = ToExpression[StringTrim[StringTake[Find[optimizedData, "Intermolecular Energy :"], {25, 50}]]];
					Close[optimizedData];
				];
				outputData[[f, 9]] = "Global Minimum intermolecular Energy = " <> ToString[minEnergyGlobal] <> "\nOptimized minimum intermolecular Energy [kcal/mole] = " <> ToString[optimizedMinimum];
				outputData[[f, 10]] = "Energy difference [kcal/mole] = " <> ToString[minEnergyGlobal - optimizedMinimum];		
				minEnergyGlobal = optimizedMinimum;
				(* Boltzmann averaging *)
				If[
					boltzmannFraction == 0.0,
					(* TRUE *)
					(* No averaging, min energy value of each configuration is taken *)					
					minEnergy = minEnergyGlobal;
					minDistance = minDistanceGlobal,
					(* FALSE *)	
					(* If fractionForAverage = 1.0 all configurational E(nonbonded) values are used for "Boltzmann average" calculation *)
					(* 0.0 < fractionForAverage < 1.0: All configurational E(nonbonded) values are sorted ascending and *)
					(* the lower "numberOfValues*fractionForAverage" E(nonbonded) values are used for "Boltzmann average" calculation *)
		     		(* Example: For 144x144x16 = 331776 E(nonbonded) values for a specific molecule distance r and *)
		     		(* a fractionForAverage of 0.25 the lowest Round(331776x0.25) = 82944 E(nonbonded) values are used for *)
		     		(* "Boltzmann average" calculation only *)
				   	xyDataWeighted = Table[{}, Length[energyDataTable]];
				   	nXyDataWeighted = Length[xyDataWeighted];
					(* Dulicated comment for better readability *) 
					(* energyDataTable is a list of lists:  *)
					(* {
						{distances1, {n*n*m energy values for distances1}},
						{distances2, {n*n*m energy values for distances2}},
						{distances3, {n*n*m energy values for distances3}},
						....
					} *)
					(* So energyDataTable[[i, 1]] returns a distance and energyDataTable[[i, 2]] returns a list of energy vales *)
		     		Do[
		     			energyDataSorted = Sort[energyDataTable[[i,2]] ];		     			
		     			fractionToMax = Ceiling[Length[energyDataSorted] * boltzmannFraction];
		     			energyDataFractioned = energyDataSorted[[;;fractionToMax]];
		     			weight = 
		     				Table[
		     					Exp[-(energyDataFractioned[[j]] - minEnergyGlobal)/(temperature * gasConstant)],(*energyData in kcal/mol;kB in kcal/(mol*K)*)
		     					{j, fractionToMax}
		     				];
		     			xyDataWeighted[[i]] = {energyDataTable[[i,1]], Total[weight * energyDataFractioned] / Total[weight]},	
		     			{i, nXyDataWeighted}
		     		];
		     		minEnergy = Min[xyDataWeighted[[All, 2]] ];
		     		minDistance = xyDataWeighted[[FirstPosition[xyDataWeighted, minEnergy][[1]] ]][[1]];
				];	
				(* --- *)			
				WriteString[particleLogFile, StringForm["Eij\n"]];
				WriteString[particleLogFile, StringForm["Raw time for determining the intermolecular minimal energy: ``\n", ToString[Now - energyCalculationTime]]];
				outputData[[f, 1]] = "Fraction for Boltzmann average = " <> ToString[boltzmannFraction] <> "\nMinimum Distance [\[CapitalARing]] = " <> ToString[minDistance];
				outputData[[f, 2]] = weightedMinEnergyTextString <> ToString[minEnergy];
				outputData[[f, 3]] = "temperature [K] = " <> ToString[temperature];
				WriteString[particleLogFile, StringForm["(Weighted) Minimum IntermolecularEnergy [kcal/mole] = ``\n", ToString[minEnergy]]];				
				If[
					!DirectoryQ[FileNameJoin[{outputPath, "Charts", particlePair, "Eij_histograms"}]],
					CreateDirectory[FileNameJoin[{outputPath, "Charts",particlePair, "Eij_histograms"}]]
				];
				energyDataHistogramTable = 
					Table[
						Histogram[
							energyDataTable[[i,2]], 
							20, (* number of bins *)
							ChartLegends -> {StringJoin[{ToString[energyDataTable[[i,1]] ], "\[CapitalARing]"}]},
		          			PlotRange -> All,
		          			PlotRangePadding -> Scaled[.05],
		          			Frame -> True,
		          			FrameLabel -> {"Eij", None},
							PlotTheme -> "Scientific",
							ChartStyle -> RGBColor["#3277a8"],
							ImageSize -> 1000
						],
						{i, Length[energyDataTable]}						
					];
				Export[
					FileNameJoin[{outputPath, "Charts", particlePair, "Eij_histograms", "Eij_histogram.nb"}],
					energyDataHistogramTable
				];					
				SetOptions[ListLinePlot, ImageSize -> {425, Automatic}];
				If[
					boltzmannFraction == 0.0,
					(* TRUE *)
					energyListLinePlots[[f]] = 
						ListLinePlot[
							Transpose@{xyData[[All, 1]], xyData[[All, 2]]},
							Frame -> True,
							FrameLabel -> {"Distance [\[CapitalARing]]", "Intermolecular Energy [kcal/mole]"},
							PlotLabel -> particlePair <> "_" <> ToString[temperature] <> "K_boltzmannFraction" <> ToString[boltzmannFraction], 
							PlotMarkers -> Automatic,
							PlotRangePadding -> Scaled[.05]
						],
					(* FALSE *)
					energyListLinePlots[[f]] = 
						ListLinePlot[
							{
								Transpose@{xyData[[All, 1]], xyData[[All, 2]]},
								Transpose@{xyDataWeighted[[All, 1]], xyDataWeighted[[All, 2]]}
							},
							Frame -> True,
							FrameLabel -> {"Distance [\[CapitalARing]]", "Intermolecular Energy [kcal/mole]"},
							PlotLabel -> particlePair <> "_" <> ToString[temperature] <> "K_boltzmannFraction" <> ToString[boltzmannFraction], 
							PlotMarkers -> Automatic,
							PlotRangePadding -> Scaled[.05]
						];	
				];
				If[
					!DirectoryQ[FileNameJoin[{outputPath, "Charts", particlePair, "Eij_ListLines"}]],
					CreateDirectory[FileNameJoin[{outputPath, "Charts", particlePair, "Eij_ListLines"}]]	
				];	
				Export[
					FileNameJoin[{outputPath, "Charts", particlePair, "Eij_ListLines", "Eij_ListLines.nb"}],
					energyListLinePlots[[f]]
				];
				Export[
					FileNameJoin[{outputPath, "Charts", particlePair, "xyData.txt"}],
					xyData
				];
				(*energyPreciseListLinePlots[[f]] = 
					ListLinePlot[
						Transpose@{preciseXyData[[All, 1]], preciseXyData[[All, 2]]},
						Frame -> True,
						FrameLabel -> {"Distance [\[CapitalARing]]", "Intermolecular Energy [kcal/mole]"},
						PlotLabel -> particlePair <> "_" <> ToString[temperature] <> "K",
						PlotMarkers -> Automatic,
						PlotRangePadding -> Scaled[.05]
					];	
				Export[
					FileNameJoin[{outputPath, "Charts", particlePair, "Eij_ListLines", "Eij_precise_ListLines.nb"}],
					energyPreciseListLinePlots[[f]]
				];	*)	

				(*Close[FileNameJoin[{scratchDirectory, "output0.out"}]];*)
				(* Copy results *)
				SetDirectory[resultDirectory];
				If[
					!DirectoryQ[particlePair],
					CreateDirectory[particlePair];
					ReadList["!attrib " <> FileNameJoin[{inputDirectory, particlePair}] <> " -r"];
				];
				SetDirectory[particlePair];
				CopyFile[
					FileNameJoin[{scratchDirectory, "output.0"}],
					FileNameJoin[{outputPath, "output.0"}],
					OverwriteTarget -> True
				];
				CopyFile[
					FileNameJoin[{scratchDirectory, "output0.out"}],
					FileNameJoin[{outputPath, "output0.out"}],
					OverwriteTarget -> True
				];
				CopyFile[
					FileNameJoin[{scratchDirectory, "output.xyz"}],
					FileNameJoin[{outputPath, "output_optimized.0"}], 
					OverwriteTarget -> True
				];
				Export[
					FileNameJoin[{resultDirectory, particlePair,"output0_optimized.out"}],
					StringJoin[{"Intermolecular Energy : ", ToString[optimizedMinimum]}], 
					"Table", 
					OverwriteTarget -> True
				];
				DeleteFile[FileNameJoin[{scratchDirectory, "output.xyz"}]];
				ConvertTinkerXYZToXYZ["output.0", "output.xyz"];
				ConvertTinkerXYZToXYZ["output_optimized.0","output_optimized.xyz"];
				ReadList[
					"!" <>
					tinkerXYZpdb <> 
					"output_optimized.0" <> 
					" -k " <> 
					keyName,
					Record
				];
				If[
					!DirectoryQ[FileNameJoin[{calculationDirectory, "OptDist"}] ],
					CreateDirectory[FileNameJoin[{calculationDirectory, "OptDist"}] ]
				];
				SetDirectory[FileNameJoin[{calculationDirectory, "OptDist"}]];
				If[
					isSameParticle == True && !FileExistsQ[particleName <> ".txt"],
					Export[particleName <> ".txt", minDistance, "Text"];
				];
				SetDirectory[scratchDirectory];
				DeleteFile[FileNames["*", "", Infinity]];
				energyCalculationTime = Now - energyCalculationTime;
				WriteString[particleLogFile, StringForm["Time for determining the intermolecular minimal energy: ``\n", ToString[energyCalculationTime]]];
				Close[particleLogFile],
				{f,Length[inputdirectoryName]}
			];
			WriteString[gblLogFile, StringForm["Energy calculation time: ``\n\n", ToString[Now - gblEnergyCalculationTime]] ];
			(*TODO: end of energy calculation*)
			gblCoordinationTime = Now;
			(* Zij calculation *)	
			allParticles = Union[oldParticles, newParticles];
            allParticlePairs = Tuples[allParticles, 2];
            oldParticlePairs = Tuples[oldParticles, 2];
            allParticlePairs = Complement[allParticlePairs, oldParticlePairs];            
			totalHistogram = Table[{}, Length[inputdirectoryName]];
			averageCharts = Table[{}, Length[inputdirectoryName]];
			neighborParticleNumbers = Table[{}, Length[allParticlePairs], Length[catchRadius]];
			neighborParticleMean = Table[{}, Length[allParticlePairs], Length[catchRadius]];
			SetSharedVariable[neighborParticleNumbers, neighborParticleMean, outputData, totalHistogram, averageCharts];
			DistributeDefinitions[allParticlePairs, allParticles, catchRadius,
				cdkJarFullPath, cellI, cellLength, cellListCellsPerRow, cellListHead, cellListIndexCube, cellListList,
				cellListMaxDistances, cellListNeighborCells, cellListNeighborCellsItems, 
				cellListNeighborCellsSolventItems, cellListResult, cellListSolubleCellsList,
				elementSymbolToAtomicNumberTable, endIndex, exportCharts, inputDirectory, inputdirectoryName, keyFileStringOrigin,
				minimizeMaxIteration, moleculeLookupTable, nDynamicIteration, nDynamicSteps, nDynamicWarmupIteration, nStepsPerRound, printInterval, randomseed, resultDirectory,
				rmsMinimizeGradient, scratchDirectory, simulationType, smiles, solventMoleculeNumber, temperature, timeStep, tinkerDynamic, 
				tinkerMinimize, tinkerXYZEdit, useCellList, vdWRadii, warmUpNeighborIndices, warmUpPrintInterval, warmUpStepNumber, warmUpTimeStep,
				waterVolumeRatio,
				ConvertXYZDataToXYZ, CreateCellList, CreateIndexCube, GetCellItems, GetMoleculeIndexByAtomIndex, GetNeighborCells, 
				GetNeighborParticleNumbers, GetSolubleCoordinates, GetSolventCoordinates, GetVdwVolume, gblLogFile
			];
    		LaunchKernels[cpuCoreNumber];
	    	x = withModifiedMemberQ@ParallelTable[
	    		particlePairTime = Now;
               	particlePair = allParticlePairs[[j]];
               	particlePairString = particlePair[[1]] <> "_" <> particlePair[[2]];
               	inputdirectoryNamePosition = FirstPosition[inputdirectoryName, particlePair[[1]] <> "_" <> particlePair[[2]], -1];
               	If[
	               	particlePair[[1]] == particlePair[[2]],
	                isSameParticle = True;
	                nRun = 1,
	                isSameParticle = False;
	                nRun = 2
	            ];	          
	            nRunIteration = 1;         	                
                If[
	                inputdirectoryNamePosition == -1,
	                inputdirectoryNamePosition = FirstPosition[inputdirectoryName, particlePair[[2]] <> "_" <> particlePair[[1]]];
	                nRunIteration = 2
	            ]; 
	            inputdirectoryNamePosition = inputdirectoryNamePosition[[1]];                        
                inputdirectoryNameString = ToString[inputdirectoryName[[inputdirectoryNamePosition]] ]; 
               	outputPath = FileNameJoin[{resultDirectory, inputdirectoryNameString}];
	            particleLogFile = OpenAppend[FileNameJoin[{outputPath, inputdirectoryNameString <> "_log.txt"}]];
	            totalHistogram[[inputdirectoryNamePosition]] = Table[{}, nRun];
	            averageCharts[[inputdirectoryNamePosition]] = Table[{}, nRun];
	            xyzNames = # <> ".xyz"& /@ particlePair;
	            xyzNames = DeleteDuplicates[xyzNames];
	            scratchDirectoryParticlePair = FileNameJoin[{scratchDirectory, particlePairString}];
	            If[
	             	!DirectoryQ[scratchDirectoryParticlePair],
	                CreateDirectory[scratchDirectoryParticlePair]
	            ];
                coordinationCalculationTime = Now;
                inputPath = FileNameJoin[{inputDirectory, inputdirectoryNameString }];
                SetDirectory[inputPath];
                CopyFile[FileNameJoin[{inputPath, #}], FileNameJoin[{scratchDirectoryParticlePair, #}]]& /@ xyzNames;
                SetDirectory[scratchDirectoryParticlePair];
                solubleName = particlePair[[1]];
                solventName = particlePair[[2]];
                solubleData = Import[solubleName <> ".xyz", "Table"];                    
                nSolubleAtoms = ToExpression[solubleData[[1, 1]]];
                If[
                	ToExpression[solubleData[[2, 1]]] == 1,
                    solubleData = Drop[solubleData, 1],
                    solubleData = Drop[solubleData, 2]
                ];
                Do[
                    solubleData[[k, 2]] = elementSymbolToAtomicNumberTable[
                           solubleData[[k, 2]]], 
                    {k, Length[solubleData]}
                ];
                solubleAtomNames = solubleData[[All, 2]];                
                solventData = Import[solventName <> ".xyz", "Table"];
                nSolventAtoms = ToExpression[solventData[[1, 1]]];
                If[
                	ToExpression[solventData[[2, 1]]] == 1,
                    solventData = Drop[solventData, 1],
                    solventData = Drop[solventData, 2]
                ];
                Do[
					solventData[[k, 2]] = 
						elementSymbolToAtomicNumberTable[
                       		solventData[[k, 2]]	
						],
                    {k, Length[solventData]}
                ];
                solventAtomNames = solventData[[All, 2]];
                (* Determine box length *)
                vdwSolventVolume = GetVdwVolume[smiles[solventName], cdkJarFullPath];
                (* Calculated water volume ratio - ratio of Vparticle and Vvdw of water = 1.7297 *)
                boxLength = (waterVolumeRatio * solventMoleculeNumber * vdwSolventVolume)^(1.0/3.0);
                (* Key file for Tinker dynamic *)
                SetDirectory[scratchDirectoryParticlePair];
                keyFileString = 
                    keyFileStringOrigin <> "a-axis " <> ToString[boxLength] <>
                    "\nRANDOMSEED " <> ToString[randomseed] <>
                    "\nEWALD" <> (* using Ewald summation for the partial charges in the periodic box *)
                    "\nTHERMOSTAT ANDERSEN" <> (* Andersen stochastic collision method *)
                    "\nSTEEPEST-DESCENT\n"; (* forces the L-BFGS optimization routine used by the MINIMIZE *)
					(*	program and other programs to perform steepest descent minimization. This option can be useful in *)
					(* conjunction with small step sizes for following minimum energy paths, but is generally inferior to the *)
					(* L-BFGS default for most optimization purposes. *)
				(* NEIGHBOR-LIST is buggy for boxes with less then 400 solvent particles *)
                (*If[
					solventMoleculeNumber >= 400,
                    keyFileString = keyFileString <> "NEIGHBOR-LIST\n" 
					(* pairwise neighbor lists for partial charge electrostatics, *) 
					(* polarize multipole electrostatics and any of the van der Waals potentials. *) 
					(* This method will yield identical energetic results to the standard double loop method. *)
                ];*)
                Export[particlePairString <>".key", keyFileString, "Text"];
                (* Build solvent box *)
                ReadList[
                	"!" <> 
                	tinkerXYZEdit <>
	                solventName <>
	                ".xyz -k " <>
	                ToString[particlePairString] <>
	                ".key " <>
	                "23 " <> (* Option 23: Create and Fill a Periodic Boundary Box *)
	                ToString[solventMoleculeNumber] <>
	                " " <>
	                ToString[boxLength] <>
	                " " <>
	                ToString[boxLength] <>
	                " " <>
	                ToString[boxLength] <>
	                " Y" <>
	                " > " <> FileNameJoin[{outputPath, ToString[particlePairString] <> "build.log"}],
	                Record
	            ];
                (* Add soluble to solvent box *)
                ReadList[
                	"!" <>
                	tinkerXYZEdit <>
                   	solubleName <>
                    ".xyz -k " <>
                    ToString[particlePairString] <>
                    ".key " <>
                    "24 " <> (* Option 24: Soak Current Molecule in Box of Solvent *)
                    solventName <>
                    ".xyz_2 " <>
                    " > " <> FileNameJoin[{outputPath, ToString[particlePairString] <> "soak.log"}],
                    Record
                ];
                dataName = particlePairString <> ".xyz";
                If[
                   solubleName == solventName,
                   RenameFile[solubleName <> ".xyz_3", dataName],
                   RenameFile[solubleName <> ".xyz_2", dataName]
                ];                
                (* Minimize solvent box *)
                boxMinimizationTime = Now;
                WriteString[
                    OpenAppend[particlePairString <> ".key"], 
                    "MAXITER " <> ToString[minimizeMaxIteration] (* Sets the maximum number of minimization iterations *)
                ];
                Close[particlePairString <> ".key"];
                ReadList[
                    "!" <> 
                    tinkerMinimize <>
                    dataName <>
                    " " <>
                    ToString[rmsMinimizeGradient] <>
                    " > " <> FileNameJoin[{outputPath, FileBaseName[dataName] <> "minimize.log"}],
                    Record
                ];
                RenameFile[
                    FileBaseName[dataName] <> ".xyz_2",
                    dataName,
                    OverwriteTarget->True
                ];
                WriteString[particleLogFile, StringForm["Coordination Number\n"]];
                WriteString[
                     particleLogFile, 
                     StringForm["Time for solvent box minimization via Tinker \"Minimize\": ``\n", ToString[Now - boxMinimizationTime]]
                ];
                boxMinimizationTime =.;
              	(* Simulation warm up *)
                dynamicWarmUpTime = Now;
                ReadList[
                	"!" <>
                    tinkerDynamic <>
                    dataName <>
                    " " <>
                    ToString[warmUpStepNumber] <>
                    " " <>
                    ToString[warmUpTimeStep] <>
                    " " <>
                    ToString[warmUpPrintInterval] <>
                    " " <>
                    ToString[simulationType] <>
                    " " <>
                    ToString[temperature] <>
                    " > " <> FileNameJoin[{outputPath, FileBaseName[dataName] <> "warmUp.log"}],
                    Record
                ];
                WriteString[
                	particleLogFile, 
                    StringForm["Time for box warm up via Tinker \"Dynamic\": ``\n", ToString[Now - dynamicWarmUpTime]]
                ];
                dynamicWarmUpTime =.;
                stream = OpenRead[particlePairString <> ".arc"];
                nAllAtoms = Read[stream, Number];
                Close[stream];
                (* Import warm up arc file *)      
				solventAtomsStartIndex = nSolubleAtoms + 1;
				nAllSolventAtoms = nAllAtoms - nSolubleAtoms;
				nSolventMolecules = nAllSolventAtoms/nSolventAtoms;
				allSolventAtomsNames = Flatten[Table[solventAtomNames, nSolventMolecules]];
				allSolventAtomsIndices = Range[nSolubleAtoms + 1, nAllAtoms];
				allAtomsNames = Join[solubleAtomNames, allSolventAtomsNames];
				moleculeLookupTable = 
					Table[
						If[
				   			j <= nSolubleAtoms,
				   			{1},
				   			{IntegerPart[(j - nSolubleAtoms - 1)/nSolventAtoms] + 2}
				   		],
				  		{j, nAllAtoms}
				  	];
				
                (* preprocessing cell list method *) 
                If[
                	useCellList,
	                (* true *)
	                cellListMaxDistances = Max[{vdWRadii[[#]] & /@ solubleAtomNames, vdWRadii[[#]] & /@ solventAtomNames}] + catchRadius[[1]];
					cellLength = boxLength/cellListMaxDistances;
					cellListCellsPerRow = IntegerPart[cellLength];
					cellLength = boxLength/cellListCellsPerRow;
					cellI = cellLength^-1;
					cellListIndexCube = CreateIndexCube[cellListCellsPerRow];
					boxLengthHalf = boxLength * 0.5;
                ];
                If[
                	useCellList,
                   	(* True - cell list *)     
                   	WriteString[particleLogFile, "Used cell list method\n"];                   		
                   	stream = OpenRead[particlePairString <> ".arc"];   
                    Skip[stream, Record, 2];
                    xyzDataWarmUpCoordsInput = ReadList[stream, {Number, Word, Real, Real, Real, Record}, nAllAtoms][[All, 3 ;; 5]];                   
                    Close[stream];                    
                    xyzDataWarmUpCoords = xyzDataWarmUpCoordsInput + (boxLengthHalf);
                    cellListResult = CreateCellList[cellListCellsPerRow, cellI, xyzDataWarmUpCoords, nAllAtoms, nSolubleAtoms]; 
                    cellListHead = cellListResult[[1]];
                    cellListList = cellListResult[[2]];
                    cellListSolubleCellsList = DeleteDuplicates[cellListResult[[3]] ]; 
                    (* neighbor cells for each cell that contains an atom of the soluted molecule; including the cell itself *)
                    cellListNeighborCells = DeleteDuplicates[Flatten[GetNeighborCells[#, cellListCellsPerRow, cellListIndexCube]&/@cellListSolubleCellsList]];
                    (* atoms/items contained in the neighborCells *)
					cellListNeighborCellsItems = DeleteCases[DeleteDuplicates[Flatten[GetCellItems[cellListHead, cellListList, #]&/@cellListNeighborCells]],0];
					(* atoms/items of solvent molecules contained in the neighborCells *)
					cellListNeighborCellsSolventItems = Drop[Sort[cellListNeighborCellsItems], nSolubleAtoms];
					(* indices of neighbor molecules *)
					warmUpNeighborIndices = 
						GetNeighborParticleNumbers[
							xyzDataWarmUpCoordsInput,
							allAtomsNames,
							nSolubleAtoms,
							cellListNeighborCellsSolventItems,
							catchRadius[[1]],
							boxLength,
							vdWRadii
						], 
                   	(* False - brute force *)
                   	WriteString[particleLogFile, "Used brute force method\n"]; 
                   	stream = OpenRead[particlePairString <> ".arc"];   
                    Skip[stream, Record, 2];
                    xyzDataWarmUpCoords = ReadList[stream, {Number, Word, Real, Real, Real, Record}, nAllAtoms][[All, 3 ;; 5]];      
                   	(* simulation *)             
                    warmUpNeighborIndices =                   
                    	GetNeighborParticleNumbers[
	                   		xyzDataWarmUpCoords,
	                   		allAtomsNames,
					        nSolubleAtoms,
					        allSolventAtomsIndices,
							catchRadius[[1]],
							boxLength,
							vdWRadii
                    	];   
                    Close[stream];                                     
                ]; 
				If[
					warmUpNeighborIndices == {},
					WriteString[
						particleLogFile, 
						StringForm[
							"warmUpNeighborIndices is empty for ``\n",
							ToString[particlePair]
						]
					];	
				];
                (* Export coordinates of warmup *)   
                stream = OpenRead[particlePairString <> ".arc"]; 
                xyzDataWarmUpCoords = ReadList[stream, {Number, Word, Real, Real, Real, Record}, nAllAtoms];       
                Export[FileNameJoin[{outputPath, FileBaseName[dataName] <> "_warmUpCoords.txt"}], xyzDataWarmUpCoords, "Table"];                  
                Close[stream];                	                    	              
                warmUpNeighborIndices = DeleteCases[DeleteDuplicates[Flatten[warmUpNeighborIndices]],0];  
				If[
					warmUpNeighborIndices == {},
					WriteString[
						particleLogFile, 
						StringForm[
							"warmUpNeighborIndices is empty for `` after DeleteDuplicates\n",
							ToString[particlePair]
						]
					];	
				];
                warmUpNeighborIndices = DeleteDuplicates[Flatten[moleculeLookupTable[[#]]&/@warmUpNeighborIndices  ]];
                Export[FileNameJoin[{outputPath, FileBaseName[dataName] <> "_warmUpNeighbors.txt"}], warmUpNeighborIndices, "List"];              
                dataName = StringJoin[{particlePairString, ".xyz"}];
                RenameFile[StringJoin[{particlePairString, ".arc"}], dataName, OverwriteTarget -> True];                  
                rawDynamicTime = Now;
                If[
                    nDynamicSteps > nStepsPerRound,
                    nIterations = Ceiling[nDynamicSteps / nStepsPerRound],
                    nIterations = 1;
                    nStepsPerRound = nDynamicSteps              	
                ];             
                rawDynamicTime=.;
                (* Analyse arc file *)
                dynamicAnalysisTime = Now;                                      
                If[
                    useCellList,
                    (* true - use cell list method *)         
                    neighborIndices = 
                    	Flatten[
                    		Last[
                    			Reap[
		                    		Do[
		                    			ReadList[
						                    "!" <>
						                    tinkerDynamic <>
						                    dataName <>
						                    " " <>
						                    ToString[nStepsPerRound] <>
						                    " " <>
						                    ToString[timeStep] <>
						                    " " <>
						                    ToString[printInterval] <>
						                    " " <>
						                    ToString[simulationType] <>
						                    " " <>
						                    ToString[temperature] <>
						                    " > " <> FileNameJoin[{outputPath, FileBaseName[dataName] <> "dynamic.log"}],
						                    Record
						                ];
						                stream = OpenRead[particlePairString <> ".arc"];						                    
						                Do[						                    	
			                    			Skip[stream, Record, 2];
			                    			xyzDataListInput = ReadList[stream, {Number, Word, Real, Real, Real, Record}, nAllAtoms][[All, 3 ;; 5]];			                    				
			                    			xyzDataList = xyzDataListInput + boxLengthHalf;
			                    			cellListResult = CreateCellList[cellListCellsPerRow, cellI, xyzDataList, nAllAtoms, nSolubleAtoms]; 
											If[
												cellListResult == {},
												WriteString[
													particleLogFile, 
													StringForm[
														"cellListResult is empty for `` in step round `` and iteration ``\n",
														ToString[particlePair],	ToString[n], ToString[m]
													]
												];	
												WriteString[
													particleLogFile, 
													StringForm[
														"Memory in use [b]: ``\n",
														ToString[N[UnitConvert[Quantity[MemoryInUse[], "Bytes"], "Megabytes"]]]
													]
												];											
											];
											If[
												cellListResult[[3]] == {},
												WriteString[
													particleLogFile, 
													StringForm[
														"cellListSolubleCellsList is empty for `` in step round `` and iteration ``\n",
														ToString[particlePair],	ToString[n], ToString[m]
													]
												];	
												WriteString[
													particleLogFile, 
													StringForm[
														"Memory in use [b]: ``\n",
														ToString[N[UnitConvert[Quantity[MemoryInUse[], "Bytes"], "Megabytes"]]]
													]
												];													
											];
			                    			cellListHead = cellListResult[[1]];
						                   	cellListList = cellListResult[[2]];
						                   	cellListSolubleCellsList = DeleteDuplicates[cellListResult[[3]] ]; 
											If[
												cellListResult[[3]] == {},
												WriteString[
													particleLogFile, 
													StringForm[
														"cellListSolubleCellsList is empty for `` in step round `` and iteration `` after DeleteDuplicates\n",
														ToString[particlePair],	ToString[n], ToString[m]
													]
												];	
												WriteString[
													particleLogFile, 
													StringForm[
														"Memory in use [b]: ``\n",
														ToString[N[UnitConvert[Quantity[MemoryInUse[], "Bytes"], "Megabytes"]]]
													]
												];													
											];
						                   	(* neighbor cells for each cell that contains an atom of the soluted molecule; including the cell itself *)
						                   	cellListNeighborCells = DeleteDuplicates[Flatten[GetNeighborCells[#, cellListCellsPerRow, cellListIndexCube]&/@cellListSolubleCellsList]];
						                   	(* atoms/items contained in the neighborCells *)
											cellListNeighborCellsItems = DeleteCases[DeleteDuplicates[Flatten[GetCellItems[cellListHead, cellListList, #]&/@cellListNeighborCells]],0];
											(* atoms/items of solvent molecules contained in the neighborCells *)
											cellListNeighborCellsSolventItems = Drop[Sort[cellListNeighborCellsItems], nSolubleAtoms];
						                   	(* indices of neighbor molecules *)
			                    			Sow[
			                    				GetNeighborParticleNumbers[
													xyzDataListInput,
							                   		allAtomsNames,
							                   		nSolubleAtoms,
							                   		cellListNeighborCellsSolventItems,
													catchRadius[[1]],
													boxLength,
													vdWRadii					
												]
			                    			],
			                    			{n, nStepsPerRound}
		                    			];
		                    			Close[stream];
		                    			stream = OpenRead[particlePairString <> ".arc"];
		                    			linesToSkip = (nAllAtoms + 2) * (nStepsPerRound - 1);
		                    			Skip[stream, Record, linesToSkip];
		                    			xyzDataList = ReadList[stream, String, (nAllAtoms+2)];
		                    			Close[stream];
		                    			Export[dataName, xyzDataList, "List"];
		                    			DeleteFile[particlePairString <> ".arc"];
		                    			(* Export last step *)
		                    			If[
						                   	m == nIterations,
						                   	stream = OpenRead[dataName];
						                   	xyzDataList = ReadList[stream, String, (nAllAtoms+2)];
			                    			Export[FileNameJoin[{outputPath, FileBaseName[dataName] <> "_lastCoords.txt"}], xyzDataList, "List"];  
			                    			Close[stream];
			                    		],
		                    			{m, nIterations}
		                    		]
		                    	]
		                    ],1
		                ];
					If[
						neighborIndices == {},
						WriteString[
							particleLogFile, 
							StringForm[
								"neighborIndices is empty for ``\n",
								ToString[particlePair]
							]
						];			
						WriteString[
							particleLogFile, 
							StringForm[
								"Memory in use [b]: ``\n",
								ToString[N[UnitConvert[Quantity[MemoryInUse[], "Bytes"], "Megabytes"]]]
							]
						];											
					];					
                    neighborIndices = Sort[DeleteCases[DeleteDuplicates[Flatten[#]], 0]] & /@ neighborIndices;
					If[
						neighborIndices == {},
						WriteString[
							particleLogFile, 
							StringForm[
								"neighborIndices is empty for `` after DeleteDuplicates\n",
								ToString[particlePair]
							]
						];	
						WriteString[
							particleLogFile, 
							StringForm[
								"Memory in use [b]: ``\n",
								ToString[N[UnitConvert[Quantity[MemoryInUse[], "Bytes"], "Megabytes"]]]
							]
						];													
					];
					WriteString[
							particleLogFile, 
							StringForm[
								"Memory in use [b]: ``\n",
								ToString[N[UnitConvert[Quantity[MemoryInUse[], "Bytes"], "Megabytes"]]]
							]
						];	
					neighborIndices = 
						DeleteDuplicates[
							Flatten[moleculeLookupTable[[#]]&/@##]
						]&/@neighborIndices,
					
                    (* false - brute force *)
                    neighborIndices = 
                    	Flatten[
                    		Last[
                    			Reap[
		                   			Do[
		                   				ReadList[
					                       "!" <>
					                       	tinkerDynamic <>
					                       	dataName <>
					                       	" " <>
					                       	ToString[nStepsPerRound] <>
					                       	" " <>
						                    ToString[timeStep] <>
						                    " " <>
						                    ToString[printInterval] <>
						                    " " <>
						                    ToString[simulationType] <>
						                    " " <>
						                    ToString[temperature] <>
						                    " > " <> FileNameJoin[{outputPath, FileBaseName[dataName] <> "dynamic.log"}],
						                    Record 
						                ];
						                stream = OpenRead[particlePairString <> ".arc"];
						                Do[						                    	
			                				Skip[stream, Record, 2];
			                				xyzDataListInput = ReadList[stream, {Number, Word, Real, Real, Real, Record}, nAllAtoms][[All, 3;;5]];
			                				Sow[
			                					GetNeighborParticleNumbers[
													xyzDataListInput,
										            allAtomsNames,
										            nSolubleAtoms,
										            allSolventAtomsIndices,
													catchRadius[[1]],
													boxLength,
													vdWRadii					
												]
		                    				],
		                    				{n, nStepsPerRound}
	                    				];										
	                    				Close[stream];
	                    				stream = OpenRead[particlePairString <> ".arc"];
	                    				linesToSkip = (nAllAtoms + 2) * (nStepsPerRound - 1);
	                    				Skip[stream, Record, linesToSkip];
	                    				xyzDataList = ReadList[stream, String, (nAllAtoms+2)];
	                    				Close[stream];
	                    				Export[dataName, xyzDataList, "List"];
	                    				DeleteFile[particlePairString <> ".arc"];
		                    			(* Export last step *)
	                    				If[
					                    	m == nIterations,
					                    	stream = OpenRead[dataName];
		                    				Export[FileNameJoin[{outputPath, FileBaseName[dataName] <> "_lastCoords.txt"}], ReadList[stream, String, (nAllAtoms+2)], "List"];  
		                    				Close[stream];
		                    			],
	                    				{m, nIterations}
	                    			]
	                    		]
	                    	],1
	                    ];
					If[
						neighborIndices == {},
						WriteString[
							particleLogFile, 
							StringForm[
								"neighborIndices is empty for ``\n",
								ToString[particlePair]
							]
						];												
					];
	                neighborIndices = DeleteCases[DeleteDuplicates[Flatten[#]], 0]&/@neighborIndices; 
					If[
						neighborIndices == {},
						WriteString[
							particleLogFile, 
							StringForm[
								"neighborIndices is empty for `` after DeleteDuplicates\n",
								ToString[particlePair]
							]
						];												
					];
                   	neighborIndices = 
						DeleteDuplicates[
							Flatten[moleculeLookupTable[[#]]&/@##]
						]&/@neighborIndices;
               	];      
    	        Export[FileNameJoin[{outputPath, FileBaseName[dataName] <> "_lastStepNeighbors.txt"}], Last[neighborIndices], "List"];                                              
				WriteString[
					particleLogFile, 
					StringForm["Time for solvent box analysis (includes import stepwise): ``\n",
					ToString[Now - dynamicAnalysisTime]]
				];						
				WriteString[particleLogFile, StringForm["Total time for determining the coordination number: ``\n", ToString[Now - coordinationCalculationTime]]];                                                          
                If[
                    exportZijShow,
                    (* true *)
                    If[
                        !DirectoryQ[FileNameJoin[{outputPath, "Charts", particlePairString, ToString[temperature], "Zij_Show"}]],
                        CreateDirectory[FileNameJoin[{outputPath, "Charts", particlePairString, ToString[temperature], "Zij_Show"}]]
                    ];
                    ReinstallJava[ClassPath -> cdkJarFullPath];
                    LoadJavaClass["org.openscience.cdk.tools.periodictable.PeriodicTable"];
                    solubleAtomNamesString = Table[PeriodicTable`getSymbol[solubleAtomNames[[k]]], {k, Length[solubleAtomNames]}];
                    solventAtomNamesString = Table[PeriodicTable`getSymbol[solventAtomNames[[k]]], {k, Length[solventAtomNames]}];
                    xyzFirst = ConvertXYZDataToXYZ[xyzDataList[[1]], solubleAtomNamesString, solventAtomNamesString, nAllAtoms, nDynamicIteration];
                    xyzLast = ConvertXYZDataToXYZ[Last[xyzDataList], solubleAtomNamesString, solventAtomNamesString, nAllAtoms, nDynamicIteration];
                    solubleCoordsFirst = GetSolubleCoordinates[xyzDataList[[1]], nDynamicIteration, solubleAtomNamesString];
                    solubleCoordsFirstDetail = TableForm[solubleCoordsFirst, TableSpacing -> {0,4}]//ToString;
                    solubleCoordsFirstDetail = ImportString[solubleCoordsFirstDetail, "XYZ"];
                    solubleCoordsLast = GetSolubleCoordinates[Last[xyzDataList], nDynamicIteration, solubleAtomNamesString];
                    solubleCoordsLastDetail = TableForm[solubleCoordsLast, TableSpacing -> {0, 4}] // ToString;
                    solubleCoordsLastDetail = ImportString[solubleCoordsLastDetail, "XYZ"];
                    Do[
                        solventCoordsFirst = GetSolventCoordinates[xyzDataList[[1]], neighborIndices[[l, 1]], solubleAtomNamesString, solventAtomNamesString];
                        solventCoordsFirstDetail = TableForm[Flatten[solventCoordsFirst, 1], TableSpacing -> {0, 4}] // ToString;
                        solventCoordsFirstDetail = ImportString[solventCoordsFirstDetail, "XYZ"];
                        solventCoordsLast = GetSolventCoordinates[
                        xyzDataList[[nDynamicIteration]], neighborIndices[[l, nDynamicIteration]], solubleAtomNamesString, solventAtomNamesString];
                        solventCoordsLastDetail = TableForm[Flatten[solventCoordsLast, 1], TableSpacing -> {0, 4}] // ToString;
                        solventCoordsLastDetail = ImportString[solventCoordsLastDetail, "XYZ"];
                        notebook = 
							CreateDocument[{
								Show[
									MoleculePlot3D[xyzFirst],
									Graphics3D[{
										Yellow,
										Opacity[0.5],
										Table[
											Sphere[
												solubleCoordsFirst[[m]][[2;;4]],
												(catchRadius[[l]] * 0.5 + vdWRadii[[elementSymbolToAtomicNumberTable[solubleCoordsFirst[[m, 1]] ]]] )														
											],
											{m, Length[solubleAtomNames]}
										]													
									}],
									Graphics3D[{
										LightGray,
										Opacity[0.5],
										Table[
											Sphere[
												solventCoordsFirst[[m, n]][[2;;4]], 
												(catchRadius[[l]] * 0.5 + vdWRadii[[elementSymbolToAtomicNumberTable[solventCoordsFirst[[m, n]][[1]] ]]] )	
											],
											{m, Length[solventCoordsFirst]}, {n, Length[solventCoordsFirst[[m]] ]}
										]													
									}],
									Boxed -> True												
								],
								Show[
									MoleculePlot3D[solubleCoordsFirstDetail],
									MoleculePlot3D[solventCoordsFirstDetail],
									Graphics3D[{
										Yellow,
										Opacity[0.5],
										Table[
											Sphere[
												solubleCoordsFirst[[m]][[2;;4]],
												(catchRadius[[l]] * 0.5 + vdWRadii[[elementSymbolToAtomicNumberTable[solubleCoordsFirst[[m]][[1]] ]]] )	
											],
											{m, Length[solubleAtomNames]}
										]													
									}],
									Graphics3D[{
										LightGray,
										Opacity[0.5],
										Table[
											Sphere[
												solventCoordsFirst[[m, n]][[2;;4]],
												(catchRadius[[l]] * 0.5 + vdWRadii[[elementSymbolToAtomicNumberTable[solventCoordsFirst[[m, n]][[1]] ]]] )	
											],
											{m, Length[solventCoordsFirst]}, {n, Length[solventCoordsFirst[[m]] ]}
										]													
									}],
									Boxed -> True
								],											
								Show[
									MoleculePlot3D[xyzLast],
									Graphics3D[{
										Yellow,
										Opacity[0.5],
										Table[
											Sphere[
												solubleCoordsLast[[m]][[2;;4]],
												(catchRadius[[l]] * 0.5 + vdWRadii[[elementSymbolToAtomicNumberTable[solubleCoordsLast[[m]][[1]] ]]] )	
											],
											{m, Length[solubleAtomNames]}
										]													
									}],
									Graphics3D[{
										LightGray,
										Opacity[0.5],
										Table[
											Sphere[
												solventCoordsLast[[m, n]][[2;;4]],
												(catchRadius[[l]] * 0.5 + vdWRadii[[elementSymbolToAtomicNumberTable[solventCoordsLast[[m, n]][[1]] ]]] )	
											],
											{m, Length[solventCoordsLast]}, {n, Length[solventCoordsLast[[m]] ]}
										]													
									}],
									Boxed -> True												
								],
								Show[
									MoleculePlot3D[solubleCoordsLastDetail],
									MoleculePlot3D[solventCoordsLastDetail],
									Graphics3D[{
										Yellow,
										Opacity[0.5],
										Table[
											Sphere[
												solubleCoordsLast[[m]][[2;;4]],
												(catchRadius[[l]] * 0.5 + vdWRadii[[elementSymbolToAtomicNumberTable[solubleCoordsLast[[m]][[1]] ]]] )	
											],
											{m, Length[solubleAtomNames]}
										]													
									}],
									Graphics3D[{
										LightGray,
										Opacity[0.5],
										Table[
											Sphere[
												solventCoordsLast[[m, n]][[2;;4]],
												(catchRadius[[l]] * 0.5 + vdWRadii[[elementSymbolToAtomicNumberTable[solventCoordsLast[[m, n]][[1]] ]]] )	
											],
											{m, Length[solventCoordsLast]}, {n, Length[solventCoordsLast[[m]] ]}
										]													
									}],
									Boxed -> True												
								],
								Visible->False
							}];	
		 					NotebookEvaluate[notebook, InsertResults -> True];
							NotebookSave[notebook, FileNameJoin[{outputPath, "Charts",particlePairString, "_", ToString[temperature], "Zij_Show", StringJoin[{"Zij_show", ToString[catchRadius[[l]] ], ".nb"}]}]];
							NotebookClose[notebook],
                        {l, Length[catchRadius]}
                    ];                       
                ];     
                WriteString[particleLogFile, "Before mean loop"];  
                Export[FileNameJoin[{outputPath, "warmupNeighborIndices.txt"}], warmUpNeighborIndices];    
                Export[FileNameJoin[{outputPath, "neighborIndices.txt"}], neighborIndices];             
                Do[
                  	neighborParticleNumbers[[j, k]] = 
						Flatten[
							Reap[
								Sow[
									Length[warmUpNeighborIndices ]
								];
								Do[
									Sow[
										Length[neighborIndices[[l]] ]
									],
									{l, nDynamicIteration}
								]
							][[2]]
						];
					neighborParticleMean[[j, k]] = Accumulate[neighborParticleNumbers[[j, k]] ]/Range[1, Length[neighborParticleNumbers[[j, k]] ] ],
					{k, Length[catchRadius]}													
				];	 
				WriteString[particleLogFile, "After mean loop"];  	                  
                If[
                   nRunIteration == 1,
                   outputData[[inputdirectoryNamePosition, 7]] = "vdwVolume(" <> solventName <> ") [\[CapitalARing]^3] = " <> ToString[vdwSolventVolume],
                   outputData[[inputdirectoryNamePosition, 8]] = "vdwVolume(" <> solventName <> ") [\[CapitalARing]^3] = " <> ToString[vdwSolventVolume]
                ];
                indexShift = (nRunIteration - 1) * 7;          
                outputData[[inputdirectoryNamePosition, 11 + indexShift]] = 
                    "Mean(" <> solubleName <> "/" <> solventName <> ") = " <> 
                    TextString[
                        Table[
                            N[Mean[neighborParticleNumbers[[j, l]]], 5], 
                            {l, Length[neighborParticleNumbers[[j]]]}
                        ], 
                        ListFormat -> {"{", ", ", "}"}
                    ];
                WriteString[particleLogFile, StringForm["Mean neighbor: ``\n", outputData[[inputdirectoryNamePosition, 11 + indexShift]]]];               
                outputData[[inputdirectoryNamePosition, 12 + indexShift]] = "boxLength [\[CapitalARing]] = " <> ToString[boxLength];
                outputData[[inputdirectoryNamePosition, 13 + indexShift]] = 
                    "\[Sigma] (" <> solubleName <> "/" <> solventName <> ") = " <>
                    TextString[
                        Table[
                            N[StandardDeviation[neighborParticleNumbers[[j, l]]], 5], 
                            {l, Length[neighborParticleNumbers[[j]]]}
                        ], 
                        ListFormat -> {"{", ", ", "}"}
                    ];
                outputData[[inputdirectoryNamePosition, 14 + indexShift]] = 
					"Min (" <> solubleName <> "/" <> solventName <> ") = " <>
                    TextString[
                        Table[
                            N[Min[neighborParticleNumbers[[j, l]]], 5], 
                            {l, Length[neighborParticleNumbers[[j]]]}
                        ], 
                        ListFormat -> {"{", ", ", "}"}
                    ];
                outputData[[inputdirectoryNamePosition, 15 + indexShift]] = 
                    "Max (" <> solubleName <> "/" <> solventName <> ") = " <>
                    TextString[
                        Table[
                            N[Max[neighborParticleNumbers[[j, l]]], 5], 
                            {l, Length[neighborParticleNumbers[[j]]]}
                        ], 
                        ListFormat -> {"{", ", ", "}"}
					];
                totalHistogram[[inputdirectoryNamePosition, nRunIteration]] = 
                    Table[
                        Histogram[
                            neighborParticleNumbers[[j, l]], 
                            ChartLegends -> StringJoin[{ToString[temperature], "K"}], 
                            PlotRange -> All, 
                            PlotRangePadding -> Scaled[.05], 
                            Frame -> True, 
                            FrameLabel -> {"Zij", None}, 
                            PlotTheme -> "Scientific", 
                            ChartStyle -> RGBColor["#3277a8"], 
                            ImageSize -> 1000
                        ], 
                        {l, Length[neighborParticleNumbers[[j]]]}
                    ];
                averageCharts[[inputdirectoryNamePosition, nRunIteration]] = 
                    Table[
                        ListLinePlot[
                            Transpose[{
                                Table[
                                    IntegerPart[warmUpStepNumber + o], 
                                    {o, 0, Length[neighborParticleMean[[j, m]]] - 1}
                                ], #}
                            ]& /@ {neighborParticleMean[[j, m]]}, 
                            Frame -> True, 
                            PlotLegends -> Placed[StringJoin[{ToString[temperature], "K"}], Below], 
                            PlotLabel -> StringForm["Average value"], 
                            FrameLabel -> {"steps", "<Z>"}, 
                            PlotRange -> All, 
                            PlotRangePadding -> Scaled[.05], 
                            PlotTheme -> "Scientific", 
                        	PlotStyle -> RGBColor["#3277a8"], 
                            ImageSize -> 1000
                        ], 
                        {m, Length[catchRadius]}
                    ];   
                (* Export Charts *)                    
                If[
                    !DirectoryQ[FileNameJoin[{outputPath, "Charts",StringJoin[{solubleName, "_", solventName}], ToString[temperature]}] ],
                    CreateDirectory[FileNameJoin[{outputPath, "Charts", StringJoin[{solubleName, "_", solventName}], ToString[temperature]}]]
                ];
                Export[
                    FileNameJoin[{outputPath, "Charts", StringJoin[{solubleName, "_", solventName}], ToString[temperature], StringJoin[{"Histogram_", ToString[temperature], "K.pdf"}]}], 
                    totalHistogram[[inputdirectoryNamePosition, nRunIteration]]
				];                   
                If[
                    exportCharts,                       
                    histograms = 
                        Table[
                        	Histogram[
                                neighborParticleMean[[j, m]], 
                                ChartLegends -> Placed[StringJoin[{ToString[temperature], "K"}], Below], 
                                PlotRange -> All, 
                                PlotRangePadding -> Scaled[.05], 
                                Frame -> True, 
                                FrameLabel -> {"Z", None},
                                PlotTheme -> "Scientific",
                                ChartStyle -> RGBColor["#3277a8"], 
                                ImageSize -> 1000
                            ], 
                            {m, Length[neighborParticleNumbers[[j]]]}
                        ];                       
                    Do[
                        Export[
                        	FileNameJoin[{
                        		outputPath, "Charts", StringJoin[{solubleName, "_", solventName}], ToString[temperature], 
                                StringJoin[{"Histogram_", ToString[catchRadius[[m]] ], "_", ToString[temperature], "K.pdf"}]
                            }], 
                        	histograms[[m]]                              	
                        ],
                        {m, Length[catchRadius]}
                    ];
                	Do[
                        Export[
                        	FileNameJoin[{
                        		outputPath, "Charts", StringJoin[{solubleName, "_", solventName}], ToString[temperature], 
                                StringJoin[{"Average_", ToString[temperature], "K.pdf"}]
                            }], 
                            averageCharts[[m]]
                        ],
                        {m, Length[averageCharts[[inputdirectoryNamePosition, nRunIteration]]]}	
					]
                ];
                (* Delete scratch directory *)
                SetDirectory[scratchDirectory];
                Run["del *.* /Q"];
		        WriteString[particleLogFile, StringForm["Total time for particle pair: ``\n", ToString[Now - particlePairTime]]];
                particlePairTime =.;
                Close[particleLogFile],
                {j, Length[allParticlePairs]}
			];	
			CloseKernels[];
			WriteString[gblLogFile, StringForm["Coordination number calculation time: ``\n\n", ToString[Now - gblCoordinationTime]] ];
			(* Generates mathematica notebook *)
			Table[
				particlePair = inputdirectoryName[[i]];
				If[
					StringSplit[particlePair, "_"][[1]] == StringSplit[particlePair, "_"][[2]],
					isSameParticle = True,
					isSameParticle = False									
				];
				outputPath = FileNameJoin[{resultDirectory, particlePair}];
				outputData[[i, 4]] = "Time to calculate minimum intermolecular energy = " <> ToString[energyCalculationTime];
				outputData[[i, 5]] = "sphereNodeNumber = " <> ToString[sphereNodeNumber];
				outputData[[i, 6]] = "circularRotationNumber = " <> ToString[rotationNumber];
				SetDirectory[FileNameJoin[{resultDirectory, particlePair}]];
					(*SetOptions[ListLinePlot, ImageSize -> {425, Automatic}];										
						xyLinePlots = 
							Table[
								ListLinePlot[
											Transpose@{xyData[[i, 1]], xyData[[i, 2, j, 1]]},
											Frame -> True,
											FrameLabel -> {"Distance [\[CapitalARing]]", "Intermolecular Energy [kcal/mole]"},
											PlotLabel -> particlePair <> "_" <> ToString[temperatures[[j]] ] <> "K",
											PlotMarkers -> Automatic,
											PlotRangePadding -> Scaled[.05]
										],
								{j, Length[temperatures]}
							];
						preciseXyLinePlots = 
							Table[
								ListLinePlot[
											Transpose@{preciseXyData[[i, 1]], preciseXyData[[i, 2, j, 1]]},
											Frame -> True,
											FrameLabel -> {"Distance [\[CapitalARing]]", "Intermolecular Energy [kcal/mole]"},
											PlotLabel -> particlePair <> "_" <> ToString[temperatures[[j]] ] <> "K",
											PlotMarkers -> Automatic,
											PlotRangePadding -> Scaled[.05]
										],
								{j, Length[temperatures]}
							];	*)
				If[
					isSameParticle == True,
					(* true *)
					SetOptions[ListLinePlot, ImageSize -> {425, Automatic}];
					notebook = 
						CreateDocument[{
							TextCell["Intermolecular interaction", "Section"],
							TextCell["Particle configuration of lowest intermolecular energy", "Subsection"],
							MoleculePlot3D[
								Import[FileNameJoin[{outputPath, "output.xyz"}]]
							],
							TableForm[energyListLinePlots[[i]], TableDirections -> Row],
							(*TableForm[energyPreciseListLinePlots[[i]], TableDirections -> Row],*)
							outputData[[i, 1]],
							outputData[[i, 2]],
							outputData[[i, 3]],
							outputData[[i, 4]],
							outputData[[i, 5]],
							outputData[[i, 6]],
							outputData[[i, 7]],
							outputData[[i, 8]],
							TextCell["Optimized minimum energy configuration","Section"],
							Import[
								FileNameJoin[{outputPath, "output_optimized.pdb"}],
								ImageResolution -> 1000
							],
							outputData[[i, 9]],
							outputData[[i, 10]],
							TextCell["Coordination number", "Section"],
							outputData[[i, 11]],
							outputData[[i, 12]],
							outputData[[i, 13]],
							outputData[[i, 14]],
							outputData[[i, 15]],
							(*totalHistogram[[i]]*)
							averageCharts[[i]]
						},
						Visible->False
					],
					(* false *)
					notebook = 
						CreateDocument[{
							TextCell["Intermolecular interaction", "Section"],
							TextCell["Particle configuration of lowest intermolecular energy", "Subsection"],
							MoleculePlot3D[
								Import[FileNameJoin[{outputPath, "output.xyz"}]]
							],
							TableForm[energyListLinePlots[[i]], TableDirections -> Row],
							(*TableForm[energyPreciseListLinePlots[[i]], TableDirections -> Row],*)
							outputData[[i, 1]],
							outputData[[i, 2]],
							outputData[[i, 3]],
							outputData[[i, 4]],
							outputData[[i, 5]],
							outputData[[i, 6]],
							outputData[[i, 7]],
							outputData[[i, 8]],
							TextCell["Optimized minimum energy configuration","Section"],
							Import[
								FileNameJoin[{outputPath, "output_optimized.pdb"}],
								ImageResolution -> 1000
							],
							outputData[[i, 9]],
							outputData[[i, 10]],
							TextCell["Coordination number", "Section"],
							TextCell[
								StringTake[
									outputData[[i, 11]],
									{
										StringPosition[outputData[[i, 11]], "("][[1, 1]] + 1, 
										StringPosition[outputData[[i, 11]], ")"][[1, 1]] - 1
									}
								], 
								"Text"
							],
							outputData[[i, 11]],
							outputData[[i, 12]],
							outputData[[i, 13]],
							outputData[[i, 14]],
							outputData[[i, 15]],
							(*totalHistogram[[i, 1]],*)
							averageCharts[[i,1]],
							outputData[[i, 16]],
							outputData[[i, 17]],
							TextCell[
								StringTake[
									outputData[[i, 18]],
									{
										StringPosition[outputData[[i, 18]], "("][[1, 1]] + 1, 
										StringPosition[outputData[[i, 18]], ")"][[1, 1]] - 1
									}
								], 
								"Text"
							],
							outputData[[i, 18]],
							outputData[[i, 19]],
							outputData[[i, 20]],
							outputData[[i, 21]],
							outputData[[i, 22]],
							outputData[[i, 23]],
							outputData[[i, 24]],
							(*totalHistogram[[i, 2]]*)
							averageCharts[[i,2]]
						},
						Visible->False
					]
				];
				NotebookEvaluate[notebook, InsertResults -> True];
				NotebookSave[notebook, FileNameJoin[{outputPath, particlePair <> ".nb"}]];
				Export[particlePair <> ".jpeg", notebook, "JPEG", "CompressionLevel" -> 0];
				NotebookClose[notebook];
				Export[particlePair <> ".dat", outputData[[i]], "Table"],					
				{i, Length[inputdirectoryName]}
			];
			(* Create overview notebook *)
			nbFilesList = Select[FileNames["*.nb", resultDirectory, Infinity] , # != "Overview.nb" &];
			If[
				createOverview && Length[nbFilesList] > 1,
				combinedNb = 
					CreateDocument[
						{},
						ShowPageBreaks -> True, 
						WindowTitle -> "Overview", 
						Visible -> False
					];
				Do[
					tempNb = NotebookOpen[nbFile, Visible -> False];
					SelectionMove[tempNb, All, Notebook];
					NotebookWrite[combinedNb, NotebookRead /@SelectedCells[tempNb]];
					SelectionMove[combinedNb, After, Notebook];
					SelectionMove[combinedNb, Previous, Cell];
					SetOptions[
						NotebookSelection[combinedNb], 
						PageBreakBelow -> True
					];
					NotebookClose[tempNb],
					{nbFile, nbFilesList}
				];
				NotebookSave[combinedNb, FileNameJoin[{resultDirectory, "Overview.nb"}]];
				Export[
					FileNameJoin[{resultDirectory, "Overview.jpeg"}],
					combinedNb, 
					"JPEG", 
					"CompressionLevel" -> 0
				];
				NotebookClose[combinedNb];
			];
			ExportParameterSet[
				parameterFileTitle,
				parameterFileTitleAbr, 
				forceField, 
				resultDirectory, 
				smiles, 
				cdkJarFullPath,
				waterVolumeRatio, 
				boltzmannFraction, 
				catchRadius, 
				weightedMinEnergyTextString,
				{cpuCoreNumber, useFibonacciSphereAlgorithm, sphereNodeNumber, rotationNumber, 
					boltzmannFraction, useCellList, solventMoleculeNumber, warmUpStepNumber, nDynamicSteps, catchRadius}
				];		
			ExportParticleSetSphereApproximantion[parameterFileTitle, parameterFileTitleAbr, resultDirectory, smiles, waterVolumeRatio, cdkJarFullPath, boltzmannFraction, weightedMinEnergyTextString];
			WriteString[gblLogFile, StringForm["Entire calculation Time: ``\n\n", ToString[Now - totalTime]] ];
			Close[gblLogFile];
			(* Clear scratch dir *)
			SetDirectory[scratchDirectory];
			dirs = Select[FileNames["*", "", Infinity], DirectoryQ];
			Do[
				SetDirectory[dirs[[1]] ];
				files = FileNames["*", "", Infinity];
				DeleteFile[#]&/@files, 	
				{i, Length[dirs]}
			];
			SetDirectory[scratchDirectory];
			DeleteDirectory[#]&/@dirs;
			(* Delete input directory *)
			SetDirectory[inputDirectory];
			dirs = Select[FileNames["*", "", Infinity], DirectoryQ];
			Do[
				SetDirectory[dirs[[1]] ];
				files = FileNames["*", "", Infinity];
				DeleteFile[#]&/@files, 	
				{i, Length[dirs]}
			];
			SetDirectory[inputDirectory];
			DeleteDirectory[#]&/@dirs;
		    (*SetDirectory[inputDirectory];
		     DeleteDirectory[#, DeleteContents -> True]&/@inputdirectoryName;	       *) 			 
			Return[];				
		];  (*End*)		
(* ::Section:: *)
(*End of Package*)
End[]
EndPackage[]