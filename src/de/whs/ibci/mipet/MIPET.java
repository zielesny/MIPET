/**
 * MIPET - Mesoscopic Interaction Parameter Estimation with Tinker
 * Copyright (C) 2025  Achim Zielesny (achim.zielesny@googlemail.com)
 * 
 * Source code is available at <https://github.com/zielesny/MIPET>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.whs.ibci.mipet;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Stream;

/**
 * Class MIPET
 * @author Mirco Daniel, Flix Bänsch
 */
public class MIPET {
    
    // <editor-fold defaultstate="collapsed" desc="Final private class variables">
    /**
     * Operating system
     */
    private static final String OS = System.getProperty("os.name");
    
    /**
     * Angstrom character
     */
    private static final char ANGSTROM = 'A';
    
    /**
     * Superscript 3 character
     */
    private static final String SUPERSCRIPTTHREE = "^3";
    
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger
            .getLogger(MIPET.class.getName());
    
    /**
     * Line separator
     */
    private static final String LINESEPARATOR = System.lineSeparator();
    
    /**
     * File separator
     */
    private static final String FILESEPARATOR = System
            .getProperty("file.separator");
    
    /**
     * Gas constant R [kcal/(mol*K)]
     */
    private static final double GASCONST = 1.98720425864E-3;
    
    /**
     * Instance object of MIPETUtility
     */
    private static final MIPETUtility MIPETUTIL = new MIPETUtility();
    
    /**
     * Water volume ratio
     */
    private static final double WATERVOLUMERATIO = 30.0036 / 
            MIPETUTIL.getVdwVolume("O");
    
     /**
     * Double for the DPD mass
     */
    private static final int MASS_DPD = 1;
    
    /**
     * Double for the value of the graphical radius of a particle
     */
    private static final double GRAPHICS_RADIUS = 0.5;
    
    /**
     * Double for the default charge of a particle
     */
    private static final int DEFAULT_CHARGE = 0;
    
    /**
     * String for the default color of a particle
     */
    private static final String DEFAULT_COLOR = "black";
    
    /**
     * String for version number of exported parameter set file
     */
    private static final String VERSION_NUMBER = "1.0.0.0";
    
    /**
     * Weighted intermolecular differential pair interaction energies
     */
    private static final ArrayList<ResultEnergyRecord> energyList = 
            new ArrayList<>();
    
    /**
     * Coordination number (average)
     */
    private static final ArrayList<ResultCNRecord> cnList = new ArrayList<>();
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    
    /**
     * Using Tinker9 flag
     */
    private static boolean isTinker9;
    
    /**
     * Using Tinker flag
     */
    private static boolean isTinkerOn;
    
    /**
     * Job file name
     */
    private static String jobFileName;
    
    /**
     * Particlenames
     */
    private static LinkedList<String> particleNames;
    
    /** 
     * MoleculeRecord
     */
    private static LinkedList<MoleculeRecord> molecules;
    
    /**
     * New particles for calculation
     */
    private static LinkedList<String> newParticles;
    
    /**
     * Old particles already calculated
     */
    private static LinkedList<String> oldParticles;
    
    /**
     * .xyz content
     */
    private static String[] xyzContent1;
    
    /**
     * .xyz alternate content
     */
    private static String[]xyzContent2;
    
    /**
     * .prm content
     */
    private static String[] prmContent1;
    
    /**
     * .prm alternate content
     */
    private static String[] prmContent2;
    
    /**
     * CPUcore number
     */
    private static int cpuCoreNumber;
    
    /**
     * Use Conformational Analysis
     */
    private static boolean isConformationalAnalysis;
    
    /**
     * Use fibonaccisphere algorithm
     */
    private static boolean isFibonacciSphereAlgorithm;

    /**
     * Use optimize sampled E(min) configuration
     */
    private static boolean isOptEmin;
    
    /**
     * Temperature
     */
    private static double temperature;
    
    /**
     * Sphere node number 
     */
    private static int sphereNodeNumber1;
    private static int sphereNodeNumber2;
    private static int sphereNodeNumber3;
    private static int sphereNodeNumber4;
    
    /**
     * Rotation configuration number
     */
    private static int rotationNumber1;
    private static int rotationNumber2;
    private static int rotationNumber3;
    private static int rotationNumber4;
    
    /**
     * Force field names for calculation of intermolecular interaction energy
     */
    private static String forcefield_IE;
    
    /**
     * Force field names for calculation of coordination number
     */
    private static String forcefield_CN;
    
    /**
     * Water model name
     */
    private static String watermodel;
    
    /**
     * ConfigurationNumber of particle 1
     */
    private static int confNumber1;
    
    /**
     * ConfigurationNumber of particle 2
     */
    private static int confNumber2;
    
    /**
     * Minimum atom distance to prevent tinker error output
     */
    private static double minAtomDistance;
    
    /**
     * Tinker's scan program
     * 0: automatic selection
     * 1: manual selection of angles to rotate
     * 2: manual selection of angles to freeze
     */
    private static String scanProgram;
    
    /**
     * Number of search directions in tinker's scan.exe
     */
    private static String nSearchDirection;
    
    /**
     * Energy threshold value in tinker's scan.exe
     */
    private static String energyThreshold;
    
    /**
     * Rms gradient value in tinker's scan.exe
     */
    private static String rmsGradient;
    
    /**
     * Rms minimize gradient value in tinker's minimize.exe
     */
    private static double rmsMinimizeGradient;
    
    /**
     * Maximum iteration value for MAXITER in .key-file
     */
    private static int minimizeMaxIteration;
    
    /**
     * Standard .key file content
     */
    private static String keyFileStringOrigin;
    
    /**
     * Dielectric constant
     */
    private static double dielectricConstant;
    
    /**
     * Directory of tinker files
     */
    private static String tinkerDirectory;
    
    /**
     * Directory of .prm files
     */
    private static String parameterDirectory;
    
    /**
     * Directory of scratch files
     */
    private static String scratchDirectory;
    
    /**
     * Directory of SMILES files
     */
    private static String smilesDirectory;
    
    /**
     * SPIECE and SMILES lookup table
     */
    private static HashMap<String, String> smiles;
    
    /**
     * Directory of calculation related datas
     */
    private static String calculationDirectory;
    
    /**
     * Directory name of particle datas
     */
    private static String sourceDirectory;
    
    /**
     * Directory name of optimized xyz datas
     */
    private static String optXYZDirectory;
    
    /**
     * Directory name of optimal distance of two particles
     */
    private static String optDistDirectory;
    
    /**
     * Directory name of output datas
     */
    private static String resultDirectory;
    
    /**
     * BufferedWriter object for global log file
     */
    private static BufferedWriter BFGblLog;
    
    /**
     * BufferedWriter for log contents
     */
    private static BufferedWriter BWParticleLog;
    
    /**
     * BufferedWriter for dat contents
     */
    private static BufferedWriter BWParticleDat;
    
    /**
     * Lower boundary of distances of two particles
     */
    private static double lowerBoundary;
    
    /**
     * Upper boundary of distances of two particles
     */
    private static double upperBoundary;
    
    /**
     * Stepsize of Prescan
     */
    private static double prescanStepSize;
    
    /**
     * Fraction of energy values used for the Boltzmann distribution
     */
    private static double boltzmannFraction;
    
    /**
     * Solvent particle number
     */
    private static int solventMoleculeNumber;
    
    /**
     * Step number for tinker's tool dynamic.exe
     */
    private static int stepNumber;
    
    /**
     * Step number per CPU core number
     */
    private static int dynamicStepsPerCore;
    
    /**
     * Timestep in fs for tinker's dynamic.exe
     */
    private static double timeStep;
    
    /**
     * Print interval in ps for tinker's dynamic.exe
     */
    private static double printInterval;
    
    /**
     * Tinker Dynamic Simulation Type
     * Available Statistical Mechanical Ensembles:
     * (1) Microcanonical (NVE)
     * (2) Canonical (NVT)
     * (3) Isoenthalpic-Isobaric (NPH)
     * (4) Isothermal-Isobaric (NPT)
     */
    private static int simulationType;
    
    /**
     * Catch radius in Angstrom
     */
    private static double catchRadius;
    
    /**
     * RMS Gradient for Tinker Optimize - to optimize particles xyz file
     */
    private static double optimizeRmsGradient;
    
    /**
     * WarmUp step number
     */
    private static int warmUpStepNumber;
    
    /**
     * Warmup time step for the simulation warmup
     */
    private static double warmUpTimeStep;
    
    /**
     * Warmup print interval 
     */
    private static double warmUpPrintInterval;
    
    /**
     * Warmup iteration timeStep in fs and printIntervall in ps
     */
    private static int nDynamicWarmupIteration;
    
    /**
     * Iteration for Tinker's dynamic.exe
     */
    private static int nDynamicIteration;
    
    /**
     * Iteration number per core for Tinker's dynamic.exe
     */
    private static int nDynamicIterationPerCore;
    
    /**
     * Path and program name of Tinker's analyze.exe
     */
    private static String tinkerAnalyze;
    
    /**
     * Path and program name of Tinker's scan.exe
     */
    private static String tinkerScan;
    
    /**
     * Path and program name of Tinker's xyzedit.exe
     */
    private static String tinkerXYZEdit;
    
    /**
     * Path and program name of Tinker's minimize.exe
     */
    private static String tinkerMinimize;
    
    /**
     * Path and program name of Tinker's dynamic.exe
     */
    private static String tinkerDynamic;
    
    /**
     * Path and program name of Tinker's optimize.exe
     */
    private static String tinkerOptimize;
    
    /**
     * Path and program name of Tinker's optrigid.exe
     */
    private static String tinkerOptrigid;
    
    /**
     * Path and program name of Tinker's xyzpdb.exe
     */
    private static String tinkerXYZPdb;
    
    /**
     * Parameterset title
     */
    private static String parameterSetTitle;
    
    /**
     * Parameterset title abbreviation
     */
    private static String parameterSetTitleAbr;
    
    /**
     * Flag for write wamup log file
     */
    private static boolean isLogWarmup;
    
    /**
     * Flag for write dynamic log file
     */
    private static boolean isLogDynamic;
    
    /**
     * Flag for write build simulation box log file
     */
    private static boolean isLogBuildBox;
    
    /**
     * Flag for write minimize simulation box log file
     */
    private static boolean isLogMinimizeBox;
    
    /**
     * Flag for write soak simulation box log file
     */
    private static boolean isLogSoakBox;
            
    /**
     * Flag for write RATTLE key word to the .key file
     */
    private static boolean isRATTLE;
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    
    /**
     * Main method
     * @param args
     */
    public static void main(String[] args) {
        //<editor-fold defaultstate="collapsed" desc="Variables declaration">
        boolean tmpIsExitCondition;
        boolean tmpIsSameParticle;
        boolean tmpHasEnergieJob;
        int tmpCurrentIndex;
        String tmpForcefield;
        String tmpParticlePair;
        String tmpKeyFileString;
        String tmpKeyFileName;
        String tmpKeyPathName;
        String tmpIEResultDirName;
        String tmpCNResultDirName;
        String tmpParticleDatFileName;
        Path tmpParticleDatFile;
        String tmpParticleName1;
        String tmpParticleName2;
        double[][] tmpXyzData1; // coordination of 1. fragment atoms
        double[][] tmpXyzData2; // coordination of 2. fragment atoms
        double[] tmpCentre1;
        double[] tmpCentre2;
        double[][][] tmpXyzRotData1; // i: configuration, j: atom, k: xyz
        double[][][] tmpXyzRotData2; // i: configuration, j: atom, k: xyz
        String[] tmpElementList1;
        String[] tmpElementList2;
        String[] tmpCmdList;
        TinkerXYZ tmpTinkerXYZ1;
        TinkerXYZ tmpTinkerXYZ2;
        Process tmpProcess;
        ProcessBuilder tmpPB;
        ArrayIndexComparator tmpComparator;
        ArrayList<String> tmpParticlePairs = new ArrayList<>(500);
        ArrayList<JobTaskRecord> tmpJobTaskRecordList = new ArrayList<>(500);
        List<double[]> tmpSphereNodeCoord;

        //</editor-fold>
        
        // <editor-fold defaultstate="collapsed" desc="Initialize and read .job file">
        Locale.setDefault(Locale.ENGLISH);
        long tmpTotalTime = System.currentTimeMillis();
        System.out.println("Initializing...");
        initialize();
        System.out.println("Reading job file...");
        readJobFile();
        
        //</editor-fold>
        
        // <editor-fold defaultstate="collapsed" desc="Set original .key file">
        if (forcefield_IE.equals("OPLSAALIGPARGEN")) {
            keyFileStringOrigin =
                "# Force Field Selection"
                + LINESEPARATOR
                + "PARAMETERS\t\""
                + parameterDirectory
                + "/"
                + "oplsaa.prm\""
                + LINESEPARATOR;
        } else {
            keyFileStringOrigin =
                "# Force Field Selection"
                + LINESEPARATOR
                + "PARAMETERS\t\""
                + parameterDirectory
                + "/"
                + forcefield_IE.toLowerCase() + ".prm\""
                + LINESEPARATOR;
        }
        
        // </editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Create log file">
        String tmpGblLogDirName = resultDirectory
                + FILESEPARATOR
                + "IE"
                + FILESEPARATOR
                + forcefield_IE;
        Path tmpGblLogDirFile = Paths.get(tmpGblLogDirName);
        if (!Files.exists(tmpGblLogDirFile)) {
            try {
                Files.createDirectories(tmpGblLogDirFile);
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE, 
                    "IOException during creating global log file.", ex);
            }
        }
        String tmpGblLogFileName = resultDirectory 
                + FILESEPARATOR
                + "IE"
                + FILESEPARATOR
                + forcefield_IE
                + FILESEPARATOR
                + "log.txt";
        Path tmpLogFile = Paths.get(tmpGblLogFileName);
        try {
            BFGblLog = new BufferedWriter(new FileWriter(tmpGblLogFileName));
            if (!Files.exists(tmpLogFile)){
                Files.createFile(tmpLogFile);
            } else {
                BFGblLog.append(LINESEPARATOR);
            }
            BFGblLog.append("CPU cores: "
                    + cpuCoreNumber 
                    + LINESEPARATOR);
            BFGblLog.append("Warm up steps for dynamic simulation: " 
                    + warmUpStepNumber
                    + LINESEPARATOR);
            BFGblLog.append("Steps for dynamic simulation: " 
                    + stepNumber
                    + LINESEPARATOR);
            BFGblLog.flush();
        } catch(IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during writing in log file.", ex);
        }
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Load surface coordinates">
        /* The coordinates for equidistantly distributed points on a sphere 
           from Technical University of Dortmund are used, thanks to 
           J. Fliege and U. Maier
	   http://www.mathematik.uni-dortmund.de/lsx/research/projects/fliege/nodes/nodes.html */
        long tmpRotCalcTime = System.currentTimeMillis();
        
        try {
            BFGblLog.append ("Time for calculation of rotation matrices: " 
                + (System.currentTimeMillis() - tmpRotCalcTime)
                + LINESEPARATOR);
            BFGblLog.flush();
        } catch(IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during writing in log file.", ex);
        }
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Prepair input">
        tmpParticlePairs = getParticlePairs();
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Make Job task record">
        String tmpIEDatFileName;
        String tmpCNDatFileName;
        int tmpParticlePairLength;
        boolean tmpHasCNJob;
        
        tmpParticlePairLength = tmpParticlePairs.size();
        
        for (int i = 0; i < tmpParticlePairLength; i++) {
            tmpParticlePair = tmpParticlePairs.get(i);
            tmpParticleName1 = tmpParticlePair.split("_")[0];
            tmpParticleName2 = tmpParticlePair.split("_")[1];
            tmpIsSameParticle = tmpParticleName1.equals(tmpParticleName2);
            tmpIEResultDirName = resultDirectory
                    + FILESEPARATOR
                    + "IE"
                    + FILESEPARATOR
                    + forcefield_IE
                    + FILESEPARATOR
                    + tmpParticlePair;
            tmpIEDatFileName = tmpIEResultDirName
                    + FILESEPARATOR
                    + tmpParticlePair
                    + ".dat";
            tmpCNResultDirName = resultDirectory
                    + FILESEPARATOR
                    + "CN"
                    + FILESEPARATOR
                    + forcefield_CN
                    + FILESEPARATOR
                    + tmpParticlePair;
            tmpCNDatFileName = tmpCNResultDirName
                    + FILESEPARATOR
                    + tmpParticlePair
                    + ".dat";
            if (Files.exists(Paths.get(tmpIEDatFileName))) {
                tmpHasEnergieJob = false;
            } else {
                tmpHasEnergieJob = true;
                try {
                    Files.createDirectories(Paths.get(tmpIEResultDirName));
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during creating result directory.", ex);
                }
            }
            if (Files.exists(Paths.get(tmpCNDatFileName))) {
                tmpHasCNJob = false;
            } else {
                tmpHasCNJob = true;
                try {
                    if (!forcefield_CN.isEmpty()) {
                        Files.createDirectories(Paths.get(tmpCNResultDirName));
                    }
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during creating CN result directory.", ex);
                }
            }
            tmpJobTaskRecordList.add(new JobTaskRecord(
                    forcefield_IE,
                    forcefield_CN,
                    tmpParticleName1, 
                    tmpParticleName2, 
                    tmpIEResultDirName,
                    tmpCNResultDirName,
                    tmpHasEnergieJob,
                    tmpHasCNJob,
                    tmpIsSameParticle,
                    false));
            if (!tmpIsSameParticle) {
                tmpHasEnergieJob = false;
                tmpJobTaskRecordList.add(new JobTaskRecord(
                        forcefield_IE,
                        forcefield_CN,
                        tmpParticleName2, 
                        tmpParticleName1, 
                        tmpIEResultDirName,
                        tmpCNResultDirName,
                        tmpHasEnergieJob,
                        tmpHasCNJob,
                        tmpIsSameParticle,
                        true));
            }
        }
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Read .xyz and .prm">
        readXyz(true, forcefield_IE);
        readPrm(forcefield_IE);
        makeMoleculeRecord();
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Conformational analysis">
        String tmpSourceName;
        String tmpTargetName;
        String tmpFileName;
        Path tmpSourcePath;
        Path tmpTargetPath;
        Path tmpOptXyzPath;
        
        if (isConformationalAnalysis) {
            System.out.println("Conformational Analysis...");
            scanParticle();
            
            // Overwrite xyzContent with new .xyz
            readXyz(false, forcefield_IE);
        } else {
            System.out.println("Conformational Analysis skipped.");
            tmpForcefield = tmpJobTaskRecordList.get(0).forcefield_IE_Name();
            
            for (String tmpParticleName : particleNames) {
                tmpOptXyzPath = Paths.get(optXYZDirectory, 
                        tmpForcefield, 
                        tmpParticleName);
                if (!Files.exists(tmpOptXyzPath)) {
                    try {
                        Files.createDirectories(tmpOptXyzPath);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during creating OptXyz directory.", 
                                ex);
                    }
                }
            }
            
            // Copy .xyz file to optXYZ directory
            for (String tmpParticleName : particleNames) {
                tmpSourceName = sourceDirectory
                        + FILESEPARATOR 
                        + tmpForcefield 
                        + FILESEPARATOR 
                        + tmpParticleName
                        + ".xyz";
                tmpFileName = FILESEPARATOR 
                        + tmpForcefield 
                        + FILESEPARATOR 
                        + tmpParticleName
                        + FILESEPARATOR 
                        + tmpParticleName
                        + ".xyz";
                tmpTargetName = optXYZDirectory + tmpFileName;

                // Check whether already scanned
                tmpSourcePath = Paths.get(tmpSourceName);
                tmpTargetPath = Paths.get(tmpTargetName);
                if (!Files.exists(tmpTargetPath)) {
                    try{
                        Files.copy(tmpSourcePath, tmpTargetPath, 
                                StandardCopyOption.REPLACE_EXISTING);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during copying .xyz file.", ex);
                    }
                }
            }
            
        }
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Variable declarations">
        byte tmpH2OPos;
        int tmpXyz1ID;
        int tmpXyz2ID;
        int tmpDistSize;
        int tmpAtomSize1;
        int tmpPrmID1;
        int tmpPrmID2;
        long tmpEnergyCalcTime;
        double tmpEqDist;
        double tmpGlbEmin;
        double tmpGlbWgtEmin;
        double tmpDistanceCandidate;
        String tmpOutputName;
        String tmpLine;
        String tmpKeyContent;
        Path tmpKeyFile;
        Path tmpOptDistFile;
        DecimalFormat decimal2;
        DecimalFormat decimal3;
        DecimalFormat decimal4;
        EnergyRecord[] tmpEnergyRecords;
        LinkedList<Double> tmpAllDistances;
        LinkedList<Double> tmpDistanceList;
        LinkedList<Distance_EnergyRecord> tmpDistWgtEmins;
        LinkedList<Distance_EnergyRecord> tmpDistEminRecords;
        LinkedList<double[][][]> tmpRotCoords;
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Variable initializations">
        decimal2 = (DecimalFormat)NumberFormat.getNumberInstance();
        decimal3 = (DecimalFormat)NumberFormat.getNumberInstance();
        decimal4 = (DecimalFormat)NumberFormat.getNumberInstance();
        decimal2.applyPattern("#0.00");
        decimal3.applyPattern("#0.000");
        decimal4.applyPattern("#0.0000");
        tmpCurrentIndex = 0;
        tmpIsExitCondition = false;
        tmpKeyFileString = keyFileStringOrigin 
                + "DIELECTRIC\t" 
                + dielectricConstant
                + LINESEPARATOR;
        tmpForcefield = forcefield_IE;
        System.out.println("Calculating intermolecular energy...");
        
        //</editor-fold>
        
        while (tmpIsExitCondition == false) {
            // Exit condition is true when all particle pair combinations
            //   were calculated.
            if (tmpJobTaskRecordList.get(tmpCurrentIndex).hasEnergieJob()) {
                //<editor-fold defaultstate="collapsed" desc="Job task record">
                tmpParticleName1 = tmpJobTaskRecordList.get(tmpCurrentIndex)
                        .particleName1();
                tmpParticleName2 = tmpJobTaskRecordList.get(tmpCurrentIndex)
                        .particleName2();
                tmpParticlePair = tmpParticleName1 + "_" + tmpParticleName2;
                tmpIsSameParticle = tmpParticleName1.equals(tmpParticleName2);
                tmpIEResultDirName = tmpJobTaskRecordList.get(tmpCurrentIndex)
                        .result_IE_PathName();
                tmpH2OPos = 0;
                if (tmpParticleName1.equals("H2O") 
                        ^ tmpParticleName2.equals("H2O")) {
                    if (tmpParticleName1.equals("H2O")) {
                        tmpH2OPos = 1;
                    } else {
                        tmpH2OPos = 2;
                    }
                }
            
                //</editor-fold>
            
                //<editor-fold defaultstate="collapsed" desc="Create log file">
                tmpParticleDatFileName = tmpIEResultDirName
                        + FILESEPARATOR
                        + tmpParticlePair
                        + ".dat";
                tmpParticleDatFile = Paths.get(tmpParticleDatFileName);
                try {
                    if (!Files.exists(tmpParticleDatFile)) {
                        Files.createFile(tmpParticleDatFile);
                    }
                    BWParticleDat = new BufferedWriter(
                            new FileWriter(tmpParticleDatFileName));
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during writing in log file.", ex);
                }
                
                //</editor-fold>

                //<editor-fold defaultstate="collapsed" desc="Centre first fragment">
                // Centre the fragments and move to the centre
                tmpXyz1ID = particleNames.indexOf(tmpParticleName1);
                tmpTinkerXYZ1 = new TinkerXYZ(tmpForcefield,
                        tmpParticleName1,
                        xyzContent1[tmpXyz1ID]);
                tmpXyzData1 = tmpTinkerXYZ1.getCoordinateList1()[0];
                tmpCentre1 = tmpTinkerXYZ1.findCentreCoordinate();
                tmpXyzData1 = tmpTinkerXYZ1.moveCoordinates(tmpXyzData1, 
                        tmpCentre1);
                if (tmpIsSameParticle) {
                    tmpTinkerXYZ2  = new TinkerXYZ(tmpForcefield,
                        tmpParticleName2,
                        xyzContent1[tmpXyz1ID]);
                } else {
                    tmpXyz2ID = particleNames.indexOf(tmpParticleName2);
                    tmpTinkerXYZ2 = new TinkerXYZ(tmpForcefield,
                            tmpParticleName2,
                            xyzContent2[tmpXyz2ID]);
                }
                tmpXyzData2 = tmpTinkerXYZ2.getCoordinateList1()[0];
                tmpCentre2 = tmpTinkerXYZ2.findCentreCoordinate();
                tmpXyzData2 = tmpTinkerXYZ2
                        .moveCoordinates(tmpXyzData2, tmpCentre2);

                //</editor-fold>
                
                //<editor-fold defaultstate="collapsed" desc="Calculate Intermolecular Energy">

                //<editor-fold defaultstate="collapsed" desc="Calculate rotated coordinates">
                tmpRotCoords = RotationUtil
                        .getRotationsCoords(sphereNodeNumber1, 
                                rotationNumber1,
                                tmpXyzData1, 
                                tmpXyzData2,
                                isFibonacciSphereAlgorithm);
                tmpXyzRotData1 = tmpRotCoords.get(0);
                tmpXyzRotData2 = tmpRotCoords.get(1);
                
                //</editor-fold>

                //<editor-fold defaultstate="collapsed" desc="Prescan">
                tmpEnergyCalcTime = System.currentTimeMillis();
                tmpDistSize = (int)Math.ceil((upperBoundary - lowerBoundary) 
                        / prescanStepSize) + 1;
                tmpAllDistances = new LinkedList<>();
                tmpDistanceList = new LinkedList<>();
                tmpDistWgtEmins = new LinkedList<>();
                tmpDistEminRecords = new LinkedList<>();
                tmpPrmID1 = particleNames.indexOf(tmpParticleName1);
                if (tmpIsSameParticle) {
                    tmpPrmID2 = tmpPrmID1;
                } else {
                    tmpPrmID2 = particleNames.indexOf(tmpParticleName2);
                }
                
                // Write .key file
                tmpKeyContent = keyFileStringOrigin 
                        + "DIELECTRIC\t" 
                        + dielectricConstant 
                        + LINESEPARATOR
                        + "NONBONDTERM ONLY"
                        + LINESEPARATOR
                        + "OPENMP-THREADS 1"
                        + LINESEPARATOR;
                tmpKeyPathName = scratchDirectory
                        + FILESEPARATOR 
                        + tmpParticlePair 
                        + ".key";
                if (tmpForcefield.equals("OPLSAALIGPARGEN")) {
                    tmpKeyContent += prmContent1[tmpPrmID1];
                    if (!tmpIsSameParticle) {
                        tmpKeyContent += prmContent2[tmpPrmID2];
                    }
                }
                MIPETUTIL.writeKeyFile(tmpKeyPathName, tmpKeyContent);

                for (int i = 0; i < tmpDistSize; i++) {
                    tmpDistanceCandidate = lowerBoundary + i * prescanStepSize;
                    tmpDistanceList.add(tmpDistanceCandidate);
                    if (tmpDistanceCandidate > upperBoundary) {
                        tmpDistanceList.add(tmpDistanceCandidate);
                        break;
                    }
                }

                tmpAllDistances.addAll(tmpDistanceList);
                double[] tmpDistances = new double[tmpDistanceList.size()];

                for (int i = 0; i < tmpDistances.length; i++) {
                    tmpDistances[i] = tmpDistanceList.get(i);
                }

                tmpEnergyRecords = new EnergyRecord[4];
                tmpEnergyRecords[0] = getInterMolecularEnergy(
                        tmpParticlePair,
                        tmpDistances, 
                        tmpTinkerXYZ1, 
                        tmpTinkerXYZ2, 
                        tmpXyzRotData1, 
                        tmpXyzRotData2,
                        1E10);
                tmpEqDist = tmpEnergyRecords[0].eqDistance();
                tmpGlbEmin = tmpEnergyRecords[0].Emin();
                tmpGlbWgtEmin = tmpEnergyRecords[0].wgtEmin();
                
                for (int i = 0; i < tmpDistances.length; i++) {
                    tmpDistEminRecords.add(new Distance_EnergyRecord(
                            tmpEnergyRecords[0].distances()[i], 
                            tmpEnergyRecords[0].Emins()[i],
                            tmpEnergyRecords[0].wgtEmins()[i]));
                }
                
                //</editor-fold>
                
                //<editor-fold defaultstate="collapsed" desc="Calculate rotated coordinates">
                tmpRotCoords = RotationUtil
                        .getRotationsCoords(sphereNodeNumber2, 
                                rotationNumber2,
                                tmpXyzData1,
                                tmpXyzData2,
                                isFibonacciSphereAlgorithm);
                tmpXyzRotData1 = tmpRotCoords.get(0);
                tmpXyzRotData2 = tmpRotCoords.get(1);
                
                //</editor-fold>

                //<editor-fold defaultstate="collapsed" desc="Precise scan">
                tmpDistanceList = new LinkedList<>();
                tmpDistSize = 9;

                for (int i = 0; i < tmpDistSize; i++) {
                    tmpDistanceCandidate = (10 * tmpEqDist - 4 + i) / 10;
                    tmpDistanceList.add(tmpDistanceCandidate);
                }

                tmpAllDistances.addAll(tmpDistanceList);
                tmpDistances = new double[tmpDistanceList.size()];

                for (int i = 0; i < tmpDistances.length; i++) {
                    tmpDistances[i] = tmpDistanceList.get(i);
                }

                tmpEnergyRecords[1] = getInterMolecularEnergy(tmpParticlePair,
                        tmpDistances, 
                        tmpTinkerXYZ1, 
                        tmpTinkerXYZ2, 
                        tmpXyzRotData1, 
                        tmpXyzRotData2,
                        tmpGlbEmin);
                if (tmpEnergyRecords[1].Emin() < tmpGlbEmin) {
                    tmpEqDist = tmpEnergyRecords[1].eqDistance();
                    tmpGlbEmin = tmpEnergyRecords[1].Emin();
                }
                if (tmpEnergyRecords[1].wgtEmin() < tmpGlbWgtEmin) {
                    tmpGlbWgtEmin = tmpEnergyRecords[1].wgtEmin();
                }
                
                for (int i = 0; i < tmpDistances.length; i++) {
                    tmpDistEminRecords.add(new Distance_EnergyRecord(
                            tmpEnergyRecords[1].distances()[i], 
                            tmpEnergyRecords[1].Emins()[i],
                            tmpEnergyRecords[1].wgtEmins()[i]));
                }
                
                //</editor-fold>
                
                //<editor-fold defaultstate="collapsed" desc="Calculate rotated coordinates">
                tmpRotCoords = RotationUtil
                        .getRotationsCoords(sphereNodeNumber3, 
                                rotationNumber3,
                                tmpXyzData1,
                                tmpXyzData2,
                                isFibonacciSphereAlgorithm);
                tmpXyzRotData1 = tmpRotCoords.get(0);
                tmpXyzRotData2 = tmpRotCoords.get(1);
                
                //</editor-fold>

                //<editor-fold defaultstate="collapsed" desc="More precise scan">
                tmpDistanceList = new LinkedList<>();
                tmpDistSize = 19;

                for (int i = 0; i < tmpDistSize; i++) {
                    tmpDistanceCandidate = (100 * tmpEqDist - 9 + i) / 100;
                    tmpDistanceList.add(tmpDistanceCandidate);
                }

                tmpAllDistances.addAll(tmpDistanceList);
                tmpDistances = new double[tmpDistanceList.size()];

                for (int i = 0; i < tmpDistances.length; i++) {
                    tmpDistances[i] = tmpDistanceList.get(i);
                }

                tmpEnergyRecords[2] = getInterMolecularEnergy(tmpParticlePair,
                        tmpDistances, 
                        tmpTinkerXYZ1, 
                        tmpTinkerXYZ2, 
                        tmpXyzRotData1, 
                        tmpXyzRotData2,
                        tmpGlbEmin);
                if (tmpEnergyRecords[2].Emin() < tmpGlbEmin) {
                    tmpEqDist = tmpEnergyRecords[2].eqDistance();
                    tmpGlbEmin = tmpEnergyRecords[2].Emin();
                }
                if (tmpEnergyRecords[2].wgtEmin() < tmpGlbWgtEmin) {
                    tmpGlbWgtEmin = tmpEnergyRecords[2].wgtEmin();
                }
                
                for (int i = 0; i < tmpDistances.length; i++) {
                    tmpDistEminRecords.add(new Distance_EnergyRecord(
                            tmpEnergyRecords[2].distances()[i], 
                            tmpEnergyRecords[2].Emins()[i],
                            tmpEnergyRecords[2].wgtEmins()[i]));
                }
                
                //</editor-fold>
                
                //<editor-fold defaultstate="collapsed" desc="Calculate rotated coordinates">
                tmpRotCoords = RotationUtil
                        .getRotationsCoords(sphereNodeNumber4, 
                                rotationNumber4,
                                tmpXyzData1,
                                tmpXyzData2,
                                isFibonacciSphereAlgorithm);
                
                //</editor-fold>
                
                //<editor-fold defaultstate="collapsed" desc="Last scan">
                tmpAllDistances.add(tmpEnergyRecords[2].eqDistance());
                
                tmpDistances = new double[1];
                tmpDistances[0] = tmpEqDist; 
                tmpEnergyRecords[3] = getInterMolecularEnergy(tmpParticlePair,
                        tmpDistances, 
                        tmpTinkerXYZ1, 
                        tmpTinkerXYZ2, 
                        tmpXyzRotData1, 
                        tmpXyzRotData2,
                        tmpGlbEmin);
                if (tmpEnergyRecords[3].Emin() < tmpGlbEmin) {
                    tmpEqDist = tmpEnergyRecords[3].eqDistance();
                    tmpGlbEmin = tmpEnergyRecords[3].Emin();
                }
                if (tmpEnergyRecords[3].wgtEmin() < tmpGlbWgtEmin) {
                    tmpGlbWgtEmin = tmpEnergyRecords[3].wgtEmin();
                }
                
                for (int i = 0; i < tmpDistances.length; i++) {
                    tmpDistEminRecords.add(new Distance_EnergyRecord(
                            tmpEnergyRecords[3].distances()[i], 
                            tmpEnergyRecords[3].Emins()[i],
                            tmpEnergyRecords[3].wgtEmins()[i]));
                }
                
                //</editor-fold>
                
                //<editor-fold defaultstate="collapsed" desc="Sort datas">
                Double[] tmpDistanceObj;
                double[][] tmpEnergySorted;
                Integer[] tmpDistanceIndices;
                
                tmpDistSize = tmpAllDistances.size();
                tmpDistanceObj = tmpAllDistances.toArray(Double[]::new);
                tmpComparator = new ArrayIndexComparator(tmpDistanceObj);
                tmpDistanceIndices = tmpComparator.createIndexArray();
                Arrays.sort(tmpDistanceIndices, tmpComparator);
                tmpEnergySorted = new double[tmpDistSize][3];

                for (int i = 0; i < tmpDistSize; i++) {
                    tmpEnergySorted[i][0] = tmpDistEminRecords
                            .get(tmpDistanceIndices[i]).distance();
                    tmpEnergySorted[i][1] = tmpDistEminRecords
                            .get(tmpDistanceIndices[i]).Emin();
                    tmpEnergySorted[i][2] = tmpDistEminRecords
                            .get(tmpDistanceIndices[i]).wgtEmin();
                }

                //</editor-fold>

                //<editor-fold defaultstate="collapsed" desc="Determining opt. Emin">
                // Delete old .key file and make new one
                tmpKeyFile = Paths.get(scratchDirectory 
                        + FILESEPARATOR
                        + tmpParticlePair
                        + ".key");
                try {
                    Files.deleteIfExists(tmpKeyFile);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE,
                            "IOException during deleting .key file in scratch.",
                            ex);
                }
                
                // Determining opt. Emin
                double tmpOptMinEnergy;
                double tmpRgdMinEnergy;
                File tmpOptFile;
                
                tmpOptMinEnergy = 0.0;
                tmpRgdMinEnergy = 0.0;
                tmpKeyFileName = tmpParticlePair + ".key";
                tmpKeyPathName = scratchDirectory 
                        + FILESEPARATOR 
                        + tmpKeyFileName;
                tmpKeyContent = tmpKeyFileString;
                if (tmpForcefield.equals("OPLSAALIGPARGEN")) {
                    tmpKeyContent += prmContent1[tmpPrmID1];
                    if (!tmpIsSameParticle) {
                        tmpKeyContent += prmContent2[tmpPrmID2];
                    }
                }
                MIPETUTIL.writeKeyFile(tmpKeyPathName, tmpKeyContent);
                
                for (int i = 0; i < 2; i++) {
                    if (i == 0) {
                        tmpCmdList = new String[] {tinkerOptimize, 
                            scratchDirectory 
                            + FILESEPARATOR 
                            + tmpParticlePair
                            + ".0",        
                            Double.toString(optimizeRmsGradient)};
                    } else {
                        tmpCmdList = new String[] {tinkerOptrigid, 
                            scratchDirectory 
                            + FILESEPARATOR 
                            + tmpParticlePair
                            + ".0",
                            Double.toString(optimizeRmsGradient)};
                    }
                    
                    // Start optimize
                    tmpAtomSize1 =tmpTinkerXYZ1.getAtomSize1();
                    tmpProcess = null;
                    tmpPB = new ProcessBuilder();
                    tmpPB.redirectErrorStream(true);
                    tmpPB.command(tmpCmdList);
                    try {
                        tmpProcess = tmpPB.start();
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during process starting.",
                                ex);
                    }

                    // This is necessary because .waitFor() will hang otherwise
                    if (tmpProcess != null) {
                        try (BufferedReader tmpBR = new BufferedReader(
                                new InputStreamReader(tmpProcess
                                        .getInputStream()))) {
                            while (tmpBR.readLine() != null ) {}
                        } catch (IOException ex) {
                            LOGGER.log(Level.SEVERE, 
                                    "IOException during writing .0 file in scratch.",
                                    ex);
                        }
                        try {
                            tmpProcess.waitFor();
                        } catch (InterruptedException ex) {
                            LOGGER.log(Level.SEVERE, 
                                    "InterruptException during processing optimize.exe",
                                    ex);
                        }
                        tmpProcess.destroy();
                    }
                    
                    // Fix .xyz file if there is H2O 
                    //  this is necessary because of a bug in tinker's optimize
                    if (forcefield_IE.equals("OPLSAALIGPARGEN") 
                            && tmpH2OPos > 0) {
                        tmpFileName = scratchDirectory
                            + FILESEPARATOR
                            + tmpParticleName1 
                            + "_"
                            + tmpParticleName2
                            + ".xyz";
                        int tmpIndex;
                        if (tmpParticleName1.equals("H2O")) {
                            tmpIndex = 1;
                        } else {
                            tmpIndex = tmpAtomSize1 + 1;
                        }
                        MIPETUTIL.fixTinkerXYZ_H2O(tmpFileName, tmpIndex);
                    }
                    
                    // Use tinker's analyze to determine intermolecular energy
                    tmpProcess = null;
                    if (i == 0)  {
                        tmpOutputName = scratchDirectory
                            + FILESEPARATOR
                            + "output0_opt.txt";
                    } else {
                        tmpOutputName = scratchDirectory
                            + FILESEPARATOR
                            + "output0_rgd.txt";
                    }
                    tmpOptFile = new File(tmpOutputName);
                    tmpCmdList = new String[] {tinkerAnalyze, 
                        scratchDirectory
                        + FILESEPARATOR 
                        + tmpParticlePair
                        + ".xyz",
                        "E"};
                    tmpPB = new ProcessBuilder();
                    tmpPB.redirectErrorStream(true);
                    tmpPB.command(tmpCmdList);
                    tmpPB.redirectOutput(tmpOptFile);
                    try {
                        tmpProcess = tmpPB.start();
                        tmpProcess.waitFor();
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during process start.",
                                ex);
                    } catch (InterruptedException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "InterruptException during reading output0_opt.txt",
                                ex);
                    }
                    tmpProcess.destroy();
                    tmpSourceName = scratchDirectory
                            + FILESEPARATOR 
                            + tmpParticlePair
                            + ".xyz";
                    if (i == 0) {
                        tmpTargetName = scratchDirectory
                                + FILESEPARATOR 
                                + tmpParticlePair
                                + "_opt"
                                + ".xyz";
                    } else {
                        tmpTargetName = scratchDirectory
                                + FILESEPARATOR 
                                + tmpParticlePair
                                + "_rgd"
                                + ".xyz";
                    }
                    try {
                        Files.move(Paths.get(tmpSourceName), 
                                Paths.get(tmpTargetName), 
                                StandardCopyOption.REPLACE_EXISTING);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException renaming .xyz file.", ex);
                    }
                    
                   // Read the intermolecular energies from .txt files
                    String tmpSearch = "Intermolecular Energy :";
                    Path tmpPath = Paths.get(tmpOutputName);

                    try (BufferedReader tmpBR = Files.newBufferedReader(
                            tmpPath, StandardCharsets.UTF_8)) {

                        while ((tmpLine = tmpBR.readLine()) != null ) {
                            if (tmpLine.contains(tmpSearch)) {
                                if (i == 0) {
                                    tmpOptMinEnergy = Double.parseDouble(tmpLine
                                        .substring(25, 50));
                                } else {
                                    tmpRgdMinEnergy = Double.parseDouble(tmpLine
                                        .substring(25, 50));
                                }
                                break;
                            }
                        }

                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during reading output0_opt.txt.",
                                ex);
                    }
                }

                //</editor-fold>
                
                energyList.add(new ResultEnergyRecord(
                        tmpParticleName1, 
                        tmpParticleName2, 
                        tmpGlbWgtEmin,
                        tmpOptMinEnergy,
                        tmpRgdMinEnergy,
                        tmpGlbEmin));

                //</editor-fold>
                
                //<editor-fold defaultstate="collapsed" desc="Write dist vs. energy datas">
                tmpOutputName = tmpJobTaskRecordList.get(tmpCurrentIndex)
                        .result_IE_PathName()
                        + FILESEPARATOR
                        + tmpParticlePair + "_dist_vs_energy.dat";
                try (BufferedWriter tmpBW = new BufferedWriter(
                        new FileWriter(tmpOutputName))) {
                    tmpBW.append("distance [" + ANGSTROM 
                            +"]  Emin(Cmin,r) [kcal/mole] <E>(r) [kcal/mole]" 
                            + LINESEPARATOR);
                    
                    for (int i = 0; i < tmpDistSize; i++) {
                        tmpBW.append(MIPETUTIL.padLeft(decimal2.format(
                                tmpEnergySorted[i][0]), 8));
                        tmpBW.append(MIPETUTIL.padLeft(decimal3.format( 
                                tmpEnergySorted[i][1]), 20));
                        tmpBW.append(MIPETUTIL.padLeft(decimal3.format( 
                                tmpEnergySorted[i][2]), 20));
                        tmpBW.append(LINESEPARATOR);
                    }
                    
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during writing output0.out in scratch."
                            , ex);
                }
                
                //</editor-fold>
                
                //<editor-fold defaultstate="collapsed" desc="Make dist vs. energy diagram">
                String tmpEnergyDataPathName = tmpJobTaskRecordList
                        .get(tmpCurrentIndex).result_IE_PathName()
                        + FILESEPARATOR
                        + tmpParticlePair + "_dist_vs_energy.dat";
                String tmpEnergyGraphicsPrefix = tmpJobTaskRecordList
                        .get(tmpCurrentIndex).particleName1()
                        + "_"
                        + tmpJobTaskRecordList
                                .get(tmpCurrentIndex).particleName2()
                        + "_";
                ChartUtil tmpChartUtil = new ChartUtil();
                boolean tmpIsOperationSuccessful = 
                        tmpChartUtil.createEnergyGraphics(
                                tmpEnergyDataPathName, 
                                tmpEnergyGraphicsPrefix);
                
                //</editor-fold>

                //<editor-fold defaultstate="collapsed" desc="Write log file">
                try {
                    BFGblLog.append(tmpParticlePair);
                    BFGblLog.append(LINESEPARATOR);
                    BWParticleDat.append("Force field: " + tmpForcefield);
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("Conformational analysis: " + isConformationalAnalysis);
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("equilibriumDistances [" + ANGSTROM + "] = "); 
                    BWParticleDat.append(decimal2.format(tmpEqDist));
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.flush();
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during writing log file.", ex);
                }

                //</editor-fold>

                //<editor-fold defaultstate="collapsed" desc="Copy results">
                boolean tmpHasH2O;
                Path tmpOriginal;
                Path tmpTarget;
                
                // Write ouput.0 file
                tmpOriginal = Paths.get(scratchDirectory, tmpParticlePair + ".0");
                tmpTarget = Paths.get(tmpIEResultDirName,"output.0");
                try {
                    Files.copy(tmpOriginal, tmpTarget,
                        StandardCopyOption.REPLACE_EXISTING);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during copying output.0", ex);
                }
                String tmpOutput0 = "Intermolecular Energy: " 
                        + tmpEnergyRecords[3].wgtEmin() + " kcal/mol";
                String tmpTargetDir = tmpIEResultDirName 
                        + FILESEPARATOR 
                        + "output0.out";        
                try (BufferedWriter tmpBW = Files.newBufferedWriter(
                        Paths.get(tmpTargetDir), StandardCharsets.UTF_8)) {
                    tmpBW.append(tmpOutput0);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during writing output0.out in scratch."
                            , ex);
                }
                
                // Write output_opt.out
                tmpOriginal = Paths.get(scratchDirectory, 
                        tmpParticlePair + "_opt.xyz");
                tmpTarget = Paths.get(tmpIEResultDirName,"output_opt.0");
                try {
                    Files.copy(tmpOriginal, tmpTarget,
                        StandardCopyOption.REPLACE_EXISTING);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during copying output.0", ex);
                }
                tmpOutput0 = "Intermolecular Energy: " 
                        + tmpOptMinEnergy + " kcal/mol";
                tmpTargetDir = tmpIEResultDirName 
                        + FILESEPARATOR 
                        + "output0_opt.out";        
                try (BufferedWriter tmpBW = Files.newBufferedWriter(
                        Paths.get(tmpTargetDir), StandardCharsets.UTF_8)) {
                    tmpBW.append(tmpOutput0);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during writing output0_opt.out "
                            + "in scratch.", ex);
                }

                // Write output_rgd.out
                tmpOriginal = Paths.get(scratchDirectory, 
                        tmpParticlePair + "_rgd.xyz");
                tmpTarget = Paths.get(tmpIEResultDirName,"output_rgd.0");
                try {
                    Files.copy(tmpOriginal, tmpTarget,
                        StandardCopyOption.REPLACE_EXISTING);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during copying output.0", ex);
                }
                tmpOutput0 = "Intermolecular Energy: " 
                        + tmpRgdMinEnergy + " kcal/mol";
                tmpTargetDir = tmpIEResultDirName 
                        + FILESEPARATOR 
                        + "output0_rgd.out";        
                try (BufferedWriter tmpBW = Files.newBufferedWriter(
                        Paths.get(tmpTargetDir), StandardCharsets.UTF_8)) {
                    tmpBW.append(tmpOutput0);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during writing output0_rgd.out "
                            + "in scratch.", ex);
                }
                
                // Generate output.xyz
                String tmpOutput0FileName = tmpIEResultDirName 
                        + FILESEPARATOR 
                        + "output.0";
                String tmpOptFileName;
                String tmpRgdFileName;
                int tmpAtomNumber1 = tmpTinkerXYZ1.getAtomSize1();
                int tmpAtomNumber2 = tmpTinkerXYZ2.getAtomSize1();
                tmpElementList1 = tmpTinkerXYZ1.getElementList1();
                tmpElementList2 = tmpTinkerXYZ2.getElementList1();
                TinkerXYZ tmpTinkerXyz = new TinkerXYZ(tmpOutput0FileName, 1, 
                        tmpAtomNumber1, tmpAtomNumber2);
                tmpFileName = tmpIEResultDirName 
                        + FILESEPARATOR 
                        + "output.xyz";
                tmpTinkerXyz.writeToXyzFile(tmpFileName);
                
                // Generate .pdb file of output.0
                tmpHasH2O = tmpParticleName1.equals("H2O") || 
                        tmpParticleName2.equals("H2O");
                String tmpKeyFileContent = "\"" 
                        + scratchDirectory
                        + FILESEPARATOR 
                        + tmpParticlePair
                        + ".key"
                        + "\"";
                MIPETUTIL.callXYZPDB(tinkerXYZPdb, 
                        tmpOutput0FileName, 
                        tmpKeyFileContent, 
                        tmpHasH2O);
                
                // Generate output_opt.xyz
                tmpOptFileName = tmpIEResultDirName 
                        + FILESEPARATOR 
                        + "output_opt.0";
                tmpTinkerXyz = new TinkerXYZ(tmpOptFileName, 1, tmpAtomNumber1, 
                        tmpAtomNumber2);
                tmpTinkerXyz.setElementList1(tmpElementList1);
                tmpTinkerXyz.setElementList2(tmpElementList2);
                tmpFileName = tmpIEResultDirName 
                        + FILESEPARATOR 
                        + "output_opt.xyz";
                tmpTinkerXyz.writeToXyzFile(tmpFileName);

                // Generate output_rgd.xyz
                tmpRgdFileName = tmpIEResultDirName 
                        + FILESEPARATOR 
                        + "output_rgd.0";
                tmpTinkerXyz = new TinkerXYZ(tmpRgdFileName, 1, tmpAtomNumber1, 
                        tmpAtomNumber2);
                tmpTinkerXyz.setElementList1(tmpElementList1);
                tmpTinkerXyz.setElementList2(tmpElementList2);
                tmpFileName = tmpIEResultDirName 
                        + FILESEPARATOR 
                        + "output_rgd.xyz";
                tmpTinkerXyz.writeToXyzFile(tmpFileName);

                // Generate .pdb file of output_opt.0
                MIPETUTIL.callXYZPDB(tinkerXYZPdb, 
                        tmpOptFileName, 
                        tmpKeyFileContent, 
                        tmpHasH2O);

                // Generate .pdb file of output_rgd.0
                MIPETUTIL.callXYZPDB(tinkerXYZPdb, 
                        tmpRgdFileName, 
                        tmpKeyFileContent, 
                        tmpHasH2O);
                
                // Further outputs
                Path tmpOptDistDir = Paths.get(optDistDirectory 
                        + FILESEPARATOR 
                        + tmpForcefield);
                if(!Files.exists(tmpOptDistDir)) {
                    try {
                        Files.createDirectories(tmpOptDistDir);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                            "IOException during creating OutDist directory.",
                            ex);
                    }
                }
                tmpFileName = optDistDirectory
                        + FILESEPARATOR
                        + tmpForcefield
                        + FILESEPARATOR
                        + tmpParticleName1
                        + ".txt";
                tmpOptDistFile = Paths.get(tmpFileName);
                if (tmpIsSameParticle && !Files.exists(tmpOptDistFile)) {
                    try (BufferedWriter tmpBW = new BufferedWriter(
                            new FileWriter(tmpFileName))) {
                        tmpBW.append(decimal4.format(tmpEqDist));
                    } catch(IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                            "IOException during writing file in OptDist directory.", 
                            ex);
                    }
                }
                tmpEnergyCalcTime  = (System.currentTimeMillis() 
                        - tmpEnergyCalcTime) / 1000;
                try {
                    BWParticleDat.append("SphereNodeNumber1: "
                            + Integer.toString(sphereNodeNumber1)
                            + LINESEPARATOR);
                    BWParticleDat.append("SphereNodeNumber2: "
                            + Integer.toString(sphereNodeNumber2)
                            + LINESEPARATOR);
                    BWParticleDat.append("SphereNodeNumber3: "
                            + Integer.toString(sphereNodeNumber3)
                            + LINESEPARATOR);
                    BWParticleDat.append("SphereNodeNumber4: "
                            + Integer.toString(sphereNodeNumber4)
                            + LINESEPARATOR);
                    BWParticleDat.append("RotationNumber1: "
                            + Integer.toString(rotationNumber1)
                            + LINESEPARATOR);
                    BWParticleDat.append("RotationNumber2: "
                            + Integer.toString(rotationNumber2)
                            + LINESEPARATOR);
                    BWParticleDat.append("RotationNumber3: "
                            + Integer.toString(rotationNumber3)
                            + LINESEPARATOR);
                    BWParticleDat.append("RotationNumber4: "
                            + Integer.toString(rotationNumber4)
                            + LINESEPARATOR);
                    BWParticleDat.append("Temperature [K]: ");
                    BWParticleDat.append(Double.toString(temperature));
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("Fraction of energy values used for the Boltzmann distribution: ");
                    BWParticleDat.append(Double.toString(boltzmannFraction));
                    BWParticleDat.append(LINESEPARATOR);
                    if (isOptEmin) {
                        BWParticleDat.append("Optimize sampled E(min) configuration: "
                                + isOptEmin);
                        BWParticleDat.append(LINESEPARATOR);
                    }
                    BWParticleDat.append("Weighted (Emin = glbMin) MinimumIntermolecularEnergy [kcal/mole]: ");
                    BWParticleDat.append(decimal4.format(tmpGlbWgtEmin));
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("""
                                         Weighted (Emin = glbMin) MinimumIntermolecularEnergy: 
                                         Weighted differential pair interaction energy with
                                          Emin = lowest differential pair interaction energy.""");
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("GlobalMinimumIntermolecularEnergy [kcal/mole]: ");
                    BWParticleDat.append(decimal4.format(tmpGlbEmin));
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("""
                                         GlobalMinimumIntermolecularEnergy:
                                           Lowest differential pair interaction energy of all dimer configurations.""");
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("Optimized minimumIntermolecularEnergy [kcal/mole]: ");
                    BWParticleDat.append(decimal4.format(tmpOptMinEnergy));
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("""
                                         Optimized minimumIntermolecularEnergy: 
                                           Differential pair interaction energy from the dimer configuration
                                           with lowest differential pair interaction energy after optimize.""");
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("Rigid-optimized minimumIntermolecularEnergy [kcal/mole]: ");
                    BWParticleDat.append(decimal4.format(tmpRgdMinEnergy));
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("""
                                         Rigid-optimized minimumIntermolecularEnergy: 
                                           Differential pair interaction energy from the dimer configuration
                                           with lowest differential pair interaction energy after optrigid.""");
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("Time to calculate minimum intermolecular energy [s]: "
                            + tmpEnergyCalcTime
                            + LINESEPARATOR);
                    BWParticleDat.close();
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during writing log file.", ex);
                }
                
                //</editor-fold>
            
                //<editor-fold defaultstate="collapsed" desc="Clean scratch directory">
                try (Stream<Path> tmpWalk = Files.walk(Paths
                        .get(scratchDirectory))) {
                    tmpWalk.sorted(Comparator.reverseOrder())
                            .filter(Files::isRegularFile)
                            .map(Path::toFile)
                            .forEach(File::delete);
                } catch(IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during deleting scratch directory.", ex);
                }
            
                //</editor-fold>
            }
            tmpCurrentIndex++;
            if (tmpCurrentIndex >= tmpJobTaskRecordList.size()) {
                tmpIsExitCondition = true;
            }
        }
        
        //<editor-fold defaultstate="collapsed" desc="Write Readme.txt">
        String tmpReadme =
                """
                *** Notice for files in subfolders ***
                *************************************
                Tinker's tools use tinker's .xyz coordinate file. This file format
                is different than the original .xyz file format. It has additional
                information like atom numbers, atom type number, information about 
                connection to other atoms. Unfortunately Tinker's tools use same 
                extension .xyz.

                1particle_2particle.dat:
                Summary information like equilibrium distance, sphere node number,
                circle node number, temperature, weighted and unweighted minimum
                differential pair interaction energy.

                1particle_2particle_dist_vs_energy.dat
                Distance vs minimum differential pair interaction energy datas.

                1particle_2particle_dist_vs_energy.svg
                Diagram plot of distance vs minimum differential pair interaction 
                energy datas.

                1particle_2particle_log.txt
                Some informations for developer.

                1particle_2particle.out
                Tinker's .xyz coordinate file of the dimer configuration with lowest
                differential pair interaction energy.

                1particle_2particle.pdb
                Coordinate file of the dimer configuration with lowest differential 
                pair interaction energy as .pdb file.

                1particle_2particle.xyz
                (Original) .xyz coordinate file of the dimer configuration with 
                lowest differential pair interaction energy.

                1particle_2particle_opt.0
                1particle_2particle_opt.pdb
                1particle_2particle_opt.xyz
                Coordinate files from the dimer configuration with lowest 
                differential pair interaction energy after optimize.

                1particle_2particle_rgd.0
                1particle_2particle_rgd.pdb
                1particle_2particle_rgd.xyz
                Coordinate files from the dimer configuration with lowest 
                differential pair interaction energy after optrigid.

                output0.out
                output0_opt.out
                output0_rgd.out
                Output files from tinker's analyze with differential pair interaction 
                energy data. 
                """;
        String ReadmeFileName = resultDirectory 
                + FILESEPARATOR
                + "IE"
                + FILESEPARATOR
                + forcefield_IE
                + FILESEPARATOR
                + "Readme.txt";
        try (BufferedWriter tmpBW = new BufferedWriter(
                new FileWriter(ReadmeFileName))) {
            tmpBW.append(tmpReadme);
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during writing Readme.txt.", ex);
        }
                
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Coordination numbers">
        System.out.println("Calculating coordination numbers...");
        if (!forcefield_CN.isEmpty()) {
            getCoordinationNumbers(tmpJobTaskRecordList);
        }
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Export parameterset">
        System.out.println("Exporting parameterset...");
        
        float tmpTotalTimeSec;
        tmpTotalTimeSec = (System.currentTimeMillis() - tmpTotalTime) / 1000;
        try {
            BFGblLog.append(LINESEPARATOR);
            BFGblLog.append("Entire calculation Time: ");
            BFGblLog.append(tmpTotalTimeSec + " s ");
            if (tmpTotalTimeSec > 3600) {
                BFGblLog.append("(" + decimal2
                        .format(tmpTotalTimeSec / 3600) + " h)");
            } else if (tmpTotalTimeSec > 60) {
                BFGblLog.append("(" + decimal2
                        .format(tmpTotalTimeSec / 60) + " min)");
            }       
            BFGblLog.close();
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, "IOException during writing log file.", 
                    ex);
        }
        
        exportParticleSetForMFSim(tmpJobTaskRecordList, 
                parameterSetTitle, 
                parameterSetTitleAbr);
        System.out.println("Ready.");
        
        //</editor-fold>
       
    } 
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Initialize method
     */
    private static void initialize() {
        setParameters();
        
        // <editor-fold defaultstate="collapsed" desc="Check and creat directories">
        if (!Files.exists(Paths.get(scratchDirectory))) {
            try {
                Files.createDirectories(Paths.get(scratchDirectory));
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE, 
                            "IOException during creating scratch directory.",
                            ex);
            }
        }
        if (!new File(resultDirectory).isFile()) {
            new File(resultDirectory).mkdirs();
        }
        if (!new File(optXYZDirectory).isFile()) {
            new File(optXYZDirectory).mkdirs();
        }
        if (!new File(optDistDirectory).isFile()) {
            new File(optDistDirectory).mkdirs();
        }
        smiles = MIPETUTIL
                .getSmilesData(smilesDirectory + "/Smiles.dat");
        
        // </editor-fold>
    }
    
    /**
     * Read the Job file
     */
    private static void readJobFile() {
        String tmpLine;
        String tmpRestString;
        Character tmpFirstChar;
        
        particleNames = new LinkedList<>();
        newParticles = new LinkedList<>();
        oldParticles = new LinkedList<>();
        
        try (BufferedReader tmpBR = new BufferedReader(
                new FileReader(jobFileName))) {
            // read jobs
            while ((tmpLine = tmpBR.readLine()) != null ) {
                if (tmpLine.isEmpty()) {
                    continue;
                }
                tmpFirstChar = tmpLine.trim().charAt(0);
                tmpRestString = tmpLine.trim().substring(1).trim();
                switch (tmpFirstChar) {
                    case '#' -> {
                        continue; // ignore comment line
                    }
                    case '*' -> {
                        String[] tmpString = tmpRestString.split("\\s+");
                        forcefield_IE = tmpString[0];
                        continue;
                    }
                    case '$' -> {
                        String[] tmpString = tmpRestString.split("\\s+");
                        forcefield_CN = tmpString[0];
                        continue;
                    }
                    case '-' -> {
                        if (!oldParticles.contains(tmpRestString)) {
                            oldParticles.add(tmpRestString);
                        }
                        if (!particleNames.contains(tmpRestString)) {
                            particleNames.add(tmpRestString);
                        }
                        continue;
                    }
                    default -> {
                        tmpLine = tmpLine.trim();
                        if (!newParticles.contains(tmpLine)) {
                            newParticles.add(tmpLine);
                        }
                        if (!particleNames.contains(tmpLine)) {
                            particleNames.add(tmpLine);
                        }
                    }
                }
            }
            if (newParticles.isEmpty()) {
                LOGGER.log(Level.SEVERE, 
                    "No new particles in the job file listed.");
            }
            tmpBR.close();
        } catch (FileNotFoundException ex) {
            throw new IllegalArgumentException("No jobfile found.");
        } catch (IOException ex) {
            throw new IllegalArgumentException("IOException during"
                    + "reading job file.");
        }
    }
    
    /**
     * SetParameters method
     */
    private static void setParameters() {
        String tmpIsFibonacciSphereAlgorithm;
        String tmpIsOptEmin;
        
        isTinker9 = MIPETUTIL.getResourceString("MIPET.Tinker9").toLowerCase()
                .equals("true");
        isTinkerOn = MIPETUTIL.getResourceString("MIPET.TinkerOn")
                .toLowerCase().equals("true");
        jobFileName = MIPETUTIL.getResourceString("MIPET.File.jobfile");
        cpuCoreNumber = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETCPUCoreNumber"));
        tmpIsFibonacciSphereAlgorithm = MIPETUTIL.getResourceString(
                "MIPETJavaUseFibonacciSphereAlgorithm");
        isFibonacciSphereAlgorithm = tmpIsFibonacciSphereAlgorithm
                .equals("true");
        tmpIsOptEmin = MIPETUTIL.getResourceString("MIPETOptEmin");
        isOptEmin = tmpIsOptEmin.equals("true");
        isConformationalAnalysis = MIPETUTIL
                .getResourceString("MIPETConformationalAnalysis")
                .equals("true");
        temperature = Double.parseDouble(MIPETUTIL.getResourceString(
                "MIPETTemperature"));
        boltzmannFraction = Double.parseDouble(MIPETUTIL.getResourceString(
                "MIPETBoltzmannFraction"));
        sphereNodeNumber1 = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETSphereNodeNumber1"));
        sphereNodeNumber2 = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETSphereNodeNumber2"));
        sphereNodeNumber3 = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETSphereNodeNumber3"));
        sphereNodeNumber4 = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETSphereNodeNumber4"));
        rotationNumber1 = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETRotationNumber1"));
        rotationNumber2 = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETRotationNumber2"));
        rotationNumber3 = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETRotationNumber3"));
        rotationNumber4 = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETRotationNumber4"));
        minAtomDistance = Double.parseDouble(MIPETUTIL.getResourceString(
                "MIPETMinAtomDistance"));
        nSearchDirection = MIPETUTIL.getResourceString(
                "MIPETNumberSearchDirection");
        energyThreshold = MIPETUTIL.getResourceString(
                "MIPETEnergyThreshold");
        rmsGradient = MIPETUTIL.getResourceString(
                "MIPETRmsGradient");
        rmsMinimizeGradient = Double.parseDouble(MIPETUTIL.getResourceString(
                "MIPETBfgsRmsGradient"));
        minimizeMaxIteration = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETBfgsMaxIteration"));
        dielectricConstant = Double.parseDouble(MIPETUTIL.getResourceString(
                "MIPETDielectricConstant"));
        parameterDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.paramDirectory");
        tinkerDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.tinker");
        if (OS.matches("Windows [0-9]*")) {
            tinkerAnalyze = tinkerDirectory + "/analyze.exe";
            tinkerScan = tinkerDirectory + "/scan.exe";
            tinkerXYZEdit = tinkerDirectory + "/xyzedit.exe";
            tinkerMinimize = tinkerDirectory + "/minimize.exe";
            tinkerDynamic = tinkerDirectory + "/dynamic.exe";
            tinkerOptimize = tinkerDirectory + "/optimize.exe";
            tinkerOptrigid = tinkerDirectory + "/optrigid.exe";
            tinkerXYZPdb = tinkerDirectory + "/xyzpdb.exe";
        } else if (OS.equals("Linux")) {
            if (isTinker9) {
                tinkerDynamic = tinkerDirectory + "/tinker9";
            } else {
                tinkerDynamic = tinkerDirectory + "/dynamic";
            }
            tinkerAnalyze = tinkerDirectory + "/analyze";
            tinkerScan = tinkerDirectory + "/scan";
            tinkerXYZEdit = tinkerDirectory + "/xyzedit";
            tinkerMinimize = tinkerDirectory + "/minimize";
            tinkerOptimize = tinkerDirectory + "/optimize";
            tinkerOptrigid = tinkerDirectory + "/optrigid";
            tinkerXYZPdb = tinkerDirectory + "/xyzpdb";
        }
        scratchDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.scratch");
        smilesDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.smiles");
        calculationDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.calculation");
        resultDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.result");
        sourceDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.source");
        optXYZDirectory = calculationDirectory + "/OptXYZ";
        optDistDirectory = calculationDirectory + "/OptDist";
        scanProgram = MIPETUTIL.getResourceString(
                "MIPETScanProgram");
        lowerBoundary = Double.parseDouble(MIPETUTIL.getResourceString(
                "MIPETLowerBoundary"));
        upperBoundary = Double.parseDouble(MIPETUTIL.getResourceString(
                "MIPETUpperBoundary"));
        prescanStepSize = Double.parseDouble(MIPETUTIL.getResourceString(
                "MIPETPrescanStepSize"));
        solventMoleculeNumber = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETSolventMoleculeNumber"));
        stepNumber = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETStepNumber"));
        dynamicStepsPerCore = stepNumber / cpuCoreNumber;
        timeStep = Double.parseDouble(MIPETUTIL.getResourceString(
                "MIPETTimeStep"));
        printInterval = Double.parseDouble(MIPETUTIL.getResourceString(
                "MIPETPrintinterval"));
        simulationType = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETSimulationType"));
        catchRadius = Double.parseDouble(MIPETUTIL.getResourceString(
                "MIPETCatchRadius"));
        optimizeRmsGradient = Double.parseDouble(MIPETUTIL.getResourceString(
                "MIPETOptimizeRmsGradient"));
        warmUpStepNumber = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETWarmUpStepNumber"));
        warmUpTimeStep = Double.parseDouble(MIPETUTIL.getResourceString(
                "MIPETWarmUpTimeStep"));
        warmUpPrintInterval = warmUpStepNumber * 1E-3;
        // timeStep in fs and printIntervall in ps
        nDynamicWarmupIteration = 
                (int) (warmUpStepNumber * warmUpTimeStep * 1E-15 /
                (warmUpPrintInterval * 1E-12));
        // timeStep in fs and printIntervall in ps
        nDynamicIteration = (int)(stepNumber * timeStep * 1E-15 /
                (printInterval * 1E-12));
        // timeStep in fs and printIntervall in ps
        nDynamicIterationPerCore = 
                (int) Math.ceil((dynamicStepsPerCore * timeStep * 1E-15) /
                (printInterval * 1E-12));
        // Write warmup log file
        isLogWarmup = MIPETUTIL.getResourceString("MIPET.Log.Warmup")
                .toLowerCase().equals("true");
        // Write dynamic log file
        isLogDynamic = MIPETUTIL.getResourceString("MIPET.Log.Dynamic")
                .toLowerCase().equals("true");
        // Write the RATTLE key word in the .key file
        isRATTLE = MIPETUTIL.getResourceString("MIPETRATTLE")
                .toLowerCase().equals("true");
        // Write build simulation box log file
        isLogBuildBox = MIPETUTIL.getResourceString("MIPET.Log.Build")
                .toLowerCase().equals("true");
        // Write minimize simulation box log file
        isLogMinimizeBox = MIPETUTIL.getResourceString("MIPET.Log.MinimizeBox")
                .toLowerCase().equals("true");
        // Write soeak simulation box log file
        isLogSoakBox = MIPETUTIL.getResourceString("MIPET.Log.Soak")
                .toLowerCase().equals("true");
        // Parameterset title
        parameterSetTitle = MIPETUTIL.getResourceString("MIPETParticleSetTitle");
        // Parameterset title abbreviation
        parameterSetTitleAbr = MIPETUTIL.getResourceString("MIPETParticleSetTitleAbr");
    }
    
    /**
     * Method getParticlePairs
     * @return particle pair names
     */
    private static ArrayList<String> getParticlePairs() {
        
        int tmpNewParticleLength;
        int tmpOldParticleLength;
        ArrayList<String> tmpParticlePairs = new ArrayList<>();
        
        tmpNewParticleLength = newParticles.size();
        tmpOldParticleLength = oldParticles.size();
       
        // Same particles
        Collections.sort(newParticles);
        Collections.sort(oldParticles);
        
        for (int i = 0; i < tmpNewParticleLength; i++) {
            tmpParticlePairs.add(
                    newParticles.get(i) + "_" + newParticles.get(i));
        }
        
        // New with old particles
        for (int i = 0; i < tmpNewParticleLength; i++) {

            for (int j = 0; j < tmpOldParticleLength; j++) {
                tmpParticlePairs.add(
                        newParticles.get(i) +"_" + oldParticles.get(j));
            }
            
        }

        // Different particles
        for (int i = 0; i < tmpNewParticleLength; i++) {

            for (int j = 0; j < i; j++) {
                tmpParticlePairs.add(
                        newParticles.get(i) + "_" + newParticles.get(j));
            }

        }    
        
        return tmpParticlePairs;
    }
    
    /**
     * Method readXyz
     * Read .xyz files and store in xyzContent.
     * If conformational analysis is set xyzContent is overwritten with optimized data.
     */
    private static void readXyz(boolean anIsOriginal, String aForcefield) {
        int tmpParticlesLength;
        String tmpParticleName;
        String tmpXyzName1;
        String tmpOldAtomType;
        String tmpNewAtomType;
        String[] tmpLines;
        String[] tmpTokens;
        StringBuilder tmpXyz2;
        
        tmpParticlesLength = particleNames.size();
        xyzContent1 = new String[tmpParticlesLength];
        xyzContent2 = new String[tmpParticlesLength];
        tmpXyz2 = new StringBuilder(2000);
        
        for (int i = 0; i < tmpParticlesLength; i++) {
            tmpParticleName = particleNames.get(i);
            
            // Read .xyz files
            if (anIsOriginal) {
                tmpXyzName1 = sourceDirectory
                        + FILESEPARATOR
                        + aForcefield
                        + FILESEPARATOR
                        + tmpParticleName
                        + ".xyz";
            } else {
                tmpXyzName1 = optXYZDirectory
                        + FILESEPARATOR
                        + aForcefield
                        + FILESEPARATOR
                        + tmpParticleName
                        + FILESEPARATOR
                        + tmpParticleName
                        + ".xyz";
            }
            try {
                if (Files.exists(Paths.get(tmpXyzName1))) {
                    xyzContent1[i] = Files.readString(Paths
                            .get(tmpXyzName1));
                } 
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE, 
                        "IOException during reading .xyz file.", ex);
            }
            
            // Change atomtype number of 2. particle to avoid redundancy
            tmpXyz2.setLength(0);
            if (aForcefield.equals("OPLSAALIGPARGEN")) {
                tmpLines = xyzContent1[i].split(LINESEPARATOR);
                tmpXyz2.append(tmpLines[0]);

                for (int j = 1; j < tmpLines.length; j++) {
                    tmpTokens = tmpLines[j].trim().split("\\s+");
                    tmpOldAtomType = " " + tmpTokens[5] + " ";
                    tmpNewAtomType = " " + Integer.toString(Integer
                            .parseInt(tmpTokens[5]) + 100) + " ";
                    tmpLines[j] = tmpLines[j].replace(tmpOldAtomType, 
                            tmpNewAtomType);
                    tmpXyz2.append(LINESEPARATOR);
                    tmpXyz2.append(tmpLines[j]);
                }

                xyzContent2[i] = tmpXyz2.toString();
            } else {
                xyzContent2[i] = xyzContent1[i];
            }
        }
        
    }
    
    /** 
     * Method readPrm
     * Read .prm files for OPLSAALIGPARGEN and convert atomtype number of 
     *  second particle
     * @param aForcefield Force field name
     */
    private static void readPrm(String aForcefield) {
        int tmpParticlesLength;
        String tmpPrmName1;
        String tmpParticleName;
        String[] tmpLines;
        String[] tmpTokens;
        StringBuilder tmpPrm2;
        
        tmpParticlesLength = particleNames.size();
        prmContent1 = new String[tmpParticlesLength];
        prmContent2 = new String[tmpParticlesLength];
        tmpPrm2 = new StringBuilder(100000);
                
        for (int i = 0; i < tmpParticlesLength; i++) {
            tmpParticleName = particleNames.get(i);
            
            // Read .prm files
            if (aForcefield.equals("OPLSAALIGPARGEN")) {
                tmpPrmName1 = parameterDirectory
                    + FILESEPARATOR
                    + aForcefield
                    + FILESEPARATOR
                    + tmpParticleName
                    + ".prm";
                try {
                    if (Files.exists(Paths.get(tmpPrmName1))) {
                        prmContent1[i] = Files.readString(Paths
                                .get(tmpPrmName1));
                    } else {
                        prmContent1[i] = "";
                    }
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during reading .prm file.", ex);
                }
            }
            
            // Change atomtype number of 2. particle to avoid redundancy
            if (aForcefield.equals("OPLSAALIGPARGEN") && 
                    !prmContent1[i].equals("")) {
                tmpPrm2.setLength(0);
                tmpLines = prmContent1[i].split("\\n");
                
                for (int j = 25; j < tmpLines.length; j++) {
                    tmpTokens = tmpLines[j].trim().split("\\s+");
                    switch (tmpTokens[0]) {
                        case    "atom", 
                                "vdw", 
                                "charge",
                                "bond",
                                "angle",
                                "torsion",
                                "imptors"-> {
                            tmpLines[j] = MIPETUTIL.changeAtomType(tmpTokens);
                        }
                    }
                    if (!tmpLines[j].equals("")) {
                        tmpPrm2.append(tmpLines[j]);
                        tmpPrm2.append(LINESEPARATOR);
                    }
                }
                prmContent2[i] = tmpPrm2.toString();
            } else {
                prmContent2[i] = "";
            }
        }
    }


    /** 
     * Method makeMoleculeRecord
     * Make molecules record for intermolecular energy calculation without tinker
     *  (only for OPLSAALIGPARGEN)
     */
    private static void makeMoleculeRecord() {
        int tmpParticleSize;
        int tmpSigmaIndex;
        int tmpChargeIndex;
        int[] tmpAtomNumber;
        int[][] tmpAtomTypes;
        double[][] tmpEpsilsons;
        double[][] tmpSigmas;
        double[][] tmpCharges;
        String tmpLine;
        String[] tmpLines;
        String[][] tmpElements;
        
        molecules = new LinkedList<>();
        tmpParticleSize = particleNames.size();
        tmpAtomNumber = new int[tmpParticleSize];
        tmpElements = new String[tmpParticleSize][];
        tmpAtomTypes = new int[tmpParticleSize][];
        tmpSigmas = new double[tmpParticleSize][];
        tmpEpsilsons = new double[tmpParticleSize][];
        tmpCharges = new double[tmpParticleSize][];
        
        // Read elements and atomTypes
        for (int i = 0; i < tmpParticleSize; i++) {
            tmpLines = xyzContent1[i].split(LINESEPARATOR);
            tmpLine = tmpLines[0].substring(0, 6).trim();
            tmpAtomNumber[i] = Integer.parseInt(tmpLine);
            tmpElements[i] = new String[tmpAtomNumber[i]];
            tmpAtomTypes[i] = new int[tmpAtomNumber[i]];
            
            for (int j = 1; j < tmpLines.length; j++) {
                tmpLine = tmpLines[j].substring(7, 11).trim();
                tmpElements[i][j - 1] = tmpLine;
                tmpLine = tmpLines[j].substring(47, 51).trim();
                tmpAtomTypes[i][j-1] = Integer.parseInt(tmpLine);
            }
            
        }
        
        // Read epsilons, sigmas and charges
        for (int i = 0; i < tmpParticleSize; i++) {
            tmpSigmaIndex = 0;
            tmpChargeIndex = 0;
            tmpLines = prmContent1[i].split(LINESEPARATOR);
            tmpSigmas[i] = new double[tmpAtomNumber[i]];
            tmpEpsilsons[i] = new double[tmpAtomNumber[i]];
            tmpCharges[i] = new double[tmpAtomNumber[i]];
            
            for (int j = 0; j < tmpLines.length; j++) {
                if (tmpLines[j].startsWith("vdw ")) {
                    tmpLine = tmpLines[j].substring(20, 32).trim();
                    tmpSigmas[i][tmpSigmaIndex] = Double.parseDouble(tmpLine);
                    tmpLine = tmpLines[j].substring(32).trim();
                    tmpEpsilsons[i][tmpSigmaIndex] = Double
                            .parseDouble(tmpLine);
                    tmpSigmaIndex++;
                }
                if (tmpLines[j].startsWith("charge ")) {
                    tmpLine = tmpLines[j].substring(20);
                    tmpCharges[i][tmpChargeIndex] = Double.parseDouble(tmpLine);
                    tmpChargeIndex++;
                }
            }
            
            molecules.add(new MoleculeRecord(
                    particleNames.get(i),
                    tmpAtomNumber[i],
                    tmpElements[i],
                    tmpAtomTypes[i],
                    tmpEpsilsons[i],
                    tmpSigmas[i],
                    tmpCharges[i]
            ));
        }
        
    }
    
    /**
     * scanParticle method
     *  Optimize, scan and optimize the input particle
     * @param aJobTaskRecordList
     *  Job task records
     */
    private static void scanParticle() {
        String tmpForcefieldName;
        String tmpKeyContent;
        String tmpKeyPathName;
        String tmpParticleXyzName;
        String tmpOptXyzDirName;
        String tmpOptXyzLogName;
        String tmpOptXyzName;
        String tmpOptArcName;
        String tmpFileName;
        String tmpFileTxyzName;
        String tmpXyzFileName;
        String tmpOutputName;
        String tmpSearchString;
        String tmpSourceName;
        String tmpTargetName;
        String[] tmpCmdList;
        Path tmpOptXyzDir;
        Path tmpOptXyzFile;
        Path tmOptArcFile;
        Path tmpSource;
        Path tmpTarget;
        TinkerXYZ tmpTinkerXYZ;
        TinkerXYZ tmpTinkerXYZ0;
        TinkerXYZ tmpAfterScan;
        Process tmpProcess;
        ProcessBuilder tmpPB;
        List<Double> tmpEnergyValues;
        LinkedList<String> tmpForcefieldList;
        List<String> tmpMinimumList;
        Double[] tmpEnergyList;
        Integer[] tmpEnergyIndices;
        double[][][] tmpCoords;
        int tmpKeyIndex;
        ArrayIndexComparator tmpComparator;
        
        tmpProcess = null;
        tmpForcefieldList = new LinkedList<>();
        
        if (!forcefield_IE.equals("")) {
            tmpForcefieldList.add(forcefield_IE);
        }
        if (!forcefield_CN.equals("")) {
            tmpForcefieldList.add(forcefield_CN);
        }
        
        for (String tmpForcefield : tmpForcefieldList) {
            
            for (String tmpParticle : particleNames) {
                tmpParticleXyzName = tmpParticle + ".xyz";
                tmpOptXyzDirName = optXYZDirectory
                        + FILESEPARATOR 
                        + tmpForcefield 
                        + FILESEPARATOR 
                        + tmpParticle;
                tmpOptXyzDir = Paths.get(tmpOptXyzDirName);
                if (!Files.exists(tmpOptXyzDir)) {
                    try {
                        Files.createDirectories(tmpOptXyzDir);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during creating OptXyz directory.", 
                                ex);
                    }
                }
                
                // Check whether already scanned
                tmpOptXyzFile = Paths.get(tmpOptXyzDirName
                        + FILESEPARATOR
                        + tmpParticleXyzName);
                if (!Files.exists(tmpOptXyzFile)) {
                    tmpSourceName = sourceDirectory
                            + FILESEPARATOR
                            + tmpForcefield
                            + FILESEPARATOR
                            + tmpParticle
                            + ".xyz";
                    tmpTargetName = tmpOptXyzDirName
                            + FILESEPARATOR
                            + tmpParticleXyzName;
                    tmpSource = Paths.get(tmpSourceName);
                    tmpTarget = Paths.get(tmpTargetName);
                    try {
                        Files.copy(tmpSource, tmpTarget, StandardCopyOption
                                .REPLACE_EXISTING);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during copying .xyz file.", ex);
                    }
                    tmpTinkerXYZ = new TinkerXYZ(tmpForcefield,
                            tmpParticle,
                            tmpTargetName);
                    tmpKeyPathName = tmpOptXyzDirName 
                            + FILESEPARATOR 
                            + tmpParticle 
                            + ".key";
                    if (tmpForcefield.equals("OPLSAALIGPARGEN")) {
                        tmpForcefieldName = "oplsaa";
                    } else {
                        tmpForcefieldName = tmpForcefield;
                    }
                    tmpKeyContent = "# Force Field Selection"
                            + LINESEPARATOR 
                            + "PARAMETERS\t\""
                            + parameterDirectory
                            + "/"
                            + tmpForcefieldName.toLowerCase() + ".prm\""
                            + LINESEPARATOR
                            + "DIELECTRIC\t" 
                            + dielectricConstant
                            + LINESEPARATOR;
                    if (tmpForcefield.equals("OPLSAALIGPARGEN")) {
                        tmpKeyIndex = particleNames.indexOf(tmpParticle);
                        if (prmContent1[tmpKeyIndex] == null || 
                                prmContent1[tmpKeyIndex].isEmpty()) {
                            readPrm(tmpForcefield);
                        }
                        tmpKeyContent += prmContent1[tmpKeyIndex];
                    }
                    
                    // Write .key file
                    MIPETUTIL.writeKeyFile(tmpKeyPathName, tmpKeyContent);
                    
                    // Create optimized .xyz files for the particles via Tinker optimize.
                    try {
                        tmpProcess = new ProcessBuilder(tinkerOptimize,
                            tmpOptXyzDirName
                            + FILESEPARATOR 
                            + tmpParticleXyzName,
                            Double.toString(optimizeRmsGradient))
                            .start();
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE,
                                "IOException during process start.", ex);
                    } finally {
                        if (tmpProcess != null) {
                            OutputStream tmpOutput;
                            try (InputStream tmpInput = tmpProcess
                                    .getInputStream()) {
                                tmpOptXyzLogName = tmpOptXyzDirName
                                        + FILESEPARATOR
                                        + tmpParticle
                                        + "_preoptimize.log";
                                tmpOutput = new FileOutputStream(
                                        tmpOptXyzLogName, false);
                                tmpInput.transferTo(tmpOutput);
                                tmpOutput.close();
                            } catch (IOException ex) {
                                LOGGER.log(Level.SEVERE, 
                                        "IOException during writing .log file.",
                                        ex);
                            }
                            try {
                                tmpProcess.waitFor();
                            } catch (InterruptedException ex) {
                                    LOGGER.log(Level.SEVERE, 
                                            "InterruptException during processing optimize.exe",
                                            ex);
                            }
                            tmpProcess.destroy();
                        }
                    }
                    tmpTinkerXYZ0 = new TinkerXYZ(
                            tmpForcefield,
                            tmpParticle,
                            tmpOptXyzDirName 
                            + FILESEPARATOR 
                            + tmpParticle 
                            + ".xyz_2");
                    tmpTinkerXYZ.setForcefieldName(forcefield_IE);
                    tmpTinkerXYZ.setCoordinateList1(tmpTinkerXYZ0
                            .getCoordinateList1(), isTinkerOn);
                    tmpSource = Paths.get(tmpOptXyzDirName, tmpParticle 
                            + ".xyz_2");
                    try {
                        Files.deleteIfExists(tmpSource);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during deleting .xyz_2 file", ex);
                    }
                    tmpXyzFileName = tmpOptXyzDirName 
                            + FILESEPARATOR 
                            + tmpParticleXyzName;
                    try (PrintWriter tmpOut = new PrintWriter((tmpXyzFileName))) {
                        tmpOut.print(tmpTinkerXYZ.getFileContent());
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE,
                                "IOException during writing .xyz file.", ex);
                    }
                    
                    // Run Tinker Scan.exe with optimized particles
                    tmpPB = new ProcessBuilder();
                    tmpPB.redirectErrorStream(true);
                    tmpCmdList = new String[]{tinkerScan,
                        tmpXyzFileName,
                        scanProgram,
                        nSearchDirection,
                        energyThreshold,
                        rmsGradient};
                    tmpOutputName = tmpOptXyzDirName
                            + FILESEPARATOR
                            + tmpParticle
                            + ".out";
                    tmpPB.redirectOutput(new File(tmpOutputName));
                    tmpPB.command(tmpCmdList);
                    try {
                        tmpProcess = tmpPB.start();
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE,
                                "IOException during tinker's scan.exe", ex);
                    } finally {
                        if (tmpProcess != null) {
                            // This is necessary because .waitFor() will hang otherwise
                            try (BufferedReader tmpBR = new BufferedReader(
                                    new InputStreamReader(tmpProcess.getInputStream()))) {
                                while (tmpBR.readLine() != null ) {
                                }
                            } catch (IOException ex) {
                            }
                            try {
                                tmpProcess.waitFor();
                            } catch (InterruptedException ex) {
                                LOGGER.log(Level.SEVERE,
                                        "IOException during tinker's scan.exe", ex);
                            }
                            tmpProcess.destroy();
                        }
                    }
                    tmpOptArcName = tmpOptXyzDirName 
                            + FILESEPARATOR 
                            + tmpParticle 
                            + ".arc";
                    tmOptArcFile = Paths.get(tmpOptArcName);
                    if (Files.exists(tmOptArcFile) && 
                            !Files.isDirectory(tmOptArcFile)) { 
                        tmpCoords = tmpTinkerXYZ.readCoordFromArc(tmpOptArcName);

                        for (int j = 0; j < tmpCoords.length; j++) {
                            tmpTinkerXYZ.setCoordinateList1(tmpCoords[j], 
                                    isTinkerOn);
                            tmpXyzFileName = tmpOptXyzDirName 
                                + FILESEPARATOR
                                + tmpParticle
                                + "_"
                                + j
                                + ".txyz";
                            tmpTinkerXYZ.writeToTxyzFile(tmpXyzFileName);
                        }

                        tmpSearchString = MIPETUTIL.getResourceString(
                                "MIPETScanMinimumEnergyText");
                        tmpMinimumList = MIPETUTIL.findList(tmpOutputName, 
                                tmpSearchString);
                        if (!tmpMinimumList.isEmpty()) {
                            tmpEnergyValues = new LinkedList<>();

                            for (int j = 0; j < tmpMinimumList.size(); j++) {
                                tmpEnergyValues.add(Double
                                        .valueOf(tmpMinimumList
                                        .get(j)
                                        .substring(58, 68)));
                            }

                            tmpEnergyList = new Double[tmpEnergyValues.size()];

                            for (int j = 0; j < tmpEnergyValues.size(); j++) {
                                tmpEnergyList[j] = tmpEnergyValues.get(j);
                            }

                            tmpComparator = new ArrayIndexComparator(
                                    tmpEnergyList);
                            tmpEnergyIndices = tmpComparator.createIndexArray();
                            Arrays.sort(tmpEnergyIndices, tmpComparator);
                            tmpOptXyzName = tmpOptXyzDirName 
                                    + FILESEPARATOR
                                    + tmpParticle 
                                    + "_0.txyz";
                            tmpAfterScan = new TinkerXYZ(tmpForcefield,
                                    tmpParticle,
                                    tmpOptXyzName); 
                            tmpTinkerXYZ.setCoordinateList1(tmpAfterScan
                                    .getCoordinateList1(), isTinkerOn);

                            // Convert tinker xyz files to .xyz file
                            for (int j = 0; j < tmpEnergyValues.size(); j++) {
                                tmpFileTxyzName = tmpOptXyzDirName 
                                    + FILESEPARATOR
                                    + tmpParticle
                                    + "_" + j +".txyz";
                                tmpTinkerXYZ0 = new TinkerXYZ(tmpForcefield,
                                        tmpParticle,
                                        tmpFileTxyzName);
                                tmpFileName = tmpOptXyzDirName 
                                    + FILESEPARATOR
                                    + tmpParticle
                                    + "_o" + j + ".xyz";
                                tmpTinkerXYZ0.writeToXyzFile(tmpFileName);
                            }

                            // Copy configuration data with lowest energy to Particle.xyz
                            tmpSourceName = tmpOptXyzDirName 
                                    + FILESEPARATOR
                                    + tmpParticle
                                    + "_0.txyz";
                            tmpTargetName = tmpOptXyzDirName 
                                    + FILESEPARATOR
                                    + tmpParticle
                                    + ".xyz";
                            tmpSource = Paths.get(tmpSourceName);
                            tmpTarget = Paths.get(tmpTargetName);
                            try {
                                if (Files.exists(tmpSource)) {Files
                                        .copy(tmpSource, tmpTarget,
                                            StandardCopyOption
                                                    .REPLACE_EXISTING);
                                }
                            } catch (IOException ex) {
                                LOGGER.log(Level.SEVERE, 
                                        "IOException during copying .xyz file", 
                                        ex);
                            }

                            // Optimize after scan.exe
                            try {
                                tmpProcess = new ProcessBuilder(tinkerOptimize,
                                        tmpOptXyzDirName
                                        + FILESEPARATOR 
                                        + tmpParticleXyzName,
                                        Double.toString(optimizeRmsGradient))
                                        .start();
                                try (InputStream tmpInput = tmpProcess
                                        .getInputStream()) {
                                    tmpOptXyzLogName = tmpOptXyzDirName
                                        + FILESEPARATOR 
                                        + tmpParticle
                                        + "_afteroptimize.log";
                                    try (OutputStream tmpOutput = 
                                            new FileOutputStream (tmpOptXyzLogName, 
                                                    false)) {
                                        tmpInput.transferTo(tmpOutput);
                                    }
                                }
                            } catch (IOException ex) {
                                LOGGER.log(Level.SEVERE, "IOException during writing .log file", ex);
                            } finally {
                                if (tmpProcess != null) {
                                    try {
                                        tmpProcess.waitFor();
                                    } catch (InterruptedException ex) {
                                        LOGGER.log(Level.SEVERE, "InterruptedException during writing optimize.exe", ex);
                                    } 
                                    tmpProcess.destroy();
                                }
                            }
                            tmpTinkerXYZ0 = new TinkerXYZ(tmpForcefield,
                                    tmpParticle,
                                    tmpOptXyzDirName 
                                    + FILESEPARATOR 
                                    + tmpParticle 
                                    + ".xyz_2");
                            tmpTinkerXYZ.setCoordinateList1(tmpTinkerXYZ0
                                    .getCoordinateList1(), isTinkerOn);
                            tmpSource = Paths.get(tmpOptXyzDirName, 
                                    tmpParticle + ".xyz_2");
                            try {
                                Files.deleteIfExists(tmpSource);
                            } catch (IOException ex) {
                                LOGGER.log(Level.SEVERE, "IOException during deleting .xyz_2 file.", ex);
                            }
                            tmpXyzFileName = tmpOptXyzDirName 
                                    + FILESEPARATOR 
                                    + tmpParticleXyzName;
                            try (PrintWriter tmpOut = new PrintWriter(
                                    tmpXyzFileName)) {
                                tmpOut.print(tmpTinkerXYZ.getFileContent());
                            } catch(IOException ex) {
                                LOGGER.log(Level.SEVERE,
                                        "IOException during writing .xyz file.", ex);
                            }
                        } 
                    }
                }
            }
        }
        
        if (forcefield_CN.equals("OPLSAALIGPARGEN")) {
            readPrm(forcefield_IE);
        }
    }
    
    /**
     * Calculate the intermolecular energies for different distances.
     * 
     * @param aParticlePair Names of both particles
     * @param aDistances Distances of both particles (centre to centre)
     * @param aTinkerXYZ1 Coordinate data of first particle
     * @param aTinkerXYZ2 Coordinate data of second particle
     * @param aRotData1 Coordinates of first particle after rotations
     * @param aRotData2 Coordinates of second particle after rotations
     * @param aMinEnergy Old minimum energy value  
     * @return Intermolecular energies for different distances, temperatures 
     *           and fractions.
     */
    private static EnergyRecord getInterMolecularEnergy(
            String aParticlePair,
            double[] aDistances,
            TinkerXYZ aTinkerXYZ1, 
            TinkerXYZ aTinkerXYZ2,
            double[][][] aRotData1,
            double[][][] aRotData2,
            double aMinEnergy) {
        int tmpAtomNumber;
        int tmpDistanceNumber;
        int tmpConfigNumber;
        int tmpChunkNumber;
        int tmpChunkSize;
        int tmpChunkRemainder;
        int tmpChunkIndex;
        int tmpRot2StartIndex;
        int tmpRot2EndIndex;
        double[][] tmpEnergyDatas;
        double[][][] tmpRotData1;
        double[][][] tmpRotData2;
        String tmpPath;
        String[] tmpCmdList;
        ArrayList<MIPETAnalyze> tmpTaskList;
        TinkerXYZ tmpTinkerXyz1;
        TinkerXYZ tmpTinkerXyz2;
        TinkerXYZ tmpTinkerXYZ;
        ExecutorService executor;
        
        // Calculate chunk size
        tmpDistanceNumber = aDistances.length;
        tmpConfigNumber = aRotData1.length * aRotData2.length;
        if (tmpConfigNumber < 1000) {
            tmpChunkNumber = 1;
            tmpChunkSize = tmpConfigNumber;
            tmpChunkRemainder = 0;
        } else {
            tmpChunkNumber = cpuCoreNumber;
            tmpChunkSize = (int)Math.ceil(aRotData2.length / tmpChunkNumber);
            tmpChunkRemainder = aRotData2.length % tmpChunkSize;
        }

        // Calculate intermolecular energy using TINKER analyze
        tmpTinkerXyz1 = aTinkerXYZ1.clone();
        tmpTinkerXyz2 = aTinkerXYZ2.clone();
        tmpRotData1 = aRotData1;
        tmpEnergyDatas = new double[tmpDistanceNumber][];
        executor = Executors.newFixedThreadPool(cpuCoreNumber);
        tmpTaskList = new ArrayList<>(2000);
        tmpPath = scratchDirectory 
                + FILESEPARATOR 
                + aParticlePair 
                + ".arc";
        
        for (int i = 0; i < tmpDistanceNumber; i++) {
            tmpRotData2 = VectorUtil.moveX(aRotData2, aDistances[i]);
            tmpChunkIndex = 0;
            tmpRot2StartIndex = 0;
            
            while (tmpChunkIndex < tmpChunkNumber) {
                tmpTinkerXYZ = new TinkerXYZ(tmpTinkerXyz1, 
                        tmpTinkerXyz2,
                        isTinkerOn);
                tmpAtomNumber = tmpTinkerXYZ.getAtomNumber();
                tmpRot2EndIndex = tmpRot2StartIndex + tmpChunkSize;
                if (tmpChunkIndex == tmpChunkNumber - 1) {
                    tmpRot2EndIndex += tmpChunkRemainder;
                } 
                double[][][] tmpRotPart2;
                tmpRotPart2 = Arrays.copyOfRange(tmpRotData2,
                        tmpRot2StartIndex,
                        tmpRot2EndIndex);
                tmpCmdList = new String[]{tinkerAnalyze, tmpPath + i + "_"
                        + tmpChunkIndex, "E"};
                tmpTaskList.add(new MIPETAnalyze(
                        tmpTinkerXYZ,
                        isTinkerOn,
                        i,
                        tmpChunkIndex,
                        tmpAtomNumber,
                        minAtomDistance,
                        tmpRotData1,
                        tmpRotPart2,
                        scratchDirectory,
                        tmpCmdList,
                        molecules));
                tmpChunkIndex++;
                tmpRot2StartIndex = tmpRot2EndIndex;
            }
            
        }
        
        int tmpDistMinIndex;
        int tmpChunkMinIndex;
        int tmpTaskIndex;
        Double tmpDistMinEnergy;
        double tmpPartMinEnergy;
        
        tmpTaskIndex = 0;
        tmpDistMinIndex = 0;
        tmpChunkMinIndex = 0;
        tmpPartMinEnergy = 1E10;
        List<Future<ArrayList<Double>>> tmpFutures = null;
        Future<ArrayList<Double>> tmpFuture;
        ArrayList<Double> tmpDistEnergies = new ArrayList<>(tmpConfigNumber);
        
        try {            
            tmpFutures = executor.invokeAll(tmpTaskList);
        } catch (InterruptedException ex) {
            Thread.currentThread().interrupt();       
        }
        executor.shutdown();
            
        for (int i = 0; i < tmpDistanceNumber; i++) {
            
            for (int j = 0; j < tmpChunkNumber; j++) {
                tmpFuture = tmpFutures.get(tmpTaskIndex);
                try {
                    tmpDistEnergies.addAll(tmpFuture.get());
                    tmpDistMinEnergy = tmpFuture.get().get(0);
                    if (tmpDistMinEnergy < tmpPartMinEnergy) {
                        tmpPartMinEnergy = tmpDistMinEnergy;
                        tmpDistMinIndex = i;
                        tmpChunkMinIndex = j;
                    }
                    tmpTaskIndex++;
                    
                    // Store all energies at same distance to tmpEnergyDatas[i]
                    if (j == tmpChunkNumber - 1) {
                        tmpEnergyDatas[i] = MIPETUTIL
                                .toPrimitive(tmpDistEnergies);
                        Arrays.sort(tmpEnergyDatas[i]);
                        tmpDistEnergies.clear();
                    }
                } catch (InterruptedException | ExecutionException ex) {
                    LOGGER.log(Level.SEVERE,
                            "InterruptException during handling tmpFuture object.",
                            ex);
                } 
            }
            
        }
        
        // Export .xyz file with lowest intermolecular energy
        if (tmpPartMinEnergy < aMinEnergy) {
            String tmpSourceFileName;
            String tmpExportFileName;
            Path tmpSource;
            Path tmpTarget;

            tmpSourceFileName = scratchDirectory
                    + FILESEPARATOR
                    + aParticlePair
                    + "_"
                    + tmpDistMinIndex
                    + "_"
                    + tmpChunkMinIndex
                    + ".0";
            tmpExportFileName = scratchDirectory 
                    + FILESEPARATOR 
                    + aParticlePair
                    + ".0";
            tmpSource = Paths.get(tmpSourceFileName);
            tmpTarget = Paths.get(tmpExportFileName);
            try {
                Files.copy(tmpSource, tmpTarget, 
                        StandardCopyOption.REPLACE_EXISTING);
            } catch(IOException ex) {
                LOGGER.log(Level.SEVERE,
                        "IOException during copying .0 file.", ex);
            }
        }
        
        //<editor-fold defaultstate="collapsed" desc="Calculate intermolecular energy of all configurations">
        // If boltzmannFraction == 0.0, no averaging, min energy value of each configuration is taken
        // If fractionForAverage = 1.0 all configurational E(nonbonded) values are used for "Boltzmann average" calculation
        // 0.0 < fractionForAverage < 1.0: All configurational E(nonbonded) values are sorted ascending and
        // the lower "numberOfValues*fractionForAverage" E(nonbonded) values are used for "Boltzmann average" calculation
        // Example: For 144x144x16 = 331776 E(nonbonded) values for a specific molecule distance r and
        // a fractionForAverage of 0.25 the lowest Round(331776x0.25) = 82944 E(nonbonded) values are used for
        // "Boltzmann average" calculation only
        int tmpFractionToMax;
        double tmpEmin;
        double tmpTempGasconst;
        double tmpWgtEmin;
        double[] tmpWeights;
        double[] tmpEnergyDataFraction;
        double[] tmpEmins;
        double[] tmpWgtEmins;

        // Find tmpMinEnergy
        tmpEmins = new double[tmpDistanceNumber];
        tmpWgtEmins = new double[tmpDistanceNumber];
        tmpTempGasconst = temperature * GASCONST;
        tmpEmin = tmpPartMinEnergy;
        tmpWgtEmin = 100.;

        for (int i = 0; i < tmpDistanceNumber; i++) {
            tmpEmins[i] = tmpEnergyDatas[i][0];
            tmpFractionToMax = (int)(tmpEnergyDatas[i].length 
                    * boltzmannFraction);
            tmpEnergyDataFraction = new double[tmpFractionToMax];
            tmpWeights = new double[tmpFractionToMax];

            for (int j = 0; j < tmpFractionToMax; j++) {
                tmpEnergyDataFraction[j] = tmpEnergyDatas[i][j];
                tmpWeights[j] = Math.exp(-(tmpEnergyDatas[i][j] - tmpEmin)
                        / tmpTempGasconst);
            }

            tmpWgtEmins[i] = MIPETUTIL.productSum(tmpWeights, 
                    tmpEnergyDataFraction) / MIPETUTIL.sum(tmpWeights);
            if (tmpWgtEmin > tmpWgtEmins[i]) {
                tmpWgtEmin = tmpWgtEmins[i];
            }
        }

        //</editor-fold>
        
        return new EnergyRecord(aDistances,
                tmpEmins,
                tmpWgtEmins, 
                aDistances[tmpDistMinIndex],
                tmpEmin,
                tmpWgtEmin);
    }
    
    /**
     * Determine coordination numbers
     * @param aJobTaskRecordList: Job task records
     */
    private static void getCoordinationNumbers(
            ArrayList<JobTaskRecord> aJobTaskRecordList) {
        
        boolean tmpIsSameParticle;
        int tmpJobTaskLength;
        int tmpJobNumber;
        int tmpJobIndex;
        int tmpXyzID1;
        int tmpXyzID2;
        long tmpTotalCNTime;
        long boxMinimizationTime;
        long tmpDynamicWarmUpTime;
        int[] tmpCNMax;
        int[] tmpCNMin;
        int[][] tmpCNs;
        double[] tmpCNMeans;
        double[] tmpVdWSolventVolumes;
        double[] tmpVdWSolutVolumes;
        double[] tmpBoxLengths;
        double[] tmpStdDeviation;
        String tmpCurrentDir;
        String tmpParticle1;
        String tmpParticle2;
        String tmpParticlePair;
        String tmpParticleLogFileName;
        String tmpDatFileName;
        String tmpKeyContent;
        String tmpKeyFixContent;
        String tmpKeyMaxiter;
        String tmpKeyFileName;
        String[] tmpCmdList;
        Path tmpSourceFile;
        Path tmpSourceFile2;
        Path tmpTargetDir;
        Path tmpParticleLogFile;
        Path tmpResultPath;
        Path tmpSource;
        Path tmpTarget;
        ProcessBuilder tmpProcessBuilder;
        Process tmpProcess;
        ExecutorService tmpExecutor;
        LinkedList<MIPETCN> tmpTaskList = new LinkedList<>();
        
        tmpJobTaskLength = aJobTaskRecordList.size();
        tmpTotalCNTime = System.currentTimeMillis();
        tmpJobNumber = 0;
        
        // <editor-fold defaultstate="collapsed" desc="Make directories and copy .xyz datas">
        for (int i = 0; i < tmpJobTaskLength; i++) {
            if (aJobTaskRecordList.get(i).hasCNJob()) {
                tmpJobNumber++;
            }
        }
        
        tmpVdWSolventVolumes = new double[tmpJobNumber];
        tmpVdWSolutVolumes = new double[tmpJobNumber];
        tmpBoxLengths = new double[tmpJobNumber];
        tmpJobIndex = 0;
        
        // Read xyzContent and prmContent
        if (!forcefield_IE.equals(forcefield_CN)) {
            if (isConformationalAnalysis) {
                readXyz(false, forcefield_CN);
            } else {
                readXyz(true, forcefield_CN);
            }
            if (forcefield_CN.equals("OPLSAALIGPARGEN")) {
                readPrm(forcefield_CN);
            }
        }
        
        for (int i = 0; i < tmpJobTaskLength; i++) {
            if (aJobTaskRecordList.get(i).hasCNJob()) {
                tmpParticle1 = aJobTaskRecordList.get(i).particleName1();
                tmpParticle2 = aJobTaskRecordList.get(i).particleName2();
                tmpParticlePair = tmpParticle1 + "_" + tmpParticle2;
                tmpIsSameParticle = aJobTaskRecordList.get(i)
                        .isSameParticle();
                tmpXyzID1 = particleNames.indexOf(tmpParticle1);
                if (tmpIsSameParticle) {
                    tmpXyzID2 = tmpXyzID1;
                } else {
                    tmpXyzID2 = particleNames.indexOf(tmpParticle2);
                }

                // Create directory in scratch directory
                tmpTargetDir = Paths.get(scratchDirectory 
                        + FILESEPARATOR 
                        + forcefield_CN
                        + FILESEPARATOR 
                        + tmpParticlePair);
                if (!Files.exists(tmpTargetDir)) {
                    try {
                        Files.createDirectories(tmpTargetDir);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE,
                            "IOException during creating scratch/forcefield/particlepair directory.",
                            ex);
                    }
                }
                
                // Copy .xyz files
                tmpTarget = Paths.get(scratchDirectory
                        + FILESEPARATOR
                        + forcefield_CN
                        + FILESEPARATOR 
                        + tmpParticlePair
                        + FILESEPARATOR 
                        + tmpParticle1
                        + ".xyz");
                try {
                    Files.writeString(tmpTarget, xyzContent1[tmpXyzID1]);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during copying files to scratch.", ex);
                }
                if (!tmpIsSameParticle) {
                    tmpTarget = Paths.get(scratchDirectory
                            + FILESEPARATOR
                            + forcefield_CN
                            + FILESEPARATOR 
                            + tmpParticlePair
                            + FILESEPARATOR 
                            + tmpParticle2
                            + ".xyz");
                    try {
                        Files.writeString(tmpTarget, xyzContent1[tmpXyzID2]);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during copying files to scratch.", ex);
                    }
                }
                tmpVdWSolventVolumes[tmpJobIndex] = MIPETUTIL
                        .getVdwVolume(smiles.get(tmpParticle2));
                if (tmpIsSameParticle) {
                    tmpVdWSolutVolumes[tmpJobIndex] = 
                            tmpVdWSolventVolumes[tmpJobIndex];
                } else {
                    tmpVdWSolutVolumes[tmpJobIndex] = MIPETUTIL.getVdwVolume(
                            smiles.get(tmpParticle1));
                }
                
                // Calculate water volume ratio - ratio of Vparticle and Vvdw 
                //   of water 1.7297
                tmpBoxLengths[tmpJobIndex] =  Math.pow(WATERVOLUMERATIO 
                        * solventMoleculeNumber 
                        * tmpVdWSolventVolumes[tmpJobIndex], 1.0/3.0);
                tmpJobIndex++;
            }
        }

        // </editor-fold>
        
        // <editor-fold defaultstate="collapsed" desc="Generate key file for dynamic">
        // OPENMP-THREADS is automatically set to core number when it is not set.
        // MIPET also works parallely. Tests showed best performance is 
        //   achievable when OPENMP-THREADS is set to 1.
        tmpKeyFixContent = 
            """
            EWALD
            OPENMP-THREADS 1
            THERMOSTAT ANDERSEN
            STEEPEST-DESCENT
            RANDOMSEED 123456
            """;
        if (solventMoleculeNumber >= 400) {
            tmpKeyFixContent += "NEIGHBOR-LIST" + LINESEPARATOR;
        }

        // </editor-fold>

        // <editor-fold defaultstate="collapsed" desc="Make simulation box">
        Boolean tmpHasCNHeadLine;
        String tmpOldXYZFileName;
        String tmpNewXYZFileName;
        String tmpResultPathName;
        String tmpBuildLogName;
        String tmpSoakLogName;
        String tmpMinimizeLogName;
        String tmpForcefieldName;
        int tmpPrm1ID;
        int tmpPrm2ID;
        
        tmpJobIndex = 0;

        for (int i = 0; i < tmpJobTaskLength; i++) {
            if (aJobTaskRecordList.get(i).hasCNJob()) {
                tmpHasCNHeadLine = false;
                tmpParticle1 = aJobTaskRecordList.get(i).particleName1(); // solute
                tmpParticle2 = aJobTaskRecordList.get(i).particleName2(); // solvent
                tmpPrm1ID = particleNames.indexOf(tmpParticle1);
                tmpPrm2ID = particleNames.indexOf(tmpParticle2);
                tmpParticlePair = tmpParticle1 + "_" + tmpParticle2;
                tmpIsSameParticle = aJobTaskRecordList.get(i)
                        .isSameParticle();
                tmpResultPathName = aJobTaskRecordList.get(i)
                        .result_CN_PathName();
                tmpCurrentDir = scratchDirectory
                        + FILESEPARATOR
                        + forcefield_CN
                        + FILESEPARATOR
                        + tmpParticlePair
                        + FILESEPARATOR;
                tmpResultPath = Paths.get(aJobTaskRecordList.get(i)
                        .result_CN_PathName());
                tmpSource = Paths.get(tmpCurrentDir 
                        + tmpParticlePair 
                        + ".xyz_2");
                tmpSourceFile2 = Paths.get(tmpCurrentDir
                        + tmpParticlePair 
                        + ".xyz_2");
                tmpTarget = Paths.get(tmpCurrentDir
                        + tmpParticlePair 
                        + ".xyz");

                // <editor-fold defaultstate="collapsed" desc="Make .key file for dynamic">
                if (forcefield_CN.equals("OPLSAALIGPARGEN")) {
                    tmpForcefieldName = "oplsaa";
                    
                } else {
                    tmpForcefieldName = forcefield_CN.toLowerCase();
                }
                tmpKeyContent = tmpKeyFixContent
                        + "PARAMETERS\t\""
                        + parameterDirectory
                        + "/"
                        + tmpForcefieldName
                        + ".prm\""
                        + LINESEPARATOR    
                        + "a-axis "
                        + tmpBoxLengths[tmpJobIndex]
                        + LINESEPARATOR;
                if (isRATTLE) {
                    tmpKeyContent += "RATTLE" + LINESEPARATOR;
                }
                tmpKeyFileName = tmpCurrentDir
                        + tmpParticlePair
                        + ".key";
                if (forcefield_CN.equals("OPLSAALIGPARGEN")) {
                    tmpKeyContent += prmContent1[tmpPrm1ID];
                    if (!tmpIsSameParticle) {
                        tmpKeyContent += prmContent2[tmpPrm2ID];
                    } 
                }
                MIPETUTIL.writeKeyFile(tmpKeyFileName, tmpKeyContent);

                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="Build solvent box">
                tmpProcessBuilder = new ProcessBuilder(tinkerXYZEdit,
                        tmpCurrentDir + tmpParticle2 + ".xyz",
                        "-k",
                        tmpKeyFileName,
                        "23",  // Option 23: Create and Fill a Periodic Boundary Box
                        Integer.toString(solventMoleculeNumber),
                        Double.toString(tmpBoxLengths[tmpJobIndex]),
                        Double.toString(tmpBoxLengths[tmpJobIndex]),
                        Double.toString(tmpBoxLengths[tmpJobIndex]),
                        "Y");
                try {
                    tmpProcess = tmpProcessBuilder.start();
                    try (BufferedReader tmpBR = new BufferedReader(
                            new InputStreamReader(tmpProcess.getInputStream()))) {
                        if (isLogBuildBox) {
                            tmpBuildLogName = tmpResultPathName
                                + FILESEPARATOR
                                + tmpParticle1
                                + "_"
                                + tmpParticle2
                                + "_build.log";
                            BufferedWriter tmpBW = new BufferedWriter(
                                    new  FileWriter(tmpBuildLogName));
                            String tmpLine;
                            
                            while ((tmpLine = tmpBR.readLine()) != null) {
                                tmpBW.append(tmpLine);
                            }
                            
                        } else {
                            // This loop is necessary for linux version.
                            while (tmpBR.readLine() != null) {
                            }
                            
                        } 
                    }
                    tmpProcess.waitFor();
                    tmpProcess.destroy();
                } catch (InterruptedException ex) {
                    Thread.currentThread().interrupt();
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, "IOException during XYZEdit.", ex);
                }

                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="Add solute to the simulation box">
                if (!tmpIsSameParticle) {
                    tmpProcessBuilder = new ProcessBuilder(tinkerXYZEdit,
                        tmpCurrentDir + tmpParticle1 + ".xyz",
                        "-k",
                        tmpKeyFileName,
                        "24", // Option 24: Soak Current Molecule in Box of Solvent
                        tmpCurrentDir + tmpParticle2 + ".xyz_2");
                    try {
                        tmpProcess = tmpProcessBuilder.start();
                        try (BufferedReader tmpBR = new BufferedReader(
                                new InputStreamReader(tmpProcess
                                        .getInputStream()))) {
                            if (isLogSoakBox) {
                                tmpSoakLogName = tmpResultPathName
                                        + FILESEPARATOR
                                        + tmpParticle1
                                        + "_"
                                        + tmpParticle2
                                        + "_soak.log";
                                BufferedWriter tmpBW = new BufferedWriter(
                                        new FileWriter(tmpSoakLogName));
                                String tmpLine;
                                
                                while ((tmpLine = tmpBR.readLine()) != null) {
                                    tmpBW.append(tmpLine);
                                }
                                
                            } else {
                                // This loop is necessary for linux version
                                while (tmpBR.readLine() != null) {
                                }
                                
                            } 
                        }
                        tmpProcess.waitFor();
                        tmpProcess.destroy();
                    } catch (InterruptedException ex) {
                        Thread.currentThread().interrupt();
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during XYZEdit.", ex);
                    }
                    
                }
                tmpOldXYZFileName = tmpCurrentDir 
                        + tmpParticle1 
                        + ".xyz_2";
                tmpNewXYZFileName = tmpCurrentDir 
                        + tmpParticlePair 
                        + ".xyz";
                try {
                    Files.move(Paths.get(tmpOldXYZFileName), 
                        Paths.get(tmpNewXYZFileName), 
                        StandardCopyOption.REPLACE_EXISTING);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                                "IOException during rename .xyz_2 to .xyz.",
                                ex);
                }

                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="Add MAXITER to .key file">
                boxMinimizationTime = System.currentTimeMillis();
                tmpKeyMaxiter = "MAXITER " + minimizeMaxIteration;
                try (BufferedWriter tmpBW = new BufferedWriter(
                        new FileWriter(tmpKeyFileName, true))) {
                            tmpBW.write(tmpKeyMaxiter);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during writing .key file.", ex);
                }

                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="Minimize solvent box">
                tmpProcessBuilder = new ProcessBuilder(tinkerMinimize,
                        tmpCurrentDir + tmpParticlePair + ".xyz",
                        Double.toString(rmsMinimizeGradient));
                try {
                    tmpProcess = tmpProcessBuilder.start();
                    try (BufferedReader tmpBR = new BufferedReader(
                            new InputStreamReader(tmpProcess.getInputStream()))) {
                        if (isLogMinimizeBox) {
                            tmpMinimizeLogName = tmpResultPathName
                                    + FILESEPARATOR 
                                    + tmpParticlePair
                                    + "_minimize.log";
                            BufferedWriter tmpBW = new BufferedWriter(
                                    new FileWriter(tmpMinimizeLogName));
                            String tmpLine;
                                
                            while ((tmpLine = tmpBR.readLine()) != null) {
                                tmpBW.append(tmpLine);
                            }
                                
                        } else {
                            // This loop is necessary for linux version.
                            while (tmpBR.readLine() != null) {
                            }
                            
                        }
                    } 
                    tmpProcess.waitFor();
                    tmpProcess.destroy();
                } catch (InterruptedException ex) {
                    Thread.currentThread().interrupt();
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during XYZEdit.", ex);
                }
                if (!aJobTaskRecordList.get(i).isReverse()) {
                    tmpParticleLogFileName = tmpResultPath
                            + FILESEPARATOR
                            + tmpParticlePair
                            + "_log.txt";
                    tmpDatFileName = tmpResultPath
                            + FILESEPARATOR
                            + tmpParticlePair
                            + ".dat";
                } else {
                    tmpParticleLogFileName = tmpResultPath
                            + FILESEPARATOR
                            + tmpParticle2
                            + "_"
                            + tmpParticle1
                            + "_log.txt";
                    tmpDatFileName = tmpResultPath
                            + FILESEPARATOR
                            + tmpParticle2
                            + "_"
                            + tmpParticle1
                            + ".dat";
                }

                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="Write .log file">
                tmpParticleLogFile = Paths.get(tmpParticleLogFileName);
                if (!tmpIsSameParticle && !Files.exists(tmpParticleLogFile)) {
                    tmpHasCNHeadLine = true;
                }
                try {
                    BWParticleDat = new BufferedWriter(
                            new FileWriter(tmpDatFileName));
                    BWParticleLog = new BufferedWriter(
                            new FileWriter(tmpParticleLogFileName));
                    if (Files.exists(tmpSourceFile2) && !tmpHasCNHeadLine) {
                        Files.move(tmpSource, tmpTarget, 
                                StandardCopyOption.ATOMIC_MOVE);
                        BWParticleLog.append(LINESEPARATOR);
                        BWParticleLog.append("*** Coordination number ***");
                        BWParticleLog.append(LINESEPARATOR);
                    }
                    BWParticleLog.append("Time for solvent box minimization via Tinker \"Minimize\" [s]: ");
                    BWParticleLog.append(Double.toString( 
                            (System.currentTimeMillis() - boxMinimizationTime) 
                                    / 1000));
                    BWParticleLog.append(LINESEPARATOR);
                    BWParticleLog.close();
                    if (!tmpIsSameParticle) {
                        BWParticleDat.append("VdWVolume(");
                        BWParticleDat.append(tmpParticle1);
                        BWParticleDat.append(") [" + ANGSTROM + SUPERSCRIPTTHREE 
                                +"] = ");
                        BWParticleDat.append(Double
                                .toString(tmpVdWSolutVolumes[tmpJobIndex]));
                        BWParticleDat.append(LINESEPARATOR);
                    }
                    BWParticleDat.append("VdWVolume(");
                    BWParticleDat.append(tmpParticle2);
                    BWParticleDat.append(") [" + ANGSTROM + SUPERSCRIPTTHREE 
                                +"] = ");
                    BWParticleDat.append(String.format("%.4f", 
                            tmpVdWSolventVolumes[tmpJobIndex]));
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.close();
                } catch(IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during Minimize.", ex);
                }

                // Make backup of particle1_particle2.xyz
                tmpSource = tmpTarget;
                tmpTarget = Paths.get(tmpCurrentDir
                        + tmpParticlePair 
                        + ".bak");
                try{
                    Files.copy(tmpSource, tmpTarget, 
                            StandardCopyOption.REPLACE_EXISTING);
                } catch(IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during Copying to .bak file", ex);
                }
                tmpJobIndex++;

                // </editor-fold>
            }
        }

        // </editor-fold>
        
        // <editor-fold defaultstate="collapsed" desc="Simulation warm up">
        for (int i = 0; i < tmpJobTaskLength; i++) {
            if (aJobTaskRecordList.get(i).hasCNJob()) {
                tmpParticle1 = aJobTaskRecordList.get(i)
                        .particleName1(); // solute
                tmpParticle2 = aJobTaskRecordList.get(i)
                        .particleName2(); // solvent
                tmpParticlePair = tmpParticle1 + "_" + tmpParticle2;
                tmpCurrentDir = scratchDirectory
                        + FILESEPARATOR
                        + forcefield_CN
                        + FILESEPARATOR
                        + tmpParticlePair
                        + FILESEPARATOR;
                if (isTinker9) {
                    tmpCmdList = new String[] {tinkerDynamic,
                        "dynamic",
                        tmpCurrentDir + tmpParticlePair + ".xyz",
                        Integer.toString(warmUpStepNumber),
                        Double.toString(warmUpTimeStep),
                        Double.toString(warmUpPrintInterval),
                        Integer.toString(simulationType),
                        Double.toString(temperature)};
                } else {
                    tmpCmdList = new String[] {tinkerDynamic,
                        tmpCurrentDir + tmpParticlePair + ".xyz",
                        Integer.toString(warmUpStepNumber),
                        Double.toString(warmUpTimeStep),
                        Double.toString(warmUpPrintInterval),
                        Integer.toString(simulationType),
                        Double.toString(temperature)};
                }
                tmpTaskList.add(new MIPETCN(tmpCmdList,
                        aJobTaskRecordList.get(i),
                        scratchDirectory,
                        catchRadius,
                        true,
                        isLogWarmup,
                        isLogDynamic,
                        isTinker9));    
            }
        }

        tmpDynamicWarmUpTime = System.currentTimeMillis();
        tmpExecutor = Executors.newFixedThreadPool(cpuCoreNumber);
        try {
            tmpExecutor.invokeAll(tmpTaskList);
        } catch (InterruptedException ex) {
            Thread.currentThread().interrupt();       
        }
        tmpExecutor.shutdown();
        String[][] tmpLabelValues = new String[1][2];
        tmpLabelValues[0][0] = "Time for solvent box warm up via Tinker \"Dynamic\" [s]: ";
        tmpLabelValues[0][1] = Double.toString((System.currentTimeMillis() 
                - tmpDynamicWarmUpTime) / 1000);
        MIPETUTIL.writeParticleLog(aJobTaskRecordList, tmpLabelValues);
        
        // </editor-fold>
        
        // <editor-fold defaultstate="collapsed" desc="Start simulation & determine coordination number">
        LinkedList<MIPETCN> tmpCNTaskList = new LinkedList<>();
        List<Future<int[]>> tmpCN = new LinkedList<>();

        for (int i = 0; i < tmpJobTaskLength; i++) {
            if (aJobTaskRecordList.get(i).hasCNJob()) {
                tmpParticle1 = aJobTaskRecordList.get(i)
                        .particleName1(); // solute
                tmpParticle2 = aJobTaskRecordList.get(i)
                        .particleName2(); // solvent
                tmpParticlePair = tmpParticle1 + "_" + tmpParticle2;
                tmpCurrentDir = scratchDirectory
                        + FILESEPARATOR
                        + forcefield_CN
                        + FILESEPARATOR
                        + tmpParticlePair
                        + FILESEPARATOR;
                if (isTinker9) {
                    tmpCmdList = new String[] {tinkerDynamic,
                        "dynamic",
                        tmpCurrentDir + tmpParticlePair + ".xyz",
                        Integer.toString(stepNumber),
                        Double.toString(timeStep),
                        Double.toString(printInterval),
                        Integer.toString(simulationType),
                        Double.toString(temperature)};
                } else {
                    tmpCmdList = new String[] {tinkerDynamic,
                        tmpCurrentDir + tmpParticlePair + ".xyz",
                        Integer.toString(stepNumber),
                        Double.toString(timeStep),
                        Double.toString(printInterval),
                        Integer.toString(simulationType),
                        Double.toString(temperature)};
                }
                tmpCNTaskList.add(new MIPETCN(tmpCmdList,
                        aJobTaskRecordList.get(i),
                        scratchDirectory,
                        catchRadius,
                        false,
                        isLogWarmup,
                        isLogDynamic,
                        isTinker9));
            }
        }
        
        long tmpArcEvaluationTime = System.currentTimeMillis();
        if (tmpJobNumber > 0) {
            tmpExecutor = Executors.newFixedThreadPool(cpuCoreNumber);
            try {
                tmpCN = tmpExecutor.invokeAll(tmpCNTaskList);
            } catch (InterruptedException ex) {
                Thread.currentThread().interrupt();       
            }
            tmpExecutor.shutdown();
            tmpLabelValues = new String[2][2];
            tmpLabelValues[0][0] = "Time for solvent box analysis (includes \".arc-evaluation\" and \"MIPETCN\") [s]: ";
            tmpLabelValues[0][1] = Double.toString((System.currentTimeMillis() 
                    - tmpArcEvaluationTime) / 1000);
            tmpLabelValues[1][0] = "Total time for determining the coordination number [s]: ";
            tmpLabelValues[1][1] = String.valueOf((System
                    .currentTimeMillis() - tmpTotalCNTime) / 1000);
            MIPETUTIL.writeParticleLog(aJobTaskRecordList, tmpLabelValues);
            String tmpLabel = "Mean neighbor (" + temperature + " K): ";
            String[] tmpValues = new String[tmpCN.size()];
            tmpCNs = new int[tmpCN.size()][];
            tmpCNMeans = new double[tmpCN.size()];
            tmpStdDeviation = new double[tmpCN.size()];
            tmpCNMax = new int[tmpCN.size()];
            tmpCNMin = new int[tmpCN.size()];
            tmpJobIndex = 0;

            for (int i = 0; i < tmpJobTaskLength; i++) {
                if (aJobTaskRecordList.get(i).hasCNJob()) {
                    tmpParticle1 = aJobTaskRecordList.get(i)
                            .particleName1(); // solute
                    tmpParticle2 = aJobTaskRecordList.get(i)
                            .particleName2(); // solvent
                    tmpParticlePair = tmpParticle1 + "_" + tmpParticle2;
                    try {
                        tmpCNs[tmpJobIndex] = tmpCN.get(tmpJobIndex).get();
                    } catch (ExecutionException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "ExecutionException during get tmpCN", ex);
                    } catch (InterruptedException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "InterruptedException during get tmpCN", ex);
                    }
                    tmpCNMeans[tmpJobIndex] = (double)MIPETUTIL.sum(
                            tmpCNs[tmpJobIndex]) / tmpCNs[tmpJobIndex].length;
                    cnList.add(new ResultCNRecord(
                            tmpParticle1, 
                            tmpParticle2, 
                            tmpCNMeans[tmpJobIndex]));
                    tmpStdDeviation[tmpJobIndex] = MIPETUTIL
                            .standarddeviation(tmpCNs[tmpJobIndex], 
                                    tmpCNMeans[tmpJobIndex]);
                    tmpCNMin[tmpJobIndex] = MIPETUTIL
                            .getMin(tmpCNs[tmpJobIndex]);
                    tmpCNMax[tmpJobIndex] = MIPETUTIL
                            .getMax(tmpCNs[tmpJobIndex]);
                    tmpValues[tmpJobIndex] = Double.toString(
                            tmpCNMeans[tmpJobIndex]);
                    tmpCurrentDir = scratchDirectory
                            + FILESEPARATOR
                            + forcefield_CN
                            + FILESEPARATOR
                            + tmpParticlePair
                            + FILESEPARATOR;
                    tmpSourceFile = Paths.get(tmpCurrentDir 
                            + tmpParticlePair
                            + ".arc");
                    try {
                        Files.deleteIfExists(tmpSourceFile);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during deleting .arc file", ex);
                    }

                    // Restore of particle1_particle2.xyz
                    tmpSource = Paths.get(tmpCurrentDir
                            + tmpParticlePair
                            + ".bak");
                    tmpTarget = Paths.get(tmpCurrentDir
                            + tmpParticlePair 
                            + ".xyz");
                    try {
                        Files.copy(tmpSource, tmpTarget, 
                        StandardCopyOption.REPLACE_EXISTING);
                    } catch(IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                        "IOException during Copying to .bak file", ex);
                    }
                    tmpJobIndex++;
                }
            }

            // </editor-fold>

            // <editor-fold defaultstate="collapsed" desc="Write data">

            MIPETUTIL.writeParticleLog(aJobTaskRecordList, tmpLabel, tmpValues);
            MIPETUTIL.writeZij_Table(aJobTaskRecordList, tmpCNs, temperature);
            tmpJobIndex = 0;

            for (int i = 0; i < tmpJobTaskLength; i++) {
                if (aJobTaskRecordList.get(i).hasCNJob()) {
                    tmpParticle1 = aJobTaskRecordList.get(i)
                            .particleName1(); // solute
                    tmpParticle2 = aJobTaskRecordList.get(i)
                            .particleName2(); // solvent
                    if (!aJobTaskRecordList.get(i).isReverse()) {
                        tmpParticlePair = tmpParticle1 + "_" + tmpParticle2;
                    } else {
                        tmpParticlePair = tmpParticle2 + "_" + tmpParticle1;
                    }
                    tmpDatFileName = aJobTaskRecordList.get(i)
                            .result_CN_PathName()
                            + FILESEPARATOR
                            + tmpParticlePair
                            + ".dat";
                    try (BufferedWriter tmpBW = new BufferedWriter(
                                new FileWriter(tmpDatFileName, true))) {
                        // Coordination number mean
                        tmpBW.append("CNmean(");
                        tmpBW.append(tmpParticle1);
                        tmpBW.append("/");
                        tmpBW.append(tmpParticle2);
                        tmpBW.append(") = ");
                        tmpBW.append(String.format("%.2f", 
                                tmpCNMeans[tmpJobIndex]));
                        tmpBW.append("    ");
                        tmpBW.append(LINESEPARATOR);

                        // Standard deviation
                        tmpBW.append("Standard deviation(");
                        tmpBW.append(tmpParticle1);
                        tmpBW.append("/");
                        tmpBW.append(tmpParticle2);
                        tmpBW.append(") = ");
                        tmpBW.append(String.format("%.2f", 
                                tmpStdDeviation[tmpJobIndex]));
                        tmpBW.append("    ");
                        tmpBW.append(LINESEPARATOR);

                        // CN Min
                        tmpBW.append("Min(");
                        tmpBW.append(tmpParticle1);
                        tmpBW.append("/");
                        tmpBW.append(tmpParticle2);
                        tmpBW.append(") = ");
                        tmpBW.append(Integer.toString(tmpCNMin[tmpJobIndex]));
                        tmpBW.append("    ");
                        tmpBW.append(LINESEPARATOR);

                        // CN Max
                        tmpBW.append("Max(");
                        tmpBW.append(tmpParticle1);
                        tmpBW.append("/");
                        tmpBW.append(tmpParticle2);
                        tmpBW.append(") = ");
                        tmpBW.append(Integer.toString(tmpCNMax[tmpJobIndex]));
                        tmpBW.append("    ");
                        tmpBW.append(LINESEPARATOR);

                        // Box length
                        tmpBW.append("BoxLength [" + ANGSTROM + "] = ");
                        tmpBW.append(String.format("%.4f", 
                                tmpBoxLengths[tmpJobIndex]));
                        tmpBW.append(LINESEPARATOR);
                        tmpBW.close();
                        tmpJobIndex++;
                    } catch(IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                        "IOException during Writing .dat file.", ex);
                    }
                }
            }

            // </editor-fold>

            // <editor-fold defaultstate="collapsed" desc="Clear scratch directory">
            try (Stream<Path> tmpWalk = Files.walk(Paths
                    .get(scratchDirectory))) {
                tmpWalk.sorted(Comparator.reverseOrder())
                        .map(Path::toFile)
                        .forEach(File::delete);
            } catch(IOException ex) {
                LOGGER.log(Level.SEVERE, 
                        "IOException during deleting scratch directory.", ex);
            }
            
            // </editor-fold>
        }
    }
    
    /**
     * Export parameter set for MFSim
     * @param aJobTaskRecords: Job task records
     * @param aTitle: Title name
     * @param aTitleAbbreviation: Title abbreviation 
     */
    private static void exportParticleSetForMFSim(
            ArrayList<JobTaskRecord> aJobTaskRecords,
            String aTitle,
            String aTitleAbbreviation) {
        
        int tmpJobLength;
        int tmpOutputIteration;
        // 0: CN = 1, Wgt = true, Opt = true, Rgd = false
        // 1: CN = 1, Wgt = false, Opt = true, Rgd = false
        // 2: CN = 1, Wgt = false, Opt = true, Rgd = true
        // 3: CN = 1, Wgt = false, Opt = false
        // 4: CN != 1, Wgt = true, Opt = true, Rgd = false
        // 5: CN != 1, Wgt = false, Opt = true, Rgd = false
        // 6: CN != 1, Wgt = false, Opt = true, Rgd = true
        // 7: CN != 1, Wgt = false, Opt = false
        String tmpParticle;
        String tmpResultsDirectory;
        
        tmpOutputIteration = 0;
        tmpJobLength = aJobTaskRecords.size();
        ArrayList<String> tmpParticleNames = new ArrayList<>(tmpJobLength);
        tmpResultsDirectory = resultDirectory 
                + FILESEPARATOR
                + "IE"
                + FILESEPARATOR
                + forcefield_IE;
        
        // Get particle names
        for (int i = 0; i < tmpJobLength; i++) {
            tmpParticle = aJobTaskRecords.get(i).particleName1();
            if (!tmpParticleNames.contains(tmpParticle)) {
                tmpParticleNames.add(tmpParticle);
            }
        }
        
        Collections.sort(tmpParticleNames);

        // particle description
        List<String> tmpParticleDescriptionString = Arrays.asList(
                "# Particle", 
                "Name", 
                "Mass[DPD]", 
                "Charge", 
                "Mass[g/mol]", 
                "Volume[A^3]", 
                "Graphics-Radius", 
                "Standard-Color");
        ArrayList<ArrayList<String>> tmpParticleDescriptionsList = 
                new ArrayList<>(tmpParticleNames.size());
        
        for(String tmpParticleName : tmpParticleNames){
            ArrayList<String> tmpParticleDescriptions = new ArrayList<>(8);
            tmpParticleDescriptions.add(tmpParticleName); //particle
            tmpParticleDescriptions.add(tmpParticleName); //full name TODO
            tmpParticleDescriptions.add(Integer.toString(MASS_DPD)); //massDPD
            tmpParticleDescriptions.add(Integer.toString(DEFAULT_CHARGE)); // charge
            try {
                tmpParticleDescriptions.add(String.format("%.4f", 
                        MIPETUTIL.getAtomicMass(smiles.get(
                                tmpParticleName), true))); // mass [g/mol]
            } catch (Exception ex) {
                tmpParticleDescriptions.add("-1");
                LOGGER.log(Level.SEVERE, ex.toString());
            } 
            tmpParticleDescriptions.add(String.format("%.4f", 
                    MIPETUTIL.getVdwVolume(smiles.get(tmpParticleName)) 
                            * WATERVOLUMERATIO)); // volume [A³]
            tmpParticleDescriptions.add(Double.toString(GRAPHICS_RADIUS)); // graphics-radius
            tmpParticleDescriptions.add(DEFAULT_COLOR); // color
            tmpParticleDescriptionsList.add(tmpParticleDescriptions);
        }
        
        // particle interactions
        int tmpEnergyListLength;
        int tmpCNListLength;
        double tmpChiNumerator;
        double tmpE12;
        double tmpE11;
        double tmpE22;
        double tmpZ11;
        double tmpZ22;
        double tmpZ12;
        double tmpZ21;
        double tmpAij;
        Double tmpEnergy;
        Double tmpCN;
        String tmpParticleName1;
        String tmpParticleName2;
        String tmpParticlePair;
        String tmpFileName;
        Set<String> tmpKeySet;
        HashMap<String, Double> tmpAijMap;
        HashMap<String, Double> tmpAijMap1;
        
        tmpEnergyListLength = energyList.size();
        tmpCNListLength = cnList.size();
        tmpAijMap = new HashMap<>(tmpJobLength);
        tmpAijMap1 = new HashMap<>(tmpJobLength); // for CN =1
        HashMap<String, Double> tmpEnergieMap = 
                new HashMap<>(tmpEnergyListLength);
        HashMap<String, Double> tmpCNMap = new HashMap<>(tmpCNListLength);
        tmpKeySet = new HashSet<>();
        tmpEnergy = 0.;
        
        while (tmpOutputIteration <= 7) {
            // Read energy data
            for (int i = 0; i < tmpEnergyListLength; i++) {
                tmpParticleName1 = energyList.get(i).particleName1();
                tmpParticleName2 = energyList.get(i).particleName2();
                tmpParticlePair = tmpParticleName1 + "_" + tmpParticleName2;
                switch (tmpOutputIteration) {
                    case 0, 4 -> tmpEnergy = energyList.get(i)
                            .wgt_Emin();
                    case 1, 5 -> tmpEnergy = energyList.get(i)
                            .wgt0_Opt_Emin();
                    case 2, 6 -> tmpEnergy = energyList.get(i)
                            .wgt0_Rgd_Emin();
                    case 3, 7 -> tmpEnergy = energyList.get(i)
                            .wgt0_Opt0_Emin();
                }
                tmpEnergieMap.put(tmpParticlePair, tmpEnergy);
            }

            // Read coordination number
            for (int i = 0; i < tmpCNListLength; i++) {
                tmpParticleName1 = cnList.get(i).particleName1();
                tmpParticleName2 = cnList.get(i).particleName2();
                tmpParticlePair = tmpParticleName1 + "_" + tmpParticleName2;
                tmpCN = cnList.get(i).cnValue();
                tmpCNMap.put(tmpParticlePair, tmpCN);
            }
            
            // Calculate aij parameters (see DPD theory pdf)
            tmpAij = 0.;

            for(int i = 0; i < tmpJobLength; i++) {
                tmpParticleName1 = aJobTaskRecords.get(i).particleName1();
                tmpParticleName2 = aJobTaskRecords.get(i).particleName2();
                if (!aJobTaskRecords.get(i).isReverse()) {
                    tmpE12 = tmpEnergieMap.get(tmpParticleName1 + "_" 
                        + tmpParticleName2);
                    tmpE11 = tmpEnergieMap.get(tmpParticleName1 + "_" 
                            + tmpParticleName1);
                    tmpE22 = tmpEnergieMap.get(tmpParticleName2 + "_" 
                            + tmpParticleName2);
                    if (!forcefield_CN.isEmpty()) {
                        tmpZ11 = tmpCNMap.get(tmpParticleName1 + "_" 
                                + tmpParticleName1);
                        tmpZ22 = tmpCNMap.get(tmpParticleName2 + "_" 
                                + tmpParticleName2);
                        tmpZ12 = tmpCNMap.get(tmpParticleName1 + "_" 
                                + tmpParticleName2);
                        tmpZ21 = tmpCNMap.get(tmpParticleName2 + "_" 
                                + tmpParticleName1);
                        tmpChiNumerator =
                                tmpZ12 * tmpE12 +
                                tmpZ21 * tmpE12 -
                                tmpZ11 * tmpE11 -
                                tmpZ22 * tmpE22;
                        tmpAij = (double)temperature / 12 + 1.7483 
                                * tmpChiNumerator;
                        tmpAijMap.put(tmpParticleName1 + "_" + tmpParticleName2, 
                                tmpAij);
                    }

                    // Calculation for CN = 1
                    tmpChiNumerator = tmpE12 + tmpE12 - tmpE11 - tmpE22;
                    tmpAij = (double)temperature/12 + 1.7483 * tmpChiNumerator;
                    tmpAijMap1.put(tmpParticleName1 + "_" + tmpParticleName2, 
                            tmpAij);
                }
            }
            
            /* Write file */
            tmpFileName = tmpResultsDirectory 
                + FILESEPARATOR 
                + aTitleAbbreviation + "_"
                +"EijFraction_" 
                + boltzmannFraction
                + "_catchRadius_" 
                + catchRadius;
            switch (tmpOutputIteration) {
                case 0 -> tmpFileName += "_CN1_Wgt.txt";
                case 1 -> tmpFileName += "_CN1_Wgt0_Opt.txt";
                case 2 -> tmpFileName += "_CN1_Wgt0_Rgd.txt";
                case 3 -> tmpFileName += "_CN1_Wgt0_Opt0.txt";
                case 4 -> tmpFileName += "_CN_Wgt.txt";
                case 5 -> tmpFileName += "_CN_Wgt0_Opt.txt";
                case 6 -> tmpFileName += "_CN_Wgt0_Rgd.txt";
                case 7 -> tmpFileName += "_CN_Wgt0_Opt0.txt";
            }
            if (!forcefield_CN.isEmpty() ||
                    (forcefield_CN.isEmpty() && tmpOutputIteration <= 3)) {
                try (BufferedWriter tmpBW = Files.newBufferedWriter(
                        Paths.get(tmpFileName))) {
                    tmpBW.append("# Particle set for MFSim created by MIPET\n");
                    tmpBW.append("# Force Field for energy calculation: ");
                    tmpBW.append(forcefield_IE);
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Force Field for coordination number calculation: ");
                    if (forcefield_CN.isEmpty()) {
                        tmpBW.append("Not calculated");
                    } else {
                        tmpBW.append(forcefield_CN);
                    }
                    tmpBW.append(forcefield_CN);
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Water model: " + watermodel);
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# CPU cores: ");
                    tmpBW.append(Integer.toString(cpuCoreNumber));
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Temperature: ");
                    tmpBW.append(Double.toString(temperature));
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Sphere nodes calculated with Fibonacci algorithm: ");
                    tmpBW.append(Boolean.toString(isFibonacciSphereAlgorithm));
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Sphere node number 1: ");
                    tmpBW.append(Integer.toString(sphereNodeNumber1));
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Sphere node number 2: ");
                    tmpBW.append(Integer.toString(sphereNodeNumber2));
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Sphere node number 3: ");
                    tmpBW.append(Integer.toString(sphereNodeNumber3));
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Sphere node number 4: ");
                    tmpBW.append(Integer.toString(sphereNodeNumber4));
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Rotation number 1: ");
                    tmpBW.append(Integer.toString(rotationNumber1));
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Rotation number 2: ");
                    tmpBW.append(Integer.toString(rotationNumber2));
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Rotation number 3: ");
                    tmpBW.append(Integer.toString(rotationNumber3));
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Rotation number 4: ");
                    tmpBW.append(Integer.toString(rotationNumber4));
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Conformational analysis: ");
                    tmpBW.append(Boolean.toString(isConformationalAnalysis));
                    tmpBW.append(LINESEPARATOR);
                    switch (tmpOutputIteration) {
                        case 0, 4 -> {
                            tmpBW.append("# Boltzmann averaging: Yes");
                            tmpBW.append(LINESEPARATOR);}
                        case 1, 2, 3, 5, 6, 7 -> {
                            tmpBW.append("# Boltzmann averaging: No");
                            tmpBW.append(LINESEPARATOR);}
                    }
                    switch (tmpOutputIteration) {
                        case 1, 6 -> {
                            tmpBW.append("# Optimize sampled E(min) configuration: Yes");
                            tmpBW.append(LINESEPARATOR);}
                        case 0, 2, 3, 4, 5, 7 -> {
                            tmpBW.append("# Optimize sampled E(min) configuration: No");
                            tmpBW.append(LINESEPARATOR);}
                    }
                    switch (tmpOutputIteration) {
                        case 2, 6 -> {
                            tmpBW.append("# Tinker's 'optrigid' used: Yes");        
                            tmpBW.append(LINESEPARATOR);}
                    }
                    switch (tmpOutputIteration) {
                        case 0, 1, 2, 3 -> {
                            tmpBW.append("# CN = 1 for all particle pairs");
                            tmpBW.append(LINESEPARATOR);}
                    }
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("[Title]\n" + aTitle + "\n[/Title]\n\n");
                    tmpBW.append("[Version]\n" 
                            + VERSION_NUMBER 
                            + "\n[/Version]\n\n");

                    // Particle description
                    tmpBW.append("[Particle Description]");
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("""
                                 # Columns:
                                 # 1. Particle (abbreviation)
                                 # 2. Molecule name of particle
                                 # 3. Mass[DPD]
                                 # 4. Charge
                                 # 5. Mass[g/mol]
                                 # 6. Volume[A^3]
                                 # 7. Graphics-Radius 
                                 # 8. Standard-Color
                                 """);
                    tmpBW.append(LINESEPARATOR);

                    for(String tmpHeader : tmpParticleDescriptionString){
                        tmpBW.append(tmpHeader + "\t");
                    }

                    for(List<String> tmpParticleDescList : 
                            tmpParticleDescriptionsList){
                        tmpBW.append(LINESEPARATOR);

                        for(String tmpDescription : tmpParticleDescList){
                            tmpBW.append(tmpDescription);
                            tmpBW.append(" ");
                        }

                    }

                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("[/Particle Description]");
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append(LINESEPARATOR);

                     // Particle interactions
                    tmpBW.append("[Particle interactions]");
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Repulsion parameters a(ij) for particle pairs for different temperatures (in K)");
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("Pair");
                    tmpBW.append(" " + Integer.toString((int)temperature));
                    switch (tmpOutputIteration) {
                        case 0, 1, 2, 3 -> {
                            tmpKeySet = tmpAijMap1.keySet();}
                        case 4, 5, 6, 7 -> {
                            tmpKeySet = tmpAijMap.keySet();}
                    }

                    for(String tmpKey : tmpKeySet) {
                        tmpBW.append(LINESEPARATOR);
                        tmpBW.append(tmpKey);
                        switch (tmpOutputIteration) {
                            case 0, 1, 2, 3 -> {
                                tmpAij = tmpAijMap1.get(tmpKey);}
                            case 4, 5, 6, 7 -> {
                                tmpAij = tmpAijMap.get(tmpKey);}
                        }
                        tmpBW.append(" " +  String.format("%.2f",tmpAij));
                    }

                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("[/Particle interactions]");
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append(LINESEPARATOR);

                    // Coordination numbers
                    switch (tmpOutputIteration) {
                        case 4, 5, 6, 7 -> {
                            tmpBW.append("[Coordination numbers]");
                            tmpCNListLength = cnList.size();

                            for (int j = 0; j < tmpCNListLength; j++) {
                                tmpBW.append(LINESEPARATOR);
                                tmpBW.append(cnList.get(j).particleName1() + "_"
                                        + cnList.get(j).particleName2());
                                tmpBW.append(String.format(" %.2f", 
                                        cnList.get(j).cnValue()));
                            }

                            tmpBW.append(LINESEPARATOR);
                            tmpBW.append("[/Coordination numbers]");
                            tmpBW.append(LINESEPARATOR);}
                    }

                    // SMILES
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("[SMILES]");
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("# Particle and corresponding SMILES of fragment molecule");
                    tmpKeySet = smiles.keySet();

                    for(String tmpKey : tmpKeySet){
                        if (tmpParticleNames.contains(tmpKey)) {
                            tmpBW.append(LINESEPARATOR);
                            tmpBW.append(tmpKey); 
                            tmpBW.append(" ");
                            tmpBW.append(smiles.get(tmpKey));
                        }
                    }

                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append("[/SMILES]");
                } catch (IOException anException) {
                    LOGGER.log(Level.SEVERE, anException.toString());
                }
            }
            tmpOutputIteration++;
        }
            
    }
}
// </editor-fold>