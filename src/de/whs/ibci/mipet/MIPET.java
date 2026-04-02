/**
 * MIPET - Mesoscopic Interaction Parameter Estimation with Tinker
 * Copyright (C) 2026  Achim Zielesny (achim.zielesny@googlemail.com)
 * <p>
 * Source code is available at <https://github.com/zielesny/MIPET>
 * <p>
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
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
import java.nio.file.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.time.Duration;
import java.time.Instant;
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
import java.util.StringTokenizer;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
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
    private static final String FILESEPARATOR = FileSystems.getDefault().getSeparator();
    
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
     * Particle names
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
     * Partial charge correction factor
     */
    private static double[] chargeCorr;
    
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
    private static String forceField_IE;
    
    /**
     * Force field names for calculation of coordination number
     */
    private static String forceField_CN;
    
    /**
     * Water model name
     */
    private static String waterModel;
    
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
     * Dielectric constant
     */
    private static double dielectricConstant;

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
     * Directory name of particle data
     */
    private static String moleculeDirectory;
    
    /**
     * Directory name of optimized xyz data
     */
    private static String optXYZDirectory;
    
    /**
     * Directory name of optimal distance of two particles
     */
    private static String optDistDirectory;
    
    /**
     * Directory name of output data
     */
    private static String resultDirectory;
    
    /**
     * BufferedWriter object for global log file
     */
    private static BufferedWriter BFGblLog;

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
     * Parameter title
     */
    private static String parameterSetTitle;
    
    /**
     * Parameter title abbreviation
     */
    private static String parameterSetTitleAbr;
    
    /**
     * Flag for write warmup log file
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

    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     */    
    public MIPET() {
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Main method
     */
    public static void main() {
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

        //</editor-fold>
        
        // <editor-fold defaultstate="collapsed" desc="Initialize and read .job file">
        Locale.setDefault(Locale.ENGLISH);
        Instant startTime = Instant.now();
        System.out.println("Initializing...");
        initialize();
        System.out.println("Reading job file...");
        readJobFile();
        
        //</editor-fold>
        
        // <editor-fold defaultstate="collapsed" desc="Set original .key file">
        /*
         * Standard .key file content
         */
        String keyFileStringOrigin;
        if (forceField_IE.equals("OPLSAALIGPARGEN")) {
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
                + forceField_IE.toLowerCase() + ".prm\""
                + LINESEPARATOR;
        }
        
        // </editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Create log file">
        String tmpGblLogDirName = resultDirectory
                + FILESEPARATOR
                + "IE"
                + FILESEPARATOR
                + forceField_IE;
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
                + forceField_IE
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
            BFGblLog.append("CPU cores: ")
                    .append(String.valueOf(cpuCoreNumber))
                    .append(LINESEPARATOR);
            BFGblLog.append("Warm up steps for dynamic simulation: ")
                    .append(String.valueOf(warmUpStepNumber))
                    .append(LINESEPARATOR);
            BFGblLog.append("Steps for dynamic simulation: ")
                    .append(String.valueOf(stepNumber))
                    .append(LINESEPARATOR);
            BFGblLog.flush();
        } catch(IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during writing in log file.", ex);
        }
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Prepare input">
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
                    + forceField_IE
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
                    + forceField_CN
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
                    if (!forceField_CN.isEmpty()) {
                        Files.createDirectories(Paths.get(tmpCNResultDirName));
                    }
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during creating CN result directory.", ex);
                }
            }
            tmpJobTaskRecordList.add(new JobTaskRecord(
                    forceField_IE,
                    forceField_CN,
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
                        forceField_IE,
                        forceField_CN,
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
        readXyz(true, forceField_IE);
        readPrm(forceField_IE);
        if (forceField_IE.equals("OPLSAALIGPARGEN")) {
            makeMoleculeRecord();
        }
        
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
            readXyz(false, forceField_IE);
        } else {
            System.out.println("Conformational Analysis skipped.");
            tmpForcefield = tmpJobTaskRecordList.getFirst()
                    .forcefield_IE_Name();
            
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
                tmpSourceName = moleculeDirectory
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
        int tmpN_atom1;
        int tmpPrmID1;
        int tmpPrmID2;
        long tmpEnergyCalcTime;
        double tmpGlbEminDist;
        double tmpGlbEmin;
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
        LinkedList<Distance_EnergyRecord> tmpDistEminRecords;
        List<double[][][]> tmpRotCoords;
        
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
        tmpForcefield = forceField_IE;
        System.out.println("Calculating intermolecular energy...");
        
        //</editor-fold>
        
        while (!tmpIsExitCondition) {
            
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
                MIPETUtility.updateStatus("Calculating " + tmpParticlePair);
            
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
                tmpRotCoords = RotationUtil.getRotationsCoords(
                        sphereNodeNumber1, 
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
                tmpDistanceList = new LinkedList<>();
                tmpDistEminRecords = new LinkedList<>();
                tmpPrmID1 = particleNames.indexOf(tmpParticleName1);
                if (tmpIsSameParticle) {
                    tmpPrmID2 = tmpPrmID1;
                } else {
                    tmpPrmID2 = particleNames.indexOf(tmpParticleName2);
                }

                /* Write .key file */
                Path keyPath = Paths.get(scratchDirectory,
                        tmpParticlePair + ".key");
                StringBuilder keyContent = new StringBuilder(
                        keyFileStringOrigin)
                        .append("DIELECTRIC\t")
                        .append(dielectricConstant)
                        .append(LINESEPARATOR)
                        .append("NONBONDTERM ONLY")
                        .append(LINESEPARATOR)
                        .append("OPENMP-THREADS 1")
                        .append(LINESEPARATOR);
                
                if ("OPLSAALIGPARGEN".equals(tmpForcefield)) {
                    keyContent.append(prmContent1[tmpPrmID1]);
                    if (!tmpIsSameParticle) {
                        keyContent.append(prmContent2[tmpPrmID2]);
                    }
                }
                MIPETUTIL.writeKeyFile(keyPath, keyContent.toString());

                for (int i = 0; i < tmpDistSize; i++) {
                    tmpDistanceCandidate = lowerBoundary + i * prescanStepSize;
                    tmpDistanceList.add(tmpDistanceCandidate);
                    if (tmpDistanceCandidate > upperBoundary) {
                        tmpDistanceList.add(tmpDistanceCandidate);
                        break;
                    }
                }

                tmpAllDistances = new LinkedList<>(tmpDistanceList);
                double[] tmpDistances = tmpDistanceList.stream()
                        .mapToDouble(d -> d).toArray();
                
                //</editor-fold>
                
                if (tmpParticlePair.equals("Na+_Na+") || 
                        tmpParticlePair.equals("Ac-_Ac-")) {
                    // Check whether particle pair is Na+_Na+ or Ac-_Ac-
                    energyList.add(new ResultEnergyRecord(
                            tmpParticleName1, 
                            tmpParticleName2, 
                            0,
                            0,
                            0,
                            0));
                } else {
                    tmpEnergyRecords = new EnergyRecord[5];
                    tmpEnergyRecords[0] = getInterMolecularEnergy(
                            tmpParticlePair,
                            tmpDistances, 
                            tmpTinkerXYZ1, 
                            tmpTinkerXYZ2, 
                            tmpXyzRotData1, 
                            tmpXyzRotData2,
                            1E10);
                    tmpGlbEminDist = tmpEnergyRecords[0].eqDistance();
                    tmpGlbEmin = tmpEnergyRecords[0].Emin();

                    for (int i = 0; i < tmpDistances.length; i++) {
                        tmpDistEminRecords.add(new Distance_EnergyRecord(
                                tmpEnergyRecords[0].distances()[i], 
                                tmpEnergyRecords[0].Emins()[i],
                                tmpEnergyRecords[0].wgtEmins()[i]));
                    }

                    //<editor-fold defaultstate="collapsed" desc="Calculate rotated coordinates">
                    if (sphereNodeNumber2 != sphereNodeNumber1 
                            || rotationNumber2 != rotationNumber1) {
                        tmpRotCoords = RotationUtil.getRotationsCoords(
                                sphereNodeNumber2, 
                                rotationNumber2,
                                tmpXyzData1,
                                tmpXyzData2,
                                isFibonacciSphereAlgorithm);
                        tmpXyzRotData1 = tmpRotCoords.get(0);
                        tmpXyzRotData2 = tmpRotCoords.get(1);
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Precise scan">
                    tmpDistanceList = new LinkedList<>();
                    tmpDistSize = 9;

                    for (int i = 0; i < tmpDistSize; i++) {
                        tmpDistanceCandidate = (10 * tmpGlbEminDist - 4 + i) 
                                / 10;
                        tmpDistanceList.add(tmpDistanceCandidate);
                    }

                    tmpAllDistances.addAll(tmpDistanceList);
                    tmpDistances = tmpDistanceList.stream()
                            .mapToDouble(d -> d)
                            .toArray();
                    tmpEnergyRecords[1] = getInterMolecularEnergy(
                            tmpParticlePair,
                            tmpDistances, 
                            tmpTinkerXYZ1, 
                            tmpTinkerXYZ2, 
                            tmpXyzRotData1, 
                            tmpXyzRotData2,
                            tmpGlbEmin);
                    if (tmpEnergyRecords[1].Emin() < tmpGlbEmin) {
                        tmpGlbEminDist = tmpEnergyRecords[1].eqDistance();
                        tmpGlbEmin = tmpEnergyRecords[1].Emin();
                    }

                    for (int i = 0; i < tmpDistances.length; i++) {
                        tmpDistEminRecords.add(new Distance_EnergyRecord(
                                tmpEnergyRecords[1].distances()[i], 
                                tmpEnergyRecords[1].Emins()[i],
                                tmpEnergyRecords[1].wgtEmins()[i]));
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Calculate rotated coordinates">
                    if (sphereNodeNumber3 != sphereNodeNumber2 
                            || rotationNumber3 != rotationNumber2) {
                        tmpRotCoords = RotationUtil.getRotationsCoords(
                                sphereNodeNumber3, 
                                rotationNumber3,
                                tmpXyzData1,
                                tmpXyzData2,
                                isFibonacciSphereAlgorithm);
                        tmpXyzRotData1 = tmpRotCoords.get(0);
                        tmpXyzRotData2 = tmpRotCoords.get(1);
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="More precise scan">
                    tmpDistanceList = new LinkedList<>();
                    tmpDistSize = 19;

                    for (int i = 0; i < tmpDistSize; i++) {
                        tmpDistanceCandidate = (100 * tmpGlbEminDist - 9 + i) 
                                / 100;
                        tmpDistanceList.add(tmpDistanceCandidate);
                    }

                    tmpAllDistances.addAll(tmpDistanceList);
                    tmpDistances = tmpDistanceList.stream()
                            .mapToDouble(d -> d)
                            .toArray();
                    tmpEnergyRecords[2] = getInterMolecularEnergy(
                            tmpParticlePair,
                            tmpDistances, 
                            tmpTinkerXYZ1, 
                            tmpTinkerXYZ2, 
                            tmpXyzRotData1, 
                            tmpXyzRotData2,
                            tmpGlbEmin);
                    if (tmpEnergyRecords[2].Emin() < tmpGlbEmin) {
                        tmpGlbEminDist = tmpEnergyRecords[2].eqDistance();
                        tmpGlbEmin = tmpEnergyRecords[2].Emin();
                    }
                    
                    for (int i = 0; i < tmpDistances.length; i++) {
                        tmpDistEminRecords.add(new Distance_EnergyRecord(
                                tmpEnergyRecords[2].distances()[i], 
                                tmpEnergyRecords[2].Emins()[i],
                                tmpEnergyRecords[2].wgtEmins()[i]));
                    }
                    
                    tmpAllDistances.add(tmpEnergyRecords[2].eqDistance());

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Calculate rotated coordinates">
                    if (sphereNodeNumber4 != sphereNodeNumber3 
                            || rotationNumber4 != rotationNumber3) {
                        tmpRotCoords = RotationUtil.getRotationsCoords(
                                sphereNodeNumber4, 
                                rotationNumber4,
                                tmpXyzData1,
                                tmpXyzData2,
                                isFibonacciSphereAlgorithm);
                        tmpXyzRotData1 = tmpRotCoords.get(0);
                        tmpXyzRotData2 = tmpRotCoords.get(1);
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Scan at equilibrium distance">
                    tmpDistances = new double[1];
                    tmpDistances[0] = tmpGlbEminDist; 
                    tmpEnergyRecords[3] = getInterMolecularEnergy(
                            tmpParticlePair,
                            tmpDistances, 
                            tmpTinkerXYZ1, 
                            tmpTinkerXYZ2, 
                            tmpXyzRotData1, 
                            tmpXyzRotData2,
                            tmpGlbEmin);
                    if (tmpEnergyRecords[3].Emin() < tmpGlbEmin) {
                        tmpGlbEminDist = tmpEnergyRecords[3].eqDistance();
                        tmpGlbEmin = tmpEnergyRecords[3].Emin();
                    }

                    for (int i = 0; i < tmpDistances.length; i++) {
                        tmpDistEminRecords.add(new Distance_EnergyRecord(
                                tmpEnergyRecords[3].distances()[i], 
                                tmpEnergyRecords[3].Emins()[i],
                                tmpEnergyRecords[3].wgtEmins()[i]));
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Find distance by min. wgtEmin">
                    int tmpDistEminRecordsSize;
                    double tmpWgtEminDist;
                    double tmpGlbWgtEmin;
                    double tmpGlbWgtEminCand;

                    tmpDistEminRecordsSize = tmpDistEminRecords.size();
                    tmpWgtEminDist = tmpDistEminRecords.getFirst().distance();
                    tmpGlbWgtEmin = tmpDistEminRecords.getFirst().wgtEmin();
                    
                    for (int i = 1; i < tmpDistEminRecordsSize; i++) {
                        tmpGlbWgtEminCand = tmpDistEminRecords.get(i).wgtEmin();
                        if (tmpGlbWgtEminCand < tmpGlbWgtEmin) {
                            tmpGlbWgtEmin = tmpGlbWgtEminCand;
                            tmpWgtEminDist = tmpDistEminRecords.get(i)
                                    .distance();
                        }
                    }
                    
                    // If tmpWgtEminDist is too far away from tmpGlbEminDist new calculations around tmpWgtEminDist will occur
                    // This is only for the output particle1_particle2_DistanceVsEaverage.jpg
                    if (Math.abs(tmpWgtEminDist - tmpGlbEminDist) > 0.1) {
                        tmpDistanceList = new LinkedList<>();
                        tmpDistSize = 19;

                        for (int i = 0; i < tmpDistSize; i++) {
                            tmpDistanceCandidate = (100 * tmpWgtEminDist - 9 
                                    + i) / 100;
                            if (!MIPETUTIL.contains(tmpAllDistances, 
                                    tmpDistanceCandidate)) {
                                tmpDistanceList.add(tmpDistanceCandidate);
                            }
                        }

                        tmpAllDistances.addAll(tmpDistanceList);
                        tmpDistances = tmpDistanceList.stream()
                                .mapToDouble(d -> d).toArray();
                        tmpEnergyRecords[4] = getInterMolecularEnergy(
                                tmpParticlePair,
                                tmpDistances, 
                                tmpTinkerXYZ1, 
                                tmpTinkerXYZ2, 
                                tmpXyzRotData1, 
                                tmpXyzRotData2,
                                tmpGlbEmin);
                        
                        for (int i = 0; i < tmpDistances.length; i++) {
                            tmpDistEminRecords.add(new Distance_EnergyRecord(
                                    tmpEnergyRecords[4].distances()[i], 
                                    tmpEnergyRecords[4].Emins()[i],
                                    tmpEnergyRecords[4].wgtEmins()[i]));
                            tmpGlbWgtEminCand = tmpEnergyRecords[4]
                                    .wgtEmins()[i];
                            if (tmpGlbWgtEmin > tmpGlbWgtEminCand) {
                                tmpGlbWgtEmin = tmpGlbWgtEminCand;
                            }
                        }
                        
                    }

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
                                .get(tmpDistanceIndices[i])
                                .distance();
                        tmpEnergySorted[i][1] = tmpDistEminRecords
                                .get(tmpDistanceIndices[i])
                                .Emin();
                        tmpEnergySorted[i][2] = tmpDistEminRecords
                                .get(tmpDistanceIndices[i])
                                .wgtEmin();
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Determining opt. Emin">
                    double tmpOptMinEnergy;
                    double tmpRgdMinEnergy;
                    File tmpOptFile;
                    
                    tmpOptMinEnergy = 0.0;
                    tmpRgdMinEnergy = 0.0;
                    
                    if (isTinkerOn) {
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
                        tmpKeyFileName = tmpParticlePair + ".key";
                        keyPath = Paths.get(scratchDirectory, tmpKeyFileName);
                        tmpKeyContent = tmpKeyFileString;
                        if (tmpForcefield.equals("OPLSAALIGPARGEN")) {
                            tmpKeyContent += prmContent1[tmpPrmID1];
                            if (!tmpIsSameParticle) {
                                tmpKeyContent += prmContent2[tmpPrmID2];
                            }
                        }
                        MIPETUTIL.writeKeyFile(keyPath, tmpKeyContent);

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
                            tmpN_atom1 =tmpTinkerXYZ1.getN_atom1();
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
                            if (forceField_IE.equals("OPLSAALIGPARGEN")
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
                                    tmpIndex = tmpN_atom1 + 1;
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
                            if (tmpProcess != null) {
                                tmpProcess.destroy();
                            }
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
                                            tmpOptMinEnergy = Double.parseDouble(tmpLine.substring(25, 50));
                                        } else {
                                            tmpRgdMinEnergy = Double.parseDouble(tmpLine.substring(25, 50));
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
                    }
                    
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
                        tmpBW.append("distance [" + ANGSTROM +
                                "]  Emin(Cmin,r) [kcal/mole] <E>(r) [kcal/mole]")
                                .append(LINESEPARATOR);

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
                    tmpChartUtil.createEnergyGraphics(
                            tmpEnergyDataPathName,
                            tmpEnergyGraphicsPrefix);

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Write log file">
                    try {
                        BFGblLog.append(tmpParticlePair);
                        BFGblLog.append(LINESEPARATOR);
                        BWParticleDat.append("Force field: ")
                                .append(tmpForcefield);
                        BWParticleDat.append(LINESEPARATOR);
                        BWParticleDat.append("Conformational analysis: "); 
                        if (isConformationalAnalysis) {
                            BWParticleDat.append("Yes");
                        } else {
                            BWParticleDat.append("No");
                        }
                        BWParticleDat.append(LINESEPARATOR);
                        BWParticleDat.append("equilibriumDistances [" + ANGSTROM + "] = "); 
                        BWParticleDat.append(decimal2.format(tmpGlbEminDist));
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

                    //<editor-fold defaultstate="collapsed" desc="Write ouput.0 file">
                    tmpOriginal = Paths.get(scratchDirectory, 
                            tmpParticlePair + ".0");
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
                    
                    //</editor-fold>
                    if (isTinkerOn) {
                        //<editor-fold defaultstate="collapsed" desc="Write output_opt.out">
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
                        
                        //</editor-fold>
                        //<editor-fold defaultstate="collapsed" desc="Write output_rgd.out">
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
                        
                        //</editor-fold>
                    }

                    //<editor-fold defaultstate="collapsed" desc="Generate output.xyz">
                    String tmpOutput0FileName = tmpIEResultDirName 
                            + FILESEPARATOR 
                            + "output.0";
                    String tmpOptFileName;
                    String tmpRgdFileName;
                    int tmpAtomNumber1 = tmpTinkerXYZ1.getN_atom1();
                    int tmpAtomNumber2 = tmpTinkerXYZ2.getN_atom1();
                    tmpElementList1 = tmpTinkerXYZ1.getElementList1();
                    tmpElementList2 = tmpTinkerXYZ2.getElementList1();
                    TinkerXYZ tmpTinkerXyz = 
                            new TinkerXYZ(tmpOutput0FileName, 1, 
                            tmpAtomNumber1, tmpAtomNumber2);
                    tmpFileName = tmpIEResultDirName 
                            + FILESEPARATOR 
                            + "output.xyz";
                    tmpTinkerXyz.writeToXyzFile(tmpFileName);
                    
                    //</editor-fold>
                    
                    // Generate .pdb file of output.0
                    TinkerToPdbConverter tmpTinkerToPdb = 
                            new TinkerToPdbConverter(); 
                    String tmpOutputFileName;
                    tmpOutputFileName = tmpOutput0FileName.substring(0, 
                            tmpOutput0FileName.length() - 2) + ".pdb";
                    tmpTinkerToPdb.convert(tmpOutput0FileName, 
                            tmpOutputFileName);
                    
                    if (isTinkerOn) {
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
                        tmpTinkerXyz = new TinkerXYZ(tmpOptFileName, 1, 
                                tmpAtomNumber1, 
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
                        tmpTinkerXyz = new TinkerXYZ(tmpRgdFileName, 1, 
                                tmpAtomNumber1, 
                                tmpAtomNumber2);
                        tmpTinkerXyz.setElementList1(tmpElementList1);
                        tmpTinkerXyz.setElementList2(tmpElementList2);
                        tmpFileName = tmpIEResultDirName 
                                + FILESEPARATOR 
                                + "output_rgd.xyz";
                        tmpTinkerXyz.writeToXyzFile(tmpFileName);

                        //<editor-fold defaultstate="collapsed" desc="Generate .pdb files">
                        MIPETUTIL.callXYZPDB(tinkerXYZPdb, 
                                tmpOptFileName, 
                                tmpKeyFileContent, 
                                tmpHasH2O);
                        MIPETUTIL.callXYZPDB(tinkerXYZPdb, 
                                tmpRgdFileName, 
                                tmpKeyFileContent, 
                                tmpHasH2O);
                        
                        //</editor-fold>
                    }

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
                            tmpBW.append(decimal4.format(tmpGlbEminDist));
                        } catch(IOException ex) {
                            LOGGER.log(Level.SEVERE, 
                                "IOException during writing file in OptDist directory.", 
                                ex);
                        }
                    }
                    tmpEnergyCalcTime  = (System.currentTimeMillis() 
                            - tmpEnergyCalcTime) / 1000;
                    try {
                        BWParticleDat.append("SphereNodeNumber1: ")
                                .append(Integer.toString(sphereNodeNumber1))
                                .append(LINESEPARATOR);
                        BWParticleDat.append("SphereNodeNumber2: ")
                                .append(Integer.toString(sphereNodeNumber2))
                                .append(LINESEPARATOR);
                        BWParticleDat.append("SphereNodeNumber3: ")
                                .append(Integer.toString(sphereNodeNumber3))
                                .append(LINESEPARATOR);
                        BWParticleDat.append("SphereNodeNumber4: ")
                                .append(Integer.toString(sphereNodeNumber4))
                                .append(LINESEPARATOR);
                        BWParticleDat.append("RotationNumber1: ")
                                .append(Integer.toString(rotationNumber1))
                                .append(LINESEPARATOR);
                        BWParticleDat.append("RotationNumber2: ")
                                .append(Integer.toString(rotationNumber2))
                                .append(LINESEPARATOR);
                        BWParticleDat.append("RotationNumber3: ")
                                .append(Integer.toString(rotationNumber3))
                                .append(LINESEPARATOR);
                        BWParticleDat.append("RotationNumber4: ")
                                .append(Integer.toString(rotationNumber4))
                                .append(LINESEPARATOR);
                        BWParticleDat.append("Temperature [K]: ");
                        BWParticleDat.append(Double.toString(temperature));
                        BWParticleDat.append(LINESEPARATOR);
                        BWParticleDat.append("Fraction of energy values used for the Boltzmann distribution: ");
                        BWParticleDat.append(Double.toString(boltzmannFraction));
                        BWParticleDat.append(LINESEPARATOR);
                        if (isOptEmin) {
                            BWParticleDat.append("Optimize sampled E(min) configuration: ")
                                    .append(String.valueOf(isOptEmin));
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
                        BWParticleDat.append("Time to calculate minimum intermolecular energy [s]: ")
                                .append(String.valueOf(tmpEnergyCalcTime))
                                .append(LINESEPARATOR);
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
                Distance vs minimum differential pair interaction energy data.

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
                + forceField_IE
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
        if (!forceField_CN.isEmpty()) {
            System.out.println("Calculating coordination numbers...");
            getCoordinationNumbers(tmpJobTaskRecordList);
        }
        
        //</editor-fold>
        //<editor-fold defaultstate="collapsed" desc="Export parameterset">
        System.out.print("\rExporting parameterset...");
        Instant endTime = Instant.now();
        Duration duration = Duration.between(startTime, endTime);
        float totalTimeSeconds = duration.toMillis() / 1000f;
        try {
            BFGblLog.append(LINESEPARATOR)
                    .append("Entire calculation Time: ")
                    .append(String.valueOf(totalTimeSeconds)).append(" s ");
            if (totalTimeSeconds > 3600) {
                BFGblLog.append("(")
                        .append(decimal2.format(totalTimeSeconds / 3600f))
                        .append(" h)");
            } else if (totalTimeSeconds > 60) {
                BFGblLog.append("(")
                        .append(decimal2.format(totalTimeSeconds / 60f))
                        .append(" min)");
            }       
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, "IOException during writing log file.", 
                    ex);
        } finally {
            if (BFGblLog != null) {
                try {
                    BFGblLog.close();
                } catch (IOException e) {
                    LOGGER.log(Level.SEVERE, "Failed to close log file.", e);
                }
            }
        }
        exportParticleSetForMFSim(tmpJobTaskRecordList, 
                parameterSetTitle, 
                parameterSetTitleAbr);
        System.out.println("\r\nReady.");
        
        //</editor-fold>
    }
    
    //</editor-fold>
    //</editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Initialize method
     */
    private static void initialize() {
        
        setParameters();
        
        // <editor-fold defaultstate="collapsed" desc="Check and create directories">
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
        
        // Checks if directory is empty
        if (!MIPETUTIL.isDirectoryEmpty(scratchDirectory)) {
            System.out.println(scratchDirectory + " is not empty.");
            System.exit(-1);
        }
        if (!MIPETUTIL.isDirectoryEmpty(resultDirectory)) {
            System.out.println(resultDirectory + " is not empty.");
            System.exit(-1);
        }
        if (!MIPETUTIL.isDirectoryEmpty(optXYZDirectory)) {
            System.out.println(optXYZDirectory + " is not empty.");
            System.exit(-1);
        }
        if (!MIPETUTIL.isDirectoryEmpty(optDistDirectory)) {
            System.out.println(optDistDirectory + " is not empty.");
            System.exit(-1);
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
        char tmpFirstChar;
        
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
                     // ignore comment line
                    }
                    case '*' -> {
                        String[] tmpString = tmpRestString.split("\\s+");
                        forceField_IE = tmpString[0];
                    }
                    case '$' -> {
                        String[] tmpString = tmpRestString.split("\\s+");
                        forceField_CN = tmpString[0];
                    }
                    case '-' -> {
                        if (!oldParticles.contains(tmpRestString)) {
                            oldParticles.add(tmpRestString);
                        }
                        if (!particleNames.contains(tmpRestString)) {
                            particleNames.add(tmpRestString);
                        }
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
        } catch (FileNotFoundException ex) {
            throw new IllegalArgumentException("No jobfile found.");
        } catch (IOException ex) {
            throw new IllegalArgumentException("IOException during"
                    + "reading job file.");
        }
        if (forceField_IE.equals("OPLSAALIGPARGEN")) {
            for (String particle : particleNames) {
                if (!smiles.containsKey(particle)) {
                    System.out.println(particle + " was not found in Smiles.dat");
                    System.exit(1);
                }

                // Check whether .xyz file exists
                Path xyzPath = Path.of(moleculeDirectory, 
                         FILESEPARATOR,
                         "OPLSAALIGPARGEN",
                         particle + ".xyz");
                if (!Files.exists(xyzPath)) {
                    System.out.println(particle + ".xyz was not found.");
                    System.exit(1);
                }
                
                // Check whether .prm file exists
                Path prmPath = Path.of(parameterDirectory, 
                         FILESEPARATOR, 
                         "OPLSAALIGPARGEN",
                         particle + ".prm");
                if (!Files.exists(prmPath)) {
                    System.out.println(particle + ".prm was not found.");
                    System.exit(1);
                }
            }
        }
    }
    
    /**
     * SetParameters method
     */
    private static void setParameters() {
        String tmpIsFibonacciSphereAlgorithm;
        String tmpIsOptEmin;
        String tmpTinkerDirectory;
        String tmpCalcDirectory;
        int tmpDynamicStepsPerCore;
        int tmpDynamicWarmupIteration;
        
        isTinker9 = MIPETUTIL.getResourceString("MIPET.Tinker9")
                .equalsIgnoreCase("true");
        isTinkerOn = MIPETUTIL.getResourceString("MIPET.TinkerOn")
                .equalsIgnoreCase("true");
        jobFileName = MIPETUTIL.getResourceString("MIPET.File.jobfile");
        cpuCoreNumber = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETCPUCoreNumber"));
        tmpIsFibonacciSphereAlgorithm = MIPETUTIL.getResourceString(
                "MIPETJavaUseFibonacciSphereAlgorithm");
        isFibonacciSphereAlgorithm = tmpIsFibonacciSphereAlgorithm
                .equals("true");
        tmpIsOptEmin = MIPETUTIL.getResourceString("MIPETOptEmin");
        isOptEmin = tmpIsOptEmin.equals("true");
        if (!isTinkerOn) {
            isConformationalAnalysis = false;
        } else {
            isConformationalAnalysis = MIPETUTIL
                .getResourceString("MIPETConformationalAnalysis")
                .equals("true");
        }
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
        /*
         * Directory of tinker files
         */
        tmpTinkerDirectory = MIPETUTIL
                .getResourceString("MIPET.Directory.tinker");
        if (OS.matches("Windows [0-9]*")) {
            tinkerAnalyze = tmpTinkerDirectory + "/analyze.exe";
            tinkerScan = tmpTinkerDirectory + "/scan.exe";
            tinkerXYZEdit = tmpTinkerDirectory + "/xyzedit.exe";
            tinkerMinimize = tmpTinkerDirectory + "/minimize.exe";
            tinkerDynamic = tmpTinkerDirectory + "/dynamic.exe";
            tinkerOptimize = tmpTinkerDirectory + "/optimize.exe";
            tinkerOptrigid = tmpTinkerDirectory + "/optrigid.exe";
            tinkerXYZPdb = tmpTinkerDirectory + "/xyzpdb.exe";
        } else if (OS.equals("Linux")) {
            if (isTinker9) {
                tinkerDynamic = tmpTinkerDirectory + "/tinker9";
            } else {
                tinkerDynamic = tmpTinkerDirectory + "/dynamic";
            }
            tinkerAnalyze = tmpTinkerDirectory + "/analyze";
            tinkerScan = tmpTinkerDirectory + "/scan";
            tinkerXYZEdit = tmpTinkerDirectory + "/xyzedit";
            tinkerMinimize = tmpTinkerDirectory + "/minimize";
            tinkerOptimize = tmpTinkerDirectory + "/optimize";
            tinkerOptrigid = tmpTinkerDirectory + "/optrigid";
            tinkerXYZPdb = tmpTinkerDirectory + "/xyzpdb";
        }
        scratchDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.scratch");
        smilesDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.smiles");
        /*
          Directory of calculation related datas
         */
        tmpCalcDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.calculation");
        resultDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.result");
        moleculeDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.source");
        optXYZDirectory = tmpCalcDirectory + "/OptXYZ";
        optDistDirectory = tmpCalcDirectory + "/OptDist";
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
        /*
          Step number per CPU core number
         */
        tmpDynamicStepsPerCore = stepNumber / cpuCoreNumber;
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
        // timeStep in fs and printInterval in ps
        /*
         * Warmup iteration timeStep in fs and printInterval in ps
         */
        tmpDynamicWarmupIteration = (int) (warmUpStepNumber * warmUpTimeStep
                * 1E-15 / (warmUpPrintInterval * 1E-12));
        // timeStep in fs and printInterval in ps
        nDynamicIteration = (int)(stepNumber * timeStep * 1E-15 /
                (printInterval * 1E-12));
        // timeStep in fs and printInterval in ps
        nDynamicIterationPerCore = 
                (int) Math.ceil((tmpDynamicStepsPerCore * timeStep * 1E-15) /
                (printInterval * 1E-12));
        // Write warmup log file
        isLogWarmup = MIPETUTIL.getResourceString("MIPET.Log.Warmup")
                .equalsIgnoreCase("true");
        // Write dynamic log file
        isLogDynamic = MIPETUTIL.getResourceString("MIPET.Log.Dynamic")
                .equalsIgnoreCase("true");
        // Write the RATTLE key word in the .key file
        isRATTLE = MIPETUTIL.getResourceString("MIPETRATTLE")
                .equalsIgnoreCase("true");
        // Write build simulation box log file
        isLogBuildBox = MIPETUTIL.getResourceString("MIPET.Log.Build")
                .equalsIgnoreCase("true");
        // Write minimize simulation box log file
        isLogMinimizeBox = MIPETUTIL.getResourceString("MIPET.Log.MinimizeBox")
                .equalsIgnoreCase("true");
        // Write soak simulation box log file
        isLogSoakBox = MIPETUTIL.getResourceString("MIPET.Log.Soak")
                .equalsIgnoreCase("true");
        // Parameter title
        parameterSetTitle = MIPETUTIL.getResourceString("MIPETParticleSetTitle");
        // Parameter title abbreviation
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
                tmpXyzName1 = moleculeDirectory
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
            
            // Change atom type number of 2. particle to avoid redundancy
            // Yoda condition protects against NullPointerException if aForcefield is null
            if ("OPLSAALIGPARGEN".equals(aForcefield) && 
                    xyzContent1[i] != null && !xyzContent1[i].isEmpty()) {
                StringBuilder tmpParsedXyz = new StringBuilder();
                String[] lines = xyzContent1[i].lines().toArray(String[]::new);

                // Append the first line (the header) if the file isn't empty
                if (lines.length > 0) {
                    tmpParsedXyz.append(lines[0]); 
                }

                for (int j = 1; j < lines.length; j++) {
                    String tmpCurrentLine = lines[j];
                    String[] tokens = tmpCurrentLine.trim().split("\\s+");

                    // Safety check: Does the line actually have at least 6 columns?
                    if (tokens.length >= 6) {
                        try {
                            // Safely parse the old atom type (column 6 / index 5) as an integer
                            int tmpOldType = Integer.parseInt(tokens[5]);
                            int tmpNewType = tmpOldType + 100;

                            String oldTypeStr = " " + tmpOldType + " ";
                            String newTypeStr = " " + tmpNewType + " ";

                            // Replace the text 
                            // (WARNING: As mentioned before, this might accidentally replace the atom ID if it matches the type!)
                            tmpCurrentLine = tmpCurrentLine
                                    .replace(oldTypeStr, newTypeStr);

                        } catch (NumberFormatException e) {
                            // Ignore gracefully if the 6th column is not a valid number (e.g., malformed lines)
                        }
                    }

                    tmpParsedXyz.append(LINESEPARATOR).append(tmpCurrentLine);
                }

                xyzContent2[i] = tmpParsedXyz.toString();

            } else {
                // Fallback if the forcefield doesn't match or the string is empty
                xyzContent2[i] = (xyzContent1[i] != null) ? xyzContent1[i] : "";
            }
            
            
            
            
            
            
            
            tmpXyz2.setLength(0);
            if (aForcefield.equals("OPLSAALIGPARGEN")) {
                tmpLines = xyzContent1[i].lines().toArray(String[]::new);
                tmpXyz2.append(tmpLines[0]);

                for (int j = 1; j < tmpLines.length; j++) {
                    tmpTokens = tmpLines[j].trim().split("\\s+");
                    tmpOldAtomType = " " + tmpTokens[5] + " ";
                    tmpNewAtomType = " "
                            + (Integer.parseInt(tmpTokens[5]) + 100) 
                            + " ";
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
        int tmpChargeCorrStartPos;
        int tmpChargeCorrEndPos;
        String tmpPrmName1;
        String tmpParticleName;
        String tmpSearch;
        
        tmpParticlesLength = particleNames.size();
        tmpSearch = "Partial charge correction factor:";
        prmContent1 = new String[tmpParticlesLength];
        prmContent2 = new String[tmpParticlesLength];
        chargeCorr = new double[tmpParticlesLength];
                
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
                
                // Read partial charge correction factor
                if (prmContent1[i].contains(tmpSearch)) {
                    tmpChargeCorrStartPos = prmContent1[i]
                            .indexOf(tmpSearch) + 33;
                    tmpChargeCorrEndPos = prmContent1[i]
                            .indexOf("\n", tmpChargeCorrStartPos);
                    chargeCorr[i] = Double.parseDouble(prmContent1[i]
                            .substring(tmpChargeCorrStartPos, 
                                    tmpChargeCorrEndPos));
                } else {
                    chargeCorr[i] = 1.;
                }
                
            }
            
            // Change atomtype number of 2. particle to avoid redundancy
            if (aForcefield.equals("OPLSAALIGPARGEN") &&
                    !prmContent1[i].isEmpty()) {
                prmContent2[i] = prmContent1[i].lines()
                    .dropWhile(line -> !line.contains("Atom Type Definitions"))
                    .map(line -> {
                        String[] tmpTokens1 = line.trim().split("\\s+");
                        if (tmpTokens1.length > 0) {
                            return switch (tmpTokens1[0]) {
                                case 
                                    "atom", 
                                    "vdw", 
                                    "charge", 
                                    "bond", 
                                    "angle", 
                                    "torsion", 
                                    "imptors" 
                                     -> MIPETUTIL.changeAtomType(tmpTokens1);
                                default 
                                     -> line;
                            };
                        }
                        return line;
                    })
                    .filter(line -> !line.isEmpty())
                    .map(line -> line + LINESEPARATOR) 
                    .collect(Collectors.joining());
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
        double[][] tmpEpsilons;
        double[][] tmpSigmas;
        double[][] tmpCharges;
        StringTokenizer tmpTokenizer;
        String[] tmpLines;
        String[] tmpWords;
        String[][] tmpElements;
        
        molecules = new LinkedList<>();
        tmpParticleSize = particleNames.size();
        tmpAtomNumber = new int[tmpParticleSize];
        tmpElements = new String[tmpParticleSize][];
        tmpAtomTypes = new int[tmpParticleSize][];
        tmpSigmas = new double[tmpParticleSize][];
        tmpEpsilons = new double[tmpParticleSize][];
        tmpCharges = new double[tmpParticleSize][];
        
        // Read elements and atomTypes
        for (int i = 0; i < tmpParticleSize; i++) {
            tmpLines = xyzContent1[i].lines().toArray(String[]::new);
            tmpWords = tmpLines[0].trim().split("\\s+");
            tmpAtomNumber[i] = Integer.parseInt(tmpWords[0]);
            tmpElements[i] = new String[tmpAtomNumber[i]];
            tmpAtomTypes[i] = new int[tmpAtomNumber[i]];
            
            for (int j = 1; j < tmpLines.length; j++) {
                tmpTokenizer = new StringTokenizer(tmpLines[j]);
                tmpTokenizer.nextToken();
                tmpElements[i][j - 1] = tmpTokenizer.nextToken();
                tmpTokenizer.nextToken();
                tmpTokenizer.nextToken();
                tmpTokenizer.nextToken();
                tmpAtomTypes[i][j-1] = Integer
                        .parseInt(tmpTokenizer.nextToken());
            }
            
        }
        
        // Read epsilons, sigmas and charges
        for (int i = 0; i < tmpParticleSize; i++) {
            tmpSigmaIndex = 0;
            tmpChargeIndex = 0;
            tmpLines = prmContent1[i].lines().toArray(String[]::new);
            tmpSigmas[i] = new double[tmpAtomNumber[i]];
            tmpEpsilons[i] = new double[tmpAtomNumber[i]];
            tmpCharges[i] = new double[tmpAtomNumber[i]];

            for (String tmpLine : tmpLines) {
                if (tmpLine.startsWith("vdw ")) {
                    tmpTokenizer = new StringTokenizer(tmpLine);
                    tmpTokenizer.nextToken();
                    tmpTokenizer.nextToken();
                    tmpSigmas[i][tmpSigmaIndex] =
                            Double.parseDouble(tmpTokenizer.nextToken());
                    tmpEpsilons[i][tmpSigmaIndex] =
                            Double.parseDouble(tmpTokenizer.nextToken());
                    tmpSigmaIndex++;
                } else if (tmpLine.startsWith("charge ")) {
                    tmpTokenizer = new StringTokenizer(tmpLine);
                    tmpTokenizer.nextToken();
                    tmpTokenizer.nextToken();
                    tmpCharges[i][tmpChargeIndex] = 
                            Double.parseDouble(tmpTokenizer.nextToken()) 
                            * chargeCorr[i];
                    tmpChargeIndex++;
                }
            }
            
            molecules.add(new MoleculeRecord(
                    particleNames.get(i),
                    tmpAtomNumber[i],
                    tmpElements[i],
                    tmpAtomTypes[i],
                    tmpEpsilons[i],
                    tmpSigmas[i],
                    tmpCharges[i]
            ));
        }
        
    }
    
    /**
     * scanParticle method
     *  Optimize, scan and optimize the input particle
     */
    private static void scanParticle() {
        String tmpForcefieldName;
        String tmpKeyContent;
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
        
        if (!forceField_IE.isEmpty()) {
            tmpForcefieldList.add(forceField_IE);
        }
        if (!forceField_CN.isEmpty()) {
            tmpForcefieldList.add(forceField_CN);
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
                    tmpSourceName = moleculeDirectory
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
                    Path keyPath = Paths.get(tmpOptXyzDirName, 
                            tmpParticle + ".key");
                    if ("OPLSAALIGPARGEN".equals(tmpForcefield)) {
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
                    MIPETUTIL.writeKeyFile(keyPath, tmpKeyContent);
                    
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
                    tmpTinkerXYZ.setForcefieldName(forceField_IE);
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

                            for (String s : tmpMinimumList) {
                                tmpEnergyValues.add(Double
                                        .valueOf(s.substring(58, 68)));
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

                            // Copy configuration data with the lowest energy to Particle.xyz
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
        
        if (forceField_CN.equals("OPLSAALIGPARGEN")) {
            readPrm(forceField_IE);
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
        
        //<editor-fold defaultstate="collapsed" desc="Calculate chunk size">
        boolean tmpIsFractionOne;
        int tmpAtomNumber;
        int tmpDistanceNumber;
        int tmpChunkNumber; // Number of chunks
        int tmpChunkSize; // Number of aRotData2 configurations
        double[][][] tmpRotData2;
        String tmpPath;
        String[] tmpCmdList;
        ArrayList<MIPETAnalyze> tmpTaskList;
        TinkerXYZ tmpTinkerXyz1;
        TinkerXYZ tmpTinkerXyz2;
        TinkerXYZ tmpTinkerXYZ;
        ExecutorService executor;
        
        String particleName1 = aParticlePair.split("_")[0];
        String particleName2 = aParticlePair.split("_")[1];
        tmpIsFractionOne = boltzmannFraction == 1.0;
        tmpDistanceNumber = aDistances.length;
        int totalElements = aRotData2.length;
        int maxChunkSize = 1000;
        int minChunks = (totalElements + maxChunkSize - 1) / maxChunkSize;
        tmpChunkNumber = ((minChunks + cpuCoreNumber - 1) / cpuCoreNumber) 
                * cpuCoreNumber;
        tmpChunkSize = (totalElements + tmpChunkNumber - 1) / tmpChunkNumber;
        
        //</editor-fold>

        //<editor-fold defaultstate="collapsed" desc="Calculate intermolecular energy using TINKER analyze">
        MoleculeRecord rec1 = null;
        MoleculeRecord rec2 = null;
        
        for (MoleculeRecord m : molecules) {
            if (m.name().equals(particleName1)) {
                rec1 = m;
            }
            if (m.name().equals(particleName2)) {
                rec2 = m;
            }
        }
        
        tmpTinkerXyz1 = aTinkerXYZ1.clone();
        tmpTinkerXyz2 = aTinkerXYZ2.clone();
        Flat3DArray sharedFlatRot1 = Flat3DArray.createFrom(aRotData1);
        
        executor = Executors.newFixedThreadPool(cpuCoreNumber);
        tmpTaskList = new ArrayList<>(3000);
        tmpPath = scratchDirectory 
                + FILESEPARATOR 
                + aParticlePair 
                + ".arc";
        int chunkIdx = 0;
        
        for (int i = 0; i < tmpDistanceNumber; i++) {
            tmpRotData2 = VectorUtil.moveX(aRotData2, aDistances[i]);
            
            for (int tmpRot2StartIdx = 0; tmpRot2StartIdx < totalElements; 
                    tmpRot2StartIdx += tmpChunkSize) {
                int tmpRot2EndIndex = Math.min(tmpRot2StartIdx + tmpChunkSize, 
                        totalElements);
                tmpTinkerXYZ = new TinkerXYZ(tmpTinkerXyz1, tmpTinkerXyz2, 
                        isTinkerOn);
                tmpAtomNumber = tmpTinkerXYZ.getN_atom();
                // This is for chunking tmpRotData2 to avoid memory issues within the thread
                double[][][] tmpRotPart2 = Arrays.copyOfRange(tmpRotData2, 
                        tmpRot2StartIdx, tmpRot2EndIndex);
                Flat3DArray sharedFlatRot2 = Flat3DArray
                        .createFrom(tmpRotPart2);
                tmpCmdList = new String[]{tinkerAnalyze, tmpPath + i + "_" 
                        + chunkIdx, "E"};
                tmpTaskList.add(new MIPETAnalyze(
                        tmpTinkerXYZ, 
                        isTinkerOn, 
                        i, 
                        chunkIdx, 
                        tmpAtomNumber,
                        minAtomDistance, 
                        sharedFlatRot1, 
                        sharedFlatRot2, 
                        scratchDirectory,
                        tmpCmdList, 
                        rec1,
                        rec2,
                        tmpIsFractionOne, 
                        temperature
                ));

                chunkIdx++;
            }
            
        }
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Calculate minimum of weighted intermolecular energy">
        /* If boltzmannFraction == 0.0, no averaging, min energy value of each configuration is taken
           If fractionForAverage = 1.0 all configurational E(nonbonded) values are used for "Boltzmann average" calculation
           0.0 < fractionForAverage < 1.0: All configurational E(nonbonded) values are sorted ascending and
           the lower "numberOfValues*fractionForAverage" E(nonbonded) values are used for "Boltzmann average" calculation
           Example: For 144x144x16 = 331776 E(nonbonded) values for a specific molecule distance r and
           a fractionForAverage of 0.25 the lowest Round(331776x0.25) = 82944 E(nonbonded) values are used for
           "Boltzmann average" calculation only */
        int tmpDistMinIndex;
        int tmpChunkMinIndex;
        int tmpTaskIndex;
        int tmpFractionToMax;
        double tmpDistMinEnergy;
        double allDistMinEnergy;
        double tmpEmin;
        double RezipRT;
        double tmpWgtEmin;
        double tmpSumWgt;
        double tmpSumWgtxE;
        double[] tmpWeights;
        double[] tmpEnergyDataFraction;
        double[] tmpEmins;
        double[] tmpWgtEmins;
        
        tmpTaskIndex = 0;
        tmpDistMinIndex = 0;
        tmpChunkMinIndex = 0;
        tmpSumWgt = 0;
        tmpSumWgtxE = 0;
        tmpDistMinEnergy = 1E10;
        allDistMinEnergy = 1E10;
        List<Future<WgtEnergyRecord>> tmpFutures = null;
        Future<WgtEnergyRecord> tmpFuture;
        List<double[]> allResultArrays = new ArrayList<>();
        tmpEmins = new double[tmpDistanceNumber];
        tmpWgtEmins = new double[tmpDistanceNumber];
        RezipRT = 1 / (temperature * GASCONST);
        tmpWgtEmin = 100.;
        
        try {            
            tmpFutures = executor.invokeAll(tmpTaskList);
        } catch (InterruptedException ex) {
            Thread.currentThread().interrupt();       
        }
        executor.shutdown();
        tmpFuture = null;
            
        for (int i = 0; i < tmpDistanceNumber; i++) {
            int totalSize = 0; // Size of energy values
                        
            for (int j = 0; j < tmpChunkNumber; j++) {
                if (tmpFutures != null) {
                    tmpFuture = tmpFutures.get(tmpTaskIndex);
                }
                try {
                    if (tmpFuture != null) {
                        // Thread results
                        double[] res = tmpFuture.get().energys();
                        allResultArrays.add(res);
                        totalSize += res.length;
                        tmpDistMinEnergy = tmpFuture.get().energys()[0];
                        if (tmpIsFractionOne) {
                            tmpSumWgt += tmpFuture.get().sumWgt(); 
                            tmpSumWgtxE += tmpFuture.get().sumWgtxE();
                        }
                    }
                    if (tmpDistMinEnergy < allDistMinEnergy) {
                        allDistMinEnergy = tmpDistMinEnergy;
                        tmpDistMinIndex = i;
                        tmpChunkMinIndex = j;
                    }
                    tmpTaskIndex++;
                } catch (InterruptedException | ExecutionException ex) {
                    LOGGER.log(Level.SEVERE,
                            "InterruptException during handling tmpFuture object.",
                            ex);
                } 
            }
            
            // Generate big target array for energy values from threads
            double[] globalEnergies = new double[totalSize];
            int currPosition = 0;
            
            for (double[] res : allResultArrays) {
                System.arraycopy(res, 0, 
                        globalEnergies, currPosition, res.length);
                currPosition += res.length;
            }
            
            allResultArrays.clear();
            Arrays.sort(globalEnergies);
            tmpEmins[i] = globalEnergies[0];
            
            // Store all energies at same distance to tmpEnergyDatas[i]
            if (tmpIsFractionOne) {
                tmpWgtEmins[i] = tmpSumWgtxE / tmpSumWgt;
                if (Double.isNaN(tmpWgtEmins[i]) || 
                        Double.isInfinite(tmpWgtEmins[i])) {
                    tmpWgtEmins[i] = 100.;
                } 
                tmpSumWgt = 0;
                tmpSumWgtxE = 0;
            } else {
                tmpFractionToMax = (int)(globalEnergies.length 
                        * boltzmannFraction);
                tmpEnergyDataFraction = new double[tmpFractionToMax];
                tmpWeights = new double[tmpFractionToMax];

                for (int j = 0; j < tmpFractionToMax; j++) {
                    tmpEnergyDataFraction[j] = globalEnergies[j];
                    tmpWeights[j] = Math.exp(-globalEnergies[j] * RezipRT);
                }

                tmpWgtEmins[i] = MIPETUTIL.productSum(tmpWeights, 
                        tmpEnergyDataFraction) / MIPETUTIL.sum(tmpWeights);
                if (Double.isNaN(tmpWgtEmins[i]) 
                        || Double.isInfinite(tmpWgtEmins[i])) {
                    tmpWgtEmins[i] = 100.;
                } 
            }
            if (tmpWgtEmin > tmpWgtEmins[i]) {
                tmpWgtEmin = tmpWgtEmins[i];
            }
        }
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Export .xyz file with the lowest intermolecular energy">
        tmpEmin = allDistMinEnergy;
        if (allDistMinEnergy < aMinEnergy) {
            int minTaskIndex = tmpChunkNumber * tmpDistMinIndex 
                    + tmpChunkMinIndex;
            TinkerXYZ minTinkerXYZ = null;
            
            if (tmpFutures != null && minTaskIndex >= 0 
                    && minTaskIndex < tmpFutures.size()) {
                try {
                    minTinkerXYZ = tmpFutures.get(minTaskIndex).get()
                            .minTinkerXYZ();
                } catch (InterruptedException ex) {
                    LOGGER.log(Level.SEVERE,
                        "Thread interrupted during future retrieval.", ex);
                } catch (ExecutionException ex) {
                    LOGGER.log(Level.SEVERE,
                        "Calculation failed in future task.", ex);
                } 
            } 
            if (minTinkerXYZ != null) {
                Path targetPath = Paths.get(scratchDirectory, 
                        aParticlePair + ".0");
                minTinkerXYZ.makeArcFile(targetPath.toString());
            }
        }
        tmpFutures = null;
        
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
        Instant startTime = Instant.now();
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
        if (!forceField_IE.equals(forceField_CN)) {
            readXyz(!isConformationalAnalysis, forceField_CN);
            if (forceField_CN.equals("OPLSAALIGPARGEN")) {
                readPrm(forceField_CN);
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
                tmpTargetDir = Paths.get(scratchDirectory, 
                        forceField_CN,
                        tmpParticlePair);
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
                tmpTarget = Paths.get(scratchDirectory,
                        forceField_CN,
                        tmpParticlePair,
                        tmpParticle1 + ".xyz");
                try {
                    Files.writeString(tmpTarget, xyzContent1[tmpXyzID1]);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during copying files to scratch.", ex);
                }
                if (!tmpIsSameParticle) {
                    tmpTarget = Paths.get(scratchDirectory,
                            forceField_CN,
                            tmpParticlePair,
                            tmpParticle2 + ".xyz");
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
        // MIPET also works parallel. Tests showed, the best performance is
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
        boolean tmpHasCNHeadLine;
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
                        + forceField_CN
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
                if (forceField_CN.equals("OPLSAALIGPARGEN")) {
                    tmpForcefieldName = "oplsaa";
                    
                } else {
                    tmpForcefieldName = forceField_CN.toLowerCase();
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
                if (forceField_CN.equals("OPLSAALIGPARGEN")) {
                    tmpKeyContent += prmContent1[tmpPrm1ID];
                    if (!tmpIsSameParticle) {
                        tmpKeyContent += prmContent2[tmpPrm2ID];
                    } 
                }
                MIPETUTIL.writeKeyFile(Paths.get(tmpKeyFileName), 
                        tmpKeyContent);

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
                                // This loop is necessary for Linux version
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
                boxMinimizationTime = System.nanoTime();
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
                    // BufferedWriter for log contents
                    try (BufferedWriter BWParticleLog = new BufferedWriter(
                             new FileWriter(tmpParticleLogFileName))) {
                        if (Files.exists(tmpSourceFile2) && !tmpHasCNHeadLine) {
                            Files.move(tmpSource, tmpTarget,
                                    StandardCopyOption.ATOMIC_MOVE);
                            BWParticleLog.append(LINESEPARATOR)
                                    .append("*** Coordination number ***")
                                    .append(LINESEPARATOR);
                        }
                        BWParticleLog.append("Time for solvent box minimization via Tinker \"Minimize\" [s]: ")
                                .append(Double.toString(
                                        (System.nanoTime()
                                                - boxMinimizationTime) / 1000d))
                                .append(LINESEPARATOR);
                    }
                    if (!tmpIsSameParticle) {
                        BWParticleDat.append("VdWVolume(")
                                .append(tmpParticle1)
                                .append(") [" + ANGSTROM + SUPERSCRIPTTHREE 
                                        +"] = ")
                                .append(Double.toString(tmpVdWSolutVolumes[tmpJobIndex]))
                                .append(LINESEPARATOR);
                    }
                    BWParticleDat.append("VdWVolume(")
                            .append(tmpParticle2)
                            .append(") [" + ANGSTROM + SUPERSCRIPTTHREE 
                                    +"] = ")
                            .append(String.format("%.4f", 
                                    tmpVdWSolventVolumes[tmpJobIndex]))
                            .append(LINESEPARATOR);
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
                        + forceField_CN
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

        tmpDynamicWarmUpTime = System.nanoTime();
        tmpExecutor = Executors.newFixedThreadPool(cpuCoreNumber);
        try {
            tmpExecutor.invokeAll(tmpTaskList);
        } catch (InterruptedException ex) {
            Thread.currentThread().interrupt();       
        }
        tmpExecutor.shutdown();
        String[][] tmpLabelValues = new String[1][2];
        tmpLabelValues[0][0] = "Time for solvent box warm up via Tinker \"Dynamic\" [s]: ";
        tmpLabelValues[0][1] = Double.toString((System.nanoTime() 
                - tmpDynamicWarmUpTime) / 1000d);
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
                        + forceField_CN
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
        
        long tmpArcEvaluationTime = System.nanoTime();
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
            tmpLabelValues[0][1] = Double.toString(
                    (System.nanoTime() - tmpArcEvaluationTime) / 1000d);
            tmpLabelValues[1][0] = "Total time for determining the coordination number [s]: ";
            Instant endTime = Instant.now();
            Duration duration = Duration.between(startTime, endTime);
            float totalTimeSeconds = duration.toMillis() /1000f;
            tmpLabelValues[1][1] = String.valueOf(totalTimeSeconds);
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
                            + forceField_CN
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
                + forceField_IE;
        
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
            tmpParticleDescriptions.add(tmpParticleName); //full name
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
            if (!tmpParticleName.equals("[Na+]")) {
                tmpParticleDescriptions.add(String.format("%.4f", 
                    MIPETUTIL.getVdwVolume(smiles.get(tmpParticleName)) 
                            * WATERVOLUMERATIO)); // volume [A³]
            }
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
        double tmpEnergy;
        double tmpCN;
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
        HashMap<String, Double> tmpEnergiesMap =
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
                tmpEnergiesMap.put(tmpParticlePair, tmpEnergy);
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
                    tmpE12 = tmpEnergiesMap.get(tmpParticleName1 + "_"
                        + tmpParticleName2);
                    tmpE11 = tmpEnergiesMap.get(tmpParticleName1 + "_"
                            + tmpParticleName1);
                    tmpE22 = tmpEnergiesMap.get(tmpParticleName2 + "_"
                            + tmpParticleName2);
                    if (!forceField_CN.isEmpty()) {
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
                        tmpAij = temperature / 12 + 1.7483
                                * tmpChiNumerator;
                        tmpAijMap.put(tmpParticleName1 + "_" + tmpParticleName2, 
                                tmpAij);
                    }

                    // Calculation for CN = 1
                    tmpChiNumerator = tmpE12 + tmpE12 - tmpE11 - tmpE22;
                    tmpAij = temperature / 12 + 1.7483 * tmpChiNumerator;
                    tmpAijMap1.put(tmpParticleName1 + "_" + tmpParticleName2, 
                            tmpAij);
                }
            }
            
            /* Write file */
            tmpFileName = tmpResultsDirectory 
                + FILESEPARATOR 
                + aTitleAbbreviation + "_"
                + "EijFraction_"
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
            if (!forceField_CN.isEmpty() || tmpOutputIteration <= 3) {
                try (BufferedWriter tmpBW = Files.newBufferedWriter(
                        Paths.get(tmpFileName))) {
                    tmpBW.append("# Particle set for MFSim created by MIPET\n")
                            .append("# Force Field for energy calculation: ")
                            .append(forceField_IE)
                            .append(LINESEPARATOR)
                            .append("# Force Field for coordination number calculation: ");
                    if (forceField_CN.isEmpty()) {
                        tmpBW.append("Not calculated");
                    } else {
                        tmpBW.append(forceField_CN);
                    }
                    tmpBW.append(forceField_CN)
                            .append(LINESEPARATOR)
                            .append("# Water model: ")
                            .append(waterModel)
                            .append(LINESEPARATOR)
                            .append("# CPU cores: ")
                            .append(Integer.toString(cpuCoreNumber))
                            .append(LINESEPARATOR)
                            .append("# Temperature: ")
                            .append(Double.toString(temperature))
                            .append(LINESEPARATOR)
                            .append("# Sphere nodes calculated with Fibonacci algorithm: ")
                            .append(Boolean.toString(isFibonacciSphereAlgorithm))
                            .append(LINESEPARATOR)
                            .append("# Sphere node number 1: ")
                            .append(Integer.toString(sphereNodeNumber1))
                            .append(LINESEPARATOR)
                            .append("# Sphere node number 2: ")
                            .append(Integer.toString(sphereNodeNumber2))
                            .append(LINESEPARATOR)
                            .append("# Sphere node number 3: ")
                            .append(Integer.toString(sphereNodeNumber3))
                            .append(LINESEPARATOR)
                            .append("# Sphere node number 4: ")
                            .append(Integer.toString(sphereNodeNumber4))
                            .append(LINESEPARATOR)
                            .append("# Rotation number 1: ")
                            .append(Integer.toString(rotationNumber1))
                            .append(LINESEPARATOR)
                            .append("# Rotation number 2: ")
                            .append(Integer.toString(rotationNumber2))
                            .append(LINESEPARATOR)
                            .append("# Rotation number 3: ")
                            .append(Integer.toString(rotationNumber3))
                            .append(LINESEPARATOR)
                            .append("# Rotation number 4: ")
                            .append(Integer.toString(rotationNumber4))
                            .append(LINESEPARATOR)
                            .append("# Conformational analysis: ")
                            .append(Boolean.toString(isConformationalAnalysis))
                            .append(LINESEPARATOR);
                    switch (tmpOutputIteration) {
                        case 0, 4 -> {
                            tmpBW.append("# Boltzmann averaging: Yes")
                                    .append(LINESEPARATOR);}
                        case 1, 2, 3, 5, 6, 7 -> {
                            tmpBW.append("# Boltzmann averaging: No")
                                    .append(LINESEPARATOR);}
                    }
                    switch (tmpOutputIteration) {
                        case 1, 6 -> {
                            tmpBW.append("# Optimize sampled E(min) configuration: Yes")
                                    .append(LINESEPARATOR);}
                        case 0, 2, 3, 4, 5, 7 -> {
                            tmpBW.append("# Optimize sampled E(min) configuration: No")
                                    .append(LINESEPARATOR);}
                    }
                    switch (tmpOutputIteration) {
                        case 2, 6 -> {
                            tmpBW.append("# Tinker's 'optrigid' used: Yes")
                                    .append(LINESEPARATOR);}
                    }
                    switch (tmpOutputIteration) {
                        case 0, 1, 2, 3 -> {
                            tmpBW.append("# CN = 1 for all particle pairs")
                                    .append(LINESEPARATOR);}
                    }
                    tmpBW.append(LINESEPARATOR)
                            .append("[Title]\n")
                            .append(aTitle)
                            .append("\n[/Title]\n\n")
                            .append("[Version]\n" 
                                    + VERSION_NUMBER 
                                    + "\n[/Version]\n\n");

                    // Particle description
                    tmpBW.append("[Particle Description]")
                            .append(LINESEPARATOR);
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
                        tmpBW.append(tmpHeader)
                                .append("\t");
                    }

                    for(List<String> tmpParticleDescList : 
                            tmpParticleDescriptionsList){
                        tmpBW.append(LINESEPARATOR);

                        for(String tmpDescription : tmpParticleDescList){
                            tmpBW.append(tmpDescription)
                                    .append(" ");
                        }

                    }

                    tmpBW.append(LINESEPARATOR)
                            .append("[/Particle Description]")
                            .append(LINESEPARATOR)
                            .append(LINESEPARATOR);

                     // Particle interactions
                    tmpBW.append("[Particle interactions]")
                            .append(LINESEPARATOR)
                            .append("# Repulsion parameters a(ij) for particle pairs for different temperatures (in K)")
                            .append(LINESEPARATOR)
                            .append("Pair")
                            .append(" ")
                            .append(Integer.toString((int) temperature));
                    switch (tmpOutputIteration) {
                        case 0, 1, 2, 3 -> {
                            tmpKeySet = tmpAijMap1.keySet();}
                        case 4, 5, 6, 7 -> {
                            tmpKeySet = tmpAijMap.keySet();}
                    }

                    for(String tmpKey : tmpKeySet) {
                        tmpBW.append(LINESEPARATOR)
                                .append(tmpKey);
                        switch (tmpOutputIteration) {
                            case 0, 1, 2, 3 -> {
                                tmpAij = tmpAijMap1.get(tmpKey);}
                            case 4, 5, 6, 7 -> {
                                tmpAij = tmpAijMap.get(tmpKey);}
                        }
                        tmpBW.append(" ").append(String.format("%.2f", tmpAij));
                    }

                    tmpBW.append(LINESEPARATOR)
                            .append("[/Particle interactions]")
                            .append(LINESEPARATOR)
                            .append(LINESEPARATOR);

                    // Coordination numbers
                    switch (tmpOutputIteration) {
                        case 4, 5, 6, 7 -> {
                            tmpBW.append("[Coordination numbers]");
                            tmpCNListLength = cnList.size();

                            for (int j = 0; j < tmpCNListLength; j++) {
                                tmpBW.append(LINESEPARATOR);
                                tmpBW.append(cnList.get(j).particleName1())
                                        .append("_")
                                        .append(cnList.get(j).particleName2());
                                tmpBW.append(String.format(" %.2f", 
                                        cnList.get(j).cnValue()));
                            }

                            tmpBW.append(LINESEPARATOR)
                                    .append("[/Coordination numbers]")
                                    .append(LINESEPARATOR);}
                    }

                    // SMILES
                    tmpBW.append(LINESEPARATOR)
                            .append("[SMILES]")
                            .append(LINESEPARATOR)
                            .append("# Particle and corresponding SMILES of fragment molecule");
                    tmpKeySet = smiles.keySet();

                    for(String tmpKey : tmpKeySet){
                        if (tmpParticleNames.contains(tmpKey)) {
                            tmpBW.append(LINESEPARATOR)
                                    .append(tmpKey)
                                    .append(" ")
                                    .append(smiles.get(tmpKey));
                        }
                    }

                    tmpBW.append(LINESEPARATOR)
                            .append("[/SMILES]");
                } catch (IOException anException) {
                    LOGGER.log(Level.SEVERE, anException.toString());
                }
            }
            if (isTinkerOn) {
                tmpOutputIteration++;
            } else {
                tmpOutputIteration = 8;
            }
        }
            
    }
    
    // </editor-fold>
}