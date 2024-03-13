/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
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
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
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
import org.apache.commons.lang3.ArrayUtils;

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
    private static final char ANGSTROM = '\u212B';
    
    /**
     * Superscript 3 character
     */
    private static final char SUPERSCRIPTTHREE = '\u00B3';
    
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger
            .getLogger(MIPET.class.getName());
    
    /**
     * Line separator
     */
    final private static String LINESEPARATOR = System.lineSeparator();
    
    /**
     * File separator
     */
    final private static String FILESEPARATOR = System
            .getProperty("file.separator");
    
    /**
     * Gas constant R [kcal/(mol*K)]
     */
    final private static double GASCONST = 1.98720425864E-3;
    
    /**
     * Instance object of MIPETUtility
     */
    final private static MIPETUtility MIPETUTIL = new MIPETUtility();
    
    /**
     * Water volume ratio
     */
    final private static double WATERVOLUMERATIO = 30.0036 / 
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
    private static final String VERSION_NUMBER = "1.0.0";
    
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
     * Job file name
     */
    private static String jobFileName;
    
    /**
     * Extra force field name for intermolecular energy calculation
     */
    final private static HashSet<String> extraForcefield_IE = new HashSet<>();
    
    /**
     * Extra force field name for coordination number calculation
     */
    final private static HashSet<String> extraForcefield_CN = new HashSet<>();
    
    /**
     * New particles for calculation
     */
    private static LinkedList<String> newParticles;
    
    /**
     * Old particles already calculated
     */
    private static LinkedList<String> oldParticles;
    
    /**
     * CPUcore number
     */
    private static int cpuCoreNumber;
    
    /**
     * Use fibonaccisphere algorithm
     */
    private static boolean isFibonacciSphereAlgorithm ;
    
    /**
     * Temperature
     */
    private static int temperature;
    
    /**
     * Sphere node number
     */
    private static int sphereNodeNumber;

    /**
     * Force field names for calculation of intermolecular interaction energy
     */
    private static String forcefield_IE;
    
    /**
     * Force field names for calculation of coordination number
     */
    private static String forcefield_CN;
    
    /**
     * ConfigurationNumber of particle 1
     */
    private static int confNumber1;
    
    /**
     * ConfigurationNumber of particle 2
     */
    private static int confNumber2;
    
    /**
     * Rotation configuration number
     */
    private static int rotationNumber;
    
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
        // Please set Your IDE output to UTF-8, otherwise it will cause issue due to symbol like Angstrom
        //<editor-fold defaultstate="collapsed" desc="Variables declaration">
        ArrayList<String> tmpParticlePairs = new ArrayList<>(500);
        ArrayList<JobTaskRecord> tmpJobTaskRecordList = new ArrayList<>(500);
        List<double[]> tmpSphereNodeCoord;
        List<double[][]> tmpRotMatrices1;
        List<double[][]> tmpRotMatrices2;
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
        String tmpParticleLogFileName;
        String tmpParticleDatFileName;
        Path tmpParticleLogFile;
        Path tmpParticleDatFile;
        String tmpParticleName1;
        String tmpParticleName2;
        String[] tmpElementList1;
        String[] tmpElementList2;
        TinkerXYZ tmpTinkerXYZ1;
        TinkerXYZ tmpTinkerXYZ2;
        Process tmpProcess;
        ProcessBuilder tmpPB;
        String[] tmpCmdList;
        ArrayIndexComparator tmpComparator;
        double[][] tmpXyzData1; // coordination of 1. fragment atoms
        double[][] tmpXyzData2; // coordination of 2. fragment atoms
        double[] tmpCentre1;
        double[] tmpCentre2;
        double[][][] tmpXyzRotData1; // i: configuration, j: atom, k: xyz
        double[][][] tmpXyzRotData2; // i: configuration, j: atom, k: xyz

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
        keyFileStringOrigin =
                "# Force Field Selection"
                + LINESEPARATOR
                + "PARAMETERS\t\""
                + parameterDirectory
                + "/"
                + forcefield_IE.toLowerCase() + ".prm\""
                + LINESEPARATOR;
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
        confNumber1 = sphereNodeNumber;
	confNumber2 = sphereNodeNumber * rotationNumber;
        // Development version
        String tmpFileNameSphereNode = 
                "de/whs/ibci/mipet/sphereNodes/SphereNodes"
                + sphereNodeNumber + ".txt";
        
        // Distribution version
//        String tmpFileNameSphereNode = 
//                "/de/whs/ibci/mipet/sphereNodes/SphereNodes"
//                + sphereNodeNumber + ".txt";
        
        // Determine rotation matrices used to rotate 
        //   the particle/atom coordinates
        if (isFibonacciSphereAlgorithm) {
            tmpSphereNodeCoord = FibonacciSphere
                    .getSphereNodes(sphereNodeNumber);
        } else {
            tmpSphereNodeCoord = RotationUtil
                .readSphereNodes (tmpFileNameSphereNode);
        }
        tmpRotMatrices1 = RotationUtil
                .getRotationMatrices1 (tmpSphereNodeCoord, 
                        new double[]{1., 0., 0.});
        tmpRotMatrices2 = RotationUtil.
                getRotationMatrices2 (tmpSphereNodeCoord, 
                        new double[] {-1.0, 0.0, 0.0}, rotationNumber);
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
                    Files.createDirectories(Paths.get(tmpCNResultDirName));
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
        
        //<editor-fold defaultstate="collapsed" desc="Optimize and scan particle">
        System.out.println("Scan particle...");
        scanParticle(tmpJobTaskRecordList);
        
        //</editor-fold>
        
        Path tmpKeyFile;
        Path tmpOptDistFile;
        String tmpOutputName;
        String tmpLine;
        
        tmpCurrentIndex = 0;
        tmpIsExitCondition = false;
        tmpKeyFileString = keyFileStringOrigin 
                + "DIELECTRIC\t" 
                + dielectricConstant;	
        System.out.println("Calculating intermolecular energy...");
        
        while (tmpIsExitCondition == false) {
            // Exit condition is true when all particle pair combinations
            //   were calculated.
            if (tmpJobTaskRecordList.get(tmpCurrentIndex).hasEnergieJob()) {
                //<editor-fold defaultstate="collapsed" desc="Job task record">
                tmpForcefield = tmpJobTaskRecordList.get(tmpCurrentIndex)
                        .forcefield_IE_Name();
                tmpParticleName1 = tmpJobTaskRecordList.get(tmpCurrentIndex)
                        .particleName1();
                tmpParticleName2 = tmpJobTaskRecordList.get(tmpCurrentIndex)
                        .particleName2();
                tmpParticlePair = tmpParticleName1 + "_" + tmpParticleName2;
                tmpIsSameParticle = tmpParticleName1.equals(tmpParticleName2);
                tmpIEResultDirName = tmpJobTaskRecordList.get(tmpCurrentIndex)
                        .result_IE_PathName();
            
                //</editor-fold>
            
                //<editor-fold defaultstate="collapsed" desc="Create log file">
                tmpParticleLogFileName = tmpIEResultDirName
                        + FILESEPARATOR
                        + tmpParticlePair 
                        + "_log.txt";
                tmpParticleDatFileName = tmpIEResultDirName
                        + FILESEPARATOR
                        + tmpParticlePair
                        + ".dat";
                tmpParticleLogFile = Paths.get(tmpParticleLogFileName);
                tmpParticleDatFile = Paths.get(tmpParticleDatFileName);
                try {
                    BWParticleLog = new BufferedWriter(
                            new FileWriter(tmpParticleLogFileName));
                    if (!Files.exists(tmpParticleLogFile)) {
                        Files.createFile(tmpParticleLogFile);
                    } else {
                        BWParticleLog.append(LINESEPARATOR);
                    }
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
                tmpTinkerXYZ1 = new TinkerXYZ(optXYZDirectory
                        + FILESEPARATOR
                        + forcefield_IE
                        + FILESEPARATOR
                        + tmpParticleName1
                        + FILESEPARATOR
                        + tmpParticleName1 + ".xyz");
                tmpXyzData1 = tmpTinkerXYZ1.getCoordinateList1()[0];
                tmpCentre1 = tmpTinkerXYZ1.findCentreCoordinate();
                tmpXyzData1 = tmpTinkerXYZ1
                        .moveCoordinates(tmpXyzData1, tmpCentre1);
                if (tmpIsSameParticle) {
                    tmpTinkerXYZ2 = tmpTinkerXYZ1.clone();
                    tmpXyzData2 = tmpXyzData1.clone();
                } else {
                    tmpTinkerXYZ2 = new TinkerXYZ(optXYZDirectory
                            + FILESEPARATOR
                            + forcefield_IE
                            + FILESEPARATOR
                            + tmpParticleName2
                            + FILESEPARATOR
                            + tmpParticleName2 + ".xyz");
                    tmpXyzData2 = tmpTinkerXYZ2.getCoordinateList1()[0];
                    tmpCentre2 = tmpTinkerXYZ2.findCentreCoordinate();
                    tmpXyzData2 = tmpTinkerXYZ2
                            .moveCoordinates(tmpXyzData2, tmpCentre2);
                }

                //</editor-fold>

                //<editor-fold defaultstate="collapsed" desc="Calculate rotated coordinates">
                // Calculates the rotated atom coordinates 
                //   using the rotation matrices
                tmpXyzRotData1 = 
                        new double[confNumber1][tmpXyzData1.length][3];
                tmpXyzRotData2 = 
                        new double[confNumber2][tmpXyzData2.length][3];

                for (int i = 0; i < confNumber1; i++) {

                    for (int j = 0; j < tmpXyzData1.length; j++) {
                            tmpXyzRotData1[i][j] = 
                                MatrixUtil.multiply(tmpRotMatrices1.get(i),
                                        tmpXyzData1[j]);
                    }

                }

                for (int i = 0; i < confNumber2; i++) {

                    for (int j = 0; j < tmpXyzData2.length; j++) {
                        tmpXyzRotData2[i][j] = 
                                MatrixUtil.multiply(tmpRotMatrices2.get(i),
                                        tmpXyzData2[j]);
                    }

                }
                //</editor-fold>

                //<editor-fold defaultstate="collapsed" desc="Calculate Intermolecular Energy">

                // Scans for the lowest intermolecular energy between particles
                // Prescan
                long tmpEnergyCalcTime = System.currentTimeMillis();
                int tmpDistSize = (int)Math.ceil(
                        (upperBoundary - lowerBoundary) / prescanStepSize) + 1;
                LinkedList<Double> tmpAllDistances = new LinkedList<>();
                LinkedList<Double> tmpDistanceList = new LinkedList<>();
                double[][] tmpEnergyDatas;
                double[][] tmpDistMinEnergyDatas;
                double tmpMinDistance;
                double tmpMinEnergy;
                double tmpDistanceCandidate;
                String tmpKeyContent;

                // Write .key file
                tmpKeyContent = keyFileStringOrigin 
                        + "DIELECTRIC\t" 
                        + dielectricConstant 
                        + LINESEPARATOR
                        + "NONBONDTERM ONLY"
                        + LINESEPARATOR
                        + "OPENMP-THREADS 1";
                tmpKeyPathName = scratchDirectory
                        + FILESEPARATOR 
                        + tmpParticlePair 
                        + ".key";
                MIPETUTIL.saveKeyFile(parameterDirectory,
                        tmpForcefield, 
                        tmpParticleName1,
                        tmpParticleName2,
                        tmpKeyPathName,
                        tmpKeyContent);

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

                EnergyRecord[] tmpEnergyRecords = new EnergyRecord[3];
                tmpEnergyRecords[0] = getInterMolecularEnergy(
                        tmpParticlePair,
                        tmpDistances, 
                        tmpTinkerXYZ1, 
                        tmpTinkerXYZ2, 
                        tmpXyzRotData1, 
                        tmpXyzRotData2,
                        1E10);
                tmpMinEnergy = tmpEnergyRecords[0].minEnergy();
                tmpMinDistance = tmpEnergyRecords[0].minDistance();

                // Precise scan
                tmpDistanceList = new LinkedList<>();
                tmpDistSize = 9;

                for (int i = 0; i < tmpDistSize; i++) {
                    tmpDistanceCandidate = tmpMinDistance - 0.4 + i * 0.1;
                    if(!MIPETUTIL
                            .contains(tmpAllDistances, tmpDistanceCandidate)) {
                        tmpDistanceList.add(tmpDistanceCandidate);
                    }
                }

                tmpAllDistances.addAll(tmpDistanceList);
                tmpDistances = new double[tmpDistanceList.size()];

                for (int i = 0; i < tmpDistances.length; i++) {
                    tmpDistances[i] = tmpDistanceList.get(i);
                }

                tmpEnergyRecords[1] = getInterMolecularEnergy(
                        tmpParticlePair,
                        tmpDistances, 
                        tmpTinkerXYZ1, 
                        tmpTinkerXYZ2, 
                        tmpXyzRotData1, 
                        tmpXyzRotData2,
                        tmpMinEnergy);
                if (tmpEnergyRecords[1].minEnergy() < tmpMinEnergy) {
                    tmpMinEnergy = tmpEnergyRecords[1].minEnergy();
                    tmpMinDistance = tmpEnergyRecords[1].minDistance();
                }

                // More precise scan
                tmpDistanceList = new LinkedList<>();
                tmpDistSize = 19;

                for (int i = 0; i < tmpDistSize; i++) {
                    tmpDistanceCandidate = tmpMinDistance - 0.09 + i * 0.01;
                    if(!MIPETUTIL
                            .contains(tmpAllDistances, tmpDistanceCandidate)) {
                        tmpDistanceList.add(tmpDistanceCandidate);
                    }
                }

                tmpAllDistances.addAll(tmpDistanceList);
                tmpDistances = new double[tmpDistanceList.size()];

                for (int i = 0; i < tmpDistances.length; i++) {
                    tmpDistances[i] = tmpDistanceList.get(i);
                }

                tmpEnergyRecords[2] = getInterMolecularEnergy(
                        tmpParticlePair,
                        tmpDistances, 
                        tmpTinkerXYZ1, 
                        tmpTinkerXYZ2, 
                        tmpXyzRotData1, 
                        tmpXyzRotData2,
                        tmpMinEnergy);
                if (tmpEnergyRecords[2].minEnergy() < tmpMinEnergy) {
                    //tmpMinEnergy = tmpEnergyRecords[2].minEnergy();
                    tmpMinDistance = tmpEnergyRecords[2].minDistance();
                }
                tmpEnergyDatas = ArrayUtils
                        .addAll(tmpEnergyRecords[0].energyDatas(),
                                tmpEnergyRecords[1].energyDatas());
                tmpEnergyDatas = ArrayUtils.addAll(tmpEnergyDatas,
                        tmpEnergyRecords[2].energyDatas());
                
                // Copy distance and minenergy datas
                tmpDistSize = tmpAllDistances.size();
                tmpDistMinEnergyDatas = new double[tmpDistSize][3];
                int tmpIteration;
                int tmpIndex = 0;
               
                for (int i = 0; i < 3; i++) {
                    tmpIteration = tmpEnergyRecords[i].distances().length;

                    for (int j = 0; j < tmpIteration; j++) {
                        tmpDistMinEnergyDatas[tmpIndex][0] = 
                                tmpEnergyRecords[i].distances()[j];
                        tmpDistMinEnergyDatas[tmpIndex][1] = 
                                tmpEnergyRecords[i].energyDatas()[j][0];
                        tmpIndex++;
                    }
                    
                }
                
                Arrays.sort(tmpDistMinEnergyDatas, 
                        (a, b) -> Double.compare(a[0], b[0]));
                
                //</editor-fold>

                //<editor-fold defaultstate="collapsed" desc="Sort datas">
                Double[] tmpDistanceObj;
                double[][] tmpEnergySorted;
                Integer[] tmpDistanceIndices;

                tmpDistanceObj = tmpAllDistances.toArray(Double[]::new);
                tmpComparator = new ArrayIndexComparator(tmpDistanceObj);
                tmpDistanceIndices = tmpComparator.createIndexArray();
                Arrays.sort(tmpDistanceIndices, tmpComparator);
                tmpEnergySorted = new double[tmpDistSize][];

                for (int i = 0; i < tmpDistSize; i++) {
                  tmpIndex = tmpDistanceIndices[i];
                  tmpEnergySorted[i] = tmpEnergyDatas[tmpIndex].clone();
                }

                //</editor-fold>

                //<editor-fold defaultstate="collapsed" desc="Minimize and analyze">

                // Delete old .key file and make new one
                tmpKeyFile = Paths.get(scratchDirectory 
                        + FILESEPARATOR
                        + tmpParticlePair+ ".key");
                try {
                    Files.deleteIfExists(tmpKeyFile);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE,
                            "IOException during deleting .key file in scratch.",
                            ex);
                }
                tmpKeyFileName = tmpParticlePair + ".key";
                tmpKeyPathName = scratchDirectory 
                        + FILESEPARATOR 
                        + tmpKeyFileName;
                MIPETUTIL.saveKeyFile(parameterDirectory, 
                        tmpForcefield,
                        tmpParticleName1,
                        tmpParticleName2,
                        tmpKeyPathName, 
                        tmpKeyFileString);

                // Optimize lowest energy conformation
                tmpCmdList = new String[] {tinkerOptimize, 
                        scratchDirectory 
                        + FILESEPARATOR 
                        + tmpParticlePair
                        + ".0",        
                        String.valueOf(optimizeRmsGradient)};
                tmpPB = new ProcessBuilder();
                tmpPB.command(tmpCmdList);
                tmpProcess = null;
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
                            new InputStreamReader(tmpProcess.getInputStream()))) {
                        while (tmpBR.readLine() != null ) {
                        }
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

                // Use tinker's analyze.exe to determine intermolecular energy
                tmpOutputName = scratchDirectory
                        + FILESEPARATOR
                        + "output0_opt.txt";
                File tmpOptFile = new File(tmpOutputName);
                tmpCmdList = new String[] {tinkerAnalyze, 
                    scratchDirectory
                    + FILESEPARATOR 
                    + tmpParticlePair
                    + ".xyz",
                    "E"};
                tmpPB = new ProcessBuilder();
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

                // Read the intermolecular energy
                String tmpSearch = "Intermolecular Energy :";
                double tmpOptMinEnergy = 0.0;
                try (BufferedReader tmpBR = new BufferedReader(
                        new FileReader(tmpOptFile), 65536)) {
                    while ((tmpLine = tmpBR.readLine()) != null ) {
                        if (tmpLine.contains(tmpSearch)) {
                            tmpOptMinEnergy = Double.parseDouble(tmpLine
                                    .substring(25, 50));
                            break;
                        }
                    }
                    tmpBR.close();
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during reading output0_opt.txt.",
                            ex);
                }

                //</editor-fold>

                //<editor-fold defaultstate="collapsed" desc="Calculate intermolecular energy of all configurations">
                int tmpFractionToMax;
                double[] tmpWeights;
                double[] tmpEnergyDataFraction;

                // If boltzmannFraction == 0.0, no averaging, min energy value of each configuration is taken
                // If fractionForAverage = 1.0 all configurational E(nonbonded) values are used for "Boltzmann average" calculation
                // 0.0 < fractionForAverage < 1.0: All configurational E(nonbonded) values are sorted ascending and
                // the lower "numberOfValues*fractionForAverage" E(nonbonded) values are used for "Boltzmann average" calculation
                // Example: For 144x144x16 = 331776 E(nonbonded) values for a specific molecule distance r and
                // a fractionForAverage of 0.25 the lowest Round(331776x0.25) = 82944 E(nonbonded) values are used for
                // "Boltzmann average" calculation only

                if (boltzmannFraction == 0.0) {
                    tmpMinEnergy = tmpOptMinEnergy;
                } else {

                    for (int i = 0; i < tmpDistSize; i++) {
                        tmpFractionToMax = (int)(tmpEnergySorted[i].length
                                * boltzmannFraction);
                        tmpEnergyDataFraction = new double[tmpFractionToMax];
                        tmpWeights = new double[tmpFractionToMax];

                        for (int j = 0; j < tmpFractionToMax; j++) {
                            tmpEnergyDataFraction[j] = tmpEnergySorted[i][j];
                            tmpWeights[j] = Math.exp(-(tmpEnergySorted[i][j] 
                                    - tmpOptMinEnergy) / (temperature * GASCONST));
                        }

                        tmpDistMinEnergyDatas[i][2] = MIPETUTIL
                                .productSum(tmpWeights, tmpEnergyDataFraction) 
                                / MIPETUTIL.sum(tmpWeights);
                    }

                    // Find minimum
                    tmpMinEnergy= 100.0;

                    for (int i = 0; i < tmpDistSize; i++) {
                        if(tmpDistMinEnergyDatas[i][2] < tmpMinEnergy) {
                            tmpMinEnergy = tmpDistMinEnergyDatas[i][2];
                            tmpMinDistance = tmpDistMinEnergyDatas[i][0];
                        }
                    }

                }
                
                energyList.add(new ResultEnergyRecord(
                        tmpParticleName1, 
                        tmpParticleName2, 
                        tmpMinEnergy));

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
                        tmpBW.append(String.format("%12.2f", 
                                tmpDistMinEnergyDatas[i][0]));
                        tmpBW.append("  ");
                        tmpBW.append(String.format("%12.3f", 
                                tmpDistMinEnergyDatas[i][1]));
                        tmpBW.append("  ");
                        tmpBW.append(String.format("%20.3f", 
                                tmpDistMinEnergyDatas[i][2]));
                        tmpBW.append(LINESEPARATOR);
                    }
                    
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during writing output0.out in scratch."
                            , ex);
                }
                
                //</editor-fold>

                //<editor-fold defaultstate="collapsed" desc="Write log file">
                try {
                    BFGblLog.append(tmpParticlePair);
                    BFGblLog.append(LINESEPARATOR);
                    BWParticleLog.append("*** Eij ***");
                    BWParticleLog.append(LINESEPARATOR);
                    BWParticleLog.append("Raw time for determining the intermolecular minimal energy [s]: ");
                    BWParticleLog.append(String.valueOf((System.currentTimeMillis() 
                                    - tmpEnergyCalcTime) / 1000));
                    BWParticleLog.append(LINESEPARATOR);
                    BWParticleLog.append("New particle pair: ");
                    BWParticleLog.append("Minimum IntermolecularEnergy [kcal/mole] = "); 
                    BWParticleLog.append(String.format("%.4f", tmpMinEnergy));
                    BWParticleLog.append(LINESEPARATOR);
                    BWParticleLog.append("Fraction of energy values used for the Boltzmann distribution: ");
                    BWParticleLog.append(Double.toString(boltzmannFraction));
                    BWParticleLog.append("    ");
                    BWParticleLog.append(LINESEPARATOR);
                    BWParticleLog.append(Integer.toString(temperature));
                    BWParticleLog.append("    ");
                    BWParticleLog.append(String.format("%.4f", tmpMinEnergy));
                    BWParticleLog.append("    ");
                    BWParticleLog.append(LINESEPARATOR);
                    BWParticleLog.flush();
                    BWParticleDat.append("equilibriumDistances [" + ANGSTROM + "] = "); 
                    BWParticleDat.append(String.format("%.2f", tmpMinDistance));
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("Fraction of energy values used for the Boltzmann distribution: ");
                    BWParticleDat.append(Double.toString(boltzmannFraction));
                    BWParticleDat.append("    ");
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("    ");
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("Temperature [K]: ");
                    BWParticleDat.append(Integer.toString(temperature));
                    BWParticleDat.append("    ");
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.flush();
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during writing log file.", ex);
                }

                //</editor-fold>

                //<editor-fold defaultstate="collapsed" desc="Copy results">
                Path tmpOriginal;
                Path tmpTarget;
                String tmpFileName;
                Boolean tmpHasH2O;
                
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
                        + tmpEnergyRecords[2].minEnergy() + " kcal/mol";
                String tmpTargetDir = tmpIEResultDirName 
                        + FILESEPARATOR 
                        + "output0.out";        
                try (BufferedWriter tmpBW = new BufferedWriter(
                        new FileWriter(tmpTargetDir))) {
                    tmpBW.append(tmpOutput0);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during writing output0.out in scratch."
                            , ex);
                }
                tmpOriginal = Paths.get(scratchDirectory, tmpParticlePair + ".xyz");
                tmpTarget = Paths.get(tmpIEResultDirName,"output_optimized.0");
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
                        + "output0_optimized.out";        
                try (BufferedWriter tmpBW = new BufferedWriter(
                        new FileWriter(tmpTargetDir))) {
                    tmpBW.append(tmpOutput0);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during writing output0_optimized.out "
                            + "in scratch.", ex);
                }
                tmpFileName = tmpIEResultDirName + FILESEPARATOR + "output.0";
                String tmpOptFileName;
                int tmpAtomNumber1 = tmpTinkerXYZ1.getAtomSize1();
                int tmpAtomNumber2 = tmpTinkerXYZ2.getAtomSize1();
                tmpElementList1 = tmpTinkerXYZ1.getElementList1();
                tmpElementList2 = tmpTinkerXYZ2.getElementList1();
                TinkerXYZ tmpTinkerXyz = new TinkerXYZ(tmpFileName, 1, 
                        tmpAtomNumber1, tmpAtomNumber2);
                tmpFileName = tmpIEResultDirName + FILESEPARATOR + "output.xyz";
                tmpTinkerXyz.writeToXyzFile(tmpFileName);
                tmpOptFileName = tmpIEResultDirName 
                        + FILESEPARATOR 
                        + "output_optimized.0";
                tmpTinkerXyz = new TinkerXYZ(tmpOptFileName, 1, tmpAtomNumber1, 
                        tmpAtomNumber2);
                tmpTinkerXyz.setElementList1(tmpElementList1);
                tmpTinkerXyz.setElementList2(tmpElementList2);
                tmpFileName = tmpIEResultDirName 
                        + FILESEPARATOR 
                        + "output_optimized.xyz";
                tmpTinkerXyz.writeToXyzFile(tmpFileName);
                String tmpKeyName = "\"" 
                        + scratchDirectory
                        + FILESEPARATOR 
                        + tmpParticlePair
                        + ".key"
                        + "\"";
                
                if (tmpParticleName1.equals("H2O") || 
                        tmpParticleName2.equals("H2O")){
                    tmpHasH2O = true;
                } else {
                    tmpHasH2O = false;
                }
                MIPETUTIL.callXYZPDB(tinkerXYZPdb, 
                        tmpOptFileName, 
                        tmpKeyName, 
                        tmpHasH2O);
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
                        tmpBW.append(String.format("%.4f", tmpMinDistance));
                    } catch(IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                            "IOException during writing file in OptDist directory.", 
                            ex);
                    }
                }
                tmpEnergyCalcTime  = (System.currentTimeMillis() 
                        - tmpEnergyCalcTime) / 1000;
                try {
                    BWParticleLog.append("Time to calculate minimum intermolecular energy [s]: " 
                            + tmpEnergyCalcTime);        
                    BWParticleLog.append(LINESEPARATOR); 
                    BWParticleLog.close();
                    BWParticleDat.append("Weighted MinimumIntermolecularEnergy [kcal/mole]: ");
                    BWParticleDat.append(String.format("%.4f", tmpMinEnergy));
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("GlobalMinimumIntermolecularEnergy [kcal/mole]: ");
                    BWParticleDat.append(String.format("%.4f", 
                            tmpEnergyRecords[2].minEnergy()));
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("Optimized minimumIntermolecularEnergy [kcal/mole]: ");
                    BWParticleDat.append(String.format("%.4f", tmpOptMinEnergy));
                    BWParticleDat.append(LINESEPARATOR);
                    BWParticleDat.append("Time to calculate minimum intermolecular energy [s]: "
                            + tmpEnergyCalcTime
                            + LINESEPARATOR);
                    BWParticleDat.append("SphereNodeNumber: "
                            + Integer.toString(sphereNodeNumber)
                            + LINESEPARATOR);
                    BWParticleDat.append("CircularRotationNumber: "
                            + Integer.toString(rotationNumber)
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
        
        System.out.println("Calculating coordination numbers...");
        //<editor-fold defaultstate="collapsed" desc="Coordination numbers">
        if (!forcefield_CN.isEmpty()) {
            getCoordinationNumbers(tmpJobTaskRecordList);
        }
        
        //</editor-fold>
        
        System.out.println("Exporting parameterset...");
        //<editor-fold defaultstate="collapsed" desc="Export parameterset">
        long tmpTotalTimeMins;
        double tmpTotalTimeHours;
        tmpTotalTimeMins = (System.currentTimeMillis() - tmpTotalTime) / 60000;
        tmpTotalTimeHours = (double)tmpTotalTimeMins / 60;
        try {
            BFGblLog.append(LINESEPARATOR);
            BFGblLog.append("Entire calculation Time: " + tmpTotalTimeMins 
                    + String.format("min (%.2f h)" , tmpTotalTimeHours));
            BFGblLog.close();
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, "IOException during writing log file.", 
                    ex);
        }
        
        exportParticleSetForMFSim(tmpJobTaskRecordList, 
                parameterSetTitle, 
                parameterSetTitleAbr);
        
        //</editor-fold>
        System.out.println("Ready.");
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
        
        LinkedList<String> tmpNewParticles = new LinkedList<>();
        LinkedList<String> tmpOldParticles = new LinkedList<>();
        
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
                        
                        for (int i = 0; i < tmpString.length - 1; i++) {
                            extraForcefield_IE.add(tmpString[i + 1]);
                        }
                        
                        continue;
                    }
                    case '$' -> {
                        String[] tmpString = tmpRestString.split("\\s+");
                        forcefield_CN = tmpString[0];
                        
                        for (int i = 0; i < tmpString.length - 1; i++) {
                            extraForcefield_CN.add(tmpString[i + 1]);
                        }
                        
                        continue;
                    }
                    case '-' -> {
                        tmpOldParticles.add(tmpRestString);
                        continue;
                    }
                    default -> tmpNewParticles.add(tmpLine.trim());
                }
            }
            newParticles = tmpNewParticles;
            if (tmpNewParticles.isEmpty()) {
                LOGGER.log(Level.SEVERE, 
                    "No new particles in the job file listed.");
            }
            oldParticles = tmpOldParticles;
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
        String tmpTinkerDirectory;
        String tmpIsFibonacciSphereAlgorithm;
        
        if (MIPETUTIL.getResourceString("MIPET.Tinker9").toLowerCase().
                equals("true")) {
            isTinker9 = true;
        } else {
            isTinker9 = false;
        }
        jobFileName = MIPETUTIL.getResourceString("MIPET.File.jobfile");
        cpuCoreNumber = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETCPUCoreNumber"));
        tmpIsFibonacciSphereAlgorithm = MIPETUTIL.getResourceString(
                "MIPETJavaUseFibonacciSphereAlgorithm");
        isFibonacciSphereAlgorithm = tmpIsFibonacciSphereAlgorithm
                .equals("true");
        temperature = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETTemperature"));
        boltzmannFraction = Double.parseDouble(MIPETUTIL.getResourceString(
                "MIPETBoltzmannFraction"));
        sphereNodeNumber = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETSphereNodeNumber"));
        rotationNumber = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETRotationNumber"));
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
        tmpTinkerDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.tinker");
        if (OS.matches("Windows [0-9]*")) {
            tinkerAnalyze = tmpTinkerDirectory + "/analyze.exe";
            tinkerScan = tmpTinkerDirectory + "/scan.exe";
            tinkerXYZEdit = tmpTinkerDirectory + "/xyzedit.exe";
            tinkerMinimize = tmpTinkerDirectory + "/minimize.exe";
            tinkerDynamic = tmpTinkerDirectory + "/dynamic.exe";
            tinkerOptimize = tmpTinkerDirectory + "/optimize.exe";
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
            tinkerXYZPdb = tmpTinkerDirectory + "/xyzpdb";
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
     * scanParticle method
     *  Optimize, scan and optimize the input particle
     */
    private static void scanParticle(
            ArrayList<JobTaskRecord> aJobTaskRecordList) {
        
        String tmpForcefield;
        String tmpParticle;
        String tmpTaskCandidate;
        String tmpKeyFileString;
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
        HashSet<String> tmpTaskName;
        Path tmpOptXyzDir;
        Path tmpOptXyzFile;
        Path tmOptArcFile;
        Path tmpOriginal;
        Path tmpTarget;
        TinkerXYZ tmpTinkerXYZ;
        TinkerXYZ tmpTinkerXYZ0;
        TinkerXYZ tmpAfterScan;
        Process tmpProcess;
        ProcessBuilder tmpPB;
        List<Double> tmpEnergyValues;
        List<String> tmpMinimumList;
        Double[] tmpEnergyList;
        Integer[] tmpEnergyIndices;
        double[][][] tmpCoords;
        ArrayIndexComparator tmpComparator;
        int tmpJobTaskRecordListSize;
        
        tmpProcess = null;
        tmpTaskName = new HashSet<>();
        tmpJobTaskRecordListSize = aJobTaskRecordList.size();
        
        for (int i = 0; i < tmpJobTaskRecordListSize; i++) {
            
            for (int j = 0; j < 2; j++) {
                if (j == 0) {
                    tmpForcefield = aJobTaskRecordList.get(i)
                            .forcefield_IE_Name();
                } else {
                    tmpForcefield = aJobTaskRecordList.get(i)
                            .forcefield_CN_Name();
                }
                
                if (!tmpForcefield.isEmpty()) {
                    
                    for (int k = 0; k < 2; k++) {
                        if ( k == 0) {
                            tmpParticle = aJobTaskRecordList.get(i)
                                    .particleName1();
                        } else {
                            tmpParticle = aJobTaskRecordList.get(i)
                                    .particleName2();
                        }
                        tmpTaskCandidate = tmpForcefield + "_" + tmpParticle;
                        if (!tmpTaskName.contains(tmpTaskCandidate)) {
                            tmpTaskName.add(tmpTaskCandidate);
                        }
                    }
                    
                }
            }
            
        }
        
        for (String tmpForcefield_Particle : tmpTaskName) {
            tmpForcefield = tmpForcefield_Particle.split("_")[0];
            tmpParticle = tmpForcefield_Particle.split("_")[1];
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
                tmpOriginal = Paths.get(sourceDirectory, tmpForcefield, 
                        tmpParticleXyzName);
                tmpTarget = Paths.get(tmpOptXyzDirName, tmpParticleXyzName);
                try {
                    Files.copy(tmpOriginal, tmpTarget, StandardCopyOption
                            .REPLACE_EXISTING);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during writing .xyz file.", ex);
                }
                tmpTinkerXYZ = new TinkerXYZ(tmpOptXyzDirName 
                        + FILESEPARATOR 
                        + tmpParticleXyzName);
                tmpKeyPathName = tmpOptXyzDirName 
                        + FILESEPARATOR 
                        + tmpParticle 
                        + ".key";
                tmpKeyFileString = "# Force Field Selection"
                + LINESEPARATOR 
                + "PARAMETERS\t\""
                + parameterDirectory
                + "/"
                + tmpForcefield.toLowerCase() + ".prm\""
                + LINESEPARATOR
                + "DIELECTRIC\t" 
                + dielectricConstant;
                        
                // Write .key file
                MIPETUTIL.saveKeyFile(parameterDirectory,
                        tmpForcefield,
                        tmpParticle,
                        tmpKeyPathName, 
                        tmpKeyFileString);
                
                // Create optimized xyz files for the particles via Tinker optimize.
                try {
                    tmpProcess = new ProcessBuilder(tinkerOptimize,
                        tmpOptXyzDirName
                        + FILESEPARATOR 
                        + tmpParticleXyzName,
                        String.valueOf(optimizeRmsGradient))
                        .start();
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE,
                            "IOException during process start.", ex);
                } finally {
                    if (tmpProcess != null) {
                        try {
                            InputStream tmpInput = tmpProcess.getInputStream();
                            tmpOptXyzLogName = tmpOptXyzDirName
                                    + FILESEPARATOR 
                                    + tmpParticle
                                    + "_preoptimize.log";
                            OutputStream tmpOutput = new FileOutputStream(
                                    tmpOptXyzLogName, false);
                            tmpInput.transferTo(tmpOutput);
                            tmpInput.close();
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
                tmpTinkerXYZ0 = new TinkerXYZ(tmpOptXyzDirName 
                        + FILESEPARATOR 
                        + tmpParticle 
                        + ".xyz_2");
                tmpTinkerXYZ.setCoordinateList1(tmpTinkerXYZ0
                        .getCoordinateList1());
                tmpOriginal = Paths.get(tmpOptXyzDirName, tmpParticle 
                        + ".xyz_2");
                try {
                    Files.deleteIfExists(tmpOriginal);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during deleting .xyz_2 file", ex);
                }
                tmpXyzFileName = tmpOptXyzDirName 
                        + FILESEPARATOR 
                        + tmpParticleXyzName;
                try (BufferedWriter tmpBW = new BufferedWriter(
                            new FileWriter(tmpXyzFileName))) {
                    tmpBW.append(tmpTinkerXYZ.getFileContent());
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE,
                            "IOException during writing .xyz file.", ex);
                }
                
                // Run Tinker Scan.exe with optimized particles
                tmpPB = new ProcessBuilder();
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
                        tmpTinkerXYZ.setCoordinateList1(tmpCoords[j]);
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
                            tmpEnergyValues.add(Double.valueOf(tmpMinimumList
                                    .get(j)
                                    .substring(58, 68)));
                        }

                        tmpEnergyList = new Double[tmpEnergyValues.size()];

                        for (int j = 0; j < tmpEnergyValues.size(); j++) {
                            tmpEnergyList[j] = tmpEnergyValues.get(j);
                        }

                        tmpComparator = new ArrayIndexComparator(tmpEnergyList);
                        tmpEnergyIndices = tmpComparator.createIndexArray();
                        Arrays.sort(tmpEnergyIndices, tmpComparator);
                        tmpOptXyzName = tmpOptXyzDirName 
                                + FILESEPARATOR
                                + tmpParticle 
                                + "_0.txyz";
                        tmpAfterScan = new TinkerXYZ(tmpOptXyzName); 
                        tmpTinkerXYZ.setCoordinateList1(tmpAfterScan
                                .getCoordinateList1());

                        // Convert tinker xyz files to .xyz file
                        for (int j = 0; j < tmpEnergyValues.size(); j++) {
                            tmpFileTxyzName = tmpOptXyzDirName 
                                + FILESEPARATOR
                                + tmpParticle
                                + "_" + j +".txyz";
                            tmpTinkerXYZ0 = new TinkerXYZ(tmpFileTxyzName);
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
                        tmpOriginal = Paths.get(tmpSourceName);
                        tmpTarget = Paths.get(tmpTargetName);
                        try {
                            if (Files.exists(tmpOriginal)) {
                                Files.copy(tmpOriginal, tmpTarget,
                                        StandardCopyOption.REPLACE_EXISTING);
                            }
                        } catch (IOException ex) {
                            LOGGER.log(Level.SEVERE, "IOException during copying .xyz file", ex);
                        }

                        // Optimize after scan.exe
                        try {
                            tmpProcess = new ProcessBuilder(tinkerOptimize,
                                    tmpOptXyzDirName
                                    + FILESEPARATOR 
                                    + tmpParticleXyzName,
                                    String.valueOf(optimizeRmsGradient))
                                    .start();
                            try (InputStream tmpInput = tmpProcess
                                    .getInputStream()) {
                                tmpOptXyzLogName = tmpOptXyzDirName
                                    + FILESEPARATOR 
                                    + tmpParticle
                                    + "_afteroptimize.log";
                                OutputStream tmpOutput = new FileOutputStream(
                                        tmpOptXyzLogName, false);
                                tmpInput.transferTo(tmpOutput);
                                tmpOutput.close();
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
                        tmpTinkerXYZ0 = new TinkerXYZ(tmpOptXyzDirName 
                                + FILESEPARATOR 
                                + tmpParticle 
                                + ".xyz_2");
                        tmpTinkerXYZ.setCoordinateList1(tmpTinkerXYZ0
                                .getCoordinateList1());
                        tmpOriginal = Paths.get(tmpOptXyzDirName, 
                                tmpParticle + ".xyz_2");
                        try {
                            Files.deleteIfExists(tmpOriginal);
                        } catch (IOException ex) {
                            LOGGER.log(Level.SEVERE, "IOException during deleting .xyz_2 file.", ex);
                        }
                        tmpXyzFileName = tmpOptXyzDirName 
                                + FILESEPARATOR 
                                + tmpParticleXyzName;
                        try (BufferedWriter tmpBW = new BufferedWriter(
                                    new FileWriter(tmpXyzFileName))) {
                            tmpBW.append(tmpTinkerXYZ.getFileContent());
                        } catch(IOException ex) {
                            LOGGER.log(Level.SEVERE,
                                    "IOException during writing .xyz file.", ex);
                        }
                    } 
                }
            }
        }

    }
    
    /**
     * Calculate the intermolecular energies for different distances.
     * 
     * @param aParticlePair Names of both particles
     * @param aDistance Distances of both particles (centre to centre)
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
            double[] aDistance,
            TinkerXYZ aTinkerXYZ1, 
            TinkerXYZ aTinkerXYZ2,
            double[][][] aRotData1,
            double[][][] aRotData2,
            double aMinEnergy) {
        
        int tmpDistanceNumber;
        int tmpConfigNumber;
        int tmpChunkNumber;
        int tmpChunkSize;
        int tmpChunkIndex;
        int tmpPart1StartIndex;
        int tmpConfigIndex;
        int tmpRotData1Index;
        int tmpRotData2Index;
        double[][][] tmpPartRotData1;
        double[][][] tmpRotData2;
        String tmpPath;
        String[] tmpCmdList;
        TinkerXYZ tmpTinkerXYZ;

        tmpDistanceNumber = aDistance.length;
        tmpChunkNumber = cpuCoreNumber;
        tmpConfigNumber = aRotData1.length * aRotData2.length;
        if (tmpDistanceNumber >= tmpChunkNumber || tmpConfigNumber < 1000) {
            tmpChunkNumber = 1;
            tmpChunkSize = tmpConfigNumber;
        } else {
            tmpChunkNumber = (int)Math.ceil((double)tmpChunkNumber / 
                    tmpDistanceNumber);
            tmpChunkSize = (int)Math.ceil((double)tmpConfigNumber 
                    / tmpChunkNumber);
            if(tmpChunkSize % aRotData2.length != 0) {
                int tmpResidual = Math.floorMod(tmpChunkSize, aRotData2.length);
                tmpChunkSize += (aRotData2.length - tmpResidual);
            }
        }
        
        // return values
        double[][] tmpEnergyDatas = new double[tmpDistanceNumber][];
        ExecutorService executor = Executors.newFixedThreadPool(cpuCoreNumber);
        LinkedList<MIPETMakeArc> tmpTaskListMakeArc = new LinkedList<>();
        tmpTinkerXYZ = new TinkerXYZ(aTinkerXYZ1, aTinkerXYZ2);
        tmpPath = scratchDirectory 
                + FILESEPARATOR 
                + aParticlePair 
                + ".arc";
        
        for (int i = 0; i < tmpDistanceNumber; i++) {
            tmpRotData2 = VectorUtil.moveX(aRotData2, aDistance[i]);
            tmpChunkIndex = 0;
            tmpConfigIndex = 0;
            tmpRotData1Index = -1;
            tmpPart1StartIndex = 0;
            
            while (tmpRotData1Index < aRotData1.length - 1) {
                tmpRotData1Index++;
                tmpRotData2Index = -1;
                
                while (tmpRotData2Index < aRotData2.length - 1) {
                    tmpRotData2Index++;
                    if ((tmpConfigIndex + 1) % tmpChunkSize == 0 
                            || tmpConfigIndex + 1 == tmpConfigNumber) {
                        tmpPartRotData1 = Arrays.copyOfRange(aRotData1, 
                                tmpPart1StartIndex,
                                tmpRotData1Index + 1);
                        tmpTaskListMakeArc.add(new MIPETMakeArc(
                                aParticlePair,
                                tmpPath + i + "_" + tmpChunkIndex,
                                aTinkerXYZ1,
                                aTinkerXYZ2,
                                tmpPartRotData1,
                                tmpRotData2,
                                tmpChunkSize,
                                minAtomDistance));
                        tmpPart1StartIndex = tmpRotData1Index + 1;
                        tmpChunkIndex++;
                    }
                    tmpConfigIndex++;
                }
            }
            
        }
        
        try {
            executor.invokeAll(tmpTaskListMakeArc);
        } catch (InterruptedException ex) {
            Thread.currentThread().interrupt();       
        }
        executor.shutdown();
        
        // Calculate intermolecular energy using TINKER analyze
        executor = Executors.newFixedThreadPool(cpuCoreNumber);
        LinkedList<MIPETAnalyze> tmpTaskList = new LinkedList<>();
        
        for (int i = 0; i < tmpDistanceNumber; i++) {

            for (int j = 0; j < tmpChunkNumber; j++) {
                tmpCmdList = new String[]{tinkerAnalyze, 
                    tmpPath + i + "_"+ j, "E"};
                tmpTaskList.add (new MIPETAnalyze (i,
                        j,
                        scratchDirectory, 
                        aParticlePair, 
                        tmpCmdList));
            }
            
        }
        
        try {
            executor.invokeAll(tmpTaskList);
        } catch (InterruptedException ex) {
            Thread.currentThread().interrupt();       
        }
        executor.shutdown();
        
        // Read intermolecular energy
        ArrayList<Double> tmpEnergyList;
        double[] tmpEnergySorted;
        double tmpValue;
        String tmpLine;
        String tmpSearch;
        String tmpValueCandidate;
        
        tmpPath = scratchDirectory 
                + FILESEPARATOR 
                + aParticlePair;
        tmpSearch = "Intermolecular Energy";
        double tmpMinEnergy = 1E10;
        int tmpDistMinIndex = -1;
        int tmpChunkMinIndex = -1;
        int tmpMinIndex = -1;
        
        for (int i = 0; i < tmpDistanceNumber; i++) {
            tmpEnergyList = new ArrayList<>(tmpConfigNumber);

            for (int j = 0; j < tmpChunkNumber; j++) {
                tmpConfigIndex = -1;
                try (BufferedReader tmpBR = new BufferedReader(new FileReader(
                        tmpPath + ".out" + i + "_" + j))) {
                    while ((tmpLine = tmpBR.readLine()) != null) {
                        if (tmpLine.contains(tmpSearch)) {
                            tmpValueCandidate = tmpLine.substring(25, 50);
                            if (!tmpValueCandidate.contains("D")) {
                                tmpValue = Double
                                       .parseDouble(tmpValueCandidate);
                                if (tmpValue != Double.NaN && 
                                        tmpValue != Double.POSITIVE_INFINITY) {
                                    tmpEnergyList.add(tmpValue);
                                    tmpConfigIndex++;
                                    if (tmpValue < tmpMinEnergy) {
                                        tmpMinEnergy = tmpValue;
                                        tmpMinIndex = tmpConfigIndex;
                                        tmpChunkMinIndex = j;
                                        tmpDistMinIndex = i;
                                    }
                                }
                            }
                        }
                    }
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE,
                            "IOException during reading .out file." , ex);
                }
            }
            
            tmpEnergyList.sort(Comparator.naturalOrder());
            tmpEnergySorted = new double[tmpEnergyList.size()];
            
            for (int j = 0; j < tmpEnergySorted.length; j++) {
                tmpEnergySorted[j] = tmpEnergyList.get(j);
            }
                    
            tmpEnergyDatas[i] = tmpEnergySorted.clone();
        }
       
        // Export .xyz file with lowest intermolecular energy
        if (aMinEnergy > tmpMinEnergy) {
            int tmpStartIndex = tmpMinIndex * (tmpTinkerXYZ.getAtomNumber() 
                    + 1);
            int tmpEndIndex = tmpStartIndex + tmpTinkerXYZ.getAtomNumber();
            String tmpFileName = tmpPath + ".arc" + tmpDistMinIndex + "_" 
                    + tmpChunkMinIndex;
            StringBuilder tmpPartArc = MIPETUTIL.readPartArcFile(tmpFileName, 
                    tmpStartIndex, tmpEndIndex);
            tmpFileName = tmpPath + ".0";
            tmpTinkerXYZ.writeToXyzFile(tmpFileName, tmpPartArc);
        }

        // Clean Scratch directory
        try (Stream<Path> tmpList = Files.list(Paths.get(scratchDirectory))) {
            tmpList.map(Path::toFile)
                    .filter(file -> file.toString().contains(".arc") ||
                            file.toString().contains(".out"))
                    .forEach(File::delete);
        } catch(IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during deleting scratch directory.", ex);
        }
        return new EnergyRecord(aDistance,
                tmpEnergyDatas, 
                aDistance[tmpDistMinIndex], 
                tmpDistMinIndex, tmpMinEnergy);
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
        
        for (int i = 0; i < tmpJobTaskLength; i++) {
            if (aJobTaskRecordList.get(i).hasCNJob()) {
                tmpParticle1 = aJobTaskRecordList.get(i).particleName1();
                tmpParticle2 = aJobTaskRecordList.get(i).particleName2();
                tmpParticlePair = tmpParticle1 + "_" + tmpParticle2;
                tmpIsSameParticle = aJobTaskRecordList.get(i)
                        .isSameParticle();

                // Copy .xyz files to scratch directory
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
                tmpSource = Paths.get(optXYZDirectory
                        + FILESEPARATOR
                        + forcefield_CN
                        + FILESEPARATOR 
                        + tmpParticle1
                        + FILESEPARATOR 
                        + tmpParticle1
                        + ".xyz");
                tmpTarget = Paths.get(scratchDirectory
                        + FILESEPARATOR
                        + forcefield_CN
                        + FILESEPARATOR 
                        + tmpParticlePair
                        + FILESEPARATOR 
                        + tmpParticle1
                        + ".xyz");
                try {
                    Files.copy(tmpSource, tmpTarget, 
                            StandardCopyOption.REPLACE_EXISTING);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during copying files to scratch.", ex);
                }
                if (!tmpIsSameParticle) {
                    tmpSource = Paths.get(optXYZDirectory
                            + FILESEPARATOR
                            + forcefield_CN
                            + FILESEPARATOR 
                            + tmpParticle2
                            + FILESEPARATOR 
                            + tmpParticle2
                            + ".xyz");
                    tmpTarget = Paths.get(scratchDirectory
                            + FILESEPARATOR
                            + forcefield_CN
                            + FILESEPARATOR 
                            + tmpParticlePair
                            + FILESEPARATOR 
                            + tmpParticle2
                            + ".xyz");
                    try {
                        Files.copy(tmpSource, tmpTarget, 
                                StandardCopyOption.REPLACE_EXISTING);
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
            "EWALD"
            + LINESEPARATOR
            + "OPENMP-THREADS 1"
            + LINESEPARATOR
            + "THERMOSTAT ANDERSEN"
            + LINESEPARATOR
            + "STEEPEST-DESCENT"
            + LINESEPARATOR
            + "RANDOMSEED 123456"
            + LINESEPARATOR;
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
        
        tmpJobIndex = 0;

        for (int i = 0; i < tmpJobTaskLength; i++) {
            if (aJobTaskRecordList.get(i).hasCNJob()) {
                tmpHasCNHeadLine = false;
                tmpParticle1 = aJobTaskRecordList.get(i)
                        .particleName1(); // solute
                tmpParticle2 = aJobTaskRecordList.get(i)
                        .particleName2(); // solvent
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
                tmpKeyContent = tmpKeyFixContent
                        + "PARAMETERS\t\""
                        + parameterDirectory
                        + "/"
                        + forcefield_CN.toLowerCase()
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
                MIPETUTIL.saveKeyFile(parameterDirectory,
                        forcefield_CN,
                        tmpParticle1,
                        tmpParticle2,
                        tmpKeyFileName,
                        tmpKeyContent);

                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="Build solvent box">
                tmpProcessBuilder = new ProcessBuilder(tinkerXYZEdit,
                        tmpCurrentDir + tmpParticle2 + ".xyz",
                        "-k",
                        tmpKeyFileName,
                        "23",  // Option 23: Create and Fill a Periodic Boundary Box
                        String.valueOf(solventMoleculeNumber),
                        String.valueOf(tmpBoxLengths[tmpJobIndex]),
                        String.valueOf(tmpBoxLengths[tmpJobIndex]),
                        String.valueOf(tmpBoxLengths[tmpJobIndex]),
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
                                new InputStreamReader(tmpProcess.getInputStream()))) {
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
                        String.valueOf(rmsMinimizeGradient));
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
                    BWParticleLog.append(String
                            .valueOf((System.currentTimeMillis() 
                                    - boxMinimizationTime) / 1000));
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
                        String.valueOf(warmUpStepNumber),
                        String.valueOf(warmUpTimeStep),
                        String.valueOf(warmUpPrintInterval),
                        String.valueOf(simulationType),
                        String.valueOf(temperature)};
                } else {
                    tmpCmdList = new String[] {tinkerDynamic,
                        tmpCurrentDir + tmpParticlePair + ".xyz",
                        String.valueOf(warmUpStepNumber),
                        String.valueOf(warmUpTimeStep),
                        String.valueOf(warmUpPrintInterval),
                        String.valueOf(simulationType),
                        String.valueOf(temperature)};
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
        tmpLabelValues[0][1] = String.valueOf((System
                .currentTimeMillis() - tmpDynamicWarmUpTime) / 1000);
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
                        String.valueOf(stepNumber),
                        String.valueOf(timeStep),
                        String.valueOf(printInterval),
                        String.valueOf(simulationType),
                        String.valueOf(temperature)};
                } else {
                    tmpCmdList = new String[] {tinkerDynamic,
                        tmpCurrentDir + tmpParticlePair + ".xyz",
                        String.valueOf(stepNumber),
                        String.valueOf(timeStep),
                        String.valueOf(printInterval),
                        String.valueOf(simulationType),
                        String.valueOf(temperature)};
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
            tmpLabelValues[0][1] = String.valueOf((System
                    .currentTimeMillis() - tmpArcEvaluationTime) / 1000);
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
                    tmpValues[tmpJobIndex] = String
                            .valueOf(tmpCNMeans[tmpJobIndex]);
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
        String tmpParticle;
        String tmpResultsDirectory;
        
        tmpJobLength = aJobTaskRecords.size();
        MIPETUtility tmpUtility = new MIPETUtility();
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
            tmpParticleDescriptions.add(String.valueOf(MASS_DPD)); //massDPD
            tmpParticleDescriptions.add(String.valueOf(DEFAULT_CHARGE)); // charge
            try {
                tmpParticleDescriptions.add(String.format("%.4f", 
                        tmpUtility.getAtomicMass(smiles.get(
                                tmpParticleName)))); // mass [g/mol]
            } catch (Exception ex) {
                tmpParticleDescriptions.add("-1");
                LOGGER.log(Level.SEVERE, ex.toString());
            } 
            tmpParticleDescriptions.add(String.format("%.4f", 
                    tmpUtility.getVdwVolume(smiles.get(tmpParticleName)) 
                            * WATERVOLUMERATIO)); // volume [A³]
            tmpParticleDescriptions.add(String.valueOf(GRAPHICS_RADIUS)); // graphics-radius
            tmpParticleDescriptions.add(DEFAULT_COLOR); // color
            tmpParticleDescriptionsList.add(tmpParticleDescriptions);
        }
        
        // particle interactions
        Double tmpEnergy;
        Double tmpCN;
        String tmpParticleName1;
        String tmpParticleName2;
        String tmpParticlePair;
        String tmpFileName;
        int tmpEnergyListLength;
        int tmpCNListLength;
        
        tmpEnergyListLength = energyList.size();
        tmpCNListLength = cnList.size();
        HashMap<String, Double> tmpEnergieMap = 
                new HashMap<>(tmpEnergyListLength);
        HashMap<String, Double> tmpCNMap = new HashMap<>(tmpCNListLength);
        
        // Read energy data
        for (int i = 0; i < tmpEnergyListLength; i++) {
            tmpParticleName1 = energyList.get(i).particleName1();
            tmpParticleName2 = energyList.get(i).particleName2();
            tmpParticlePair = tmpParticleName1 + "_" + tmpParticleName2;
            tmpEnergy = energyList.get(i).energyValue();
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
        double tmpChiNumerator;
        double tmpE12;
        double tmpE11;
        double tmpE22;
        double tmpZ11;
        double tmpZ22;
        double tmpZ12;
        double tmpZ21;
        double tmpAij;
        
        Set<String> tmpKeySet;
        HashMap<String, Double> tmpAijMap = new HashMap<>(tmpJobLength);
        HashMap<String, Double> tmpAijMap1 = new HashMap<>(tmpJobLength);
        
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
                tmpAij = (double)temperature / 12 + 1.7483 * tmpChiNumerator;
                tmpAijMap1.put(tmpParticleName1 + "_" + tmpParticleName2, 
                        tmpAij);
            }
        }

        /* Write file */
        for (int i = 0; i < 2; i++) {
            if (i == 0) {
                tmpFileName = tmpResultsDirectory 
                + FILESEPARATOR 
                + aTitleAbbreviation 
                + "_"
                +"EijFraction_" 
                + boltzmannFraction
                + "_catchRadius_" 
                + catchRadius 
                + ".txt";
            } else {
                tmpFileName = tmpResultsDirectory 
                + FILESEPARATOR 
                + aTitleAbbreviation 
                + "_"
                +"EijFraction_" 
                + boltzmannFraction
                + "_catchRadius_" 
                + catchRadius 
                + "_CN1.txt";
            }
            try (BufferedWriter tmpBW = new BufferedWriter(
                new FileWriter(tmpFileName))) {
                tmpBW.append("# Particle set for MFSim created by MIPET\n");
                tmpBW.append("# Force Field for energy calculation: ");
                tmpBW.append(forcefield_IE);
                tmpBW.append(LINESEPARATOR);
                if (i == 0) {
                    tmpBW.append("# Force Field for coordination number calculation: ");
                    tmpBW.append(forcefield_CN);
                    tmpBW.append(LINESEPARATOR);
                }
                tmpBW.append("# Water model: TIP5P");
                tmpBW.append(LINESEPARATOR);
                tmpBW.append("# CPU cores: ");
                tmpBW.append(String.valueOf(cpuCoreNumber));
                tmpBW.append(LINESEPARATOR);
                tmpBW.append("# Temperature: ");
                tmpBW.append(String.valueOf(temperature));
                tmpBW.append(LINESEPARATOR);
                tmpBW.append("# Sphere nodes calculated with Fibonacci algorithm: ");
                tmpBW.append(String.valueOf(isFibonacciSphereAlgorithm));
                tmpBW.append(LINESEPARATOR);
                tmpBW.append("# Sphere node number: ");
                tmpBW.append(String.valueOf(sphereNodeNumber));
                tmpBW.append(LINESEPARATOR);
                tmpBW.append("# Sphere rotation number: ");
                tmpBW.append(String.valueOf(rotationNumber));
                tmpBW.append(LINESEPARATOR);
                tmpBW.append("# Fraction for Boltzmann averaging: ");
                tmpBW.append(String.valueOf(boltzmannFraction));
                tmpBW.append(LINESEPARATOR);
                if (i == 1) {
                    tmpBW.append("# CN = 1 for all particle pairs");
                    tmpBW.append(LINESEPARATOR);
                }
                tmpBW.append(LINESEPARATOR);
                tmpBW.append("[Title]\n" + aTitle + "\n[/Title]\n\n");
                tmpBW.append("[Version]\n" 
                        + VERSION_NUMBER 
                        + "\n[/Version]\n\n");
            
                // Particle description
                tmpBW.append("[Particle Description]");
                tmpBW.append(LINESEPARATOR);

                for(String tmpHeader : tmpParticleDescriptionString){
                    tmpBW.append(tmpHeader + "\t");
                }

                for(List<String> tmpParticleDescList : tmpParticleDescriptionsList){
                    tmpBW.append(LINESEPARATOR);

                    for(String tmpDescription : tmpParticleDescList){
                        tmpBW.append(String.valueOf(tmpDescription));
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
                tmpBW.append("# Repulsion parameters a(ij) for the temperature (in K)");
                tmpBW.append(LINESEPARATOR);
                tmpBW.append("Pair");
                tmpBW.append(" " + String.valueOf((int)temperature));
                if (i == 0) {
                    tmpKeySet = tmpAijMap.keySet();
                } else {
                    tmpKeySet = tmpAijMap1.keySet();
                }

                for(String tmpKey : tmpKeySet) {
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append(tmpKey);
                    if (i == 0) {
                        tmpAij = tmpAijMap.get(tmpKey);
                    } else {
                        tmpAij = tmpAijMap1.get(tmpKey);
                    }
                    tmpBW.append(" " +  String.format("%.2f",tmpAij));
                }

                tmpBW.append(LINESEPARATOR);
                tmpBW.append("[/Particle interactions]");
                tmpBW.append(LINESEPARATOR);
                tmpBW.append(LINESEPARATOR);
                
                // Coordination numbers
                if (!forcefield_CN.isEmpty()) {
                    if (i == 0) {
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
                        tmpBW.append(LINESEPARATOR);
                    }
                }

                // Particle SMILES
                tmpBW.append(LINESEPARATOR);
                tmpBW.append("[Particle SMILES]");
                tmpBW.append(LINESEPARATOR);
                tmpBW.append("# Particle\tSMILES");
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
                tmpBW.append("[/Particle SMILES]");
            } catch (IOException anException) {
                LOGGER.log(Level.SEVERE, anException.toString());
            }
        }
    }
}
// </editor-fold>