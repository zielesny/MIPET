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
import java.io.IOException;
import java.io.InputStream;
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
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.apache.commons.io.FileUtils;

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
    private static List<String> particleNames;
    
    /**
     * New particles for calculation
     */
    private static List<String> newParticles;
    
    /**
     * Old particles already calculated
     */
    private static List<String> oldParticles;
    
    /** 
     * MoleculeRecord
     */
    private static List<MoleculeRecord> molecules;
    
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
     * 
     * @param args Input arguments
     */
    public static void main(String[] args) {
        //<editor-fold defaultstate="collapsed" desc="Variables declaration">
        ProcessBuilder builder;
        ArrayList<JobTaskRecord> jobTaskRecordList = new ArrayList<>(500);

        //</editor-fold>
        
        // <editor-fold defaultstate="collapsed" desc="Initialize and read .job file">
        Locale.setDefault(Locale.ENGLISH);
        Instant startTime = Instant.now();
        System.out.println("Initializing...");
        initialize(args);
        System.out.println("Reading job file...");
        try {
            readJobFile();
            System.out.println("  The job file is OK."  );
        } catch (IllegalArgumentException | IllegalStateException ex) {
            LOGGER.log(Level.SEVERE, "Critical error during reading job file: " 
                    + ex.getMessage(), ex);
            System.err.println("Please verify that all .xyz and .prm files listed in the job file exist.");
        }
        
        //</editor-fold>
        
        // <editor-fold defaultstate="collapsed" desc="Set original .key file">
        /*
         * Standard .key file content
         */
        String keyFileOrigin =
                "# Force Field Selection"
                + LINESEPARATOR
                + "PARAMETERS\t\""
                + parameterDirectory
                + "/";
        if (forceField_IE != null) {
            if ("OPLSAALIGPARGEN".equals(forceField_IE)) {
                keyFileOrigin +=
                        "oplsaa.prm\""
                        + LINESEPARATOR;
            } else {
                keyFileOrigin +=
                        forceField_IE + ".prm\""
                        + LINESEPARATOR;
            }
        }
        
        // </editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Create log file">
        Path gblLogPath;
        Path logPath;
        if (forceField_IE == null) {
            gblLogPath = Paths.get(resultDirectory, "CN", forceField_CN);
            logPath = Paths.get(resultDirectory, "CN", forceField_CN, 
                    "log.txt");
        } else {
            gblLogPath = Paths.get(resultDirectory, "IE", forceField_IE);
            logPath = Paths.get(resultDirectory, "IE", forceField_IE, 
                    "log.txt");
        }
        try {
            Files.createDirectories(gblLogPath);
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                "IOException during creating global log file.", ex);
        }
        boolean fileAlreadyExists = Files.exists(logPath);
        try {
            BFGblLog = Files.newBufferedWriter(logPath, StandardCharsets.UTF_8);
            if (fileAlreadyExists) {
                BFGblLog.newLine();
            }
            BFGblLog.append("CPU cores: ")
                    .append(String.valueOf(cpuCoreNumber));
            BFGblLog.newLine();
            BFGblLog.append("Warm up steps for dynamic simulation: ")
                    .append(String.valueOf(warmUpStepNumber));
            BFGblLog.newLine();
            BFGblLog.append("Steps for dynamic simulation: ")
                    .append(String.valueOf(stepNumber));
            BFGblLog.newLine();
            BFGblLog.flush();
        } catch(IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during writing in log file.", ex);
        }
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Prepare input">
        ArrayList<String> particlePairs = getParticlePairs();
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Make Job task record">
        String particle1;
        String particle2;
        String particlePair;
        boolean isSameParticle;
        int particlePairLength = particlePairs.size();
        
        for (int i = 0; i < particlePairLength; i++) {
            particlePair = particlePairs.get(i);
            particle1 = particlePair.split("_")[0];
            particle2 = particlePair.split("_")[1];
            isSameParticle = particle1.equals(particle2);
            String ieResultDirName = resultDirectory
                    + FILESEPARATOR
                    + "IE"
                    + FILESEPARATOR
                    + forceField_IE
                    + FILESEPARATOR
                    + particlePair;
            String cnResultDirName = resultDirectory
                    + FILESEPARATOR
                    + "CN"
                    + FILESEPARATOR
                    + forceField_CN
                    + FILESEPARATOR
                    + particlePair;
            boolean hasEnergieJob = forceField_IE != null;
            if (hasEnergieJob) {
                try {
                    Files.createDirectories(Paths.get(ieResultDirName));
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE,
                            "IOException during creating result directory.", ex);
                }
            }
            boolean hasCNJob = forceField_CN != null;
            if (hasCNJob) {
                try {
                    Files.createDirectories(Paths.get(cnResultDirName));
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during creating CN result directory.", ex);
                }
            }
            jobTaskRecordList.add(new JobTaskRecord(
                    forceField_IE,
                    forceField_CN,
                    particle1, 
                    particle2, 
                    ieResultDirName,
                    cnResultDirName,
                    hasEnergieJob,
                    hasCNJob,
                    isSameParticle,
                    false));
            if (!isSameParticle) {
                hasEnergieJob = false;
                jobTaskRecordList.add(new JobTaskRecord(
                        forceField_IE,
                        forceField_CN,
                        particle2, 
                        particle1, 
                        ieResultDirName,
                        cnResultDirName,
                        hasEnergieJob,
                        hasCNJob,
                        isSameParticle,
                        true));
            }
        }
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Read .xyz and .prm">
        if (forceField_IE != null) {
            readXyz(true, forceField_IE);
            readPrm(forceField_IE);
            makeMoleculeRecord();
        }
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Conformational analysis">
        String forcefield;
        Path optXyzPath;
        
        if (isConformationalAnalysis) {
            System.out.println("Conformational Analysis...");
            scanParticle();
            
            // Overwrite xyzContent with new .xyz
            readXyz(false, forceField_IE);
        } else {
            System.out.println("Conformational Analysis skipped.");
            if (forceField_IE != null) {
                forcefield = jobTaskRecordList.getFirst().forcefield_IE_Name();

                for (String particle : particleNames) {
                    optXyzPath = Paths.get(optXYZDirectory, forcefield,
                            particle);
                    if (!Files.exists(optXyzPath)) {
                        try {
                            Files.createDirectories(optXyzPath);
                        } catch (IOException ex) {
                            LOGGER.log(Level.SEVERE,
                                    "IOException during creating OptXyz directory.",
                                    ex);
                        }
                    }
                }

                // Copy .xyz file to optXYZ directory
                for (String particleName : particleNames) {
                    // Check whether already scanned
                    Path sourcePath = Paths.get(
                            moleculeDirectory,
                            forcefield,
                            particleName + ".xyz");
                    Path targetPath = Paths.get(
                            optXYZDirectory,
                            forcefield,
                            particleName,
                            particleName + ".xyz");
                    if (!Files.exists(targetPath)) {
                        try {
                            Files.copy(sourcePath, targetPath);
                        } catch (IOException ex) {
                            LOGGER.log(Level.SEVERE,
                                    "IOException during copying .xyz file.", ex);
                        }
                    }
                }
                
            }
        }
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Variable initializations">
        DecimalFormat decimal2 = (DecimalFormat)NumberFormat
                .getNumberInstance();
        DecimalFormat decimal3 = (DecimalFormat)NumberFormat
                .getNumberInstance();
        DecimalFormat decimal4 = (DecimalFormat)NumberFormat
                .getNumberInstance();
        decimal2.applyPattern("#0.00");
        decimal3.applyPattern("#0.000");
        decimal4.applyPattern("#0.0000");
        boolean isExitCondition = false;
        int currIndex = 0;
        String keyFileString = keyFileOrigin 
                + "DIELECTRIC\t" 
                + dielectricConstant
                + LINESEPARATOR;
        forcefield = forceField_IE;
        System.out.println("Calculating intermolecular energy...");
        
        //</editor-fold>
        
        while (!isExitCondition) {
            // Exit condition is true when all particle pair combinations
            //   were calculated.
            if (jobTaskRecordList.get(currIndex).hasEnergieJob()) {
                //<editor-fold defaultstate="collapsed" desc="Job task record">
                particle1 = jobTaskRecordList.get(currIndex).particleName1();
                particle2 = jobTaskRecordList.get(currIndex).particleName2();
                particlePair = particle1 + "_" + particle2;
                isSameParticle = particle1.equals(particle2);
                String ieResultDirName = jobTaskRecordList.get(currIndex)
                        .result_IE_PathName();
                byte pos_H2O = 0;
                if (particle1.equals("H2O") ^ particle2.equals("H2O")) {
                    if (particle1.equals("H2O")) {
                        pos_H2O = 1;
                    } else {
                        pos_H2O = 2;
                    }
                }
                MIPETUtility.updateStatus("  Calculating " + particlePair);
            
                //</editor-fold>
                //<editor-fold defaultstate="collapsed" desc="Create log file">
                Path particleDatPath = Paths.get(ieResultDirName,
                        particlePair + ".dat");
                try {
                    BWParticleDat = Files.newBufferedWriter(particleDatPath, 
                            StandardCharsets.UTF_8);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during writing in log file.", ex);
                }
                
                //</editor-fold>
                //<editor-fold defaultstate="collapsed" desc="Centre first fragment">
                // Centre the fragments and move to the centre
                int xyz1ID = particleNames.indexOf(particle1);
                TinkerXYZ tinkerXYZ1 = new TinkerXYZ(forcefield, particle1,
                        xyzContent1[xyz1ID]);
                TinkerXYZ tinkerXYZ2;
                double[][] data1 = tinkerXYZ1.getCoordinateList1()[0];
                double[] centre1 = tinkerXYZ1.findCenterCoordinate();
                data1 = tinkerXYZ1.moveCoordinates(data1, centre1);
                if (isSameParticle) {
                    tinkerXYZ2  = new TinkerXYZ(forcefield,
                        particle2,
                        xyzContent1[xyz1ID]);
                } else {
                    int xyz2ID = particleNames.indexOf(particle2);
                    tinkerXYZ2 = new TinkerXYZ(forcefield,
                            particle2,
                            xyzContent2[xyz2ID]);
                }
                double[][] xyzData2 = tinkerXYZ2.getCoordinateList1()[0];
                double[] centre2 = tinkerXYZ2.findCenterCoordinate();
                xyzData2 = tinkerXYZ2.moveCoordinates(xyzData2, centre2);

                //</editor-fold>
                //<editor-fold defaultstate="collapsed" desc="Calculate Intermolecular Energy">

                //<editor-fold defaultstate="collapsed" desc="Calculate rotated coordinates">
                RotationRecord rotRecord = RotationUtil.getRotationsCoords(
                        sphereNodeNumber1, 
                        rotationNumber1,
                        data1, 
                        xyzData2,
                        isFibonacciSphereAlgorithm);
                
                // xyz-configurations [i]: configuration, [j]: atom, [k]: xyz
                double[][][] xyzRotData1 = rotRecord.rotData1();
                double[][][] xyzRotData2 = rotRecord.rotData2();
                
                //</editor-fold>
                //<editor-fold defaultstate="collapsed" desc="Prescan">
                long energyCalcTime = System.currentTimeMillis();
                int distSize = (int)Math.ceil((upperBoundary - lowerBoundary) 
                        / prescanStepSize) + 1;
                List<Double> distList = new ArrayList<>();
                List<Distance_EnergyRecord> distEminRecords = new ArrayList<>();
                int prmID1 = particleNames.indexOf(particle1);
                int prmID2;
                if (isSameParticle) {
                    prmID2 = prmID1;
                } else {
                    prmID2 = particleNames.indexOf(particle2);
                }

                /* Write .key file */
                Path keyPath = Paths.get(scratchDirectory,
                        particlePair + ".key");
                StringBuilder keyContent = new StringBuilder(
                        keyFileOrigin)
                        .append("DIELECTRIC\t")
                        .append(dielectricConstant)
                        .append(LINESEPARATOR)
                        .append("NONBONDTERM ONLY")
                        .append(LINESEPARATOR)
                        .append("OPENMP-THREADS 1")
                        .append(LINESEPARATOR);
                
                if ("OPLSAALIGPARGEN".equals(forcefield)) {
                    keyContent.append(prmContent1[prmID1]);
                    if (!isSameParticle) {
                        keyContent.append(LINESEPARATOR);
                        keyContent.append(prmContent2[prmID2]);
                    }
                }
                MIPETUTIL.writeKeyFile(keyPath, keyContent.toString());
                double distCandidate;

                for (int i = 0; i < distSize; i++) {
                    distCandidate = lowerBoundary + i * prescanStepSize;
                    distList.add(distCandidate);
                    if (distCandidate > upperBoundary) {
                        distList.add(distCandidate);
                        break;
                    }
                }

                List<Double> allDistances = new ArrayList<>(distList);
                double[] distances = distList.stream()
                        .mapToDouble(d -> d)
                        .toArray();
                
                //</editor-fold>
                EnergyRecord[] energyRecords;
                double glbEmin;
                if (particlePair.equals("Na+_Na+") || 
                        particlePair.equals("Ac-_Ac-")) {
                    
                    // Check whether particle pair is Na+_Na+ or Ac-_Ac-
                    energyList.add(new ResultEnergyRecord(
                            particle1, 
                            particle2, 
                            0,
                            0,
                            0,
                            0));
                } else {
                    energyRecords = new EnergyRecord[5];
                    energyRecords[0] = getInterMolecularEnergy(
                            particlePair,
                            distances, 
                            tinkerXYZ1, 
                            tinkerXYZ2, 
                            xyzRotData1, 
                            xyzRotData2,
                            1E10);
                    double glbEminDist = energyRecords[0].eqDistance();
                    glbEmin = energyRecords[0].eMin();

                    for (int i = 0; i < distances.length; i++) {
                        distEminRecords.add(new Distance_EnergyRecord(
                                energyRecords[0].distances()[i], 
                                energyRecords[0].eMins()[i],
                                energyRecords[0].wgtEmins()[i]));
                    }

                    //<editor-fold defaultstate="collapsed" desc="Calculate rotated coordinates">
                    if (sphereNodeNumber2 != sphereNodeNumber1 
                            || rotationNumber2 != rotationNumber1) {
                        rotRecord = RotationUtil.getRotationsCoords(
                                sphereNodeNumber2, 
                                rotationNumber2,
                                data1,
                                xyzData2,
                                isFibonacciSphereAlgorithm);
                        xyzRotData1 = rotRecord.rotData1();
                        xyzRotData2 = rotRecord.rotData2();
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Precise scan">
                    distList = new LinkedList<>();
                    distSize = 9;

                    for (int i = 0; i < distSize; i++) {
                        distCandidate = (10 * glbEminDist - 4 + i) 
                                / 10;
                        distList.add(distCandidate);
                    }

                    allDistances.addAll(distList);
                    distances = distList.stream()
                            .mapToDouble(d -> d)
                            .toArray();
                    energyRecords[1] = getInterMolecularEnergy(
                            particlePair,
                            distances, 
                            tinkerXYZ1, 
                            tinkerXYZ2, 
                            xyzRotData1, 
                            xyzRotData2,
                            glbEmin);
                    if (energyRecords[1].eMin() < glbEmin) {
                        glbEminDist = energyRecords[1].eqDistance();
                        glbEmin = energyRecords[1].eMin();
                    }

                    for (int i = 0; i < distances.length; i++) {
                        distEminRecords.add(new Distance_EnergyRecord(
                                energyRecords[1].distances()[i], 
                                energyRecords[1].eMins()[i],
                                energyRecords[1].wgtEmins()[i]));
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Calculate rotated coordinates">
                    if (sphereNodeNumber3 != sphereNodeNumber2 
                            || rotationNumber3 != rotationNumber2) {
                        rotRecord = RotationUtil.getRotationsCoords(
                                sphereNodeNumber3, 
                                rotationNumber3,
                                data1,
                                xyzData2,
                                isFibonacciSphereAlgorithm);
                        xyzRotData1 = rotRecord.rotData1();
                        xyzRotData2 = rotRecord.rotData2();
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="More precise scan">
                    distList = new LinkedList<>();
                    distSize = 19;

                    for (int i = 0; i < distSize; i++) {
                        distCandidate = (100 * glbEminDist - 9 + i) / 100;
                        distList.add(distCandidate);
                    }

                    allDistances.addAll(distList);
                    distances = distList.stream()
                            .mapToDouble(d -> d)
                            .toArray();
                    energyRecords[2] = getInterMolecularEnergy(
                            particlePair,
                            distances, 
                            tinkerXYZ1, 
                            tinkerXYZ2, 
                            xyzRotData1, 
                            xyzRotData2,
                            glbEmin);
                    if (energyRecords[2].eMin() < glbEmin) {
                        glbEminDist = energyRecords[2].eqDistance();
                        glbEmin = energyRecords[2].eMin();
                    }
                    
                    for (int i = 0; i < distances.length; i++) {
                        distEminRecords.add(new Distance_EnergyRecord(
                                energyRecords[2].distances()[i], 
                                energyRecords[2].eMins()[i],
                                energyRecords[2].wgtEmins()[i]));
                    }
                    
                    allDistances.add(energyRecords[2].eqDistance());

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Calculate rotated coordinates">
                    if (sphereNodeNumber4 != sphereNodeNumber3 
                            || rotationNumber4 != rotationNumber3) {
                        rotRecord = RotationUtil.getRotationsCoords(
                                sphereNodeNumber4, 
                                rotationNumber4,
                                data1,
                                xyzData2,
                                isFibonacciSphereAlgorithm);
                        xyzRotData1 = rotRecord.rotData1();
                        xyzRotData2 = rotRecord.rotData2();
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Scan at equilibrium distance">
                    distances = new double[1];
                    distances[0] = glbEminDist; 
                    energyRecords[3] = getInterMolecularEnergy(
                            particlePair,
                            distances, 
                            tinkerXYZ1, 
                            tinkerXYZ2, 
                            xyzRotData1, 
                            xyzRotData2,
                            glbEmin);
                    if (energyRecords[3].eMin() < glbEmin) {
                        glbEminDist = energyRecords[3].eqDistance();
                        glbEmin = energyRecords[3].eMin();
                    }

                    for (int i = 0; i < distances.length; i++) {
                        distEminRecords.add(new Distance_EnergyRecord(
                                energyRecords[3].distances()[i], 
                                energyRecords[3].eMins()[i],
                                energyRecords[3].wgtEmins()[i]));
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Find distance by min. wgtEmin">
                    int distEminRecordsSize = distEminRecords.size();
                    double wgtEminDist = distEminRecords.getFirst().distance();
                    double glbWgtEmin = distEminRecords.getFirst().wgtEmin();
                    double glbWgtEminCand;
                    
                    for (int i = 1; i < distEminRecordsSize; i++) {
                        glbWgtEminCand = distEminRecords.get(i).wgtEmin();
                        if (glbWgtEminCand < glbWgtEmin) {
                            glbWgtEmin = glbWgtEminCand;
                            wgtEminDist = distEminRecords.get(i).distance();
                        }
                    }
                    
                    // If wgtEminDist is too far away from glbEminDist new calculations around wgtEminDist will occur
                    // This is only for the output particle1_particle2_DistanceVsEaverage.jpg
                    if (Math.abs(wgtEminDist - glbEminDist) > 0.1) {
                        distList = new LinkedList<>();
                        distSize = 19;

                        for (int i = 0; i < distSize; i++) {
                            distCandidate = (100 * wgtEminDist - 9 
                                    + i) / 100;
                            if (!MIPETUTIL.contains(allDistances, 
                                    distCandidate)) {
                                distList.add(distCandidate);
                            }
                        }

                        allDistances.addAll(distList);
                        distances = distList.stream()
                                .mapToDouble(d -> d).toArray();
                        energyRecords[4] = getInterMolecularEnergy(
                                particlePair,
                                distances, 
                                tinkerXYZ1, 
                                tinkerXYZ2, 
                                xyzRotData1, 
                                xyzRotData2,
                                glbEmin);
                        
                        for (int i = 0; i < distances.length; i++) {
                            distEminRecords.add(new Distance_EnergyRecord(
                                    energyRecords[4].distances()[i], 
                                    energyRecords[4].eMins()[i],
                                    energyRecords[4].wgtEmins()[i]));
                            glbWgtEminCand = energyRecords[4]
                                    .wgtEmins()[i];
                            if (glbWgtEmin > glbWgtEminCand) {
                                glbWgtEmin = glbWgtEminCand;
                            }
                        }
                        
                    }

                    //<editor-fold defaultstate="collapsed" desc="Sort datas">
                    distSize = allDistances.size();
                    Double[] distObj = allDistances
                            .toArray(Double[]::new);
                    ArrayIndexComparator comparator = 
                            new ArrayIndexComparator(distObj);
                    Integer[] distIndices = comparator.createIndexArray();
                    Arrays.sort(distIndices, comparator);
                    double[][] sortedResults = new double[distSize][3];

                    for (int i = 0; i < distSize; i++) {
                        sortedResults[i][0] = distEminRecords
                                .get(distIndices[i])
                                .distance();
                        sortedResults[i][1] = distEminRecords
                                .get(distIndices[i])
                                .eMin();
                        sortedResults[i][2] = distEminRecords
                                .get(distIndices[i])
                                .wgtEmin();
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Determining opt. Emin">
                    double optMinEnergy = 0.0;
                    double rgdMinEnergy = 0.0;
                    
                    if (isTinkerOn) {
                        // Delete old .key file and make new one
                        Path keyOldPath = Paths.get(scratchDirectory, 
                                particlePair + ".key");
                        try {
                            Files.deleteIfExists(keyOldPath);
                        } catch (IOException ex) {
                            LOGGER.log(Level.SEVERE,
                                    "IOException during deleting .key file in scratch.",
                                    ex);
                        }

                        // Determining opt. Emin
                        Path keyNewPath = Paths.get(scratchDirectory, 
                                particlePair + ".key");
                        String optKeyContent = keyFileString;
                        if (forcefield.equals("OPLSAALIGPARGEN")) {
                            optKeyContent += prmContent1[prmID1];
                            if (!isSameParticle) {
                                optKeyContent += prmContent2[prmID2];
                            }
                        }
                        MIPETUTIL.writeKeyFile(keyNewPath, optKeyContent);
                        String[] cmdList;

                        for (int i = 0; i < 2; i++) {
                            if (i == 0) {
                                cmdList = new String[] {tinkerOptimize, 
                                    scratchDirectory 
                                    + FILESEPARATOR 
                                    + particlePair
                                    + ".0",        
                                    Double.toString(optimizeRmsGradient)};
                            } else {
                                cmdList = new String[] {tinkerOptrigid, 
                                    scratchDirectory 
                                    + FILESEPARATOR 
                                    + particlePair
                                    + ".0",
                                    Double.toString(optimizeRmsGradient)};
                            }

                            // Start optimize
                            builder = new ProcessBuilder();
                            builder.command(cmdList);
                            builder.redirectErrorStream(true);
                            builder.redirectOutput(ProcessBuilder
                                    .Redirect.DISCARD); 
                            Process process = null;
                            try {
                                process = builder.start();
                                if (process != null) {
                                    try {
                                        process.waitFor();
                                    } catch (InterruptedException ex) {
                                        Thread.currentThread().interrupt();
                                        LOGGER.log(Level.SEVERE,
                                                "InterruptException during processing optimize.exe", ex);
                                        process.destroy();
                                    }
                                }
                            } catch (IOException ex) {
                                LOGGER.log(Level.SEVERE, 
                                        "IOException during process starting.",
                                        ex);
                            }

                            // Fix .xyz file if there is H2O 
                            //  this is necessary because of a bug in tinker's optimize
                            if (forceField_IE.equals("OPLSAALIGPARGEN")
                                    && pos_H2O > 0) {
                                String fileName = scratchDirectory
                                    + FILESEPARATOR
                                    + particle1 
                                    + "_"
                                    + particle2
                                    + ".xyz";
                                int idx;
                                if (particle1.equals("H2O")) {
                                    idx = 1;
                                } else {
                                    idx = tinkerXYZ1.getN_atom1() + 1;
                                }
                                MIPETUTIL.fixTinkerXYZ_H2O(fileName, idx);
                            }

                            // Use tinker's analyze to determine intermolecular energy
                            String outputName;
                            if (i == 0)  {
                                outputName = scratchDirectory
                                    + FILESEPARATOR
                                    + "output0_opt.txt";
                            } else {
                                outputName = scratchDirectory
                                    + FILESEPARATOR
                                    + "output0_rgd.txt";
                            }
                            File optFile = new File(outputName);
                            cmdList = new String[] {tinkerAnalyze, 
                                scratchDirectory
                                + FILESEPARATOR 
                                + particlePair
                                + ".xyz",
                                "E"};
                            builder = new ProcessBuilder();
                            builder.redirectErrorStream(true);
                            builder.command(cmdList);
                            builder.redirectOutput(optFile);
                            try {
                                process = builder.start();
                                process.waitFor();
                            } catch (IOException ex) {
                                LOGGER.log(Level.SEVERE, 
                                        "IOException during process start.",
                                        ex);
                            } catch (InterruptedException ex) {
                                LOGGER.log(Level.SEVERE, 
                                        "InterruptException during reading output0_opt.txt",
                                        ex);
                            }
                            if (process != null) {
                                process.destroy();
                            }
                            Path sourcePath = Paths.get(scratchDirectory,
                                    particlePair + ".xyz");
                            Path targetPath;
                            if (i == 0) {
                                targetPath = Paths.get(scratchDirectory,
                                        particlePair + "_opt" + ".xyz");
                            } else {
                                targetPath = Paths.get(scratchDirectory,
                                        particlePair + "_rgd" + ".xyz");
                            }
                            try {
                                Files.move(sourcePath, targetPath, 
                                        StandardCopyOption.REPLACE_EXISTING);
                            } catch (IOException ex) {
                                LOGGER.log(Level.SEVERE, 
                                        "IOException renaming .xyz file.", ex);
                            }

                           // Read the intermolecular energies from .txt files
                            String searchStr = "Intermolecular Energy :";
                            Path path = Paths.get(outputName);

                            try (BufferedReader reader = Files
                                    .newBufferedReader(path, 
                                            StandardCharsets.UTF_8)) {
                                String line;

                                while ((line = reader.readLine()) != null ) {
                                    if (line.contains(searchStr)) {
                                        try {
                                            String energyValueStr = line.substring(25).trim();
                                            double extractedEnergy = Double.parseDouble(energyValueStr);
                                            if (i == 0) {
                                                optMinEnergy = extractedEnergy;
                                            } else {
                                                rgdMinEnergy = extractedEnergy;
                                            }
                                        } catch (NumberFormatException | StringIndexOutOfBoundsException parseEx) {
                                            LOGGER.log(Level.WARNING, "Error during parsing of energy value in file: " + outputName, parseEx);
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
                        particle1, 
                        particle2, 
                        glbWgtEmin,
                        optMinEnergy,
                        rgdMinEnergy,
                        glbEmin));
                    
                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Write dist vs. energy datas">
                    Path outputPath = Paths.get(
                            jobTaskRecordList.get(currIndex)
                                    .result_IE_PathName(),
                                    particlePair + "_dist_vs_energy.dat");
                    try (PrintWriter writer = new PrintWriter(
                            Files.newBufferedWriter(outputPath))) {
                        writer.append("distance [" + ANGSTROM +
                                "]  Emin(Cmin,r) [kcal/mole] <E>(r) [kcal/mole]");
                        writer.println();

                        for (int i = 0; i < distSize; i++) {
                            writer.printf("%8s%20s%20s%n",
                                    decimal2.format(sortedResults[i][0]),
                                    decimal3.format(sortedResults[i][1]),
                                    decimal3.format(sortedResults[i][2]));
                        }

                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during writing output0.out in scratch."
                                , ex);
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Make dist vs. energy diagram">
                    String energyDataPathName = jobTaskRecordList
                            .get(currIndex).result_IE_PathName()
                            + FILESEPARATOR
                            + particlePair + "_dist_vs_energy.dat";
                    String energyGraphicsPrefix = jobTaskRecordList
                            .get(currIndex).particleName1()
                            + "_"
                            + jobTaskRecordList
                                    .get(currIndex).particleName2()
                            + "_";
                    ChartUtil chartUtil = new ChartUtil();
                    chartUtil.createEnergyGraphics(
                            energyDataPathName,
                            energyGraphicsPrefix);

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Write log file">
                    try {
                        BFGblLog.append(particlePair);
                        BFGblLog.append(LINESEPARATOR);
                        BWParticleDat.append("Force field: ")
                                .append(forcefield);
                        BWParticleDat.append(LINESEPARATOR);
                        BWParticleDat.append("Conformational analysis: "); 
                        if (isConformationalAnalysis) {
                            BWParticleDat.append("Yes");
                        } else {
                            BWParticleDat.append("No");
                        }
                        BWParticleDat.append(LINESEPARATOR);
                        BWParticleDat.append("equilibriumDistances [" + ANGSTROM + "] = "); 
                        BWParticleDat.append(decimal2.format(glbEminDist));
                        BWParticleDat.append(LINESEPARATOR);
                        BWParticleDat.flush();
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during writing log file.", ex);
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Copy results">
                    //<editor-fold defaultstate="collapsed" desc="Write ouput.0 file">
                    Path original = Paths.get(scratchDirectory, 
                            particlePair + ".0");
                    Path target = Paths.get(ieResultDirName, "output.0");
                    try {
                        Files.copy(original, target, 
                                StandardCopyOption.REPLACE_EXISTING);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during copying output.0", ex);
                    }
                    String output0 = "Intermolecular Energy: " 
                            + energyRecords[3].wgtEmin() + " kcal/mol";
                    Path targetPath = Paths.get(ieResultDirName, 
                            "output0.out");        
                    try (BufferedWriter writer = Files.newBufferedWriter(
                            targetPath, StandardCharsets.UTF_8)) {
                        writer.append(output0);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during writing output0.out in scratch."
                                , ex);
                    }
                    
                    //</editor-fold>
                    if (isTinkerOn) {
                        //<editor-fold defaultstate="collapsed" desc="Write output_opt.out">
                        original = Paths.get(scratchDirectory, 
                                particlePair + "_opt.xyz");
                        target = Paths.get(ieResultDirName,"output_opt.0");
                        try {
                            Files.copy(original, target,
                                StandardCopyOption.REPLACE_EXISTING);
                        } catch (IOException ex) {
                            LOGGER.log(Level.SEVERE, 
                                    "IOException during copying output.0", ex);
                        }
                        output0 = "Intermolecular Energy: " 
                                + optMinEnergy + " kcal/mol";
                        targetPath = Paths.get(ieResultDirName, 
                                "output0_opt.out");        
                        try (BufferedWriter writer = Files.newBufferedWriter(
                                targetPath, StandardCharsets.UTF_8)) {
                            writer.append(output0);
                        } catch (IOException ex) {
                            LOGGER.log(Level.SEVERE, 
                                    "IOException during writing output0_opt.out "
                                    + "in scratch.", ex);
                        }
                        
                        //</editor-fold>
                        //<editor-fold defaultstate="collapsed" desc="Write output_rgd.out">
                        original = Paths.get(scratchDirectory, 
                                particlePair + "_rgd.xyz");
                        target = Paths.get(ieResultDirName, "output_rgd.0");
                        try {
                            Files.copy(original, target,
                                StandardCopyOption.REPLACE_EXISTING);
                        } catch (IOException ex) {
                            LOGGER.log(Level.SEVERE, 
                                    "IOException during copying output.0", ex);
                        }
                        output0 = "Intermolecular Energy: " 
                                + rgdMinEnergy + " kcal/mol";
                        targetPath = Paths.get(ieResultDirName, 
                                "output0_rgd.out");        
                        try (BufferedWriter writer = Files.newBufferedWriter(
                                targetPath, StandardCharsets.UTF_8)) {
                            writer.append(output0);
                        } catch (IOException ex) {
                            LOGGER.log(Level.SEVERE, 
                                    "IOException during writing output0_rgd.out "
                                    + "in scratch.", ex);
                        }
                        
                        //</editor-fold>
                    }

                    //<editor-fold defaultstate="collapsed" desc="Generate output.xyz">
                    Path output0Path = Paths.get(ieResultDirName, "output.0");
                    int atomNumber1 = tinkerXYZ1.getN_atom1();
                    int atomNumber2 = tinkerXYZ2.getN_atom1();
                    String[] elementList1 = tinkerXYZ1.getElementList1();
                    String[] elementList2 = tinkerXYZ2.getElementList1();
                    TinkerXYZ tinkerXyz = new TinkerXYZ(output0Path, 1, 
                            atomNumber1, atomNumber2);
                    String fileName = ieResultDirName 
                            + FILESEPARATOR 
                            + "output.xyz";
                    tinkerXyz.writeToXyzFile(fileName);
                    
                    //</editor-fold>
                    
                    // Generate .pdb file of output.0
                    Path pdbPath = output0Path.resolveSibling("output.pdb");
                    TinkerToPdbConverter.convert(output0Path, pdbPath);
                    
                    if (isTinkerOn) {
                        // Generate output_opt.xyz
                        Path optPath = Paths.get(ieResultDirName,
                                "output_opt.0");
                        tinkerXyz = new TinkerXYZ(optPath, 1, 
                                atomNumber1, atomNumber2);
                        tinkerXyz.setElementList1(elementList1);
                        tinkerXyz.setElementList2(elementList2);
                        fileName = ieResultDirName 
                                + FILESEPARATOR 
                                + "output_opt.xyz";
                        tinkerXyz.writeToXyzFile(fileName);

                        // Generate output_rgd.xyz
                        Path rgdPath = Paths.get(ieResultDirName, 
                                "output_rgd.0");
                        tinkerXyz = new TinkerXYZ(rgdPath, 1, 
                                atomNumber1, atomNumber2);
                        tinkerXyz.setElementList1(elementList1);
                        tinkerXyz.setElementList2(elementList2);
                        fileName = ieResultDirName 
                                + FILESEPARATOR 
                                + "output_rgd.xyz";
                        tinkerXyz.writeToXyzFile(fileName);

                        //<editor-fold defaultstate="collapsed" desc="Generate .pdb files">
                        pdbPath = optPath.resolveSibling("output_opt.pdb");
                        TinkerToPdbConverter.convert(optPath, pdbPath);
                        pdbPath = rgdPath.resolveSibling("output_rgd.pdb");
                        TinkerToPdbConverter.convert(rgdPath, pdbPath);
                        
                        //</editor-fold>
                    }

                    // Further outputs
                    Path optDistDir = Paths.get(optDistDirectory 
                            + FILESEPARATOR 
                            + forcefield);
                    try {
                        Files.createDirectories(optDistDir);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                            "IOException during creating OutDist directory.",
                            ex);
                    }
                    Path optDistPath = Paths.get(optDistDirectory, forcefield,
                            particle1 + ".txt");
                    if (isSameParticle && !Files.exists(optDistPath)) {
                        try (BufferedWriter writer = Files.newBufferedWriter(
                                optDistPath, StandardCharsets.UTF_8)) {
                            writer.append(decimal4.format(glbEminDist));
                        } catch(IOException ex) {
                            LOGGER.log(Level.SEVERE, 
                                "IOException during writing file in OptDist directory.", 
                                ex);
                        }
                    }
                    energyCalcTime  = (System.currentTimeMillis() 
                            - energyCalcTime) / 1000;
                    try {
                        BWParticleDat.append("SphereNodeNumber1: ")
                                .append(Integer.toString(sphereNodeNumber1));
                        BWParticleDat.newLine();
                        BWParticleDat.append("SphereNodeNumber2: ")
                                .append(Integer.toString(sphereNodeNumber2));
                        BWParticleDat.newLine();
                        BWParticleDat.append("SphereNodeNumber3: ")
                                .append(Integer.toString(sphereNodeNumber3));
                        BWParticleDat.newLine();
                        BWParticleDat.append("SphereNodeNumber4: ")
                                .append(Integer.toString(sphereNodeNumber4));
                        BWParticleDat.newLine();
                        BWParticleDat.append("RotationNumber1: ")
                                .append(Integer.toString(rotationNumber1));
                        BWParticleDat.newLine();
                        BWParticleDat.append("RotationNumber2: ")
                                .append(Integer.toString(rotationNumber2));
                        BWParticleDat.newLine();
                        BWParticleDat.append("RotationNumber3: ")
                                .append(Integer.toString(rotationNumber3));
                        BWParticleDat.newLine();
                        BWParticleDat.append("RotationNumber4: ")
                                .append(Integer.toString(rotationNumber4));
                        BWParticleDat.newLine();
                        BWParticleDat.append("Temperature [K]: ");
                        BWParticleDat.append(Double.toString(temperature));
                        BWParticleDat.append(LINESEPARATOR);
                        BWParticleDat.append("Fraction of energy values used for the Boltzmann distribution: ");
                        BWParticleDat.append(Double.toString(boltzmannFraction));
                        BWParticleDat.newLine();
                        if (isOptEmin) {
                            BWParticleDat.append("Optimize sampled E(min) configuration: ")
                                    .append(String.valueOf(isOptEmin));
                            BWParticleDat.newLine();
                        }
                        BWParticleDat.append("Weighted (Emin = glbMin) MinimumIntermolecularEnergy [kcal/mole]: ");
                        BWParticleDat.append(decimal4.format(glbWgtEmin));
                        BWParticleDat.newLine();
                        BWParticleDat.append("""
                                             Weighted (Emin = glbMin) MinimumIntermolecularEnergy:
                                             Weighted differential pair interaction energy with
                                              Emin = lowest differential pair interaction energy.""");
                        BWParticleDat.newLine();
                        BWParticleDat.append("GlobalMinimumIntermolecularEnergy [kcal/mole]: ");
                        BWParticleDat.append(decimal4.format(glbEmin));
                        BWParticleDat.newLine();
                        BWParticleDat.append("""
                                             GlobalMinimumIntermolecularEnergy:
                                               Lowest differential pair interaction energy of all dimer configurations.""");
                        BWParticleDat.newLine();
                        if (isTinkerOn) {
                            BWParticleDat.append("Optimized minimumIntermolecularEnergy [kcal/mole]: ");
                            BWParticleDat.append(decimal4.format(optMinEnergy));
                            BWParticleDat.newLine();
                            BWParticleDat.append("""
                                             Optimized minimumIntermolecularEnergy:
                                               Differential pair interaction energy from the dimer configuration
                                               with lowest differential pair interaction energy after optimize.""");
                            BWParticleDat.newLine();
                            BWParticleDat.append("Rigid-optimized minimumIntermolecularEnergy [kcal/mole]: ");
                            BWParticleDat.append(decimal4.format(rgdMinEnergy));
                            BWParticleDat.newLine();
                            BWParticleDat.append("""
                                             Rigid-optimized minimumIntermolecularEnergy:
                                               Differential pair interaction energy from the dimer configuration
                                               with lowest differential pair interaction energy after optrigid.""");
                            BWParticleDat.newLine();
                        }
                        BWParticleDat.append("Time to calculate minimum intermolecular energy [s]: ")
                                .append(String.valueOf(energyCalcTime));
                        BWParticleDat.newLine();
                        BWParticleDat.close();
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during writing log file.", ex);
                    }

                    //</editor-fold>
                    //<editor-fold defaultstate="collapsed" desc="Clean scratch directory">
                    try {
                        FileUtils.deleteDirectory(new File(scratchDirectory));
                    } catch(IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during deleting scratch directory.", ex);
                    }

                    //</editor-fold>
                }
                
            }
            currIndex++;
            if (currIndex >= jobTaskRecordList.size()) {
                isExitCondition = true;
            }
        }
        
        //<editor-fold defaultstate="collapsed" desc="Write Readme.txt">
        String readmeStr =
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
        if (forceField_IE != null) {
            Path readmePath = Paths.get(resultDirectory, "IE", forceField_IE,
                    "Readme.txt");
            try (BufferedWriter writer = Files.newBufferedWriter(readmePath,
                    StandardCharsets.UTF_8)) {
                writer.append(readmeStr);
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE,
                        "IOException during writing Readme.txt.", ex);
            }
        }
                
        //</editor-fold>
        //<editor-fold defaultstate="collapsed" desc="Coordination numbers">
        if (forceField_CN != null) {
            System.out.println("Calculating coordination numbers...");
            getCoordinationNumbers(jobTaskRecordList);
        }
        
        //</editor-fold>
        //<editor-fold defaultstate="collapsed" desc="Export parameterset">
        System.out.print("\rExporting parameterset...          ");
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
        if (forceField_IE != null) {
            exportParticleSetForMFSim(jobTaskRecordList,
                    parameterSetTitle,
                    parameterSetTitleAbr);
        }
        System.out.println("\r\nReady.");
        
        //</editor-fold>
    }
    
    //</editor-fold>
    //</editor-fold>
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Initialize method
     * 
     * @param args Input arguments
     */
    private static void initialize(String args[]) {
        // Read title/version
        String appTitle = MIPET.class.getPackage().getImplementationVersion();

        // Fallback for test in NetBeans
        if (appTitle == null) {
            appTitle = "MIPET (Version 1.0.6)";
        }

        // Ask for versionnumber 
        if (args.length > 0 && (args[0].equals("-v") || args[0]
                .equals("--version"))) {
            System.out.println(appTitle);
            return; // Quit the program
        }
        setParameters();
        
        // <editor-fold defaultstate="collapsed" desc="Check and create directories">
        try {
            Files.createDirectories(Paths.get(scratchDirectory));
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                        "IOException during creating scratch directory.", ex);
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
        } else if (!MIPETUTIL.isDirectoryEmpty(resultDirectory)) {
            System.out.println(resultDirectory + " is not empty.");
            System.exit(-1);
        } else if (!MIPETUTIL.isDirectoryEmpty(optXYZDirectory)) {
            System.out.println(optXYZDirectory + " is not empty.");
            System.exit(-1);
        }else  if (!MIPETUTIL.isDirectoryEmpty(optDistDirectory)) {
            System.out.println(optDistDirectory + " is not empty.");
            System.exit(-1);
        }
        
        // Path finding for SMILES (instead using System.getProperty) ---
        String smilesPath = "Molecules" + FILESEPARATOR 
                + "SMILES" + FILESEPARATOR + "Smiles.dat";
        try {
            // Get path of mipet.jar
            File jarFile = new java.io.File(MIPET.class.getProtectionDomain()
                    .getCodeSource().getLocation().toURI());
            if (jarFile.isFile()) {
                // Path "app" in installation directory
                String appDir = jarFile.getParent();
                File testPath = new java.io.File(appDir, smilesPath);
                if (testPath.exists()) {
                    smilesPath = testPath.getAbsolutePath();
                }
            }
        } catch (Exception e) {
            // Ignorie and Fallback using relative path (for development)
        }
        smiles = MIPETUTIL.getSmilesData(smilesPath);
        
        // </editor-fold>
    }
    
    /**
     * Read the Job file
     */
    private static void readJobFile() {
        particleNames = new ArrayList<>();
        newParticles = new ArrayList<>();
        oldParticles = new ArrayList<>();
        Path jobPath = Paths.get(jobFileName);
        try (BufferedReader reader = Files.newBufferedReader(jobPath, 
                StandardCharsets.UTF_8)) {
            String line;
            
            // read jobs
            while ((line = reader.readLine()) != null ) {
                String trimmedLine = line.trim();
                if (line.isEmpty()) {
                    continue;
                }
                char firstChar = trimmedLine.charAt(0);
                String restString = trimmedLine.substring(1).trim();
                switch (firstChar) {
                    case '#' -> {
                     // ignore comment line
                    }
                    case '*' -> {
                        String[] str = restString.split("\\s+");
                        forceField_IE = str[0];
                    }
                    case '$' -> {
                        String[] str = restString.split("\\s+");
                        forceField_CN = str[0];
                    }
                    case '-' -> {
                        if (!oldParticles.contains(restString)) {
                            oldParticles.add(restString);
                        }
                        if (!particleNames.contains(restString)) {
                            particleNames.add(restString);
                        }
                    }
                    default -> {
                        line = line.trim();
                        if (!newParticles.contains(line)) {
                            newParticles.add(line);
                        }
                        if (!particleNames.contains(line)) {
                            particleNames.add(line);
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
        if (forceField_IE == null && forceField_CN == null) {
            System.out.println("No forcefield was defined in .job file.");
            throw new IllegalStateException("No forcefield was defined in .job file.");
        }
        if ("OPLSAALIGPARGEN".equals(forceField_IE)) {
            
            for (String particle : particleNames) {
                if (!smiles.containsKey(particle)) {
                    System.out.println(particle + " was not found in Smiles.dat");
                    throw new IllegalStateException(particle 
                            + " was not found in Smiles.dat");
                }

                // Check whether .xyz file exists
                Path xyzPath = Path.of(moleculeDirectory, 
                         "OPLSAALIGPARGEN",
                         particle + ".xyz");
                if (!Files.exists(xyzPath)) {
                    System.out.println(particle + ".xyz was not found.");
                    throw new IllegalStateException(particle 
                            + ".xyz was not found.");
                }
                
                // Check whether .prm file exists
                Path prmPath = Path.of(parameterDirectory, 
                         "OPLSAALIGPARGEN",
                         particle + ".prm");
                if (!Files.exists(prmPath)) {
                    System.out.println(particle + ".prm was not found.");
                    throw new IllegalStateException(particle 
                            + ".prm was not found.");
                }
            }
            
        }
    }
    
    /**
     * SetParameters method
     */
    private static void setParameters() {
        // Determine the main path of the installation
        String appBaseDir = "";
        try {
            File jarFile = new java.io.File(MIPET.class.getProtectionDomain()
                    .getCodeSource().getLocation().toURI());
            if (jarFile.isFile()) {
                appBaseDir = jarFile.getParent() + File.separator; // e.g. "C:\programs\MIPET 1.0.6\app\"
            }
        } catch (Exception e) {
            // Ignore
        }
        isTinker9 = MIPETUTIL.getResourceString("MIPET.Tinker9")
                .equalsIgnoreCase("true");
        isTinkerOn = MIPETUTIL.getResourceString("MIPET.TinkerOn")
                .equalsIgnoreCase("true");
        jobFileName = MIPETUTIL.getResourceString("MIPET.File.jobfile");
        // Smart fallback for the job file
        // 1. Try to find the file in current working directory
        // 2. If not found then use the app-directory of the installation
        if (!Files.exists(Paths.get(jobFileName)) && !appBaseDir.isEmpty()) {
            Path fallbackJobPath = Paths.get(appBaseDir, jobFileName);
            if (Files.exists(fallbackJobPath)) {
                jobFileName = fallbackJobPath.toString();
            }
        }
        cpuCoreNumber = Integer.parseInt(MIPETUTIL.getResourceString(
                "MIPETCPUCoreNumber"));
        isFibonacciSphereAlgorithm = MIPETUTIL.getResourceString(
                "MIPETJavaUseFibonacciSphereAlgorithm").equals("true");
        isOptEmin = MIPETUTIL.getResourceString("MIPETOptEmin").equals("true");
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
        parameterDirectory = appBaseDir + MIPETUTIL.getResourceString(
                "MIPET.Directory.paramDirectory");
        /*
         * Directory of tinker files
         */
        String tinkerDirectory = MIPETUTIL
                .getResourceString("MIPET.Directory.tinker");
        if (OS.matches("Windows [0-9]*")) {
            tinkerAnalyze = tinkerDirectory + "/analyze.exe";
            tinkerScan = tinkerDirectory + "/scan.exe";
            tinkerXYZEdit = tinkerDirectory + "/xyzedit.exe";
            tinkerMinimize = tinkerDirectory + "/minimize.exe";
            tinkerDynamic = tinkerDirectory + "/dynamic.exe";
            tinkerOptimize = tinkerDirectory + "/optimize.exe";
            tinkerOptrigid = tinkerDirectory + "/optrigid.exe";
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
        }
        scratchDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.scratch");
        smilesDirectory = appBaseDir + MIPETUTIL.getResourceString(
                "MIPET.Directory.smiles");
        /*
          Directory of calculation related datas
         */
        String calcDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.calculation");
        resultDirectory = MIPETUTIL.getResourceString(
                "MIPET.Directory.result");
        moleculeDirectory = appBaseDir + MIPETUTIL.getResourceString(
                "MIPET.Directory.source");
        optXYZDirectory = calcDirectory + "/OptXYZ";
        optDistDirectory = calcDirectory + "/OptDist";
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
        int dynamicStepsPerCore = stepNumber / cpuCoreNumber;
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
        int dynamicWarmupIteration = (int) (warmUpStepNumber * warmUpTimeStep
                * 1E-15 / (warmUpPrintInterval * 1E-12));
        // timeStep in fs and printInterval in ps
        nDynamicIteration = (int)(stepNumber * timeStep * 1E-15 /
                (printInterval * 1E-12));
        // timeStep in fs and printInterval in ps
        nDynamicIterationPerCore = 
                (int) Math.ceil((dynamicStepsPerCore * timeStep * 1E-15) /
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
        int newParticleSize = newParticles.size();
        int oldParticleSize = oldParticles.size();
        int initialCapacity = newParticleSize + 
                (newParticleSize * oldParticleSize) +
                (newParticleSize * (newParticleSize - 1) / 2);
        ArrayList<String> particlePairs = new ArrayList<>(initialCapacity);
        Collections.sort(newParticles);
        Collections.sort(oldParticles);
       
        // Same particles
        for (String p: newParticles) {
            particlePairs.add(p + "_" + p);
        }
        
        // New with old particles
        for (String nP : newParticles) {

            for (String oP : oldParticles) {
                particlePairs.add(nP + "_" + oP);
            }
            
        }

        // Different particles
        for (int i = 0; i < newParticleSize; i++) {

            for (int j = 0; j < i; j++) {
                particlePairs.add(
                        newParticles.get(i) + "_" + newParticles.get(j));
            }

        }    
        
        return particlePairs;
    }
    
    /**
     * Method readXyz
     * Read .xyz files and store in xyzContent.
     * If conformational analysis is set xyzContent is overwritten with optimized data.
     * 
     * @param anIsOriginal Flag for whether original xyz-file from moleculeDirectory will be used (true)
     *                    or optimized (wih same chosen forcefield) xyz- file (false)
     * @param aForcefield Forcefield name
     */
    private static void readXyz(boolean anIsOriginal, String aForcefield) {
        int size = particleNames.size();
        xyzContent1 = new String[size];
        xyzContent2 = new String[size];
        
        for (int i = 0; i < size; i++) {
            String particleName = particleNames.get(i);
            Path filePath = buildPath(anIsOriginal, aForcefield, particleName);
            
            try{
                String content = Files.readString(filePath);
                xyzContent1[i] = content;
                if ("OPLSAALIGPARGEN".equals(aForcefield) && content != null &&
                        !content.isEmpty()) {
                    xyzContent2[i] = processXyzContent(content);
                } else {
                    xyzContent2[i] = content;
                }
            } catch(IOException ex) {
                LOGGER.log(Level.SEVERE, "Error during reading of the file: " + 
                        filePath, ex);
            }
        }
        
    }
    
    /**
     * Helping method for readXyz() to change atomnumber if second particle is different
     * 
     * @param content xyz-file content
     * @return Changed xyz-file content
     */
    private static String processXyzContent(String content) {
        StringBuilder sb = new StringBuilder();
        String[] lines = content.lines().toArray(String[]::new);

        if (lines.length > 0) {
            sb.append(lines[0]); // Header
        }

        for (int j = 1; j < lines.length; j++) {
            String[] tokens = lines[j].trim().split("\\s+");
            if (tokens.length >= 6) {
                try {
                    // Change of the sixth token (index 5)
                    int newType = Integer.parseInt(tokens[5]) + 100;
                    tokens[5] = String.valueOf(newType);

                    // Build the line together again 
                    sb.append(LINESEPARATOR)
                      .append(String.format("%-3s %12s %12s %12s %5s %5s", 
                              (Object[]) tokens)); 
                } catch (NumberFormatException e) {
                    sb.append(LINESEPARATOR)
                            .append(lines[j]);
                }
            } else {
                sb.append(System.lineSeparator()).append(lines[j]);
            }
        }
        
        return sb.toString();
    }

    /**
     * Build path from strings
     * 
     * @param isOriginal Flag for whether original xyz-file from moleculeDirectory will be used (true)
     *                    or optimized (wih same chosen forcefield) xyz- file (false)
     * @param forcefield Chosen forcefield
     * @param particleName Particle name
     * @return Path Path name
     */
    private static Path buildPath(boolean isOriginal, String forcefield, 
            String particleName) {
        if (isOriginal) {
            return Paths.get(moleculeDirectory, forcefield, 
                    particleName + ".xyz");
        } else {
            return Paths.get(optXYZDirectory, forcefield, particleName, 
                    particleName + ".xyz");
        }
    }
    
    /** 
     * Method readPrm
     * Read .prm files for OPLSAALIGPARGEN and convert atomtype number of 
     *  second particle
     * @param aForcefield Force field name
     */
    private static void readPrm(String aForcefield) {
        final String SEARCH = "Partial charge correction factor:";
        int particlesSize = particleNames.size();
        
        prmContent1 = new String[particlesSize];
        prmContent2 = new String[particlesSize];
        chargeCorr = new double[particlesSize];
                
        for (int i = 0; i < particlesSize; i++) {
            String particleName = particleNames.get(i);
            
            // Initialize
            prmContent1[i] = "";
            prmContent2[i] = "";
            chargeCorr[i] = 0.0;
            
            if (!"OPLSAALIGPARGEN".equals(aForcefield)) {
                continue;
            }
            
            // Read .prm file
            Path path = Paths.get(parameterDirectory, 
                    aForcefield, particleName + ".prm");
            try {
                if (Files.exists(path)) {
                    prmContent1[i] = Files.readString(path);
                }
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE, 
                        "IOException during reading .prm file.", ex);
            }
            
            // Read partial charge correction factor
            if (prmContent1[i].contains(SEARCH)) {
                int chargeCorrStartPos = prmContent1[i].indexOf(SEARCH) + 33;
                int chargeCorrEndPos = prmContent1[i].indexOf("\n", 
                        chargeCorrStartPos);
                
                String valueStr = prmContent1[i].substring(chargeCorrStartPos, 
                        chargeCorrEndPos).trim();
                chargeCorr[i] = Double.parseDouble(valueStr);
            } else {
                chargeCorr[i] = 1.0;
            }
            
            // Atom type of second particle
            if (!prmContent1[i].isEmpty()) {
                prmContent2[i] = prmContent1[i].lines()
                    .dropWhile(line -> !line.contains("Atom Type Definitions"))
                    .map(line -> {
                        String[] tokens1 = line.trim().split("\\s+");
                        if (tokens1.length > 0) {
                            return switch (tokens1[0]) {
                                case "atom", "vdw", "charge", "bond", "angle", "torsion", "imptors" 
                                     -> MIPETUTIL.changeAtomType(tokens1);
                                default -> line;
                            };
                        }
                        return line;
                    })
                    .filter(line -> !line.isEmpty())
                    .map(line -> line + LINESEPARATOR) 
                    .collect(Collectors.joining());
            }
        }
        
    }

    /** 
     * Method makeMoleculeRecord
     * Make molecules record for intermolecular energy calculation without tinker
     *  (only for OPLSAALIGPARGEN - with assignment check in prm-file)
     */
    private static void makeMoleculeRecord() {
        molecules = new LinkedList<>();
        int particleSize = particleNames.size();
        
        for (int i = 0; i < particleSize; i++) {
            // Read xyz-file (coordinates and atom types)
            String[] xyzLines = xyzContent1[i].lines().toArray(String[]::new);
            String[] xyzHeader = xyzLines[0].trim().split("\\s+");
            int numAtoms = Integer.parseInt(xyzHeader[0]);
            String[] elements = new String[numAtoms];
            int[] atomTypes = new int[numAtoms];
            
            // Begin with after header -> j = 1
            for (int j = 1; j < xyzLines.length; j++) {
                String[] tokens = xyzLines[j].trim().split("\\s+");
                elements[j - 1] = tokens[1]; // Element (2. Column)
                atomTypes[j - 1] = Integer.parseInt(tokens[5]); // Atom type (6. Column)
            }
            
            if (!"OPLSAALIGPARGEN".equals(forceField_IE)) {
                molecules.add(new MoleculeRecord(
                        particleNames.get(i),
                        numAtoms,
                        elements,
                        atomTypes,
                        null,
                        null,
                        null
                ));
                return;
            }
            
            // Read prm-file and copy in maps
            Map<Integer, Double> sigmaMap = new HashMap<>();
            Map<Integer, Double> epsilonMap = new HashMap<>();
            Map<Integer, Double> chargeMap = new HashMap<>();
            String[] prmLines = prmContent1[i].lines().toArray(String[]::new);
            
            for (String line : prmLines) {
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 3) {
                    continue; // Ignore empty line
                }
                if (tokens[0].equals("vdw")) {
                    int type = Integer.parseInt(tokens[1]); // Atom type nummer
                    sigmaMap.put(type, Double.valueOf(tokens[2]));
                    epsilonMap.put(type, Double.valueOf(tokens[3]));
                } else if (tokens[0].equals("charge")) {
                    int type = Integer.parseInt(tokens[1]); // Atom type nummer
                    chargeMap.put(type, Double.parseDouble(tokens[2]) 
                            * chargeCorr[i]);
                }
            }
            
            // Assign parameters to the atoms based on their type number
            double[] sigmas = new double[numAtoms];
            double[] epsilons = new double[numAtoms];
            double[] charges = new double[numAtoms];

            for (int j = 0; j < numAtoms; j++) {
                int currentAtomType = atomTypes[j];
                sigmas[j] = sigmaMap.getOrDefault(currentAtomType, 0.0);
                epsilons[j] = epsilonMap.getOrDefault(currentAtomType, 0.0);
                charges[j] = chargeMap.getOrDefault(currentAtomType, 0.0);
            }
            
            molecules.add(new MoleculeRecord(
                    particleNames.get(i),
                    numAtoms,
                    elements,
                    atomTypes,
                    epsilons,
                    sigmas,
                    charges
            ));
        }
        
    }
    
    /**
     * scanParticle method
     *  Optimize, scan and optimize the input particle
     */
    private static void scanParticle() {
        LinkedList<String> forcefieldList = new LinkedList<>();
        
        if (!forceField_IE.isEmpty()) {
            forcefieldList.add(forceField_IE);
        }
        if (!forceField_CN.isEmpty()) {
            forcefieldList.add(forceField_CN);
        }
        
        for (String forcefield : forcefieldList) {
            
            for (String particle : particleNames) {
                String particleXyzName = particle + ".xyz";
                Path optXyzDir = Paths.get(optXYZDirectory, forcefield, 
                        particle);
                try {
                    Files.createDirectories(optXyzDir);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during creating OptXyz directory.", 
                            ex);
                    continue;
                }
                
                // Check whether already scanned
                Path targetXyzPath = Paths.get(optXyzDir.toString(), 
                        particleXyzName);
                if (!Files.exists(targetXyzPath)) {
                    Path source = Paths.get(moleculeDirectory, forcefield,
                            particle + ".xyz");
                    try {
                        Files.copy(source, targetXyzPath);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during copying .xyz file.", ex);
                        continue;
                    }
                    TinkerXYZ tinkerXYZ = new TinkerXYZ(forcefield, particle,
                            targetXyzPath);
                    Path keyPath = Paths.get(optXyzDir.toString(), 
                            particle + ".key");
                    String forcefieldName = "OPLSAALIGPARGEN".equals(forcefield) 
                            ? "oplsaa" : forcefield;
                    String keyContent = "# Force Field Selection"
                            + LINESEPARATOR 
                            + "PARAMETERS\t\""
                            + parameterDirectory
                            + "/"
                            + forcefieldName.toLowerCase() + ".prm\""
                            + LINESEPARATOR
                            + "DIELECTRIC\t" 
                            + dielectricConstant
                            + LINESEPARATOR;
                    if (forcefield.equals("OPLSAALIGPARGEN")) {
                        int keyIndex = particleNames.indexOf(particle);
                        if (prmContent1[keyIndex] == null || 
                                prmContent1[keyIndex].isEmpty()) {
                            readPrm(forcefield);
                        }
                        keyContent += prmContent1[keyIndex];
                    }
                    
                    // Write .key file
                    MIPETUTIL.writeKeyFile(keyPath, keyContent);
                    
                    // Run Tinker Optimize (PRE)
                    Path optXyzLogName = optXyzDir.resolve(particle 
                            + "_preoptimize.log");
                    ProcessBuilder pbOpt = new ProcessBuilder(tinkerOptimize, 
                            targetXyzPath.toString(), 
                            Double.toString(optimizeRmsGradient));
                    pbOpt.redirectErrorStream(true);
                    pbOpt.redirectOutput(optXyzLogName.toFile());

                    try {
                        Process process = pbOpt.start();
                        process.waitFor();
                    } catch (IOException | InterruptedException ex) {
                        LOGGER.log(Level.SEVERE, "Exception during tinker's optimize.exe (pre)", ex);
                    }

                    TinkerXYZ tinkerXYZ0 = new TinkerXYZ(forcefield, particle, 
                            optXyzDir.resolve(particle + ".xyz_2"));
                    tinkerXYZ.setForcefieldName(forceField_IE);
                    tinkerXYZ.setCoordinateList1(tinkerXYZ0
                            .getCoordinateList1(), isTinkerOn);

                    try {
                        Files.deleteIfExists(optXyzDir.resolve(particle 
                                + ".xyz_2"));
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, "IOException during deleting .xyz_2 file", ex);
                    }

                    try (PrintWriter writer = new PrintWriter(targetXyzPath
                            .toFile())) {
                        writer.print(tinkerXYZ.getFileContent());
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, "IOException during writing .xyz file.", ex);
                    }
                    
                    // Run Tinker Scan
                    Path scanOutputName = optXyzDir.resolve(particle + ".out");
                    ProcessBuilder pbScan = new ProcessBuilder(tinkerScan, 
                            targetXyzPath.toString(), 
                            scanProgram, 
                            nSearchDirection, 
                            energyThreshold, 
                            rmsGradient);
                    pbScan.redirectErrorStream(true);
                    pbScan.redirectOutput(scanOutputName.toFile());

                    try {
                        Process process = pbScan.start();
                        process.waitFor();
                    } catch (IOException | InterruptedException ex) {
                        LOGGER.log(Level.SEVERE, "Exception during tinker's scan.exe", ex);
                    }

                    Path optArcFile = optXyzDir.resolve(particle + ".arc");
                    if (Files.exists(optArcFile) && 
                            !Files.isDirectory(optArcFile)) {
                        double[][][] coords = tinkerXYZ.readCoordFromArc(
                                optArcFile.toString());

                        for (int j = 0; j < coords.length; j++) {
                            tinkerXYZ.setCoordinateList1(coords[j], isTinkerOn);
                            Path txyzFileName = optXyzDir.resolve(particle 
                                    + "_" + j + ".txyz");
                            tinkerXYZ.writeToTxyzFile(txyzFileName.toString());
                        }

                        String searchString = MIPETUTIL.getResourceString(
                                "MIPETScanMinimumEnergyText");
                        List<String> minimumList = MIPETUTIL.findList(
                                scanOutputName.toString(), searchString);

                        if (!minimumList.isEmpty()) {
                            List<Double> energyValues = new LinkedList<>();

                            for (String s : minimumList) {
                                energyValues.add(Double.valueOf(s
                                        .substring(58, 68).trim()));
                            }

                            Double[] tmpEnergyList = energyValues.toArray(
                                    Double[]::new);
                            ArrayIndexComparator comparator = 
                                    new ArrayIndexComparator(tmpEnergyList);
                            Integer[] energyIndices = comparator
                                    .createIndexArray();
                            Arrays.sort(energyIndices, comparator);
                            Path optXyzPath = optXyzDir.resolve(particle 
                                    + "_0.txyz");
                            TinkerXYZ tinkerAfterScan = new TinkerXYZ(
                                    forcefield, particle, optXyzPath);
                            tinkerXYZ.setCoordinateList1(tinkerAfterScan
                                    .getCoordinateList1(), isTinkerOn);

                            // Convert tinker xyz files to .xyz file
                            for (int j = 0; j < energyValues.size(); j++) {
                                Path txyzPath = optXyzDir.resolve(
                                        particle + "_" + j + ".txyz");
                                tinkerXYZ0 = new TinkerXYZ(forcefield, particle,
                                        txyzPath);
                                Path particlePath = optXyzDir.resolve(particle 
                                        + "_o" + j + ".xyz");
                                tinkerXYZ0.writeToXyzFile(
                                        particlePath.toString());
                            }

                            // Copy configuration data with the lowest energy to Particle.xyz
                            Path bestEnergeySource = optXyzDir
                                    .resolve(particle + "_0.txyz");
                            try {
                                if (Files.exists(bestEnergeySource)) {
                                    Files.copy(bestEnergeySource, 
                                            targetXyzPath, 
                                            StandardCopyOption.REPLACE_EXISTING);
                                }
                            } catch (IOException ex) {
                                LOGGER.log(Level.SEVERE, "IOException during copying lowest energy .txyz file", ex);
                            }
                            
                            // --- RUN TINKER OPTIMIZE (AFTER) ---
                            Path afterOptLogName = optXyzDir.resolve(particle 
                                    + "_afteroptimize.log");
                            ProcessBuilder pbOptAfter = new ProcessBuilder(
                                    tinkerOptimize, 
                                    targetXyzPath.toString(), 
                                    Double.toString(optimizeRmsGradient));
                            pbOptAfter.redirectErrorStream(true);
                            pbOptAfter.redirectOutput(afterOptLogName.toFile());
                            try {
                                Process process = pbOptAfter.start();
                                process.waitFor();
                            } catch (IOException | InterruptedException ex) {
                                LOGGER.log(Level.SEVERE, "Exception during tinker's optimize.exe (after)", ex);
                            }
                            tinkerXYZ0 = new TinkerXYZ(forcefield, 
                                    particle, 
                                    optXyzDir.resolve(particle + ".xyz_2"));
                            tinkerXYZ.setCoordinateList1(tinkerXYZ0
                                    .getCoordinateList1(), isTinkerOn);
                            try {
                                Files.deleteIfExists(optXyzDir.resolve(
                                        particle + ".xyz_2"));
                            } catch (IOException ex) {
                                LOGGER.log(Level.SEVERE, "IOException during deleting .xyz_2 file.", ex);
                            }
                            try (PrintWriter writer = new PrintWriter(
                                    targetXyzPath.toFile())) {
                                writer.print(tinkerXYZ.getFileContent());
                            } catch (IOException ex) {
                                LOGGER.log(Level.SEVERE, "IOException during writing .xyz file.", ex);
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
        String particleName1 = aParticlePair.split("_")[0];
        String particleName2 = aParticlePair.split("_")[1];
        boolean isFractionOne = boltzmannFraction == 1.0;
        int distNumber = aDistances.length;
        int totalElements = aRotData2.length;
        int maxChunkSize = 1000;
        int minChunks = (totalElements + maxChunkSize - 1) / maxChunkSize;
        int chunkNumber = ((minChunks + cpuCoreNumber - 1) / cpuCoreNumber) 
                * cpuCoreNumber;
        int chunkSize = (totalElements + chunkNumber - 1) / chunkNumber; // Number of aRotData2 configurations
        
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
        
        TinkerXYZ tinkerXyz1 = aTinkerXYZ1.clone();
        TinkerXYZ tinkerXyz2 = aTinkerXYZ2.clone();
        Flat3DArray sharedFlatRot1 = Flat3DArray.createFrom(aRotData1);
        ExecutorService executor = Executors.newFixedThreadPool(cpuCoreNumber);
        List<MIPETAnalyze> taskList = new ArrayList<>(3000);
        String path = scratchDirectory 
                + FILESEPARATOR 
                + aParticlePair 
                + ".arc";
        int chunkIdx = 0;
        
        for (int i = 0; i < distNumber; i++) {
            double[][][] rotData2 = VectorUtil.moveX(aRotData2, aDistances[i]);
            
            for (int rot2StartIdx = 0; rot2StartIdx < totalElements; 
                    rot2StartIdx += chunkSize) {
                int rot2EndIndex = Math.min(rot2StartIdx + chunkSize, 
                        totalElements);
                TinkerXYZ tinkerXYZ = new TinkerXYZ(tinkerXyz1, tinkerXyz2, 
                        isTinkerOn);
                int atomNumber = tinkerXYZ.getN_atom();
                
                // This is for chunking rotData2 to avoid memory issues within the thread
                double[][][] rotPart2 = Arrays.copyOfRange(rotData2, 
                        rot2StartIdx, rot2EndIndex);
                Flat3DArray sharedFlatRot2 = Flat3DArray
                        .createFrom(rotPart2);
                String[] cmdList = new String[]{tinkerAnalyze, 
                    path + i + "_" + chunkIdx, "E"};
                taskList.add(new MIPETAnalyze(
                        tinkerXYZ, 
                        isTinkerOn, 
                        i, 
                        chunkIdx, 
                        atomNumber,
                        minAtomDistance, 
                        sharedFlatRot1, 
                        sharedFlatRot2, 
                        scratchDirectory,
                        cmdList, 
                        rec1,
                        rec2,
                        isFractionOne, 
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
        int taskIdx = 0;
        int distMinIdx = 0;
        int chunkMinIdx = 0;
        final double REZIP_RT = 1 / (temperature * GASCONST);
        double sumWgt = 0; // Sum of weights
        double sumWgtxE = 0; // Sum of (weight * energy)
        double distMinEnergy = 1E10;
        double allDistMinEnergy = 1E10;
        double wgtEmin = 100.;
        List<Future<WgtEnergyRecord>> futures = null;
        List<double[]> allResultArrays = new ArrayList<>();
        double[] energyMins = new double[distNumber];
        double[] wgtEmins = new double[distNumber];
        
        try {            
            futures = executor.invokeAll(taskList);
        } catch (InterruptedException ex) {
            Thread.currentThread().interrupt();       
        }
        executor.shutdown();
        Future<WgtEnergyRecord> future = null;
            
        for (int i = 0; i < distNumber; i++) {
            int totalSize = 0; // Size of energy values
                        
            for (int j = 0; j < chunkNumber; j++) {
                if (futures != null) {
                    future = futures.get(taskIdx);
                }
                try {
                    if (future != null) {
                        // Thread results
                        double[] res = future.get().energys();
                        allResultArrays.add(res);
                        totalSize += res.length;
                        distMinEnergy = future.get().energys()[0];
                        if (isFractionOne) {
                            sumWgt += future.get().sumWgt(); 
                            sumWgtxE += future.get().sumWgtxE();
                        }
                    }
                    if (distMinEnergy < allDistMinEnergy) {
                        allDistMinEnergy = distMinEnergy;
                        distMinIdx = i;
                        chunkMinIdx = j;
                    }
                    taskIdx++;
                } catch (InterruptedException | ExecutionException ex) {
                    LOGGER.log(Level.SEVERE,
                            "InterruptException during handling future object.",
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
            energyMins[i] = globalEnergies[0];
            
            // Store all energies at same distance to energyDatas[i]
            if (isFractionOne) {
                wgtEmins[i] = sumWgtxE / sumWgt;
                if (Double.isNaN(wgtEmins[i]) || 
                        Double.isInfinite(wgtEmins[i])) {
                    wgtEmins[i] = 100.;
                } 
                sumWgt = 0;
                sumWgtxE = 0;
            } else {
                int fractionToMax = (int)(globalEnergies.length 
                        * boltzmannFraction);
                double[] energyDataFraction = new double[fractionToMax];
                double[] weights = new double[fractionToMax];

                for (int j = 0; j < fractionToMax; j++) {
                    energyDataFraction[j] = globalEnergies[j];
                    weights[j] = Math.exp(-globalEnergies[j] * REZIP_RT);
                }

                wgtEmins[i] = MIPETUTIL.productSum(weights, 
                        energyDataFraction) / MIPETUTIL.sum(weights);
                if (Double.isNaN(wgtEmins[i]) 
                        || Double.isInfinite(wgtEmins[i])) {
                    wgtEmins[i] = 100.;
                } 
            }
            if (wgtEmin > wgtEmins[i]) {
                wgtEmin = wgtEmins[i];
            }
        }
        
        //</editor-fold>
        
        //<editor-fold defaultstate="collapsed" desc="Export .xyz file with the lowest intermolecular energy">
        double energyMin = allDistMinEnergy;
        if (allDistMinEnergy < aMinEnergy) {
            int minTaskIndex = chunkNumber * distMinIdx 
                    + chunkMinIdx;
            TinkerXYZ minTinkerXYZ = null;
            
            if (futures != null && minTaskIndex >= 0 
                    && minTaskIndex < futures.size()) {
                try {
                    minTinkerXYZ = futures.get(minTaskIndex).get()
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
        futures = null;
        
        //</editor-fold>

        return new EnergyRecord(aDistances,
                energyMins,
                wgtEmins, 
                aDistances[distMinIdx],
                energyMin,
                wgtEmin);
    }
    
    /**
     * Determine coordination numbers
     * @param aJobTaskRecordList: Job task records
     */
    private static void getCoordinationNumbers(
            ArrayList<JobTaskRecord> aJobTaskRecordList) {
        int jobTaskLength = aJobTaskRecordList.size();
        int jobSize = 0;
        Instant startTime = Instant.now();
        
        // <editor-fold defaultstate="collapsed" desc="Make directories and copy .xyz datas">
        for (int i = 0; i < jobTaskLength; i++) {
            if (aJobTaskRecordList.get(i).hasCNJob()) {
                jobSize++;
            }
        }
        
        double[] vdWSolventVolumes = new double[jobSize];
        double[] vdWSolutVolumes = new double[jobSize];
        double[] boxLengths = new double[jobSize];
        int jobIndex = 0;
        
        // Read xyzContent and prmContent
        if (!forceField_CN.equals(forceField_IE)) {
            readXyz(!isConformationalAnalysis, forceField_CN);
            if (forceField_CN.equals("OPLSAALIGPARGEN")) {
                readPrm(forceField_CN);
            }
        }
        
        boolean isSameParticle;
        String particle1;
        String particle2;        
        String particlePair;
        
        for (int i = 0; i < jobTaskLength; i++) {
            if (aJobTaskRecordList.get(i).hasCNJob()) {
                particle1 = aJobTaskRecordList.get(i).particleName1();
                particle2 = aJobTaskRecordList.get(i).particleName2();
                particlePair = particle1 + "_" + particle2;
                isSameParticle = aJobTaskRecordList.get(i).isSameParticle();
                int xyzID1 = particleNames.indexOf(particle1);
                int xyzID2;
                if (isSameParticle) {
                    xyzID2 = xyzID1;
                } else {
                    xyzID2 = particleNames.indexOf(particle2);
                }

                // Create directory in scratch directory
                Path scratchPath = Paths.get(scratchDirectory, 
                        forceField_CN,
                        particlePair);
                try {
                    Files.createDirectories(scratchPath);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE,
                        "IOException during creating scratch/forcefield/particlepair directory.",
                        ex);
                }
                
                // Copy .xyz files
                Path targetPath = Paths.get(scratchDirectory,
                        forceField_CN,
                        particlePair,
                        particle1 + ".xyz");
                try {
                    Files.writeString(targetPath, xyzContent1[xyzID1]);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during copying files to scratch.", ex);
                }
                if (!isSameParticle) {
                    targetPath = Paths.get(scratchDirectory,
                            forceField_CN,
                            particlePair,
                            particle2 + ".xyz");
                    try {
                        Files.writeString(targetPath, xyzContent1[xyzID2]);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during copying files to scratch.", ex);
                    }
                }
                vdWSolventVolumes[jobIndex] = MIPETUTIL
                        .getVdwVolume(smiles.get(particle2));
                if (isSameParticle) {
                    vdWSolutVolumes[jobIndex] = 
                            vdWSolventVolumes[jobIndex];
                } else {
                    vdWSolutVolumes[jobIndex] = MIPETUTIL.getVdwVolume(
                            smiles.get(particle1));
                }
                
                // Calculate water volume ratio - ratio of Vparticle and Vvdw 
                //   of water 1.7297
                boxLengths[jobIndex] =  Math.pow(WATERVOLUMERATIO 
                        * solventMoleculeNumber 
                        * vdWSolventVolumes[jobIndex], 1.0/3.0);
                jobIndex++;
            }
        }

        // </editor-fold>
        
        // <editor-fold defaultstate="collapsed" desc="Generate key file for dynamic">
        // OPENMP-THREADS is automatically set to core number when it is not set.
        // MIPET also works parallel. Tests showed, the best performance is
        //   achievable when OPENMP-THREADS is set to 1.
        String keyFixContent = 
            """
            EWALD
            OPENMP-THREADS 1
            THERMOSTAT ANDERSEN
            STEEPEST-DESCENT
            RANDOMSEED 123456
            """;
        if (solventMoleculeNumber >= 400) {
            keyFixContent += "NEIGHBOR-LIST" + LINESEPARATOR;
        }

        // </editor-fold>

        // <editor-fold defaultstate="collapsed" desc="Make simulation box">
        jobIndex = 0;

        for (int i = 0; i < jobTaskLength; i++) {
            if (aJobTaskRecordList.get(i).hasCNJob()) {
                boolean hasCNHeadLine = false;
                particle1 = aJobTaskRecordList.get(i).particleName1(); // solute
                particle2 = aJobTaskRecordList.get(i).particleName2(); // solvent
                int prm1ID = particleNames.indexOf(particle1);
                int prm2ID = particleNames.indexOf(particle2);
                particlePair = particle1 + "_" + particle2;
                isSameParticle = aJobTaskRecordList.get(i)
                        .isSameParticle();
                Path currPath = Paths.get(scratchDirectory, forceField_CN,
                        particlePair);
                Path resultPath = Paths.get(aJobTaskRecordList.get(i)
                        .result_CN_PathName());
                Path sourcePath = currPath.resolve(particlePair + ".xyz_2");
                Path sourcePath2 = currPath.resolve(particlePair + ".xyz_2");
                Path targetPath = currPath.resolve(particlePair + ".xyz");

                // <editor-fold defaultstate="collapsed" desc="Make .key file for dynamic">
                String forcefieldName;
                if ("OPLSAALIGPARGEN".equals(forceField_CN)) {
                    forcefieldName = "oplsaa";
                } else {
                    forcefieldName = forceField_CN.toLowerCase();
                }
                String keyContent = keyFixContent
                        + "PARAMETERS\t\""
                        + parameterDirectory
                        + "/"
                        + forcefieldName
                        + ".prm\""
                        + LINESEPARATOR    
                        + "a-axis "
                        + boxLengths[jobIndex]
                        + LINESEPARATOR;
                if (isRATTLE) {
                    keyContent += "RATTLE" + LINESEPARATOR;
                }
                Path keyPath = currPath.resolve(particlePair + ".key");
                if (forceField_CN.equals("OPLSAALIGPARGEN")) {
                    keyContent += prmContent1[prm1ID];
                    if (!isSameParticle) {
                        keyContent += prmContent2[prm2ID];
                    } 
                }
                MIPETUTIL.writeKeyFile(keyPath, keyContent);

                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="Build solvent box">
                ProcessBuilder pBuilder = new ProcessBuilder(tinkerXYZEdit, 
                        currPath.resolve(particle2 + ".xyz").toString(),
                        "-k",
                        keyPath.toString(),
                        "23",  // Option 23: Create and Fill a Periodic Boundary Box
                        Integer.toString(solventMoleculeNumber),
                        Double.toString(boxLengths[jobIndex]),
                        Double.toString(boxLengths[jobIndex]),
                        Double.toString(boxLengths[jobIndex]),
                        "Y");
                pBuilder.redirectErrorStream(true);
                if (isLogBuildBox) {
                    Path buildLogPath = resultPath.resolve(particle1 + "_" 
                            + particle2 + "_build.log");
                    pBuilder.redirectOutput(buildLogPath.toFile());
                } else {
                    // This loop is necessary for linux version.
                    pBuilder.redirectOutput(ProcessBuilder.Redirect.DISCARD);
                }
                Process process;
                try {
                    process = pBuilder.start();
                    process.waitFor();
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, "IOException during XYZEdit.", ex);
                } catch (InterruptedException ex) {
                    Thread.currentThread().interrupt();
                    LOGGER.log(Level.SEVERE, "Process was interrupted.", ex);
                } 

                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="Add solute to the simulation box">
                if (!isSameParticle) {
                    pBuilder = new ProcessBuilder(tinkerXYZEdit,
                        currPath.resolve(particle1 + ".xyz").toString(),
                        "-k",
                        keyPath.toString(),
                        "24", // Option 24: Soak Current Molecule in Box of Solvent
                        currPath + particle2 + ".xyz_2");
                    try {
                        process = pBuilder.start();
                        try (InputStream inStream = 
                                process.getInputStream()) {
                            if (isLogSoakBox) {
                                Path soakLogPath = resultPath.resolve(particle1
                                        + "_" + particle2 + "_soak.log");
                                Files.copy(inStream, soakLogPath, 
                                        StandardCopyOption.REPLACE_EXISTING);
                            } else {
                                // This loop is necessary for Linux version
                                inStream.transferTo(OutputStream.nullOutputStream());
                            } 
                        }
                        process.waitFor();
                    } catch (InterruptedException ex) {
                        Thread.currentThread().interrupt();
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during XYZEdit.", ex);
                    }
                    
                }
                Path oldXYZPath = currPath.resolve(particle1 + ".xyz_2");
                Path newXYZPath = currPath.resolve(particlePair + ".xyz");
                try {
                    Files.move(oldXYZPath, newXYZPath, 
                            StandardCopyOption.REPLACE_EXISTING);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                                "IOException during rename .xyz_2 to .xyz.",
                                ex);
                }

                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="Add MAXITER to .key file">
                long boxMinimizationTime = System.nanoTime();
                String keyMaxiter = "MAXITER " + minimizeMaxIteration;
                try (BufferedWriter writer = Files.newBufferedWriter(
                        keyPath, StandardCharsets.UTF_8, 
                        StandardOpenOption.APPEND)) {
                            writer.write(keyMaxiter);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during writing .key file.", ex);
                }

                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="Minimize solvent box">
                pBuilder = new ProcessBuilder(tinkerMinimize,
                        currPath.resolve(particlePair + ".xyz").toString(),
                        Double.toString(rmsMinimizeGradient));
                pBuilder.redirectErrorStream(true);
                if (isLogMinimizeBox) {
                    Path minimizeLogPath = resultPath.resolve(particlePair 
                            + "_minimize.log");
                    pBuilder.redirectOutput(minimizeLogPath.toFile());
                } else {
                    pBuilder.redirectOutput(ProcessBuilder.Redirect.DISCARD);
                }
                
                try {
                    process = pBuilder.start();
                    process.waitFor();
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, "IOException during XYZEdit.", ex);
                }catch (InterruptedException ex) {
                    Thread.currentThread().interrupt();
                    LOGGER.log(Level.SEVERE, "Process was interrupted", ex);
                } 
                String particleLogName;
                String datFileName;
                if (!aJobTaskRecordList.get(i).isReverse()) {
                    particleLogName = resultPath
                            + FILESEPARATOR
                            + particlePair
                            + "_log.txt";
                    datFileName = resultPath
                            + FILESEPARATOR
                            + particlePair
                            + ".dat";
                } else {
                    particleLogName = resultPath
                            + FILESEPARATOR
                            + particle2
                            + "_"
                            + particle1
                            + "_log.txt";
                    datFileName = resultPath
                            + FILESEPARATOR
                            + particle2
                            + "_"
                            + particle1
                            + ".dat";
                }

                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="Write .log file">
                Path particleLogPath = Paths.get(particleLogName);
                if (!isSameParticle && !Files.exists(particleLogPath)) {
                    hasCNHeadLine = true;
                }
                try {
                    Path datFilePath = Paths.get(datFileName);
                    BWParticleDat = Files.newBufferedWriter(datFilePath, 
                            StandardCharsets.UTF_8);
                    
                    // BufferedWriter for log contents
                    particleLogPath = Paths.get(particleLogName);
                    try (BufferedWriter BWParticleLog = Files.newBufferedWriter(
                            particleLogPath, StandardCharsets.UTF_8)) {
                        if (Files.exists(sourcePath2) && !hasCNHeadLine) {
                            Files.move(sourcePath, targetPath,
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
                    if (!isSameParticle) {
                        BWParticleDat.append("VdWVolume(")
                                .append(particle1)
                                .append(") [" + ANGSTROM + SUPERSCRIPTTHREE 
                                        +"] = ")
                                .append(Double.toString(vdWSolutVolumes[jobIndex]))
                                .append(LINESEPARATOR);
                    }
                    BWParticleDat.append("VdWVolume(")
                            .append(particle2)
                            .append(") [" + ANGSTROM + SUPERSCRIPTTHREE 
                                    +"] = ")
                            .append(String.format("%.4f", 
                                    vdWSolventVolumes[jobIndex]))
                            .append(LINESEPARATOR);
                    BWParticleDat.close();
                } catch(IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during Minimize.", ex);
                }

                // Make backup of particle1_particle2.xyz
                sourcePath = targetPath;
                targetPath = currPath.resolve(particlePair + ".bak");
                try{
                    Files.copy(sourcePath, targetPath, 
                            StandardCopyOption.REPLACE_EXISTING);
                } catch(IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during Copying to .bak file", ex);
                }
                jobIndex++;

                // </editor-fold>
            }
        }

        // </editor-fold>
        
        // <editor-fold defaultstate="collapsed" desc="Simulation warm up">
        String[] cmdList;
        LinkedList<MIPETCN> taskList = new LinkedList<>();
        
        for (int i = 0; i < jobTaskLength; i++) {
            if (aJobTaskRecordList.get(i).hasCNJob()) {
                particle1 = aJobTaskRecordList.get(i).particleName1(); // solute
                particle2 = aJobTaskRecordList.get(i).particleName2(); // solvent
                particlePair = particle1 + "_" + particle2;
                Path keyPath = Paths.get(scratchDirectory, forceField_CN, 
                        particlePair);
                if (isTinker9) {
                    cmdList = new String[] {tinkerDynamic,
                        "dynamic",
                        keyPath.resolve(particlePair + ".xyz").toString(),
                        Integer.toString(warmUpStepNumber),
                        Double.toString(warmUpTimeStep),
                        Double.toString(warmUpPrintInterval),
                        Integer.toString(simulationType),
                        Double.toString(temperature)};
                } else {
                    cmdList = new String[] {tinkerDynamic,
                        keyPath.resolve(particlePair + ".xyz").toString(),
                        Integer.toString(warmUpStepNumber),
                        Double.toString(warmUpTimeStep),
                        Double.toString(warmUpPrintInterval),
                        Integer.toString(simulationType),
                        Double.toString(temperature)};
                }
                taskList.add(new MIPETCN(cmdList,
                        aJobTaskRecordList.get(i),
                        scratchDirectory,
                        catchRadius,
                        true,
                        isLogWarmup,
                        isLogDynamic,
                        isTinker9));    
            }
        }

        long dynamicWarmUpTime = System.nanoTime();
        ExecutorService executor = Executors.newFixedThreadPool(cpuCoreNumber);
        try {
            executor.invokeAll(taskList);
        } catch (InterruptedException ex) {
            Thread.currentThread().interrupt();       
        }
        executor.shutdown();
        String[][] labelValues = new String[1][2];
        labelValues[0][0] = "Time for solvent box warm up via Tinker \"Dynamic\" [s]: ";
        labelValues[0][1] = Double.toString((System.nanoTime() 
                - dynamicWarmUpTime) / 1000d);
        MIPETUtility.writeParticleLog(aJobTaskRecordList, labelValues);
        
        // </editor-fold>
        
        // <editor-fold defaultstate="collapsed" desc="Start simulation & determine coordination number">
        LinkedList<MIPETCN> cnTaskList = new LinkedList<>();

        for (int i = 0; i < jobTaskLength; i++) {
            if (aJobTaskRecordList.get(i).hasCNJob()) {
                particle1 = aJobTaskRecordList.get(i).particleName1(); // solute
                particle2 = aJobTaskRecordList.get(i).particleName2(); // solvent
                particlePair = particle1 + "_" + particle2;
                Path currPath = Paths.get(scratchDirectory, forceField_CN,
                        particlePair);
                if (isTinker9) {
                    cmdList = new String[] {tinkerDynamic,
                        "dynamic",
                        currPath.resolve(particlePair + ".xyz").toString(),
                        Integer.toString(stepNumber),
                        Double.toString(timeStep),
                        Double.toString(printInterval),
                        Integer.toString(simulationType),
                        Double.toString(temperature)};
                } else {
                    cmdList = new String[] {tinkerDynamic,
                        currPath.resolve(particlePair + ".xyz").toString(),
                        Integer.toString(stepNumber),
                        Double.toString(timeStep),
                        Double.toString(printInterval),
                        Integer.toString(simulationType),
                        Double.toString(temperature)};
                }
                cnTaskList.add(new MIPETCN(cmdList,
                        aJobTaskRecordList.get(i),
                        scratchDirectory,
                        catchRadius,
                        false,
                        isLogWarmup,
                        isLogDynamic,
                        isTinker9));
            }
        }
        
        long arcEvaluationTime = System.nanoTime();
        List<Future<int[]>> cnFutures = new LinkedList<>();
        if (jobSize > 0) {
            executor = Executors.newFixedThreadPool(cpuCoreNumber);
            try {
                cnFutures = executor.invokeAll(cnTaskList);
            } catch (InterruptedException ex) {
                Thread.currentThread().interrupt();       
            }
            executor.shutdown();
            labelValues = new String[2][2];
            labelValues[0][0] = "Time for solvent box analysis (includes \".arc-evaluation\" and \"MIPETCN\") [s]: ";
            labelValues[0][1] = Double.toString(
                    (System.nanoTime() - arcEvaluationTime) / 1000d);
            labelValues[1][0] = "Total time for determining the coordination number [s]: ";
            Instant endTime = Instant.now();
            Duration duration = Duration.between(startTime, endTime);
            float totalTimeSeconds = duration.toMillis() /1000f;
            labelValues[1][1] = String.valueOf(totalTimeSeconds);
            MIPETUtility.writeParticleLog(aJobTaskRecordList, labelValues);
            String labelStr = "Mean neighbor (" + temperature + " K): ";
            String[] valueStrs = new String[cnFutures.size()];
            int[][] cns = new int[cnFutures.size()][];
            double[] cnMeans = new double[cnFutures.size()];
            double[] stdDeviation = new double[cnFutures.size()];
            int[] cnMaxes = new int[cnFutures.size()];
            int[] cnMins = new int[cnFutures.size()];
            jobIndex = 0;

            for (int i = 0; i < jobTaskLength; i++) {
                if (aJobTaskRecordList.get(i).hasCNJob()) {
                    particle1 = aJobTaskRecordList.get(i).particleName1(); // solute
                    particle2 = aJobTaskRecordList.get(i).particleName2(); // solvent
                    particlePair = particle1 + "_" + particle2;
                    try {
                        cns[jobIndex] = cnFutures.get(jobIndex).get();
                    } catch (ExecutionException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "ExecutionException during get CN", ex);
                    } catch (InterruptedException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "InterruptedException during get CN", ex);
                    }
                    cnMeans[jobIndex] = (double)MIPETUTIL.sum(
                            cns[jobIndex]) / cns[jobIndex].length;
                    cnList.add(new ResultCNRecord(
                            particle1, 
                            particle2, 
                            cnMeans[jobIndex]));
                    stdDeviation[jobIndex] = MIPETUTIL
                            .standarddeviation(cns[jobIndex], 
                                    cnMeans[jobIndex]);
                    cnMins[jobIndex] = MIPETUTIL.getMin(cns[jobIndex]);
                    cnMaxes[jobIndex] = MIPETUTIL.getMax(cns[jobIndex]);
                    valueStrs[jobIndex] = Double.toString(cnMeans[jobIndex]);
                    Path currPath = Paths.get(scratchDirectory, forceField_CN,
                            particlePair);
                    Path sourceFile = currPath.resolve(particlePair + ".arc");
                    try {
                        Files.deleteIfExists(sourceFile);
                    } catch (IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                                "IOException during deleting .arc file", ex);
                    }

                    // Restore of particle1_particle2.xyz
                    Path sourcePath = currPath.resolve(particlePair + ".bak");
                    Path targetPath = currPath.resolve(particlePair + ".xyz");
                    try {
                        Files.copy(sourcePath, targetPath, 
                        StandardCopyOption.REPLACE_EXISTING);
                    } catch(IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                        "IOException during Copying to .bak file", ex);
                    }
                    jobIndex++;
                }
            }

            // </editor-fold>

            // <editor-fold defaultstate="collapsed" desc="Write data">

            MIPETUTIL.writeParticleLog(aJobTaskRecordList, labelStr, valueStrs);
            MIPETUTIL.writeZij_Table(aJobTaskRecordList, cns, temperature);
            jobIndex = 0;

            for (int i = 0; i < jobTaskLength; i++) {
                if (aJobTaskRecordList.get(i).hasCNJob()) {
                    particle1 = aJobTaskRecordList.get(i).particleName1(); // solute
                    particle2 = aJobTaskRecordList.get(i).particleName2(); // solvent
                    if (!aJobTaskRecordList.get(i).isReverse()) {
                        particlePair = particle1 + "_" + particle2;
                    } else {
                        particlePair = particle2 + "_" + particle1;
                    }
                    Path datPath = Paths.get(aJobTaskRecordList.get(i)
                            .result_CN_PathName(), particlePair + ".dat");
                    try (BufferedWriter writer = Files.newBufferedWriter(
                            datPath, StandardCharsets.UTF_8)) {
                        // Coordination number mean
                        writer.append("CNmean(")
                                .append(particle1)
                                .append("/")
                                .append(particle2)
                                .append(") = ")
                                .append(String.format("%.2f", 
                                        cnMeans[jobIndex]));
                        writer.append("    ");
                        writer.newLine();

                        // Standard deviation
                        writer.append("Standard deviation(")
                                .append(particle1)
                                .append("/")
                                .append(particle2)
                                .append(") = ")
                                .append(String.format("%.2f", 
                                        stdDeviation[jobIndex]));
                        writer.append("    ");
                        writer.newLine();

                        // CN Min
                        writer.append("Min(")
                                .append(particle1)
                                .append("/")
                                .append(particle2)
                                .append(") = ")
                                .append(Integer.toString(cnMins[jobIndex]))
                                .append("    ");
                        writer.newLine();

                        // CN Max
                        writer.append("Max(")
                                .append(particle1)
                                .append("/")
                                .append(particle2)
                                .append(") = ")
                                .append(Integer.toString(cnMaxes[jobIndex]));
                        writer.append("    ");
                        writer.newLine();

                        // Box length
                        writer.append("BoxLength [" + ANGSTROM + "] = ")
                                .append(String.format("%.4f", 
                                        boxLengths[jobIndex]));
                        writer.newLine();
                        writer.close();
                        jobIndex++;
                    } catch(IOException ex) {
                        LOGGER.log(Level.SEVERE, 
                        "IOException during Writing .dat file.", ex);
                    }
                }
            }

            // </editor-fold>

            // <editor-fold defaultstate="collapsed" desc="Clear scratch directory">
            File folder = new File(scratchDirectory);
            try {
                FileUtils.deleteDirectory(folder);
            } catch(IOException ex) {
                LOGGER.log(Level.SEVERE, 
                        "IOException during deleting scratch directory.", ex);
            }Path scratchPath = Paths.get(scratchDirectory);
            try {
                if (Files.exists(scratchPath)) {
                    try (Stream<Path> walk = Files.walk(scratchPath)) {
                        walk.sorted(Comparator.reverseOrder()).forEach(path -> {
                            try {
                                Files.delete(path);
                            } catch (IOException e) {
                                LOGGER.log(Level.WARNING, () -> 
                                        "Couldn't delete: " 
                                                + path 
                                                + " - Reason: " 
                                                + e.getMessage());
                            }
                        });
                    }
                }
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE, "IOException beim Durchsuchen des Scratch-Verzeichnisses.", ex);
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
        
        // 0: CN = 1, Wgt = true, Opt = true, Rgd = false
        // 1: CN = 1, Wgt = false, Opt = true, Rgd = false
        // 2: CN = 1, Wgt = false, Opt = true, Rgd = true
        // 3: CN = 1, Wgt = false, Opt = false
        // 4: CN != 1, Wgt = true, Opt = true, Rgd = false
        // 5: CN != 1, Wgt = false, Opt = true, Rgd = false
        // 6: CN != 1, Wgt = false, Opt = true, Rgd = true
        // 7: CN != 1, Wgt = false, Opt = false
        int jobSize = aJobTaskRecords.size();

        // particle description
        List<String> particleDescriptionString = Arrays.asList(
                "# Particle", 
                "Name", 
                "Mass[DPD]", 
                "Charge", 
                "Mass[g/mol]", 
                "Volume[A^3]", 
                "Graphics-Radius", 
                "Standard-Color");
        ArrayList<ArrayList<String>> particleDescriptionsList = 
                new ArrayList<>(particleNames.size());
        
        for(String particleName : particleNames){
            ArrayList<String> particleDescriptions = new ArrayList<>(8);
            particleDescriptions.add(particleName); //particle
            particleDescriptions.add(particleName); //full name
            particleDescriptions.add(Integer.toString(MASS_DPD)); //massDPD
            particleDescriptions.add(Integer.toString(DEFAULT_CHARGE)); // charge
            try {
                particleDescriptions.add(String.format("%.4f", 
                        MIPETUTIL.getAtomicMass(smiles.get(
                                particleName), true))); // mass [g/mol]
            } catch (Exception ex) {
                particleDescriptions.add("-1");
                LOGGER.log(Level.SEVERE, ex.toString());
            } 
            if (!particleName.equals("[Na+]")) {
                particleDescriptions.add(String.format("%.4f", 
                    MIPETUTIL.getVdwVolume(smiles.get(particleName)) 
                            * WATERVOLUMERATIO)); // volume [A³]
            }
            particleDescriptions.add(Double.toString(GRAPHICS_RADIUS)); // graphics-radius
            particleDescriptions.add(DEFAULT_COLOR); // color
            particleDescriptionsList.add(particleDescriptions);
        }
        
        // particle interactions
        
        String particleName1;
        String particleName2;
        String particlePair;
        
        int energyListSize = energyList.size();
        int cnListSize = cnList.size();
        HashMap<String, Double> aijMap = new HashMap<>(jobSize);
        HashMap<String, Double> aijMap1 = new HashMap<>(jobSize); // for CN =1
        HashMap<String, Double> energiesMap = new HashMap<>(energyListSize);
        HashMap<String, Double> cnMap = new HashMap<>(cnListSize);
        Set<String> keySet = new HashSet<>();
        double energy = 0.;
        int outputIteration = 0;
        
        while (outputIteration <= 7) {
            // Read energy data
            for (int i = 0; i < energyListSize; i++) {
                particleName1 = energyList.get(i).particleName1();
                particleName2 = energyList.get(i).particleName2();
                particlePair = particleName1 + "_" + particleName2;
                switch (outputIteration) {
                    case 0, 4 -> energy = energyList.get(i)
                            .wgt_Emin();
                    case 1, 5 -> energy = energyList.get(i)
                            .wgt0_Opt_Emin();
                    case 2, 6 -> energy = energyList.get(i)
                            .wgt0_Rgd_Emin();
                    case 3, 7 -> energy = energyList.get(i)
                            .wgt0_Opt0_Emin();
                }
                energiesMap.put(particlePair, energy);
            }

            // Read coordination number
            for (int i = 0; i < cnListSize; i++) {
                particleName1 = cnList.get(i).particleName1();
                particleName2 = cnList.get(i).particleName2();
                particlePair = particleName1 + "_" + particleName2;
                double cn = cnList.get(i).cnValue();
                cnMap.put(particlePair, cn);
            }
            
            // Calculate aij parameters (see DPD theory pdf)
            double aij = 0.;

            if (forceField_IE != null) {
                
                for (int i = 0; i < jobSize; i++) {
                    particleName1 = aJobTaskRecords.get(i).particleName1();
                    particleName2 = aJobTaskRecords.get(i).particleName2();
                    if (!aJobTaskRecords.get(i).isReverse()) {
                        double e12 = energiesMap.get(particleName1 + "_"
                                + particleName2);
                        double e11 = energiesMap.get(particleName1 + "_"
                                + particleName1);
                        double e22 = energiesMap.get(particleName2 + "_"
                                + particleName2);
                        double chiNumerator;
                        if (forceField_CN != null) {
                            double z11 = cnMap.get(particleName1 + "_"
                                    + particleName1);
                            double z22 = cnMap.get(particleName2 + "_"
                                    + particleName2);
                            double z12 = cnMap.get(particleName1 + "_"
                                    + particleName2);
                            double z21 = cnMap.get(particleName2 + "_"
                                    + particleName1);
                            chiNumerator
                                    = z12 * e12
                                    + z21 * e12
                                    - z11 * e11
                                    - z22 * e22;
                            aij = temperature / 12 + 1.7483
                                    * chiNumerator;
                            aijMap.put(particleName1 + "_" + particleName2, 
                                    aij);
                        }

                        // Calculation for CN = 1
                        chiNumerator = e12 + e12 - e11 - e22;
                        aij = temperature / 12 + 1.7483 * chiNumerator;
                        aijMap1.put(particleName1 + "_" + particleName2, aij);
                    }
                }
                
            }
            
            /* Write file */
            String resultsDirectory = resultDirectory 
                + FILESEPARATOR
                + "IE"
                + FILESEPARATOR
                + forceField_IE;
            String fileName = resultsDirectory 
                + FILESEPARATOR 
                + aTitleAbbreviation + "_"
                + "EijFraction_"
                + boltzmannFraction
                + "_catchRadius_" 
                + catchRadius;
            switch (outputIteration) {
                case 0 -> fileName += "_CN1_Wgt.txt";
                case 1 -> fileName += "_CN1_Wgt0_Opt.txt";
                case 2 -> fileName += "_CN1_Wgt0_Rgd.txt";
                case 3 -> fileName += "_CN1_Wgt0_Opt0.txt";
                case 4 -> fileName += "_CN_Wgt.txt";
                case 5 -> fileName += "_CN_Wgt0_Opt.txt";
                case 6 -> fileName += "_CN_Wgt0_Rgd.txt";
                case 7 -> fileName += "_CN_Wgt0_Opt0.txt";
            }
            if (forceField_CN != null || outputIteration <= 3) {
                try (BufferedWriter writer = Files.newBufferedWriter(
                        Paths.get(fileName), StandardCharsets.UTF_8)) {
                    writer.append("# Particle set for MFSim created by ")
                            .append(MIPET.class.getPackage()
                                    .getImplementationTitle())
                            .append(" ")
                            .append(MIPET.class.getPackage()
                                    .getImplementationVersion())
                            .append(LINESEPARATOR)
                            .append("# Force Field for energy calculation: ")
                            .append(forceField_IE)
                            .append(LINESEPARATOR)
                            .append("# Force Field for coordination number calculation: ");
                    if (forceField_CN == null) {
                        writer.append("Not calculated");
                    } else {
                        writer.append(forceField_CN);
                    }
                    writer.append(forceField_CN)
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
                    switch (outputIteration) {
                        case 0, 4 -> {
                            writer.append("# Boltzmann averaging: Yes")
                                    .append(LINESEPARATOR);}
                        case 1, 2, 3, 5, 6, 7 -> {
                            writer.append("# Boltzmann averaging: No")
                                    .append(LINESEPARATOR);}
                    }
                    switch (outputIteration) {
                        case 1, 6 -> {
                            writer.append("# Optimize sampled E(min) configuration: Yes")
                                    .append(LINESEPARATOR);}
                        case 0, 2, 3, 4, 5, 7 -> {
                            writer.append("# Optimize sampled E(min) configuration: No")
                                    .append(LINESEPARATOR);}
                    }
                    switch (outputIteration) {
                        case 2, 6 -> {
                            writer.append("# Tinker's 'optrigid' used: Yes")
                                    .append(LINESEPARATOR);}
                    }
                    switch (outputIteration) {
                        case 0, 1, 2, 3 -> {
                            writer.append("# CN = 1 for all particle pairs")
                                    .append(LINESEPARATOR);}
                    }
                    writer.append(LINESEPARATOR)
                            .append("[Title]\n")
                            .append(aTitle)
                            .append("\n[/Title]\n\n")
                            .append("[Version]\n" 
                                    + VERSION_NUMBER 
                                    + "\n[/Version]\n\n");

                    // Particle description
                    writer.append("[Particle Description]")
                            .append(LINESEPARATOR);
                    writer.append("""
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
                    writer.append(LINESEPARATOR);

                    for(String header : particleDescriptionString){
                        writer.append(header)
                                .append("\t");
                    }

                    for(List<String> particleDescList : 
                            particleDescriptionsList){
                        writer.append(LINESEPARATOR);

                        for(String description : particleDescList){
                            writer.append(description)
                                    .append(" ");
                        }

                    }

                    writer.append(LINESEPARATOR)
                            .append("[/Particle Description]")
                            .append(LINESEPARATOR)
                            .append(LINESEPARATOR);

                     // Particle interactions
                    writer.append("[Particle interactions]")
                            .append(LINESEPARATOR)
                            .append("# Repulsion parameters a(ij) for particle pairs for different temperatures (in K)")
                            .append(LINESEPARATOR)
                            .append("Pair")
                            .append(" ")
                            .append(Integer.toString((int) temperature));
                    switch (outputIteration) {
                        case 0, 1, 2, 3 -> {
                            keySet = aijMap1.keySet();}
                        case 4, 5, 6, 7 -> {
                            keySet = aijMap.keySet();}
                    }

                    for(String key : keySet) {
                        writer.append(LINESEPARATOR)
                                .append(key);
                        switch (outputIteration) {
                            case 0, 1, 2, 3 -> {
                                aij = aijMap1.get(key);}
                            case 4, 5, 6, 7 -> {
                                aij = aijMap.get(key);}
                        }
                        writer.append(" ").append(String.format("%.2f", aij));
                    }

                    writer.append(LINESEPARATOR)
                            .append("[/Particle interactions]")
                            .append(LINESEPARATOR)
                            .append(LINESEPARATOR);

                    // Coordination numbers
                    switch (outputIteration) {
                        case 4, 5, 6, 7 -> {
                            writer.append("[Coordination numbers]");
                            cnListSize = cnList.size();

                            for (int j = 0; j < cnListSize; j++) {
                                writer.append(LINESEPARATOR);
                                writer.append(cnList.get(j).particleName1())
                                        .append("_")
                                        .append(cnList.get(j).particleName2());
                                writer.append(String.format(" %.2f", 
                                        cnList.get(j).cnValue()));
                            }

                            writer.append(LINESEPARATOR)
                                    .append("[/Coordination numbers]")
                                    .append(LINESEPARATOR);}
                    }

                    // SMILES
                    writer.append(LINESEPARATOR)
                            .append("[SMILES]")
                            .append(LINESEPARATOR)
                            .append("# Particle and corresponding SMILES of fragment molecule");
                    keySet = smiles.keySet();

                    for(String key : keySet){
                        if (particleNames.contains(key)) {
                            writer.append(LINESEPARATOR)
                                    .append(key)
                                    .append(" ")
                                    .append(smiles.get(key));
                        }
                    }

                    writer.append(LINESEPARATOR)
                            .append("[/SMILES]");
                } catch (IOException anException) {
                    LOGGER.log(Level.SEVERE, anException.toString());
                }
            }
            if (isTinkerOn) {
                outputIteration++;
            } else {
                outputIteration = 8;
            }
        }
            
    }
    
    // </editor-fold>
}