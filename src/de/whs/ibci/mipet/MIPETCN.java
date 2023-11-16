/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package de.whs.ibci.mipet;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.LinkedList;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Stream;

/**
 * MIPETCN class
   This class helps to process determine coordination number (CN) parallaly.
 * 
 * @author Mirco Daniel
 */
public class MIPETCN implements Callable<int[]> {
    // <editor-fold defaultstate="collapsed" desc="Final Class variables">
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
     * Flag for warmup
     */
    private final Boolean ISWARMUP;
    
    /**
     * Command list for tinker's analyze.exe
     */
    private final String[] COMMAND_LIST;
    
    /**
     * Parameter and particle names 
     */
    private final JobTaskRecord JOBTASK_RECORD;
    
    /**
     * Scratch directory name
     */
    private final String SCRATCH_DIR;
    
    /**
     * File separator
     */
    private final String FILESEPARATOR = System.getProperty("file.separator");
    
    /**
     * Catch radius in Angstrom
     */
    private final double CATCH_RADIUS;
    
    /**
     * Flag for write warmup log file
     */
    private final boolean ISLOGWARMUP;
    
    /**
     * Flag for write dynamic log file
     */
    private final boolean ISLOGDYNAMIC;
    
    /**
     * Flag if Tinker9 is used or not
     */
    private final boolean ISTINKER9;
    
    // </editor-fold> 
    
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    
    public MIPETCN(String[] aCommandList,
            JobTaskRecord aJobTaskRecord, 
            String aScratchDir,
            double aCatchRadius,
            boolean anIsWarmup,
            boolean anIsLogWarmup,
            boolean anIsLogDynamic,
            boolean anIsTinker9) {
        this.COMMAND_LIST = aCommandList;
        this.JOBTASK_RECORD = aJobTaskRecord;
        this.SCRATCH_DIR = aScratchDir;
        this.CATCH_RADIUS = aCatchRadius;
        this.ISWARMUP = anIsWarmup;
        this.ISLOGWARMUP = anIsLogWarmup;
        this.ISLOGDYNAMIC = anIsLogDynamic;
        this.ISTINKER9 = anIsTinker9;
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    @Override
    public int[] call() {
        boolean tmpIsSameParticle;
        int tmpAtomNumber1;
        int tmpAtomNumber2;
        int tmpIteration;
        int tmpStep;
        int tmpStepsPerRound;
        int tmpRestSteps;
        int tmpStepNumber;
        int tmpStepNumberOfLastIteration;
        int[] tmpResultNeighborNumber;
        LinkedList<int[]> tmpCurrentNeighbors;
        LinkedList<Integer> tmpNeighborNumbers;
        double tmpBoxLength;
        String tmpForcefield;
        String tmpParticle1;
        String tmpParticle2;
        String tmpParticlePair;
        String tmpCurrentDir;
        String tmpSourceName;
        String tmpTargetName;
        String tmpTargetXyzName;
        String tmpResultPathName;
        String tmpXyzFileName2;
        String tmpWarmUpLogName;
        String tmpDynamicLogName;
        String tmpArcName;
        String tmpXyzName;
        String[] tmpCommandList;
        String[] tmpElements1;
        String[] tmpElements2;
        TinkerXYZ tmpTxyz1;
        TinkerXYZ tmpTxyz2;
        Path tmpSource;
        Path tmpTarget;
        ProcessBuilder tmpPBuilder;
        Process tmpProcess;
        CoordinatesRecord tmpCoordRecord;
        
        tmpNeighborNumbers = new LinkedList<>();
        MIPETUtility MIPET4JUtil = new MIPETUtility();
        tmpStepsPerRound = 2000;
        tmpForcefield = this.JOBTASK_RECORD.forcefield_CN_Name();
        tmpParticle1 = this.JOBTASK_RECORD.particleName1(); // solute
        tmpParticle2 = this.JOBTASK_RECORD.particleName2(); // solvent
        tmpIsSameParticle = tmpParticle1.equals(tmpParticle2);
        tmpParticlePair = tmpParticle1 + "_" + tmpParticle2;
        tmpResultPathName = this.JOBTASK_RECORD.result_CN_PathName();
        tmpCurrentDir = this.SCRATCH_DIR
                    + FILESEPARATOR
                    + tmpForcefield
                    + FILESEPARATOR
                    + tmpParticlePair
                    + FILESEPARATOR;
        tmpTxyz1 = new TinkerXYZ(tmpCurrentDir + tmpParticle1 + ".xyz");
        tmpBoxLength = MIPET4JUtil.getBoxLength(tmpCurrentDir 
                + tmpParticlePair + ".xyz");
        tmpAtomNumber1 = tmpTxyz1.getAtomNumber();
        tmpElements1 = tmpTxyz1.getElementList1();
        if (tmpIsSameParticle) {
            tmpAtomNumber2 = tmpAtomNumber1;
            tmpElements2 = tmpElements1.clone();
        } else {
            tmpXyzFileName2 = tmpCurrentDir + tmpParticle2 + ".xyz";
            tmpTxyz2 = new TinkerXYZ(tmpXyzFileName2);
            tmpAtomNumber2 = tmpTxyz2.getAtomNumber();
            tmpElements2 = tmpTxyz2.getElementList1();
        }
        if (this.ISTINKER9) {
            tmpStepNumber = Integer.parseInt(this.COMMAND_LIST[3]);
        } else {
            tmpStepNumber = Integer.parseInt(this.COMMAND_LIST[2]);
        }
        tmpRestSteps = tmpStepNumber % tmpStepsPerRound;
        tmpCommandList = this.COMMAND_LIST.clone();
        if (tmpStepNumber > tmpStepsPerRound) {
            tmpIteration = (int)Math.ceil((double)tmpStepNumber 
                    / tmpStepsPerRound);
            tmpStep = tmpStepsPerRound;
            if (this.ISTINKER9) {
                tmpCommandList[3] = Integer.toString(tmpStep);
            } else {
                tmpCommandList[2] = Integer.toString(tmpStep);
            }
            if ( tmpRestSteps == 0) {
                tmpStepNumberOfLastIteration = tmpStep;
            } else {
                tmpStepNumberOfLastIteration = tmpRestSteps;
            }
        } else {
            tmpIteration = 1;
            tmpStep = tmpStepNumber;
            tmpStepNumberOfLastIteration = tmpStep;
        }
        if(ISWARMUP) {
            tmpPBuilder = new ProcessBuilder();
            tmpPBuilder.command(this.COMMAND_LIST);
            try {
                tmpProcess = tmpPBuilder.start();
                try (BufferedReader tmpBR = new BufferedReader(
                        new InputStreamReader(tmpProcess.getInputStream()))) {
                    if(ISLOGWARMUP) {
                        tmpWarmUpLogName = tmpResultPathName
                        + FILESEPARATOR
                        + tmpParticle1
                        + "_"
                        + tmpParticle2
                        + "_warmUp.log";
                        BufferedWriter tmpBW = new BufferedWriter(
                                new FileWriter(tmpWarmUpLogName));
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
            } catch(IOException | InterruptedException ex) {
                LOGGER.log(Level.SEVERE,
                        "Exception during tinker's dynamic.exe", ex);
            }
            
            // Rename .arc to .xyz
            tmpSourceName = tmpCurrentDir + tmpParticlePair + ".arc";
            tmpSource = Paths.get(tmpSourceName);
            tmpTargetXyzName = tmpCurrentDir + tmpParticlePair + ".xyz";
            tmpTarget = Paths.get(tmpTargetXyzName);
            try {
                Files.move(tmpSource, tmpTarget, 
                        StandardCopyOption.ATOMIC_MOVE);
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE, 
                        "IOException during remove .arc file to .xyz file.",
                        ex);
            }
            
            // Write _lastStepNeighbors.txt file
            tmpCoordRecord = MIPET4JUtil.getCoordinatesFromArcFile(
                    tmpTargetXyzName, tmpAtomNumber1, tmpAtomNumber2);
            tmpCurrentNeighbors = MIPET4JUtil.getNeighborNumbersBruteForce(
                    tmpCoordRecord, 
                    tmpElements1, 
                    tmpElements2,
                    tmpBoxLength,
                    this.CATCH_RADIUS);
            tmpTargetName = tmpResultPathName 
                    + FILESEPARATOR
                    + tmpParticlePair 
                    + "_warmUpNeighbors.txt";
            int tmpIterationSize = tmpCurrentNeighbors.size();
            int tmpNeighborsSize = tmpCurrentNeighbors
                    .get(tmpIterationSize - 1).length;
            
            try (BufferedWriter tmpBW = new BufferedWriter(
                    new FileWriter(tmpTargetName))) {
                for (int i = 0; i < tmpNeighborsSize; i++) {
                    tmpBW.append(String.valueOf(tmpCurrentNeighbors
                            .get(tmpIterationSize - 1)[i]));
                    tmpBW.append(LINESEPARATOR);
                }
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE, 
                        "IOException during writing _warmUpNeighbors.txt file.",
                        ex);
            }
            
            // Copy .xyz file after warmup to result directory
            tmpSource = Paths.get(tmpTargetXyzName);
            tmpTarget = Paths.get(tmpResultPathName 
                    + FILESEPARATOR
                    + tmpParticlePair 
                    + "_warmUpCoords.xyz");
            try {
                Files.copy(tmpSource, tmpTarget, StandardCopyOption
                        .REPLACE_EXISTING);
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE, 
                        "IOException during copying .xyz file.", ex);
            }
            
            // Delete useless files
            try (Stream<Path> tmpList = Files.list(Paths.get(tmpCurrentDir))) {
                tmpList.filter(file -> !Files.isDirectory(file))
                        .filter(file -> file.endsWith(".xyz_2"))
                        .map(Path::toFile)
                        .forEach(File::delete);
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE, 
                        "IOException during deleting .xyz_2 file.", ex);
            }
            return null;
        }
        int tmpIterationSize;
        int tmpNeighborSize;
        int tmpCurrentNeighborsSize;
        
        tmpDynamicLogName = tmpResultPathName
                        + FILESEPARATOR
                        + tmpParticle1
                        + "_"
                        + tmpParticle2
                        + "_dynamic.log";
        
        for (int i = 0; i < tmpIteration; i++) {
            tmpPBuilder = new ProcessBuilder();
            if (i == tmpIteration - 1 && tmpRestSteps != 0) {
                if (this.ISTINKER9) {
                    tmpCommandList[3] = Integer
                        .toString(tmpStepNumberOfLastIteration);
                } else {
                    tmpCommandList[2] = Integer
                        .toString(tmpStepNumberOfLastIteration);
                }
            }
            tmpPBuilder.command(tmpCommandList);
            try {
                tmpProcess = tmpPBuilder.start();
                try (BufferedReader tmpBR = new BufferedReader(
                        new InputStreamReader(tmpProcess.getInputStream()))) {
                    if (ISLOGDYNAMIC) {
                        BufferedWriter tmpBW = new BufferedWriter(
                                new FileWriter(tmpDynamicLogName));
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
            } catch(IOException | InterruptedException ex) {
                LOGGER.log(Level.SEVERE,
                        "Exception during tinker's dynamic.exe", ex);
            }
            
            tmpArcName = tmpCurrentDir + tmpParticlePair + ".arc";
            tmpCoordRecord = MIPET4JUtil
                    .getCoordinatesFromArcFile(tmpArcName, tmpAtomNumber1, 
                            tmpAtomNumber2);
            tmpCurrentNeighbors = MIPET4JUtil
                    .getNeighborNumbersBruteForce(tmpCoordRecord, 
                              tmpElements1, 
                              tmpElements2,
                              tmpBoxLength,
                              this.CATCH_RADIUS);
            //      Cell index method
            /*tmpNeighborNumber = MIPET4JUtil
                  .getNeighborNumbers(tmpCoordRecord, 
                          tmpElements1, 
                          tmpElements2,
                          tmpBoxLength,
                          this.CATCH_RADIUS); */
            
            // Copy .xyz file after the last simulation step to result directory
            if (i == tmpIteration - 1) {
                tmpArcName = tmpCurrentDir 
                        + tmpParticlePair 
                        + ".arc";
                tmpXyzName = tmpResultPathName
                    + FILESEPARATOR
                    + tmpParticlePair
                    + "_lastCoords.xyz";
                MIPET4JUtil.writeLastPartToXYZ(tmpArcName, 
                        tmpXyzName, tmpStepNumberOfLastIteration);
                
                // Write _lastStepNeighbors.txt file
                tmpTargetName = tmpResultPathName 
                        + FILESEPARATOR
                        + tmpParticlePair 
                        + "_lastStepNeighbors.txt";
                tmpIterationSize = tmpCurrentNeighbors.size();
                tmpNeighborSize = tmpCurrentNeighbors
                        .get(tmpIterationSize - 1).length;
                try (BufferedWriter tmpBW = new BufferedWriter(
                        new FileWriter(tmpTargetName))) {
                    for (int j = 0; j < tmpNeighborSize; j++) {
                        tmpBW.append(String.valueOf(tmpCurrentNeighbors
                                        .get(tmpIterationSize - 1)[j]));
                        tmpBW.append(LINESEPARATOR);
                    }
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during writing _lastStepNeighbors.txt file.",
                            ex);
                }
            }
            
            if (tmpIteration > 1 && i < tmpIteration - 1) {
                try {
                    Files.deleteIfExists(Paths.get(tmpCurrentDir 
                            + tmpParticlePair 
                            + ".xyz"));
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                        "IOException during deleting .xyz file.", ex);
                }
                tmpSourceName = tmpCurrentDir + tmpParticlePair + ".arc";
                tmpTargetName = tmpCurrentDir + tmpParticlePair + ".xyz";
                MIPET4JUtil.writeLastPartToXYZ(tmpSourceName, 
                        tmpTargetName, tmpStep);
                try {
                    Files.deleteIfExists(Paths.get(tmpCurrentDir 
                            + tmpParticlePair 
                            + ".arc"));
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during deleting .arc file.", ex);
                }
            }
            tmpCurrentNeighborsSize = tmpCurrentNeighbors.size();
            
            for (int j = 0; j < tmpCurrentNeighborsSize; j++) {
                tmpNeighborNumbers.add(tmpCurrentNeighbors.get(j).length);
            }
            
        }
        
        tmpResultNeighborNumber = tmpNeighborNumbers.stream()
                .mapToInt(Integer::intValue).toArray();
        return tmpResultNeighborNumber;
    }
    
    // </editor-fold>
}
