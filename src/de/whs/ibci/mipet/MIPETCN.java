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

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.*;
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
    private final String FILESEPARATOR = FileSystems.getDefault()
            .getSeparator();
    
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
    
    /**
     * Constructor of MIPETCN
     * 
     * @param aCommandList: Command list for tinker's dynamic.exe
     * @param aJobTaskRecord: Jobtasks
     * @param aScratchDir: Scratch directory name
     * @param aCatchRadius: Catch radius
     * @param anIsWarmup: Whether warmup will be done
     * @param anIsLogWarmup: Whether log data for warmup will be generated or not
     * @param anIsLogDynamic: Whether log data for dynamic will be generated or not
     * @param anIsTinker9: Whether tinker9 (GPU supported calculation) is used or not
     */
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
        LinkedList<int[]> currNeighbors;
        ProcessBuilder pBuilder;
        CoordinatesRecord coordRecord;
        
        LinkedList<Integer> neighborNumbers = new LinkedList<>();
        MIPETUtility MIPET4JUtil = new MIPETUtility();
        int stepsPerRound = 2000;
        String forcefield = this.JOBTASK_RECORD.forcefield_CN_Name();
        String particle1 = this.JOBTASK_RECORD.particleName1(); // solute
        String particle2 = this.JOBTASK_RECORD.particleName2(); // solvent
        boolean isSameParticle = particle1.equals(particle2);
        String particlePair = particle1 + "_" + particle2;
        String resultPathName = this.JOBTASK_RECORD.result_CN_PathName();
        Path currPath = Paths.get(this.SCRATCH_DIR, forcefield, particlePair);
        TinkerXYZ tXyz1 = new TinkerXYZ(forcefield, particle1,
                currPath.resolve(particle1 + ".xyz"));
        TinkerXYZ tXyz2;
        double boxLength = MIPET4JUtil.getBoxLength(currPath.resolve(
                particlePair + ".xyz"));
        int atomNumber1 = tXyz1.getN_atom();
        int atomNumber2;
        String[] elements1 = tXyz1.getElementList1();
        String[] elements2;
        if (isSameParticle) {
            atomNumber2 = atomNumber1;
            elements2 = elements1.clone();
        } else {
            tXyz2 = new TinkerXYZ(forcefield, particle2,
                    currPath.resolve(particle2 + ".xyz"));
            atomNumber2 = tXyz2.getN_atom();
            elements2 = tXyz2.getElementList1();
        }
        int stepNumber;
        if (this.ISTINKER9) {
            stepNumber = Integer.parseInt(this.COMMAND_LIST[3]);
        } else {
            stepNumber = Integer.parseInt(this.COMMAND_LIST[2]);
        }
        int step;
        int restStep = stepNumber % stepsPerRound;
        int stepOfLastIteration;
        String[] commandList = this.COMMAND_LIST.clone();
        int iteration;
        if (stepNumber > stepsPerRound) {
            iteration = (int)Math.ceil((double)stepNumber 
                    / stepsPerRound);
            step = stepsPerRound;
            if (this.ISTINKER9) {
                commandList[3] = Integer.toString(step);
            } else {
                commandList[2] = Integer.toString(step);
            }
            if ( restStep == 0) {
                stepOfLastIteration = step;
            } else {
                stepOfLastIteration = restStep;
            }
        } else {
            iteration = 1;
            step = stepNumber;
            stepOfLastIteration = step;
        }
        Process process;
        if(ISWARMUP) {
            pBuilder = new ProcessBuilder();
            pBuilder.command(this.COMMAND_LIST);
            pBuilder.redirectErrorStream(true);
            if (ISLOGWARMUP) {
                Path warmUpLogPath = Paths.get(resultPathName,
                        particle1 + "_" + particle2 + "_warmUp.log");
                pBuilder.redirectOutput(warmUpLogPath.toFile());
            } else {
                // This loop is necessary for linux version
                pBuilder.redirectOutput(ProcessBuilder.Redirect.DISCARD);
            }
            try {
                process = pBuilder.start();
                process.waitFor();
            } catch(IOException ex) {
                LOGGER.log(Level.SEVERE,
                        "IOException during tinker's dynamic.exe", ex);
            } catch(InterruptedException ex) {
                LOGGER.log(Level.SEVERE,
                        "InterruptException during tinker's dynamic.exe", ex);
            }
            
            // Rename .arc to .xyz
            Path source = currPath.resolve(particlePair + ".arc");
            Path target = currPath.resolve(particlePair + ".xyz");
            try {
                Files.move(source, target, StandardCopyOption.ATOMIC_MOVE);
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE, 
                        "IOException during remove .arc file to .xyz file.",
                        ex);
            }
            
            // Write _lastStepNeighbors.txt file
            coordRecord = MIPET4JUtil.getCoordinatesFromArcFile(
                    target, atomNumber1, atomNumber2);
            currNeighbors = MIPET4JUtil.getNeighborNumbersBruteForce(
                    coordRecord, 
                    elements1, 
                    elements2,
                    boxLength,
                    this.CATCH_RADIUS);
            Path  targetPath = Paths.get(resultPathName, 
                    particlePair + "_warmUpNeighbors.txt");
            int iterationSize = currNeighbors.size();
            int neighborsSize = currNeighbors.get(iterationSize - 1).length;
            
            try (BufferedWriter writer = Files.newBufferedWriter(targetPath)) {
                
                for (int i = 0; i < neighborsSize; i++) {
                    writer.append(String.valueOf(currNeighbors
                            .get(iterationSize - 1)[i]));
                    writer.newLine();
                }
                
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE, 
                        "IOException during writing _warmUpNeighbors.txt file.",
                        ex);
            }
            
            // Copy .xyz file after warmup to result directory
            source = target;
            target = Paths.get(resultPathName, 
                    particlePair + "_warmUpCoords.xyz");
            try {
                Files.copy(source, target, StandardCopyOption
                        .REPLACE_EXISTING);
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE, 
                        "IOException during copying .xyz file.", ex);
            }
            
            // Delete useless files
            try (Stream<Path> pathList = Files.list(currPath)) {
                pathList.filter(file -> !Files.isDirectory(file))
                        .filter(file -> file.endsWith(".xyz_2"))
                        .map(Path::toFile)
                        .forEach(File::delete);
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE, 
                        "IOException during deleting .xyz_2 file.", ex);
            }
            return null;
        }
        Path dynamicLogPath = Paths.get(resultPathName,
                particle1 + "_" + particle2 + "_dynamic.log");
        
        for (int i = 0; i < iteration; i++) {
            pBuilder = new ProcessBuilder();
            if (i == iteration - 1 && restStep != 0) {
                if (this.ISTINKER9) {
                    commandList[3] = Integer
                        .toString(stepOfLastIteration);
                } else {
                    commandList[2] = Integer
                        .toString(stepOfLastIteration);
                }
            }
            pBuilder.command(commandList);
            pBuilder.redirectErrorStream(true);
            if (ISLOGDYNAMIC) {
                pBuilder.redirectOutput(dynamicLogPath.toFile());
            } else {
                // This loop is necessary for linux version
                pBuilder.redirectOutput(ProcessBuilder.Redirect.DISCARD);
            }
            try {
                process = pBuilder.start();
                process.waitFor();
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE,
                        "IOException during tinker's dynamic.exe", ex);
            } catch ( InterruptedException ex) {
                LOGGER.log(Level.SEVERE,
                        "InterruptedException during tinker's dynamic.exe", ex);
            }
            Path arcPath = currPath.resolve(particlePair + ".arc");
            coordRecord = MIPET4JUtil.getCoordinatesFromArcFile(arcPath, 
                    atomNumber1, atomNumber2);
            currNeighbors = MIPET4JUtil
                    .getNeighborNumbersBruteForce(coordRecord, 
                              elements1, 
                              elements2,
                              boxLength,
                              this.CATCH_RADIUS);
            //      Cell index method
            /*tmpNeighborNumber = MIPET4JUtil
                  .getNeighborNumbers(tmpCoordRecord, 
                          tmpElements1, 
                          tmpElements2,
                          tmpBoxLength,
                          this.CATCH_RADIUS); */
            
            // Copy .xyz file after the last simulation step to result directory
            if (i == iteration - 1) {
                arcPath = currPath.resolve(particlePair + ".arc");
                Path xyzPath = Paths.get(resultPathName, particlePair 
                        + "_lastCoords.xyz");
                MIPET4JUtil.writeLastPartToXYZ(arcPath, 
                        xyzPath, stepOfLastIteration);
                
                // Write _lastStepNeighbors.txt file
                Path targetPath = Paths.get(resultPathName,
                        particlePair + "_lastStepNeighbors.txt");
                int iterationSize = currNeighbors.size();
                int neighborSize = currNeighbors
                        .get(iterationSize - 1).length;
                try (BufferedWriter writer = Files.newBufferedWriter(
                        targetPath)) {
                    
                    for (int j = 0; j < neighborSize; j++) {
                        writer.append(String.valueOf(currNeighbors
                                        .get(iterationSize - 1)[j]));
                        writer.append(LINESEPARATOR);
                    }
                    
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during writing _lastStepNeighbors.txt file.",
                            ex);
                }
            }
            
            if (iteration > 1 && i < iteration - 1) {
                try {
                    Files.deleteIfExists(Paths.get(currPath 
                            + particlePair 
                            + ".xyz"));
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                        "IOException during deleting .xyz file.", ex);
                }
                Path sourcePath = currPath.resolve(particlePair + ".arc");
                Path targetPath = currPath.resolve(particlePair + ".xyz");
                MIPET4JUtil.writeLastPartToXYZ(sourcePath, targetPath, step);
                try {
                    Files.deleteIfExists(Paths.get(currPath 
                            + particlePair 
                            + ".arc"));
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                            "IOException during deleting .arc file.", ex);
                }
            }
            int currentNeighborsSize = currNeighbors.size();
            
            for (int j = 0; j < currentNeighborsSize; j++) {
                neighborNumbers.add(currNeighbors.get(j).length);
            }
            
        }
        
        int[] resultNeighborNumber = neighborNumbers.stream()
                .mapToInt(Integer::intValue).toArray();
        return resultNeighborNumber;
    }
    
    // </editor-fold>
}