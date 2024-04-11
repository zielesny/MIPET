/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package de.whs.ibci.mipet;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * MIPET4JavaAnalyse class
 *   This class helps to process tinker's analyze.exe parallely
 * @author Mirco Daniel
 */
public class MIPETAnalyze implements Callable<ArrayList<Double>> {
    
    // <editor-fold defaultstate="collapsed" desc="Final Class variables">
    /**
     * Instance object of MIPETUtility
     */
    final private static MIPETUtility MIPETUTIL = new MIPETUtility();
    
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger
            .getLogger(MIPET.class.getName());
    
    /**
     * File separator
     */
    private final String FILESEPARATOR = System
            .getProperty("file.separator");
    
    /**
     * Tinkerxyz object
     */
    private final TinkerXYZ TINKERXYZ;
    private final double MINATOMDISTANCE;
    private final double[][][] ROTDATA1;
    private final double[][][] ROTDATA2;
    private final int CHUNKSIZE;
    
    /**
     * Index number (zero-based) of current distance
     */
    private final int DISTANCEINDEX;
    
    /**
     * Index number (zero-based) of current chunk
     */
    private final int CHUNKINDEX;
    
    /**
     * Atom number of both molecules
     */
    private final int ATOMNUMBER;
    
    /**
     * Scratch directory name
     */
    private final String SCRATCH_DIR;
    
    /**
     * Both particle names
     */
    private final String PARTICLE_PAIR;
    
    /**
     * Command list for tinker's analyze.exe
     */
    private final String[] COMMAND_LIST;
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    
    /**
     * Constructor of MIPET4JavaAnalyze class
     * 
     * @param aTinkerXYZ
     * @param aDistanceIndex
     * @param aChunkIndex
     * @param anAtomNumber
     * @param aChunkSize
     * @param aScratchDir
     *   Scratch directory name
     * @param aMinAtomDistance
     * @param aParticlePair
     *   Particle pair name
     * @param aRotData1
     * @param aCommandList 
     *   Command list for tinker's analyze.exe
     * @param aRotData2
     */
    public MIPETAnalyze(
            TinkerXYZ aTinkerXYZ,
            int aDistanceIndex,
            int aChunkIndex,
            int anAtomNumber,
            int aChunkSize,
            double aMinAtomDistance,
            double[][][] aRotData1,
            double[][][] aRotData2,
            String aScratchDir,
            String aParticlePair,
            String[] aCommandList) {
        this.ROTDATA1 = aRotData1;
        this.ROTDATA2 = aRotData2;
        this.CHUNKSIZE = aChunkSize;
        this.MINATOMDISTANCE = aMinAtomDistance;
        this.TINKERXYZ = aTinkerXYZ;
        this.DISTANCEINDEX = aDistanceIndex;
        this.CHUNKINDEX = aChunkIndex;
        this.ATOMNUMBER = anAtomNumber;
        this.SCRATCH_DIR = aScratchDir;
        this.PARTICLE_PAIR = aParticlePair;
        this.COMMAND_LIST = aCommandList;
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    
    @Override
    public ArrayList<Double> call() {
        boolean tmpIs2Close;
        int tmpRot1Size;
        int tmpRot2Size;
        int tmpChunkIndex;
        ArrayList<Double> tmpEnergyList;
        String tmpArcFileName;
        String tmpOutFileName;
        String tmpMinFileName;
        String tmpValueCandidate;
        Path tmpArcPath;

        tmpEnergyList = new ArrayList<>(this.CHUNKSIZE);
        tmpOutFileName = this.SCRATCH_DIR
                + this.FILESEPARATOR
                + this.PARTICLE_PAIR + ".out" 
                + this.DISTANCEINDEX + "_"
                + this.CHUNKINDEX;
        
        // Check if the particles are not too close together
        // Save .arc file in scratch directory
        tmpRot1Size = this.ROTDATA1.length;
        tmpRot2Size = this.ROTDATA2.length;
        tmpChunkIndex = 0;
        tmpArcFileName = this.SCRATCH_DIR
                + FILESEPARATOR 
                + this.PARTICLE_PAIR 
                + ".arc"
                + this.DISTANCEINDEX
                + "_"
                + this.CHUNKINDEX;
        tmpArcPath = Paths.get(tmpArcFileName);
        try (BufferedWriter tmpBW = Files.newBufferedWriter(tmpArcPath, 
                StandardCharsets.UTF_8)) {
            
            for (int i = 0; i < tmpRot1Size; i++) {
                this.TINKERXYZ.setCoordinateList1(this.ROTDATA1[i]);
                
                for (int j = 0; j < tmpRot2Size; j++) {
                    tmpIs2Close = MIPETUTIL.isTooClose(
                            this.ROTDATA1[i], 
                            this.ROTDATA2[j], 
                            this.MINATOMDISTANCE);
                    if (!tmpIs2Close) {
                        this.TINKERXYZ.setHeader(this.PARTICLE_PAIR);
                        this.TINKERXYZ.setCoordinateList2(this.ROTDATA2[j]);
                        tmpBW.append(this.TINKERXYZ.getFileContent());
                    }
                    tmpChunkIndex++;
                    if (tmpChunkIndex >= this.CHUNKSIZE) {
                        break;
                    }
                }
                
            }
            
        } catch(IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during making .arc file.", ex);
        }
        
        // Start analyze
        StringBuilder tmpPartArc;
        ProcessBuilder tmpPBuilder;
        Process tmpProcess;
        int tmpConfigIndex;
        int tmpMinIndex;
        double tmpMinEnergy;
        double tmpValue;
        String tmpSearch;
        String tmpLine;
        Path tmpOutPath;
        
        tmpPBuilder = new ProcessBuilder();
        tmpPBuilder.command(this.COMMAND_LIST);
        tmpPBuilder.redirectOutput(new File(tmpOutFileName));
        tmpSearch = "Intermolecular Energy";
        tmpMinEnergy = 1E10;
        tmpMinIndex = -1;
        tmpConfigIndex = -1;
        
        try {
            tmpProcess = tmpPBuilder.start();
            tmpProcess.waitFor();
            tmpProcess.destroy();
        } catch(IOException | InterruptedException ex) {
            LOGGER.log(Level.SEVERE,
                    "Exception during tinker's analyze.exe", ex);
        }
        
        // Read .arc files and find intermolecular energy
        tmpOutPath = Paths.get(tmpOutFileName);
        try (BufferedReader tmpBR = Files.newBufferedReader(tmpOutPath)) {
            while ((tmpLine = tmpBR.readLine()) != null) {
                if (tmpLine.contains(tmpSearch)) {
                    tmpValueCandidate = tmpLine.substring(25, 50);
                    if (!tmpValueCandidate.contains("D")) {
                        tmpValue = Double.parseDouble(tmpValueCandidate);
                        tmpEnergyList.add(tmpValue);
                        tmpConfigIndex++;
                        if (tmpValue < tmpMinEnergy) {
                            tmpMinEnergy = tmpValue;
                            tmpMinIndex = tmpConfigIndex;
                        }
                    }
                }
            }
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE,
                    "IOException during reading .out file." , ex);
        }
            
        // Export .xyz file with lowest intermolecular energy
        int tmpStartIndex;
        int tmpEndIndex;
        
        tmpStartIndex = tmpMinIndex * (this.ATOMNUMBER + 1);
        tmpEndIndex = tmpStartIndex + this.ATOMNUMBER;
        tmpArcFileName = this.SCRATCH_DIR
                + this.FILESEPARATOR
                + this.PARTICLE_PAIR
                + ".arc" 
                + this.DISTANCEINDEX 
                + "_" 
                + this.CHUNKINDEX;
        tmpOutFileName = this.SCRATCH_DIR
                + this.FILESEPARATOR
                + this.PARTICLE_PAIR
                + ".out" 
                + this.DISTANCEINDEX 
                + "_" 
                + this.CHUNKINDEX;
        tmpPartArc = MIPETUTIL.readPartArcFile(tmpArcFileName, 
                tmpStartIndex, tmpEndIndex);
        tmpMinFileName = this.SCRATCH_DIR
                + this.FILESEPARATOR
                + this.PARTICLE_PAIR
                +"_"
                + this.DISTANCEINDEX 
                + "_" 
                + this.CHUNKINDEX
                + ".0";
        TINKERXYZ.writeToXyzFile(tmpMinFileName, tmpPartArc);
        tmpEnergyList.sort(Comparator.naturalOrder());

        // Delete .arc file and .out file
        try {
            Files.delete(Paths.get(tmpArcFileName));
            Files.delete(Paths.get(tmpOutFileName));
        } catch(IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during deleting files in scratch directory.", ex);
        }
        
        return tmpEnergyList;
    }
    
    // </editor-fold>
}
