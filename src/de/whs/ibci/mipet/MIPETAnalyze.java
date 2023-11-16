/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package de.whs.ibci.mipet;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * MIPET4JavaAnalyse class
 *   This class helps to process tinker's analyze.exe parallely
 * @author Mirco Daniel
 */
public class MIPETAnalyze implements Callable<String> {
    
    // <editor-fold defaultstate="collapsed" desc="Final Class variables">
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger
            .getLogger(MIPET.class.getName());
    
    /**
     * Index number (zero-based) of current distance
     */
    private final int DISTANCEINDEX;
    
    /**
     * Index number (zero-based) of current chunk
     */
    private final int CHUNKINDEX;
    
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
     * @param aIndex
     *   Index number of current distance.
     * @param aScratchDir
     *   Scratch directory name
     * @param aParticlePair
     *   Particle pair name
     * @param aCommandList 
     *   Command list for tinker's analyze.exe
     */
    public MIPETAnalyze(
            int aDistanceIndex,
            int aChunkIndex,
            String aScratchDir,
            String aParticlePair,
            String[] aCommandList) {
        this.DISTANCEINDEX = aDistanceIndex;
        this.CHUNKINDEX = aChunkIndex;
        this.SCRATCH_DIR = aScratchDir;
        this.PARTICLE_PAIR = aParticlePair;
        this.COMMAND_LIST = aCommandList;
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    
    @Override
    public String call() {
        ProcessBuilder tmpPBuilder = new ProcessBuilder();
        Process tmpProcess;
        tmpPBuilder.command(this.COMMAND_LIST);
        tmpPBuilder.redirectOutput(new File(this.SCRATCH_DIR 
                + "/" + this.PARTICLE_PAIR + ".out" + this.DISTANCEINDEX
                + "_" + this.CHUNKINDEX));
        try {
            tmpProcess = tmpPBuilder.start();
            tmpProcess.waitFor();
            tmpProcess.destroy();
        } catch(IOException | InterruptedException ex) {
            LOGGER.log(Level.SEVERE,
                    "Exception during tinker's analyze.exe", ex);
        }
        return "Job done";
    }
    
    // </editor-fold>
}
