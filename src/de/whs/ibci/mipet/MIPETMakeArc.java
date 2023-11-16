/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package de.whs.ibci.mipet;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *  MIPETMakeArc class
    This class helps to process make .arc files parallely
 * @author Mirco Daniel
 */
public class MIPETMakeArc implements Callable<String>{
    // <editor-fold defaultstate="collapsed" desc="Final Class variables">
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger
            .getLogger(MIPET.class.getName());
    private static final MIPETUtility MIPET4JUTIL = 
            new MIPETUtility();
    private final String PARTICLEPAIR;
    private final String PATH;
    private final TinkerXYZ TINKERXYZ1;
    private final TinkerXYZ TINKERXYZ2;
    private final double[][][] ROTDATA1;
    private final double[][][] ROTDATA2;
    private final int CHUNKSIZE;
    private final double MINATOMDISTANCE;
    
    // </editor-fold> 
    
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    /**
     * Constructor MIPET4JavaMakeArc
     * 
     * @param aParticlePair
     *   Particle pair name
     * @param aPath
     *   Path name
     * @param aTinkerXYZ1
     *   TinkerXYZ object of first particle
     * @param aTinkerXYZ2
     *   TinkerXYZ object of second particle
     * @param aRotData1
     *   Rotated coordinates of first particle
     * @param aRotData2
     *   Rotated coordinates of second particle
     * @param aChunkSize
     *   Chunk size
     * @param aMinAtomDistance 
     *   Minium atom to atom distance is taken account
     */
    public MIPETMakeArc(
            String aParticlePair,
            String aPath,
            TinkerXYZ aTinkerXYZ1,
            TinkerXYZ aTinkerXYZ2,
            double[][][] aRotData1,
            double[][][] aRotData2,
            int aChunkSize,
            double aMinAtomDistance) {
        this.PARTICLEPAIR = aParticlePair;
        this.PATH = aPath;
        this.TINKERXYZ1 = aTinkerXYZ1;
        this.TINKERXYZ2 = aTinkerXYZ2;
        this.ROTDATA1 = aRotData1;
        this.ROTDATA2 = aRotData2;
        this.CHUNKSIZE = aChunkSize;
        this.MINATOMDISTANCE = aMinAtomDistance;
    }
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    @Override
    public String call() {
        FileWriter tmpOutputFile = null;
        int tmpRot1Size;
        int tmpRot2Size;
        int tmpChunkIndex;
        boolean tmpIs2Close;
        TinkerXYZ tmpTinkerXYZ;
        
        tmpTinkerXYZ = new TinkerXYZ(this.TINKERXYZ1, this.TINKERXYZ2);

        // Check if the particles are not too close together
        // Save tmpOutput in scratch directory
        tmpRot1Size = this.ROTDATA1.length;
        tmpRot2Size = this.ROTDATA2.length;
        tmpChunkIndex = 0;
        
        try (BufferedWriter tmpBW = new BufferedWriter(
                new FileWriter(this.PATH), 65536)) {
            
            for (int i = 0; i < tmpRot1Size; i++) {
                tmpTinkerXYZ.setCoordinateList1(this.ROTDATA1[i]);
                
                for (int j = 0; j < tmpRot2Size; j++) {
                    tmpIs2Close = MIPET4JUTIL.isTooClose(
                            this.ROTDATA1[i], 
                            this.ROTDATA2[j], 
                            this.MINATOMDISTANCE);
                    if (!tmpIs2Close) {
                        tmpTinkerXYZ.setHeader(this.PARTICLEPAIR);
                        tmpTinkerXYZ.setCoordinateList2(this.ROTDATA2[j]);
                        tmpBW.append(tmpTinkerXYZ.getFileContent());
                    }
                    tmpChunkIndex++;
                    if (tmpChunkIndex >= this.CHUNKSIZE) {
                        break;
                    }
                }
                
            }

        } catch(IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException in MIPET4JavaMakeArc.", ex);
        }
        return "";
    }
    
    // </editor-fold>
}
