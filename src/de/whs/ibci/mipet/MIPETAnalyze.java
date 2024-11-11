/**
 * MIPET - Mesoscopic Interaction Parameter Estimation with Tinker
 * Copyright (C) 2024  Achim Zielesny (achim.zielesny@googlemail.com)
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
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
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
    private final String FILESEPARATOR = System.getProperty("file.separator");
    
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
        String tmpMinFileName;
        String tmpValueCandidate;
        TinkerXYZ tmpTinkerXYZ;

        tmpEnergyList = new ArrayList<>(this.CHUNKSIZE);
        tmpTinkerXYZ = this.TINKERXYZ;
        
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
        Path tmpPath = Paths.get(tmpArcFileName);
        try (BufferedWriter tmpBW = Files.newBufferedWriter(tmpPath, 
                StandardCharsets.UTF_8)) {
            
            for (int i = 0; i < tmpRot1Size; i++) {
                tmpTinkerXYZ.setCoordinateList1(this.ROTDATA1[i]);
                
                for (int j = 0; j < tmpRot2Size; j++) {
                    tmpIs2Close = MIPETUTIL.isTooClose(
                            this.ROTDATA1[i], 
                            this.ROTDATA2[j], 
                            this.MINATOMDISTANCE);
                    if (!tmpIs2Close) {
                        tmpTinkerXYZ.setHeader(this.PARTICLE_PAIR);
                        tmpTinkerXYZ.setCoordinateList2(this.ROTDATA2[j]);
                        // call getDistance
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
                    "IOException during making .arc file.", ex);
        }
        
        // Start analyze.exe
        //  read .arc files and find intermolecular energy
        ProcessBuilder tmpPBuilder;
        Process tmpProcess;
        int tmpConfigIndex;
        int tmpMinIndex;
        double tmpMinEnergy;
        double tmpValue;
        String tmpSearch;
        String tmpLine;
        
        tmpPBuilder = new ProcessBuilder();
        tmpPBuilder.redirectErrorStream(true);
        tmpPBuilder.command(this.COMMAND_LIST);
        tmpSearch = "Intermolecular Energy";
        tmpMinEnergy = 1E10;
        tmpMinIndex = -1;
        tmpConfigIndex = -1;
        
        try {
            tmpProcess = tmpPBuilder.start();
            try (InputStream tmpInStream = tmpProcess.getInputStream();
                    BufferedReader tmpBR = new BufferedReader(
                            new InputStreamReader(tmpInStream))){

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
                
                tmpProcess.waitFor();
                tmpProcess.destroy();
            }
        } catch(IOException | InterruptedException ex) {
            LOGGER.log(Level.SEVERE,
                    "Exception during tinker's analyze.exe", ex);
        }
            
        // Export .xyz file with lowest intermolecular energy
        int tmpStartIndex;
        int tmpEndIndex;
        StringBuilder tmpPartArc;
        
        tmpStartIndex = tmpMinIndex * (this.ATOMNUMBER + 1);
        tmpEndIndex = tmpStartIndex + this.ATOMNUMBER;
        tmpArcFileName = this.SCRATCH_DIR
                + this.FILESEPARATOR
                + this.PARTICLE_PAIR
                + ".arc" 
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

        // Delete .arc files
        try {
            Files.delete(Paths.get(tmpArcFileName));
        } catch(IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during deleting files in scratch directory.", ex);
        }
        
        return tmpEnergyList;
    }
    
    // </editor-fold>
}
