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
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * MIPET4JavaAnalyse class
 *   This class helps to process tinker's analyze.exe parallely.
 * @author Mirco Daniel
 */
public class MIPETAnalyze implements Callable<ArrayList<Double>> {
    
    // <editor-fold defaultstate="collapsed" desc="Final Class variables">
    /**
     * Instance object of MIPETUtility
     */
    private static final MIPETUtility MIPETUTIL = new MIPETUtility();
    
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
     * Conversion factor Joule to Calorie (therm.)
     */
    private final double J_CAL = 0.2390057361; 
    
    /**
     * Coulomb constant in N*m^2/C^2
     */
    private final double COULOMB = 8.9875517862E9; 
    
    /**
     * Elementary charge in C
     */
    private final double ELEMENTCHARGE = 1.602176634E-19;
    
    /**
     * Avogadro constant in 1/mole
     */
    private final double AVOGADRO = 6.02214076E23;
    
    /**
     * Tinkerxyz object
     */
    private final TinkerXYZ TINKERXYZ;
    
    private final boolean ISTINKERON;
    private final double MINATOMDISTANCE;
    private final double[][][] ROTDATA1;
    private final double[][][] ROTDATA2;
    
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
     * Command list for tinker's analyze.exe
     */
    private final String[] COMMAND_LIST;
    
    /**
     * MoleculeRecord
     */
    private final LinkedList<MoleculeRecord> MOLECULES;
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    
    /**
     * Constructor of MIPET4JavaAnalyze class
     * 
     * @param aTinkerXYZ
     * @param aTinkerOn Flag for whether tinker is used or not
     * @param aDistanceIndex
     * @param aChunkIndex
     * @param anAtomNumber
     * @param aScratchDir Scratch directory name
     * @param aMinAtomDistance
     * @param aRotData1 Coordinates of first particle
     * @param aCommandList Command list for tinker's analyze.exe
     * @param aRotData2 Coordinatdes of sendond particle
     * @param aMolecules Molecule datas
     */
    public MIPETAnalyze(
            TinkerXYZ aTinkerXYZ,
            boolean aTinkerOn,
            int aDistanceIndex,
            int aChunkIndex,
            int anAtomNumber,
            double aMinAtomDistance,
            double[][][] aRotData1,
            double[][][] aRotData2,
            String aScratchDir,
            String[] aCommandList,
            LinkedList<MoleculeRecord> aMolecules) {
        this.TINKERXYZ = aTinkerXYZ;
        this.ISTINKERON = aTinkerOn;
        this.DISTANCEINDEX = aDistanceIndex;
        this.CHUNKINDEX = aChunkIndex;
        this.ATOMNUMBER = anAtomNumber;
        this.MINATOMDISTANCE = aMinAtomDistance;
        this.ROTDATA1 = aRotData1;
        this.ROTDATA2 = aRotData2;
        this.SCRATCH_DIR = aScratchDir;
        this.COMMAND_LIST = aCommandList;
        this.MOLECULES = aMolecules;
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    
    @Override
    public ArrayList<Double> call() {
        boolean tmpIs2Close;
        int tmpRot1Size;
        int tmpRot2Size;
        int tmpChunkSize;
        int tmpChunkIndex;
        int tmpAtomSize1;
        int tmpAtomSize2;
        int tmpMoleculeSize;
        int tmpID1;
        int tmpID2;
        int tmpConfigIndex;
        int tmpMinIndex;
        double tmpChargeQ; // in C^2
        double tmpFactor; // in kcal * m / (mole * C^2)
        double tmpSigma; // in Angstrom
        double tmpEpsilon; // in kcal/mole
        double tmpDistances[][]; // in Angstrom
        double tmpEpsilons1[];
        double tmpEpsilons2[];
        double tmpSigmas1[];
        double tmpSigmas2[];
        double tmpCharges1[];
        double tmpCharges2[];
        double tmpCoulombEnergy; // in kcal/mole
        double tmpLJEnergy; // in kcal/mole
        double tmpRatio;
        double tmpRatio2; // tmpRatio^2
        double tmpRatio6; // tmpRatio^6
        double tmpValue;
        double tmpMinEnergy; // in kcal/mole
        ArrayList<Double> tmpEnergyList;
        String tmpForcefield;
        String tmpParticlePair;
        String tmpParticle1;
        String tmpParticle2;
        String tmpArcFileName;
        String tmpMinFileName;
        String tmpValueCandidate;
        TinkerXYZ tmpTinkerXYZ;
        TinkerXYZ tmpTinkerXYZMin;
        
        tmpChargeQ = ELEMENTCHARGE * ELEMENTCHARGE;
        tmpFactor = AVOGADRO * J_CAL * tmpChargeQ * COULOMB * 1E7; 
        tmpTinkerXYZ = this.TINKERXYZ;
        tmpTinkerXYZMin = new TinkerXYZ();
        tmpForcefield = tmpTinkerXYZ.getForcefieldName();
        tmpParticle1 = this.TINKERXYZ.getParticleName1();
        tmpParticle2 = this.TINKERXYZ.getParticleName2();
        tmpParticlePair =  tmpParticle1 + "_" + tmpParticle2;
        tmpAtomSize1 = tmpTinkerXYZ.getAtomSize1();
        tmpAtomSize2 = tmpTinkerXYZ.getAtomSize2();
        tmpID1 = 0;
        tmpID2 = 0;
        tmpMinEnergy = 1E10;
        tmpMinIndex = -1;
        tmpConfigIndex = -1;
        tmpSigmas1 = null;
        tmpSigmas2 = null;
        tmpEpsilons1 = null;
        tmpEpsilons2 = null;
        tmpCharges1 = null;
        tmpCharges2 = null;
        
        if (tmpForcefield.equals("OPLSAALIGPARGEN")) {
            tmpMoleculeSize = MOLECULES.size();
            
            for (int i = 0; i < tmpMoleculeSize; i++) {
                if (tmpParticle1.equals(MOLECULES.get(i).name())) {
                    tmpID1 = i;
                    break;
                }
            }
        
            for (int i = 0; i < tmpMoleculeSize; i++) {
                if (tmpParticle2.equals(MOLECULES.get(i).name())) {
                    tmpID2 = i;
                    break;
                }
            }
            
            tmpSigmas1 = MOLECULES.get(tmpID1).sigmas();
            tmpSigmas2 = MOLECULES.get(tmpID2).sigmas();
            tmpEpsilons1 = MOLECULES.get(tmpID1).epsilons();
            tmpEpsilons2 = MOLECULES.get(tmpID2).epsilons();
            tmpCharges1 = MOLECULES.get(tmpID1).charges();
            tmpCharges2 = MOLECULES.get(tmpID2).charges();
        } 
        tmpMinFileName = this.SCRATCH_DIR
                    + this.FILESEPARATOR
                    + tmpParticlePair
                    +"_"
                    + this.DISTANCEINDEX 
                    + "_" 
                    + this.CHUNKINDEX
                    + ".0";
        
        // Check if the particles are not too close together
        // Save .arc file in scratch directory
        tmpRot1Size = this.ROTDATA1.length;
        tmpRot2Size = this.ROTDATA2.length;
        tmpChunkSize = tmpRot1Size * tmpRot2Size;
        tmpEnergyList = new ArrayList<Double>(tmpChunkSize);
        tmpChunkIndex = 0;
        tmpTinkerXYZ.setHeader(tmpParticlePair, ISTINKERON);
        if (!tmpForcefield.equals("OPLSAALIGPARGEN") || ISTINKERON) {
            tmpArcFileName = this.SCRATCH_DIR
                    + FILESEPARATOR 
                    + tmpParticlePair 
                    + ".arc"
                    + this.DISTANCEINDEX
                    + "_"
                    + this.CHUNKINDEX;
            Path tmpPath = Paths.get(tmpArcFileName);
            try (BufferedWriter tmpBW = Files.newBufferedWriter(tmpPath, 
                    StandardCharsets.UTF_8)) {

                for (int i = 0; i < tmpRot1Size; i++) {
                    tmpTinkerXYZ.setCoordinateList1(this.ROTDATA1[i], 
                            ISTINKERON);

                    for (int j = 0; j < tmpRot2Size; j++) {
                        tmpIs2Close = MIPETUTIL.isTooClose(
                                this.ROTDATA1[i], 
                                this.ROTDATA2[j], 
                                this.MINATOMDISTANCE);
                        if (!tmpIs2Close) {
                            tmpTinkerXYZ.setCoordinateList2(this.ROTDATA2[j], 
                                    ISTINKERON);
                            tmpBW.append(tmpTinkerXYZ.getFileContent());
                        } 
                        tmpChunkIndex++;
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
            String tmpSearch;
            String tmpLine;

            tmpPBuilder = new ProcessBuilder();
            tmpPBuilder.redirectErrorStream(true);
            tmpPBuilder.command(this.COMMAND_LIST);
            if (tmpAtomSize1 == 1 && tmpAtomSize2 == 1) {
                tmpSearch = "Total Potential Energy";
            } else {
                tmpSearch = "Intermolecular Energy";
            }
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
                    + tmpParticlePair
                    + ".arc" 
                    + this.DISTANCEINDEX 
                    + "_" 
                    + this.CHUNKINDEX;
            tmpPartArc = MIPETUTIL.readPartArcFile(tmpArcFileName, 
                    tmpStartIndex, tmpEndIndex);
            TINKERXYZ.writeToXyzFile(tmpMinFileName, tmpPartArc);

            // Delete .arc files
            try {
                Files.delete(Paths.get(tmpArcFileName));
            } catch(IOException ex) {
                LOGGER.log(Level.SEVERE, 
                        "IOException during deleting files in scratch directory.", ex);
            }
        } else {
            for (int i = 0; i < tmpRot1Size; i++) {
                tmpTinkerXYZ.setCoordinateList1(this.ROTDATA1[i], 
                        ISTINKERON);

                for (int j = 0; j < tmpRot2Size; j++) {
                    tmpIs2Close = MIPETUTIL.isTooClose(
                            this.ROTDATA1[i], 
                            this.ROTDATA2[j], 
                            this.MINATOMDISTANCE);
                    if (!tmpIs2Close) {
                        tmpTinkerXYZ.setCoordinateList2(this.ROTDATA2[j], 
                                ISTINKERON);
                        tmpTinkerXYZ.setDistances();
                        tmpDistances = tmpTinkerXYZ.getDistances();
                        tmpValue = 0;
                        
                        for (int k = 0; k < tmpAtomSize1; k++) {
                            
                            for (int m = 0; m < tmpAtomSize2; m++) {
                                tmpSigma = Math.sqrt(tmpSigmas1[k] 
                                        * tmpSigmas2[m]);
                                tmpEpsilon = Math.sqrt(tmpEpsilons1[k] 
                                        * tmpEpsilons2[m]);
                                tmpRatio =  tmpSigma
                                        / tmpDistances[k][m];
                                tmpRatio2 = tmpRatio * tmpRatio;
                                tmpRatio6 = tmpRatio2 * tmpRatio2 * tmpRatio2;
                                tmpLJEnergy = 4 * tmpEpsilon * tmpRatio6 
                                        * (tmpRatio6 - 1);
                                tmpCoulombEnergy = tmpFactor 
                                        * tmpCharges1[k] 
                                        * tmpCharges2[m] 
                                        / tmpDistances[k][m];
                                tmpValue += tmpCoulombEnergy + tmpLJEnergy;
                            }
                                
                        }
                        
                        tmpEnergyList.add(tmpValue);
                        if (tmpValue < tmpMinEnergy) {
                            tmpMinEnergy = tmpValue;
                            tmpTinkerXYZMin = tmpTinkerXYZ.clone(); 
                        }
                    }
                    tmpChunkIndex++;
                }

            }
            tmpTinkerXYZMin.makeArcFile(tmpMinFileName);
            
        }
        tmpEnergyList.sort(Comparator.naturalOrder());
        return tmpEnergyList;
    }
    
    // </editor-fold>
}
