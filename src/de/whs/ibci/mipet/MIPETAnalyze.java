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
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * MIPET4JavaAnalyse class
 *   This class helps to process tinker's analyze.exe parallely.
 * @author Mirco Daniel
 */
public class MIPETAnalyze implements Callable<WgtEnergyRecord> {
    
    // <editor-fold defaultstate="collapsed" desc="Final Class variables">
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger
            .getLogger(MIPETAnalyze.class.getName());
    
    /**
     * File separator
     */
    private final String FILESEPARATOR = FileSystems.getDefault()
            .getSeparator();

    /**
     * Tinkerxyz object
     */
    private final TinkerXYZ TINKERXYZ;
    
    /**
     * Whether Tinker is used (true) or not (false)
     */
    private final boolean ISTINKERON;
    
    /**
     * Whether Fraction is used (true) or not (false, all data will be used)
     */
    private final boolean ISFRACTIONONE;
    
    /**
     * Minimum of atom to atom distance
     */
    private final double MINATOMDISTANCE;
    
    /**
     * Flattened coordinates of molecule 1
     */
    private final Flat3DArray FLAT_ROT1;
    
    /**
     * Flattened coordinates of molecule 2
     */
    private final Flat3DArray FLAT_ROT2;
    
    /**
     * Number of Configurations (atoms) from molecule 1 and 2
     */
    private final int NUM_ROT1, NUM_ATOMS1;
    private final int NUM_ROT2, NUM_ATOMS2;
    
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
     * MoleculeRecords
     */
    private final MoleculeRecord MOLECULE_1;
    private final MoleculeRecord MOLECULE_2;
    
    /**
     * Temperature in K
     */
    private final double TEMP;

    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    
    /**
     * Constructor of MIPET4JavaAnalyze class
     * 
     * @param aTinkerXYZ: Tinker xyz-data
     * @param aTinkerOn: Flag for whether tinker is used or not
     * @param aDistanceIndex: Distance index
     * @param aChunkIndex: Chunk index
     * @param anAtomNumber: Atom number
     * @param aScratchDir: Scratch directory name
     * @param aFlatRot1: Flattened coordinates of molecule 1
     * @param aFlatRot2: Flattened coordinates of molecule 2
     * @param aMinAtomDistance: Minimum atom to atom distance (if smaller it will not be calculated)
     * @param aCommandList: Command list for tinker's analyze.exe
     * @param aMolRec1: Molecule datas of molecule1
     * @param aMolRec2: Molecule datas of molecule2
     * @param aIsFractionOne: Whether Boltzmann fraction = 1.0 (true) or not (false)
     * @param aTemperature: Temperature [K]
     */
    public MIPETAnalyze(
            TinkerXYZ aTinkerXYZ,
            boolean aTinkerOn,
            int aDistanceIndex,
            int aChunkIndex,
            int anAtomNumber,
            double aMinAtomDistance,
            Flat3DArray aFlatRot1,
            Flat3DArray aFlatRot2,
            String aScratchDir,
            String[] aCommandList,
            MoleculeRecord aMolRec1,
            MoleculeRecord aMolRec2,
            boolean aIsFractionOne,
            double aTemperature) {
        this.TINKERXYZ = aTinkerXYZ;
        this.ISTINKERON = aTinkerOn;
        this.DISTANCEINDEX = aDistanceIndex;
        this.CHUNKINDEX = aChunkIndex;
        this.ATOMNUMBER = anAtomNumber;
        this.MINATOMDISTANCE = aMinAtomDistance;
        this.FLAT_ROT1 = aFlatRot1;
        this.FLAT_ROT2 = aFlatRot2;
        this.SCRATCH_DIR = aScratchDir;
        this.COMMAND_LIST = aCommandList;
        this.MOLECULE_1 = aMolRec1;
        this.MOLECULE_2 = aMolRec2;
        this.ISFRACTIONONE = aIsFractionOne;
        this.TEMP = aTemperature;
        this.NUM_ATOMS1 = aTinkerXYZ.getN_atom1();
        this.NUM_ATOMS2 = aTinkerXYZ.getN_atom2();
        this.NUM_ROT1 = aFlatRot1.dim1;
        this.NUM_ROT2 = aFlatRot2.dim1;
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    
    @Override
    public WgtEnergyRecord call() {
        /*
         * Gas constant R [kcal/(mol*K)]
         */
        final double GASCONST = 1.98720425864E-3;

        /*
         * Conversion factor Joule to Calorie (therm.)
         */
        final double j_CAL = 0.2390057361;

        /*
         * Coulomb constant in N*m^2/C^2
         */
        final double COULOMB = 8.9875517862E9;

        /*
         * Elementary charge in C
         */
        final double ELEMENTCHARGE = 1.602176634E-19;

        /*
         * Avogadro constant in 1/mole
         */
        double AVOGADRO = 6.02214076E23;
        
        /*
        * Inverse of R*T
        */
        final double INV_RT = 1 / (GASCONST * this.TEMP);
        
        boolean tmpIs2Close;
        int tmpChunkSize;
        int tmpChunkIndex;
        int tmpConfigIndex;
        int tmpMinIndex;
        int energyCount = 0; // Counter for the valid calculations of energy
        double tmpChargeQ; // in C^2
        double tmpFactor; // in kcal * m / (mole * C^2)
        double[] tmpEpsilons1;
        double[] tmpEpsilons2;
        double[] tmpSigmas1;
        double[] tmpSigmas2;
        double[] tmpCharges1;
        double[] tmpCharges2;
        double tmpValue;
        double localMinEnergy; // in kcal/mole
        double tmpWgt;
        double tmpWgtxE;
        double tmpSumWgt = 0;
        double tmpSumWgtxE = 0;
        
        double[] energyArray;
        String tmpForcefield;
        String tmpParticlePair;
        String tmpParticle1;
        String tmpParticle2;
        String tmpArcFileName;
        String tmpMinFileName;
        String tmpValueCandidate;
        TinkerXYZ tinkerXYZ;
        TinkerXYZ tmpTinkerXYZMin;

        tmpChargeQ = ELEMENTCHARGE * ELEMENTCHARGE;
        tmpFactor = AVOGADRO * j_CAL * tmpChargeQ * COULOMB * 1E7;
        tinkerXYZ = this.TINKERXYZ;
        tmpTinkerXYZMin = new TinkerXYZ();
        tmpForcefield = tinkerXYZ.getForcefieldName();
        tmpParticle1 = this.TINKERXYZ.getParticleName1();
        tmpParticle2 = this.TINKERXYZ.getParticleName2();
        tmpParticlePair =  tmpParticle1 + "_" + tmpParticle2;
        localMinEnergy = Double.MAX_VALUE;
        tmpMinIndex = -1;
        tmpConfigIndex = -1;
        tmpSigmas1 = null;
        tmpSigmas2 = null;
        tmpEpsilons1 = null;
        tmpEpsilons2 = null;
        tmpCharges1 = null;
        tmpCharges2 = null;
        
        if (tmpForcefield.equals("OPLSAALIGPARGEN")) {
            tmpSigmas1 = MOLECULE_1.sigmas();
            tmpSigmas2 = MOLECULE_2.sigmas();
            tmpEpsilons1 = MOLECULE_1.epsilons();
            tmpEpsilons2 = MOLECULE_2.epsilons();
            tmpCharges1 = MOLECULE_1.charges();
            tmpCharges2 = MOLECULE_2.charges();
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
        tmpChunkSize = NUM_ROT1 * NUM_ROT2;
        energyArray = new double[tmpChunkSize];
        tmpChunkIndex = 0;
        tinkerXYZ.setHeader(tmpParticlePair, ISTINKERON);
        if (!tmpForcefield.equals("OPLSAALIGPARGEN") || ISTINKERON) {
            // <editor-fold defaultstate="collapsed" desc="Use Tinker">
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
                
                for (int i = 0; i < NUM_ROT1; i++) {
                    double[][] currentCoords1 = extractSingleConf(FLAT_ROT1, i, 
                            NUM_ATOMS1);
                    tinkerXYZ.setCoordinateList1(currentCoords1, ISTINKERON);

                    for (int j = 0; j < NUM_ROT2; j++) {
                        tmpIs2Close = isTooCloseFlat(FLAT_ROT1, i, FLAT_ROT2, j, 
                                NUM_ATOMS1, NUM_ATOMS2, MINATOMDISTANCE);

                        if (!tmpIs2Close) {
                            double[][] currentCoords2 = 
                                    extractSingleConf(FLAT_ROT2, j, NUM_ATOMS2);
                            tinkerXYZ.setCoordinateList2(currentCoords2, 
                                    ISTINKERON);
                            tmpBW.append(tinkerXYZ.getFileContent());
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
            if (NUM_ATOMS1 == 1 && NUM_ATOMS2 == 1) {
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
                                tmpValue = Double
                                       .parseDouble(tmpValueCandidate);
                                if (ISFRACTIONONE) {
                                    tmpWgt = Math.exp(-tmpValue * INV_RT);
                                    tmpWgtxE = tmpWgt * tmpValue;
                                    tmpSumWgt += tmpWgt;
                                    tmpSumWgtxE += tmpWgtxE;
                                } else {
                                    energyArray[energyCount] = tmpValue;
                                    energyCount++;
                                }
                                tmpConfigIndex++;
                                if (tmpValue < localMinEnergy) {
                                    localMinEnergy = tmpValue;
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
            tmpPartArc = MIPETUtility.readPartArcFile(tmpArcFileName, 
                    tmpStartIndex, tmpEndIndex);
            TINKERXYZ.writeToXyzFile(tmpMinFileName, tmpPartArc);

            // Delete .arc files
            try {
                Files.delete(Paths.get(tmpArcFileName));
            } catch(IOException ex) {
                LOGGER.log(Level.SEVERE, 
                        "IOException during deleting files in scratch directory.", ex);
            }
            
            // </editor-fold>
        } else {
            // <editor-fold defaultstate="collapsed" desc="Don't use Tinker">
            // Pre-computation of sigma and epsilon
            double[] tmpSigma3_1 = new double[NUM_ATOMS1];
            double[] tmpSqrtEps1 = new double[NUM_ATOMS1];
            
            for (int i = 0; i < NUM_ATOMS1; i++) {
                tmpSigma3_1[i] = tmpSigmas1[i] * tmpSigmas1[i] *tmpSigmas1[i];
                tmpSqrtEps1[i] = Math.sqrt(tmpEpsilons1[i]);
            }
            
            double[] tmpSigma3_2 = new double[NUM_ATOMS2];
            double[] tmpSqrtEps2 = new double[NUM_ATOMS2];
            
            for (int i = 0; i < NUM_ATOMS2; i++) {
                tmpSigma3_2[i] = tmpSigmas2[i] * tmpSigmas2[i] *tmpSigmas2[i];
                tmpSqrtEps2[i] = Math.sqrt(tmpEpsilons2[i]);
            }
            
            int localBestRot1 = -1;
            int localBestRot2 = -1;

            for (int i = 0; i < NUM_ROT1; i++) {
                
                for (int j = 0; j < NUM_ROT2; j++) {

                    // Fast distance check
                    if (isTooCloseFlat(this.FLAT_ROT1, i, this.FLAT_ROT2, j, 
                            NUM_ATOMS1, NUM_ATOMS2, MINATOMDISTANCE)) {
                        continue;
                    }

                    double energyValue = 0.0;

                    // Main calculation: double loop over the atoms
                    for (int k = 0; k < NUM_ATOMS2; k++) {
                        // Coordinates of the molecule 1 (dyretly from the flattend array)
                        double ax = this.FLAT_ROT1.get(i, k, 0);
                        double ay = this.FLAT_ROT1.get(i, k, 1);
                        double az = this.FLAT_ROT1.get(i, k, 2);

                        double k_Sigma3 = tmpSigma3_1[k];
                        double k_SqrtEps = tmpSqrtEps1[k];
                        double k_CoulombBase = tmpFactor * tmpCharges1[k];

                        for (int m = 0; m < NUM_ATOMS2; m++) {
                            // Distance calculation
                            double dx = ax - this.FLAT_ROT2.get(j, m, 0);
                            double dy = ay - this.FLAT_ROT2.get(j, m, 1);
                            double dz = az - this.FLAT_ROT2.get(j, m, 2);

                            double r2 = dx * dx + dy * dy + dz * dz;
                            double invR2 = 1.0 / r2;

                            // Lennard-Jones Optimierung (r6 aus invR2)
                            double invR6 = invR2 * invR2 * invR2;
                            double sigma6 = k_Sigma3 * tmpSigma3_2[m];
                            double ratio6 = sigma6 * invR6;
                            double epsProd = k_SqrtEps * tmpSqrtEps2[m];
                            double lj = 4.0 * epsProd * ratio6 * (ratio6 - 1.0);
                            double coulomb = k_CoulombBase * tmpCharges2[m] 
                                    * Math.sqrt(invR2);
                            energyValue += lj + coulomb;
                        }
                        
                    }

                    if (ISFRACTIONONE) {
                        // Ignore very small values
                        if (energyValue <= 20.0) {
                            double wgt = Math.exp(-energyValue * INV_RT);
                            tmpSumWgt += wgt;
                            tmpSumWgtxE += wgt * energyValue;
                        }
                    } else {
                        energyArray[energyCount] =energyValue;
                        energyCount++;
                    }

                    if (energyValue < localMinEnergy) {
                        localMinEnergy = energyValue;
                        localBestRot1 = i;
                        localBestRot2 = j;
                    }
                }
                
            }
            
            if (localBestRot1 != -1 && localBestRot2 != -1) {
                double[][] bestCoords1 = extractSingleConf(FLAT_ROT1, 
                        localBestRot1, NUM_ATOMS1);
                double[][] bestCoords2 = extractSingleConf(FLAT_ROT2, 
                        localBestRot2, NUM_ATOMS2);
                tinkerXYZ.setCoordinateList1(bestCoords1, ISTINKERON);
                tinkerXYZ.setCoordinateList2(bestCoords2, ISTINKERON);
                tmpTinkerXYZMin = tinkerXYZ.clone(); 
            }
            
            //tmpTinkerXYZMin.makeArcFile(tmpMinFileName);
        }
        
        // </editor-fold>
        if (ISFRACTIONONE) {
            double[] resultList = new double[1];
            resultList[0] = localMinEnergy;
            return new WgtEnergyRecord(resultList, tmpSumWgt, tmpSumWgtxE, 
                    tmpTinkerXYZMin);
        } else {
            Arrays.sort(energyArray, 0, energyCount);
            double[] finalEnergies = Arrays.copyOf(energyArray, energyCount);
            return new WgtEnergyRecord(finalEnergies, 0, 0, tmpTinkerXYZMin);
        }
    }
    // </editor-fold>
            
    // </editor-fold>
    
    /**
     * Checks whether two molecules are too close each other
     * 
     * @param c1 Flatarray of molecule A
     * @param i Configuration index of molecule A
     * @param c2 Flatarray of molecule B
     * @param j Configuration index of molecule B
     * @param size1 Atom number of molecule A
     * @param size2 Atom number of molecule B
     * @param minDist Minimum distance of atom of molecule A and atom of molecule B
     * @return Whether two molecules are too close (true) or not (false)
     */
    private boolean isTooCloseFlat(Flat3DArray c1, int i, 
            Flat3DArray c2, int j, int size1, int size2, double minDist) {
        double minDist2 = minDist * minDist;
        
        for (int k = 0; k < size1; k++) {
            double ax = c1.get(i, k, 0);
            double ay = c1.get(i, k, 1);
            double az = c1.get(i, k, 2);
            
            for (int m = 0; m < size2; m++) {
                double dx = ax - c2.get(j, m, 0);
                double dy = ay - c2.get(j, m, 1);
                double dz = az - c2.get(j, m, 2);
                if ((dx*dx + dy*dy + dz*dz) < minDist2) {
                    return true;
                }
            }
            
        }
        
        return false;
    }
    
    /**
     * Converts flat array to double[][][] array
     * 
     * @param flat Flat array
     * @param numRot Number of conformation
     * @param numAtoms Number of atoms
     * @return Coordinate datas in double[][][] array
     */
    private double[][][] inflate(double[] flat, int numRot, int numAtoms) {
        double[][][] inflated = new double[numRot][numAtoms][3];
        
        for (int r = 0; r < numRot; r++) {
            
            for (int a = 0; a < numAtoms; a++) {
                int base = (r * numAtoms + a) * 3;
                inflated[r][a][0] = flat[base];
                inflated[r][a][1] = flat[base + 1];
                inflated[r][a][2] = flat[base + 2];
            }
            
        }
        
        return inflated;
    }
    
    /**
     * Converts best configuration (in flat array) to double[][]
     * 
     * @param flat Flat array
     * @param bestConfIdx Index of the best configuration
     * @param numAtoms Number of atoms
     * @return Coodinates in double[][]
     */
    private double[][] extractSingleConf(Flat3DArray flat, int bestConfIdx, 
            int numAtoms) {
        double [][] result = new double[numAtoms][3];
        
        for (int a = 0; a < numAtoms; a++) {
            result[a][0] = flat.get(bestConfIdx, a, 0);
            result[a][1] = flat.get(bestConfIdx, a, 1);
            result[a][2] = flat.get(bestConfIdx, a, 2);
        }
    
        return result;
    }
}