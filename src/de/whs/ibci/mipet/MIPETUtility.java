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
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.periodictable.PeriodicTable;
import org.openscience.cdk.geometry.volume.VABCVolume;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;

/**
 * Class MIPETUtility
 * @author Mirco Daniel
 */
public class MIPETUtility{
    
    //<editor-fold defaultstate="collapsed" desc="Final class variables">
    
    /**
     * Resource bundle name
     */
    private static final String BUNDLE_NAME_INTERN = "de.whs.ibci.mipet.MIPET";
    private static final String BUNDLE_NAME_EXTERN = "MIPET.properties";
    
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
    private static final String FILESEPARATOR = FileSystems
            .getDefault().getSeparator();
    
    /**
     * Cache for calculated VdW volumes (SMILES -> Volume)
     */
    private final Map<String, Double> vdwVolumeCache = new HashMap<>();

    /**
     * Cache for calculated atomic masses (SMILES -> Mass)
     */
    private final Map<String, Double> atomicMassCache = new HashMap<>();
    
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Private class variables">
    
    /**
     * Resource bundle
     */
    private ResourceBundle RESOURCE_BUNDLE;
    
    /**
     * Smiles parser
     */
    private SmilesParser smilesParser;
    
    /**
     * Van der Waals radius of elements
     */
    private double[] vdWRadii = new double[93];
    
    /**
     * Atomic numbers of elements
     */
    private HashMap<String, Integer> atomicNumber;
    
    //</editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    
    /** 
     * Constructor
     */
    public MIPETUtility(){
        initialize();
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    
    /**
     * Return resource string for key
     *
     * @param aKey Key
     * @return Resource string for key
     */
    public String getResourceString(String aKey) {
        try {
            return RESOURCE_BUNDLE.getString(aKey).trim();
        } catch (MissingResourceException e) {
            return "Key '" + aKey + "' not found.";
        }
    }
    
    /**
     * Get method for Van der Waals volume by using VABCVolume of CDK
     * 
     * @param smilesString
     *   input fragment as SMILES 
     * @return 
     *   Van der Waals volume by using VAB Volume of CDK
     */
    public double getVdwVolume(String smilesString) {
        if (smilesString == null || smilesString.isEmpty()) {
            throw new IllegalArgumentException("Smiles was null or empty.");
        }
        // Check Cache
        if (vdwVolumeCache.containsKey(smilesString)) {
            return vdwVolumeCache.get(smilesString);
        }
        
        double vabcVolume = 0.0;
        try {
            if (!smilesString.equals("[Na+]")) {
                IAtomContainer particle =smilesParser.parseSmiles(smilesString);
                AtomContainerManipulator.
                        percieveAtomTypesAndConfigureAtoms(particle);
                vabcVolume = VABCVolume.calculate(particle);
                
                // Save to Cache
                vdwVolumeCache.put(smilesString, vabcVolume);
            }
        } catch (CDKException ex) {
            LOGGER.log(Level.SEVERE,"CDKException was thrown.", ex);
        }
        return vabcVolume;
    }
    
    /**
     * Get method for calculation of atomic mass by using CDK 
     * @param aSmilesString
     *   input fragment as SMILES
     * @param anIsSmiles
     *   Flag if aSmilesString is a SMILES or not (element)
     * @return
     *   Atomic mass of the fragment
     */
    public double getAtomicMass(String aSmilesString, boolean anIsSmiles) {
        if (aSmilesString == null || aSmilesString.trim().isEmpty()) {
            return 0;
        }
        if (aSmilesString.trim().equalsIgnoreCase("lp")) {
            return 0;
        }
        String cacheKey = (anIsSmiles ? "S_" : "E_") + aSmilesString.trim();
        
        // Check cache
        if (atomicMassCache.containsKey(cacheKey)) {
            return atomicMassCache.get(cacheKey);
        }
        
        double atomicMass = 0;
        try {
            IAtomContainer particle;
            if (anIsSmiles) {
                particle = smilesParser.parseSmiles(aSmilesString);
            } else {
                particle = smilesParser.parseSmiles("[" + aSmilesString + "]");
            }
            AtomContainerManipulator
                    .percieveAtomTypesAndConfigureAtoms(particle);
            atomicMass  = AtomContainerManipulator.getMass(particle);
            
            // Save to cach
            atomicMassCache.put(cacheKey, atomicMass);
        } catch (CDKException ex) {
            LOGGER.log(Level.SEVERE, 
                    "CDKException was thrown in getAtomicMass.", ex);
        }
        return atomicMass;
    }
    
    /**
     * Get method for atomic number table by using CDK
     * @return 
     * atomic number table (key: chemical symbol, value: atomic number)
     */
    public HashMap<String, Integer> getAtomicNumberTable() {
        
        HashMap<String, Integer> tmpAtomicNumberTable = new HashMap<>();
        
        for (int i = 1; i < 93; i++) {
            tmpAtomicNumberTable.put(PeriodicTable.getSymbol(i), i); 
        }
        
        return tmpAtomicNumberTable;
    }
    
    /**
     * Get method for vdW-radii of atoms by using CDK
     *   (remarks: there are some atoms without data)
     * 
     * @return 
     * vdW-radii in Angstrom, index number = atomic number
     */
    public double[] getVdWRadii() {
        
        for (int i = 1; i < 93; i++) {
            if (PeriodicTable
                    .getVdwRadius(PeriodicTable.getSymbol(i)) != null) {
                vdWRadii[i] = PeriodicTable
                    .getVdwRadius(PeriodicTable.getSymbol(i));
            }
        }
        
        return vdWRadii;
    }
    
    /**
     * Read the box length from a given .arc file
     * @param aPath Pathname of .arc file
     * @return Box length in Angstrom
     */
    public double getBoxLength(Path aPath) {
        double boxLength = 0.0;
        try (BufferedReader reader = Files.newBufferedReader(aPath, 
                StandardCharsets.UTF_8)) {
            reader.readLine();
            String line = reader.readLine();
            if (line != null) {
                String[] lineArray = line.trim().split("\\s+", 2);
                boxLength = Double.parseDouble(lineArray[0]);
            }
        } catch (FileNotFoundException ex) {
            throw new IllegalArgumentException("File not found in getBoxLength.");
        } catch (IOException ex) {
            throw new IllegalArgumentException("IOException in getBoxLength.");            
        }
        return boxLength;
    }
    
    /**
     * Checks if a given double value is in a list
     * @param aList
     *   A List with double values
     * @param aValue
     *   A value
     * @return 
     *   True if the value is in the list (otherwise false)
     */
    public boolean contains(List<Double> aList, double aValue) {
        double threshold = 0.001;
        
        for (double value : aList) {
            if (Math.abs(value - aValue) < threshold) {
                return true;
            }
        }
    
        return false;
    }
    
    /**
     * Get method for the neighbor particle number of solute
     * 
     * @param aCoordRecord
     *   Coordinate data of solute and solvent particles
     * @param anElements1
     *   Element names of solute particle
     * @param anElements2
     *   Element names of solvent particles
     * @param aBoxLength
     *   Length of simulation box
     * @param aCatchRadius
     *   Catch radius
     * @return 
     *   Neighbor particle number of solute of each simulation iteration
     */
    public LinkedList<int[]> getNeighborNumbersBruteForce(
            CoordinatesRecord aCoordRecord,
            String[] anElements1,
            String[] anElements2,
            double aBoxLength,
            double aCatchRadius) {
        
        final double ONEHALF = 0.5;
        LinkedList<int[]> tmpResult;
        double[][][] tmpCoord1 = aCoordRecord.coord1();
        double[][][][] tmpCoord2 = aCoordRecord.coord2();
        double tmpDeltaX;
        double tmpDeltaY;
        double tmpDeltaZ;
        double tmpDistQ;
        double tmpMinDist;
        double tmpMinDistQ;
        int tmpIndex;
        int tmpAtomSize1;
        int tmpAtomSize2;
        int tmpPartSize2;
        int tmpIterSize;
        int[] tmpNeighborIndices;
        
        tmpResult = new LinkedList<>();
        tmpAtomSize1 = aCoordRecord.coord1()[0].length;
        tmpAtomSize2 = aCoordRecord.coord2()[0][0].length;
        tmpPartSize2 = aCoordRecord.coord2()[0].length;
        tmpIterSize = aCoordRecord.coord1().length;
        double[] tmpVdWRadii1 = new double[tmpAtomSize1];
        double[] tmpVdWRadii2 = new double[tmpAtomSize2];
        HashSet<String> tmpIgnoreElements = new HashSet<>();
        HashSet<Integer> tmpNeighborPart = new HashSet<>();
        HashSet<Integer> tmpIgnoreIndex1 = new HashSet<>();
        HashSet<Integer> tmpIgnoreIndex2 = new HashSet<>();
        
        tmpIgnoreElements.add("M");
        tmpIgnoreElements.add("Lp");
        for (int i = 0; i < tmpAtomSize1; i++) {
            if (! tmpIgnoreElements.contains(anElements1[i])) {
                tmpVdWRadii1[i] = vdWRadii[atomicNumber.get(anElements1[i])];
            } else {
                tmpIgnoreIndex1.add(i);
            }
        }
        
        for (int i = 0; i < tmpAtomSize2; i++) {
            if (! tmpIgnoreElements.contains(anElements2[i])) {
                tmpVdWRadii2[i] = vdWRadii[atomicNumber.get(anElements2[i])];
            } else {
                tmpIgnoreIndex2.add(i);
            }
        }
        
        for (int i = 0; i < tmpIterSize; i++) {

            for (int j = 0; j < tmpAtomSize1; j++) {
                if (!tmpIgnoreIndex1.contains(j)) {
                    for (int k = 0; k < tmpPartSize2; k++) {
                        if (tmpNeighborPart.contains(k)) {
                            continue;
                        }

                        for (int l = 0; l < tmpAtomSize2; l++) {
                            if (!tmpIgnoreIndex2.contains(l)) {
                                tmpDeltaX = tmpCoord1[i][j][0] 
                                        - tmpCoord2[i][k][l][0];
                                tmpDeltaY = tmpCoord1[i][j][1] 
                                        - tmpCoord2[i][k][l][1];
                                tmpDeltaZ = tmpCoord1[i][j][2] 
                                        - tmpCoord2[i][k][l][2];
                                if(tmpDeltaX > aBoxLength * ONEHALF)
                                    tmpDeltaX -= aBoxLength;
                                if(tmpDeltaX <= -aBoxLength * ONEHALF)
                                    tmpDeltaX += aBoxLength;
                                if(tmpDeltaY > aBoxLength * ONEHALF)
                                    tmpDeltaY -= aBoxLength;
                                if(tmpDeltaY <= -aBoxLength * ONEHALF)
                                    tmpDeltaY += aBoxLength;
                                if(tmpDeltaZ > aBoxLength * ONEHALF)
                                    tmpDeltaZ -= aBoxLength;
                                if(tmpDeltaZ <= -aBoxLength * ONEHALF)		
                                    tmpDeltaZ += aBoxLength;
                                tmpDistQ = tmpDeltaX * tmpDeltaX 
                                        + tmpDeltaY * tmpDeltaY 
                                        + tmpDeltaZ * tmpDeltaZ;
                                tmpMinDist = aCatchRadius + tmpVdWRadii1[j] 
                                        + tmpVdWRadii2[l];
                                tmpMinDistQ = tmpMinDist * tmpMinDist;
                                if (tmpDistQ <= tmpMinDistQ) {
                                    tmpNeighborPart.add(k);
                                    break;
                                } 
                            }
                        }
                    
                    }
                    
                }
            }
            
            tmpIndex = 0;
            tmpNeighborIndices = new int[tmpNeighborPart.size()];
            
            for (int tmpCN : tmpNeighborPart) {
                tmpNeighborIndices[tmpIndex++] = tmpCN;
            }
            
            Arrays.sort(tmpNeighborIndices);
            tmpResult.add(tmpNeighborIndices);
            tmpNeighborPart.clear();
        }
        
        return tmpResult;
    }
    
    /**
     * Get method for the neighbor particle number of solute
     * 
     * @param aCoordRecord
     *   Coordinate data of solute and solvent particles
     * @param anElements1
     *   Element names of solute particle
     * @param anElements2
     *   Element names of solvent particles
     * @param aBoxDist
     *   Length of simulation box
     * @param aCatchRadius
     *   Catch radius
     * @return 
     *   Neighbor particle number of solute of each simulation iteration
     */
    public LinkedList<int[]> getNeighborNumbers(
            CoordinatesRecord aCoordRecord,
            String[] anElements1,
            String[] anElements2,
            double aBoxDist,
            double aCatchRadius) {
        final double ONEHALF = 0.5;
        LinkedList<int[]> tmpResult = new LinkedList<>();
        int tmpIndex;
        int tmpCellsInRow;
        int tmpCellIndex;
        int tmpCellX;
        int tmpCellY;
        int tmpCellZ;
        int tmpAtomNumber1;
        int tmpAtomNumber2;
        int tmpPartNumber2;
        int tmpIterNumber;
        int[] tmpNeighborIndices;
        double tmpDeltaX;
        double tmpDeltaY;
        double tmpDeltaZ;
        double tmpDistQ;
        double tmpMinDist;
        double tmpMinDistQ;
        double tmpMaxRadii1;
        double tmpMaxRadii2;
        double tmpCellDist; // cell edge length in Angstrom
        double[] tmpVdWRadii1;
        double[] tmpVdWRadii2;
        double[][][] tmpCoord1;
        double[][][][] tmpCoord2;
        
        tmpAtomNumber1 = aCoordRecord.coord1()[0].length;
        tmpAtomNumber2 = aCoordRecord.coord2()[0][0].length;
        tmpPartNumber2 = aCoordRecord.coord2()[0].length;
        tmpIterNumber = aCoordRecord.coord1().length;
        tmpCoord1 = aCoordRecord.coord1();
        tmpCoord2 = aCoordRecord.coord2();
        tmpVdWRadii1 = new double[tmpAtomNumber1];
        tmpVdWRadii2 = new double[tmpAtomNumber2];
        HashSet<Integer> tmpNeighborCells = new HashSet<>();
        HashSet<Integer> tmpNeighborPartCand = new HashSet<>();
        HashSet<Integer> tmpNeighborPart = new HashSet<>();
        
        for (int i = 0; i < tmpAtomNumber1; i++) {
            tmpVdWRadii1[i] = vdWRadii[atomicNumber.get(anElements1[i])];
        }
        
        for (int i = 0; i < tmpAtomNumber2; i++) {
            tmpVdWRadii2[i] = vdWRadii[atomicNumber.get(anElements2[i])];
        }
        
        // Determine cell edge length and cell number
        tmpMaxRadii1 = this.getMax(tmpVdWRadii1);
        tmpMaxRadii2 = this.getMax(tmpVdWRadii2);
        tmpCellDist = aBoxDist / 
                (tmpMaxRadii1 + tmpMaxRadii2 + aCatchRadius);
        tmpCellsInRow = (int) tmpCellDist;
        tmpCellDist = aBoxDist / tmpCellsInRow;
        
        for (int i = 0; i < tmpIterNumber; i++) {
            tmpCellIndex = 1;
            
            // Determine the solute particle cells
            for (int j = 0; j < tmpAtomNumber1; j++) {
                tmpCellX = (int) ((tmpCoord1[i][j][0] + aBoxDist / 2) 
                        / tmpCellDist) + 1;
                if(tmpCellX > tmpCellsInRow) {
                    tmpCellX = tmpCellsInRow;
                }
                tmpCellY = (int) ((tmpCoord1[i][j][1] + aBoxDist / 2) 
                        / tmpCellDist) + 1;
                if(tmpCellY > tmpCellsInRow) {
                    tmpCellY = tmpCellsInRow;
                }
                tmpCellZ = (int) ((tmpCoord1[i][j][2] + aBoxDist / 2) 
                        / tmpCellDist) + 1;
                if(tmpCellZ > tmpCellsInRow) {
                    tmpCellZ = tmpCellsInRow;
                }
                tmpCellIndex = tmpCellX + (tmpCellY - 1) * tmpCellsInRow +
                        (tmpCellZ - 1) * tmpCellsInRow * tmpCellsInRow;
                tmpNeighborCells.add(tmpCellIndex);
            }
            
            // Determine the neighbor cells
            tmpNeighborCells = this
                    .getNeighborCellNumbers(tmpCellsInRow, tmpNeighborCells);
            
            // Determine the neighbor particles
            for (int j = 0; j < tmpPartNumber2; j++) {
                
                for (int k = 0; k < tmpAtomNumber2; k++) {
                    tmpCellX = (int) ((tmpCoord2[i][j][k][0] + aBoxDist / 2)
                            / tmpCellDist) + 1;
                    if(tmpCellX > tmpCellsInRow) {
                        tmpCellX = tmpCellsInRow;
                    }
                    tmpCellY = (int) ((tmpCoord2[i][j][k][1] + aBoxDist / 2) 
                            / tmpCellDist) + 1;
                    if(tmpCellY > tmpCellsInRow) {
                        tmpCellY = tmpCellsInRow;
                    }
                    tmpCellZ = (int) ((tmpCoord2[i][j][k][2] + aBoxDist / 2) 
                            / tmpCellDist) + 1;
                    if(tmpCellZ > tmpCellsInRow) {
                        tmpCellZ = tmpCellsInRow;
                    }
                    tmpCellIndex = tmpCellX + (tmpCellY - 1) * tmpCellsInRow +
                        (tmpCellZ - 1) * tmpCellsInRow * tmpCellsInRow;
                    if (tmpNeighborCells.contains(tmpCellIndex)) {
                        if(tmpNeighborPartCand.contains(j)) {
                            break;
                        } else {
                            tmpNeighborPartCand.add(j);
                            break;
                        }
                    }
                }
            }
            
            for (int j = 0; j < tmpAtomNumber1; j++) {
                
                for (int k = 0; k < tmpPartNumber2; k++) {
                    if (tmpNeighborPartCand.contains(k)) {
                        for (int l = 0; l < tmpAtomNumber2; l++) {
                            tmpDeltaX = tmpCoord1[i][j][0] 
                                    - tmpCoord2[i][k][l][0];
                            tmpDeltaY = tmpCoord1[i][j][1] 
                                    - tmpCoord2[i][k][l][1];
                            tmpDeltaZ = tmpCoord1[i][j][2] 
                                    - tmpCoord2[i][k][l][2];
                            if(tmpDeltaX > aBoxDist * ONEHALF)
                                tmpDeltaX -= aBoxDist;
                            if(tmpDeltaX <= -aBoxDist * ONEHALF)
                                tmpDeltaX += aBoxDist;
                            if(tmpDeltaY > aBoxDist * ONEHALF)
                                tmpDeltaY -= aBoxDist;
                            if(tmpDeltaY <= -aBoxDist * ONEHALF)
                                tmpDeltaY += aBoxDist;
                            if(tmpDeltaZ > aBoxDist * ONEHALF)
                                tmpDeltaZ -= aBoxDist;
                            if(tmpDeltaZ <= -aBoxDist * ONEHALF)		
                                tmpDeltaZ += aBoxDist;
                            tmpDistQ = tmpDeltaX * tmpDeltaX + tmpDeltaY 
                                    * tmpDeltaY + tmpDeltaZ * tmpDeltaZ;
                            tmpMinDist = aCatchRadius + tmpVdWRadii1[j] 
                                    + tmpVdWRadii2[l];
                            tmpMinDistQ = tmpMinDist * tmpMinDist;
                            if (tmpDistQ <= tmpMinDistQ) {
                                tmpNeighborPart.add(k);
                                break;
                            }
                        }
                        
                    }
                }
                
            }
            
            tmpIndex = 0;
            tmpNeighborIndices = new int[tmpNeighborPart.size()];
            
            for (int tmpCN : tmpNeighborPart) {
                tmpNeighborIndices[tmpIndex++] = tmpCN;
            }
            
            Arrays.sort(tmpNeighborIndices);
            tmpResult.add(tmpNeighborIndices);
            tmpNeighborPart.clear();
        }
        
        return tmpResult;
    }
    
    /**
     * Determine cell numbers solute particle and its neighbor cells
     * @param aCellsInRow Numbers of cells per edge
     * @param aSoluteCellIndices All cell numbers of of located solute particle
     * @return Cell numbers of solute particle
     */
    public HashSet<Integer> getNeighborCellNumbers(int aCellsInRow, 
            HashSet<Integer> aSoluteCellIndices) {
        HashSet<Integer> tmpNeighborIndices;
        int[][][] tmpIndexCube;
        int tmpCellIndex;
        
        tmpNeighborIndices = new HashSet<>();
        tmpIndexCube = 
                new int[aCellsInRow + 2][aCellsInRow + 2][aCellsInRow + 2];
        tmpCellIndex = 1;
        
        for (int i = 1; i <= aCellsInRow; i++) {
            
            for (int j = 1; j <= aCellsInRow; j++) {
                
                for (int k = 1; k <= aCellsInRow; k++) {
                    tmpIndexCube[i][j][k] = tmpCellIndex;
                    tmpCellIndex++;
                }
                
            }
            
        }
        
        // bottom + top face
        for (int j = 1; j <= aCellsInRow; j++) {
            
            for (int k = 1; k <= aCellsInRow; k++) {
                tmpIndexCube[0][j][k] = tmpIndexCube[aCellsInRow][j][k];
                tmpIndexCube[aCellsInRow + 1][j][k] = tmpIndexCube[1][j][k];
            }
            
        }
        
        // left + right face
        for (int i = 0; i <= aCellsInRow + 1; i++) {
            
            for (int j = 1; j <= aCellsInRow; j++) {
                tmpIndexCube[i][j][0] = tmpIndexCube[i][j][aCellsInRow];
                tmpIndexCube[i][j][aCellsInRow + 1] = tmpIndexCube[i][j][1];
            }
            
        }
        
        // front + back face
        for (int i = 0; i <= aCellsInRow + 1; i++) {
            
            for (int k = 0; k <= aCellsInRow + 1; k++) {
                tmpIndexCube[i][0][k] = tmpIndexCube[i][aCellsInRow][k];
                tmpIndexCube[i][aCellsInRow + 1][k] = tmpIndexCube[i][1][k];
            }
        }
        
        for (int tmpSoluteCellIndex : aSoluteCellIndices) {
            tmpCellIndex = 1;
            
            for (int i = 1; i <= aCellsInRow; i++) {
            
                for (int j = 1; j <= aCellsInRow; j++) {
                
                    for (int k = 1; k <= aCellsInRow; k++) {
                        if (tmpSoluteCellIndex == tmpCellIndex) {
                            tmpNeighborIndices.add(tmpIndexCube[i][j][k]);
                            tmpNeighborIndices.add(tmpIndexCube[i][j][k + 1]);
                            tmpNeighborIndices.add(tmpIndexCube[i][j][k - 1]);
                            tmpNeighborIndices.add(tmpIndexCube[i][j + 1][k]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i][j + 1][k + 1]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i][j + 1][k - 1]);
                            tmpNeighborIndices.add(tmpIndexCube[i][j - 1][k]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i][j - 1][k + 1]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i][j - 1][k - 1]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i + 1][j][k]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i + 1][j][k + 1]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i + 1][j][k - 1]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i + 1][j + 1][k]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i + 1][j + 1][k + 1]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i + 1][j + 1][k - 1]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i + 1][j - 1][k]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i + 1][j - 1][k + 1]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i + 1][j - 1][k - 1]);
                            tmpNeighborIndices.add(tmpIndexCube[i - 1][j][k]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i - 1][j][k + 1]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i - 1][j][k - 1]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i - 1][j + 1][k]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i - 1][j + 1][k + 1]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i - 1][j + 1][k - 1]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i - 1][j - 1][k]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i - 1][j - 1][k + 1]);
                            tmpNeighborIndices
                                    .add(tmpIndexCube[i - 1][j - 1][k - 1]);
                        }
                        tmpCellIndex++;
                    }
                
                }
            
            }
        }
        
        return tmpNeighborIndices;
    }
    
    /**
     * Read the coordinates from an .arc file
     * 
     * @param path File name
     * @param anAtomNumber1 Atom number of solute particle
     * @param anAtomNumber2 Atom number of solvent particle
     * @return All coordinates of .arc file
     */
    public CoordinatesRecord getCoordinatesFromArcFile(Path path,
            int anAtomNumber1, int anAtomNumber2) {
        // Check parameters
        if (path == null) {
            throw new IllegalArgumentException("Null was passed to the " 
                    + "setFileContent method.");
        }
        if (anAtomNumber1 <= 0 || anAtomNumber2 <= 0) {
            throw new IllegalArgumentException("Illegal numbers" 
                    + "passed in getCoordinatesFromArcFile.");
        }
        
        int tmpLength1;
        int tmpLength2;
        int tmpLineCounter; 
        int tmpIndex;
        int tmpIndex1;
        int tmpIndex2;
        int tmpParticleIndex2;
        int tmpParticleNumber2;
        int tmpAtomNumber;
        double[][] tmpCoord1;
        double[][][] tmpCoord2;
        double[][][] tmpCoordinate1;
        double[][][][] tmpCoordinate2;
        String tmpLine;
        String[] tmpLineArray;
        LinkedList<double[][]> tmpCoordList1 = new LinkedList<>();
        LinkedList<double[][][]> tmpCoordList2 = new LinkedList<>();
        CoordinatesRecord tmpCoordRecord;
        
        tmpLineCounter = 1; 
        tmpIndex = 0;
        tmpIndex1 = 0;
        tmpIndex2 = 0;
        tmpParticleIndex2 = 0;
        tmpParticleNumber2 = 0;
        tmpCoord1 = new double[anAtomNumber1][3];
        try (BufferedReader tmpBR = Files.newBufferedReader(path,
                StandardCharsets.UTF_8)) {
            tmpBR.mark(80);
            tmpLine = tmpBR.readLine();
            tmpAtomNumber = Integer.parseInt(tmpLine.substring(0, 6).trim());
            tmpBR.reset();
            tmpParticleNumber2 = (tmpAtomNumber - anAtomNumber1) 
                    / anAtomNumber2;
            tmpCoord2 = new double[tmpParticleNumber2][anAtomNumber2][3];
            
            while((tmpLine = tmpBR.readLine()) != null) {
                if (tmpIndex >= 2 && tmpIndex <= anAtomNumber1 + 1) {
                    tmpLineArray = this.split(tmpLine.trim());
                    tmpCoord1[tmpIndex1][0] = Double.parseDouble(
                            tmpLineArray[0]);
                    tmpCoord1[tmpIndex1][1] = Double.parseDouble(
                            tmpLineArray[1]);
                    tmpCoord1[tmpIndex1][2] = Double.parseDouble(
                            tmpLineArray[2]);
                    tmpLineCounter++;
                    tmpIndex1++;
                    if (tmpIndex1 >= anAtomNumber1) {
                        tmpIndex1 = 0;
                    }
                } else if (tmpIndex > anAtomNumber1 + 1) {
                    tmpLineArray = this.split(tmpLine.trim());
                    tmpCoord2[tmpParticleIndex2][tmpIndex2][0] = Double
                            .parseDouble(tmpLineArray[0]);
                    tmpCoord2[tmpParticleIndex2][tmpIndex2][1] = Double
                            .parseDouble(tmpLineArray[1]);
                    tmpCoord2[tmpParticleIndex2][tmpIndex2][2] = Double
                            .parseDouble(tmpLineArray[2]);
                    tmpLineCounter++;
                    tmpIndex2++;
                    if (tmpIndex2 >= anAtomNumber2) {
                        tmpIndex2 = 0;
                        tmpParticleIndex2++;
                        if (tmpParticleIndex2 >= tmpParticleNumber2) {
                            tmpParticleIndex2 = 0;
                            tmpCoordList1.add(tmpCoord1);
                            tmpCoordList2.add(tmpCoord2);
                            tmpCoord1 = new double[anAtomNumber1][3];
                            tmpCoord2 = new double[tmpParticleNumber2]
                                    [anAtomNumber2][3];
                        }
                    }
                } else {
                    tmpLineCounter++;
                }
                tmpIndex = (tmpLineCounter - 1) % (tmpAtomNumber + 2);
            }
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException in getCoordinatesFromArcFile.", ex);
        }
        tmpLength1 = tmpCoordList1.size();
        tmpLength2 = tmpCoordList2.size();
        tmpCoordinate1 = new double[tmpLength1][anAtomNumber1][3];
        tmpCoordinate2 =
                new double[tmpLength2][tmpParticleNumber2][anAtomNumber1][3];
        
        for (int i = 0; i < tmpLength1; i++) {
            tmpCoordinate1[i] = tmpCoordList1.get(i);
            tmpCoordinate2[i] = tmpCoordList2.get(i);
        }
        
        tmpCoordRecord = new CoordinatesRecord(tmpCoordinate1, tmpCoordinate2);
        return tmpCoordRecord;
    }
    
    /**
     * Determine (widest) molecular diameter 
     * 
     * @param aTinkerXYZ TinkerXYZ object
     * @return (Widest) molecular diameter
     */
    public double getMolecularDiameter(TinkerXYZ aTinkerXYZ) {
        int tmpAtomSize;
        final int tmpAtomicNumber1;
        int tmpAtomicNumber2;
        int tmpIndex1;
        int tmpIndex2;
        final double tmpVdW1;
        final double tmpVdW2;
        double tmpDeltaX;
        double tmpDeltaY;
        double tmpDeltaZ;
        double tmpMaxDistanceQ;
        double tmpMaxDistanceQCandidate;
        double tmpReturn;
        double[][] tmpCoordinates;
        String[] tmpElements;
        
        tmpAtomSize = aTinkerXYZ.getN_atom1();
        tmpElements = aTinkerXYZ.getElementList1();
        tmpCoordinates = aTinkerXYZ.getCoordinateList1()[0];
        tmpMaxDistanceQ = 0.0;
        tmpIndex1 = 0;
        tmpIndex2 = 0;
        
        for (int i = 0; i < tmpAtomSize; i++) {
            
            for (int j = i + 1; j < tmpAtomSize; j++) {
                tmpDeltaX = tmpCoordinates[i][0] - tmpCoordinates[j][0];
                tmpDeltaY = tmpCoordinates[i][1] - tmpCoordinates[j][1];
                tmpDeltaZ = tmpCoordinates[i][2] - tmpCoordinates[j][2];
                tmpMaxDistanceQCandidate = tmpDeltaX * tmpDeltaX + 
                        tmpDeltaY * tmpDeltaY +
                        tmpDeltaZ * tmpDeltaZ;
                if (tmpMaxDistanceQCandidate > tmpMaxDistanceQ) {
                    tmpMaxDistanceQ = tmpMaxDistanceQCandidate;
                    tmpIndex1 = i;
                    tmpIndex2 = j;
                }
            }
            
        }
        
        tmpReturn = Math.sqrt(tmpMaxDistanceQ);
        tmpAtomicNumber1 = PeriodicTable
                .getAtomicNumber(tmpElements[tmpIndex1]);
        tmpAtomicNumber2 = PeriodicTable
                .getAtomicNumber(tmpElements[tmpIndex2]);
        tmpVdW1 = this.getVdWRadii()[tmpAtomicNumber1];
        tmpVdW2 = this.getVdWRadii()[tmpAtomicNumber2];
        tmpReturn += tmpVdW1 + tmpVdW2;
        return tmpReturn;
    }
    
    /**
     * Write the last part of .arc file to .xyz file
     * @param anArcPath
     *   Path name of the .arc file
     * @param aXyzPath
     *   Path name of the .xyz file
     * @param aSteps 
     *   Step number in the .arc file
     */
    public void writeLastPartToXYZ(Path anArcPath, Path aXyzPath, int aSteps) {
        try (BufferedReader tmpBR = Files.newBufferedReader(anArcPath, 
                StandardCharsets.UTF_8);
                BufferedWriter tmpBW = Files.newBufferedWriter(aXyzPath, 
                        StandardCharsets.UTF_8)) {
            int headerSize = 2;
            String line = tmpBR.readLine();
            int atomSize = Integer.parseInt(line.substring(0, 6).trim());
            int ignoreLineSize = (atomSize + headerSize) * (aSteps - 1) - 1;
            
            for (int i = 0; i < ignoreLineSize; i++) {
                tmpBR.readLine();
            }
            
            while((line = tmpBR.readLine()) != null) {
                tmpBW.append(line);
                tmpBW.append(LINESEPARATOR);
            }
                
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException in writeLastPartToXYZ.", ex);
        }
    }
    
    /**
     * Get SMILES data method
     *   Loads the SMILES code of the "known" molecules
     * @param aSmilesFileName: File name of SMILES
     * @return Molecule name and SMILES as hashmap
     */
    public HashMap<String, String> getSmilesData(String aSmilesFileName) {
        HashMap <String, String> tmpSmiles = new HashMap<>();
        Path path = Paths.get(aSmilesFileName);
        try (BufferedReader tmpBR = Files.newBufferedReader(path,
                StandardCharsets.UTF_8)) {
            String line;
            
            while((line = tmpBR.readLine()) != null) {
                String[] lineArray = line.trim().split("\\s+");
                tmpSmiles.put(lineArray[0], lineArray[1]);
            }
            
        } catch (FileNotFoundException ex) {
           LOGGER.log(Level.SEVERE, 
                    "FileNotFoundException in getSmilesData.", ex);
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException in getSmilesData.", ex);
        }
        return tmpSmiles;
    }
    
    /**
     * Determines if the smallest distance between two particles 
     *   is less than minDistance or not. This method is used to exclude
     *   two particle configurations, which are too close. This would produce
     *   too high value of intermolecular energy which will be expressed 
     *   by asterisks in tinker's output file.
     * 
     * @param aCoord1
     *   Coordinates of atoms of first particle
     * @param aCoord2
     *   Coordinates of atoms of second particle
     * @param aMinDist
     *   Minimum distance
     * @return 
     *   true if the distance of two particles is too close
     *   false otherwise
     */
    public Boolean isTooClose(double[][] aCoord1, double[][] aCoord2, 
            double aMinDist) {
        double minDistQ = aMinDist * aMinDist;
        
        for (double[] p1 : aCoord1) {
            
            for (double[] p2 : aCoord2) {
                double dx = p1[0] - p2[0];
                double dy = p1[1] - p2[1];
                double dz = p1[2] - p2[2];
                double distanceQ = dx * dx + dy * dy + dz * dz;
                if(distanceQ < minDistQ){
                    return true;
                }
            }
            
        }
        
        return false;
    }
    
    /**
     * Checks if a path is empty or not
     * 
     * @param path Directory name
     * @return true: is empty, false: isn't empty
     */
    public boolean isDirectoryEmpty(String path){
        Path directory = Paths.get(path);
        if (Files.isDirectory(directory)) {
            try (Stream<Path> entries = Files.list(directory)) {
                return !entries.findFirst().isPresent();
            } catch (IOException ex) {
                System.err.println("Error during access to the directory: " 
                    + ex.getMessage());
            }
        }
        return false;
    }
    
    /**
     * Returns all Lines with aSearchString
     * @param fileName
     *   A Text file
     * @param searchString
     *   A search string
     * @return 
     *   All lines with aSearchString in it as a list.
     */
    public List<String> findList(String fileName, String searchString) {
        Path filePath = Paths.get(fileName);
        try (Stream<String> lines = Files.lines(filePath)) {
            return lines.filter(line -> line.contains(searchString))
                    .collect(Collectors.toList());

        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, "IOException in findList.", ex);
            return new ArrayList<>();
        }
    }

    /**
     * Calculates the sum of the array elements.
     * 
     * @param aDoubleArray A double array
     * @return Sum of the array elements
     */
    public double sum(double[] aDoubleArray) {
        double tmpSum = 0;

        for (double tmpValue : aDoubleArray) {
            tmpSum += tmpValue;
        }
        
        return tmpSum;
    }
    
    /**
     * Calculates the sum of the array elements.
     * 
     * @param anIntArray A double array
     * @return Sum of the array elements
     */
    public long sum(int[] anIntArray) {
        int tmpSum = 0;

        for (int tmpValue : anIntArray) {
            tmpSum += tmpValue;
        }
        
        return tmpSum;
    }
    
    /**
     * Calculates the sum of element products of two arrays.
     * 
     * @param aDoubleArray1 A double array
     * @param aDoubleArray2 Another double array
     * @return Sum of element products of two arrays
     */
    public double productSum(double[] aDoubleArray1, double[] aDoubleArray2) {
        double tmpSum = 0;
        
        for (int i = 0; i < aDoubleArray1.length; i++) {
            tmpSum += aDoubleArray1[i] * aDoubleArray2[i];
        }
        return tmpSum;
    }
    
    /**
     * Calculates the standard deviation of values
     * @param aValues values
     * @param aMean mean of values
     * @return standard deviation
     */
    public double standarddeviation(int[] aValues, double aMean) {
        double tmpSum = 0;
        double tmpReturnValue;

        for (int tmpValue : aValues) {
            tmpSum += (tmpValue - aMean) * (tmpValue - aMean);
        }
        
        tmpReturnValue = Math.sqrt(1.0 / (aValues.length - 1) * tmpSum);
        return tmpReturnValue;
    }
    
    /**
     * Determine the biggest value from an integer array
     * @param aValues integer array
     * @return biggest value of an integer array
     */
    public int getMax(int[] aValues) {
        int tmpMax = aValues[0];
        
        for (int i = 1; i < aValues.length; i++) {
            if (tmpMax < aValues[i]) {
                tmpMax = aValues[i];
            }
        }
        
        return tmpMax;
    }
    
    /**
     * Determine the biggest value from a double array
     * @param aValues double array
     * @return biggest value of a double array
     */
    public double getMax(double[] aValues) {
        double tmpMax = aValues[0];
        
        for (int i = 1; i < aValues.length; i++) {
            if (tmpMax < aValues[i]) {
                tmpMax = aValues[i];
            }
        }
        
        return tmpMax;
    }
    
    
    /**
     * Determine the smallest value of an integer array
     * @param aValues integer array
     * @return smallest value of an integer array
     */
    public int getMin(int[] aValues) {
        int tmpMin = aValues[0];
        
        for (int i = 1; i < aValues.length; i++) {
            if (tmpMin > aValues[i]) {
                tmpMin = aValues[i];
            }
        }
        
        return tmpMin;
    }
    
    /**
     * Determine the distance of two points
     * @param aCoord1 xyz coordinates of point1
     * @param aCoord2 xyz coordinates of point2
     * @return distance of point1 and point2
     */
    public double getDistance(double[] aCoord1, double[] aCoord2) {
        double tmpDeltaX;
        double tmpDeltaY;
        double tmpDeltaZ;
        double tmpResult;
        
        tmpDeltaX = aCoord2[0] - aCoord1[0];
        tmpDeltaY = aCoord2[1] - aCoord1[1];
        tmpDeltaZ = aCoord2[2] - aCoord1[2];
        tmpResult = Math.sqrt(tmpDeltaX * tmpDeltaX + 
                tmpDeltaY * tmpDeltaY +
                tmpDeltaZ * tmpDeltaZ);
        return tmpResult;
    }
    
    /**
     * Determine the distances to shortest neighbor
     *   This is only for test purpose
     * @return Distances to shortest neighbor
     */
    public double[] getNextDistance() {
        String sphereNodeFileName = "resources/de/whs/ibci/mipet/sphereNodes/SphereNodes";
        int[] sphereNodeNumber = new int[] {4, 9, 16, 25, 36, 49, 64, 81, 100, 
            121, 144, 169, 196, 225, 256, 289, 324, 361, 400, 441, 484, 529, 
            576, 625, 676, 729, 784, 841, 900};
        double[] result = new double[29];
        double shortestDistQ;
        
        for (int i = 0; i < 29; i++) {
            int tmpMolNumber = sphereNodeNumber[i];
            String tmpFileName = sphereNodeFileName + tmpMolNumber + ".txt";
            Path path = Paths.get(tmpFileName);
            try (BufferedReader tmpBR = Files.newBufferedReader(path,
                    StandardCharsets.UTF_8)) {
                double[][] matrix = new double[sphereNodeNumber[i]][3];
                int atomIndex = 0;
                shortestDistQ = 1000.;
                String tmpLine;
                
                while ((tmpLine = tmpBR.readLine()) != null) {
                    int startIndex = tmpLine.indexOf("{");
                    int endIndex = tmpLine.lastIndexOf("}");
                    String subString = tmpLine
                            .substring(startIndex + 1, endIndex);
                    String[] tmpTokens = subString.split(",");
                    matrix[atomIndex][0] = Double.parseDouble(
                            tmpTokens[0].trim());
                    matrix[atomIndex][1] = Double.parseDouble(
                            tmpTokens[1].trim());
                    matrix[atomIndex][2] = Double.parseDouble(
                            tmpTokens[2].trim());
                    atomIndex++;
                }
                
                for (int j = 1; j < tmpMolNumber; j++) {
                    double distQX = matrix[j][0] - matrix[0][0];
                    double distQY = matrix[j][1] - matrix[0][1];
                    double distQZ = matrix[j][2] - matrix[0][2];
                    double distQ = distQX * distQX
                            + distQY * distQY
                            + distQZ * distQZ;
                    if (distQ < shortestDistQ) {
                        shortestDistQ = distQ;
                    }
                }
                
                result[i] = Math.sqrt(shortestDistQ);
            } catch(IOException ex) {
            }
        }
        
        return result;
    }
    
    /**
     * Write status message in console
     * 
     * @param message Status message
     */
    public static void updateStatus(String message) {
        // \r Goto start position, \u001b[K Clear the rest
        System.out.print("\r\u001b[K" + message);
        System.out.flush(); // Force immediatly output
    }
    
    /**
     * Read a part of .arc file and returns the content as StringBuilder object.
     * 
     * @param aPath
     *   A Path name
     * @param aStartIndex
     *   Zero-based start line index
     * @param aEndIndex
     *   Zero-based end line index
     * @return 
     * Content of the selected part of .arc file
     */
    public static StringBuilder readPartArcFile(Path aPath, 
            int aStartIndex, int aEndIndex) {
        String line;
        int index = 0;
        StringBuilder strBuilder = new StringBuilder();
        try (BufferedReader reader = Files.newBufferedReader(aPath, 
                StandardCharsets.UTF_8)) {
            
            while ((line = reader.readLine()) != null) {
                if (aStartIndex <= index && index <= aEndIndex) {
                    strBuilder.append(line);
                    strBuilder.append(LINESEPARATOR);
                }
                index++;
            }
            
        } catch (IOException ex) {
                LOGGER.log(Level.SEVERE,
                        "IOException during reading" + aPath.toString() , ex);
        }
        return strBuilder;
    }
    
    /**
     * Write Particle logfiles
     * 
     * @param aJobTaskRecords
     *   First Particle name, second Particle name, result directory name etc.
     * @param aLabelValues
     *   0: Label 1: Value
     */
    public static void writeParticleLog(ArrayList<JobTaskRecord> aJobTaskRecords,
            String[][] aLabelValues) {
        int jobTaskLength = aJobTaskRecords.size();
        
        for (int i = 0; i < jobTaskLength; i++) {
            if (aJobTaskRecords.get(i).hasCNJob()) {
                String particle1 = aJobTaskRecords.get(i).particleName1(); // solute
                String particle2 = aJobTaskRecords.get(i).particleName2(); // solvent
                String particlePair = particle1 + "_" + particle2;
                String resultPathName = aJobTaskRecords.get(i)
                        .result_CN_PathName();
                Path logPath = Paths.get(resultPathName, 
                        particlePair + "_log.txt");
                try (BufferedWriter writer = Files.newBufferedWriter(logPath, 
                        StandardCharsets.UTF_8,
                        StandardOpenOption.CREATE,
                        StandardOpenOption.APPEND)) {

                    for (String[] aLabelValue1 : aLabelValues) {
                        writer.append(aLabelValue1[0]);
                        writer.append(aLabelValue1[1]);
                        writer.newLine();
                    }
                    
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during Writing logfile.", ex);
                }
            }
        }
        
    }
    
    /**
     * Write Particlelogfiles
     * 
     * @param aJobTaskRecords
     *   Job task record list
     * @param aLabel
     *   Label name
     * @param aValues
     *   Value names
     */
    public void writeParticleLog(ArrayList<JobTaskRecord> aJobTaskRecords,
            String aLabel, String[] aValues) {
        int jobTaskLength = aJobTaskRecords.size();
        int index = 0;
        
        for (int i = 0; i < jobTaskLength; i++) {
            if (aJobTaskRecords.get(i).hasCNJob()) {
                String particle1 = aJobTaskRecords.get(i).particleName1(); // solute
                String particle2 = aJobTaskRecords.get(i).particleName2(); // solvent
                String particlePair = particle1 + "_" + particle2;
                String resultPathName = aJobTaskRecords.get(i).result_CN_PathName();
                Path logPath = Paths.get(resultPathName, 
                        particlePair + "_log.txt");

                try (BufferedWriter writer = Files.newBufferedWriter(logPath, 
                        StandardCharsets.UTF_8,
                        StandardOpenOption.CREATE,
                        StandardOpenOption.APPEND)) {
                    writer.append(aLabel)
                            .append("(")
                            .append(particle1)
                            .append("/")
                            .append(particle2)
                            .append("): ")
                            .append(aValues[index]);
                    writer.newLine();
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during Writing logfile.", ex);
                }
                index++;
            }
        }
        
    }
    
    /**
     * Write Zij data
     * 
     * @param aJobTaskRecords
     *   Parameter and particle names etc.
     * @param aCNs
     *   Coordination numbers
     * @param aTemperature 
     *   Temperature
     */
    public void writeZij_Table(ArrayList<JobTaskRecord> aJobTaskRecords,
            int[][] aCNs,
            double aTemperature) {
        int jobTaskLength = aJobTaskRecords.size();
        int index = 0;
        
        for (int i = 0; i < jobTaskLength; i++) {
            if (aJobTaskRecords.get(i).hasCNJob()) {
                String particle1 = aJobTaskRecords.get(i).particleName1(); // solute
                String particle2 = aJobTaskRecords.get(i).particleName2(); // solvent
                String particlePair = particle1 + "_" + particle2;
                String targetDirName = aJobTaskRecords.get(i).result_CN_PathName();
                Path zijPath = Paths.get(targetDirName,
                        particlePair + "_ZijTable.dat");
                try (BufferedWriter tmpBW = Files.newBufferedWriter(zijPath, 
                        StandardCharsets.UTF_8,
                        StandardOpenOption.CREATE,
                        StandardOpenOption.APPEND)) {
                    tmpBW.append("temperature [K]: ")
                            .append(Double.toString(aTemperature));
                    tmpBW.newLine();
                    tmpBW.append(Integer.toString(aCNs[index][0]));

                    for (int j = 1; j < aCNs[index].length; j++) {
                        tmpBW.newLine();
                        tmpBW.append(Integer.toString(aCNs[index][j]));
                    }

                    tmpBW.newLine();
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during Writing logfile.", ex);
                }
                index++;
            }
        }
        
    }
    
    /**
     * Get Parameter name and particle pair names
     * @param aPathName
     *   Path name list
     * @return 
     *   0: Parameter name
     *   1: Particle name of first particle
     *   2: Particle name of second particle
     */
    public String[][] getParameterParticleNameList(
            LinkedList<String> aPathName) {
        if (aPathName == null || aPathName.isEmpty())  {
            throw new IllegalArgumentException("Illegal argument was used in getParameterParticleName().");
        }
        LinkedList<String[]> tmpParameterParticleList = new LinkedList<>();

        for (String s : aPathName) {
            int paramStartIndex = -1;
            int particleStartIndex = -1;
            int separateIndex = -1;


            for (int j = s.length() - 1; j >= 0; j--) {
                if (s.charAt(j) == '_') {
                    separateIndex = j;
                } else if (s.charAt(j) ==
                        FILESEPARATOR.charAt(0) &&
                        particleStartIndex == -1) {
                    particleStartIndex = j + 1;
                } else {
                    if (s.charAt(j) == FILESEPARATOR.charAt(0)) {
                        paramStartIndex = j + 1;
                        break;
                    } else if (j == 0) {
                        paramStartIndex = 0;
                    }
                }
            }

            String[] parameterParticle = new String[3];
            String parameter = s.substring(paramStartIndex, 
                    particleStartIndex - 1);
            parameterParticle[0] = parameter;
            String particle1 = s.substring(particleStartIndex, separateIndex);
            parameterParticle[1] = particle1;
            String particle2 = s.substring(separateIndex + 1);
            parameterParticle[2] = particle2;
            tmpParameterParticleList.add(parameterParticle);
            if (!parameterParticle[1].equals(parameterParticle[2])) {
                parameterParticle = new String[3];
                parameterParticle[0] = parameter;
                parameterParticle[1] = particle2;
                parameterParticle[2] = particle1;
                tmpParameterParticleList.add(parameterParticle);
            }
        }
        
        String[][] result = tmpParameterParticleList
                .toArray(new String[tmpParameterParticleList.size()][3]);
        return result;
    }
    
    /**
     * Save .key file
     * @param keyPath
     *   Key path name
     * @param aContent
     *   Content of .key file
     */
    public void writeKeyFile(Path keyPath, String aContent) {
        try (PrintWriter tmpOut = new PrintWriter(Files
                .newBufferedWriter(keyPath, StandardCharsets.UTF_8))) {
            tmpOut.print(aContent);
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during writing .key file.", ex);
        }
    }
    
    /**
     * Write distance and energy values to create a diagram
     * 
     * @param aFileName: File name
     * @param aDistances: Distances
     * @param aDistanceIndices: Indices of distances
     * @param aEnergySorted: Sorted energy values
     */
    public void writeDistance_Energy(String aFileName, Double[] aDistances, 
            Integer[] aDistanceIndices, double[][] aEnergySorted) {
        Path path = Paths.get(aFileName);
        try (BufferedWriter writer = Files.newBufferedWriter(path,
                StandardCharsets.UTF_8)) {
            writer.append("[distances]");
            writer.newLine();
            
            for (int i = 0; i < aDistances.length; i++) {
                writer.append(String.valueOf(aDistances[aDistanceIndices[i]]));
                writer.newLine();
            }
            
            writer.append("[/distances]");
            writer.newLine();
            writer.append("[energies]");
            writer.newLine();
            int tmpEnergiesNumber;
            
            for (double[] aEnergySorted1 : aEnergySorted) {
                tmpEnergiesNumber = aEnergySorted1.length;
                
                for (int j = 0; j < tmpEnergiesNumber; j++) {
                    writer.append(String.valueOf(aEnergySorted1[j]));
                    writer.append("    ");
                }
                
                writer.newLine();
            }
            
            writer.append("[/energies]");
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during writing distanceEnergy file.", ex);
        }
    }
    
    /**
     * Take a String and change an atom type number to another one
     *
     * @param aTokens Tokens of a Line
     * @return String with changed atom type number
     */
    public String changeAtomType(String[] aTokens) {
        String tmpResult;
        
        tmpResult = null;
        switch (aTokens[0]) {
            case "atom" -> {
                aTokens[1] = String.valueOf((Integer.
                        parseInt(aTokens[1]) + 100));
                aTokens[2] = aTokens[1];
                tmpResult = String.join("",
                        aTokens[0], 
                        " ".repeat(8),
                        aTokens[1],
                        " ".repeat(2),
                        aTokens[2],
                        " ".repeat(4),
                        this.padLeft(aTokens[3], 2),
                        " ".repeat(4),
                        aTokens[4],
                        " ".repeat(9),
                        aTokens[5],
                        " ".repeat(5),
                        this.padLeft(aTokens[6], 6),
                        " ".repeat(5),
                        aTokens[7]);
            }
            case "vdw" -> {
                aTokens[1] = String.valueOf((Integer.
                        parseInt(aTokens[1]) + 100));
                tmpResult = String.join("",
                        aTokens[0], 
                        " ".repeat(9),
                        aTokens[1],
                        " ".repeat(11),
                        aTokens[2],
                        " ".repeat(3),
                        aTokens[3]);
            }
            case "bond" -> {
                aTokens[1] = String.valueOf((Integer.
                        parseInt(aTokens[1]) + 100));
                aTokens[2] = String.valueOf((Integer.
                        parseInt(aTokens[2]) + 100));
                tmpResult = String.join("",
                        aTokens[0], 
                        " ".repeat(8),
                        aTokens[1],
                        " ".repeat(2),
                        aTokens[2],
                        " ".repeat(11),
                        aTokens[3],
                        " ".repeat(3),
                        aTokens[4]);
            }
            case "angle" -> {
                aTokens[1] = String.valueOf((Integer.
                        parseInt(aTokens[1]) + 100));
                aTokens[2] = String.valueOf((Integer.
                        parseInt(aTokens[2]) + 100));
                aTokens[3] = String.valueOf((Integer.
                        parseInt(aTokens[3]) + 100));
                tmpResult = String.join("",
                        aTokens[0], 
                        " ".repeat(7),
                        aTokens[1],
                        " ".repeat(2),
                        aTokens[2],
                        " ".repeat(2),
                        aTokens[3],
                        " ".repeat(4),
                        aTokens[4],
                        " ".repeat(3),
                        aTokens[5]);
            }
            case "torsion" -> {
                aTokens[1] = String.valueOf((Integer.
                        parseInt(aTokens[1]) + 100));
                aTokens[2] = String.valueOf((Integer.
                        parseInt(aTokens[2]) + 100));
                aTokens[3] = String.valueOf((Integer.
                        parseInt(aTokens[3]) + 100));
                aTokens[4] = String.valueOf((Integer.
                        parseInt(aTokens[4]) + 100));
                tmpResult = String.join("",
                        aTokens[0], 
                        " ".repeat(5),
                        aTokens[1],
                        " ".repeat(2),
                        aTokens[2],
                        " ".repeat(2),
                        aTokens[3],
                        " ".repeat(2),
                        aTokens[4],
                        " ".repeat(8),
                        aTokens[5],
                        " ".repeat(2),
                        aTokens[6],
                        " ".repeat(2),
                        aTokens[7],
                        " ".repeat(2),
                        aTokens[8],
                        " ".repeat(1),
                        aTokens[9],
                        " ".repeat(2),
                        aTokens[10],
                        " ".repeat(2),
                        aTokens[11],
                        " ".repeat(2),
                        aTokens[12],
                        " ".repeat(2),
                        aTokens[13]);
            }
            case "imptors" -> {
                aTokens[1] = String.valueOf((Integer.
                        parseInt(aTokens[1]) + 100));
                aTokens[2] = String.valueOf((Integer.
                        parseInt(aTokens[2]) + 100));
                aTokens[3] = String.valueOf((Integer.
                        parseInt(aTokens[3]) + 100));
                aTokens[4] = String.valueOf((Integer.
                        parseInt(aTokens[4]) + 100));
                tmpResult = String.join("",
                        aTokens[0], 
                        " ".repeat(5),
                        aTokens[1],
                        " ".repeat(2),
                        aTokens[2],
                        " ".repeat(2),
                        aTokens[3],
                        " ".repeat(2),
                        aTokens[4],
                        " ".repeat(8),
                        aTokens[5],
                        " ".repeat(2),
                        aTokens[6],
                        " ".repeat(2),
                        aTokens[7]);
            }
            case "charge" -> {
                aTokens[1] = String.valueOf((Integer.
                        parseInt(aTokens[1]) + 100));
                tmpResult = String.join("",
                        aTokens[0], 
                        " ".repeat(9),
                        aTokens[1],
                        " ".repeat(11),
                        aTokens[2]);
            }
        }
        return tmpResult; 
    }
    
    /**
     * Fills left side of the string with spaces so the string is right aligned
     *   It is much faster version of String.format()
     * @param aInput
     *   Input string
     * @param aPadUpTo
     *   Total length of Input string and spaces
     * @return 
     *   String with left filled with spaces
     */
    public String padLeft(String aInput, int aPadUpTo) {
        String tmpSB;

        // Check parameters
        if (aInput == null || aInput.isEmpty()) {
            throw new IllegalArgumentException("aInput is null or empty.");
        } else if (aPadUpTo <= 0) {
            throw new IllegalArgumentException("aPadUpTo should be positive.");
        }

        char tmpPadChar = ' ';
        tmpSB = String.valueOf(tmpPadChar)
                .repeat(Math.max(0, aPadUpTo - aInput.length())) + aInput;
        return tmpSB;
    }
    
    /**
     * Fills right side of string with spaces so the string is left aligned
     *   Much faster version of String.format()
     * @param aInput
     *   Input string
     * @param aPadUpTo
     *   Total length of Input string and spaces
     * @return 
     *   String with right filled with spaces
     */
    public String padRight(String aInput, int aPadUpTo) {
        String tmpSb;

        if (aInput == null || aInput.isEmpty()) {
            throw new IllegalArgumentException("aInput is null or empty.");
        } else if (aPadUpTo <= 0) {
            throw new IllegalArgumentException("aPadUpTo should be positive.");
        }

        char tmpPadChar = ' ';
        tmpSb = aInput + String.valueOf(tmpPadChar)
                .repeat(Math.max(0, aPadUpTo - aInput.length()));
        return tmpSb;
    }
    
    /**
     * Workaround method to fix TIP5P due to Lone pairs
     * @param aFileName: .xyz-file name with H2O in it
     * @param aStartIndex: Number in .xyz-file
     */
    public void fixTinkerXYZ_H2O(String aFileName, int aStartIndex) {
        String line;
        String newLine;
        StringBuilder strBuilder = new StringBuilder(5000);
        int startIndex = aStartIndex;
        int index = 0;
        Path path = Paths.get(aFileName);
        
        try (BufferedReader tmpBR = Files.newBufferedReader(path, 
                StandardCharsets.UTF_8)) {
            
            while ((line = tmpBR.readLine()) != null ) {
                if (index == startIndex) {
                    newLine = line.substring(0, 8)
                            + " O"
                            + line.substring(11);
                    line = newLine;
                } else if (index == startIndex + 1 
                        || index == startIndex + 2) {
                    newLine = line.substring(0, 8)
                            + " H"
                            + line.substring(11);
                    line = newLine;
                } else if (index == startIndex + 3 
                        || index == startIndex + 4) {
                    newLine = line.substring(0, 8)
                            + "Lp"
                            + line.substring(11);
                    line = newLine;
                }
                strBuilder.append(line);
                strBuilder.append(LINESEPARATOR);
                index++;
            }
            
        }catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during reading .xyz file in fixTinerXYZ_H2O",
                    ex);
        }
        path.toFile().delete();
        try (BufferedWriter tmpBW = Files.newBufferedWriter(path,
                StandardCharsets.UTF_8)) {
            tmpBW.append(strBuilder);
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during writing .xyz file in fixTinkerXYZ_H2O."
                            , ex);
        }
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    
    private void initialize() {
        ResourceBundle tmpBundle;
        Path tmpExternalPath = Paths.get(BUNDLE_NAME_EXTERN);

        if (Files.exists(tmpExternalPath)) {
            // 1. Try: Load the external file (for the production)
            try (InputStream is = Files.newInputStream(tmpExternalPath)) {
                tmpBundle = new PropertyResourceBundle(is);
            } catch (IOException ex) {
            // Fallback during reading the file 
                tmpBundle = ResourceBundle.getBundle(BUNDLE_NAME_INTERN, 
                        Locale.getDefault(), this.getClass().getClassLoader());
            }
        } else {
            // 2. Fallback: Read internal file (for the development)
            tmpBundle = ResourceBundle.getBundle(BUNDLE_NAME_INTERN, 
                    Locale.getDefault(), this.getClass().getClassLoader());
        }
        RESOURCE_BUNDLE = tmpBundle;
        smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        atomicNumber = this.getAtomicNumberTable();
        vdWRadii = this.getVdWRadii();
    }
    
    /**
     * Split a string line and return the coordinates only
     * @param aString string line
     * @return  coordinates
     */
    private String[] split(String aString) {
        boolean tmpIsSpace;
        int tmpStart;
        int tmpEnd;
        int tmpTokenIndex;
        int tmpIndex;
        
        String[] tmpSplitedString = new String[3];
        tmpStart = 0;
        tmpTokenIndex = 0;
        tmpIndex = 0;
        while ((tmpEnd = aString.indexOf(' ', tmpStart)) >= 0) {
            if (tmpTokenIndex > 1) {
                tmpSplitedString[tmpIndex] = aString.substring(tmpStart, 
                        tmpEnd);
                tmpIndex++;
            }
            tmpTokenIndex++;
            if(tmpTokenIndex > 4) {
                break;
            }
            tmpIsSpace = true;
            tmpStart = tmpEnd + 1;
            while (tmpIsSpace) {
                if (aString.charAt(tmpStart) != ' ') {
                    tmpIsSpace = false;
                } else {
                    tmpStart++;
                }
            }
        }
        return tmpSplitedString;
    }
    
    // </editor-fold>

}