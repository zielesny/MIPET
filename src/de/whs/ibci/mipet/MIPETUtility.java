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
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.function.Consumer;
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
public class MIPETUtility {
    
    //<editor-fold defaultstate="collapsed" desc="Final class variables">
    
    /**
     * Resource bundle name
     */
    private static final String BUNDLE_NAME_INTERN = "de.whs.ibci.mipet.MIPET";
    private static final String BUNDLE_NAME_EXTERN = "MIPET.properties";
    
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(
            MIPET.class.getName());
    
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
    
    /** * Constructor
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
     * @param smilesStr input fragment as SMILES 
     * @return Van der Waals volume by using VAB Volume of CDK
     */
    public double getVdwVolume(String smilesStr) {
        if (smilesStr == null || smilesStr.isEmpty()) {
            throw new IllegalArgumentException("Smiles was null or empty.");
        }
        // Check Cache
        if (vdwVolumeCache.containsKey(smilesStr)) {
            return vdwVolumeCache.get(smilesStr);
        }
        
        double vabcVolume = 0.0;
        try {
            if (!smilesStr.equals("[Na+]")) {
                IAtomContainer particle = smilesParser.parseSmiles(smilesStr);
                AtomContainerManipulator.
                        percieveAtomTypesAndConfigureAtoms(particle);
                vabcVolume = VABCVolume.calculate(particle);
                
                // Save to Cache
                vdwVolumeCache.put(smilesStr, vabcVolume);
            }
        } catch (CDKException ex) {
            LOGGER.log(Level.SEVERE,"CDKException was thrown.", ex);
        }
        return vabcVolume;
    }
    
    /**
     * Get method for calculation of atomic mass by using CDK 
     * @param aSmilesStr Input fragment as SMILES
     * @param anIsSmiles Flag if aSmilesString is a SMILES or not (element)
     * @return
     * Atomic mass of the fragment
     */
    public double getAtomicMass(String aSmilesStr, boolean anIsSmiles) {
        if (aSmilesStr == null || aSmilesStr.trim().isEmpty()) {
            return 0;
        }
        if (aSmilesStr.trim().equalsIgnoreCase("lp")) {
            return 0;
        }
        String cacheKey = (anIsSmiles ? "S_" : "E_") + aSmilesStr.trim();
        
        // Check cache
        if (atomicMassCache.containsKey(cacheKey)) {
            return atomicMassCache.get(cacheKey);
        }
        
        double atomicMass = 0;
        try {
            IAtomContainer particle;
            if (anIsSmiles) {
                particle = smilesParser.parseSmiles(aSmilesStr);
            } else {
                particle = smilesParser.parseSmiles("[" + aSmilesStr + "]");
            }
            AtomContainerManipulator
                    .percieveAtomTypesAndConfigureAtoms(particle);
            atomicMass  = AtomContainerManipulator.getMass(particle);
            
            // Save to cache
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
        
        HashMap<String, Integer> atomicNumberTable = new HashMap<>();
        
        for (int i = 1; i < 93; i++) {
            atomicNumberTable.put(PeriodicTable.getSymbol(i), i); 
        }
        
        return atomicNumberTable;
    }
    
    /**
     * Get method for vdW-radii of atoms by using CDK
     *  (remarks: there are some atoms without data)
     * @return  vdW-radii in Angstrom, index number = atomic number
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
     * A List with double values
     * @param aValue
     * A value
     * @return 
     * True if the value is in the list (otherwise false)
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
     * @param aCoordRecord Coordinate data of solute and solvent particles
     * @param anElements1 Element names of solute particle
     * @param anElements2 Element names of solvent particles
     * @param aBoxLength Length of simulation box
     * @param aCatchRadius Catch radius
     * @return Neighbor numbers for each iteration
     * Neighbor particle number of solute of each simulation iteration
     */
    public LinkedList<int[]> getNeighborNumbersBruteForce(
            CoordinatesRecord aCoordRecord,
            String[] anElements1,
            String[] anElements2,
            double aBoxLength,
            double aCatchRadius) {
        
        final double ONEHALF = 0.5;
        
        int atomSize1 = aCoordRecord.coord1()[0].length;
        int atomSize2 = aCoordRecord.coord2()[0][0].length;
        int partSize2 = aCoordRecord.coord2()[0].length;
        HashSet<String> ignoreElements = new HashSet<>();
        HashSet<Integer> neighborPart = new HashSet<>();
        HashSet<Integer> ignoreIndex1 = new HashSet<>();
        HashSet<Integer> ignoreIndex2 = new HashSet<>();
        ignoreElements.add("M");
        ignoreElements.add("Lp");
        double[] vdWRadii1 = new double[atomSize1];
        double[] vdWRadii2 = new double[atomSize2];
        
        for (int i = 0; i < atomSize1; i++) {
            if (! ignoreElements.contains(anElements1[i])) {
                vdWRadii1[i] = vdWRadii[atomicNumber.get(anElements1[i])];
            } else {
                ignoreIndex1.add(i);
            }
        }
        
        for (int i = 0; i < atomSize2; i++) {
            if (! ignoreElements.contains(anElements2[i])) {
                vdWRadii2[i] = vdWRadii[atomicNumber.get(anElements2[i])];
            } else {
                ignoreIndex2.add(i);
            }
        }
        
        int iterSize = aCoordRecord.coord1().length;
        double[][][] coord1 = aCoordRecord.coord1();
        double[][][][] coord2 = aCoordRecord.coord2();
        LinkedList<int[]> resultList = new LinkedList<>();
        
        for (int i = 0; i < iterSize; i++) {

            for (int j = 0; j < atomSize1; j++) {
                if (!ignoreIndex1.contains(j)) {
                    
                    for (int k = 0; k < partSize2; k++) {
                        if (neighborPart.contains(k)) {
                            continue;
                        }

                        for (int l = 0; l < atomSize2; l++) {
                            if (!ignoreIndex2.contains(l)) {
                                double deltaX = coord1[i][j][0] 
                                        - coord2[i][k][l][0];
                                double deltaY = coord1[i][j][1] 
                                        - coord2[i][k][l][1];
                                double deltaZ = coord1[i][j][2] 
                                        - coord2[i][k][l][2];
                                if(deltaX > aBoxLength * ONEHALF)
                                    deltaX -= aBoxLength;
                                if(deltaX <= -aBoxLength * ONEHALF)
                                    deltaX += aBoxLength;
                                if(deltaY > aBoxLength * ONEHALF)
                                    deltaY -= aBoxLength;
                                if(deltaY <= -aBoxLength * ONEHALF)
                                    deltaY += aBoxLength;
                                if(deltaZ > aBoxLength * ONEHALF)
                                    deltaZ -= aBoxLength;
                                if(deltaZ <= -aBoxLength * ONEHALF)		
                                    deltaZ += aBoxLength;
                                double distQ = deltaX * deltaX 
                                        + deltaY * deltaY 
                                        + deltaZ * deltaZ;
                                double minDist = aCatchRadius + vdWRadii1[j] 
                                        + vdWRadii2[l];
                                double minDistQ = minDist * minDist;
                                if (distQ <= minDistQ) {
                                    neighborPart.add(k);
                                    break;
                                } 
                            }
                        }
                    
                    }
                    
                }
            }
            
            int index = 0;
            int[] neighborIndices = new int[neighborPart.size()];
            
            for (int cn : neighborPart) {
                neighborIndices[index++] = cn;
            }
            
            Arrays.sort(neighborIndices);
            resultList.add(neighborIndices);
            neighborPart.clear();
        }
        
        return resultList;
    }
    
    /**
     * Get method for neighbornumbers
     *  It uses block-wise algorithm to minimize the output file size 
     * 
     * @param block Coodinates of one block
     * @param anElements1 Element names of solut molecule
     * @param anElements2 Element names of solvent molecule
     * @param aBoxLength Box lenghth in Angstrom
     * @param aCatchRadius Catch radius in Angstrom
     * @return Neighbornumbers as an array
     */
    public int[] getNeighborNumbersBruteForceForBlock(
            CoordinatesBlock block,
            String[] anElements1,
            String[] anElements2,
            double aBoxLength,
            double aCatchRadius) {

        final double ONEHALF = 0.5;

        // Dimensions from the block
        double[][] coord1 = block.getCoord1();
        double[][][] coord2 = block.getCoord2();
        int atomSize1 = coord1.length;
        int partSize2 = coord2.length;
        int atomSize2 = coord2[0].length;
        HashSet<String> ignoreElements = new HashSet<>();
        HashSet<Integer> neighborPart = new HashSet<>();
        HashSet<Integer> ignoreIndex1 = new HashSet<>();
        HashSet<Integer> ignoreIndex2 = new HashSet<>();
        ignoreElements.add("M");
        ignoreElements.add("Lp");
        double[] vdWRadii1 = new double[atomSize1];
        double[] vdWRadii2 = new double[atomSize2];

        for (int i = 0; i < atomSize1; i++) {
            if (!ignoreElements.contains(anElements1[i])) {
                vdWRadii1[i] = vdWRadii[atomicNumber.get(anElements1[i])];
            } else {
                ignoreIndex1.add(i);
            }
        }

        for (int i = 0; i < atomSize2; i++) {
            if (!ignoreElements.contains(anElements2[i])) {
                vdWRadii2[i] = vdWRadii[atomicNumber.get(anElements2[i])];
            } else {
                ignoreIndex2.add(i);
            }
        }

        // Calculation block
        for (int j = 0; j < atomSize1; j++) {
            if (!ignoreIndex1.contains(j)) {

                for (int k = 0; k < partSize2; k++) {
                    if (neighborPart.contains(k)) {
                        continue;
                    }

                    for (int l = 0; l < atomSize2; l++) {
                        if (!ignoreIndex2.contains(l)) {
                            // Direct access via 2D and 3D array of the active block
                            double deltaX = coord1[j][0] - coord2[k][l][0];
                            double deltaY = coord1[j][1] - coord2[k][l][1];
                            double deltaZ = coord1[j][2] - coord2[k][l][2];

                            if (deltaX > aBoxLength * ONEHALF) {
                                deltaX -= aBoxLength;
                            }
                            if (deltaX <= -aBoxLength * ONEHALF) {
                                deltaX += aBoxLength;
                            }
                            if (deltaY > aBoxLength * ONEHALF) {
                                deltaY -= aBoxLength;
                            }
                            if (deltaY <= -aBoxLength * ONEHALF) {
                                deltaY += aBoxLength;
                            }
                            if (deltaZ > aBoxLength * ONEHALF) {
                                deltaZ -= aBoxLength;
                            }
                            if (deltaZ <= -aBoxLength * ONEHALF) {
                                deltaZ += aBoxLength;
                            }

                            double distQ = deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ;
                            double minDist = aCatchRadius + vdWRadii1[j] + vdWRadii2[l];
                            double minDistQ = minDist * minDist;

                            if (distQ <= minDistQ) {
                                neighborPart.add(k);
                                break;
                            }
                        }
                    }
                }
            }
        }

        // Return block
        int index = 0;
        int[] neighborIndices = new int[neighborPart.size()];
        
        for (int cn : neighborPart) {
            neighborIndices[index++] = cn;
        }

        java.util.Arrays.sort(neighborIndices);
        return neighborIndices;
    }
    
    /**
     * Get method for the neighbor particle number of solute (Using cell index method)
     * 
     * @param aCoordRecord Coordinate data of solute and solvent particles
     * @param anElements1 Element names of solute particle
     * @param anElements2 Element names of solvent particles
     * @param aBoxDist Length of simulation box
     * @param aCatchRadius Catch radius
     * @return Neighbor particle number of solute of each simulation iteration
     */
    public LinkedList<int[]> getNeighborNumbers(
            CoordinatesRecord aCoordRecord,
            String[] anElements1,
            String[] anElements2,
            double aBoxDist,
            double aCatchRadius) {
        final double ONEHALF = 0.5;
        int[] neighborIndices;
        int atomNumber1 = aCoordRecord.coord1()[0].length;
        int atomNumber2 = aCoordRecord.coord2()[0][0].length;
        int partNumber2 = aCoordRecord.coord2()[0].length;
        int iterNumber = aCoordRecord.coord1().length;
        double[][][] coord1 = aCoordRecord.coord1();
        double[][][][] coord2 = aCoordRecord.coord2();
        double[] vdWRadii1 = new double[atomNumber1];
        double[] vdWRadii2 = new double[atomNumber2];
        
        for (int i = 0; i < atomNumber1; i++) {
            vdWRadii1[i] = vdWRadii[atomicNumber.get(anElements1[i])];
        }
        
        for (int i = 0; i < atomNumber2; i++) {
            vdWRadii2[i] = vdWRadii[atomicNumber.get(anElements2[i])];
        }
        
        // Determine cell edge length and cell number
        double maxRadii1 = this.getMax(vdWRadii1);
        double maxRadii2 = this.getMax(vdWRadii2);
        double cellDist = aBoxDist / (maxRadii1 + maxRadii2 + aCatchRadius); // cell edge length in Angstrom
        int cellX;
        int cellY;
        int cellZ;
        int cellIndex;
        int cellsInRow = (int) cellDist;
        cellDist = aBoxDist / cellsInRow;
        LinkedList<int[]> resultList = new LinkedList<>();
        HashSet<Integer> neighborCells = new HashSet<>();
        HashSet<Integer> neighborPartCand = new HashSet<>();
        HashSet<Integer> neighborPart = new HashSet<>();
        
        for (int i = 0; i < iterNumber; i++) {
            cellIndex = 1;
            
            // Determine the solute particle cells
            for (int j = 0; j < atomNumber1; j++) {
                cellX = (int) ((coord1[i][j][0] + aBoxDist / 2) / cellDist) + 1;
                if(cellX > cellsInRow) {
                    cellX = cellsInRow;
                }
                cellY = (int) ((coord1[i][j][1] + aBoxDist / 2) / cellDist) + 1;
                if(cellY > cellsInRow) {
                    cellY = cellsInRow;
                }
                cellZ = (int) ((coord1[i][j][2] + aBoxDist / 2) / cellDist) + 1;
                if(cellZ > cellsInRow) {
                    cellZ = cellsInRow;
                }
                cellIndex = cellX + (cellY - 1) * cellsInRow +
                        (cellZ - 1) * cellsInRow * cellsInRow;
                neighborCells.add(cellIndex);
            }
            
            // Determine the neighbor cells
            neighborCells = this.getNeighborCellNumbers(cellsInRow, 
                    neighborCells);
            
            // Determine the neighbor particles
            for (int j = 0; j < partNumber2; j++) {
                
                for (int k = 0; k < atomNumber2; k++) {
                    cellX = (int) ((coord2[i][j][k][0] + aBoxDist / 2)
                            / cellDist) + 1;
                    if(cellX > cellsInRow) {
                        cellX = cellsInRow;
                    }
                    cellY = (int) ((coord2[i][j][k][1] + aBoxDist / 2) 
                            / cellDist) + 1;
                    if(cellY > cellsInRow) {
                        cellY = cellsInRow;
                    }
                    cellZ = (int) ((coord2[i][j][k][2] + aBoxDist / 2) 
                            / cellDist) + 1;
                    if(cellZ > cellsInRow) {
                        cellZ = cellsInRow;
                    }
                    cellIndex = cellX + (cellY - 1) * cellsInRow +
                        (cellZ - 1) * cellsInRow * cellsInRow;
                    if (neighborCells.contains(cellIndex)) {
                        if(neighborPartCand.contains(j)) {
                            break;
                        } else {
                            neighborPartCand.add(j);
                            break;
                        }
                    }
                }
            }
            
            for (int j = 0; j < atomNumber1; j++) {
                
                for (int k = 0; k < partNumber2; k++) {
                    if (neighborPartCand.contains(k)) {
                        
                        for (int l = 0; l < atomNumber2; l++) {
                            double deltaX = coord1[i][j][0] 
                                    - coord2[i][k][l][0];
                            double deltaY = coord1[i][j][1] 
                                    - coord2[i][k][l][1];
                            double deltaZ = coord1[i][j][2] 
                                    - coord2[i][k][l][2];
                            if(deltaX > aBoxDist * ONEHALF)
                                deltaX -= aBoxDist;
                            if(deltaX <= -aBoxDist * ONEHALF)
                                deltaX += aBoxDist;
                            if(deltaY > aBoxDist * ONEHALF)
                                deltaY -= aBoxDist;
                            if(deltaY <= -aBoxDist * ONEHALF)
                                deltaY += aBoxDist;
                            if(deltaZ > aBoxDist * ONEHALF)
                                deltaZ -= aBoxDist;
                            if(deltaZ <= -aBoxDist * ONEHALF)		
                                deltaZ += aBoxDist;
                            double distQ = deltaX * deltaX + deltaY * deltaY 
                                    + deltaZ * deltaZ;
                            double minDist = aCatchRadius + vdWRadii1[j] 
                                    + vdWRadii2[l];
                            double minDistQ = minDist * minDist;
                            if (distQ <= minDistQ) {
                                neighborPart.add(k);
                                break;
                            }
                        }
                        
                    }
                }
                
            }
            
            int index = 0;
            neighborIndices = new int[neighborPart.size()];
            
            for (int cn : neighborPart) {
                neighborIndices[index++] = cn;
            }
            
            Arrays.sort(neighborIndices);
            resultList.add(neighborIndices);
            neighborPart.clear();
        }
        
        return resultList;
    }
    
    /**
     * Determine cell numbers solute particle and its neighbor cells
     * @param aCellsInRow Numbers of cells per edge
     * @param aSoluteCellIndices All cell numbers of of located solute particle
     * @return Cell numbers of solute particle
     */
    public HashSet<Integer> getNeighborCellNumbers(int aCellsInRow, 
            HashSet<Integer> aSoluteCellIndices) {
        HashSet<Integer> neighborIndices = new HashSet<>();
        int[][][] indexCube = 
                new int[aCellsInRow + 2][aCellsInRow + 2][aCellsInRow + 2];
        int cellIdx = 1;
        
        for (int i = 1; i <= aCellsInRow; i++) {
            
            for (int j = 1; j <= aCellsInRow; j++) {
                
                for (int k = 1; k <= aCellsInRow; k++) {
                    indexCube[i][j][k] = cellIdx;
                    cellIdx++;
                }
                
            }
            
        }
        
        // bottom + top face
        for (int j = 1; j <= aCellsInRow; j++) {
            
            for (int k = 1; k <= aCellsInRow; k++) {
                indexCube[0][j][k] = indexCube[aCellsInRow][j][k];
                indexCube[aCellsInRow + 1][j][k] = indexCube[1][j][k];
            }
            
        }
        
        // left + right face
        for (int i = 0; i <= aCellsInRow + 1; i++) {
            
            for (int j = 1; j <= aCellsInRow; j++) {
                indexCube[i][j][0] = indexCube[i][j][aCellsInRow];
                indexCube[i][j][aCellsInRow + 1] = indexCube[i][j][1];
            }
            
        }
        
        // front + back face
        for (int i = 0; i <= aCellsInRow + 1; i++) {
            
            for (int k = 0; k <= aCellsInRow + 1; k++) {
                indexCube[i][0][k] = indexCube[i][aCellsInRow][k];
                indexCube[i][aCellsInRow + 1][k] = indexCube[i][1][k];
            }
        }
        
        for (int soluteCellIdx : aSoluteCellIndices) {
            cellIdx = 1;
            
            for (int i = 1; i <= aCellsInRow; i++) {
            
                for (int j = 1; j <= aCellsInRow; j++) {
                
                    for (int k = 1; k <= aCellsInRow; k++) {
                        if (soluteCellIdx == cellIdx) {
                            neighborIndices.add(indexCube[i][j][k]);
                            neighborIndices.add(indexCube[i][j][k + 1]);
                            neighborIndices.add(indexCube[i][j][k - 1]);
                            neighborIndices.add(indexCube[i][j + 1][k]);
                            neighborIndices
                                    .add(indexCube[i][j + 1][k + 1]);
                            neighborIndices
                                    .add(indexCube[i][j + 1][k - 1]);
                            neighborIndices.add(indexCube[i][j - 1][k]);
                            neighborIndices
                                    .add(indexCube[i][j - 1][k + 1]);
                            neighborIndices
                                    .add(indexCube[i][j - 1][k - 1]);
                            neighborIndices
                                    .add(indexCube[i + 1][j][k]);
                            neighborIndices
                                    .add(indexCube[i + 1][j][k + 1]);
                            neighborIndices
                                    .add(indexCube[i + 1][j][k - 1]);
                            neighborIndices
                                    .add(indexCube[i + 1][j + 1][k]);
                            neighborIndices
                                    .add(indexCube[i + 1][j + 1][k + 1]);
                            neighborIndices
                                    .add(indexCube[i + 1][j + 1][k - 1]);
                            neighborIndices
                                    .add(indexCube[i + 1][j - 1][k]);
                            neighborIndices
                                    .add(indexCube[i + 1][j - 1][k + 1]);
                            neighborIndices
                                    .add(indexCube[i + 1][j - 1][k - 1]);
                            neighborIndices.add(indexCube[i - 1][j][k]);
                            neighborIndices
                                    .add(indexCube[i - 1][j][k + 1]);
                            neighborIndices
                                    .add(indexCube[i - 1][j][k - 1]);
                            neighborIndices
                                    .add(indexCube[i - 1][j + 1][k]);
                            neighborIndices
                                    .add(indexCube[i - 1][j + 1][k + 1]);
                            neighborIndices
                                    .add(indexCube[i - 1][j + 1][k - 1]);
                            neighborIndices
                                    .add(indexCube[i - 1][j - 1][k]);
                            neighborIndices
                                    .add(indexCube[i - 1][j - 1][k + 1]);
                            neighborIndices
                                    .add(indexCube[i - 1][j - 1][k - 1]);
                        }
                        cellIdx++;
                    }
                
                }
            
            }
        }
        
        return neighborIndices;
    }
    
    /**
     * processArcFileBlockByBlock method
     * 
     * @param sourcePath Source path name
     * @param anAtomNumber1 Atom number of first molecule
     * @param anAtomNumber2 Atom number of second molecule
     * @param elements1 Element names of first molecule
     * @param elements2 Element names of second molecule
     * @param boxLength Box lengths in Angstrom
     * @param catchRadius Catch radius in Angsrom
     * @param tinkerProcess Tinker process
     * @param resultConsumer Result consumer
     */
    
    public void processArcFileBlockByBlock(Path sourcePath,
            int anAtomNumber1, int anAtomNumber2,
            String[] elements1, String[] elements2,
            double boxLength, double catchRadius,
            Process tinkerProcess, 
            Consumer<int[]> resultConsumer) {
        if (sourcePath == null || resultConsumer == null) {
            throw new IllegalArgumentException("Paths and Consumer must not be null.");
        }
        
        // Wait for a short time until Tinker creates the file
        int waitAttempts = 0;
        while (!Files.exists(sourcePath) && tinkerProcess.isAlive() 
                && waitAttempts < 100) {
            try {
                Thread.sleep(50);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }
            waitAttempts++;
        }
        try (BufferedReader reader = Files.newBufferedReader(sourcePath, 
                StandardCharsets.UTF_8)) {
            String firstLine = reader.readLine();
            
            // If Tinker starts extremely slow, wait a short time for the header
            while (firstLine == null && tinkerProcess.isAlive()) {
                try {
                    Thread.sleep(10);
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                }
                firstLine = reader.readLine();
            }

            if (firstLine == null) {
                return;
            }
            int atomNumber = Integer.parseInt(firstLine.substring(0, 6).trim());
            int particleNumber2 = (atomNumber - anAtomNumber1) / anAtomNumber2;

            while (true) {
                // Read next block
                CoordinatesBlock block = readNextBlock(reader, anAtomNumber1, 
                        anAtomNumber2, particleNumber2, tinkerProcess);

                if (block == null) {
                    // If there is no block but tinker is still running -> wait and retry (Live-Streaming)
                    if (tinkerProcess.isAlive()) {
                        try {
                            Thread.sleep(5); // Wait 5 ms
                        } catch (InterruptedException ex) {
                            Thread.currentThread().interrupt();
                            break;
                        }
                        continue;
                    } else {
                        // Tinker completed and there are no more blocks left
                        break;
                    }
                }

                // Calculate the neighbors for the parsed block
                int[] currNeighbors = getNeighborNumbersBruteForceForBlock(
                        block, elements1, elements2, boxLength, catchRadius);

                resultConsumer.accept(currNeighbors);
            }

        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, "Error processing ARC file block by block", ex);
        }
    }
    
    private static CoordinatesBlock readNextBlock(
            BufferedReader reader,
            int anAtomNumber1,
            int anAtomNumber2,
            int particleNumber2,
            Process tinkerProcess) throws IOException {

        String blockHeader = reader.readLine();
        if (blockHeader == null) {
            return null;
        }

        double[][] coord1 = new double[anAtomNumber1][3];
        double[][][] coord2 = new double[particleNumber2][anAtomNumber2][3];
        String[] lineArray;

        // 2. Read coordinates of the first particle (Solute)
        for (int i = 0; i < anAtomNumber1; i++) {
            String line = reader.readLine();

            while (line == null && tinkerProcess.isAlive()) {
                try {
                    Thread.sleep(5); // Give the system I/O breathing room
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                }
                line = reader.readLine();
            }

            // STRICT SHIELD: Immediate safety exit if line remains null
            if (line == null) {
                if (!tinkerProcess.isAlive()) {
                    return null; // Silent closure, Tinker terminated normally
                }
                LOGGER.log(Level.WARNING, "ARC file ended unexpectedly while Tinker is still running (Atom1).");
                return null;
            }

            lineArray = split(line.trim());
            if (lineArray == null || lineArray[0] == null) {
                return null;
            }
            coord1[i][0] = Double.parseDouble(lineArray[0]);
            coord1[i][1] = Double.parseDouble(lineArray[1]);
            coord1[i][2] = Double.parseDouble(lineArray[2]);
        }

        // 3. Read coordinates of the second particle (Solvent)
        for (int p = 0; p < particleNumber2; p++) {
            for (int a = 0; a < anAtomNumber2; a++) {
                String line = reader.readLine();

                while (line == null && tinkerProcess.isAlive()) {
                    try {
                        Thread.sleep(5);
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                    }
                    line = reader.readLine();
                }

                // STRICT SHIELD:
                if (line == null) {
                    if (!tinkerProcess.isAlive()) {
                        return null; // Silent closure, Tinker terminated normally
                    }
                    LOGGER.log(Level.WARNING, "ARC file ended unexpectedly while Tinker is still running (Atom2).");
                    return null;
                }

                lineArray = split(line.trim());
                if (lineArray == null || lineArray[0] == null) {
                    return null;
                }
                coord2[p][a][0] = Double.parseDouble(lineArray[0]);
                coord2[p][a][1] = Double.parseDouble(lineArray[1]);
                coord2[p][a][2] = Double.parseDouble(lineArray[2]);
            }
        }

        return new CoordinatesBlock(coord1, coord2);
    }
    
    /**
     * Determine (widest) molecular diameter 
     * 
     * @param aTinkerXYZ TinkerXYZ object
     * @return (Widest) molecular diameter
     */
    public double getMolecularDiameter(TinkerXYZ aTinkerXYZ) {
        int atomSize = aTinkerXYZ.getN_atom1();
        String[] elements = aTinkerXYZ.getElementList1();
        double[][] coord = aTinkerXYZ.getCoordinateList1()[0];
        double maxDistanceQ = 0.0;
        int idx1 = 0;
        int idx2 = 0;
        
        for (int i = 0; i < atomSize; i++) {
            
            for (int j = i + 1; j < atomSize; j++) {
                double deltaX = coord[i][0] - coord[j][0];
                double deltaY = coord[i][1] - coord[j][1];
                double deltaZ = coord[i][2] - coord[j][2];
                double maxDistQCandidate = deltaX * deltaX 
                        + deltaY * deltaY +deltaZ * deltaZ;
                if (maxDistQCandidate > maxDistanceQ) {
                    maxDistanceQ = maxDistQCandidate;
                    idx1 = i;
                    idx2 = j;
                }
            }
            
        }
        
        double diameter = Math.sqrt(maxDistanceQ);
        final int atomicNumber1 = PeriodicTable.getAtomicNumber(
                elements[idx1]);
        final int atomicNumber2 = PeriodicTable.getAtomicNumber(
                elements[idx2]);
        final double vdWRadii1 = this.getVdWRadii()[atomicNumber1];
        final double vdWRadii2 = this.getVdWRadii()[atomicNumber2];
        diameter += vdWRadii1 + vdWRadii2;
        return diameter;
    }
    
    /**
     * Write the last part of .arc file to .xyz file
     * @param anArcPath
     * Path name of the .arc file
     * @param aXyzPath
     * Path name of the .xyz file
     * @param aSteps 
     * Step number in the .arc file
     */
    public void writeLastPartToXYZ(Path anArcPath, Path aXyzPath, int aSteps) {
        try (BufferedReader reader = Files.newBufferedReader(anArcPath, 
                StandardCharsets.UTF_8);
                BufferedWriter writer = Files.newBufferedWriter(aXyzPath, 
                        StandardCharsets.UTF_8)) {
            int headerSize = 2;
            String line = reader.readLine();
            int atomSize = Integer.parseInt(line.substring(0, 6).trim());
            int ignoreLineSize = (atomSize + headerSize) * (aSteps - 1) - 1;
            
            for (int i = 0; i < ignoreLineSize; i++) {
                reader.readLine();
            }
            
            while((line = reader.readLine()) != null) {
                writer.append(line);
                writer.newLine();
            }
                
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException in writeLastPartToXYZ.", ex);
        }
    }
    
    /**
     * Get SMILES data method
     * Loads the SMILES code of the "known" molecules
     * @param aSmilesFileName: File name of SMILES
     * @return Molecule name and SMILES as hashmap
     */
    public HashMap<String, String> getSmilesData(String aSmilesFileName) {
        HashMap <String, String> smiles = new HashMap<>();
        Path path = Paths.get(aSmilesFileName);
        try (BufferedReader reader = Files.newBufferedReader(path,
                StandardCharsets.UTF_8)) {
            String line;
            
            while((line = reader.readLine()) != null) {
                String[] lineArray = line.trim().split("\\s+");
                smiles.put(lineArray[0], lineArray[1]);
            }
            
        } catch (FileNotFoundException ex) {
           LOGGER.log(Level.SEVERE, 
                    "FileNotFoundException in getSmilesData.", ex);
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException in getSmilesData.", ex);
        }
        return smiles;
    }
    
    /**
     * Determines if the smallest distance between two particles 
     * is less than minDistance or not. This method is used to exclude
     * two particle configurations, which are too close. This would produce
     * too high value of intermolecular energy which will be expressed 
     * by asterisks in tinker's output file.
     * 
     * @param aCoord1 Coordinates of atoms of first particle
     * @param aCoord2 Coordinates of atoms of second particle
     * @param aMinDist Minimum distance
     * @return true if the distance of two particles is too close false otherwise
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
     * @param fileName A Text file
     * @param searchStr A search string
     * @return All lines with aSearchString in it as a list.
     */
    public List<String> findList(String fileName, String searchStr) {
        Path path = Paths.get(fileName);
        try (Stream<String> lines = Files.lines(path)) {
            return lines.filter(line -> line.contains(searchStr))
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
        double sum = 0;

        for (double value : aDoubleArray) {
            sum += value;
        }
        
        return sum;
    }
    
    /**
     * Calculates the sum of the array elements.
     * 
     * @param anIntArray A double array
     * @return Sum of the array elements
     */
    public long sum(int[] anIntArray) {
        int sum = 0;

        for (int value : anIntArray) {
            sum += value;
        }
        
        return sum;
    }
    
    /**
     * Calculates the sum of element products of two arrays.
     * 
     * @param aDoubleArray1 A double array
     * @param aDoubleArray2 Another double array
     * @return Sum of element products of two arrays
     */
    public double productSum(double[] aDoubleArray1, double[] aDoubleArray2) {
        double sum = 0;
        
        for (int i = 0; i < aDoubleArray1.length; i++) {
            sum += aDoubleArray1[i] * aDoubleArray2[i];
        }
        return sum;
    }
    
    /**
     * Calculates the standard deviation of values
     * @param aValues values
     * @param aMean mean of values
     * @return standard deviation
     */
    public double standarddeviation(int[] aValues, double aMean) {
        double sum = 0;
        double returnValue;

        for (int value : aValues) {
            sum += (value - aMean) * (value - aMean);
        }
        
        returnValue = Math.sqrt(1.0 / (aValues.length - 1) * sum);
        return returnValue;
    }
    
    /**
     * Determine the biggest value from an integer array
     * @param aValues integer array
     * @return biggest value of an integer array
     */
    public int getMax(int[] aValues) {
        int max = aValues[0];
        
        for (int i = 1; i < aValues.length; i++) {
            if (max < aValues[i]) {
                max = aValues[i];
            }
        }
        
        return max;
    }
    
    /**
     * Determine the biggest value from a double array
     * @param aValues double array
     * @return biggest value of a double array
     */
    public double getMax(double[] aValues) {
        double max = aValues[0];
        
        for (int i = 1; i < aValues.length; i++) {
            if (max < aValues[i]) {
                max = aValues[i];
            }
        }
        
        return max;
    }
    
    /**
     * Determine the smallest value of an integer array
     * @param aValues integer array
     * @return smallest value of an integer array
     */
    public int getMin(int[] aValues) {
        int min = aValues[0];
        
        for (int i = 1; i < aValues.length; i++) {
            if (min > aValues[i]) {
                min = aValues[i];
            }
        }
        
        return min;
    }
    
    /**
     * Determine the distance of two points
     * @param aCoord1 xyz coordinates of point1
     * @param aCoord2 xyz coordinates of point2
     * @return distance of point1 and point2
     */
    public double getDistance(double[] aCoord1, double[] aCoord2) {
        double deltaX = aCoord2[0] - aCoord1[0];
        double deltaY = aCoord2[1] - aCoord1[1];
        double deltaZ = aCoord2[2] - aCoord1[2];
        double result = Math.sqrt(deltaX * deltaX + deltaY * deltaY 
                + deltaZ * deltaZ);
        return result;
    }
    
    /**
     * Determine the distances to shortest neighbor
     * This is only for test purpose
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
            int molNumber = sphereNodeNumber[i];
            Path path = Paths.get(sphereNodeFileName + molNumber + ".txt");
            try (BufferedReader reader = Files.newBufferedReader(path,
                    StandardCharsets.UTF_8)) {
                double[][] matrix = new double[sphereNodeNumber[i]][3];
                int atomIdx = 0;
                shortestDistQ = 1000.;
                String line;
                
                while ((line = reader.readLine()) != null) {
                    int startIndex = line.indexOf("{");
                    int endIndex = line.lastIndexOf("}");
                    String subString = line
                            .substring(startIndex + 1, endIndex);
                    String[] tokens = subString.split(",");
                    matrix[atomIdx][0] = Double.parseDouble(tokens[0].trim());
                    matrix[atomIdx][1] = Double.parseDouble(tokens[1].trim());
                    matrix[atomIdx][2] = Double.parseDouble(tokens[2].trim());
                    atomIdx++;
                }
                
                for (int j = 1; j < molNumber; j++) {
                    double distQX = matrix[j][0] - matrix[0][0];
                    double distQY = matrix[j][1] - matrix[0][1];
                    double distQZ = matrix[j][2] - matrix[0][2];
                    double distQ = distQX * distQX + distQY * distQY
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
        System.out.print("\r\u001b[K" + message);
        System.out.flush(); 
    }
    
    /**
     * Read a part of .arc file and returns the content as StringBuilder object.
     * 
     * @param aPath A Path name
     * @param aStartIdx Zero-based start line index
     * @param aEndIdx Zero-based end line index
     * @return Content of the selected part of .arc file
     */
    public static StringBuilder readPartArcFile(Path aPath, int aStartIdx, 
            int aEndIdx) {
        String line;
        int idx = 0;
        StringBuilder strBuilder = new StringBuilder();
        try (BufferedReader reader = Files.newBufferedReader(aPath, 
                StandardCharsets.UTF_8)) {
            
            while ((line = reader.readLine()) != null) {
                if (aStartIdx <= idx && idx <= aEndIdx) {
                    strBuilder.append(line);
                    strBuilder.append(LINESEPARATOR);
                }
                idx++;
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
     * @param aJobTaskRecords First Particle name, second Particle name, result directory name etc.
     * @param aLabelValues 0: Label 1: Value
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
     * @param aJobTaskRecords Job task record list
     * @param aLabel Label name
     * @param aValues Value names
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
     * @param aJobTaskRecords Parameter and particle names etc.
     * @param aCNs Coordination numbers
     * @param aTemperature Temperature
     */
    public void writeZij_Table(ArrayList<JobTaskRecord> aJobTaskRecords,
            int[][] aCNs,
            double aTemperature) {
        int jobTaskLength = aJobTaskRecords.size();
        int idx = 0;
        
        for (int i = 0; i < jobTaskLength; i++) {
            if (aJobTaskRecords.get(i).hasCNJob()) {
                String particle1 = aJobTaskRecords.get(i).particleName1(); // solute
                String particle2 = aJobTaskRecords.get(i).particleName2(); // solvent
                String particlePair = particle1 + "_" + particle2;
                String targetDirName = aJobTaskRecords.get(i).result_CN_PathName();
                Path zijPath = Paths.get(targetDirName,
                        particlePair + "_ZijTable.dat");
                try (BufferedWriter writer = Files.newBufferedWriter(zijPath, 
                        StandardCharsets.UTF_8,
                        StandardOpenOption.CREATE,
                        StandardOpenOption.APPEND)) {
                    writer.append("temperature [K]: ")
                            .append(Double.toString(aTemperature));
                    writer.newLine();
                    writer.append(Integer.toString(aCNs[idx][0]));

                    for (int j = 1; j < aCNs[idx].length; j++) {
                        writer.newLine();
                        writer.append(Integer.toString(aCNs[idx][j]));
                    }

                    writer.newLine();
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during Writing logfile.", ex);
                }
                idx++;
            }
        }
        
    }
    
    /**
     * Get Parameter name and particle pair names
     * 
     * @param aPath Path name list
     * @return 0: Parameter name
     *         1: Particle name of first particle
     *         2: Particle name of second particle
     */
    public String[][] getParameterParticleNameList(LinkedList<String> aPath) {
        if (aPath == null || aPath.isEmpty())  {
            throw new IllegalArgumentException("Illegal argument was used in getParameterParticleName().");
        }
        LinkedList<String[]> parameterParticleList = new LinkedList<>();

        for (String s : aPath) {
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
            parameterParticleList.add(parameterParticle);
            if (!parameterParticle[1].equals(parameterParticle[2])) {
                parameterParticle = new String[3];
                parameterParticle[0] = parameter;
                parameterParticle[1] = particle2;
                parameterParticle[2] = particle1;
                parameterParticleList.add(parameterParticle);
            }
        }
        
        String[][] result = parameterParticleList.toArray(
                new String[parameterParticleList.size()][3]);
        return result;
    }
    
    /**
     * Save .key file
     * 
     * @param keyPath Key path name
     * @param aContent Content of .key file
     */
    public void writeKeyFile(Path keyPath, String aContent) {
        try {
            Files.writeString(keyPath, aContent);
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
            int energySize;
            
            for (double[] aEnergySorted1 : aEnergySorted) {
                energySize = aEnergySorted1.length;
                
                for (int j = 0; j < energySize; j++) {
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
        String result = null;
        switch (aTokens[0]) {
            case "atom" -> {
                aTokens[1] = String.valueOf((Integer.
                        parseInt(aTokens[1]) + 100));
                aTokens[2] = aTokens[1];
                result = String.join("",
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
                result = String.join("",
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
                result = String.join("",
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
                result = String.join("",
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
                result = String.join("",
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
                result = String.join("",
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
                result = String.join("",
                        aTokens[0], 
                        " ".repeat(9),
                        aTokens[1],
                        " ".repeat(11),
                        aTokens[2]);
            }
        }
        return result; 
    }
    
    /**
     * Fills left side of the string with spaces so the string is right aligned
     * It is much faster version of String.format()
     * @param aInput
     * Input string
     * @param aPadUpTo
     * Total length of Input string and spaces
     * @return 
     * String with left filled with spaces
     */
    public String padLeft(String aInput, int aPadUpTo) {
        if (aInput == null || aInput.isEmpty()) {
            throw new IllegalArgumentException("aInput is null or empty.");
        } else if (aPadUpTo <= 0) {
            throw new IllegalArgumentException("aPadUpTo should be positive.");
        }

        char padChar = ' ';
        String resultStr = String.valueOf(padChar)
                .repeat(Math.max(0, aPadUpTo - aInput.length())) + aInput;
        return resultStr;
    }
    
    /**
     * Fills right side of string with spaces so the string is left aligned
     * Much faster version of String.format()
     * @param aInput
     * Input string
     * @param aPadUpTo
     * Total length of Input string and spaces
     * @return 
     * String with right filled with spaces
     */
    public String padRight(String aInput, int aPadUpTo) {
        if (aInput == null || aInput.isEmpty()) {
            throw new IllegalArgumentException("aInput is null or empty.");
        } else if (aPadUpTo <= 0) {
            throw new IllegalArgumentException("aPadUpTo should be positive.");
        }

        char padChar = ' ';
        String resultStr = aInput + String.valueOf(padChar)
                .repeat(Math.max(0, aPadUpTo - aInput.length()));
        return resultStr;
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
        
        try (BufferedReader reader = Files.newBufferedReader(path, 
                StandardCharsets.UTF_8)) {
            
            while ((line = reader.readLine()) != null ) {
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
        try (BufferedWriter writer = Files.newBufferedWriter(path,
                StandardCharsets.UTF_8)) {
            writer.append(strBuilder);
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during writing .xyz file in fixTinkerXYZ_H2O."
                            , ex);
        }
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    private void initialize() {
        ResourceBundle bundle = null;
        
        // Determine the real terminal path
        String realWorkingDir = System.getenv("PWD");
        Path externalPath;
        if (realWorkingDir != null) {
            externalPath = Paths.get(realWorkingDir, BUNDLE_NAME_EXTERN);
        } else {
            // Fallback for Windows / NetBeans
            externalPath = Paths.get(BUNDLE_NAME_EXTERN);
        }
        if (Files.exists(externalPath)) {
            try (InputStream is = Files.newInputStream(externalPath)) {
                bundle = new PropertyResourceBundle(is);
            } catch (IOException ex) {
                // Fallback see bottom
            }
        }

        // Fallback: There is no external file or it couldn't read
        // Then: Read standard bundle.
        if (bundle == null) {
            bundle = ResourceBundle.getBundle(BUNDLE_NAME_INTERN,
                    Locale.getDefault(), this.getClass().getClassLoader());
        }
        RESOURCE_BUNDLE = bundle;
        smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        atomicNumber = this.getAtomicNumberTable();
        vdWRadii = this.getVdWRadii();
    }
    
    /**
     * Split a string line and return the coordinates only
     * @param aString string line
     * @return  coordinates
     */
    private static String[] split(String aString) {
        int endIdx;
        String[] splitedStrs = new String[3];
        int startIdx = 0;
        int tokenIdx = 0;
        int idx = 0;
        while ((endIdx = aString.indexOf(' ', startIdx)) >= 0) {
            if (tokenIdx > 1) {
                splitedStrs[idx] = aString.substring(startIdx, endIdx);
                idx++;
            }
            tokenIdx++;
            if(tokenIdx > 4) {
                break;
            }
            boolean isSpace = true;
            startIdx = endIdx + 1;
            while (isSpace) {
                if (aString.charAt(startIdx) != ' ') {
                    isSpace = false;
                } else {
                    startIdx++;
                }
            }
        }
        return splitedStrs;
    }
    // </editor-fold>
}