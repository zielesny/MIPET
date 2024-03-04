/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package de.whs.ibci.mipet;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.MissingResourceException;
import java.util.PropertyResourceBundle;
import java.util.ResourceBundle;
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
    private final static Logger LOGGER = Logger
            .getLogger(MIPET.class.getName());
    
     /**
     * Line separator
     */
    final private String LINESEPARATOR = System.getProperty("line.separator");
    
    /**
     * File separator
     */
    final private String FILESEPARATOR = System.getProperty("file.separator");
    
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
     *   Van der Waals volume by using VABVolume of CDK
     */
    public double getVdwVolume(String smilesString) {
        if (smilesString == null || smilesString.isEmpty()) {
            throw new IllegalArgumentException("Smiles was null or empty.");
        }
        
        double vabcVolume = 0.0;
        try {
            IAtomContainer particle =smilesParser.parseSmiles(smilesString);
            AtomContainerManipulator.
                    percieveAtomTypesAndConfigureAtoms(particle);
            vabcVolume = VABCVolume.calculate(particle);
        } catch (CDKException ex) {
            LOGGER.log(Level.SEVERE, 
                    "CDKException was thrown.", ex);
        }
        
        return vabcVolume;
    }
    
    /**
     * Get method for calculation of atomic mass by using CDK 
     * @param smilesString
     *   input fragment as SMILES
     * @return
     *   Atomic mass of the fragment
     */
    public double getAtomicMass(String smilesString) {
        double atomicMass = 0;
        try {
            IAtomContainer particle = smilesParser.parseSmiles(smilesString);
            AtomContainerManipulator
                    .percieveAtomTypesAndConfigureAtoms(particle);
            atomicMass  = AtomContainerManipulator.getMass(particle);
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
        double[] tmpVdWRadii = new double[93];
        
        for (int i = 1; i < 93; i++) {
            if (PeriodicTable
                    .getVdwRadius(PeriodicTable.getSymbol(i)) != null) {
                tmpVdWRadii[i] = PeriodicTable
                    .getVdwRadius(PeriodicTable.getSymbol(i));
            }
        }
        
        return tmpVdWRadii;
    }
    
    /**
     * Read the boxlength from a given .arc file
     * @param aFileName Filename of .arc file
     * @return Boxlength in Angstrom
     */
    public double getBoxLength(String aFileName) {
        String tmpReadLine;
        String[] tmpReadLineArray;
        double tmpBoxLength = 0.0;
        
        try (BufferedReader tmpBR = new BufferedReader(
                new FileReader(aFileName))) {
            tmpBR.readLine();
            tmpReadLine = tmpBR.readLine();
            if (tmpReadLine != null) {
                tmpReadLineArray = tmpReadLine.trim().split("\\s+");
                tmpBoxLength = Double.parseDouble(tmpReadLineArray[0]);
            }
        } catch (FileNotFoundException ex) {
            throw new IllegalArgumentException("File not found in getBoxLength.");
        } catch (IOException ex) {
            throw new IllegalArgumentException("IOException in getBoxLength.");            
        }
        return tmpBoxLength;
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
    public boolean contains(LinkedList<Double> aList, double aValue) {
        double tmpThreshold = 0.001;
        
        for (double tmpList : aList) {
            if (Math.abs(tmpList - aValue) < tmpThreshold) {
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
            tmpNeighborPart.clear();
            
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
     * @param aFileName File name
     * @param anAtomNumber1 Atom number of solute particle
     * @param anAtomNumber2 Atom number of slvent particle
     * @return All coordinates of .arc file
     */
    public CoordinatesRecord getCoordinatesFromArcFile(String aFileName,
            int anAtomNumber1, int anAtomNumber2) {
        // Check parameters
        if (aFileName == null || aFileName.isEmpty()) {
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
        try (BufferedReader tmpBR = new BufferedReader(
                new FileReader(aFileName), 65536)) {
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
            
            tmpBR.close();
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
     * Write the last part of .arc file to .xyz file
     * @param anArcFileName
     *   File name of the .arc file
     * @param aXyzFileName
     *   File name of the .xyz file
     * @param aSteps 
     *   Step number in the .arc file
     */
    public void writeLastPartToXYZ(String anArcFileName, String aXyzFileName,
            int aSteps) {
        int tmpAtomNumber;
        int tmpHeader;
        int tmpIgnoreLineNumber;
        String tmpLine;

        tmpHeader = 2;
        
        try (BufferedReader tmpBR = new BufferedReader(new FileReader(
                anArcFileName), 65536);
                BufferedWriter tmpBW = new BufferedWriter(new FileWriter(
                        aXyzFileName))) {
            tmpLine = tmpBR.readLine();
            tmpAtomNumber = Integer.parseInt(tmpLine.substring(0, 6).trim());
            tmpIgnoreLineNumber = (tmpAtomNumber + tmpHeader) 
                    * (aSteps - 1) - 1;
            
            for (int i = 0; i < tmpIgnoreLineNumber; i++) {
                tmpBR.readLine();
            }
            
            while((tmpLine = tmpBR.readLine()) != null) {
                tmpBW.append(tmpLine);
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
     * @param aSmilesFileName
     * @return 
     */
    public HashMap<String, String> getSmilesData(String aSmilesFileName) {
        HashMap <String, String> tmpSmiles = new HashMap<>();
        String tmpReadLine;
        String[] tmpReadLineArray;
        
        try (BufferedReader tmpBR = new BufferedReader(
                new FileReader(aSmilesFileName), 65536)) {
            
            while((tmpReadLine = tmpBR.readLine()) != null) {
                tmpReadLineArray = tmpReadLine.trim().split("\\s+");
                tmpSmiles.put(tmpReadLineArray[0], tmpReadLineArray[1]);
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
     *   two particle configurations, wich are too close. This would produce
     *   too high value of intermolecular energy which will be expressed 
     *   by asterisks in tinker's outputfile.
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
        
        double tmpXDelta;
        double tmpYDelta;
        double tmpZDelta;
        double tmpDistanceQ;
        double tmpMinDistQ;
        boolean tmpIsTooClose = false;
        
        tmpMinDistQ = aMinDist * aMinDist;
        
        for (double[] coord11 : aCoord1) {
            for (double[] coord21 : aCoord2) {
                tmpXDelta = coord11[0] - coord21[0];
                tmpYDelta = coord11[1] - coord21[1];
                tmpZDelta = coord11[2] - coord21[2];
                tmpDistanceQ = tmpXDelta * tmpXDelta 
                        + tmpYDelta * tmpYDelta 
                        + tmpZDelta * tmpZDelta;
                if(tmpDistanceQ < tmpMinDistQ){
                    tmpIsTooClose = true;
                    break;
                }
            }
        }
        return tmpIsTooClose;
    }
    
    /**
     * Returns all Lines with aSearchString
     * @param aFileName
     *   A Text file
     * @param aSearchString
     *   A search string
     * @return 
     *   All lines with aSearchString in it as a list.
     */
    public List<String> findList(String aFileName, String aSearchString) {
        int tmpLineNum;
        Path tmpFile = Paths.get(aFileName);
        List<String> tmpList = new LinkedList<>();
        String tmpLine;
            
        tmpLineNum = 0;
        try{
            Scanner tmpScanner = new Scanner(tmpFile);
            
            while(tmpScanner.hasNextLine()){
                tmpLine = tmpScanner.nextLine();
                tmpLineNum++;
                if(tmpLine.contains(aSearchString)){
                   tmpList.add(tmpLine);
                }
            }
            
        } catch (FileNotFoundException ex) {
            LOGGER.log(Level.SEVERE, 
                    "FileNotFoundException in findList.", ex);
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException in findList.", ex);
        }
        return tmpList;
    }
    
    /**
     * Calculates the sum of the array elements.
     * 
     * @param aDoubleArray A double array
     * @return Sum of the array elements
     */
    public double sum(double[] aDoubleArray) {
        double tmpSum = 0;
        
        for (int i = 0; i < aDoubleArray.length; i++) {
            tmpSum += aDoubleArray[i];
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
        
        for (int i = 0; i < anIntArray.length; i++) {
            tmpSum += anIntArray[i];
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
        
        for (int i = 0; i < aValues.length; i++) {
            tmpSum += (aValues[i] - aMean) * (aValues[i] - aMean);
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
     * Read a part of .arc file and returns the content as StringBuilder object.
     * 
     * @param aFileName
     *   A File name
     * @param aStartIndex
     *   Zero-based start line index
     * @param aEndIndex
     *   Zero-based end line index
     * @return 
     * Content of the selected part of .arc file
     */
    public StringBuilder readPartArcFile(String aFileName, int aStartIndex, 
            int aEndIndex) {
        String tmpLine;
        int tmpIndex = 0;
        StringBuilder tmpSB = new StringBuilder();
        try (BufferedReader tmpBR = new BufferedReader(
                new FileReader(aFileName), 65536)) {
            
            while ((tmpLine = tmpBR.readLine()) != null) {
                if (aStartIndex <= tmpIndex && tmpIndex <= aEndIndex) {
                    tmpSB.append(tmpLine);
                    tmpSB.append(LINESEPARATOR);
                }
                tmpIndex++;
            }
            
        } catch (IOException ex) {
                LOGGER.log(Level.SEVERE,
                        "IOException during reading" + aFileName , ex);
        }
        return tmpSB;
    }
    
    /**
     * Write Particlelogfiles
     * 
     * @param aJobTaskRecords
     *   First Particlename, second Particlename, result directory name etc.
     * @param aLabelValues
     *   0: Label 1: Value
     */
    public void writeParticleLog(ArrayList<JobTaskRecord> aJobTaskRecords,
            String[][] aLabelValues) {
        int tmpJobTaskLength;
        String tmpParticle1;
        String tmpParticle2;
        String tmpParticlePair;
        String tmpParticleLogFileName;
        String tmpResultPathName;
        
        tmpJobTaskLength = aJobTaskRecords.size();
        
        for (int i = 0; i < tmpJobTaskLength; i++) {
            if (aJobTaskRecords.get(i).hasCNJob()) {
                tmpParticle1 = aJobTaskRecords.get(i).particleName1(); // solute
                tmpParticle2 = aJobTaskRecords.get(i).particleName2(); // solvent
                tmpParticlePair = tmpParticle1 + "_" + tmpParticle2;
                tmpResultPathName = aJobTaskRecords.get(i).result_CN_PathName();
                tmpParticleLogFileName = tmpResultPathName
                        + FILESEPARATOR
                        + tmpParticlePair
                        + "_log.txt";

                try (BufferedWriter tmpBW = new BufferedWriter(new FileWriter(
                        tmpParticleLogFileName, true))) {

                    for (String[] aLabelValue1 : aLabelValues) {
                        tmpBW.append(aLabelValue1[0]);
                        tmpBW.append(aLabelValue1[1]);
                        tmpBW.append(LINESEPARATOR);
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
        
        int tmpJobTaskLength;
        int tmpIndex;
        String tmpParticle1;
        String tmpParticle2;
        String tmpParticlePair;
        String tmpParticleLogFileName;
        String tmpResultPathName;
        
        tmpJobTaskLength = aJobTaskRecords.size();
        tmpIndex = 0;
        
        for (int i = 0; i < tmpJobTaskLength; i++) {
            if (aJobTaskRecords.get(i).hasCNJob()) {
                tmpParticle1 = aJobTaskRecords.get(i).particleName1(); // solute
                tmpParticle2 = aJobTaskRecords.get(i).particleName2(); // solvent
                tmpParticlePair = tmpParticle1 + "_" + tmpParticle2;
                tmpResultPathName = aJobTaskRecords.get(i).result_CN_PathName();
                tmpParticleLogFileName = tmpResultPathName
                        + FILESEPARATOR
                        + tmpParticlePair
                        + "_log.txt";

                try (BufferedWriter tmpBW = new BufferedWriter(new FileWriter(
                        tmpParticleLogFileName, true))) {
                    tmpBW.append(aLabel);
                    tmpBW.append("(");
                    tmpBW.append(tmpParticle1);
                    tmpBW.append("/");
                    tmpBW.append(tmpParticle2);
                    tmpBW.append("): ");
                    tmpBW.append(aValues[tmpIndex]);
                    tmpBW.append(LINESEPARATOR);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during Writing logfile.", ex);
                }
                tmpIndex++;
            }
        }
        
    }
    
    /**
     * Write Zij datas
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
            int aTemperature) {
        
        int tmpJobTaskLength;
        int tmpIndex;
        String tmpParticle1;
        String tmpParticle2;
        String tmpParticlePair;
        String tmpZij_TableFileName;
        String tmpTargetDirName;
        
        tmpJobTaskLength = aJobTaskRecords.size();
        tmpIndex = 0;
        
        for (int i = 0; i < tmpJobTaskLength; i++) {
            if (aJobTaskRecords.get(i).hasCNJob()) {
                tmpParticle1 = aJobTaskRecords.get(i).particleName1(); // solute
                tmpParticle2 = aJobTaskRecords.get(i).particleName2(); // solvent
                tmpParticlePair = tmpParticle1 + "_" + tmpParticle2;
                tmpTargetDirName = aJobTaskRecords.get(i).result_CN_PathName();
                tmpZij_TableFileName = tmpTargetDirName
                        + FILESEPARATOR
                        + tmpParticlePair
                        + "_ZijTable.dat";
                try (BufferedWriter tmpBW = new BufferedWriter(new FileWriter(
                        tmpZij_TableFileName, true))) {
                    tmpBW.append("temperature [K]: ");
                    tmpBW.append(Integer.toString(aTemperature));
                    tmpBW.append(LINESEPARATOR);
                    tmpBW.append(Integer.toString(aCNs[tmpIndex][0]));

                    for (int j = 1; j < aCNs[tmpIndex].length; j++) {
                        tmpBW.append(LINESEPARATOR);
                        tmpBW.append(Integer.toString(aCNs[tmpIndex][j]));
                    }

                    tmpBW.append(LINESEPARATOR);
                } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE, 
                    "IOException during Writing logfile.", ex);
                }
                tmpIndex++;
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
        
        int tmpParamStartIndex;
        int tmpParticleStartIndex;
        int tmpSeparateIndex;
        int tmpPathNameLength = aPathName.size();
        LinkedList<String[]> tmpParameterParticleList = new LinkedList<>();
        String tmpParameter;
        String tmpParticle1;
        String tmpParticle2;
        String[] tmpParameterParticle;
        String[][] tmpReturn;
        
        for (int i = 0; i < tmpPathNameLength; i++) {
            tmpParamStartIndex = -1;
            tmpParticleStartIndex = -1;
            tmpSeparateIndex = -1;
            
            
            for (int j = aPathName.get(i).length() - 1; j >= 0; j--) {
                if (aPathName.get(i).charAt(j) == '_') {
                    tmpSeparateIndex = j;
                } else if (aPathName.get(i).charAt(j) == 
                        FILESEPARATOR.charAt(0) && 
                        tmpParticleStartIndex == -1) {
                                tmpParticleStartIndex = j + 1;
                } else if (tmpParamStartIndex == -1) {
                    if (aPathName.get(i).charAt(j) == FILESEPARATOR.charAt(0)) {
                        tmpParamStartIndex = j + 1;
                        break;
                    } else if (j == 0) {
                        tmpParamStartIndex = 0;
                    }
                }
            }
            
            tmpParameterParticle = new String[3];
            tmpParameter = aPathName.get(i)
                    .substring(tmpParamStartIndex, tmpParticleStartIndex - 1);
            tmpParameterParticle[0] = tmpParameter;
            tmpParticle1 = aPathName.get(i)
                    .substring(tmpParticleStartIndex, tmpSeparateIndex);
            tmpParameterParticle[1] = tmpParticle1; 
            tmpParticle2 = aPathName.get(i)
                    .substring(tmpSeparateIndex + 1, aPathName.get(i).length());
            tmpParameterParticle[2] = tmpParticle2;
            tmpParameterParticleList.add(tmpParameterParticle);
            if(!tmpParameterParticle[1].equals(tmpParameterParticle[2])) {
                tmpParameterParticle = new String[3];
                tmpParameterParticle[0] = tmpParameter;
                tmpParameterParticle[1] = tmpParticle2;
                tmpParameterParticle[2] = tmpParticle1;
                tmpParameterParticleList.add(tmpParameterParticle);
            }
        }
        
        tmpReturn = tmpParameterParticleList
                .toArray(new String[tmpParameterParticleList.size()][3]);
        return tmpReturn;
    }
    
    /**
     * Save .key file
     * @param aForcefieldDir
     *   Directory name of .prm files
     * @param aForcefield
     *   Force field name
     * @param aParticle
     *   Particle name
     * @param aKeyFileName
     *   Key file name
     * @param aContent
     *   Content of .key file
     */
    public void saveKeyFile(String aForcefieldDir, 
            String aForcefield,
            String aParticle,
            String aKeyFileName, 
            String aContent) {
        this.saveKeyFile(aForcefieldDir, 
                aForcefield,
                aParticle, 
                "", 
                aKeyFileName, 
                aContent);
    }
    
    /**
     * Save .key file
     * @param aForcefieldDir
     *   Directory name of .prm files
     * @param aForcefield
     *   Force field name
     * @param aParticle1
     *   Particle name of first particle
     * @param aParticle2
     *   Particle name of second particle
     * @param aKeyFileName
     *   Key file name
     * @param aContent
     *   Content of .key file
     */
    public void saveKeyFile(String aForcefieldDir, 
            String aForcefield,
            String aParticle1,
            String aParticle2,
            String aKeyFileName, 
            String aContent) {
        
        String tmpExtraPrmName;
        String tmpKeyFileString;
        String tmpParticle;
        int tmpRoundtripSize;
        
        if (aParticle2.isEmpty()) {
            tmpRoundtripSize = 1;
        } else {
            tmpRoundtripSize = 2;
        }
        tmpParticle = aParticle1;
        tmpKeyFileString = aContent;
        
        for (int i = 0; i < tmpRoundtripSize; i++) {
            if (i == 1) {
                tmpParticle = aParticle2;
            } 
            tmpExtraPrmName = aForcefieldDir
                    + FILESEPARATOR
                    + aForcefield
                    + FILESEPARATOR
                    + tmpParticle
                    + ".prm";
            // Check if there is a extra .prm file
            try {
                if (Files.exists(Paths.get(tmpExtraPrmName))) {
                    tmpKeyFileString += LINESEPARATOR;
                    tmpKeyFileString += Files.readString(Path
                            .of(tmpExtraPrmName));
                }
            } catch (IOException ex) {
                    LOGGER.log(Level.SEVERE,
                        "IOException during read extra .prm file.", ex);
            }
        }
        
        // Write .key file
        try (BufferedWriter tmpBW = new BufferedWriter(
                new FileWriter(aKeyFileName))) {
            tmpBW.append(tmpKeyFileString);
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during writing .key file.", ex);
        }
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    
    private void initialize() {
        // For the development
        RESOURCE_BUNDLE = ResourceBundle.getBundle(BUNDLE_NAME_INTERN, 
                Locale.getDefault(), this.getClass().getClassLoader());
        // For the distribution
//        try {
//            RESOURCE_BUNDLE = new PropertyResourceBundle(Files
//                    .newInputStream(Paths.get(BUNDLE_NAME_EXTERN)));
//        } catch (IOException ex) {
//            LOGGER.log(Level.SEVERE, 
//                    "IOException during initialize().", ex);
//        }
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