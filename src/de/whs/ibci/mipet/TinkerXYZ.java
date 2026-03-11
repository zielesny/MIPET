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
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.StringReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.ArrayUtils;

/**
 * Class TinkerXYZ
 * @author Mirco Daniel
 */
public class TinkerXYZ implements Cloneable {
    
    //<editor-fold defaultstate="collapsed" desc="Final class variables">

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
     * Standard capacity for stringbuilder object
     */
    final int STRINGBUILDER_CAPACITY = 4096;
    
    /**
     * Standard BufferdReader buffersize 
     */
    final int READER_BUFFERSIZE = 65536;;
    
    /**
     * Line separator
     */
    final private String LINESEPARATOR = System.lineSeparator();
    
    /**
     * Atomic mass of elements
     */
    final private HashMap<String, Double> atomicMassMap = new HashMap<>();

    //</editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    
    /**
     * Name of first particle
     */
    private String particleName1;
    
    /**
     * Name of second particle
     */
    private String particleName2;
    
    /**
     * Name of force field
     */
    private String forcefieldName;
    
    /**
     * Atom number of all particles
     */
    private int atomNumber;
    
    /**
     * Header (String behind atomnumber)
     */
    private String header;
    
    /**
     * Atom number of first particle
     */
    private int N_atom1;
    
    /**
     * Atom number of second particle
     */
    private int N_atom2;
    
    /**
     * Second particle number
     */
    private int N_particle2;
    
    /**
     * Comment line
     */
    private String comment;
    
    /**
     * Element names of first particle
     */
    private String[] elementList1;
    
    /**
     * Element names of second particle
     */
    private String[] elementList2;
    
    /**
     * Atomic masses of first particle
     */
    private double[] atomicMassList1;
    
    /**
     * Atomic masses of second particle
     */
    private double[] atomicMassList2;
    
    /**
     * Atom Coordinates of first particle
     * [i, j, k] i: sim.iteration 
     *           j: atom id 
     *           k: xyz
     */
    private double[][][] coordinateList1;
    
    /**
     * Atom Coordinates of second particle
     * [i, j, k, L] i: sim.iteration j: particle k: atom, L: xyz
     */
    private double[][][][] coordinateList2;
    
    /**
     * Parameters of specific parameterset defined in *.prm of first particle
     */
    private int[] atomTypeList1;
    
    /**
     * Parameters of specific parameterset defined in *.prm of second particle
     */
    private int[] atomTypeList2;
    
    /**
     * Connectionlist of first particle 
     */
    private int[][] connectionList1;
    
    /**
     * Connectionlist of second particle 
     */
    private int[][] connectionList2;
    
    /**
     * Distances of atom1 to atom2
     *  i: atom index of particle1
     *  j: atom index of particle2
     */
    private double[][] distances;
    
    /**
     * File content
     */
    private StringBuilder fileContent;
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    /**
     * Constructor TinkerXYZ
     */
    public TinkerXYZ() {
    }
    
    /**
     * Constructor TinkerXYZ
     * 
     * @param aForcefieldName Name of the force field
     * @param aParticleName Name of the particle
     * @param aTxyzFileName Tinker xyz filename or Tinker xyz as a string
     */
    public TinkerXYZ(String aForcefieldName, 
            String aParticleName, 
            String aTxyzFileName) {
        // Check parameters
        if(aTxyzFileName == null || aTxyzFileName.isEmpty()) {
            throw new IllegalArgumentException("Null was passed to the " 
                    + "constructor TinkerXYZ.");
        }
        
        String tmpTxyz;
        Path tmpPath;
        
        if (new File(aTxyzFileName).isFile()) {
            tmpPath = Path.of(aTxyzFileName);
            try {
                tmpTxyz = Files.readString(tmpPath);
            } catch (IOException ex) {
                throw new IllegalArgumentException("IOException in TinkerXYZ.");
            }
            this.initialize(tmpTxyz, 1, 0, 0);
        } else {
            this.initialize(aTxyzFileName, 1, 0, 0);
        }
        this.forcefieldName = aForcefieldName;
        this.particleName1 = aParticleName;
    }
    
    /**
     * Constructor TinkerXYZ
     * 
     * @param aTxyzContent Tinker xyz filename
     * @param atomSize Atom number of particle
     */
    public TinkerXYZ(String aTxyzContent, int atomSize) {
        this.initialize(aTxyzContent, 1, atomSize, 0);
    }
    
    /**
     * Constructor TinkerXYZ for getCoordinationNumber only
     * 
     * @param aTxyzFileName Tinker xyz filename
     * @param anIterationSize Number of iterations
     * @param anAtomSize1 Number of atoms from first fragment
     * @param anAtomSize2 Number of atoms from second fragment
     */
    public TinkerXYZ(String aTxyzFileName, int anIterationSize, 
            int anAtomSize1, int anAtomSize2) {
        String tmpTxyz;
        Path tmpPath;
        
        tmpPath = Path.of(aTxyzFileName);
        try {
            tmpTxyz = Files.readString(tmpPath);
        } catch (IOException ex) {
            throw new IllegalArgumentException("IOException in TinkerXYZ.");
        }
        this.initialize(tmpTxyz, anIterationSize, anAtomSize1, anAtomSize2);
    }
    /**
     * Constructor TinkerXYZ
     * 
     * @param aTinkerXYZ1 TinkerXYZ object of first fragment
     * @param aTinkerXYZ2 TinkerXYZ object of second fragment
     * @param aTinkerOn Flag for whether tinker is used or not
     */
    public TinkerXYZ(TinkerXYZ aTinkerXYZ1, TinkerXYZ aTinkerXYZ2, 
            boolean aTinkerOn) {
        this.initialize2(aTinkerXYZ1, aTinkerXYZ2, aTinkerOn);
        
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Public properties">
    /**
     * Returns the name of first particle
     * @return name of first particle
     */
    public String getParticleName1() {
        return this.particleName1;
    }
    
    /**
     * Returns the name of second particle
     * @return name of second particle
     */
    public String getParticleName2() {
        return this.particleName2;
    }
    
    /**
     * Returns the forcefield name
     * 
     * @return forcefield name
     */
    public String getForcefieldName() {
        return this.forcefieldName;
    }
    
    /**
     * Returns the atom number
     * 
     * @return the atom number
     */
    public int getN_atom() {
        return this.atomNumber;
    }
    
    /**
     * Returns the header (String behind atomnumber in the 1st line)
     * 
     * @return Header
     */
    public String getHeader() {
        return this.header;
    }
    
    /**
     * Returns the number of atoms in the first molecule
     * 
     * @return Atom number of first particle
     */
    public int getN_atom1() {
        return this.N_atom1;
    }
    
    /**
     * Returns the number of atoms in the second molecule
     * 
     * @return Atom number of second particle
     */
    public int getN_atom2() {
        return this.N_atom2;
    }
    
    /**
     * Returns the particle number of the second particle
     * 
     * @return Number of solvent particles in a simulation box
     */
    public int getN_particle2() {
        return this.N_particle2;
    }
    
    /**
     * Returns the comment line
     * 
     * @return comment
     */
    public String getComment() {
        return this.comment;
    }
            
    /**
     * Returns element names of first particle
     * 
     * @return Elementlist of first praticle
     */
    public String[] getElementList1() {
        if (this.elementList1 != null) {
            return this.elementList1;
        } else {
            throw new NullPointerException("Return value of getElementList1"
                    + "is null.");
        }
    }
    
    /**
     * Returns element names of second particle
     * 
     * @return Elementlist of second praticle
     */
    public String[] getElementList2() {
        if (this.elementList2 != null) {
            return this.elementList2;
        } else {
            throw new NullPointerException("Return value of getElementList2"
                    + "is null.");
        }
    }
    
    /** 
     * Returns atomic masses of the first particle
     * 
     * @return Atomic masses of first particle
     */
    public double[] getAtomicMassList1() {
        HashSet<String> tmpIgnoreElements = new HashSet<>();
        
        tmpIgnoreElements.add("M");
        tmpIgnoreElements.add("Lp");
        
        if (this.elementList1 != null) {
            this.atomicMassList1 = new double[this.N_atom1];
            
            for (int i = 0; i < this.N_atom1; i++) {
                if (!tmpIgnoreElements.contains(this.elementList1[i])) {
                    this.atomicMassList1[i] = this.atomicMassMap
                            .get(this.elementList1[i]);
                }
            }
            
            return this.atomicMassList1;
        } else {
            throw new NullPointerException("Return value of getElementList1 is null.");
        }
    }
    
    /** 
     * Returns atomic masses of the second particle
     * 
     * @return Atomic masses of second particle
     */
    public double[] getAtomicMassList2() {
        HashSet<String> tmpIgnoreElements = new HashSet<>();
        
        tmpIgnoreElements.add("M");
        tmpIgnoreElements.add("Lp");
        
        if (this.elementList2 != null) {
            this.atomicMassList2 = new double[this.N_atom2];
            
            for (int i = 0; i < this.N_atom2; i++) {
                if (!tmpIgnoreElements.contains(this.elementList2[i])) {
                    this.atomicMassList2[i] = this.atomicMassMap
                            .get(this.elementList2[i]);                    
                }
            }
            
            return this.atomicMassList2;
        } else {
            throw new NullPointerException("Return value of getElementList2 is null.");
        }
    }
    
    /**
     * Returns the atom coordinates of first particle
     * 
     * @return coordinates of the atoms[i, j, k]
     * i: sim.iteration j: atom k: xyz
     */
    public double[][][] getCoordinateList1() {
        if (this.coordinateList1 != null) {
            return this.coordinateList1;
        } else {
            throw new NullPointerException("Return value of getCoordinateList1"
                    + "is null.");
        }
    }
    
    /**
     * Returns the atom coordinates of second particle
     * @return coordinates[i, j, k, l] 
     * i: sim.iteration j: particle k: atom l: xyz
     */
    public double[][][][] getCoordinateList2() {
         if (this.coordinateList1 != null) {
            return this.coordinateList2;
        } else {
            throw new NullPointerException("Return value of getCoordinateList2"
                    + "is null.");
        }
    }
    
    /** 
     * Returns the coordinate of centre of mass from first particle of each 
     *  simulation
     * 
     * @return coordinate of center of mass from first particle[i][j]
     * i: atom index
     * j = 0: x
     * j = 1: y
     * j = 2: z
     */
    public double[][] getCentreOfMass1() {
        int tmpIterationSize;
        double tmpSumMass;
        double tmpCentreX;
        double tmpCentreY;
        double tmpCentreZ;
        double[][] tmpCentreOfMass;

        tmpIterationSize = this.coordinateList1.length;
        tmpSumMass = 0.0;
        tmpCentreX = 0.0;
        tmpCentreY = 0.0;
        tmpCentreZ = 0.0;
        tmpCentreOfMass = new double[tmpIterationSize][3];
        
        for (int i = 0; i < this.N_atom1 ; i++) {
            tmpSumMass += this.atomicMassList1[i];
        }
        
        for (int i = 0; i < tmpIterationSize; i++) {
            
            for (int j = 0; j < this.N_atom1; j++) {
                tmpCentreX += this.coordinateList1[i][j][0] * 
                        this.atomicMassList1[j];
                tmpCentreY += this.coordinateList1[i][j][1] * 
                        this.atomicMassList1[j];
                tmpCentreZ += this.coordinateList1[i][j][2] * 
                        this.atomicMassList1[j];
            }
            
            tmpCentreOfMass[i][0] = tmpCentreX / tmpSumMass;
            tmpCentreOfMass[i][1] = tmpCentreY / tmpSumMass;
            tmpCentreOfMass[i][2] = tmpCentreZ / tmpSumMass;
        }
        
        return tmpCentreOfMass;    
    }
    
    /** 
     * Returns the coordinate of centre of mass from second particle
     * 
     * @return coordinate of center of mass from second particle
     * [i][j][k]    i: index of simulation; 
     *              j: index of partikel2
     *              k = 0: x
     *                  1: y
     *                  2: z
     */
    public double[][][] getCentreOfMass2() {
        
        int tmpIterationSize;
        double tmpSumMass;
        double tmpCentreX;
        double tmpCentreY;
        double tmpCentreZ;
        double[][][] tmpCentreOfMass;

        tmpIterationSize = this.coordinateList1.length;
        tmpSumMass = 0.0;
        
        tmpCentreOfMass = new double[tmpIterationSize][this.N_particle2][3];
        
        for (int i = 0; i < this.N_atom2 ; i++) {
            tmpSumMass += this.atomicMassList2[i];
        }
        
        for (int i = 0; i < tmpIterationSize; i++) {

            for (int j = 0; j < this.N_particle2; j++) {
                tmpCentreX = 0.0;
                tmpCentreY = 0.0;
                tmpCentreZ = 0.0;

                for (int k = 0; k < this.N_atom2; k++) {
                    tmpCentreX += this.coordinateList2[i][j][k][0] * 
                            this.atomicMassList2[k];
                    tmpCentreY += this.coordinateList2[i][j][k][1] * 
                            this.atomicMassList2[k];
                    tmpCentreZ += this.coordinateList2[i][j][k][2] * 
                            this.atomicMassList2[k];
                }
                
                tmpCentreOfMass[i][j][0] = tmpCentreX / tmpSumMass;
                tmpCentreOfMass[i][j][1] = tmpCentreY / tmpSumMass;
                tmpCentreOfMass[i][j][2] = tmpCentreZ / tmpSumMass;
            }
            
        }
        
        return tmpCentreOfMass;    
    }
    
    
    /**
     * Returns the distances of atom1 to atom2
     *  i: atom index of particle1
     *  j: iteration index
     *  k: atom index of particle2
     * @return distances of atom1 to atom2
     */
    public double[][] getDistances() {
        if (this.distances != null) {
            return this.distances;
        } else {
            throw new NullPointerException("Return value of getDistances"
                    + "is null.");
        }
    }
    
    /** 
     * Returns the distances of centre of mass from first particle to 
     *  centre of mass from second particle(s)
     * 
     * @param aBoxLength Box length in Angstrom
     * @return distances between particle1 to particle2 
     * [i][j][k] i: index of simulation; j: atom index of partikel2
     * k = 0: x
     * k = 1: y
     * k = 2: z
     */
    public double[][] getPBCDistances(double aBoxLength) {
        final double tmpHalfBoxLength = 0.5 * aBoxLength;
        int tmpIterationSize;
        double[][] tmpDistance;
        double[][] tmpCentre1;
        double[][][] tmpCentre2;
        
        tmpIterationSize = this.coordinateList1.length;
        tmpDistance = new double[tmpIterationSize][this.N_particle2];
        tmpCentre1 = this.getCentreOfMass1();
        tmpCentre2 = this.getCentreOfMass2();
        
        for (int i = 0; i < tmpIterationSize; i++) {
            double c1x = tmpCentre1[i][0];
            double c1y = tmpCentre1[i][1];
            double c1z = tmpCentre1[i][2];
            double[][] currCentre2I = tmpCentre2[i];
            double[] currDistances = tmpDistance[i];
            
            for (int j = 0; j < this.N_particle2; j++) {
                double[] c2 = currCentre2I[j];
                double tmpDeltaX = c1x - c2[0];
                double tmpDeltaY = c1y - c2[1];
                double tmpDeltaZ = c1z - c2[2];
                if(tmpDeltaX > tmpHalfBoxLength) {
                    tmpDeltaX -= aBoxLength;
                } else if(tmpDeltaX <= -tmpHalfBoxLength) {
                    tmpDeltaX += aBoxLength;
                }
                if(tmpDeltaY > tmpHalfBoxLength) {
                    tmpDeltaY -= aBoxLength;
                } else if(tmpDeltaY <= -tmpHalfBoxLength) {
                    tmpDeltaY += aBoxLength;
                }
                if(tmpDeltaZ > tmpHalfBoxLength) {
                    tmpDeltaZ -= aBoxLength;
                } else if(tmpDeltaZ <= -tmpHalfBoxLength) {
                    tmpDeltaZ += aBoxLength;
                }
                currDistances[j] = Math.sqrt(tmpDeltaX * tmpDeltaX +
                        tmpDeltaY * tmpDeltaY +
                        tmpDeltaZ * tmpDeltaZ);
            }
            
        }
        
        return tmpDistance;
    }
    
    /**
     * Returns parameters of first particle (6. column of .txyz-file)
     * 
     * @return Parameters of first particle
     */
    public int[] getAtomTypeList1() {
        return this.atomTypeList1;
    }
    
    /**
     * Returns parameters of second particle (6. column of .txyz-file)
     * 
     * @return Parameters of second particle
     */
    public int[] getAtomTypeList2() {
        return this.atomTypeList2;
    }
    
    /**
     * Returns connection list of first particle
     * 
     * @return Connection list of first particle
     */
    public int[][] getConnectionList1() {
        return this.connectionList1;
    }
    
    /**
     * Returns connection list of second particle
     * 
     * @return Connection list of second particle
     */
    public int[][] getConnectionList2() {
        return this.connectionList2;
    }
    
    /**
     * Returns .txyz file content as stringbuilder
     * 
     * @return FileContent 
     */
    public StringBuilder getFileContent() {
        return this.fileContent;
    }
    
    /**
     * Set the header (also in fileContent)
     * 
     * @param aHeader A new Header
     * @param aTinkerOn Flag for whether tinker is used or not
     */
    public void setHeader(String aHeader, boolean aTinkerOn) {
        this.header = aHeader;
        
        if (!this.forcefieldName.equals("OPLSAALIGPARGEN") || aTinkerOn) {
            int tmpStartIndex;
            int tmpEndIndex;
            Boolean tmpIsFirstNonSpace;
            tmpStartIndex = 0;
            tmpIsFirstNonSpace = false;
            tmpEndIndex = this.fileContent.indexOf(LINESEPARATOR);
            String tmpHeader = this.fileContent.substring(0, tmpEndIndex);

            // Determine first space after non-space character
            for (int i = 0; i < tmpHeader.length(); i++) {
                if (tmpHeader.charAt(i) != ' ') {
                    tmpIsFirstNonSpace = true;
                } else {
                    if (tmpIsFirstNonSpace) {
                        tmpStartIndex = i + 1;
                        break;
                    }
                }
            }

            this.fileContent.replace(tmpStartIndex, tmpEndIndex, aHeader);
        }
    }
    
    /**
     * Set the name of force field
     * @param aParticleName Name of force field
     */
    public void setForcefieldName(String aParticleName) {
        // Check parameters
        if(aParticleName == null) {
            throw new IllegalArgumentException("Null was passed to the " 
                    + "setForcefieldName method.");
        }
        
        this.forcefieldName = aParticleName;
    }
    
    /**
     * Set the name of particle1
     * @param aParticleName Name of particle1
     */
    public void setParticleName1(String aParticleName) {
        // Check parameters
        if(aParticleName == null) {
            throw new IllegalArgumentException("Null was passed to the " 
                    + "setParticleName1 method.");
        }
        
        this.particleName1 = aParticleName;
    }
    
    /**
     * Set the name of particle2
     * @param aParticleName Name of particle2
     */
    public void setParticleName2(String aParticleName) {
        // Check parameters
        if(aParticleName == null) {
            throw new IllegalArgumentException("Null was passed to the " 
                    + "setParticleName2 method.");
        }
        
        this.particleName2 = aParticleName;
    }
    
    /**
     * Sets aFileContent attribute
     * 
     * @param aFileContent TinkerXYZ file content as StringBuilder object
     */
    public void setFileContent(StringBuilder aFileContent) {
        // Check parameters
        if(aFileContent == null) {
            throw new IllegalArgumentException("Null was passed to the " 
                    + "setFileContent method.");
        }
        
        this.fileContent = aFileContent;
    }

    /**
     * Overwrite the element list
     * 
     * @param anElementList A new element list
     */
    public void setElementList1 (String[] anElementList) {
        // Check parameters
        if (anElementList == null || anElementList.length == 0) {
            throw new IllegalArgumentException("Null or empty object was passed"
                    + "to the setElementList1 method.");
        }
        
        this.elementList1 = anElementList;
        int tmpSkipLines;
        if (this.comment == null || this.comment.isEmpty()) {
            tmpSkipLines = 1;
        }else{
            tmpSkipLines = 2;
        }
        int tmpStartIndex = 0;

        for (int i = 0; i < tmpSkipLines; i++) {
            tmpStartIndex = this.fileContent
                .indexOf(this.LINESEPARATOR, tmpStartIndex);
            tmpStartIndex += this.LINESEPARATOR.length();
        }
        
        String tmpElement;
        
        for (int i = 0; i < this.N_atom1; i++) {
            tmpStartIndex += 8;
            tmpElement = anElementList[i];
            this.fileContent.replace(tmpStartIndex, tmpStartIndex + 4,
                    MIPETUTIL.padRight(tmpElement, 3));
            tmpStartIndex = this.fileContent
                .indexOf(this.LINESEPARATOR, tmpStartIndex);
            tmpStartIndex += this.LINESEPARATOR.length();
        }
    }
    
    /**
     * Overwrite the element list
     * 
     * @param anElementList A new element list
     */
    public void setElementList2 (String[] anElementList) {
        
        // Check parameters
        if (anElementList == null || anElementList.length == 0) {
            throw new IllegalArgumentException("Null or empty object was passed"
                    + "to the setElementList1 method.");
        }
        
        this.elementList2 = anElementList;
        int tmpSkipLines;
        if (this.comment == null || this.comment.isEmpty()) {
            tmpSkipLines = 1;
        } else{
            tmpSkipLines = 2;
        }
        int tmpStartIndex = 0;
        tmpSkipLines += this.elementList1.length;

        for (int i = 0; i < tmpSkipLines; i++) {
            tmpStartIndex = this.fileContent
                .indexOf(this.LINESEPARATOR, tmpStartIndex);
            tmpStartIndex += this.LINESEPARATOR.length();
        }
        String tmpElement;
        
        for (int i = 0; i < this.N_atom2; i++) {
            tmpStartIndex += 8;
            tmpElement = anElementList[i];
            this.fileContent.replace(tmpStartIndex, tmpStartIndex + 4,
                    MIPETUTIL.padRight(tmpElement, 3));
            tmpStartIndex = this.fileContent
                .indexOf(this.LINESEPARATOR, tmpStartIndex);
            tmpStartIndex += this.LINESEPARATOR.length();
        }
    }
    
    /**
     * Set first coordinate list and overwrite stringbuilder
     * 
     * @param aCoordinates1 
     *   [i, j, k] i: simulation/config. index, j: atom, k: xyz
     * @param aTinkerOn
     *   Flag for whether tinker is used or not.
     */
    public void setCoordinateList1(double[][][] aCoordinates1, 
            boolean aTinkerOn) {
        
        // Check parameters
        if (aCoordinates1 == null || aCoordinates1.length == 0) {
            throw new IllegalArgumentException("Null or empty object was passed"
                    + "to the setCoordinateList1 method.");
        }
        
        this.coordinateList1 = aCoordinates1.clone();
        
        if (!this.forcefieldName.equals("OPLSAALIGPARGEN") || aTinkerOn) {
            int tmpSkipLines;
            DecimalFormat tmpDF = new DecimalFormat("0.000000", 
                    DecimalFormatSymbols.getInstance(Locale.ENGLISH));
            if (this.comment == null || this.comment.isEmpty()){
                tmpSkipLines = 1;
            }else{
                tmpSkipLines = 2;
            }
            int tmpStartIndex = 0;
            String tmpX;
            String tmpY;
            String tmpZ;

            for (int i = 0; i < tmpSkipLines; i++) {
                tmpStartIndex = this.fileContent
                        .indexOf(this.LINESEPARATOR, tmpStartIndex)
                        + this.LINESEPARATOR.length();
            }

            for (int i = 0; i < this.N_atom1; i++) {
                tmpStartIndex += 11;
                tmpX = tmpDF.format(aCoordinates1[0][i][0]);
                tmpX = MIPETUTIL.padLeft(tmpX, 10);
                this.fileContent.replace(tmpStartIndex, tmpStartIndex + 10, tmpX);
                tmpStartIndex += 11;
                tmpY = tmpDF.format(aCoordinates1[0][i][1]);
                tmpY = MIPETUTIL.padLeft(tmpY, 11);
                this.fileContent.replace(tmpStartIndex, tmpStartIndex + 11, tmpY);
                tmpStartIndex += 12;
                tmpZ = tmpDF.format(aCoordinates1[0][i][2]);
                tmpZ = MIPETUTIL.padLeft(tmpZ, 11);
                this.fileContent.replace(tmpStartIndex, tmpStartIndex + 11, tmpZ);
                tmpStartIndex = this.fileContent
                    .indexOf(this.LINESEPARATOR, tmpStartIndex)
                    + this.LINESEPARATOR.length();
            }
            
        }
    }
    
    /**
     * Overloaded method
     * Set first coordinate list and overwrite stringbuilder
     * 
     * @param aCoord1  [i][j] i: atomid, j: xyz
     * @param aTinkerOn Flag for whether tinker is used or not
     */
    public void setCoordinateList1(double[][] aCoord1, boolean aTinkerOn) {
        
        // Check parameters
        if (aCoord1 == null || aCoord1.length == 0) {
            throw new IllegalArgumentException("Null or empty object was passed"
                    + "to the setCoordinateList1 method.");
        }
        
        double[][][] tmpCoord1 = new double[1][][];
        tmpCoord1[0] = aCoord1;
        setCoordinateList1(tmpCoord1, aTinkerOn);
    }
    
    /**
     * Set second coordinate list and overwrite stringbuilder
     * 
     * @param aCoord2[i][j][k][l]
     *   i: sim.iteration j: particle k: atom l: xyz
     * @param aTinkerOn Flag for whether tinker is used or not
     */
    public void setCoordinateList2(double[][][][] aCoord2, boolean aTinkerOn) {
        
        // Check parameters
        if (aCoord2 == null || aCoord2.length == 0) {
            throw new IllegalArgumentException("Null or empty object was passed"
                    + "to the setCoordinateList2 method.");
        }
        
        this.coordinateList2 = aCoord2;
        if (!this.forcefieldName.equals("OPLSAALIGPARGEN") || aTinkerOn) {
            DecimalFormat tmpDF = new DecimalFormat("0.000000", 
            DecimalFormatSymbols.getInstance(Locale.ENGLISH));
            Boolean tmpHasComment = true;
            if (this.comment == null || this.comment.isEmpty()){
                tmpHasComment = false;
            }
            int tmpStartIndex = 0;
            int tmpSkip = this.N_atom1 + 1;
            if (tmpHasComment){
                tmpSkip++;
            }
            String tmpX;
            String tmpY;
            String tmpZ;

            for (int i = 0; i < tmpSkip; i++) {
                tmpStartIndex = this.fileContent
                        .indexOf(this.LINESEPARATOR, tmpStartIndex)
                        + this.LINESEPARATOR.length();
            }

            for (int i = 0; i < this.N_atom2; i++) {
                tmpStartIndex += 11;
                tmpX = tmpDF.format(aCoord2[0][0][i][0]);
                tmpX = MIPETUTIL.padLeft(tmpX, 10);
                this.fileContent.replace(tmpStartIndex, tmpStartIndex + 10, tmpX);
                tmpStartIndex += 11;
                tmpY = tmpDF.format(aCoord2[0][0][i][1]);
                tmpY = MIPETUTIL.padLeft(tmpY, 11);
                this.fileContent.replace(tmpStartIndex, tmpStartIndex + 11, tmpY);
                tmpStartIndex += 12;
                tmpZ = tmpDF.format(aCoord2[0][0][i][2]);
                tmpZ = MIPETUTIL.padLeft(tmpZ, 11);
                this.fileContent.replace(tmpStartIndex, tmpStartIndex + 11, tmpZ);
                tmpStartIndex = this.fileContent
                    .indexOf(this.LINESEPARATOR, tmpStartIndex)
                        + this.LINESEPARATOR.length();
            }
            
        }
    }
    
    /**
     * Overloaded method
     * Set second coordinate list and overwrite stringbuilder
     * 
     * @param aCoord2 
     * [i, j] i: atomid, j: xyz
     * @param aTinkerOn Flag for whether tinker is used or not
     */
    public void setCoordinateList2(double[][] aCoord2, boolean aTinkerOn) {
        
        // Check parameters
        if (aCoord2 == null || aCoord2.length == 0) {
            throw new IllegalArgumentException("Null or empty object was passed"
                    + "to the setCoordinateList2 method.");
        }
        
        double[][][][] tmpCoord2;
                
        tmpCoord2 = new double[1][1][][];
        tmpCoord2[0][0] = aCoord2;
        setCoordinateList2(tmpCoord2, aTinkerOn);
    }
    
    /**
     * Set atom types of particle 1
     */
    public void setAtomTypeList1() {
        int tmpAtomSize1;
        int[] tmpAtomTypes;
        
        tmpAtomSize1 = this.N_atom1;
        tmpAtomTypes = new int[tmpAtomSize1];
        
        for (int i = 0; i < tmpAtomSize1; i++) {
            tmpAtomTypes[i] = this.atomTypeList1[i] + 100;
        }
        
        this.atomTypeList1 = tmpAtomTypes;
    }
    
    /**
     * Set distances
     */
    public void setDistances() {
        
        // Check coordinates
        if (this.coordinateList1 == null || this.coordinateList1.length == 0) {
            throw new IllegalArgumentException("Null or empty coordinateList1"
                    + "found in setDistances method.");
        }
        if (this.coordinateList2 == null || this.coordinateList2.length == 0) {
            throw new IllegalArgumentException("Null or empty coordinateList2"
                    + "found in setDistances method.");
        }
        
        this.distances = 
                new double[this.N_atom1][this.N_atom2];
        double[][] coords1 = this.coordinateList1[0];
        double[][] coords2 = this.coordinateList2[0][0];
        
        for (int i = 0; i < this.N_atom1; i++) {
            double[] distRow = this.distances[i];
            double[] coord1 = coords1[i];
            double x1 = coord1[0];
            double y1 = coord1[1];
            double z1 = coord1[2];
                
            for (int j = 0; j < this.N_atom2; j++) {
                double[] coord2 = coords2[j];
                double tmpDeltaX = coord2[0] - x1;
                double tmpDeltaY = coord2[1] - y1;
                double tmpDeltaZ = coord2[2] - z1;
                distRow[j] = Math.sqrt(tmpDeltaX * tmpDeltaX
                        + tmpDeltaY * tmpDeltaY 
                        + tmpDeltaZ * tmpDeltaZ);
            }
            
        }
        
    }
    
    // </editor-fold>
        
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * initialize mothod
     * 
     * @param aTinkerXyz
     *   Tinker xyz content
     * @param itrerationSize
     *   number of iteration
     * @param anAtomSize1
     *   atom number of first particle
     * @param anAtomSize2
     *   atom number of second particle
     */
    private void initialize(String aTinkerXyz, int anIterationSize, 
            int anAtomSize1, int anAtomSize2) {
        
        // Check parameters
        if (aTinkerXyz == null || aTinkerXyz.isEmpty()) {
            throw new IllegalArgumentException(".xyz file is null or empty.");
        } else if (anIterationSize <= 0) {
            throw new IllegalArgumentException("iterationSize should be not "
                    + "negative or zero.");
        } else if (anAtomSize1 < 0) {
            throw new IllegalArgumentException("anAtomSize1 should be not "
                    + "negative.");
        } else if (anAtomSize2 < 0) {
            throw new IllegalArgumentException("anAtomSize2 should be not "
                    + "negative.");
        }
        
        int tmpAtomNumber;
        int tmpConnectionSize;
        boolean tmpHasCommentLine;
        String tmpReadLine;
        String[] tmpLines;
        Reader tmpStringReader;
        
        tmpStringReader = new StringReader(aTinkerXyz);
        tmpHasCommentLine = false;
        this.fileContent = new StringBuilder(STRINGBUILDER_CAPACITY);
        
        try (BufferedReader tmpBR = new BufferedReader(tmpStringReader)) {
            
            // read first line
            tmpReadLine = tmpBR.readLine();
            if (tmpReadLine != null) {
                this.fileContent.append(tmpReadLine);
                this.fileContent.append(this.LINESEPARATOR);
                tmpLines = tmpReadLine.trim().split("\\s+");
                tmpAtomNumber = Integer.parseInt(tmpLines[0]);
                this.atomNumber = tmpAtomNumber;
                this.header = "";
                if (tmpLines != null && tmpLines.length > 1) {

                    for (int i = 1; i < tmpLines.length; i++) {
                        this.header +=  tmpLines[i];
                    }

                }
                if (anAtomSize2 == 0) {
                    this.N_atom1 = tmpAtomNumber;
                    this.N_atom2 = 0;
                    this.N_particle2 = 0;
                } else {
                    this.N_atom1 = anAtomSize1;
                    this.N_atom2 = anAtomSize2;
                    this.N_particle2 = (tmpAtomNumber - anAtomSize1) 
                            / anAtomSize2;
                    this.coordinateList2 = new double[anIterationSize]
                            [this.N_particle2][this.N_atom2][3];
                }
                this.coordinateList1 = 
                        new double[anIterationSize][this.N_atom1][3];
            
                // read second line
                tmpBR.mark(80);
                tmpReadLine = tmpBR.readLine();
                tmpLines = tmpReadLine.trim().split("\\s+");
                if (!tmpLines[0].equals("1")) {
                    tmpHasCommentLine = true;
                    this.comment = tmpReadLine;
                    this.fileContent
                            .append(tmpReadLine)
                            .append(this.LINESEPARATOR);
                } else {
                    tmpBR.reset();
                }
                this.elementList1 = new String[this.N_atom1];
                this.elementList2 = new String[anAtomSize2];
                this.atomTypeList1 = new int[this.N_atom1];
                this.atomTypeList2 = new int[anAtomSize2]; 
                this.connectionList1 = new int[this.N_atom1][];
                this.connectionList2 = new int[anAtomSize2][];

                for (int i = 0; i < anIterationSize; i++) {

                    for (int j = 0; j < this.N_atom1; j++) {
                        String tmpLine = tmpBR.readLine();
                        this.fileContent
                                .append(tmpLine)
                                .append(this.LINESEPARATOR);
                        String[] tmpTokens = tmpLine.trim().split("\\s+");
                        double x = Double.parseDouble(tmpTokens[2]);
                        double y = Double.parseDouble(tmpTokens[3]);
                        double z = Double.parseDouble(tmpTokens[4]);
                        this.coordinateList1[i][j] = new double[]{x, y, z};
                        if (i == 0) {
                            this.elementList1[j] = tmpTokens[1];
                            this.atomTypeList1[j] = Integer
                                    .parseInt(tmpTokens[5]);
                            int connectionSize = tmpTokens.length - 6;
                            this.connectionList1[j] = new int[connectionSize];

                            for (int k = 0; k < connectionSize; k++) {
                                this.connectionList1[j][k] = Integer
                                        .parseInt(tmpTokens[k + 6]);
                            }
                            
                        }
                    }

                    // read following lines of second particle
                    for (int j = 0; j < this.N_particle2; j++) {

                        for (int k = 0; k < this.N_atom2; k++) {
                            tmpReadLine = tmpBR.readLine();
                            this.fileContent.append(tmpReadLine)
                                    .append(this.LINESEPARATOR);
                            tmpLines = tmpReadLine.trim().split("\\s+");
                            double x = Double.parseDouble(tmpLines[2]);
                            double y = Double.parseDouble(tmpLines[3]);
                            double z = Double.parseDouble(tmpLines[4]);
                            this.coordinateList2[i][j][k] = 
                                    new double[]{x, y, z};
                            tmpConnectionSize = tmpLines.length - 6;
                            if(i == 0 && j == 0) {
                                this.elementList2[k] = tmpLines[1];
                                this.atomTypeList2[k] = Integer.
                                        parseInt(tmpLines[5]);
                                this.connectionList2[k] = 
                                        new int[tmpConnectionSize];
                                for (int l = 0; l < tmpConnectionSize; l++) {
                                        this.connectionList2[k][l] = Integer
                                                .parseInt(tmpLines[l + 6]);
                                }
                                
                            }
                        }
                        
                    }

                    if (i < anIterationSize) {
                        tmpBR.readLine();
                        if (tmpHasCommentLine) {
                            tmpBR.readLine();
                        }
                    }
                }
            }
        } catch (FileNotFoundException ex) {
            throw new IllegalArgumentException("File not found in TinkerXYZ.");
        } catch (IOException ex) {
            throw new IllegalArgumentException("IOException in TinkerXYZ.");            
        }
        
        // Set atomic mass map
        double tmpAtomicMass;
        
        for (int i = 0; i < this.N_atom1; i++) {
            if (!this.atomicMassMap.containsKey(this.elementList1[i])) {
                tmpAtomicMass = MIPETUTIL.getAtomicMass(
                        this.elementList1[i], false);
                this.atomicMassMap.put(this.elementList1[i], tmpAtomicMass);
            }
        }
        
        for (int i = 0; i < this.N_atom2; i++) {
            if (!this.atomicMassMap.containsKey(this.elementList2[i])) {
                tmpAtomicMass = MIPETUTIL.getAtomicMass(
                        this.elementList2[i], false);
                this.atomicMassMap.put(this.elementList2[i], tmpAtomicMass);
            }
        }
        
        // Set atomic mass list
        this.atomicMassList1 = this.getAtomicMassList1();
        if (this.N_atom2 > 0) {
            this.atomicMassList2 = this.getAtomicMassList2();
        }
    }
   
    /**
     * Build a new TinkerXYZ object out of two TinkerXYZ objects
     * 
     * @param aTinkerXYZ1
     *   First TinkerXYZ object
     * @param aTinkerXYZ2 
     *   Second TinkerXYZ object
     * @param aTinkerOn Flag for whether tinker is used or not
     */
    private void initialize2(TinkerXYZ aTinkerXYZ1, TinkerXYZ aTinkerXYZ2,
            boolean aTinkerOn) {
        
        // Check parameters
        if (aTinkerXYZ1 == null) {
            throw new IllegalArgumentException("aTinkerXYZ1 is null.");
        } else if (aTinkerXYZ2 == null) {
            throw new IllegalArgumentException("aTinkerXYZ2 is null.");
        }
        
        int tmpAtomSize1;
        int tmpAtomSize2;
        int tmpNConnection;
        String tmpAtomNumber;
        DecimalFormat tmpDF;
        String tmpX;
        String tmpY;
        String tmpZ;
        String tmpParameter;
        String tmpConnection;
        TinkerXYZ tmpTinkerXyz1;
        TinkerXYZ tmpTinkerXyz2;

        tmpTinkerXyz1 = aTinkerXYZ1;
        tmpTinkerXyz2 = aTinkerXYZ2;
        tmpDF = new DecimalFormat("0.000000", DecimalFormatSymbols
                .getInstance(Locale.ENGLISH));
        tmpAtomSize1 = tmpTinkerXyz1.getN_atom();
        tmpAtomSize2 = tmpTinkerXyz2.getN_atom();
        this.forcefieldName = tmpTinkerXyz1.getForcefieldName();
        this.particleName1 = tmpTinkerXyz1.getParticleName1();
        this.particleName2 = tmpTinkerXyz2.getParticleName1();
        this.atomNumber = tmpAtomSize1 + tmpAtomSize2;
        this.N_atom1 = tmpAtomSize1;
        this.N_atom2 = tmpAtomSize2;
        this.N_particle2 = (this.atomNumber - this.N_atom1) 
                            / this.N_atom2;
        this.elementList1 = tmpTinkerXyz1.getElementList1().clone();
        this.elementList2 = tmpTinkerXyz2.getElementList1().clone();
        this.atomicMassList1 = tmpTinkerXyz1.getAtomicMassList1().clone();
        this.atomicMassList2 = tmpTinkerXyz2.getAtomicMassList2().clone();
        this.atomTypeList1 = tmpTinkerXyz1.getAtomTypeList1().clone();
        this.atomTypeList2 = tmpTinkerXyz2.getAtomTypeList1().clone();
        
        this.connectionList1 = tmpTinkerXyz1.getConnectionList1().clone();
        this.connectionList2 = tmpTinkerXyz2.getConnectionList1().clone();
        this.connectionList2 = this.correctConnectionList(this.connectionList2);
        this.coordinateList1 = new double[1][tmpAtomSize1][3];
        this.coordinateList1 = tmpTinkerXyz1.coordinateList1.clone();
        this.coordinateList2 = new double[1][1][tmpAtomSize2][3];
        this.coordinateList2[0] = tmpTinkerXyz2.coordinateList1.clone();
        if (!this.forcefieldName.equals("OPLSAALIGPARGEN") || aTinkerOn) {
            this.fileContent = new StringBuilder(STRINGBUILDER_CAPACITY);
            this.fileContent.append(tmpTinkerXyz1.getFileContent());
            this.fileContent.replace(0, 6, MIPETUTIL.padLeft(Integer
                    .toString(this.atomNumber) + "  ", 6));
            
            for (int i = 0; i < tmpAtomSize2; i++) {
                tmpAtomNumber = Integer.toString(i + tmpAtomSize1 + 1);
                this.fileContent.append(MIPETUTIL.padLeft(tmpAtomNumber, 6));
                this.fileContent.append("  ");
                this.fileContent.append(MIPETUTIL.padRight(tmpTinkerXyz2
                        .getElementList1()[i], 2));
                tmpX = tmpDF.format(tmpTinkerXyz2
                        .getCoordinateList1()[0][i][0]);
                this.fileContent.append(MIPETUTIL.padLeft(tmpX, 11));
                tmpY = tmpDF.format(tmpTinkerXyz2
                        .getCoordinateList1()[0][i][1]);
                this.fileContent.append(MIPETUTIL.padLeft(tmpY, 12));
                tmpZ = tmpDF.format(tmpTinkerXyz2
                        .getCoordinateList1()[0][i][2]);
                this.fileContent.append(MIPETUTIL.padLeft(tmpZ, 12));
                tmpParameter = Integer.toString(tmpTinkerXyz2
                        .getAtomTypeList1()[i]);
                this.fileContent.append(MIPETUTIL.padLeft(tmpParameter, 6));
                tmpNConnection = tmpTinkerXyz2
                        .getConnectionList1()[i].length;

                for (int j = 0; j < tmpNConnection; j++) {
                    tmpConnection = Integer.toString( 
                            tmpTinkerXyz2.getConnectionList1()[i][j] + tmpAtomSize1);
                    this.fileContent.append(MIPETUTIL.padLeft(tmpConnection, 6));
                }

                this.fileContent.append(this.LINESEPARATOR);
            }
            
        }
    }
    
    /**
     * Correct connection list for the second particle
     * 
     * @param aIntegerList Old connection list
     * @return Corrected connection list
     */
    private int[][] correctConnectionList(int[][] aIntegerList) {
        int[][] tmpResult = new int[aIntegerList.length][];
        
        for (int i = 0; i < aIntegerList.length; i++) {
            tmpResult[i] = new int[aIntegerList[i].length];
            
            for (int j = 0; j < aIntegerList[i].length; j++) {
                tmpResult[i][j] = aIntegerList[i][j] + this.N_atom1;
            }
            
        }
        
        return tmpResult;
    }
    
    /**
     * Clone method for 2d-array
     * 
     * @param AnOriginal Original 2d-array
     * @return Cloned deep copied 2d-array
     */
    private double[][] clone2DArray(double[][] AnOriginal) {
        if (AnOriginal == null) {
            return null;
        }
        double[][] tmpCopy = new double[AnOriginal.length][];

        for (int i = 0; i < AnOriginal.length; i++) {
            if (AnOriginal[i] != null) {
                tmpCopy[i] = AnOriginal[i].clone();
            }
        }

        return tmpCopy;
    }
    
    /**
     * Clone method for 3d-array
     * 
     * @param AnOriginal Original 3d-array
     * @return Cloned deep copied 3d-array
     */
    private double[][][] clone3DArray(double[][][] AnOriginal) {
        if (AnOriginal == null) {
            return null;
        }

        double[][][] tmpCopy = new double[AnOriginal.length][][];

        for (int i = 0; i < AnOriginal.length; i++) {
            if (AnOriginal[i] != null) {
                tmpCopy[i] = new double[AnOriginal[i].length][];

                for (int j = 0; j < AnOriginal[i].length; j++) {
                    if (AnOriginal[i][j] != null) {
                        tmpCopy[i][j] = AnOriginal[i][j].clone();
                    }
                }
                
            }
        }

        return tmpCopy;
    }
    
    /**
     * Clone method for 4d-array
     * 
     * @param AnOriginal Original 4d-array
     * @return Cloned deep copied 4d-array
     */
    private double[][][][] clone4DArray(double[][][][] AnOriginal) {
        if (AnOriginal == null) {
            return null;
        }

        double[][][][] tmpCopy = new double[AnOriginal.length][][][];

        for (int i = 0; i < AnOriginal.length; i++) {
            if (AnOriginal[i] != null) {
                tmpCopy[i] = new double[AnOriginal[i].length][][];

                for (int j = 0; j < AnOriginal[i].length; j++) {
                    if (AnOriginal[i][j] != null) {
                        tmpCopy[i][j] = AnOriginal[i][j].clone();
                    }
                    
                    for (int k = 0; k < AnOriginal[i][j].length; k++) {
                        if (AnOriginal[i][j][k] != null) {
                            tmpCopy[i][j][k] = AnOriginal[i][j][k].clone();
                        }
                    }
                    
                }
                
            }
        }

        return tmpCopy;
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    
    /**
     * Find the mass independent centre coordinate
     * 
     * @return Mass independent centre coordinate
     */
    public double[] findCentreCoordinate() {
        return this.findCentreCoordinate(this.coordinateList1[0]);
    }
    
    /**
     * Find the mass independent centre coordinate
     * 
     * @param aCoords
     *   Coordinates of atoms
     * @return Mass independent centre coordinate
     */
    public double[] findCentreCoordinate(double[][] aCoords) {
        
        // Check parameters
        if (aCoords == null || aCoords.length == 0) {
            throw new IllegalArgumentException("Null or empty object was passed"
                    + "to the findCentreCoordinate method.");
        }
        
        double tmpX = 0.0;
        double tmpY = 0.0;
        double tmpZ = 0.0;
        double[] tmpCentreCoord = new double[3];
        
        for (double[] tmpCoord : aCoords) {
            tmpX += tmpCoord[0];
            tmpY += tmpCoord[1];
            tmpZ += tmpCoord[2];
        }
        tmpCentreCoord[0] = tmpX / aCoords.length; 
        tmpCentreCoord[1] = tmpY / aCoords.length;
        tmpCentreCoord[2] = tmpZ / aCoords.length;
        return tmpCentreCoord;
    }
    
    /**
     * Move coordinates of all fragment atoms by the given vector
     * 
     * @param aCoords Coordinates of atoms
     * @param aVector Vector coordinate
     * @return  Moved coordinates of atoms
     */
    public double[][] moveCoordinates(double[][] aCoords, 
            double[] aVector) {
        
        // Check parameters
        if (aCoords == null || aCoords.length == 0) {
            throw new IllegalArgumentException("aCoords is null or empty.");
        } else if (aVector == null || aVector.length == 0) {
            throw new IllegalArgumentException("aVector is null or empty.");
        }
                
        double[][] tmpCoordinates = new double[aCoords.length][3];

        for (int i = 0; i < aCoords.length; i++) {
            tmpCoordinates[i][0] = aCoords[i][0] - aVector[0];
            tmpCoordinates[i][1] = aCoords[i][1] - aVector[1];
            tmpCoordinates[i][2] = aCoords[i][2] - aVector[2];
        }
        
        return tmpCoordinates;
    }
        
    /**
     * Converts tinker object to .xyz file
     * 
     * @param aXyzFileName File name of .xyz file
     */
    public void writeToXyzFile(String aXyzFileName) {
        
        // Check parameters
        if (aXyzFileName == null || aXyzFileName.isEmpty()) {
            throw new IllegalArgumentException("aXyzFileName in "
                    + "writeToXyzFile() is null.");
        }
        
        File tmpTargetFile = new File(aXyzFileName);
        int tmpCoordSize1 = this.coordinateList1[0].length;
        int tmpCoordSize2;
        
        if (this.coordinateList2 != null) {
            tmpCoordSize2 = this.coordinateList2[0][0].length;
        } else {
            tmpCoordSize2 = 0;
        }
        int tmpCoordSize = tmpCoordSize1 + tmpCoordSize2;
        String[] tmpElements;
        double[][] tmpCoordList = new double[tmpCoordSize][3]; 
        
        tmpElements = ArrayUtils.addAll(this.elementList1, this.elementList2);
        
        for (int i = 0; i < this.coordinateList1[0].length; i++) {
            tmpCoordList[i][0] = this.coordinateList1[0][i][0];
            tmpCoordList[i][1] = this.coordinateList1[0][i][1];
            tmpCoordList[i][2] = this.coordinateList1[0][i][2];
        }
        
        if (this.coordinateList2 != null) {
            for (int i = 0; i < this.coordinateList2[0][0].length; i++) {
                tmpCoordList[i + tmpCoordSize1][0] = this.
                        coordinateList2[0][0][i][0];
                tmpCoordList[i + tmpCoordSize1][1] = this.
                        coordinateList2[0][0][i][1];
                tmpCoordList[i + tmpCoordSize1][2] = this.
                        coordinateList2[0][0][i][2];
            }
        }
                
        try (BufferedWriter tmpBW = new BufferedWriter(
                new FileWriter(tmpTargetFile))) {
            // first line
            tmpBW.append(Integer.toString(this.atomNumber)
                    +"    "
                    + this.header 
                    + this.LINESEPARATOR);
        
            // comment line
            if (this.comment != null) {
                tmpBW.append(this.comment + this.LINESEPARATOR);
            }
        
            // body lines
            new Locale.Builder().setLanguage("en")
                    .setRegion("UK")
                    .build();
            DecimalFormat decimalFormat = (DecimalFormat)NumberFormat
                .getNumberInstance();
            decimalFormat.applyPattern("#0.000000");
            
            for (int i = 0; i < this.atomNumber; i++) {
                tmpBW.append("  ");
                tmpBW.append(MIPETUTIL.padRight(tmpElements[i], 2));
                tmpBW.append(MIPETUTIL.padLeft(decimalFormat
                    .format(tmpCoordList[i][0]), 12));
                tmpBW.append(MIPETUTIL.padLeft(decimalFormat
                    .format(tmpCoordList[i][1]), 12));
                tmpBW.append(MIPETUTIL.padLeft(decimalFormat
                    .format(tmpCoordList[i][2]), 12));
                tmpBW.append(LINESEPARATOR);
            }
            tmpBW.close();
        } catch(IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException in writeToXyzFile().", ex);
        }
    }
    
    /**
     * Converts tinker xyz object to .xyz file
     * 
     * @param aXyzFileName .xyz file name
     * @param aData Stringbuilder object with data
     */
    public void writeToXyzFile(String aXyzFileName, StringBuilder aData) {
        
        // Check parameters
        if (aXyzFileName == null || aXyzFileName.isEmpty()) {
            throw new IllegalArgumentException("aXyzFileName in " 
                    + "writeToXyzFile() is null.");
        } else if (aData == null || aData.isEmpty()) {
            throw new IllegalArgumentException("aData in " 
                    + "writeToXyzFile() is null or empty.");
        }
        
        try (BufferedWriter tmpBW = new BufferedWriter(
                new FileWriter(aXyzFileName))) {
            tmpBW.append(aData);
        } catch(IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException in writeToXyzFile().", ex);
        }
    }
    
    /**
     * Converts tinker xyz object to .txyz file
     * 
     * @param aTxyzFileName .txyz file
     */
    public void writeToTxyzFile(String aTxyzFileName) {
        
        // Check parameters
        if (aTxyzFileName == null || aTxyzFileName.isEmpty()) {
            throw new IllegalArgumentException("aTxyzFileName in " 
                    + "writeToTxyzFile() is null or empty.");
        }
        
        this.writeToXyzFile(aTxyzFileName, this.fileContent);
    }
        
    /**
     * Read the coordinates of a .arc file and give it as double[][][]
     *   [i,j,k]: i: configuration j: atom k: xyz
     * 
     * @param aFileName A file name
     * @return Coordinates data of a .arc file
     */
    public double[][][] readCoordFromArc(String aFileName) {
        
        // Check parameters
        if (aFileName == null || aFileName.isEmpty()) {
            throw new IllegalArgumentException("aFileName in readCoordFromArc is null or empty.");
        }
        
        File tmpSourceFile = new File(aFileName);
        int tmpAtomSize = 0;
        int tmpLineCounter = 0;
        int tmpAtomIndex = 0;
        String tmpLine;
        ArrayList<double[][]> tmpCoordsList = new ArrayList<>();
        double[][] tmpCoords;
        double[][][] tmpResult;
        
        try (BufferedReader tmpBR = new BufferedReader (
                new FileReader (tmpSourceFile), this.READER_BUFFERSIZE)) {
            // Read first line
            tmpBR.mark(80);
            tmpLine = tmpBR.readLine();
            tmpAtomSize = Integer.parseInt(tmpLine.substring(0, 6).trim());
            tmpLineCounter++;
            tmpBR.reset();
            tmpCoords = new double [tmpAtomSize][3];
            
            while((tmpLine = tmpBR.readLine()) != null) {
                if (tmpLineCounter % (tmpAtomSize + 1) != 1) {
                    tmpCoords[tmpAtomIndex][0] = Double.parseDouble(tmpLine
                            .substring(12, 23));
                    tmpCoords[tmpAtomIndex][1] = Double.parseDouble(tmpLine
                            .substring(24, 35));
                    tmpCoords[tmpAtomIndex][2] = Double.parseDouble(tmpLine
                            .substring(36, 47));
                    tmpAtomIndex++;
                    if (tmpLineCounter % (tmpAtomSize + 1) == 0) {
                        tmpCoordsList.add(tmpCoords.clone());
                        tmpCoords = new double [tmpAtomSize][3];
                        tmpAtomIndex = 0;
                    }
                }
                tmpLineCounter++;
            }
            
            tmpBR.close();
        } catch(IOException ex) {
            LOGGER.log(Level.SEVERE, 
                    "IOException in readCoordFromArc().", ex);
        }
        tmpResult = new double[tmpCoordsList.size()][tmpAtomSize][3];
        
        for (int i = 0; i < tmpCoordsList.size(); i++) {
            tmpResult[i] = tmpCoordsList.get(i);
        }
        
        return tmpResult;
    }
    
    /**
     * Make .arc file
     * 
     * @param aFileName Output file name
     */
    public void makeArcFile(String aFileName) {
        // Check parameters
        if (aFileName == null || aFileName.isEmpty()) {
            throw new IllegalArgumentException("aFileName in makeArcFile is null or empty.");
        }
        
        int tmpIndex;
        int tmpLineNumber;
        DecimalFormat decimal6;
        StringBuilder tmpContent;
        
        tmpLineNumber = this.N_atom1 + this.N_atom2 + 5;
        tmpContent = new StringBuilder(tmpLineNumber * 80);
        decimal6 = (DecimalFormat)NumberFormat.getNumberInstance();
        decimal6.applyPattern("#0.000000");
        
        // Fill the content from TinkerXYZ object
        tmpContent.append(MIPETUTIL
                .padLeft(String.valueOf(this.atomNumber), 9))
                .append(" ")
                .append(this.header)
                .append(LINESEPARATOR);
        tmpIndex = 1;
        
        for (int i = 0; i < this.N_atom1; i++) {
            double[] tmpCurrCoords = this.coordinateList1[0][i];
            int[] tmpCurrConnections = this.connectionList1[i];
            
            tmpContent.append(MIPETUTIL.padLeft(String.valueOf(tmpIndex), 9))
                    .append("   ")
                    .append(MIPETUTIL.padRight(this.elementList1[i], 3))
                    .append(MIPETUTIL.padLeft(decimal6.format(
                            tmpCurrCoords[0]), 12))
                    .append(MIPETUTIL.padLeft(decimal6.format(
                            tmpCurrCoords[1]), 12))
                    .append(MIPETUTIL.padLeft(decimal6.format(
                            tmpCurrCoords[2]), 12))
                    .append(MIPETUTIL.padLeft(
                            Integer.toString(this.atomTypeList1[i]), 6));
            
            tmpIndex ++;
            int tmpConnectionSize = tmpCurrConnections.length;
            
            for (int j = 0; j < tmpConnectionSize; j++) {
                tmpContent.append(MIPETUTIL.padLeft(
                    Integer.toString(tmpCurrConnections[j]), 6));
            }
            
            if (i < this.N_atom1) {
                tmpContent.append(LINESEPARATOR);
            }
        }
        
        for (int i = 0; i < N_atom2; i++) {
            double[] tmpCurrCoords = this.coordinateList2[0][0][i];
            int[] tmpCurrConnections = this.connectionList2[i];
            
            tmpContent.append(MIPETUTIL.padLeft(String.valueOf(tmpIndex), 9))
                    .append("   ")
                    .append(MIPETUTIL.padRight(this.elementList2[i], 3))
                    .append(MIPETUTIL.padLeft(decimal6
                            .format(tmpCurrCoords[0]), 12))
                    .append(MIPETUTIL.padLeft(decimal6
                            .format(tmpCurrCoords[1]), 12))
                    .append(MIPETUTIL.padLeft(decimal6
                            .format(tmpCurrCoords[2]), 12))
                    .append(MIPETUTIL.padLeft(Integer
                            .toString(this.atomTypeList2[i]), 6));
            tmpIndex++;
            int tmpConnectionSize = tmpCurrConnections.length;
            
            for (int j = 0; j < tmpConnectionSize; j++) {
                tmpContent.append(MIPETUTIL.padLeft(
                        Integer.toString(tmpCurrConnections[j]), 6));
            }
            
            if (i < this.N_atom2 - 1) {
                tmpContent.append(LINESEPARATOR);
            }
        }
        
        try (PrintWriter tmpOut = new PrintWriter((aFileName))) {
            tmpOut.print(tmpContent);
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE,
                    "IOException during writing .arc file in makeArcFile.", ex);
        }
    }
       
    /**
     * Clone method for TinkerXYZ class
     */
    @Override
    public TinkerXYZ clone() {
        try {
            TinkerXYZ cloned = (TinkerXYZ) super.clone();
            
            if (this.atomicMassList1 != null) {
                cloned.atomicMassList1 = this.atomicMassList1.clone();
            } 
            if (this.atomTypeList1 != null) {
                cloned.atomTypeList1 = this.atomTypeList1.clone();
            }
            if (this.coordinateList1 != null) {
                cloned.coordinateList1 = clone3DArray(this.coordinateList1);
            } 
            if (this.N_atom2 > 0) {
                if (this.atomicMassList2 != null) 
                    cloned.atomicMassList2 = this.atomicMassList2.clone();
                if (this.atomTypeList2 != null) 
                    cloned.atomTypeList2 = this.atomTypeList2.clone();
                if (this.coordinateList2 != null) 
                    cloned.coordinateList2 = clone4DArray(this.coordinateList2);
                if (this.distances != null) 
                    cloned.distances = this.clone2DArray(this.distances);
            }
            return cloned;
        } catch (final CloneNotSupportedException ex) {
            throw new AssertionError("Cloning failed.", ex);
        }
    }
    
    // </editor-fold>
}