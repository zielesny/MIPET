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
    private int atomSize1;
    
    /**
     * Atom number of second particle
     */
    private int atomSize2;
    
    /**
     * Second particle number
     */
    private int particleSize2;
    
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
     * [i, j, k] i: sim.iteration j: atomid k: xyz
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
     * @param aTxyzFileName Tinker xyz filename
     * @param atomSize Atom number of particle
     */
    public TinkerXYZ(String aTxyzFileName, int atomSize) {
        this.initialize(aTxyzFileName, 1, atomSize, 0);
    }
    
    /**
     * Constructor TinkerXYZ
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
    
    public String getForcefieldName() {
        return this.forcefieldName;
    }
    
    /**
     * Returns the atom number
     * 
     * @return the atom number
     */
    public int getAtomNumber() {
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
     * Returns the atom number of the first particle
     * 
     * @return Atom number of first particle
     */
    public int getAtomSize1() {
        return this.atomSize1;
    }
    
    /**
     * Returns the atom number of the second particle
     * 
     * @return Atom number of second particle
     */
    public int getAtomSize2() {
        return this.atomSize2;
    }
    
    /**
     * Returns the particle number of the second particle
     * 
     * @return Number of solvent particles in a simulation box
     */
    public int getParticleSize2() {
        return this.particleSize2;
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
            this.atomicMassList1 = new double[this.atomSize1];
            
            for (int i = 0; i < this.atomSize1; i++) {
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
            this.atomicMassList2 = new double[this.atomSize2];
            
            for (int i = 0; i < this.atomSize2; i++) {
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
        
        for (int i = 0; i < this.atomSize1 ; i++) {
            tmpSumMass += this.atomicMassList1[i];
        }
        
        for (int i = 0; i < tmpIterationSize; i++) {
            
            for (int j = 0; j < this.atomSize1; j++) {
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
     * [i][j][k] i: index of simulation; j: index of partikel2
     * k = 0: x
     * k = 0: y
     * k = 0: z
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
        
        tmpCentreOfMass = new double[tmpIterationSize][this.particleSize2][3];
        
        for (int i = 0; i < this.atomSize2 ; i++) {
            tmpSumMass += this.atomicMassList2[i];
        }
        
        for (int i = 0; i < tmpIterationSize; i++) {

            for (int j = 0; j < this.particleSize2; j++) {
                tmpCentreX = 0.0;
                tmpCentreY = 0.0;
                tmpCentreZ = 0.0;

                for (int k = 0; k < this.atomSize2; k++) {
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
        final double ONEHALF = 0.5;
        int tmpIterationSize;
        double[][] tmpCentre1;
        double[][][] tmpCentre2;
        double tmpDeltaX;
        double tmpDeltaY;
        double tmpDeltaZ;
        
        tmpIterationSize = this.coordinateList1.length;
        double[][] tmpDistance = new double[tmpIterationSize][this.particleSize2];
        tmpCentre1 = this.getCentreOfMass1();
        tmpCentre2 = this.getCentreOfMass2();
        
        for (int i = 0; i < tmpIterationSize; i++) {
            
            for (int j = 0; j < this.particleSize2; j++) {
                tmpDeltaX = tmpCentre1[i][0] - tmpCentre2[i][j][0];
                tmpDeltaY = tmpCentre1[i][1] - tmpCentre2[i][j][1];
                tmpDeltaZ = tmpCentre1[i][2] - tmpCentre2[i][j][2];
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
                tmpDistance[i][j] = Math.sqrt(tmpDeltaX * tmpDeltaX +
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
        
        for (int i = 0; i < this.atomSize1; i++) {
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
        
        for (int i = 0; i < this.atomSize2; i++) {
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

            for (int i = 0; i < this.atomSize1; i++) {
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
            int tmpSkip = this.atomSize1 + 1;
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

            for (int i = 0; i < this.atomSize2; i++) {
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
        
        double tmpDeltaX;
        double tmpDeltaY;
        double tmpDeltaZ;
        
        this.distances = 
                new double[this.atomSize1][this.atomSize2];
        
        for (int i = 0; i < this.atomSize1; i++) {
                
            for (int j = 0; j < this.atomSize2; j++) {
                tmpDeltaX = this.coordinateList2[0][0][j][0] 
                        - this.coordinateList1[0][i][0];
                tmpDeltaY = this.coordinateList2[0][0][j][1] 
                        - this.coordinateList1[0][i][1];
                tmpDeltaZ = this.coordinateList2[0][0][j][2] 
                        - this.coordinateList1[0][i][2];
                this.distances[i][j] = Math.sqrt(tmpDeltaX * tmpDeltaX
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
     *   Tinker xyz filename
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
        double tmpX;
        double tmpY;
        double tmpZ;
        String tmpReadLine;
        String[] tmpReadLineArray;
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
                tmpReadLineArray = tmpReadLine.trim().split("\\s+");
                tmpAtomNumber = Integer.parseInt(tmpReadLineArray[0]);
                this.atomNumber = tmpAtomNumber;
                this.header = "";
                if (tmpReadLineArray != null && tmpReadLineArray.length > 1) {

                    for (int i = 1; i < tmpReadLineArray.length; i++) {
                        this.header +=  tmpReadLineArray[i];
                    }

                }
                if (anAtomSize2 == 0) {
                    this.atomSize1 = tmpAtomNumber;
                    this.atomSize2 = 0;
                    this.particleSize2 = 0;
                } else {
                    this.atomSize1 = anAtomSize1;
                    this.atomSize2 = anAtomSize2;
                    this.particleSize2 = (tmpAtomNumber - anAtomSize1) 
                            / anAtomSize2;
                    this.coordinateList2 = new double[anIterationSize]
                            [this.particleSize2][this.atomSize2][3];
                }
                this.coordinateList1 = 
                        new double[anIterationSize][this.atomSize1][3];
            
                // read second line
                tmpBR.mark(80);
                tmpReadLine = tmpBR.readLine();
                tmpReadLineArray = tmpReadLine.trim().split("\\s+");
                if (!tmpReadLineArray[0].equals("1")) {
                    tmpHasCommentLine = true;
                    this.comment = tmpReadLine;
                    this.fileContent.append(tmpReadLine);
                    this.fileContent.append(this.LINESEPARATOR);
                } else {
                    tmpBR.reset();
                }
                this.elementList1 = new String[this.atomSize1];
                this.elementList2 = new String[anAtomSize2];
                this.atomTypeList1 = new int[this.atomSize1];
                this.atomTypeList2 = new int[anAtomSize2]; 
                this.connectionList1 = new int[this.atomSize1][];
                this.connectionList2 = new int[anAtomSize2][];

                for (int i = 0; i < anIterationSize; i++) {

                    for (int j = 0; j < this.atomSize1; j++) {
                        tmpReadLine = tmpBR.readLine();
                        this.fileContent.append(tmpReadLine);
                        this.fileContent.append(this.LINESEPARATOR);
                        tmpReadLineArray = tmpReadLine.trim().split("\\s+");
                        tmpX = Double.parseDouble(tmpReadLineArray[2]);
                        tmpY = Double.parseDouble(tmpReadLineArray[3]);
                        tmpZ = Double.parseDouble(tmpReadLineArray[4]);
                        this.coordinateList1[i][j] = new double[3];
                        this.coordinateList1[i][j][0] = tmpX;
                        this.coordinateList1[i][j][1] = tmpY;
                        this.coordinateList1[i][j][2] = tmpZ;
                        tmpConnectionSize = tmpReadLineArray.length - 6;
                        if (i == 0) {
                            this.elementList1[j] = tmpReadLineArray[1];
                            this.atomTypeList1[j] = Integer.
                                    parseInt(tmpReadLineArray[5]);
                            this.connectionList1[j] = new int[tmpConnectionSize];

                            for (int k = 0; k < tmpConnectionSize; k++) {
                                this.connectionList1[j][k] = Integer
                                        .parseInt(tmpReadLineArray[k + 6]);
                            }
                        }
                    }

                    // read following lines of second particle
                    for (int j = 0; j < this.particleSize2; j++) {

                        for (int k = 0; k < this.atomSize2; k++) {
                            tmpReadLine = tmpBR.readLine();
                            this.fileContent.append(tmpReadLine);
                            this.fileContent.append(this.LINESEPARATOR);
                            tmpReadLineArray = tmpReadLine.trim().split("\\s+");
                            tmpX = Double.parseDouble(tmpReadLineArray[2]);
                            tmpY = Double.parseDouble(tmpReadLineArray[3]);
                            tmpZ = Double.parseDouble(tmpReadLineArray[4]);
                            this.coordinateList2[i][j][k] = new double[3];
                            this.coordinateList2[i][j][k][0] = tmpX;
                            this.coordinateList2[i][j][k][1] = tmpY;
                            this.coordinateList2[i][j][k][2] = tmpZ;
                            tmpConnectionSize = tmpReadLineArray.length - 6;
                            if(i == 0 && j == 0) {
                                this.elementList2[k] = tmpReadLineArray[1];
                                this.atomTypeList2[k] = Integer.
                                        parseInt(tmpReadLineArray[5]);
                                this.connectionList2[k] = 
                                        new int[tmpConnectionSize];
                                for (int l = 0; l < tmpConnectionSize; l++) {
                                        this.connectionList2[k][l] = Integer
                                                .parseInt(tmpReadLineArray[l + 6]);
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
        
        for (int i = 0; i < this.atomSize1; i++) {
            if (!this.atomicMassMap.containsKey(this.elementList1[i])) {
                tmpAtomicMass = MIPETUTIL.getAtomicMass(
                        this.elementList1[i], false);
                this.atomicMassMap.put(this.elementList1[i], tmpAtomicMass);
            }
        }
        
        for (int i = 0; i < this.atomSize2; i++) {
            if (!this.atomicMassMap.containsKey(this.elementList2[i])) {
                tmpAtomicMass = MIPETUTIL.getAtomicMass(
                        this.elementList2[i], false);
                this.atomicMassMap.put(this.elementList2[i], tmpAtomicMass);
            }
        }
        
        // Set atomic mass list
        this.atomicMassList1 = this.getAtomicMassList1();
        if (this.atomSize2 > 0) {
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
        tmpAtomSize1 = tmpTinkerXyz1.getAtomNumber();
        tmpAtomSize2 = tmpTinkerXyz2.getAtomNumber();
        this.forcefieldName = tmpTinkerXyz1.getForcefieldName();
        this.particleName1 = tmpTinkerXyz1.getParticleName1();
        this.particleName2 = tmpTinkerXyz2.getParticleName1();
        this.atomNumber = tmpAtomSize1 + tmpAtomSize2;
        this.atomSize1 = tmpAtomSize1;
        this.atomSize2 = tmpAtomSize2;
        this.particleSize2 = (this.atomNumber - this.atomSize1) 
                            / this.atomSize2;
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
        int[][] tmpIntegerList = new int[aIntegerList.length][];
        
        for (int i = 0; i < aIntegerList.length; i++) {
            tmpIntegerList[i] = new int[aIntegerList[i].length];
            
            for (int j = 0; j < aIntegerList[i].length; j++) {
                tmpIntegerList[i][j] = aIntegerList[i][j] + this.atomSize1;
            }
            
        }
        
        return tmpIntegerList;
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
        
        StringBuilder tmpContent;
        int tmpIndex;
        int tmpConnectionSize;
        DecimalFormat decimal6;
        
        tmpContent = new StringBuilder();
        decimal6 = (DecimalFormat)NumberFormat.getNumberInstance();
        decimal6.applyPattern("#0.000000");
        
        // Fill the content from TinkerXYZ object
        tmpContent.append(MIPETUTIL.padLeft(
                String.valueOf(this.atomNumber), 9));
        tmpContent.append(" ");
        tmpContent.append(this.header);
        tmpContent.append(LINESEPARATOR);
        tmpIndex = 1;
        
        for (int i = 0; i < this.atomSize1; i++) {
            tmpContent.append(MIPETUTIL.padLeft(String.valueOf(tmpIndex), 9));
            tmpIndex ++;
            tmpContent.append("   ");
            tmpContent.append(MIPETUTIL.padRight(this.elementList1[i], 3));
            tmpContent.append(MIPETUTIL.padLeft(decimal6.format(
                    this.coordinateList1[0][i][0]), 12));
            tmpContent.append(MIPETUTIL.padLeft(decimal6.format(
                    this.coordinateList1[0][i][1]), 12));
            tmpContent.append(MIPETUTIL.padLeft(decimal6.format(
                    this.coordinateList1[0][i][2]), 12));
            tmpContent.append(MIPETUTIL.padLeft(Integer.toString(this.atomTypeList1[i]), 6));
            tmpConnectionSize = this.connectionList1[i].length;
            
            for (int j = 0; j < tmpConnectionSize; j++) {
                tmpContent.append(MIPETUTIL.padLeft(
                    Integer.toString(this.connectionList1[i][j]), 6));
            }
            
            if (i < this.atomSize1) {
                tmpContent.append(LINESEPARATOR);
            }
        }
        
        for (int i = 0; i < atomSize2; i++) {
            tmpContent.append(MIPETUTIL.padLeft(String.valueOf(tmpIndex), 9));
            tmpIndex++;
            tmpContent.append("   ");
            tmpContent.append(MIPETUTIL.padRight(this.elementList2[i], 3));
            tmpContent.append(MIPETUTIL.padLeft(decimal6.format(
                    this.coordinateList2[0][0][i][0]), 12));
            tmpContent.append(MIPETUTIL.padLeft(decimal6.format(
                    this.coordinateList2[0][0][i][1]), 12));
            tmpContent.append(MIPETUTIL.padLeft(decimal6.format(
                    this.coordinateList2[0][0][i][2]), 12));
            tmpContent.append(MIPETUTIL.padLeft(Integer.toString(this.atomTypeList2[i]), 6));
            tmpConnectionSize = this.connectionList2[i].length;
            
            for (int j = 0; j < tmpConnectionSize; j++) {
                tmpContent.append(MIPETUTIL.padLeft(
                    Integer.toString(this.connectionList2[i][j]), 6));
            }
            
            if (i < this.atomSize2 - 1) {
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
        TinkerXYZ tmpTinkerXYZ;
        
        try {
            tmpTinkerXYZ = (TinkerXYZ)super.clone();
            if (tmpTinkerXYZ == null) {
                throw new InternalError("TinkerXYZ is null.");
            }
            tmpTinkerXYZ.header = this.header;
            tmpTinkerXYZ.forcefieldName = this.forcefieldName;
            tmpTinkerXYZ.particleName1 = this.particleName1;
            tmpTinkerXYZ.atomNumber = this.atomNumber;
            tmpTinkerXYZ.atomSize1 = this.atomSize1;
            tmpTinkerXYZ.atomicMassList1 = this.atomicMassList1.clone();
            tmpTinkerXYZ.coordinateList1 = this.coordinateList1.clone();
            tmpTinkerXYZ.atomTypeList1 = this.atomTypeList1.clone();
            tmpTinkerXYZ.connectionList1 = this.connectionList1.clone();
            if (atomSize2 > 0) {
                tmpTinkerXYZ.particleName2 = this.particleName2;
                tmpTinkerXYZ.atomSize2 = this.atomSize2;
                tmpTinkerXYZ.atomicMassList2 = this.atomicMassList2.clone();
                tmpTinkerXYZ.coordinateList2 = this.coordinateList2.clone();
                tmpTinkerXYZ.atomTypeList2 = this.atomTypeList2.clone();
                tmpTinkerXYZ.connectionList2 = this.connectionList2.clone();
                tmpTinkerXYZ.distances = this.distances.clone();
            }
            return tmpTinkerXYZ;
        } catch (final CloneNotSupportedException ex) {
            throw new InternalError(ex.getMessage());
        }
    }
    
    // </editor-fold>
}
