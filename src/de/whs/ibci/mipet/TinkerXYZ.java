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
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.StringReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

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
    private final int STRINGBUILDER_CAPACITY = 4096;
    
    /**
     * Line separator
     */
    private final String LINESEPARATOR = System.lineSeparator();
    
    /**
     * Atomic mass of elements
     */
    private final Map<String, Double> atomicMassMap = new HashMap<>();

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
    private int N_atom;
    
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
        this.initialize2(aTinkerXYZ1, aTinkerXYZ2);
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
        return this.N_atom;
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
        if (this.atomicMassList1 != null) {
            return this.atomicMassList1;
        }

        this.atomicMassList1 = new double[this.N_atom1];
        
        for (int i = 0; i < this.N_atom1; i++) {
            String element = this.elementList1[i];

            // Check of pseudo atoms
            if (element.equalsIgnoreCase("Lp") || 
                    element.equalsIgnoreCase("M")) {
                this.atomicMassList1[i] = 0.0;
            } else {
                Double mass = this.atomicMassMap.get(element);
                this.atomicMassList1[i] = (mass != null) ? mass : 0.0;
            }
        }
        return this.atomicMassList1;
    }
    
    /** 
     * Returns atomic masses of the second particle
     * 
     * @return Atomic masses of second particle
     */
    public double[] getAtomicMassList2() {
        if (this.atomicMassList2 != null) {
            return this.atomicMassList2;
        }
        
        this.atomicMassList2 = new double[this.N_atom2];
        
        for (int i = 0; i < this.N_atom2; i++) {
            String element = this.elementList2[i];

            // Check of pseudo atoms
            if (element.equalsIgnoreCase("Lp") || 
                    element.equalsIgnoreCase("M")) {
                this.atomicMassList2[i] = 0.0;
            } else {
                Double mass = this.atomicMassMap.get(element);
                this.atomicMassList2[i] = (mass != null) ? mass : 0.0;
            }
        }
        return this.atomicMassList2;
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
     * Returns .txyz file content as String
     *  It prevents manipulation from other classes.
     * 
     * @return FileContent as String
     */
    public synchronized String getFileContent() {
        StringBuilder sb = new StringBuilder(STRINGBUILDER_CAPACITY);

        // Header line (number of atoms + header)
        sb.append(MIPETUTIL.padLeft(Integer.toString(this.N_atom), 6))
          .append("  ")
          .append(this.header != null ? this.header : "")
          .append(LINESEPARATOR);

        // Comment line
        if (this.comment != null && !this.comment.isEmpty()) {
            sb.append(this.comment).append(LINESEPARATOR);
        }

        // Atoms of the first particle
        for (int i = 0; i < this.N_atom1; i++) {
            appendAtomLine(sb, i, this.elementList1[i], 
                           this.coordinateList1[0][i], 
                           this.atomTypeList1[i], 
                           this.connectionList1[i]);
        }
        
        // Atoms of the rest particles
        for (int p = 0; p < this.N_particle2; p++) {
            
            for (int i = 0; i < this.N_atom2; i++) {
                appendAtomLine(sb, i + this.N_atom1 + (p * this.N_atom2), 
                               this.elementList2[i], 
                               this.coordinateList2[0][p][i], 
                               this.atomTypeList2[i], 
                               this.connectionList2[i]);
            }
            
        }
        
        return sb.toString();
    }
    
    
    
    /**
     * Set the header (also in fileContent)
     *  Synchronized prevents race conditions in replace operations
     * 
     * @param aHeader A new Header
     */
    public synchronized void setHeader(String aHeader) {
        this.header = aHeader;
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
        // Check
        if (this.atomTypeList1 == null) {
            return;
        }
        
        int tmpAtomSize1 = this.N_atom1;
        int[] tmpAtomTypes =new int[tmpAtomSize1];
        
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
        
        int atomNumber;
        Reader stringReader = new StringReader(aTinkerXyz);
        boolean hasCommentLine = false;
        
        try (BufferedReader tmpBR = new BufferedReader(stringReader)) {
            
            // read first line
            String tmpReadLine = tmpBR.readLine();
            if (tmpReadLine != null) {
                String[] words = tmpReadLine.trim().split("\\s+", 2);
                atomNumber = Integer.parseInt(words[0]);
                this.N_atom = atomNumber;
                this.header = "";
                if (words != null && words.length > 1) {
                    this.header = words[1];
                }
                if (anAtomSize2 == 0) {
                    this.N_atom1 = atomNumber;
                    this.N_atom2 = 0;
                    this.N_particle2 = 0;
                } else {
                    this.N_atom1 = anAtomSize1;
                    this.N_atom2 = anAtomSize2;
                    this.N_particle2 = (atomNumber - anAtomSize1) 
                            / anAtomSize2;
                    this.coordinateList2 = new double[anIterationSize]
                            [this.N_particle2][this.N_atom2][3];
                }
                this.coordinateList1 = 
                        new double[anIterationSize][this.N_atom1][3];
            
                // read second line
                tmpBR.mark(80);
                tmpReadLine = tmpBR.readLine();
                words = tmpReadLine.trim().split("\\s+");
                if (!words[0].equals("1")) {
                    hasCommentLine = true;
                    this.comment = tmpReadLine;
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
                            words = tmpReadLine.trim().split("\\s+");
                            double x = Double.parseDouble(words[2]);
                            double y = Double.parseDouble(words[3]);
                            double z = Double.parseDouble(words[4]);
                            this.coordinateList2[i][j][k] = 
                                    new double[]{x, y, z};
                            int connectionSize = words.length - 6;
                            if(i == 0 && j == 0) {
                                this.elementList2[k] = words[1];
                                this.atomTypeList2[k] = Integer.
                                        parseInt(words[5]);
                                this.connectionList2[k] = 
                                        new int[connectionSize];
                                for (int l = 0; l < connectionSize; l++) {
                                        this.connectionList2[k][l] = Integer
                                                .parseInt(words[l + 6]);
                                }
                                
                            }
                        }
                        
                    }

                    if (i < anIterationSize) {
                        tmpBR.readLine();
                        if (hasCommentLine) {
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
    private void initialize2(TinkerXYZ aTinkerXYZ1, TinkerXYZ aTinkerXYZ2) {
        
        // Check parameters
        if (aTinkerXYZ1 == null) {
            throw new IllegalArgumentException("aTinkerXYZ1 is null.");
        } else if (aTinkerXYZ2 == null) {
            throw new IllegalArgumentException("aTinkerXYZ2 is null.");
        }
        
        TinkerXYZ tinkerXyz1 = aTinkerXYZ1;
        TinkerXYZ tinkerXyz2 = aTinkerXYZ2;
        int atomSize1 = tinkerXyz1.getN_atom();
        int atomSize2 = tinkerXyz2.getN_atom();
        this.forcefieldName = tinkerXyz1.getForcefieldName();
        this.particleName1 = tinkerXyz1.getParticleName1();
        this.particleName2 = tinkerXyz2.getParticleName1();
        this.N_atom = atomSize1 + atomSize2;
        this.N_atom1 = atomSize1;
        this.N_atom2 = atomSize2;
        this.N_particle2 = (this.N_atom - this.N_atom1) 
                            / this.N_atom2;
        this.elementList1 = tinkerXyz1.getElementList1().clone();
        this.elementList2 = tinkerXyz2.getElementList1().clone();
        this.atomicMassList1 = tinkerXyz1.getAtomicMassList1().clone();
        this.atomicMassList2 = tinkerXyz2.getAtomicMassList2().clone();
        this.atomTypeList1 = tinkerXyz1.getAtomTypeList1().clone();
        this.atomTypeList2 = tinkerXyz2.getAtomTypeList1().clone();
        this.connectionList1 = tinkerXyz1.getConnectionList1().clone();
        this.connectionList2 = tinkerXyz2.getConnectionList1().clone();
        this.connectionList2 = this.correctConnectionList(this.connectionList2);
        this.coordinateList1 = new double[1][atomSize1][3];
        this.coordinateList1 = tinkerXyz1.coordinateList1.clone();
        this.coordinateList2 = new double[1][1][atomSize2][3];
        this.coordinateList2[0] = tinkerXyz2.coordinateList1.clone();
    }
    
    /**
     * Correct connection list for the second particle
     * 
     * @param aIntegerList Old connection list
     * @return Corrected connection list
     */
    private int[][] correctConnectionList(int[][] aIntegerList) {
        int[][] result = new int[aIntegerList.length][];
        
        for (int i = 0; i < aIntegerList.length; i++) {
            result[i] = new int[aIntegerList[i].length];
            
            for (int j = 0; j < aIntegerList[i].length; j++) {
                result[i][j] = aIntegerList[i][j] + this.N_atom1;
            }
            
        }
        
        return result;
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
        double[][] copy = new double[AnOriginal.length][];

        for (int i = 0; i < AnOriginal.length; i++) {
            if (AnOriginal[i] != null) {
                copy[i] = AnOriginal[i].clone();
            }
        }

        return copy;
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

        double[][][] copy = new double[AnOriginal.length][][];

        for (int i = 0; i < AnOriginal.length; i++) {
            if (AnOriginal[i] != null) {
                copy[i] = new double[AnOriginal[i].length][];

                for (int j = 0; j < AnOriginal[i].length; j++) {
                    if (AnOriginal[i][j] != null) {
                        copy[i][j] = AnOriginal[i][j].clone();
                    }
                }
                
            }
        }

        return copy;
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

        double[][][][] copy = new double[AnOriginal.length][][][];

        for (int i = 0; i < AnOriginal.length; i++) {
            if (AnOriginal[i] != null) {
                copy[i] = new double[AnOriginal[i].length][][];

                for (int j = 0; j < AnOriginal[i].length; j++) {
                    if (AnOriginal[i][j] != null) {
                        copy[i][j] = AnOriginal[i][j].clone();
                    }
                    
                    for (int k = 0; k < AnOriginal[i][j].length; k++) {
                        if (AnOriginal[i][j][k] != null) {
                            copy[i][j][k] = AnOriginal[i][j][k].clone();
                        }
                    }
                    
                }
                
            }
        }

        return copy;
    }
    
    /**
     * Helping method to formating of tinkerXYZ lines.
     * 
     * @param sb StringBuilder object
     * @param index index number of elements
     * @param element Element name
     * @param coords xyz coordinates
     * @param type Atom type number
     * @param connections Connection index
     */
    private void appendAtomLine(StringBuilder sb, int index, String element,
            double[] coords, int type, int[] connections) {
        sb.append(MIPETUTIL.padLeft(Integer.toString(index + 1), 6))
                .append("  ")
                .append(MIPETUTIL.padRight(element, 3))
                .append(String.format(Locale.ENGLISH, "%12.6f%12.6f%12.6f", 
                        coords[0], coords[1], coords[2]))
                .append(MIPETUTIL.padLeft(Integer.toString(type), 6));
        if (connections != null) {
            for (int conn : connections) {
                sb.append(MIPETUTIL.padLeft(Integer.toString(conn), 6));
            }
        }
        sb.append(LINESEPARATOR);
    }
    
    /**
     * Helping method to formating of xyz-file lines.
     * 
     * @param sb StringBuilder object
     * @param element Element name
     * @param coords xyz coordinates
     */
    private void appendXyzLine(StringBuilder sb, String element, 
            double[] coords) {
        sb.append(MIPETUTIL.padRight(element, 3))
                .append(String.format(Locale.ENGLISH, "%12.6f%12.6f%12.6f", 
                        coords[0], coords[1], coords[2]));
        sb.append(LINESEPARATOR);
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
        
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        double[] centreCoord = new double[3];
        
        for (double[] coord : aCoords) {
            x += coord[0];
            y += coord[1];
            z += coord[2];
        }
        centreCoord[0] = x / aCoords.length; 
        centreCoord[1] = y / aCoords.length;
        centreCoord[2] = z / aCoords.length;
        return centreCoord;
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
                
        double[][] coordinates = new double[aCoords.length][3];

        for (int i = 0; i < aCoords.length; i++) {
            coordinates[i][0] = aCoords[i][0] - aVector[0];
            coordinates[i][1] = aCoords[i][1] - aVector[1];
            coordinates[i][2] = aCoords[i][2] - aVector[2];
        }
        
        return coordinates;
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
        
        // Make a local StringBuilder only for writing
        StringBuilder sb = new StringBuilder(STRINGBUILDER_CAPACITY);

        // Header for .xyz
        sb.append(this.N_atom).append("    ")
                .append(this.header)
                .append(LINESEPARATOR);

        // Comment line
        if (this.comment != null) {
            sb.append(this.comment).append(LINESEPARATOR);
        }
        
        // Body
        for (int i = 0; i < this.N_atom1; i++) {
            this.appendXyzLine(sb, this.elementList1[i], 
                    this.coordinateList1[0][i]);
        }
        
        // Particle 2
        for (int p = 0; p < this.N_particle2; p++) {
            
            for (int i = 0; i < this.N_atom2; i++) {
                this.appendXyzLine(sb, this.elementList2[i], 
                        this.coordinateList2[0][p][i]);
            }
            
        }
        
        // Write to file
        try (PrintWriter out = new PrintWriter(new FileWriter(aXyzFileName))) {
            out.print(sb.toString());
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, "Error writing .xyz file", ex);
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
        Path xyzPath = Paths.get(aXyzFileName);
        try (BufferedWriter writer = Files.newBufferedWriter(xyzPath, 
                StandardCharsets.UTF_8)) {
            writer.append(aData);
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
        
        //this.writeToXyzFile(aTxyzFileName, this.fileContent);
        try (PrintWriter writer = new PrintWriter(aTxyzFileName)) {
            writer.print(this.getFileContent());
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, "IOException during writing .txyz file.", ex);
        }
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
        ArrayList<double[][]> coordsList = new ArrayList<>();
        Path path = Paths.get(aFileName);
        try (BufferedReader reader = Files.newBufferedReader(path, 
                StandardCharsets.UTF_8)){
            String line = reader.readLine();
            if (line == null) {
                return new double[0][0][0];
            }
            int tmpAtomSize = Integer.parseInt(line.trim()
                    .split("\\s+")[0]);

            while (line != null) {
                // Header line
                double[][] currSnapshot = new double[tmpAtomSize][3];

                for (int i = 0; i < tmpAtomSize; i++) {
                    line = reader.readLine();
                    if (line == null) break;

                    String[] tokens = line.trim().split("\\s+");
                    // XYZ block Columns: 0:index, 1:element, 2:X, 3:Y, 4:Z
                    currSnapshot[i][0] = Double.parseDouble(tokens[2]);
                    currSnapshot[i][1] = Double.parseDouble(tokens[3]);
                    currSnapshot[i][2] = Double.parseDouble(tokens[4]);
                }

                coordsList.add(currSnapshot);

                // Read the next Header line for the next snapshot
                line = reader.readLine(); 
            }
        } catch (IOException | NumberFormatException ex) {
            LOGGER.log(Level.SEVERE, "Error reading .arc file: " + aFileName, ex);
        }

        // Converting of the list to the final 3d-array
        return coordsList.toArray(double[][][]::new);
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
        int lineNumber = this.N_atom1 + this.N_atom2 + 5;
        StringBuilder content = new StringBuilder(lineNumber * 80);
        DecimalFormat decimal6 = (DecimalFormat)NumberFormat
                .getNumberInstance();
        decimal6.applyPattern("#0.000000");
        
        // Fill the content from TinkerXYZ object
        content.append(MIPETUTIL
                .padLeft(String.valueOf(this.N_atom), 9))
                .append(" ")
                .append(this.header)
                .append(LINESEPARATOR);
        int tmpIndex = 1;
        
        for (int i = 0; i < this.N_atom1; i++) {
            double[] currCoords = this.coordinateList1[0][i];
            int[] currConnections = this.connectionList1[i];
            
            content.append(MIPETUTIL.padLeft(String.valueOf(tmpIndex), 9))
                    .append("   ")
                    .append(MIPETUTIL.padRight(this.elementList1[i], 3))
                    .append(MIPETUTIL.padLeft(decimal6.format(
                            currCoords[0]), 12))
                    .append(MIPETUTIL.padLeft(decimal6.format(
                            currCoords[1]), 12))
                    .append(MIPETUTIL.padLeft(decimal6.format(
                            currCoords[2]), 12))
                    .append(MIPETUTIL.padLeft(
                            Integer.toString(this.atomTypeList1[i]), 6));
            
            tmpIndex ++;
            int connectionSize = currConnections.length;
            
            for (int j = 0; j < connectionSize; j++) {
                content.append(MIPETUTIL.padLeft(
                    Integer.toString(currConnections[j]), 6));
            }
            
            if (i < this.N_atom1) {
                content.append(LINESEPARATOR);
            }
        }
        
        for (int i = 0; i < N_atom2; i++) {
            double[] currCoords = this.coordinateList2[0][0][i];
            int[] currConnections = this.connectionList2[i];
            
            content.append(MIPETUTIL.padLeft(String.valueOf(tmpIndex), 9))
                    .append("   ")
                    .append(MIPETUTIL.padRight(this.elementList2[i], 3))
                    .append(MIPETUTIL.padLeft(decimal6
                            .format(currCoords[0]), 12))
                    .append(MIPETUTIL.padLeft(decimal6
                            .format(currCoords[1]), 12))
                    .append(MIPETUTIL.padLeft(decimal6
                            .format(currCoords[2]), 12))
                    .append(MIPETUTIL.padLeft(Integer
                            .toString(this.atomTypeList2[i]), 6));
            tmpIndex++;
            int connectionSize = currConnections.length;
            
            for (int j = 0; j < connectionSize; j++) {
                content.append(MIPETUTIL.padLeft(
                        Integer.toString(currConnections[j]), 6));
            }
            
            if (i < this.N_atom2 - 1) {
                content.append(LINESEPARATOR);
            }
        }
        
        try (PrintWriter writer = new PrintWriter((aFileName))) {
            writer.print(content);
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
        TinkerXYZ cloned;
        try {
            cloned = (TinkerXYZ) super.clone();
        } catch (final CloneNotSupportedException ex) {
            throw new AssertionError("Cloning failed.", ex);
        }
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
    }
    
    // </editor-fold>
}