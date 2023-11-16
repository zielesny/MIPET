/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package de.whs.ibci.mipet;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.ArrayList;
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
     * Logger of this class
     */
    private static final Logger LOGGER = Logger
            .getLogger(MIPET.class.getName());
    /**
     * Standard capacity for stringbuilder object
     */
    final int STRINGBUILDER_CAPACITY;
    
    /**
     * Standard BufferdReader buffersize 
     */
    final int READER_BUFFERSIZE;
    
    /**
     * Line separator
     */
    final private String LINESEPARATOR;

    //</editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    
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
     * second particle number
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
    private int[] parameterList1;
    
    /**
     * Parameters of specific parameterset defined in *.prm of second particle
     */
    private int[] parameterList2;
    
    /**
     * Connectionlist of first particle 
     */
    private int[][] connectionList1;
    
    /**
     * Connectionlist of second particle 
     */
    private int[][] connectionList2;
    
    /**
     * File content
     */
    private StringBuilder fileContent;
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    
    /**
     * Constructor TinkerXYZ
     * 
     * @param aTxyzFileName Tinker xyz filename
     */
    public TinkerXYZ(String aTxyzFileName) {
        this.STRINGBUILDER_CAPACITY = 4096;
        this.READER_BUFFERSIZE = 65536;
        this.LINESEPARATOR = System.getProperty("line.separator");
        this.initialize(aTxyzFileName, 1, 0, 0);
    }
    
    /**
     * Constructor TinkerXYZ
     * 
     * @param aTxyzFileName Tinker xyz filename
     * @param atomSize Atom number of particle
     */
    public TinkerXYZ(String aTxyzFileName, int atomSize) {
        this.STRINGBUILDER_CAPACITY = 4096;
        this.READER_BUFFERSIZE = 65536;
        this.LINESEPARATOR = System.getProperty("line.separator");
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
        this.STRINGBUILDER_CAPACITY = 4096;
        this.READER_BUFFERSIZE = 65536;
        this.LINESEPARATOR = System.getProperty("line.separator");
        this.initialize(aTxyzFileName, anIterationSize, anAtomSize1, 
                anAtomSize2);
    }
    /**
     * Constructor TinkerXYZ
     * 
     * @param aTinkerXYZ1 TinkerXYZ object of first fragment
     * @param aTinkerXYZ2 TinkerXYZ object of second fragment
     */
    public TinkerXYZ(TinkerXYZ aTinkerXYZ1, TinkerXYZ aTinkerXYZ2) {
        this.STRINGBUILDER_CAPACITY = 4096;
        this.READER_BUFFERSIZE = 65536;
        this.LINESEPARATOR = System.getProperty("line.separator");
        this.initialize2(aTinkerXYZ1, aTinkerXYZ2);
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Public properties">
    
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
     * Set the header (also in fileContent)
     * 
     * @param aHeader A new Header
     */
    public void setHeader(String aHeader) {
        this.header = aHeader;
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
     * Returns parameters of first particle (6. column of .txyz-file)
     * 
     * @return Parameters of first particle
     */
    public int[] getParameterList1() {
        return this.parameterList1;
    }
    
    /**
     * Returns parameters of second particle (6. column of .txyz-file)
     * 
     * @return Parameters of second particle
     */
    public int[] getParameterList2() {
        return this.parameterList2;
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
                    padRightLoopSB(tmpElement, 3));
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
                    padRightLoopSB(tmpElement, 3));
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
     */
    public void setCoordinateList1(double[][][] aCoordinates1) {
        
        // Check parameters
        if (aCoordinates1 == null || aCoordinates1.length == 0) {
            throw new IllegalArgumentException("Null or empty object was passed"
                    + "to the setCoordinateList1 method.");
        }
        
        this.coordinateList1 = aCoordinates1.clone();
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
            tmpStartIndex += 12;
            tmpX = tmpDF.format(aCoordinates1[0][i][0]);
            this.fileContent.replace(tmpStartIndex, tmpStartIndex + 11,
                    padLeftLoopSB(tmpX, 11));
            tmpStartIndex += 12;
            tmpY = tmpDF.format(aCoordinates1[0][i][1]);
            this.fileContent.replace(tmpStartIndex, tmpStartIndex + 11,
                    padLeftLoopSB(tmpY, 11));
            tmpStartIndex += 12;
            tmpZ = tmpDF.format(aCoordinates1[0][i][2]);
            this.fileContent.replace(tmpStartIndex, tmpStartIndex + 11,
                    padLeftLoopSB(tmpZ, 11));
            tmpStartIndex = this.fileContent
                .indexOf(this.LINESEPARATOR, tmpStartIndex);
            tmpStartIndex += this.LINESEPARATOR.length();
        }
    }
    
    /**
     * Overloaded method
     * Set first coordinate list and overwrite stringbuilder
     * 
     * @param aCoord1  [i][j] i: atomid, j: xyz
     */
    public void setCoordinateList1(double[][] aCoord1) {
        
        // Check parameters
        if (aCoord1 == null || aCoord1.length == 0) {
            throw new IllegalArgumentException("Null or empty object was passed"
                    + "to the setCoordinateList1 method.");
        }
        
        double[][][] tmpCoord1 = new double[1][][];
        tmpCoord1[0] = aCoord1;
        setCoordinateList1(tmpCoord1);
    }
    
    /**
     * Set second coordinate list and overwrite stringbuilder
     * 
     * @param aCoord2[i][j][k][l]
     *   i: sim.iteration j: particle k: atom l: xyz
     */
    public void setCoordinateList2(double[][][][] aCoord2) {
        
        // Check parameters
        if (aCoord2 == null || aCoord2.length == 0) {
            throw new IllegalArgumentException("Null or empty object was passed"
                    + "to the setCoordinateList2 method.");
        }
        
        this.coordinateList2 = aCoord2;
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
            tmpStartIndex += 12;
            tmpX = tmpDF.format(aCoord2[0][0][i][0]);
            this.fileContent.replace(tmpStartIndex, tmpStartIndex + 11,
                    padLeftLoopSB(tmpX, 11));
            tmpStartIndex += 12;
            tmpY = tmpDF.format(aCoord2[0][0][i][1]);
            this.fileContent.replace(tmpStartIndex, tmpStartIndex + 11,
                    padLeftLoopSB(tmpY, 11));
            tmpStartIndex += 12;
            tmpZ = tmpDF.format(aCoord2[0][0][i][2]);
            this.fileContent.replace(tmpStartIndex, tmpStartIndex + 11,
                    padLeftLoopSB(tmpZ, 11));
            tmpStartIndex = this.fileContent
                .indexOf(this.LINESEPARATOR, tmpStartIndex);
            tmpStartIndex += this.LINESEPARATOR.length();
        }
    }
    
    /**
     * Overloaded method
     * Set second coordinate list and overwrite stringbuilder
     * 
     * @param aCoord2 
     * [i, j] i: atomid, j: xyz
     */
    public void setCoordinateList2(double[][] aCoord2) {
        
        // Check parameters
        if (aCoord2 == null || aCoord2.length == 0) {
            throw new IllegalArgumentException("Null or empty object was passed"
                    + "to the setCoordinateList2 method.");
        }
        
        double[][][][] tmpCoord2 = new double[1][1][][];
        tmpCoord2[0][0] = aCoord2;
        setCoordinateList2(tmpCoord2);
    }
    
    // </editor-fold>
        
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * initialize mothod
     * 
     * @param aTxyzFileName
     *   Tinker xyz filename
     * @param itrerationSize
     *   number of iteration
     * @param anAtomSize1
     *   atom number of first particle
     * @param anAtomSize2
     *   atom number of second particle
     */
    private void initialize(String aTxyzFileName, int anIterationSize, 
            int anAtomSize1, int anAtomSize2) {
        
        // Check parameters
        if (aTxyzFileName == null || aTxyzFileName.isEmpty()) {
            throw new IllegalArgumentException("Filename is null or empty.");
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
        
        this.fileContent = new StringBuilder(STRINGBUILDER_CAPACITY);
        File tmpSourceFile = new File(aTxyzFileName);
        String tmpReadLine;
        String[] tmpReadLineArray;
        int tmpAtomNumber;
        int tmpConnectionSize;
        Boolean tmpHasCommentLine = false;
        double tmpX;
        double tmpY;
        double tmpZ;
        
        try (BufferedReader tmpBR = new BufferedReader(
                new FileReader(tmpSourceFile), this.READER_BUFFERSIZE)) {
            
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
                this.parameterList1 = new int[this.atomSize1];
                this.parameterList2 = new int[anAtomSize2]; 
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
                            this.parameterList1[j] = Integer.
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
                                this.parameterList2[k] = Integer.
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
            tmpBR.close();
        } catch (FileNotFoundException ex) {
            throw new IllegalArgumentException("File not found in TinkerXYZ.");
        } catch (IOException ex) {
            throw new IllegalArgumentException("IOException in TinkerXYZ.");            
        }
    }
   
    /**
     * Build a new TinkerXYZ object out of two TinkerXYZ objects
     * 
     * @param aTinkerXYZ1
     *   First TinkerXYZ object
     * @param aTinkerXYZ2 
     *   Second TinkerXYZ object
     */
    private void initialize2(TinkerXYZ aTinkerXYZ1, TinkerXYZ aTinkerXYZ2) {
        
        // Check parameters
        if (aTinkerXYZ1 == null) {
            throw new IllegalArgumentException("aTinkerXYZ1 is null.");
        } else if (aTinkerXYZ2 == null) {
            throw new IllegalArgumentException("aTinkerXYZ2 is null.");
        }
        
        String tmpAtomNumber;
        DecimalFormat tmpDF = new DecimalFormat("0.000000", 
        DecimalFormatSymbols.getInstance(Locale.ENGLISH));
        String tmpX;
        String tmpY;
        String tmpZ;
        String tmpParameter;
        int tmpAtomSize1 = aTinkerXYZ1.getAtomNumber();
        int tmpAtomSize2 = aTinkerXYZ2.getAtomNumber();
        int tmpNConnection;
        String tmpConnection;
        
        this.atomNumber = tmpAtomSize1 + tmpAtomSize2;
        this.atomSize1 = aTinkerXYZ1.getAtomNumber();
        this.atomSize2 = aTinkerXYZ2.getAtomNumber();
        this.elementList1 = aTinkerXYZ1.getElementList1();
        this.elementList2 = aTinkerXYZ2.getElementList1();
        this.parameterList1 = aTinkerXYZ1.getParameterList1();
        this.parameterList2 = aTinkerXYZ2.getParameterList2();
        this.connectionList1 = aTinkerXYZ1.getConnectionList1();
        this.connectionList2 = aTinkerXYZ2.getConnectionList2();
        this.fileContent = new StringBuilder(STRINGBUILDER_CAPACITY);
        this.fileContent.append(aTinkerXYZ1.getFileContent());
        this.fileContent.replace(0, 6, 
                padLeftLoopSB(Integer.toString(this.atomNumber) + "  ", 6));
        this.coordinateList1 = new double[1][tmpAtomSize1][3];
        this.coordinateList1 = aTinkerXYZ1.coordinateList1;
        this.coordinateList2 = new double[1][1][tmpAtomSize2][3];
        this.coordinateList2[0] = aTinkerXYZ2.coordinateList1;
        
        for (int i = 0; i < tmpAtomSize2; i++) {
            tmpAtomNumber = Integer.toString(i + tmpAtomSize1 + 1);
            this.fileContent.append(padLeftLoopSB(tmpAtomNumber, 6));
            this.fileContent.append("  ");
            this.fileContent.append(padRightLoopSB(aTinkerXYZ2
                    .getElementList1()[i], 3));
            tmpX = tmpDF.format(aTinkerXYZ2.getCoordinateList1()[0][i][0]);
            this.fileContent.append(padLeftLoopSB(tmpX, 12));
            tmpY = tmpDF.format(aTinkerXYZ2.getCoordinateList1()[0][i][1]);
            this.fileContent.append(padLeftLoopSB(tmpY, 12));
            tmpZ = tmpDF.format(aTinkerXYZ2.getCoordinateList1()[0][i][2]);
            this.fileContent.append(padLeftLoopSB(tmpZ, 12));
            tmpParameter = Integer.toString(aTinkerXYZ2.getParameterList1()[i]);
            this.fileContent.append(padLeftLoopSB(tmpParameter, 7));
            tmpNConnection = aTinkerXYZ2.getConnectionList1()[i].length;
            
            for (int j = 0; j < tmpNConnection; j++) {
                tmpConnection = Integer.toString( 
                        aTinkerXYZ2.getConnectionList1()[i][j] + tmpAtomSize1);
                this.fileContent.append(padLeftLoopSB(tmpConnection, 6));
            }
            
            this.fileContent.append(this.LINESEPARATOR);
        }
        
    }
    
    /**
     * Fills left side of the string with spaces so the string is right aligned
     *   It is much faster version of String.format()
     * @param aInput
     *   Inputstring
     * @param aPadUpTo
     *   Total length of Inputstring and spaces
     * @return 
     *   String with left filled with spaces
     */
    private static String padLeftLoopSB(String aInput, int aPadUpTo) {
        
        // Check parameters
        if (aInput == null || aInput.isEmpty()) {
            throw new IllegalArgumentException("aInput is null or empty.");
        } else if (aPadUpTo <= 0) {
            throw new IllegalArgumentException("aPadUpTo should be positive.");
        }
        
        StringBuilder tmpSB = new StringBuilder();
        char tmpPadChar = ' ';
        
        for (int toPrepend = aPadUpTo - aInput.length(); 
                toPrepend > 0; toPrepend--) {
            tmpSB.append(tmpPadChar);
        }
        
        tmpSB.append(aInput);
        return tmpSB.toString();
    }
    
    /**
     * Fills right side of string with spaces so the string is left aligned
     *   Much faster version of String.format()
     * @param aInput
     *   Inputstring
     * @param aPadUpTo
     *   Total length of Inputstring and spaces
     * @return 
     *   String with right filled with spaces
     */
    private static String padRightLoopSB(String aInput, int aPadUpTo) {
        
        if (aInput == null || aInput.isEmpty()) {
            throw new IllegalArgumentException("aInput is null or empty.");
        } else if (aPadUpTo <= 0) {
            throw new IllegalArgumentException("aPadUpTo should be positive.");
        }
        
        StringBuilder tmpSb = new StringBuilder();
        char tmpPadChar = ' ';
        tmpSb.append(aInput);
        
        for (int toAppend = aInput.length(); 
                toAppend < aPadUpTo; toAppend++) {
            tmpSb.append(tmpPadChar);
        }
        
        return tmpSb.toString();
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
            tmpBW.append(String.format("%d", this.atomNumber)
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
                tmpBW.append(padRightLoopSB(tmpElements[i], 5));
                tmpBW.append(padLeftLoopSB(decimalFormat
                    .format(tmpCoordList[i][0]), 12));
                tmpBW.append(padLeftLoopSB(decimalFormat
                    .format(tmpCoordList[i][1]), 12));
                tmpBW.append(padLeftLoopSB(decimalFormat
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
            throw new IllegalArgumentException("aFileName is null or empty.");
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
            return tmpTinkerXYZ;
        } catch (final CloneNotSupportedException ex) {
            throw new InternalError(ex.getMessage());
        }
    }
    
    // </editor-fold>
}
