/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package de.whs.ibci.mipet;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The RotationModel class provides methods to calculate rotation matrices
 * and import 3D coordinates from text files in the resource folder.
 *
 * @author Veit Hucke, Felix Baensch
 */
public class RotationUtil {

    //<editor-fold desc="Private final class constants" defaultstate="collapsed">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(RotationUtil.class.getName());
    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Reads the input file (only .txt) and saves the coordinates as a list of arrays
     *
     * @param aFilePath Path of the file
     * @return A List of all sphere node coordinates, saved as double arrays
     */
    public static List<double[]> readSphereNodes(String aFilePath) {
        if (aFilePath == null || aFilePath.isEmpty()) {
            throw new NullPointerException("Null was passed to readSphereNodes");
        }
        
        List<double[]> tmpCoordinateList = new ArrayList<>();
        // Use this in the development phase
        ClassLoader tmpClassLoader = RotationUtil.class.getClassLoader();
        File tmpFile = new File(tmpClassLoader.getResource(aFilePath)
                .getFile());
        try (BufferedReader tmpReader = new BufferedReader(
                new FileReader(tmpFile))) {
         
        // Use this for the deployment
//        InputStream tmpInputStream = RotationUtil.class
//                .getResourceAsStream(aFilePath);
//        try (BufferedReader tmpReader = new BufferedReader(
//                new InputStreamReader(tmpInputStream))) {
            
            String tmpNewSubString;
            String tmpLine;
            String[] tmpStringTokens;
            int tmpStartIndex;
            int tmpEndIndex;
            int tmpIndex;
            double tmpTokenAsDouble;
            double[] tmpParsedTokens;
            
            while ((tmpLine = tmpReader.readLine()) != null) {
                tmpStartIndex = tmpLine.indexOf("{");
                tmpEndIndex = tmpLine.indexOf("}");
                tmpNewSubString = tmpLine.substring(tmpStartIndex + 1, 
                        tmpEndIndex);
                tmpStringTokens = tmpNewSubString.split(",\\s+");
                tmpParsedTokens = new double[3];
                tmpIndex = 0;
                for (String tmpToken : tmpStringTokens) {
                    try {
                        String tmpNewToken = tmpToken
                                .replaceAll("\\*10\\^", "E");
                        tmpTokenAsDouble = Double.parseDouble(tmpNewToken);
                        tmpParsedTokens[tmpIndex] = tmpTokenAsDouble;
                        tmpIndex++;
                    } catch (NullPointerException | NumberFormatException 
                            anException) {
                        RotationUtil.LOGGER.log(Level.SEVERE,
                                "String is null or could not be parsed.",
                                anException);
                    }
                }
                tmpCoordinateList.add(tmpParsedTokens);
            }
            tmpReader.close();
            return tmpCoordinateList;
        } catch(IOException | NullPointerException | IllegalArgumentException 
                anException) {
            RotationUtil.LOGGER.log(Level.SEVERE,
                    "File could not be read, FilePath is Null or array could " 
                    + "not be added to the list.",
                    anException);
            return null;
        }
    }
    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="Properties">
    /**
     * Calculates the rotation matrix to rotate around an axis by given radians.
     *
     * @param anAngle Angle to rotate by (in radians).
     * @param aRotationAxis Axis to rotate around.
     * @return The rotation matrix.
     */
    public static double[][] getRotationMatrix(double anAngle, double[] aRotationAxis) {
        try {
            double[] tmpUnitAxis = VectorUtil.normalizeVector(aRotationAxis);
            double[][] tmpDyadic = VectorUtil.dyadicProduct(tmpUnitAxis);
            double[][] tmpSkewSymmetricMatrix = VectorUtil.generateSkewSymmetricMatrix(tmpUnitAxis);
            return MatrixUtil.add(
                    (MatrixUtil.multiplyWithScalar(Math.cos(anAngle), MatrixUtil.UNIT_MATRIX_3D)),
                    (MatrixUtil.multiplyWithScalar((1 - Math.cos(anAngle)), tmpDyadic)),
                    (MatrixUtil.multiplyWithScalar(Math.sin(anAngle), tmpSkewSymmetricMatrix))
            );
        } catch(IllegalArgumentException anException) {
            RotationUtil.LOGGER.log(Level.SEVERE, "Matrices could not be added.", anException);
            return null;
        }
    }

    /**
     * Calculates the rotation matrix to map vector1 on vector2.
     *
     * @param aVector1 Vector to map onto another vector.
     * @param aVector2 Vector to be mapped on.
     * @return The rotation matrix.
     */
    public static double[][] getRotationMatrix(double[] aVector1, double[] aVector2) {
        try {
            double[] tmpUnitVector1 = VectorUtil.normalizeVector(aVector1);
            double[] tmpUnitVector2 = VectorUtil.normalizeVector(aVector2);
            double tmpAngle = Math.acos(VectorUtil.dotProduct(tmpUnitVector1, tmpUnitVector2)); //Angle between v1 and v2, only for unit vectors!
            double[] tmpRotationAxis = VectorUtil.crossProduct(aVector1, aVector2); //Axis to rotate around (vector n is perpendicular to v1 and v2)
            return RotationUtil.getRotationMatrix(tmpAngle, tmpRotationAxis);
        } catch(IllegalArgumentException anException) {
            RotationUtil.LOGGER.log(Level.SEVERE, "Vectors may not have the same dimension.", anException);
            return null;
        }
    }

    /**
     * Gets the matrices to perform all necessary rotations.
     *
     * @param sphereNodeCoordinates Coordinates of the sphere nodes.
     * @param aVector Vector to map the sphere coordinates on.
     * @return List of all rotation matrices.
     */
    public static List<double[][]> getRotationMatrices1(List<double[]> sphereNodeCoordinates, double[] aVector) {
        List<double[][]> tmpRotationMatrixList = new ArrayList<>();
        for(double[] anElement : sphereNodeCoordinates) {
            tmpRotationMatrixList.add(RotationUtil.getRotationMatrix(anElement, aVector));
        }
        return tmpRotationMatrixList;
    }

    /**
     * Gets the matrices to perform all necessary rotations. Here, the sphere additionally rotates around the input vector.
     *
     * @param sphereNodeCoordinates Coordinates of the sphere nodes.
     * @param aVector Vector to map the sphere coordinates on.
     * @param aRotationNumber Number of equal rotations around itself.
     * @return List of all rotation matrices.
     */
    public static List<double[][]> getRotationMatrices2(List<double[]> sphereNodeCoordinates, double[] aVector, int aRotationNumber){
        List<double[][]> tmpRotationMatrixList = new ArrayList<>();
        try {
            int tmpSphereNodeNumber = sphereNodeCoordinates.size();
            for (int i = 0; i < tmpSphereNodeNumber; i++) {
                for (int j = 0; j < aRotationNumber; j++) {
                    tmpRotationMatrixList.add(
                            MatrixUtil.multiply(
                                    RotationUtil.getRotationMatrix(((2 * Math.PI * (j + 1)) / aRotationNumber), aVector),
                                    RotationUtil.getRotationMatrix(sphereNodeCoordinates.get(i), aVector)
                            )
                    );
                }
            }
            return tmpRotationMatrixList;
        } catch(IllegalArgumentException | NullPointerException anException){
            RotationUtil.LOGGER.log(Level.SEVERE,"Matrices could not be multiplied." ,anException);
            return null;
        }
    }
    //</editor-fold>
}

