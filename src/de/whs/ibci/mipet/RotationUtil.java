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
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The RotationModel class provides methods to calculate rotation matrices
 * and import 3D coordinates from text files in the resource folder.
 *
 * @author Veit Hucke, Felix Baensch, Mirco Daniel
 */
public class RotationUtil {

    //<editor-fold desc="Private final class constants" defaultstate="collapsed">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(RotationUtil.class.getName());
    //</editor-fold>

    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor RotationUtil
     */
    public RotationUtil() {
    }
    
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Reads the input file (only .txt) and saves the coordinates as a list of arrays
     *
     * @param aFilePath Path of the file
     * @return A List of all sphere node coordinates, saved as double arrays
     */
    public static LinkedList<double[]> readSphereNodes(String aFilePath) {
        if (aFilePath == null || aFilePath.isEmpty()) {
            throw new NullPointerException("Null was passed to readSphereNodes");
        }
        
        LinkedList<double[]> tmpCoordinateList = new LinkedList<>();
        // Use this in the development phase
        File tmpFile = new File(aFilePath);
        try (BufferedReader tmpReader = new BufferedReader(
                new FileReader(tmpFile))) {
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
     *  It uses direct calculation (Rodrigues formular)
     * @param angle Angle to rotate by (in radians).
     * @param rotationAxis Axis to rotate around.
     * @return The rotation matrix.
     */
    public static double[][] getRotationMatrix(double angle, 
            double[] rotationAxis) {
        try {
            // Normalize axis
            double[] u = VectorUtil.normalizeVector(rotationAxis);
            double x = u[0];
            double y = u[1];
            double z = u[2];
            
            // Calculate trigonometry
            double cos = Math.cos(angle);
            double sin = Math.sin(angle);
            double oneMinusCos = 1.0 - cos;
            
            // Calculate often used values
            double xyC = x * y * oneMinusCos;
            double xzC = x * z * oneMinusCos;
            double yzC = y * z * oneMinusCos;
            double xSin = x * sin;
            double ySin = y * sin;
            double zSin = z * sin;
            
            return new double[][] {
                {cos + x * x * oneMinusCos, xyC - zSin, xzC + ySin},
                {xyC + zSin, cos + y * y * oneMinusCos, yzC - xSin},
                {xzC - ySin, yzC + xSin, cos + z * z * oneMinusCos}
            };
        } catch(IllegalArgumentException anException) {
            RotationUtil.LOGGER.log(Level.SEVERE, 
                    "Matrices could not be calculated.", anException);
            return null;
        }
    }

    /**
     * Calculates the rotation matrix to map vector1 on vector2.
     *
     * @param vector1 Vector to map onto another vector.
     * @param vector2 Vector to be mapped on.
     * @return The rotation matrix.
     */
    public static double[][] getRotationMatrix(double[] vector1, double[] vector2) {
        if (vector1 == null || vector2 == null || 
                vector1.length != vector2.length) {
            RotationUtil.LOGGER.log(Level.SEVERE, 
                    "Vectors must be non-null and have the same size");
            return null;
        }
        if (vector1.length != 3) {
            RotationUtil.LOGGER.log(Level.SEVERE, 
                    "Vectors must be 3-dimensional for cross product.");
        return null;
        }
        
        double[] unitV1 = VectorUtil.normalizeVector(vector1);
        double[] unitV2 = VectorUtil.normalizeVector(vector2);
        
        // Calculate skalarproduct and clipping (prevent NaN-Bug)
        double dotProduct = VectorUtil.dotProduct(unitV1, unitV2);
        double clampedDot = Math.max(-1.0, Math.min(1.0, dotProduct)); 
        double angle = Math.acos(clampedDot);
        
        // Cross produkt
        double[] rotationAxis = VectorUtil.crossProduct(unitV1, unitV2); 

        return RotationUtil.getRotationMatrix(angle, rotationAxis);
    }

    /**
     * Gets the matrices to perform all necessary rotations.
     *
     * @param sphereNodeCoordinates Coordinates of the sphere nodes.
     * @param aVector Vector to map the sphere coordinates on.
     * @return List of all rotation matrices.
     */
    public static List<double[][]> getRotationMatrices1(
            List<double[]> sphereNodeCoordinates, double[] aVector) {
        LinkedList<double[][]> tmpRotationMatrixList = new LinkedList<>();
        
        for(double[] anElement : sphereNodeCoordinates) {
            tmpRotationMatrixList.add(RotationUtil
                    .getRotationMatrix(anElement, aVector));
        }
        
        return tmpRotationMatrixList;
    }
    
    /**
     * Get the coordinatates of particle1 and particle2
     * 
     * @param sphereNodeNumber: Sphere node number
     * @param rotNumber: Rotation node number
     * @param xyzData1: Original coordinates of particle1
     * @param xyzData2: Original coordinates of particle2
     * @param isFibonacciSphereAlgorithm: Flag whether fibonaccis phere algorithm is used
     * @return Coordinates of particle1 after rotations
     */
    public static List<double[][][]> getRotationsCoords(
            int sphereNodeNumber,
            int rotNumber, 
            double[][] xyzData1,
            double[][] xyzData2,
            boolean isFibonacciSphereAlgorithm) {
        int confCount1;
        int confCount2;
        int atomSize1;
        int atomSize2;
        double[][][] rotatedXyzData1;
        double[][][] rotatedXyzData2;
        List<double[][]> rotationMatrices1;
        List<double[][]> rotationMatrices2;
        List<double[]> sphereNodeCoord;
        List<double[][][]> result;
        
        atomSize1 = xyzData1.length;
        atomSize2 = xyzData2.length;
        confCount1 = sphereNodeNumber;
        confCount2 = sphereNodeNumber * rotNumber;
        rotatedXyzData1 = new double[confCount1][xyzData1.length][3];
        rotatedXyzData2 = new double[confCount2][xyzData2.length][3];
        
        // Determine rotation matrices used to rotate 
        //   the particle/atom coordinates
        if (isFibonacciSphereAlgorithm) {
            sphereNodeCoord = FibonacciSphere
                    .getSphereNodes(sphereNodeNumber);
        } else {
            //<editor-fold defaultstate="collapsed" desc="Load surface coordinates">
            /* The coordinates for equidistantly distributed points on a sphere 
                from Technical University of Dortmund are used, thanks to 
                J. Fliege and U. Maier
                http://www.mathematik.uni-dortmund.de/lsx/research/projects/fliege/nodes/nodes.html */
            String fileNameSphereNode = 
                  "/de/whs/ibci/mipet/sphereNodes/SphereNodes"
                  + sphereNodeNumber + ".txt";
            Path nodePath = Paths.get(fileNameSphereNode);
            if (!Files.exists(nodePath)) {
                fileNameSphereNode = 
                    "de/whs/ibci/mipet/sphereNodes/SphereNodes"
                    + sphereNodeNumber + ".txt";
            } 
            sphereNodeCoord = RotationUtil
                    .readSphereNodes (fileNameSphereNode);

            //</editor-fold>
        }
        if (atomSize1 == 1) {
            rotatedXyzData1 = new double[][][] {{{0., 0., 0.}}};
        } else {
            final double[] xAxisVector = {1.0, 0.0, 0.0};
            rotationMatrices1 = RotationUtil
                    .getRotationMatrices1 (sphereNodeCoord, xAxisVector);
            rotatedXyzData1 = 
                    new double[confCount1][xyzData1.length][];

            for (int i = 0; i < confCount1; i++) {
                double[][] currRotationMatrix = rotationMatrices1.get(i);

                for (int j = 0; j < xyzData1.length; j++) {
                    rotatedXyzData1[i][j] = MatrixUtil
                            .multiply(currRotationMatrix, xyzData1[j]);
                }

            }

        }
        if (atomSize2 == 1) {
            rotatedXyzData2 = new double[][][] {{{0., 0., 0.}}};
        } else {
            final double[] xAxisVector = {-1.0, 0.0, 0.0};
            rotationMatrices2 = RotationUtil.getRotationMatrices2 (
                    sphereNodeCoord, xAxisVector, rotNumber);

            for (int i = 0; i < confCount2; i++) {
                double[][] currRotationMatrix = rotationMatrices2.get(i);

                for (int j = 0; j < xyzData2.length; j++) {
                    rotatedXyzData2[i][j] = MatrixUtil
                            .multiply(currRotationMatrix, xyzData2[j]);
                }

            }

        }
        result = new LinkedList<>();
        result.add(rotatedXyzData1);
        result.add(rotatedXyzData2);
        return result;
    }
    
    /**
     * Gets the matrices to perform all necessary rotations. Here, the sphere additionally rotates around the input vector.
     *
     * @param sphereNodeCoordinates Coordinates of the sphere nodes.
     * @param vector Vector to map the sphere coordinates on.
     * @param rotationNumber Number of equal rotations around itself.
     * @return List of all rotation matrices.
     */
    public static List<double[][]> getRotationMatrices2(
            List<double[]> sphereNodeCoordinates, 
            double[] vector, 
            int rotationNumber){
        // Caculate capacity
        int nodeCount = sphereNodeCoordinates.size();
        int totalMatrices = nodeCount * rotationNumber;
        List<double[][]> rotationMatrixList = new ArrayList<>(totalMatrices);
        try {
            // Precalculation of angle rotations
            double[][][] angleRotations = new double[rotationNumber][][];
            
            for (int i = 0; i < rotationNumber; i++) {
                double angle = (2.0 * Math.PI * (i + 1)) / rotationNumber;
                angleRotations[i] = RotationUtil
                        .getRotationMatrix(angle, vector);
            }
            
            for (int i = 0; i < nodeCount; i++) {
                double[][] nodeRotation = RotationUtil.getRotationMatrix(
                        sphereNodeCoordinates.get(i), vector);
                
                for (int j = 0; j < rotationNumber; j++) {
                    rotationMatrixList.add(MatrixUtil.multiply(
                            angleRotations[j], nodeRotation));
                }
                
            }
            
            return rotationMatrixList;
        } catch(IllegalArgumentException | NullPointerException anException){
            RotationUtil.LOGGER.log(Level.SEVERE,"Matrices could not be multiplied." ,anException);
            return null;
        }
    }
    //</editor-fold>
}

