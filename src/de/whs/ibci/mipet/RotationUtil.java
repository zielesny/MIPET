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
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
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
     * @param inputStream Inputstream
     * @return A List of all sphere node coordinates, saved as double arrays
     */
    public static List<double[]> readSphereNodes(InputStream inputStream) {
        if (inputStream == null) {
            throw new IllegalArgumentException("InputStream cannot be null (File not found in classpath).");
        }
        List<double[]> coordList = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(
                new InputStreamReader(inputStream, StandardCharsets.UTF_8))) {
            String line;
            
            while ((line = reader.readLine()) != null) {
                int startIdx = line.indexOf('{');
                int endIdx = line.indexOf('}');
                if (startIdx == -1 || endIdx == -1 || startIdx >= endIdx) {
                    continue; 
                }
                String newSubString = line.substring(startIdx + 1, endIdx);
                String[] tokens = newSubString.split(",\\s+");
                if (tokens.length != 3) {
                    continue;
                }
                double[] parsedTokens = new double[3];
                boolean parseError = false;

                for (int i = 0; i < tokens.length; i++) {
                    try {
                        parsedTokens[i] = Double.parseDouble(tokens[i].replace("*10^", "E"));
                    } catch (NumberFormatException ex) {
                        LOGGER.log(Level.SEVERE, "String could not be parsed: " + tokens[i], ex);
                        parseError = true;
                        break;
                    }
                }

                if (!parseError) {
                    coordList.add(parsedTokens);
                }
            }
            return coordList;

        } catch (Exception ex) {
            LOGGER.log(Level.SEVERE, "Error reading the input stream.", ex);
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
        List<double[][]> rotMatrixList = new LinkedList<>();
        
        for(double[] anElement : sphereNodeCoordinates) {
            rotMatrixList.add(RotationUtil
                    .getRotationMatrix(anElement, aVector));
        }
        
        return rotMatrixList;
    }
    
    /**
     * Get the coordinatates of molecule1 and molecule2
     * 
     * @param sphereNodeNumber: Sphere node number
     * @param rotNumber: Rotation node number
     * @param xyzData1: Original coordinates of particle1
     * @param xyzData2: Original coordinates of particle2
     * @param isFibonacciSphereAlgorithm: Flag whether fibonaccis phere algorithm is used
     * @return Coordinates of particle1 after rotations
     */
    public static RotationRecord getRotationsCoords(
            int sphereNodeNumber,
            int rotNumber, 
            double[][] xyzData1,
            double[][] xyzData2,
            boolean isFibonacciSphereAlgorithm) {
        int confSize1 = sphereNodeNumber;
        int confSize2 = sphereNodeNumber * rotNumber;
        
        // Load sphere nodes
        List<double[]> sphereNodeCoord;
        if (isFibonacciSphereAlgorithm) {
            sphereNodeCoord = FibonacciSphere.getSphereNodes(sphereNodeNumber);
        } else {
            //<editor-fold defaultstate="collapsed" desc="Load surface coordinates">
            /* The coordinates for equidistantly distributed points on a sphere 
                from Technical University of Dortmund are used, thanks to 
                J. Fliege and U. Maier
                http://www.mathematik.uni-dortmund.de/lsx/research/projects/fliege/nodes/nodes.html */
            String resourcePath = "/de/whs/ibci/mipet/sphereNodes/SphereNodes"
                    + sphereNodeNumber + ".txt";
            try (InputStream is = RotationUtil.class.getResourceAsStream(
                    resourcePath)) {
                if (is == null) {
                    throw new RuntimeException("Could'nt find the file in the JAR: " 
                            + resourcePath);
                }
                sphereNodeCoord = RotationUtil.readSphereNodes(is);
            } catch (Exception e) {
                throw new RuntimeException("Error during load sphereNodes.", e);
            }

            //</editor-fold>
        }
        
        // Data of molecule 1
        double[][][] rotatedXyzData1;
        if (xyzData1.length == 1) {
            rotatedXyzData1 = new double[confSize1][1][3];
        } else {
            final double[] xAxisVector = {1.0, 0.0, 0.0};
            List<double[][]> rotationMatrices1 = RotationUtil
                    .getRotationMatrices1(sphereNodeCoord, xAxisVector);
            rotatedXyzData1 = new double[confSize1][xyzData1.length][3];

            for (int i = 0; i < confSize1; i++) {
                double[][] currRotationMatrix = rotationMatrices1.get(i);

                for (int j = 0; j < xyzData1.length; j++) {
                    rotatedXyzData1[i][j] = MatrixUtil
                            .multiply(currRotationMatrix, xyzData1[j]);
                }

            }

        }
        
        // Data of molecule 2
        double[][][] rotatedXyzData2;
        if (xyzData2.length == 1) {
            rotatedXyzData2 = new double[confSize2][1][3];
        } else {
            final double[] xAxisVector = {-1.0, 0.0, 0.0};
            List<double[][]> rotationMatrices2 = RotationUtil
                    .getRotationMatrices2(sphereNodeCoord, xAxisVector, 
                            rotNumber);
            rotatedXyzData2 = new double[confSize2][xyzData2.length][3];

            for (int i = 0; i < confSize2; i++) {
                double[][] currRotationMatrix = rotationMatrices2.get(i);

                for (int j = 0; j < xyzData2.length; j++) {
                    rotatedXyzData2[i][j] = MatrixUtil
                            .multiply(currRotationMatrix, xyzData2[j]);
                }

            }

        }
        return new RotationRecord(rotatedXyzData1, rotatedXyzData2);
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
        
        // Checks 
        if (sphereNodeCoordinates == null || vector == null) {
            LOGGER.log(Level.SEVERE, 
                    "Input parameters (coordinates or vector) must not be null.");
            return Collections.emptyList();
        }
        if (rotationNumber <= 0) {
            LOGGER.log(Level.SEVERE, "rotationNumber must be strictly positive.");
            return Collections.emptyList();
        }
        
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
        } catch(IllegalArgumentException e){
            LOGGER.log(Level.SEVERE, "Math/Matrix operation failed: " 
                    + e.getMessage(), e);
            return Collections.emptyList();
        }
    }
    
    //</editor-fold>
}