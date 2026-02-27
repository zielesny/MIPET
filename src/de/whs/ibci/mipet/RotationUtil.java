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
import java.util.LinkedList;
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
    public static double[][] getRotationMatrix(double anAngle, 
            double[] aRotationAxis) {
        try {
            double[] tmpUnitAxis = VectorUtil.normalizeVector(aRotationAxis);
            double[][] tmpDyadic = VectorUtil.dyadicProduct(tmpUnitAxis);
            double[][] tmpSkewSymmetricMatrix = VectorUtil
                    .generateSkewSymmetricMatrix(tmpUnitAxis);
            return MatrixUtil.add(
                    (MatrixUtil.multiplyWithScalar(Math.cos(anAngle), 
                            MatrixUtil.UNIT_MATRIX_3D)),
                    (MatrixUtil.multiplyWithScalar((1 - Math.cos(anAngle)), 
                            tmpDyadic)),
                    (MatrixUtil.multiplyWithScalar(Math.sin(anAngle), 
                            tmpSkewSymmetricMatrix)));
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
            double tmpAngle = Math.acos(VectorUtil.dotProduct(tmpUnitVector1, 
                    tmpUnitVector2)); //Angle between v1 and v2, only for unit vectors!
            double[] tmpRotationAxis = VectorUtil
                    .crossProduct(aVector1, aVector2); //Axis to rotate around (vector n is perpendicular to v1 and v2)
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
    public static LinkedList<double[][]> getRotationMatrices1(
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
     * @param aSphereNodeNumber: Sphere node number
     * @param aRotNumber: Rotation node number
     * @param aXyzData1: Original coordinates of particle1
     * @param aXyzData2: Original coordinates of particle2
     * @param aIsFibonacciSphereAlgorithm: Flag whether fibonaccis phere algorithm is used
     * @return Coordinates of particle1 after rotations
     */
    public static LinkedList<double[][][]> getRotationsCoords(
            int aSphereNodeNumber,
            int aRotNumber, 
            double[][] aXyzData1,
            double[][] aXyzData2,
            boolean aIsFibonacciSphereAlgorithm) {
        int tmpConfNumber1;
        int tmpConfNumber2;
        int tmpAtomSize1;
        int tmpAtomSize2;
        double[][][] tmpXyzRotData1;
        double[][][] tmpXyzRotData2;
        LinkedList<double[][]> tmpRotMatrices1;
        LinkedList<double[][]> tmpRotMatrices2;
        LinkedList<double[]> tmpSphereNodeCoord;
        LinkedList<double[][][]> tmpResult;
        
        tmpAtomSize1 = aXyzData1.length;
        tmpAtomSize2 = aXyzData2.length;
        tmpConfNumber1 = aSphereNodeNumber;
        tmpConfNumber2 = aSphereNodeNumber * aRotNumber;
        tmpXyzRotData1 = new double[tmpConfNumber1][aXyzData1.length][3];
        tmpXyzRotData2 = new double[tmpConfNumber2][aXyzData2.length][3];
        
        // Determine rotation matrices used to rotate 
        //   the particle/atom coordinates
        if (aIsFibonacciSphereAlgorithm) {
            tmpSphereNodeCoord = FibonacciSphere
                    .getSphereNodes(aSphereNodeNumber);
        } else {
            //<editor-fold defaultstate="collapsed" desc="Load surface coordinates">
            /* The coordinates for equidistantly distributed points on a sphere 
                from Technical University of Dortmund are used, thanks to 
                J. Fliege and U. Maier
                http://www.mathematik.uni-dortmund.de/lsx/research/projects/fliege/nodes/nodes.html */


            // Development version
//            String tmpFileNameSphereNode = 
//                    "de/whs/ibci/mipet/sphereNodes/SphereNodes"
//                    + aSphereNodeNumber + ".txt";

            // Distribution version
          String tmpFileNameSphereNode = 
                  "/de/whs/ibci/mipet/sphereNodes/SphereNodes"
                  + aSphereNodeNumber + ".txt";
            Path tmpNodeFile = Paths.get(tmpFileNameSphereNode);
            if (!Files.exists(tmpNodeFile)) {
                LOGGER.log(Level.SEVERE, "No NodeFile found.");
            } 
            tmpSphereNodeCoord = RotationUtil
                    .readSphereNodes (tmpFileNameSphereNode);

            //</editor-fold>
        }
        if (tmpAtomSize1 == 1) {
            tmpXyzRotData1 = new double[1][1][];
            tmpXyzRotData1[0][0] = new double[] {0., 0., 0.};
        } else {
            tmpRotMatrices1 = RotationUtil
                    .getRotationMatrices1 (tmpSphereNodeCoord, 
                            new double[] {1., 0., 0.});

            for (int i = 0; i < tmpConfNumber1; i++) {

                for (int j = 0; j < aXyzData1.length; j++) {
                    tmpXyzRotData1[i][j] = MatrixUtil
                            .multiply(tmpRotMatrices1.get(i), aXyzData1[j]);

                }

            }

        }
        if (tmpAtomSize2 == 1) {
            tmpXyzRotData2 = new double[1][1][];
            tmpXyzRotData2[0][0] = new double[] {0., 0., 0.};
        } else {
            tmpRotMatrices2 = RotationUtil
                    .getRotationMatrices2 (tmpSphereNodeCoord, 
                            new double[] {-1., 0., 0.}, aRotNumber);

            for (int i = 0; i < tmpConfNumber2; i++) {

                for (int j = 0; j < aXyzData2.length; j++) {
                    tmpXyzRotData2[i][j] = MatrixUtil
                            .multiply(tmpRotMatrices2.get(i), aXyzData2[j]);
                }

            }

        }
        tmpResult = new LinkedList<>();
        tmpResult.add(tmpXyzRotData1);
        tmpResult.add(tmpXyzRotData2);
        return tmpResult;
    }
    
    /**
     * Gets the matrices to perform all necessary rotations. Here, the sphere additionally rotates around the input vector.
     *
     * @param sphereNodeCoordinates Coordinates of the sphere nodes.
     * @param aVector Vector to map the sphere coordinates on.
     * @param aRotationNumber Number of equal rotations around itself.
     * @return List of all rotation matrices.
     */
    public static LinkedList<double[][]> getRotationMatrices2(
            List<double[]> sphereNodeCoordinates, 
            double[] aVector, 
            int aRotationNumber){
        LinkedList<double[][]> tmpRotationMatrixList = new LinkedList<>();
        int tmpSphereNodeNumber;
        
        tmpSphereNodeNumber = sphereNodeCoordinates.size();
        try {
            
            for (int i = 0; i < tmpSphereNodeNumber; i++) {
                for (int j = 0; j < aRotationNumber; j++) {
                    tmpRotationMatrixList.add(MatrixUtil.multiply(
                            RotationUtil.getRotationMatrix(((2 * Math.PI * (j + 1)) / aRotationNumber), aVector),
                            RotationUtil.getRotationMatrix(sphereNodeCoordinates.get(i), aVector)));
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

