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

import java.util.logging.Logger;

/**
 * This class provides basic matrix algebra functionalities.
 *
 * @author Veit Hucke
 */
public class MatrixUtil {
    //<editor-fold desc="private final class constants" defaultstate="collapsed">
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(MatrixUtil.class.getName());
    //</editor-fold>

    //<editor-fold desc="Public static class constants" defaultstate="collapsed">
    
    /**
     * 3D unit matrix
     */
    public static double[][] UNIT_MATRIX_3D = {{1.0, 0.0, 0.0}, 
        {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="Public static methods">
    /**
     * Calculates the product of a scalar with a matrix.
     *
     * @param aScalar Input scalar.
     * @param aMatrix Input Matrix.
     * @return Product of the scalar multiplication.
     */
    public static double[][] multiplyWithScalar(double aScalar, 
            double[][] aMatrix) {
        double[][] tmpNewMatrix = new double[aMatrix[0]
                .length][aMatrix[1].length];
        
        for(int i = 0; i < aMatrix[0].length; i++) {
            
            for (int j = 0; j < aMatrix[1].length; j++) {
                tmpNewMatrix[i][j] = aMatrix[i][j] * aScalar;
            }
        }
        return tmpNewMatrix;
    }

    /**
     * Calculates the sum of N matrices.
     *
     * @param aMatrixList Variable number of input matrices.
     * @return The sum of the matrices.
     * @throws IllegalArgumentException If the matrices do not have the same amount of rows and columns.
     */
    public static double[][] add(double[][] ... aMatrixList) 
            throws IllegalArgumentException {
        int tmpRows = aMatrixList[0].length;
        int tmpColumns = aMatrixList[0][0].length;
        
        for(double[][] aMatrix : aMatrixList) {
        
            for(int i = 0; i < tmpRows; i++) {
                if (aMatrix.length != tmpRows || aMatrix[i].length != tmpColumns) {
                    throw new IllegalArgumentException("Columns or rows of one or more matrices do not match.");
                }
            }
        }
        double[][] tmpNewMatrix = new double[tmpRows][tmpColumns];
        
        for(int i = 0; i < tmpRows; i++) {
            
            for(int j = 0; j < tmpColumns; j++) {
                
                for(double[][] aMatrix : aMatrixList) {
                    tmpNewMatrix[i][j] += aMatrix[i][j];
                }
            }
        }
        return tmpNewMatrix;
    }

    /**
     * Calculates the product of two matrices.
     *
     * @param aMatrix1 First Matrix.
     * @param aMatrix2 Second Matrix.
     * @return Product of two matrices.
     * @throws IllegalArgumentException If matrices cannot be multiplied 
     *   (rows and columns do not match).
     */
    public static double[][] multiply(double[][] aMatrix1, double[][] aMatrix2) 
            throws IllegalArgumentException {
        int tmpMatrix1Rows = aMatrix1.length;
        int tmpMatrix1Columns = aMatrix1[0].length;
        int tmpMatrix2Rows = aMatrix2.length;
        int tmpMatrix2Columns = aMatrix2[0].length;
        if(tmpMatrix1Columns != tmpMatrix2Rows) {
            throw new IllegalArgumentException("Rows of Matrix1 do not match columns of Matrix2.");
        }
        double[][] tmpNewMatrix = new double[tmpMatrix1Rows][tmpMatrix2Columns];
        for(int i = 0; i < tmpMatrix1Rows; i++) {
            for(int j = 0; j < tmpMatrix2Columns; j++) {
                for(int k = 0; k < tmpMatrix1Columns; k++) {
                    tmpNewMatrix[i][j] += aMatrix1[i][k] * aMatrix2[k][j];
                }
            }
        }
        return tmpNewMatrix;
    }
    
    /**
     * Overloaded method for multiply (dot product) a matrix and a vector
     *   in the order mentioned.
     * @param aMatrix
     *   Rotation matrix
     * @param aVector
     *   A vector
     * @return 
     *   Dot product of rotation matrix and a vector.
     */
    public static double[] multiply(double[][] aMatrix, double[] aVector) {
        double[] tmpResult = new double[3];
        
        for(int i = 0; i < 3; i++) {
            
            for(int j = 0; j < 3; j++) {
                tmpResult[i] += aMatrix[i][j] * aVector[j];
            }
        }
        return tmpResult;
    }
    //</editor-fold>
}
