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

import java.util.logging.Logger;

/**
 * This class provides basic matrix algebra functionalities.
 *
 * @author Veit Hucke, Mirco Daniel
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

    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /** 
     * Constructor MatrixUtil
     */
    public MatrixUtil() {
    }
    
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
     * @param matrixA First Matrix.
     * @param matrixB Second Matrix.
     * @return Product of two matrices.
     */
    public static double[][] multiply(double[][] matrixA, double[][] matrixB) {
        int rowsA = matrixA.length;
        int colsA = matrixA[0].length;
        int rowsB = matrixB.length;
        int colsB = matrixB[0].length;
        
        if(colsA != rowsB) {
            throw new IllegalArgumentException("Columns of Matrix A do not match rows of Matrix B.");
        }
        
        double[][] resultMatrix = new double[rowsA][colsB];
        
        // Optimized loop order: i->k->j
        for(int i = 0; i < rowsA; i++) {
            double[] rowA = matrixA[i]; // Cache current row from A
            double[] rowResult = resultMatrix[i];
            
            for(int k = 0; k < colsB; k++) {
                double aVal = rowA[k];
                double[] rowB = matrixB[k]; // Cache current row from B
                
                // Innerst loop accesses memory sequentially
                for(int j = 0; j < colsB; j++) {
                    rowResult[j] += aVal * rowB[j];
                }
                
            }
            
        }
        return resultMatrix;
    }
    
    /**
     * Overloaded method for multiply (dot product) a matrix and a vector
     *   in the order mentioned.
     * @param matrix
     *   Rotation matrix
     * @param vector
     *   A vector
     * @return 
     *   Dot product of rotation matrix and a vector.
     */
    public static double[] multiply(double[][] matrix, double[] vector) {
        int rows = matrix.length;
        int cols = matrix[0].length;
        
        // Security check
        if (cols != vector.length) {
            throw new IllegalArgumentException("Matrix columns (" + cols 
                    + ") must match vector length (" + vector.length + ").");
        }
        
        double[] result = new double[rows];
        
        for(int i = 0; i < rows; i++) {
            double[] row = matrix[i]; // Cache current row
            double sum = 0.0; // Local accumulator
            
            for(int j = 0; j < cols; j++) {
                sum += row[j] * vector[j];
            }
            
            result[i] = sum;
        }
        
        return result;
    }
    
    //</editor-fold>
}
