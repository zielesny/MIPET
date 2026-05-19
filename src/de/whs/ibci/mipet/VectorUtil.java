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
 * This class provides basic vector algebra functionalities.
 *
 * @author Veit Hucke, Mirco Daniel
 */
public class VectorUtil {
    //<editor-fold desc="Private final class constants" defaultstate="collapsed">
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(VectorUtil.class.getName());
    //</editor-fold>

    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor
     */
    public VectorUtil() {
    }
    
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Public static methods">
    /**
     * Method to calculate the Dot product of two vectors.
     *
     * @param aVector1 First vector.
     * @param aVector2 Second vector.
     * @throws IllegalArgumentException If vectors do not have the same dimension.
     * @return The scalar product of the first and second vector.
     */
    public static double dotProduct(double[] aVector1, double[] aVector2) throws IllegalArgumentException {
        int vector1Dim = aVector1.length;
        int vector2Dim = aVector2.length;
        if (vector1Dim != vector2Dim) {
            throw new IllegalArgumentException("Vectors do not have the same dimension.");
        }
        double product = 0;
        
        for (int i = 0; i < aVector1.length; i++) {
            product += aVector1[i] * aVector2[i];
        }
        
        return product;
    }

    /**
     * Method to calculate the Cross product of two vectors.
     *
     * @param v1 First 3D vector.
     * @param v2 Second 3D vector.
     * @return A new array representing the normal vector (cross product).
     */
    public static double[] crossProduct(double[] v1, double[] v2) {
        if (v1.length != 3 || v2.length != 3) {
            throw new IllegalArgumentException("Both vectors must have exactly 3 dimensions.");
        }
        return new double[] {
            (v1[1] * v2[2]) - (v1[2] * v2[1]),
            (v1[2] * v2[0]) - (v1[0] * v2[2]),
            (v1[0] * v2[1]) - (v1[1] * v2[0])
        };
    }

    /**
     * Calculates the norm (or magnitude) of a vector.
     *
     * @param aVector Input vector.
     * @return The Norm of the input vector.
     */
    public static double normOfVector(double[] aVector) {
        double normQ = 0;
        for(double element : aVector) {
            normQ += element * element;
        }
        return Math.sqrt(normQ);
    }

    /**
     * Normalizes a Vector.
     *
     * @param aVector Vector to normalize.
     * @return Normalized vector.
     */
    public static double[] normalizeVector(double[] aVector) {
        double normOfVector = normOfVector(aVector);
        return multiplyWithScalar((1 / normOfVector), aVector);
    }

    /**
     * Calculates the product of a scalar with a vector.
     *
     * @param aScalar Input scalar.
     * @param aVector Input vector.
     * @return Product of the scalar multiplication.
     */
    public static double[] multiplyWithScalar(double aScalar, double[] aVector) {
        double[] newVector = new double[aVector.length];
        
        for(int i = 0; i < aVector.length; i++) {
            newVector[i] = aVector[i] * aScalar;
        }
        
        return newVector;
    }

    /**
     * Calculates the sum of N vectors.
     *
     * @param aVectorList Variable number of input vectors
     * @throws IllegalArgumentException If the input vector do not have the same dimension
     * @return The sum of the Vectors
     */
    public static double[] add(double[] ...aVectorList) throws IllegalArgumentException {
        int dimension = aVectorList[0].length;
        
        for(double[] aVector : aVectorList) {
            if(aVector.length != dimension) {
                throw new IllegalArgumentException("One or more vectors do not have the same dimension.");
            }
        }
        
        double[] newVector = new double[dimension];
        
        for(int i = 0; i < dimension; i++) {
            
            for(double[] aVector : aVectorList) {
                newVector[i] += aVector[i];
            }
            
        }
        
        return newVector;
    }

    /**
     * Calculates the new vector which is moved in x-direction by aXLength
     * 
     * @param vectors A 2D array containing the vectors to be moved (format: n x 3).
     * @param xLength The distance by which the vectors should be moved along the X-axis.
     * @return A new 2D array containing the moved vectors.
     */
    public static double[][] moveX(double[][] vectors, double xLength) {
        double[][] newVectors = new double[vectors.length][3];
        
        for (int i = 0; i < vectors.length; i++) {
            newVectors[i][0] = vectors[i][0] + xLength;
            newVectors[i][1] = vectors[i][1];
            newVectors[i][2] = vectors[i][2];
        }
        
        return newVectors;
    } 

    /**
     * Overloaded method moveX
     * 
     * @param aVectorMoveTo i: particleId, j: atomId k: xyz
     * @param aXLength Length in x-direction
     * @return New vectors which are moved in x-direction by aXLength
     */
    public static double[][][] moveX(double[][][] aVectorMoveTo, 
            double aXLength) {
        int dim1 = aVectorMoveTo.length;
        int dim2 = aVectorMoveTo[0].length;  
        double[][][] newVecotList = new double[dim1][dim2][3];
        
        for (int i = 0; i < dim1; i++) {
            newVecotList[i] = moveX(aVectorMoveTo[i], aXLength);
        }
        
        return newVecotList;
    }

    /**
     * Calculates the dyadic product (outer product) of a vector v and its transposed version v^T.
     *
     * @param aVector Input vector.
     * @return A nxn Matrix.
     */
    public static double[][] dyadicProduct(double[] aVector) {
        double[][] dyadic = new double[3][3];
        
        for(int i = 0; i < aVector.length; i++) {
            
            for(int j = 0; j < aVector.length; j++) {
                dyadic[i][j] = aVector[i] * aVector[j];
            }
            
        }
        
        return dyadic;
    }

    /**
     * Calculates the cross product matrix of the input vector.
     * This matrix has the same effect as calculating the cross product with this input vector.
     *
     * @param aVector Input vector.
     * @return A nxn matrix.
     */
    public static double[][] generateSkewSymmetricMatrix(double[] aVector) {
        return new double[][] {
                {0, -aVector[2], aVector[1]},
                {aVector[2], 0, -aVector[0]},
                {-aVector[1], aVector[0], 0}
        };
    }
    
    //</editor-fold>
}