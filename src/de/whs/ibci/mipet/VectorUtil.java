/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package de.whs.ibci.mipet;

import java.util.logging.Logger;

/**
 * This class provides basic vector algebra functionalities.
 *
 * @author Veit Hucke
 */
public class VectorUtil {
    //<editor-fold desc="Private final class constants" defaultstate="collapsed">
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(VectorUtil.class.getName());
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
        int tmpVector1Dim = aVector1.length;
        int tmpVector2Dim = aVector2.length;
        if (tmpVector1Dim != tmpVector2Dim) {
            throw new IllegalArgumentException("Vectors do not have the same dimension.");
        }
        double tmpProduct = 0;
        for (int i = 0; i < aVector1.length; i++) {
            tmpProduct += aVector1[i] * aVector2[i];
        }
        return tmpProduct;
    }

    /**
     * Method to calculate the Cross product of two vectors.
     *
     * @param aVector1 First vector.
     * @param aVector2 Second vector.
     * @throws IllegalArgumentException If vectors do not have the same dimension.
     * @return A new vector orthogonal to the first and second vector.
     */
    public static double[] crossProduct(double[] aVector1, double[] aVector2) throws IllegalArgumentException {
        //cave: magnitude of n equals the sine of the angle between v1 and v2!
        int tmpVector1Dim = aVector1.length;
        int tmpVector2Dim = aVector2.length;
        if (tmpVector1Dim != tmpVector2Dim) {
            throw new IllegalArgumentException("Vectors do not have the same dimension.");
        }
        double[] tmpNormalVector = new double[3];
        tmpNormalVector[0] = (aVector1[1] * aVector2[2]) - (aVector1[2] * aVector2[1]);
        tmpNormalVector[1] = (aVector1[2] * aVector2[0]) - (aVector1[0] * aVector2[2]);
        tmpNormalVector[2] = (aVector1[0] * aVector2[1]) - (aVector1[1] * aVector2[0]);
        return tmpNormalVector;
    }

    /**
     * Calculates the norm (or magnitude) of a vector.
     *
     * @param aVector Input vector.
     * @return The Norm of the input vector.
     */
    public static double normOfVector(double[] aVector) {
        double tmpNormSquared = 0;
        for(double element : aVector) {
            tmpNormSquared += Math.pow(element, 2);
        }
        return Math.sqrt(tmpNormSquared);
    }

    /**
     * Normalizes a Vector.
     *
     * @param aVector Vector to normalize.
     * @return Normalized vector.
     */
    public static double[] normalizeVector(double[] aVector) {
        double tmpNormOfVector = normOfVector(aVector);
        return multiplyWithScalar((1/tmpNormOfVector), aVector);
    }

    /**
     * Calculates the product of a scalar with a vector.
     *
     * @param aScalar Input scalar.
     * @param aVector Input vector.
     * @return Product of the scalar multiplication.
     */
    public static double[] multiplyWithScalar(double aScalar, double[] aVector) {
        double[] tmpNewVector = new double[aVector.length];
        for(int i = 0; i < aVector.length; i++) {
            tmpNewVector[i] = aVector[i] * aScalar;
        }
        return tmpNewVector;
    }

    /**
     * Calculates the sum of N vectors.
     *
     * @param aVectorList Variable number of input vectors
     * @throws IllegalArgumentException If the input vector do not have the same dimension
     * @return The sum of the Vectors
     */
    public static double[] add(double[] ...aVectorList) throws IllegalArgumentException {
        int tmpDimension = aVectorList[0].length;
        for(double[] aVector : aVectorList) {
            if(aVector.length != tmpDimension) {
                throw new IllegalArgumentException("One or more vectors do not have the same dimension.");
            }
        }
        double[] tmpNewVector = new double[tmpDimension];
        for(int i = 0; i < tmpDimension; i++) {
            for(double[] aVector : aVectorList) {
                tmpNewVector[i] += aVector[i];
            }
        }
        return tmpNewVector;
    }

    /**
     * Calculates the new vector which is moved in x-direction by aXLength
     * 
     * @param aVectorMoveTo A vector to be moved
     * @param aXLength Length in x-direction
     * @return New vectors which are moved in x-direction by aXLength
     */
    public static double[][] moveX(double[][] aVectorMoveTo, double aXLength) {
        int tmpDim = aVectorMoveTo.length;
        double[][] tmpNewVectorList = new double[tmpDim][3];
        
        for (int i = 0; i < tmpDim; i++) {
            
            for (int j = 0; j < 3; j++) {
                tmpNewVectorList[i][0] = aVectorMoveTo[i][0] + aXLength;
                tmpNewVectorList[i][1] = aVectorMoveTo[i][1];
                tmpNewVectorList[i][2] = aVectorMoveTo[i][2];
            }
        }
        return tmpNewVectorList;
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
        int tmpDim1 = aVectorMoveTo.length;
        int tmpDim2 = aVectorMoveTo[0].length;  
        double[][][] tmpNewVecotList = new double[tmpDim1][tmpDim2][3];
        
        for (int i = 0; i < tmpDim1; i++) {
            tmpNewVecotList[i] = moveX(aVectorMoveTo[i], aXLength);
        }
        return tmpNewVecotList;
    }

    /**
     * Calculates the dyadic product (outer product) of a vector v and its transposed version v^T.
     *
     * @param aVector Input vector.
     * @return A nxn Matrix.
     */
    public static double[][] dyadicProduct(double[] aVector) {
        double[][] tmpDyadic = new double[3][3];
        for(int i = 0; i < aVector.length; i++) {
            for(int j = 0; j < aVector.length; j++) {
                tmpDyadic[i][j] = aVector[i] * aVector[j];
            }
        }
        return tmpDyadic;
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