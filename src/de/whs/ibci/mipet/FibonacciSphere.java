/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package de.whs.ibci.mipet;

import java.util.LinkedList;
import java.util.List;

/**
 * Class FibonacciSphere
 * @author Daniel
 * Saff, E.B. and Kuijlaars, A.B.J. (1997) Distributing Many Points on a Sphere. The Mathematical Intelligencer, 19, 5-11.
 * https://doi.org/10.1007/BF03024331
 * 
 */
public class FibonacciSphere {
    
    /**
     * Determine the coordinates of nodes of fibonaccisphere 
     * @param aNodeNumber 
     *   Number of nodes
     * @return 
     *   xyz-Coordinates of nodes
     */
    public static List<double[]> getSphereNodes(int aNodeNumber) {
        int tmpNodeNumber = aNodeNumber;
        double tmpPhi; // Golden angle in radians
        double tmpTheta; // Golden angle increment
        double tmpX; // x-coordinate [-1; 1]
        double tmpY; // y-coordinate
        double tmpZ; // z-coordinate
        double tmpRadius; // Radius at tmpY
        List<double[]> tmpCoordsList = new LinkedList<>();
        double[] tmpCoords;
        
        tmpPhi = (3.0 - Math.sqrt(5.0)) * Math.PI;
         
        for (int i = 0; i < tmpNodeNumber; i++) {
            tmpY = 1 - (i / (double)(tmpNodeNumber - 1)) * 2;
            tmpRadius = Math.sqrt(1 - tmpY * tmpY);
            tmpTheta = i * tmpPhi;
            tmpX = Math.cos(tmpTheta) * tmpRadius;
            tmpZ = Math.sin(tmpTheta) * tmpRadius;
            tmpCoords = new double[3];
            tmpCoords[0]= tmpX;
            tmpCoords[1]= tmpY;
            tmpCoords[2]= tmpZ;
            tmpCoordsList.add(tmpCoords);
        }
        
        return tmpCoordsList;
    }
    
}
