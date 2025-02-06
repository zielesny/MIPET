/**
 * MIPET - Mesoscopic Interaction Parameter Estimation with Tinker
 * Copyright (C) 2025  Achim Zielesny (achim.zielesny@googlemail.com)
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

import java.util.LinkedList;

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
    public static LinkedList<double[]> getSphereNodes(int aNodeNumber) {
        int tmpNodeNumber = aNodeNumber;
        double tmpPhi; // Golden angle in radians
        double tmpTheta; // Golden angle increment
        double tmpX; // x-coordinate [-1; 1]
        double tmpY; // y-coordinate
        double tmpZ; // z-coordinate
        double tmpRadius; // Radius at tmpY
        double[] tmpCoords;
        LinkedList<double[]> tmpCoordsList = new LinkedList<>();
        
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