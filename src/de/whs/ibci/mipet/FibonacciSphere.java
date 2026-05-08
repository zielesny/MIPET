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

import java.util.LinkedList;

/**
 * Class FibonacciSphere
 * @author Daniel
 * Saff, E.B. and Kuijlaars, A.B.J. (1997) Distributing Many Points on a Sphere. The Mathematical Intelligencer, 19, 5-11.
 * <a href="https://doi.org/10.1007/BF03024331">...</a>
 *
 */
public class FibonacciSphere {
    
    /**
     * Constructor FibonacciSphere
     */
    public FibonacciSphere() {
    }
    
    /**
     * Determine the coordinates of nodes of fibonaccisphere 
     * @param aNodeNumber 
     *   Number of nodes
     * @return 
     *   xyz-Coordinates of nodes
     */
    public static LinkedList<double[]> getSphereNodes(int aNodeNumber) {
        double phi; // Golden angle in radians
        double theta; // Golden angle increment
        double x; // x-coordinate [-1; 1]
        double y; // y-coordinate
        double z; // z-coordinate
        double radius; // Radius at tmpY
        double[] coords;
        LinkedList<double[]> coordsList = new LinkedList<>();
        
        phi = (3.0 - Math.sqrt(5.0)) * Math.PI;
         
        for (int i = 0; i < aNodeNumber; i++) {
            y = 1 - (i / (double)(aNodeNumber - 1)) * 2;
            radius = Math.sqrt(1 - y * y);
            theta = i * phi;
            x = Math.cos(theta) * radius;
            z = Math.sin(theta) * radius;
            coords = new double[3];
            coords[0]= x;
            coords[1]= y;
            coords[2]= z;
            coordsList.add(coords);
        }
        
        return coordsList;
    }
}