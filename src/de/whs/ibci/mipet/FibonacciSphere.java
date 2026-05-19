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

import java.util.ArrayList;
import java.util.List;

/**
 * Class FibonacciSphere
 * @author Daniel
 * Saff, E.B. and Kuijlaars, A.B.J. (1997) Distributing Many Points on a Sphere. The Mathematical Intelligencer, 19, 5-11.
 * <a href="https://doi.org/10.1007/BF03024331">...</a>
 *
 */
public class FibonacciSphere {
    /**
     * Golden angle
     */
    private static final double GOLDEN_ANGLE = (3.0 - Math.sqrt(5.0)) * Math.PI;
    /**
     * Constructor FibonacciSphere
     */
    public FibonacciSphere() {
    }
    
    /**
     * Determine the coordinates of nodes of fibonaccisphere 
     * @param aNodeNumber Number of nodes
     * @return xyz-Coordinates of nodes
     */
    public static List<double[]> getSphereNodes(int aNodeNumber) {
        if (aNodeNumber < 1) {
            throw new IllegalArgumentException("Node number must be at least 1.");
        }
        List<double[]> coordsList = new ArrayList<>(aNodeNumber);
        if (aNodeNumber == 1) {
            coordsList.add(new double[]{0.0, 1.0, 0.0});
            return coordsList;
        }
         
        for (int i = 0; i < aNodeNumber; i++) {
            double y = 1.0 - (i / (double)(aNodeNumber - 1)) * 2.0;
            double radius = Math.sqrt(1.0 - y * y);
            double theta = i * GOLDEN_ANGLE;
            double x = Math.cos(theta) * radius;
            double z = Math.sin(theta) * radius;
            coordsList.add(new double[]{x, y, z});
        }
        
        return coordsList;
    }
}