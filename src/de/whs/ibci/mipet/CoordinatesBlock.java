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

/**
 *
 * @author Mirco Daniel
 */
public class CoordinatesBlock {
    private final double[][] coord1;
    private final double[][][] coord2;

    public CoordinatesBlock(double[][] coord1, double[][][] coord2) {
        this.coord1 = coord1;
        this.coord2 = coord2;
    }

    public double[][] getCoord1() { 
        return coord1; 
    }
    
    public double[][][] getCoord2() { 
        return coord2; 
    }
}