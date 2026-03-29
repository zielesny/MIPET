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
 * Class Flat3DArray
 * 
 * @author Mirco Daniel
 */
public class Flat3DArray {
    public final int dim1, dim2, dim3;
    private final double[] data;
    private final int stride1, stride2;
   
    /**
     * Constructor of the class
     * 
     * @param d1 dimension of the configuration
     * @param d2 dimension of the atom
     * @param d3 dimension of the coordinate
     */
    public Flat3DArray(int d1, int d2, int d3) {
        this.dim1 = d1;
        this.dim2 = d2;
        this.dim3 = d3;
        this.stride2 = d3;
        this.stride1 = d2 * d3;
        this.data = new double[d1 * stride1];
    }

    /**
     * Sets a value of flattend array
     * 
     * @param i index of the configuration
     * @param j index of the atom
     * @param k index of the coordinate (x=0, y=1, z=2)
     * @param value The new value
     */
    public void set(int i, int j, int k, double value) {
        data[i * stride1 + j * stride2 + k] = value;
    }

    /**
     * Return the flattened array
     * @param i index of the configuration
     * @param j index of the atom
     * @param k index of the coordinate (x=0, y=1, z=2)
     * @return Flat data array
     */
    public double get(int i, int j, int k) {
        return data[i * stride1 + j * stride2 + k];
    }

    /**
     * Converts a normal double[][][] into a flat structure
     * 
     * @param source A source 3d-array
     * @return Flattened 1d-array
     */
    public static Flat3DArray createFrom(double[][][] source) {
        int d1 = source.length;
        int d2 = source[0].length;
        int d3 = source[0][0].length;
        
        Flat3DArray flat = new Flat3DArray(d1, d2, d3);
        
        for (int i = 0; i < d1; i++) {
            
            for (int j = 0; j < d2; j++) {
                System.arraycopy(source[i][j], 0, flat.data,
                        i * flat.stride1 + j * flat.stride2, d3);
            }
            
        }
        
        return flat;
    }

    /**
     * Getter for flatterned 1d-array
     * 
     * @return Flatterned 1d-array
     */
    public double[] getRawData() {
        return data;
    }
}