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

import java.util.Comparator;

/**
 * Class ArrayIndexComparator
 * @author Mirco Daniel
 */
public class ArrayIndexComparator implements Comparator<Integer>{
    private final Double[] array;

    /**
     * Constructor ArrayIndexComparator
     * @param anArray: A double array
     */
    public ArrayIndexComparator(Double[] anArray){
        this.array = anArray;
    }

    /**
     * Creates an index array object
     * @return Index array as an object
     */
    public Integer[] createIndexArray(){
        Integer[] indices = new Integer[array.length];
        
        for (int i = 0; i < array.length; i++)
        {
            indices[i] = i;
        }
        
        return indices;
    }

    @Override
    public int compare(Integer aIndex1, Integer aIndex2){
        return this.array[aIndex1].compareTo(this.array[aIndex2]);
    }
}