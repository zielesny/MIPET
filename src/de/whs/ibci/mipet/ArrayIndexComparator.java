/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package de.whs.ibci.mipet;

import java.util.Comparator;

/**
 *
 * @author Mirco Daniel
 */
public class ArrayIndexComparator implements Comparator<Integer>{
    private final Double[] array;

    public ArrayIndexComparator(Double[] anArray)
    {
        this.array = anArray;
    }

    public Integer[] createIndexArray()
    {
        Integer[] tmpIndexes = new Integer[array.length];
        
        for (int i = 0; i < array.length; i++)
        {
            tmpIndexes[i] = i;
        }
        
        return tmpIndexes;
    }

    @Override
    public int compare(Integer aIndex1, Integer aIndex2)
    {
        return array[aIndex1].compareTo(array[aIndex2]);
    }
    
    
}
