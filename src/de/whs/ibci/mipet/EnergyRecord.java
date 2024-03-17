/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package de.whs.ibci.mipet;

/**
 * Record for storage energy datas
 * 
 * @author Mirco Daniel
 */

/**
 * Energy record
 * @param distances 
 *   Distances of particle centre i and particle centre j
 * @param energyDatas
 *   sorted energy values of all configurations
 *   [i][j]: i distance index, j configuration index
 * @param minDistance 
 *   equilibrium distance in Angstrom 
 * @param minDistanceIndex 
 *   index of equilibrium distance array
 * @param minEnergy 
 *   minimum of intermolecular energy
 */
public record EnergyRecord(
        double[] distances,
        double[][] energyDatas,
        double minDistance,
        int minDistanceIndex,
        double minEnergy) {
}
