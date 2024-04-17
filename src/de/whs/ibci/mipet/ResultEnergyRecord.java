/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Record.java to edit this template
 */
package de.whs.ibci.mipet;

/**
 * Record for storage weighted intermolecular differential pair interaction 
 *   energies
 * @author Mirco Daniel
 * @param particleName1
 *   Name of particle 1
 * @param particleName2
 *   Name of particle 2
 * @param wgtMinEnergy
 *   Lowest weighted differential pair interaction energy
 * @param minEnergy
 *   Lowest differential pair interaction energy
 * @param optMinEnergy
 *   Optimized lowest differential pair interaction energy
 */
public record ResultEnergyRecord(
        String particleName1,    
        String particleName2,
        double wgtMinEnergy,
        double minEnergy,
        double optMinEnergy) {
}
