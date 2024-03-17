/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Record.java to edit this template
 */
package de.whs.ibci.mipet;

/**
 * Record for storage weighted intermolecular differential pair interaction 
 *   energies
 * @param particleName1
 *   Name of particle 1
 * @param particleName2
 *   Name of particle 2
 * @param minWgtEnergy
 *   Lowest weighted differential pair interaction energy
 * @param minGlbEnergy
 *   Lowest differential pair interaction energy
 * @author Mirco Daniel
 */
public record ResultEnergyRecord(
        String particleName1,    
        String particleName2,
        double minWgtEnergy,
        double minGlbEnergy) {
}
