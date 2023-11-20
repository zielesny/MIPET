/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Record.java to edit this template
 */
package de.whs.ibci.mipet;

/**
 * Record for storage weighted intermolecular differential pair interaction 
 *   energies
 * 
 * @author Mirco Daniel
 */
public record ResultEnergyRecord(
        String particleName1,    
        String particleName2,
        double energyValue) {
}
