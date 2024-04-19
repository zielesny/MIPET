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
 * @param wgt_Opt_MinEnergy
 *   Weighted differential pair interaction energy (Lowest = energy from optimized configuration)
 * @param wgt_Rgd_MinEnergy
 *   Weighted differential pair interaction energy (Lowest = energy from rigid-optimized configuration)
 * @param wgt_Opt0_MinEnergy
 *   Weighted differential pair interaction energy (Lowest = lowest energy from all calculated configurations)
 * @param wgt0_Opt_MinEnergy
 *   Differential pair interaction energy from optimized configuration
 * @param wgt0_Rgd_MinEnergy
 *   Differential pair interaction energy from rigid-optimized configuration
 * @param wgt0_Opt0_MinEnergy
 *   Differential pair interaction energy from lowest energy from all calculated configurations
 */
public record ResultEnergyRecord(
        String particleName1,    
        String particleName2,
        double wgt_Opt_MinEnergy,
        double wgt_Rgd_MinEnergy,
        double wgt_Opt0_MinEnergy,
        double wgt0_Opt_MinEnergy,
        double wgt0_Rgd_MinEnergy,
        double wgt0_Opt0_MinEnergy) {
}
