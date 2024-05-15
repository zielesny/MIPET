/**
 * MIPET - Mesoscopic Interaction Parameter Estimation with Tinker
 * Copyright (C) 2024  Achim Zielesny (achim.zielesny@googlemail.com)
 * 
 * Source code is available at <https://github.com/zielesny/MIPET>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
