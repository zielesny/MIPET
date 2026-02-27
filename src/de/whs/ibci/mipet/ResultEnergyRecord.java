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
 * Record for storage weighted intermolecular differential pair interaction 
 *   energies
 * @author Mirco Daniel
 * @param particleName1:
 *   Name of particle 1
 * @param particleName2:
 *   Name of particle 2
 * @param wgt_Emin:
 *   Weighted differential pair interaction energy (Lowest = Lowest energy configuration at same stage)
 * @param wgt0_Opt_Emin:
 *   Differential pair interaction energy from optimized configuration
 * @param wgt0_Rgd_Emin:
 *   Differential pair interaction energy from rigid-optimized configuration
 * @param wgt0_Opt0_Emin:
 *   Differential pair interaction energy from lowest energy from all calculated configurations
 */
public record ResultEnergyRecord(
        String particleName1,    
        String particleName2,
        double wgt_Emin,
        double wgt0_Opt_Emin,
        double wgt0_Rgd_Emin,
        double wgt0_Opt0_Emin) {
}