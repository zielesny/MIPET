/**
 * MIPET - Mesoscopic Interaction Parameter Estimation with Tinker
 * Copyright (C) 2025  Achim Zielesny (achim.zielesny@googlemail.com)
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
 * Record for storage energy datas
 * 
 * @author Mirco Daniel
 */

/**
 * Energy record
 * @param distances 
 *   Distances of particle centre i and particle centre j
 * @param Eminsintermolecular energies of all configurations
 *   [i]: i distance index
 *   Minimum of 
 * @param wgtEmins
 *   Minimum of weighted intermolecular energies of all configurations
 *   [i]: i distance index
 * @param eqDistance 
 *   Equilibrium distance in Angstrom with wgtEmin
 * @param wgtEmin 
 *   Minimum of weighted intermolecular energy
 */
public record EnergyRecord(
        double[] distances,
        double[] Emins,
        double[] wgtEmins,
        double eqDistance,
        double Emin,
        double wgtEmin) {
}
