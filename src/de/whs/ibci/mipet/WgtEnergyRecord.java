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
 *   energies, Sum of Weights and Sum of weight * energy
 * @author Mirco Daniel
 * @param energys
 *   Weighted energy values (sorted)
 * @param sumWgt
 *   Sum of weights
 * @param sumWgtxE
 *   Sum of weight * energy value
 * @param tinkerXYZ
 *   TinkerXYZ object with local lowest energy
 */
public record WgtEnergyRecord(
        double[] energys,
        double sumWgt,
        double sumWgtxE,
        TinkerXYZ minTinkerXYZ) {
}
