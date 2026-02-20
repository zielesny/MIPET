/**
 * MIPET - Mesoscopic Interaction Parameter Estimation with Tinker
 * Copyright (C) 2026  Achim Zielesny (achim.zielesny@googlemail.com)
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
 * Record for job tasks
 * 
 * @author Mirco Daniel
 * @param forcefield_IE_Name: Name of the force field for intermolecular energy calculation
 * @param forcefield_CN_Name: Name of the force field for coordinations number calculation
 * @param particleName1: Name of the first particle
 * @param particleName2: Name of the second particle
 * @param result_IE_PathName: Name of the energy result directory
 * @param result_CN_PathName: Name of the coordination number result directory
 * @param hasEnergieJob: Flag of energy calculation job
 * @param hasCNJob: Flag of coordination number calculation job
 * @param isSameParticle: Flag of wether both particles are same
 * @param isReverse: Flag of wether sequence of particles is reversed 
 */
public record JobTaskRecord(
        String forcefield_IE_Name,
        String forcefield_CN_Name,
        String particleName1,    
        String particleName2,
        String result_IE_PathName,
        String result_CN_PathName,
        boolean hasEnergieJob,
        boolean hasCNJob,
        boolean isSameParticle,
        boolean isReverse) {
}