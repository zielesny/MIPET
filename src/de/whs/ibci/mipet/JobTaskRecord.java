/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Record.java to edit this template
 */
package de.whs.ibci.mipet;

import java.util.HashSet;

/**
 * Record for storage job tasks
 * 
 * @author Mirco Daniel
 */

/**
 * @param forcefieldName: Name of the force field
 * @param particleName1: Name of the first particle
 * @param particleName2: Name of the second particle
 * @param result_IE_PathName: Name of the energy result directory
 * @param result_CN_PathName: Name of the coordination number result directory
 * @param hasSpecialForcefield_IE: Flag if there is a special force field file for intermolecular energy calculation
 * @param hasSpecialForcefield_CN: Flag if there is a special force field file for coordination number calculation
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
