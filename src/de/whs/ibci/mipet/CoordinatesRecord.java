/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package de.whs.ibci.mipet;

/**
 * Record for storage coordinates datas
 * 
 * @author Mirco Daniel
 */

/**
 * @param coord1: Coordinates of solute particle
 *   [i][j][k]: i: simulation step, j: atom, k: xyz
 * @param coord2: Coordinates of solvent particles
 *   [i][j][k][l]: i: simulation step, j: particle k: atom l: xyz
 */
public record CoordinatesRecord(
        double[][][] coord1,
        double[][][][] coord2) {
}
