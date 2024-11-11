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
package de.whs.ibci.mipet.test;

import de.whs.ibci.mipet.JobTaskRecord;
import de.whs.ibci.mipet.MIPETUtility;
import de.whs.ibci.mipet.TinkerXYZ;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test TinkerXYZ class
 * @author Mirco Daniel
 */
public class TestTinkerXYZ {
    // <editor-fold defaultstate="collapsed" desc="Final private class variables">
    /**
     * Instance object of MIPETUtility
     */
    final private static MIPETUtility MIPETUTIL = new MIPETUtility();
    
    /**
     * Water volume ratio
     */
    final private static double WATERVOLUMERATIO = 30.0036 / 
            MIPETUTIL.getVdwVolume("O");

    // </editor-fold>
    
    /**
     * SPIECE and SMILES lookup table
     */
    private static HashMap<String, String> smiles;
    
    @Test
    @SuppressWarnings("empty-statement")
    public void TestTinkerXYZ() {
        // Test with only one particle
        TinkerXYZ txyz;
        txyz = new TinkerXYZ("./testdata/Me.txyz", 5);
        int atomNumber = txyz.getAtomNumber();
        Assert.assertEquals(5, atomNumber);
        
        double[][][] coord = txyz.getCoordinateList1();
        Assert.assertEquals(-0.859, coord[0][0][0], 0.001);
        Assert.assertEquals(0.589, coord[0][0][1], 0.001);
        Assert.assertEquals(0.000, coord[0][0][2], 0.001);
        Assert.assertEquals(-1.929, coord[0][4][0], 0.001);
        Assert.assertEquals(0.589, coord[0][4][1], 0.001);
        Assert.assertEquals(0.000, coord[0][4][2], 0.001);
        
        String[] tmpExpectedElementList = new String[]{"C", "H", "H", "H", "H"};
        String[] tmpElementList = txyz.getElementList1();
        Assert.assertArrayEquals(tmpExpectedElementList, tmpElementList);
        
        int[] tmpExpectedParameterList = new int[]{1, 5, 5, 5, 5};
        int[] tmpParameterList = txyz.getParameterList1();
        Assert.assertArrayEquals(tmpExpectedParameterList, tmpParameterList);
        
        int[][] tmpExpectedConnectList = new int[][]{{2,3,4,5},{1},{1},{1},{1}};
        int[][] tmpConnectList = txyz.getConnectionList1();
        Assert.assertArrayEquals(tmpExpectedConnectList, tmpConnectList);
        
        // Test with two sorts of particles
        txyz = new TinkerXYZ("./testdata/H2O_Et.txyz", 1, 3, 8);
        atomNumber = txyz.getAtomNumber();
        Assert.assertEquals(75, atomNumber);
        Assert.assertEquals(3, txyz.getAtomSize1());
        Assert.assertEquals(8, txyz.getAtomSize2());
        double[] expectedCoord11 = new double[] {0.0, 0.0, 0.110440};
        double[] expectedCoord12 = new double[]{0.0, 0.777553, -0.448220};
        double[] expectedCoord13 = new double[]{0.0, -0.777553, -0.448220};
        Assert.assertArrayEquals(expectedCoord11, 
                txyz.getCoordinateList1()[0][0], 0.000001);
        Assert.assertArrayEquals(expectedCoord12,
                txyz.getCoordinateList1()[0][1], 0.000001);
        Assert.assertArrayEquals(expectedCoord13, 
                txyz.getCoordinateList1()[0][2], 0.000001);
        
        double[] expectedCoord2_9_8 = 
                new double[]{1.558169, -0.563318, 3.581470};
        Assert.assertArrayEquals(expectedCoord2_9_8, 
                txyz.getCoordinateList2()[0][8][7], 0.000001);
                
        tmpExpectedElementList = new String[]{"O", "H", "H"};
        tmpElementList = txyz.getElementList1();
        Assert.assertArrayEquals(tmpExpectedElementList, tmpElementList);
        tmpExpectedElementList = 
                new String[]{"C", "H", "H", "H", "C", "H", "H", "H"};
        tmpElementList = txyz.getElementList2();
        Assert.assertArrayEquals(tmpExpectedElementList, tmpElementList);
        
        tmpExpectedParameterList = new int[]{36, 37, 37};
        tmpParameterList = txyz.getParameterList1();
        Assert.assertArrayEquals(tmpExpectedParameterList, tmpParameterList);
        tmpExpectedParameterList = new int[]{27, 28, 28, 28, 27, 28, 28, 28};
        tmpParameterList = txyz.getParameterList2();
        Assert.assertArrayEquals(tmpExpectedParameterList, tmpParameterList);
        
        tmpExpectedConnectList = 
                new int[][]{{5, 6, 7, 8}, {4}, {4}, {4},
                    {4, 9, 10, 11}, {8}, {8}, {8}};
        tmpConnectList = txyz.getConnectionList2();
        Assert.assertArrayEquals(tmpExpectedConnectList, tmpConnectList);
        
        // Test make one txyz out of two
        String txyzFileName1 = "./testdata/H2O.txyz";
        String txyzFileName2 = "./testdata/Me.txyz";
        TinkerXYZ tmpTxyz1 = new TinkerXYZ(txyzFileName1);
        TinkerXYZ tmpTxyz2 = new TinkerXYZ(txyzFileName2);
        
        txyz = new TinkerXYZ(tmpTxyz1, tmpTxyz2);
        int tmpExpectedAtomSize = 8;
        int tmpAtomSize = txyz.getAtomNumber();
        Assert.assertEquals(tmpExpectedAtomSize, tmpAtomSize);
        
        tmpExpectedElementList = new String[]{"O", "H", "H"};
        tmpElementList = txyz.getElementList1();
        Assert.assertArrayEquals(tmpExpectedElementList, tmpElementList);
        tmpExpectedElementList = new String[]{"C", "H", "H", "H", "H"};
        tmpElementList = txyz.getElementList2();
        Assert.assertArrayEquals(tmpExpectedElementList, tmpElementList);
        
        // Test .arc file
        txyz = new TinkerXYZ("./testdata/H2O_Et_simulated.arc", 100, 3, 8);
        double[] expectedCoord100_2_8 = 
                new double[]{1.899158, -0.746083, 4.3549006};
        Assert.assertArrayEquals(expectedCoord100_2_8,
                txyz.getCoordinateList2()[99][8][7], 0.00001);
    }
    
    @Test
    public void TestSetCoordinateList1() throws IOException {
        TinkerXYZ tmpTxyz1 = new TinkerXYZ("./testdata/H2O.txyz", 3);
        TinkerXYZ tmpTxyz2 = new TinkerXYZ("./testdata/H2Omod1.txyz", 3);
        double[][][] tmpCoords = tmpTxyz2.getCoordinateList1();
        tmpTxyz1.setCoordinateList1(tmpCoords);
    }
    
    @Test
    public void TestFindCentreCoordinate() {
        double[] tmpCentre;
        
        TinkerXYZ tmpMethane = new TinkerXYZ("./testdata/Me.txyz");
        tmpCentre= tmpMethane.findCentreCoordinate(tmpMethane.
                getCoordinateList1()[0]);
        Assert.assertEquals(-0.8588, tmpCentre[0], 0.000001);
        Assert.assertEquals(0.5892, tmpCentre[1], 0.000001);
        Assert.assertEquals(0.000, tmpCentre[2], 0.000001);
    }
    
    @Test
    public void TestMoveCoordinates() {
        double[] tmpCentre;
        
        TinkerXYZ tmpMethane = new TinkerXYZ("./testdata/Me.txyz");
        tmpCentre = tmpMethane.findCentreCoordinate(tmpMethane.
                getCoordinateList1()[0]);
        double[][] tmpMoved;
        tmpMoved = tmpMethane.moveCoordinates(tmpMethane.getCoordinateList1()[0], tmpCentre);
        Assert.assertEquals(0.0, tmpMoved[0][0], 0.001);
        Assert.assertEquals(0.0, tmpMoved[0][1], 0.001);
        Assert.assertEquals(0.0, tmpMoved[0][2], 0.001);
    }
    
    @Test
    public void TestGetAtomicMassList1() {
        TinkerXYZ tmpMethane = new TinkerXYZ("./testdata/Me.txyz");
        double[] tmpMass = new double[5];
        tmpMass[0] = tmpMethane.getAtomicMassList1()[0];
        tmpMass[1] = tmpMethane.getAtomicMassList1()[1];
        tmpMass[2] = tmpMethane.getAtomicMassList1()[2];
        tmpMass[3] = tmpMethane.getAtomicMassList1()[3];
        tmpMass[4] = tmpMethane.getAtomicMassList1()[4];
        Assert.assertEquals(12.0107358, tmpMass[0], 0.000001);
        Assert.assertEquals(1.0079407, tmpMass[1], 0.000001);
        Assert.assertEquals(1.0079407, tmpMass[2], 0.000001);
        Assert.assertEquals(1.0079407, tmpMass[3], 0.000001);
        Assert.assertEquals(1.0079407, tmpMass[4], 0.000001);
    }
    
    @Test
    public void TestGetCentreOfMass1() {
        TinkerXYZ tmpMethane = new TinkerXYZ("./testdata/Me.txyz");
        double[][] tmpCentreCoordinate;
        tmpCentreCoordinate = tmpMethane.getCentreOfMass1();
        Assert.assertEquals(-0.8589372, tmpCentreCoordinate[0][0], 0.000001);
        Assert.assertEquals(0.58906283, tmpCentreCoordinate[0][1], 0.000001);
        Assert.assertEquals(0.000000, tmpCentreCoordinate[0][2], 0.000001);
    }
    
    @Test
    public void TestGetCentreOfMass2() {
        TinkerXYZ tmpH2O_Et = new TinkerXYZ("./testdata/H2O_Et.txyz", 1, 3, 8);
        double[][] tmpCentreCoord1;
        double[][][] tmpCentreCoords2;
        
        tmpCentreCoord1 = tmpH2O_Et.getCentreOfMass1();
        tmpCentreCoords2 = tmpH2O_Et.getCentreOfMass2();
        Assert.assertEquals(0.0, tmpCentreCoord1[0][0], 0.000001);
        Assert.assertEquals(0.0, tmpCentreCoord1[0][1], 0.000001);
        Assert.assertEquals(0.0479268, tmpCentreCoord1[0][2], 0.000001);
        Assert.assertEquals(-1.32664569, tmpCentreCoords2[0][0][0], 0.000001);
        Assert.assertEquals(-4.419935481, tmpCentreCoords2[0][0][1], 0.000001);
        Assert.assertEquals(2.9358169, tmpCentreCoords2[0][0][2], 0.000001);
        Assert.assertEquals(0.309953149, tmpCentreCoords2[0][8][0], 0.000001);
        Assert.assertEquals(0.313282164, tmpCentreCoords2[0][8][1], 0.000001);
        Assert.assertEquals(3.860684465, tmpCentreCoords2[0][8][2], 0.000001);
    }
    
    @Test
    public void TestgetPBCDistances() {
        TinkerXYZ tmpH2O_Et = new TinkerXYZ("./testdata/H2O_Et.txyz", 1, 3, 8);
        double[][] tmpDistances;
        double tmpBoxLength;
        
        tmpBoxLength = 9.070700;
        tmpDistances = tmpH2O_Et.getPBCDistances(tmpBoxLength);
        Assert.assertEquals(5.443870652, tmpDistances[0][0], 0.000001);
        Assert.assertEquals(3.838142478, tmpDistances[0][8], 0.000001);
    }
    
    @Test
    // This is only for experimental purpose
    public void TestRDF() {
        String tmpParticle2;
        String tmpDatFileName;
        Path tmpTargetDir;
        int tmpSolventNumber;
        int tmpIntervalSize;
        int tmpIntervalIndex;
        int tmpDistancesSize;
        int tmpDistanceIndex;
        int[] tmpDensity;
        double tmpVdWSolventVolume;
        double tmpBoxLength;
        double tmpBulkDensity;
        double tmpIntervalWidth;
        double tmpBoundary;
        double tmpCubic1;
        double tmpCubic2;
        double tmpShellVolume;
        TinkerXYZ tmpMeOH = new TinkerXYZ("Z:/Scratch/MeOH_MeOH.xyz", 1, 6, 6);
        double[][] tmpDistances;
        double[] tmpRDFs;
        
        
        Locale.setDefault(Locale.ENGLISH);
        tmpParticle2 = "MeOH";
        tmpSolventNumber = 5000;
        tmpIntervalIndex = 0;
        tmpDistanceIndex = 0;
        tmpTargetDir = Paths.get("Z:\\Scratch\\");
        tmpDatFileName = "Z:/Scratch/distances.dat";
        tmpIntervalWidth = 1.0;
        
        if (!Files.exists(tmpTargetDir)) {
            try {
                Files.createDirectories(tmpTargetDir);
            } catch (IOException ex) {
            }
        }
        smiles = MIPETUTIL
                .getSmilesData("Molecules/SMILES/" + "Smiles.dat");
        tmpVdWSolventVolume = MIPETUTIL.getVdwVolume(smiles.get(tmpParticle2));
        
        // Calculate water volume ratio - ratio of Vparticle and Vvdw 
        //   of water 1.7297
        tmpBoxLength =  Math.pow(WATERVOLUMERATIO 
                * tmpSolventNumber 
                * tmpVdWSolventVolume, 1.0/3.0);
        
        // Bulk particle density
        tmpBulkDensity = (double)tmpSolventNumber / 
                (tmpBoxLength * tmpBoxLength * tmpBoxLength);
        
        // Determine distances before warmup
        tmpDistances = tmpMeOH.getPBCDistances(tmpBoxLength);
        tmpDistancesSize = tmpDistances[0].length - 1;
        Arrays.sort(tmpDistances[0]);
        tmpIntervalSize = (int)Math.ceil(tmpDistances[0][tmpDistancesSize]);
        tmpDensity = new int[tmpIntervalSize];
        tmpRDFs = new double[tmpIntervalSize];
        tmpBoundary = tmpIntervalWidth;
        
        while (tmpIntervalIndex < tmpIntervalSize && 
                tmpDistanceIndex < tmpDistancesSize) {
            if (tmpDistances[0][tmpDistanceIndex] < tmpBoundary) {
                tmpDensity[tmpIntervalIndex]++;
                tmpDistanceIndex++;
            } else {
                tmpBoundary += tmpIntervalWidth;
                tmpIntervalIndex++;
            }
        }
        
        for (int i = 0; i < tmpIntervalSize; i++) {
            tmpCubic1 = Math.pow(i * tmpIntervalWidth, 3);
            tmpCubic2 = Math.pow((i + 1) * tmpIntervalWidth, 3);
            tmpShellVolume = Math.PI * 4 / 3 * (tmpCubic2 - tmpCubic1);
            tmpRDFs[i] = tmpDensity[i] / tmpShellVolume / tmpBulkDensity;
        }
        
        try (BufferedWriter tmpBW = new BufferedWriter(
                new FileWriter(tmpDatFileName, true))) {
            for (double tmpRDF: tmpRDFs) {
                tmpBW.append(String.format("%.2f", tmpRDF));
                tmpBW.append("\n");
            }
            
        } catch (IOException ex) {
        }
        
        
        // Warmup
        
        // 
    }
}
