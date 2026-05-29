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

import de.whs.ibci.mipet.CoordinatesRecord;
import java.io.IOException;
import java.util.HashMap;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.exception.CDKException;
import de.whs.ibci.mipet.MIPETUtility;
import de.whs.ibci.mipet.TinkerXYZ;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.LinkedList;

/**
 *
 * @author Mirco Daniel
 */
public class TestMIPET4JavaUtility {
    
    private MIPETUtility mipetUtil;
    
    @Before
    public void beforeTest(){
        mipetUtil = new MIPETUtility();
    }

    @Test
    public void testGetVdwVolume() throws InvalidSmilesException, CDKException {
        String smilesString = "";
        double volume;

        try {
            volume = mipetUtil.getVdwVolume(smilesString);
            Assert.fail("IllegalArgumentException not thrown");
        } catch (IllegalArgumentException expected) {
        }
        
        smilesString = "C";        
        volume = mipetUtil.getVdwVolume(smilesString);
        Assert.assertEquals(25.8524433266667, volume, 0.01);
        
        smilesString = "c1ccccc1";
        volume = mipetUtil.getVdwVolume(smilesString);
        Assert.assertEquals(81.1665316528, volume, 0.01);
        
        smilesString = "CC(=O)O";
        volume = mipetUtil.getVdwVolume(smilesString);
        Assert.assertEquals(58.0924226528555, volume, 0.01);        
    }
    
    @Test
    public void testGetAtomicMass() throws InvalidSmilesException, 
            CDKException {
        String smilesString = "";
        double atomicMass;
        
        atomicMass = mipetUtil.getAtomicMass(smilesString, true);
        Assert.assertEquals(0.0, atomicMass, 0.01);
        
        smilesString = "C";        
        atomicMass = mipetUtil.getAtomicMass(smilesString, true);
        Assert.assertEquals(16.04, atomicMass, 0.01);
        
        smilesString = "c1ccccc1";
        atomicMass = mipetUtil.getAtomicMass(smilesString, true);
        Assert.assertEquals(78.11, atomicMass, 0.01);
        
        smilesString = "CC(=O)O";
        atomicMass = mipetUtil.getAtomicMass(smilesString, true);
        Assert.assertEquals(60.05, atomicMass, 0.01);
    }
    
    @Test
    public void testGetAtomicNumberTable() {
        Assert.assertEquals(3, mipetUtil.getAtomicNumberTable()
                .get("Li").intValue());
        Assert.assertEquals(92, mipetUtil.getAtomicNumberTable()
                .get("U").intValue());
    }
    
    @Test
    public void testGetVdWRadii() {
        Assert.assertEquals(1.7, mipetUtil.getVdWRadii()[6], 0.1);
        Assert.assertEquals(1.55, mipetUtil.getVdWRadii()[8], 0.01);
    }
    
    @Test
    public void testMIPET4GetResourceString(){
        String tmpReturnString = 
                mipetUtil.getResourceString(
                        "MIPET.Directory.tinker");
        Assert.assertEquals("tinker", tmpReturnString);
    }
    
    @Test
    public void getSmilesData() throws IOException {
        String tmpSmilesFileName = 
                "Molecules/SMILES/" + "Smiles.dat";
        HashMap<String, String> tmpSmiles = 
                mipetUtil.getSmilesData(tmpSmilesFileName);
        Assert.assertEquals("N1C=CN=C1", tmpSmiles.get("Imidazole"));
    }
    
    @Test
    public void testIsTooClose() {
        double minDistance = 1.0;
        double[][] p1 = new double[3][3];
        p1[0] = new double[]{1.6513234, -0.0589104, 0.4708854};
        p1[1] = new double[]{2.6113234, -0.0589104, 0.4708854};
        p1[2] = new double[]{1.3308688, 0.8460254, 0.4708854};
        double[][] p2 = new double[3][3];
        p2[0] = new double[]{-0.5213109, -0.9148672, -0.4292483};
        p2[1] = new double[]{-0.8066268, -1.6238965, -1.0101719};
        p2[2] = new double[]{-1.1393048, -0.1827312, -0.897485};
        Boolean tmpResult = mipetUtil.isTooClose(p1, p2, minDistance);
        Assert.assertFalse(tmpResult);
        
        minDistance = 2.51;
        tmpResult = mipetUtil.isTooClose(p1, p2, minDistance);
        Assert.assertTrue(tmpResult);
    }
    
    @Test
    public void testSum() {
        double[] tmpArray = new double[]{1.1, 2.2, 3.3, 4.4, 5.5};
        double tmpResult = mipetUtil.sum(tmpArray);
        double tmpExpect = 16.5;
        Assert.assertEquals(tmpExpect, tmpResult, 0.01);
    }
    
    @Test
    public void testProductSum() {
        double[] tmpArray1 = new double[]{1.1, 2.2, 3.3, 4.4, 5.5};
        double[] tmpArray2 = new double[]{1.0, 2.0, 3.0, 4.0, 5.0};
        double tmpResult = mipetUtil.productSum(tmpArray1, tmpArray2);
        double tmpExpect = 60.5;
        Assert.assertEquals(tmpExpect, tmpResult, 0.01);
    }
    
    @Test
    public void testStandardDeviation() {
        int[] tmpValues = new int[]{1, 3, 4, 1, 5, 7};
        double tmpResult = mipetUtil.standarddeviation(tmpValues, 3.5);
        double tmpExpect = 2.3452;
        Assert.assertEquals(tmpExpect, tmpResult, 0.0001);
    }
    
    @Test
    public void testMax() {
        int[] tmpValues = new int[]{1, 3, 4, 1, 5, 7};
        int tmpResult = mipetUtil.getMax(tmpValues);
        int tmpExpect = 7;
        Assert.assertEquals(tmpExpect, tmpResult);
    }
    
    @Test
    public void testMin() {
        int[] tmpValues = new int[]{1, 3, 4, 1, 5, 7};
        int tmpResult = mipetUtil.getMin(tmpValues);
        int tmpExpect = 1;
        Assert.assertEquals(tmpExpect, tmpResult);
    }
    
    @Test
    public void testGetParameterParticleNameList() {
        LinkedList<String> tmpPathName = new LinkedList<>();
        tmpPathName.add("MMFF\\H2O_Me");
        String[][] tmpResult = mipetUtil
                .getParameterParticleNameList(tmpPathName);
        Assert.assertEquals("MMFF", tmpResult[0][0]);
        Assert.assertEquals("H2O", tmpResult[0][1]);
        Assert.assertEquals("Me", tmpResult[0][2]);
        
        tmpPathName = new LinkedList<>();
        tmpPathName.add("input\\MMFF\\H2O_Me");
        tmpResult = mipetUtil
                .getParameterParticleNameList(tmpPathName);
        Assert.assertEquals("MMFF", tmpResult[0][0]);
        Assert.assertEquals("H2O", tmpResult[0][1]);
        Assert.assertEquals("Me", tmpResult[0][2]);
    }
    
    @Test
    public void testGetBoxLength() {
        Path path = Paths.get("testdata/H2O_Et_simulated.arc");
        double tmpBoxLength = mipetUtil.getBoxLength(path);
        Assert.assertEquals(9.070700, tmpBoxLength, 0.000001);
    }
    
    @Test
    public void testGetNeighborCellNumber() {
        int tmpCellsInRow = 6;
        HashSet tmpNeighborIndices;
        HashSet<Integer> tmpCellIndices;
        
        tmpCellIndices = new HashSet<>();
        tmpCellIndices.add(8);
        tmpNeighborIndices = mipetUtil
                .getNeighborCellNumbers(tmpCellsInRow, tmpCellIndices);
        Assert.assertEquals(27, tmpNeighborIndices.size());
        Assert.assertTrue(tmpNeighborIndices.contains(1));
        Assert.assertTrue(tmpNeighborIndices.contains(2));
        Assert.assertTrue(tmpNeighborIndices.contains(3));
        Assert.assertTrue(tmpNeighborIndices.contains(7));
        Assert.assertTrue(tmpNeighborIndices.contains(8));
        Assert.assertTrue(tmpNeighborIndices.contains(9));
        Assert.assertTrue(tmpNeighborIndices.contains(13));
        Assert.assertTrue(tmpNeighborIndices.contains(14));
        Assert.assertTrue(tmpNeighborIndices.contains(15));
        Assert.assertTrue(tmpNeighborIndices.contains(37));
        Assert.assertTrue(tmpNeighborIndices.contains(38));
        Assert.assertTrue(tmpNeighborIndices.contains(39));
        Assert.assertTrue(tmpNeighborIndices.contains(43));
        Assert.assertTrue(tmpNeighborIndices.contains(44));
        Assert.assertTrue(tmpNeighborIndices.contains(45));
        Assert.assertTrue(tmpNeighborIndices.contains(49));
        Assert.assertTrue(tmpNeighborIndices.contains(50));
        Assert.assertTrue(tmpNeighborIndices.contains(51));
        Assert.assertTrue(tmpNeighborIndices.contains(181));
        Assert.assertTrue(tmpNeighborIndices.contains(182));
        Assert.assertTrue(tmpNeighborIndices.contains(183));
        Assert.assertTrue(tmpNeighborIndices.contains(187));
        Assert.assertTrue(tmpNeighborIndices.contains(188));
        Assert.assertTrue(tmpNeighborIndices.contains(189));
        Assert.assertTrue(tmpNeighborIndices.contains(193));
        Assert.assertTrue(tmpNeighborIndices.contains(194));
        Assert.assertTrue(tmpNeighborIndices.contains(195));
        
        tmpCellIndices = new HashSet<>();
        tmpCellIndices.add(1);
        tmpCellIndices.add(2);
        tmpNeighborIndices = mipetUtil
                .getNeighborCellNumbers(tmpCellsInRow, tmpCellIndices);
        Assert.assertEquals(36, tmpNeighborIndices.size());
        Assert.assertTrue(tmpNeighborIndices.contains(12));
        Assert.assertTrue(tmpNeighborIndices.contains(6));
        Assert.assertTrue(tmpNeighborIndices.contains(36));
        Assert.assertTrue(tmpNeighborIndices.contains(31));
        Assert.assertTrue(tmpNeighborIndices.contains(32));
        Assert.assertTrue(tmpNeighborIndices.contains(33));
        Assert.assertTrue(tmpNeighborIndices.contains(38));
        Assert.assertTrue(tmpNeighborIndices.contains(42));
        Assert.assertTrue(tmpNeighborIndices.contains(72));
        Assert.assertTrue(tmpNeighborIndices.contains(67));
        Assert.assertTrue(tmpNeighborIndices.contains(68));
        Assert.assertTrue(tmpNeighborIndices.contains(69));
        Assert.assertTrue(tmpNeighborIndices.contains(216));
        Assert.assertTrue(tmpNeighborIndices.contains(211));
        Assert.assertTrue(tmpNeighborIndices.contains(212));
        Assert.assertTrue(tmpNeighborIndices.contains(213));
    }
    
    @Test
    public void testGetMolecularDiameter() {
        double tmpResult;
        String tmpForceField;
        String tmpParticle;
        String tmpFileName;
        TinkerXYZ tmpTinkerXYZ;
        
        tmpForceField = "MMFF";
        tmpParticle = "Me";
        tmpFileName = "testdata/Me.txyz";
        tmpTinkerXYZ = new TinkerXYZ(tmpForceField, tmpParticle, tmpFileName);
        tmpResult = mipetUtil.getMolecularDiameter(tmpTinkerXYZ);
        Assert.assertEquals(4.148, tmpResult, 0.001);
    }
    
}
