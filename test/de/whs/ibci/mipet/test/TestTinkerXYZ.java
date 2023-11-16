/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/UnitTests/JUnit4TestClass.java to edit this template
 */
package de.whs.ibci.mipet.test;

import de.whs.ibci.mipet.TinkerXYZ;
import java.io.IOException;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test TinkerXYZ class
 * @author Mirco Daniel
 */
public class TestTinkerXYZ {
    
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
        TinkerXYZ tmpMethane = null;
        
        tmpMethane = new TinkerXYZ("./testdata/Me.txyz");
        double[] tmpCentre;
        tmpCentre= tmpMethane.findCentreCoordinate(tmpMethane.
                getCoordinateList1()[0]);
        Assert.assertEquals(-0.8588, tmpCentre[0], 0.000001);
        Assert.assertEquals(0.5892, tmpCentre[1], 0.000001);
        Assert.assertEquals(0.000, tmpCentre[2], 0.000001);
    }
    
    @Test
    public void TestMoveCoordinates() {
        TinkerXYZ tmpMethane = null;
        tmpMethane = new TinkerXYZ("./testdata/Me.txyz");
        double[] tmpCentre;
        tmpCentre = tmpMethane.findCentreCoordinate(tmpMethane.
                getCoordinateList1()[0]);
        double[][] tmpMoved;
        tmpMoved = tmpMethane.moveCoordinates(tmpMethane.getCoordinateList1()[0], tmpCentre);
        Assert.assertEquals(0.0, tmpMoved[0][0], 0.001);
        Assert.assertEquals(0.0, tmpMoved[0][1], 0.001);
        Assert.assertEquals(0.0, tmpMoved[0][2], 0.001);
    }
}
