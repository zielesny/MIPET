package de.whs.ibci.matrixrotation;

import de.whs.ibci.mipet.MatrixUtil;
import de.whs.ibci.mipet.RotationUtil;
import de.whs.ibci.mipet.VectorUtil;
import org.junit.Assert;
import org.junit.Test;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Test class for the rotation methods
 */
public class RotationModelTest {

    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(RotationModelTest.class.getName());

    //<editor-fold defaultstate="collapsed" desc="NxN test methods">
    @Test
    public void testMatrices4x4() {
        this.testMatricesNxN(4, 1e-15);
    }
    @Test
    public void testMatrices9x9() {
        this.testMatricesNxN(9, 1e-14);
    }
    @Test
    public void testMatrices100x100() {
        this.testMatricesNxN(100, 1e-15);
    }
    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="NxNxM test methods">
    @Test
    public void testMatrices4x4x4() {
        this.testMatricesNxNxM(4, 4, 1e-15);
    }
    @Test
    public void testMatrices4x4x8() {
        this.testMatricesNxNxM(4, 8, 1e-15);
    }
    @Test
    public void testMatrices4x4x16() {
        this.testMatricesNxNxM(4, 16, 1e-15);
    }
    @Test
    public void testMatrices9x9x4() {
        this.testMatricesNxNxM(9, 4, 1e-14);
    }
    @Test
    public void testMatrices9x9x8() {
        this.testMatricesNxNxM(9, 8, 1e-14);
    }
    @Test
    public void testMatrices9x9x16() {
        this.testMatricesNxNxM(9, 16, 1e-14);
    }
    @Test
    public void testMatrices100x100x4() {
        this.testMatricesNxNxM(100, 4, 1e-14);
    }
    @Test
    public void testMatrices100x100x8() {
        this.testMatricesNxNxM(100, 8, 1e-14);
    }
    @Test
    public void testMatrices100x100x16() {
        this.testMatricesNxNxM(100, 16, 1e-14);
    }
    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="test other methods">
    @Test
    public void testMultiply() {
        double[][] tmpRotMatrix = new double[3][3];
        double[] tmpVector = new double[3];
        double[] tmpExpect;
        double[] tmpResult = new double[3];
        
        tmpRotMatrix[0][0] = 0.788675;
        tmpRotMatrix[0][1] = -0.211325;
        tmpRotMatrix[0][2] = -0.57735;
        tmpRotMatrix[1][0] = -0.211325;
        tmpRotMatrix[1][1] = 0.788675;
        tmpRotMatrix[1][2] = -0.57735;
        tmpRotMatrix[2][0] = 0.57735;
        tmpRotMatrix[2][1] = 0.57735;
        tmpRotMatrix[2][2] = 0.57735;
        tmpVector[0] = 1.0;
        tmpVector[1] = 2.0;
        tmpVector[2] = 3.0;
        tmpExpect = new double[]{-1.36603, -0.366025, 3.46410};
        tmpResult = MatrixUtil.multiply(tmpRotMatrix, tmpVector);
        
        Assert.assertEquals(tmpExpect[0], tmpResult[0], 0.00001);
        Assert.assertEquals(tmpExpect[1], tmpResult[1], 0.00001);
        Assert.assertEquals(tmpExpect[2], tmpResult[2], 0.00001);
    }
    
    @Test
    public void testMoveX() {
        double[][] tmpVector = new double[1][3];
        double[][] tmpExpect = new double[1][3];
        double[][] tmpResult = new double[1][3];
        
        tmpVector[0][0] = 0.788675;
        tmpVector[0][1] = -0.211325;
        tmpVector[0][2] = -0.57735;
        tmpExpect[0][0] = 1.788675;
        tmpExpect[0][1] = -0.211325;
        tmpExpect[0][2] = -0.57735;
        tmpResult = VectorUtil.moveX(tmpVector, 1);
        
        Assert.assertEquals(tmpExpect[0][0], tmpResult[0][0], 0.000001);
        Assert.assertEquals(tmpExpect[0][1], tmpResult[0][1], 0.000001);
        Assert.assertEquals(tmpExpect[0][2], tmpResult[0][2], 0.000001);
    }
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="private methods">
    /**
     * Method which test a rotation matrix calculated with given sphere nodes and rotation number
     *
     * @param aSphereNodeNumber Number of sphere nodes
     * @param aRotationNumber Number of rotations
     * @param aDelta Accepted difference between expected and actual value
     */
    private void testMatricesNxNxM(int aSphereNodeNumber, int aRotationNumber, 
            double aDelta) {
        String tmpFileNameSphereNode = 
                "de/whs/ibci/mipet4java/sphereNodes/SphereNodes" 
                + aSphereNodeNumber + ".txt";
        String tmpFileNameRotationMatrices = 
                "de/whs/ibci/mipet4java/expectedRotationMatrices/expectedRotationMatrices" 
                + aSphereNodeNumber + "x" 
                + aSphereNodeNumber + "x" 
                + aRotationNumber + ".txt";
        List<double[]> tmpSphereNodeList = RotationUtil.readSphereNodes(tmpFileNameSphereNode);
        List<double[][]> tmpActualRotationMatrices = RotationUtil.getRotationMatrices2(tmpSphereNodeList, new double[] {-1.0, 0.0, 0.0}, aRotationNumber);
        List<double[][]> tmpExpectedRotationMatrices = new ArrayList<>();
        try {
            ClassLoader tmpClassLoader = this.getClass().getClassLoader(); //new
            File tmpFile = new File(tmpClassLoader.getResource(tmpFileNameRotationMatrices).getFile()); //new
            BufferedReader tmpReader = new BufferedReader(new FileReader(tmpFile));
            String tmpLine = tmpReader.readLine();
            while(tmpLine != null) {
                int tmpStartOfSubString = tmpLine.indexOf("{");
                int tmpEndOfSubString = tmpLine.lastIndexOf("}");
                String tmpNewSubString = tmpLine.substring(tmpStartOfSubString + 2, tmpEndOfSubString - 1);
                double[][] tmpDoubleRowTokens = new double[3][3];
                String[] tmpRowTokens = tmpNewSubString.split("\\},\\{");
                int i = 0;
                for (String tmpRowToken : tmpRowTokens) {
                    String[] tmpColumnTokens = tmpRowToken.split(",");
                    double[] tmpDoubleColumnTokens = new double[tmpColumnTokens.length];
                    int j = 0;
                    for (String tmpColumnToken : tmpColumnTokens) {
                        tmpDoubleColumnTokens[j] = Double.parseDouble(tmpColumnToken);
                        j++;
                    }
                    tmpDoubleRowTokens[i] = tmpDoubleColumnTokens;
                    i++;
                }
                tmpExpectedRotationMatrices.add(tmpDoubleRowTokens);
                tmpLine = tmpReader.readLine();
            }
        } catch (Exception anException) {
            RotationModelTest.LOGGER.log(Level.SEVERE, anException.toString());
        }
        this.checkMatrices(tmpExpectedRotationMatrices, tmpActualRotationMatrices, aDelta);
    }

    /**
     * Method which tests a rotation matrix calculated with given sphere nodes
     *
     * @param aSphereNodeNumber Number of sphere nodes
     * @param aDelta Accepted difference between expected and actual value
     */
    private void testMatricesNxN(int aSphereNodeNumber, double aDelta) {
        String tmpFileNameSphereNode;
        String tmpFileNameRotationMatrices;
        List<double[]> tmpSphereNodeList;
        List<double[][]> tmpActualRotationMatrices;
        List<double[][]> tmpExpectedRotationMatrices;
        
        tmpFileNameSphereNode = "de/whs/ibci/mipet4java/sphereNodes/SphereNodes"
                + aSphereNodeNumber + ".txt";
        tmpFileNameRotationMatrices = "de/whs/ibci/mipet4java/expectedRotationMatrices/expectedRotationMatrices" 
                + aSphereNodeNumber + "x" + aSphereNodeNumber + ".txt";
        tmpSphereNodeList = RotationUtil.readSphereNodes(tmpFileNameSphereNode);
        tmpActualRotationMatrices = RotationUtil
                .getRotationMatrices1(tmpSphereNodeList, 
                        new double[] {1.0, 0.0, 0.0});
        tmpExpectedRotationMatrices = new ArrayList<>();
        
        try {
            ClassLoader tmpClassLoader = this.getClass().getClassLoader(); //new
            File tmpFile = new File(tmpClassLoader.getResource(tmpFileNameRotationMatrices).getFile()); //new
            BufferedReader tmpReader = new BufferedReader(new FileReader(tmpFile));
            String tmpLine = tmpReader.readLine();
            while(tmpLine != null) {
                int tmpStartOfSubString = tmpLine.indexOf("{");
                int tmpEndOfSubString = tmpLine.lastIndexOf("}");
                String tmpNewSubString = tmpLine.substring(tmpStartOfSubString + 2, tmpEndOfSubString - 1);
                double[][] tmpDoubleRowTokens = new double[3][3];
                String[] tmpRowTokens = tmpNewSubString.split("\\},\\{");
                int i = 0;
                for (String tmpRowToken : tmpRowTokens) {
                    String[] tmpColumnTokens = tmpRowToken.split(",");
                    double[] tmpDoubleColumnTokens = new double[tmpColumnTokens.length];
                    int j = 0;
                    for (String tmpColumnToken : tmpColumnTokens) {
                        tmpDoubleColumnTokens[j] = Double.parseDouble(tmpColumnToken);
                        j++;
                    }
                    tmpDoubleRowTokens[i] = tmpDoubleColumnTokens;
                    i++;
                }
                tmpExpectedRotationMatrices.add(tmpDoubleRowTokens);
                tmpLine = tmpReader.readLine();
            }
        } catch (Exception anException) {
            RotationModelTest.LOGGER.log(Level.SEVERE, anException.toString());
        }
        this.checkMatrices(tmpExpectedRotationMatrices, tmpActualRotationMatrices, aDelta);
    }

    /**
     * Checks an expected rotation matrix with the actual matrix
     *
     * @param anExpectedMatrix Expected rotation matrix
     * @param anActualMatrix Actual rotation matrix
     * @param aDelta Accepted difference between actual and expected value
     */
    private void checkMatrices(List<double[][]> anExpectedMatrix, List<double[][]> anActualMatrix , double aDelta) {
        for(int i = 0; i < anExpectedMatrix.size(); i++) {
            for (int j = 0; j < anExpectedMatrix.get(i).length; j++) {
                double[] tmpExpectedRow  = anExpectedMatrix.get(i)[j];
                double[] tmpActualRow = anActualMatrix.get(i)[j];
                Assert.assertArrayEquals(tmpExpectedRow, tmpActualRow, aDelta);
            }
        }
    }
    //</editor-fold>
}