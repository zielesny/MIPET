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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.stream.ImageOutputStream;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Utility class for chart generation
 * 
 * @author Achim Zielesny
 */
public class ChartUtil {
    
    /**
     * Tiny threshold
     */
    private final static double TINY_THRESHOLD = 1E-6;
    /**
     * Name of data series
     */
    private final static String DATA_SERIES_NAME = "XyData";
    /**
     * Whitespace pattern
     */
    private final static Pattern WHITESPACE_PATTERN = Pattern.compile("\\s+");
    /**
     * Full HD image size
     */
    private final static int[] FULL_HD_IMAGE_SIZE = new int[] {1920, 1080};
    
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger
            .getLogger(MIPET.class.getName());
    
    /**
     * Constructor
     */
    public ChartUtil() {
        // Nothing to initialize
    }

    /**
     * Creates energy graphics files in Full HD resolution in folder of 
     * aEnergyDataPathName.
     * 
     * @param aEnergyDataPathName File pathname with energy data
     * @param anEnergyGraphicsPrefix Prefix of energy graphics filenames
     * @return True: Operation was successful, false: Otherwise
     */
    public boolean createEnergyGraphics(
        String aEnergyDataPathName, 
        String anEnergyGraphicsPrefix
    ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aEnergyDataPathName == null || aEnergyDataPathName.isEmpty()) {
            return false;
        }
        if (anEnergyGraphicsPrefix == null || anEnergyGraphicsPrefix.isEmpty()) {
            return false;
        }
        if (!(new File(aEnergyDataPathName)).isFile()) {
            return false;
        }
        
        // </editor-fold>
        
        try {
            // <editor-fold defaultstate="collapsed" desc="Read raw energy data">
            String[][] energyData = this
                    .readJaggedStringArrayFromFile(aEnergyDataPathName);
            if (energyData == null || energyData.length < 2) {
                return false;
            }
            // </editor-fold>
            
            // <editor-fold defaultstate="collapsed" desc="Set diagram labels">
            if (energyData[0].length != 6) {
                return false;
            }
            String xAxisLabel = energyData[0][0] + " " + energyData[0][1];
            String yAxisLabel1 = energyData[0][2] + " " + energyData[0][3];
            String yAxisLabel2 = energyData[0][4] + " " + energyData[0][5];
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Convert energy data">
            // Convert energy data
            double[] distances = new double[energyData.length - 1];
            double[] eMinValues = new double[energyData.length - 1];
            double[] eAverageValues = new double[energyData.length - 1];
            for (int i = 1; i < energyData.length; i++) {
                distances[i - 1] = Double.parseDouble(energyData[i][0]);
                eMinValues[i - 1] = Double.parseDouble(energyData[i][1]);
                eAverageValues[i - 1] = Double
                        .parseDouble(energyData[i][2]);
            }
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Create and save chart 1: Distance vs. Emin">
            String diagramLabel1 = yAxisLabel1 + " vs. " + xAxisLabel;
            JFreeChart eMinVsDistanceChart = 
                this.createXyChart(
                    distances, 
                    eMinValues, 
                    diagramLabel1, 
                    xAxisLabel, 
                    yAxisLabel1, 
                    true,  // anIsThickLines
                    true,  // anIsShapePaint
                    true,  // anIsFillColorWhite
                    false  // anIsOutlinePaintWhite
                );
            if (eMinVsDistanceChart == null) {
                return false;
            }
            BufferedImage eMinVsDistanceImage = 
                this.getImage(
                    eMinVsDistanceChart, 
                    FULL_HD_IMAGE_SIZE[0], 
                    FULL_HD_IMAGE_SIZE[1]
                );
            if (eMinVsDistanceImage == null) {
                return false;
            }
            String eMinVsDistanceImagePathname = 
                new File(aEnergyDataPathName).getParent() + 
                File.separatorChar + 
                anEnergyGraphicsPrefix + 
                "DistanceVsEmin.jpg";
            if (!ChartUtil.writeJpegImageToFile(
                    eMinVsDistanceImage, 
                    new File(eMinVsDistanceImagePathname)
                )
            ) {
                return false;
            }
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Create and save chart 1: Distance vs. Eaverage">
            String diagramLabel2 = yAxisLabel2 + " vs. " + xAxisLabel;
            JFreeChart eAverageVsDistanceChart = 
                this.createXyChart(
                    distances, 
                    eAverageValues, 
                    diagramLabel2, 
                    xAxisLabel, 
                    yAxisLabel2, 
                    true,  // anIsThickLines
                    true,  // anIsShapePaint
                    true,  // anIsFillColorWhite
                    false  // anIsOutlinePaintWhite
                );
            if (eAverageVsDistanceChart == null) {
                return false;
            }
            BufferedImage eAverageVsDistanceImage = 
                this.getImage(
                    eAverageVsDistanceChart, 
                    FULL_HD_IMAGE_SIZE[0], 
                    FULL_HD_IMAGE_SIZE[1]
                );
            if (eAverageVsDistanceImage == null) {
                return false;
            }
            String eAverageVsDistanceImagePathname = 
                new File(aEnergyDataPathName).getParent() + 
                File.separatorChar + 
                anEnergyGraphicsPrefix + 
                "DistanceVsEaverage.jpg";
            // </editor-fold>
            return ChartUtil.writeJpegImageToFile(
                    eAverageVsDistanceImage,
                    new File(eAverageVsDistanceImagePathname));
        } catch (NumberFormatException anException) {
            LOGGER.log(Level.SEVERE, 
                    "IOException during creating chart.", anException);
            return false;
        }
    }
    
    /**
     * Creates XY chart from x and y values.
     *
     * @param aXValues X values
     * @param aYValues Y values
     * @param aTitle Title of diagram
     * @param aXAxisLabel Label of x-axis
     * @param aYAxisLabel Label of y-axis
     * @param anIsThickLines True: Thick lines, false: Thin lines
     * @param anIsShapePaint True: Shapes are painted, false: Shapes are not
     *   painted
     * @param anIsFillColorWhite If shapes are painted: True: Fill color of
     *   shapes is white, false: Fill color of shapes is default color
     * @param anIsOutlinePaintWhite If shapes are painted: True: Outline of
     *   shapes is white, false: Outline of shapes is default color
     * @return Chart or null if chart could not be created
     */
    private JFreeChart createXyChart(
        double[] aXValues, 
        double[] aYValues,
        String aTitle,
        String aXAxisLabel,
        String aYAxisLabel,
        boolean anIsThickLines,
        boolean anIsShapePaint, 
        boolean anIsFillColorWhite, 
        boolean anIsOutlinePaintWhite
    ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aXValues == null || aYValues == null || aXValues.length != aYValues.length) {
            return null;
        }
        // </editor-fold>
        try {
            // <editor-fold defaultstate="collapsed" desc="- Create data series">
            XYSeries xySeries = new XYSeries(DATA_SERIES_NAME);

            double xMin;
            double xMax;
            double yMin;
            double yMax;
            xMin = Double.MAX_VALUE;
            xMax = -Double.MAX_VALUE;
            yMin = Double.MAX_VALUE;
            yMax = -Double.MAX_VALUE;
            for (int i = 0; i < aXValues.length; i++) {
                // Safeguard for NaN values in value item
                try {
                    if (Double.isNaN(aXValues[i]) || Double.isNaN(aYValues[i])) {
                        continue;
                    }
                } catch (Exception anException) {
                    continue;
                }
                xMin = Math.min(xMin, aXValues[i]);
                xMax = Math.max(xMax, aXValues[i]);
                yMin = Math.min(yMin, aYValues[i]);
                yMax = Math.max(yMax, aYValues[i]);
                xySeries.add(aXValues[i], aYValues[i]);
            }
            double xOffset;
            if (Math.abs(xMax - xMin) < xMax * TINY_THRESHOLD) {
                xOffset = (xMax + xMin) * 0.5 * 0.05;
            } else {
                xOffset = (xMax - xMin) * 0.05;
            }
            double x1 = xMin - xOffset;
            double x2 = xMax + xOffset;
            double yOffset;
            if (Math.abs(yMax - yMin) < yMax * TINY_THRESHOLD) {
                yOffset = (yMax + yMin) * 0.5 * 0.05;
            } else {
                yOffset = (yMax - yMin) * 0.05;
            }
            double y1 = yMin - yOffset;
            double y2 = yMax + yOffset;

            XYDataset xyDataset = new XYSeriesCollection(xySeries);
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="- Create chart">
            JFreeChart xyChart = 
                ChartFactory.createXYLineChart(
                    aTitle, // Title
                    aXAxisLabel, // xAxisLabel
                    aYAxisLabel, // yAxisLabel
                    xyDataset, // dataset
                    PlotOrientation.VERTICAL, // orientation
                    false, // legend flag
                    false, // tooltips flag
                    false  // URLs flag
                ); 
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="- Set axis ranges">
            XYPlot xyPlot = (XYPlot) xyChart.getPlot();
            ValueAxis xAxis = xyPlot.getDomainAxis();
            xAxis.setRange(x1, x2);
            ValueAxis yAxis = xyPlot.getRangeAxis();
            yAxis.setRange(y1, y2);
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="- Set outline, shapes and thickness">
            XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) xyPlot.getRenderer();
            renderer.setPaint(Color.BLACK);
            if (anIsThickLines) {
                renderer.setStroke(new BasicStroke(3f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL));
            }
            if (anIsShapePaint) {
                renderer.setShape(new Ellipse2D.Float(-5.0f, -5.0f, 10.0f, 10.0f));
                renderer.setShapesVisible(true);
                renderer.setShapesFilled(true);
                if (anIsFillColorWhite) {
                    renderer.setUseFillPaint(true);
                    renderer.setFillPaint(Color.white);
                } else if (anIsOutlinePaintWhite) {
                    renderer.setOutlinePaint(Color.white);
                    renderer.setUseOutlinePaint(true);
                }
            }
            // </editor-fold>
            return xyChart;
        } catch (Exception anException) {
            return null;
        }
    }

    /**
     * Returns an image of the chart
     * 
     * @param aChart Chart
     * @param aWidth Width of the image
     * @param aHeight Height of the image
     * @return An image of the chart or null if image can not be created
     */
    private BufferedImage getImage(
        JFreeChart aChart, 
        int aWidth, 
        int aHeight
    ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aChart == null || aWidth <= 0 || aHeight <= 0) {
            return null;
        }
        // </editor-fold>
        try {
            return aChart.createBufferedImage(aWidth, aHeight);
        } catch (Exception anException) {
            return null;
        }
    }

    /**
     * Reads jagged string array from file. Each line is splitted after one or
     * more whitespace characters.
     *
     * @param aSourceFilePathname Full pathname of source (may be null then null
     * is returned)
     * @return Jagged string array or null if jagged string array could not be
     * read
     */
    private String[][] readJaggedStringArrayFromFile(
            String aSourceFilePathname) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSourceFilePathname == null || aSourceFilePathname.isEmpty()) {
            return null;
        }
        // </editor-fold>
        
        Path filePath = Paths.get(aSourceFilePathname);
        if (!Files.isRegularFile(filePath)) {
            return null;
        }
        try (BufferedReader reader = Files.newBufferedReader(filePath)){
            LinkedList<String[]> linkedList = new LinkedList<>();
            String line;
            
            while ((line = reader.readLine()) != null) {
                String[] items = this.splitAndTrim(line.trim());
                if (items != null) {
                    linkedList.add(items);
                }
            }
            
            return linkedList.toArray(String[][]::new);
        } catch (IOException anException) {
            return null;
        } 
    }
    
    /**
     * Splits a string into a string array after one or more whitespace
     * characters.
     *
     * @param aString String to split
     * @return A string array containing the trimmed split strings or null if
     * aString was null or empty
     */
    private String[] splitAndTrim(
        String aString
    ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aString == null || aString.isEmpty()) {
            return null;
        }
        // </editor-fold>
        String[] items = WHITESPACE_PATTERN.split(aString.trim());
        if (items == null || items.length == 0) {
            return null;
        }
        // NOTE: Trim-operation is not necessary since all whitespace characters 
        // are removed
        return items;
    }

    /**
     * Writes JPEG image with defined quality to file. NOTE: ImageIO is NOT
     * thread safe! (thus synchronized)
     *
     * @param anImage Image
     * @param aFile File object
     * @return True: Operation was successful, false: Otherwise
     */
    private static synchronized boolean writeJpegImageToFile(
        BufferedImage anImage, 
        File aFile
    ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (anImage == null) {
            return false;
        }
        if (aFile == null) {
            return false;
        }
        try {
            if (aFile.isFile()) {
                if (!aFile.delete()) {
                    return false;
                }
            }
        } catch (Exception anException) {
            return false;
        }
        // </editor-fold>
        try {
            Iterator<ImageWriter> iterator = ImageIO
                    .getImageWritersByFormatName("jpg");
            if (iterator.hasNext()) {
                ImageWriter imageWriter = iterator.next();
                try (ImageOutputStream ioStream = ImageIO
                        .createImageOutputStream(aFile)) {
                    imageWriter.setOutput(ioStream);
                    ImageWriteParam parameters = imageWriter
                            .getDefaultWriteParam();
                    parameters.setCompressionMode(
                            ImageWriteParam.MODE_EXPLICIT);
                    // Parameter 1.0f: High image quality
                    parameters.setCompressionQuality(1.0f);
                    imageWriter.write(null, new IIOImage(anImage, null, null), 
                            parameters);
                    ioStream.flush();
                    imageWriter.dispose();
                }
                return true;
            } else {
                return false;
            }
        } catch (IOException anException) {
            return false;
        }
    }
    
}