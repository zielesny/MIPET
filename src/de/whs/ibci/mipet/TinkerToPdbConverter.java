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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * TinkerToPdbConverter Class
 * 
 * @author Mirco Daniel
 */
public class TinkerToPdbConverter {

    public static void convert(String inputFilePath, String outputFilePath) {
        try (BufferedReader reader = 
                new BufferedReader(new FileReader(inputFilePath));
            BufferedWriter writer =
                    new BufferedWriter(new FileWriter(outputFilePath))) {

            String tmpLine = reader.readLine();
            if (tmpLine == null) {
                System.err.println("The tinker xyz input file is empty.");
                return;
            }

            int tmpNumAtoms; // Unnecessary, but maybe for the future implementation
            int tmpAtomCounter;
            int tmpAtomId;
            String tmpAtomName;
            String tmpElement;
            String tmpPdbLine;
            String tmpRecordType;
            String tmpPdbElementColumn;
            String[] tmpFirstLineTokens;
            String[] tmpTokens;
            StringBuilder tmpConectBuilder;
            
            // First line: Number of atoms 
            tmpFirstLineTokens = tmpLine.trim().split("\\s+");
            tmpNumAtoms = Integer.parseInt(tmpFirstLineTokens[0]);
            tmpAtomCounter = 1;
            List<String> tmpConectLines = new ArrayList<>();
            List<Integer> tmpBondedAtoms;

            while ((tmpLine = reader.readLine()) != null) {
                tmpLine = tmpLine.trim();
                if (tmpLine.isEmpty()) continue;

                tmpTokens = tmpLine.split("\\s+");
                
                // Tinker Format: [ID] [Name] [X] [Y] [Z] [Typ] [Connected atom1] [Connected atom2] ...
                if (tmpTokens.length >= 6) { 
                    tmpAtomId = Integer.parseInt(tmpTokens[0]); // The ID of current atom
                    tmpAtomName = tmpTokens[1];
                    double x = Double.parseDouble(tmpTokens[2]);
                    double y = Double.parseDouble(tmpTokens[3]);
                    double z = Double.parseDouble(tmpTokens[4]);

                    // Isolate element symbol
                    tmpElement = tmpAtomName.replaceAll("[^A-Za-z]", "");
                    if (tmpElement.length() > 2) {
                        tmpElement = tmpElement.substring(0, 2);
                    }

                    // --- For lone pairs and dummy atoms ---
                    tmpRecordType = "ATOM  ";

                    if (tmpElement.equals("LP") || tmpElement.equals("M")) {
                        tmpRecordType = "HETATM"; // Mark as not-standard-atom
                        tmpPdbElementColumn = " X"; // Generic element "X", so that the viewer not crashes
                    } else {
                        tmpPdbElementColumn = String.format("%2s", tmpElement); // normal element
                    }
                    
                    // PDB-format for the atom (Strict column width)
                    String pdbLine = String.format(Locale.US,
                            "ATOM  %5d %-4s MOL A   1    %8.3f%8.3f%8.3f  1.00  0.00          %2s",
                            tmpAtomId, // Use the ID from the Tinker file
                            formatPdbAtomName(tmpAtomName),
                            x,
                            y,
                            z,
                            String.format("%2s", tmpElement).toUpperCase()
                    );

                    writer.write(pdbLine);
                    writer.newLine();

                    // --- Handle CONECT ---
                    tmpBondedAtoms = new ArrayList<>();
                    
                    for (int i = 6; i < tmpTokens.length; i++) {
                        tmpBondedAtoms.add(Integer.valueOf(tmpTokens[i]));
                    }

                    // PDB allows max. 4 bonds pro CONECT-line
                    // If an atom has 5 bonds, we need two lines.
                    for (int i = 0; i < tmpBondedAtoms.size(); i += 4) {
                        tmpConectBuilder = new StringBuilder("CONECT");
                        
                        // The start atom (column 7-11)
                        tmpConectBuilder.append(String
                                .format(Locale.US, "%5d", tmpAtomId));
                        
                        // Insert up to 4 atoms (column 12-16, 17-21, 22-26, 27-31)
                        for (int j = 0; j < 4 
                                && (i + j) < tmpBondedAtoms.size(); j++) {
                            tmpConectBuilder.append(
                                    String.format(Locale.US, "%5d", 
                                            tmpBondedAtoms.get(i + j)));
                        }
                        
                        tmpConectLines.add(tmpConectBuilder.toString());
                    }

                    tmpAtomCounter++;
                }
            }
            
            // Write all saved CONECT-lines at the end of the file
            for (String tmpConectLine : tmpConectLines) {
                writer.write(tmpConectLine);
                writer.newLine();
            }

            // End-tag for pdb file
            writer.write("END");
            writer.newLine();
        } catch (IOException e) {
            System.err.println("Error during reading/writing of the file: " 
                    + e.getMessage());
        } catch (NumberFormatException e) {
            System.err.println("Format error. Please, check if the file is in tinker-xyz format. " + e.getMessage());
        }
    }
      
    /**
     * Helping method for the alignment of PDB-atom name
     * @param aName: atom name
     * @return Formated atom name
     */
    private static String formatPdbAtomName(String aName) {
        if (aName.length() == 1 || aName.length() == 2 || aName.length() == 3) {
            return " " + String.format("%-3s", aName); 
        }
        return aName.substring(0, Math.min(aName.length(), 4));
    }
}