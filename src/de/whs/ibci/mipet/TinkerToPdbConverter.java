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
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * TinkerToPdbConverter Class
 * 
 * @author Mirco Daniel
 */
public class TinkerToPdbConverter {
    
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    /**
     * Constructor TinkerToPdbConverter
     */
    public TinkerToPdbConverter() {
    }

    // </editor-fold> 
    
    /**
     * Converts tinker's xyz format to pdb format
     * @param inputPath Tinker's xyz-file
     * @param outputPath Pdb-file
     */
    public static void convert(Path inputPath, Path outputPath) {
        try (BufferedReader reader = Files.newBufferedReader(inputPath);
            BufferedWriter writer =Files.newBufferedWriter(outputPath)) {

            String line = reader.readLine();
            if (line == null) {
                System.err.println("The tinker xyz input file is empty.");
                return;
            }

            int numAtoms; // Unnecessary, but maybe for the future implementation
            int atomCounter;
            int atomId;
            String atomName;
            String element;
            String[] firstLineTokens;
            String[] tokens;
            StringBuilder connectBuilder;
            
            // First line: Number of atoms 
            firstLineTokens = line.trim().split("\\s+");
            numAtoms = Integer.parseInt(firstLineTokens[0]);
            atomCounter = 1;
            List<String> connectLines = new ArrayList<>();
            List<Integer> bondedAtoms;

            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                tokens = line.split("\\s+");
                
                // Tinker Format: [ID] [Name] [X] [Y] [Z] [Typ] [Connected atom1] [Connected atom2] ...
                if (tokens.length >= 6) { 
                    atomId = Integer.parseInt(tokens[0]); // The ID of current atom
                    atomName = tokens[1];
                    double x = Double.parseDouble(tokens[2]);
                    double y = Double.parseDouble(tokens[3]);
                    double z = Double.parseDouble(tokens[4]);

                    // Isolate element symbol
                    element = atomName.replaceAll("[^A-Za-z]", "");
                    if (element.length() > 2) {
                        element = element.substring(0, 2);
                    }

                    // --- For lone pairs and dummy atoms ---
                    String recordType = "ATOM  ";
                    String pdbElementColumn;
                    if (element.equals("LP") || element.equals("M")) {
                        recordType = "HETATM"; // Mark as not-standard-atom
                        pdbElementColumn = " X"; // Generic element "X", so that the viewer not crashes
                    } else {
                        pdbElementColumn = String.format("%2s", element); // normal element
                    }
                    
                    // PDB-format for the atom (Strict column width)
                    String pdbLine = String.format(Locale.US,
                            "ATOM  %5d %-4s MOL A   1    %8.3f%8.3f%8.3f  1.00  0.00          %2s",
                            atomId, // Use the ID from the Tinker file
                            formatPdbAtomName(atomName),
                            x,
                            y,
                            z,
                            String.format("%2s", element).toUpperCase()
                    );
                    writer.write(pdbLine);
                    writer.newLine();

                    // --- Handle CONECT ---
                    bondedAtoms = new ArrayList<>();
                    
                    for (int i = 6; i < tokens.length; i++) {
                        bondedAtoms.add(Integer.valueOf(tokens[i]));
                    }

                    // PDB allows max. 4 bonds pro CONECT-line
                    // If an atom has 5 bonds, we need two lines.
                    for (int i = 0; i < bondedAtoms.size(); i += 4) {
                        connectBuilder = new StringBuilder("CONECT");
                        
                        // The start atom (column 7-11)
                        connectBuilder.append(String
                                .format(Locale.US, "%5d", atomId));
                        
                        // Insert up to 4 atoms (column 12-16, 17-21, 22-26, 27-31)
                        for (int j = 0; j < 4 
                                && (i + j) < bondedAtoms.size(); j++) {
                            connectBuilder.append(
                                    String.format(Locale.US, "%5d", 
                                            bondedAtoms.get(i + j)));
                        }
                        
                        connectLines.add(connectBuilder.toString());
                    }

                    atomCounter++;
                }
            }
            
            // Write all saved CONECT-lines at the end of the file
            for (String conectLine : connectLines) {
                writer.write(conectLine);
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