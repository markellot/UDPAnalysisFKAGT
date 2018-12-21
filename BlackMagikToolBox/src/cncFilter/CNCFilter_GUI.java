package cncFilter;


import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.NoSuchElementException;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import general.BamChrChanger;
import general.FindBam;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

/*
 * A filter designed specifically for the CNC variants called in the MakeBamROC program.
 * It looks at both coverage across multiple individuals and by locus/position to
 * determine if the variant is a truly extreme novel exon deletion (true Call no Call).
 */
public class CNCFilter_GUI {
	// Initialize global variables
	private int Iind;
	private int Igene;
	private int chrIndex;
	private int posIndex;
	private int posendIndex;
	// private static int naIndex;

	// private static ArrayList<Integer> family;
	private static ArrayList<SamReader> familyBam;
	

	private ArrayList<ArrayList<String>> bamChrHeaders  = new ArrayList<ArrayList<String>>();
	private HashMap<String, ArrayList<String>> bamChrMap = new HashMap<String, ArrayList<String>>();
	
	private File vsFile;
	private String dest;
	private File pedFile;
	private File bamlo;
	private int clusterSize;
	private String exonBoundConfig;
	private int CNCExon;
	private File output;
	
	public CNCFilter_GUI(File vsFile, String dest, File ped, File bammy, int cluster, HashMap<String, ArrayList<String>> bamChrMap, String exonBoundConfig, int CNCExon) throws IOException {
		
		this.vsFile = vsFile;
		this.dest = dest;
		this.pedFile = ped;
		this.bamlo = bammy;
		this.clusterSize = cluster;
		this.bamChrMap = bamChrMap;
		this.exonBoundConfig = exonBoundConfig;
		this.CNCExon = CNCExon;
		
		CNCer();
	}
	
	public void CNCer() throws IOException {
		
		/// Get input VarSifter file
		//File vsFile = TXTFile.getVSFile();
		// String vsFileString = vsFile.getPath();

		/// Get input pedigree file
		//File pedFile = TXTFile.getPedigreeFile();

		/// Get BAM directory file
		//File bamlo = TXTFile.getBAMDirectoryFile();
		String bamloPath = bamlo.getPath();

		/// Allow user to choose what minimum size CNC cluster they want.
		//int clusterSize = clusterMin();
		// int clusterSize = Integer.parseInt(JOptionPane.showInputDialog(new
		// Frame("Choose cluster size minimum"),
		// "Please write how many CNC variants you want to be considered a CNC cluster
		// size minimum. \n"
		// + " Please keep your input as a positive integer. If you choose 1, then the
		// bounds \n "
		// + "of the exon will be +/- 75 bp of the CNC variant's position in the
		// genome.",
		// JOptionPane.PLAIN_MESSAGE));
		// System.out.println("clusterQuant: " + clusterSize);

		/// Prompt user to input output file name
		
		BufferedReader pedData = new BufferedReader(new FileReader(pedFile));
		String Line = pedData.readLine();

		JFrame jframe = canceler();
		jframe.setLocationRelativeTo(null);
		jframe.setVisible(true);

		/// Read through the pedigree file
		while (Line != null) {
			String[] curLine = Line.split("\t");
			// For the proband
			initializer(vsFile, dest, curLine, bamloPath, clusterSize);

			/// Leaving out the portion in analyzing the sibling since we're not doing an FB
			/// analysis
			Line = pedData.readLine();

		}
		
		pedData.close();
		jframe.dispose();
		
		
	}
	
	public File getOutput() {
		
		return output;
		
	}

	public void initializer(File vsFile, String outputFile, String[] trio, String bamDirectory, int clusterQuant)
			throws IOException {

		/// Read in inputs and outputs
		// System.out.println("outputFile: " + outputFile);
		// System.out.println("trio: " + Arrays.toString(trio));

		BufferedReader vsData = new BufferedReader(new FileReader(vsFile));
		
		output = new File(outputFile + ".vs");
		
		PrintWriter filterWriter = new PrintWriter(output);

		String Line = vsData.readLine();
		String[] curLine = Line.split("\t");
		filterWriter.println(Line); // Print headers to output

		/// Define headers
		ArrayList<String> headers = new ArrayList<String>();
		for (int i = 0; i < curLine.length; i++) {
			headers.add(curLine[i]);
		}

		Iind = headers.indexOf("Index");
		Igene = headers.indexOf("Gene_name");
		chrIndex = headers.indexOf("Chr");
		posIndex = headers.indexOf("LeftFlank");
		posendIndex = headers.indexOf("RightFlank");

		int mendInconsisIndex = headers.indexOf("MendInconsis");
		int mendFlagIndex = headers.indexOf("MendFlags");

		// Map of index-->Line information
		HashMap<String, String> ind2LineMap = new HashMap<String, String>();
		// Map of gene name --> indices with that gene name
		//HashMap<String, ArrayList<String>> gene2IndsMap = new HashMap<String, ArrayList<String>>();

		// Map of chrom --> inds
		HashMap<String, ArrayList<String>> chrom2IndsMap = new HashMap<String, ArrayList<String>>();
		/*
		ArrayList<String> sureSelectFams = new ArrayList<String>();
		sureSelectFams.add("UDP369");
		sureSelectFams.add("UDP3350");
		sureSelectFams.add("UDP984");
		sureSelectFams.add("UDP2439");
		sureSelectFams.add("UDP1418");
		sureSelectFams.add("UDP2993");
		sureSelectFams.add("UDP2992");
		sureSelectFams.add("UDP2296");
		sureSelectFams.add("UDP2600");
		sureSelectFams.add("UDP2601");
		*/
		//String gene = ""; // Gene name TO DELETE
		
		// For each variant
		Line = vsData.readLine();
		while (Line != null) {
			curLine = Line.split("\t");
			String mendFlag = curLine[mendFlagIndex]; // Identifier
			int mendInconsis = Integer.parseInt(curLine[mendInconsisIndex]);
			if (mendFlag.contains("CNC") || mendInconsis == 4) { // If variant is a
				// potential CNC
				// if (cncEval(curLine)) {
				String ind = curLine[Iind]; // index
				// String gene = curLine[Igene]; // Gene name
				//gene = curLine[Igene];
				String chrom = curLine[chrIndex]; // Chromosome
				//System.out.println("potential CNC gene: " + gene);
				
				// Put each CNC candidate in map for index-->line info
				ind2LineMap.put(ind, Line);

				/// Fill out the maps from above
				// By gene name
				/*
				if (!gene2IndsMap.containsKey(gene)) {
					ArrayList<String> indArr = new ArrayList<String>();
					indArr.add(ind);
					gene2IndsMap.put(gene, indArr);
				} else {
					ArrayList<String> indArr = gene2IndsMap.get(gene);
					indArr.add(ind);
					gene2IndsMap.put(gene, indArr);
				}
				*/
				// By position in a chromosome
				if (!chrom2IndsMap.containsKey(chrom)) {
					ArrayList<String> indArr = new ArrayList<String>();
					indArr.add(ind);
					chrom2IndsMap.put(chrom, indArr);
				} else {
					ArrayList<String> indArr = chrom2IndsMap.get(chrom);
					indArr.add(ind);
					chrom2IndsMap.put(chrom, indArr);
				}
			} else {
				filterWriter.println(Line); // Print line to output as normal if not CNC
			}
			Line = vsData.readLine();
		}
		vsData.close();
		
		
		/// Get the boundary to search with the bamCompare program by the locus
		/// Map bounds of region {chrom #, start, end} --> indices contained within
		/// Set that contains potential CNC clusters
		HashMap<String[], ArrayList<String>> boundMap = new HashMap<String[], ArrayList<String>>();

		for (String cncChrom : chrom2IndsMap.keySet()) {
			ArrayList<String> chromIndsList = chrom2IndsMap.get(cncChrom);
			// System.out.println("cncChrom: " + cncChrom);
			// System.out.println("chromIndsList: " + list2StrStr(chromIndsList));
			HashMap<String[], ArrayList<String>> chromBoundsMap = cncBoundsByPosition(cncChrom, chromIndsList,
					ind2LineMap, clusterQuant);
			// ArrayList<String[]> posMapStore = new ArrayList<String[]>();
			
			
			for (String[] chromBoundsArr : chromBoundsMap.keySet()) {
				// System.out.println("chromBoundsArr: " + Arrays.toString(chromBoundsArr));
				if (!boundMap.containsKey(chromBoundsArr) && chromBoundsArr != null) {
					boundMap.put(chromBoundsArr, chromBoundsMap.get(chromBoundsArr));
				}
			}
		}
		
		
		/// Evaluate the read depth for each potential CNC for people
		// For every potential CNC region
		HashMap<String, String> BamMap = new HashMap<String, String>();
		FindBam.InitializeBam(bamDirectory, BamMap);

		if (!boundMap.isEmpty()) {
			ArrayList<String[]> cncLines = new ArrayList<String[]>();
			for (String[] boundArr : boundMap.keySet()) {
				if (boundArr != null) {
					//System.out.println("array: " + Arrays.toString(boundArr));
					boolean covEnough = bamReadCompare(trio, boundArr, BamMap);
					// If the read depth in the family members is enough, print to output
					if (covEnough) {
						ArrayList<String> cncIndList = boundMap.get(boundArr);
						/// If family was sequenced with Sure Select Human All Exon capture kit,
						/// we need to further evaluate if the CNC is a true CNC or a product of
						/// poor capture.
						boolean probandRangeNovel = probandRangeCoverage(trio, boundArr, BamMap);
						if (probandRangeNovel) {
							for (String cncInd : cncIndList) {
								// Retrieve the appropriate variant line given the index
								String cncLine = ind2LineMap.get(cncInd);
								String[] cncLineSplit = cncLine.split("\t");
								cncLines.add(cncLineSplit);
							}
							// } else {
							// System.out.println("False CNC! Not novel in the proband range");
						}
						// } else {
						// System.out.println("False CNC! Not enough coverage");
					}
				}
			}
			
			//CNC Exon Filtering Addin which removes any CNCs that are not within an exon's bounds
			if(CNCExon == 0 && cncLines.size() > 0) {
				
				//System.out.println("\n\n");
				
				//Sorts the CNCs
				Collections.sort(cncLines, new Comparator<String[]>() {
	
					@Override
					public int compare(String[] line1, String[] line2) {
						if (!line1[chrIndex].equals(line2[chrIndex])) {
							return line1[chrIndex].compareTo(line2[chrIndex]);
						} else if (Integer.parseInt(line1[posIndex]) == Integer.parseInt(line2[posIndex])) {
							return 0;
						} else if (Integer.parseInt(line1[posIndex]) < Integer.parseInt(line2[posIndex])) {
							return -1;
						} else {
							return 1;
						}
					}
					
				});
				
				//Gets Exon Bound Config
				BufferedReader exonReder = new BufferedReader(new FileReader(new File(exonBoundConfig)));
				
				String Linez;
				int i = 0;
				while((Linez = exonReder.readLine()) != null) {
					
					if(i >= cncLines.size()) { //if done looking at all CNC's from VS file
						
						break;
						
					}
					
					//Removes the chr from the VS lines
					String vsChr = cncLines.get(i)[chrIndex].substring(3);
					
					//Retrieves each exon bound line
					String [] lineSplit = Linez.split("\t");
					
					String chr = lineSplit[0];
					int startPos = Integer.parseInt(lineSplit[1]);
					int endPos = Integer.parseInt(lineSplit[2]);
					
					if (!vsChr.equals(chr)) { //if Chr in VS file and Exon bounds are not the same
						
						if(vsChr.compareTo(chr) < 0) { //And the VS file chr is before the Exon Bound chr (aka it moved on)
							
							//System.out.println("1: Chr: " + chr + " Pos: " + endPos + " i: " + i);
							//System.out.println("1-2: Chr: " + vsChr + " Pos: " + cncLines.get(i)[posIndex] + " i: " + i);
							i++;
							while(i < cncLines.size()) { //Checks the next position in the VS list to determine if a CNC is found for the current bound but another section of list
								
								if((Integer.parseInt(cncLines.get(i)[posIndex]) < startPos)) { //if the start position bound is already past the next cnc line, go to next cnc line
									
									//System.out.println("3: Chr: " + chr + " Pos: " +  startPos + " - " + endPos + " i: " + i);
									//System.out.println("3-2: Chr: " + vsChr + " Pos: " + cncLines.get(i)[posIndex] + " i: " + i);
									i++;
									
								} else if (Integer.parseInt(cncLines.get(i)[posIndex]) >= startPos && Integer.parseInt(cncLines.get(i)[posIndex]) <= endPos) { //if it matches the cnc bounds
									
									String output = cncLines.get(i)[0];
									
									for(int j = 1; j < cncLines.get(i).length; j++) {
										
										output += "\t" + cncLines.get(i)[j];
										
									}
									
									filterWriter.println(output);
									//System.out.println("2: Chr: " + chr + " Pos: " +  startPos + " - " + endPos + " i: " + i);
									//System.out.println("2-2: Chr: " + vsChr + " Pos: " + cncLines.get(i)[posIndex] + " i: " + i);
									
									i++;
									
								} else if(Integer.parseInt(cncLines.get(i)[posIndex]) > endPos) { //If not past, can resume normal iterations from reading config lines
									
									break;
								}
								
							}
							
						}
						
					} else if (Integer.parseInt(cncLines.get(i)[posIndex]) >= startPos && Integer.parseInt(cncLines.get(i)[posIndex]) <= endPos) { //if it matches the cnc bounds
						
						String output = cncLines.get(i)[0];
						
						for(int j = 1; j < cncLines.get(i).length; j++) {
							
							output += "\t" + cncLines.get(i)[j];
							
						}
						
						filterWriter.println(output);
						i++;
						while(i < cncLines.size()) { //Checks the next position in the VS list to determine if a CNC is found for the current bound but another section of list
							
							if((Integer.parseInt(cncLines.get(i)[posIndex]) < startPos)) { //if the start position bound is already past the next cnc line, go to next cnc line
								
								//System.out.println("3: Chr: " + chr + " Pos: " +  startPos + " - " + endPos + " i: " + i);
								//System.out.println("3-2: Chr: " + vsChr + " Pos: " + cncLines.get(i)[posIndex] + " i: " + i);
								i++;
								
							} else if (Integer.parseInt(cncLines.get(i)[posIndex]) >= startPos && Integer.parseInt(cncLines.get(i)[posIndex]) <= endPos) { //if it matches the cnc bounds
								
								String output2 = cncLines.get(i)[0];
								
								for(int j = 1; j < cncLines.get(i).length; j++) {
									
									output2 += "\t" + cncLines.get(i)[j];
									
								}
								
								filterWriter.println(output2);
								//System.out.println("2: Chr: " + chr + " Pos: " +  startPos + " - " + endPos + " i: " + i);
								//System.out.println("2-2: Chr: " + vsChr + " Pos: " + cncLines.get(i)[posIndex] + " i: " + i);
								
								i++;
								
							} else if(Integer.parseInt(cncLines.get(i)[posIndex]) > endPos) { //If not past, can resume normal iterations from reading config lines
								
								break;
							}
							
						}
						//System.out.println("2: Chr: " + chr + " Pos: " +  startPos + " - " + endPos + " i: " + i);
						//System.out.println("2-2: Chr: " + vsChr + " Pos: " + cncLines.get(i)[posIndex] + " i: " + i);
						
					} else if (Integer.parseInt(cncLines.get(i)[posIndex]) < endPos) { //If the config lines has read past the cnc lines, then move on to the next position in CNC list
						
						//System.out.println("3: Chr: " + chr + " Pos: " +  startPos + " - " + endPos + " i: " + i);
						//System.out.println("3-2: Chr: " + vsChr + " Pos: " + cncLines.get(i)[posIndex] + " i: " + i);
						i++;
						while(i < cncLines.size()) { //Checks the next position in the VS list to determine if a CNC is found for the current bound but another section of list
							
							if((Integer.parseInt(cncLines.get(i)[posIndex]) < startPos)) { //if the start position bound is already past the next cnc line, go to next cnc line
								
								//System.out.println("3: Chr: " + chr + " Pos: " +  startPos + " - " + endPos + " i: " + i);
								//System.out.println("3-2: Chr: " + vsChr + " Pos: " + cncLines.get(i)[posIndex] + " i: " + i);
								i++;
								
							} else if (Integer.parseInt(cncLines.get(i)[posIndex]) >= startPos && Integer.parseInt(cncLines.get(i)[posIndex]) <= endPos) { //if it matches the cnc bounds
								
								String output = cncLines.get(i)[0];
								
								for(int j = 1; j < cncLines.get(i).length; j++) {
									
									output += "\t" + cncLines.get(i)[j];
									
								}
								
								filterWriter.println(output);
								//System.out.println("2: Chr: " + chr + " Pos: " +  startPos + " - " + endPos + " i: " + i);
								//System.out.println("2-2: Chr: " + vsChr + " Pos: " + cncLines.get(i)[posIndex] + " i: " + i);
								
								i++;
								
							} else if(Integer.parseInt(cncLines.get(i)[posIndex]) > endPos) { //If not past, can resume normal iterations from reading config lines
								
								break;
							}
							
						}
						
					}
					
				}
				exonReder.close();
				
				//System.out.println(cncLines.size());
				//System.out.println(k);
				//System.out.println(cncLines.size()-k);
				//filterWriter.println(cncLine);
				
			} else if(CNCExon == 1) { //If not filtering by Exon Boundaries (Prints all CNC's)
				
				for (String[] cncLins : cncLines) {

					String output = cncLins[0];
					
					for(int j = 1; j < cncLins.length; j++) {
						
						output += "\t" + cncLins[j];
						
					}
					
					filterWriter.println(output);
				}
				
			}
		}
		vsData.close();
		filterWriter.close();
	}

	/*
	 * Evaluate variant to see if it's a potential CNC variant. Uses the same logic
	 * to ID CNC's in the KaylaKode.
	 */
	/*
	 * public static boolean cncEval(String[] line) { int deletion = 0;
	 * 
	 * if (Integer.parseInt(line[family.get(2) + 2]) == 0 &&
	 * Integer.parseInt(line[family.get(2) + 1]) == 0) { if
	 * (Integer.parseInt(line[family.get(0) + 2]) >= 5 &&
	 * Integer.parseInt(line[family.get(1) + 2]) >= 5) {
	 * 
	 * for (int i = naIndex; i < line.length; i += 3) { if (Integer.parseInt(line[i
	 * + 2]) < 5 && Integer.parseInt(line[i + 1]) < 5) { deletion++; } }
	 * 
	 * if (deletion < 3) { return true; } } } return false; }
	 */

	/*
	 * If multiple potential CNC variants within one gene (ID'd by name), returns
	 * int array with start position of first CNC variant and end position of last
	 * CNC variant in that locus. Input is a map for index-->Line
	 */
	/*
	public String[] cncBoundsofGene(ArrayList<String> cncInds, HashMap<String, String> cncMap, int clusterQuant) {
		/// Ask user how many variants they want in a cluster to qualify for a CNC

		if (!cncInds.isEmpty() && cncQuantByGene(cncInds, clusterQuant)) { // if there are multiple potential CnC
																			// variants in a locus
			/// For each index in the locus
			ArrayList<String> cncStartArr = new ArrayList<String>();
			ArrayList<String> cncEndArr = new ArrayList<String>();
			String chrom = "";
			for (String cncIndex : cncInds) {
				String[] cncLine = cncMap.get(cncIndex).split("\t");
				chrom = cncLine[chrIndex];
				String startPos = cncLine[posIndex];
				String endPos = cncLine[posendIndex];
				cncStartArr.add(startPos);
				cncEndArr.add(endPos);
			}
			/// Sort through start and end position arrays to find first and last CNC
			/// variant
			Collections.sort(cncStartArr);
			Collections.sort(cncEndArr, Collections.reverseOrder());

			String cncStart = "";
			String cncEnd = "";
			/// Get start and end of potential deleted exon
			if (clusterQuant == 1) {
				cncStart = Integer.toString(Integer.parseInt(cncStartArr.get(0)) - 75);
				cncEnd = Integer.toString(Integer.parseInt(cncStartArr.get(0)) + 75);
			} else {
				cncStart = cncStartArr.get(0);
				cncEnd = cncEndArr.get(0);
			}
			String[] cncPosBounds = { chrom, cncStart, cncEnd };
			// System.out.println("cncPosBounds: " + Arrays.toString(cncPosBounds));

			return cncPosBounds;
		} else {
			return null;
		}
	}
*/	
	/*
	 * Returns true if there are multiple potential CNC variants within one gene
	 * (ID'd by name).
	 */
	/*
	public boolean cncQuantByGene(ArrayList<String> cncMap, int clusterMin) {
		if (cncMap.size() >= clusterMin) {
			return true;
		} else {
			return false;
		}
	}
*/
	// /*
	// * Evaluates how many short reads are in the bam files of a group of people
	// * (i.e. family members of the proband, person with next lowest coverage)
	// using
	// * the bamCompare program. If there are too few in anyone (beyond the
	// proband),
	// * then return false.
	// */
	// public boolean cncCoverageCheck(String person, int[] positions) {
	//
	// return true;
	// }

	/*
	 * Inputs an array of all the indices within a chromosome that were labeled CNC.
	 * Returns a list of indices grouped together by position.
	 */
	public HashMap<String[], ArrayList<String>> cncBoundsByPosition(String chrom, ArrayList<String> indexArr,
			HashMap<String, String> cncMap, int clusterQuant) {
		/// Sort in ascending order for the indices of the CNC variants (assuming
		/// indices are in order by position)
		if (!indexArr.isEmpty() && indexArr != null) {
			Collections.sort(indexArr);
			// System.out.println("indexArr: " + list2StrStr(indexArr));

			// Ask if the distance between variants is within 1000 bp
			HashMap<String[], ArrayList<String>> posArrMap = new HashMap<String[], ArrayList<String>>();
			ArrayList<String> posIndList = new ArrayList<String>();
			int count = 1; // Number of variants within 250 bp of each other
			int startBound = 0;
			int endBound = 0;
			
			if(indexArr.size() == 1) {
				
				clusterQuant = 1;
				
			}
			
			if (clusterQuant == 1) { /// Only 1 variant needed for cluster
				for (String ind : indexArr) {
					String[] line1 = cncMap.get(ind).split("\t");
					//System.out.println("line1: " + Arrays.toString(line1));
					String ind1 = line1[Iind];
					int pos1 = Integer.parseInt(line1[posIndex]);
					startBound = pos1 - 75;
					endBound = pos1 + 75;
					if (!posIndList.contains(ind1)) {
						posIndList.add(ind1);
					}
				}
				
				String[] posArr = { chrom, Integer.toString(startBound), Integer.toString(endBound) };
				posArrMap.put(posArr, posIndList);
				
			} else {
				
				for (int i = 0; i < indexArr.size() - 1; i++) {
					String[] line1 = cncMap.get(indexArr.get(i)).split("\t");
					String[] line2 = cncMap.get(indexArr.get(i + 1)).split("\t");

					//System.out.println("line1: " + Arrays.toString(line1));
					//System.out.println("line2: " + Arrays.toString(line2));

					String ind1 = line1[Iind];
					String ind2 = line2[Iind];
					int pos1 = Integer.parseInt(line1[posIndex]);
					int pos2 = Integer.parseInt(line2[posIndex]);
					int endPos2 = Integer.parseInt(line2[posendIndex]);

					/// Set initial positions as the two starting pieces.
					if (startBound == 0) {
						startBound = pos1 + 1;
					}

					if (endBound == 0) {
						endBound = pos2 - 1;
					}
					// System.out.println("Difference between pos2-pos1: " + (pos2 - pos1));
					/// Evaluate distance between positions
					if (Math.abs(pos2 - pos1) <= 250) {
						if (pos1 < startBound) {
							startBound = pos1;
						}
						if (endPos2 > endBound) {
							endBound = endPos2;
						}
						count++;
						if (!posIndList.contains(ind1)) {
							posIndList.add(ind1);
						}
						if (!posIndList.contains(ind2)) {
							posIndList.add(ind2);
						}
						// If at end and all variants were within 250 bp, then add to the absCount
						if (i == indexArr.size() - 2 && count >= clusterQuant) {

							String[] posArr = { chrom, Integer.toString(startBound), Integer.toString(endBound) };

							//System.out.println("chrom w clusterQuant: " + chrom + ", count: " + count + ", posArr"
							//		+ Arrays.toString(posArr) + ", posIndList: " + list2StrStr(posIndList));
							// System.out.println(
							// "posArr" + Arrays.toString(posArr) + ", posIndList: " +
							// list2StrStr(posIndList));
							posArrMap.put(posArr, posIndList);
						}
					} else {
						// If there are more than 250 bp between variants, then check if there are
						// multiple (>=clusterQuant) variants already found within 250 bp of each other.
						// If yes, then
						// add to absolute count of CNC variant groups in that chromosome. Also, set
						// count back to 1 and loop again
						if (count >= clusterQuant) {
							String[] posArr = new String[3];
							posArr[0] = chrom;
							posArr[1] = Integer.toString(startBound);
							posArr[2] = Integer.toString(endBound);
							// System.out.println(
							// "posArr" + Arrays.toString(posArr) + ", posIndList: " +
							// list2StrStr(posIndList));
							//System.out.println("chrom w clusterQuant: " + chrom + ", count: " + count + ", posArr"
							//		+ Arrays.toString(posArr) + ", posIndList: " + list2StrStr(posIndList));

							posArrMap.put(posArr, posIndList);
						}
						// Reset the counters and boundaries
						count = 1;
						startBound = 0;
						endBound = 0;
						posIndList = new ArrayList<String>();
					}
				}
			}
			return posArrMap;
		} else {
			return null;
		}
	}

	/*
	 * Compares the number of reads per person (i.e. family, proband, person with
	 * next lowest coverage)
	 */
	public boolean bamReadCompare(String[] people, String[] cncLocs, HashMap<String, String> BamMap)
			throws IOException {
		/// Determine bam file locations
		// HashMap<String, String> BamMap = new HashMap<String, String>();
		// String bamlo =
		// "U:/KaylaEmily/ForwardBackwardsPerfamily/full_path_to_bamfile_auto.txt";
		// FindBam.InitializeBam(bamlo, BamMap);

		/// Locate where the appropriate bam files are
		String[] bamLocation = FindBam.MakeBamString(people, BamMap);
		AddBam(bamLocation, bamChrMap);

		/// For each family member's bam file
		for (int i = 0; i < familyBam.size(); i++) {
			if (i != 2) { // If not the proband (b/c proband was shown as 0 all the way)
				SamReader sR = familyBam.get(i);
				
				ArrayList<Integer> readCountList = bamReadCount(cncLocs, sR, i);
				int poorSpotCount = 0;
				for (int readCount : readCountList) {
					// If there are less than 5 reads at 5+ spots within the region, then assume
					// the region overall is not adequately sequenced
					if (readCount < 5) {
						poorSpotCount++;
						if (readCountList.size() >= 10 && poorSpotCount >= 5) {
							//System.out.println("False CNC!");
							return false;
						} else if (readCountList.size() < 10 && poorSpotCount >= 4) {
							//System.out.println("False CNC!");
							return false;
						} else if (readCountList.size() < 7 && poorSpotCount >= 3) {
							//System.out.println("False CNC!");
							return false;
						} else if (readCountList.size() < 5 && poorSpotCount >= 2) {
							//System.out.println("False CNC!");
							return false;
						} else if (readCountList.size() < 3 && poorSpotCount >= 1) {
							//System.out.println("False CNC!");
							return false;
						} else if (readCountList.size() < 2 && poorSpotCount >= 1) {
							//System.out.println("False CNC!");
							return false;
						}
					}
				}
			}
		}
		
		//System.out.println("Chrom: " + cncLocs[0] + "Positions: " + cncLocs[1] + " & " + cncLocs[2]);
		
		return true;
	}

	/*
	 * Evaluates if the proband has reads in the regions immediately flanking the
	 * presumed CNC position (bounds determined by parents' pileups). If it's a true
	 * CNC, then there should be no reads in that region. Otherwise, it's a false
	 * positive.
	 */
	public boolean probandRangeCoverage(String[] people, String[] cncLocs, HashMap<String, String> BamMap)
			throws IOException {
		// HashMap<String, String> BamMap = new HashMap<String, String>();
		// String bamlo =
		// "U:/KaylaEmily/ForwardBackwardsPerfamily/full_path_to_bamfile_auto.txt";
		// FindBam.InitializeBam(bamlo, BamMap);

		/// Locate where the appropriate bam files are
		String[] bamLocation = FindBam.MakeBamString(people, BamMap);
		AddBam(bamLocation, bamChrMap);
		
		/// Extract bounds of flanking region from parents
		ArrayList<String[]> flankRegionList = parentBounds(cncLocs);
		
		if(flankRegionList.size() < 2) { //if parents do not exist at one of the flanking region sites
			
			return false;
			
		}
		
		/// NOTE: Worry here is what if there's both a true positive exon and a false
		/// positive exon in the same gene? This only evaluates the region overall
		for (String[] flankRegionBounds : flankRegionList) {
			//System.out.println("parentBounds: " + Arrays.toString(flankRegionBounds));
			/// Get number of reads at each position in the region for the proband
			SamReader probandSR = familyBam.get(2);
			ArrayList<Integer> readCountList = bamReadCount(flankRegionBounds, probandSR, 2);
			for (int readCount : readCountList) {
				if (readCount > 0) { // If any reads in proband, is a false CNC
					//System.out.println("Reads in the proband! False CNC!");
					return false;
				}
			}
		}
		return true;
	}
	
	/*
	 * Returns the number of reads within a given region of a person's bam file.
	 */
	public ArrayList<Integer> bamReadCount(String[] bamLocs, SamReader samReader, int i) throws IOException {
		SAMRecordIterator iterator;
		SAMRecord samRecord;

		ArrayList<Integer> readCountsFull = new ArrayList<Integer>();

		/// Establish bounds of the region of interest
		// System.out.println("bamLocs: " + Arrays.toString(bamLocs));
		String chrom = bamLocs[0];
		int startLoc = Integer.parseInt(bamLocs[1]);
		int endLoc = Integer.parseInt(bamLocs[2]);
		int curLoc = startLoc;
		
		//Translator - changes chromosome input depending on the type of BAM chrs the BAMs take in
		if(!bamChrHeaders.get(i).contains(chrom)) {
			
			chrom = BamChrChanger.translate(chrom);
			
		}
		
		try {
			
			while (curLoc <= endLoc) {
				iterator = samReader.queryOverlapping(chrom, curLoc, curLoc);
				// Store read names so that the same read doesn't get counted twice
				HashSet<String> readNameSet = new HashSet<String>();
				while (iterator.hasNext()) {
					samRecord = iterator.next();
					String qName = samRecord.getReadName(); // Extract read QNAME
					readNameSet.add(qName);
				}
				/// Count number of unique reads within the region
				readCountsFull.add(readNameSet.size());
				iterator.close();
				curLoc = curLoc + 10; // NOTE: This is a metric that may need changing. Done for sake of computational
										// efficiency
			}
		}catch(IllegalArgumentException j) { //Sends Error Message if translation does not work
				System.gc();
				JOptionPane.showMessageDialog(new Frame("Error"), "The Sam Reader passed an invalid Chromosome input to"
						+ " the BAM file. The chromosome sent was: " +  chrom + ". The BAM header format will be outputted as \"BAM_Header_Error.txt\"."
								+ "Please edit the \"BAM_Chr_Config\" accordingly. The Program will now exit");
				
				BamChrChanger.writeHeaderErrors(bamChrHeaders.get(i));
				
				System.exit(0);
				
				
				
		}
		return readCountsFull;
	}
	
	/*
	 * Retrieve the BAM file one at a time for the whole family
	 */
	public void AddBam(String[] bamlocation, HashMap<String, ArrayList<String>> bamChrMap) throws IOException {
		// store all the bam files in an arraylist

		familyBam = new ArrayList<SamReader>();
		for (int i = 0; i < bamlocation.length; i++) {
			File bamFile = new File(bamlocation[i]);
			// create a samreader factory
			SamReaderFactory srf = SamReaderFactory.make();
			srf.validationStringency(ValidationStringency.LENIENT);
			SamReader samR = srf.open(bamFile);
			familyBam.add(samR);
			
			//retrieves BAM chr list for each person
			bamChrHeaders.add(bamChrMap.get(bamlocation[i]));
		}
	}

	/*
	 * Extract bounds of parent's pileup relative to proband's at presumed CNC
	 * position. (NOTE: This isn't the bounds of the full pileup for a nearby exon.
	 * This takes only the reads that overlap with the CNC position of interest
	 * (POI) because it is computationally less taxing and it's unnecessary to go to
	 * the full bounds of the exon.)
	 */
	public ArrayList<String[]> parentBounds(String[] bamLocs) throws IOException {
		SAMRecordIterator iterator;
		SAMRecord samRecord;
		ArrayList<String[]> parentBoundList = new ArrayList<String[]>();

		String chrom = bamLocs[0];
		/// Initialize boundaries
		String startBound = bamLocs[1];
		String endBound = bamLocs[2];
		
		for (int i = 1; i < bamLocs.length; i++) { // For each ID'd CNC position (start + end of a region)
			int bamLocInt = Integer.parseInt(bamLocs[i]);
			ArrayList<Integer> startList = new ArrayList<Integer>();
			ArrayList<Integer> endList = new ArrayList<Integer>();
			for (int j = 0; j < 2; j++) { // For each parent
				SamReader sR = familyBam.get(j);
				
				//Translator - changes chromosome input depending on the type of BAM chrs the BAMs take in
				if(!bamChrHeaders.get(j).contains(chrom)) {
					
					chrom = BamChrChanger.translate(chrom);
					
				}
				
				try {
						
					
					/// Extract reads that overlap with CNC position of interest
					iterator = sR.queryOverlapping(chrom, bamLocInt, bamLocInt);
	
					while(iterator.hasNext()) {
						samRecord = iterator.next();
						startList.add(samRecord.getAlignmentStart());
						endList.add(samRecord.getAlignmentEnd());
						
						//System.out.println("Alignment Start: " + samRecord.getAlignmentStart());
						//System.out.println("Alignment End: " + samRecord.getAlignmentEnd());
						//System.out.println(iterator.hasNext());
					}
					iterator.close();
				}catch(IllegalArgumentException k) { //Sends Error Message if translation does not work
					System.gc();
					JOptionPane.showMessageDialog(new Frame("Error"), "The Sam Reader passed an invalid Chromosome input to"
							+ " the BAM file. The chromosome sent was: " +  chrom + ". The BAM header format will be outputted as \"BAM_Header_Error.txt\"."
									+ "Please edit the \"BAM_Chr_Config\" accordingly. The Program will now exit");
					
					BamChrChanger.writeHeaderErrors(bamChrHeaders.get(j));
					
					System.exit(0);
				}	
				
				
			}
				
			// System.out.println("startBoundList: " + list2StrInt(startList));
			// System.out.println("endBoundList: " + list2StrInt(endList));
			/// Get min of startList and max of endList
			//System.out.println("StartList size: " + startList.size());
			//System.out.println("EndList size: " + endList.size());
			
			//System.out.println("\n\n\n");
			
			//The parents may not always exist at each of the flanking regions of the CNC, so this is put into place
			if(startList.size() > 0 && endList.size() > 0) {
				startBound = Integer.toString(Collections.min(startList));
				endBound = Integer.toString(Collections.max(endList));
				String[] parentBoundArr = { chrom, startBound, endBound };
				parentBoundList.add(parentBoundArr);
			}
			// System.out.println("Original parentBounds: " +
			// Arrays.toString(parentBoundArr));
		}
		return parentBoundList;
	}
	
	public String list2StrIntArr(ArrayList<Integer[]> list) {
		StringBuilder sb = new StringBuilder();
		if (!list.isEmpty()) {
			for (Integer[] s : list) {
				if (s != null) {
					String sStr = Arrays.toString(s);
					sb.append(sStr);
					sb.append(";");
				}
			}
		}
		return sb.toString();
	}

	public String list2StrInt(ArrayList<Integer> list) {
		StringBuilder sb = new StringBuilder();
		if (!list.isEmpty()) {
			for (Integer s : list) {
				if (s != null) {
					String sStr = Integer.toString(s);
					sb.append(sStr);
					sb.append(";");
				}
			}
		}
		return sb.toString();
	}

	public String list2StrStr(ArrayList<String> list) {
		StringBuilder sb = new StringBuilder();
		if (!list.isEmpty()) {
			for (String s : list) {
				if (s != null && !s.isEmpty()) {
					sb.append(s);
					sb.append(";");
				}
			}
		}
		return sb.toString();
	}

	public String arr2TabStr(String[] arr) {
		StringBuilder sb = new StringBuilder();
		if (arr != null) {
			for (String s : arr) {
				if (s != null) {
					sb.append(s);
					sb.append("\t");
				}
			}
		}
		return sb.toString();
	}
	
	
	public JFrame canceler() {
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("FB CNC Filter");
		jframe.setSize(500, 100);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);

		JLabel jtext = new JLabel("Generating the files with CNC's filter. Use the cancel botton below to abort.",
				SwingConstants.CENTER);
		gbc.gridx = 0;
		gbc.gridy = 0;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(5, 5, 5, 5);
		gbc.weightx = 1.0;
		jpanel.add(jtext, gbc);

		JButton abort = new JButton("Cancel");
		abort.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				jframe.dispose();

				// filterWriter2.close();
				System.gc();
				System.exit(0);
			}
		});

		gbc.gridx = 0;
		gbc.gridy = 21;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(5, 5, 5, 5);
		gbc.weightx = 1.0;
		jpanel.add(abort, gbc);

		jframe.getContentPane().add(jpanel);
		return jframe;
	}
}
