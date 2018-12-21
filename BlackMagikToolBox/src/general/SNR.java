package general;
import java.awt.Frame;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import javax.swing.JOptionPane;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class SNR {
	/*
	 * Determines mismatch density in the region +/- 20bp of the variant position
	 */
	public static double variantdensity(String chrom, int pos, ArrayList<SamReader> samReader,
			ArrayList<HashMap<String, Integer>> genotypemap, ArrayList<HashMap<String, Integer>> goodgenotypemap,
			int length, ArrayList<ArrayList<String>> bamChrHeaders) throws IOException {
		SAMRecordIterator iterator;
		SAMRecord samrecord;

		/// Establish boundaries of query
		int posprior = pos - 20;
		int posposterior = pos + 20;
		
		/// Count number of mismatches at each location
		HashMap<Integer, Integer> locations = new HashMap<Integer, Integer>(); // position in genome--># of indels/SNPs
																				// called at that position
		ArrayList<Integer> linearMismatch = new ArrayList<Integer>(); // list of mismatch positions in the pileup within
																		// +/-20bp of the variant of interest
		ArrayList<Integer> mismatchPerRead = new ArrayList<Integer>(); // list of # mismatches for each read in the
																		// pileup
		int BaseMaped = 0; /// Track total # bases read in the pileup

		for (int i = 0; i < samReader.size(); i++) {
			
			//Translator - changes chromosome input depending on the type of BAM chrs the BAMs take in
			if(!bamChrHeaders.get(i).contains(chrom)) {
				
				chrom = BamChrChanger.translate(chrom);
				
			}
			
			try {
				iterator = samReader.get(i).queryOverlapping(chrom, posprior, posposterior);
				while (iterator.hasNext()) {
					samrecord = iterator.next();
	
					/// Get read metadata
					String MDtag = (String) samrecord.getAttribute("MD");
					String cigar = samrecord.getCigar().toString();
					int start = samrecord.getAlignmentStart();
					int end = samrecord.getAlignmentEnd();
	
					/// If cigar isn't empty, get the amount of "noise" and count variants within
					/// the query region
					if (!cigar.equals("*")) {
						mismatchPerRead.add(mismatchposition(MDtag, start, cigar, locations));
	
						if (pos >= start && pos <= end) {
							String read = samrecord.getReadString();
							int mapquality = samrecord.getMappingQuality();
							byte[] bases = samrecord.getBaseQualities();
	
							CountVThere(genotypemap.get(i), goodgenotypemap.get(i), cigar, read, pos, start, length,
									mapquality, bases);
						}
						
						/// Count the number of bases mapped under the region
						BaseMaped += (end - start + 1);
					}
				}
				iterator.close();
			}  catch(IllegalArgumentException j) { //Sends Error Message if translation does not work
				
				System.gc();
				JOptionPane.showMessageDialog(new Frame("Error"), "The Sam Reader passed an invalid Chromosome input to"
						+ " the BAM file. The chromosome sent was: " +  chrom + ". The BAM header format will be outputted as \"BAM_Header_Error.txt\"."
								+ "Please edit the \"BAM_Chr_Config\" accordingly. The Program will now exit");
				
				BamChrChanger.writeHeaderErrors(bamChrHeaders.get(i));
				
				System.exit(0);
				
				
				
			}
		}

		int Mismatch = 0; /// initialize # mismatches in the pileup

		for (Integer key : locations.keySet()) {
			Mismatch += locations.get(key);
			if (key >= posprior && key <= posposterior) {
				/// Group together mismatches within +/- 20bp of the variant position
				linearMismatch.add(locations.get(key));
			}
		}

		/// Calculate SNR
		double SNR = -7;

		if (BaseMaped >= 1) {
			/// Scale down to 300*50 square
			double scaling = (double) 15000 / (double) BaseMaped;

			/// General mismatch rate of the region
			double mismatchrate = (double) Mismatch / (double) BaseMaped;

			/// Calculate a metric ("covariance") for distribution of mismatches in the
			/// region
			double verticov = Histo(linearMismatch, scaling);

			/// Handles when cov is 0 (no mismatch or no read mapped there)
			/// if no read mapped or low read map (<6) high error rate
			/// if no mismatch, high depth, low error rate
			if (verticov > 0) {
				if (verticov < 1) {
					verticov = 1;
				}
				SNR = 100 * Math.log10(1 / mismatchrate) / verticov;
			}
		}
		return SNR;
	}

	/*
	 * Returns heuristic covariance value (measure of spatial distribution of
	 * mismatches in the pileup)
	 */
	public static double Histo(ArrayList<Integer> mismatch, double scaling) {

		Collections.sort(mismatch);
		Collections.reverse(mismatch);

		double cov = 0;
		int loc = 1;
		for (int i = 1; i < mismatch.size(); i += 2) {
			cov += Math.pow(loc, 2) * (mismatch.get(i)) * scaling;
			if ((i + 1) < mismatch.size()) {
				cov += Math.pow(loc, 2) * (mismatch.get(i + 1)) * scaling;
			}
			loc++;
		}

		if (mismatch.size() > 1) {
			cov = Math.sqrt(cov / ((double) mismatch.size() - 1));
		}
		return cov;
	}

	/*
	 * Map unique read sequences to genotypemap and goodgenotypemaps to the number
	 * of times each of those sequences appears in the pileup. This subroutine
	 * parses through the CIGAR string for each read. "Good quality" reads are
	 * defined as those with (map quality)>30 and (base quality for the variant of
	 * interest)>20
	 */
	public static void CountVThere(HashMap<String, Integer> genotypemap, HashMap<String, Integer> goodgenotypemap,
			String cigar, String read, int pos, int start, int length, int mapquality, byte[] bases) {

		// quality check:
		if (read.length() != bases.length) {
			//System.out.println("Base quality length does not match, something is kinda wrong");
		}
		// Cut the soft and hard clipping
		String cigarnumber = cigar.replaceAll("\\D+", ";");
		cigarnumber = cigarnumber.substring(0, cigarnumber.length() - 1);
		String[] cigarN = cigarnumber.split(";");

		// Replace all the number in a cigar string as ; and get rid of the heading ;
		String cigarletter = cigar.replaceAll("[^A-Za-z]+", ";").substring(1);
		String[] cigarL = cigarletter.split(";");

		/// Extract all the numbers
		int whereread = pos - start; /// Tracks position in the read string
		int wheregenome = start - 1; /// Tracks position in the reference genome
		String reads = ""; /// Initializes full read string sequence
		int basequality = 35;

		/// Parse through CIGAR string
		/// Track genomic position and what the base quality is at the current position.
		/// Also start building the sequence of the read.
		for (int i = 0; i < cigarL.length; i++) {

			/// Get rid of S and H clipping
			if (cigarL[i].equals("S")) {
				whereread += Integer.parseInt(cigarN[i]);
			} else if (cigarL[i].equals("H")) {

			} else if (cigarL[i].equals("M")) { /// Match (not indel or clipped)

				wheregenome += Integer.parseInt(cigarN[i]);
				
				if (pos <= wheregenome) {
					reads = Character.toString(read.charAt(whereread));
					
					if ((int) bases[whereread] < basequality) {
						basequality = (int) bases[whereread];
					}

					if ((i + 1) < cigarL.length && pos == wheregenome) {
						/// If next segment of the read being evaluated is an indel
						if (cigarL[i + 1].equals("I") && length > 0) {
							/// Insertion: read and move through read string like usual
							for (int ii = 0; ii < Integer.parseInt(cigarN[i + 1]); ii++) {
								reads += Character.toString(read.charAt(whereread + ii + 1));

								if ((int) bases[whereread + ii + 1] < basequality) {
									basequality = (int) bases[whereread + ii + 1];
								}
							}
						} else if (cigarL[i + 1].equals("D") && length > 0) {
							/// Deletion: add "D" characters for however long deleted sequence is
							/// and reduce supposed length of variant of interest (if it's an indel)
							for (int ii = 0; ii < Integer.parseInt(cigarN[i + 1]); ii++) {
								reads += "D";
								length--;
							}
						}
					}

					/// Parse through the rest of the segment of the read until the variant of
					/// interest (length of the reference minus 1)
					for (int ii = 0; ii < length; ii++) {
						whereread += 1;
						if (whereread < read.length()) {
							reads += Character.toString(read.charAt(whereread));
							if ((int) bases[whereread] < basequality) {
								basequality = (int) bases[whereread];
							}
						}
					}

					/// Put read sequence and # reads with that sequence in genotypemap
					if (!genotypemap.containsKey(reads)) {
						genotypemap.put(reads, 1);
					} else {
						genotypemap.put(reads, genotypemap.get(reads) + 1);
					}
					
					/// Put good quality read sequence and # reads with that sequence in
					/// goodgenotypemap
					if (basequality > 20 && mapquality > 30) {
						if (!goodgenotypemap.containsKey(reads)) {
							goodgenotypemap.put(reads, 1);
						} else {
							goodgenotypemap.put(reads, goodgenotypemap.get(reads) + 1);
						}
					}
					break;
				}
			} else if (cigarL[i].equals("D")) { // Deletion
				wheregenome += Integer.parseInt(cigarN[i]);
				whereread -= Integer.parseInt(cigarN[i]);
				if (pos <= wheregenome) {
					break;
				}
			} else if (cigarL[i].equals("I")) { // Insertion
				whereread += Integer.parseInt(cigarN[i]);
			}
		}
	}

	public static String replaceAll(String MDtag, String x) {
		StringBuilder MD = new StringBuilder(MDtag);
		int index = MD.indexOf(x);
		while (index != -1) {
			MD.setCharAt(index, 'K');
			index = MD.indexOf(x);
		}
		return MD.toString();

	}

	/*
	 * Counts mismatch locations and deletion locations. Deletion and insertion
	 * position is where it starts exactly with n conservation of the right most
	 * reference. Returns measurement of "noise," or number of indels and SNPs
	 * called in the read being evaluated.
	 */
	public static Integer mismatchposition(String MDtag, int start, String cigar, HashMap<Integer, Integer> locations) {
		int noise = 0;

		// Separate letters and numbers in MD tag
		MDtag = replaceAll(MDtag, "^");
		String[] MDnumber = MDtag.replaceAll("\\D+", ";").split(";");
		String[] MDletter = MDtag.replaceAll("[^A-Za-z]+", ";").split(";");

		int trackmismatch = start; // Tracks position in genome being evaluated

		/// Read through MD tag
		for (int i = 1; i < MDletter.length; i++) {
			if (Integer.parseInt(MDnumber[i - 1]) != 0) { /// If there are indels/SNPs present
				trackmismatch += Integer.parseInt(MDnumber[i - 1]);

				if (MDletter[i].length() == 2) {
					noise++;
				} else {
					noise++;
				}

				/// Add to locations map: position in genome of indel/SNP--># of indels/SNPs
				/// called at that position
				if (locations.containsKey(trackmismatch)) {
					locations.put(trackmismatch, locations.get(trackmismatch) + 1);
				} else {
					locations.put(trackmismatch, 1);
				}

				trackmismatch += 1;

			} else {
				if (MDletter[i].length() != 2) {
					noise++;
					/// Add to locations map: position in genome of indel/SNP--># of indels/SNPs in
					/// read
					if (locations.containsKey(trackmismatch)) {
						locations.put(trackmismatch, locations.get(trackmismatch) + 1);
					} else {
						locations.put(trackmismatch, 1);
					}
				}
				trackmismatch += 1;
			}
		}

		/// Eliminate the soft and hard clipping using the cigar string
		// Replace the alphabet with ";" so only numbers remain
		String cigarnumber = cigar.replaceAll("\\D+", ";");
		cigarnumber = cigarnumber.substring(0, cigarnumber.length() - 1);
		String[] cigarN = cigarnumber.split(";");

		// Replace all the number in a cigar string as ";" and get rid of the heading
		// ";"
		String cigarletter = cigar.replaceAll("[^A-Za-z]+", ";").substring(1);
		// Extract all the alphabets
		String[] cigarL = cigarletter.split(";");

		int trackmismatch2 = start - 1;

		/// Read through the cigar string. This is necessary because MD tags cannot
		/// contain information about insertions. For more information about this,
		/// please see
		/// https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files
		for (int i = 0; i < cigarL.length; i++) {
			if (cigarL[i].equals("S") || cigarL[i].equals("H")) { // Hard and soft clipping

			} else {
				if (cigarL[i].equals("I")) { // Insertion
					noise++;
					if (locations.containsKey(trackmismatch2 + 1)) {
						locations.put(trackmismatch2 + 1, locations.get(trackmismatch2 + 1) + 1);
					} else {
						locations.put(trackmismatch2 + 1, 1);
					}
				} else {
					trackmismatch2 += Integer.parseInt(cigarN[i]);
				}
			}
		}
		return noise;
	}

	/*
	 * Returns an array of two arrays with number of ref and alt reads: {(ref count,
	 * alt count), (good ref count, good alt count)}
	 */
	public static String[] callgenotype(HashMap<String, Integer> genotypemap, HashMap<String, Integer> goodgenotypemap,
			String ref, String var) {
		int refcount = 0;
		int altcount = 0;

		if (ref.length() == 1 && var.length() == 1) { // SNP
			// If this is a SNP place, some of the indel positions might contain the
			// reference or variant read
			for (String key : genotypemap.keySet()) {
				if (Character.toString(key.charAt(0)).equals(ref)) {
					refcount += genotypemap.get(key);
				} else if (Character.toString(key.charAt(0)).equals(var)) {
					altcount += genotypemap.get(key);
				}
			}

		} else {
			// If it's an indel, just get the count
			for (String key : genotypemap.keySet()) {
				// boolean isDeletion = false;
				if (key.contains("D") && (ref.length() > var.length())) { // Confirms deletion
					// isDeletion = true;
					String oldKey = key;
					key = key.replace("D", "");
					if (key.equals(ref)) {
						refcount = genotypemap.get(ref);
					}
					if (key.equals(var)) {
						altcount = genotypemap.get(oldKey);
					}
				} else if (ref.length() < var.length()) { // Insertion
					String oldKey = key;
					if (key.equals(ref)) {
						refcount = genotypemap.get(ref);
					}
					if (key.equals(var)) {
						altcount = genotypemap.get(oldKey);
					}
				}
			}
		}

		int goodrefcount = 0;
		int goodaltcount = 0;

		if (ref.length() == 1 && var.length() == 1) {
			/// If this is a SNP, the call at the position of interest might contain the
			/// reference or variant read
			for (String key : goodgenotypemap.keySet()) {
				if (Character.toString(key.charAt(0)).equals(ref)) {
					goodrefcount += goodgenotypemap.get(key);
				} else if (Character.toString(key.charAt(0)).equals(var)) {
					goodaltcount += goodgenotypemap.get(key);
				}
			}

		}

		/// If it's an indel, just get the count
		else {
			if (goodgenotypemap.containsKey(ref)) {
				goodrefcount = goodgenotypemap.get(ref);
			}
			if (goodgenotypemap.containsKey(var)) {
				goodaltcount = goodgenotypemap.get(var);
			}
		}
		return new String[] { refcount + ";" + altcount, goodrefcount + ";" + goodaltcount };
	}
}
