package confettiFilter;



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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Set;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import general.BamChrChanger;
import general.FindBam;
import general.SNR;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class PPDNFilter {
	private ArrayList<SamReader> familyBam;
	private int minorfield; /// Index of field denoting if reference is minor
	private int homVarIndex;
	private int homRefIndex;
	private int genotypeIndex;
	private int genoThreshHold;
	private int baseFilter;
	private int baseQualThresh;
	private double badReadRatioMax;	
	
	private File vsFile;
	private String destination;
	private File pedFile;
	private File bamlo;
	private ArrayList<ArrayList<String>> bamChrHeaders = new ArrayList<ArrayList<String>>();
	private HashMap<String, ArrayList<String>> bamChrMap;
	
	private File output;
	//private PrintWriter errorWriter;
	
	public PPDNFilter(File vsFile, String dest, File ped, File bammy, int minor, int homRef, int homVar, int geno, HashMap<String, ArrayList<String>> bamChrMap, int genoThresh, int baseFilter, int baseQualThresh, double badReadRatio) throws IOException {
		
		this.vsFile = vsFile;
		this.destination = dest;
		this.pedFile = ped;
		this.bamlo = bammy;
		this.minorfield = minor;
		this.homRefIndex = homRef;
		this.homVarIndex = homVar;
		this.genotypeIndex = geno;
		this.bamChrMap = bamChrMap;
		this.genoThreshHold = genoThresh;
		this.baseFilter = baseFilter;
		this.baseQualThresh = baseQualThresh;
		this.badReadRatioMax = badReadRatio;
		
		confettiParty();
		
	}
	
	
	public void confettiParty() throws IOException {
		
		/// Get input VarSifter file
		//File vsFile = TXTFile.getVSFile();
		
		/// Get input pedigree file
		//File pedFile = TXTFile.getPedigreeFile();
		
		/// Get BAM directory file
		//File bamlo = TXTFile.getBAMDirectoryFile();
		String bamloPath = bamlo.getPath();

		HashMap<String, String> BamMap = new HashMap<String, String>();
		FindBam.InitializeBam(bamloPath, BamMap);
		
		BufferedReader pedData = new BufferedReader(new FileReader(pedFile));
		String Line = pedData.readLine();

		JFrame jframe = canceler(pedData);
		jframe.setLocationRelativeTo(null);
		jframe.setVisible(true);

		/// Read through the pedigree file
		while (Line != null) {
			String[] curLine = Line.split("\t");
			// String proband = curLine[2];
			String[] bamlocation = FindBam.MakeBamString(curLine, BamMap);
			initializer(vsFile, destination, curLine, bamlocation, pedData);

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
	
	
	public void initializer(File vsFile, String outputFile, String[] trio, String[] bamlocation,
			BufferedReader pedigree) throws IOException {
		// File vsFile = TXTFile.getVSFilefromstring2(VarSifterFile, trio[2]);
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
		
		
		int Iind = headers.indexOf("Index");
		int Igene = headers.indexOf("Gene_name");
		int chrIndex = headers.indexOf("Chr");
		int posIndex = headers.indexOf("LeftFlank");
		int posendIndex = headers.indexOf("RightFlank");
		int refIndex = headers.indexOf("ref_allele");
		int varIndex = headers.indexOf("var_allele");
		// int refalleleminor = headers.indexOf("All_Gahlref_is_minor");
		// int refalleleminor = headers.indexOf("Gahl_UDPNref_is_minor");
		int mendFlagIndex = headers.indexOf("MendFlags");
		int miscFlagIndex = headers.indexOf("MiscFlags");

		int vcfPosIndex = headers.indexOf("VCF_Position");
		int vcfRefIndex = headers.indexOf("VCF_Ref");
		int vcfVarIndex = headers.indexOf("VCF_Var");

		// int allGahlHomVarIndex = headers.indexOf("All_Gahlc_homvar");
		// int allGahlHomRefIndex = headers.indexOf("All_Gahlc_homref");

		// int allGahlHomVarIndex = headers.indexOf("Gahl_UDPNc_homvar");
		// int allGahlHomRefIndex = headers.indexOf("Gahl_UDPNc_homref");
		//
		// int csHomVarIndex = headers.indexOf("hg19_CS_938c_homvar");
		//
		// int allGahlGenotypeCovIndex = headers.indexOf("Gahl_UDPNc_genotypes");
		// int clinSeqGenotypeCovIndex = headers.indexOf("hg19_CS_938c_genotypes");

		AddBam(bamlocation, bamChrMap);
		
		
		//errorWriter = new PrintWriter(new File(vsFile.getParent() + "\\Error_Info.txt"));
		
		Line = vsData.readLine();
		while (Line != null) {
			curLine = Line.split("\t");
			String mendFlag = curLine[mendFlagIndex];
			String miscFlag = curLine[miscFlagIndex];

			/// **TEMPORARILY CHANGED***
			if (miscFlag.contains("UD") || mendFlag.contains("DN")) { // If variant is an unsalvaged DN
				String ind = curLine[Iind];
				String gene = curLine[Igene];
				String chrom = curLine[chrIndex];
				String leftflank = curLine[posIndex];
				String rightflank = curLine[posendIndex];
				String refprint = curLine[refIndex];
				String altprint = curLine[varIndex];
				int refminor = Integer.parseInt(curLine[minorfield]);
				int homVarCount = Integer.parseInt(curLine[homVarIndex]);
				int homRefCount = Integer.parseInt(curLine[homRefIndex]);
				// int csHomVarCount = Integer.parseInt(curLine[csHomVarIndex]);
				int genoCount = Integer.parseInt(curLine[genotypeIndex]);
				// int csGenoCount = Integer.parseInt(curLine[clinSeqGenotypeCovIndex]);

				// int vcfPos = Integer.parseInt(curLine[vcfPosIndex]);
				// String vcfRef = curLine[vcfRefIndex];
				// String vcfVar = curLine[vcfVarIndex];

				int place = Integer.parseInt(leftflank) + 1;
				int endplace = Integer.parseInt(rightflank) - 1;
				int length = endplace - place;

				// System.out.println("DN gene: " + gene);

				/// General filter tests
				boolean readDepthCheck = readDepthEval(chrom, place, familyBam, length);
				boolean extremeNovelCheck = true;
				if (refminor == 1) {
					extremeNovelCheck = extremePopNovel(homRefCount);
				} else {
					extremeNovelCheck = extremePopNovel(homVarCount);
				}
				boolean genotypeCovCheck = genotypeCoverageEval(genoCount);
				
				/// Confetti filter tests
				boolean horizontalCheck = horizontalPPDN(chrom, place, familyBam, length, trio);
				//errorWriter.println("\n");
				boolean verticalCheck = verticalPPDN(chrom, place, familyBam, length, trio, refprint, altprint);
				//errorWriter.println("\n");
				boolean haploCheck = haploCountEngine(chrom, place, familyBam, length, trio, refprint);
				//errorWriter.println("\n\n");
				
				if (readDepthCheck && extremeNovelCheck && genotypeCovCheck && horizontalCheck && verticalCheck
						&& haploCheck) { // Passed the PP DN filters
					filterWriter.println(Line);
				} else {
					//System.out.println("Removed a gene: " + gene);
					if (!haploCheck) {
						//System.out.println("failed haplo");
						//errorWriter.println("Failed Haplo");
					}
					if (!verticalCheck) {
						//System.out.println("failed vertical");
						//errorWriter.println("Failed Vertical");
					}
					if (!horizontalCheck) {
						//System.out.println("failed horizontal");
						//errorWriter.println("Failed Horisontal");
					}
					if (!readDepthCheck) {
						//System.out.println("failed read depth");
						//errorWriter.println("Failed Read Depth");
					}
					if (!extremeNovelCheck) {
						//System.out.println("failed extreme novel");
						//errorWriter.println("Failed Extreme Novel");
					}
					if (!genotypeCovCheck) {
						//System.out.println("failed geno check");
						//errorWriter.println("Failed Genotype Check");
					}
				}
			} else {
				filterWriter.println(Line); // Print line to output file as normal
			}
			Line = vsData.readLine();
		}
		//errorWriter.close();
		vsData.close();
		filterWriter.close();
	}

	/*
	 * Assesses the quality of each read in the bam file within the position
	 * indicated by the variant. Adapted from Kayla Gu's SNR.java code. This is
	 * called "horizontal" because it assesses each read in a "horizontal" direction
	 * in the bam file. In other words, it detects the bits of "confetti" in a set
	 * of reads in the horizontal direction for each read. If there are too many
	 * bits of confetti in the read, the read is determined to be bad (a poor read).
	 * If there are enough bad reads in the pileup (at least 25% of the pileup),
	 * then the entire pileup is determined to be bad quality (False). Otherwise,
	 * it's a good quality pileup (True), horizontally speaking.
	 */
	public  boolean horizontalPPDN(String chrom, int pos, ArrayList<SamReader> samReader, int length,
			String[] family) throws IOException {
		SAMRecordIterator iterator;
		SAMRecord samRecord;

		/// Positions of interest include +/-10 of variant of interest
		int startBound = pos;
		int endBound = pos + length;

		/// Iterate through the samReader
		for (int i = 0; i < samReader.size(); i++) {
			
			//Translator - changes chromosome input depending on the type of BAM chrs the BAMs take in
			if(!bamChrHeaders.get(i).contains(chrom)) {
				
				chrom = BamChrChanger.translate(chrom);
				
			}
			//errorWriter.println("Hor PPDN");
			try {
				iterator = samReader.get(i).queryOverlapping(chrom, startBound, endBound);
				// System.out.println("family member: " + family[i]);
				int readCount = 0;
				int poorReadCount = 0;
				int l = 0;
				while (iterator.hasNext()) {
					samRecord = iterator.next();
					String cigar = samRecord.getCigar().toString();
					int alignStart = samRecord.getAlignmentStart();
					int alignEnd = samRecord.getAlignmentEnd();
					if (!cigar.equals("*")) {
						if (pos >= alignStart && pos <= alignEnd) {
							int mapQuality = samRecord.getMappingQuality();
							byte[] bases = samRecord.getBaseQualities();
							int poorBaseCount = 0;
							for (byte base : bases) {
								if ((int) base < baseFilter) { // poor base quality score if less than input cutoff 
									poorBaseCount++;
								}
							} //# of poor bases is > the input threshold, map quality is < than 30 or number of reads is less than 10 adds one to bad read count
							if (poorBaseCount > baseQualThresh || mapQuality < 30 || bases.length < 10) { // heuristic
								poorReadCount++;
								if(i == 1) {
									//errorWriter.println("Read #: " + l + ": PBC: " + poorBaseCount + " MQ: " + mapQuality + " Base Length: " + bases.length);
								}
							}
						}
						readCount++;
						l++;
					}
				}
				iterator.close();
				double badReadRatio = (double) poorReadCount / readCount;
				
				//errorWriter.println(badReadRatio);
				//errorWriter.println(poorReadCount);
				//errorWriter.println(readCount);
				
				//If one person's ratio is > the input bad read ratio cutoff, the line is not printed since it fails
				if (badReadRatio >= badReadRatioMax) {
					return false; // This set of reads fails the horizontal base quality check.
				}
			}  catch(IllegalArgumentException j) { //Sends Error Message if translation does not work
				System.gc();
				JOptionPane.showMessageDialog(new Frame("Error"), "The Sam Reader passed an invalid Chromosome input to"
						+ " the BAM file. The chromosome sent was: " +  chrom + ". The BAM header format will be outputted as \"BAM_Header_Error.txt\"."
								+ "Please edit the \"BAM_Chr_Config\" accordingly. The Program will now exit");
				
				BamChrChanger.writeHeaderErrors(bamChrHeaders.get(i));
				
				System.exit(0);
				
				
				
			}
		}
		return true;
	}

	/*
	 * In contrast to the horizontal filter above (horizontalPPDNFilter), the
	 * vertical filter looks at the number of variants called at the variant
	 * position of interest across the family's pileups. If the parents' and
	 * siblings' variant genotype calls are together greater than the proband's,
	 * then it is likely that someone was mis-called and/or that the quality of the
	 * reads was relatively poor. The same goes for if EVERYONE has greater than 3
	 * variant calls in their pileup. Given that this filter is meant for a supposed
	 * de novo, variant calls in every family member (including the apparent
	 * unaffected sib and the parents) throws the likelihood of a true positive DN
	 * call into doubt.
	 */
	public boolean verticalPPDN(String chrom, int pos, ArrayList<SamReader> familySamReader, int length,
			String[] family, String ref, String var) throws IOException {

		// int startBound = pos - 1;
		// int endBound = pos + 1;

		/// Make genotype maps for the family. One for genotype calls in general and
		/// another for good calls.
		ArrayList<HashMap<String, Integer>> genotypemap = new ArrayList<HashMap<String, Integer>>();
		for (int i = 0; i < familySamReader.size(); i++) {
			genotypemap.add(new HashMap<String, Integer>());
		}
		ArrayList<HashMap<String, Integer>> goodgenotypemap = new ArrayList<HashMap<String, Integer>>();
		for (int i = 0; i < familySamReader.size(); i++) {
			goodgenotypemap.add(new HashMap<String, Integer>());
		}
		/// Get array of alt reads for each family member
		/// We're using SNR.variantdensity here because the genotype mapping and calling
		/// was already written by Kayla in her SNR program. No sense in reinventing the
		/// wheel.
		double Snrscore = SNR.variantdensity(chrom, pos, familyBam, genotypemap, goodgenotypemap, length, bamChrHeaders);

		/// Iterate through each family member to get good alt counts
		ArrayList<String[]> arrayReadCount = new ArrayList<String[]>();
		ArrayList<Integer> arrayAltCount = new ArrayList<Integer>();
		int probVarCount = 0;
		int sibParVarCount = 0;

		int altMultiCount = 0;
		int superAltMultiCount = 0;
		for (int i = 0; i < family.length; i++) {

			String[] readCount = SNR.callgenotype(genotypemap.get(i), goodgenotypemap.get(i), ref, var);
			int varCount = Integer.parseInt(readCount[0].split(";")[1]);
			if (i == 2) {
				probVarCount = varCount;
			} else {
				sibParVarCount += varCount;
			}
			arrayReadCount.add(readCount);
			arrayAltCount.add(varCount);
			/// If everyone has more than 3 good alt reads, something's wrong.
			// This isn't the most elegant way to do it, but it works!
			if (varCount >= 2) {
				altMultiCount++;
			}
			if (varCount >= 3 && i != 2) {
				superAltMultiCount++;
			}
		}
		//errorWriter.println("Vertical PPDN");
		//errorWriter.println(sibParVarCount);
		//errorWriter.println(probVarCount);
		//errorWriter.println(superAltMultiCount);
		//errorWriter.println(altMultiCount);
		
		/// If variants in the siblings and parents are together greater than the
		/// probands', something's wrong.
		if (sibParVarCount > probVarCount) {
			return false;
		}

		if (superAltMultiCount > 0) { // For ppl with more than 3 variant counts (outside of proband)
			return false;
		}
		if (altMultiCount >= 4) { // Parents, proband, and at least 1 sib have var calls
			return false;
		}
		return true;
	}

	/*
	 * Because each person has only two of each chromosome (assuming no mosaicism or
	 * deletions/duplications), there should only be at most 2 haplotypes in a given
	 * region, one for each chromosome--ref and alt/variant. If there are more (i.e.
	 * multiple reads don't match up with each other), then it is likely the read
	 * was mis-aligned/mis-mapped, which can cause confetti regions. position, the
	 * call is ref (refMap) or not (altMap). The key is the read sequence itself
	 * (readString) and the value(s) are the VOI position's coordinates within the
	 * readString ([posVarStart, posVarEnd]). This gets a bit tricky given the
	 * potential of soft clipping and/or indels within the raw read sequence. In
	 * order to mark where indels have occurred in the sequence, the readString has
	 * been padded with "D"'s for deletions and "I"'s for insertions. These do not
	 * change the readString sequence overall, but rather mark where in the
	 * readString those indels may be.
	 * 
	 * Once all the reads have been mapped to their respective maps (refMap or
	 * altMap), then all reads with variant reads at the VOI position (altMap) are
	 * compared with one another. If the two reads differ from each other within the
	 * positions that they align (as anchored by the VOI position), then they
	 * represent different haplotypes. On the other hand, if they agree within all
	 * those positions, then they are part of the same haplotype.
	 * 
	 * If a read differs from every other read, then there's a potential that the
	 * former is just a badly aligned read and the others form their own. If there
	 * are too many of those "one off reads" (here determined to be either 10+ OR at
	 * least the total number of reads - 2 (provided there are at least 5 reads)),
	 * then the region is just overall poorly aligned and mismapped. If the number
	 * of disagreeing query reads is at least twice the number of agreeing query
	 * reads for a given subject read, then the subject read forms a unique
	 * haplotype. If there are at least 3 different haplotypes within one region for
	 * a person, then the region is considered confetti.
	 */
	/*
	 * public static boolean haplotypeCount(String chrom, int pos,
	 * ArrayList<SamReader> samReader, int length, String[] family, String ref) {
	 * SAMRecordIterator iterator; SAMRecord samRecord;
	 * 
	 * 
	 * /// Positions of interest include +/-50 of variant's of interest int
	 * startBound = pos; int endBound = pos + length;
	 * 
	 * /// Iterate through each family member for (int i = 0; i < samReader.size();
	 * i++) { // System.out.println("Family member: " + family[i]); // Extract the
	 * reads that are contained within the startBound, endBound range. // Note: This
	 * DOESN'T limit the readString to within startBound, endBound. It // gives the
	 * full string if a portion of it is contained within the range. iterator =
	 * samReader.get(i).queryOverlapping(chrom, startBound, endBound);
	 * 
	 * HashMap<String, ArrayList<Integer[]>> altMap = new HashMap<String,
	 * ArrayList<Integer[]>>(); HashMap<String, ArrayList<Integer[]>> refMap = new
	 * HashMap<String, ArrayList<Integer[]>>();
	 * 
	 * // Iterate through each read int readCount = 0; samReaderBreak: while
	 * (iterator.hasNext()) { samRecord = iterator.next(); // First, assess map
	 * quality. If fail, reject the whole read and move to the // next. int mapQ =
	 * samRecord.getMappingQuality(); if (mapQ < 30) { if (iterator.hasNext()) {
	 * samRecord = iterator.next(); } else { iterator.close(); return false; } } ///
	 * Determine where the variant of interest (VOI) is in the read int start =
	 * samRecord.getAlignmentStart(); int end = samRecord.getAlignmentEnd();
	 * 
	 * // if the read doesn't contain the variant, move to the next read while
	 * (((end < pos) || (pos < start))) { if (iterator.hasNext()) { samRecord =
	 * iterator.next(); start = samRecord.getAlignmentStart(); end =
	 * samRecord.getAlignmentEnd(); } else { break samReaderBreak; } }
	 * 
	 * /// Extract cigar string String cigar = samRecord.getCigar().toString();
	 * String read = samRecord.getReadString();
	 * 
	 * /// Eliminate the soft and hard clipping using the cigar string // Replace
	 * the alphabet with ";" so only numbers remain String cigarnumber =
	 * cigar.replaceAll("\\D+", ";"); // Get rid of the last ";" cigarnumber =
	 * cigarnumber.substring(0, cigarnumber.length() - 1); // Extract all the
	 * numbers String[] cigarN = cigarnumber.split(";"); // Replace all the number
	 * in a cigar string as ";," get rid of the heading ";" String cigarletter =
	 * cigar.replaceAll("[^A-Za-z]+", ";").substring(1); // Extract all the letters
	 * String[] cigarL = cigarletter.split(";");
	 * 
	 * // Anchor where the VOI's position is in the readString int posVarStart = pos
	 * - start; int posVarEnd = posVarStart + length;
	 * 
	 * int whereread = 0; // where in the readString (relative) int wheregenome =
	 * start - 1; // where in the genome (actual) String reads = ""; /// Iterate
	 * through the cigar string for (int ii = 0; ii < cigarL.length; ii++) { if
	 * (cigarL[ii].equals("S")) { // Soft clipping if (ii == 0) { // Soft clipped at
	 * the beginning of the read whereread += Integer.parseInt(cigarN[ii]); } } else
	 * if (cigarL[ii].equals("H")) { // Hard clipping...no adjustment needed
	 * 
	 * } else if (cigarL[ii].equals("M")) { // Matching wheregenome +=
	 * Integer.parseInt(cigarN[ii]); // Read in the rest of the substring from
	 * whereread for (int iii = 0; iii < Integer.parseInt(cigarN[ii]); iii++) { if
	 * (whereread < read.length()) { reads +=
	 * Character.toString(read.charAt(whereread)); } whereread++; }
	 * 
	 * if ((ii + 1) < cigarL.length) { if (cigarL[ii + 1].equals("I")) { // Next in
	 * cigar is insertion for (int iii = 0; iii < Integer.parseInt(cigarN[ii + 1]);
	 * iii++) { // Add inserted sequence to final readString reads += "I"; // Mark
	 * where the insert is in the read } reads +=
	 * Character.toString(read.charAt(whereread));
	 * 
	 * // If the insert is part of the VOI, then need to update VOI end whereread +=
	 * Integer.parseInt(cigarN[ii + 1]); } else if (cigarL[ii + 1].equals("D")) { //
	 * Next in cigar is deletion for (int iii = 0; iii < Integer.parseInt(cigarN[ii
	 * + 1]); iii++) { reads += "D"; // Pad deletion area with D's } // If D padding
	 * occurs before VOI, then need to update VOI's position in // the readString //
	 * If the D padding is part of the VOI, then need to update VOI end if
	 * (whereread >= posVarStart && whereread <= posVarEnd) { posVarEnd +=
	 * Integer.parseInt(cigarN[ii + 1]); } whereread += Integer.parseInt(cigarN[ii +
	 * 1]); } } /// Determine if the base(s) at posVar is ref or alt // TODO: These
	 * posVarStart and posVarEnd values are gonna be updated thanks to // the indels
	 * and soft clipping....this is gonna be done above } else if
	 * (cigarL[ii].equals("D")) { // deletion wheregenome +=
	 * Integer.parseInt(cigarN[ii]); whereread -= Integer.parseInt(cigarN[ii]); //
	 * bc we padded with D's earlier? } else if (cigarL[ii].equals("I")) { //
	 * insertion whereread += Integer.parseInt(cigarN[ii]); }
	 * 
	 * // Reached end of cigar string if (ii == (cigarL.length - 1)) { Integer[]
	 * posVarArr = { posVarStart, posVarEnd }; if ((posVarEnd + 1) < reads.length())
	 * { // if end of read doesn't include end of variant // Add 1 to posVarEnd
	 * because end index in substring is exclusive if (reads.substring(posVarStart,
	 * posVarEnd + 1).equals(ref)) { if (!refMap.containsKey(reads)) {
	 * ArrayList<Integer[]> initPosVarArr = new ArrayList<Integer[]>();
	 * initPosVarArr.add(posVarArr); refMap.put(reads, initPosVarArr); } else { ///
	 * It's okay if updatePosVarArr already contains posVarArr; /// it means
	 * multiple reads were sequenced with the same start /// and end positions.
	 * ArrayList<Integer[]> updatePosVarArr = refMap.get(reads);
	 * updatePosVarArr.add(posVarArr); refMap.put(reads, updatePosVarArr); } } else
	 * { if (!altMap.containsKey(reads)) { ArrayList<Integer[]> initPosVarArr = new
	 * ArrayList<Integer[]>(); initPosVarArr.add(posVarArr); altMap.put(reads,
	 * initPosVarArr); } else { ArrayList<Integer[]> updatePosVarArr =
	 * altMap.get(reads); updatePosVarArr.add(posVarArr); altMap.put(reads,
	 * updatePosVarArr); } } } } } readCount++; } iterator.close();
	 * 
	 * int altTotalCount = 0; // Total number of variants at VOI for (String read :
	 * altMap.keySet()) { // System.out.println("read: " + read); altTotalCount +=
	 * altMap.get(read).size(); }
	 * 
	 * ArrayList<String> altReadList = new ArrayList<String>(altMap.keySet());
	 * ArrayList<Integer[]> squabbleList = new ArrayList<Integer[]>();
	 * ArrayList<Double> squabbleRatioList = new ArrayList<Double>(); int haploCount
	 * = 0; int oneOffReads = 0;
	 * 
	 * // System.out.println("Number of reads: " + readCount +
	 * ", number alt reads: " + // altTotalCount + ", altReadList size: " +
	 * altReadList.size());
	 * 
	 * /// Iterate through the altMap to determine mismatches for (int j = 0; j <
	 * altReadList.size(); j++) { int disagreeCount = 0; // Number of times subject
	 * read differs from a query read int agreeCount = 0; // Number of times no
	 * mismatches between subject and query reads // System.out.println("j: " + j);
	 * String subjectRead = altReadList.get(j);
	 * 
	 * // Combine with the others in the altMap, but don't have to repeat combos.
	 * ArrayList<Integer[]> posV1ArrList = altMap.get(subjectRead); //
	 * System.out.println("altMap arr1: " + list2StrIntArr(posV1ArrList)); for
	 * (Integer[] posV1Arr : posV1ArrList) { //
	 * System.out.println("Disagree count pre: " + disagreeCount + ", Agree count: "
	 * // + agreeCount);
	 * 
	 * // System.out.println("posV1Arr: " + Arrays.toString(posV1Arr));
	 * 
	 * for (int jj = j + 1; jj < altReadList.size(); jj++) { //
	 * System.out.println("jj: " + jj);
	 * 
	 * String queryRead = altReadList.get(jj); ArrayList<Integer[]> posV2ArrList =
	 * altMap.get(queryRead); // query read // System.out.println("altMap arr2: " +
	 * list2StrIntArr(posV2ArrList));
	 * 
	 * compareBreak: for (Integer[] posV2Arr : posV2ArrList) { /// Establish start
	 * and end boundaries of the comparison // Start int posV1Start = posV1Arr[0];
	 * int posV2Start = posV2Arr[0]; // End // int posV1End = posV1Arr[1]; // int
	 * posV2End = posV2Arr[1];
	 * 
	 * /// Iterate through comparing each base radiating out from the VOI in each
	 * /// string. However, if you hit the end of a string, break.
	 * 
	 * // Backwards direction int disagreeCheck = 0; for (int jjj = 0; jjj <= 50;
	 * jjj++) { String readV1 = null; String readV2 = null;
	 * 
	 * if (((posV1Start - jjj) > 0) && ((posV2Start - jjj) > 0)) { readV1 =
	 * Character.toString(subjectRead.charAt(posV1Start - jjj)); readV2 =
	 * Character.toString(queryRead.charAt(posV2Start - jjj)); } // Check if
	 * starting bound of string < 50 away from VOI // If it is, then break out of
	 * the loop/stop comparing this set of reads. if (readV1 != null && readV2 !=
	 * null) { if ((!readV1.equals(readV2))) { //
	 * System.out.println("back disagree! jjj: " + jjj); disagreeCount++; //
	 * System.out.println("DISAGREE!"); break compareBreak; } } }
	 * 
	 * /// Forwards direction // System.out.println("Moving to the fwd");
	 * forwardBreak: for (int jjjj = 0; jjjj <= 50; jjjj++) { String readV1 = null;
	 * String readV2 = null; // If hit the end of one readString, then end it. if
	 * (((posV1Start + jjjj) < subjectRead.length()) && ((posV2Start + jjjj) <
	 * queryRead.length())) { readV1 =
	 * Character.toString(subjectRead.charAt(posV1Start + jjjj)); readV2 =
	 * Character.toString(queryRead.charAt(posV2Start + jjjj)); } else { // if
	 * ((posV1Start + jjjj) >= subjectRead.length()) { //
	 * System.out.println("greater than subject length of " + subjectRead.length() +
	 * // ": " + (posV1Start + jjjj)); // } else if ((posV2Start + jjjj) >=
	 * queryRead.length()) { // System.out.println("greater than query length of " +
	 * queryRead.length() + ": // " + (posV2Start + jjjj)); // } //
	 * System.out.println("wtf , jjjj: " + jjjj); // // if (jjjj == 50 ||
	 * readV1.equals("") || readV2.equals("") // || readV1.equals(null) ||
	 * readV2.equals(null) // || ((posV1Start + jjjj) == subjectRead.length()) // ||
	 * ((posV2Start + jjjj) == queryRead.length())) { //
	 * System.out.println("readV1: " + readV1 + "; readV2: " + readV2); //
	 * System.out.println("end! jjjj: " + jjjj); // // If hit the end bound without
	 * mismatches, then assume it means reads match // agreeCount++; //
	 * System.out.println("agree!"); // break forwardBreak; // } } // Check if
	 * starting bound of string < 50 away from VOI // If it is, then break out of
	 * the loop. if (readV1 != null && readV2 != null) { if
	 * ((!readV1.equals(readV2))) { // System.out.println("fwd disagree! jjjj: " +
	 * jjjj); disagreeCount++; // System.out.println("DISAGREE!"); break
	 * compareBreak; } else if (jjjj == 50) { agreeCount++; //
	 * System.out.println("agreed bc reached the end!"); //
	 * System.out.println("agree!"); break forwardBreak; } } else if (readV1 == null
	 * || readV2 == null || ((posV1Start + jjjj) == subjectRead.length()) ||
	 * ((posV2Start + jjjj) == queryRead.length())) { //
	 * System.out.println("subject length of " + subjectRead.length() + ": " + //
	 * (posV1Start + jjjj)); // System.out.println("query length of " +
	 * queryRead.length() + ": " + // (posV2Start + jjjj)); //
	 * System.out.println("readV1: " + readV1 + "; readV2: " + readV2); //
	 * System.out.println("end! jjjj: " + jjjj); // If hit the end bound without
	 * mismatches, then assume it means reads match agreeCount++; //
	 * System.out.println("agree!"); break forwardBreak; // // Remove from altMap
	 * because would match with others that matched // subjectRead // //
	 * altMap.remove(queryRead); // // altReadList.remove(queryRead); } //
	 * System.out.println("jjjj: " + jjjj + ", readV1: " + readV1 + ", readV2: " +
	 * // readV2); } } } // Integer[] squabbleArr = { disagreeCount, agreeCount };
	 * 
	 * Double squabbleRatio = (double) disagreeCount / agreeCount; //
	 * System.out.println("disagreeCount: " + disagreeCount + " agreeCount: " + //
	 * agreeCount // + " squabble ratio: " + squabbleRatio);
	 * squabbleRatioList.add(squabbleRatio);
	 * 
	 * // If every other read differs from the subject read, then likely subject is
	 * // bad. // if (squabbleRatio == Double.POSITIVE_INFINITY && disagreeCount >
	 * 1) { if (squabbleRatio == Double.POSITIVE_INFINITY) { // If there are at
	 * least 10 reads that don't match with others or all reads in // the altMap
	 * differ from each other, it's a bad pile up oneOffReads++; } }
	 * 
	 * } /// *** MOVE THIS AND SEE WHERE IT FITS*** //
	 * System.out.println("one off read count: " + oneOffReads); if (oneOffReads >=
	 * 10 || (oneOffReads >= (altTotalCount - 2) && altTotalCount > 5)) {
	 * System.out.println("too many one off reads!"); return false; }
	 * 
	 * /// Determine how "disagreeable" the pileup is (aka how bad the confetti is).
	 * // If the ratio is about 50%, then that means about 2 haplotypes (half //
	 * disagreed, half agreed?). One off reads don't count though. for (double ratio
	 * : squabbleRatioList) { // DisagreeCount : agreeCount for each read if (ratio
	 * >= 2 && (ratio != Double.POSITIVE_INFINITY)) { haploCount++; } }
	 * 
	 * if (haploCount >= 3) { // More than ref and VOI haplotypes
	 * System.out.println("haplotype count: " + haploCount); return false; } }
	 * return true;
	 * 
	 * }
	 */

	/*
	 * Because each person has only two of each chromosome (assuming no mosaicism or
	 * deletions/duplications), there should only be at most 2 haplotypes in a given
	 * region, one for each chromosome--ref and alt/variant. If there are more (i.e.
	 * multiple reads don't match up with each other), then it is likely the read
	 * was mis-aligned/mis-mapped, which can cause confetti regions. position, the
	 * call is ref (refMap) or not (altMap). The key is the read sequence itself
	 * (readString) and the value(s) are the VOI position's coordinates within the
	 * readString ([posVarStart, posVarEnd]). This gets a bit tricky given the
	 * potential of soft clipping and/or indels within the raw read sequence. In
	 * order to mark where indels have occurred in the sequence, the readString has
	 * been padded with "D"'s for deletions and "I"'s for insertions. These do not
	 * change the readString sequence overall, but rather mark where in the
	 * readString those indels may be.
	 * 
	 * Once all the reads have been mapped to their respective maps (refMap or
	 * altMap), then all reads with variant reads at the VOI position (altMap) are
	 * compared with one another. If the two reads differ from each other within the
	 * positions that they align (as anchored by the VOI position), then they
	 * represent different haplotypes. On the other hand, if they agree within all
	 * those positions, then they are part of the same haplotype.
	 * 
	 * If a read differs from every other read, then there's a potential that the
	 * former is just a badly aligned read and the others form their own. If there
	 * are too many of those "one off reads" (here determined to be either 10+ OR at
	 * least the total number of reads - 2 (provided there are at least 5 reads)),
	 * then the region is just overall poorly aligned and mismapped. If the number
	 * of disagreeing query reads is at least twice the number of agreeing query
	 * reads for a given subject read, then the subject read forms a unique
	 * haplotype. If there are at least 3 different haplotypes within one region for
	 * a person, then the region is considered confetti.
	 */
	public boolean haploCountEngine(String chrom, int pos, ArrayList<SamReader> samReader, int length,
			String[] family, String ref) throws IOException {
		SAMRecordIterator iterator;
		SAMRecord samRecord;

		/// Positions of interest include +/-50 of variant's of interest
		int startBound = pos;
		int endBound = pos + length;

		/// Iterate through each family member
		for (int i = 0; i < samReader.size(); i++) {
			
			//Translator - changes chromosome input depending on the type of BAM chrs the BAMs take in
			if(!bamChrHeaders.get(i).contains(chrom)) {
				
				chrom = BamChrChanger.translate(chrom);
				
			}
			
			try {
				// Extract the reads that are contained within the startBound, endBound range.
				iterator = samReader.get(i).queryOverlapping(chrom, startBound, endBound);
	
				HashMap<String, HashSet<Integer[]>> altMap = new HashMap<String, HashSet<Integer[]>>();
				HashMap<String, HashSet<Integer[]>> refMap = new HashMap<String, HashSet<Integer[]>>();
	
				HashMap<String, Integer[]> readPosMap = new HashMap<String, Integer[]>();
				HashMap<String, ArrayList<Integer>> readQualMap = new HashMap<String, ArrayList<Integer>>();
	
				/// Iterate through each read and extract sequence using CIGAR string
				int readCount = 0;
				//errorWriter.println("Haplo check: ");
				samReaderBreak: while (iterator.hasNext()) {
					samRecord = iterator.next();
					// First, assess map quality. If too poor, reject the whole read and move to the
					// next.
					int mapQ = samRecord.getMappingQuality();
					if (mapQ < 30) {
						if (iterator.hasNext()) {
							samRecord = iterator.next();
						} else {
							//errorWriter.println("MAPQ 2 smol: " + mapQ);
							iterator.close();
							return false;
						}
					}
	
					/// Extract base quality array of the string
					byte[] baseQArr = samRecord.getBaseQualities();
	
					/// Determine where the variant of interest (VOI) is in the read
					int start = samRecord.getAlignmentStart();
					int end = samRecord.getAlignmentEnd();
	
					// If the read doesn't contain the variant position, move to the next read
					while (((end < pos) || (pos < start))) {
						if (iterator.hasNext()) {
							samRecord = iterator.next();
							start = samRecord.getAlignmentStart();
							end = samRecord.getAlignmentEnd();
						} else {
							break samReaderBreak;
						}
					}
	
					/// Extract cigar string
					String cigar = samRecord.getCigar().toString();
					String read = samRecord.getReadString();
	
					/// Eliminate the soft and hard clipping using the cigar string
					// Replace the alphabet with ";" so only numbers remain
					String cigarnumber = cigar.replaceAll("\\D+", ";");
					// Get rid of the last ";"
					cigarnumber = cigarnumber.substring(0, cigarnumber.length() - 1);
					// Extract all the numbers
					String[] cigarN = cigarnumber.split(";");
					// Replace all the number in a cigar string as ";," get rid of the heading ";"
					String cigarletter = cigar.replaceAll("[^A-Za-z]+", ";").substring(1);
					// Extract all the letters
					String[] cigarL = cigarletter.split(";");
	
					/// Anchor where the VOI's position is in the readString
					int posVarStart = pos - start;
					int posVarEnd = posVarStart + length;
	
					/// CIGAR string interpretation
					int whereread = 0; // where in the readString (relative)
					String reads = "";
					ArrayList<Integer> baseQualList = new ArrayList<Integer>();
	
					if (baseQArr.length == read.length()) {
						for (int ii = 0; ii < cigarL.length; ii++) {
							if (cigarL[ii].equals("S")) { // Soft clipping
								if (ii == 0) { // Soft clipped at the beginning of the read
									whereread += Integer.parseInt(cigarN[ii]);
								}
							} else if (cigarL[ii].equals("H")) { // Hard clipping...no adjustment needed
	
							} else if (cigarL[ii].equals("M")) { // Matching
								// Read in the rest of the substring from whereread
								for (int iii = 0; iii < Integer.parseInt(cigarN[ii]); iii++) {
									if (whereread < read.length()) {
										reads += Character.toString(read.charAt(whereread));
	
										byte baseQByte = baseQArr[whereread];
										baseQualList.add((int) baseQByte);
									}
									whereread++;
								}
	
								if ((ii + 1) < cigarL.length) {
									if (cigarL[ii + 1].equals("I")) { // Next in cigar is insertion
										for (int iii = 0; iii < Integer.parseInt(cigarN[ii + 1]); iii++) {
											reads += "I"; // Mark where the insert is in the read
	
											baseQualList.add(20); /// Default value for inserted bases
	
										}
										reads += Character.toString(read.charAt(whereread));
	
										byte baseQByte = baseQArr[whereread];
										baseQualList.add((int) baseQByte);
	
										whereread += Integer.parseInt(cigarN[ii + 1]);
									} else if (cigarL[ii + 1].equals("D")) { // Next in cigar is deletion
										for (int iii = 0; iii < Integer.parseInt(cigarN[ii + 1]); iii++) {
											reads += "D"; // Pad deletion area with D's
	
											baseQualList.add(20); /// Default value for deleted bases
										}
										/// If D padding occurs before VOI, then need to update VOI's position in
										/// the readString
										if (whereread >= posVarStart && whereread <= posVarEnd) {
											posVarEnd += Integer.parseInt(cigarN[ii + 1]);
										}
										whereread += Integer.parseInt(cigarN[ii + 1]);
									}
								}
							} else if (cigarL[ii].equals("D")) { // deletion
								whereread -= Integer.parseInt(cigarN[ii]); // bc we padded with D's earlier
							} else if (cigarL[ii].equals("I")) { // insertion
								whereread += Integer.parseInt(cigarN[ii]);
							}
	
							/// Reached end of cigar string
							if (ii == (cigarL.length - 1)) {
								Integer[] posVarArr = { posVarStart, posVarEnd };
								if ((posVarEnd + 1) < reads.length()) { // if end of read doesn't include end of variant
									/// Add 1 to posVarEnd because end index in substring is exclusive
									if (reads.substring(posVarStart, posVarEnd + 1).equals(ref)) {
										if (!refMap.containsKey(reads)) {
											HashSet<Integer[]> initPosVarArr = new HashSet<Integer[]>();
											initPosVarArr.add(posVarArr);
											refMap.put(reads, initPosVarArr);
										} else {
											/// It's okay if updatePosVarArr already contains posVarArr;
											/// it means multiple reads were sequenced with the same start
											/// and end positions.
											HashSet<Integer[]> updatePosVarArr = refMap.get(reads);
											updatePosVarArr.add(posVarArr);
											refMap.put(reads, updatePosVarArr);
										}
									} else {
										if (!altMap.containsKey(reads)) {
											HashSet<Integer[]> initPosVarArr = new HashSet<Integer[]>();
											initPosVarArr.add(posVarArr);
											altMap.put(reads, initPosVarArr);
										} else {
											HashSet<Integer[]> updatePosVarArr = altMap.get(reads);
											updatePosVarArr.add(posVarArr);
											altMap.put(reads, updatePosVarArr);
										}
									}
								}
								/// Maintain a map of each read and their alignment start/end (for later use)
								Integer[] readBounds = { start, end };
	
								readPosMap.put(reads, readBounds);
								readQualMap.put(reads, baseQualList);
							}
						}
						readCount++;
					}
				}
				iterator.close();
	
				/// Get haploCount and one off read count
				int[] altHaploArr = haploCounter(altMap, pos, readPosMap, readQualMap, readCount);
				int altHaploCount = altHaploArr[0];
				int altOneOffCount = altHaploArr[1];
				
				int[] refHaploArr = haploCounter(refMap, pos, readPosMap, readQualMap, readCount);
				int refHaploCount = refHaploArr[0];
				int refOneOffCount = refHaploArr[1];
	
				int combinedOneOffCount = altOneOffCount + refOneOffCount;
				int combinedHaploCount = altHaploCount + refHaploCount;
				
				//errorWriter.println(pos);
				//errorWriter.println(readCount);
				//errorWriter.println(altOneOffCount);
				//errorWriter.println(refOneOffCount);
				//errorWriter.println(altHaploCount);
				//errorWriter.println(refHaploCount);
				
				// System.out.println("combined haplocount: " + combinedHaploCount);
				// System.out.println("combined one off count: " + combinedOneOffCount);
				
				/// NOTE: TALK ABT THIS W TOM
				if (combinedOneOffCount >= 10 && combinedOneOffCount >= (0.2 * readCount)) {
					// System.out.println("too many one offs: " + combinedOneOffCount);
					return false;
				}
	
				if (combinedHaploCount >= 3) {
					// System.out.println("too many haplotypes: " + "alt: " + altHaploCount + ",
					// ref: " + refHaploCount);
					return false;
				}
			}  catch(IllegalArgumentException j) { //Sends Error Message if translation does not work
				System.gc();
				JOptionPane.showMessageDialog(new Frame("Error"), "The Sam Reader passed an invalid Chromosome input to"
						+ " the BAM file. The chromosome sent was: " +  chrom + ". The BAM header format will be outputted as \"BAM_Header_Error.txt\"."
								+ "Please edit the \"BAM_Chr_Config\" accordingly. The Program will now exit");
				
				BamChrChanger.writeHeaderErrors(bamChrHeaders.get(i));
				
				System.exit(0);
				
				
				
			}
		}

		return true;
	}

	public int[] haploCounter(HashMap<String, HashSet<Integer[]>> readMap, int varPos,
			HashMap<String, Integer[]> readPosMap, HashMap<String, ArrayList<Integer>> readQualMap, int readCount) {
		int haploCount = 0;

		ArrayList<String> readList = new ArrayList<String>(readMap.keySet());

		ArrayList<HashSet<String>> haploReadList = new ArrayList<HashSet<String>>();
		HashSet<String> prevReadSet = new HashSet<String>();

		/// Initialize map of previously read reads.
		/// Genomic position of variant -->List{Variant calls}
		HashMap<Integer, ArrayList<PutativeVarCall>> potVarMap = new HashMap<Integer, ArrayList<PutativeVarCall>>();

		HashMap<Integer, HashSet<String>> potVarCallMap = new HashMap<Integer, HashSet<String>>();

		/// Iterate through the altMap to determine mismatches
		for (int j = 0; j < readList.size(); j++) {
			// System.out.println("j: " + j);
			HashSet<String> prevReadList = new HashSet<String>();
			String subjectRead = readList.get(j);

			// Combine with the others in the altMap, but don't have to repeat combos.
			HashSet<Integer[]> posV1ArrList = readMap.get(subjectRead);
			ArrayList<Integer> posV1QualList = readQualMap.get(subjectRead);

			for (Integer[] posV1Arr : posV1ArrList) {
				for (int jj = 0; jj < readList.size(); jj++) {
					// System.out.println("jj: " + jj);
					String queryRead = readList.get(jj);
					HashSet<Integer[]> posV2ArrList = readMap.get(queryRead); // query read
					ArrayList<Integer> posV2QualList = readQualMap.get(queryRead);

					for (Integer[] posV2Arr : posV2ArrList) {
						/// Establish start and end boundaries of the comparison
						// Start
						int posV1Start = posV1Arr[0];
						int posV2Start = posV2Arr[0];

						boolean seedBack = false;
						boolean seedFront = false;
						boolean seedAgrees = true;

						/// Iterate through comparing each base radiating out from the VOI in each
						/// string. However, if you hit the end of a string, break.

						// Backwards direction
						for (int jjj = 0; jjj <= 50; jjj++) {

							String readV1 = null;
							String readV2 = null;

							int qualV1 = (int) Double.POSITIVE_INFINITY;
							int qualV2 = (int) Double.POSITIVE_INFINITY;

							if (((posV1Start - jjj) > 0) && ((posV2Start - jjj) > 0)) {
								readV1 = Character.toString(subjectRead.charAt(posV1Start - jjj));
								readV2 = Character.toString(queryRead.charAt(posV2Start - jjj));

								qualV1 = posV1QualList.get(posV1Start - jjj);
								qualV2 = posV2QualList.get(posV2Start - jjj);
							}
							// Check if starting bound of string < 30 away from VOI
							// If it is, then break out of the loop/stop comparing this set of reads.
							if (!(readV1 == null) && !(readV2 == null)) {
								if ((!readV1.equals(readV2)) && (qualV1 >= 20 && qualV2 >= 20)) {

									seedAgrees = false;
									int varPosition = varPos - jjj;

									/// Add to variant count map
									if ((!potVarMap.containsKey(varPosition) || potVarMap.isEmpty())) { /// New position

										PutativeVarCall varCall1 = new PutativeVarCall(readV1, varPosition,
												subjectRead);

										PutativeVarCall varCall2 = new PutativeVarCall(readV2, varPosition, queryRead);

										ArrayList<PutativeVarCall> varArrList = new ArrayList<PutativeVarCall>();
										varArrList.add(varCall1);
										varArrList.add(varCall2);

										potVarMap.put(varPosition, varArrList);

										HashSet<String> varCallSet = new HashSet<String>();
										varCallSet.add(readV1);
										varCallSet.add(readV2);
										potVarCallMap.put(varPosition, varCallSet);

									} else { // Variant position previously seen
										ArrayList<PutativeVarCall> varArrList = potVarMap.get(varPosition);
										HashSet<String> varCallSet = potVarCallMap.get(varPosition);

										for (ListIterator<PutativeVarCall> iter = varArrList.listIterator(); iter
												.hasNext();) {

											PutativeVarCall potVar = iter.next();

											// if (varPosition == 123171558 && potVar.getCall().equals("A")) {
											// System.out.println("readV1: " + readV1 + "; readV2: " + readV2
											// + "; potVar: " + potVar.getCall() + "; varPos: " + varPosition
											// + ", other Calls: " + set2Str(varCallSet) + ", matches readV1: "
											// + (readV1.equals(potVar.getCall())) + ", matches readV2: "
											// + (readV2.equals(potVar.getCall()))
											// + ", subject matches potVar PRS: "
											// + (potVar.prevReadSet.contains(subjectRead))
											// + ", query matches potVar PRS: "
											// + (potVar.prevReadSet.contains(queryRead))
											// + "; varArrList.size():" + varArrList.size() + ", setSize: "
											// + potVar.getCount() + ": " + potVar.printSet() + ", subject: "
											// + subjectRead + ", query: " + queryRead);
											// }
											// System.out.println("varPos: " + varPosition + "; " + "varArrList.size():
											// "
											// + varArrList.size());
											// if (!prevReadList.contains(subjectRead)
											// || !prevReadList.contains(queryRead)) { /// Either read hasn't been
											/// eval'd yet
											if (readV1.equals(potVar.getCall())) {
												// ReadV1 variant call previously seen at this pos
												potVar.addRead(subjectRead);
												prevReadList.add(subjectRead);

											} else if (readV2.equals(potVar.getCall())) {
												// ReadV2 variant call previously seen at this pos
												potVar.addRead(queryRead);
												prevReadList.add(subjectRead);

											} else {
												/// readV1 != readV2 != potVar call

												if (!potVar.prevReadSet.contains(subjectRead)) {

													if (!varCallSet.contains(readV1)) {
														PutativeVarCall varCall1 = new PutativeVarCall(readV1,
																varPosition, subjectRead);
														iter.add(varCall1);
														varCallSet.add(readV1);
														prevReadList.add(subjectRead);

													}
												} else if (!potVar.prevReadSet.contains(queryRead)) {

													if (!varCallSet.contains(readV2)) {
														PutativeVarCall varCall2 = new PutativeVarCall(readV1,
																varPosition, subjectRead);
														iter.add(varCall2);
														varCallSet.add(readV2);
														prevReadList.add(queryRead);
													}
												}
											}
										}
									}
								}
								if (jjj == 15) { /// Hit end of seed region
									seedBack = true;
								}
							}
						}

						// Forwards direction
						forwardBreak: for (int jjjj = 0; jjjj <= 50; jjjj++) {

							String readV1 = null;
							String readV2 = null;
							// If hit the end of one readString, then end it.
							if (((posV1Start + jjjj) < subjectRead.length())
									&& ((posV2Start + jjjj) < queryRead.length())) {
								readV1 = Character.toString(subjectRead.charAt(posV1Start + jjjj));
								readV2 = Character.toString(queryRead.charAt(posV2Start + jjjj));
							}
							// Check if end bound of string > 25 away from VOI
							// If it is, then break out of the loop.
							if (!(readV1 == null) && !(readV2 == null)) {
								if ((!readV1.equals(readV2))) {

									seedAgrees = false;
									int varPosition = varPos + jjjj;
									/// Add to variant count map

									if ((!potVarMap.containsKey(varPosition) || potVarMap.isEmpty())) { /// New position

										PutativeVarCall varCall1 = new PutativeVarCall(readV1, varPosition,
												subjectRead);
										varCall1.addCall(readV1);
										varCall1.addCall(readV2);

										PutativeVarCall varCall2 = new PutativeVarCall(readV2, varPosition, queryRead);
										varCall2.addCall(readV1);
										varCall2.addCall(readV2);

										ArrayList<PutativeVarCall> varArrList = new ArrayList<PutativeVarCall>();
										varArrList.add(varCall1);
										varArrList.add(varCall2);

										potVarMap.put(varPosition, varArrList);

										HashSet<String> varCallSet = new HashSet<String>();
										varCallSet.add(readV1);
										varCallSet.add(readV2);

										potVarCallMap.put(varPosition, varCallSet);

									} else { // Variant position previously seen
										ArrayList<PutativeVarCall> varArrList = potVarMap.get(varPosition);
										HashSet<String> varCallSet = potVarCallMap.get(varPosition);
										
										/// THIS IS JUST A BUG WAITING TO HAPPEN. 3/12: WHAT DID I SAY
										
										for (ListIterator<PutativeVarCall> iter = varArrList.listIterator(); iter
												.hasNext();) {
											
											PutativeVarCall potVar = iter.next();

											// if (varPosition == 123171604 && potVar.getCall().equals("C")) {
											// System.out.println("readV1: " + readV1 + "; readV2: " + readV2
											// + "; potVar: " + potVar.getCall() + "; varPos: " + varPosition
											// + ", other Calls: " + set2Str(varCallSet) + ", matches readV1: "
											// + (readV1.equals(potVar.getCall())) + ", matches readV2: "
											// + (readV2.equals(potVar.getCall()))
											// + ", subject matches potVar PRS: "
											// + (potVar.prevReadSet.contains(subjectRead))
											// + ", query matches potVar PRS: "
											// + (potVar.prevReadSet.contains(queryRead))
											// + "; varArrList.size():" + varArrList.size() + ", setSize: "
											// + potVar.getCount() + ": " + potVar.printSet() + ", subject: "
											// + subjectRead + ", query: " + queryRead);
											// }
											// System.out.println("varPos: " + varPosition + "; " + "varArrList.size():
											// "
											// + varArrList.size());
											// if (!prevReadList.contains(subjectRead)
											// || !prevReadList.contains(queryRead)) { /// Either read hasn't been
											/// eval'd yet
											if (readV1.equals(potVar.getCall())) {
												// ReadV1 variant call previously seen at this pos
												potVar.addRead(subjectRead);
												prevReadList.add(subjectRead);

											} else if (readV2.equals(potVar.getCall())) {
												// ReadV2 variant call previously seen at this pos
												potVar.addRead(queryRead);
												prevReadList.add(subjectRead);

											} else {

												if (!potVar.prevReadSet.contains(subjectRead)) {

													if (!varCallSet.contains(readV1)) {
														PutativeVarCall varCall1 = new PutativeVarCall(readV1,
																varPosition, subjectRead);
														iter.add(varCall1);
														varCallSet.add(readV1);
														prevReadList.add(subjectRead);

													}
												} else if (!potVar.prevReadSet.contains(queryRead)) {

													if (!varCallSet.contains(readV2)) {
														PutativeVarCall varCall2 = new PutativeVarCall(readV1,
																varPosition, subjectRead);
														iter.add(varCall2);
														varCallSet.add(readV2);
														prevReadList.add(queryRead);
													}
												}
											}
										}
									}
								}
								if (jjjj == 15 && seedAgrees) { /// Front end of seed agrees
									seedFront = true;
								}
							} else if (readV1 == null || readV2 == null || ((posV1Start + jjjj) == subjectRead.length())
									|| ((posV2Start + jjjj) == queryRead.length())) {
								/// Seed matches and extended region still agrees
								if (seedFront && seedBack && seedAgrees) {
									if (haploReadList.isEmpty() || (haploReadList.size() == 0)) {
										HashSet<String> haploReads = new HashSet<String>();
										haploReads.add(subjectRead);
										haploReads.add(queryRead);
										haploReadList.add(haploReads);

										prevReadSet.add(subjectRead);
										prevReadSet.add(queryRead);
									} else {
										/// Quick and dirty list to avoid ConcurrentModificationException
										ArrayList<HashSet<String>> newHaplos = new ArrayList<HashSet<String>>();
										/// NOTE: This can definitely be done better with an iterator to also avoid
										/// ConcurrentModificationException.
										for (HashSet<String> haploReads : haploReadList) {
											if (prevReadSet.contains(subjectRead) && !prevReadSet.contains(queryRead)) {
												haploReads.add(queryRead);
												prevReadSet.add(queryRead);
											} else if (prevReadSet.contains(queryRead)
													&& !prevReadSet.contains(subjectRead)) {
												haploReads.add(subjectRead);
												prevReadSet.add(queryRead);
											} else if (!prevReadSet.contains(subjectRead)
													&& !prevReadSet.contains(queryRead)) {
												HashSet<String> haploReadNew = new HashSet<String>();
												haploReadNew.add(subjectRead);
												haploReadNew.add(queryRead);
												newHaplos.add(haploReadNew);
												// haploReadList.add(haploReadNew);

												prevReadSet.add(subjectRead);
												prevReadSet.add(queryRead);
											}
										}
										haploReadList.addAll(newHaplos);
									}

								}
								/// Reached end of either subject or query readString
								/// If hit the end bound without mismatches, then assume it means reads match

								break forwardBreak;
							}
						}
					}

				}
			}
		}

		// for (int potVarPos : potVarMap.keySet()) {
		// System.out.println(
		// "putative variant position: " + potVarPos + "; size: " +
		// potVarMap.get(potVarPos).size() + ": ");
		// ArrayList<PutativeVarCall> pvL = potVarMap.get(potVarPos);
		// for (PutativeVarCall pvc : pvL) {
		// if (pvc.getCount() == 1) {
		// System.out.print("putative var pos: " + potVarPos + "--> " + pvc.getCall() +
		// ", " + pvc.getCount()
		// + ": " + (pvc.printSet()) + " \n");
		// } else {
		// System.out.print(
		// "putative var pos: " + potVarPos + "--> " + pvc.getCall() + ", " +
		// pvc.getCount() + "; \n");
		// }
		// }
		// }

		/// Count and eliminate one off reads from evaluation
		HashSet<String> oneOffReadList = new HashSet<String>(); /// List of one off bad reads
		ArrayList<String> copyAltReadList = new ArrayList<String>(readList);

		/// *** WARNING: MUCH BUGS POTENTIALLY PRESENT
		boolean oneOffReadsPresent = true;
		while (oneOffReadsPresent) {
			oneOffReadsPresent = false; /// Reset boolean each time the potVarMap gets evaluated
			for (int potVarPos : potVarMap.keySet()) {
				ArrayList<PutativeVarCall> pvL = potVarMap.get(potVarPos);
				for (PutativeVarCall pvc : pvL) {
					if (pvc.getCount() <= 1 && copyAltReadList.contains(pvc.printSet().replaceAll("[^a-zA-Z]", ""))) {
						String oneOffRead = pvc.printSet().replaceAll("[^a-zA-Z]", "");
						Integer[] oneOffReadBound = readPosMap.get(oneOffRead);
						int oneOffReadStart = oneOffReadBound[0];
						int oneOffReadEnd = oneOffReadBound[1];
						oneOffReadList.add(oneOffRead);

						/// If the varCall is a one off bad read, then remove its read from the
						/// altReadList
						copyAltReadList.removeAll(oneOffReadList);
						/// Also remove that read from potVarMap
						for (int potVarPos2 : potVarMap.keySet()) {
							ArrayList<PutativeVarCall> pvL2 = potVarMap.get(potVarPos2);

							for (PutativeVarCall pvc2 : pvL2) {
								/// Remove from pvc2 if it's prevReadSet contains the one off and the genomic
								/// position of the pvc2 is within the boundaries of the one off read
								if (pvc2.prevReadSet.contains(oneOffRead)
										&& ((pvc.getPos() <= oneOffReadEnd) && (pvc.getPos() >= oneOffReadStart))) {
									pvc2.prevReadSet.remove(oneOffRead);
								}
							}
						}
						oneOffReadsPresent = true;
					}
				}
			}
		}

		// System.out.println("one off read count: " + oneOffReadList.size());
		// for (int potVarPos : potVarMap.keySet()) {
		//// System.out
		//// .println("putative variant position: " + potVarPos + "; size: " +
		// potVarMap.get(potVarPos).size());
		// ArrayList<PutativeVarCall> pvL = potVarMap.get(potVarPos);
		// for (PutativeVarCall pvc : pvL) {
		//// System.out.print(
		//// "putative var pos: " + potVarPos + "--> " + pvc.getCall() + ", " +
		// pvc.getCount() + "; \n");
		// // ": " + pvc.printSet() +
		// }
		// }

		for (Iterator<HashSet<String>> it = haploReadList.iterator(); it.hasNext();) {
			HashSet<String> next = it.next();
			/// Remove one off reads and resulting empty sets from the haploblock list
			if (!Collections.disjoint(next, oneOffReadList)) {
				next.removeAll(oneOffReadList);
			}
			if (next.isEmpty()) {
				it.remove();
			}
			if (next.size() >= 0.2 * readCount) {
				haploCount++;
			}
		}
		// System.out.println("haploCount: " + haploCount);
		// System.out.println("oneoffcount: " + oneOffReadList.size());
		int[] haploArr = { haploCount, oneOffReadList.size() };
		return haploArr;
	}

	/*
	 * Read depth filter. If the read depth is too low, then the likelihood that the
	 * correct genotype call was just missed during sequencing goes up. Example,
	 * sequence coverage is 5 reads vs. 30 reads. If there were 2 calls in 5 vs. 2
	 * calls in 30...which do you think is more likely that the genotype call of ref
	 * was more correct? In order to avoid issues of allele skew due to abundance of
	 * relatively bad reads from one family member, should also include in the
	 * filter something regarding bam quality/noise..."Error" or "SNR" should be
	 * sufficient.
	 */
	public boolean readDepthEval(String chrom, int pos, ArrayList<SamReader> samReader, int length) throws IOException {
		// TODO: Write out the pseudocode first and think about allele skew
		SAMRecordIterator iterator;
		SAMRecord samRecord;

		/// Positions of interest (of the VOI)
		int startBound = pos;
		int endBound = pos + length;

		/// Iterate through each family member
		for (int i = 0; i < samReader.size(); i++) {
			// Extract the reads that are contained within the startBound, endBound range.
			// Note: This DOESN'T limit the readString to within startBound, endBound. It
			// gives the full string if a portion of it is contained within the range.
			
			//Translator - changes chromosome input depending on the type of BAM chrs the BAMs take in
			if(!bamChrHeaders.get(i).contains(chrom)) {
				
				chrom = BamChrChanger.translate(chrom);
				
			}
			
			try {
				
				iterator = samReader.get(i).queryOverlapping(chrom, startBound, endBound);
				int readDepth = 0;
				while (iterator.hasNext()) {
					readDepth++;
					samRecord = iterator.next();
				}
				iterator.close();
				
				// If poor read depth, then remove
				/// TODO: The read depth min of 8 is a heuristic. We may need to take into
				// account how often the variant has appeared (or not) in ExAC and/or gnomAD.
				if (readDepth < 8) {
					return false;
				}
				
			}  catch(IllegalArgumentException j) { //Sends Error Message if translation does not work
				System.gc();
				JOptionPane.showMessageDialog(new Frame("Error"), "The Sam Reader passed an invalid Chromosome input to"
						+ " the BAM file. The chromosome sent was: " +  chrom + ". The BAM header format will be outputted as \"BAM_Header_Error.txt\"."
								+ "Please edit the \"BAM_Chr_Config\" accordingly. The Program will now exit");
				
				BamChrChanger.writeHeaderErrors(bamChrHeaders.get(i));
				
				System.exit(0);
				
				
				
			}
		}
		return true;
	}

	/*
	 * Extreme novel filter. If there are homs in the population, that means there
	 * are hets, which means other people in the population will have the variant(s)
	 * as hets. This indicates that the population freq of this supposed DN gene is
	 * too high because 2 people with the same DN would have had to get together and
	 * make a child that's hom.
	 */
	public boolean extremePopNovel(int homVarCount) {
		if (homVarCount > 0) {
			return false;
		}
		return true;
	}

	/*
	 * Genotype coverage. If not enough people were genotyped at the position of the
	 * VOI, then our power for whether or not the population frequency filters are
	 * large enough goes down. This is a Poisson problem.
	 */
	public boolean genotypeCoverageEval(int genotypeCount) {
		if (genotypeCount < genoThreshHold) {
			return false;
		}
		return true;
	}

	/*
	 * Stores all the bam files in the familyBam arrayList.
	 */
	public void AddBam(String[] bamlocation, HashMap<String, ArrayList<String>> bamChrMap) throws IOException {
		familyBam = new ArrayList<SamReader>();
		for (int i = 0; i < bamlocation.length; i++) {
			File bamFile = new File(bamlocation[i]);
			// Create a Samreader factory
			SamReaderFactory srf = SamReaderFactory.make();
			srf.validationStringency(ValidationStringency.LENIENT);
			SamReader samR = srf.open(bamFile);
			familyBam.add(samR);
			
			//retrieves BAM chr list for each person
			bamChrHeaders.add(bamChrMap.get(bamlocation[i]));
		}
	}

	public String list2StrStr(List<String> list) {
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

	public String set2Str(Set<String> set) {
		StringBuilder sb = new StringBuilder();
		for (String s : set) {
			sb.append(s);
			sb.append("; ");
		}
		return sb.toString();
	}

	public JFrame canceler(BufferedReader vsData) {
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("FB Confetti Filter");
		jframe.setSize(500, 100);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);

		JLabel jtext = new JLabel("Generating the files with DN's filtered. Use the cancel botton below to abort.",
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
				try {
					vsData.close();
				} catch (IOException e1) {
				}

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
