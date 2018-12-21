package kaylaKode;

import java.awt.Component;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Comparator;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.filechooser.FileNameExtensionFilter;

//this is the most current 
//De novo recessive
//different kind of parental inconsistent: mom uniparental disomy or dad uniparental disomy

/*
 *KaylaKode is a pipeline that enhances upon the original compound het filter developed by Tom,
 * Lukas and Shane. Gracie and Emily have also contributed to its completion along the way.
 * It first employs multi-threads to reduces the burden on memory and CPU usage
 * It completes the analysis by adding the "Watson is mutant" part of the analysis
 * It also adds the Xlink analysis
 *
 * It is more user-friendly, the user input panel supports scroll-down, check and uncheck
 * 
 * The input file is in varsifter format. The header column needs to include "Index", "type", "chr"
 * "LeftFlank", "ref_allele", "var_allele", "HGMDtages", "EPHRED", "EPCPHRED" and "PHRED", all other fields are 
 * from user input
 * 
 * The gene_config file and the varsifter file need to be sorted by chromosome and chrM or any chromosome that 
 * has underscore in them are not accounted for due to currently lack of CADD scores for those lines.
 * 
 * 
 * 9/17/2016
 * @author Kayla Gu
 * 
 * This is the version of KaylaKode that was sent to Tom and me in August 2018. I've reconfigured it for a GUI version.
 * 4/22/2018
 * @author Anchi Wu
 * 
 * 
 * 
 */

public class KaylaKodeNew_GUI2 {

	private class halfHetNode {
		public String[] lineSplit; // the actual linesplit
		public int strength; // the strength based on the type 3 is high, 2 is intermediate, 1 is weak
		public String id; // the id field for this line
		public String ref; // the reference allele on this line
		public String var; // the variant allele on this line
		private double phred; // the phred score for this line
		private double ephred; // the eigen phred score for this line
		private double epcphred;
		public int parent; // 0 is father, 1 is mother
		public ArrayList<Integer> unAffSibs; // tracks if the unaffected siblings are homs (0 for var, 1 for ref) or
												// hets (3)
		
		public int bestStrength; // this stores the bestStrenght - set to 0 at init, overwritten by score method
		public double bestVMMC; // stores the bestVMM CADD score of any pairs - set to 0 at init, overwritten by
								// score method
		public double bestVMME; // stores the bestVMM EIGEN score of any pairs - set to 0 at init, overwritten
								// by score method
		public double bestVMMEPC; // stores the bestVMM EIGEN score of any pairs - set to 0 at init, overwritten
									// by score method
		public int HomRec;// Stores the homozygous recessive value
		public int Pinconsist;// Stores the Mendelian Inconsistency value
		public int Xlinkk; // Stores the Xlink value
		// public int Menddom; // Stores the Mendelian dominant vlaue
		public PriorityQueue<Object[]> pairs; // creates a priority queue showing the other lines that pair with this
												// one, sorted on bestStrength then VMM CADD

		public halfHetNode(String[] inputLine) {
			lineSplit = inputLine;
			id = inputLine[indIndex];
			ref = inputLine[refIndex];
			var = inputLine[varIndex];
			phred = Double.parseDouble(inputLine[ICADD]);
			ephred = Double.parseDouble(inputLine[EphredIndex]);
			epcphred = Double.parseDouble(inputLine[EPCphredIndex]);

			if (homHet(lineSplit[family.get(0)], ref, var) == 3) { // if the father is a het
				parent = 0; // father is het
			} else {
				parent = 1; // mother is het
			}

			String temp = inputLine[typeIndex].toLowerCase().trim();
			if (temp.contains(";")) {
				String[] strArray = temp.split(";");
				strength = 0;
				for (String s : strArray) {
					s = s.toLowerCase().trim();
					if (strong.contains(s)) {
						strength = 3;
					} else if (strong.contains(s) && strength < 2) {
						strength = 2;
					} else if (strength < 1) {
						strength = 1;
					}
				}
			} else if (strong.contains(temp)) {
				strength = 3;
			} else if (inter.contains(temp)) {
				strength = 2;
			} else {
				strength = 1;
			}

			unAffSibs = new ArrayList<Integer>();
			for (int i = 3; i < family.size(); i++) {
				if (!sibsAffected.get(i - 3)) { // if the sibling is unaffected....
					unAffSibs.add(homHet(lineSplit[family.get(i)], ref, var)); // add their homhet score
				}
			}

			bestStrength = 0;
			bestVMMC = 0.0;
			HomRec = 0;
			Pinconsist = 0;
			Xlinkk = 0;

			pairs = new PriorityQueue<Object[]>(10, new Comparator<Object[]>() {
				public int compare(Object[] y, Object[] x) {
					int score1 = (int) Math.signum((int) x[1] - (int) y[1]);

					if (score1 == 0) {
						return (int) Math.signum((double) x[2] - (double) y[2]);
					}
					return score1;
				}
			});
			
		}

		/*
		 * Does the scoring when a halfHet is matched
		 */
		public void score(halfHetNode h) {
			Object[] temp1 = new Object[5]; // goes into this object
			Object[] temp2 = new Object[5]; // goes into the input node

			temp1[0] = h.id; // first index has the id number
			temp2[0] = id;
			temp1[1] = strength * h.strength; // second number as the strength, which is the multiplication of the two
												// individual strengths
			temp2[1] = temp1[1];
			temp1[2] = Math.pow(11,
					(Math.log10(phred + 1) / Math.log10(11)) * (Math.log10(h.phred + 1) / Math.log10(11))); // calculates
																											// the VMM
																											// based on
																											// the two
																											// PRHED
																											// scores
			temp2[2] = temp1[2];
			temp1[3] = Math.pow(11,
					(Math.log10(ephred + 1) / Math.log10(11)) * (Math.log10(h.ephred + 1) / Math.log10(11))); // calculates
																												// the
																												// VMM
																												// based
																												// on
																												// the
																												// two
																												// PRHED
																												// scores
			temp2[3] = temp1[3];
			temp1[4] = Math.pow(11,
					(Math.log10(epcphred + 1) / Math.log10(11)) * (Math.log10(h.epcphred + 1) / Math.log10(11))); // calculates
																													// the
																													// VMM
																													// based
																													// on
																													// the
																													// two
																													// PRHED
																													// scores
			temp2[4] = temp1[4];

			pairs.add(temp1);
			h.pairs.add(temp2);

			if ((int) temp1[1] > bestStrength) { // adds the new strengths to the appropriate nodes if they are higher
				bestStrength = (int) temp1[1];
			}
			if ((int) temp2[1] > h.bestStrength) {
				h.bestStrength = (int) temp2[1];
			}

			if ((double) temp1[2] > bestVMMC) { // adds the new bestVMM to the appropriate nodes if they are higher
				bestVMMC = (double) temp1[2];
			}
			if ((double) temp2[2] > h.bestVMMC) {
				h.bestVMMC = (double) temp2[2];
			}
			if ((double) temp1[3] > bestVMME) { // adds the new bestVMM to the appropriate nodes if they are higher
				bestVMME = (double) temp1[3];
			}
			if ((double) temp2[3] > h.bestVMME) {
				h.bestVMME = (double) temp2[3];
			}
			if ((double) temp1[4] > bestVMMEPC) { // adds the new bestVMM to the appropriate nodes if they are higher
				bestVMMEPC = (double) temp1[4];
			}
			if ((double) temp2[4] > h.bestVMMEPC) {
				h.bestVMMEPC = (double) temp2[4];
			}
		}

		/*
		 * Returns a string to be printed out for this halfHet
		 */
		public String getLine() {
			String output = "";
			long index = Long.parseLong(lineSplit[indIndex].replaceAll("\\D+", ""));
			if (ifR && indexer.containsKey(index)) {
				double[] scorevalue = indexer.get(index);

				if (scorevalue[0] >= bestVMMC && scorevalue[1] >= bestVMME && scorevalue[2] >= bestVMMEPC) {

					return output;
				}

			}

			for (int i = 0; i < naIndex; i++) { // all fields up the family data
				output += lineSplit[i] + "\t";
			}
			output += Integer.toString(bestStrength) + "\t"; // add the beststrength field
			output += Double.toString(bestVMMC) + "\t"; // add the bestvmm CADD field
			output += Double.toString(bestVMME) + "\t"; // add the bestvmm EIGEN field
			output += Double.toString(bestVMMEPC) + "\t"; // add the bestvmm EIGENPC field
			// if (ifM) {
			output += Integer.toString(Pinconsist) + "\t";
			// } // add mendalian inconsistent
			// output += Integer.toString(Menddom) + "\t";// add dominant denovo
			output += Integer.toString(HomRec) + "\t";
			// if (probandgender) {
			output += Integer.toString(Xlinkk) + "\t";
			// }
			// add the bestvmm EIGENPC field
			while (!pairs.isEmpty()) { // take each of pairs from the Pairs queue and create a string, as the
										// MendHetRec field
				output += ((String) pairs.poll()[0]) + ",";
			}

			for (int i = naIndex; i < lineSplit.length; i++) { // add family data
				output += "\t" + lineSplit[i];
			}

			double[] scores = { bestVMMC, bestVMME, bestVMMEPC };
			indexer.put(index, scores);

			return output; // return the string

		}
	}

	/*
	 * Standard chromosomes that's in the config list
	 */

	private class IntervalNode {
		private String name; // name of this gene
		private String transcript;
		private String chr; // the chromosome of this node
		private long start1; // start
		private long start2;
		private long stop1; // stop
		private long stop2;
		private IntervalNode nextNode; // next node in the interval linked list
		private ArrayList<halfHetNode> halfHets;

		/*
		 * Constructor for IntervalNode class
		 */
		public IntervalNode(String[] configLine) {
			name = configLine[geneConfigInd];
			chr = configLine[chrConfigInd];
			start1 = Long.parseLong(configLine[startConfigInd]) - Long.parseLong(configLine[alphConfigInd]);
			start2 = Long.parseLong(configLine[startConfigInd]) + Long.parseLong(configLine[gamaConfigInd]);
			stop1 = Long.parseLong(configLine[endConfigInd]) - Long.parseLong(configLine[deltConfigInd]);
			stop2 = Long.parseLong(configLine[endConfigInd]) + Long.parseLong(configLine[betaConfigInd]);
			transcript = configLine[transConfigInd] + " (" + Long.toString(start1) + ", " + Long.toString(stop2) + ")";
			nextNode = null;
			halfHets = new ArrayList<halfHetNode>();
		}

		public void add(String[] str) {
			halfHets.add(new halfHetNode(str));
		}

		/*
		 * Returns true if position is within boundaries of node.
		 */
		public boolean inInterval(long pos) {
			if ((pos >= start1) && (pos <= start2)) {
				return true;
			} else if ((pos >= stop1) && (pos <= stop2)) {
				return true;
			}

			return false;
		}

		/*
		 * Print all the halfHets in an IntervalNode.
		 */
		public void printLines(PrintWriter writer) {
			if (halfHets.size() != 0) {
				for (halfHetNode h : halfHets) {
					String tmp = h.getLine();
					if (!tmp.equals("")) {
						writer.println(tmp);
					}

				}
			}
		}
	}

	private BufferedReader geneConfigs; // reader for the gene boundary config file
	private int geneConfigInd; // name of the gene
	private int transConfigInd;
	private int chrConfigInd; // chromosome of this interval
	private int startConfigInd; // start transcription position
	private int endConfigInd; // end transcription position
	private int alphConfigInd; // alpha (range in front of a gene) position
	private int betaConfigInd; // beta (range at the end of a gen) position
	private int gamaConfigInd;
	private int deltConfigInd;
	private BufferedReader strengthConfigs; // reader for the strength list config file
	private ArrayList<String> strong; // holds labels that get a strong score (3)
	private ArrayList<String> inter; // holds labels that get an intermediate score (2)
	private boolean ifP;// whether to include population frequency in mendelian frequency
	private boolean ifM;// whetehr to include mendelian inconsistency at all
	private boolean ifunique;
	private boolean ifR;// whether to print unduplicate;
	private int minorfield;
	// private ArrayList<Integer> vartot;

	private File vsFile; // input file that is being filtered
	private String fileName; // destination of filtered file
	private BufferedReader vsData; // used to read the varSifter
	private PrintWriter filterWriter; // used to write out the filter result

	private String Line;
	private String[] curLine; // the current line
	private String[] scroller;
	private List<String> headers; // header line

	private ArrayList<int[]> populations; // indices of the ref allele, var allele, coverage counts
	private ArrayList<int[]> xfilters; // Xlink filters
	private ArrayList<int[]> xfilterc;// Xlink ref and var counts
	private ArrayList<Integer> xfilterExt;
	private int extData;

	private ArrayList<Integer> family; // parents at 0 (father) and 1 (mother), proband at 2, siblings at 3 and onwards

	private ArrayList<Boolean> sibsAffected; // boolean values of whether siblings are affected
	private ArrayList<Boolean> sibsGender;// boolean values of the siblings gender, 1 is male
	private HashMap<Long, double[]> indexer;
	private boolean probandgender;// probandgender==1 is male

	private double rejectionThreshold; // maximum tolerated minor allele frequency, it is halves for xlink
	private int cutcnt2; // top number of hemvar for xlink

	private int indIndex; // index of the Index field
	private int typeIndex; // index of the type field
	private int Itype; // index of locus type
	private int chrIndex; // index of the chromosome field
	private int refIndex; // index of the reference allele
	private int varIndex; // index of the variant allele
	private int posIndex; // index of the left_flank
	private int ICADD; // index of the phred field
	private int EphredIndex; // index of the phred field
	private int EPCphredIndex; // index of the phred field
	private int hgmdIndex; // index of the HGMDtags field
	private int naIndex; // index of the start of the family data
	private boolean chrM;
	private int polyphen;
	private List<IntervalNode> intervals;
	private int threads;
	private int familytotal;
	
	private String strengthPath;
	private String genePath;
	
	private File output;
	
	public KaylaKodeNew_GUI2(File vsFile, String dest, ArrayList<Integer> fam, ArrayList<Boolean> sibsAffect, ArrayList<Boolean> sibGend, boolean proGend, boolean[] popBools, ArrayList<int[]> populs, ArrayList<int[]> xFiltS,
			ArrayList<int[]> xFiltC, double rejectThresh, int thredz, int cutCNT, String strengthPath, String genePath, int minor, ArrayList<Integer> xFiltExt, int extData) throws IOException, InterruptedException {
		
		this.vsFile = vsFile;
		this.fileName = dest;
		this.family = fam;
		this.sibsAffected = sibsAffect;
		this.sibsGender = sibGend;
		this.probandgender = proGend;
		this.ifP = popBools[0];
		this.ifM = popBools[1];
		this.ifunique = popBools[2];
		this.ifR = popBools[3];
		this.chrM = popBools[4];
		this.minorfield = minor;
		this.populations = populs;
		this.xfilters = xFiltS;
		this.xfilterc = xFiltC;
		this.xfilterExt = xFiltExt;
		this.extData = extData;
		this.rejectionThreshold = rejectThresh;
		this.threads = thredz;
		this.cutcnt2 = cutCNT;
		this.strengthPath = strengthPath;
		this.genePath = genePath;
		
		
		
		initializer();
		Core();
	}
	
	
	/*
	 * User input for VarSifter (VS) file and output file name.
	 *
	
	public void uI() {
		//vsFile = getVSFile();
		//fileName = getDestination();
	}
	 */
	
	public File getOutput() {
		
		return output;
		
	}
	
	
	/// We are only counting for chr1 to chr22 chrX chrY chrUn
	/*
	 * Starts thread, and reads in gene configuration and strength type
	 * configuration files.
	 */
	public void Core() throws IOException, InterruptedException {
		File f = new File(System.getProperty("java.class.path")); // finds the folder that the JAR is stored in
		File dir = f.getAbsoluteFile().getParentFile();
		String path = dir.toString(); // the path to the folder

		try {
			strengthConfigs = new BufferedReader(new FileReader(strengthPath));
			// strengthConfigs = new BufferedReader(new FileReader(path +
			// "/Strength_Config.txt"));
			// strengthConfigs = new BufferedReader(new FileReader("Strength_Config.txt"));

		} catch (FileNotFoundException e) {
			JOptionPane.showMessageDialog(new JPanel(),
					"The strength config file could not be found. Please ensure the file is in the same directory as the .jar.");
			vsData.close();
			filterWriter.close();
			System.gc();
			System.exit(0);
		}

		strong = new ArrayList<String>(); // an ArrayList of types that get a strong score
		inter = new ArrayList<String>(); // an ArrayList of types that get a weak score
		ArrayList<String> tempArray = strong; // set the tempArray to point to the strong array
		while (strengthConfigs.ready()) {
			String line = strengthConfigs.readLine();
			if (line.substring(0, 2).equals("##")) { // if it's a header line, ignore
				continue;
			}

			if (line.equals("#Strong")) { // if the line says strong, set the pointer array to Strong and continue
				tempArray = strong;
				continue;
			}

			if (line.equals("#Intermediate")) { // if the line is the intermediate header, reset the pointer array and
												// continue
				tempArray = inter;
				continue;
			}
			tempArray.add(line.toLowerCase().trim()); // add the string to the appropriate arraylist
		}
		strengthConfigs.close();

		try {
			// geneConfigs = new BufferedReader(new
			// FileReader("Gene_Config_Splits_USE.txt"));
			geneConfigs = new BufferedReader(new FileReader(genePath));
			// geneConfigs = new BufferedReader(new FileReader(path +
			// "/Gene_Config_Splits_USE.txt")); // finds the config file; must be stored in
			// the
			// same folder as the JAR
		} catch (FileNotFoundException e) {
			JOptionPane.showMessageDialog(new JPanel(),
					"The gene config file could not be found. Please ensure the file is in the same directory as the .jar.");
			vsData.close();
			filterWriter.close();
			System.gc();
			System.exit(0);
		}

		geneConfigs.readLine(); // gets rid of the info line
		String[] temp = geneConfigs.readLine().split("\t");
		List<String> configHeaders = java.util.Arrays.asList(temp);

		// get the headers within the config file
		geneConfigInd = configHeaders.indexOf("#Name");
		transConfigInd = configHeaders.indexOf("Transcript");
		chrConfigInd = configHeaders.indexOf("Chromosome");
		startConfigInd = configHeaders.indexOf("Transcription_Start");
		endConfigInd = configHeaders.indexOf("Transcription_End");
		alphConfigInd = configHeaders.indexOf("Alpha");
		betaConfigInd = configHeaders.indexOf("Beta");
		gamaConfigInd = configHeaders.indexOf("Gamma");
		deltConfigInd = configHeaders.indexOf("Delta");

		intervals = Collections.synchronizedList(new ArrayList<IntervalNode>()); // make it threadsafe
		intervals.add(new IntervalNode(geneConfigs.readLine().split("\t"))); // adds the first interval as an interval
																				// node (see above)
		while (geneConfigs.ready()) { // adds the rest of the config file
			intervals.add(new IntervalNode(geneConfigs.readLine().split("\t")));
			intervals.get(intervals.size() - 2).nextNode = intervals.get(intervals.size() - 1);
		}
		geneConfigs.close(); // close the reader

		JFrame jframe = canceler();
		jframe.setLocationRelativeTo(null);
		jframe.setVisible(true);

		String[] Split;
		String lineChr;
		long linePos;
		Line = vsData.readLine();

		// go to the first node
		IntervalNode curNode = intervals.get(0); // the current node
		IntervalNode slider;
		ExecutorService execService = Executors.newFixedThreadPool(threads);
		
		while (Line != null) {
			Split = Line.split("\t");
			lineChr = Split[chrIndex];
			linePos = Long.parseLong(Split[posIndex]);
			
			//If Loctype column does not exist
			String loctype = "";
			
			if(Itype != -1) {
				loctype = Split[Itype];
			} else if(typeIndex != -1) {
				
				if((Split[typeIndex].contains("synonymous") || Split[typeIndex].contains("SYNONYMOUS")) || (Split[typeIndex].contains("stopgain") || Split[typeIndex].contains("STOPGAIN"))
						|| (Split[typeIndex].contains("frameshift") || Split[typeIndex].contains("FRAMESHIFT"))) {
					
					loctype = "exonic";
					
				}
				
			} else {
				
				JOptionPane.showMessageDialog(new Frame("Error"),
						"There is no 'loc_type' or 'type' header in the input file. Please close and try again.");
				vsData.close();
				filterWriter.close();
				System.exit(0);
				
			}
			
			double caddscore = Double.parseDouble(Split[ICADD]);
			int caddcutoff = 2;

			/// If variant is exonic, set default CADD cutoff is 7.
			if (loctype.contains("exonic") && !loctype.contains("ncRNA")) {
				caddcutoff = 7;
			}

			if (!lineChr.contains("_") && !lineChr.contains("Y")) {
				if (halfHetDetection(Split)) { // if this is a halfHet
					// System.out.println("lineChr: " + lineChr);
					// System.out.println("name: " + curNode.name + "; chr: " + curNode.chr);
					while (!curNode.chr.equals(lineChr)) { /// halfHets must be on same chromosome
						curNode = curNode.nextNode;
					}
					slider = curNode.nextNode;
					if (curNode.inInterval(linePos)) {
						curNode.add(Split);
					}
					while (slider != null && (linePos >= slider.start1) && slider.chr.equals(lineChr)) {
						if (slider.inInterval(linePos)) {
							slider.add(Split);
						}
						slider = slider.nextNode;
					}
				} else {
					if (caddscore >= caddcutoff) {
						nonPar(Split, lineChr, linePos); /// Evaluate and print line if non-compound het.
					}
					// execService.execute(new Lineprocessor(Split,lineChr,linePos));
				}
			}
			Line = vsData.readLine();
		}

		for (IntervalNode n : intervals) {
			execService.execute(new Cmpd(n.halfHets));
		}

		execService.shutdown();

		for (IntervalNode n : intervals) {
			n.printLines(filterWriter);
		}
		// line is at the end
		filterWriter.close();
		vsData.close();
		jframe.dispose();
		
	}

	public class Cmpd extends Thread {
		private ArrayList<halfHetNode> interv;

		public Cmpd(ArrayList<halfHetNode> interv) {
			this.interv = interv;
		}

		public void run() {
			interv = findCHetPairs(interv);
		}
	}

	/*
	 * only run X-linked filter on chromosome X, nonautosomal regions
	 * 
	 * @param n - position in the genome
	 */
	public boolean autosome(Long n) {

		if (n < 2699520) {
			return false;
		}
		if (n > 154931043) {
			return false;
		}
		if (88400000 < n && 92000000 > n) {
			return false;
		} else {
			return true;
		}

	}

	/*
	 * Boolean. Returns true if variant is mitochondrial.
	 */
	public boolean chrMDetection(String[] line) {
		if (line[polyphen].equalsIgnoreCase("possibly damaging")
				|| line[polyphen].equalsIgnoreCase("probably damaging")) {

			String ref;
			String var;

			if (Integer.parseInt(line[minorfield]) != 1) {
				ref = line[refIndex];
				var = line[varIndex];
			} else {// if all gahl reference is minor, switch ref and var
				ref = line[varIndex];
				var = line[refIndex];
			}
			int geno = homHet(line[family.get(2)], ref, var);
			if (!(geno == 1)) { // if geno is not ref return true
				return Mcheck(line);
			}
		}
		return false;
	}

	/*
	 * Lineprocessor takes in the config of particular chromosome it runs through
	 * varsifter file for the corresponding chromosome and runs several filters on
	 * them it doesn't not account for chrM or any chromosome with "_" in them
	 */
	// public class Lineprocessor extends Thread {
	// private String[] liny;
	// private String chr;
	// private long pos;
	//
	// public Lineprocessor(String[] line, String chr, long pos) {
	// this.liny = line;
	// this.chr = chr;
	// this.pos = pos;
	// }
	//
	// public void run() {
	//
	// if (chr.equals("chrM") && chrM) {
	// if (chrMDetection(liny)) {
	// filterWriter.println(PrintchrM(liny));
	// }
	// }
	//
	// else if (probandgender && chr.equals("chrX") && autosome(pos + 1)) {
	// if (Xlinkcheck(liny)) {
	// filterWriter.println(Xlinkprint(liny));
	// }
	// } else if (HomRecDetection(liny)) {
	// filterWriter.println(PrintHomo(liny));
	// } else if (CnCtest(liny)) {
	// filterWriter.println(PrintCnC(liny));
	// } else {
	// Integer[] TF = Pinconsist(liny, ifP);
	// if (ifM & TF[0] == 1) {
	// filterWriter.println(PrintInconsist(liny, TF));
	// } else if ((!ifM) & TF[1] == 1) {
	// filterWriter.println(PrintInconsist(liny, TF));
	// }
	//
	// }
	//
	// }
	//
	// }

	public void initializer() throws IOException {
		try {
			// vsData = new BufferedReader(new FileReader("LYSTALL.vs"));//input VarSifter
			// file - file provided by user
			
			vsData = new BufferedReader(new FileReader(vsFile));
			
			output = new File(fileName + ".vs");
			
			filterWriter = new PrintWriter(output);

		} catch (FileNotFoundException error) {
			JOptionPane.showMessageDialog(new Frame("Error"),
					"The VarSifter file could not be found, possibly because it is open elsewhere. Please close and try again.");
			System.exit(0);
		}

		Line = vsData.readLine();
		curLine = Line.split("\t");
		headers = java.util.Arrays.asList(curLine);
		indexer = new HashMap<Long, double[]>();

		// finds the headers within the VarSifter file being filtered
		// need Index, type, chr, LeftFlank, ref_allele, var_allele, HGMDtags, PHRED
		// columns
		// running through the annotation pipeline will ensure all necessary columns are
		// present and correct
		indIndex = headers.indexOf("Index");
		typeIndex = headers.indexOf("type");
		chrIndex = headers.indexOf("Chr");
		posIndex = headers.indexOf("LeftFlank");
		refIndex = headers.indexOf("ref_allele");
		varIndex = headers.indexOf("var_allele");
		Itype = headers.indexOf("loc_type");
		polyphen = headers.indexOf("Prediction");
		if (polyphen == -1) {
			polyphen = 0;
		}
		ICADD = headers.indexOf("PHRED");
		hgmdIndex = headers.indexOf("HGMDtags");
		EphredIndex = headers.indexOf("EPHRED");
		EPCphredIndex = headers.indexOf("EPCPHRED");
		
		genotypeIndices(curLine, curLine.length);
		
		// rejectionThreshold = getThreshold(); //prompt the user to enter the rejection
		// threshold
		// rejectionThreshold=0.04;

		String output = ""; // create the new header line

		for (int i = 0; i < naIndex; i++) { // add up to the family data
			output += curLine[i] + "\t";
		}
		
		//scroller = output.split("\t");
		//getPedigree(); // prompt the user to enter the pedigree
		//getPopulations();

		output += "Best_Strength" + "\t"; // add the four new fields used to view comphets
		output += "Best_VMMC" + "\t";
		output += "Best_VMME" + "\t";
		output += "Best_VMMEPC" + "\t";
		// if (ifM) {
		output += "MendInconsis" + "\t";
		// }
		// output += "DeNovoDom" + "\t";
		output += "MendHomRec" + "\t";
		// if (probandgender) {
		output += "Xlink" + "\t";
		// }
		output += "MendHetRec";

		for (int i = naIndex; i < curLine.length; i++) { // add the family data curLine is the current line container
			output += "\t" + curLine[i];
		}

		filterWriter.println(output);

	}
	

	//gets the start of the metadata (naPos) and the distance between each of the genotypes
	public void genotypeIndices(String[] lineSplit, int columns) {
		int counter = 0;
		int tempVal = 0;
		while (tempVal <2 && counter < columns) {								//should capture initial genotype index
			int val = lineSplit[counter].length();								//and distance between genotype indices - if the headder is at least 3 char long. Otherwise, will bypass.
			
			Pattern r = Pattern.compile("\\.NA$");
			Matcher m = r.matcher(lineSplit[counter]);
			
			
			if (val < 3 ) {
				counter++;														// don't go looking for sample data that is 3 characters long if the headder is less than 3 characters
			}																	// this is for idiots that make headers with less than 3 characters - that is you Lukas!
			//System.out.println(val);
			else {
				if (m.find()) {		
					if (tempVal == 0) {
						naIndex = counter;
					}
					else {
						//naDistanceApart = counter - naIndex;
					}
				tempVal++;			
				}
				counter++;
			}
		}
	}
	
	
	/*
	 * Print line for non-compound het variant to output.
	 */
	public void nonPar(String[] liny, String chr, long pos) {
		if (chr.equals("chrM")) {
			if (chrMDetection(liny)) { /// Mitochondrial
				filterWriter.println(PrintchrM(liny));
			}
		} else if (HomRecDetection(liny)) { /// Homozygous recessve
			filterWriter.println(PrintHomo(liny));
		} else {
			/// TODO: This is where we'd put in the filter for X-linked DN in females
			if (probandgender) { // if proband is male
				if (chr.equals("chrX") && autosome(pos + 1)
						&& MendelianState.hemy(liny[family.get(2)], liny[refIndex], liny[varIndex])) {
					if (Xlinkcheck(liny)) { /// X-linked
						filterWriter.println(Xlinkprint(liny));
					}
				} else {
					int inconsistent = Pinconsist(liny);
					if (inconsistent == 2 || inconsistent == 3 || inconsistent == 4) {
						filterWriter.println(PrintInconsist(liny, inconsistent));	
					}
				}
			} else { // if female
				int inconsistent = Pinconsist(liny);
				if (inconsistent == 2 || inconsistent == 3 || inconsistent == 4) {
					filterWriter.println(PrintInconsist(liny, inconsistent));
				}
			}
		}
		// else if (chr.equals("chrX") && autosome(pos + 1)
		// && MendelianState.hemy(liny[family.get(2)], liny[refIndex], liny[varIndex]))
		// {
		// if (Xlinkcheck(liny)) { /// X-linked
		// filterWriter.println(Xlinkprint(liny));
		// }
		// }

		// else { /// Mendelian inconsistent
		// int inconsistent = Pinconsist(liny);
		// if (inconsistent == 2 || inconsistent == 3 || inconsistent == 4) {
		// filterWriter.println(PrintInconsist(liny, inconsistent));
		// }
		// }
	}

	public boolean CnCtest(String[] liny) {
		int deletion = 0;

		if (Integer.parseInt(liny[family.get(2) + 2]) == 0 && Integer.parseInt(liny[family.get(2) + 1]) == 0) {
			if (Integer.parseInt(liny[family.get(0) + 2]) >= 5 && Integer.parseInt(liny[family.get(1) + 2]) >= 5) {

				for (int i = naIndex; i < liny.length; i += 3) {
					if (Integer.parseInt(liny[i + 2]) < 5 && Integer.parseInt(liny[i + 1]) < 5) {
						deletion++;
					}
				}

				if (deletion < 3) {
					return true;
				}
			}
		}
		return false;
	}

	// tests if the fields entered by the user in the getPedigree function are
	// actually in the header
	public boolean inHeader(String word) {
		if (headers.indexOf(word) == -1) {
			return false;
		}
		return true;
	}

	public JFrame canceler() {
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("KaylaKode");
		jframe.setSize(500, 100);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);

		JLabel jtext = new JLabel("Generating the Kode output file. Use the cancel botton below to abort.",
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
				filterWriter.close();
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

	// this will test raw data for the population criteria
	// returns false if the line is to be discarded and true if not
	public boolean frequencyChecker(String alleleCount, String coverage, Double nn) {
		//System.out.println("alleleCount: " + alleleCount);
		int altCount = Integer.parseInt(alleleCount);
		int totCount = Integer.parseInt(coverage);

		if (totCount > 700) { // if the coverage is this high, we can take the measured frequency at face
								// value
			double rate = ((double) altCount) / ((double) totCount);
			double thred = 1.2 * nn;
			return rate < thred; // always hedge above the threshold by a factor of 1.2

		}

		// for smaller populations, the full calculation must be performed
		// for more on the mathematics of this, consult the manuel
		double fact = 1;
		double prob = Math.pow(2.718, -altCount) * Math.pow(altCount, 0) / fact; // adds the 0 case
		for (int i = 1; i <= Math.ceil(totCount * nn); i++) { // loops through cases 1->rejectionthreshold
			fact = fact * i;
			prob += Math.pow(2.718, -altCount) * Math.pow(altCount, i) / fact;
		}
		if (prob < 0.05) { // if the probability that the value is below the rejection threshold is less
							// then 5%, reject
			return false;
		}
		return true;
	}

	/*
	 * Detects halfHets, checking both Mendelian consistency and frequency data.
	 * Detects if a line qualifies as a halfHet
	 */
	public boolean halfHetDetection(String[] line) {
		String ref;
		String var;

		if (Integer.parseInt(line[minorfield]) != 1) {
			ref = line[refIndex];
			var = line[varIndex];
		} else {// if all gahl reference is minor, switch ref and var
			ref = line[varIndex];
			var = line[refIndex];
		}

		int parentVal = homHet(line[family.get(0)], ref, var) + homHet(line[family.get(1)], ref, var); // stores the sum
																										// of the
																										// parents
																										// genotypes

		if ((parentVal >= 3) && (parentVal < 6)) { // only one parent can be het, so either het(3) + homref(1) = 4 or
													// het(3) + homvar(0) = 3 will work
			if (homHet(line[family.get(2)], ref, var) == 3) { // proband must be het
				for (int i = 3; i < family.size(); i++) { // all affected sibs must be het, unaffected sibs cannot be
															// homvar
					int sibgeno = homHet(line[family.get(i)], ref, var);
					if (sibsAffected.get(i - 3) && sibgeno != 3 && sibgeno != -5) { // if sib is affected and not het,
																					// reject
						return false; // fail
					} else if (!sibsAffected.get(i - 3) && homHet(line[family.get(i)], ref, var) == 0) { // if sib is
																											// unaffected
																											// and
																											// homvar,
																											// reject
						return false; // fail
					}
				}
				return Fcheck(line);
			}
		}
		return false; // default to false
	}

	// This section defines homRec variants
	public boolean HomRecDetection(String[] line) {
		String ref;
		String var;

		if (Integer.parseInt(line[minorfield]) != 1) {
			ref = line[refIndex];
			var = line[varIndex];
		} else {// if all gahl reference is minor, switch ref and var
			ref = line[varIndex];
			var = line[refIndex];
		}
		int parentVal = homHet(line[family.get(0)], ref, var) + homHet(line[family.get(1)], ref, var); // stores the sum
																										// of the
																										// parents
																										// genotypes
		if (parentVal == 6) { // both parents are het
			if (homHet(line[family.get(2)], ref, var) == 0) { // proband cannot be het
				for (int i = 3; i < family.size(); i++) { // all affected sibs must be het, unaffected sibs cannot be
															// homvar
					int sibgeno = homHet(line[family.get(i)], ref, var);
					if (sibsAffected.get(i - 3) && sibgeno != 0 && sibgeno != -5) { // if sib is affected and not homvar
																					// nor NA, reject
						return false; // fail
					} else if (!sibsAffected.get(i - 3) && sibgeno == 0) { // if sib is unaffected and homvar, reject
						return false; // fail
					}
				}
				return Fcheck(line);
			}
		}
		return false; // default to false
	}

	// Xlink filter
	public boolean Xlinkcheck(String[] line) {
		String ref;
		String var;

		if (Integer.parseInt(line[minorfield]) != 1) {
			ref = line[refIndex];
			var = line[varIndex];
		} else {// if all gahl reference is minor, switch ref and var
			ref = line[varIndex];
			var = line[refIndex];
		}

		// mother is het and father is not hemvar and proband is hemref
		// AW: mother is het and father is not hemvar and proband is HEMVAR?
		if (homHet(line[family.get(1)], ref, var) == 3 && homHet(line[family.get(0)], ref, var) != 0
				&& homHet(line[family.get(2)], ref, var) == 0) {
			for (int i = 3; i < family.size(); i++) { // all male affected is hemvar
				int sibgeno = homHet(line[family.get(i)], ref, var);
				// if gender is male
				if (sibsGender.get(i - 3)) {
					if (sibsAffected.get(i - 3) && sibgeno != 0 && sibgeno != -5) { // if sib is affected and not NA or
																					// homvar, reject
						return false; // fail
					} else if (!sibsAffected.get(i - 3) && sibgeno == 0) { // if sib is unaffected and homvar, reject
						return false; // fail
					}
				}
			}
			return Xcheck(line);
		}
		return false;
	}

	/*
	 * public boolean Uniquevariant(String[] line){ String ref; String var;
	 * 
	 * if(Integer.parseInt(line[minorfield])!=1){ ref = line[refIndex]; var =
	 * line[varIndex]; } else{//if all gahl reference is minor, switch ref and var
	 * ref=line[varIndex]; var=line[refIndex]; } int varcount=0; //mother is het and
	 * father is not hemvar and proband is hemref
	 * 
	 * if( Integer.parseInt(line[family.get(2)+1])<15 ||homHet(line[family.get(2)],
	 * ref, var)==1 ||homHet(line[family.get(2)], ref, var)==-5){ return false; }
	 * for (int i = 0; i < family.size(); i++) { //all male affected is hemvar int
	 * gen=homHet(line[family.get(i)], ref, var);
	 * 
	 * if(gen==3){ varcount++; } if(gen==0){ varcount+=2; } }
	 * 
	 * return Integer.parseInt(line[vartot.get(0)])<(varcount+3); }
	 */

	/*
	 * 
	 * /////////////////////////////////////////////////////////////// Check for
	 * /////////////////////////////////////////////////////////////// Mendalian
	 * /////////////////////////////////////////////////////////////// inconsistency
	 * // returns 1 in the first element in a string[] if any of the following is
	 * true: // Both parents homref(1+1) the proband is not homref (!1, !-5) // Both
	 * parents homvar (0+0) the proband is not homvar (!0, !-5) // One parent
	 * homvar, one parent homref (0+1), the proband is not het (!3, !-5) // One
	 * parent is het, one homref (3+1), the proband is homvar (0) // One parent is
	 * het, one homvar (3+0), the proband is homref (1) // Double check if it
	 * doesn't need to pass population check
	 * 
	 * // return 1 in the second element in the a string[] if de novo dominate is
	 * true: // pass the population frequency check // parents are homref(1+1), the
	 * affected is het (3) The unaffected is homref 1 // or NA // parents are
	 * homvar(0+0), the affected is het (3) the unaffected is homvar 0 // or NA
	 * 
	 * public Integer[] Pinconsist(String[] line, boolean ifP) { Integer[] Pboolean
	 * = { 0, 0 }; String ref; String var;
	 * 
	 * if (Integer.parseInt(line[minorfield]) != 1) { ref = line[refIndex]; var =
	 * line[varIndex]; } else {// if all gahl reference is minor, switch ref and var
	 * ref = line[varIndex]; var = line[refIndex]; }
	 * 
	 * int parentVal = homHet(line[family.get(0)], ref, var) +
	 * homHet(line[family.get(1)], ref, var); // stores the sum // of the // parents
	 * // genotypes int probandVal = homHet(line[family.get(2)], ref, var);// store
	 * the proband genotype
	 * 
	 * if (parentVal == 2) { if (probandVal != 1 && probandVal != -5) { Pboolean[0]
	 * = 1; // check denovo dominant if (probandVal == 3) {// proband is a het int
	 * result2 = 0; for (int ii = 3; ii < family.size(); ii++) { // all affected
	 * sibs must be het, unaffected sibs has // to be homref int sibgeno =
	 * homHet(line[family.get(ii)], ref, var); if (sibsAffected.get(ii - 3) &&
	 * sibgeno != 3 && sibgeno != -5) { // condition for not deno result2++; } if
	 * (!sibsAffected.get(ii - 3) && sibgeno != 1 && sibgeno != -5) {// condition
	 * for not deno result2++; } }
	 * 
	 * if (result2 == 0) { Pboolean[1] = 1; }
	 * 
	 * } } else { return Pboolean; } }
	 * 
	 * if (parentVal == 0) { if (probandVal != 0 && probandVal != -5) { Pboolean[0]
	 * = 1; // then check for de novo if (probandVal == 3) {// proband is a het int
	 * result2 = 0; for (int ii = 3; ii < family.size(); ii++) { // all affected
	 * sibs must be het, unaffected sibs has // to be homref int sibgeno =
	 * homHet(line[family.get(ii)], ref, var); if (sibsAffected.get(ii - 3) &&
	 * sibgeno != 3 && sibgeno != -5) { result2++; }
	 * 
	 * if (!sibsAffected.get(ii - 3) && sibgeno != 0 && sibgeno != -5) { result2++;
	 * } } if (result2 == 0) { Pboolean[1] = 1; } } } else { return Pboolean; } } //
	 * One parent is het, one homref (3+1), the proband is homvar (0) if (parentVal
	 * == 4) { if (probandVal == 0) { Pboolean[0] = 1; } else { return Pboolean; } }
	 * // One parent homvar, one parent homref (0+1), children are not het (!3, !-5)
	 * if (parentVal == 1) { if (probandVal != 3 && probandVal != 5) { Pboolean[0] =
	 * 1; } else { return Pboolean; } } // One parent is het, one homvar (3+0), the
	 * proband is homref (1) if (parentVal == 3) { if (probandVal == 1) {
	 * Pboolean[0] = 1; } else { return Pboolean; } } if (ifP && !Fcheck(line)) {
	 * Pboolean[0] = 0; }
	 * 
	 * if ((!ifP) && !Fcheck(line)) { Pboolean[1] = 0; }
	 * 
	 * return Pboolean; }
	 */

	/*
	 * Determine type of Mendelian inconsistency of the variant. Returns a value
	 * based on apparent model: 0 - consistent 1 - inconsistent 2 - de novo 3 -
	 * hemizygous 4 - call no call (extreme novel deleted exon) 5 - non-novel
	 * deletion
	 * 
	 * Returns 1 in the first element in a string[] if any of the following is true:
	 * Both parents homref(1+1) the proband is not homref (!1, !-5) Both parents
	 * homvar (0+0) the proband is not homvar (!0, !-5) One parent homvar, one
	 * parent homref (0+1), the proband is not het (!3, !-5) One parent is het, one
	 * homref (3+1), the proband is homvar (0) One parent is het, one homvar (3+0),
	 * the proband is homref (1) Double check if it doesn't need to pass population
	 * check
	 * 
	 * Return 1 in the second element in the a string[] if de novo dominant is true:
	 * Pass the population frequency check parents are homref(1+1), the affected is
	 * het (3) The unaffected is homref 1 or NA parents are homvar(0+0), the
	 * affected is het (3) the unaffected is homvar 0 or NA
	 */
	public Integer Pinconsist(String[] line) {

		int inconsistent = 0;

		String ref;
		String var;

		if (Integer.parseInt(line[minorfield]) != 1) {
			ref = line[refIndex];
			var = line[varIndex];

		} else { // if all gahl reference is minor, switch ref and var
			ref = line[varIndex];
			var = line[refIndex];
		}

		boolean mut = true;
		if (ref.length() == 1 && var.length() == 1) {
			mut = false;
		}
		int probandVal = homHet(line[family.get(2)], ref, var);// store the proband genotype

		/// Proband is hemizygous or homozygous for the major allele
		if (probandVal != 1) {

			int inconsistentfirstcheck = MendelianState.trueconsist(line[family.get(0)], line[family.get(1)],
					line[family.get(2)], mut);

			if (inconsistentfirstcheck == 1) {
				inconsistent = 1;
				/// Store the sum of the parents' genotypes
				int parentVal = homHet(line[family.get(0)], ref, var) + homHet(line[family.get(1)], ref, var);

				// Check for de novo
				/// TODO
				if (probandVal == 3 && (parentVal == 2 || parentVal == 0)) {
					inconsistent = 2;
					for (int ii = 3; ii < family.size(); ii++) { // all affected sibs must be het, unaffected sibs has
						// to be homref
						int sibgeno = homHet(line[family.get(ii)], ref, var);
						if (sibsAffected.get(ii - 3) && sibgeno == (parentVal / 2)) {
							inconsistent = 1;
							break;
						}
						if (!sibsAffected.get(ii - 3) && sibgeno != (parentVal / 2) && sibgeno != -5) {// condition for
							// not deno
							inconsistent = 1;
							break;
						}
					}
				} else if (probandVal + parentVal == 4 && probandVal != 3) {
					inconsistent = 3;
					for (int ii = 3; ii < family.size(); ii++) { // all affected sibs must be het, unaffected sibs has
						// to be homref
						int sibgeno = homHet(line[family.get(ii)], ref, var);
						if (sibsAffected.get(ii - 3) && sibgeno != probandVal && sibgeno != -5) {
							inconsistent = 1;
							break;
						}
						if (!sibsAffected.get(ii - 3) && sibgeno == probandVal) {// condition for not deno
							inconsistent = 1;
							break;
						}
					}
				}
			} else if (inconsistentfirstcheck == 2) {
				if (CnCtest(line)) {
					inconsistent = 4;
				} else {
					inconsistent = 5;
				}
			}
		}

		if (inconsistent != 0 && inconsistent != 4) {
			if (!inconsistFcheck(line)) {
				inconsistent = 0;
			}
		}

		return inconsistent;
	}

	/*
	 * Determine if variant is X-linked.
	 */
	public boolean Xcheck(String[] line) {
		if (line[hgmdIndex].equals("DM")) { // if it's HGMD DM, return true regardless of population statistics
			return true;
		}
		// if all gahl reference is minor, filter all the reference fields
		if (Integer.parseInt(line[minorfield]) == 1) {
			for (int[] a : xfilters) {
				if ((a.length == 3) && (!frequencyChecker(line[a[0]], line[a[2]], rejectionThreshold / 4))) {
					/// if there is raw data, use the frequencychecker method (above)
					return false; // fail
				} else if ((a.length == 2) && Boolean.parseBoolean(line[a[0]])) {
					/// if the length is one, it is a filter field - just parse the boolean.
					return false; // fail
				}
			}
			for (int[] aa : xfilterc) {
				if (Integer.parseInt(line[aa[0]]) > cutcnt2) {
					return false;
				}
			}
		}
		/// if all gahl reference is not minor, filter all the variance field
		else {
			for (int[] aa : xfilterc) {
				if (Integer.parseInt(line[aa[1]]) > cutcnt2) {
					return false;
				}
			}
			for (int[] a : xfilters) {
				if ((a.length == 3) && (!frequencyChecker(line[a[1]], line[a[2]], rejectionThreshold / 4))) {
					/// if there is raw data, use the frequencychecker method (above)
					return false; // fail
				} else if ((a.length == 2) && Boolean.parseBoolean(line[a[1]])) {
					/// if the length is one, it is a filter field - just parse the boolean.
					return false; // fail
				}
			}
		}
		if(extData == 0) {
			
			//Checks External Databases for a count greater than the cutoff count, if so, then it rejects variants
			for (int a : xfilterExt) {
				if (Integer.parseInt(line[a]) > cutcnt2) {
					return false;
				}
			}
		}
		
		return true;
	}

	/*
	 * Returns true if variant is not too frequent, globally (not just raw data).
	 * However, if HGMD tag or PolyPhen tags indicate that the variant is
	 * pathogenic, return true regardless of frequency.
	 */
	public boolean Fcheck(String[] line) {
		if (line[hgmdIndex].equals("DM") || line[polyphen].equals("probably damaging")
				|| line[polyphen].equals("possibly damaging") || Double.parseDouble(line[ICADD]) >= 15
				|| Double.parseDouble(line[EphredIndex]) >= 20 || Double.parseDouble(line[EPCphredIndex]) >= 20) {

			if (Integer.parseInt(line[minorfield]) != 1) {
				for (int[] a : populations) {
					if ((a.length == 3) && (!frequencyChecker(line[a[1]], line[a[2]], 0.1))) {
						/// if there is raw data, use the frequencychecker method (above)
						return false; // fail
					}
				}

			} else {
				for (int[] a : populations) {
					if ((a.length == 3) && (!frequencyChecker(line[a[0]], line[a[2]], 0.1))) {
						/// if there is raw data, use the frequencychecker method (above)
						return false; // fail
					}
				}
			}
			// if it's HGMD DM, return true regardless of population statistics
			return true;
		}
		// if all gahl reference is not minor
		// filter all the variance fields

		if (Integer.parseInt(line[minorfield]) != 1) {
			for (int[] a : populations) {
				if ((a.length == 3) && (!frequencyChecker(line[a[1]], line[a[2]], rejectionThreshold))) {
					/// if there is raw data, use the frequencychecker method (above)
					return false; // fail

				} else if ((a.length == 2) && Boolean.parseBoolean(line[a[1]])) {
					/// if the length is one, it is a filter field - just parse the boolean.
					return false; // fail
				}
			}
		}
		// if all gahl reference is minor
		else {
			for (int[] a : populations) {
				if ((a.length == 3) && (!frequencyChecker(line[a[0]], line[a[2]], rejectionThreshold))) {
					/// if there is raw data, use the frequencychecker method (above)
					return false; // fail
				}
				// else if ((a.length == 2) && Boolean.parseBoolean(line[a[0]])) { //if the
				// length is one, it is a filter field - just parse the boolean
				// return false; //fail
				// }
			}
		}
		return true;
	}

	/*
	 * Returns true if Mendelian inconsistent variant is not too frequent, globally
	 * (not just raw data). However, if HGMD tag or PolyPhen tags indicate that the
	 * variant is pathogenic, return true regardless of frequency.
	 */
	public boolean inconsistFcheck(String[] line) {
		if (line[hgmdIndex].equals("DM") || line[polyphen].equals("probably damaging")
				|| line[polyphen].equals("possibly damaging") || Double.parseDouble(line[ICADD]) >= 15
				|| Double.parseDouble(line[EphredIndex]) >= 20 || Double.parseDouble(line[EPCphredIndex]) >= 20) {

			if (Integer.parseInt(line[minorfield]) != 1) {
				for (int[] a : populations) {
					if ((a.length == 3) && (!frequencyChecker(line[a[1]], line[a[2]], 0.1))) {
						/// if there is raw data, use the frequencychecker method (above)
						return false; // fail
					}
				}
			} else {
				for (int[] a : populations) {
					if ((a.length == 3) && (!frequencyChecker(line[a[0]], line[a[2]], 0.1))) {
						/// if there is raw data, use the frequencychecker method (above)
						return false; // fail
					}
				}
			}
			/// If it's HGMD DM, return true regardless of population statistics
			return true;
		}

		/// If all gahl reference is not minor, filter all the variant's fields
		if (Integer.parseInt(line[minorfield]) != 1) {
			for (int[] a : populations) {
				if ((a.length == 3) && (!frequencyChecker(line[a[1]], line[a[2]], rejectionThreshold))) {
					/// if there is raw data, use the frequencychecker method (above)
					return false; // fail

				} else if ((a.length == 2) && Boolean.parseBoolean(line[a[1]])) {
					/// if the length is one, it is a filter field - just parse the boolean.
					return false; // fail
				}
			}
		}
		// if all gahl reference is minor
		else {
			for (int[] a : populations) {
				if ((a.length == 3) && (!frequencyChecker(line[a[0]], line[a[2]], rejectionThreshold))) {
					/// if there is raw data, use the frequencychecker method (above)
					return false; // fail
				}

				// else if ((a.length == 2) && Boolean.parseBoolean(line[a[0]])) { //if the
				// length is one, it is a filter field - just parse the boolean
				// return false; //fail
				// }
			}
		}
		return true;
	}

	/*
	 * Returns true if mitochondrial variant isn't too common.
	 */

	public boolean Mcheck(String[] line) {
		if (line[hgmdIndex].equals("DM")) { // if it's HGMD DM, return true regardless of population statistics
			return true;
		}
		// if all gahl reference is not minor
		// filter all the variance fields

		if (Integer.parseInt(line[minorfield]) != 1) {
			for (int[] a : populations) {
				if ((a.length == 3) && (!frequencyChecker(line[a[1]], line[a[2]], rejectionThreshold / 4))) {
					/// if there is raw data, use the frequencychecker method (above)
					return false; // fail

				} else if ((a.length == 2) && Boolean.parseBoolean(line[a[1]])) {
					/// if the length is one, it is a filter field - just parse the boolean
					return false; // fail
				}
			}
		}
		/// if all gahl reference is minor
		else {
			for (int[] a : populations) {
				if ((a.length == 3) && (!frequencyChecker(line[a[0]], line[a[2]], rejectionThreshold / 4))) {
					/// If there is raw data, use the frequencychecker method (above)
					return false; // fail

				} else if ((a.length == 2) && Boolean.parseBoolean(line[a[0]])) {
					/// if the length is one, it is a filter field - just parse the boolean
					return false; // fail
				}
			}
		}
		return true;
	}

	/*
	 * Tests for the genotype of an input string - returns 0 for homVar, 1 for
	 * homRef, 3 for heterozygous, -5 for other.
	 */
	public int homHet(String genotype, String ref, String var) {
		String homRef;
		String homVar;
		String het1;
		String het2;

		if ((ref.length() > 1) || (var.length() > 1)) {
			homRef = ref + ":" + ref;
			homVar = var + ":" + var;
			het1 = ref + ":" + var;
			het2 = var + ":" + ref;
			// NA="NA";
		} else {
			homRef = ref + ref;
			homVar = var + var;
			het1 = ref + var;
			het2 = var + ref;
			// NA="NA";
		}

		// homvar or hemvar
		if (genotype.equals(homVar) || genotype.equals(var)) {
			return 0;
			// homref or hemref
		} else if (genotype.equals(homRef) || genotype.equals(ref)) {
			return 1;
			// het
		} else if ((genotype.equals(het1)) || (genotype.equals(het2))) {
			return 3;
		}

		else {
			return -5;
		}
	}

	/*
	 * Takes an ArrayList of halfHets, tests what phases takes the halfHets for a
	 * gene and processes it down to the compound heterozygous recessive matches.
	 */
	public ArrayList<halfHetNode> findCHetPairs(ArrayList<halfHetNode> halfHets) {
		boolean[] matches = new boolean[halfHets.size()]; // this array tracks which halfHets are phased; if not phased,
															// the corresponding position is left as false

		for (int i = 0; i < matches.length; i++) {
			matches[i] = false;
		}

		for (int i = 0; i < halfHets.size(); i++) { // loop through the entire array of halfHets
			if (halfHets.get(i).parent == 0) { // skip paternal half hets, everything that passes is maternal
				continue;
			}

			for (int j = 0; j < halfHets.size(); j++) { // reloop through the entire array of halfHets and attempt to
														// phase
				if (halfHets.get(j).parent == 0) { // if the halfHet we're attempting to phase with is paternal
					boolean unaffs = true; // boolean to track for this pair if unaffected siblings phase, which they
											// should not

					for (int k = 0; k < halfHets.get(i).unAffSibs.size(); k++) {
						int sibsVal = halfHets.get(i).unAffSibs.get(k) + halfHets.get(j).unAffSibs.get(k);

						if (sibsVal == 6) { // if the same unaffected sib is ever het in both lines (which would be
											// 3+3)...
							unaffs = false; // make the temp variable false
						}
					}

					if (unaffs) { // if one of the sibs is false, this will not pass. Otherwise
						halfHets.get(i).score(halfHets.get(j)); // score the halfHets (see the halfHetNode object at
																// the top)
						matches[i] = true; // both positions in the matches array are now true, since these two both
											// phase
						matches[j] = true;
					}
				}
			}
		}

		for (int i = halfHets.size() - 1; i >= 0; i--) { // loop through the halfHet array backwards
			if (!matches[i]) { // if any of the positions in matches is false, that means there was no phasing
				halfHets.remove(i); // since no phasing, remove
			}
		}

		return halfHets; // Return the remaining halfHets, all of which phase
	}

	/*
	 * Returns a string to be printed in the output for this variant if homozygous
	 * recessive.
	 */
	public String PrintHomo(String[] splitline) {
		String output = "";
		for (int i = 0; i < naIndex; i++) { // all fields up the family data
			output += splitline[i] + "\t";
		}
		output += "0" + "\t"; // add the beststrength field
		output += "0" + "\t"; // add the bestvmm CADD field
		output += "0" + "\t"; // add the bestvmm EIGEN field
		output += "0" + "\t"; // add the bestvmm EIGENPC field
		output += Integer.toString(0) + "\t"; // MendInconsist
		output += Integer.toString(1) + "\t"; // MendHomRec
		output += Integer.toString(0) + "\t"; // Xlink
		output += Integer.toString(0) + ","; // MendHetRec
		// add the bestvmm EIGENPC field

		for (int i = naIndex; i < splitline.length; i++) { // add family data
			output += "\t" + splitline[i];
		}
		return output; // return the string
	}

	/*
	 * Returns a string to be printed in the output for this variant if Call No Call
	 * (CNC).
	 */
	public String PrintCnC(String[] splitline) {
		String output = "";
		for (int i = 0; i < naIndex; i++) { // all fields up the family data
			output += splitline[i] + "\t";
		}
		output += "0" + "\t"; // add the beststrength field
		output += "0" + "\t"; // add the bestvmm CADD field
		output += "0" + "\t"; // add the bestvmm EIGEN field
		output += "0" + "\t"; // add the bestvmm EIGENPC field
		output += Integer.toString(4) + "\t"; // MendInconsist
		output += Integer.toString(0) + "\t"; // MendHomRec
		output += Integer.toString(0) + "\t"; // Xlink
		output += Integer.toString(0) + ","; // MendHetRec
		// add the bestvmm EIGENPC field

		for (int i = naIndex; i < splitline.length; i++) { // add family data
			output += "\t" + splitline[i];
		}
		return output; // return the string
	}

	/*
	 * Returns a string to be printed in the output for this variant if
	 * mitochondrial.
	 */
	public String PrintchrM(String[] splitline) {
		String output = "";
		for (int i = 0; i < naIndex; i++) { // all fields up the family data
			output += splitline[i] + "\t";
		}
		output += "0" + "\t"; // add the beststrength field
		output += "0" + "\t"; // add the bestvmm CADD field
		output += "0" + "\t"; // add the bestvmm EIGEN field
		output += "0" + "\t"; // add the bestvmm EIGENPC field
		output += Integer.toString(0) + "\t"; // MendInconsist
		output += Integer.toString(0) + "\t"; // MendHomRec
		output += Integer.toString(0) + "\t"; // Xlink
		output += Integer.toString(0) + ","; // MendHetRec
		// add the bestvmm EIGENPC field

		for (int i = naIndex; i < splitline.length; i++) { // add family data
			output += "\t" + splitline[i];
		}
		return output; // return the string
	}

	/*
	 * Returns a string to be printed in the output for this variant if X-linked.
	 */
	public String Xlinkprint(String[] splitline) {
		String output = "";
		for (int i = 0; i < naIndex; i++) { // all fields up the family data
			output += splitline[i] + "\t";
		}
		output += "0" + "\t"; // add the beststrength field
		output += "0" + "\t"; // add the bestvmm CADD field
		output += "0" + "\t"; // add the bestvmm EIGEN field
		output += "0" + "\t"; // add the bestvmm EIGENPC field
		output += Integer.toString(0) + "\t"; // MendInconsist
		output += Integer.toString(0) + "\t"; // MendHomRec
		output += Integer.toString(1) + "\t"; // Xlink
		output += Integer.toString(0) + ","; // MendHetRec
		// add the bestvmm EIGENPC field

		for (int i = naIndex; i < splitline.length; i++) { // add family data
			output += "\t" + splitline[i];
		}
		return output; // return the string
	}

	/*
	 * Returns a string to be printed in the output for this variant if Mendelian
	 * inconsistent.
	 */
	public String PrintInconsist(String[] splitline, int inconsistent) {
		// 1 inconsistent, 2 de novo, 3, hemi
		String output = "";
		for (int i = 0; i < naIndex; i++) { // all fields up the family data
			output += splitline[i] + "\t";
		}
		output += "0" + "\t"; // Add the beststrength field
		output += "0" + "\t"; // Add the bestvmm CADD field
		output += "0" + "\t"; // Add the bestvmm EIGEN field
		output += "0" + "\t"; // Add the bestvmm EIGENPC field
		output += Integer.toString(inconsistent) + "\t"; // MendInconsis
		output += Integer.toString(0) + "\t"; // MendHomRec
		output += Integer.toString(0) + "\t"; // Xlink
		output += Integer.toString(0) + ","; // MendHetRec
		// Add the bestVMM EIGENPC field
		for (int i = naIndex; i < splitline.length; i++) { // Add family data
			output += "\t" + splitline[i];
		}
		return output; // Return the string
	}

	// public void descriptions() {
	// String intro = "Hello, thank you for using the KaylaKode.\n" + "This program
	// includes filters for: \n"
	// + "Xlink \n" + "Mendelian Inconsistent \n" + "Homozygous Recessive \n" +
	// "Compound Heterozygous \n"
	// + "DeNovo Dominant \n"
	// + "Please make sure Strength_Config.txt and Gene_Config_Splits_USE.txt are
	// saved in the same path as the .jar file.";
	// JOptionPane.showMessageDialog(new Frame("Introduction"), intro);
	// }
	
	
}
