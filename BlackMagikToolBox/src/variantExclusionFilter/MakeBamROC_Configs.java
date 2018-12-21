package variantExclusionFilter;


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
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

import general.Pedigree;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class MakeBamROC_Configs {
	/// Initialize global variables
	private BufferedReader vsData;
	private File vsFile;

	private String Line;
	private String[] curLine;
	private static ArrayList<String> headers;
	
	private int chrIndex;
	private int posIndex;
	private int posendIndex;
	private int refIndex;
	private int varIndex;

	//private ArrayList<Integer> family;

	private int minorField;
	private ArrayList<Integer> totACInds;
	private ArrayList<Integer> totANInds;
	private ArrayList<Integer> totAFInds;
	private ArrayList<Integer> homRefInds;
	private ArrayList<Integer> homVarInds;
	private ArrayList<Integer> hetInds;
	private ArrayList<Integer> othInds;
	private ArrayList<Integer> genotypeInds;
	private ArrayList<Integer> hemRefInds;
	private ArrayList<Integer> hemVarInds;
	private ArrayList<Integer> maxMAF1Inds;
	private ArrayList<Integer> maxMAF2Inds;
	private ArrayList<Integer> afrAFInds;
	private ArrayList<Integer> easAFInds;
	private ArrayList<Integer> amrAFInds;
	private ArrayList<Integer> eurAFInds;
	private ArrayList<Integer> sasAFInds;
	private ArrayList<Integer> othAFInds;
	private ArrayList<Integer> afrACInds;
	private ArrayList<Integer> easACInds;
	private ArrayList<Integer> amrACInds;
	private ArrayList<Integer> eurACInds;
	private ArrayList<Integer> sasACInds;
	private ArrayList<Integer> othACInds;
	private ArrayList<Integer> afrANInds;
	private ArrayList<Integer> easANInds;
	private ArrayList<Integer> amrANInds;
	private ArrayList<Integer> eurANInds;
	private ArrayList<Integer> sasANInds;
	private ArrayList<Integer> othANInds;

	private int homVarThreshold;
	private int hemVarThreshold;
	private int hetThreshold;
	
	private String destination;
	private File pedFile;
	private ArrayList<Boolean> sibsAffected; // boolean values of whether siblings are affected
	private ArrayList<Boolean> sibsGender;// boolean values of the siblings gender, 1 is male
	
	private Double[] cutArr;
	private File headerFile;
	private ArrayList<Integer> genotypeThreshList;
	
	private File ROCOutput;
	private File BAMOutput;
	
	
	/*
	 * The purpose of this class is to control the inputs and outputs of the
	 * MakeBamROC.java program. The inputs should include (1) a bam curated VS file
	 * (after ideally having been run through the FBKaylaScriptBam/BamCuration.java
	 * program) and (2) a pedigree file with all the families of the cohort of
	 * analysis (in this case, the UDP cohort). There will be 2 tab-delimited text
	 * outputs: (1) one output with the same variants as the VS file but annotated
	 * with several flags (the main ones being pop, bam curation quality, Mendelian
	 * model, and misc.) and the salvage status of each pedigree member's allele for
	 * a given variant, and (2) an output with the same information as (1) BUT with
	 * variants filtered out based on the flags annotated in (1) and a given
	 * CADD/VMM cutoff score, as determined by an array of CADD cutoff values.
	 * -Anchi Wu, 9.29.17
	 * 
	 * This program is similar to MakeBamROC.java but is different in that it requests
	 * more user input (e.g., the header config files) to allow more flexibility for
	 * the user. 
	 * -Anchi Wu, 5.16.2018
	 */
	
	public MakeBamROC_Configs(File vsFile, String dest, File ped, Double [] caddThresh, int homVarThresh, int hetThresh, int hemVarThresh, File headerConfig, ArrayList<Integer> genoThreshList,
			ArrayList<Boolean> sibsAffect, ArrayList<Boolean> sibGend) throws IOException {
		
		this.vsFile = vsFile;
		this.destination = dest;
		this.pedFile = ped;
		this.cutArr = caddThresh;
		this.homVarThreshold = homVarThresh;
		this.hetThreshold = hetThresh;
		this.hemVarThreshold = hemVarThresh;
		this.headerFile = headerConfig;
		this.genotypeThreshList = genoThreshList;
		this.sibsAffected = sibsAffect;
		this.sibsGender = sibGend;
		
		variantExclusionFilter();
	}
	
	
	public void variantExclusionFilter() throws IOException {
		
		
		/// Input VarSifter file
		//File vsFile = TXTFile.getVSFile();
		// String vsFileString = vsFile.getPath();

		/// Input pedigree file
		//File pedFile = TXTFile.getPedigreeFile();
		BufferedReader pedData = new BufferedReader(new FileReader(pedFile));

		/// Input BAM file directory
		//File bamlo = TXTFile.getBAMDirectoryFile();
		
		/*
		String bamloPath = bamlo.getPath();
		HashMap<String, String> BamMap = new HashMap<String, String>();
		FindBam.InitializeBam(bamloPath, BamMap);
		*/
		
		/// Input CADD threshold
		//Double[] cutArr = delThresholdCalc();
		
		/// Ask user to input population varaiant genotype (homvar, het, hemvar)
		/// maximums. (Default is 3)
		//getPopulationGenotypes();
		

		String outputLocROC = destination + "_ROCFiltered.vs";
		String outputLocBam = destination + "_BamCurated.vs";
		
		//MakeBamROC_Configs curateFilter = new MakeBamROC_Configs();
		
		JFrame jframe = canceler();
		jframe.setLocationRelativeTo(null);
		jframe.setVisible(true);

		/// Iterate through each family being analyzed
		String Line = pedData.readLine();
		while (Line != null) {
			String[] pedLine = Line.split("\t");
			initializer(vsFile, outputLocROC, outputLocBam, pedLine, cutArr);
			Line = pedData.readLine();
		}
		
		pedData.close();
		jframe.dispose();
		
	}
	
	public File getOutput() {
		
		return ROCOutput;
		
	}
	
	/*
	 * This is the main driver of the combined bam curation and ROC filtration steps
	 * of the Forwwards Backwards Analysis. This takes inputs from the given bam
	 * curated VS file and outputs two files: (1) one output with the same variants
	 * as the VS file but annotated with several flags (the main ones being pop, bam
	 * curation quality, Mendelian model, and misc.) and the salvage status of each
	 * pedigree member's allele for a given variant, and (2) an output with the same
	 * information as (1) BUT with variants filtered out based on the flags
	 * annotated in (1) and a given CADD/VMM cutoff score, as determined by
	 * caddcutofffeed.
	 */
	public void initializer(File vsFile, String outputROCFile, String outputBamFile, String[] trio, Double[] caddcutofffeed) throws IOException {
		/// Initialize variant info variables
		int Igenename;
		int Iind;
		int Igene;
		int Itype;
		int ICADD;
		int naIndex;

		int IMendInconsis;
		int IMendHomRec;
		int IXlinkd;
		int IMendHetRec; // Index of the other part of a potential cmpd het pairing
		int Ivmmc;

		// Misc
		int IDMtag;
		int IMCAP;

		// Bam Quality
		int IerrorVal;
		int Isnr;
		
		/// Open VS file and initialize output files
		vsData = new BufferedReader(new FileReader(vsFile));
		
		ROCOutput = new File(outputROCFile);
		BAMOutput = new File(outputBamFile);
		
		PrintWriter bamFilterWriter = new PrintWriter(BAMOutput);
		PrintWriter ROCFilterWriter = new PrintWriter(ROCOutput);

		/// Initialize and index variant headers
		Line = vsData.readLine();
		curLine = Line.split("\t");

		headers = new ArrayList<String>();
		for (int i = 0; i < curLine.length; i++) {
			headers.add(curLine[i]);
		}

		ArrayList<String> mygene = new ArrayList<String>();

		/// Initialize and index variant information variables
		Iind = headers.indexOf("Index");
		Igene = headers.indexOf("Gene_name");
		chrIndex = headers.indexOf("Chr");
		posIndex = headers.indexOf("LeftFlank");
		posendIndex = headers.indexOf("RightFlank");
		refIndex = headers.indexOf("ref_allele");
		varIndex = headers.indexOf("var_allele");
		
		int typeIndex = headers.indexOf("type");
		Itype = headers.indexOf("loc_type");
		Igenename = headers.indexOf("Gene_name");
		ICADD = headers.indexOf("PHRED");

		IMendInconsis = headers.indexOf("MendInconsis");
		IMendHomRec = headers.indexOf("MendHomRec");
		IXlinkd = headers.indexOf("Xlink");
		IMendHetRec = headers.indexOf("MendHetRec");
		Ivmmc = headers.indexOf("Best_VMMC");

		// Misc
		IDMtag = headers.indexOf("HGMDtags");
		IMCAP = headers.indexOf("M_CAPv1.0");

		// Bam Quality
		IerrorVal = headers.indexOf("Error");
		Isnr = headers.indexOf("SNR");

		/// Initialize user-input fields
		minorField = 0;
		totACInds = new ArrayList<Integer>();
		totANInds = new ArrayList<Integer>();
		totAFInds = new ArrayList<Integer>();
		homRefInds = new ArrayList<Integer>();
		homVarInds = new ArrayList<Integer>();
		hetInds = new ArrayList<Integer>();
		othInds = new ArrayList<Integer>();
		genotypeInds = new ArrayList<Integer>();
		hemRefInds = new ArrayList<Integer>();
		hemVarInds = new ArrayList<Integer>();
		maxMAF1Inds = new ArrayList<Integer>();
		maxMAF2Inds = new ArrayList<Integer>();
		afrAFInds = new ArrayList<Integer>();
		easAFInds = new ArrayList<Integer>();
		amrAFInds = new ArrayList<Integer>();
		eurAFInds = new ArrayList<Integer>();
		sasAFInds = new ArrayList<Integer>();
		othAFInds = new ArrayList<Integer>();
		afrACInds = new ArrayList<Integer>();
		easACInds = new ArrayList<Integer>();
		amrACInds = new ArrayList<Integer>();
		eurACInds = new ArrayList<Integer>();
		sasACInds = new ArrayList<Integer>();
		othACInds = new ArrayList<Integer>();
		afrANInds = new ArrayList<Integer>();
		easANInds = new ArrayList<Integer>();
		amrANInds = new ArrayList<Integer>();
		eurANInds = new ArrayList<Integer>();
		sasANInds = new ArrayList<Integer>();
		othANInds = new ArrayList<Integer>();

		/// Read config file for remaining headers
		//File headerFile = TXTFile.getHeaderConfigFile();
		readHeaders(headerFile);
		
		/// Define family relationship based on pedigree
		ArrayList<Integer> family = Pedigree.getPedigreefromString(trio, headers);

		/// Establish where patient variant information begins in the columns
		naIndex = genotypeIndices(curLine, curLine.length);

		/// Establish AllCount and GoodCount indices for each family member
		HashMap<Integer, Integer[]> familyCountIndMap = new HashMap<Integer, Integer[]>();
		for (int i = naIndex; i < headers.size(); i++) {
			if (headers.get(i).contains(".NA")) {
				Integer[] countArr = { new Integer(i + 3), new Integer(i + 4) }; // [Index of AllCount,Index of
																					// GoodCount]
				familyCountIndMap.put(i, countArr);
			}
		}

		/// Write headers to the output files
		// Add flag headers to the beginning of the output files
		String[] flagHeaders = { "PopFlags", "BamFlags", "MendFlags", "MiscFlags" };

		String headerAndFlags = curateEnd(curLine, naIndex, family, flagHeaders); // Concatenate flag headers and other
																					// headers
		String fullHeader = headerAndFlags;

		// Print headers for the output files
		bamFilterWriter.println(fullHeader);
		ROCFilterWriter.println(fullHeader);

		/// Prepare for the compound het
		ArrayList<CmphPair_Configs> cmpdhetpair = new ArrayList<CmphPair_Configs>();

		/// List of bad flags for later ROC filtering
		ArrayList<String> badFlags = new ArrayList<String>();
		badFlags.add("HP"); // high population stats of the variant
		badFlags.add("LC"); // low sampling of the position/variant in the population
		badFlags.add("HF"); // too high pop freq
		badFlags.add("HC"); // too high absolute pop count
		// badFlags.add("GP"); // too high in the all gahl cohort
		badFlags.add("BG"); // in the big gene list
		badFlags.add("BR"); // bad read
		badFlags.add("OU"); // bad bam curation quality

		/// Big gene list (for BG flag)
		ArrayList<String> bigGene = new ArrayList<String>();
		bigGene.add("NEB");
		bigGene.add("TTN");
		bigGene.add("OBSCN");
		bigGene.add("KRT");
		bigGene.add("OR");
		bigGene.add("TAS");
		
		/*
		/// Input minimum genotype count threshold
		ArrayList<Integer> genotypeThreshList = new ArrayList<Integer>();
		for (int genoHeaderInd : genotypeInds) {
			genotypeThreshList.add(genoCountThreshold(genoHeaderInd));
		}
		*/
		
		/// Read through each variant of the VS file
		while ((Line = vsData.readLine()) != null) {
			curLine = Line.split("\t");

			/// Define variant information variables for current variant
			/// Required columns for all input VS files
			String chrom = curLine[chrIndex];
			String leftflank = curLine[posIndex];
			String rightflank = curLine[posendIndex];
			String refprint = curLine[refIndex];
			String altprint = curLine[varIndex];
			int place = Integer.parseInt(leftflank) + 1;
			int endplace = Integer.parseInt(rightflank) - 1;

			int index = Integer.parseInt(curLine[Iind]);
			String gene = curLine[Igenename];
			int inconsist = Integer.parseInt(curLine[IMendInconsis]);
			int homrec = Integer.parseInt(curLine[IMendHomRec]);
			int Xlink = Integer.parseInt(curLine[IXlinkd]);
			String cmpd = curLine[IMendHetRec];
			
			//If loctype column does not exist
			String loctype = "";
			
			if(Itype != -1) {
				loctype = curLine[Itype];
			} else if(typeIndex != -1) {
				
				if((curLine[typeIndex].contains("synonymous") || curLine[typeIndex].contains("SYNONYMOUS")) || (curLine[typeIndex].contains("stopgain") || curLine[typeIndex].contains("STOPGAIN"))
						|| (curLine[typeIndex].contains("frameshift") || curLine[typeIndex].contains("FRAMESHIFT"))) {
					
					loctype = "exonic";
					
				}
				
			} else {
				
				JOptionPane.showMessageDialog(new Frame("Error"),
						"There is no 'loc_type' or 'type' header in the input file. Please close and try again.");
				vsData.close();
				bamFilterWriter.close();
				ROCFilterWriter.close();
				System.exit(0);
				
			}
			
			
			double caddscore = Double.parseDouble(curLine[ICADD]);
			double vmmcscore = Double.parseDouble(curLine[Ivmmc]);

			/// Retrieve values for each of the user-specified headers
			int refMinor = Integer.parseInt(curLine[minorField]); /// Denotes if reference is minor
			
			/// Extract values for each parameter from variant's line
			ArrayList<Integer> totACList = retrieveValInt(totACInds); // Population variant allele count (AC)
			ArrayList<Integer> totANList = retrieveValInt(totANInds); // Pop total allele count (AN)
			ArrayList<Double> totAFList = retrieveValDouble(totAFInds); // Pop variant allele frequency (AF)
			ArrayList<Integer> homRefList = retrieveValInt(homRefInds); // Pop homref count
			ArrayList<Integer> homVarList = retrieveValInt(homVarInds); // Pop homvar count
			ArrayList<Integer> hetList = retrieveValInt(hetInds); // Pop het count
			ArrayList<Integer> othList = retrieveValInt(othInds); // Pop other genotype count
			ArrayList<Integer> genotypeList = retrieveValInt(genotypeInds); // Pop total genotype count
			ArrayList<Integer> hemRefList = retrieveValInt(hemRefInds); // Pop hemref count
			ArrayList<Integer> hemVarList = retrieveValInt(hemVarInds); // Pop hemvar count
			ArrayList<Boolean> maxMAF1List = retrieveValBoolean(maxMAF1Inds); // MaxMAF>1%
			ArrayList<Boolean> maxMAF2List = retrieveValBoolean(maxMAF2Inds); // MaxMAF>2%
			ArrayList<Double> afrAFList = retrieveValDouble(afrAFInds); // AFR AF
			ArrayList<Double> easAFList = retrieveValDouble(easAFInds); // EAS AF
			ArrayList<Double> amrAFList = retrieveValDouble(amrAFInds); // AMR AF
			ArrayList<Double> eurAFList = retrieveValDouble(eurAFInds); // EUR AF
			ArrayList<Double> sasAFList = retrieveValDouble(sasAFInds); // SAS AF
			ArrayList<Double> othAFList = retrieveValDouble(othAFInds); // other ethnicity (OTH) AF
			ArrayList<Integer> afrACList = retrieveValInt(afrACInds); // AFR AC
			ArrayList<Integer> easACList = retrieveValInt(easACInds); // EAS AC
			ArrayList<Integer> amrACList = retrieveValInt(amrACInds); // AMR AC
			ArrayList<Integer> eurACList = retrieveValInt(eurACInds); // EUR AC
			ArrayList<Integer> sasACList = retrieveValInt(sasACInds); // SAS AC
			ArrayList<Integer> othACList = retrieveValInt(othACInds); // OTH AC
			ArrayList<Integer> afrANList = retrieveValInt(afrANInds); // AFR AN
			ArrayList<Integer> easANList = retrieveValInt(easANInds); // EAS AN
			ArrayList<Integer> amrANList = retrieveValInt(amrANInds); // AMR AN
			ArrayList<Integer> eurANList = retrieveValInt(eurANInds); // EUR AN
			ArrayList<Integer> sasANList = retrieveValInt(sasANInds); // SAS AN
			ArrayList<Integer> othANList = retrieveValInt(othANInds); // OTH AN

			// Future TODO: Allow user to choose which ethnicities they'd like to use,
			// between afr, eas, amr, eur, sas, and oth

			// Bam Quality
			double errorVal = Double.parseDouble(curLine[IerrorVal]);

			String DM = curLine[IDMtag];
			Double mcapscore = Double.parseDouble(curLine[IMCAP]);

			/// Initialize flag columns
			HashSet<String> popFlag = new HashSet<String>();
			HashSet<String> bamFlag = new HashSet<String>();
			HashSet<String> mendFlag = new HashSet<String>();
			HashSet<String> miscFlag = new HashSet<String>();

			// Flag based on the MCAP score or HGMD annotation (supersedes other flags)
			boolean special = false;
			if (DM.equals("DM") || mcapscore > 0.5) {
				special = true;
			}
			/// To be added in later once Tom's finished annotating in ClinVar
			// if (clinVarVal.equals("pathogenic")) {
			// special = true;
			// }

			/// Flag based on if gene is a big gene
			Boolean count = true;
			String[] geneList = gene.split(",");
			for (int i = 0; i < geneList.length; i++) {
				String tmpGene = geneList[i].replaceAll("[^A-Za-z0-9]+", "");
				if (bigGene.indexOf(tmpGene) != -1) { // Works for all but the OR genes
					count = false;
					miscFlag.add("BG");
				}
				for (String bigGeneOne : bigGene) {
					if (tmpGene.matches(bigGeneOne + "\\d+(.*)")) {
						count = false;
						miscFlag.add("BG");
					}
				}

			}

			/// Flag based on genotype population counts and freqs (# of ppl who've been
			/// sequenced at this variant position)
			boolean MAF = true;	// True means maxMAF < 1%
			boolean MAF2 = true;	// True means maxMAF < 2%
			boolean looskmaf = true; // True means ethnic subpop MAF < threshold

			if (refMinor == 0) { /// Reference is major allele
				/// General pop maxMAF filter
				for (boolean maxMAF1 : maxMAF1List) {
					if (maxMAF1) {
						MAF = false;
					}
				}

				for (boolean maxMAF2 : maxMAF2List) {
					if (maxMAF2) {
						MAF2 = false;
					}
				}

				/// Ethnic-specific maxMAF filter
				/// AFR
				for (double afrAF : afrAFList) {
					if (afrAF > 0.02) {
						looskmaf = false;
					}
				}

				/// EAS
				for (double easAF : easAFList) {
					if (easAF > 0.02) {
						looskmaf = false;
					}
				}

				/// AMR
				for (double amrAF : amrAFList) {
					if (amrAF > 0.02) {
						looskmaf = false;
					}
				}

				/// EUR
				for (double eurAF : eurAFList) {
					if (eurAF > 0.02) {
						looskmaf = false;
					}
				}

				/// SAS
				for (double sasAF : sasAFList) {
					if (sasAF > 0.02) {
						looskmaf = false;
					}
				}

				/// Other
				for (double othAF : othAFList) {
					if (othAF > 0.02) {
						looskmaf = false;
					}
				}

			} else { /// Reference is the minor allele
				boolean isOth = false;
				for (int othCount : othList) {
					if (othCount != 0) {
						miscFlag.add("Oth");
						isOth = true;
						break;
					}
				}
				/// Assuming locus is bi-allelic (only the variant and the ref), do inverse
				/// Poisson calculation
				if (!isOth) {
					/// Note: This is the step that requires each AC entry to have an accompanying
					/// AN entry.
					MAF = invPoissonCalcExec(totACList, totANList, 0.01); // MaxMAF < 1% Filter
					MAF2 = invPoissonCalcExec(totACList, totANList, 0.02); // MaxMAF < 2% Filter
					
					boolean afrMAF2 = invPoissonCalcExec(afrACList, afrANList, 0.02);
					boolean easMAF2 = invPoissonCalcExec(easACList, easANList, 0.02);
					boolean amrMAF2 = invPoissonCalcExec(amrACList, amrANList, 0.02);
					boolean eurMAF2 = invPoissonCalcExec(eurACList, eurANList, 0.02);
					boolean sasMAF2 = invPoissonCalcExec(sasACList, sasANList, 0.02);
					boolean othMAF2 = invPoissonCalcExec(othACList, othANList, 0.02);

					if (afrMAF2 || easMAF2 || amrMAF2 || eurMAF2 || sasMAF2 || othMAF2) {
						looskmaf = false; /// Ethnic subpops > 1% and 2%
					}
				}
			}

			/// Low genotype count filter
			for (int i = 0; i < genotypeThreshList.size(); i++) {
				int genoMin = genotypeThreshList.get(i);
				int genoCount = genotypeList.get(i);
				if (genoCount < genoMin) {
					miscFlag.add("LC");
				}
			}

			// Confirmed disease marker
			if (special) {
				miscFlag.add("DM");
				count = true;
			}

			// maxMAF filter
			if (!MAF || !MAF2 || !looskmaf) {
				popFlag.add("HP");
				popFlag.add("HF");
			}

			// Write string with variant-specific information (before family)
			// String outMeta = curLine[0];
			//
			// for (int i = 1; i < naIndex; i++) {
			// outMeta += "\t" + curLine[i];
			// }

			/// Determine total reads in the family (for later avgreads calculation)
			ArrayList<String[]> Arrayreadcount = new ArrayList<String[]>();
			int totalReads = 0;
			for (int familyInd : family) {
				String allCountStr = curLine[familyInd + 3];
				String goodCountStr = curLine[familyInd + 4];

				String[] allCount = allCountStr.split(";");
				String[] goodCount = goodCountStr.split(";");

				// Extract and flag salvage status
				bamFlag.add(allCount[0]);
				if (allCount[0].contains("NN") || allCount[0].contains("NP")) {
					bamFlag.add("SA"); // Mark as salvaged
				}

				String[] countArr = { allCountStr, goodCountStr };
				Arrayreadcount.add(countArr);

				// Calculate total reads of the variant across the family
				totalReads += Integer.parseInt(goodCount[0]) + Integer.parseInt(goodCount[1]);
				// outMeta += "\t" + curLine[familyInd] + "\t" + allCountStr + "\t" +
				// goodCountStr;
			}

			/// Calculate VMM cutoff
			Double caddcutoff = caddcutofffeed[0];
			if (loctype.contains("exonic") && !loctype.contains("ncRNA")) {
				caddcutoff = caddcutofffeed[1];
			}
			Double vmmcutoff = 30 + 2 * (caddcutofffeed[1] - 20);

			/// Filter based on VMM and CADD score cutoffs
			// Compound het filter
			if (vmmcscore > vmmcutoff && MAF2 && count) {
				cmpdhetpair.add(new CmphPair_Configs(Line, IMendHetRec));
				boolean findamatch = false;
				for (int i = 0; i < cmpdhetpair.size(); i++) {
					String indextm = curLine[Iind].replaceAll("\\D+", "");
					if (cmpdhetpair.get(i).pair(indextm, Line)) {
						findamatch = true;
					}
				}
			} else if (inconsist == 4) { // CnC
				// Label as CnC and add to the CnC map
				mendFlag.add("CNC");

				/// Bam curation!!
				// Calculate average number of reads per person in a family
				double avgread = (double) totalReads / (double) family.size();
				if (errorVal > 0.01 && avgread < 4) {
					bamFlag.add("BR"); // Bad read
				} else if (errorVal > 0.02) {
					bamFlag.add("BR"); // Bad read
				}

				// Set the flags
				String popStr = set2Str(popFlag);
				String bamStr = set2Str(bamFlag);
				String mendStr = set2Str(mendFlag);
				String miscStr = set2Str(miscFlag);

				String[] fullFlagStr = { popStr, bamStr, mendStr, miscStr };

				/// Write to the ROC output file
				// Filter based on flags
				Boolean filterPass = true;
				if (mendFlag.contains("CNC")) {
					if (miscFlag.contains("LC") || bamFlag.contains("BR")) {
						filterPass = false;
					}
				}

				if (filterPass) {
					ROCFilterWriter.println(curateEnd(curLine, naIndex, family, fullFlagStr));
				} else if (special && !bamFlag.contains("OU") && !bamFlag.contains("BR")) {
					// If was marked as DM by HGMD and decent quality bam region
					ROCFilterWriter.println(curateEnd(curLine, naIndex, family, fullFlagStr));
				}

			} else if (vmmcscore == 0 && caddscore >= caddcutoff) { // All other models
				// De novo
				if (inconsist == 2) {
					mendFlag.add("DN");
					// if (Ghet > 1) {
					// miscFlag.add("GP");
					// }
					if (refMinor == 1) {
						for (int hetCount : hetList) {
							if (hetCount > hetThreshold) {
								popFlag.add("HP");
								popFlag.add("HC");
							}
						}

						/// TODO: ADD ref allele exac threshold for de novo
						// for (int refCount : )
						// if (hghetc > 2 || refalleleexac > 3) {
						// popFlag.add("HP");
						// popFlag.add("HC");
						// }
						// } else if (hghetc > 2 || alleleexac > 3 || !MAF) {
						// popFlag.add("HP");
						// popFlag.add("HC");
					} else {
						for (int hetCount : hetList) {
							if (hetCount > hetThreshold) {
								popFlag.add("HP");
								popFlag.add("HC");
							}
						}
						/// TODO: Add in the total allele filter for de novo
						// for (int varCount : totACList) {
						// if (varCount > homVarThreshold)
						// }
					}
					if (!bamFlag.contains("SA")) {
						miscFlag.add("UD"); // unsalvaged de novo
					}
				} else if (inconsist == 3) { // Hemizygous
					mendFlag.add("HM");
					if (refMinor == 1) {
						for (int homRefCount : homRefList) {
							if (homRefCount > homVarThreshold) {
								popFlag.add("HP");
								popFlag.add("HC");
							}
						}
					} else {
						for (int homVarCount : homVarList) {
							if (homVarCount > homVarThreshold) {
								popFlag.add("HP");
								popFlag.add("HC");
							}
						}
					}
				} else if (Xlink == 1) {
					// Xlink
					mendFlag.add("XL");
					if (refMinor == 1) {
						for (int hemRefCount : hemRefList) {
							if (hemRefCount > hemVarThreshold) {
								popFlag.add("HP");
								popFlag.add("HC");
							}
						}
					} else {
						for (int hemVarCount : hemVarList) {
							if (hemVarCount > hemVarThreshold) {
								popFlag.add("HP");
								popFlag.add("HC");
							}
						}
					}
				} else if (homrec == 1) { // Homrec
					mendFlag.add("HR");
					if (refMinor == 1) {
						boolean salvageStatus = true;
						for (int homRefCount : homRefList) {
							if (homRefCount > homVarThreshold) {
								popFlag.add("HP");
								popFlag.add("HC");
								salvageStatus = false;
							}
						}
						if (salvageStatus && bamFlag.contains("SA")) {
							bamFlag.add("SR");
						}
					} else {
						boolean salvageStatus = true;
						for (int homVarCount : homVarList) {
							if (homVarCount > homVarThreshold) {
								popFlag.add("HP");
								popFlag.add("HC");
								salvageStatus = false;
							}
						}
						if (salvageStatus && bamFlag.contains("SA")) {
							bamFlag.add("SR");
						}
					}
				}

				/// Bam curation!!
				// Calculate average number of reads per person in a family
				double avgread = (double) totalReads / (double) family.size();
				if (errorVal > 0.01 && avgread < 4) {
					bamFlag.add("BR"); // Bad read
				} else if (errorVal > 0.02) {
					bamFlag.add("BR"); // Bad read
				}

				// Begin annotation based on bam quality
				if (avgread < 20) {
					if (mendFlag.contains("DN")) { // If variant is flagged as de novo, varcount (+/-)=
						int varindex = 1 - refMinor;
						int varcount = 0;

						for (int i = 0; i < family.size(); i++) {
							if(i > 2 && sibsAffected.get(i - 3)) { //Checks if siblings are affected
								varcount += Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[varindex]);
							} else if (i != 2) { // If have at least one parent and proband?? (AW)
								varcount -= Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[varindex]);
							} else {
								varcount += Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[varindex]);
							}
						}

						if (varcount < 3) {
							bamFlag.add("OU"); // Poor bam quality
						}
						
					} else if (mendFlag.contains("CM")) { // If flag is cmpd het
						int hetcount = 0;
						for (int i = 0; i < 2; i++) {
							if (Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[0]) > 2
									&& Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[1]) > 2) {
								hetcount++;
							}
						}

						if (hetcount == 2 || !(Integer.parseInt(Arrayreadcount.get(2)[1].split(";")[0]) > 0
								&& Integer.parseInt(Arrayreadcount.get(2)[1].split(";")[1]) > 0)) {
							bamFlag.add("OU");
						}

					} else if (mendFlag.contains("HR")) { // If flag is homozygous recessive
						int varindex = 1 - refMinor; // 1 - varindex = refMinor
						boolean homrecB = false;
						if (Integer.parseInt(Arrayreadcount.get(0)[1].split(";")[1 - varindex]) < 2) { //father
							homrecB = true;
						}
						if (Integer.parseInt(Arrayreadcount.get(1)[1].split(";")[1 - varindex]) < 2) { //mother
							homrecB = true;
						}
						if (Integer.parseInt(Arrayreadcount.get(2)[1].split(";")[1 - varindex]) > 2) { //proband
							homrecB = true;
						}
						for (int i = 3; i < family.size(); i++) {
							if(!sibsAffected.get(i - 3)) { //Checks if siblings are affected
								if (Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[1 - varindex]) < 2) {
									homrecB = true;
								}
							} else {
								
								if (Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[1 - varindex]) > 2) {
									homrecB = true;
								}
								
							}
						}
						if (homrecB) {
							bamFlag.add("OU");
						}
					} else if (mendFlag.contains("XL")) { // If flag is x-linked
						int varindex = 1 - refMinor;
						boolean xl = false;

						if (Integer.parseInt(Arrayreadcount.get(2)[1].split(";")[1 - varindex]) > 2) {
							xl = true;
						}
						for (int i = 0; i < family.size(); i++) {

							if(sibsGender.get(i - 3)) { //Check if siblings' genders are male
							
								if (!sibsAffected.get(i - 3)) { //Checks if siblings are affected
									if (Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[1 - varindex]) < 2) {
										xl = true;
									}
								} else if (sibsAffected.get(i - 3)) {
									if (Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[1 - varindex]) > 2) {
										xl = true;
									}
								}
							}
						}
						if (xl) {
							bamFlag.add("OU");
						}
					}

				} else {
					if (mendFlag.contains("DN")) { // If flag is de novo
						int varindex = 1 - refMinor;
						int varcount = 0;

						for (int i = 0; i < family.size(); i++) {
							if(i > 2 && sibsAffected.get(i - 3)) { //Checks if siblings are affected
								varcount += Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[varindex]);
							} else if (i != 2) {
								varcount -= Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[varindex]);
							} else {
								varcount += Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[varindex]);
							}
						}

						if (varcount < 3) {
							bamFlag.add("OU");
						}

					} else if (mendFlag.contains("HR")) { // If flag is homozygous recessive
						int varindex = 1 - refMinor;
						boolean homrecB = false;
						if (Integer.parseInt(Arrayreadcount.get(0)[1].split(";")[1 - varindex]) < 5) { //father
							homrecB = true;
						}
						if (Integer.parseInt(Arrayreadcount.get(1)[1].split(";")[1 - varindex]) < 5) { //mother
							homrecB = true;
						}
						if (Integer.parseInt(Arrayreadcount.get(2)[1].split(";")[1 - varindex]) > 5) { //proband
							homrecB = true;
						}
						for (int i = 3; i < family.size(); i++) {
							if(!sibsAffected.get(i - 3)) { //Checks if siblings are affected
								if (Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[1 - varindex]) < 5) {
									homrecB = true;
								}
							} else {
								
								if (Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[1 - varindex]) > 5) {
									homrecB = true;
								}
								
							}
						}

						if (homrecB) {
							bamFlag.add("OU");
						}

					} else if (mendFlag.contains("XL")) { // If flag is x-linked
						int varindex = 1 - refMinor;
						boolean xl = false;
						
						if (Integer.parseInt(Arrayreadcount.get(2)[1].split(";")[1 - varindex]) > 5) {
							xl = true;
						}
						for (int i = 3; i < family.size(); i++) {
							
							if(sibsGender.get(i - 3)) { //Check if siblings' genders are male
							
								if (!sibsAffected.get(i - 3)) { //Checks if siblings are affected
									if (Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[1 - varindex]) < 5) {
										xl = true;
									}
								} else if (sibsAffected.get(i - 3)) {
									if (Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[1 - varindex]) > 5) {
										xl = true;
									}
								}
							}
						}
						if (xl) {
							bamFlag.add("OU");
						}
					}
				}
				
				/// Write to the bam curation output file
				String popStr = set2Str(popFlag);
				String bamStr = set2Str(bamFlag);
				String mendStr = set2Str(mendFlag);
				String miscStr = set2Str(miscFlag);

				String[] fullFlagStr = { popStr, bamStr, mendStr, miscStr };
				bamFilterWriter.println(curateEnd(curLine, naIndex, family, fullFlagStr));

				/// Write to the ROC output file
				// Filter based on flags
				Boolean filterPass = true;

				if (!Collections.disjoint(badFlags, popFlag)) {
					filterPass = false;
				} else if (!Collections.disjoint(badFlags, bamFlag)) {
					filterPass = false;
				} else if (!Collections.disjoint(badFlags, mendFlag)) {
					filterPass = false;
				} else if (!Collections.disjoint(badFlags, miscFlag)) {
					filterPass = false;
				}
				
				if (filterPass) {
					ROCFilterWriter.println(curateEnd(curLine, naIndex, family, fullFlagStr));
				} else if (special && !bamFlag.contains("OU") && !bamFlag.contains("BR")) {
					// if was marked as DM by HGMD and decent quality bam
					ROCFilterWriter.println(curateEnd(curLine, naIndex, family, fullFlagStr));
				}
			}
		}

		/// Compound het filter
		// For each compound het pairing
		for (int j = 0; j < cmpdhetpair.size(); j++) {
			String[] getPair = cmpdhetpair.get(j).printOut(homVarInds, homRefInds, minorField, maxMAF1Inds, Igene,
					mygene, homVarThreshold);

			// If there is a pairing (b/t 2 halfHets), then write to the file)
			if (getPair.length == 2) {

				ArrayList<String[]> cmFlagList = new ArrayList<String[]>();
				String mendCMFlagStr = "CM";

				ArrayList<Boolean> filterPassArr = new ArrayList<Boolean>();

				for (String cmLineStr : getPair) {
					String[] cmLine = cmLineStr.split("\t");

					/// Initialize flag columns
					HashSet<String> bamCMFlags = new HashSet<String>();
					HashSet<String> miscCMFlags = new HashSet<String>();

					ArrayList<String[]> Arrayreadcount = new ArrayList<String[]>();
					ArrayList<String> bamStrArr = new ArrayList<String>();
					int totalCMReads = 0;
					for (int familyInd : family) {
						String allCountStr = cmLine[familyInd + 3];
						String goodCountStr = cmLine[familyInd + 4];

						String[] allCount = allCountStr.split(";");
						String[] goodCount = goodCountStr.split(";");

						/// Extract and flag salvage status
						bamCMFlags.add(allCount[0]);
						bamStrArr.add(allCountStr);

						if (allCount[0].contains("NN") || allCount[0].contains("NP")) {
							bamCMFlags.add("SA"); // Salvaged flag
							bamCMFlags.add("SH"); // Salvaged het flag
						}

						String[] countArr = { allCountStr, goodCountStr };
						Arrayreadcount.add(countArr);

						// Calculate total reads of the variant across the family
						totalCMReads += Integer.parseInt(goodCount[0]) + Integer.parseInt(goodCount[1]);
					}

					/// CM-specific bam curation
					int hetcount = 0;
					double avgread = (double) totalCMReads / (double) family.size();
					double errorVal = Double.parseDouble(cmLine[IerrorVal]);

					// Calculate average number of reads per person in a family
					if (errorVal > 0.01 && avgread < 4) {
						bamCMFlags.add("BR"); // Bad read
					} else if (errorVal > 0.02) {
						bamCMFlags.add("BR"); // Bad read
					}

					// For the parents
					if (avgread < 20) {
						for (int i = 0; i < 2; i++) {
							if (Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[0]) > 2
									&& Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[1]) > 2) {
								hetcount++;
							}
						}

						if (hetcount == 2 || !(Integer.parseInt(Arrayreadcount.get(2)[1].split(";")[0]) > 0
								&& Integer.parseInt(Arrayreadcount.get(2)[1].split(";")[1]) > 0)) {
							bamCMFlags.add("OU");
						}
					} else {
						for (int i = 0; i < 2; i++) {
							if (Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[0]) > 5
									&& Integer.parseInt(Arrayreadcount.get(i)[1].split(";")[1]) > 5) {
								hetcount++;
							}
						}

						if (hetcount == 2 || !(Integer.parseInt(Arrayreadcount.get(2)[1].split(";")[0]) > 0
								&& Integer.parseInt(Arrayreadcount.get(2)[1].split(";")[1]) > 0)) {
							bamCMFlags.add("OU");
						}
					}

					/// Evaluate each halfHet's flags
					String popCMFlagStr = cmpdhetpair.get(j).returnCMPopFlags();
					String bamCMFlagStr = set2Str(bamCMFlags);
					String miscCMFlagStr = set2Str(miscCMFlags);

					String[] CMFlagArr = { popCMFlagStr, bamCMFlagStr, mendCMFlagStr, miscCMFlagStr };
					cmFlagList.add(CMFlagArr);

					Boolean filterPass = true;
					if (badFlags.contains(popCMFlagStr)) {
						filterPass = false;
					} else if (!Collections.disjoint(badFlags, bamCMFlags)) {
						filterPass = false;
					}

					filterPassArr.add(filterPass);
				}

				String[] halfHetLine1 = getPair[0].split("\t");
				String[] halfHetLine2 = getPair[1].split("\t");

				// Re-calculate VMM scores between the pairs (in case the previously matching
				// halfHet was rejected b/c of things like pop.
				double caddHalfHet1 = Double.parseDouble(halfHetLine1[ICADD]);
				double caddHalfHet2 = Double.parseDouble(halfHetLine2[ICADD]);
				double vmmHalfHet1 = Double.parseDouble(halfHetLine1[Ivmmc]);
				double vmmHalfHet2 = Double.parseDouble(halfHetLine2[Ivmmc]);

				if (vmmHalfHet1 != vmmHalfHet2) {
					double newVMM = VMMCalculation(caddHalfHet1, caddHalfHet2);
					// Replace the old VMM score with the new VMM score
					halfHetLine1[Ivmmc] = Double.toString(newVMM);
					halfHetLine2[Ivmmc] = Double.toString(newVMM);
					getPair[0] = String.join("\t", halfHetLine1);
					getPair[1] = String.join("\t", halfHetLine2);
				}

				// Print the halfHets to the output ROC and bam files
				// If no bad flags in either halfHet, then print to ROC output
				for (int i = 0; i < getPair.length; i++) {
					String[] halfHetLine = getPair[i].split("\t");
					/// Print to bam curate file
					bamFilterWriter.println(curateEnd(halfHetLine, naIndex, family, cmFlagList.get(i)));
					/// Print to ROC files
					if (!filterPassArr.contains(false)) {
						ROCFilterWriter.println(curateEnd(halfHetLine, naIndex, family, cmFlagList.get(i)));
					}
				}
			}

		}
		vsData.close();
		bamFilterWriter.close();
		ROCFilterWriter.close();
	}

	/*
	 * Read in headers from user-input header configuration file.
	 */
	public void readHeaders(File headerFile) throws IOException {
		BufferedReader headerData = new BufferedReader(new FileReader(headerFile));
		String Line = headerData.readLine();

		ArrayList<String> unknownHeaders = new ArrayList<String>();

		while (Line != null) {
			/// If Line is not empty or a comment
			if (!Line.startsWith("#") && !Line.isEmpty() && !Line.trim().equals("") && !Line.trim().equals("\n")) {
				String[] headerLine = Line.split("=");
				String categ = headerLine[0]; // Category
				String vsCol = headerLine[1]; // Column name in VS file
				if (categ.equals("Tot_AC") && (headers.indexOf(vsCol) != -1)) {
					totACInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("Tot_AN") && (headers.indexOf(vsCol) != -1)) {
					totANInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("Tot_AF") && (headers.indexOf(vsCol) != -1)) {
					totAFInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("HomRef_Count") && (headers.indexOf(vsCol) != -1)) {
					homRefInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("HomVar_Count") && (headers.indexOf(vsCol) != -1)) {
					homVarInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("Het_Count") && (headers.indexOf(vsCol) != -1)) {
					hetInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("Other_Count") && (headers.indexOf(vsCol) != -1)) {
					othInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("Genotype_Count") && (headers.indexOf(vsCol) != -1)) {
					genotypeInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("HemRef_Count") && (headers.indexOf(vsCol) != -1)) {
					hemRefInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("HemVar_Count") && (headers.indexOf(vsCol) != -1)) {
					hemVarInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("MaxMAF1") && (headers.indexOf(vsCol) != -1)) {
					maxMAF1Inds.add(headers.indexOf(vsCol));
				} else if (categ.equals("MaxMAF2") && (headers.indexOf(vsCol) != -1)) {
					maxMAF2Inds.add(headers.indexOf(vsCol));
				} else if (categ.equals("AFR_AF") && (headers.indexOf(vsCol) != -1)) {
					afrAFInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("EAS_AF") && (headers.indexOf(vsCol) != -1)) {
					easAFInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("AMR_AF") && (headers.indexOf(vsCol) != -1)) {
					amrAFInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("EUR_AF") && headers.indexOf(vsCol) != -1) {
					eurAFInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("SAS_AF") && (headers.indexOf(vsCol) != -1)) {
					sasAFInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("OTH_AF") && (headers.indexOf(vsCol) != -1)) {
					othAFInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("AFR_AC") && (headers.indexOf(vsCol) != -1)) {
					afrACInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("EAS_AC") && (headers.indexOf(vsCol) != -1)) {
					easACInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("AMR_AC") && (headers.indexOf(vsCol) != -1)) {
					amrACInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("EUR_AC") && (headers.indexOf(vsCol) != -1)) {
					eurACInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("SAS_AC") && (headers.indexOf(vsCol) != -1)) {
					sasACInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("OTH_AC") && (headers.indexOf(vsCol) != -1)) {
					othACInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("AFR_AN") && (headers.indexOf(vsCol) != -1)) {
					afrANInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("EAS_AN") && (headers.indexOf(vsCol) != -1)) {
					easANInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("AMR_AN") && (headers.indexOf(vsCol) != -1)) {
					amrANInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("EUR_AN") && (headers.indexOf(vsCol) != -1)) {
					eurANInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("SAS_AN") && (headers.indexOf(vsCol) != -1)) {
					sasANInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("OTH_AN") && (headers.indexOf(vsCol) != -1)) {
					othANInds.add(headers.indexOf(vsCol));
				} else if (categ.equals("Reference_is_minor") && (headers.indexOf(vsCol) != -1)) {
					
					minorField = headers.indexOf(vsCol); /// There should be only one "reference is minor" field
															/// selected.
				} else {
					unknownHeaders.add(Line);
				}
			}
			Line = headerData.readLine();
		}

		/// If there are unknown header categories or unknown headers for the VarSifter
		/// file within the header config file, returns error message.
		if (!unknownHeaders.isEmpty()) {
			JOptionPane.showMessageDialog(new JPanel(),
					"Headers or categories not found for \n" + list2Str(unknownHeaders, "\n"));
			System.gc();
			System.exit(0);
		}

		headerData.close();
	}
	

	//gets the start of the metadata (naPos) and the distance between each of the genotypes
	public int genotypeIndices(String[] lineSplit, int columns) {
		int counter = 0;
		int tempVal = 0;
		int naIndex = -1;
		
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
		
		return naIndex;
		
	}
	

	/*
	 * Retrieve values for the variant based on the header category and the headers
	 * specified by the user in the header configuration file. This is for
	 * integer-specific values.
	 */
	public ArrayList<Integer> retrieveValInt(ArrayList<Integer> indexList) {
		ArrayList<Integer> valList = new ArrayList<Integer>();
		for (int ind : indexList) {
			valList.add(Integer.parseInt(curLine[ind]));
		}
		return valList;
	}

	/*
	 * Retrieve values for the variant based on the header category and the headers
	 * specified by the user in the header configuration file. This is for
	 * boolean-specific values.
	 */
	public ArrayList<Boolean> retrieveValBoolean(ArrayList<Integer> indexList) {
		ArrayList<Boolean> valList = new ArrayList<Boolean>();
		for (int ind : indexList) {
			valList.add(Boolean.parseBoolean(curLine[ind]));
		}
		return valList;
	}

	/*
	 * Retrieve values for the variant based on the header category and the headers
	 * specified by the user in the header configuration file. This is for
	 * double-specific values.
	 */
	public ArrayList<Double> retrieveValDouble(ArrayList<Integer> indexList) {
		ArrayList<Double> valList = new ArrayList<Double>();
		for (int ind : indexList) {
			valList.add(Double.parseDouble(curLine[ind]));
		}
		return valList;
	}
	
	public static void addgene(ArrayList<String> genelist, String gene) {
		if (genelist.indexOf(gene) == -1) {
			genelist.add(gene);
		}
	}

	/*
	 * Append strings of flags (e.g. BamFlags) to the end of each variant's metadata
	 * and before family's genotype metadata.
	 */
	public static String curateEnd(String[] curLine, int naIndex, ArrayList<Integer> family, String[] flagHeaders) {
		String out = curLine[0];
		for (int i = 1; i < naIndex; i++) {
			out += "\t" + curLine[i];
		}
		String flagString = flagHeaders[0];
		for (int i = 1; i < flagHeaders.length; i++) {
			flagString += "\t" + flagHeaders[i];
		}

		out = out + "\t" + flagString;
		for (int i = 0; i < family.size(); i++) {
			out += "\t" + curLine[family.get(i)] + "\t" + curLine[family.get(i) + 1] + "\t"
					+ curLine[family.get(i) + 2];
		}
		return out;
	}

	/*
	 * public static String curateStart(String[] curLine, int naIndex,
	 * ArrayList<Integer> family, String[] flagHeaders) { String out = curLine[0];
	 * for (int i = 1; i < naIndex; i++) { out += "\t" + curLine[i]; } String
	 * flagString = flagHeaders[0]; for (int i = 1; i < flagHeaders.length; i++) {
	 * flagString += "\t" + flagHeaders[i]; }
	 * 
	 * out = flagString + "\t" + out; for (int i = 0; i < family.size(); i++) { out
	 * += "\t" + curLine[family.get(i)] + "\t" + curLine[family.get(i) + 1] + "\t" +
	 * curLine[family.get(i) + 2]; } return out; }
	 */

	/*
	 * Convert string array to tab-delimited string.
	 */
	public static String makestring(String[] curLine) {
		String first = curLine[0];
		for (int i = 1; i < curLine.length; i++) {
			first += "\t" + curLine[i];
		}
		return first;
	}

	/*
	 * Map gene's frequency to genelist hashMap (gene-->count)
	 */
	public static void addhashgene(HashMap<String, Integer> genelist, String gene) {
		if (genelist.containsKey(gene)) {
			genelist.put(gene, genelist.get(gene) + 1);
		} else {
			genelist.put(gene, 1);
		}
	}

	/*
	 * Convert list to semicolon-delimited string
	 */
	public static String list2Str(ArrayList<String> list, String delimiter) {
		StringBuilder sb = new StringBuilder();
		if (!list.isEmpty()) {
			for (String s : list) {
				if (s != null && !s.isEmpty()) {
					sb.append(s);
					sb.append(delimiter);
				}
			}
		}
		return sb.toString();
	}

	/*
	 * Convert list to semicolon-delimited string
	 */
	public static String list2StrInt(ArrayList<Integer> list, String delimiter) {
		StringBuilder sb = new StringBuilder();
		if (!list.isEmpty()) {
			for (Integer i : list) {
				String s = i.toString();
				if (s != null && !s.isEmpty()) {
					sb.append(s);
					sb.append(delimiter);
				}
			}
		}
		return sb.toString();
	}

	/*
	 * Convert set to tab-delimited string
	 */
	public static String set2Str(Set<String> set) {
		StringBuilder sb = new StringBuilder();
		for (String s : set) {
			sb.append(s);
			sb.append(";");
		}
		return sb.toString();
	}

	/*
	 * Calculations for inverse Poisson distribution probability of a given group's
	 * allele frequency, for comparison with assumed MAF threshold of either 1% or
	 * 2%. This is used specifically in the above code in the case that there is a
	 * biallelic distribution (reference and variant) AND that the reference allele
	 * is the minor allele.
	 * 
	 * alleleCount = expected number of occurrences of variant in a pop (lambda), nn
	 * = threshold (ex. 1% or 2%), coverage = total number of events possible.
	 * nn*coverage = k, actual number of events that happened. prob is the
	 * cumulative inverse Poisson distribution chance of seeing less variant alleles
	 * than the average. Thus, 1-prob is the chance of seeing more variant alleles
	 * than the average. If (1-prob)<0.05, then there is a p<0.05 chance of seeing
	 * more variant alleles in the sample than the average.
	 * 
	 * Returning false means there is a p<0.05 chance of seeing more variant alleles
	 * in the sample than the average. Returning true means there is a p>=0.05
	 * chance of seeing more variant alleles in the sample than the average.
	 */
	public boolean invPoissonCalc(double alleleCount, double coverage, double nn) {

		if (coverage > 700) { // if the coverage is this high, we can take the
			// measured frequency at face value
			double rate = ((double) alleleCount) / ((double) coverage);
			double thred = 1.2 * nn;
			return rate < thred; // always hedge above the threshold by a factor of 1.2
		}

		/// For smaller populations, the full inverse Poisson distribution calculation
		/// must be performed. Adapted from Lukas Vlahos's compound het filter.
		double fact = 1; // initialize factorial calculation
		double prob = Math.pow(2.718, -alleleCount) * Math.pow(alleleCount, 0) / fact; // adds the 0 case (k=0)
		for (int i = 1; i <= Math.ceil(coverage * nn); i++) { // loops through cases 1->rejectionthreshold
			fact = fact * i;
			prob += Math.pow(2.718, -alleleCount) * Math.pow(alleleCount, i) / fact; // cumulative inverse Poisson
																						// distribution
		}

		if ((1 - prob) < 0.05) {
			return false;
		}
		return true;
	}

	/*
	 * Implements the inverse Poisson calculation (invPoissonCalc) for variants
	 * where the reference denotes the minor allele. Returning true means
	 * maxMAF<threshold. Returning false means maxMAF>=threshold for 95% confidence.
	 */
	public boolean invPoissonCalcExec(ArrayList<Integer> acList, ArrayList<Integer> anList, double threshold) {
		boolean j = true;
		if (acList.size() != anList.size()) { // Not every AC has a corresponding AN value
			// System.out.println("Uh oh: " + list2StrInt(acList, ",") + ", " +
			// list2StrInt(anList, ";") + "; size: " + acList.size()
			// + ", " + anList.size());
			JOptionPane.showMessageDialog(new JPanel(),
					"Each variant allele count header should have a corresponding \n"
							+ "total allele count header. Please ensure this in the header\n" + "configuration file.");
			System.gc();
			System.exit(0);
		} else {
			for (int i = 0; i < acList.size(); i++) {
				double maj_ac = (double) acList.get(i);
				double an = (double) anList.get(i);
				double min_ac = an - maj_ac;
				j = invPoissonCalc(min_ac, an, threshold);
				if (!j) {
					j = false;
					break;
				} else {
					j = true;
				}
			}
		}
		return j;
	}

	/*
	 * Calculates VMM score. Is mainly used because if the CM pair has non-matching
	 * VMM scores, then it's because the other halfHets with matching VMM scores got
	 * rejected due to filters (e.g. pop, count, etc.). Without this, if the VMM
	 * scores aren't matching, then the assumed VMM score of the two is the
	 * lower...which isn't entirely accurate.
	 */
	public double VMMCalculation(double cadd1, double cadd2) {
		double vmmScore = Math.pow(11, (Math.log(cadd1 + 1) / Math.log(11)) * (Math.log(cadd2 + 1) / Math.log(11)));
		return vmmScore;
	}
	
	public static JFrame canceler() {
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("FBBamCurationROC");
		jframe.setSize(500, 100);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);

		JLabel jtext = new JLabel(
				"Generating the bam curate and ROC filtered files. Use the cancel botton below to abort.",
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
