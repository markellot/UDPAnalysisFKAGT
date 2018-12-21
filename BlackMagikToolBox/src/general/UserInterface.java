package general;

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
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.filechooser.FileNameExtensionFilter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class UserInterface {
	
	//All Variables to be used for the User Interface, most of them are passed into the different modules
	
	//Variables for looking into the header of the VS file
	private int naIndex;
	private List<String> headzList;
	private String[] scroller;
	
	//General Variables
	private File vsFile;
	private String destination;
	private File pedFile;
	private String configPath; 
	private File BamDirectoryFile;
	private String outputDir;
	private HashMap<String, ArrayList<String>> bamChrMap;
	
	//Ethnicity Matcher Specific Variables
	private File sampleFile;
	private String cutoffCriteria;
	private String[] markerPaths;
	
	//Salvage Pathway Specific Variables
	private int minor;
	private int yesSetPop;
	private int[] populationHeaders = new int[7];
	
	//Kayla Kode Specific Variables
	private String strengthPath;
	private String genePath;
	
	//Pedigree
	private int numSibs;
	private ArrayList<Integer> family;
	private ArrayList<Boolean> sibsAffected;
	private ArrayList<Boolean> sibsGender;
	private boolean probandgender;
	private ArrayList<String> familyNames;
	
	//Populations
	private boolean[] popBools;
	private int[] popInts;
	private ArrayList<int[]> populations;
	private ArrayList<int[]> xfilters;
	private ArrayList<int[]> xfilterc;
	private ArrayList<Integer> xfilterExt = new ArrayList<Integer>();
	private int extData;
	private int popExtInt;
	private double rejectionThreshold;
	private int threads;
	private int cutcnt2;
	
	//Nothing needed for Broad Bam File Curator :)
	
	//Variant Exclusion Filter
	private Double[] caddThresh;
	private int homVarThreshold;
	private int hetThreshold;
	private int hemVarThreshold;
	private File ROCHeaderConfig;
	private ArrayList<Integer> genotypeThreshList;
	
	//Confetti Filter
	private int homRefIndex;
	private int homVarIndex;
	private int genotypeIndex;
	private int genotypeThresh;
	private int baseFilter;
	private int baseQualThresh;
	private double badReadRatio;
	
	//CNC Filter
	private int clusterSize;
	private String exonBoundConfig;
	private int CNCExon;
	
	//Constructor basically just calls the entire interface so that it asks everything moment object
	//is initialized
	public UserInterface() throws IOException, InterruptedException {
		
		constructUI();
		
	}
	
	//The main "runner" UI method - this is where all the UI is, and split apart by general variables
	//and modular specific variables
	public void constructUI() throws IOException, InterruptedException {
		
		//User Introduction
		
		//Skipping Descriptions
		int k = JOptionPane.showConfirmDialog(new Frame("Skip Intro"), "Would you like to skip the introduction?", "Skip Intro?",
				JOptionPane.YES_NO_OPTION);
		
		if (k == 1) {
			
			TXTFile.descriptions(); //Description of pipeline
			
		}
		
		//General Inputs
		vsFile = TXTFile.getVSFile(); //Initial VarSifter File
		
		//Checks if VS file is a file and readable
		if(!(vsFile.isFile() && vsFile.canRead())) {
			
			JOptionPane.showMessageDialog(new Frame("Error"),
					"Error\nInput File is not a file/ cannot be read, system exiting.");
			System.exit(1);
			
		}
		
		System.out.println("VSFile Path: " + vsFile.getAbsolutePath());
		
		
		destination = TXTFile.getDestination(); //Gets String for the output folder name 
		//pedFile = TXTFile.getPedigreeFile();
		System.out.println("Output Name: " + destination);
		
		//Gets Config Directory
		File configz = new File("Configs");
		if(!configz.exists()) { //Checks if Config Directory is in same folder as jar
			
			configz = getConfigDirec(); //if not, asks user to locate the directory
			
		}  else { //if so, confirms the config found in directory and allows user to change in case they do not want to use that config direc
			
			int j = JOptionPane.showConfirmDialog(new Frame("Config Confirmer"), "Do you want to use this config Directory?: " + configz.getAbsolutePath(), "Config Confirmer",
					JOptionPane.YES_NO_OPTION);
			
			if (j == 1) {
				
				configz = getConfigDirec();
				
			}
		}
		configPath = configz.getAbsolutePath();
		
		System.out.println("Config Directory Path: " + configPath);
		
		//Gets BAM directory file
		BamDirectoryFile = new File(configPath + "\\BAM_Directory_Config.txt");
		if(!BamDirectoryFile.exists()) {
			
			BamDirectoryFile = TXTFile.getBAMDirectoryFile();
			
		}
		
		//Checks if BAM directory file is a file and readable
		if(!(BamDirectoryFile.isFile() && BamDirectoryFile.canRead())) {
			
			JOptionPane.showMessageDialog(new Frame("Error"),
					"Error\nBam Directory File is not a file / cannot be read, system exiting.");
			System.exit(1);
			
		}
		
		
		System.out.println("Bam Directory File: " + BamDirectoryFile.getAbsolutePath());
		
		
		//For reading the headers and developing the variables needed for some of the User Inputs
		BufferedReader buffy = new BufferedReader(new FileReader(vsFile));
		
		String headLine = buffy.readLine();
		
		String[] headz = (headLine.split("\t"));
 		
		headzList = java.util.Arrays.asList(headLine.split("\t"));
		
		genotypeIndices(headz, headz.length);
		
		String scrollHeadLine = "";
		
		for(int i = 0; i < naIndex; i++) {
			
			scrollHeadLine += headz[i] + "\t";
			
		}
		
		scroller = scrollHeadLine.split("\t");
		
		buffy.close();
		
		VSFileChecker(headzList);
		
		//Start of the modular specific UI
		
		
		//Ethnicity Matcher
		ethnicityMatcherUI();
		
		System.out.println("\nEthnicity Matcher UI");
		System.out.println("Cutoff Criteria: " + cutoffCriteria);
		for(int i = 0; i < markerPaths.length; i++) {
			System.out.println("Marker Paths #" + (i+1) + ": " +  markerPaths[i]);
		}
		
		
		//Salvage Pathway
		salvagePathwayUI();
		
		System.out.println("\nSalvage Pathway UI");
		System.out.println("Minor Allele Index: " + minor);
		System.out.println("Yes Set Population Number: " + yesSetPop);
		if(yesSetPop == 1) {
			for(int i = 0; i < populationHeaders.length; i++) {
				
				System.out.println("Population Headers #" + (i+1) + ": " + populationHeaders[i]);
				
			}
		}
		
		
		//Kayla Kode
		kaylaKodeUI();
		
		//This is the initialization for the supporting method that essentially translates the chr inputs to the
		//BAM files based on the BAMs' formatting
		bamChrMapping();
		
		
		System.out.println("\nKayla Kode UI");
		System.out.println("Strength Config Path: " + strengthPath);
		System.out.println("Gene Config Path: " + genePath);
		System.out.println("Number of siblings: " + numSibs);
		for(int i = 0; i < family.size(); i++) {
			System.out.println("Family Member #" + (i+1) + "'s index: " + family.get(i));
		}
		for(int i = 0; i < sibsAffected.size(); i++) {
			System.out.println("Affected Sib #" + (i+1) + ": " + sibsAffected.get(i));
		}
		for(int i = 0; i < sibsGender.size(); i++) {
			System.out.println("Gender for Sib #" + (i+1) + ": " + sibsGender.get(i));
		}
		System.out.println("Proband Gender: " + probandgender);
		for(int i = 0; i < familyNames.size(); i++) {
			System.out.println("Family Member #" + (i+1) + "'s name: " + familyNames.get(i));
		}
		System.out.println("ifP: " + popBools[0]);
		System.out.println("ifM: " + popBools[1]);
		System.out.println("ifunique: " + popBools[2]);
		System.out.println("ifR: " + popBools[3]);
		System.out.println("chrM: " + popBools[4]);
		System.out.println("NumPops1: " + popInts[0]);
		System.out.println("NumPops2: " + popInts[1]);
		System.out.println("NumPops3: " + popInts[2]);
		System.out.println("NumPops4: " + popInts[3]);
		System.out.println("NumPops5: " + popInts[4]);
		System.out.println("Minor Field Index: " + minor);
		
		for(int i = 0; i < populations.size(); i++) {
			for(int j =0; j < populations.get(i).length; j++) {
				System.out.println("Populations #" + (i+1) + " and index #" + (j+1) + ": " + populations.get(i)[j]);
			}
		}
		
		for(int i = 0; i < xfilters.size(); i++) {
			for(int j =0; j < xfilters.get(i).length; j++) {
				System.out.println("XfilterS #" + (i+1) + "and index #" + (j+1) + ": " + xfilters.get(i)[j]);
			}
		}
		for(int i = 0; i < xfilterc.size(); i++) {
			for(int j =0; j < xfilterc.get(i).length; j++) {
				System.out.println("XfilterC #" + (i+1) + "and index #" + (j+1) + ": " + xfilterc.get(i)[j]);
			}
		}
		System.out.println("Rejection Threshold: " + rejectionThreshold);
		System.out.println("Threads: " + threads);
		System.out.println("CutCNT2: " + cutcnt2);
		
		
		//Nothing needed for BroadBam File Curation :)
		
		//Variant Exclusion Filter
		varExclusionUI();
		
		System.out.println("\nVariant Exclusion UI");
		for(int i = 0; i < caddThresh.length; i++) {
			
			System.out.println("Cadd Threshold #" + i + ": " + caddThresh[i]);
			
		}
		System.out.println("Hom Var Threshold: " + homVarThreshold);
		System.out.println("Het Threshold: " + hetThreshold);
		System.out.println("Hem Var Threshold: " + hemVarThreshold);
		System.out.println("ROC Header Config File: " + ROCHeaderConfig);
		for(int i = 0; i < genotypeThreshList.size(); i++) {
			System.out.println("Genotype Threshold List #" + (i+1) + ": " + genotypeThreshList.get(i));
		}
		
		//Confetti Filter
		confettiUI();
		
		System.out.println("\nConfetti UI");
		System.out.println("Hom Ref Index: " + homRefIndex);
		System.out.println("Hom Var Index: " + homVarIndex);
		System.out.println("Genotype Index: " + genotypeIndex);
		
		//CNC Filter
		cncUI();
		
		System.out.println("\nCNC Filter UI");
		System.out.println("Cluster Size: " + clusterSize);
		
		
	}
	
//////////////////////////////////////////////Modular Specific UI Methods///////////////////////////////////////////////////////////////////////////////////////////////////////
	
	public void ethnicityMatcherUI() {
		
		//sampleFile = TXTFile.getSampleFile(); - No Sample file needed as Sample file is generated by Kayla Kode
		
		cutoffCriteria = TXTFile.getCoverageCutoff();
		
		/// Check if entered criteria cutoff fits requirement
		if (!cutoffCriteria.contains(";") || !cutoffCriteria.replaceAll("\\d", "").equals(";")) {
			JOptionPane.showMessageDialog(new JPanel(), "Cutoff input does not fit requirement. System exiting.");
			System.gc();
			System.exit(0);
		}

		String[] criteria = cutoffCriteria.split(";");
		if (Integer.parseInt(criteria[0]) > 200 || Integer.parseInt(criteria[1]) > 500) {
			JOptionPane.showMessageDialog(new JPanel(), "Cutoff input does not fit requirement. System exiting.");
			System.gc();
			System.exit(0);
		}
		
		
		//Gets all the Marker Configs found in the Configs folder
		markerPaths = new String[13];
		
		markerPaths[0] = configChecker("EASMarker");
		markerPaths[1] = configChecker("AMRMarker");
		markerPaths[2] = configChecker("EURMarker");
		markerPaths[3] = configChecker("AFRMarker");
		markerPaths[4] = configChecker("PairWiseMarker/AMRneu_EUR.txt");
		markerPaths[5] = configChecker("PairWiseMarker/AMRneu_SAS.txt");
		markerPaths[6] = configChecker("PairWiseMarker/EURneu_AMR.txt");
		markerPaths[7] = configChecker("PairWiseMarker/EURneu_SAS.txt");
		markerPaths[8] = configChecker("PairWiseMarker/SASneu_AMR.txt");
		markerPaths[9] = configChecker("PairWiseMarker/SASneu_EUR.txt");
		markerPaths[10] = configChecker("PairWiseMarker/EUR_AMR.txt");
		markerPaths[11] = configChecker("PairWiseMarker/SAS_AMR.txt");
		markerPaths[12] = configChecker("PairWiseMarker/SAS_EUR.txt");
		
		
		
	}
	
	public void salvagePathwayUI() throws IOException {
		
		
		getMinorFieldHeader(scroller); //Gets the header dictating ref is minor in VarSifter file
		
		//Asks user if they want to set the input population fields different from 1K genomes
		yesSetPop = askUserPopulation();
		if (yesSetPop == 1)	{
			getPopulationHeaders(scroller);
		} else {
			
			int popAFFieldIndex = headzList.indexOf("1KGen_Allele_freq");
			int popACFieldIndex = headzList.indexOf("1KGen_Total_Allele_c");
			int afrAFFieldIndex = headzList.indexOf("1KGen_AFR_freq");
			int easAFFieldIndex = headzList.indexOf("1KGen_EAS_freq");
			int amrAFFieldIndex = headzList.indexOf("1KGen_AMR_freq");
			int eurAFFieldIndex = headzList.indexOf("1KGen_EUR_freq");
			int sasAFFieldIndex = headzList.indexOf("1KGen_SAS_freq");
			
			if(popAFFieldIndex == -1 || popACFieldIndex == -1 || afrAFFieldIndex == -1 || easAFFieldIndex == -1 || 
					amrAFFieldIndex == -1 || eurAFFieldIndex == -1 || sasAFFieldIndex == -1) {
				
				JOptionPane.showMessageDialog(new JPanel(),
						"Headers not found for 1K Genomes, please input the population headers manually.");
				yesSetPop = 0;
				getPopulationHeaders(scroller);
				
			}
			
		}
		
	}
	
	public void kaylaKodeUI() throws IOException {
		
		//Gets the pedigree of the family, this will be used for generating the pedfiles and sample ids for the entire pipeline
		getPedigree(headzList);
		
		//Gets other information required for KaylaKode including cutoff numbers, raw, count, and filter data
		//This also shifts based on if proband is a male or a female for X-Linked conditions
		getPopulations(scroller);
		
		
		if(probandgender) { //If Proband is a male, ask if a user wants to use external databases to filter for X-linked conditions
			
			extData = askExtDatabase();
			
			if(extData == 0) {
				
				getExternalHemiHeader(scroller); //Gets headers from these external databases
				
			}
		}
		
		//Gets configs required for KaylaKode
		strengthPath = configChecker("Strength_Config.txt");
		genePath = configChecker("Gene_Config.txt");
		
	}
	
	public void varExclusionUI() throws IOException {
		
		caddThresh = delThresholdCalc(); //Inputs the CADD threshold for filtering variants
		getPopulationGenotypes(); //Input for thresholds on Hom Var, Het, and Hem Vars found in other people in the respective column found in the ROC Header Config
		
		//Gets Header Config for Var Exclusion Filter
		ROCHeaderConfig = new File(configChecker("BamROC_HeaderConfig.txt"));
		
		BufferedReader headerData = new BufferedReader(new FileReader(ROCHeaderConfig));
		String Line = headerData.readLine();
		
		ArrayList<Integer> genoInd = new ArrayList<Integer>();
		
		ArrayList<String> unknownHeaders = new ArrayList<String>();
		
		//Prechecks Config file for genotype headers or unknown headers
		while (Line != null) {
			/// If Line is not empty or a comment
			if (!Line.startsWith("#") && !Line.isEmpty() && !Line.trim().equals("") && !Line.trim().equals("\n")) {
				String[] headerLine = Line.split("=");
				String categ = headerLine[0]; // Category
				String vsCol = headerLine[1]; // Column name in VS file
				if (checkROCConfigHeaders(categ, vsCol) == 1) { //If Config line is a genotype header needed for UI genotyper thresholds
					
					genoInd.add(headzList.indexOf(vsCol));
				} else if (checkROCConfigHeaders(categ, vsCol) == 2) { //If Config Line is not found in VS file
					
					unknownHeaders.add(Line);
					
				}
				
				
			}
			
			Line = headerData.readLine();
		}
		
		
		/// If there are unknown header categories or unknown headers for the VarSifter
		/// file within the header config file, returns error message.
		if (!unknownHeaders.isEmpty()) {
			
			StringBuilder sb = new StringBuilder();
			if (!unknownHeaders.isEmpty()) {
				for (String s : unknownHeaders) {
					if (s != null && !s.isEmpty()) {
						sb.append(s);
						sb.append("\n");
					}
				}
			}
			
			JOptionPane.showMessageDialog(new JPanel(),
					"Headers or categories not found for \n" + sb.toString());
			System.gc();
			System.exit(0);
		}
		
		/// Input minimum genotype count threshold for each genotype found in ROC config list
		genotypeThreshList = new ArrayList<Integer>();
		for (int genoHeaderInd : genoInd) {
			genotypeThreshList.add(genoCountThreshold(genoHeaderInd, headzList));
		}
				
		
		headerData.close();
		
		
	}
	
	public void confettiUI() throws IOException {
		
		getHeaders(scroller); //Gets all information required for confetti filter
		
		
	}
	
	public void cncUI() {
		
		clusterSize = clusterMin(); //Inputs a Cluster Size input for the cnc filter - this will affect the bounds
									//of the CNC the filter looks for
		
		CNCExon = askCNCExonOnly(); //Asks if we want to filter based on if CNC is in an exonic region
		
		if(CNCExon == 0) {
			exonBoundConfig = configChecker("Exon_Boundary_Config.txt");
		}
	}
	
	
	//////////////////////////////////////////////Supporting UI Methods///////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Checks VarSifter file headers to see if the headers used throughout the pipeline are missing
	public void VSFileChecker(List<String> headers) {
		int chrIndex = headers.indexOf("Chr"); // finds indices of the five headers
		// to the left
		int posIndex = headers.indexOf("VCF_Position"); // if header does not exist,
				// String.indexOf returns
		int refIndex = headers.indexOf("VCF_Ref"); // a -1
		int varIndex = headers.indexOf("VCF_Var");
		int mutTypeIndex = headers.indexOf("muttype");
		int typeIndex = headers.indexOf("type");
		int indexIndex = headers.indexOf("Index");
		int LeftIndex = headers.indexOf("LeftFlank");
		int RightIndex = headers.indexOf("RightFlank");
		int refAlleleIndex = headers.indexOf("ref_allele");
		int varAlleleIndex = headers.indexOf("var_allele");
		int geneIndex = headers.indexOf("Gene_name");
		
		int predIndex = headers.indexOf("Prediction");
		int phredIndex = headers.indexOf("PHRED");
		int hgmdIndex = headers.indexOf("HGMDtags");
		int ePhredIndex = headers.indexOf("EPHRED");
		int epcPhredIndex = headers.indexOf("EPCPHRED");
		int mcapIndex = headers.indexOf("M_CAPv1.0");
		
		
		boolean indexError = false;
		
		String errorMessage = "";
		if (chrIndex == -1) {
			errorMessage += "The header 'Chr' appears no where in the header line of the VarSifter file.\n"; // missing 'Chr' header
			indexError = true;
		}
		if (posIndex == -1) {
			errorMessage += "The header 'VCF_Position' appears no where in the header line of the VarSifter file.\n"; // missing 'VCF_Position' header
			indexError = true;
		}
		if (refIndex == -1) {
			errorMessage += "The header 'VCF_Ref' appears no where in the header line of the VarSifter file.\n"; // missing 'VCF_Ref' header
			indexError = true;
		}
		if (varIndex == -1) {
			errorMessage += "The header 'VCF_Var' appears no where in the header line of the VarSifter file.\n"; // missing 'VCF_Var' header
			indexError = true;
		}
		if (mutTypeIndex == -1) {
			errorMessage += "The header 'muttype' appears no where in the header line of the VarSifter file.\n"; // missing 'muttype' header
			indexError = true;
		}
		if (typeIndex == -1) {
			errorMessage += "The header 'type appears no where in the header line of the VarSifter file.\n"; // missing 'type' header
			indexError = true;
		}
		if (indexIndex == -1) {
			errorMessage += "The header 'Index' appears no where in the header line of the VarSifter file.\n"; // missing 'Index' header
			indexError = true;
		}
		if (LeftIndex == -1) {
			errorMessage += "The header 'LeftFlank' appears no where in the header line of the VarSifter file.\n"; // missing 'LeftFlank' header
			indexError = true;
		}
		if (RightIndex == -1) {
			errorMessage += "The header 'RightFlank' appears no where in the header line of the VarSifter file.\n"; // missing 'RightFlank' header
			indexError = true;
		}
		if (refAlleleIndex == -1) {
			errorMessage += "The header 'ref_allele' appears no where in the header line of the VarSifter file.\n"; // missing 'ref_allele' header
			indexError = true;
		}
		if (predIndex == -1) {
			errorMessage += "The header 'Prediction' appears no where in the header line of the VarSifter file.\nPlease edit the VarSifter file accordingly or edit in pseudomaster edit files as referenced in manual.\n"; // missing 'Prediction' header
			indexError = true;
		}
		if (varAlleleIndex == -1) {
			errorMessage += "The header 'var_allele' appears no where in the header line of the VarSifter file.\n"; // missing 'var_allele' header
			indexError = true;
		}
		if (geneIndex == -1) {
			errorMessage += "The header 'Gene_name' appears no where in the header line of the VarSifter file.\n"; // missing 'Gene_name' header
			indexError = true;
		}
		if (phredIndex == -1) {
			errorMessage += "The header 'PHRED' appears no where in the header line of the VarSifter file.\nPlease edit the VarSifter file accordingly or edit in pseudomaster edit files as referenced in manual.\n"; // missing 'PHRED' header
			indexError = true;
		}
		if (hgmdIndex == -1) {
			errorMessage += "The header 'HGMDtags' appears no where in the header line of the VarSifter file.\nPlease edit the VarSifter file accordingly or edit in pseudomaster edit files as referenced in manual.\n"; // missing 'HGMDtags' header
			indexError = true;
		}
		if (ePhredIndex == -1) {
			errorMessage += "The header 'EPHRED' appears no where in the header line of the VarSifter file.\nPlease edit the VarSifter file accordingly or edit in pseudomaster edit files as referenced in manual.\n"; // missing 'EPHRED' header
			indexError = true;
		}
		if (epcPhredIndex == -1) {
			errorMessage += "The header 'EPCPHRED' appears no where in the header line of the VarSifter file.\nPlease edit the VarSifter file accordingly or edit in pseudomaster edit files as referenced in manual.\n"; // missing 'EPCPHRED' header
			indexError = true;
		}
		if (mcapIndex == -1) {
			errorMessage += "The header 'M_CAPv1.0' appears no where in the header line of the VarSifter file.\nPlease edit the VarSifter file accordingly or edit in pseudomaster edit files as referenced in manual.\n"; // missing 'M_CAPv1.0' header
			indexError = true;
		}
		
		if (indexError) { // if any of the above if statements are true,
			errorMessage += "As stated above. One or more of the necessary headers do not appear" // then
								 /*the
								 header
								 does
								 not
								 appear
								 within
								 the*/
					+ " in your VarSifter file. \nPlease ensure these headers are in " // header
					 /*line.
					 Thus,
					 we
					 print
					 the
					 missing
					 header(s)*/
					+ "your VarSifter file prior to running this program. These headers " // along
						 /*with
						 the
						 statement
						 below
						 and
						 the
						 program
						 exits.*/
					+ "can be added by running the annotation utility.";
			JOptionPane.showMessageDialog(new Frame("Error Message"),
					errorMessage);
			System.exit(0);
		}
		
	}
	
	//creates a FileChooser that allow the selection of  directory containing all the config files
	public File getConfigDirec() {
		JFileChooser browseTo = new JFileChooser(); //creates a file chooser object
		browseTo.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY); //limits selectable objects to directories
		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please select the directory containing the " +
				"config files.");
		int returnVal = browseTo.showOpenDialog(new JPanel()); //opens a window for user to browse to directory
		
	    if(returnVal == JFileChooser.APPROVE_OPTION) {
	    	return browseTo.getSelectedFile(); //returns selected directory
	    } else {
			int i = TXTFile.wantToContinue("No directory selected."); 
			if (i == 1) {
				System.exit(0);
			}
			return this.getConfigDirec();
			
		}	
	}
	
	//Checks config directory for the config based on the string inputed, if it does not exist in file, it will look for it
	//Changes mode of looking for config based on if it is a directory or a .txt file
	public String configChecker(String config) {
		
		File configz = new File(configPath + "\\" + config);
		
		if(!configz.exists()) {
			
			configz = getConfigFile(configz.getAbsolutePath(), config);
			
		}
		
		if(configz.getAbsolutePath().endsWith(".txt") && !(configz.isFile() && configz.canRead())) {
			
			JOptionPane.showMessageDialog(new Frame("Error"),
					"Error\n" + config + " config file is not a file/ cannot be read, system exiting.");
			System.exit(1);
			
		}
		
		return configz.getAbsolutePath();
		
	}
	

	//creates a FileChooser that allows the selection of an input marker config File
	//Changes 'mode' based on if config is a .txt file or a directory file
	public File getConfigFile(String path, String name) {
		JFileChooser browseTo = new JFileChooser(); //creates a file chooser object
		
		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Config with the following path has not been found:\n\n" + path + "\n\nPlease select the edit config with the name:\n\n\"" + name + "\"");
		
		if(!path.endsWith(".txt")) {
			
			browseTo.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY); //limits selectable objects to directories
			
		} else {
			
			FileNameExtensionFilter filter = new FileNameExtensionFilter(
					"Text Files", "txt");
			browseTo.setFileFilter(filter);
		}
		
		int returnVal = browseTo.showOpenDialog(new JPanel()); //opens a window for user to browse to config
		
		if(returnVal == JFileChooser.APPROVE_OPTION) {
			return browseTo.getSelectedFile(); //gets the path of the selected file
			
		} else {
			int i = TXTFile.wantToContinue("No Config selected."); 
				if (i == 1) {
					System.exit(0);
				}
				return this.getConfigFile(path, name);
		}
	}
	
	//Generates a ped file based on output from family names from KaylaKode getPedigree() method
	public void generatePedFile(String dest) throws IOException {
		
		File ped = new File(dest +  "\\" + "PedFile" + familyNames.get(2) + ".txt");
		
		PrintWriter printy = new PrintWriter(ped);
		
		String printedLion = "";
		
		for(int i = 0; i < familyNames.size()-1; i++) {
			
			printedLion += familyNames.get(i) + "\t";
			
		}
		
		printedLion += familyNames.get(familyNames.size() - 1);
		
		printy.println(printedLion);
		
		printy.close();
		
		pedFile = ped;
		

		System.out.println("Ped File Path: " + pedFile.getAbsolutePath());
		
		
	}
	
	//Generates a sample file based on output from family names from KaylaKode getPedigree() method
	public void generateSampleFile(String dest) throws IOException {
		
		File samp = new File(dest + "\\" + "SampleFile" + familyNames.get(2) + ".txt");
		
		PrintWriter printy = new PrintWriter(samp);
		
		
		for(int i = 0; i < familyNames.size(); i++) {
			
			printy.println(familyNames.get(i));			
		}
		
		
		printy.close();
		
		sampleFile = samp;
		System.out.println("Sample File Path: " + sampleFile.getAbsolutePath());
	}
	
	//Method that returns a hashmap of <Location of Bam file, list of headers found from BAM file>
	//Also initializes the lists used to translate from chr formatting found within the config file
	public HashMap<String, ArrayList<String>> getbamChrMap() throws IOException {
		
		bamChrMapping();
		BamChrChanger.initializeConfig(configPath, outputDir);
		
		return bamChrMap;
		
	}
	
	//Generates HashMap of <Location of Bam file, list of headers found from BAM file>
	public void bamChrMapping() throws IOException {
		
		HashMap<String, String> BamMap = new HashMap<String, String>();
		/// Path to the BAM file directory file
		//String bamlo = path + "/EthnicityMatching_Bam.txt"; This is hard-coded BAM Dir File
		FindBam.InitializeBam(BamDirectoryFile.getAbsolutePath(), BamMap);
		
		bamChrMap = new HashMap<String, ArrayList<String>>();
		
		for(int i = 0; i < familyNames.size(); i++) {
			
			String bamLoc = FindBam.MakeBamString(familyNames.get(i), BamMap);
			
			bamChrMap.put(bamLoc, getBamChrHeaders(bamLoc));
			
		}
		
		
		
		
	}
	
	//Reads BAM file to get the chrs
	public ArrayList<String> getBamChrHeaders(String bamLoc) throws IOException {
		
		//Creates SAMReader to look for headers
		File bamFile = new File(bamLoc);
		
		SamReader samTheReader = 
		SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile);     
		
		SAMRecordIterator idioterator = samTheReader.iterator();
		
		String headz = "";
		
		if(idioterator.hasNext()) {
			
			SAMRecord recorder = idioterator.next();    
			
			headz = (recorder.getHeader().getTextHeader());
			
		} else {
			

			JOptionPane.showMessageDialog(new Frame("Error"), "BAM Header could not be read. Program will now exit");
			System.exit(0);
			
			
		}
		
		
		
		String[] headzArr = headz.split("\n");
		
		ArrayList<String> headerz = new ArrayList<String>();
		
		ArrayList<String> SN = new ArrayList<String>();
		
		int counter = 0;
		
		//@SQ lines contain the information needed for the chr names
		//SN: dictates the chr names
		for(int i = 0; i < headzArr.length; i++) {
			
			if(headzArr[i].startsWith("@SQ")) {
				
				headerz.add(headzArr[i]);
				
				
				String[] sqLine = headerz.get(counter).split("\t");
				
				for(int j = 0; j < sqLine.length; j++) {
					
					if(sqLine[j].startsWith("SN:")) {
						
						SN.add(sqLine[j].substring(sqLine[j].indexOf(":") + 1));
						
					}
					
				}
				
				
				
				counter++;
			}
			
		}
		
		
		//r.close();
		idioterator.close();
		
		samTheReader.close();
		//sr.close();
		
		return(SN);
		
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
						//naDistanceApart = counter - naPos;
					}
				tempVal++;			
				}
				counter++;
			}
		}
	}
	
	//Sets the output directory so that Ped and sample files can be located into there
	public void setOutputDir(String dir) {
		
		outputDir = dir;
		
		
	}
	

	/*
	 * Asks for user input on which field they want to denote if the reference
	 * represents the minor allele. This is derived from the KaylaKode GUI code. The
	 * output file, input VS file, and pedigree files are included as inputs for the
	 * cancellation option, so that the PrintWriter and BufferedReader streams are
	 * successfully closed if the user chooses to abort the job.
	 */
	public void getMinorFieldHeader(String[] scroller) {
		/// Create an initial frame, panel, dialog, and gridBagConstraints
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("Input Panel");
		jframe.setSize(2000, 20000);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		final JDialog diag = new JDialog(jframe, "Input Panel", true);

		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);

		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(2, 2, 2, 2);
		gbc.weightx = 1.0;
		gbc.gridwidth = 2;

		/// Request user input for field denoting if reference is minor allele.
		JLabel refIsMinorAsk = new JLabel("Enter minor allele filter (1 denotes reference being minor): ");
		gbc.gridx = 0;
		gbc.gridy = 0;
		gbc.gridwidth = 1;
		jpanel.add(refIsMinorAsk, gbc);

		JComboBox msglist10 = new JComboBox(scroller);
		msglist10.setEditable(true);

		msglist10.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				minor = msglist10.getSelectedIndex();
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 0;
		jpanel.add(msglist10, gbc);

		/// Makes the submit button and ensures all fields are full.
		JButton okButton = new JButton("Submit");
		okButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				boolean temp = true;
				if (minor == 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please input minor allele filter.");
					temp = false;
				}

				if (temp) { /// if fields are full, get rid of the panel and proceed
					diag.dispose();
				}
			}
		});

		gbc.gridx = 2;
		gbc.gridy = 0;
		jpanel.add(okButton, gbc);
		
		/// Gives user option to cancel inputting info.
		JButton cancelButton = new JButton("Cancel");
		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				jframe.dispose();
				System.gc();
				System.exit(0);

			}
		});

		gbc.gridx = 2;
		gbc.gridy = 1;
		jpanel.add(cancelButton, gbc);

		/// Makes the panel appear
		diag.getContentPane().add(jpanel);
		diag.pack();
		diag.setLocationRelativeTo(null);
		diag.setVisible(true);
	}

	/*
	 * Ask user if they want to use the default population allele frequency fields
	 * (set in the initializer), or choose their own.
	 */
	public int askUserPopulation() {
		int j = JOptionPane.showConfirmDialog(new Frame("Input prompt"),
				"Would you like to use the default allele frequency fields (1000 Genomes), \n"
						+ "or set your own? There must be fields for the total population \n "
						+ "allele frequency and count, and allele frequency fields for AFR, \n" 
						+ "EAS, AMR, EUR, and SAS ethnic subpopulations.",
				"User Prompt", JOptionPane.YES_NO_OPTION);
		/// 0 = yes. 1 = no.
		return j; // User wants to use default
	}

	/*
	 * Ask user for input in population frequency headers (AFR, EAS, AMR< EUR, SAS).
	 * This is derived from the KaylaKode GUI code. The output file, input VS file,
	 * and pedigree files are included as inputs for the cancellation option, so
	 * that the PrintWriter and BufferedReader streams are successfully closed if
	 * the user chooses to abort the job.
	 */
	public void getPopulationHeaders(String[] scroller) {
		// Create an initial frame, panel, dialog, and gridBagConstraints
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("Input Panel");
		jframe.setSize(2000, 20000);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		final JDialog diag = new JDialog(jframe, "Input Panel", true);

		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);

		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(2, 2, 2, 2);
		gbc.weightx = 1.0;
		gbc.gridwidth = 2;
		
		
		/// Reset indices to 0
		/*
		popAFFieldIndex = 0;
		popACFieldIndex = 0;
		afrAFFieldIndex = 0;
		easAFFieldIndex = 0;
		amrAFFieldIndex = 0;
		eurAFFieldIndex = 0;
		sasAFFieldIndex = 0;
		 */
		
		populationHeaders[0] = 0; //popAFFieldIndex
		populationHeaders[1] = 0; //popACFieldIndex
		populationHeaders[2] = 0; //afrAFFieldIndex
		populationHeaders[3] = 0; //easAFFieldIndex
		populationHeaders[4] = 0; //amrAFFieldIndex
		populationHeaders[5] = 0; //eurAFFieldIndex
		populationHeaders[6] = 0; //sasAFFieldIndex
		
		/// Requests user input for field denoting population allele frequency
		JLabel popAlleleFreqAsk = new JLabel(
				"Enter field for population variant allele frequency (e.g. 1KGen_Allele_freq): ");
		gbc.gridx = 0;
		gbc.gridy = 0;
		gbc.gridwidth = 1;
		jpanel.add(popAlleleFreqAsk, gbc);

		JComboBox msglist10 = new JComboBox(scroller);
		msglist10.setEditable(true);

		msglist10.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				//popAFFieldIndex = msglist10.getSelectedIndex();
				populationHeaders[0] = msglist10.getSelectedIndex();
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 0;
		jpanel.add(msglist10, gbc);

		/// Requests user input for field denoting population allele count
		JLabel popAlleleCountAsk = new JLabel("Enter field for population variant allele count (e.g. 1KGen_Alt_Allele_c): ");
		gbc.gridx = 0;
		gbc.gridy = 1;
		gbc.gridwidth = 1;
		jpanel.add(popAlleleCountAsk, gbc);

		JComboBox msglist11 = new JComboBox(scroller);
		msglist11.setEditable(true);

		msglist11.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				//popACFieldIndex = msglist11.getSelectedIndex();
				populationHeaders[1] = msglist11.getSelectedIndex();
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 1;
		jpanel.add(msglist11, gbc);

		/// Requests user input for field denoting African allele frequency
		JLabel afrAFAsk = new JLabel("Enter field for African (AFR) allele frequency: ");
		gbc.gridx = 0;
		gbc.gridy = 2;
		gbc.gridwidth = 1;
		jpanel.add(afrAFAsk, gbc);

		JComboBox msglist12 = new JComboBox(scroller);
		msglist12.setEditable(true);

		msglist12.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				//afrAFFieldIndex = msglist12.getSelectedIndex();
				populationHeaders[2] = msglist12.getSelectedIndex();
				
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 2;
		jpanel.add(msglist12, gbc);
		
		/// Requests user input for field denoting East Asian allele frequency
		JLabel easAFAsk = new JLabel("Enter field for East Asian (EAS) allele frequency: ");
		gbc.gridx = 0;
		gbc.gridy = 3;
		gbc.gridwidth = 1;
		jpanel.add(easAFAsk, gbc);

		JComboBox msglist13 = new JComboBox(scroller);
		msglist13.setEditable(true);

		msglist13.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				//easAFFieldIndex = msglist13.getSelectedIndex();
				populationHeaders[3] = msglist13.getSelectedIndex();
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 3;
		jpanel.add(msglist13, gbc);
		
		/// Requests user input for field denoting American/Hispanic allele frequency
		JLabel amrAFAsk = new JLabel("Enter field for American/Hispanic (AMR) allele frequency: ");
		gbc.gridx = 0;
		gbc.gridy = 4;
		gbc.gridwidth = 1;
		jpanel.add(amrAFAsk, gbc);
		
		JComboBox msglist14 = new JComboBox(scroller);
		msglist14.setEditable(true);

		msglist14.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				//amrAFFieldIndex = msglist14.getSelectedIndex();
				populationHeaders[4] = msglist14.getSelectedIndex();
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 4;
		jpanel.add(msglist14, gbc);

		/// Requests user input for field denoting European allele frequency
		JLabel eurAFAsk = new JLabel("Enter field for European (EUR) allele frequency: ");
		gbc.gridx = 0;
		gbc.gridy = 5;
		gbc.gridwidth = 1;
		jpanel.add(eurAFAsk, gbc);

		JComboBox msglist15 = new JComboBox(scroller);
		msglist15.setEditable(true);

		msglist15.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				//eurAFFieldIndex = msglist15.getSelectedIndex();
				populationHeaders[5] = msglist15.getSelectedIndex();
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 5;
		jpanel.add(msglist15, gbc);

		/// Requests user input for field denoting South Asian allele frequency
		JLabel sasAFAsk = new JLabel("Enter field for South Asian (SAS) allele frequency: ");
		gbc.gridx = 0;
		gbc.gridy = 6;
		gbc.gridwidth = 1;
		jpanel.add(sasAFAsk, gbc);

		JComboBox msglist16 = new JComboBox(scroller);
		msglist16.setEditable(true);

		msglist16.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				//sasAFFieldIndex = msglist16.getSelectedIndex();
				populationHeaders[6] = msglist16.getSelectedIndex();
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 6;
		jpanel.add(msglist16, gbc);

		/// Makes the submit button and ensures all fields are full.
		JButton okButton = new JButton("Submit");
		okButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				boolean temp = true;
				if (populationHeaders[0] == 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"),
							"Please input population variant allele frequency.");
					temp = false;
				}
				if (populationHeaders[1] == 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"),
							"Please input population variant allele count.");
					temp = false;
				}
				//if (afrAFFieldIndex == 0) {
				if (populationHeaders[2] == 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"),
							"Please input African allele frequency field.");
					temp = false;
				}

				//if (easAFFieldIndex == 0) {
				if (populationHeaders[3] == 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"),
							"Please input East Asian allele frequency field.");
					temp = false;
				}

				//if (amrAFFieldIndex == 0) {
				if (populationHeaders[4] == 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"),
							"Please input American/Hispanic allele frequency field.");
					temp = false;
				}

				//if (eurAFFieldIndex == 0) {
				if (populationHeaders[5] == 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"),
							"Please input European allele frequency field.");
					temp = false;
				}

				//if (sasAFFieldIndex == 0) {
				if (populationHeaders[6] == 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"),
							"Please input South Asian allele frequency field.");
					temp = false;
				}

				if (temp) { /// if fields are full, get rid of the panel and proceed
					diag.dispose();
				}
			}
		});

		gbc.gridx = 4;
		gbc.gridy = 0;
		jpanel.add(okButton, gbc);

		/// Gives user option to cancel inputting info.
		JButton cancelButton = new JButton("Cancel");
		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
//				try {
//					output.close();
//					inputVS.close();
//					pedReader.close();
//				} catch (IOException e1) {
//
//				}
				jframe.dispose();
				System.gc();
				System.exit(0);

			}
		});
		gbc.gridx = 4;
		gbc.gridy = 1;
		jpanel.add(cancelButton, gbc);
		
		// makes the panel appear
		diag.getContentPane().add(jpanel);
		diag.pack();
		diag.setLocationRelativeTo(null);
		diag.setVisible(true);
		diag.setAlwaysOnTop(true);
	}
	
	//This gets the Pedigree information of the family (Father, Mother, Proband, Siblings)
	//Also gets gender of proband, and siblings as well as affected status of siblings
	public void getPedigree(List<String> headers) {
		// this is almost entirely UI garbage - you can ignore JLabels in particular
		// during bug fixes

		// initializes the global variables
		numSibs = 0;
		family = new ArrayList<Integer>(); // this will store the index of the genotypes of family members
		sibsAffected = new ArrayList<Boolean>(); // for each sib, stores whether the sib is affected (true) or
		familyNames = new ArrayList<String>();
		
		// unaffected (false)
		sibsGender = new ArrayList<Boolean>();
		for (int i = 0; i < 3; i++) {
			family.add(null);       // adds null in the father/mother/proband positions to start - must be
			familyNames.add("");	// overwritten by user
		}

		// creates the initial frame, panel, dialog, and gbc
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("Input Panel");
		jframe.setSize(300, 400);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		final JDialog diag = new JDialog(jframe, "Input Panel", true);

		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);

		// each family data field will let the user type in a string - if that string is
		// not found in the header...
		// the program will notify the user and wipe the field - the user must then
		// reenter

		// adding the original pedigree field labels and text boxes
		JLabel pedDesc = new JLabel("Please enter the full ID number of the family members.");
		gbc.gridx = 0;
		gbc.gridy = 0;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(2, 2, 2, 2);
		gbc.weightx = 1.0;
		gbc.gridwidth = 5;
		jpanel.add(pedDesc, gbc);

		JLabel father = new JLabel("Father:");
		gbc.gridy = 1;
		gbc.gridwidth = 1;
		jpanel.add(father, gbc);

		// asks for the father's label
		JTextField fatherField = new JTextField(20);
		fatherField.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent e) {
			}

			public void focusLost(FocusEvent e) {
				if (headers.indexOf(fatherField.getText() + ".NA") == -1) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"),
							"The string entered could not be found. Please try again.");
					fatherField.setText("");
				} else {
					family.set(0, headers.indexOf(fatherField.getText() + ".NA"));
					if(fatherField.getText() != null || fatherField.getText() != "") {
						familyNames.set(0, fatherField.getText());
					}
				}
			}
		});
		gbc.gridx = 1;
		jpanel.add(fatherField, gbc);

		JLabel mother = new JLabel("Mother:");
		gbc.gridx = 0;
		gbc.gridy = 2;
		jpanel.add(mother, gbc);

		// asks for the mother's label
		JTextField motherField = new JTextField(20);
		motherField.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent e) {
			}

			public void focusLost(FocusEvent e) {
				if (headers.indexOf(motherField.getText() + ".NA") == -1) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"),
							"The string entered could not be found. Please try again.");
					motherField.setText("");
				} else {
					family.set(1, headers.indexOf(motherField.getText() + ".NA"));
					if(motherField.getText() != null || motherField.getText() != "") {
						familyNames.set(1, motherField.getText());
					}
				}
			}
		});
		gbc.gridx = 1;
		jpanel.add(motherField, gbc);

		JLabel proband = new JLabel("ProBand:");
		gbc.gridx = 0;
		gbc.gridy = 3;
		jpanel.add(proband, gbc);

		// asks for the proband's label
		JTextField pbField = new JTextField(20);
		pbField.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent e) {
			}

			public void focusLost(FocusEvent e) {
				if (headers.indexOf(pbField.getText() + ".NA") == -1) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"),
							"The string entered could not be found. Please try again.");
					pbField.setText("");
				} else {
					family.set(2, headers.indexOf(pbField.getText() + ".NA"));
					if(pbField.getText() != null || pbField.getText() != "") {
						familyNames.set(2, pbField.getText());
					}
				}
			}
		});
		gbc.gridx = 1;
		jpanel.add(pbField, gbc);

		probandgender = false;
		JLabel Mlabel = new JLabel("Male:");
		gbc.gridx = 2;
		gbc.gridy = 3;
		jpanel.add(Mlabel, gbc);

		// check box for the user to mark if the sbling is affected
		JCheckBox malabel = new JCheckBox();
		malabel.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				probandgender = !probandgender;
			}
		});
		gbc.gridx = 3;
		gbc.gridy = 3;
		jpanel.add(malabel, gbc);
		
		// make the adder button - this will let the user add sib
		JButton addButton = new JButton("Add Sib");
		addButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				family.add(null); // need to add a new position to the family array list
				familyNames.add("");
				sibsAffected.add(false); // and a new position to the sib affected array list
				sibsGender.add(false);
				
				// just a label
				JLabel tempLabel = new JLabel("Sib " + Integer.toString(numSibs + 1) + ":");
				gbc.gridx = 0;
				gbc.gridy = numSibs + 4;
				jpanel.add(tempLabel, gbc);

				// functions in the same way as the fahter/mother/proband fields, and will also
				// check the header
				JTextField tempField = new JTextField(20);
				tempField.addFocusListener(new FocusListener() {
					final int x = numSibs;

					public void focusGained(FocusEvent e) {
					}

					public void focusLost(FocusEvent e) {
						if (headers.indexOf(tempField.getText() + ".NA") == -1) {
							JOptionPane.showMessageDialog(new Frame("Input prompt"),
									"The string entered could not be found. Please try again.");
							tempField.setText("");
						} else {
							family.set(x + 3, headers.indexOf(tempField.getText() + ".NA"));
							
							if(tempField.getText() != null || tempField.getText() != "") {
								familyNames.set((x+3), tempField.getText());
							}
						}
					}
				});
				gbc.gridx = 1;
				gbc.gridy = numSibs + 4;
				jpanel.add(tempField, gbc);

				JLabel nlabel = new JLabel("Male:");
				gbc.gridx = 4;
				jpanel.add(nlabel, gbc);

				// check box for the user to mark if the sbling is affected
				JCheckBox nnlabel = new JCheckBox();
				nnlabel.addActionListener(new ActionListener() {
					final int z = numSibs;

					public void actionPerformed(ActionEvent e) {
						sibsGender.set(z, !sibsGender.get(z));
					}
				});
				gbc.gridx = 5;
				gbc.gridy = numSibs + 4;
				jpanel.add(nnlabel, gbc);
				// just a label

				JLabel tempLabel2 = new JLabel("Affected:");
				gbc.gridx = 2;
				jpanel.add(tempLabel2, gbc);

				// check box for the user to mark if the sibling is affected
				JCheckBox tempBox = new JCheckBox();
				tempBox.addActionListener(new ActionListener() {
					final int y = numSibs;

					public void actionPerformed(ActionEvent e) {
						sibsAffected.set(y, !sibsAffected.get(y)); // every time it is (un)checked, the corresponding
																	// position in sibsaffected is reveresed
					}
				});
				gbc.gridx = 3;
				jpanel.add(tempBox, gbc);

				numSibs++; // numSibs is increased

				diag.pack();
				diag.validate();
			}
		});
		gbc.gridx = 4;
		gbc.gridy = 1;
		jpanel.add(addButton, gbc);

		// makes the submit button
		// will make sure that all fields are full before allowing submission
		JButton okButton = new JButton("Submit");
		okButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				boolean temp = true;
				for (Integer i : family) {
					if (i == null) { // if any of the positions is still null, make the user refill out
						JOptionPane.showMessageDialog(new Frame("Input prompt"),
								"Please fill all fields before submitting.");
						temp = false;
						break;
					}
				}

				if (temp) { // if fields are full, get rid of the panel and proceed
					diag.dispose();
				}
			}
		});
		gbc.gridy = 2;
		jpanel.add(okButton, gbc);

		/// Gives user option to cancel inputting info.
		JButton cancelButton = new JButton("Cancel");
		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
//				try {
//					output.close();
//					inputVS.close();
//					pedReader.close();
//				} catch (IOException e1) {
//
//				}
				jframe.dispose();
				System.gc();
				System.exit(0);

			}
		});
		gbc.gridy = 3;
		jpanel.add(cancelButton, gbc);
		
		// makes the panel appear
		diag.getContentPane().add(jpanel);
		diag.pack();
		diag.setLocationRelativeTo(null);
		diag.setVisible(true);
	}
	

	// prompts the user for the headers of the population fields
	// REMEMBER - the user must click the fields after scrolling to the desired
	// header location
	public void getPopulations(String [] scroller) {
		
		/*
		ifP = true;
		ifM = true;
		ifunique = true;
		ifR = true;
		chrM = true;
		*/
		
		popBools = new boolean[5];
		popBools[0] = true; //ifP
		popBools[1] = true; //ifM
		popBools[2] = true; //ifunique
		popBools[3] = true; //ifR
		popBools[4] = true; //chrM
		
		
		// initializes the global variables
		// counts the number of population sets used for filtration
		/*
		numPops = 0;
		numPops2 = 0;
		numPops3 = 0;
		numPops4 = 0;
		numPops5 = 0;
		minorfield = 0;
		*/
		popInts = new int[5];
		popInts[0] = 0; //numPops
		popInts[1] = 0; //numPops2
		popInts[2] = 0; //numPops3
		popInts[3] = 0; //numPops4
		popInts[4] = 0; //numPops5
		
		
		// vartot=new ArrayList<Integer>();
		
		
		
		populations = new ArrayList<int[]>();
		xfilters = new ArrayList<int[]>();
		xfilterc = new ArrayList<int[]>();
		// populations is an array list of interger arrays containing the indexes of the
		// fields being filtered on
		// the sub arrays can be either length 1, containing a filter field, or length
		// 2, containing raw data (allele and chromosome counts)

		// creates the initial frame, panel, dialog, and gbc
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("Input Panel");
		jframe.setSize(2000, 20000);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		final JDialog diag = new JDialog(jframe, "Input Panel", true);
		
		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);
		
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(2, 2, 2, 2);
		gbc.weightx = 1.0;
		gbc.gridwidth = 2;
		// jpanel.add(pedDesc, gbc);

		// just a label

		/*
		 * JLabel Plabel1 = new
		 * JLabel("Would you like to include mendelian inconsistency?");
		 * Plabel1.setAlignmentX(Component.LEFT_ALIGNMENT); gbc.gridx = 0; gbc.gridy =
		 * 0; gbc.gridwidth=1; jpanel.add(Plabel1,gbc);;
		 * 
		 * JCheckBox Pbox1 = new JCheckBox();
		 * 
		 * Pbox1.addActionListener(new ActionListener() {
		 * 
		 * 
		 * public void actionPerformed(ActionEvent e) { popBools[1]=!popBools[1]; } }); gbc.gridx =1;
		 * gbc.gridy=0; jpanel.add(Pbox1, gbc);
		 * 
		 * 
		 * 
		 * JLabel Plabel2 = new
		 * JLabel("Would you like to filter mendelian inconsistency with population frequency?"
		 * ); Plabel2.setAlignmentX(Component.LEFT_ALIGNMENT); gbc.gridx = 0; gbc.gridy
		 * = 1; gbc.gridwidth=1; jpanel.add(Plabel2,gbc);;
		 * 
		 * JCheckBox Pbox2 = new JCheckBox(); Pbox2.addActionListener(new
		 * ActionListener() {
		 * 
		 * 
		 * public void actionPerformed(ActionEvent e) { popBools[0]=!popBools[0]; } }); gbc.gridx =1;
		 * gbc.gridy=1; jpanel.add(Pbox2, gbc);
		 * 
		 * JLabel Plabel3 = new
		 * JLabel("Would you like to print only rows with nonrepeated index?");
		 * Plabel3.setAlignmentX(Component.LEFT_ALIGNMENT); gbc.gridx = 0; gbc.gridy =
		 * 2; gbc.gridwidth=1; jpanel.add(Plabel3,gbc);;
		 * 
		 * JCheckBox Pbox3 = new JCheckBox();
		 * 
		 * Pbox3.addActionListener(new ActionListener() {
		 * 
		 * 
		 * public void actionPerformed(ActionEvent e) { popBools[3]=!popBools[3]; //if checked, only
		 * print nonduplicate } }); gbc.gridx =1; gbc.gridy=2; jpanel.add(Pbox3, gbc);
		 * 
		 * JLabel Plabel300 = new
		 * JLabel("Would you like to return deleterious popBools[4] scored by polyphen Prediction?"
		 * ); Plabel300.setAlignmentX(Component.LEFT_ALIGNMENT); gbc.gridx = 0;
		 * gbc.gridy = 3; gbc.gridwidth=1; jpanel.add(Plabel300,gbc);;
		 * 
		 * JCheckBox Pbox300 = new JCheckBox();
		 * 
		 * Pbox300.addActionListener(new ActionListener() {
		 * 
		 * 
		 * public void actionPerformed(ActionEvent e) { popBools[4]=!popBools[4]; //if checked, only
		 * print nonduplicate } }); gbc.gridx =1; gbc.gridy=3; jpanel.add(Pbox300, gbc);
		 */
		
		JLabel cutoff = new JLabel(
				"Please enter a maximum frequency in decimal format for the minor allele in population data. The default is 0.04. ",
				SwingConstants.LEFT);
		gbc.gridx = 0;
		gbc.gridy = 0;
		gbc.gridwidth = 1;
		jpanel.add(cutoff, gbc);

		rejectionThreshold = 0.04;
		JTextField cuto = new JTextField(10);
		cuto.setEditable(true);
		cuto.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent e) {
			}

			public void focusLost(FocusEvent e) {

				if (cuto.getText().equals("")) {
					cuto.setText("0.04");
				} else {
					try {
						Double.parseDouble(cuto.getText());
					} catch (NumberFormatException e1) {
						JOptionPane.showMessageDialog(new Frame("Population Statistics"),
								"Please enter a valid number");
					}
				}
				rejectionThreshold = Double.parseDouble(cuto.getText());
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 0;
		jpanel.add(cuto, gbc);

		JLabel cutcount00 = new JLabel("Please enter the number of threads. The default is 60. ", SwingConstants.LEFT);
		gbc.gridx = 0;
		gbc.gridy = 1;
		gbc.gridwidth = 1;

		jpanel.add(cutcount00, gbc);

		threads = 60;
		JTextField cutcnt00 = new JTextField(10);
		cutcnt00.setEditable(true);
		cutcnt00.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent e) {
			}

			public void focusLost(FocusEvent e) {

				if (cutcnt00.getText().equals("")) {
					cutcnt00.setText("60");
				} else {
					try {
						Integer.parseInt(cutcnt00.getText());
					} catch (NumberFormatException e1) {
						JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please enter a valid number");
					}
				}
				threads = Integer.parseInt(cutcnt00.getText());
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 1;
		jpanel.add(cutcnt00, gbc);

		/*
		 * familytotal=7777777; JLabel cutcountdd=new
		 * JLabel("Please enter total number of samples. ", SwingConstants.LEFT);
		 * gbc.gridx=0; gbc.gridy=6; gbc.gridwidth=1;
		 * 
		 * jpanel.add(cutcountdd,gbc);
		 * 
		 * JTextField cutcntdd=new JTextField(10); cutcntdd.setEditable(true);
		 * cutcntdd.addFocusListener(new FocusListener() { public void
		 * focusGained(FocusEvent e) { }
		 * 
		 * public void focusLost(FocusEvent e) {
		 * 
		 * if(cutcntdd.getText().equals("")){ cutcntdd.setText("777777"); } else{ try {
		 * Integer.parseInt(cutcntdd.getText()); } catch (NumberFormatException e1) {
		 * JOptionPane.showMessageDialog(new Frame("Input prompt"),
		 * "Please enter a valid number"); } }
		 * familytotal=Integer.parseInt(cutcntdd.getText()); } } ); gbc.gridx = 1;
		 * gbc.gridy = 6; jpanel.add(cutcntdd, gbc);
		 */

		if (probandgender) {
			JLabel cutcount = new JLabel("Please enter the maximum number of hemizygous variant count for Xlink. Default is 5.",
					SwingConstants.LEFT);
			gbc.gridx = 0;
			gbc.gridy = 2;
			gbc.gridwidth = 1;

			jpanel.add(cutcount, gbc);
			
			cutcnt2 = 5;
			JTextField cutcntn = new JTextField(10);
			cutcntn.setEditable(true);
			cutcntn.addFocusListener(new FocusListener() {
				public void focusGained(FocusEvent e) {
				}

				public void focusLost(FocusEvent e) {
					
					if (cutcntn.getText().equals("")) {
						cutcntn.setText(Integer.toString(cutcnt2));
					} else {
						try {
							Integer.parseInt(cutcntn.getText());
						} catch (NumberFormatException e1) {
							JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please enter a valid number");
						}
					}
					cutcnt2 = Integer.parseInt(cutcntn.getText());
				}
			});
			gbc.gridx = 1;
			gbc.gridy = 2;
			jpanel.add(cutcntn, gbc);
		}

		JComboBox msglist13 = new JComboBox(scroller);
		msglist13.setEditable(true);
		
		final int fixed = 3;

		// make the adder button for raw data, adding a new sub array of length 2 to the
		// populations arraylist
		JButton addButton2 = new JButton("Add Raw Data");
		addButton2.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				populations.add(new int[] { 1, 1, 1 });

				final int x = popInts[0];
				final int x2 = popInts[1];
				final int x3 = popInts[2];
				final int x4 = popInts[3];
				final int x5 = popInts[4];

				JLabel AC = new JLabel(
						"Population " + Integer.toString(x + x2 + x3 + x4 + x5 + 1) + " Reference Allele Count:");
				gbc.gridx = 0;
				gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed;
				gbc.gridwidth = 1;
				jpanel.add(AC, gbc);

				int[] temp = { 1, 1, 1 };
				JComboBox msglist = new JComboBox(scroller);
				msglist.setEditable(true);

				msglist.addActionListener(new ActionListener() {

					public void actionPerformed(ActionEvent e) {
						temp[0] = msglist.getSelectedIndex();
					}
				});

				gbc.gridx = 1;
				gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed;
				jpanel.add(msglist, gbc);

				JLabel AC2 = new JLabel(
						"Population " + Integer.toString(x + x2 + x3 + x4 + x5 + 1) + " Alternate Allele Count:");
				gbc.gridx = 0;
				gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed + 1;
				gbc.gridwidth = 1;
				jpanel.add(AC2, gbc);

				JComboBox msglist3 = new JComboBox(scroller);
				msglist3.setEditable(true);

				msglist3.addActionListener(new ActionListener() {

					public void actionPerformed(ActionEvent e) {

						temp[1] = msglist3.getSelectedIndex();
						// vartot.add(temp[1]);
					}
				});

				gbc.gridx = 1;
				gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed + 1;
				jpanel.add(msglist3, gbc);

				JLabel CC = new JLabel(
						"Population " + Integer.toString(x + x2 + x3 + x4 + x5 + 1) + " Genotype Count:");
				gbc.gridx = 0;
				gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed + 2;
				gbc.gridwidth = 1;
				jpanel.add(CC, gbc);

				JComboBox msglist2 = new JComboBox(scroller);
				msglist2.setEditable(true);
				msglist2.addActionListener(new ActionListener() {

					public void actionPerformed(ActionEvent e) {
						temp[2] = msglist2.getSelectedIndex();
						
						popBools[2] = true;
					}
				});

				gbc.gridx = 1;
				gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed + 2;
				jpanel.add(msglist2, gbc);

				popInts[0]++;
				// populations.add(temp);
				populations.set(x + x2, temp);

				diag.pack();
				diag.validate();
				diag.setLocationRelativeTo(null);
				diag.setVisible(true);
			}
		});
		gbc.gridx = 4;
		gbc.gridy = 0;
		jpanel.add(addButton2, gbc);

		// makes the adder button for filter fields, adding a new sub array of length 1
		// to the populations arraylist
		JButton addButton = new JButton("Add Filter Data");
		addButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				populations.add(new int[] { 1, 1 });

				final int x = popInts[0];
				final int x2 = popInts[1];// alt and ref filter field
				final int x3 = popInts[2];
				final int x4 = popInts[3];
				final int x5 = popInts[4];
				JLabel AC = new JLabel(
						"Population " + Integer.toString(x + x2 + x3 + x4 + x5 + 1) + " Filter Field (reference):");
				gbc.gridx = 0;
				gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed;
				gbc.gridwidth = 1;
				jpanel.add(AC, gbc);

				int[] temp = { 1, 1 };
				JComboBox msglist = new JComboBox(scroller);
				msglist.setEditable(true);
				msglist.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						temp[0] = msglist.getSelectedIndex();
					}
				});

				gbc.gridx = 1;
				gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed;
				jpanel.add(msglist, gbc);

				JLabel CC = new JLabel(
						"Population " + Integer.toString(x + x2 + x3 + x4 + x5 + 1) + " Filter Field (alternate):");
				gbc.gridx = 0;
				gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed + 1;
				gbc.gridwidth = 1;
				jpanel.add(CC, gbc);

				JComboBox msglist2 = new JComboBox(scroller);
				msglist2.setEditable(true);
				msglist2.addActionListener(new ActionListener() {

					public void actionPerformed(ActionEvent e) {
						temp[1] = msglist2.getSelectedIndex();
					}
				});

				gbc.gridx = 1;
				gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed + 1;
				jpanel.add(msglist2, gbc);

				popInts[1]++;
				populations.set(x + x2, temp);
				// populations.add(temp);

				diag.pack();
				diag.validate();
				diag.setLocationRelativeTo(null);
				diag.setVisible(true);
			}
		});
		gbc.gridx = 4;
		gbc.gridy = 1;
		jpanel.add(addButton, gbc);

		if (probandgender) {
			JButton addButton40 = new JButton("Add Xlink Raw Data");
			addButton40.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					xfilters.add(new int[] { 1, 1, 1 });

					final int x = popInts[0];
					final int x2 = popInts[1];
					final int x3 = popInts[2];
					final int x4 = popInts[3];// Hemi ref and var and total
					final int x5 = popInts[4];

					JLabel AC5 = new JLabel("Population " + Integer.toString(x + x2 + x3 + x4 + x5 + 1)
							+ " Hemizygous Reference Allele Count:");
					gbc.gridx = 0;
					gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed;
					gbc.gridwidth = 1;
					jpanel.add(AC5, gbc);

					int[] temp = { 1, 1, 1 };
					JComboBox msglist5 = new JComboBox(scroller);
					msglist5.setEditable(true);

					msglist5.addActionListener(new ActionListener() {

						public void actionPerformed(ActionEvent e) {

							temp[0] = msglist5.getSelectedIndex();
						}
					});

					gbc.gridx = 1;
					gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed;
					jpanel.add(msglist5, gbc);

					JLabel AC6 = new JLabel("Population " + Integer.toString(x + x2 + x3 + x4 + x5 + 1)
							+ " Hemizygous Alternate Allele Count:");
					gbc.gridx = 0;
					gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed + 1;
					gbc.gridwidth = 1;
					jpanel.add(AC6, gbc);

					JComboBox msglist7 = new JComboBox(scroller);
					msglist7.setEditable(true);

					msglist7.addActionListener(new ActionListener() {

						public void actionPerformed(ActionEvent e) {
							temp[1] = msglist7.getSelectedIndex();
						}
					});

					gbc.gridx = 1;
					gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed + 1;
					jpanel.add(msglist7, gbc);

					JLabel CC8 = new JLabel("Population " + Integer.toString(x + x2 + x3 + x4 + x5 + 1)
							+ " Hemizygous Genotype Count:");
					gbc.gridx = 0;
					gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed + 2;
					gbc.gridwidth = 1;
					jpanel.add(CC8, gbc);

					JComboBox msglist9 = new JComboBox(scroller);
					msglist9.setEditable(true);
					msglist9.addActionListener(new ActionListener() {

						public void actionPerformed(ActionEvent e) {
							temp[2] = msglist9.getSelectedIndex();
						}
					});

					gbc.gridx = 1;
					gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed + 2;
					jpanel.add(msglist9, gbc);

					popInts[3]++;
					xfilters.set(x4 + x3, temp);
					// xfilters.add(temp);

					diag.pack();
					diag.validate();
					diag.setLocationRelativeTo(null);
					diag.setVisible(true);
				}
			});
			gbc.gridx = 4;
			gbc.gridy = 3;
			jpanel.add(addButton40, gbc);

			JButton addButton30 = new JButton("Add Xlink Filter Data");
			addButton30.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					xfilters.add(new int[] { 1, 1 });

					final int x = popInts[0];
					final int x2 = popInts[1];
					final int x3 = popInts[2];// hemizygoug MAF
					final int x4 = popInts[3];
					final int x5 = popInts[4];
					JLabel AC3 = new JLabel("Population " + Integer.toString(x + x2 + x3 + x4 + x5 + 1)
							+ " Filter Field (Hemizygous Reference MAF):");
					gbc.gridx = 0;
					gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed;
					gbc.gridwidth = 1;
					jpanel.add(AC3, gbc);

					int[] temp3 = { 1, 1 };
					// JList<String> tempJL1 = new JList<String>(curLine);
					JComboBox msglist3 = new JComboBox(scroller);
					msglist3.setEditable(true);

					msglist3.addActionListener(new ActionListener() {

						public void actionPerformed(ActionEvent e) {
							temp3[0] = msglist3.getSelectedIndex();
						}
					});

					gbc.gridx = 1;
					gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed;
					jpanel.add(msglist3, gbc);

					JLabel CC3 = new JLabel("Population " + Integer.toString(x + x2 + x3 + x4 + x5 + 1)
							+ " Filter Field (Hemizygous Alternate MAF):");
					gbc.gridx = 0;
					gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed + 1;
					gbc.gridwidth = 1;
					jpanel.add(CC3, gbc);

					JComboBox msglist4 = new JComboBox(scroller);
					msglist4.setEditable(true);
					msglist4.addActionListener(new ActionListener() {

						public void actionPerformed(ActionEvent e) {
							temp3[1] = msglist4.getSelectedIndex();
						}
					});

					gbc.gridx = 1;
					gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed + 1;
					jpanel.add(msglist4, gbc);

					popInts[2]++;
					xfilters.set(x3 + x4, temp3);

					diag.pack();
					diag.validate();
					diag.setLocationRelativeTo(null);
					diag.setVisible(true);
				}
			});
			gbc.gridx = 4;
			gbc.gridy = 2;
			jpanel.add(addButton30, gbc);

			JButton addButton10 = new JButton("Add Xlink Count Data");
			addButton10.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					xfilterc.add(new int[] { 1, 1 });

					final int x = popInts[0];
					final int x2 = popInts[1];
					final int x3 = popInts[2];
					final int x4 = popInts[3];
					final int x5 = popInts[4];// hemizygous count
					JLabel AC3 = new JLabel("Population " + Integer.toString(x + x2 + x3 + x4 + x5 + 1)
							+ " Hemizygous Reference Count:");
					gbc.gridx = 0;
					gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed;
					gbc.gridwidth = 1;
					jpanel.add(AC3, gbc);

					int[] temp3 = { 1, 1 };
					JComboBox msglist3 = new JComboBox(scroller);
					msglist3.setEditable(true);

					msglist3.addActionListener(new ActionListener() {

						public void actionPerformed(ActionEvent e) {
							temp3[0] = msglist3.getSelectedIndex();
						}
					});

					gbc.gridx = 1;
					gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed;
					jpanel.add(msglist3, gbc);

					JLabel CC3 = new JLabel("Population " + Integer.toString(x + x2 + x3 + x4 + x5 + 1)
							+ " Hemizygous Alternate Count:");
					gbc.gridx = 0;
					gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed + 1;
					gbc.gridwidth = 1;
					jpanel.add(CC3, gbc);

					JComboBox msglist4 = new JComboBox(scroller);
					msglist4.setEditable(true);
					msglist4.addActionListener(new ActionListener() {

						public void actionPerformed(ActionEvent e) {
							temp3[1] = msglist4.getSelectedIndex();
						}
					});

					gbc.gridx = 1;
					gbc.gridy = x * 3 + x2 * 2 + x3 * 2 + x4 * 3 + x5 * 2 + fixed + 1;
					jpanel.add(msglist4, gbc);

					popInts[4]++;
					xfilterc.set(x5, temp3);

					diag.pack();
					diag.validate();
					diag.setLocationRelativeTo(null);
					diag.setVisible(true);
				}
			});
			gbc.gridx = 4;
			gbc.gridy = 4;
			jpanel.add(addButton10, gbc);

		}
		// makes the submit button
		JButton okButton = new JButton("Submit");
		okButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				boolean temp = true;
				if (populations.size() != 0) {
					for (int[] i : populations) {
						for (int n = 0; n < i.length; n++) {
							if (i[n] == 1) { // if any of the positions is still null, make the user refill out
								JOptionPane.showMessageDialog(new Frame("Input prompt"),
										"Please fill all population fields before submitting.");
								temp = false;
								break;
							}
						}
					}
				}
				if (xfilters.size() != 0) {
					for (int[] i : xfilters) {
						for (int n = 0; n < i.length; n++) {
							if (i[n] == 1) { // if any of the positions is still null, make the user refill out
								JOptionPane.showMessageDialog(new Frame("Input prompt"),
										"Please fill all xlink fields before submitting.");
								temp = false;
								break;
							}
						}
					}
				}
				if (xfilterc.size() != 0) {
					for (int[] i : xfilterc) {
						for (int n = 0; n < i.length; n++) {
							if (i[n] == 1) { // if any of the positions is still null, make the user refill out
								JOptionPane.showMessageDialog(new Frame("Input prompt"),
										"Please fill all xlink count fields before submitting.");
								temp = false;
								break;
							}
						}
					}
				}
				if (temp) { // if fields are full, get rid of the panel and proceed
					diag.dispose();
				}
			}
		});
		gbc.gridx = 5;
		gbc.gridy = 0;
		jpanel.add(okButton, gbc);
		

		/// Gives user option to cancel inputting info.
		JButton cancelButton = new JButton("Cancel");
		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
//				try {
//					output.close();
//					inputVS.close();
//					pedReader.close();
//				} catch (IOException e1) {
//
//				}
				jframe.dispose();
				System.gc();
				System.exit(0);

			}
		});
		gbc.gridy = 1;
		jpanel.add(cancelButton, gbc);
		
		
		// makes the panel appear
		diag.getContentPane().add(jpanel);
		diag.pack();
		diag.setLocationRelativeTo(null);
		diag.setVisible(true);
	}
	
	//Asks if user would like to filter the Xlinked diseases based on external hemizygous columns
	public int askExtDatabase() {
		int j = JOptionPane.showConfirmDialog(new Frame("Input prompt"),
				"Would you like to use External Databases to filter through the Hemizygous Columns\n"
				+ "of databases such as ExAC or GnomAD",
				"User Prompt", JOptionPane.YES_NO_OPTION);
		/// 0 = yes. 1 = no.
		return j; // User wants to use default
	}
	
	
	//Gets the External Hemizygous filtering columns
	public void getExternalHemiHeader(String[] scroller) {
		
		
		/// Create an initial frame, panel, dialog, and gridBagConstraints
		JFrame jframe = new JFrame();				
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("Input Panel");
		jframe.setSize(2000, 20000);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		final JDialog diag = new JDialog(jframe, "Input Panel", true);
		
		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);
		
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(2, 2, 2, 2);
		gbc.weightx = 1.0;
		gbc.gridwidth = 2;
		
		/// Request user input for field of hemizygous columns
		JLabel extDesc = new JLabel("Please enter the hemizygous count columns from external databases such as ExAC and GnomAD");
		gbc.gridx = 0;
		gbc.gridy = 0;
		gbc.gridwidth = 1;
		jpanel.add(extDesc, gbc);
		
		
		xfilterExt.add(new Integer(0));
		popExtInt = 0;
		
		final int x = popExtInt;
		
		/// Request user input for header of hemizygous columns
		JLabel extAsk = new JLabel("Hemizygous External Database #" + Integer.toString(x + 1));
		gbc.gridx = 0;
		gbc.gridy = x + 1;
		gbc.gridwidth = 1;
		jpanel.add(extAsk, gbc);
		
		
		JComboBox msglist = new JComboBox(scroller);
		msglist.setEditable(true);
		
		msglist.addActionListener(new ActionListener() {
			
			public void actionPerformed(ActionEvent e) {
				xfilterExt.set(x, msglist.getSelectedIndex());
			}
		});
		gbc.gridx = 1;
		gbc.gridy = (x+1);
		jpanel.add(msglist, gbc);
		
		// make the adder button - this will let the user add databases
		JButton addButton = new JButton("Add Database");
		addButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				
				
				popExtInt++;
				
				final int x = popExtInt;
				
				xfilterExt.add(new Integer(0));
				
				
				/// Request user input for header of hemizygous columns
				JLabel extAsk = new JLabel("Hemizygous External Database #" + Integer.toString(x + 1));
				gbc.gridx = 0;
				gbc.gridy = x + 1;
				gbc.gridwidth = 1;
				jpanel.add(extAsk, gbc);
				
				
				JComboBox msglist = new JComboBox(scroller);
				msglist.setEditable(true);

				msglist.addActionListener(new ActionListener() {
					
					public void actionPerformed(ActionEvent e) {
						xfilterExt.set(x, msglist.getSelectedIndex());
					}
				});
				gbc.gridx = 1;
				gbc.gridy = (x+1);
				jpanel.add(msglist, gbc);
				
				diag.pack();
				diag.validate();
			}
		});
		gbc.gridx = 2;
		gbc.gridy = 0;
		jpanel.add(addButton, gbc);
		
		
		
		/// Makes the submit button and ensures all fields are full.
		JButton okButton = new JButton("Submit");
		okButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				boolean temp = true;
				for (int i : xfilterExt) {
					if (i == 0) { // if any of the positions is still null, make the user refill out
						JOptionPane.showMessageDialog(new Frame("Input prompt"),
								"Please fill all hemizygous external database fields before submitting.");
						temp = false;
						break;
					}
				}
				if (temp) { /// if fields are full, get rid of the panel and proceed
					diag.dispose();
				}
			}
		});

		gbc.gridx = 2;
		gbc.gridy = 1;
		jpanel.add(okButton, gbc);
		
		/// Gives user option to cancel inputting info.
		JButton cancelButton = new JButton("Cancel");
		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				jframe.dispose();
				System.gc();
				System.exit(0);

			}
		});

		gbc.gridx = 2;
		gbc.gridy = 2;
		jpanel.add(cancelButton, gbc);

		/// Makes the panel appear
		diag.getContentPane().add(jpanel);
		diag.pack();
		diag.validate();
		diag.setLocationRelativeTo(null);
		diag.setVisible(true);
	}
	
	
	
	/*
	 * User input of exonic CADD threshold
	 */
	public Double[] delThresholdCalc() {
		
		String thresh = JOptionPane.showInputDialog(new Frame("Choose threshold"),
				"Please write an exonic CADD threshold. 19 or 20 are recommended.", JOptionPane.PLAIN_MESSAGE);
		
		//Determines if they want to cancel or inputed an empty value accidentally
		if (thresh == null || thresh.equals("")) {
			int i = TXTFile.wantToContinue("You did not enter a valid CADD Threshold.");
			if (i == 1) {
				System.exit(0);
			}
			return delThresholdCalc();
		}
		
		Double exonCut;
		//Determines if input is a double value
		try {
			exonCut = Double.parseDouble(thresh);
		} catch (NumberFormatException e1) {
			int i = TXTFile.wantToContinue("You did not enter a valid CADD Threshold.");
			if (i == 1) {
				System.exit(0);
			}
			return delThresholdCalc();
		}
		
		//Calculates CADD threshold for intronic variants
		double intronCut = 0.75 * exonCut;
		Double[] cutArr = { intronCut, exonCut };

		/// Ask user if they're okay with this
		int confirmed = TXTFile.wantToContinue("Your intronic CADD threshold is " + fix(Double.toString(intronCut))
				+ " and exonic CADD threshold is " + exonCut + ". Is this okay?");

		if (confirmed == JOptionPane.YES_OPTION) {
			return cutArr;
		} else { // If user not satisfied, then ask user to input again.
			delThresholdCalc();
		}
		return cutArr;
	}
	

	//Truncator Method that truncates the intronic CADD number to 3 decimal places
	public String fix(String tmp) {
		double tmp2 = Double.parseDouble(tmp);
		String form = "#0.";
		for (int x = 0; x < 3; x++) {
			form = form + "#";
		}
		
		NumberFormat formatter1 = new DecimalFormat(form);
		String out = formatter1.format(tmp2);
		
		if (Double.parseDouble(out) == 0) {
			out = "0";
		}
			
		return out;
	}
	
	/*
	 * User input for minimum cutoff counts in datasets for each genotype (homvar, hemvar,
	 * het...).
	 */
	public void getPopulationGenotypes() {
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("Input Panel");
		jframe.setSize(2000, 20000);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		final JDialog diag = new JDialog(jframe, "Input Panel", true);

		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);

		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(2, 2, 2, 2);
		gbc.weightx = 1.0;
		gbc.gridwidth = 3;

		/// Request user input for thresholds
		/// HomVar
		JLabel cutoff00 = new JLabel(
				"Please enter a maximum number of people in a dataset who can have a homozygous variant genotype. Default is 3.",
				SwingConstants.LEFT);
		gbc.gridx = 0;
		gbc.gridy = 0;
		gbc.gridwidth = 1;
		jpanel.add(cutoff00, gbc);
		
		homVarThreshold = 3;
		JTextField homVarInput = new JTextField(10);
		homVarInput.setEditable(true);
		homVarInput.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent e) {
			}
			
			public void focusLost(FocusEvent e) {

				if (homVarInput.getText().equals("")) {
					homVarInput.setText("3");
				} else {
					try {
						Integer.parseInt(homVarInput.getText());
					} catch (NumberFormatException e1) {
						JOptionPane.showMessageDialog(new Frame("Population Statistics"),
								"Please enter a valid number");
					}
				}
				homVarThreshold = Integer.parseInt(homVarInput.getText());
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 0;
		jpanel.add(homVarInput, gbc);

		/// HomVar
		JLabel cutoff01 = new JLabel(
				"Please enter a maximum number of people in a dataset who can have a heterozygous genotype. Default is 3.",
				SwingConstants.LEFT);
		gbc.gridx = 0;
		gbc.gridy = 1;
		gbc.gridwidth = 1;
		jpanel.add(cutoff01, gbc);

		hetThreshold = 3;
		JTextField hetInput = new JTextField(10);
		hetInput.setEditable(true);
		hetInput.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent e) {
			}

			public void focusLost(FocusEvent e) {

				if (hetInput.getText().equals("")) {
					hetInput.setText("3");
				} else {
					try {
						Integer.parseInt(hetInput.getText());
					} catch (NumberFormatException e1) {
						JOptionPane.showMessageDialog(new Frame("Population Statistics"),
								"Please enter a valid number");
					}
				}
				hetThreshold = Integer.parseInt(hetInput.getText());
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 1;
		jpanel.add(hetInput, gbc);

		/// HemVar
		/// HomVar
		JLabel cutoff02 = new JLabel(
				"Please enter a maximum number of people in a dataset who can have a hemizygous variant genotype. Default is 3.",
				SwingConstants.LEFT);
		gbc.gridx = 0;
		gbc.gridy = 2;
		gbc.gridwidth = 1;
		jpanel.add(cutoff02, gbc);

		hemVarThreshold = 3;
		JTextField hemVarInput = new JTextField(10);
		hemVarInput.setEditable(true);
		hemVarInput.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent e) {
			}

			public void focusLost(FocusEvent e) {

				if (hemVarInput.getText().equals("")) {
					hemVarInput.setText("3");
				} else {
					try {
						Integer.parseInt(hemVarInput.getText());
					} catch (NumberFormatException e1) {
						JOptionPane.showMessageDialog(new Frame("Population Statistics"),
								"Please enter a valid number");
					}
				}
				hemVarThreshold = Integer.parseInt(hemVarInput.getText());
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 2;
		jpanel.add(hemVarInput, gbc);

		/// Submit button
		JButton okButton = new JButton("Submit");
		okButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				boolean temp = true;
				if (homVarThreshold <= 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"),
							"Please input homozygous variant count minimum.");
					temp = false;
				}
				if (hetThreshold <= 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"),
							"Please input heterozygous count minimum.");
					temp = false;
				}
				if (hemVarThreshold <= 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"),
							"Please input hemizygous variant count minimum.");
					temp = false;
				}

				if (temp) { // if fields are full, get rid of the panel and proceed
					diag.dispose();
				}
			}
		});
		gbc.gridx = 4;
		gbc.gridy = 0;
		jpanel.add(okButton, gbc);

		/// Cancel button
		JButton cancelButton = new JButton("Cancel");
		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				jframe.dispose();
				System.gc();
				System.exit(0);
			}
		});
		gbc.gridx = 4;
		gbc.gridy = 1;
		jpanel.add(cancelButton, gbc);

		/// Makes the panel appear
		diag.getContentPane().add(jpanel);
		diag.pack();
		diag.setLocationRelativeTo(null);
		diag.setVisible(true);
	}
	
	//Prechecks Config file for genotype headers or unknown headers
	//Returns 1 if it is a genotype header, 2 if it is an unknown header, and 0 if anything else
	public int checkROCConfigHeaders(String categ, String vsCol) {
		
		if (categ.equals("Tot_AC") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("Tot_AN") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("Tot_AF") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("HomRef_Count") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("HomVar_Count") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("Het_Count") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("Other_Count") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("Genotype_Count") && (headzList.indexOf(vsCol) != -1)) {
			return 1;
		} else if (categ.equals("HemRef_Count") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("HemVar_Count") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("MaxMAF1") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("MaxMAF2") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("AFR_AF") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("EAS_AF") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("AMR_AF") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("EUR_AF") && headzList.indexOf(vsCol) != -1) {
			return 0;
		} else if (categ.equals("SAS_AF") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("OTH_AF") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("AFR_AC") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("EAS_AC") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("AMR_AC") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("EUR_AC") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("SAS_AC") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("OTH_AC") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("AFR_AN") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("EAS_AN") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("AMR_AN") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("EUR_AN") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("SAS_AN") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("OTH_AN") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else if (categ.equals("Reference_is_minor") && (headzList.indexOf(vsCol) != -1)) {
			return 0;
		} else {
			
			return 2;
			
		}
		
	}
	

	/*
	 * User input of minimum genotype count cutoff of all genotypes that are found in the ROC config
	 */
	public int genoCountThreshold(int headerIndex, List<String> headers) {
		String genoHeader = headers.get(headerIndex);
		
		String genoThreshStr = JOptionPane.showInputDialog(new Frame("Choose threshold"),
				"Please write an integer-value minimum genotype count for the dataset " + genoHeader + ". \n",
				JOptionPane.PLAIN_MESSAGE);
		try {
			int genoThresh = Integer.parseInt(genoThreshStr);
			return genoThresh;
		} catch (NumberFormatException e) {
			int confirmed = TXTFile.wantToContinue("The value you entered is not an integer. Please try again.");
			if (confirmed == JOptionPane.NO_OPTION) {
				System.gc();
				System.exit(0);

			}
			return genoCountThreshold(headerIndex, headers);
		}
	}
	

	
	/*
	 * Ask user for input in headers. This is derived from the KaylaKode GUI code.
	 * The output file, input VS file, and pedigree files are included as inputs for
	 * the cancellation option, so that the PrintWriter and BufferedReader streams
	 * are successfully closed if the user chooses to abort the job.
	 */
	public void getHeaders(String[] scroller) {
		// Create an initial frame, panel, dialog, and gridBagConstraints
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("Input Panel");
		jframe.setSize(2000, 20000);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		final JDialog diag = new JDialog(jframe, "Input Panel", true);

		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);

		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(2, 2, 2, 2);
		gbc.weightx = 1.0;
		gbc.gridwidth = 2;
		
		/// Requests user input for field denoting population homref count
		JLabel homRefAsk = new JLabel("Enter field for population homozygous reference count: ");
		gbc.gridx = 0;
		gbc.gridy = 0;
		gbc.gridwidth = 1;
		jpanel.add(homRefAsk, gbc);

		JComboBox msglist11 = new JComboBox(scroller);
		msglist11.setEditable(true);

		msglist11.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				homRefIndex = msglist11.getSelectedIndex();
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 0;
		jpanel.add(msglist11, gbc);

		/// Requests user input for field denoting population homvar count
		JLabel homVarAsk = new JLabel("Enter field for population homozygous variant count: ");
		gbc.gridx = 0;
		gbc.gridy = 1;
		gbc.gridwidth = 1;
		jpanel.add(homVarAsk, gbc);

		JComboBox msglist12 = new JComboBox(scroller);
		msglist12.setEditable(true);

		msglist12.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				homVarIndex = msglist12.getSelectedIndex();
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 1;
		jpanel.add(msglist12, gbc);

		/// Requests user input for field denoting population genotype count
		JLabel genotypeAsk = new JLabel("Enter field for population genotype count: ");
		gbc.gridx = 0;
		gbc.gridy = 2;
		gbc.gridwidth = 1;
		jpanel.add(genotypeAsk, gbc);

		JComboBox msglist13 = new JComboBox(scroller);
		msglist13.setEditable(true);

		msglist13.addActionListener(new ActionListener() {
			
			public void actionPerformed(ActionEvent e) {
				genotypeIndex = msglist13.getSelectedIndex();
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 2;
		jpanel.add(msglist13, gbc);
		
		// Requests user input for field denoting population genotype threshold
		JLabel genoThresh = new JLabel(
				"Please input an integer-value minimum genotype threshold: ",
				SwingConstants.LEFT);
		gbc.gridx = 0;
		gbc.gridy = 3;
		gbc.gridwidth = 1;
		jpanel.add(genoThresh, gbc);

		genotypeThresh = 0;
		JTextField geno = new JTextField(10);
		geno.setEditable(true);
		geno.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent e) {
			}

			public void focusLost(FocusEvent e) {

				if (geno.getText().equals("")) {
					geno.setText(Integer.toString(genotypeThresh));
				} else {
					try {
						Integer.parseInt(geno.getText());
					} catch (NumberFormatException e1) {
						JOptionPane.showMessageDialog(new Frame("Genotype Threshold"),
								"Please enter a valid number");
					}
				}
				genotypeThresh = Integer.parseInt(geno.getText());
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 3;
		jpanel.add(geno, gbc);
		
		// Requests user input for integer that determines if a base has a low Base Quality Score
		JLabel baseFilterText = new JLabel(
				"Please input an integer-value minimum base quality score: ",
				SwingConstants.LEFT);
		gbc.gridx = 0;
		gbc.gridy = 4;
		gbc.gridwidth = 1;
		jpanel.add(baseFilterText, gbc);

		baseFilter = 20;
		JTextField baseFilterAsk = new JTextField(10);
		baseFilterAsk.setEditable(true);
		baseFilterAsk.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent e) {
			}

			public void focusLost(FocusEvent e) {

				if (baseFilterAsk.getText().equals("")) {
					baseFilterAsk.setText(Integer.toString(baseFilter));
				} else {
					try {
						Integer.parseInt(baseFilterAsk.getText());
					} catch (NumberFormatException e1) {
						JOptionPane.showMessageDialog(new Frame("Base Filter Threshold"),
								"Please enter a valid number");
					}
				}
				baseFilter = Integer.parseInt(baseFilterAsk.getText());
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 4;
		jpanel.add(baseFilterAsk, gbc);
		
		// Requests user input for integer that determines the number of bad base qualities to make a read considered a bad read
		JLabel poorBaseThresh = new JLabel(
				"Please input an integer-value minimum for number of bad base qualities to be considered a bad read: ",
				SwingConstants.LEFT);
		gbc.gridx = 0;
		gbc.gridy = 5;
		gbc.gridwidth = 1;
		jpanel.add(poorBaseThresh, gbc);
		
		baseQualThresh = 10;
		JTextField poorBaseThreshAsk = new JTextField(10);
		poorBaseThreshAsk.setEditable(true);
		poorBaseThreshAsk.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent e) {
			}

			public void focusLost(FocusEvent e) {

				if (poorBaseThreshAsk.getText().equals("")) {
					poorBaseThreshAsk.setText(Integer.toString(baseQualThresh));
				} else {
					try {
						Integer.parseInt(poorBaseThreshAsk.getText());
					} catch (NumberFormatException e1) {
						JOptionPane.showMessageDialog(new Frame("Bad Base Reads Threshold"),
								"Please enter a valid number");
					}
				}
				baseQualThresh = Integer.parseInt(poorBaseThreshAsk.getText());
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 5;
		jpanel.add(poorBaseThreshAsk, gbc);
		
		// Requests user input for double that determines the ratio of bad reads/ total reads in a region
		JLabel badReadRatioText = new JLabel(
				"Please input an decimal for the percent of bad reads accepted per region: ",
				SwingConstants.LEFT);
		gbc.gridx = 0;
		gbc.gridy = 6;
		gbc.gridwidth = 1;
		jpanel.add(badReadRatioText, gbc);
		
		badReadRatio = 0.25;
		JTextField badReadRatioAsk = new JTextField(10);
		badReadRatioAsk.setEditable(true);
		badReadRatioAsk.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent e) {
			}

			public void focusLost(FocusEvent e) {

				if (badReadRatioAsk.getText().equals("")) {
					badReadRatioAsk.setText(Double.toString(badReadRatio));
				} else {
					try {
						Double.parseDouble(badReadRatioAsk.getText());
					} catch (NumberFormatException e1) {
						JOptionPane.showMessageDialog(new Frame("Bad Read Ratio"),
								"Please enter a valid number");
					}
				}
				badReadRatio = Double.parseDouble(badReadRatioAsk.getText());
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 6;
		jpanel.add(badReadRatioAsk, gbc);
		
		
		/// Makes the submit button and ensures all fields are full.
		JButton okButton = new JButton("Submit");
		okButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				boolean temp = true;
				if (homVarIndex == 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please input homozygous variant field.");
					temp = false;
				}

				if (homRefIndex == 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"),
							"Please input homozygous reference field.");
					temp = false;
				}

				if (genotypeIndex == 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please input genotype count field.");
					temp = false;
				}
				
				if (genotypeThresh <= 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please input a valid genotype threshold.");
					temp = false;
				}
				
				if (baseFilter <= 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please input a valid base quality filter.");
					temp = false;
				}
				
				if (baseQualThresh <= 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please input a valid number for base quality bad read threshold.");
					temp = false;
				}
				
				if (badReadRatio <= 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please input a valid ratio for the % of bad reads accepted by the confetti filter.");
					temp = false;
				}
				
				if (temp) { /// if fields are full, get rid of the panel and proceed
					diag.dispose();
				}
			}
		});

		gbc.gridx = 4;
		gbc.gridy = 0;
		jpanel.add(okButton, gbc);

		/// Gives user option to cancel inputting info.
		JButton cancelButton = new JButton("Cancel");
		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				
				jframe.dispose();
				System.gc();
				System.exit(0);

			}
		});
		gbc.gridx = 4;
		gbc.gridy = 1;
		jpanel.add(cancelButton, gbc);

		// makes the panel appear
		diag.getContentPane().add(jpanel);
		diag.pack();
		diag.setLocationRelativeTo(null);
		diag.setVisible(true);
	}
	


	/*
	 * Requests user input for CNC cluster minimum size.
	 */
	public int clusterMin() {
		String clusterQuant = JOptionPane.showInputDialog(new Frame("Choose cluster size minimum"),
				"Please write how many CNC variants you want to be considered a CNC cluster size minimum. \n"
						+ " Please keep your input as a positive integer. If you choose 1, then the bounds \n "
						+ "of the exon will be +/- 75 bp of the CNC variant's position in the genome.",
				JOptionPane.PLAIN_MESSAGE);
		if (isInteger(clusterQuant)) {
			int clusterInt = Integer.parseInt(clusterQuant);
			if (clusterInt > 0) {
				return Integer.parseInt(clusterQuant);
			} else {
				int i = TXTFile.wantToContinue("Please input a positive integer.");
				if (i == 1) {
					System.exit(0);
				}
				return clusterMin();
			}
		} else {
			int i = TXTFile.wantToContinue("Please input a positive integer.");
			if (i == 1) {
				System.exit(0);
			}
			return clusterMin();
		}
	}
	
	//Asks if user wants to filter out all non-exonic CNC variants
	public int askCNCExonOnly() {
		int j = JOptionPane.showConfirmDialog(new Frame("Input prompt"),
				"Would you like to filter CNCs by removing all CNCs that do not lie within an exon boundary?",
				"User Prompt", JOptionPane.YES_NO_OPTION);
		/// 0 = yes. 1 = no.
		return j; // User wants to use default
	}
	
	
	/*
	 * Returns true if input is an integer.Derived from
	 * https://www.quora.com/How-do-you-check-if-a-value-in-Java-is-of-an-integer-
	 * type
	 */
	public boolean isInteger(Object object) {
		if (object instanceof Integer) {
			return true;
		} else {
			String str = object.toString();
			try {
				Integer.parseInt(str);
			} catch (Exception e) {
				return false;
			}
		}
		return true;
	}
	
	
	////////////////////////////////////////////////////////////////////////////GETTER METHODS///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	public File getVSFile() {
		
		return vsFile;
		
	}
	
	public String getDest() {
		
		return destination;
		
	}
	
	public File getPedFile() {
		
		return pedFile;
		
	}
	
	public File getBamDir() {
		
		return BamDirectoryFile;
		
	}
 	
	public String getConfigPath() {
		
		return configPath;
		
	}
	
	///////////////////////////////////////////////////////////////Ethnicity Matcher///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	public File getSampleFile() {
		
		return sampleFile;
		
	}
	
	public String getCriteria() {
		
		return cutoffCriteria;
		
	}
	
	public String[] getMarkers() {
		
		return markerPaths;
		
	}
	
	///////////////////////////////////////////////Salvage Pathway/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	public int getMinor() {
		
		return minor;
		
	}
	
	public int getYesPop() {
		
		return yesSetPop;
		
	}
	
	public int[] getPopIndexes() {
		
		return populationHeaders;
		
	}
	
	///////////////////////////////////////////////Kayla Kode/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	public ArrayList<Integer> getFam() {
		
		return family;
		
	}
	
	public ArrayList<Boolean> getSibsAffected() {
		
		return sibsAffected;
		
	}
	
	public ArrayList<Boolean> getSibsGender() {
		
		return sibsGender;
		
	}
	
	public boolean getProbandGend() {
		
		return probandgender;
	}
	
	public boolean[] getPopBools() {
		
		return popBools;
		
	}
	
	public ArrayList<int[]> getPopulations() {
		
		return populations;
		
	}
	
	public ArrayList<int[]> getXFilterS() {
		
		return xfilters;
	}
	
	public ArrayList<int[]> getXFilterC() {
		
		return xfilterc;
	}
	

	public ArrayList<Integer> getXFilterExt() {
		
		return xfilterExt;
		
	}
	
	public int getExtData() {
		
		return extData;
		
	}
	
	public double getRejectionThreshold() {
		
		return rejectionThreshold;
		
	}
	
	public int getThreads() {
		
		return threads;
		
	}
	
	public int getCutCnt2() {
		
		return cutcnt2;
	}
	
	public String getStrengthPath() {
		
		return strengthPath;
	}
	
	public String getGenePath() {
		
		return genePath;
		
	}
	
	
	///////////////////////////////////////////////Variant Exclusion Filter/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public Double[] getCaddThresh() {
		
		return caddThresh;
		
	}
	
	public int getHomVarThreshold() {
		
		return homVarThreshold;
		
	}
	

	public int getHetThreshold() {
		
		return hetThreshold;
		
	}

	public int getHemVarThreshold() {
		
		return hemVarThreshold;
		
	}
	
	public File getROCHeaderConfig() {
		
		return ROCHeaderConfig;
	}
	
	public ArrayList<Integer> getGenoThreshList() {
		
		return genotypeThreshList;
		
	}
	
	
	///////////////////////////////////////////////Confetti Filter/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public int getHomRefIndex() {
		
		return homRefIndex;
		
	}
	
	public int getHomVarIndex() {
		
		return homVarIndex;
		
	}
	
	public int getGenotypeIndex() {
		
		return genotypeIndex;
		
	}
	
	public int getGenotypeThreshold() {
		
		return genotypeThresh;
		
	}
	
	public int getBaseFilter() {
		
		return baseFilter;
		
	}
	
	public int getBaseQualThresh() {
		
		return baseQualThresh;
		
	}
	

	public double getBadReadRatio() {
		
		return badReadRatio;
		
	}
	
	
	///////////////////////////////////////////////CNC Filter/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public int getClusterMin() {
		
		return clusterSize;
		
	}
	
	public String getExonBoundConfig() {
		
		return exonBoundConfig;
		
	}
	
	public int getCNCExon() {
		
		return CNCExon;
		
	}
	
	//Prints out file with all of user inputs
	public void generateUserFile(String dest) throws IOException {
		
		File userOutput = new File(dest);
		
		PrintWriter userWriter = new PrintWriter(userOutput);
		
		userWriter.println("General UI");
		userWriter.println("VSFile Path: " + vsFile.getAbsolutePath());
		
		userWriter.println("Output Name: " + destination);
		
		userWriter.println("Ped File Path: " + pedFile.getAbsolutePath());
		
		userWriter.println("Bam Directory File: " + BamDirectoryFile.getAbsolutePath());
		
		userWriter.println();
		userWriter.println("Ethnicity Matcher UI");
		userWriter.println("Sample File Path: " + sampleFile.getAbsolutePath());
		userWriter.println("Cutoff Criteria: " + cutoffCriteria);
		for(int i = 0; i < markerPaths.length; i++) {
			userWriter.println("Marker Paths #" + (i+1) + ": " +  markerPaths[i]);
		}
		
		userWriter.println();
		userWriter.println("Salvage Pathway UI");
		userWriter.println("Minor Allele Index + Header: " + minor + ": " + scroller[minor]);
		userWriter.println("Using 1K Genomes for Salvage Population Data (0 for yes, 1 for no): " + yesSetPop);
		if(yesSetPop == 1) {
			for(int i = 0; i < populationHeaders.length; i++) {
				
				userWriter.println("Population Headers #" + (i+1) + ": " + populationHeaders[i]);
				
			}
		}
		
		userWriter.println();
		userWriter.println("Kayla Kode UI");
		userWriter.println("Strength Config Path: " + strengthPath);
		userWriter.println("Gene Config Path: " + genePath);
		userWriter.println("Number of siblings: " + numSibs);
		for(int i = 0; i < family.size(); i++) {
			userWriter.println("Family Member #" + (i+1) + "'s index: " + family.get(i));
		}
		for(int i = 0; i < familyNames.size(); i++) {
			userWriter.println("Family Member #" + (i+1) + "'s name: " + familyNames.get(i));
		}
		for(int i = 0; i < sibsAffected.size(); i++) {
			userWriter.println("Affected Sib #" + (i+1) + ": " + sibsAffected.get(i));
		}
		for(int i = 0; i < sibsGender.size(); i++) {
			userWriter.println("Gender for Sib #" + (i+1) + ": " + sibsGender.get(i));
		}
		userWriter.println("Proband Gender: " + probandgender);
		userWriter.println("ifP: " + popBools[0]);
		userWriter.println("ifM: " + popBools[1]);
		userWriter.println("ifunique: " + popBools[2]);
		userWriter.println("ifR: " + popBools[3]);
		userWriter.println("chrM: " + popBools[4]);
		userWriter.println("Number of Raw Data Inputs: " + popInts[0]);
		userWriter.println("Number of Filter Data Inputs: " + popInts[1]);
		userWriter.println("Number of Hemizygous Filter Data Inputs: " + popInts[2]);
		userWriter.println("Number of Hemizygous Raw Data Inputs: " + popInts[3]);
		userWriter.println("Number of Hemizygous Count Data Inputs: " + popInts[4]);
		userWriter.println("Minor Field Index: " + minor);
		for(int i = 0; i < populations.size(); i++) {
			for(int j =0; j < populations.get(i).length; j++) {
				userWriter.println("Populations #" + (i+1) + " and header index #" + (j+1) + ": " + populations.get(i)[j] + ": " + scroller[populations.get(i)[j]]);
			}
		}
		
		for(int i = 0; i < xfilters.size(); i++) {
			for(int j =0; j < xfilters.get(i).length; j++) {
				userWriter.println("Xfilters #" + (i+1) + " and header index #" + (j+1) + ": " + xfilters.get(i)[j] + ": " + scroller[xfilters.get(i)[j]]);
			}
		}
		for(int i = 0; i < xfilterc.size(); i++) {
			for(int j =0; j < xfilterc.get(i).length; j++) {
				userWriter.println("Xfilter Count Column #" + (i+1) + " and header index #" + (j+1) + ": " + xfilterc.get(i)[j] + ": " + scroller[xfilterc.get(i)[j]]);
			}
		}
		userWriter.println("Use External Databases?: " + extData);
		for(int i = 0; i < xfilterExt.size(); i++) {
			userWriter.println("Xfilter External Databases" + (i+1) + ": " + xfilterExt.get(i) + ": " + scroller[xfilterExt.get(i)]);
		}
		userWriter.println("Rejection Threshold: " + rejectionThreshold);
		userWriter.println("Threads: " + threads);
		userWriter.println("X-link Cutoff Count: " + cutcnt2);
		
		userWriter.println();
		userWriter.println("Variant Exclusion UI");
		for(int i = 0; i < caddThresh.length; i++) {
			
			userWriter.println("Cadd Threshold #" + i + ": " + caddThresh[i]);
			
		}
		
		userWriter.println("Hom Var Threshold: " + homVarThreshold);
		userWriter.println("Het Threshold: " + hetThreshold);
		userWriter.println("Hem Var Threshold: " + hemVarThreshold);
		userWriter.println("ROC Header Config File: " + ROCHeaderConfig);
		for(int i = 0; i < genotypeThreshList.size(); i++) {
			userWriter.println("Genotype Threshold List #" + (i+1) + ": " + genotypeThreshList.get(i));
		}
		userWriter.println("Genotype Threshold: " + genotypeThresh);
		userWriter.println("Base Quality score cutoff: " + baseFilter);
		userWriter.println("Number of Bad Base Qualities for a Bad Read: " + baseQualThresh);
		userWriter.println("Bad Read Ratio: " + badReadRatio);
		
		
		userWriter.println();
		userWriter.println("Confetti UI");
		userWriter.println("Hom Ref Index: " + homRefIndex + ": " + scroller[homRefIndex]);
		userWriter.println("Hom Var Index: " + homVarIndex + ": " + scroller[homVarIndex]);
		userWriter.println("Genotype Index: " + genotypeIndex + ": " + scroller[genotypeIndex]);
		
		userWriter.println();
		userWriter.println("CNC Filter UI");
		userWriter.println("Cluster Size: " + clusterSize);
		userWriter.println("Use CNC Exon Filtering?: " + CNCExon);
		userWriter.println("Exon Boundaries Config Path: " + exonBoundConfig);
		
		userWriter.close();
		
	}
	
}
