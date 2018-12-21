import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import javax.swing.filechooser.FileNameExtensionFilter;

public class UniversalUtility5v6SnpEffPipeLine {
	
	
	//names of the various programs being run, will only affect UI
	private String prog1 = "VarSifter Cleanup";
	//private String prog2 = "HGMD Fixer";
	private String prog3 = "Annotation Primer";
	//private String prog4 = "CADD Input Generator";
	private String prog5 = "Column Deleter";
	//private String prog6 = "PolyPhen Input Generator";
	private String prog7 = "Index Adder";
	
	//indexes of fields that may be indexed depending on the utility being run
	private int chrIndex;
	private int leftFlankIndex;
	private int refAlleleIndex;
	private int varAlleleIndex;
	private int mutTypeIndex;
	private int HGMDtags;
	private int HGMDhgvs;
	private int HGMDyn;
	private int vcfPos;
	private int vcfRef;
	private int vcfVar;  

	//values within fields
	private String refAllele;
	private String varAllele;
	private String leftFlank;
	private String currentChrom = "";
	
	//variables changed by UI
	private boolean[] selections; 
	private File vsFile;
	private String newFileName; 
	private String genomeDirectory; 
	private ArrayList<Integer> toDelete;
	
	//variables changed by initializer
	private BufferedReader vsData;
	private PrintWriter writer;
	private PrintWriter caddWriter;
	private PrintWriter pphenWriter;
	private int caddNum; //the number of the currnt CADD chunk being written
	private int pphenNum; //the number of the current PolyPhen chunk being writtn
	private int columns; 
	private int caddCounter; //counts the number of CADD lines put into chunks
	private int pphenCounter; //counts the number of Polyphen lines put into chunks
	private String path;
	
	//scanner for chromosomes
	private BufferedReader refGenome;
	private String chromLine;
	private int lineLength;
	private int lineNumber; 
	
	private int naPos; //used to find the start point to edit metadata genotypes
	private int naDistanceApart; //used to find the distance between metadata genotypes

	boolean newFile = false; //stores if we need to create a new VS file
	private long newIndex;
	
	private boolean needSnpEff;
	
	//time
	private long bRunner;
	private long aRunner;
	private long aSort;
	private long aThread;
	private String sortTime;
	private String snpEffTime;
	
	//Temp File + Parallelization
	private List<File> tempIn = new ArrayList<File>();
	private List<File> tempOut = new ArrayList<File>();
	private List<BufferedReader> readIn = new ArrayList<BufferedReader>();
	private List<PrintWriter> readOut = new ArrayList<PrintWriter>();
	private ExecutorService executor;
	private final int NTHREDS = Runtime.getRuntime().availableProcessors();
	private int sepChroms = 0;
	private File dir;
	
	
	/*----------------------------------------------------------------------------Beginning of UI Methods--------------------------------
	 * --------------------------------------------------------------------------------------------------------------------------*/	
	
	
	//User Interface
	public void uI() throws IOException {
		descriptions();
		
		vsFile = getVSFile();
		
		selections = utilitySelection();
		//field 0 - Friendly-izing the VS file - var cleanup
		//field 1 - HGMD fixing
		//field 2 - adds CADD predictions/fields - annotation primer
		//field 3 - saves as 100,000 line chunks - cadd input
		//field 4 - deleting columns
		//field 5 - polyphen input
		//field 6 - adds index
		
		
		if (selections[0] || selections[1] || selections[2] || selections[4] || selections[6]) {
			newFileName = getDestination();
			newFile = true;
		}
		
		if (selections[6]) {
			newIndex = 1;
		}
		
		if (selections[0]) {
			genomeDirectory = getChromoDirec();
		}
		
		needSnpEff = needsSnpEff();
		
		if(needSnpEff) {

			SnpEffEr derpy = new SnpEffEr(vsFile, selections[4]);
			snpEffTime = derpy.getTimeMessage();
			toDelete = derpy.getToDelete();
			
			Sorter shorty = new Sorter(derpy.getOutput(), derpy.getChrIndex(), derpy.getLeftFlankIndex());
			
			vsFile = shorty.getOutput();
			sortTime = shorty.getTimeMessage();
			
		} else {
			
			
			snpEffTime = "Not applicable\n\n";
			sortTime = "Not applicable\n\n";
			
			
		}
		
		
	}
	
	//generates a message window describing the available utilities
	public void descriptions() {
		String intro = "This program offers several utilities for corrrecting VarSifter files. This window offers brief \n"
						+ "explanations of how each of these utilities works. When you are ready, click 'OK' and you will be \n"
						+ "prompted to select which utilities you'd like to run. \n\n";
		String desc1 = prog1 + ": \n"
						+ "This utility adjusts the VarSifter to match VCF format requirements, adding the previous base to \n"
						+ "the left side of the reference and variant alleles, as well as decrementing the left flank by 1. \n"
						+ "You will need to provide a folder with a reference genome and an input VarSifter file sorted by \n "
						+ "chromosome and then left flank (in that order). Your VarSifter file must also contain the following \n"
						+ "column headers: Chr, LeftFlank, ref_allele, var_allele, muttype. \n\n";
		
		String desc3 = prog3 + ": \n "
				+ "This utility prepares a VarSifter file to receive CADD scores. This is done by the addition of five \n"
				+ "new fields. The first three account for the prediction of CADD's simplification of Indels by eating \n"
				+ "into superfluous bases, which will allow remerging the files CADD returns. The remaining two fields \n"
				+ "are the fields for the raw and phred scores to be incorporated. \n\n";
		
		String desc5 = prog5 + ": \n"
						+ "This utility allows you to select columns for deletion. \n\n";

		String desc7 = prog7 + ":\n"
						+ "This utility will create an index field if there is none. Use on Appistry files.";
		
		/*
		 * Discontinued as of 6/5/18 - Faris and John

		String desc4 = prog4 + ": \n"
						+ "This utility will write the modified VarSifter file in 100,000 line chunks with only VCF columns in addition to \n"
						+ "the modified large file. These 100,000 line files can then be sent to CADD without violating size restrictions. \n\n";
		
		String desc6 = prog6 + ":\n"
				+ "This utility will generate seperate files for input into Polyphen as batch queries. \n\n";
		
		String desc2 = prog2 + ": \n"
				+ "This utility corrects ambiguities in the HGMD tag field by checking the HGMD Allele against the \n"
				+ "variant for a given line. Your input VarSifter file must have the following column headers: ref_allele \n"
				+ "var_allele, HGMDtags, HGMD_ALLELES_in_GENOME_STRAND. \n\n";
		
		*/
		
		//JOptionPane.showMessageDialog(new Frame("Introduction"), intro + desc1 + desc2 + desc3 + desc4 + desc5 + desc6 + desc7);
		JOptionPane.showMessageDialog(new Frame("Introduction"), intro + desc1 + desc3 + desc5 + desc7);
		
	}
	
	public boolean needsSnpEff() {
		
		int j = JOptionPane.showConfirmDialog(new Frame("Input Prompt"), "Is this a NISC file or a Varsifter file in which SnpEff does not need to be performed?", "Perform SnpEff?", 
				JOptionPane.YES_NO_OPTION);
		
		if(j == 0) {
			
			return false;
		}
		else {
			
			return true;
			
		}
		
		
	}
	
	//creates a checkbox window where the user can select the utilities to run - returns a boolean array where true means run
	public boolean[] utilitySelection() {
		boolean[] selections = new boolean[7];
		selections[0] = false;
		selections[2] = false;

		selections[4] = false;
		
		selections[6] = false;
		
		/*
		 * Discontinued as of 6/5/18 - Faris and John
		selections[1] = false;
		selections[3] = false;
		selections[5] = false;
		
		*/
		
		JCheckBox check1 = new JCheckBox(prog1);
		check1.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				selections[0] = !selections[0];
			}
		});
		
		/*
		 * Discontinued as of 6/5/18 - Faris and John
		 *  
		JCheckBox check2 = new JCheckBox(prog2);
		check2.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) {
				selections[1] = !selections[1];
			}
		});
		*/
		JCheckBox check3 = new JCheckBox(prog3);
		check3.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				selections[2] = !selections[2];
			}
		});
		/*
		 * Discontinued as of 6/5/18 - Faris and John
		 * 
		JCheckBox check4 = new JCheckBox(prog4);
		check4.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {`
				selections[3] = !selections[3];
			}
		});
		*/
		
		JCheckBox check5 = new JCheckBox(prog5);
		check5.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				selections[4] = !selections[4];
			}
		});
		/*
		 * Discontinued as of 6/5/18 - Faris and John
		 * 
		JCheckBox check6 = new JCheckBox(prog6);
		check6.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				selections[5] = !selections[5];
			}
		});
		*/
		JCheckBox check7 = new JCheckBox(prog7);
		check7.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				selections[6] = !selections[6];
			}
		});
		
		String message = "Please select the utilities you would like to run.";
		
		//Object[] params = {message, check1, check2, check3, check4, check5, check6, check7};
		Object[] params = {message, check1, check3, check5, check7};
		int j = JOptionPane.showConfirmDialog(new Frame("Input prompt"), params, "Disconnect Products", JOptionPane.OK_CANCEL_OPTION);
		
		if (j == 2) {
			System.exit(0);
		}
		
		return selections;
	}
	
	//creates a FileChooser that allows the selection of an input VS File
	public File getVSFile() {
		JFileChooser browseTo = new JFileChooser(); //creates a file chooser object
		FileNameExtensionFilter filter = new FileNameExtensionFilter("VARSIFTER FILES", "vs", "txt"); 
		browseTo.setFileFilter(filter); //limits the viewable files to .vs and .txt
		
		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please select an input VarSifter file.");
		int returnVal = browseTo.showOpenDialog(new JPanel()); //opens a window for user to browse to vs file to be edited
		
	    if(returnVal == JFileChooser.APPROVE_OPTION) {
	    	return browseTo.getSelectedFile(); //gets the path of the selected file
	     } else {
				int i = this.wantToContinue("No VarSifter file selected."); 
				if (i == 1) {
					System.exit(0);
				}
				return this.getVSFile();
			}	
	}
	
	//gets the name of the new file from the user
	public String getDestination() { //prompts the user to enter the name of a destination for the edited file
		String fileName = (String) JOptionPane.showInputDialog(new Frame("Input prompt"), 
				"Enter a file name for the edited VarSifter. This file will be saved in the same folder as the input file."); 
		if (fileName == null) {
			int i = this.wantToContinue("You did not enter a valid file name.");
			if (i == 1) {
				System.exit(0);
			}
			return this.getDestination();
		} else if (fileName.equals("")) {
			int i = this.wantToContinue("You did not enter a valid file name.");
			if (i == 1) {
				System.exit(0);
			}
			return this.getDestination();
		}
		return fileName;
	}
	
	//creates a FileChooser that allow the selection of  directory containing a reference genome
	public String getChromoDirec() {
		JFileChooser browseTo = new JFileChooser(); //creates a file chooser object
		browseTo.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY); //limits selectable objects to directories
		
		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please select the directory containing the " +
				"reference genome in seperate files by chromosome.");
		int returnVal = browseTo.showOpenDialog(new JPanel()); //opens a window for user to browse to directory containing reference genome
		
	    if(returnVal == JFileChooser.APPROVE_OPTION) {
	    	return browseTo.getSelectedFile().getPath(); //gets path of the selected directory
	    } else {
			int i = this.wantToContinue("No directory selected."); 
			if (i == 1) {
				System.exit(0);
			}
			return this.getChromoDirec();
		}	
	}
/*----------------------------------------------------------------------------End of UI Methods--------------------------------
 * --------------------------------------------------------------------------------------------------------------------------*/	
	
	
/*----------------------------------------------------------------------------Beginning of Initializer Methods--------------------------------
* --------------------------------------------------------------------------------------------------------------------------*/	
		
	public void initializer() throws IOException {

		try {																//input VarSifter file - file provided by user
			vsData = new BufferedReader(new FileReader(vsFile));		
			if (newFile) {
				writer = new PrintWriter(vsFile.getParent() + "\\" + newFileName + ".vs", "UTF-8");
			}
			
			//path = vsFile.getPath().substring(0,vsFile.getPath().lastIndexOf(vsFile.getName()));
			
			/*
			 * Discontinued as of 6/6/18 - Faris and John
			String temp = vsFile.getName().substring(0, vsFile.getName().lastIndexOf("."));
			
			if (selections[3]) {
				caddNum = 1;
				caddCounter = 0;
						

				(new File(path + "\\CADDInputs_" + vsFile.getName().substring(0, vsFile.getName().indexOf(".")))).mkdirs();
				caddWriter = new PrintWriter(new File(path + "\\CADDInputs_" + temp + "\\" + temp + "_result" + Integer.toString(caddNum)) + ".txt", "UTF-8");		//output VCF file - file provided by user
				caddWriter.println("##fileformat=VCFv4.0");									//print first header line
				caddWriter.println("#CHROM\tPOS\tID\tREF\tALT");
			}
			
			if (selections[5]) {
				pphenNum = 1;
				pphenCounter = 0;
				
				(new File(path + "\\PolyPhenInputs_" + vsFile.getName().substring(0, vsFile.getName().indexOf(".")))).mkdirs();
				pphenWriter = new PrintWriter(new File(path + "\\PolyPhenInputs_" + temp + "\\" + temp + "_result" + Integer.toString(pphenNum)) + ".txt", "UTF-8");
			}
			
			*/
		} catch (FileNotFoundException error) {
			JOptionPane.showMessageDialog(new Frame("Error"), "The VarSifter file could not be found, possibly because it is open elsewhere. Please close and try again.");
			System.exit(0);
		}
		
		String curLine = vsData.readLine();
		String[] curLineSplit = curLine.split("\t");
		List<String> headers = java.util.Arrays.asList(curLineSplit); 
		columns = curLineSplit.length;

		genotypeIndices(curLineSplit);									//finds the start of the metadata - necesary to change genotypes or add CADD columns
		
		
		if ((selections[2])) {
			vcfPos = naPos;
			vcfRef = naPos+1;
			vcfVar = naPos+2;
		} 
		/*
		 * Discontinued as of 6/6/18 - John and Faris
		else if (selections[3]) {
			vcfPos = headers.indexOf("VCF_Position");
			vcfRef = headers.indexOf("VCF_Ref");
			vcfVar = headers.indexOf("VCF_Var");
			
			if ((vcfPos == -1) || (vcfPos == -1) || (vcfPos == -1)) {
				JOptionPane.showMessageDialog(new Frame("Error Message"), "You are missing the necessary annotation fields. Please run the Annotation Primer Utility.");
				System.exit(0);
			}
		}
		*/
		
		
		if (newFile) {
			if (selections[6]) {
				curLine = "Index\t" + curLine; 
			}		
		}
		
		if (selections[2]) {				//adds the CADD columns if that's what the user wants
			curLine = "";
			curLineSplit = annotater(curLineSplit, "VCF_Position\tVCF_Ref\tVCF_Var\t");
			for (int i = 0; i < curLineSplit.length; i++) { //reforms curLine
				curLine = curLine.concat(curLineSplit[i]).concat("\t");
			}
			curLine.concat(curLineSplit[columns-1]);
			headers = java.util.Arrays.asList(curLineSplit); 
		}

		chrIndex = headers.indexOf("Chr");									//finds indices of the five headers to the left
		leftFlankIndex = headers.indexOf("LeftFlank");						//if header does not exist, String.indexOf returns
		refAlleleIndex = headers.indexOf("ref_allele");						//a -1
		varAlleleIndex = headers.indexOf("var_allele");						
		mutTypeIndex = headers.indexOf("muttype");
		
		/*
		HGMDtags = headers.indexOf("HGMDtags");
		HGMDhgvs = headers.indexOf("HGMD_HGVS");
		HGMDyn = headers.indexOf("REF_MATCHES_HGMD_DISEASE_ALLELE");
		*/
		
		
		boolean indexError = false;
		String errorMessage =  "";
		if (chrIndex == -1) {
			errorMessage += "The header 'Chr' appears no where in the header line.\n";   			//missing 'Chr' header
			indexError = true;
		}
		if (leftFlankIndex == -1) {
			errorMessage += "The header 'LeftFlank' appears no where in the header line.\n";		//missing 'LeftFlank' header
			indexError = true;
		}
		if (refAlleleIndex == -1) {
			errorMessage += "The header 'ref_allele' appears no where in the header line.\n";		//missing 'ref_allele' header
			indexError = true;
		}
		if (varAlleleIndex == -1) {
			errorMessage += "The header 'var_allele' appears no where in the header line.\n";		//missing 'var_allele' header
			indexError = true;
		}
		if (mutTypeIndex == -1) {
			errorMessage += "The header 'muttype' appears no where in the header line.\n";		//missing 'muttype' header
			indexError = true;	
		}
		/*
		if ((selections[1]) && (HGMDtags == -1)) {
			errorMessage += "The header 'HGMDtags' appears no where in the header line.\n";		//missing 'muttype' header
			indexError = true;	
		}
		if ((selections[1]) && (HGMDhgvs == -1)) {
			errorMessage += "The header 'HGMD_HGVS' appears no where in the header line.\n";		//missing 'muttype' header
			indexError = true;	
		}
		if ((selections[1]) && (HGMDyn == -1)) {
			errorMessage += "The header 'REF_MATCHES_HGMD_DISEASE_ALLELE' appears no where in the header line.\n";		//missing 'muttype' header
			indexError = true;	
		}
		*/
		
		if (indexError) {																	//if any of the above if statements are true, 
			errorMessage += "As stated above. One or more of the necessary headers do not appear"	//then the header does not appear within the
						  + " in your VarSifter file. \nPlease ensure these headers are in "		//header line. Thus, we print the missing header(s)
					      + "your VarSifter file prior to running this program.";					//along with the statement below and the program exits.
			JOptionPane.showMessageDialog(new Frame("Error Message"), errorMessage);
			System.exit(0);
		}	
		
		if (selections[4]) { //if the user wants to delete
			
			if(!needSnpEff) {
				
				getDeleteColumns(curLineSplit);
				
			}
			columns = columns - toDelete.size(); //adjust the number of columns as appropriate
			
			curLineSplit = deleter(curLineSplit); //deletes the rows from curLineSplit
			curLine = "";
			for (int i = 0; i < curLineSplit.length; i++) { //reforms curLine
				curLine = curLine.concat(curLineSplit[i]).concat("\t");
			}
			curLine.concat(curLineSplit[columns-1]);
			headers = java.util.Arrays.asList(curLineSplit);  
		}
		
		if (selections[6]) {
			
			curLine = "Index\t" + curLine;
			
		}
		
		if (newFile) {
			
			writer.println(curLine.trim());
			
		}
		
		
	}
	
	//gets the start of the metadata (naPos) and the distance between each of the genotypes
	public void genotypeIndices(String[] lineSplit) {
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
						naPos = counter;
					}
					else {
						naDistanceApart = counter - naPos;
					}
				tempVal++;			
				}
				counter++;
			}
		}
	}
	
	//adds the CADD columns in the appropriate place if called
	public String[] annotater(String[] lineSplit, String newFields) {		
		String newLine = "";
		for (int i=0; i < vcfPos; i++) {						//for loop concatenates all but last field, with field followed by tab			
			newLine = newLine.concat(lineSplit[i]).concat("\t");
		}
		newLine = newLine.concat(newFields);
		
		for (int j=vcfPos; j < lineSplit.length-1; j++) {
			newLine = newLine.concat(lineSplit[j]).concat("\t");
		}
		newLine = newLine.concat(lineSplit[lineSplit.length-1]);
		
		return newLine.split("\t");					//prints the header line of the VarSifter file
	}
	
	
	//prompts the user to enter the headers of columns they want to delete
	public void getDeleteColumns(String[] curLineSplit) {
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("Column Deleter");
		jframe.setSize(((curLineSplit.length - 7) / 15) * 7 + 50, 400);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		final JDialog diag = new JDialog(jframe, "Delete Columns", true);		
		
		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);
		
		JLabel jtext = new JLabel("Please select which columns you'd like to delete.", SwingConstants.CENTER);
		gbc.gridx = 0;
		gbc.gridy = 0;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(5, 5, 5, 5);
		gbc.weightx = 1.0;
		jpanel.add(jtext, gbc);
		
		JCheckBox tempBox;
		boolean[] boxSelections = new boolean[naPos];
		int xPos = 0;
		int yPos = 1;
		
		for (int i = 1; i < naPos-7; i++) {
			final int x = i+7;
			
			String temp; 
			if (curLineSplit[x].length() > 25) {
				temp = curLineSplit[x].substring(0,24).concat("...");
			} else {
				temp = curLineSplit[x];
			}
			
			tempBox = new JCheckBox(temp, false);
			tempBox.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					boxSelections[x] = !boxSelections[x];
				}
				});
			gbc.gridx = xPos;
			gbc.gridy = yPos;
			gbc.fill = GridBagConstraints.HORIZONTAL;
			gbc.anchor = GridBagConstraints.WEST;
			gbc.insets = new Insets(5, 5, 5, 5);
			gbc.weightx = 1.0;
			jpanel.add(tempBox, gbc);
			
			yPos++;
			if (yPos >= 20) {
				yPos = 1;
				xPos++;
			}
		}
		
		JButton okButton = new JButton("Ok");
		okButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				diag.dispose();
			}
		});
		
		gbc.gridx = 0;
		gbc.gridy = 21;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(5, 5, 5, 5);
		gbc.weightx = 1.0;
		jpanel.add(okButton, gbc);
		
		JButton cancelButton = new JButton("Cancel");
		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				diag.dispose(); 
				int i = wantToContinue(""); 
				if (i == 1) {
					System.exit(0);
				}
				getDeleteColumns(curLineSplit);
			}
		});
		gbc.gridx = 1;
		gbc.gridy = 21;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(5, 5, 5, 5);
		gbc.weightx = 1.0;
		jpanel.add(cancelButton, gbc);
		
		diag.getContentPane().add(jpanel);
		diag.pack();
		diag.setLocationRelativeTo(null);
		diag.setVisible(true);
		
		toDelete = new ArrayList<Integer>();
		for (int j = 0; j < boxSelections.length; j++) {
			if (boxSelections[j]) {
				toDelete.add(j);
			}
		}
	}
	
	//deletes columns based on the indices in the toDelete matrix
	public String[] deleter(String[] curLineSplit) {
		ArrayList<String> temp = new ArrayList<String>(Arrays.asList(curLineSplit));
		for (int j = toDelete.size() - 1; j >= 0; j--) {
			temp.remove((int) toDelete.get(j));
		}
		curLineSplit = new String[temp.size()];
		curLineSplit = temp.toArray(curLineSplit);
		return curLineSplit; 
	}
	
	
	
	
/*----------------------------------------------------------------------------END of Initializer Methods--------------------------------
* --------------------------------------------------------------------------------------------------------------------------*/	
		

/*----------------------------------------------------------------------------Beginning of Runner Methods--------------------------------
* --------------------------------------------------------------------------------------------------------------------------*/	
	
	
	public void runner() throws IOException, InterruptedException {
		initializer();
		
		//UI canceler
		JFrame jframe = canceler();
		jframe.setLocationRelativeTo(null);
		jframe.setVisible(true);
		
		
		//time
		bRunner = System.currentTimeMillis(); //Check 0
		

		//created executor to organize threads
		executor = Executors.newWorkStealingPool(NTHREDS);
		
		//created necessary sorts for threading
		System.out.println("B4 Sortin Temp");
		createSortIntoTemp();
		
		aSort = System.currentTimeMillis();
		
		System.out.println("After Sorting Temp + Creating Threads");
		System.out.println("Number of seperate Chromosome Threads: " + sepChroms);
		
		executor.shutdown();
		
		//waiting for threads to be finished
		try {
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		System.out.println("Threads all finished");
		
		//closing readers
		System.out.println("Closing Readers");
		for(int i = 0; i < sepChroms; i++) {
			
			readIn.get(i).close();
			readOut.get(i).close();
			
		}
		
		aThread = System.currentTimeMillis();
		
		//merging files
		System.out.println("Merging Files");
		String out = "";
		for(int i = 0; i < sepChroms; i++) {
			
			readIn.set(i, new BufferedReader(new FileReader(tempOut.get(i))));
			
			while((out = readIn.get(i).readLine()) != null) {
				
				writer.println(out);
			}
		}
		
		aRunner = System.currentTimeMillis();
		
		//Ending
		
		//deleting files and closing stuff
		for(int i = 0; i < sepChroms; i++) {
			
			readIn.get(i).close();
			readOut.get(i).close();
			tempIn.get(i).delete();
			tempOut.get(i).delete();	
			
		}
		
		dir.delete();
		
		jframe.dispose();
		
		//ending screen and time
		

		JFrame finalFrame = new JFrame("Completed");
		
		finalFrame.setUndecorated(true);
		finalFrame.setVisible(true);
		finalFrame.setLocationRelativeTo(null);
		finalFrame.setAlwaysOnTop(true);
		
		String timeMessage = ("The Line By Line Utilities Program Took " + (aRunner - bRunner) + "ms ("
				+ ((aRunner - bRunner) / 1000) + "." + (((aRunner - bRunner) % 1000) / 100) + " seconds)\n\nThe Splitter Took: " + (aSort - bRunner) + "ms (" 
				+ ((aSort - bRunner) / 1000 + "." + (((aSort - bRunner) % 1000) / 100) + " seconds)\nThe Parallelization Analyzer (after the Splitter)  Took: " + (aThread - aSort) + "ms ("
				+ ((aThread - aSort) / 1000 + "." + (((aThread - aSort) % 1000) / 100) + " seconds)\nThe Merger Took: " + (aRunner - aThread) + "ms ("
				+ (aRunner - aThread) / 1000 + "." + (((aRunner - aThread % 1000) / 100 + "seconds)")))));
		
		
		JOptionPane.showMessageDialog(finalFrame, ("Completed. \n\n") + snpEffTime + "\n\n" + sortTime + "\n\n" + timeMessage);
		vsData.close(); 
		
		if (writer != null) {
			writer.close(); 
		}
		
		if (caddWriter != null) {
			caddWriter.close();
			
		}
		
		if (pphenWriter != null) {
			pphenWriter.close();
		}
		
		if (refGenome != null) {
			refGenome.close();
		}
			
		System.gc();
		System.exit(0);
	}
	
	
	
	
	public JFrame canceler() {
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("Line-by-Line Utilities");
		jframe.setSize(500, 100);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);
		
		JLabel jtext = new JLabel("Line By Line Program running. If you'd like to abort, hit the cancel button below.", SwingConstants.CENTER);
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
				try {
					vsData.close();

					//deleting files and closing stuff
					for(int i = 0; i < sepChroms; i++) {
						
						tempIn.get(i).delete();
						tempOut.get(i).delete();
						readIn.get(i).close();
						readOut.get(i).close();
						
					}
					dir.delete();
					
				} catch (IOException e1) {
				}
				
				writer.close();
				if (caddWriter != null) {
					caddWriter.close();
				}
				if (refGenome != null) {
					try {
						refGenome.close();
					} catch (IOException e1) {
					}
				} 
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
	
	//creates files and sorts original VarSifter file by chromosome and sets readers and writers appropriately
	public void createSortIntoTemp() throws IOException {
		String line1;
		String line2;
		String[] lineArr1;
		String[] lineArr2;
		
		//reads first two lines
		line1 = vsData.readLine();
		line2 = vsData.readLine();
		
		//Directory of where temp files go

		int suffix = vsFile.getName().lastIndexOf(".");
		if(suffix == -1) {
			
			suffix = vsFile.getName().length();
			
		}
		dir = new File(vsFile.getParent() + "\\LBLTempFiles" + vsFile.getName().substring(0, suffix));
		dir.mkdir();
		dir.deleteOnExit();
		
		//loop until end of varsifter file
		while(line2 != null) {
			
			//Creating new temp files
			tempIn.add(File.createTempFile("Part " + Integer.toString(sepChroms+1), ".txt", dir));
			tempIn.get(sepChroms).deleteOnExit();
			
			tempOut.add(File.createTempFile("Part " + Integer.toString(sepChroms+1), ".txt", dir));
			tempOut.get(sepChroms).deleteOnExit();
			
			readOut.add(new PrintWriter(tempIn.get(sepChroms)));
			
			
			
			
			do {
				
				//prints to temp 
				readOut.get(sepChroms).println(line1);
				
				if(line2 == null) {
					
					break;
					
				}
				
				
				lineArr1 = line1.split("\t");
				lineArr2 = line2.split("\t");
				
				//checks for change in chromosome
				if(!(lineArr1[chrIndex].equals(lineArr2[chrIndex]))) {
					
					/* removed for time efficiency
					 * 
					if(lineArr1[chrIndex].length() >= 5 && lineArr2[chrIndex].length() >= 5) {
						
						if(!(lineArr1[chrIndex].substring(0, 5).equals(lineArr2[chrIndex].substring(0, 5)))) {
							
							line1 = line2;
							line2 = vsData.readLine();
							break;
							
						}
						
					} else if(!(lineArr1[chrIndex].substring(0, 4).equals(lineArr2[chrIndex].substring(0, 4)))) {
						
						line1 = line2;
						line2 = vsData.readLine();
						break;
						
					}
					*/
					
					//changes line and exits loop
					line1 = line2;
					line2 = vsData.readLine();
					break;
					
					
				}
				//changes line
				line1 = line2;
				line2 = vsData.readLine();
				
			} while(true);
			
			//sets buffered readers and writers and adds to new chromosome numbers
			readIn.add(new BufferedReader(new FileReader(tempIn.get(sepChroms))));
			readOut.get(sepChroms).close();
			readOut.set(sepChroms, new PrintWriter(tempOut.get(sepChroms)));
			
			//Starts Thread Execution After Temp File is made
			Runnable worker = new Analysisv2(readIn.get(sepChroms), readOut.get(sepChroms), selections, genomeDirectory, chrIndex, leftFlankIndex, refAlleleIndex, varAlleleIndex,
					mutTypeIndex, newIndex, naPos, naDistanceApart, vcfPos, toDelete);
			executor.execute(worker);
		
			
			sepChroms++;
			
		}
	}
	
	

/*----------------------------------------------------------------------------END of Runner Methods--------------------------------
* --------------------------------------------------------------------------------------------------------------------------*/	
	
	
	//prompts user if they want to continue if they've hit cancel
	public int wantToContinue(String input) {
		int j = JOptionPane.showConfirmDialog(new Frame("Input prompt"), input + " Do you want to continue?", "User Prompt", 
				JOptionPane.YES_NO_OPTION);
		return j;
	}
	
	public static void main(String[] args) throws IOException, InterruptedException {
		
		UniversalUtility5v6SnpEffPipeLine derp = new UniversalUtility5v6SnpEffPipeLine();
		derp.uI();
		derp.runner();
		

	}

}
