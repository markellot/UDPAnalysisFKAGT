package general;

import java.awt.Frame;
import java.io.File;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.filechooser.FileNameExtensionFilter;

public class TXTFile {
	
	/*
	 * For the user welcome screen for each program.
	 */
	public static void descriptions() {
		//Fix THIS LATER
		String intro = "Hello, thank you for using the Black Magik Tool Box Program. \n"
				+ "This program is a collection of the seven different modules for the Forwards Backwards Analysis Automated Programs developed by Kayla Gu, Anchi Wu, Faris Aziz, and Tom Markello. \nThis window offers brief"
				+ "explanations of how all of these programs work, and how to use them.";
		String desc1 = "\n\nEthnicity Matcher Module.\n"
				+ "This module determines the ethnicity for the given samples. \n"
				+ "You will be prompted to enter the sample file that contains the ID to be ethnically determined and the genotype quality filter criteria. \n"
				+ "Please include in the same folder as this program the following marker folders: \n" + "EASMarker \n"
				+ "AMRMarker \n" + "EURMarker \n" + "AFRMarker \n" + "PairWiseMarker \n"
				+ "Please also include under the same folder the file that contains the paths to the bam files for all the sample ID under study. \n"
				+ "The bam config file should be named EthnicityMatching_Bam.txt.\n"
				+ "The result file and the report file will be saved under the same folder of the input sample file.";
		
		String desc2 = "\n\nSalvage Pathway Module\n"
				+ "This program makes perfamily file whilst correcting false positive Mendelian Inconsistencies. \n"
				+ "You will be prompted to enter the sample file that contains the pedigrees whose perfamilies are to be made \n"
				+ "and the VarSifter files from which the perfamilies will be made. \n" + "\n"
				+ "The perfamily is automatically named as ProbandName+Perfamily.vs. A report file is named as Report_ProbandName+Perfamily.vs\n"
				+ "Both will be saved under the same folder as the input VarSifter file.\n" + "\n"
				+ "Please include under the same folder the file that contains the paths to the bam files for all the sample ID under study. \n"
				+ "The bam config file should be named Salvage_Bam.txt.\n" + "\n"
				+ "Please also include under the same folder the file that contains the ethnicities for the parents. \n"
				+ "Such ethnicity file should have no header and contains only the ID and the corresponding ethnicity separated by tab. \n"
				+ "For example: UDP984 EUR \n"
				+ "The ethnicity config file should be named EthnicityMatchingResult_Salvage.txt.";
		
		String desc3 = "\n\nKayla Kode Module\n This program analyzes the "
				+ "following inheritance models: \n" + "Compound heterozygous recessive \n"
				+ "Homozygous recessive \n" + "Hemizygous \n" + "X-Linked \n" + "De novo dominant \n"
				+ "Mendelian inconsistent \n \n"
				+ "Please make sure a strength configuration and gene configuration file for the \n"
				+ "compound het filter are available. You will be prompted to choose and enter a pedigree file \n"
				+ "that contains the families to be analyzed, a directory file containing the \n"
				+ "locations of the family's BAM files, the name of the output file, and both the \n"
				+ "aforementioned strength and gene config files. The output of this program is a VarSifter file";
		
		String part2 = "Part Two";
		
		String desc4 = "\n\nBroad BAM Curation/SNR calculator program. \n"
				+ "This program is the 5th piece of the full Forwards Backwards pipeline, and evaluates \n"
				+ "the variant density in the region of interest to calculate SNR and error values for \n"
				+ "the variant in question. \n \n"
				+ "You will be prompted to choose and enter a VarSifter file (ideally post-Kayla Kode). \n"
				+ "You will also be prompted to choose and enter a pedigree file that contains the families to be analyzed, \n"
				+ "as well as a directory file containing the locations of the family's BAM files. \n"
				+ "The output of this program is a VarSifter file, with two new columns 'SNR' and 'Error', \n"
				+ "and should have '_BamShortOutput' appended to the input file name.";
		
		String desc5 = "\n\nVariant Exclusion Filter (ROC filter)\n"
				+ "This program is the 6th piece of the full Forwards Backwards pipeline, and evaluates \n"
				+ "all variants identified earlier by the KaylaKode and filters them based on population frequency, \n"
				+ " broad SNR (calculated earlier in Step 5), and predictions of deleteriousness. " + "\n"
				+ "You will be prompted to choose and enter a VarSifter file (ideally post-ROC Filter). \n"
				+ "You will also be prompted to choose and enter a pedigree file that contains the families to be analyzed, \n"
				+ "a directory file containing the locations of the family's BAM files, a deleteriousness scoring threshold \n"
				+ "specific population genotype maximum values, the output file's name, the header configuration file, \n"
				+ "and a population genotype minimum for all genotypes. \n\n"
				+ "The outputs of this program are two VarSifter files. The first, in which all variants in the input \n"
				+ "VS file will only be annotated with flags (i.e. 'HP'= 'high population') based on the above criteria, \n"
				+ "appends '_BamCurated' to the end of the input file name. The second, in which variants are actually \n"
				+ "filtered based on those above flags, appends '_ROCFiltered' to the end of the input file name.";
		
		String part3 = "Part Three";
		String part4 = "Part Four";
		
		
		String desc6 = "\n\nThe pedigree-aware, multiparametric BAM file noise evaluator (Confetti Filter).\n"
				+ "This program is the 7th piece of the full Forwards Backwards pipeline, and evaluates \n"
				+ "potential de novo's identified earlier by the KaylaKode to determine if they are true \n"
				+ "positive de novo's. \n" + "\n"
				+ "You will be prompted to choose and enter a pedigree file that contains the families to be analyzed. \n"
				+ "You will also be prompted to choose and enter a VarSifter file (ideally post-ROC Filter), \n"
				+ "as well as a directory file containing the locations of the family's BAM files. \n"
				+ "The output of this program is a VarSifter file that should have '_DN.vs' appended to the input file name.";
		
		String desc7 = "\n\nCall no Call (CNC, extreme novel deleted exon) filter. \n"
				+ "This program is the 8th and final piece of the full Forwards Backwards pipeline, and evaluates \n"
				+ "potential CNC's identified earlier by the KaylaKode to determine if they are true \n"
				+ "positive CNC's. \n" + "\n"
				+ "You will be prompted to choose and enter a pedigree file that contains the families to be analyzed. \n"
				+ "You will also be prompted to choose and enter a VarSifter file (ideally post-Confetti Filter), \n"
				+ "as well as a directory file containing the locations of the family's BAM files. \n"
				+ "The output of this program is a VarSifter file that should have '_CNC.vs' appended to the input file name.";
		
		JOptionPane.showMessageDialog(new Frame("Introduction"), (intro + desc1 + desc2));
		JOptionPane.showMessageDialog(new Frame("Introduction"), (part2 + desc3 + desc4));
		JOptionPane.showMessageDialog(new Frame("Introduction"), (part3 + desc5 + desc6));
		JOptionPane.showMessageDialog(new Frame("Introduction"), (part4 + desc7));
	}
	

	/*
	 * User prompt to input input sample file.
	 */
	/*
	public static File getSampleFile() {
		JFileChooser browseTo = new JFileChooser(); // creates a file chooser object
		FileNameExtensionFilter filter = new FileNameExtensionFilter("Text Files", "txt");
		browseTo.setFileFilter(filter); // limits the viewable files to .vs and .txt

		JOptionPane.showMessageDialog(new Frame("Input prompt"),
				"Select an input sample file (in .txt format)." + "\n"
						+ "The file should contain the ID number of the samples to be ethnically determined.\n"
						+ "One per each line. For example:" + "\n" + "UDP984" + "\n" + "UDP2495" + "\n" + "UDP2496");
		int returnVal = browseTo.showOpenDialog(new JPanel()); // opens a window for user to browse to vs file to be
																// edited

		if (returnVal == JFileChooser.APPROVE_OPTION) {
			System.out.println(browseTo.getSelectedFile().toString());
			return browseTo.getSelectedFile(); // gets the path of the selected file
		} else {
			int i = wantToContinue("No Sample file selected.");
			if (i == 1) {
				System.exit(0);
			}
			return getSampleFile();
		}
	}
	*/

	/*
	 * Prompt user input for genotype score and coverage cutoff values.
	 */
	public static String getCoverageCutoff() { 
		String cutoff = (String) JOptionPane.showInputDialog(new Frame("Input prompt"), 
				"Enter the score cutoff value and the coverage cutoff value, separated by semicolon. \n"+
		"For example, input 9;14. 9 is the cutoff for genotype score, 14 is the cutoff for coverage."+"\n"+"The suggested cutoff is 9;14."+"\n"
		+"Score cutoff value cannot be above 200. Coverage cutoff cannot be above 500."); 
		if (cutoff == null || cutoff.equals("")) {
			int i = wantToContinue("You did not enter valid cutoff values.");
			if (i == 1) {
				System.exit(0);
			}
			return getCoverageCutoff();
		}
		return cutoff;
	}
	
	/*
	 * Requests user input for VarSifter flat file.
	 */
	public static File getVSFile() {
		JFileChooser browseTo = new JFileChooser(); // creates a file chooser object
		FileNameExtensionFilter filter = new FileNameExtensionFilter("VARSIFTER FILES", "vs", "txt");
		browseTo.setFileFilter(filter); // limits the viewable files to .vs and .txt

		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Select an input VarSifter file.");
		int returnVal = browseTo.showOpenDialog(new JPanel()); // opens a window for user to browse to vs file to be
																// edited
		
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			return browseTo.getSelectedFile(); // gets the path of the selected file
		} else {
			int i = wantToContinue("No VarSifter file selected.");
			if (i == 1) {
				System.exit(0);
			}
			return getVSFile();
		}
	}
	
	/*
	 * Requests user input for pedigree file
	 */
	/*
	public static File getPedigreeFile() {
		JFileChooser browseTo = new JFileChooser(); // creates a file chooser object
		FileNameExtensionFilter filter = new FileNameExtensionFilter("TXT FILES", "vs", "txt");
		browseTo.setFileFilter(filter); // limits the viewable files to .vs and .txt

		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Select an input pedigree file.");
		int returnVal = browseTo.showOpenDialog(new JPanel()); // opens a window for user to browse to vs file to be
																// edited

		if (returnVal == JFileChooser.APPROVE_OPTION) {
			System.out.println(browseTo.getSelectedFile().toString());
			return browseTo.getSelectedFile(); // gets the path of the selected file
		} else {
			int i = wantToContinue("No pedigree file selected.");
			if (i == 1) {
				System.exit(0);
			}
			return getPedigreeFile();
		}
	}
	*/
	
	/*
	 * This class includes exception handles and input of the output file name
	 */
	public static String getDestination() { 
		String fileName = (String) JOptionPane.showInputDialog(new Frame("Input prompt"), 
				"Enter a file name for the output file. This file will be saved in the same folder as the input file, along with a report file.\nDo not include suffixes such as .txt or .vs as those will automatically be generated."); 
		if (fileName == null || fileName.equals("")) {
			int i = wantToContinue("You did not enter a valid file name.");
			if (i == 1) {
				System.exit(0);
			}
			return getDestination();
		}
		return fileName;
	}
	
	
	/*
	 * Requests user input for BAM file directory file
	 */
	public static File getBAMDirectoryFile() {
		JFileChooser browseTo = new JFileChooser(); // creates a file chooser object
		FileNameExtensionFilter filter = new FileNameExtensionFilter("TEXT FILES", "txt", "tsv");
		browseTo.setFileFilter(filter); // limits the viewable files to .txt and .tsv

		JOptionPane.showMessageDialog(new Frame("Input prompt"), "BAM Directory File not found under the Config folder. Please select BAM file directory file.");
		int returnVal = browseTo.showOpenDialog(new JPanel()); // opens a window for user to browse to vs file to be
																// edited

		if (returnVal == JFileChooser.APPROVE_OPTION) {
			//System.out.println(browseTo.getSelectedFile().toString());
			return browseTo.getSelectedFile(); // gets the path of the selected file
		} else {
			int i = wantToContinue("No BAM directory file selected.");
			if (i == 1) {
				System.exit(0);
			}
			return getBAMDirectoryFile();
		}
	}
	
	/*
	 * This is piped so unnecessary
	 * For Salvage specifically. Requests user input for ethnicity config file.
	 *
	public static File getEthConfigFile() {
		JFileChooser browseTo = new JFileChooser(); // creates a file chooser object
		FileNameExtensionFilter filter = new FileNameExtensionFilter("TEXT FILES", "txt", "tsv");
		browseTo.setFileFilter(filter); // limits the viewable files to .txt and .tsv

		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Select ethnicity config file.");
		int returnVal = browseTo.showOpenDialog(new JPanel()); // opens a window for user to browse to vs file to be
																// edited

		if (returnVal == JFileChooser.APPROVE_OPTION) {
			System.out.println(browseTo.getSelectedFile().toString());
			return browseTo.getSelectedFile(); // gets the path of the selected file
		} else {
			int i = wantToContinue("No ethnicity config file selected.");
			if (i == 1) {
				System.exit(0);
			}
			return getEthConfigFile();
		}
	}
	
	*/
	
	/*
	 * For KaylaKode. Requests user input for gene config file.
	 */
	/*
	public static File getGeneConfig() {
		JFileChooser browseTo = new JFileChooser();
		FileNameExtensionFilter gcFilter = new FileNameExtensionFilter("GENE CONFIG FILES", "txt");
		browseTo.setFileFilter(gcFilter); // limits the viewable files to .txt

		int i = JOptionPane.YES_OPTION;

		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Select a Gene Config file.");
		int returnVal = browseTo.showOpenDialog(new JPanel()); // opens a window for user to browse to vs file
																// to be edited

		if (returnVal == JFileChooser.APPROVE_OPTION) {
			i = JOptionPane.NO_OPTION;
			return browseTo.getSelectedFile(); // gets the path of the selected file
		} else {
			i = wantToContinue("No Gene Config file selected.");
			if (i == JOptionPane.NO_OPTION) {
				System.exit(0);
			}
			return getGeneConfig();
		}

	}
	 */
	/*
	 * For KaylaKode. Requests user input for strength file.
	 */
	/*
	public static File getStrengthFile() {
		JFileChooser browseTo = new JFileChooser();
		FileNameExtensionFilter gcFilter = new FileNameExtensionFilter("STRENGTH FILES", "txt");
		browseTo.setFileFilter(gcFilter); // limits the viewable files to .txt

		int i = JOptionPane.YES_OPTION;

		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Select a strength file.");
		int returnVal = browseTo.showOpenDialog(new JPanel()); // opens a window for user to browse to vs file
																// to be edited

		if (returnVal == JFileChooser.APPROVE_OPTION) {
			i = JOptionPane.NO_OPTION;
			return browseTo.getSelectedFile(); // gets the path of the selected file
		} else {
			i = wantToContinue("No strength file selected.");
			if (i == JOptionPane.NO_OPTION) {
				System.exit(0);
			}
			return getStrengthFile();
		}
	}
	*/
	/*
	 * For MakeBamROC. Requests user input for header config file.
	 */
	/*
	public static File getHeaderConfigFile() {
		JFileChooser browseTo = new JFileChooser();
		FileNameExtensionFilter gcFilter = new FileNameExtensionFilter("HEADER CONFIG FILES", "txt");
		browseTo.setFileFilter(gcFilter); // limits the viewable files to .txt

		int i = JOptionPane.YES_OPTION;

		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Select a header config file.");
		int returnVal = browseTo.showOpenDialog(new JPanel()); // opens a window for user to browse to vs file
																// to be edited

		if (returnVal == JFileChooser.APPROVE_OPTION) {
			i = JOptionPane.NO_OPTION;
			return browseTo.getSelectedFile(); // gets the path of the selected file
		} else {
			i = wantToContinue("No header config file selected.");
			if (i == JOptionPane.NO_OPTION) {
				System.exit(0);
			}
			return getHeaderConfigFile();
		}
	}
	*/
	
	/*
	 * Ask user if they would like to continue with the run or abort.
	 */
	public static int wantToContinue(String input) {
		int j = JOptionPane.showConfirmDialog(new Frame("Input prompt"), input + " Would you like to continue?", "User Prompt", 
				JOptionPane.YES_NO_OPTION);
		return j;
	}
	
	/*
	public static File getVSFilefromstring(String loc) {

		File n = new File(loc);
		if (n.isFile() & n.canRead()) {
			return n;
		} else {
			JOptionPane.showMessageDialog(new Frame("Input prompt"),
					"VarSifter file not found at" + loc + ", system exiting.");
			System.exit(1);
		}
		return n;
	}
	*/
	/*
	public static File getVSFilefromstring2(String loc, String np) {

		File n = new File(loc);
		if (n.isFile() & n.canRead()) {
			return n;
		} else {
			JOptionPane.showMessageDialog(new Frame("Input prompt"),
					"VarSifter file not found for " + loc + ", system exiting.");
			System.exit(1);
		}
		return n;
	}
	*/


	/*
	 * public static void salvageDescriptions() { String intro =
	 * "Hello, thank you for using the Salvage Pathway program.\n" + "\n" +
	 * "This program makes perfamily file whilst correcting false positive Mendelian Inconsistencies. \n"
	 * +
	 * "You will be prompted to enter the sample file that contains the pedigrees whose perfamilies are to be made \n"
	 * + "and the VarSifter files from which the perfamilies will be made. \n" +
	 * "\n" +
	 * "The perfamily is automatically named as ProbandName+Perfamily.vs. A report file is named as Report_ProbandName+Perfamily.vs\n"
	 * + "Both will be saved under the same folder as the input VarSifter file.\n" +
	 * "\n" +
	 * "Please include under the same folder the file that contains the paths to the bam files for all the sample ID under study. \n"
	 * + "The bam config file should be named Salvage_Bam.txt.\n" + "\n" +
	 * "Please also include under the same folder the file that contains the ethnicities for the parents. \n"
	 * +
	 * "Such ethnicity file should have no header and contains only the ID and the corresponding ethnicity separated by tab. \n"
	 * + "For example: UDP984 EUR \n" +
	 * "The ethnicity config file should be named EthnicityMatchingResult_Salvage.txt.\n"
	 * ; JOptionPane.showMessageDialog(new Frame("Introduction"), intro); }
	 */

	/*
	 * public static void kaylaKodeDescriptions() { String intro =
	 * "Hello, thank you for using the KaylaKode program. This program analyzes the \n"
	 * + "following inheritance models: \n" + "Compound heterozygous recessive \n" +
	 * "Homozygous recessive \n" + "Hemizygous \n" + "X-Linked \n" +
	 * "De novo dominant \n" + "Mendelian inconsistent \n \n" +
	 * "Please make sure a strength configuration and gene configuration file for the \n"
	 * +
	 * "compound het filter are available. You will be prompted to choose and enter a pedigree file \n"
	 * +
	 * "that contains the families to be analyzed, a directory file containing the \n"
	 * +
	 * "locations of the family's BAM files, the name of the output file, and both the \n"
	 * +
	 * "aforementioned strength and gene config files. The output of this program is a VarSifter file \n"
	 * ; JOptionPane.showMessageDialog(new Frame("KaylaKode Intro"), intro); }
	 */

	/*
	 * public static void bamDescriptions() { String intro =
	 * "Hello, thank you for using the broad BAM Curation/SNR calculator program. \n \n"
	 * +
	 * "This program is the 5th piece of the full Forwards Backwards pipeline, and evaluates \n"
	 * +
	 * "the variant density in the region of interest to calculate SNR and error values for \n"
	 * + "the variant in question. \n \n" +
	 * "You will be prompted to choose and enter a VarSifter file (ideally post-Kayla Kode). \n"
	 * +
	 * "You will also be prompted to choose and enter a pedigree file that contains the families to be analyzed, \n"
	 * +
	 * "as well as a directory file containing the locations of the family's BAM files. \n"
	 * +
	 * "The output of this program is a VarSifter file, with two new columns 'SNR' and 'Error', \n"
	 * + "and should have '_BamShortOutput' appended to the input file name.";
	 * JOptionPane.showMessageDialog(new Frame("Broad BAM Curation Introduction"),
	 * intro); }
	 */

	/*
	 * public static void rocDescriptions() { String intro =
	 * "Hello, thank you for using the population frequency and deleteriousness score threshold \n "
	 * + "(ROC filter). \n" + "\n" +
	 * "This program is the 6th piece of the full Forwards Backwards pipeline, and evaluates \n"
	 * +
	 * "all variants identified earlier by the KaylaKode and filters them based on population frequency, \n"
	 * +
	 * " broad SNR (calculated earlier in Step 5), and predictions of deleteriousness. "
	 * + "\n" +
	 * "You will be prompted to choose and enter a VarSifter file (ideally post-ROC Filter). \n"
	 * +
	 * "You will also be prompted to choose and enter a pedigree file that contains the families to be analyzed, \n"
	 * +
	 * "as well as a directory file containing the locations of the family's BAM files. \n"
	 * +
	 * "The outputs of this program are two VarSifter files. The first, in which all variants in the input \n"
	 * +
	 * "VS file will only be annotated with flags (i.e. 'HP'= 'high population') based on the above criteria, \n"
	 * +
	 * "appends '_BamCurated' to the end of the input file name. The second, in which variants are actually \n"
	 * +
	 * "filtered based on those above flags, appends '_ROCFiltered' to the end of the input file name. \n"
	 * ; JOptionPane.showMessageDialog(new Frame("Confetti Filter Introduction"),
	 * intro); }
	 */

	/*
	 * public static void rocConfigDescriptions() { String intro =
	 * "Hello, thank you for using the population frequency and deleteriousness score threshold \n "
	 * + "(ROC filter). \n" + "\n" +
	 * "This program is the 6th piece of the full Forwards Backwards pipeline, and evaluates \n"
	 * +
	 * "all variants identified earlier by the KaylaKode and filters them based on population frequency, \n"
	 * +
	 * " broad SNR (calculated earlier in Step 5), and predictions of deleteriousness. "
	 * + "\n" +
	 * "You will be prompted to choose and enter a VarSifter file (ideally post-ROC Filter). \n"
	 * +
	 * "You will also be prompted to choose and enter a pedigree file that contains the families to be analyzed, \n"
	 * +
	 * "a directory file containing the locations of the family's BAM files, a deleteriousness scoring threshold \n"
	 * +
	 * "specific population genotype maximum values, the output file's name, the header configuration file, \n"
	 * + "and a population genotype minimum for all genotypes. \n\n" +
	 * "The outputs of this program are two VarSifter files. The first, in which all variants in the input \n"
	 * +
	 * "VS file will only be annotated with flags (i.e. 'HP'= 'high population') based on the above criteria, \n"
	 * +
	 * "appends '_BamCurated' to the end of the input file name. The second, in which variants are actually \n"
	 * +
	 * "filtered based on those above flags, appends '_ROCFiltered' to the end of the input file name. \n"
	 * ; JOptionPane.showMessageDialog(new Frame("Confetti Filter Introduction"),
	 * intro); }
	 */

	/*
	 * public static void ppdnDescription() { String intro =
	 * "Hello, thank you for using the pedigree-aware, multiparametric BAM file noise evaluator \n "
	 * + "(confetti filter). \n" + "\n" +
	 * "This program is the 7th piece of the full Forwards Backwards pipeline, and evaluates \n"
	 * +
	 * "potential de novo's identified earlier by the KaylaKode to determine if they are true \n"
	 * + "positive de novo's. \n" + "\n" +
	 * "You will be prompted to choose and enter a pedigree file that contains the families to be analyzed. \n"
	 * +
	 * "You will also be prompted to choose and enter a VarSifter file (ideally post-ROC Filter), \n"
	 * +
	 * "as well as a directory file containing the locations of the family's BAM files. \n"
	 * +
	 * "The output of this program is a VarSifter file that should have '_DN.vs' appended to the input file name. \n"
	 * ; JOptionPane.showMessageDialog(new Frame("Confetti Filter Introduction"),
	 * intro); }
	 */

	/*
	 * public static void CNCDescriptions() { String intro =
	 * "Hello, thank you for using the Call no Call (CNC, extreme novel deleted exon) filter. \n"
	 * + "\n" +
	 * "This program is the 8th and final piece of the full Forwards Backwards pipeline, and evaluates \n"
	 * +
	 * "potential CNC's identified earlier by the KaylaKode to determine if they are true \n"
	 * + "positive CNC's. \n" + "\n" +
	 * "You will be prompted to choose and enter a pedigree file that contains the families to be analyzed. \n"
	 * +
	 * "You will also be prompted to choose and enter a VarSifter file (ideally post-Confetti Filter), \n"
	 * +
	 * "as well as a directory file containing the locations of the family's BAM files. \n"
	 * +
	 * "The output of this program is a VarSifter file that should have '_CNC.vs' appended to the input file name. \n"
	 * ; JOptionPane.showMessageDialog(new Frame("CNC Filter Introduction"), intro);
	 * }
	 */
}
