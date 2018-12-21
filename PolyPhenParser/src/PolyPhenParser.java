import java.awt.Frame;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.List;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.filechooser.FileNameExtensionFilter;

public class PolyPhenParser {
	//generates and edit file compatable with the editor program from PolyPhen2 output
	
	private File input; //input PPhen file
	private String newFileName; //filename for resulting edit file
	
	private BufferedReader inputData; //reads the input file
	private PrintWriter writer; //writers out the edit file
	private File unsorted;
	
	private String inputLine;
	private String[] inputLineSplit;
	private String output = "";
	
	private int protPosIndex;
	private int refAAIndex;
	private int altAAIndex;
	private int predictionIndex;
	private int probIndex;
	private int sensIndex; //TPR
	private int specIndex; //1-FPR
	private int locIndex;
	
	private long aRunner;
	private long bRunner;
	private String timeMessage;
	
	public void uI() {
		input = getInputFile(); 
		newFileName = getDestination();
	}
	
	//prompts the user to give a PolyPhen input file
	public File getInputFile() {
		JFileChooser browseTo = new JFileChooser(); //creates a file chooser object
		FileNameExtensionFilter filter = new FileNameExtensionFilter("Text Files", "txt"); 
		browseTo.setFileFilter(filter); //limits the viewable files to .txt
		
		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please select a PolyPhen output file as input.");
		int returnVal = browseTo.showOpenDialog(new JPanel()); //opens a window for user to browse to vs file to be edited
		
	    if(returnVal == JFileChooser.APPROVE_OPTION) {
	    	return browseTo.getSelectedFile(); //gets the path of the selected file
	     } else {
				int i = this.wantToContinue("No input file selected."); 
				if (i == 1) {
					System.exit(0);
				}
				return this.getInputFile();
			}	
	}
	
	//prompts the user to enter the name of a destination for the edited file
	public String getDestination() { 
		String fileName = (String) JOptionPane.showInputDialog(new Frame("Input prompt"), 
				"Enter a file name for the edit file. This file will be saved in the same folder as the input file."); 
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
	
	public int wantToContinue(String input) {
		int j = JOptionPane.showConfirmDialog(new Frame("Input prompt"), input + " Do you want to continue?", "User Prompt", 
				JOptionPane.YES_NO_OPTION);
		return j;
	}
	
	public void initializer() throws IOException {
		try {																//input VarSifter file - file provided by user
			inputData = new BufferedReader(new FileReader(input));
			unsorted = new File(input.getParent() + "\\" + newFileName + "_unsorted.txt");
			writer = new PrintWriter(unsorted, "UTF-8");					
		} catch (FileNotFoundException error) {
			JOptionPane.showMessageDialog(new Frame("Error"), "The input PolyPhen file could not be found, possibly because it is open elsewhere. Please close and try again.");
			System.exit(0);
		}
		
		inputLine = inputData.readLine();
		inputLineSplit = inputLine.split("\t");
		
		for (int i = 0; i < inputLineSplit.length; i++) {
			if (inputLineSplit[i].contains("o_pos")) {
				protPosIndex = i;
			} else if (inputLineSplit[i].contains("o_aa1")) {
				refAAIndex = i;
			} else if (inputLineSplit[i].contains("o_aa2")) {
				altAAIndex = i;
			} else if (inputLineSplit[i].contains("prediction")) {
				predictionIndex = i;
			} else if (inputLineSplit[i].contains("pph2_prob")) {
				probIndex = i;
			} else if (inputLineSplit[i].contains("pph2_TPR")) {
				sensIndex = i;
			} else if (inputLineSplit[i].contains("pph2_FPR")) {
				specIndex = i;
			}
		}
		
		locIndex = inputLineSplit.length;		
		writer.println("#CHROM\tPOS\tID\tREF\tALT\tPPhen_GeneName\tPPhen_Transcript\tProtein_POS\tRef_AA\tAlt_AA\tPrediction\tProb\tSensitivity\tSpecificity");
	}
	
	public void runner() throws IOException {
		
		bRunner = System.currentTimeMillis();
		
		while (inputData.ready()) {
			inputLine = inputData.readLine();
			inputLineSplit = inputLine.split("\t");
			
			output = "";
			if (inputLine.substring(0, 1).equals("#")) {
				break;
			}
			
			output += locParser(inputLineSplit[locIndex]) + "\t";
			output += inputLineSplit[protPosIndex].trim() + "\t";
			output += inputLineSplit[refAAIndex].trim() + "\t";
			output += inputLineSplit[altAAIndex].trim() + "\t";
			output += inputLineSplit[predictionIndex].trim() + "\t";
			
			if (inputLineSplit[probIndex].contains("?") || inputLineSplit[sensIndex].contains("?") || inputLineSplit[specIndex].contains("?")) {
				output += "0" + "\t" + "0" + "\t" + "0";
			} else {
				output += fix(inputLineSplit[probIndex].trim()) + "\t";
				output += fix(inputLineSplit[sensIndex].trim()) + "\t";
				output += fix(Double.toString(1 - Double.parseDouble(inputLineSplit[specIndex])));
			}
			
			writer.println(output);
		}
		
		aRunner = System.currentTimeMillis();
		
		inputData.close();
		writer.close();
		
		PolySort pSort = new PolySort(unsorted, 0, 1, 3, 4, 11, input.getParent() + "\\" + newFileName + ".txt");
		
		timeMessage = "\nThe Parser Program Took " + (aRunner - bRunner) + "ms ("
				+ ((aRunner - bRunner) / 1000) + "." + (((aRunner - bRunner) % 1000) / 100) + " seconds)";
		
		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Completed." + timeMessage + "\n" + pSort.getTimeMessage());	
		System.gc();
		System.exit(0);
	}
	
	public String locParser(String pphenData) {
		String[] split = pphenData.split("\\|"); 
		
		String ref = split[1].substring(0,1);
		String var = split[1].substring(1);
		String pos = split[0].substring(split[0].indexOf(":") + 1);		
		String chrom = split[0].substring(split[0].indexOf("chr") + 3, split[0].indexOf(":"));
		
		String gene;
		String trans;
		
		if(split[split.length-2].startsWith("uc")) { //only for when transcript is uc
			gene = split[split.length - 1];
			trans = split[split.length - 2];
		} else { //Normal Format
			gene = split[split.length - 2];
			trans = split[split.length - 1];
		}
		return chrom + "\t" + pos + "\t-\t" + ref + "\t" + var + "\t" + gene + "\t" + trans;
	}
	
	//Truncator Method that truncates each number to 3 decimal places, removes scientific notation, and turns 0.0 into 0
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
	
	
	public static void main(String[] args) throws IOException {
		PolyPhenParser ppParser = new PolyPhenParser();
		ppParser.uI();
		ppParser.initializer();
		ppParser.runner(); 
	}
}