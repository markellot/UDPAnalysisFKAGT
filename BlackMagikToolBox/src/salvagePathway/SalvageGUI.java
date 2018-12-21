package salvagePathway;
import general.*;
import java.awt.Frame;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import javax.swing.JOptionPane;

public class SalvageGUI {
	
	private File vsFile;
	private File vsPedigree;
	private File bamlo;
	private File ethlo;
	private int minor;
	private int yesPop;
	private int[] popHeaders;
	private File output;
	private String destination;
	private HashMap<String, ArrayList<String>> bamChrMap;
	
	public SalvageGUI(File vsFile, File vsPed, File bamDir, File ethConfig, int minor, int yesPop, int[] popHeadz, String dest, HashMap<String, ArrayList<String>> bamChrMap) throws IOException, InterruptedException{
		
		this.vsFile = vsFile;
		this.vsPedigree = vsPed;
		this.bamlo = bamDir;
		this.ethlo = ethConfig;
		this.minor = minor;
		this.yesPop = yesPop;
		this.popHeaders = popHeadz;
		this.destination = dest;
		this.bamChrMap = bamChrMap;
		
		FBGUI();
	}
	
	public void FBGUI() throws IOException, InterruptedException {
		
		// get the path of this program
		// File f = new File(System.getProperty("java.class.path")); // finds the folder
		// that the JAR is stored in
		// File dir = f.getAbsoluteFile().getParentFile();
		// String path = f.getPath();
		// String path = dir.toString(); // the path to the folder
		// String path = "U:/KaylaEmily/10_SalvagePathway/SalvagePathwaySoftware";
		// System.out.println("Path: " + path);
		// JOptionPane.showMessageDialog(new Frame("Input prompt"), "JAR path is " +
		// path);

		/// Request user input the VarSifter file from which the output VS is made
		
		//File vsFile = TXTFile.getVSFile();
		//String allgahlfile = vsFile.toString();

		String outputfilepath = vsFile.getPath().substring(0, vsFile.getPath().indexOf(vsFile.getName()));

		/// Request user input for the pedgiree.
		//File vsPedigree = TXTFile.getPedigreeFile();
		// the paths for all the output files
		String perfamilyoutput = destination + ".vs";
		String reportoutput = destination + "_Salvage_Report.txt";

		/// Read the pedigree
		BufferedReader pedData = new BufferedReader(new FileReader(vsPedigree));

		// Read in each pedigree line
		String Line = pedData.readLine();

		/// Input BAM file directory
		//File bamlo = TXTFile.getBAMDirectoryFile();
		String bamloPath = bamlo.getPath();
		HashMap<String, String> BamMap = new HashMap<String, String>();
		FindBam.InitializeBam(bamloPath, BamMap);
		
		/// Input the ethnicity config file to make the ethnicity HashMap
		//File ethlo = TXTFile.getEthConfigFile();
		String ethloPath = ethlo.getPath();
		HashMap<String, String> ethnicMap = new HashMap<String, String>();
		FindBam.InitializeBam(ethloPath, ethnicMap);

		/// Generate the output file.
		SalvagePipeline salvagepipe = new SalvagePipeline();
		/*
		while (Line != null) {
			String[] curLine = Line.split("\t");

			String[] bamlocation = FindBam.MakeBamString(curLine, BamMap);

			/// Retrieve parent ethnicities from config file.
			String[] parentethnic = FindBam.MakeEthString(new String[] { curLine[0], curLine[1] }, ethnicMap);
			
			salvagepipe.initializer(allgahlfile, curLine, bamlocation,
					perfamilyoutput + curLine[2] + "Perfamily" + ".vs", reportoutput + curLine[2] + "Perfamily.vs",
					parentethnic, pedData);

			Line = pedData.readLine();
		}
		
		*/
		
		String[] curLine = Line.split("\t");
		
		String[] bamlocation = FindBam.MakeBamString(curLine, BamMap);
		
		/// Retrieve parent ethnicities from config file.
		String[] parentethnic = FindBam.MakeEthString(new String[] { curLine[0], curLine[1] }, ethnicMap);
		
		salvagepipe.initializer(vsFile, curLine, bamlocation,
				perfamilyoutput, reportoutput,
				parentethnic, pedData, minor, yesPop, popHeaders, bamChrMap);
		
		//Line = pedData.readLine();
		
		pedData.close();
		
		output = salvagepipe.getOutput();
		
	}
	
	public File getOutput() {
		
		return output;
		
	}
}
