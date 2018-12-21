package ethnicityMatcher;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URISyntaxException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;

import javax.print.DocFlavor.URL;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import general.*;

public class EthnicityCentralControl {
	
	private String bamDirectoryFile;
	private String[] ancMarkers;
	private File sampleFile;
	private String destination;
	private String cutoff;
	private HashMap<String, ArrayList<String>> bamChrMap;
	
	private File output;
	
	public EthnicityCentralControl(String bamDir, String[] ancMarkers, File sampleFile, String cutoff, String dest, HashMap<String, ArrayList<String>> bamChrMap) throws IOException {
		
		this.bamDirectoryFile = bamDir;
		this.ancMarkers = ancMarkers;
		this.sampleFile = sampleFile;
		this.cutoff = cutoff;
		this.destination = dest;
		this.bamChrMap = bamChrMap;
		
		
		ethnicityControl();
	}
	
	public void ethnicityControl() throws IOException {
		
				
				/// ***Put the ethnic config files and BAM file directory file under the same
				/// path as the .jar file
				PrintWriter report;
				PrintWriter report2;
				/// First define where all the markers are:
				// File f = new File(System.getProperty("java.class.path")); // finds the folder
				// that the JAR is stored in
//				File f = new File(System.getProperty("user.dir"));
//				System.out.println("fStr: " + f.toString());
//				File dir = f.getAbsoluteFile().getParentFile();
////				String path = dir.toString(); // the path to the folder
//				// System.out.println("path get name: " + dir.getName());
//				//
//				Path currentRelativePath = Paths.get("");
//				String path = currentRelativePath.toAbsolutePath().toString();
//				System.out.println("path: " + path);
//				System.out.println("Current relative path is: " + s);

				//
				// java.net.URL location =
				// EthnicityCentralControl.class.getProtectionDomain().getCodeSource().getLocation();
				// System.out.println(location.getFile());

				/// Retrieve the file with the ID's of the people (samples) to be analyzed for
				/// ethnicity
				//File sampleFile = TXTFile.getSampleFile();
				BufferedReader vsData = new BufferedReader(new FileReader(sampleFile));
				String Line = vsData.readLine();
				ArrayList<String> sample = new ArrayList<String>();
				while (Line != null) {
					sample.add(Line);
					Line = vsData.readLine();
				}
				vsData.close();
				/*
				/// Fill in the path for all the marker groups
				String configEAS = path + "/EASMarker";
				String configAMR = path + "/AMRMarker";
				String configEUR = path + "/EURMarker";
				String configAFR = path + "/AFRMarker";

				String configAMRneuEUR = path + "/PairWiseMarker/AMRneu_EUR.txt";
				String configAMRneuSAS = path + "/PairWiseMarker/AMRneu_SAS.txt";
				String configEURneuAMR = path + "/PairWiseMarker/EURneu_AMR.txt";
				String configEURneuSAS = path + "/PairWiseMarker/EURneu_SAS.txt";
				String configSASneuAMR = path + "/PairWiseMarker/SASneu_AMR.txt";
				String configSASneuEUR = path + "/PairWiseMarker/SASneu_EUR.txt";
				String configEURvAMR = path + "/PairWiseMarker/EUR_AMR.txt";
				String configSASvAMR = path + "/PairWiseMarker/SAS_AMR.txt";
				String configSASvEUR = path + "/PairWiseMarker/SAS_EUR.txt";
				*/
				/// Then define where all the BAM files are
				HashMap<String, String> BamMap = new HashMap<String, String>();
				/// Path to the BAM file directory file
				//String bamlo = path + "/EthnicityMatching_Bam.txt"; This is hard-coded BAM Dir File
				FindBam.InitializeBam(bamDirectoryFile, BamMap);
				
				/// Define the report output name
				String reportname = destination + ".txt";
				String reportname2 = reportname.substring(0, reportname.lastIndexOf(".")) + "_Ethnicity_Report.txt";
				
				String[] cutoffArr = cutoff.split(";");
				
				
				//String path2 = sampleFile.getPath().substring(0, sampleFile.getPath().indexOf(sampleFile.getName()));
				output = new File(reportname);
				report = new PrintWriter(output);
				report2 = new PrintWriter(new File(reportname2));
				
				// Request user input for cutoff criteria
				//String criteriacutoff = TXTFile.getCoverageCutoff();
				/*
				 * Moved to User Interface 
				/// Check if entered criteria cutoff fits requirement
				if (!criteriacutoff.contains(";") || !criteriacutoff.replaceAll("\\d", "").equals(";")) {
					JOptionPane.showMessageDialog(new JPanel(), "Cutoff input does not fit requirement. System exiting.");
					System.gc();
					System.exit(0);
				}

				String[] criteria = criteriacutoff.split(";");
				if (Integer.parseInt(criteria[0]) > 200 || Integer.parseInt(criteria[1]) > 500) {
					JOptionPane.showMessageDialog(new JPanel(), "Cutoff input does not fit requirement. System exiting.");
					System.gc();
					System.exit(0);
				}
				*/

				/// Put the marker sets address into a HashMap (name->config file path)
				HashMap<String, String> ConfigMap = new HashMap<String, String>();
				/*
				
				ConfigMap.put("configEAS", configEAS);
				ConfigMap.put("configAMR", configAMR);
				ConfigMap.put("configEUR", configEUR);
				ConfigMap.put("configAFR", configAFR);
				ConfigMap.put("configAMRneuEUR", configAMRneuEUR);
				ConfigMap.put("configAMRneuSAS", configAMRneuSAS);
				ConfigMap.put("configEURneuAMR", configEURneuAMR);
				ConfigMap.put("configEURneuSAS", configEURneuSAS);
				ConfigMap.put("configSASneuAMR", configSASneuAMR);
				ConfigMap.put("configSASneuEUR", configSASneuEUR);
				ConfigMap.put("configEURvAMR", configEURvAMR);
				ConfigMap.put("configSASvAMR", configSASvAMR);
				ConfigMap.put("configSASvEUR", configSASvEUR);
				
				*/
				
				ConfigMap.put("configEAS", ancMarkers[0]);
				ConfigMap.put("configAMR", ancMarkers[1]);
				ConfigMap.put("configEUR", ancMarkers[2]);
				ConfigMap.put("configAFR", ancMarkers[3]);
				ConfigMap.put("configAMRneuEUR", ancMarkers[4]);
				ConfigMap.put("configAMRneuSAS", ancMarkers[5]);
				ConfigMap.put("configEURneuAMR", ancMarkers[6]);
				ConfigMap.put("configEURneuSAS", ancMarkers[7]);
				ConfigMap.put("configSASneuAMR", ancMarkers[8]);
				ConfigMap.put("configSASneuEUR", ancMarkers[9]);
				ConfigMap.put("configEURvAMR", ancMarkers[10]);
				ConfigMap.put("configSASvAMR", ancMarkers[11]);
				ConfigMap.put("configSASvEUR", ancMarkers[12]);
				
				/// Make a new instance of the EthnicityMatching engine
				EthnicityFromBam ethnicityengine = new EthnicityFromBam();
				
				/// Print to the output ethnicity sample configuration file
				report.println("ID" + "\t" + "Ethnicity_matched");
				
				/// Print to the report output file.
				for (int i = 0; i < sample.size(); i += 1) {
					report2.println(sample.get(i));
					report2.println("MarkerGroup_set" + "\t" + "Marker_used" + "\t" + "MarkerScore");
					String ethnicity = ethnicityengine.EthPipeLine(FindBam.MakeBamString(sample.get(i), BamMap), sample.get(i),
							ConfigMap, Integer.parseInt(cutoffArr[0]), Integer.parseInt(cutoffArr[1]), report2, bamChrMap);
					report.println(sample.get(i) + "\t" + ethnicity);
					
				}
				report.close();
				report2.close();
				
//				/// The program is now complete
//				JOptionPane.showMessageDialog(new JPanel(),
//						"The program is now complete. Thank you for using the UDP bioinformatics service.");
//				System.gc();
//				System.exit(0);
				
				
	}
	
	public File getOutput() {
		
		return output;
		
	}
	
}
