package general;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.filechooser.FileNameExtensionFilter;

public class BamChrChanger {
	
	//ArrayLists that contain the two different chr formats and effectively translates them from the old version to the new one if needed
	private static ArrayList<String> oldChrs = new ArrayList<String>();
	private static ArrayList<String> newChrs = new ArrayList<String>();
	private static String outputDir;
	
	//Initializes the chr lists by retrieving the bam chr config and reading the two inputs
	public static void initializeConfig(String configPath, String dir) throws IOException {
		

		File configz = new File(configPath + "\\" + "BAM_Chr_Config.txt");
		
		if(!configz.exists()) {
			
			configz = getConfigFile(configz.getAbsolutePath(), "BAM_Chr_Config.txt");
			
		}
		

		if(!(configz.isFile() && configz.canRead())) {
			
			JOptionPane.showMessageDialog(new Frame("Error"),
					"Error\nBam_Chr_Config config file is not a file/ cannot be read, system exiting.");
			System.exit(1);
			
		}
		

		BufferedReader configReder = new BufferedReader(new FileReader(configz));
		
		String Line;
		
		while((Line = configReder.readLine()) != null) {
			
			if (!Line.startsWith("#") && !Line.isEmpty() && !Line.trim().equals("") && !Line.trim().equals("\n")) {
				
				String[] configline = Line.split("=");
				
				oldChrs.add(configline[0]);
				newChrs.add(configline[1]);
				
				
			}
			
		}
		
		configReder.close();
		
		outputDir = dir;
		
	}
	

	//creates a FileChooser that allows the selection of an input marker config File
	public static File getConfigFile(String path, String name) {
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
				return getConfigFile(path, name);
		}
	}
	
	
	//Translates to the new chr format if the old chr is not used in the BAMs
	public static String translate(String chr) {
		
		return newChrs.get(oldChrs.indexOf(chr));
		
		
	}

	public static void writeHeaderErrors(ArrayList<String> bamchrs) throws IOException {
		
		File errors = new File(outputDir + "\\BAM_Header_Error.txt");
		
		PrintWriter errorWriter = new PrintWriter(errors);
		
		errorWriter.println("BAM Headers:");
		
		for(int i = 0; i < bamchrs.size(); i++) {
			
			errorWriter.println(bamchrs.get(i));
			
			
		}
		
		errorWriter.close();
	}
	
	
}
