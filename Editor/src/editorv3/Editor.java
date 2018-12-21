package editorv3;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
//import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
//import java.io.FileWriter;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

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

public class Editor {

	// Object designed to hold lines - 0 for VS, 1 for edit
	private class LineHolder implements Comparable<LineHolder> {
		private String[] lineSplit;
		public String chrom;
		public int position;
		public String ref;
		public String var;
		
		
		//numFile can be anything if doing vsFiles
		public LineHolder(String inputLine, int type, int numFile) {
			lineSplit = inputLine.split("\t");
			if (lineSplit.length == 0) {
				JOptionPane
						.showMessageDialog(
								new Frame("Error"),
								"Your file contained a blank line. Please remove all blank lines and try again.");
				System.gc();
				System.exit(0);
			}
			
			if (type == 0) {
				chrom = lineSplit[chrIndex].substring(3);
				position = Integer.parseInt(lineSplit[posIndex]);
				ref = lineSplit[refIndex];
				var = lineSplit[varIndex];
			} else {
				chrom = lineSplit[editHeaderIndexes.get(numFile)[0]];
				position = Integer.parseInt(lineSplit[editHeaderIndexes.get(numFile)[1]]);
				ref = lineSplit[editHeaderIndexes.get(numFile)[2]];
				var = lineSplit[editHeaderIndexes.get(numFile)[3]];
			}
		}
		
		//numFile can be anything if doing vsFiles
		public LineHolder(String[] inputLine, int type, int numFile) {
			lineSplit = inputLine;
			if (type == 0) {
				chrom = lineSplit[chrIndex].substring(3);
				position = Integer.parseInt(lineSplit[posIndex]);
				ref = lineSplit[refIndex];
				var = lineSplit[varIndex];
			} else {
				chrom = lineSplit[editHeaderIndexes.get(numFile)[0]];
				position = Integer.parseInt(lineSplit[editHeaderIndexes.get(numFile)[1]]);
				ref = lineSplit[editHeaderIndexes.get(numFile)[2]];
				var = lineSplit[editHeaderIndexes.get(numFile)[3]];
			}
		}

		@Override
		public int compareTo(LineHolder line) {
			if (!chrom.equals(line.chrom)) {
				return chrom.compareTo(line.chrom);
			} else if (position == line.position) {
				return 0;
			} else if (position < line.position) {
				return -1;
			} else {
				return 1;
			}
		}
		
		public boolean isEqual(LineHolder line) {
			if (compareTo(line) != 0) {
				return false;
			} else if (!ref.equals(line.ref)) {
				return false;
			} else if (!var.equals(line.var)) {
				return false;
			} else {
				return true;
			}
		}

		@SuppressWarnings("unused")
		public String getLine() {
			String output = "";
			for (int i = 0; i < lineSplit.length-1; i++) {
				output = output.concat(lineSplit[i]).concat("\t");
			}
			output = output.concat(lineSplit[lineSplit.length-1]);
			return output;
		}
		
		
		public String[] getLineSplit() {
			return lineSplit;
		}
	}

	// Comparator for the Queue
	private class LineHolderComp implements Comparator<LineHolder> {
		@Override
		public int compare(LineHolder line1, LineHolder line2) {
			return line1.compareTo(line2);
		}
	}

	// The Queue of VSLines
	private PriorityQueue<LineHolder> vsQueue;

	// the Queue of Edit lines
	//private PriorityQueue<LineHolder> editQueue;
	private ArrayList<PriorityQueue<LineHolder>> editQueues = new ArrayList<PriorityQueue<LineHolder>>();

	// Input/Output
	private File vsFile;
	private String vsFilename;
	private String newFileName;
	//private File editFile;
	private BufferedReader vsReader;
	//private BufferedWriter writer;
	private PrintWriter writer;
	private ArrayList<BufferedReader> editFileReaders = new ArrayList<BufferedReader>();
	
	
	// Line Manipulation Variables
	private LineHolder vsLine;
	
	// VS file required Headers
	private int chrIndex;
	private int posIndex;
	private int refIndex;
	private int varIndex;
	private int mutTypeIndex;
	private int typeIndex;
	private int naIndex; // used to determine the index at which to add new
							// fields'

	// edit file required headers
	/*
	private int editChrIndex;
	private int editPosIndex;
	private int editRefIndex;
	private int editVarIndex;
	 */
	
	// ArrayLists of LineHolder objects
	private ArrayList<LineHolder> vsLines = new ArrayList<LineHolder>();
	//private ArrayList<LineHolder> editLines = new ArrayList<LineHolder>();
	private ArrayList<ArrayList<LineHolder>> editLinesList = new ArrayList<ArrayList<LineHolder>>();

	// ArrayList of headers that need to be added/edited
	//private ArrayList<int[]> editPairs; 
	private ArrayList<ArrayList<int[]>> editPairList = new ArrayList<ArrayList<int[]>>();// each entry has 2 ints: index 0 is the
													 // index in the vs file, index 1 is the
													 // index in the edit file
	private int newCols = 0; // counts the number of columns that have to be
								// added
	//private ArrayList<String> defaults = new ArrayList<String>();

	//private boolean editDone = false;
	private boolean[] editDoneList;
	
	//CADD input files for unscored variants
	private boolean CADDedit = false;
	private PrintWriter CADDwriter;
	private int CADDcounter = 0;
	private int CADDfileNum = 1;
	private String CADDpath;
	private int CADDLines;
	
	//PPhen input files for unscored variants
	private boolean PPHENedit = false; 
	private PrintWriter PPHENwriter;
	private int PPHENcounter = 0;
	private int PPHENfileNum = 1;
	private String PPHENpath;
	
	//Eigen Editor
	private boolean EigenEdit = false;
	private int eigenIndex;
	private int[] eigenRaws = new int[2];
	private int[] eigenPhreds = new int[2];
	private int rawIndex;
	private int phredIndex;
	
	
	//Multi MasterEdit files
	private String editFileDirec;
	private File editConfig; 
	private ArrayList<File> editFiles = new ArrayList<File>(); //list of edit files
	private ArrayList <ArrayList<String[]>> MEDefaults = new ArrayList<ArrayList<String[]>>(); //List of Default Values based on config file, First List is per file, second list is per line, array is for the split between the "=" sign (Header name [0], Default Value[1]
	private int numFiles = 0;
	private int caddIndex;
	private int polyIndex;
	private ArrayList<int[]> editHeaderIndexes = new ArrayList<int[]>(); //List of indexes for chr, pos, ref and var per ME file
	
	//time
	private long bRunner;
	private long aRunner;
	private String timeMessage;
	
	//Counting the lines
	private long lineCount = 0;
	
	
	// ///////////////////////////////// UI Methods // /////////////////////////////////////////////////
	
	public void uI() throws IOException {
		vsFile = getVSFile(); //get Varsifter File
		vsFilename = vsFile.getName(); //get Name of Varsifter File
		newFileName = getDestination(); //get Name of Output (does not need extension after since .txt will be added anyway
		editFileDirec = getEditDirec(); //gets directory of ME files
		editConfigAnalysis(); //analyzes config file for ME Files and Default Values
		MEChecker();
		if(CADDedit) {
			
			CADDLines = getCADDSize();
			
		}
	}
	
	//creates a FileChooser that allow the selection of  directory containing all the master edit files
	public String getEditDirec() {
		JFileChooser browseTo = new JFileChooser(); //creates a file chooser object
		browseTo.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY); //limits selectable objects to directories
		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please select the directory containing the " +
				"master edit files.");
		int returnVal = browseTo.showOpenDialog(new JPanel()); //opens a window for user to browse to directory
		
	    if(returnVal == JFileChooser.APPROVE_OPTION) {
	    	return browseTo.getSelectedFile().getPath(); //gets path of the selected directory
	    } else {
			int i = this.wantToContinue("No directory selected."); 
			if (i == 1) {
				System.exit(0);
			}
			return this.getEditDirec();
			
		}	
	}
	
	
	
	// gets the input VarSifterFile
	public File getVSFile() {
		JFileChooser browseTo = new JFileChooser(); // creates a file chooser
													// object
		FileNameExtensionFilter filter = new FileNameExtensionFilter(
				"VARSIFTER FILES", "vs", "txt");
		browseTo.setFileFilter(filter); // limits the viewable files to .vs and
										// .txt

		JOptionPane
				.showMessageDialog(
						new Frame("Input prompt"),
						"Please select an input VarSifter file. This file must have the 'VCF_Position',\n"
								+ "'VCF_Ref', and 'VCF_Var' fields from the Annotation utility as well \n"
								+ "as the 'Chr' field.");
		int returnVal = browseTo.showOpenDialog(new JPanel()); // opens a window
																// for user to
																// browse to vs
																// file to be
																// edited

		if (returnVal == JFileChooser.APPROVE_OPTION) {
			return browseTo.getSelectedFile(); // gets the path of the selected
												// file
		} else {
			int i = this.wantToContinue("No VarSifter file selected.");
			if (i == 1) {
				System.exit(0);
			}
			return this.getVSFile();
		}
	}
	
	// gets the name of the new file from the user
	public String getDestination() { // prompts the user to enter the name of a
										// destination for the edited file
		String fileName = (String) JOptionPane
				.showInputDialog(
						new Frame("Input prompt"),
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
	

	// gets the number of lines that will be sent to CADD per input file
	public int getCADDSize() {
		
		String size = (String) JOptionPane
				.showInputDialog(
						new Frame("Input prompt"),
						"Enter the number of lines for the input files you want to send to CADD");
		
		//Accounts for null, empty, not a number, decimal, greater than limit or less than or equal to 0
		if (size == null || size.equals("") || !size.matches("-?\\d+(\\.\\d+)?") || size.contains(".") || Integer.parseInt(size) > 999998 || Integer.parseInt(size) <= 0) { 
			int i = this.wantToContinue("You did not enter a valid input");
			if (i == 1) {
				System.exit(0);
			}
			return this.getCADDSize();
		}
		
		int caddLines = (Integer.parseInt(size) - 2); //Accounts for first two lines in each input file
		
		return caddLines;
	}
	
	// prompts user if they want to continue if they've hit cancel
	public int wantToContinue(String input) {
		int j = JOptionPane.showConfirmDialog(new Frame("Input prompt"), input
				+ " Do you want to continue?", "User Prompt",
				JOptionPane.YES_NO_OPTION);
		return j;
	}
	

	//creates a FileChooser that allows the selection of an input config File
	public File getConfigFile() {
		JFileChooser browseTo = new JFileChooser(); //creates a file chooser object
		
		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Config file not found. Please select the edit config file.");
		int returnVal = browseTo.showOpenDialog(new JPanel()); //opens a window for user to browse to vs file to be edited
		
		if(returnVal == JFileChooser.APPROVE_OPTION) {
			return browseTo.getSelectedFile(); //gets the path of the selected file
			
		} else {
			int i = this.wantToContinue("No Config file selected."); 
				if (i == 1) {
					System.exit(0);
				}
				return this.getConfigFile();
		}
	}
		
		// gets the input Master Edit File if the file was not found in the config file and ME Directory
		public File getMEFile(String path) {
			JFileChooser browseTo = new JFileChooser(); // creates a file chooser
														// object
			FileNameExtensionFilter filter = new FileNameExtensionFilter(
					"VARSIFTER FILES", "vs", "txt");
			browseTo.setFileFilter(filter); // limits the viewable files to .vs and
											// .txt

			JOptionPane
					.showMessageDialog(
							new Frame("Input prompt"), "Edit File with the path: " + path + " is not found. Please select the Master Edit file");
			int returnVal = browseTo.showOpenDialog(new JPanel()); // opens a window
																	// for user to
																	// browse to vs
																	// file to be
																	// edited

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				return browseTo.getSelectedFile(); // gets the path of the selected
													// file
			} else {
				int i = this.wantToContinue("No VarSifter file selected.");
				if (i == 1) {
					System.exit(0);
				}
				return this.getMEFile();
			}
		}
		
		
		//Receives and Analyzes config file for the ME files, Headers, and Default Values
		public void editConfigAnalysis() throws IOException {
			
			//gets config file
			editConfig = new File("edit_config.txt");
			if(!editConfig.exists()) {
				
				editConfig = getConfigFile();
				
			} else {
				
				int j = JOptionPane.showConfirmDialog(new Frame("Config Confirmer"), "Do you want to use this config?: " + editConfig.getAbsolutePath(), "Config Confirmer",
						JOptionPane.YES_NO_OPTION);
				
				if (j == 1) {
					
					editConfig = getConfigFile();
					
				}
			}
			
			BufferedReader configReder = new BufferedReader(new FileReader(editConfig));
			String line;
			
			//skips past all the initial lines til it occurs a "##"
			while((line = configReder.readLine()) != null) {
				
				if(line.trim().startsWith("##") && !line.trim().isEmpty()) {
					
					break;
					
					
				}
			}
			
			//if it never found a "##"
			if(line == null) {
				
				JOptionPane.showMessageDialog(new Frame("Error"), "Config file does not include a \"##\", the Program will now exit");
				System.exit(0);
				
			}
			
			//looping per line in the config file
			while(line != null) {
				
				//checks if there is a "##" for the next file
				if(line.trim().startsWith("##") && !line.trim().isEmpty()) {
					
					ArrayList<String[]> fileList = new ArrayList<String[]>();
					
					editFiles.add(new File(editFileDirec + "\\" + configReder.readLine()));
					
					if(!editFiles.get(numFiles).exists()) {
						
						editFiles.set(numFiles, getMEFile(editFiles.get(numFiles).getAbsolutePath()));
						
					}
					
					if(line.trim().startsWith("##CADD")) {
						

						//makes CADD file preparations
						CADDedit = true;
						caddIndex = numFiles;
						
						
						CADDpath = vsFile.getPath().substring(0,vsFile.getPath().indexOf(vsFile.getName())); 
						(new File(CADDpath + "\\" + vsFilename.substring(0, vsFilename.indexOf(".")) + "_CADD_Input_Files")).mkdirs();
						
						CADDwriter = new PrintWriter(CADDpath + "\\" + vsFilename.substring(0, vsFilename.indexOf(".")) 
								+ "_CADD_Input_Files\\" + vsFilename.substring(0, vsFilename.indexOf(".")) +  "_CADD_input_file_1.txt", "UTF-8");
						
						CADDwriter.println("##fileformat=VCFv4.0 - File 1 to be sent to CADD");
						CADDwriter.println("#CHROM\tPOS\tID\tREF\tALT");
						
						//Goes between each line searching for the default values
						
						
					} else if(line.trim().startsWith("##POLY")) { // Check if this is a poly file
						
						//makes Poly file preparations
						
						PPHENedit = true;
						polyIndex = numFiles;
						
						PPHENpath = vsFile.getPath().substring(0,vsFile.getPath().indexOf(vsFile.getName())); 
						(new File(PPHENpath + "\\" + vsFilename.substring(0, vsFilename.indexOf(".")) + "_PolyPhen_Input_Files")).mkdirs();
						
						PPHENwriter = new PrintWriter(PPHENpath + "\\" + vsFilename.substring(0, vsFilename.indexOf(".")) 
							+ "_PolyPhen_Input_Files\\" + vsFilename.substring(0, vsFilename.indexOf(".")) +  "_PolyPhen_input_file_1.txt", "UTF-8");
						
						
						
					} else if(line.trim().startsWith("##Eigen")) {
						
						EigenEdit = true;
						eigenIndex = numFiles;
					}
					

					while((line = configReder.readLine()) != null) {
						//Stop if at the next ME file
						if(line.startsWith("##")) {
							
							break;
							
						}
						//Potential Error Message if it does not contain a "=" sign
						if(!line.contains("=")) {
								
							JOptionPane.showMessageDialog(new Frame("Error"), "Lines do not include an \"=\" sign the program will now exit");
							System.exit(0);
						}
						
						//Splits each line and gets information
						String[] defaults = line.split("=");
						
						
						if(defaults.length == 1 && line.substring(line.lastIndexOf(line)).equals("=")) {
							String temp = defaults[0];
							
							String[] def = {temp, ""};
							
							fileList.add(def);
							
						}else {
							
							fileList.add(defaults);
							
						}
						
					}
					
						
					//Adds each list to the Default List per file
					MEDefaults.add(fileList);
					numFiles++; //<-- Very important for configuring each file	
					
				}
			}
			
			configReder.close();
			
			
			
		}
		
		
		public void MEChecker() {
			String message = "Do the Master Edit files shown below correspond to the correct version?\n\n";
			for(int i = 0; i < numFiles; i++) {
				
				message += editFiles.get(i).getAbsolutePath() + "\n\n";
				
			}
			
			int j = JOptionPane.showConfirmDialog(new Frame("Confirm Information"), message, "Confirm Information", JOptionPane.YES_NO_CANCEL_OPTION);
			
			if(j == 0) {
				
				return;
			}
			
			if(j==2) {
				
				int i = wantToContinue(""); 
				if (i == 1) {
					System.exit(0);
				}
				MEChecker();
				
			}
			
			boolean[] switchME = new boolean[numFiles];
			
			Object[] params = new Object[numFiles+1];
			
			
			for(int i = 0; i < numFiles; i++) {
				
				final int x = i;
				
				switchME[i] = false;
				
				JCheckBox check = new JCheckBox(editFiles.get(i).getAbsolutePath());
				
				check.addActionListener(new ActionListener() {
				
					public void actionPerformed(ActionEvent e) {
					
						switchME[x] = !switchME[x];
					
					}
				});
				
				params[i+1] = check;
				
			}				
			String massage = "Please select which Master Edit file you would like to switch out";
			
			params[0] = massage;
			
			int k = JOptionPane.showConfirmDialog(new Frame("Input Prompt"), params, "Change Master Edit file", JOptionPane.OK_CANCEL_OPTION);
				
			if(k == 2) {
				
				int i = wantToContinue(""); 
				if (i == 1) {
					System.exit(0);
				}
				MEChecker();
				
			}
			
			for(int i = 0; i < numFiles; i++) {
				
				if(switchME[i]) {
					
					showColumns(i);
					editFiles.set(i, getMEFile());
					
				}
				
				
			}
			
			
		}
		
		public void showColumns(int numFile) {
			
			JFrame jframe = new JFrame();
			GridBagLayout gbl = new GridBagLayout();
			GridBagConstraints gbc = new GridBagConstraints();
			jframe.setLayout(gbl);
			jframe.setTitle("Columns in Master Edit File");
			jframe.setSize((MEDefaults.get(numFile).size()/15) * 7 + 50, 400);
			jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			final JDialog diag = new JDialog(jframe, "Columns Revealer", true);
			
			JPanel jpanel = new JPanel();
			jpanel.setLayout(gbl);
			
			JLabel jtext = new JLabel("The following columns are shown for the MasterEdit file resolving to the following path:\n" + editFiles.get(numFile).getAbsolutePath(), SwingConstants.CENTER);
			gbc.gridx = 0;
			gbc.gridy = 0;
			
			
			gbc.fill = GridBagConstraints.HORIZONTAL;
			gbc.anchor = GridBagConstraints.WEST;
			gbc.insets = new Insets(5,5,5,5);
			gbc.weightx = 1.0;
			jpanel.add(jtext, gbc);
			
			int xPos = 0;
			int yPos = 2;
			
			gbc.gridx = xPos;
			gbc.gridy = yPos;
			gbc.fill = GridBagConstraints.HORIZONTAL;
			gbc.anchor = GridBagConstraints.WEST;
			gbc.insets = new Insets(5, 5, 5, 5);
			gbc.weightx = 1.0;
			jtext = new JLabel("#Chrom");
			jpanel.add(jtext, gbc);
			
			yPos++;

			gbc.gridx = xPos;
			gbc.gridy = yPos;
			gbc.fill = GridBagConstraints.HORIZONTAL;
			gbc.anchor = GridBagConstraints.WEST;
			gbc.insets = new Insets(5, 5, 5, 5);
			gbc.weightx = 1.0;
			jtext = new JLabel("Pos");
			jpanel.add(jtext, gbc);
			
			yPos++;

			gbc.gridx = xPos;
			gbc.gridy = yPos;
			gbc.fill = GridBagConstraints.HORIZONTAL;
			gbc.anchor = GridBagConstraints.WEST;
			gbc.insets = new Insets(5, 5, 5, 5);
			gbc.weightx = 1.0;
			jtext = new JLabel("Ref");
			jpanel.add(jtext, gbc);
			
			yPos++;

			gbc.gridx = xPos;
			gbc.gridy = yPos;
			gbc.fill = GridBagConstraints.HORIZONTAL;
			gbc.anchor = GridBagConstraints.WEST;
			gbc.insets = new Insets(5, 5, 5, 5);
			gbc.weightx = 1.0;
			jtext = new JLabel("Alt");
			jpanel.add(jtext, gbc);
			
			yPos++;
			
			
			for(int i = 0; i < MEDefaults.get(numFile).size(); i++) {
				
				String temp; 
				if (MEDefaults.get(numFile).get(i)[0].length() > 25) {
					temp = MEDefaults.get(numFile).get(i)[0].substring(0,16).concat("...").concat(MEDefaults.get(numFile).get(i)[0].substring(MEDefaults.get(numFile).get(i)[0].length()-10));
				} else {
					temp = MEDefaults.get(numFile).get(i)[0];
				}
				
				jtext = new JLabel(temp);
				gbc.gridx = xPos;
				gbc.gridy = yPos;
				gbc.fill = GridBagConstraints.HORIZONTAL;
				gbc.anchor = GridBagConstraints.WEST;
				gbc.insets = new Insets(5, 5, 5, 5);
				gbc.weightx = 1.0;
				jpanel.add(jtext, gbc);
				
				yPos++;
				if(yPos >= 20) {
					
					yPos = 1;
					xPos++;
					
				}
				
			}
			
			jtext = new JLabel("Please ensure the master edit you input next will have the above columns in the right order.", SwingConstants.CENTER);
			gbc.gridx = 0;
			gbc.gridy = 21;
			gbc.fill = GridBagConstraints.HORIZONTAL;
			gbc.anchor = GridBagConstraints.WEST;
			gbc.insets = new Insets(5,5,5,5);
			gbc.weightx = 1.0;
			jpanel.add(jtext, gbc);
			
			
			
			JButton okButton = new JButton("Ok");
			okButton.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					diag.dispose();
				}
			});
			
			gbc.gridx = 0;
			gbc.gridy = 22;
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
					MEChecker();
				}
			});
			
			gbc.gridx = 1;
			gbc.gridy = 22;
			gbc.fill = GridBagConstraints.HORIZONTAL;
			gbc.anchor = GridBagConstraints.WEST;
			gbc.insets = new Insets(5, 5, 5, 5);
			gbc.weightx = 1.0;
			jpanel.add(cancelButton, gbc);
			
			diag.getContentPane().add(jpanel);
			diag.pack();
			diag.setLocationRelativeTo(null);
			diag.setVisible(true);
			
		}
		

		// gets the input Master Edit File if the file was not found in the config file and ME Directory
		public File getMEFile() {
			JFileChooser browseTo = new JFileChooser(); // creates a file chooser
														// object
			FileNameExtensionFilter filter = new FileNameExtensionFilter(
					"VARSIFTER FILES", "vs", "txt");
			browseTo.setFileFilter(filter); // limits the viewable files to .vs and
											// .txt

			JOptionPane
					.showMessageDialog(
							new Frame("Input prompt"), "Please select the new Master Edit file");
			int returnVal = browseTo.showOpenDialog(new JPanel()); // opens a window
																	// for user to
																	// browse to vs
																	// file to be
																	// edited

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				return browseTo.getSelectedFile(); // gets the path of the selected
													// file
		} else {
			int i = this.wantToContinue("No VarSifter file selected.");
			if (i == 1) {
				System.exit(0);
			}
			return this.getMEFile();
		}
	}
		
		
	public void initializer() throws IOException {
		
		//File writingFile = new File(vsFile.getParent() + "\\" + newFileName);
		
		//writer = new BufferedWriter(new FileWriter(writingFile), 1000000000); 	//1 GB
		writer = new PrintWriter(vsFile.getParent() + "\\" + newFileName + ".vs", "UTF-8");
		
		try {																//input VarSifter file - file provided by user
			vsReader = new BufferedReader(new FileReader(vsFile), 1000000000); //VarSifter file must be presorted by chromosome then position (newLeftFlank)
		}
		catch (FileNotFoundException error) {
			JOptionPane.showMessageDialog(new Frame("Error"), "The VarSifter file could not be found, possibly because it is open elsewhere. Please close and try again.");
			writer.close();
			System.exit(0);
		}
		
		try {
			for(int i = 0; i < numFiles; i++) {   //adds all edit readers needed for reading all edit files
				
				editFileReaders.add(new BufferedReader(new FileReader(editFiles.get(i))));
				
			}
		}
		catch (FileNotFoundException error) {
			JOptionPane.showMessageDialog(new Frame("Error"), "The edit file could not be found, possibly because it is open elsewhere. Please close and try again.");			
			writer.close();
			vsReader.close();
			System.exit(0);
		}
		
		List<String> headers = java.util.Arrays.asList(vsReader.readLine().split("\t")); //finding headers from the VS File
		
		getVSHeaders(headers); //Gets relevant index info from headers
		
		
		
		//gets Headers
		String [] vsHeaders = new String[naIndex - 1 + newCols]; //Has length of naIndex(the number of indexes until the population data is found, then the newColumns
		ArrayList<String[]> editHeadz = new ArrayList<String []>();
		for(int i = 0; i < numFiles; i++) {
			String temp = editFileReaders.get(i).readLine();	//finding headers from the edit file
			while (temp.substring(0,1).equals("#") && (temp.substring(1,2).equals("#"))) { //moves through lines until an appropriately formatted header line is found
				temp = editFileReaders.get(i).readLine();
			}
			List<String> editHeaders = java.util.Arrays.asList(temp.split("\t"));
			getEditHeaders(editHeaders); //Gets relevant index information from headers
			
			editColumnDetector(headers, editHeaders, i); //finding the fields that need to be edited
			
			vsQueue = new PriorityQueue<LineHolder>(1000, new LineHolderComp()); //creates the Queue that will manage the VS lines
			
			
			//Creates editQueues that will manage the editFiles
			editQueues.add(new PriorityQueue<LineHolder>(5, new LineHolderComp()));
			editQueues.get(i).add(new LineHolder(editFileReaders.get(i).readLine(), 1, i));
			
			String[] editTemp = editHeaders.toArray(new String[editHeaders.size()]);
			
				
			String[] vsTemp = headers.toArray(new String[headers.size()]);
			vsHeaders = addColumns(vsTemp); //editing Header Array until it has the appropriate length
			
			editHeadz.add(editTemp); //adds edit file Headers for each edit file
			
		}
		
		for(int i = 0; i < numFiles; i++) {
			
			vsHeaders = editIn(vsHeaders, editHeadz.get(i), i); //edits in the information for each edit file to the header array
			
		}
		
		vsWriter(vsHeaders); //writes headers
		
		//initializes the done list which references if each edit file is finished
		editDoneList = new boolean[numFiles];
		for(int i = 0; i < numFiles; i++) {
			
			editDoneList[i] = false;
			
		}
		
		if(EigenEdit) {
			
			eigenRaws[0] = editPairList.get(eigenIndex).get(0)[0];
			eigenPhreds[0] = editPairList.get(eigenIndex).get(1)[0];
			eigenRaws[1] = editPairList.get(eigenIndex).get(2)[0];
			eigenPhreds[1] = editPairList.get(eigenIndex).get(3)[0];
			
		} else if(headers.indexOf("Eigen") != -1) {
			
			EigenEdit = true;
			
			eigenRaws[0] = headers.indexOf("Eigen");
			
			eigenPhreds[0] = headers.indexOf("EPHRED");
			
			eigenRaws[1] = headers.indexOf("EigenPC");
			
			eigenPhreds[1] = headers.indexOf("EPCPHRED");
			
		}
		
		if(EigenEdit && (!CADDedit && !(headers.indexOf("RawScore") != -1))) {
			
			EigenEdit = false;
			
		} else {
			
			if(CADDedit) {
				
				rawIndex = editPairList.get(caddIndex).get(0)[0];
				phredIndex = editPairList.get(caddIndex).get(1)[0];
				
				
			} else if(headers.indexOf("RawScore") != -1) {
				
				rawIndex = headers.indexOf("RawScore");
				phredIndex = headers.indexOf("PHRED");
				
				
			} else {
				
				System.out.println("You faked up...");
				
			}
			
		}
		
		
	}
		/*
	// gets the default values for new columns
	public String getDefault(String field) {
		String defaultStr = (String) JOptionPane.showInputDialog(new Frame(
				"Input prompt"), "Please enter a default value for the field '"
				+ field + "'. Default values cannot contain tabs.");
		if (defaultStr == null) {
			return "";
		} else if (defaultStr.equals("")) {
			return defaultStr;
		}
		return defaultStr;
	}

	 */
		
	// gets headers from the input VS line
	public void getVSHeaders(List<String> headers) {
		chrIndex = headers.indexOf("Chr"); // finds indices of the five headers
											// to the left
		posIndex = headers.indexOf("VCF_Position"); // if header does not exist,
													// String.indexOf returns
		refIndex = headers.indexOf("VCF_Ref"); // a -1
		varIndex = headers.indexOf("VCF_Var");
		mutTypeIndex = headers.indexOf("muttype");
		typeIndex = headers.indexOf("type");
		
		//This snippet finds the naindex, the column numnber in the headders list where the metadata columns ends and sample data columns begin. 
		// The identifier is that the last three charaters in the headder are ".NA"
		// A tricky issue is when the column header is less that three characters total! simply searching for ".NA" will hang the program!
		
		for (int i = 5; i < headers.size(); i++) {	
			int val = headers.get(i).length();
			// here is where you have to account for headders with less than 3 characters
			Pattern r = Pattern.compile("\\.NA$");
			Matcher m = r.matcher(headers.get(i));
			
			if(val > 3){
				if (m.find()) {
					naIndex = i;
					break;
				}
				
			}
		}

		boolean indexError = false;
		String errorMessage = "";
		if (chrIndex == -1) {
			errorMessage += "The header 'Chr' appears no where in the header line of the VarSifter file.\n"; // missing 'Chr' header
			indexError = true;
		}
		if (posIndex == -1) {
			errorMessage += "The header 'VCF_Position' appears no where in the header line of the VarSifter file.\n"; // missing 'LeftFlank' header
			indexError = true;
		}
		if (refIndex == -1) {
			errorMessage += "The header 'VCF_Ref' appears no where in the header line of the VarSifter file.\n"; // missing 'ref_allele' header
			indexError = true;
		}
		if (varIndex == -1) {
			errorMessage += "The header 'VCF_Var' appears no where in the header line of the VarSifter file.\n"; // missing 'var_allele' header
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

	// gets headers from the input edit line
	public void getEditHeaders(List<String> headers) {
		
		int[] thisEditIndex = new int[4]; //Chr is 0, Pos is 1, Ref is 2, Var is 3
		
		thisEditIndex[0] = headers.indexOf("#CHROM"); // finds indices of the five
													// headers to the left
		thisEditIndex[1] = headers.indexOf("POS"); // if header does not exist,
												// String.indexOf returns
		thisEditIndex[2] = headers.indexOf("REF"); // a -1
		
		thisEditIndex[3] = headers.indexOf("ALT");

		boolean indexError = false;
		String errorMessage = "";
		if (thisEditIndex[0] == -1) {
			errorMessage += "The header '#CHROM' appears no where in the header line of the edit file.\n"; // missing 'Chr' header
			indexError = true;
		}
		if (thisEditIndex[1] == -1) {
			errorMessage += "The header 'POS' appears no where in the header line of the edit file.\n"; // missing 'LeftFlank' header
			indexError = true;
		}
		if (thisEditIndex[2] == -1) {
			errorMessage += "The header 'REF' appears no where in the header line of the edit file.\n"; // missing 'ref_allele' header
			indexError = true;
		}
		if (thisEditIndex[3] == -1) {
			errorMessage += "The header 'VAR' appears no where in the header line of the edit file.\n"; // missing 'var_allele' header
			indexError = true;
		}
		if (indexError) { // if any of the above if statements are true,
			errorMessage += "As stated above. One or more of the necessary headers do not appear" 	// then the header does not appear within the
					+ " in your edit file. \nPlease ensure these headers are in " 				  	// header line. Thus, we print the missing header(s)
					+ "your edit file prior to running this program."; 	// along with the statement below and the program exits.
			JOptionPane.showMessageDialog(new Frame("Error Message"), errorMessage);
			System.exit(0);
		}
		
		editHeaderIndexes.add(thisEditIndex);
		
	}

	// detects the new columns of the edit file, tells the user that they're
	// going to be added
	public void editColumnDetector(List<String> headers,
			List<String> editHeaders, int numFile) throws IOException {
		int rightBound = Math.max(editHeaderIndexes.get(numFile)[0],
				Math.max(editHeaderIndexes.get(numFile)[1], Math.max(editHeaderIndexes.get(numFile)[2], editHeaderIndexes.get(numFile)[3]))); //finds the right bound most index 	
		ArrayList<int[]> editPairs = new ArrayList<int[]>();
		
		int[] temp;
		int deletions = 0;
		for (int i = rightBound + 1; i < editHeaders.size(); i++) {  //loops through the edit header line and checks if the original vs header line already has the header in its list
			if (!headers.contains(editHeaders.get(i))) {
				newCols++;
				temp = new int[2];
				temp[0] = naIndex - 1 + newCols; //The new header index for the location of where it will be in the VS file
				temp[1] = i; //the location of that same respective header located in the edit file
				editPairs.add(temp);
				
				//Does the Default Name found in the config match with the header name?
				if(!(MEDefaults.get(numFile).get(i-deletions-rightBound-1)[0].equals(editHeaders.get(i)))) {
					
					
					
					JOptionPane.showMessageDialog(new Frame("Error"), "The edit file has conflicting headers with the config file. Please edit the config file and try again.\n"
							+ "The Header in the Config File is: " + MEDefaults.get(numFile).get(i-deletions-rightBound)[0] + " and the Header in the corresponding Edit File is: " + editHeaders.get(i));
					writer.close();
						
					for(int k = 0; k < numFiles; k++) {
						
						editFileReaders.get(k).close();
						
					}
					vsReader.close();
					System.exit(0);
					
					
				}
				
			} else {
				temp = new int[2];
				temp[0] = headers.indexOf(editHeaders.get(i)); //The new header index for the location of where it will be in the VS file
				temp[1] = i; //the location of that same respective header located in the edit file
				editPairs.add(temp);
				
				//Does the Default Name found in the config match with the header name?
				if(!(MEDefaults.get(numFile).get(i-deletions-rightBound-1)[0].equals(editHeaders.get(i)))) {
					

					
					JOptionPane.showMessageDialog(new Frame("Error"), "The edit file has conflicting headers with the config file. Please edit the config file and try again.\n"
							+ "The Header in the Config File is: " + MEDefaults.get(numFile).get(i-deletions-rightBound)[0] + " and the Header in the corresponding Edit File is: " + editHeaders.get(i));
					writer.close();
					
					for(int k = 0; k < numFiles; k++) {
						
						editFileReaders.get(k).close();
						
					}
					vsReader.close();
					System.exit(0);
					
					
				}
				
				//removes the default line from the config because there is no need to add a new column for an already existing column
				MEDefaults.get(numFile).remove(i-deletions - rightBound - 1);
				deletions++;
				
			}
		}
		
		//Error Message
		if (editPairs.size() == 0) {
			JOptionPane
					.showMessageDialog(new Frame("Input prompt"),
							"There were no fields to edit or add. The program will now close.");
			System.exit(0);
		}
		
		editPairList.add(editPairs);
		
	}

	//creates the cancel window
	public JFrame canceler() {
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("Editor");
		jframe.setSize(500, 100);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);
		
		JLabel jtext = new JLabel("Editor Program running. If you'd like to abort, hit the cancel button below.", SwingConstants.CENTER);
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
					vsReader.close();
				} catch (IOException e1) {
				}
				
				writer.close();
				
				if (CADDwriter != null) {
					CADDwriter.close();
					try {
						gzipFile(CADDpath + "\\" + vsFilename.substring(0, vsFilename.indexOf(".")) + "_CADD_Input_Files\\" + vsFilename.substring(0, vsFilename.indexOf(".")) 
						+ "_CADD_input_file_" + Integer.toString(CADDfileNum) + ".txt", CADDpath + "\\" + vsFilename.substring(0, vsFilename.indexOf(".")) + "_CADD_Input_Files\\" + vsFilename.substring(0, vsFilename.indexOf(".")) 
						+ "_CADD_input_file_" + Integer.toString(CADDfileNum) + ".gz");
					} catch (IOException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}
				}
				
				try {
					
					for(int i = 0; i < numFiles; i++) {
						
						editFileReaders.get(i).close();
						
					}
				} catch (IOException e1) {
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
	
	public void runner() throws IOException {
		JFrame jframe = canceler();
		jframe.setLocationRelativeTo(null);
		jframe.setVisible(true);
		
		//Before Running 
		bRunner = System.currentTimeMillis();
		
		while (vsReader.ready()) {
			while (vsQueue.size() < 1000) { // build up the Queue
				if (!vsReader.ready()) { // if we've reached the end of the
											// file, break
					break;
				}
				vsRead(); //otherwise, read
			}

			while (vsQueue.size() > 500) { // once the queue is at capacity,
											// process it down
				merger();
			}
		}
		
		while (vsQueue.size() > 0) { // still going to be have of a queue to process (fencepost)
			merger();
		}
		
		//After Running
		aRunner = System.currentTimeMillis();
		
		//close all readers/writers
		jframe.dispose();
		
		
		
		//Final Message Screen
		timeMessage = ("The Editor Program Took " + (aRunner - bRunner) + "ms ("
				+ ((aRunner - bRunner) / 1000) + "." + (((aRunner - bRunner) % 1000) / 100) + " seconds");
		
		
		
		JFrame finalFrame = new JFrame("Completed");
		
		finalFrame.setUndecorated( true );
		finalFrame.setVisible( true );
		finalFrame.setLocationRelativeTo( null );
		finalFrame.setAlwaysOnTop( true );
		
		
		JOptionPane.showMessageDialog(finalFrame, ("Completed.\n\n") + timeMessage);
		
		
		//Closes everything
		vsReader.close();
		writer.close();
		

		for(int i = 0; i < numFiles; i++) {
			
			editFileReaders.get(i).close();
			
		}
		
		if (CADDedit) {
			CADDwriter.close();
			gzipFile(CADDpath + "\\" + vsFilename.substring(0, vsFilename.indexOf(".")) + "_CADD_Input_Files\\" + vsFilename.substring(0, vsFilename.indexOf(".")) 
			+ "_CADD_input_file_" + Integer.toString(CADDfileNum) + ".txt", CADDpath + "\\" + vsFilename.substring(0, vsFilename.indexOf(".")) + "_CADD_Input_Files\\" + vsFilename.substring(0, vsFilename.indexOf(".")) 
			+ "_CADD_input_file_" + Integer.toString(CADDfileNum) + ".txt.gz");
		}
		if (PPHENedit) {
			PPHENwriter.close();
		}
		
		System.exit(0);
	}
	
	//Merges the vsfile with the edit files
	public void merger() throws IOException {
		vsLine = vsQueue.poll(); //pull the first line from the Queue of VS lines
		vsLines.add(vsLine); //add to the ArrayList of varSifter lines
		
		//if there are still things in the queue, we want all lines with the same coordinates - peek into the queue to check coordinates, add to the array list if they match what is currently in the list
		while ((!vsQueue.isEmpty()) && (vsLine.compareTo(vsQueue.peek()) == 0)) {
			vsLines.add(vsQueue.poll()); //adds to arraylist and pulls out of queue simultaneously
		}
		
		String VCFoutput; 
		String PPHENoutput;
		
		//Checks to see if all the edit files are done reading or if the edit files are all ahead of the vsFile
		if(isEditPast()) {
			
			//if so, goes through vsLines and print all the Cadd and Polyphen lines to their respetive output file
			for(LineHolder l : vsLines) {
				

				if(l.getLineSplit()[mutTypeIndex].equals("INDEL")) {
					
					l = new LineHolder(eigenEdit(l.getLineSplit()), 0, -1);
					
				}
				
				vsWriter(l.getLineSplit());
				if(CADDedit && (!(l.chrom.contains("U") || l.chrom.contains("g") || l.chrom.contains("M") || l.chrom.contains("_") || l.ref.equals("N")))) {
					VCFoutput = "";
					VCFoutput += l.chrom + "\t" + l.position + "\t.\t" + l.ref + "\t" + l.var;
					CADDwriter.println(VCFoutput);
					CADDchunk();
					
					
				} else if ((PPHENedit) && (l.getLineSplit()[mutTypeIndex].equals("SNP")) && ((l.getLineSplit()[typeIndex].contains("NON_SYNONYMOUS")) || (l.getLineSplit()[typeIndex].contains("nonsynonymous")) || (l.getLineSplit()[typeIndex].contains("missense"))) && (!l.getLineSplit()[chrIndex].contains("g"))) {
					PPHENoutput = "";
					PPHENoutput += "chr" + l.chrom + ":" + l.position + " " + l.ref + "/" + l.var;
					PPHENwriter.println(PPHENoutput);
					PPHENchunk();
				}
				
				
			}
			
			vsLines.clear();
			return;
			
		}
		
		//loops through each edit file checking if they contain the information for each vsLine and adds them to the editLinesList
		for(int i = 0; i < numFiles; i++) {	
			
			while ((!editDoneList[i]) && (editQueues.get(i).peek().compareTo(vsLines.get(0)) < 0)) { // while the edit queue is lagging behind the current point of the vsFile, get rid of the edit lines
				editRead(i); //skips over edit lines
			}
			
			ArrayList<LineHolder> editLine = new ArrayList<LineHolder>();
			while ((!editDoneList[i]) && (editQueues.get(i).peek().compareTo(vsLines.get(0)) == 0)) {
				
				editLine.add(editRead(i));  // pulls the lines at the front of the editqueue that match the varsifter arraylist into the edit array list
				
			}
			
			editLinesList.add(editLine);
		}
		
		mergeHelper(); //will process the vs and edit line array lists
		
	}

	public void mergeHelper() throws IOException {
		
		Boolean[][] matched = new Boolean[vsLines.size()][numFiles]; //a boolean array showing which of the vslines have been paired with an edit line - default false
		for (int k = 0; k < matched.length; k++) {  //k = number vsLine and i = numFile
			for(int i = 0; i < matched[0].length; i++) {
				matched[k][i] = false; //sets defaults
			}
		}
		
		//Loops through each file and checks each edit
		for(int i = 0 ; i < numFiles; i++) {
			
			for(LineHolder j : editLinesList.get(i)) {
				int temp = 0;
				for(LineHolder k : vsLines) {
					
					if(j.isEqual(k)) {
						
						k = new LineHolder(editIn(k.getLineSplit(), j.getLineSplit(), i), 0, -1);
						matched[temp][i] = true;
						
					}
					
					temp++;
				}
				
				
			}
			
			
		}
		
		String VCFoutput;
		String PPHENoutput;
		int temp2 = 0;
		for (LineHolder j : vsLines) { //once everything is processed, print out each line in the varsifter array list
			
			if(j.getLineSplit()[mutTypeIndex].equals("INDEL")) {
				
				j = new LineHolder(eigenEdit(j.getLineSplit()), 0, -1);
				
			}
			
			vsWriter(j.getLineSplit()); //prints the lines
			
			//Loops through each file looking to see if any of the vslines did not match with the cadd or poly edit files
			for(int i = 0; i < numFiles; i++) {
				
				if (!matched[temp2][i]) { //for lines that are not matched
					if ((CADDedit) && (i == caddIndex) && (!(j.chrom.contains("U") || j.chrom.contains("g") || j.chrom.contains("M") || j.chrom.contains("_") || j.ref.equals("N")))) { //print CADD if appropriate
						VCFoutput = "";
						VCFoutput += j.chrom + "\t" + j.position + "\t.\t" + j.ref + "\t" + j.var; 
						CADDwriter.println(VCFoutput);
						CADDchunk();
					} else if ((PPHENedit) && (i == polyIndex) && (j.getLineSplit()[mutTypeIndex].equals("SNP")) && ((j.getLineSplit()[typeIndex].contains("NON_SYNONYMOUS")) || (j.getLineSplit()[typeIndex].contains("nonsynonymous")) || (j.getLineSplit()[typeIndex].contains("missense")) && (!j.getLineSplit()[chrIndex].contains("g")))) { //print PPHEN if appropriate
						PPHENoutput = "";
						PPHENoutput += "chr" + j.chrom + ":" + j.position + " " + j.ref + "/" + j.var;
						PPHENwriter.println(PPHENoutput);
						PPHENchunk();
					}
				}
			}
			
			temp2++;
		}
		
		editLinesList.clear();
		vsLines.clear();
	}
	
	
	//Checks to see if all the edit files are done reading or if the edit files are all ahead of the vsFile
	public boolean isEditPast() {
		
		boolean done = true;
		boolean past = true;
		
		for(int i = 0; i < numFiles; i++) {
			
			if(!editDoneList[i]) {
				
				done = false;
				
				if(!(vsLines.get(0).compareTo(editQueues.get(i).peek()) < 0)) {
					
					past = false;
					
					
				}
				
			} 
		
		}
		
		return (done || past);
		
	}
	
	//Makes new Cadd files if the counter is greater than or equal to 16383 lines
	public void CADDchunk() throws IOException {
		CADDcounter++;
		
		
		if (CADDcounter >= CADDLines) {
			CADDwriter.close();
			
			gzipFile(CADDpath + "\\" + vsFilename.substring(0, vsFilename.indexOf(".")) + "_CADD_Input_Files\\" + vsFilename.substring(0, vsFilename.indexOf(".")) 
			+ "_CADD_input_file_" + Integer.toString(CADDfileNum) + ".txt", CADDpath + "\\" + vsFilename.substring(0, vsFilename.indexOf(".")) + "_CADD_Input_Files\\" + vsFilename.substring(0, vsFilename.indexOf(".")) 
			+ "_CADD_input_file_" + Integer.toString(CADDfileNum) + ".txt.gz");
			
			CADDfileNum++;
			CADDcounter = 0;
			//String tempFilename = vsFile.getName().substring(0, vsFile.getName().indexOf("."));
			CADDwriter = new PrintWriter(CADDpath + "\\" + vsFilename.substring(0, vsFilename.indexOf(".")) + "_CADD_Input_Files\\" + vsFilename.substring(0, vsFilename.indexOf(".")) 
					+ "_CADD_input_file_" + Integer.toString(CADDfileNum) + ".txt", "UTF-8");
			//CADDwriter.println("##fileformat=VCFv4.1 _ File " + Integer.toString(CADDfileNum) + " to be sent to CADD");
			CADDwriter.println("##fileformat=VCFv4.1");
			CADDwriter.println("#CHROM\tPOS\tID\tREF\tALT");
		}
	}
	
	public static void gzipFile(String from, String to) throws IOException {
	    // Create stream to read from the from file
	    FileInputStream in = new FileInputStream(from);
	    // Create stream to compress data and write it to the to file.
	    GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(to));
	    // Copy bytes from one stream to the other
	    byte[] buffer = new byte[4096];
	    int bytes_read;
	    while ((bytes_read = in.read(buffer)) != -1)
	      out.write(buffer, 0, bytes_read);
	    // And close the streams
	    in.close();
	    out.close();
	}
	
	public String[] eigenEdit(String[] curlineSplit) {
		
		curlineSplit[eigenRaws[0]] = curlineSplit[rawIndex];
		curlineSplit[eigenRaws[1]] = curlineSplit[rawIndex];
		curlineSplit[eigenPhreds[0]] = curlineSplit[phredIndex];
		curlineSplit[eigenPhreds[1]] = curlineSplit[phredIndex];
		
		return curlineSplit;
		
	}
	
	
	//Makes new PolyPhen files if the counter is greater than 5000000 lines
	public void PPHENchunk() throws FileNotFoundException, UnsupportedEncodingException {
		PPHENcounter++;
		if (PPHENcounter >= 5000000) {
			PPHENwriter.close();
			PPHENfileNum++;
			PPHENcounter = 0;
			//String tempFilename = vsFile.getName().substring(0, vsFile.getName().indexOf("."));
			PPHENwriter = new PrintWriter(PPHENpath + "\\" + vsFilename.substring(0, vsFilename.indexOf(".")) + "_PolyPhen_Input_Files\\" + vsFilename.substring(0, vsFilename.indexOf(".")) 
					+ "_PolyPhen_input_file_" + Integer.toString(PPHENfileNum) + ".txt", "UTF-8");
		}
	}
	
	//Reads from the Varsifter file and adds columns to it
	public void vsRead() throws IOException {
		String line = vsReader.readLine();
		String[] lineSplit = line.split("\t", naIndex+5); //splits the line and leaves the sample data after the 5th index hanging in the last element
		lineSplit = addColumns(lineSplit);
		vsQueue.add(new LineHolder(lineSplit, 0, -1));
		
		lineCount++;
		if(lineCount % 10000 == 0) {
			
			System.out.println(lineCount);
			
		}
	}

	// pops off the top of the Edit Queue and reads in a new line behind it if
	// there is another line to be read
	//Takes parameter numFile to divide between edit readers and lists
	public LineHolder editRead(int numFile) throws IOException {
		if (editFileReaders.get(numFile).ready()) {
			editQueues.get(numFile).add(new LineHolder(editFileReaders.get(numFile).readLine(), 1, numFile));
		} else {
			editDoneList[numFile] = true;
		}
		
		return editQueues.get(numFile).poll();
	}
	

	// adds the necessary blank columns to the curLineSplit, which is given as
	// input
	public String[] addColumns(String[] curLineSplit) {
		String curLine = "";
		for (int i = 0; i < naIndex; i++) {
			curLine = curLine.concat(curLineSplit[i]).concat("\t");
		}
		
		for (int j = 0; j < numFiles; j++) {
			for(int i = 0; i < MEDefaults.get(j).size(); i++) {
				curLine = curLine.concat(MEDefaults.get(j).get(i)[1]).concat("\t");
			}
		}

		for (int k = naIndex; k < curLineSplit.length; k++) {
			curLine = curLine.concat(curLineSplit[k]).concat("\t");
		}
		curLineSplit = curLine.trim().split("\t", naIndex+5+newCols);
		return curLineSplit;
	}
	
	
	
	// edits the fields based on the input vs and edit lines
	public String[] editIn(String[] curLineSplit, String[] editLineSplit, int numFile) {
		
		//edits the original varsifter line with the line in the edit file using the editPairList to reference the locations for each data point
		for(int j = 0; j < editPairList.get(numFile).size(); j++) {
			
			curLineSplit[editPairList.get(numFile).get(j)[0]] = editLineSplit[editPairList.get(numFile).get(j)[1]];
			
		}
		
		
		return curLineSplit;
	}
	
	
	// takes a String[] of a line and concats, then prints
	public void vsWriter(String[] lineSplit) throws IOException {
		String output = lineSplit[0]; //gets the first item in the tab delimited array
		
		for (int i = 1; i < lineSplit.length; i++) { //concatenates the rest of the line
			output += "\t" + lineSplit[i];
		}
		output = output.trim();

		if (output.length() != 0) {
			writer.println(output); //writes the line; trim should remove the extra white space
		}
	}

	public static void main(String[] args) throws IOException {
		Editor derp = new Editor();
		derp.uI();
		derp.initializer();
		derp.runner();
	}

}