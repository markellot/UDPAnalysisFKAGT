import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.*;
/*
 * Sorry for improper variable name conventions ;)
 * 
 * -Faris
 */

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

public class SnpEffEr {
	
	private File vsFile;
	private BufferedReader reder;
	private File output;
	private PrintWriter righter;
	
	private int chrIndex;
	private int refIndex;
	private int varIndex;
	private int mutTypeIndex;
	private int naPos;
	private int naDistanceApart;
	private int columns;
	private int acIndex;
	private int afIndex;
	private int totIndex;
	private int leftFlankIndex;
	
	private ArrayList<Integer> toDelete;
	private boolean columnDeleter;
	
	private String line;
	private String[] lineArr;
	private List<String> headz = new ArrayList<String>();
	
	private File headerConfig;
	private File glConfig;
	
	private List<String> oldGl = new ArrayList<String>();
	private List<String> newGl = new ArrayList<String>();
	
	//Parallelize
	private File dir;
	private int sepFiles = 0;
	private List<File> tempIn = new ArrayList<File>();
	private List<File> tempOut = new ArrayList<File>();
	private List<BufferedReader> readIn = new ArrayList<BufferedReader>();
	private List<PrintWriter> readOut = new ArrayList<PrintWriter>();
	private ExecutorService executor;
	private final int NTHREDS = Runtime.getRuntime().availableProcessors();
	
	
	private long bRunner;
	private long aRunner;
	private String timeMessage;
	
	public SnpEffEr(File vsFile, boolean columnDeleter) throws IOException {
		
		this.vsFile = vsFile;
		this.columnDeleter = columnDeleter;
		
		initializer();
		runner();
		
	}
	
	
	public void initializer() throws IOException {
		
		output = new File(vsFile.getAbsolutePath().substring(0, vsFile.getAbsolutePath().lastIndexOf(".")) + "_snpEffed.vs");
		
		reder = new BufferedReader(new FileReader(vsFile));
		righter = new PrintWriter(output);
		
		line = reder.readLine();
		lineArr = line.split("\t");
		
		
		for(int i = 0; i < lineArr.length; i++) {
			
			lineArr[i] = lineArr[i].trim();
			
		}
		
		
		lineArr = truncator(lineArr);
		
		
		headz = java.util.Arrays.asList(lineArr);
		
		columns = lineArr.length;
		
		genotypeIndices(lineArr);
		glRenamer();
		
		
		chrIndex = headz.indexOf("Chr");
		leftFlankIndex = headz.indexOf("LeftFlank");	
		refIndex = headz.indexOf("ref_allele");
		varIndex = headz.indexOf("var_allele");
		mutTypeIndex = headz.indexOf("muttype");
		
		
		if(columnDeleter) {
			
			getDeleteColumns(lineArr);
			
		}
		
		
		String output = "";
		
		for(int i = 0; i < lineArr.length-1; i++) {
			
			output = output.concat(lineArr[i]).concat("\t").replaceAll(" ", "");
			
		}
		
		
		output = output.concat(lineArr[lineArr.length-1]).replaceAll(" ", "");
				
		
		
		righter.println(output);
		
	}
	
	
	public File getOutput() {
		
		return output;
		
	}
	
	public int getChrIndex() {
		
		return chrIndex;
	}
	
	public int getLeftFlankIndex() {
		
		return leftFlankIndex;
		
	}
	
	public ArrayList<Integer> getToDelete() {
		
		
		return toDelete;
		
		
	}
	
	public String[] truncator(String [] linezArr) throws IOException {
		
		headerConfig = new File("header_config.txt");
		if(!headerConfig.exists()) {
			
			headerConfig = getConfigFile();
			
		}
		
		BufferedReader configReder = new BufferedReader(new FileReader(headerConfig));
		
		String Line;
		
		while((Line = configReder.readLine()) != null) {
			
			if (!Line.startsWith("#") && !Line.isEmpty() && !Line.trim().equals("") && !Line.trim().equals("\n")) {
				
				String[] headerLine = Line.split("=");
				
				String old = headerLine[0];
				String notSoOld = headerLine[1];
				
				
				for(int i = 0; i < linezArr.length; i++) {
						
					if(linezArr[i].equals(old)) {
						
						linezArr[i] = notSoOld;
						
						if(linezArr[i].equals("AlleleCount")) {
							
							acIndex= i;
							
						} else if(linezArr[i].equals("AlleleFreq")) {							
							afIndex= i;
								
						} else if(linezArr[i].equals("TotAlleleNum")) {
							
							totIndex= i;
							
						}
							
							
						break;
						
						
					}
				}
				
			}
			
		}
		
		
		configReder.close();
		
		return linezArr;
	}
	
	
	public void glRenamer() throws IOException {
		
		glConfig = new File("glchrom_config.txt");
		if(!glConfig.exists()) {
			
			glConfig = getGLConfigFile();
			
		}
		
		BufferedReader configReder = new BufferedReader(new FileReader(glConfig));
		
		String Line;
		
		while((Line = configReder.readLine()) != null) {
			
			if (!Line.startsWith("#") && !Line.isEmpty() && !Line.trim().equals("") && !Line.trim().equals("\n")) {
				
				String[] configline = Line.split("=");
				
				oldGl.add(configline[0]);
				newGl.add(configline[1]);
				
				
			}
			
		}
		
		configReder.close();
	}
	
	
	//creates a FileChooser that allows the selection of an input VS File
	public File getConfigFile() {
		JFileChooser browseTo = new JFileChooser(); //creates a file chooser object
		
		JOptionPane.showMessageDialog(new Frame("Input prompt"), "Config file not found. Please select the header config file.");
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
	
	
	//creates a FileChooser that allows the selection of an input VS File
		public File getGLConfigFile() {
			JFileChooser browseTo = new JFileChooser(); //creates a file chooser object
			
			JOptionPane.showMessageDialog(new Frame("Input prompt"), "Config file not Found. Please select the GL Renamer config file.");
			int returnVal = browseTo.showOpenDialog(new JPanel()); //opens a window for user to browse to vs file to be edited
			
		    if(returnVal == JFileChooser.APPROVE_OPTION) {
		    	return browseTo.getSelectedFile(); //gets the path of the selected file
		     } else {
					int i = this.wantToContinue("No Config file selected."); 
					if (i == 1) {
						System.exit(0);
					}
					return this.getGLConfigFile();
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
	
	
	
	public void runner() throws IOException {
		
		
		//UI canceler
		JFrame jframe = canceler();
		jframe.setLocationRelativeTo(null);
		jframe.setVisible(true);
		
		bRunner = System.currentTimeMillis();
		
		//created executor to organize threads
		executor = Executors.newWorkStealingPool(NTHREDS);
		
		createSortIntoTemp();
		
		executor.shutdown();
		
		//waiting for threads to be finished
		try {
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		

		//closing readers
		System.out.println("Closing Readers");
		for(int i = 0; i < sepFiles; i++) {
			
			readIn.get(i).close();
			readOut.get(i).close();
			
		}

		//merging files
		System.out.println("Merging Files");
		String out = "";
		for(int i = 0; i < sepFiles; i++) {
			
			readIn.set(i, new BufferedReader(new FileReader(tempOut.get(i))));
			
			while((out = readIn.get(i).readLine()) != null) {
				
				righter.println(out);
			}
		}
		
		aRunner = System.currentTimeMillis();
		
		jframe.dispose();
		
		//ending screen and time
		
		timeMessage = ("The SnpEffEr Program Took " + (aRunner - bRunner) + "ms ("
		+ ((aRunner - bRunner) / 1000) + "." + (((aRunner - bRunner) % 1000) / 100) + " seconds)\n\n");	
		reder.close();
		righter.close();
		

		//deleting files and closing stuff
		for(int i = 0; i < sepFiles; i++) {
			readIn.get(i).close();
			readOut.get(i).close();
			tempIn.get(i).delete();
			tempOut.get(i).delete();
			
		}
		
		dir.delete();
		
		
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
		
		JLabel jtext = new JLabel("Snp Eff Program Running. If you'd like to abort, hit the cancel button below.", SwingConstants.CENTER);
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
					
					reder.close();
					
					righter.close();
					

					//deleting files and closing stuff
					for(int i = 0; i < sepFiles; i++) {
						
						readIn.get(i).close();
						readOut.get(i).close();
						tempIn.get(i).delete();
						tempOut.get(i).delete();
						
						
					}
					dir.delete();
					
					
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
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
		line1 = reder.readLine();
		line2 = reder.readLine();
		
		
		//Directory of where temp files go

		int suffix = vsFile.getName().lastIndexOf(".");
		if(suffix == -1) {
			
			suffix = vsFile.getName().length();
			
		}
		dir = new File(vsFile.getParent() + "\\SnpEfferTempFiles" + vsFile.getName().substring(0, suffix));
		dir.mkdir();
		dir.deleteOnExit();
		
		//loop until end of varsifter file
		while(line2 != null) {
			
			
			//Creating new temp files
			tempIn.add(File.createTempFile("Part " + Integer.toString(sepFiles+1), ".txt", dir));
			tempIn.get(sepFiles).deleteOnExit();
			
			tempOut.add(File.createTempFile("Part " + Integer.toString(sepFiles+1), ".txt", dir));
			tempOut.get(sepFiles).deleteOnExit();
			
			readOut.add(new PrintWriter(tempIn.get(sepFiles)));
			

			
			
			do {
				
				//prints to temp 
				readOut.get(sepFiles).println(line1);
				
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
							line2 = reder.readLine();
							break;
							
						}
						
					} else if(!(lineArr1[chrIndex].substring(0, 4).equals(lineArr2[chrIndex].substring(0, 4)))) {
						
						line1 = line2;
						line2 = reder.readLine();
						break;
						
					}
					*/
					
					//changes line and exits loop
					line1 = line2;
					line2 = reder.readLine();
					break;
					
					
				}
				//changes line
				line1 = line2;
				line2 = reder.readLine();
				
			} while(true);
			
			//sets buffered readers and writers and adds to new chromosome numbers
			readIn.add(new BufferedReader(new FileReader(tempIn.get(sepFiles))));
			readOut.get(sepFiles).close();
			readOut.set(sepFiles, new PrintWriter(tempOut.get(sepFiles)));
			
			//Starts Thread Execution After Temp File is made
			Runnable worker = new SnpEffLineChanger(readIn.get(sepFiles), readOut.get(sepFiles), naPos, naDistanceApart, chrIndex, refIndex, varIndex,
					mutTypeIndex, afIndex, totIndex, acIndex, oldGl, newGl, sepFiles);
			executor.execute(worker);
		
			
			sepFiles++;
			
		}
	}
	
	
	public String getTimeMessage() {
		
		return timeMessage;
		
	}
	
	
	//prompts user if they want to continue if they've hit cancel
	public int wantToContinue(String input) {
		int j = JOptionPane.showConfirmDialog(new Frame("Input prompt"), input + " Do you want to continue?", "User Prompt", 
				JOptionPane.YES_NO_OPTION);
		return j;
	}
	
}
