package broadLevelBamFileCuration;

import general.*;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class AutoBamCuration {
	private BufferedReader vsData;
	private File vsFile;
	private PrintWriter filterWriter;

	private String Line;
	private String[] curLine;
	private ArrayList<String> headers;
	private ArrayList<SamReader> familybam;
	private ArrayList<ArrayList<String>> bamChrHeaders = new ArrayList<ArrayList<String>>();
	private HashMap<String, ArrayList<String>> bamChrMap;

	private int chrIndex;
	private int posIndex;
	private int posendIndex;
	private int refIndex;
	private int varIndex;
	private int muttypeIndex;
	private int refalleleminor;
	private int refalleleminor2;
	private int minorfield; /// Index of field denoting if reference is minor.
	private int naIndex;

	private ArrayList<Integer> family;
	private File pedFile;
	private File bamlo;
	private String destination;
	private File output;
	
	
	public AutoBamCuration(File vsFile, File ped, File bammy, String dest, int minor, HashMap<String, ArrayList<String>> bamChrMap) throws IOException {
		
		this.vsFile = vsFile;
		this.pedFile = ped;
		this.bamlo = bammy;
		this.destination = dest;
		this.minorfield = minor;
		this.bamChrMap = bamChrMap;
		
		Curator();
	}
	
	public void Curator() throws IOException {

		/// Input VarSifter file
		//File vsFile = TXTFile.getVSFile();

		/// Input pedigree file
		//File pedFile = TXTFile.getPedigreeFile();
		BufferedReader pedData = new BufferedReader(new FileReader(pedFile));

		/// Input BAM file directory
		HashMap<String, String> BamMap = new HashMap<String, String>();
		//File bamlo = TXTFile.getBAMDirectoryFile();
		String bamloPath = bamlo.getPath();
		FindBam.InitializeBam(bamloPath, BamMap);
		
		
		/// Initialize and run broad BAM curation

		/// Iterate through the families being analyzed
		String Line = pedData.readLine();

		JFrame jframe = canceler();
		jframe.setLocationRelativeTo(null);
		jframe.setVisible(true);
		

		while (Line != null) {
			String[] curLine = Line.split("\t");
			// String proband = curLine[2];
			String[] bamlocation = FindBam.MakeBamString(curLine, BamMap);
			initializer(vsFile, destination, curLine, bamlocation, pedData);
			Line = pedData.readLine();
		}
		
		pedData.close();
		jframe.dispose();
		
		
		
	}
	
	public File getOutput() {
		
		return output;
		
	}
	
	
	public void initializer(File vsFile, String outputfile, String[] trio, String[] bamlocation, BufferedReader pedigree)
			throws IOException {

		vsData = new BufferedReader(new FileReader(vsFile));
		output = new File(outputfile + ".vs");
		filterWriter = new PrintWriter(output);

		Line = vsData.readLine();

		curLine = Line.split("\t");
		headers = new ArrayList<String>();
		for (int i = 0; i < curLine.length; i++) {
			headers.add(curLine[i]);
		}


		chrIndex = headers.indexOf("Chr");
		posIndex = headers.indexOf("LeftFlank");
		posendIndex = headers.indexOf("RightFlank");
		refIndex = headers.indexOf("ref_allele");
		varIndex = headers.indexOf("var_allele");
		muttypeIndex = headers.indexOf("muttype");
//		refalleleminor = headers.indexOf("All_Gahlref_is_minor");
//		refalleleminor2 = headers.indexOf("Gahl_UDPNref_is_minor");

		genotypeIndices(curLine, curLine.length);

		ArrayList<Integer> family = Pedigree.getPedigreefromString(trio, headers);
		String firstline = curLine[0];
		for (int i = 1; i < naIndex; i++) {
			firstline += "\t" + curLine[i];
		}

		firstline += "\t" + "Error";
		firstline += "\t" + "SNR";
		
		for (int i = 0; i < family.size(); i++) {
			firstline += "\t" + curLine[family.get(i)] + "\t" + curLine[family.get(i) + 1] + "\t"
					+ curLine[family.get(i) + 2] + "\t" + "AllCount" + "\t" + "GoodCount";
		}

		filterWriter.println(firstline);

		AddBam(bamlocation, bamChrMap);

		Line = vsData.readLine();

		/// Iterate through each variant in the VS file
		while (Line != null) {
			curLine = Line.split("\t");

			/// Assign values
			String chrom = curLine[chrIndex];
			String leftflank = curLine[posIndex];
			String rightflank = curLine[posendIndex];
			String refprint = curLine[refIndex];
			String altprint = curLine[varIndex];
			// int refminor = Integer.parseInt(curLine[refalleleminor]);
			int refminor = Integer.parseInt(curLine[minorfield]);
			
			int place = Integer.parseInt(leftflank) + 1;
			int endplace = Integer.parseInt(rightflank) - 1;
			
			// /// Test if the allele is not SNP
			// String mut = curLine[muttypeIndex];
			// Boolean ifmut = !mut.equals("SNP");
			
			/// Initialize the length of the indel
			String ref = refprint;
			String alt = altprint;
			int length = endplace - place; /// Length of variant if it's an indel
			
			/// Mark if variant is an indel
			if (refprint.length() > altprint.length()) {
				alt = altprint.substring(0, 1);
				for (int i = 0; i < (refprint.length() - altprint.length()); i++) {
					alt += "D";
				}
				alt += altprint.substring(1, altprint.length());
			}
			
			/// Initialize maps for # reads and good reads in pileup
			/// (read sequence-># reads with that sequence in the pileup)
			ArrayList<HashMap<String, Integer>> genotypemap = new ArrayList<HashMap<String, Integer>>();
			for (int i = 0; i < familybam.size(); i++) {
				genotypemap.add(new HashMap<String, Integer>());
			}
			
			ArrayList<HashMap<String, Integer>> goodgenotypemap = new ArrayList<HashMap<String, Integer>>();
			for (int i = 0; i < familybam.size(); i++) {
				goodgenotypemap.add(new HashMap<String, Integer>());
			}

			/// Determine the environmental error rate from the trio (SNR)
			double Snrscore = SNR.variantdensity(chrom, place, familybam, genotypemap, goodgenotypemap, length, bamChrHeaders);
			
			/// Determine the error rate, cap at 0.5
			double error = 0.001;
			if (Snrscore > 1) {
				error = 1 / Snrscore * 0.25;
			}
			if (error <= 0.001) {
				error = 0.001;
				int totalread = 0;
				for (int i = 0; i < genotypemap.size(); i++) {
					totalread += genotypemap.get(i).keySet().size();
				}
				if (totalread < 11) {
					error = 0.005;
				}
			}
			
			/// Read BAM files for genotype calls and # reads & good reads
			ArrayList<String[]> Arrayreadcount = new ArrayList<String[]>();
			for (int i = 0; i < family.size(); i++) {
				String[] readcount = SNR.callgenotype(genotypemap.get(i), goodgenotypemap.get(i), ref, alt);
				Arrayreadcount.add(readcount);
			}

			/// Write output lines
			String out = curLine[0];

			// Copy the metadata for the variant
			for (int i = 1; i < naIndex; i++) {
				out += "\t" + curLine[i];
			}

			/// Add the SNR and error values
			out += "\t" + Double.toString(error);
			out += "\t" + Snrscore;

			/// For each family member, annotate with salvage status, # reads (ref, alt), #
			/// good reads (ref,alt)
			for (int i = 0; i < family.size(); i++) {
				String salvagestatus = "";
				if (Integer.parseInt(curLine[family.get(i) + 1]) < 0) {
					salvagestatus += "N";
				} else {
					salvagestatus += "P";
				}

				if (Integer.parseInt(curLine[family.get(i) + 2]) < 0) {
					salvagestatus += "N;";
				} else {
					salvagestatus += "P;";
				}
				out += "\t" + curLine[family.get(i)] + "\t" + curLine[family.get(i) + 1] + "\t"
						+ curLine[family.get(i) + 2] + "\t" + salvagestatus + Arrayreadcount.get(i)[0] + "\t"
						+ Arrayreadcount.get(i)[1];
			}

			/// Print output line to output
			filterWriter.println(out);

			/// Move to the next variant
			Line = vsData.readLine();
		}

		vsData.close();
		filterWriter.close();
	}

	/*
	 * Add family member's BAM file to list of family's BAM files
	 */
	public void AddBam(String[] bamlocation, HashMap<String, ArrayList<String>> bamChrMap) throws IOException {
		/// Store all the bam files in an arraylist
		familybam = new ArrayList<SamReader>();
		for (int i = 0; i < bamlocation.length; i++) {
			File bamFile = new File(bamlocation[i]);
			
			/// Create a samreader factory
			SamReaderFactory srf = SamReaderFactory.make();
			srf.validationStringency(ValidationStringency.LENIENT);
			SamReader samR = srf.open(bamFile);
			familybam.add(samR);
			
			//retrieves BAM chr list for each person
			bamChrHeaders.add(bamChrMap.get(bamlocation[i]));
		}
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
						//naDistanceApart = counter - naIndex;
					}
				tempVal++;			
				}
				counter++;
			}
		}
	}
	

	/*
	 * Ask user for input in headers. This is derived from the KaylaKode GUI code.
	 * The output file, input VS file, and pedigree files are included as inputs for
	 * the cancellation option, so that the PrintWriter and BufferedReader streams
	 * are successfully closed if the user chooses to abort the job.
	 */
	
	/* Unneeded
	
	public static void getHeaders(String[] scroller, PrintWriter output, BufferedReader inputVS,
			BufferedReader pedReader) {
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
				minorfield = msglist10.getSelectedIndex();
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
				if (minorfield == 0) {
					JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please input minor allele filter.");
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
				try {
					output.close();
					inputVS.close();
					pedReader.close();
				} catch (IOException e1) {

				}
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
	*/
	/*
	 * Add gene to genelist.
	 */
	public static void addgene(ArrayList<String> genelist, String gene) {
		if (genelist.indexOf(gene) == -1) {
			genelist.add(gene);
		}
	}

	/*
	 * Aborts run if user clicks "Cancel."
	 */
	public JFrame canceler() {
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("FBBamCuration");
		jframe.setSize(500, 100);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);

		JLabel jtext = new JLabel("Broad Level BAM file curation/ SNR calculation module.",
				SwingConstants.CENTER);
		gbc.gridx = 0;
		gbc.gridy = 0;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(5, 5, 5, 5);
		gbc.weightx = 1.0;
		jpanel.add(jtext, gbc);
		
		JLabel jtext2 = new JLabel("Generating the bam curate files. Use the cancel botton below to abort.",
				SwingConstants.CENTER);
		gbc.gridx = 0;
		gbc.gridy = 1;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(5, 5, 5, 5);
		gbc.weightx = 1.0;
		jpanel.add(jtext2, gbc);
		
		JButton abort = new JButton("Cancel");
		abort.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				jframe.dispose();
				try {
					vsData.close();
				} catch (IOException e1) {
				}
				System.gc();
				System.exit(0);
			}
		});

		gbc.gridx = 0;
		gbc.gridy = 2;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(5, 5, 5, 5);
		gbc.weightx = 1.0;
		jpanel.add(abort, gbc);
		
		jpanel.validate();
		jframe.getContentPane().add(jpanel);
		return jframe;
	}
}
