package salvagePathway;

import general.*;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;

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
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

public class SalvagePipeline {

	private BufferedReader vsData;
	private File vsFile;
	private PrintWriter filterWriter;
	// the report format should be reconsidered and what is put in it
	private PrintWriter report;// a report is generated automatically to reflect what has changed
	private String Line;
	private String[] curLine;
	private ArrayList<String> headers;
	private ArrayList<SamReader> familybam;
	
	private ArrayList<ArrayList<String>> bamChrHeaders  = new ArrayList<ArrayList<String>>();
	
	private int chrIndex;
	private int posIndex;
	private int posendIndex;
	private int refIndex;
	private int varIndex;
	private int muttypeIndex;
	private int motherfreq;
	private int fatherfreq;

	private static int minorField;
	private static int popAFFieldIndex;
	private static int popACFieldIndex;
	private static int sasAFFieldIndex;
	private static int easAFFieldIndex;
	private static int afrAFFieldIndex;
	private static int eurAFFieldIndex;
	private static int amrAFFieldIndex;

	private JFrame cancelFrame;
	private int CADD;

	private ArrayList<Integer> family;
	
	private File output;
	
	public void initializer(File varsifterFile, String[] trio, String[] bamlocation, String filename,
			String reportlocation, String[] parentethnic, BufferedReader pedigree, int minor, int yesPop, int[] populationIndexes, HashMap<String, ArrayList<String>> bamChrMap) throws IOException {
		
		// check if the BAM files are complete
		if (bamlocation.length < trio.length) {
			JOptionPane.showMessageDialog(new JPanel(),
					"Please make sure the BAM files are included for everyone in the family because the siblings will be salvaged using their bam files. The system is now exiting.");
			System.gc();
			System.exit(0);
		}
		
		/// Initialize the VarSifter (VS) headers
		vsFile = varsifterFile;
		vsData = new BufferedReader(new FileReader(vsFile));
		
		Line = vsData.readLine();
		
		
		/// filterWriter.println(Line);
		curLine = Line.split("\t");
		headers = new ArrayList<String>();
		for (int i = 0; i < curLine.length; i++) {
			headers.add(curLine[i]);
		}

		/// Check for the essential columns in the input VS file
		chrIndex = headers.indexOf("Chr");
		posIndex = headers.indexOf("LeftFlank");
		posendIndex = headers.indexOf("RightFlank");
		refIndex = headers.indexOf("ref_allele");
		varIndex = headers.indexOf("var_allele");
		muttypeIndex = headers.indexOf("muttype");
		CADD = headers.indexOf("PHRED");
		
//		motherfreq = headers.indexOf("1KGen_Allele_freq");
//		fatherfreq = headers.indexOf("1KGen_Allele_freq");
//		popACFieldIndex = headers.indexOf("1KGen_Total_Allele_c");
		// refalleleminor = headers.indexOf("All_Gahlref_is_minor");
		// refalleleminor2 = headers.indexOf("Gahl_UDPNref_is_minor");

		/// Set field denoting if reference is minor
		//getMinorFieldHeader(headers.toArray(new String[headers.size()]));
		
		minorField = minor;
		
		/// Default allele frequency fields
		popAFFieldIndex = headers.indexOf("1KGen_Allele_freq");
		popACFieldIndex = headers.indexOf("1KGen_Total_Allele_c");
		afrAFFieldIndex = headers.indexOf("1KGen_AFR_freq");
		easAFFieldIndex = headers.indexOf("1KGen_EAS_freq");
		amrAFFieldIndex = headers.indexOf("1KGen_AMR_freq");
		eurAFFieldIndex = headers.indexOf("1KGen_EUR_freq");
		sasAFFieldIndex = headers.indexOf("1KGen_SAS_freq");

		/// Ask user if they want to input population allele frequency fields or use
		/// defaults
		/*
		int yesSetPop = askUserPopulation();
		if (yesSetPop == 1)	{
			getPopulationHeaders(headers.toArray(new String[headers.size()]));
		}
		*/
		if(yesPop == 1) {
			
			popAFFieldIndex = populationIndexes[0];
			popACFieldIndex = populationIndexes[1];
			afrAFFieldIndex = populationIndexes[2];
			easAFFieldIndex = populationIndexes[3];
			amrAFFieldIndex = populationIndexes[4];
			eurAFFieldIndex = populationIndexes[5];
			sasAFFieldIndex = populationIndexes[6];
			
		}
		
		/// Initialize the population frequencies used for the parents
		fatherfreq = popAFFieldIndex;
		motherfreq = popAFFieldIndex;
		
		if (parentethnic[0].equals("SAS")) {
			fatherfreq = sasAFFieldIndex;
		} else if (parentethnic[0].equals("EAS")) {
			fatherfreq = easAFFieldIndex;
		} else if (parentethnic[0].equals("AMR")) {
			fatherfreq = amrAFFieldIndex;
		} else if (parentethnic[0].equals("EUR")) {
			fatherfreq = eurAFFieldIndex;
		} else if (parentethnic[0].equals("AFR")) {
			fatherfreq = afrAFFieldIndex;
		}

		if (parentethnic[1].equals("SAS")) {
			motherfreq = sasAFFieldIndex;
		} else if (parentethnic[1].equals("EAS")) {
			motherfreq = easAFFieldIndex;
		} else if (parentethnic[1].equals("AMR")) {
			motherfreq = amrAFFieldIndex;
		} else if (parentethnic[1].equals("EUR")) {
			motherfreq = eurAFFieldIndex;
		} else if (parentethnic[1].equals("AFR")) {
			motherfreq = afrAFFieldIndex;
		}

		/// If any of the above columns are not included, give an error message
		if (chrIndex == -1 || posIndex == -1 || posendIndex == -1 || refIndex == -1 || varIndex == -1
				|| muttypeIndex == -1 || motherfreq == -1 || fatherfreq == -1 || popACFieldIndex == -1
				|| minorField == -1) {
			JOptionPane.showMessageDialog(new Frame("Error"),
					"Incomplete columns: please make sure the VarSifter file is annotated with Chr, LeftFlank, RightFlank, ref_allele, var_allele and muttype.");
			System.exit(0);
			vsData.close();
		}
		
		/// Retrieve the family relationships from the pedigree
		family = Pedigree.getPedigreefromString(trio, headers);
		
		/// Initialize the BAM files
		AddBam(bamlocation, bamChrMap);
		
		output = new File(filename);
		/// Open and begin writing to the output VS file.
		filterWriter = new PrintWriter(output);
		// print out the headers
		filterWriter.println(Line);
		
		/// Generate the output file
		report = new PrintWriter(new File(reportlocation));
		report.println(trio[0] + "\t" + parentethnic[0] + "\t" + trio[1] + "\t" + parentethnic[1]);
		report.println("Chr" + "\t" + "LeftFlank" + "\t" + "RightFlank" + "\t" + "ref_allele" + "\t" + "var_allele"
				+ "\t" + "muttype" + "\t" + "proband.oldcall" + "\t" + "mother.oldcall" + "\t" + "father.oldcall" + "\t"
				+ "proband.newcall" + "\t" + "probandScore" + "\t" + "mother.newcall" + "\t" + "motherScore" + "\t"
				+ "father.newcall" + "\t" + "fatherScore" + "\t" + "p.probandAA" + "\t" + "p.probandAB" + "\t"
				+ "p.probandBB" + "\t" + "p.motherAA" + "\t" + "p.motherAB" + "\t" + "p.motherBB" + "\t" + "p.fatherAA"
				+ "\t" + "p.fatherAB" + "\t" + "p.fatherBB" + "\t" + "Salvage" + "\t" + "Situation" + "\t" + "TotalCov"
				+ "\t" + "CADDScore" + "\t" + "ErrorRate");

		runline(trio);
	}
	
	public File getOutput() {
		
		return output;
		
	}
	
	/*
	 * Retrieve the BAM file one at a time for the whole family
	 */
	public void AddBam(String[] bamlocation, HashMap<String, ArrayList<String>> bamChrMap) throws IOException {
		// store all the bam files in an arraylist

		familybam = new ArrayList<SamReader>();
		for (int i = 0; i < bamlocation.length; i++) {
			File bamFile = new File(bamlocation[i]);
			// create a samreader factory
			SamReaderFactory srf = SamReaderFactory.make();
			srf.validationStringency(ValidationStringency.LENIENT);
			SamReader samR = srf.open(bamFile);
			familybam.add(samR);
			
			bamChrHeaders.add(bamChrMap.get(bamlocation[i]));
		}
	}

	/*
	 * Run the salvaging/re-genotyping process for each variant in the VS file.
	 */
	public void runline(String[] trio) throws IOException {
		// read in the second line of the varsifter file
		Line = vsData.readLine();
		String chrom;

		/// Allows the user to abort the program in the middle of the run.
		/// The canceling frame will indicate the proband's name.
		cancelFrame = canceler(trio);
		cancelFrame.setLocationRelativeTo(null);
		cancelFrame.setVisible(true);

		while (Line != null) {
			curLine = Line.split("\t");
			Boolean denovosave = false;

			/// Initialize header values
			chrom = curLine[chrIndex];
			String leftflank = curLine[posIndex];
			String rightflank = curLine[posendIndex];
			String refprint = curLine[refIndex];
			String altprint = curLine[varIndex];
			String oldfather = curLine[family.get(0)];
			String oldmother = curLine[family.get(1)];
			String oldchild = curLine[family.get(2)];
			int refisminor = Integer.parseInt(curLine[minorField]);

			/// Store the genotypes for all the families in an arraylist
			ArrayList<String> geno = new ArrayList<String>();
			for (int i = 0; i < family.size(); i++) {
				geno.add(curLine[family.get(i)]);
			}

			/// Test if the family contain non-reference allele
			if (Other.perfamily(geno, refprint, altprint, refisminor)) {

				/// Main salvaging/re-genotyping engine

				if (chrom.contains("M") || chrom.contains("_") || chrom.contains("Y") || chrom.contains("X")) {
					/// skip over chrUn, chrM and chrY and chrX
				} else {

					/// Store the literal print for the reference and alt allele
					String ref = refprint;
					String alt = altprint;

					/// Define the region of interest
					int place = Integer.parseInt(leftflank) + 1;
					int endplace = Integer.parseInt(rightflank) - 1;

					/// Test if the allele isn't a SNP
					String mut = curLine[muttypeIndex];
					Boolean ifmut = !mut.equals("SNP");

					/// Initialize the length of the indel
					int length = endplace - place;

					/// Check if the trio is Mendelian inconsistent
					if (!Other.trueconsist(oldmother, oldfather, oldchild, ifmut)) {

						/// if not everyone are diploid, then it can't be salvaged.
						if (oldchild.equals("NA") || oldmother.equals("NA") || oldfather.equals("NA")) {

						}
						// Situation where there are a third allele, in which case, skip this place
						else if (Other.thirdallele(oldmother, oldfather, oldchild, refprint, altprint, ifmut)) {

						} else {

							if (Other.Denovo(oldfather, oldmother, oldchild, refprint, altprint)) {
								denovosave = true; /// De novo
							}

							/// For deletion site, add "D" characters into the places of deletion in the
							/// alt reference
							if (refprint.length() > altprint.length()) {
								alt = altprint.substring(0, 1);
								for (int i = 0; i < (refprint.length() - altprint.length()); i++) {
									alt += "D";
								}
								alt += altprint.substring(1, altprint.length());
							}

							/// Get population frequency for the parents' calls. The default is 0.5.
							double fathervarfreq = 0.5;
							double mothervarfreq = 0.5;

							if (Double.parseDouble(curLine[popACFieldIndex]) > 0) {
								fathervarfreq = Double.parseDouble(curLine[fatherfreq]);
								mothervarfreq = Double.parseDouble(curLine[motherfreq]);

								if (fathervarfreq <= 0) {
									fathervarfreq = (double) 1 / (Double.parseDouble(curLine[popACFieldIndex]) * 1.3);
								}

								if (mothervarfreq <= 0) {
									mothervarfreq = (double) 1 / (Double.parseDouble(curLine[popACFieldIndex]) * 1.3);
								}
							}

							/// Initialize the HashMap to store all the allele information which will be
							/// extracted from the BAM file (read sequence-># times that sequence is in the
							/// pileup)
							ArrayList<HashMap<String, Integer>> genotypemap = new ArrayList<HashMap<String, Integer>>();
							genotypemap.add(new HashMap<String, Integer>());
							genotypemap.add(new HashMap<String, Integer>());
							genotypemap.add(new HashMap<String, Integer>());

							/// Store the BAM files separately for the trio
							/// NOTE: this is a relic of a version of the salvage that doesn't salvage the
							/// siblings
							ArrayList<SamReader> triobam = new ArrayList<SamReader>();
							triobam.add(familybam.get(0));
							triobam.add(familybam.get(1));
							triobam.add(familybam.get(2));

							/// Determine the environmental error rate from the trio
							double Snrscore = SNR_Salvage.variantdensity(chrom, place, triobam, genotypemap, length, bamChrHeaders);
							/// Determine the error rate, cap at 0.25
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

							} else if (error > 0.25) {
								error = 0.25;
							}

							/// Determine the genotype call, allele count and genotype probabilities

							/// Call father's most likely genotype
							int[] fathercount = { 0, 0 };
							double[] fathergen = Callgenotype.callgenotype(genotypemap.get(0), ref, alt, fathervarfreq,
									error, fathercount);

							/// Call mother's most likely genotype
							int[] mothercount = { 0, 0 };
							double[] mothergen = Callgenotype.callgenotype(genotypemap.get(1), ref, alt, mothervarfreq,
									error, mothercount);

							/// Use parent's information to call child's most likely genotype
							int[] childcount = { 0, 0 };
							double[] childgen = Callgenotype.Childgenotype(genotypemap.get(2), ref, alt, error,
									childcount, fathergen, mothergen, fathercount, mothercount);

							/// Prepare to print out the final salvaged genotype
							String refref;
							String refalt;
							String altalt;
							if (mut.equals("INDEL")) {
								refref = refprint + ":" + refprint;
								refalt = refprint + ":" + altprint;
								altalt = altprint + ":" + altprint;
							} else {
								refref = refprint + refprint;
								refalt = refprint + altprint;
								altalt = altprint + altprint;
							}
							String[] genotypes = new String[] { refref, refalt, altalt, "NA" };

							/// If nobody is NA and if the call is still inconsistent, salvage the parents
							/// using the child.
							if ((int) fathergen[3] != 3 && (int) mothergen[3] != 3 && (int) childgen[3] != 3) {
								/// Determine which parent is inconsistent
								int ifconsist = Other.trueconsistindex(genotypes[(int) mothergen[3]],
										genotypes[(int) fathergen[3]], genotypes[(int) childgen[3]], ifmut);
								/// Check if there might be a hemizygous variant
								boolean hemy = false;

								if (((int) childgen[3] + (int) mothergen[3] + (int) fathergen[3]) == 3
										&& (int) childgen[3] != 1) {
									if (fathergen[(int) fathergen[3]] > 0.98 && mothergen[(int) mothergen[3]] > 0.98
											&& childgen[(int) childgen[3]] > 0.98
											&& (fathercount[0] + fathercount[1]) >= 5
											&& (mothercount[0] + mothercount[1]) >= 5
											&& (childcount[0] + childcount[1]) >= 5) {
										hemy = true;
									}
								}
								/// If both inconsistent
								if (ifconsist == 0) {
									/// Salvage the one with the lower read depth
									if ((mothergen[(int) mothergen[3]]) > fathergen[(int) fathergen[3]]) {
										fathergen = Callgenotype.Parentgenotype(fathergen, mothergen, childgen, hemy);
										if (Other.trueconsistindex(genotypes[(int) mothergen[3]],
												genotypes[(int) fathergen[3]], genotypes[(int) childgen[3]],
												ifmut) == 1) {
											mothergen = Callgenotype.Parentgenotype(mothergen, fathergen, childgen,
													hemy);
										}
									} else {
										mothergen = Callgenotype.Parentgenotype(mothergen, fathergen, childgen, hemy);
										if (Other.trueconsistindex(genotypes[(int) mothergen[3]],
												genotypes[(int) fathergen[3]], genotypes[(int) childgen[3]],
												ifmut) == 2) {
											fathergen = Callgenotype.Parentgenotype(fathergen, mothergen, childgen,
													hemy);
										}
									}
								}
								/// If mother is inconsistent
								else if (ifconsist == 1) {
									mothergen = Callgenotype.Parentgenotype(mothergen, fathergen, childgen, hemy);
								}
								/// If father is inconsistent
								else if (ifconsist == 2) {
									fathergen = Callgenotype.Parentgenotype(fathergen, mothergen, childgen, hemy);
								}
							}

							/// Calculate the genotype scores
							int fatherq = Other.newscore((1 - fathergen[(int) fathergen[3]]));
							int motherq = Other.newscore((1 - mothergen[(int) mothergen[3]]));
							int childq = Other.newscore((1 - childgen[(int) childgen[3]]));

							// store the printing result
							String[] fatherprint = { genotypes[(int) fathergen[3]], String.valueOf(fatherq),
									String.valueOf(fathercount[0] + fathercount[1]) };
							String[] motherprint = { genotypes[(int) mothergen[3]], String.valueOf(motherq),
									String.valueOf(mothercount[0] + mothercount[1]) };
							String[] childprint = { genotypes[(int) childgen[3]], String.valueOf(childq),
									String.valueOf(childcount[0] + childcount[1]) };

							/// If it's in an excessively noisy region (error rate above 0.03),
							/// assign one negative value (so it would be blacked on varsifter)

							if (error > 0.03) {
								curLine[family.get(0)] = fatherprint[0];
								curLine[family.get(0) + 1] = "-" + fatherprint[1];
								curLine[family.get(0) + 2] = fatherprint[2];

								curLine[family.get(1)] = motherprint[0];
								curLine[family.get(1) + 1] = "-" + motherprint[1];
								curLine[family.get(1) + 2] = motherprint[2];

								curLine[family.get(2)] = childprint[0];
								curLine[family.get(2) + 1] = "-" + childprint[1];
								curLine[family.get(2) + 2] = childprint[2];
								// do the siblings

								for (int i = 3; i < family.size(); i++) {
									/// sibling doesn't have weird allele that's not reference or alternate
									if (!Other.sibthirdallele(geno.get(i), refprint, altprint, ifmut)) {
										HashMap<String, Integer> sibgenotypemap = new HashMap<String, Integer>();

										/// salvage the sibling
										SNR_Salvage.SibVD(chrom, place, familybam.get(i), sibgenotypemap, length, bamChrHeaders.get(i));
										int[] sibcount = { 0, 0 };
										double[] sibgen = Callgenotype.Childgenotype(sibgenotypemap, ref, alt, error,
												sibcount, fathergen, mothergen, fathercount, mothercount);
										int sibq = Other.newscore((1 - sibgen[(int) sibgen[3]]));
										String[] sibprint = { genotypes[(int) sibgen[3]], String.valueOf(sibq),
												String.valueOf(sibcount[0] + sibcount[1]) };

										// print the siblings
										curLine[family.get(i)] = sibprint[0];
										curLine[family.get(i) + 1] = "-" + sibprint[1];
										curLine[family.get(i) + 2] = sibprint[2];
									}
								}

							} else {
								/// If it's not in a confetti region, negative both the score and coverage
								curLine[family.get(0)] = fatherprint[0];
								curLine[family.get(0) + 1] = "-" + fatherprint[1];
								curLine[family.get(0) + 2] = "-" + fatherprint[2];

								curLine[family.get(1)] = motherprint[0];
								curLine[family.get(1) + 1] = "-" + motherprint[1];
								curLine[family.get(1) + 2] = "-" + motherprint[2];

								curLine[family.get(2)] = childprint[0];
								curLine[family.get(2) + 1] = "-" + childprint[1];
								curLine[family.get(2) + 2] = "-" + childprint[2];

								for (int i = 3; i < family.size(); i++) {
									if (!Other.sibthirdallele(geno.get(i), refprint, altprint, ifmut)) {
										HashMap<String, Integer> sibgenotypemap = new HashMap<String, Integer>();
										SNR_Salvage.SibVD(chrom, place, familybam.get(i), sibgenotypemap, length, bamChrHeaders.get(i));
										int[] sibcount = { 0, 0 };
										double[] sibgen = Callgenotype.Childgenotype(sibgenotypemap, ref, alt, error,
												sibcount, fathergen, mothergen, fathercount, mothercount);
										int sibq = Other.newscore((1 - sibgen[(int) sibgen[3]]));
										String[] sibprint = { genotypes[(int) sibgen[3]], String.valueOf(sibq),
												String.valueOf(sibcount[0] + sibcount[1]) };
										curLine[family.get(i)] = sibprint[0];
										curLine[family.get(i) + 1] = "-" + sibprint[1];
										curLine[family.get(i) + 2] = "-" + sibprint[2];
									}
								}
							}

							String saved = "Yes";
							/// If it's still inconsistent
							if (Other.Pinconsist(Other.convert((int) mothergen[3]), Other.convert((int) fathergen[3]),
									Other.convert((int) childgen[3]))) {
								saved = "No";
							}

							/// Check what the final Mendelian state is
							String situations;
							String[] situation = { "compoundhet", "de novo", "other" };
							int situationscode = Other.situations(Other.convert((int) mothergen[3]),
									Other.convert((int) fathergen[3]), Other.convert((int) childgen[3]));
							situations = situation[situationscode];

							/// Total coverage
							int totalcount = fathercount[0] + fathercount[1] + mothercount[0] + mothercount[1]
									+ childcount[0] + childcount[1];

							/// If the variant is annotated with CADD
							String caddscore = "-777";
							if (CADD != -1) {
								caddscore = curLine[CADD];
							}
							
							/// Print the report
							String reportlin = chrom + "\t" + leftflank + "\t" + rightflank + "\t" + refprint + "\t"
									+ altprint + "\t" + mut + "\t" + oldchild + "\t" + oldmother + "\t" + oldfather
									+ "\t" + genotypes[(int) childgen[3]] + "\t" + String.valueOf(childq) + "\t"
									+ genotypes[(int) mothergen[3]] + "\t" + String.valueOf(motherq) + "\t"
									+ genotypes[(int) fathergen[3]] + "\t" + String.valueOf(fatherq) + "\t"
									+ childgen[0] + "\t" + childgen[1] + "\t" + childgen[2] + "\t" + mothergen[0] + "\t"
									+ mothergen[1] + "\t" + mothergen[2] + "\t" + fathergen[0] + "\t" + fathergen[1]
									+ "\t" + fathergen[2] + "\t" + saved + "\t" + situations + "\t" + totalcount + "\t"
									+ caddscore + "\t" + error;
							
							report.println(reportlin);
							
							if (situationscode != 1 && denovosave) {
								filterWriter.println(Line);
							}

							/// Reassemble the line
							Line = curLine[0];
							for (int i = 1; i < curLine.length; i++) {
								Line += "\t" + curLine[i];
							}

						}
					}

				}
				
				/// Store the new genotypes for the family
				geno = new ArrayList<String>();
				for (int i = 0; i < family.size(); i++) {
					geno.add(curLine[family.get(i)]);
				}

				/// If they are still not all reference
				if (Other.perfamily(geno, refprint, altprint, refisminor)) {
					filterWriter.println(Line);

				}
			}
			Line = vsData.readLine();

		}

		/// Close everything
		vsData.close();
		for (int i = 0; i < familybam.size(); i++) {
			familybam.get(i).close();
		}
		filterWriter.close();
		report.close();
		cancelFrame.dispose();
	}

	/*
	 * Allows user to abort the run at any point during the run.
	 */
	public JFrame canceler(String[] trio) {
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("Salvage Pathway");
		jframe.setSize(500, 100);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);

		JLabel jtext = new JLabel("Generating salvaged VarSifter file for " + trio[2] + ". "
				+ "If you'd like to abort, hit the cancel button below.", SwingConstants.CENTER);
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
				try {
					vsData.close();
				} catch (IOException e1) {
				}
				filterWriter.close();
				report.close();
				System.gc();
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
		jframe.pack();
		return jframe;
	}

}
