package ethnicityMatcher;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class EthnicityFromBam {

	private SamReader familybam;
	private ArrayList<String> bamChrHeaders;
	
	private JFrame cancelFrame;

	/*
	 * Returns the ethnicity for one person at a time
	 */
	public String EthPipeLine(String bamfilelocation, String probandname, HashMap<String, String> ConfigMap,
			int setcoverage, int setscore, PrintWriter report2, HashMap<String, ArrayList<String>> bamChrMap) throws IOException {

		// Open the BAM file
		AddBam(bamfilelocation, bamChrMap);

		/// Allows user to abort the program in the middle
		cancelFrame = canceler(probandname);
		cancelFrame.setLocationRelativeTo(null);
		cancelFrame.setVisible(true);

		String personethnicity = "Unknown";

		/// Initialize the marker scores for each marker set
		double AFR = 0;
		double EAS = 0;
		double AMR = 0;
		double EUR = 0;

		/// Calculate the scores for each ethnicity

		/// Keep track of the number of African (AFR) marker sets used
		double AFRcount = 0;
		for (int i = 1; i < 5; i++) {
			/// Calculate the log ratio score from this marker set
			double AFRcoretmp = Markerfile.MarkerScore("AFR", ConfigMap.get("configAFR"), i, familybam, setcoverage,
					setscore, 65, report2, bamChrHeaders);
			/// If the marker set does have more than 55 markers
			if (AFRcoretmp != -7777) {
				AFR += AFRcoretmp;
				AFRcount++;
			}
		}

		/// Do the same for the unAFR marker set and the rest of the ethnicities
		for (int i = 1; i < 5; i++) {
			double AFRcoretmp2 = Markerfile.MarkerScore("unAFR", ConfigMap.get("configAFR"), i, familybam, setcoverage,
					setscore, 65, report2, bamChrHeaders);
			if (AFRcoretmp2 != -7777) {
				AFR += AFRcoretmp2;
				AFRcount++;
			}
		}
		/// Calculate the average AFR score
		AFR /= AFRcount;

		String decision = "Unknown";
		/// If the marker score is >3, then accept AFR as the ethnicity
		if (AFR > 3) {
			decision = "AFR";
		} else if (AFR < -3) {

			/// Repeat for East Asian (EAS) marker
			double EAScount = 0;
			for (int i = 1; i < 5; i++) {
				double EAStmp = Markerfile.MarkerScore("EAS", ConfigMap.get("configEAS"), i, familybam, setcoverage,
						setscore, 65, report2, bamChrHeaders);
				if (EAStmp != -7777) {
					EAS += EAStmp;
					EAScount++;
				}
			}
			for (int i = 1; i < 5; i++) {
				double EAStmp2 = Markerfile.MarkerScore("unEAS", ConfigMap.get("configEAS"), i, familybam, setcoverage,
						setscore, 65, report2, bamChrHeaders);
				if (EAStmp2 != -7777) {
					EAS += EAStmp2;
					EAScount++;
				}
			}
			EAS /= EAScount;

			if (EAS > 3) {
				decision = "EAS";

			} else if (EAS < -3) {
				// Repeat for Native American/Hispanic (AMR)
				double AMRcount = 0;
				for (int i = 1; i < 5; i++) {
					double AMRtmp = Markerfile.MarkerScore("AMR", ConfigMap.get("configAMR"), i, familybam, setcoverage,
							setscore, 65, report2, bamChrHeaders);
					if (AMRtmp != -7777) {
						AMR += AMRtmp;
						AMRcount++;
					}
				}
				for (int i = 1; i < 5; i++) {
					double AMRtmp2 = Markerfile.MarkerScore("unAMR", ConfigMap.get("configAMR"), i, familybam,
							setcoverage, setscore, 65, report2, bamChrHeaders);
					if (AMRtmp2 != -7777) {
						AMR += AMRtmp2;
						AMRcount++;
					}
				}
				AMR /= AMRcount;

				if (AMR > 3) {
					decision = "AMR";

				} else if (AMR < -3) {
					/// Repeat for European (EUR)
					double EURcount = 0;
					for (int i = 1; i < 5; i++) {
						double EURtmp = Markerfile.MarkerScore("EUR", ConfigMap.get("configEUR"), i, familybam,
								setcoverage, setscore, 65, report2, bamChrHeaders);
						if (EURtmp != -7777) {
							EUR += EURtmp;
							EURcount++;
						}
					}

					for (int i = 1; i < 5; i++) {
						double EURtmp2 = Markerfile.MarkerScore("unEUR", ConfigMap.get("configEUR"), i, familybam,
								setcoverage, setscore, 65, report2, bamChrHeaders);
						if (EURtmp2 != -7777) {
							EUR += EURtmp2;
							EURcount++;
						}
					}
					EUR /= EURcount;

					if (EUR > 3) {
						decision = "EUR";
					} else if (EUR < -3) {
						decision = "SAS"; /// South Asian (SAS)
					}
				}
			} else { // If EAS can't be accepted or rejected
				decision = "SuperUnknown";
			}

		} else { /// If AFR can't be accepted or rejected
			decision = "SuperUnknown";
		}
		
		report2.println(
				"Final_AFR_Score" + "\t" + "Final_EAS_Score" + "\t" + "Final_EUR_Score" + "\t" + "Final_AMR_Score");
		report2.println(AFR + "\t" + EAS + "\t" + EUR + "\t" + AMR);

		double EURvAMR = 0;
		double SASvAMR = 0;
		double SASvEUR = 0;

		double AMRneuEURscore = 0;
		double AMRneuSASscore = 0;
		double EURneuSASscore = 0;

		double AMRneuEURscore_c = 0;
		double AMRneuSASscore_c = 0;
		double EURneuSASscore_c = 0;

		/// First decide the four dominant ethnicities
		personethnicity = decision;
		if (personethnicity.equals("Unknown")) {
			/// Calculate the scores in each pairwise marker files
			double AMRneuEUR = MarkerfilePair.MarkerScore(ConfigMap.get("configAMRneuEUR"), familybam, setcoverage,
					setscore, 35, report2, bamChrHeaders);
			double AMRneuSAS = MarkerfilePair.MarkerScore(ConfigMap.get("configAMRneuSAS"), familybam, setcoverage,
					setscore, 35, report2, bamChrHeaders);
			double EURneuAMR = MarkerfilePair.MarkerScore(ConfigMap.get("configEURneuAMR"), familybam, setcoverage,
					setscore, 35, report2, bamChrHeaders);
			double EURneuSAS = MarkerfilePair.MarkerScore(ConfigMap.get("configEURneuSAS"), familybam, setcoverage,
					setscore, 35, report2, bamChrHeaders);
			double SASneuAMR = MarkerfilePair.MarkerScore(ConfigMap.get("configSASneuAMR"), familybam, setcoverage,
					setscore, 35, report2, bamChrHeaders);
			double SASneuEUR = MarkerfilePair.MarkerScore(ConfigMap.get("configSASneuEUR"), familybam, setcoverage,
					setscore, 35, report2, bamChrHeaders);

			if (AMRneuEUR != -7777) {
				AMRneuEURscore += AMRneuEUR;
				AMRneuEURscore_c++;
			}

			if (EURneuAMR != -7777) {
				AMRneuEURscore -= EURneuAMR;
				AMRneuEURscore_c++;
			}

			if (AMRneuSAS != -7777) {
				AMRneuSASscore += AMRneuSAS;
				AMRneuSASscore_c++;
			}

			if (SASneuAMR != -7777) {
				AMRneuSASscore -= SASneuAMR;
				AMRneuSASscore_c++;
			}

			if (EURneuSAS != -7777) {
				EURneuSASscore += EURneuSAS;
				EURneuSASscore_c++;
			}

			if (SASneuEUR != -7777) {
				EURneuSASscore -= SASneuEUR;
				EURneuSASscore_c++;
			}

			/// Calculate the average for each pairwise comparison
			if (AMRneuEURscore_c > 0) {
				AMRneuEURscore /= AMRneuEURscore_c;
			}
			if (AMRneuSASscore_c > 0) {
				AMRneuSASscore /= AMRneuSASscore_c;
			}
			if (EURneuSASscore_c > 0) {
				EURneuSASscore /= EURneuSASscore_c;
			}

			/// These sets do not need to be scaled, so the scale index is -1
			EURvAMR += MarkerfilePair.MarkerScore(ConfigMap.get("configEURvAMR"), familybam, setcoverage, setscore, -1,
					report2, bamChrHeaders);
			SASvAMR += MarkerfilePair.MarkerScore(ConfigMap.get("configSASvAMR"), familybam, setcoverage, setscore, -1,
					report2, bamChrHeaders);
			SASvEUR += MarkerfilePair.MarkerScore(ConfigMap.get("configSASvEUR"), familybam, setcoverage, setscore, -1,
					report2, bamChrHeaders);

			report2.println("EURvAMR" + "\t" + "SASvAMR" + "\t" + "SASvEUR" + "\t" + "AMRneuEURscore" + "\t"
					+ "AMRneuSASscore" + "\t" + "EURneuSASscore");
			report2.println(EURvAMR + "\t" + SASvAMR + "\t" + SASvEUR + "\t" + AMRneuEURscore + "\t" + AMRneuSASscore
					+ "\t" + EURneuSASscore);

			/// Determine the ethnicty using the pairwise scores
			personethnicity = CalculateBay.whichpairethnic(AMRneuEURscore, AMRneuSASscore, EURneuSASscore, EURvAMR,
					SASvAMR, SASvEUR);

		}

		/// Close the cancel dialog box window
		cancelFrame.dispose();
		return personethnicity;

	}

	/*
	 * Allows user to abort the run at any point during the run.
	 */
	public JFrame canceler(String trio) {
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("EthnicityMatching");
		jframe.setSize(500, 100);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);

		JLabel jtext = new JLabel("Determining the ethnicity compositions for " + trio + ". "
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

					familybam.close();

					System.gc();
					System.exit(0);
				} catch (IOException e1) {
				}
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

	/*
	 * Retrieves BAM files for the whole family.
	 */
	public void AddBam(String bamlocation, HashMap<String, ArrayList<String>> bamChrMap) throws IOException {
		
		File bamFile = new File(bamlocation);
		SamReaderFactory srf = SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.LENIENT);
		SamReader samR = srf.open(bamFile);
		familybam = samR;
		
		//retrieves BAM chr list for each person
		bamChrHeaders = bamChrMap.get(bamlocation);
		
	}

}
