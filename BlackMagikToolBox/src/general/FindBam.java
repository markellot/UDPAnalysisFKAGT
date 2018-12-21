package general;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;

import javax.swing.JOptionPane;
import javax.swing.JPanel;

public class FindBam {

	/*
	 * Finds the BAM file paths for the given samples
	 */
	public static String[] MakeBamString(String[] family, HashMap<String, String> bamfilemap) {
		String[] bamlocation = new String[family.length];

		for (int i = 0; i < bamlocation.length; i++) {
			if (bamfilemap.get(family[i]) != null) {
				bamlocation[i] = bamfilemap.get(family[i]);
			} else {
				JOptionPane.showMessageDialog(new JPanel(), "Bam file not found for " + family[i]);
				System.gc();
				System.exit(0);
			}
		}

		return bamlocation;

	}
	
	/*
	 * Find the paths for the BAM files given the UDP ID name. (For Ethnicity Matcher)
	 */
	public static String MakeBamString(String family, HashMap<String, String> bamfilemap) {

		if (bamfilemap.get(family) != null) {
			return bamfilemap.get(family);
		} else {
			JOptionPane.showMessageDialog(new JPanel(), "Bam file not found for " + family + ". \n"
					+ "Check pedigree file for excessive whitespaces.");
			System.gc();
			System.exit(0);
		}

		return "SomethingWeird";
	}

	/*
	 * Take the ethnicity config file or BAM file directory file and put each
	 * sample's information into a HashMap for quick reference later in the process.
	 */
	public static void InitializeBam(String bammaplocation, HashMap<String, String> bammap) throws IOException {
		File vsFile = new File(bammaplocation);
		BufferedReader vsData = new BufferedReader(new FileReader(vsFile));

		String Line = vsData.readLine();
		while (Line != null) {
			// filterWriter.println(Line);
			String[] curLine = Line.split("\t");
			// System.out.println(Line);
			bammap.put(curLine[0], curLine[1]);

			Line = vsData.readLine();
		}
		vsData.close();
	}

	/*
	 * Find the ethnicity for the given samples
	 */
	public static String[] MakeEthString(String[] family, HashMap<String, String> bamfilemap) {
		String[] bamlocation = new String[family.length];

		for (int i = 0; i < bamlocation.length; i++) {
			if (bamfilemap.get(family[i]) != null) {
				bamlocation[i] = bamfilemap.get(family[i]);
			} else {
				JOptionPane.showMessageDialog(new JPanel(), "Ethnicity not found for " + family[i]);
				System.gc();
				System.exit(0);
			}
		}
		return bamlocation;
	}

}
