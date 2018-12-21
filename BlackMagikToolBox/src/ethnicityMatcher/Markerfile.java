package ethnicityMatcher;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import htsjdk.samtools.SamReader;

public class Markerfile {
	private static BufferedReader markerfile;

	/*
	 * Returns score for each single ethnicity evaluated.
	 */
	public static double MarkerScore(String A, String config, int number, SamReader samReader, int coverage, int score,
			int scale, PrintWriter report2, ArrayList<String> bamChrHeadz) throws IOException {

		/// Open and read appropriate ancestry informative marker file
		String AB = config + "/" + A + number + ".txt";
		markerfile = new BufferedReader(new FileReader(AB));
		String[] firstline = markerfile.readLine().split("\t");

		/// Initialize count of the number of markers hit
		int markercount = 0;

		/// Initialize ethnicity marker score
		double ethnicityscore = 0;
		ArrayList<String> header = new ArrayList<String>();

		for (int i = 0; i < firstline.length; i++) {
			header.add(firstline[i]);
		}

		int chr = header.indexOf("Chr");
		int left = header.indexOf("LeftFlank");
		int right = header.indexOf("RightFlank");
		int ref = header.indexOf("ref_allele");
		int alt = header.indexOf("var_allele");
		int hetfreq = header.indexOf("Het_Score");
		int horfreq = header.indexOf("Hor_Score");
		int hoafreq = header.indexOf("Hoa_Score");

		/// Avg allele frequency is used as prior probability for genotyping
		int varfreq = header.indexOf("1KGen_Allele_freq");

		String Line = markerfile.readLine();
		ArrayList<String[]> genotypeline = new ArrayList<String[]>();

		while (Line != null) {
			String[] cursplit = Line.split("\t");
			if (cursplit[0].contains("Group")) {

				/// Calculate log ratio score of most likely genotype call (homVar, homRef, Het)
				double tmp = GenotypeGroupScore.variantdensity(genotypeline, samReader, chr, left, right, ref, alt,
						hetfreq, horfreq, hoafreq, varfreq, coverage, score, bamChrHeadz);

				// If no genotype is found in the group, the default tmp value is -7777
				if (tmp != -7777) {
					markercount++;
					ethnicityscore += tmp;
				}
				genotypeline = new ArrayList<String[]>();

			} else {
				genotypeline.add(cursplit);
			}

			Line = markerfile.readLine();
		}

		/// If >45 markers are used, ethnicity score is scaled.
		if (markercount > (scale - 20)) {
			ethnicityscore = ethnicityscore / (double) markercount * (double) scale;
		}
		// Otherwise, return a signal that not enough markers are hit in this marker set
		else {
			ethnicityscore = -7777;
		}

		report2.println(A + number + "\t" + markercount + "\t" + ethnicityscore);
		return ethnicityscore;
	}
}
