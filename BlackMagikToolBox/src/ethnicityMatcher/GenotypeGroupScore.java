package ethnicityMatcher;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class GenotypeGroupScore {
	
	/*
	 * Returns log ratio for the most likely genotype call (HomRef, HomVar, or Het)
	 */
	public static double variantdensity(ArrayList<String[]> grouplines, SamReader samReader, int chro, int leftposition,
			int rightposition, int ref, int alt, int het, int hor, int hoa, int varfreq, int coverage, int score, ArrayList<String> bamChrHeadz) throws IOException {
		/// Set default value if no genotype is confidently called in the group marker
		/// set
		double defaultvalue = -7777;

		/// Iterate through each line in the group set
		for (int i = 0; i < grouplines.size(); i++) {
			String[] line = grouplines.get(i);
			String chrom = line[chro];
			int place = Integer.parseInt(line[leftposition]) + 1;
			int endplace = Integer.parseInt(line[rightposition]) - 1;
			String refallele = line[ref];
			String altallele = line[alt];

			/// The variant frequency is used as the prior for genotyping
			double varfrequency = Double.parseDouble(line[varfreq]);
			int length = endplace - place;

			/// Extract precalculated log ratio score for HomRef, Het, and HomVar genotypes
			double[] scores = { Double.parseDouble(line[hor]), Double.parseDouble(line[het]),
					Double.parseDouble(line[hoa]) };
			HashMap<String, Integer> genotypemap = new HashMap<String, Integer>();
			
			/// Determine the error rate
			double Snrscore = ethSNR.variantdensity(chrom, place, samReader, genotypemap, length, bamChrHeadz);

			/// Determine the error rate, cap at 0.5
			double error = 0.001;
			if (Snrscore > 1) {
				error = 1 / Snrscore * 0.2;
			} else if (error > 0.25) {
				error = 0.25;
			}

			/// Handles deletions
			String altt = altallele;
			if (refallele.length() > altallele.length()) {
				altallele = altt.substring(0, 1);
				for (int ii = 0; ii < (refallele.length() - altt.length()); ii++) {
					altallele += "D";
				}
				altallele += altt.substring(1, altt.length());
			}

			/// Determine if genotype call is HomRef, HomVar, or Het
			int geno = genotyping(genotypemap, varfrequency, refallele, altallele, error, coverage, score);
			if (geno != -1) {
				return scores[geno];
			}
		}
		return defaultvalue;
	}

	/*
	 * Returns genotype score based on likelihoods genotype is HomRef, Het, or
	 * HomVar.
	 */
	public static int genotyping(HashMap<String, Integer> genotypemap, double varfreq, String ref, String var,
			double error, int coverage, int score) {

		int result = -1;

		/// If variant frequency is 0, a minimum variant frequency is assigned
		if (varfreq <= 0) {
			varfreq = (double) 1 / (double) 10000;
		}

		/// Extract the reference allele count and the alternate allele count
		int refcount = 0;
		int altcount = 0;

		/// Count reads in pileup with reference call and variant call
		if (ref.length() == 1 && var.length() == 1) { /// SNP

			for (String key : genotypemap.keySet()) {
				if (Character.toString(key.charAt(0)).equals(ref)) {
					refcount += genotypemap.get(key);
				} else if (Character.toString(key.charAt(0)).equals(var)) {
					altcount += genotypemap.get(key);
				}
			}
		} else {
			if (genotypemap.containsKey(ref)) {
				refcount = genotypemap.get(ref);
			}

			if (genotypemap.containsKey(var)) {
				altcount = genotypemap.get(var);
			}
		}

		double hor = 0;
		double het = 0;
		double hoa = 0;

		/// See if the position meets the coverage cutoff criteria
		if ((refcount + altcount) > coverage) {

			/// Bayesian calculator
			hor = Math.pow((1 - error), refcount) * Math.pow(error, altcount); /// HomRef
			het = Math.pow(0.5, (altcount + refcount)); /// Het
			hoa = Math.pow((1 - error), altcount) * Math.pow(error, refcount); /// HomAlt or HomVar

			double tot = hor + het + hoa;

			hor /= tot;
			het /= tot;
			hoa /= tot;

			/// Multiply with the population frequency
			hor *= Math.pow((1 - varfreq), 2);
			hoa *= Math.pow(varfreq, 2);
			het *= varfreq * (1 - varfreq) * 2;

			/// Normalize again
			tot = hor + hoa + het;

			hor /= tot;
			het /= tot;
			hoa /= tot;

			/// Calculate a new genotype score and check if it satisfies the score threshold
			if (hor >= het && hor >= hoa && newscore(1 - hor) > score) {	/// HomRef
				result = 0;
			} else if (hoa >= het && hoa >= hor && newscore(1 - hoa) > score) { /// HomVar
				result = 2;
			} else if (het >= hoa && het >= hor && newscore(1 - het) > score) {	/// Het
				result = 1;
			}
		}
		return result;
	}

	/*
	 * Returns a new genotype score based on probability call.
	 */
	public static int newscore(double n) {
		double d = 0;
		if (n == 0) {
			d = 358;
		}
		if (n > 0) {
			d = -10 * Math.log10(n);
		}

		if (d > 358) {	/// Set maximum score to 358 (heuristic)
			d = 358;
		}

		return (int) d;
	}

}
