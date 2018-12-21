package salvagePathway;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class Callgenotype {
	/*
	 * Used primarily in the Salvage Pathway Returns array: { P(homref), P(het),
	 * P(homvar), which genotype is most likely} For 4th element: 0 - homref, 1 -
	 * homvar, 2 - het
	 */
	public static double[] callgenotype(HashMap<String, Integer> genotypemap, String ref, String var, double varfreq,
			double error, int[] allelecount) {
		/// If the variant frequency is low, assign the lowest 15000
		double[] result = { 0, 0, 0, 3 };

		int refcount = 0;
		int altcount = 0;

		if (ref.length() == 1 && var.length() == 1) {
			/// If this is a SNP place, the call at the position of interest may be either
			/// reference or variant of interest.
			for (String key : genotypemap.keySet()) {
				if (Character.toString(key.charAt(0)).equals(ref)) {
					refcount += genotypemap.get(key);
				} else if (Character.toString(key.charAt(0)).equals(var)) {
					altcount += genotypemap.get(key);
				}
			}

		}
		/// If variant is an indel, just get the count
		else {
			if (genotypemap.containsKey(ref)) {
				refcount = genotypemap.get(ref);
			}

			if (genotypemap.containsKey(var)) {
				altcount = genotypemap.get(var);
			}
		}

		if ((refcount + altcount) > 60) { /// if combined read count is above 60, scale the counts
			double totscale = (double) 60 / (refcount + altcount);
			refcount *= totscale;
			altcount *= totscale;
		}

		allelecount[0] = refcount;
		allelecount[1] = altcount;
		double hor = 0;
		double het = 0;
		double hoa = 0;

		// Ehrenfest model for likelihoods genotype calls are homref (hor), het, homvar
		// (hoa)
		if ((refcount + altcount) > 0) {
			hor = Math.pow((1 - error), refcount) * Math.pow(error, altcount);
			het = Math.pow(0.5, (altcount + refcount));
			hoa = Math.pow((1 - error), altcount) * Math.pow(error, refcount);

			/// Normalization
			double tot = hor + het + hoa;

			hor /= tot;
			het /= tot;
			hoa /= tot;

			/// Multiply by population frequency
			hor *= Math.pow((1 - varfreq), 2);
			hoa *= Math.pow(varfreq, 2);
			het *= varfreq * (1 - varfreq) * 2;

			tot = hor + hoa + het;

			hor /= tot;
			het /= tot;
			hoa /= tot;

			result[0] = hor;
			result[1] = het;
			result[2] = hoa;

			/// Check which genotype has the highest probability
			if (hor >= het) {
				if (hor >= hoa) {
					result[3] = 0;
				}
			} else if (hoa >= het) {
				if (hoa >= hor) {
					result[3] = 2;
				}
			} else if (het >= hoa) {
				if (het >= hor) {
					result[3] = 1;
				}
			}
		}

		return result;
	}

	/*
	 * Use parent's most likely genotypes and Bayesian probability logic to
	 * determine most likely genotype of child. Returns an array similar to the
	 * array returned by callgenotype (above).
	 */
	public static double[] Childgenotype(HashMap<String, Integer> genotypemap, String ref, String var, double error,
			int[] childcount, double[] father, double[] mother, int[] fathercount, int[] mothercount) {
		double[] result = { 0, 0, 0, 3 };
		int refcount = 0;
		int altcount = 0;

		if (ref.length() == 1 && var.length() == 1) {

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

		if ((refcount + altcount) > 60) {
			double totscale = (double) 60 / (refcount + altcount);
			refcount *= totscale;
			altcount *= totscale;
		}

		childcount[0] = refcount;
		childcount[1] = altcount;

		double hor = 0;
		double het = 0;
		double hoa = 0;

		if ((refcount + altcount) > 0) {
			hor = Math.pow((1 - error), refcount) * Math.pow(error, altcount);
			het = Math.pow(0.5, (altcount + refcount));
			hoa = Math.pow((1 - error), altcount) * Math.pow(error, refcount);

			double tot = hor + het + hoa;

			hor /= tot;
			het /= tot;
			hoa /= tot;

			// store it temporarily in result
			result[0] = hor;
			result[1] = het;
			result[2] = hoa;

			if (hor >= het) {
				if (hor >= hoa) {
					result[3] = 0;
				}
			} else if (hoa >= het) {
				if (hoa >= hor) {
					result[3] = 2;
				}
			} else if (het >= hoa) {
				if (het >= hor) {
					result[3] = 1;
				}
			}

			if ((int) father[3] != 3 && (int) mother[3] != 3) {
				/// Check if the hemizygous model fits
				boolean hemy = false;
				if (((int) result[3] + (int) mother[3] + (int) father[3]) == 3 && (int) result[3] != 1) {
					if (father[(int) father[3]] > 0.98 && mother[(int) mother[3]] > 0.98
							&& result[(int) result[3]] > 0.98 && (fathercount[0] + fathercount[1]) >= 5
							&& (mothercount[0] + mothercount[1]) >= 5 && (childcount[0] + childcount[1]) >= 5) {
						hemy = true;

					}
				}

				double E = 2e-7;

				/// if hemizygosity is possible, use the mutation rate of deletion
				if (hemy) {
					E = 1e-2;
				}
				// else if((childcount[0]+childcount[1])>50){
				// E=1e-5;
				// }
				double E2 = Math.pow(E, 2);

				/// population prior
				/// double[] hor, het, hoa

				double co1 = 1 / (1 + E + E2);
				double co2 = 1 / (1 + 0.5 * E);
				double co3 = 1 / (1 + 2 * E);

				double PAA = mother[0] * (co1 * father[0] + co2 * 0.5 * father[1] + co3 * E * father[2])
						+ mother[1] * (co2 * 0.5 * father[0] + 0.25 * father[1] + co2 * 0.5 * E * father[2])
						+ mother[2] * (co3 * E * father[0] + co2 * 0.5 * E * father[1] + co1 * E2 * father[2]);

				double PAB = mother[0] * (co1 * E * father[0] + co2 * 0.5 * father[1] + co3 * father[2])
						+ mother[1] * (co2 * 0.5 * father[0] + 0.5 * father[1] + co2 * 0.5 * father[2])
						+ mother[2] * (co3 * father[0] + co2 * 0.5 * father[1] + co1 * E * father[2]);

				double PBB = mother[0] * (co1 * E2 * father[0] + co2 * 0.5 * E * father[1] + co3 * E * father[2])
						+ mother[1] * (co2 * 0.5 * E * father[0] + 0.25 * father[1] + co2 * 0.5 * father[2])
						+ mother[2] * (co3 * E * father[0] + co2 * 0.5 * father[1] + co1 * father[2]);

				tot = PAA + PAB + PBB;
				PAA /= tot;
				PAB /= tot;
				PBB /= tot;

				tot = PAA * hor + PAB * het + PBB * hoa;

				hor *= PAA / tot;
				het *= PAB / tot;
				hoa *= PBB / tot;

				result[0] = hor;
				result[1] = het;
				result[2] = hoa;

				if (hor >= het) {
					if (hor >= hoa) {
						result[3] = 0;
					}
				} else if (hoa >= het) {
					if (hoa >= hor) {
						result[3] = 2;
					}
				} else if (het >= hoa) {
					if (het >= hor) {
						result[3] = 1;
					}
				}
			}
		}

		return result;
	}

	/*
	 * Determine most likely genotype for a parent based off of the other parent and
	 * the child. Uses the same Bayesian logic as above (callgenotype,
	 * Childgenotype) and returns an array similar to those methods.
	 */
	public static double[] Parentgenotype(double[] parent1, double[] parent2, double[] child, boolean hemy) {
		double[] result = parent1;

		double hor = result[0];
		double het = result[1];
		double hoa = result[2];

		double E = 2e-7;

		// add hemizygous possibility for deletion mutation, 1/10000 rate of observing
		// deletion
		if (hemy) {
			E = 1e-2;
		}

		double E2 = Math.pow(E, 2);
		double co1 = 1 / (1.5 + E);
		double co2 = 1 / (1.5 * E + E2);
		double co3 = 1 / (0.75 + 0.5 * E);
		double athird = (double) 1 / (double) 3;

		// population prior
		double PAA = parent2[0] * (co1 * child[0] + E * co1 * child[1] + E2 * co2 * child[2])
				+ parent2[1] * (0.5 * co3 * child[0] + athird * child[1] + 0.5 * E * co3 * child[2])
				+ parent2[2] * (E * co2 * child[0] + co1 * child[1] + E * co1 * child[2]);

		double PAB = parent2[0] * (0.5 * co1 * child[0] + 0.5 * co1 * child[1] + 0.5 * E * co2 * child[2])
				+ parent2[1] * (0.25 * co3 * child[0] + athird * child[1] + 0.25 * co3 * child[2])
				+ parent2[2] * (0.5 * E * co2 * child[0] + 0.5 * co1 * child[1] + 0.5 * co1 * child[2]);

		double PBB = parent2[0] * (E * co1 * child[0] + co1 * child[1] + E * co2 * child[2])
				+ parent2[1] * (0.5 * co3 * child[0] + athird * child[1] + 0.5 * co3 * child[2])
				+ parent2[2] * (E2 * co2 * child[0] + E * co1 * child[1] + co1 * child[2]);

		double tot = PAA + PAB + PBB;
		PAA /= tot;
		PAB /= tot;
		PBB /= tot;

		tot = PAA * hor + PAB * het + PBB * hoa;

		hor *= PAA / tot;
		het *= PAB / tot;
		hoa *= PBB / tot;

		result[0] = hor;
		result[1] = het;
		result[2] = hoa;

		if (hor >= het) {
			if (hor >= hoa) {
				result[3] = 0;
			}
		} else if (hoa >= het) {
			if (hoa >= hor) {
				result[3] = 2;
			}
		} else if (het >= hoa) {
			if (het >= hor) {
				result[3] = 1;
			}
		}

		return result;

	}

}
