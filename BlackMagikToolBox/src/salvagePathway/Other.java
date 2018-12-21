package salvagePathway;
import java.awt.Frame;
import java.util.ArrayList;

import javax.swing.JOptionPane;

public class Other {
	/*
	 * Used primarily in the Salvage Pathway.
	 */

	/*
	 * Converts given value (genotype (homref, het, and homvar) values to values
	 * compatible with Pinconsist's value system in the SalvagePipeline.
	 */
	public static int convert(double n) {
		int result = 0;
		if (n == 0) {
			result = 1;
		}
		if (n == 1) {
			result = 3;
		}
		if (n == 2) {
			result = 0;
		}

		return result;
	}

	/*
	 * Returns true if position is in an autosomal chromosome
	 */
	public static boolean autosome(int n) {

		if (n < 2699520) {
			return false;
		}
		if (n > 154931043) {
			return false;
		}
		if (88400000 < n && 92000000 > n) {
			return false;
		} else {
			return true;
		}

	}

	/*
	 * Optimize the order with which the parents are salvaged. 0 - mother and father
	 * are equally inconsistent 1 - mother is inconsistent 2 - father is
	 * inconsistent 3 - everyone's genotypes are consistent -1 - child is not
	 * diploid
	 */
	public static int trueconsistindex(String mom, String dad, String child, Boolean mut) {
		int incon = 0;
		ArrayList<String> momalleles = alleles(mom, mut);
		ArrayList<String> dadalleles = alleles(dad, mut);
		ArrayList<String> childalleles = alleles(child, mut);

		if (childalleles.size() == 2) {
			if (dadalleles.contains(childalleles.get(0)) && momalleles.contains(childalleles.get(1))
					|| dadalleles.contains(childalleles.get(1)) && momalleles.contains(childalleles.get(0))) {
				incon = 3;
			} else if (dadalleles.contains(childalleles.get(0)) && momalleles.contains(childalleles.get(0))
					|| dadalleles.contains(childalleles.get(1)) && momalleles.contains(childalleles.get(1))
					|| !dadalleles.contains(childalleles.get(1)) && !momalleles.contains(childalleles.get(1))) {
				incon = 0;
			} else if (dadalleles.contains(childalleles.get(0)) || dadalleles.contains(childalleles.get(1))) {
				incon = 1;
			} else {
				incon = 2;
			}
		} else {
			incon = -1;
		}

		return incon;

	}

	/*
	 * Evaluate if family is truly consistent. Returns true if consistent. Returns
	 * false, otherwise.
	 */
	public static boolean trueconsist(String mom, String dad, String child, Boolean mut) {
		if (child.equals("NA") && mom.equals("NA") && dad.equals("NA")) {
			return true;
		}
		if (child.equals("NA") || mom.equals("NA") || dad.equals("NA")) {
			return false;
		}

		ArrayList<String> momalleles = alleles(mom, mut);
		ArrayList<String> dadalleles = alleles(dad, mut);
		ArrayList<String> childalleles = alleles(child, mut);

		if (childalleles.size() == 2) {
			if (dadalleles.contains(childalleles.get(0)) && momalleles.contains(childalleles.get(1))) {
				return true;
			} else if (dadalleles.contains(childalleles.get(1)) && momalleles.contains(childalleles.get(0))) {
				return true;
			}
		}

		return false;

	}

	/*
	 * Return list of child's alleles at a given position.
	 */
	public static ArrayList<String> alleles(String child, Boolean mut) {

		ArrayList<String> alleles = new ArrayList<String>();
		if (mut && !child.contains(":")) {
			alleles.add(child);
		} else if (mut && child.contains(":")) {
			alleles.add(child.split(":")[0]);
			alleles.add(child.split(":")[1]);
		} else if (!mut) {
			alleles.add(Character.toString(child.charAt(0)));
			if (child.length() > 1) {
				alleles.add(Character.toString(child.charAt(1)));
			}
		}

		return alleles;
	}

	/*
	 * Returns true (variant line will be printed to output) if no one is NA in
	 * genotype and genotype does not contain ref.
	 */
	public static boolean perfamily(ArrayList<String> genotype, String refprint, String altprint, int refisminor) {
		String ref = refprint;
		String alt = altprint;

		if (refisminor == 1) {
			alt = refprint;
			ref = altprint;
		}

		boolean yesprint = false;
		int whoseNA = 0;

		for (int i = 0; i < genotype.size(); i++) {
			String gen = genotype.get(i);
			if (!gen.equals("NA")) {
				// SNP
				if (ref.length() == 1 && alt.length() == 1) {
					if (gen.length() == 2) {
						if (!Character.toString(gen.charAt(0)).equals(ref)
								|| !Character.toString(gen.charAt(1)).equals(ref)) {
							yesprint = true;
						}
					} else if (gen.length() == 2) {
						if (!Character.toString(gen.charAt(0)).equals(ref)) {
							yesprint = true;
						}
					}
				}
				// indel
				else if (gen.contains(":")) {
					if (!gen.split(":")[0].equals(ref) || !gen.split(":")[1].equals(ref)) {
						yesprint = true;
					}
				} else if (!gen.equals(ref)) {
					yesprint = true;
				}
			} else {
				whoseNA++;
			}
		}

		if (yesprint == false) {
			// add when the proband is the only one NA
			if (whoseNA == 1 && genotype.get(2).equals("NA")) {
				yesprint = true;
			}
		}

		return yesprint;
	}

	/*
	 * Returns true if a family member has an allele that is not the ref or alt
	 * specified in the VarSifter file.
	 */
	public static boolean thirdallele(String mother, String father, String child, String ref, String alt, Boolean mut) {

		ArrayList<String> momalleles = alleles(mother, mut);
		ArrayList<String> dadalleles = alleles(father, mut);
		ArrayList<String> childalleles = alleles(child, mut);

		for (int i = 0; i < momalleles.size(); i++) {
			if (!momalleles.get(i).equals(ref) && !momalleles.get(i).equals(alt)) {
				return true;
			}
		}

		for (int i = 0; i < dadalleles.size(); i++) {
			if (!dadalleles.get(i).equals(ref) && !dadalleles.get(i).equals(alt)) {
				return true;
			}
		}

		for (int i = 0; i < childalleles.size(); i++) {
			if (!childalleles.get(i).equals(ref) && !childalleles.get(i).equals(alt)) {
				return true;
			}
		}

		return false;
	}

	/*
	 * Returns true is a sibling has an allele that is not the ref or alt specified
	 * in the VarSifter file.
	 */
	public static boolean sibthirdallele(String child, String ref, String alt, Boolean mut) {

		ArrayList<String> childalleles = alleles(child, mut);

		for (int i = 0; i < childalleles.size(); i++) {
			if (!childalleles.get(i).equals(ref) && !childalleles.get(i).equals(alt)) {
				return true;
			}
		}

		return false;
	}

	/*
	 * Returns true if genotypes among family members are inconsistent.
	 */
	public static boolean Pinconsist(int mother, int father, int child) {
		int parentVal = mother + father; // stores the sum of the parents genotypes
		int probandVal = child;// store the proband genotype

		if (parentVal == 2) {
			if (probandVal != 1 && probandVal != -5) {
				return true;
			}
		}

		if (parentVal == 0) {
			if (probandVal != 0 && probandVal != -5) {
				return true;
			}
		}
		// One parent is het, one homref (3+1), the proband is homvar (0)
		if (parentVal == 4) {
			if (probandVal == 0) {
				return true;
			}
		}
		// One parent homvar, one parent homref (0+1), children are not het (!3, !-5)
		if (parentVal == 1) {
			if (probandVal != 3 && probandVal != 5) {
				return true;
			}
		}
		// One parent is het, one homvar (3+0), the proband is homref (1)
		if (parentVal == 3) {
			if (probandVal == 1) {
				return true;
			}
		}
		return false;
	}

	// {"compoundhet","de novo","other"};
	/*
	 * Returns an integer value depending on the apparent inheritance models among
	 * family members. 0 - compound het, 1 - de novo, 2 - default
	 */
	public static int situations(int mother, int father, int child) {
		int parentVal = mother + father; // stores the sum of the parents genotypes
		int probandVal = child;// store the proband genotype
		int situations = 2;

		// compound het
		if (parentVal == 3 && probandVal == 3) {
			situations = 0;
		} else if (parentVal == 4 && probandVal == 3) {
			situations = 0;
		}
		// de novo
		else if (parentVal == 2) {
			if (probandVal == 3) {
				situations = 1;
			}
		}

		else if (parentVal == 0) {
			if (probandVal == 3) {
				situations = 1;
			}
		}

		return situations;

	}

	/*
	 * Returns true if it is apparent child has a de novo mutation (salvaged de
	 * novo).
	 */
	public static boolean Denovo(String father, String mother, String child, String ref, String var) {
		String homRef;
		String homVar;
		String het1;
		String het2;

		if ((ref.length() > 1) || (var.length() > 1)) {
			homRef = ref + ":" + ref;
			homVar = var + ":" + var;
			het1 = ref + ":" + var;
			het2 = var + ":" + ref;
			// NA="NA";
		} else {
			homRef = ref + ref;
			homVar = var + var;
			het1 = ref + var;
			het2 = var + ref;
			// NA="NA";
		}

		if (child.equals(het1) || child.equals(het2)) {
			// parents homvar, child het
			if (father.equals(homVar) && mother.equals(homVar)) {
				return true;

				// homref or hemref
			} else if (father.equals(homRef) && mother.equals(homRef)) {
				return true;
				// het
			}
		}
		return false;
	}

	/*
	 * Calculate the new genotype score
	 */
	public static int newscore(double n) {
		double d = 0;
		if (n == 0) {
			d = 358;
		}
		if (n > 0) {
			// the maximum score, why is 358 I'm not sure
			d = -10 * Math.log10(n);
		}

		if (d > 358) {
			d = 358;
		}

		return (int) d;
	}
	
}
