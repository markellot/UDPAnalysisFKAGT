package kaylaKode;
import java.util.ArrayList;

public class MendelianState {

	/*
	 * Returns 0 if the genotypes in the family members are consistent. Return 2 if
	 * someone's genotype data is missing. Return 1, otherwise.
	 */
	public static int trueconsist(String mom, String dad, String child, Boolean mut) {
		// 0 is consistent
		// 1 is somebody is inconsistent
		if (child.equals("NA") && mom.equals("NA") && dad.equals("NA")) {
			return 0;
		}
		if (child.equals("NA") || mom.equals("NA") || dad.equals("NA")) {
			return 2;
		}

		ArrayList<String> momalleles = alleles(mom, mut);
		ArrayList<String> dadalleles = alleles(dad, mut);
		ArrayList<String> childalleles = alleles(child, mut);

		if (childalleles.size() == 2) {
			if (dadalleles.contains(childalleles.get(0)) && momalleles.contains(childalleles.get(1))) {
				return 0;
			} else if (dadalleles.contains(childalleles.get(1)) && momalleles.contains(childalleles.get(0))) {
				return 0;
			}
		}
		return 1;

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
	 * Returns true if child's genotype appears hemizygous.
	 */
	public static boolean hemy(String child, String ref, String var) {

		if (ref.length() == 1 && var.length() == 1) {
			if (child.length() == 1) {
				return true;
			}
		} else if (!child.contains(":")) {
			return true;
		}

		return false;
	}

}
