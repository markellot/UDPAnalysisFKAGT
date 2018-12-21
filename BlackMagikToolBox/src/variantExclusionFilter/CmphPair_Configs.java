package variantExclusionFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

public class CmphPair_Configs {
	ArrayList<String> pairmatch;
	ArrayList<String> linepairs;
	private ArrayList<String> popFlags;

	/*
	 * Defines CmphPair (compound het pair) object. Stores the potential HalfHet
	 * (line) in the linepairs array (array of VS lines) and the potential HalfHet's
	 * other partners (pairindex) in the pairmatch array (array of indices).
	 */
	public CmphPair_Configs(String line, int pairindex) {
		linepairs = new ArrayList<String>();
		pairmatch = new ArrayList<String>();
		linepairs.add(line);
		String[] pairstores = line.split("\t")[pairindex].split(",");
		for (int i = 0; i < pairstores.length; i++) {
			pairmatch.add(pairstores[i]);
		}
	}

	/*
	 * If query gene's index matches that of the subject gene's, then they're a
	 * potential pair. Also, add the query gene's VS line into linepairs.
	 */
	public boolean pair(String index, String liny) {
		boolean result = false;
		for (int i = 0; i < pairmatch.size(); i++) {
			// some format of the file had the index changes so it prints index", this is to
			// rid of the " problem
			String tmp = pairmatch.get(i).replaceAll("\\D+", "");
			if (tmp.equals(index)) {
				result = true;
				linepairs.add(liny);
			}
		}
		return result;
	}

	/*
	 * Returns an array with the subject and query halfHets, if they're a successful
	 * compound het pairing.
	 */
	/*
	 * public String[] printOut(int Igahlhomvar, int Igahlhomref, int
	 * Irefalleleminor, int Ihomozygousexac, int Ihghomref, int Ihghomvar, int
	 * Iukmaf1, int Ikmaf1, int Iexacmaf1, int Ignommaf1, int Ignomexmaf1, int
	 * Igene, ArrayList<String> genelist) { if (linepairs.size() > 1) { for (int i =
	 * 1; i < linepairs.size(); i++) {
	 * 
	 * if (testHP(linepairs.get(0), linepairs.get(i), Igahlhomvar, Igahlhomref,
	 * Irefalleleminor, Ihomozygousexac, Ihghomref, Ihghomvar, Iukmaf1, Ikmaf1,
	 * Iexacmaf1, Ignommaf1, Ignomexmaf1, Igene, genelist)) { String genetmp1 =
	 * linepairs.get(0).split("\t")[Igene]; String genetmp2 =
	 * linepairs.get(i).split("\t")[Igene];
	 * 
	 * // System.out.println(genetmp1+" "+genetmp2); if (genetmp1.equals(genetmp2)
	 * || genetmp1.contains(genetmp2) || genetmp2.contains(genetmp1)) { if
	 * (containsGeneCmpd(genetmp2, genetmp1, genelist)) { return new String[] {
	 * linepairs.get(0), linepairs.get(i) }; } } } } } return new String[] { "" }; }
	 */

	/*
	 * Returns an array with the subject and query halfHets, if they're a successful
	 * compound het pairing.
	 */
	public String[] printOut(ArrayList<Integer> homVarIndsList, ArrayList<Integer> homRefIndsList, int Irefalleleminor,
			ArrayList<Integer> maxMAF1IndsList, int Igene, ArrayList<String> genelist, int homVarMin) {
		if (linepairs.size() > 1) {
			for (int i = 1; i < linepairs.size(); i++) {

				if (testHP(linepairs.get(0), linepairs.get(i), homVarIndsList, homRefIndsList, Irefalleleminor,
						maxMAF1IndsList, Igene, genelist, homVarMin)) {
					String genetmp1 = linepairs.get(0).split("\t")[Igene];
					String genetmp2 = linepairs.get(i).split("\t")[Igene];

					//System.out.println("Made it through testing for linepairs");
					if (genetmp1.equals(genetmp2) || genetmp1.contains(genetmp2) || genetmp2.contains(genetmp1)) {
						if (containsGeneCmpd(genetmp2, genetmp1, genelist)) {
							//System.out.println("printed something");
							return new String[] { linepairs.get(0), linepairs.get(i) };
						}
					}
				}
			}
		}
		return new String[] { "" };
	}

	/*
	 * For each gene in the genelist, if either gene1 or gene2 are in the
	 * genelist-->this particular pair has already been analyzed by the filter. If
	 * not, if gene1 and gene2 share the same name, then add it to the genelist.
	 */
	public static boolean containsGeneCmpd(String gene1, String gene2, ArrayList<String> genelist) {
		for (int i = 0; i < genelist.size(); i++) {
			if (genelist.get(i).contains(gene1) || genelist.get(i).contains(gene2) || gene1.contains(genelist.get(i))
					|| gene2.contains(genelist.get(i))) {
				return false;
			}
		}
		if (gene1.contains(gene2)) {
			genelist.add(gene2);
		} else {
			genelist.add(gene1);
		}
		return true;
	}

	/*
	 * Tests the subject and query genes (gene1, gene2) if they match the pop freq
	 * and count criteria for a potentially deleterious compound het. TODO: Add
	 * gnomAD info into this for the exact count and the MAFs.
	 */
	public boolean testHP(String line1, String line2, ArrayList<Integer> homVarIndsList,
			ArrayList<Integer> homRefIndsList, int Irefalleleminor, ArrayList<Integer> maxMAF1IndsList, int Igene,
			ArrayList<String> genelist, int homVarMin) {
		String[] line1cur = line1.split("\t");
		String[] line2cur = line2.split("\t");

		ArrayList<Integer> homVarList1 = new ArrayList<Integer>();
		ArrayList<Integer> homRefList1 = new ArrayList<Integer>();
		ArrayList<Boolean> maxMAF1List1 = new ArrayList<Boolean>();

		ArrayList<Integer> homVarList2 = new ArrayList<Integer>();
		ArrayList<Integer> homRefList2 = new ArrayList<Integer>();
		ArrayList<Boolean> maxMAF1List2 = new ArrayList<Boolean>();

		for (int hvInd : homVarIndsList) {
			homVarList1.add(Integer.parseInt(line1cur[hvInd]));
			homVarList2.add(Integer.parseInt(line2cur[hvInd]));
		}

		for (int hrInd : homRefIndsList) {
			homRefList1.add(Integer.parseInt(line1cur[hrInd]));
			homRefList2.add(Integer.parseInt(line2cur[hrInd]));
		}

		for (int maf1Ind : maxMAF1IndsList) {
			maxMAF1List1.add(Boolean.parseBoolean(line1cur[maf1Ind]));
			maxMAF1List2.add(Boolean.parseBoolean(line2cur[maf1Ind]));
		}

		int refminor = Integer.parseInt(line1cur[Irefalleleminor]);
		int refminor2 = Integer.parseInt(line2cur[Irefalleleminor]);

		/// Test if variants pass the population frequency filter
		boolean MAFv1 = false;
		for (boolean maxmaf1 : maxMAF1List1) {
			if (!maxmaf1) {
				MAFv1 = true;
			}
		}

		boolean MAFv2 = false;
		for (boolean maxmaf2 : maxMAF1List2) {
			if (!maxmaf2) {
				MAFv2 = true;
			}
		}

		if (refminor == 1) {
			MAFv1 = true;
			homVarList1 = homRefList1;
			// Ghomvar = Ghomref;
			// hghomvarc = hghomrefc;
		}

		if (refminor2 == 1) {
			MAFv2 = true;
			homVarList2 = homRefList2;
			// secondGhomvar = secondGhomref;
			// secondhghomvarc = secondhghomrefc;
		}

		popFlags = new ArrayList<String>();
		if (MAFv1 || MAFv2) {	// Passes pop freq filter
			/// Annotate and return based on absolute count
			boolean highCount = true; /// True if too many homozygous variants in population datasets
			for (int i = 0; i < homVarList1.size(); i++) {
				int homVarCount1 = homVarList1.get(i);
				int homVarCount2 = homVarList2.get(i);
				if (homVarCount1 < homVarMin || homVarCount2 < homVarMin) {
					highCount = false;
				} else {
					highCount = true;
					popFlags.add("HC");
					return false;
				}
			}

			if (highCount) {
				popFlags.add("HC");
				return false;
			}	else {
				popFlags.add("");
				return true;
			}
		} else {
			popFlags.add("HF");
			return false;
		}
	}

	/*
	 * Tests the subject and query genes (gene1, gene2) if they match the pop freq
	 * and count criteria for a potentially deleterious compound het. TODO: Add
	 * gnomAD info into this for the exact count and the MAFs.
	 */
	/*
	 * public boolean testHP(String line1, String line2, int Igahlhomvar, int
	 * Igahlhomref, int Irefalleleminor, int Ihomozygousexac, int Ihghomref, int
	 * Ihghomvar, int Iukmaf1, int Ikmaf1, int Iexacmaf1, int Ignommaf1, int
	 * Ignomexmaf1, int Igene, ArrayList<String> genelist) { String[] line1cur =
	 * line1.split("\t"); String[] line2cur = line2.split("\t");
	 * 
	 * int secondGhomref = Integer.parseInt(line2cur[Igahlhomref]); int
	 * secondGhomvar = Integer.parseInt(line2cur[Igahlhomvar]);
	 * 
	 * int secondrefminor = Integer.parseInt(line2cur[Irefalleleminor]);
	 * 
	 * int secondexac = Integer.parseInt(line2cur[Ihomozygousexac]);
	 * 
	 * int secondhghomrefc = Integer.parseInt(line2cur[Ihghomref]); int
	 * secondhghomvarc = Integer.parseInt(line2cur[Ihghomvar]);
	 * 
	 * boolean secondukmaf1b = Boolean.parseBoolean(line2cur[Iukmaf1]); boolean
	 * secondkmaf1b = Boolean.parseBoolean(line2cur[Ikmaf1]);
	 * 
	 * boolean secondexacmaf1b = Boolean.parseBoolean(line2cur[Iexacmaf1]);
	 * 
	 * boolean secondgnommaf1b = Boolean.parseBoolean(line2cur[Ignommaf1]); boolean
	 * secondgnomexmaf1b = Boolean.parseBoolean(line2cur[Ignomexmaf1]);
	 * 
	 * int Ghomref = Integer.parseInt(line1cur[Igahlhomref]); int Ghomvar =
	 * Integer.parseInt(line1cur[Igahlhomvar]);
	 * 
	 * int refminor = Integer.parseInt(line1cur[Irefalleleminor]);
	 * 
	 * int exac = Integer.parseInt(line1cur[Ihomozygousexac]);
	 * 
	 * int hghomrefc = Integer.parseInt(line1cur[Ihghomref]); int hghomvarc =
	 * Integer.parseInt(line1cur[Ihghomvar]);
	 * 
	 * boolean ukmaf1b = Boolean.parseBoolean(line1cur[Iukmaf1]); boolean kmaf1b =
	 * Boolean.parseBoolean(line1cur[Ikmaf1]);
	 * 
	 * boolean exacmaf1 = Boolean.parseBoolean(line1cur[Iexacmaf1]);
	 * 
	 * boolean gnommaf1b = Boolean.parseBoolean(line1cur[Ignommaf1]); boolean
	 * gnomexmaf1b = Boolean.parseBoolean(line1cur[Ignomexmaf1]);
	 * 
	 * boolean secondMAF = false; if (!secondukmaf1b && !secondkmaf1b &&
	 * !secondexacmaf1b && !secondgnommaf1b && !secondgnomexmaf1b) { secondMAF =
	 * true; }
	 * 
	 * boolean MAF = false; if (!ukmaf1b && !kmaf1b && !exacmaf1 && !gnommaf1b &&
	 * !gnomexmaf1b) { MAF = true; }
	 * 
	 * if (refminor == 1) { MAF = true; Ghomvar = Ghomref; hghomvarc = hghomrefc; }
	 * 
	 * if (secondrefminor == 1) { secondMAF = true; secondGhomvar = secondGhomref;
	 * secondhghomvarc = secondhghomrefc; }
	 * 
	 * // might do another carrier criteria : ghomvar<3 or second ghomvar<3 //
	 * if((MAF||secondMAF) && (Ghomvar<5 || secondGhomvar<5) && (hghomvarc<5 || //
	 * secondhghomvarc<5) && (exac<5 || secondexac<5)){ // if ((MAF || secondMAF) &&
	 * (Ghomvar < 3 || secondGhomvar < 3) && (hghomvarc < // 3 || secondhghomvarc <
	 * 3) && (exac < 3 || secondexac < 3)) {
	 * 
	 * Note to self: need to fix this somehow so that the flags are assigned for
	 * each halfHet, not for the pair as a whole.
	 * 
	 * // Annotate and return based on population frequency popFlags = new
	 * ArrayList<String>(); if (MAF || secondMAF) { // Annotate and return based on
	 * absolute count if ((Ghomvar < 3 || secondGhomvar < 3) && (hghomvarc < 3 ||
	 * secondhghomvarc < 3) && (exac < 3 || secondexac < 3)) { popFlags.add("");
	 * return true; } else { popFlags.add("HC"); return false; } } else {
	 * popFlags.add("HF"); return false; } }
	 */

	public String returnCMPopFlags() {
		if (!popFlags.isEmpty()) {
			StringBuilder sb = new StringBuilder();
			for (String s : popFlags) {
				sb.append(s);
				sb.append(";");
			}
			return sb.toString();
		} else {
			return "";
		}

	}

}
