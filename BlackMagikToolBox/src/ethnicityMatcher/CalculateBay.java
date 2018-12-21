package ethnicityMatcher;
import java.util.ArrayList;
import java.util.Collections;

public class CalculateBay {
	/*
	 * The logic to do pairwise comparison
	 */
	public static String whichpairethnic(double AMRneuEURscore, double AMRneuSASscore, double EURneuSASscore,
			double EURvAMR, double SASvAMR, double SASvEUR) {
		String output = "Unknown";

		/// The result for two sets of comparison scores
		String candidate = "Unknown";
		String candidate2 = "Unknown";

		/// Result of the first comparison set
		if (AMRneuSASscore > 3 && AMRneuEURscore > 3) {
			candidate = "AMR";
		} else if (EURneuSASscore > 3 && (-AMRneuEURscore) > 3) {
			candidate = "EUR";
		} else if ((-EURneuSASscore) > 3 && (-AMRneuSASscore > 3)) {
			candidate = "SAS";
		}

		/// Result of the second comparison set
		if ((-SASvAMR) > 3 && (-EURvAMR) > 3) {
			candidate2 = "AMR";
		}
		if ((-SASvEUR) > 3 && EURvAMR > 3) {
			candidate2 = "EUR";
		}
		if (SASvEUR > 3 && SASvAMR < 3) {
			candidate2 = "SAS";
		}

		/// Determine whether the two comparisons agree with each other.
		/// If they do, return the candidate. Otherwise, return "unknown."
		if (candidate2.equals(candidate)) {
			output = candidate;
		}

		return output;

	}

}
