package confettiFilter;
import java.util.HashSet;
import java.util.Set;

//import PPDNFilter.PutativeVarCall;

public class PutativeVarCall {
	private int pos; /// Genomic posiiton
	private String varCall; /// Variant call
	public HashSet<String> prevReadSet; /// Set of previously read reads
	public HashSet<String> callSet; /// Record of calls at this position
	/// Constructor
	/// Maintains count of times a particular variant call has been made at the
	/// position along with associated reads with that call

	public PutativeVarCall() {
		varCall = "";
		pos = 0;
		prevReadSet = new HashSet<String>();
		callSet = new HashSet<String>();
	}

	public PutativeVarCall(String s, int p, String read) {
		varCall = s;
		pos = p;
		prevReadSet = new HashSet<String>();
		prevReadSet.add(read);
		callSet = new HashSet<String>();
	}

	/// Maintain list of previously read reads
	public void addRead(String s) {
		prevReadSet.add(s);
	}

	public boolean equals(PutativeVarCall pvc) {
		if ((this.varCall).equals(pvc.varCall)) {
			return true;
		}
		return false;
	}

	public String getCall() {
		return varCall;
	}
	
	public int getPos() {
		return pos;
	}
	
	public int getCount() {
		return prevReadSet.size();
	}
	
	public String printSet() {
		return set2Str(prevReadSet);
	}
	
	public void addCall(String s) {
		callSet.add(s);
	}
	
	public String printOtherCalls() {
		return set2Str(callSet);
	}

	public static String set2Str(Set<String> set) {
		StringBuilder sb = new StringBuilder();
		for (String s : set) {
			sb.append(s);
			sb.append("; ");
		}
		return sb.toString();
	}

}
