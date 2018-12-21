import java.util.*;

import javax.swing.JOptionPane;

import java.awt.Frame;
import java.io.*;
import java.text.SimpleDateFormat;

public class Analysisv2 extends Thread{
	
	private BufferedReader red;
	private PrintWriter righter;
	private boolean[] selections;
	
	private BufferedReader refGenome;
	private String currentChrom;
	private String chromLine;
	private int lineLength;
	private int lineNumber;
	private String genomeDirectory;
	
	//values within fields
	private String refAllele;
	private String varAllele;
	private String leftFlank;
	
	private ArrayList<Integer> toDelete;
	
	//indexes of fields that may be indexed depending on the utility being run
	private int chrIndex;
	private int leftFlankIndex;
	private int refAlleleIndex;
	private int varAlleleIndex;
	private int mutTypeIndex;
	private long newIndex;
	
	private int naPos;
	private int naDistanceApart;
	private int vcfPos;
	
	
	private int init = -1;
	private String chrom = "";
	
	//time
	private long after;
	private long before;
	private long allTotal = 0;
	private long step1;
	private long oneTotal = 0;
	private long step2;
	private long twoTotal = 0;
	private long step3;
	private long threeTotal = 0;
	private long afterTotal = 0;
	
	//Constructor
	public Analysisv2(BufferedReader red, PrintWriter righter, boolean[] selections, String genomeDirectory, int chrIndex, int leftFlankIndex, 
			int refAlleleIndex, int varAlleleIndex, int mutTypeIndex, long newIndex, int naPos, int naDistApart, int vcfPos, ArrayList<Integer> toDelete) {
		
		this.red = red;
		this.righter = righter;
		this.selections = selections;
		currentChrom = "";
		this.genomeDirectory = genomeDirectory;
		this.chrIndex = chrIndex;
		this.leftFlankIndex = leftFlankIndex;
		this.refAlleleIndex = refAlleleIndex;
		this.varAlleleIndex = varAlleleIndex;
		this.mutTypeIndex = mutTypeIndex;
		this.newIndex = newIndex;
		this.naPos = naPos;
		this.naDistanceApart = naDistApart;
		this.vcfPos = vcfPos;
		
		this.toDelete = toDelete;
		
		
	}
	
	//Runner for each Thread
	@Override
	public void run() {
		
		String read;
		String[] readLineSplit;
		try {
			//while loop reading through each chromosome 
			while((read = red.readLine()) != null) {
				
				before = System.currentTimeMillis();
				
				try {
					
					readLineSplit = read.split("\t");
					
					//initializing the current chrom
					if(init == -1) {
						
						chrom = readLineSplit[chrIndex];
						init = 0;
					}
					
					analyze(readLineSplit, read);       //brute force of each line and passes each line into analysis method
				} catch (NumberFormatException | IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				//time
				after = System.currentTimeMillis();
				allTotal += (after-before);
				afterTotal += (after - step3);
				
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//thread statistics
		System.out.println("Thread Completed for Part: " + chrom);
		
		long yourmilliseconds = System.currentTimeMillis();
		SimpleDateFormat sdf = new SimpleDateFormat("MMM dd,yyyy HH:mm");    
		Date resultdate = new Date(yourmilliseconds);
		System.out.println("Current Time: " + sdf.format(resultdate));
		
		System.out.println("Total Time: " + allTotal);
		System.out.println("Step 1: " + oneTotal);
		System.out.println("Step 2: " + twoTotal);
		System.out.println("Step 3: " + threeTotal);
		System.out.println("Step 4: " + afterTotal);
		
		
	}
	
	//The methods for analysis
	public void analyze(String[] curLineSplit, String curLine) throws NumberFormatException, FileNotFoundException, IOException {
		
		step1 = System.currentTimeMillis();
		if (selections[0]) { //fixes indels in VarSifter
			
			curLineSplit = vsCleaner(curLineSplit);
			
		}
		step2 = System.currentTimeMillis();
		
		if (selections[2]) { //annotates
			
			curLineSplit = annotater(curLineSplit, vcfFields(curLineSplit));
			
		}
		
		step3 = System.currentTimeMillis();
		
		
		if (selections[4]) { //deletes columns
			curLineSplit = deleter(curLineSplit);			
			//curLine = "";
			StringBuffer temp = new StringBuffer();
			for (int i = 0; i < curLineSplit.length-1; i++) { //reforms curLine
				//curLine = curLine.concat(curLineSplit[i]).concat("\t");
				temp.append(curLineSplit[i]).append("\t");
			}
			//curLine.concat(curLineSplit[curLineSplit.length-1]);
			temp.append(curLineSplit[curLineSplit.length-1]);
			curLineSplit = temp.toString().split("\t");
		}
		
		//uses StringBuffers to more efficiently concatenate with strings (and is also thread safe)
		//output = "";
		StringBuffer temp = new StringBuffer();
		for (int i=0; i < curLineSplit.length-1; i++) {						//for loop concatenates all but last field, with field followed by tab
			//output = output.concat(curLineSplit[i]).concat("\t");
			temp.append(curLineSplit[i]).append("\t");
		}
		//output = output.concat(curLineSplit[curLineSplit.length-1]);		//concatenates final field, with no tab following it
		temp.append(curLineSplit[curLineSplit.length-1]);
		
		//index adder
		if (selections[6]) {
			//output = Long.toString(newIndex) + "\t" + output;
			temp = new StringBuffer(Long.toString(newIndex) + "\t" + temp.toString());
			newIndex++;
		}
			
		righter.println(temp.toString().trim());									//writes line to output file				
		
		//more time
		oneTotal += (step1-before);
		twoTotal += (step2 - step1);
		threeTotal += (step3 - step2);
		
	}
	
	
	//edits indels in both the ref and var locations and then metadata genotypes
	public String[] vsCleaner(String[] curLineSplit) throws NumberFormatException, FileNotFoundException, IOException {
		if ((curLineSplit[mutTypeIndex].equals("INDEL")) && ((curLineSplit[refAlleleIndex].length() > 1) || (curLineSplit[varAlleleIndex].length() > 1))) {						//if variant line is an INDEL
			
			
			String refBase = getRefGenome(genomeDirectory, curLineSplit[chrIndex], Integer.parseInt(curLineSplit[leftFlankIndex]));
			//retrieve reference base with method, getRefGenome
			if (curLineSplit[refAlleleIndex].equals("''")) { 					//if refAllele has double apostrophes, reassign double single quotes as reference base
				refAllele = refBase;								
			}
			else {													//else, just concatenate reference base as prefix to refAllele
				refAllele = refBase.concat(curLineSplit[refAlleleIndex]); 			
			}
			if (curLineSplit[varAlleleIndex].equals("''")) { 					//if varAllele has double apostrophes, reassign double single quotes as reference base
				varAllele = refBase;
			}
			else {	
				varAllele = refBase.concat(curLineSplit[varAlleleIndex]); 		//else, just concatenate reference base as prefix to varAllele
			}
			leftFlank = Integer.toString(Integer.parseInt(curLineSplit[leftFlankIndex]) - 1); //decrement left flank by one (-1)

			curLineSplit[leftFlankIndex] = leftFlank;							//reassign leftFlank to new value
			curLineSplit[refAlleleIndex] = refAllele;							//reassign refAllele with updates
			curLineSplit[varAlleleIndex] = varAllele;							//reassign varAllele with updates
			
			curLineSplit = genotypeLineFixer(curLineSplit, refBase);			//fixes the metadata genotypes
			
			
		}
		
		return curLineSplit;
		
	}
	
	
	//gets the nucleotide from the reference genome as part of the vsCleaner
	public String getRefGenome(String directory, String chromosome, int leftFlank) throws FileNotFoundException, IOException {
		/* Given a chromosome and position, getRefGenome reads in a FASTA file
		 * corresponding to the given chromosome, retrieves the base at the specified
		 * position within the given chromosome, and returns that base  
		 * 		  
		 * Parameters:
		 * 		chromosome - the abbreviation "chr" concatenated with a chrosomome number, or X, or Y as a string
		 * 		leftFlank - the leftFlank position of the variant as it appears in the VarSifter file
		 * 
		 * Returns:
		 * 		the reference genome base of a specified chromosome and position 
		 * 
		 */
		
		if (!chromosome.equals(currentChrom)) {
			if (!currentChrom.equals("")) {
				refGenome.close();										//if this is not the first file scanner, close previous scanner
			}
			
			File file; 
			try {
				file = new File(directory + "\\" + chromosome + ".fa");
				refGenome = new BufferedReader(new FileReader(file));
			} catch (FileNotFoundException error) {
				String errorMessage = "Your reference genome was missing " + chromosome + ". Please add this file and restart.";
				JOptionPane.showMessageDialog(new Frame("Error Message"), errorMessage);
				System.exit(0);
			}
			
			refGenome.readLine();											//skips over header line
			chromLine = refGenome.readLine();									//grabs first line of reference bases
			lineLength = chromLine.length();     							//assigns LINELENGTH to the number of characters in the FASTA file lines (not the header line)
			currentChrom = chromosome;									//reassigns chromosome value
			lineNumber = 1;
			
		}
				
		int numbOfLines = (leftFlank-1)/lineLength - lineNumber + 1;
		lineNumber += numbOfLines;
		
		for (int i = 0; i < numbOfLines; i++) { 					//for loops over the lines of reference bases
			chromLine = refGenome.readLine();
		}
		
		int positionInLine = (leftFlank-1) % lineLength;   					//determines the position within a given line
		
		return chromLine.substring(positionInLine, positionInLine+1).toUpperCase(); 		//indexes into line to retrieve and then return the reference base
	}

	//adjusts the genotypes in the metadata during vscleanup
	public String[] genotypeLineFixer(String[] splitLine, String refBase) {
		String[] genotypeSplit;
		
		 // Not sure how this for loop works but it somehow fixes an error...
		// now uses StringBuffers() for more memory-efficient concatenation of strings
		StringBuffer temp = new StringBuffer();
		
		//String temp = "";
		
		for (int i = 0; i < splitLine.length - 1; i++) {
			
			temp.append(splitLine[i]).append("\t");
			//temp = temp.concat(splitLine[i]).concat("\t");
		}
		
		//temp = temp.concat(splitLine[splitLine.length-1]);
		temp.append(splitLine[splitLine.length-1]);
		
		splitLine = temp.toString().split("\t");
		
		
		//StringBuffer ref = new StringBuffer(refBase);
		//tried using StringBuffers here, but maybe not necessary or worth the time to code
		for (int i=naPos; i < splitLine.length; i+= naDistanceApart) {
			
			splitLine[i] = splitLine[i].replace("''", "");
			
			if ((!splitLine[i].equals("NA")) && (splitLine[i].contains(":"))) {
				
				if (splitLine[i].indexOf(":") == 0) {
					
					splitLine[i] = refBase.concat(":").concat(refBase).concat(splitLine[i].substring(1));
					//splitLine[i] = ref.append(":").append(refBase).append(splitLine[i].substring(1)).toString();
					
				} else if (splitLine[i].indexOf(":") == splitLine[i].length()-1) {
					
					splitLine[i] = refBase.concat(splitLine[i]).concat(refBase);
					//splitLine[i] = ref.append(splitLine[i]).append(refBase).toString();
					
				} else {
					
					genotypeSplit = splitLine[i].split(":");
					
					splitLine[i] = refBase.concat(genotypeSplit[0]).concat(":").concat(refBase).concat(genotypeSplit[1]);
					//splitLine[i] = ref.append(genotypeSplit[0]).append(":").append(refBase).append(genotypeSplit[1]).toString();
					
				}
			} else if (splitLine[i].equals("")) {
				splitLine[i] = refBase;
			} else if (!splitLine[i].equals("NA")) {
				
				splitLine[i] = refBase.concat(splitLine[i]);
				//splitLine[i] = ref.append(splitLine[i]).toString();
			}
			//else do nothing - leave it alone
		}
		
		return splitLine;
	}
	
	
	//adds the CADD columns in the appropriate place if called
	public String[] annotater(String[] lineSplit, String newFields) {
		
		//changed concatenation here to make use of StringBuffers
		StringBuffer newLine = new StringBuffer();

		//String newLine = "";
		
		for (int i=0; i < vcfPos; i++) {						//for loop concatenates all but last field, with field followed by tab			
			
			//newLine = newLine.concat(lineSplit[i]).concat("\t");
			
			newLine.append(lineSplit[i]).append("\t");
		}
		
		//newLine = newLine.concat(newFields);
		
		newLine.append(newFields);
		
		for (int j=vcfPos; j < lineSplit.length-1; j++) {
			
			//newLine = newLine.concat(lineSplit[j]).concat("\t");
			
			newLine.append(lineSplit[j]).append("\t");
		}
		
		//newLine = newLine.concat(lineSplit[lineSplit.length-1]);
		
		newLine.append(lineSplit[lineSplit.length-1]);
		
		//String newString = newLine.toString();
		
		return newLine.toString().split("\t");					//prints the header line of the VarSifter file
		
	}
	

	//adds the CADD columns to each line
	public String vcfFields(String[] lineSplit) throws NumberFormatException, FileNotFoundException, IOException {

		String refNew = lineSplit[refAlleleIndex];
		String varNew = lineSplit[varAlleleIndex];
		int position = Integer.parseInt(lineSplit[leftFlankIndex]) + 1;
		
		while (refNew.length()>1 && varNew.length()>1) {
			if (refNew.substring(refNew.length()-1).equals(varNew.substring(varNew.length()-1))) {
				refNew = refNew.substring(0, refNew.length()-1);
				varNew = varNew.substring(0, varNew.length()-1);
			}
			
			else {
				break;
			}
		} //eating in
		
		while (refNew.length()>1 && varNew.length()>1) {
			if (refNew.substring(0, 1).equals(varNew.substring(0,1))) {
				refNew = refNew.substring(1);
				varNew = varNew.substring(1);
				position++;
			}
			else {
				break;
			}
		} //making new genomes for CADD
		
		return Integer.toString(position) + "\t" + refNew + "\t" + varNew + "\t";

	}
	
	
	//deletes columns based on the indices in the toDelete matrix
	public String[] deleter(String[] curLineSplit) {
		ArrayList<String> temp = new ArrayList<String>(Arrays.asList(curLineSplit));
		for (int j = toDelete.size() - 1; j >= 0; j--) {
			temp.remove((int) toDelete.get(j));
		}
		curLineSplit = new String[temp.size()];
		curLineSplit = temp.toArray(curLineSplit);
		return curLineSplit; 
	}
	
	
	
	
	
	
}
