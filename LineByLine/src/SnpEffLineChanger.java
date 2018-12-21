import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

public class SnpEffLineChanger extends Thread{
	private BufferedReader red;
	private PrintWriter righter;
	
	private int naPos;
	private int naDistanceApart;
	
	private int chrIndex;
	private int refIndex;
	private int varIndex;
	private int mutTypeIndex;
	private int afIndex;
	private int totIndex;
	private int acIndex;
	
	private List<String> oldGl;
	private List<String> newGl;
	
	private long before;
	private long after;
	private long allTotal;
	private int FileNum;
	
	public SnpEffLineChanger(BufferedReader red, PrintWriter righter, int naPos, int naDistApart, int chrIndex, int refIndex, int varIndex,
			int mutTypeIndex, int afIndex, int anIndex, int acIndex, List<String> oldGl, List<String> newGl, int FileNum) {
		
		this.red = red;
		this.righter = righter;
		this.naPos = naPos;
		this.naDistanceApart = naDistApart;
		this.chrIndex = chrIndex;
		this.refIndex = refIndex;
		this.varIndex = varIndex;
		this.mutTypeIndex = mutTypeIndex;
		this.afIndex = afIndex;
		this.acIndex = acIndex;
		this.totIndex = anIndex;
		this.oldGl = oldGl;
		this.newGl = newGl;
		this.FileNum = FileNum;
	}
	
	@Override
	public void run() {

		String read;
		String[] readLineSplit;
		try {
			//while loop reading through each chromosome 
			while((read = red.readLine()) != null) {
				
				before = System.currentTimeMillis();
				
				readLineSplit = read.split("\t");
				
				
				readLineSplit = lineChanger(readLineSplit);       //brute force of each line and passes each line into analysis method
				
				//time
				after = System.currentTimeMillis();
				allTotal += (after-before);
				
				//uses StringBuffers to more efficiently concatenate with strings (and is also thread safe)
				//output = "";
				StringBuffer temp = new StringBuffer();
				for (int i=0; i < readLineSplit.length-1; i++) {						//for loop concatenates all but last field, with field followed by tab
					//output = output.concat(readLineSplit[i]).concat("\t");
					temp.append(readLineSplit[i]).append("\t");
				}
				//output = output.concat(readLineSplit[readLineSplit.length-1]);		//concatenates final field, with no tab following it
				temp.append(readLineSplit[readLineSplit.length-1]);
				
				
				righter.println(temp.toString().trim());
				
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//thread statistics
		System.out.println("Thread Completed for Part: " + FileNum);
		
		long yourmilliseconds = System.currentTimeMillis();
		SimpleDateFormat sdf = new SimpleDateFormat("MMM dd,yyyy HH:mm");    
		Date resultdate = new Date(yourmilliseconds);
		System.out.println("Current Time: " + sdf.format(resultdate));
		
		System.out.println("Total Time: " + allTotal);
		
		
	}
	
	
	public String[] lineChanger(String[] linezArr) {

		//uses StringBuffers to more efficiently concatenate with strings (and is also thread safe)
		//output = "";
		StringBuffer temp = new StringBuffer();
		for (int i=0; i < linezArr.length-1; i++) {						//for loop concatenates all but last field, with field followed by tab
			//output = output.concat(readLineSplit[i]).concat("\t");
			temp.append(linezArr[i]).append("\t");
		}
		//output = output.concat(readLineSplit[readLineSplit.length-1]);		//concatenates final field, with no tab following it
		temp.append(linezArr[linezArr.length-1]);
		
		
		linezArr = temp.toString().trim().split("\t");
		
		for(int i = naPos; i < linezArr.length; i += naDistanceApart) {
						
			if(linezArr[i].equals("**")) {
				
				linezArr[i] = "'':''";
				
			}
			if(linezArr[i].contains("*A")) {
				
				linezArr[i] = linezArr[i].replace("*A", "'':A");
				
			}
			if(linezArr[i].contains("*C")) {
				
				linezArr[i] = linezArr[i].replace("*C", "'':C");
				
			}
			if(linezArr[i].contains("*G")) {
				
				linezArr[i] = linezArr[i].replace("*G", "'':G");
				
			}
			if(linezArr[i].contains("*T")) {
				
				linezArr[i] = linezArr[i].replace("*T", "'':T");
				
			} if(linezArr[i].contains("*:")) {
				
				linezArr[i] = linezArr[i].replace("*:", "'':");
				
			} if(linezArr[i].contains(":*")) {
				
				linezArr[i] = linezArr[i].replace(":*", ":''");
				
			}
			
		}
		
		if(linezArr[refIndex].length() > 1) {
			
			linezArr[mutTypeIndex] = "INDEL";
			
		}
		
		if(linezArr[varIndex].contains("*")) {
			
			linezArr[mutTypeIndex] = "INDEL";
			linezArr[varIndex] = "''";			
			
		}
		
		if(linezArr[afIndex].equals("NaN")) {
			
			linezArr[afIndex] = Double.toString(((double)Integer.parseInt(linezArr[acIndex])/(double)Integer.parseInt(linezArr[totIndex])));
					
		}
		
		if(linezArr[chrIndex].equals("chrMT")) {
			
			linezArr[chrIndex] = "chrM";
			
		}
		
		if(linezArr[chrIndex].contains("GL")) {
			
			for(int i = 0; i < oldGl.size(); i++) {
				
				if(oldGl.get(i).equals(linezArr[chrIndex])) {
					
					linezArr[chrIndex] = newGl.get(i);
					break;
					
				}
				
				
			}
			
		}
		
		return linezArr;
	}
	
	
}
