package ethnicityMatcher;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import htsjdk.samtools.SamReader;

public class MarkerfilePair {
private static BufferedReader markerfile;

//go through each marker set and compute the log ratio score for the marker set
public static double MarkerScore(String config,  SamReader samReader, int coverage, int score, int scale, PrintWriter report2, ArrayList<String> bamChrHeadz) throws IOException{

	//read in the marker set
	markerfile = new BufferedReader(new FileReader(config));
	String[] firstline=markerfile.readLine().split("\t");
	int markercount=0;

	ArrayList<String> header=new ArrayList<String>();
	
	for(int i=0;i<firstline.length;i++){
		header.add(firstline[i]);
	}
	
	double ethnicityscore=0;
	
	//in the marker set, extract precalculated score for each genotype type
	int chr=header.indexOf("Chr");
	int left=header.indexOf("LeftFlank");
	int right=header.indexOf("RightFlank");
	int ref=header.indexOf("ref_allele");
	int alt=header.indexOf("var_allele");
	int hetfreq=header.indexOf("Het_Score");
	int horfreq=header.indexOf("Hor_Score");
	int hoafreq=header.indexOf("Hoa_Score");
	int varfreq=header.indexOf("1KGen_Allele_freq");
	

	String Line=markerfile.readLine();
	ArrayList<String[]> genotypeline=new ArrayList<String[]>();
	
	while(Line!=null){
		String[] cursplit=Line.split("\t");
		//each different marker set is separated by "Group"
		if(cursplit[0].contains("Group")){
			//calculate the first well-sequenced genotype in the group set
			double tmp=GenotypeGroupScore.variantdensity(genotypeline, samReader, chr, left, right, ref, alt, hetfreq, horfreq, hoafreq, varfreq, coverage, score, bamChrHeadz);
			
			//if there is a well-sequenced genotype in the group set
			if(tmp!=-7777){
				markercount++;
				ethnicityscore+=tmp;
			}
			genotypeline=new ArrayList<String[]>();
		}
		else{
			//put each group together
		genotypeline.add(cursplit);	
		}	
		Line=markerfile.readLine();
	}	
	
	//if there are enough markers hit and the scale is above 0, scale the marker set
	//the scale is below zero when the marker set doesn't need to be scaled
	if(markercount>(scale-15)  && scale>0){
		ethnicityscore=ethnicityscore/(double)markercount*(double)scale;
	}
	
	if(markercount==0){
		ethnicityscore=-7777;
	}
	
	report2.println(config+"\t"+markercount+"\t"+ethnicityscore);
return ethnicityscore;
}
}
