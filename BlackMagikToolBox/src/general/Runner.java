package general;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JFrame;
import javax.swing.JOptionPane;

import broadLevelBamFileCuration.AutoBamCuration;
import cncFilter.CNCFilter_GUI;
import confettiFilter.PPDNFilter;
import salvagePathway.SalvageGUI;
import variantExclusionFilter.MakeBamROC_Configs;
import ethnicityMatcher.EthnicityCentralControl;
import kaylaKode.KaylaKodeNew_GUI2;;

public class Runner {

	public static void main(String[] args) throws IOException, InterruptedException {
		
		//This is the entire engine for the User interface System of the Black Magik Tool Box Program
		//The Constructor Automatically calls upon the rest of the user interface system, and thus all the UI is done right from the start.
		UserInterface user = new UserInterface();
		
		
		//This is the entire chain of variables extracted from the UI through the use of retrieval methods and piped into the various modules afterwards.
		//General Variables
		File vsFile = user.getVSFile(); //Varsifter File
		String dest = user.getDest(); //Output File Name
		File bamDir = user.getBamDir(); //Bam File Directory
		String configPaths = user.getConfigPath();
		
		//Module Specific Variables
		
		//Ethnicity Matcher
		String cutoffCriteria = user.getCriteria();
		String[] markers = user.getMarkers();
		
		//Salvage Pathway
		int minorField = user.getMinor();
		int yesPop = user.getYesPop();
		int [] popHeaders = user.getPopIndexes();
		
		//KaylaKode
		ArrayList<Integer> family = user.getFam();
		ArrayList<Boolean> sibsAffected = user.getSibsAffected();
		ArrayList<Boolean> sibGen = user.getSibsGender();
		boolean proGend = user.getProbandGend();
		boolean[] popBools = user.getPopBools();
		ArrayList<int[]> populations = user.getPopulations();
		ArrayList<int[]> xFilterS = user.getXFilterS();
		ArrayList<int[]> xFilterC = user.getXFilterC();
		ArrayList<Integer> xFilterExt = user.getXFilterExt();
		int extData = user.getExtData();
		double rejectionThreshold = user.getRejectionThreshold();
		int threads = user.getThreads();
		int cutCnt2 = user.getCutCnt2();
		String strengthPath = user.getStrengthPath();
		String genePath = user.getGenePath();
		
		//Variant Exclusion Filter
		Double [] caddThreshold = user.getCaddThresh();
		int homVarThresh = user.getHomVarThreshold();
		int hetThresh = user.getHetThreshold();
		int hemVarThresh = user.getHemVarThreshold();
		File ROCHeaderConfig = user.getROCHeaderConfig();
		ArrayList<Integer> genoThreshList = user.getGenoThreshList();
		
		//Confetti Filter
		int homRefIndex = user.getHomRefIndex();
		int homVarIndex = user.getHomVarIndex();
		int genoIndex = user.getGenotypeIndex();
		int genotypeThresh = user.getGenotypeThreshold();
		int baseFilter = user.getBaseFilter();
		int baseQualThresh = user.getBaseQualThresh();
		double badReadRatio = user.getBadReadRatio();
		
		//CNC Filter
		int clusterSize = user.getClusterMin();
		String exonConfig = user.getExonBoundConfig();
		int CNCExon = user.getCNCExon();
		
		//Variables used for the Time Messages
		long bRunner;
		long aRunner;
		long[] bModules = new long[7];
		long[] aModules = new long[7];
		String [] timeMessages = new String[7];
		String timeMessage;
		
		//Output Folder
		File outputDir = new File(vsFile.getParent() + "\\" + dest + "_BlackBox_Output");
		outputDir.mkdir();
		
		
		dest = outputDir.getAbsolutePath() + "\\" + dest; //Puts destination in output folder
		
		//Generates Ped and Sample files
		user.generateSampleFile(outputDir.getAbsolutePath());
		user.generatePedFile(outputDir.getAbsolutePath());
		
		File sampleFile = user.getSampleFile(); //Sample IDs
		File pedFile = user.getPedFile(); //Pedigree File
		
		//Gets BAM chr hashmap for translating chr inputs for different BAMs
		HashMap<String, ArrayList<String>> bamChrMap = user.getbamChrMap();
		
		//Generates file with all of user inputs
		user.generateUserFile(outputDir.getAbsolutePath() + "\\User_Inputs.txt");
		
		//Start of timer
		bRunner = System.currentTimeMillis();
		
		//All modules running, with time statistics for before and after each module
		//The output of each module is "piped" to the next one
		//Each module has their respective constructors to take in all the necessary inputs for each program
		
		bModules[0] = System.currentTimeMillis();
		
		EthnicityCentralControl eth = new EthnicityCentralControl(bamDir.getAbsolutePath(), markers, sampleFile, cutoffCriteria, (dest + "_EthnicityMatched"), bamChrMap);
		
		aModules[0] = System.currentTimeMillis();
		
		timeMessages[0] = ("\n\nThe Ethnicity Matcher Took " + (aModules[0] - bModules[0]) + "ms ("
				+ ((aModules[0] - bModules[0]) / 1000) + "." + (((aModules[0] - bModules[0]) % 1000) / 100) + " seconds).");
		
		File ethOutput = eth.getOutput();
		
		
		
		bModules[1] = System.currentTimeMillis();
		
		SalvageGUI salv = new SalvageGUI(vsFile, pedFile, bamDir, ethOutput, minorField, yesPop, popHeaders, (dest + "_SalvagePathwayed"), bamChrMap);
		
		aModules[1] = System.currentTimeMillis();
		
		timeMessages[1] = ("\n\nThe Salvage Pathway Took " + (aModules[1] - bModules[1]) + "ms ("
				+ ((aModules[1] - bModules[1]) / 1000) + "." + (((aModules[1] - bModules[1]) % 1000) / 100) + " seconds).");
		
		File salvOutput = salv.getOutput();
		
		
		
		bModules[2] = System.currentTimeMillis();
		
		KaylaKodeNew_GUI2 kKode = new KaylaKodeNew_GUI2(salvOutput, dest + "_KaylaKoded", family, sibsAffected, sibGen, proGend, popBools, populations, xFilterS, xFilterC, rejectionThreshold, threads, cutCnt2, strengthPath, genePath, minorField, xFilterExt, extData);
		
		aModules[2] = System.currentTimeMillis();
		
		timeMessages[2] = ("\n\nThe Kayla Kode Took " + (aModules[2] - bModules[2]) + "ms ("
				+ ((aModules[2] - bModules[2]) / 1000) + "." + (((aModules[2] - bModules[2]) % 1000) / 100) + " seconds).");
		
		File kOutput = kKode.getOutput();
		
		
		bModules[3] = System.currentTimeMillis();
		
		AutoBamCuration zeCurator = new AutoBamCuration(kOutput, pedFile, bamDir, (dest + "_AutoBamCurated"), minorField, bamChrMap);

		aModules[3] = System.currentTimeMillis();
		
		timeMessages[3] = ("\n\nThe Broad Level BAM Curator/ SNR Calculator Took " + (aModules[3] - bModules[3]) + "ms ("
				+ ((aModules[3] - bModules[3]) / 1000) + "." + (((aModules[3] - bModules[3]) % 1000) / 100) + " seconds).");
		
		File curatorOutput = zeCurator.getOutput();
		
		
		bModules[4] = System.currentTimeMillis();
		MakeBamROC_Configs varExclusFilter = new MakeBamROC_Configs(curatorOutput, (dest + "_VariantExcluded"), pedFile, caddThreshold, homVarThresh, hetThresh, hemVarThresh, ROCHeaderConfig, genoThreshList, sibsAffected, sibGen);

		aModules[4] = System.currentTimeMillis();
		
		timeMessages[4] = ("\n\nThe Variant Exclusion Filter Took " + (aModules[4] - bModules[4]) + "ms ("
				+ ((aModules[4] - bModules[4]) / 1000) + "." + (((aModules[4] - bModules[4]) % 1000) / 100) + " seconds).");
		
		File varExclusOutput = varExclusFilter.getOutput();
		
		
		
		bModules[5] = System.currentTimeMillis();
		PPDNFilter confetti = new PPDNFilter(varExclusOutput, (dest+ "_Confettied"), pedFile, bamDir, minorField, homRefIndex, homVarIndex, genoIndex, bamChrMap, genotypeThresh, baseFilter, baseQualThresh, badReadRatio);
		
		aModules[5] = System.currentTimeMillis();
		
		timeMessages[5] = ("\n\nThe Confetti Filter Took " + (aModules[5] - bModules[5]) + "ms ("
				+ ((aModules[5] - bModules[5]) / 1000) + "." + (((aModules[5] - bModules[5]) % 1000) / 100) + " seconds).");
		
		File confettiOutput = confetti.getOutput();
		
		
		
		bModules[6] = System.currentTimeMillis();
		CNCFilter_GUI cncer = new CNCFilter_GUI(confettiOutput, dest, pedFile, bamDir, clusterSize, bamChrMap, exonConfig, CNCExon);
		
		aModules[6] = System.currentTimeMillis();
		
		timeMessages[6] = ("\n\nThe CNC Filter Took " + (aModules[6] - bModules[6]) + "ms ("
				+ ((aModules[6] - bModules[6]) / 1000) + "." + (((aModules[6] - bModules[6]) % 1000) / 100) + " seconds).");
		
		File finalOutput = cncer.getOutput();
		
		aRunner = System.currentTimeMillis();
		
		//Final Message Screen
		timeMessage = ("The BlackBox Program Took " + ((aRunner - bRunner) / (1000*60)) + "." + (((aRunner - bRunner) % (1000*60)) / 100) + " minutes or " + (aRunner - bRunner) + "ms ("
				+ ((aRunner - bRunner) / 1000) + "." + (((aRunner - bRunner) % 1000) / 100) + " seconds).");
		
		String destinationMessage = ("Thank you for using the Black Magik Tool Box Program!.\n\nThe Output File is located at the path found below:\n" + finalOutput.getAbsolutePath() + "\n\nTime Statistics:\n\n");
		
		
		JFrame finalFrame = new JFrame("Completed");
		
		finalFrame.setUndecorated(true);
		finalFrame.setVisible(true);
		finalFrame.setLocationRelativeTo(null);
		finalFrame.setAlwaysOnTop(true);
		
		//Displays final output screen, puts time statistics information in output folder and exits program.
		JOptionPane.showMessageDialog(finalFrame, ("Completed.\n\n") + destinationMessage +  timeMessage + timeMessages[0] + timeMessages[1] + timeMessages[2] + timeMessages[3] + timeMessages[4] + timeMessages[5] + timeMessages[6]);
		
		File timeFile = new File(outputDir.getAbsolutePath() + "\\" + "Time_Statistics.txt");
		
		PrintWriter timeWalker = new PrintWriter(timeFile);
		
		timeWalker.println("Time Statistics");
		
		timeWalker.println();
		
		timeWalker.println(timeMessage.replaceAll("\n", ""));
		
		for(int i = 0; i < timeMessages.length; i++) {
			
			timeWalker.println();
			
			timeWalker.println(timeMessages[i].replaceAll("\n", ""));
			
		}
		
		timeWalker.close();
		
		System.gc();
		System.exit(0);
		
	}
	
	
}
