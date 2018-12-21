import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import java.util.*;
public class PolySort {
	
	private File unsorted;
	private File sorted;
	private BufferedReader reder;
	private PrintWriter righter; //"Rights" to the final file
	private String newFileName;
	
	private int chrIndex;
	private int posIndex;
	private int refIndex;
	private int altIndex;
	private int probIndex;
	
	private String header;
	
	
	private long bRunner;
	private long aRunner;
	private String timeMessage;
	
	String line;
	String line2;
	
	
	public PolySort(File unsorted, int chrIndex, int posIndex, int refIndex, int altIndex, int probIndex, String outputName) throws IOException {
		
		this.unsorted = unsorted;
		this.chrIndex = chrIndex;
		this.posIndex = posIndex;
		this.refIndex = refIndex;
		this.altIndex = altIndex;
		this.probIndex = probIndex;
		this.newFileName = outputName;
		
		sort();
	}
	
	public File getOutput() {
		
		return sorted;
		
	}
	
	public String getTimeMessage() {
		
		return timeMessage;
		
	}
	
	public void sort() throws IOException {
		
		sorted = new File(newFileName);
		reder = new BufferedReader(new FileReader(unsorted));
		righter = new PrintWriter(sorted);
		
		
		//UI canceler
		JFrame jframe = canceler();
		jframe.setLocationRelativeTo(null);
		jframe.setVisible(true);
		
		
		bRunner = System.currentTimeMillis();
		
		//Headers
		header = reder.readLine();
		righter.println(header);
		
		line = reder.readLine();
		line2 = reder.readLine();
		
		while(true) {
			
			if(duplicate(line, line2)) {
				
				probSort(line, line2);
				
			} else {
				
				righter.println(line);
				
			}
			
			line = line2;
			line2 = reder.readLine();
			
			if(line == null) {
				
				break;
				
			}
			else if(line2 == null) {
				
				righter.println(line);
				break;
			}
			
		}
		
		
		aRunner = System.currentTimeMillis();
		
		jframe.dispose();
		
		timeMessage = "\nThe Sorting Program Took " + (aRunner - bRunner) + "ms ("
				+ ((aRunner - bRunner) / 1000) + "." + (((aRunner - bRunner) % 1000) / 100) + " seconds)";
		
		righter.close();
		reder.close();
		
	}
	
	
	public JFrame canceler() {
		JFrame jframe = new JFrame();
		GridBagLayout gbl = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		jframe.setLayout(gbl);
		jframe.setTitle("Line-by-Line Utilities");
		jframe.setSize(500, 100);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		JPanel jpanel = new JPanel();
		jpanel.setLayout(gbl);
		
		JLabel jtext = new JLabel("Sorting Program running. If you'd like to abort, hit the cancel button below.", SwingConstants.CENTER);
		gbc.gridx = 0;
		gbc.gridy = 0;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(5, 5, 5, 5);
		gbc.weightx = 1.0;
		jpanel.add(jtext, gbc);
		
		JButton abort = new JButton("Cancel");
		abort.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				jframe.dispose();
				System.gc();
				
				try {
					
					//IO Closing
					righter.close();
					reder.close();
					
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
				
				
				System.exit(0);
			}
		});
		
		gbc.gridx = 0;
		gbc.gridy = 21;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.insets = new Insets(5, 5, 5, 5);
		gbc.weightx = 1.0;
		jpanel.add(abort, gbc);
		
		jframe.getContentPane().add(jpanel);
		return jframe; 
	}
	
	private boolean duplicate(String liner, String liner2) {
		
		String[] linez1 = liner.split("\t");
		String[] linez2 = liner2.split("\t");
		if((linez1[chrIndex].equals(linez2[chrIndex])) && (linez1[posIndex].equals(linez2[posIndex])) && (linez1[refIndex].equals(linez2[refIndex])) && (linez1[altIndex].equals(linez2[altIndex]))) {
			
			
			return true;
			
		} else {
			
			return false;
			
		}
	}
	
	private void probSort(String liner, String liner2) throws IOException {
		
		String[] linerArr = liner.split("\t");
		String[] liner2Arr = liner2.split("\t");
		
		ArrayList<String []> linez = new ArrayList<String[]>();
		
		linez.add(linerArr);
		linez.add(liner2Arr);
		
		String liner3 = reder.readLine();
		if(liner3 != null) {
			while(duplicate(liner3, liner2)) {
				
				String[] liner3Arr = liner3.split("\t");
				linez.add(liner3Arr);
				
				liner3 = reder.readLine();
				
				if(liner3 == null) {
					
					break;
					
				}
			}
		}
		line2 = liner3;
		
		ArrayList<String[]> sortList = mergeSort(linez, probIndex);
		
		
		/* For Non Deduping
		for(int i = 0; i < sortList.size(); i++) {
			
			String [] sortLineArr = sortList.get(i);
			
			String sortLine = "";
			
			for(int j = 0; j < sortLineArr.length - 1; j++) {
				
				sortLine += sortLineArr[j] + "\t";
				
			}
			
			sortLine += sortLineArr[sortLineArr.length-1];
			
			righter.println(sortLine);
			
		}
		
		*/
		
		String [] sortLineArr = sortList.get(0);
		
		String sortLine = "";
		
		for(int j = 0; j < sortLineArr.length - 1; j++) {
			
			sortLine += sortLineArr[j] + "\t";
			
		}
		
		sortLine += sortLineArr[sortLineArr.length-1];
		
		righter.println(sortLine);
	}
	
	private ArrayList<String[]> mergeSort(ArrayList<String[]> linezArr, int probIndex) {
		
		ArrayList<String[]> left = new ArrayList<String[]>();
		ArrayList<String[]> right = new ArrayList<String[]>();
		
		if(linezArr.size() <= 1) {
			
			return linezArr;
			
		} else {
			
			int middle = linezArr.size()/2;
			for(int i = 0; i < middle; i++) {
				
				left.add(linezArr.get(i));
				
			}
			
			for(int i = middle; i < linezArr.size(); i++) {
				
				right.add(linezArr.get(i));
				
			}
			
			left = mergeSort(left, probIndex);
			right = mergeSort(right, probIndex);
			
			return merge(left, right, probIndex);
			
			
		}
		
	}
	
	private static ArrayList<String[]> merge(ArrayList<String[]> left, ArrayList<String[]> right, int probIndex) {
		
		ArrayList<String[]> result = new ArrayList<String[]>();
		
		while(left.size() > 0 && right.size() > 0) {
			
			if(Double.parseDouble(left.get(0)[probIndex]) > Double.parseDouble(right.get(0)[probIndex])) {
				
				result.add(left.get(0));
				left.remove(0);
				
				
			} else if(Double.parseDouble(left.get(0)[probIndex]) <= Double.parseDouble(right.get(0)[probIndex])) {
				
				result.add(right.get(0));
				right.remove(0);
				
			}
			
			
		}
		
		if(left.size() > 0) {
			
			for(int i = 0; i < left.size(); i++) {
				
				result.add(left.get(i));
				
			}
			
		}
		
		if(right.size() > 0) {
			
			for(int i = 0; i < right.size(); i++) {
				
				result.add(right.get(i));
				
			}
			
		}
		
		return result;
		
	}
	
	
	
	
}
