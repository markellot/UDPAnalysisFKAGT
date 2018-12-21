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
public class Sorter {
	
	private File unsorted;
	private File sorted;
	private BufferedReader reder;
	private PrintWriter righter; //"Rights" to the final file
	private PrintWriter wronger; //"Rights" to each temp file
	
	private int chrIndex;
	private int posIndex;
	
	private String header;
	
	//temp files
	private File dir;
	private List<File> chunkies = new ArrayList<File>();
	
	private final int LINESPERFILE = 10000;
	
	private long bRunner;
	private long aRunner;
	private String timeMessage;
	
	public Sorter(File unsorted, int chrIndex, int posIndex) throws IOException {
		
		this.unsorted = unsorted;
		this.chrIndex = chrIndex;
		this.posIndex = posIndex;
		
		sort();
	}
	
	public File getOutput() {
		
		return sorted;
		
	}
	
	public void sort() throws IOException {
		
		sorted = new File(unsorted.getAbsolutePath().substring(0, unsorted.getAbsolutePath().lastIndexOf(".")) + "_sorted.txt");
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
		
		ArrayList<String[]> linez = new ArrayList<String[]>();
		
		//Directory of where temp files go
		
		int suffix = unsorted.getName().lastIndexOf(".");
		if(suffix == -1) {
			
			suffix = unsorted.getName().length();
			
		}
		dir = new File(unsorted.getParent() + "\\SortingTempFiles" + unsorted.getName().substring(0, suffix));
		dir.mkdir();
		
		int numFiles = 0;
		
		String line = "";
		String[] lineArr;
		
		while(true) {
			
			//Creating new temp files
			chunkies.add(File.createTempFile("Part " + Integer.toString(numFiles+1), ".txt", dir));
			chunkies.get(numFiles).deleteOnExit();
			wronger = new PrintWriter(chunkies.get(numFiles));
			
			
			for(int i = 0; i < LINESPERFILE; i++) {
				
				line = reder.readLine();
				
				if(line == null) {
					
					lineArr = null;
					break;
					
				}
				
				lineArr = line.split("\t");
				linez.add(lineArr);
				
			}
						
			linez = mergeSort(linez, chrIndex, posIndex);
						
			for(int i = 0; i < linez.size(); i++) {
				
				String output = "";
				
				for(int j = 0; j < linez.get(i).length-1; j++) {
					
					output = output.concat(linez.get(i)[j]).concat("\t");
					
				}
				
				output = output.concat(linez.get(i)[linez.get(i).length-1]);
				
				wronger.println(output);
				
			}
			
			wronger.close();
			numFiles++;
			linez.clear();
			
			if(line==null) {
				
				break;
				
			}
			
		}
		
		mergeFiles(chunkies, numFiles, chrIndex, posIndex);
		
		aRunner = System.currentTimeMillis();
		
		jframe.dispose();
		
		
		timeMessage = ("The Sorting Program Took " + (aRunner - bRunner) + "ms ("
		+ ((aRunner - bRunner) / 1000) + "." + (((aRunner - bRunner) % 1000) / 100) + " seconds)\n\n");	
		
		righter.close();
		wronger.close();
		reder.close();
		

		for(int i = 0; i < chunkies.size(); i++) {
			
			chunkies.get(i).delete();
			
		}
		dir.delete();
		
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
					wronger.close();
					reder.close();
					
					for(int i = 0; i < chunkies.size(); i++) {
						
						chunkies.get(i).delete();
						
					}
					
					dir.delete();
					
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
	
	
	private ArrayList<String[]> mergeSort(ArrayList<String[]> linezArr, int chrIndex, int posIndex) {
		
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
			
			left = mergeSort(left, chrIndex, posIndex);
			right = mergeSort(right, chrIndex, posIndex);
			
			return merge(left, right, chrIndex, posIndex);
			
			
		}
		
	}
	
	private static ArrayList<String[]> merge(ArrayList<String[]> left, ArrayList<String[]> right, int chrIndex, int posIndex) {
		
		ArrayList<String[]> result = new ArrayList<String[]>();
		
		while(left.size() > 0 && right.size() > 0) {
			
			if(left.get(0)[chrIndex].compareTo(right.get(0)[chrIndex]) < 0) {
				
				result.add(left.get(0));
				left.remove(0);
				
				
			} else if(left.get(0)[chrIndex].compareTo(right.get(0)[chrIndex]) > 0) {
				
				result.add(right.get(0));
				right.remove(0);
				
			} else {
				
				if(Integer.parseInt(left.get(0)[posIndex]) <= Integer.parseInt(right.get(0)[posIndex])) { //Should we have the <= here? or <
					
					result.add(left.get(0));
					left.remove(0);
					
				} else {
					
					result.add(right.get(0));
					right.remove(0);
					
				}
				
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
	
	
	
	private void mergeFiles(List<File> chunkies, int numFiles, int chrIndex, int posIndex) throws IOException {
		
		ArrayList<BufferedReader> rederz = new ArrayList<BufferedReader>();
		
		ArrayList<String[]> filerows = new ArrayList<String[]>();
		
		boolean hasRowsLeft = false;
		
		for(int i = 0; i < numFiles; i++) {
			
			rederz.add(new BufferedReader(new FileReader(chunkies.get(i))));
			
			String line;
			if((line = rederz.get(i).readLine()) != null) {
				
				filerows.add(line.split("\t"));
				hasRowsLeft = true;
				
			}
			else {
				
				filerows.add(null);
				
			}
			
		}
		
		String[] lineArr;
		while(hasRowsLeft) {
			
			String chr;
			String pos;
			int minIndex;
			
			lineArr = filerows.get(0);
			
			if(lineArr != null) {
				
				chr = lineArr[chrIndex];
				pos = lineArr[posIndex];
				minIndex = 0;
				
			}
			else {
				
				chr = null;
				pos = null;
				minIndex = -1;
				
			}
			
			for(int i = 1; i < filerows.size(); i++) {
				
				lineArr = filerows.get(i);
				if(chr != null && pos != null) {
					
					if(lineArr != null) {
						
						if(lineArr[chrIndex].compareTo(chr) < 0) {
							
							minIndex = i;
							chr = lineArr[chrIndex];
							pos = lineArr[posIndex];
							
						}
						else if(lineArr[chrIndex].compareTo(chr) == 0) {
							
							if(Integer.parseInt(lineArr[posIndex]) <= Integer.parseInt(pos)) {
								
								minIndex = i;
								
								chr = lineArr[chrIndex];
								pos = lineArr[posIndex];
								
							}
							
						}
						
					}
					
				} else {
					
					
					if(lineArr != null) {
						
						minIndex = i;
						chr = lineArr[chrIndex];
						pos = lineArr[posIndex];
						
						
					}
					
					
				}
				
			}
			
			if(minIndex < 0) {
				
				hasRowsLeft = false;
				
			} else {
				
				String output = "";
				
				for(int i = 0; i < filerows.get(minIndex).length-1; i++) {
					
					output = output.concat(filerows.get(minIndex)[i]).concat("\t");
					
				}
				
				output = output.concat(filerows.get(minIndex)[filerows.get(minIndex).length-1]);
				
				
				righter.println(output);
				
				String line = rederz.get(minIndex).readLine();
				
				if(line != null) {
					
					filerows.set(minIndex, line.split("\t"));
					
				}
				
				else {
					
					filerows.set(minIndex, null);
					
				}
			}
			
			for(int i = 0; i < filerows.size(); i++) {
				
				hasRowsLeft = false;
				if(filerows.get(i) != null) {
					
					if(minIndex<0) {
						
						JOptionPane.showMessageDialog(new Frame("ERROR"), "Minimum Index for lowest value is less than 0 and program found non-null - System exitting...");
						System.exit(-1);
						
					}
					
					hasRowsLeft = true;
					break;
					
				}

			}
			
			if(!hasRowsLeft) {
				
				for(int i = 0; i < filerows.size(); i++) {
					
					if(filerows.get(i) == null) {
						
						String line = rederz.get(i).readLine();
						
						if(line != null) {
							
							hasRowsLeft = true;
							filerows.set(i, line.split("\t"));
							
						}
						
					}
					
				}
				
			}
			
		}
		
		for(int i = 0; i < rederz.size(); i++) {
			
			rederz.get(i).close();
			
		}
		
		
		
	}
	
	public String getTimeMessage() {
		
		return timeMessage;
		
	}
	
	
	
	
	
}
