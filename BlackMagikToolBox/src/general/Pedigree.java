package general;

//take pedigree input
//or take in a config file

import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.List;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

public class Pedigree {
	
	/*
	 * Returns list denoting familial relationships based on string containing
	 * pedigree info.
	 */
	public static ArrayList<Integer> getPedigreefromString(String[] family, ArrayList<String> headers) {
		ArrayList<Integer> list = new ArrayList<Integer>();
		for (int i = 0; i < family.length; i++) {
			int indice = headers.indexOf(family[i] + ".NA");
			if (indice != -1) {
				list.add(headers.indexOf(family[i] + ".NA"));
			} else {
				JOptionPane.showMessageDialog(new Frame("Input prompt"), "Person not found.");
				System.exit(1);
			}
		}
		return list;

	}
	
}
