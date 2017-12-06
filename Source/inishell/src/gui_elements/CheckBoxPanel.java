/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/*	This file is part of INIshell.
*
*   INIshell is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   INIshell is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with INIshell.  If not, see <http://www.gnu.org/licenses/>.
*
*/

package gui_elements;

import main.*;

import java.awt.Color;
import java.awt.Font;
import java.awt.event.ItemEvent;
import java.util.HashMap;

import javax.swing.JCheckBox;
import javax.swing.JOptionPane;
import javax.swing.JLabel;

import org.w3c.dom.Element;
import org.w3c.dom.Node;

public class CheckBoxPanel extends OptionPanel {

	private static final long serialVersionUID = 552186457208697362L;
	private final JCheckBox[] boxes;
	String concat;
	private String setKey = null, setValue = null;
	private boolean doSet = false;
	private HashMap setMap;


	/**
	 * Constructor for a panel that can hold a number of checkboxes to compose
	 * the value for a key/value set.
	 *
	 * @param element
	 * @throws GUIBuildException
	 */
	public CheckBoxPanel(Element element, ControlledPanel parent) throws GUIBuildException {
		super(element, parent);

		concat = "";
		boxes = new JCheckBox[values.length];
		for (int i = 0; i < values.length; i++) {
			boxes[i] = new JCheckBox(values[i]);

			this.add(boxes[i], "cell 1 " + i + ", wrap, ay top");

			String extra = "";
			if (i > 0)
				extra = ", gapleft 12";

			this.add(ControlledPanel.createHelpPane(this.helptexts[i]), "cell 6 "+ i +", width 200:400:400, wrap" + extra);

			if (optional) {
				setNewBackground(boxes[i], Color.white);
				//boxes[i].setForeground(GUI.DFLT_VAL_COLOR);
				//boxes[i].setFont( boxes[i].getFont().deriveFont(Font.ITALIC) );
			} else {
				setNewBackground(boxes[i], GUI.WARNING_BACKGROUND_COLOR);
				//boxes[i].setFont( boxes[i].getFont().deriveFont(Font.PLAIN) );
			}
			boxes[i].setVisible(true);
			boxes[i].addItemListener(this);
			if (defaultTrues.contains(values[i])) {
				boxes[i].setSelected(true);
				concat += values[i] + " ";
			}
		}
	}

	@Override
	public void close() {
		for (final JCheckBox box : boxes)
			box.setSelected(false);
	}
	
	private void setNewBackground(JCheckBox box, final Color new_color) {
		box.setOpaque(false);
		box.setBackground(new_color);
		box.setOpaque(true);
	}

	@Override
	public String getValue() {
		if (!isOptional() && concat.equals("")) {
			JOptionPane.showMessageDialog(null, "No option was selected for "
									+ getKey() + " (section '" + getSection().toUpperCase() + "')"
									+ ".\nA value is required. Your .ini-file is probably incorrect.",
									"Problem when building .ini file",
									JOptionPane.WARNING_MESSAGE);
		}

		if (concat.equals("")) return null;
		else	return concat.trim();
	}

	public synchronized void set(HashMap hm, String key, String value) {
		if (hm==null) {
			System.out.println("\thashmap null in CheckBoxPanel. Set canceled for " + key + " value " + value);
			return;
		}
		final String val = (String)hm.get(hashKey);
		if (val==null) {
			System.out.println("**** Key "+hashKey+" not found!");
			return;
		}

		boolean oneBoxSelected = false;
		final String[] tokens = (val.toUpperCase()).split("\\s+");
		for (final JCheckBox box : boxes) {
			box.setSelected(false);
			if (!optional)
				setNewBackground(box, GUI.WARNING_BACKGROUND_COLOR);

			for (final String token : tokens) {
				if (token.equals(box.getText())) {
					//System.out.println("\t\tSelecting Checkbox: " + token);
					doSet = true;
					setKey = key;
					setValue = value;
					setMap = hm;
					box.setSelected(true);
					oneBoxSelected = true;
				}
			}
		}
		
		if (oneBoxSelected &&  !optional) { //one box has been ticked in a non-optional set -> set to normal
			for (final JCheckBox box : boxes) {
				setNewBackground(box, Color.white);
			}
		}
	}


	@Override
	public synchronized void itemStateChanged(ItemEvent arg0) {
		final String selectedOption = ((JCheckBox) arg0.getItemSelectable()).getText();

		if (arg0.getStateChange() == ItemEvent.SELECTED) {
			concat += selectedOption + " ";
		} else {
			concat = concat.replaceAll(new String(" " + selectedOption + " "), " ");
			concat = concat.replaceFirst(new String("^" + selectedOption + " "), ""); //at the beginning of the concat
		}

		if (!optional) {
			if (concat.equals("")) {
				for (final JCheckBox box : boxes) {
					setNewBackground(box, GUI.WARNING_BACKGROUND_COLOR);
				}
			} else {
				for (final JCheckBox box : boxes) {
					setNewBackground(box, Color.white);
				}
			}
		}

		//super.handleItemEvent(arg0, selectedOption);
		super.handleItemEvent(arg0, selectedOption, doSet, setMap, setKey, setValue);
		doSet = false;
		setKey = null;
		setMap = null;
		setValue = null;
	}
}
