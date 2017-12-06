/***********************************************************************************/
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

import java.awt.event.ItemEvent;

import javax.swing.JComboBox;
import javax.swing.JOptionPane;
import java.util.HashMap;
import javax.swing.*;
import java.awt.*;

import main.GUI;

import org.w3c.dom.Element;

/**
 * The AlternativePanel contains a dropdown list with alternatives to choose from.
 *
 * @author Thomas Egger
 *
 */
public class AlternativePanel extends OptionPanel {

	private static final long serialVersionUID = 4468480774697888798L;
	private JComboBox<String> box;

	private String setKey = null, setValue = null;
	private boolean doSet = false, bool_choice = false;
	private HashMap setMap = null;

	/**
	 * Constructor for the AlternativePanel
	 *
	 * @param element the <parameter ...> element for which the panel is created.
	 * @throws GUIBuildException
	 */
	public AlternativePanel(Element element, ControlledPanel parent) throws GUIBuildException {
		super(element, parent);

		box = new JComboBox<String>(labels);
		box.insertItemAt("", 0);     //Add an empty label
		box.setVisible(true);
		box.setSelectedIndex(0);
		box.addItemListener(this);

		//Issue 263: On OSX the box.setBackground method does not function without this tweak:
		if (System.getProperty("os.name").equals("Mac OS X")) {
			box.setRenderer(new ColorCellRenderer(box.getPreferredSize().width));
			//((JLabel)box.getRenderer()).setHorizontalAlignment(SwingConstants.CENTER);
		} else
			((JLabel)box.getRenderer()).setHorizontalAlignment(SwingConstants.CENTER);

		if (!optional) box.setBackground(GUI.WARNING_BACKGROUND_COLOR);

		if (defaultTrues.size() > 1) {
			throw new GUIBuildException("Element " + getHashKey() + " has multiple default values");
		} else if (defaultTrues.size() == 1) {
			box.setSelectedItem(defaultTrues.getFirst());
			if(optional) {
				box.setFont(box.getFont().deriveFont(Font.ITALIC));
				box.setForeground(GUI.DFLT_VAL_COLOR);
			}
		}
		
		//is this a TRUE/FALSE choice?
		final String first_val = values[0].toUpperCase();
		if (first_val.equals("TRUE") || first_val.equals("FALSE")) //boolean choice?
			bool_choice = true;
		
		for (int ii=1; ii<values.length; ii++) { //loop through all existant items
			//make sure that either we have only booleans or no booleans
			if (values[ii].toUpperCase().equals("TRUE") || values[ii].toUpperCase().equals("FALSE")) {
				if (!bool_choice) 
					throw new GUIBuildException("Element " + getHashKey() + " has a mix of boolean and non-boolean values");
				continue;
			}
			if (bool_choice) 
					throw new GUIBuildException("Element " + getHashKey() + " has a mix of boolean and non-boolean values");
		}

		this.add(box, "cell 1 0, growx, ay top, wrap");
	}

	@Override
	public void close() {
		this.box.setSelectedIndex(0);
	}

	@Override
	public String getValue() {

		if (box.getSelectedIndex() != 0) {
			final int index = box.getSelectedIndex() - 1;
			if (index < 0) return "";

			return values[index];
		}

		if (!isOptional()) {
			JOptionPane.showMessageDialog(null, "No option was selected for "
									+ getKey() + " (section '" + getSection().toUpperCase() + "')"
									+ ".\nA value is required. Your .ini-file is probably incorrect.",
									"Problem when building .ini file",
									JOptionPane.WARNING_MESSAGE);
		}

		return null;
	}

	@Override
	public synchronized void set(HashMap hm, String key, String value) {
		if (hm==null) {
			System.out.println("\thashmap null in AlternativePanel. Set canceled for " + key + " value " + value);
			return;
		}
		doSet = true;
		setKey = key;
		setValue = value;
		setMap = hm;
		if (!hm.containsKey(hashKey)) return;
		String val = ((String)hm.get(hashKey)).toUpperCase();
		if (bool_choice) {
			if (val.equals("0") || val.equals("F")) val="FALSE";
			if (val.equals("1") || val.equals("T")) val="TRUE";
		}

		for (int ii=0; ii<values.length; ii++) { //loop through all existant items
			if (val.equals(values[ii].toUpperCase())) { //compare to val
				this.box.setSelectedIndex(ii+1);
				break;
			}
		}
	}

	@Override
	public synchronized void itemStateChanged(ItemEvent arg0) {
		if(optional) {
			box.setFont(box.getFont().deriveFont(Font.PLAIN));
			box.setForeground(GUI.VAL_COLOR);
		}

		//Deal with help text
		help.setText("");

		if ((box.getSelectedIndex() > 0) && (helptexts[box.getSelectedIndex()-1] != null)) {
			help.setText("<html>" + helptexts[box.getSelectedIndex()-1] + "</html>");
		} else {
			if (helptext != null)
				help.setText("<html>" + helptext + "</html>");
		}

		//handle selected Item
		final String selectedOption = (String) arg0.getItem();
		if (selectedOption == "") {
			if (!optional) {
				box.setBackground(GUI.WARNING_BACKGROUND_COLOR);
			}
		} else {
			if (!optional) {
				box.setBackground(GUI.NORMAL_BACKGROUND_COLOR);
			}
			//System.out.println("\tDOSET of AlternativePanel: " + doSet);
			super.handleItemEvent(arg0, selectedOption, doSet, setMap, setKey, setValue);
			doSet = false;
			setKey = null;
			setMap = null;
			setValue = null;
		}
	}
}

class ColorCellRenderer implements ListCellRenderer {
	protected DefaultListCellRenderer defaultRenderer = new DefaultListCellRenderer();
	private Dimension preferredSize = new Dimension(0, 18);

	public ColorCellRenderer(int width) {
		width -= 50;
		if (width < 100) width = 100; //Default size on OSX
		preferredSize = new Dimension(width, 18);
	}

	public Component getListCellRendererComponent(JList list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
		JLabel renderer = (JLabel) defaultRenderer.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
		renderer.setPreferredSize(preferredSize);
		return renderer;
	}
}
