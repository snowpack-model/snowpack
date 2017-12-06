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

import java.awt.event.*;
import java.util.HashMap;

import java.awt.Font;
import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.SwingConstants;

import main.GUI;

import org.w3c.dom.Element;

/**
 * A panel that maintains a number in an editable textfield,
 * optionally with a maximum and minimum value.
 * @author korhammer
 *
 */
/**
 * @author korhammer
 *
 */
public class IntegerPanel extends TextfieldPanel implements ActionListener, FocusListener {
	/**
	 *
	 */
	private static final long serialVersionUID = -8959578535487185306L;

	private int defaultvalue;
	private int maximumvalue;
	private int minimumvalue;

	private boolean printPlus;

	private final JButton down;
	private final JButton up;

	/**
	 * Constructor for a NumberPanel.
	 *
	 * @param element
	 * @throws GUIBuildException
	 */
	public IntegerPanel(Element element, ControlledPanel parent) throws GUIBuildException {
		super(element, parent);

		setDefaultvalue(0);
		minimumvalue = Integer.MIN_VALUE;
		maximumvalue = Integer.MAX_VALUE;
		printPlus = false;

		if (element.getAttribute("type").equals("integer+")) {
			setPrintplus(true);
		} else {
			setPrintplus(false);
		}

		if (element.hasAttribute("default")) {
			setDefaultvalue(element.getAttribute("default"));
			if(optional) {
				textfield.setFont(textfield.getFont().deriveFont(Font.ITALIC));
				textfield.setForeground(GUI.DFLT_VAL_COLOR);
			}
		}
		if (element.hasAttribute("maximum")) {
			setMaximumvalue(element.getAttribute("maximum"));
		}
		if (element.hasAttribute("minimum")) {
			setMinimumvalue(element.getAttribute("minimum"));
		}

		textfield.setColumns(4);
		textfield.invalidate();
		textfield.setHorizontalAlignment(SwingConstants.CENTER);

		up = new JButton("+");
		up.addActionListener(this);
		down = new JButton("-");
		down.addActionListener(this);
		this.add(up, "split 2, growx, ay top");
		this.add(down, "growx, ay top");
		textfield.addActionListener(this);
		textfield.addFocusListener(this);

	}

	/**
	 * Parses a String to an integer. Different from the standard methods as it
	 * handles the + sign e.g. in time zones.
	 *
	 * @param text
	 * @return the parsed number
	 * @throws NumberFormatException
	 */
	private static int parse(String text) throws NumberFormatException {
		final char first = text.charAt(0);

		if (first == '+')
			return Integer.parseInt(text.substring(1));
		if (first == '-')
			return Integer.parseInt(text);
		return Integer.parseInt(text);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see
	 * java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	@Override
	public void actionPerformed(ActionEvent e) {
		if(optional) {
			textfield.setFont(textfield.getFont().deriveFont(Font.PLAIN));
			textfield.setForeground(GUI.VAL_COLOR);
		}
		if (e.getSource() == up) {
			manipulateNumber(1);
		} else if (e.getSource() == down) {
			manipulateNumber(-1);
		} else if (e.getSource() == textfield) {
			parseTextfield();
		}
	}

	@Override
	public void set(HashMap hm, String key, String value) {
		if (hm==null) {
			System.out.println("\thashmap null in IntegerPanel. Set canceled for " + key + " value " + value);
			return;
		}
		final String val = (String)hm.get(hashKey);
		if (val != null) {
			try {
				Integer.parseInt(val);
				textfield.setText(val);
			} catch (Exception e) {} //Do nothing if conversion fails
		}
	}

	@Override
	public void focusGained(FocusEvent arg0) {

	}

	@Override
	public void focusLost(FocusEvent arg0) {
		if(optional) {
			textfield.setFont(textfield.getFont().deriveFont(Font.PLAIN));
			textfield.setForeground(GUI.VAL_COLOR);
		}
		parseTextfield();
	}

	/**
	 * Manipulates the number in the text field by parsing it, adding the
	 * increment and writing it back into the textfield.
	 *
	 * @param increment
	 */
	public void manipulateNumber(int increment) {
		int number;

		try {
			number = parseTextfield();
		} catch (final NumberFormatException nfe) {
			return;
		}

		if (number + increment <= maximumvalue
				&& number + increment >= minimumvalue) {
			number += increment;
		} else if (number + increment < minimumvalue) {
			number = minimumvalue;
		} else if (number + increment > maximumvalue) {
			number = maximumvalue;
		}

		printToTextfield(number);
	}

	/**
	 * Calls the parse method on the textfield.
	 *
	 * @return
	 * @throws NumberFormatException
	 */
	public int parseTextfield() throws NumberFormatException {
		try {
			return parse(textfield.getText());
		} catch (final NumberFormatException nfe) {
			JOptionPane.showMessageDialog(null, "The textfield " + getKey()
					+ " does not contain a number.");
			throw nfe;
		}
	}

	/**
	 * Prints a number to the text field optionally with a plus sign.
	 *
	 * @param number
	 */
	public void printToTextfield(int number) {
		if (this.printPlus && number > 0) {
			textfield.setText("+" + number);
		} else {
			textfield.setText("" + number);
		}
		checkBgColor();
	}

	public void setDefaultvalue(int defaultvalue) {
		this.defaultvalue = defaultvalue;
		this.printToTextfield(defaultvalue);
	}

	public void setDefaultvalue(String defaultvalue) {
		this.defaultvalue = parse(defaultvalue);
		this.printToTextfield(this.defaultvalue);
	}

	public void setMaximumvalue(int maximumvalue) {
		this.maximumvalue = maximumvalue;
	}

	public void setMaximumvalue(String maximumvalue) {
		this.maximumvalue = parse(maximumvalue);
	}

	public void setMinimumvalue(int minimumvalue) {
		this.minimumvalue = minimumvalue;
	}

	public void setMinimumvalue(String minimumvalue) {
		this.minimumvalue = parse(minimumvalue);
	}

	public void setPrintplus(boolean printPlus) {
		this.printPlus = printPlus;
		printToTextfield(parse(textfield.getText()));
	}
}
