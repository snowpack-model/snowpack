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

import main.GUI;

import java.awt.*;
import java.awt.event.*;
import java.util.HashMap;

import javax.swing.JTextField;
import javax.swing.JOptionPane;

import org.w3c.dom.Element;

/**
 * A panel that contains a textfield.
 *
 * @author korhammer et egger
 *
 */
public class TextfieldPanel extends ControlledPanel implements KeyListener {

	private static final long serialVersionUID = 1L;
	private String defaultText;

	protected JTextField textfield;

	/**
	 * The constructor for a TextfieldPanel.
	 *
	 * @param element
	 */
	public TextfieldPanel(Element element, ControlledPanel parent) {
		super(element, parent);
		textfield = new JTextField(10);
		textfield.setVisible(true);
		this.add(textfield, "width 20:100:400, cell 1 0, growx, ay top, gaptop 3");
		textfield.addKeyListener(this);

		defaultText = "";
		if (element.hasAttribute("default")) {
			defaultText = element.getAttribute("default");
			if(optional) {
				textfield.setFont(textfield.getFont().deriveFont(Font.ITALIC));
				textfield.setForeground(GUI.DFLT_VAL_COLOR);
			}
		}
		textfield.setHorizontalAlignment(JTextField.CENTER);
		textfield.setText(defaultText);

		if (!optional)
			textfield.setBackground(GUI.WARNING_BACKGROUND_COLOR);

		this.checkBgColor();
		textfield.invalidate();
	}



	@Override
	public void close() {
		// TODO Auto-generated method stub

	}

	@Override
	public void set(HashMap hm, String key, String value) {
		if (hm==null) {
			System.out.println("\thashmap null in TextfieldPanel. Set canceled for " + key + " value " + value);
			return;
		}
		final String val = (String)hm.get(hashKey);
		if (val != null) textfield.setText(val);
		this.checkBgColor();
	}

	@Override
	public String getValue() {
		final String thetext = textfield.getText();

		if (!isOptional() && thetext.equals("")) {
			JOptionPane.showMessageDialog(null, "No value was entered for "
									+ getKey() + " (section '" + getSection().toUpperCase() + "')"
									+ ".\nA value is required. Your .ini-file is probably incorrect.",
									"Problem when building .ini file",
									JOptionPane.WARNING_MESSAGE);
		}

		if (thetext.equals(""))
			return null;

		return thetext;
	}

	/** Handle the key typed event from the text field. */
	public void keyTyped(KeyEvent e) {}

	/** Handle the key-pressed event from the text field. */
	public void keyPressed(KeyEvent e) {}

	/** Handle the key-released event from the text field. */
	public void keyReleased(KeyEvent e) {
		if(optional) {
			textfield.setFont(textfield.getFont().deriveFont(Font.PLAIN));
			textfield.setForeground(GUI.VAL_COLOR);
		}
		checkBgColor();
	}

	public void checkBgColor() {
		final String text = textfield.getText();

		if (!optional) {
			if (text.length() > 0) {
				textfield.setBackground(GUI.NORMAL_BACKGROUND_COLOR);
			} else {
				textfield.setBackground(GUI.WARNING_BACKGROUND_COLOR);
			}
			textfield.invalidate();
		}
	}
}
