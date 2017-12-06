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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.util.HashMap;

import java.awt.Font;
import java.text.DecimalFormat;
import java.text.ParseException;
import java.text.DecimalFormatSymbols;
import javax.swing.JFormattedTextField;
import javax.swing.JOptionPane;

import main.GUI;

import org.w3c.dom.Element;

public class DecimalPanel extends TextfieldPanel implements ActionListener, FocusListener{

	/**
	 *
	 */
	private static final long serialVersionUID = 1555742510677354674L;
	double defaultvalue;
	double maximumvalue;
	double minimumvalue;

	JFormattedTextField textfield;

	public DecimalPanel(Element element, ControlledPanel parent) {
		super(element, parent);
		this.remove(super.textfield);

		DecimalFormat format = new DecimalFormat( "###,##0.######" );
		DecimalFormatSymbols custom = new DecimalFormatSymbols();
		custom.setDecimalSeparator('.');
		format.setDecimalFormatSymbols(custom);
		format.setGroupingUsed(false);

		textfield = new JFormattedTextField(format);
		textfield.setHorizontalAlignment(JFormattedTextField.CENTER);
		textfield.addActionListener(this);
		textfield.addFocusListener(this);
		this.add(textfield,"cell 1 0, growx, ay top, gaptop 3");
		textfield.setVisible(true);
		this.invalidate();
		defaultvalue = 0.0;
		maximumvalue = Double.MAX_VALUE;
		minimumvalue = -Double.MAX_VALUE;


		if(element.hasAttribute("minimum")){
			minimumvalue = Double.parseDouble(element.getAttribute("minimum"));
		}
		if(element.hasAttribute("maximum")){
			maximumvalue = Double.parseDouble(element.getAttribute("maximum"));
		}
		if(element.hasAttribute("default")){
			defaultvalue = Double.parseDouble(element.getAttribute("default"));
			if(optional) {
				textfield.setFont(textfield.getFont().deriveFont(Font.ITALIC));
				textfield.setForeground(GUI.DFLT_VAL_COLOR);
			}
		}
		setToDefaultValue();
	}

	private void setToDefaultValue() {
		textfield.setValue(new Double(defaultvalue));
		this.invalidate();
	}

	private void checkValue(){
		try {
			textfield.commitEdit();
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		final double currentvalue = ((Number)textfield.getValue()).doubleValue();
		
		boolean value_ok = true;
		if(currentvalue > maximumvalue){
			textfield.setValue(new Double(maximumvalue));
			value_ok = false;
		}
		else if(currentvalue<minimumvalue){
			textfield.setValue(new Double(minimumvalue));
			value_ok = false;
		}
		
		if (!value_ok) {
			final String msg = "The value you entered is out of range for key '"+key+"'.\nIt has been reset to the closest acceptable value.";
			final String title = "Out of range";
			final int returnval = JOptionPane.showConfirmDialog(null, msg, title, JOptionPane.DEFAULT_OPTION, JOptionPane.ERROR_MESSAGE);
		}
	}

	public String getValue(){
		return textfield.getValue().toString();
	}

	@Override
	public void set(HashMap hm, String key, String value) {
		if (hm==null) {
			System.out.println("\thashmap null in DecimalPanel. Set canceled for " + key + " value " + value);
			return;
		}
		final String val = (String)hm.get(hashKey);
		if (val != null) {
			try {
				Double.parseDouble(val);
				textfield.setText(val);
				textfield.commitEdit();
			} catch (Exception e) {} //Do nothing if conversion fails
		}
	}

	@Override
	public void focusGained(FocusEvent arg0) {
		checkValue();
	}

	@Override
	public void focusLost(FocusEvent arg0) {
		if(optional) {
			textfield.setFont(textfield.getFont().deriveFont(Font.PLAIN));
			textfield.setForeground(GUI.VAL_COLOR);
		}
		checkValue();
	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		checkValue();
	}
}
