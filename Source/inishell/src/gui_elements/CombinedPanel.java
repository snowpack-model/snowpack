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
import java.text.NumberFormat;
import java.util.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;

import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFormattedTextField;
import javax.swing.JTextField;
import javax.swing.JLabel;
import javax.swing.JComboBox;

import main.XMLHelper;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import net.miginfocom.swing.MigLayout;

/**
 * @brief The CombinedPanel may contain other panels, such as check and drop down boxes and text fields
 *
 * @author Thomas Egger
 */
public class CombinedPanel extends OptionPanel {
	enum Type { CHOICE, INTEGER, DECIMAL, STRING, DROPDOWN }

	/**
	 *
	 */
	private static final long serialVersionUID = -6751282938445852408L;;

	private final LinkedList<JComponent> valueComponents;

	/**
	 * Constructor for a panel that can hold text, numbers and check boxes. It
	 * corresponds to the parameter type "combination" in the configuration file
	 * and its subcomponents correspond to the option values.
	 *
	 * @param element
	 * @throws GUIBuildException
	 */
	public CombinedPanel(Element element, ControlledPanel parent) throws GUIBuildException {
		super(element, parent, false); //doing layout by itself

		final int column1 = 200 - getHierarchy() * 5;
		this.setLayout(new MigLayout("gapx 6, wrap 8, ins n 5 n n", "["+column1+"!][78!][78!][78!][78!][78!][78!][20!][200:400:400]", ""));

		this.add(jlabel, "gapleft " + getHierarchy() * 5);
		this.add(help, "cell 8 0, wrap");

		valueComponents = new LinkedList<JComponent>();

		final Element[] optionElements = XMLHelper.getElementArray(element, "option");

		byte lblcounter = 0;

		DecimalFormat format = (DecimalFormat)DecimalFormat.getInstance();
		DecimalFormatSymbols custom = new DecimalFormatSymbols();
		custom.setDecimalSeparator('.');
		format.setDecimalFormatSymbols(custom);
		format.setGroupingUsed(false);

		DecimalFormat integerFormat = (DecimalFormat)DecimalFormat.getIntegerInstance();
		integerFormat.setGroupingUsed(false);

		int row = 0;
		int column = 1;

		for (int i=0;i<optionElements.length;i++) {
			final Type optionType = Type.valueOf(optionElements[i].getAttribute("type").toUpperCase());
			final String initialValue = optionElements[i].getAttribute("value");
			final String label = optionElements[i].getAttribute("label");

			JComponent component = null;
			JLabel lblcomponent  = null;

			NodeList nodelist = optionElements[i].getChildNodes();

			switch (optionType) {
			case INTEGER:
				component = new JFormattedTextField(integerFormat);
				((JFormattedTextField)component).setColumns(8);
				if ((initialValue != null) && (!initialValue.equals("")))
					((JFormattedTextField)component).setValue(new Integer(initialValue));
				break;
			case DECIMAL:
				component = new JFormattedTextField(format);
				((JFormattedTextField)component).setColumns(8);
				if ((initialValue != null) && (!initialValue.equals("")))
					((JFormattedTextField)component).setValue((new Double(initialValue)));
				break;
			case STRING:
				component = new JTextField(10);
				((JTextField)component).setText(initialValue);
				break;
			case CHOICE:
				component = new JCheckBox(optionElements[i].getAttribute("value"));
				((JCheckBox)component).addItemListener(this);
				break;
			case DROPDOWN:
				component = new JComboBox();

				for (int jj = 0; jj < nodelist.getLength(); jj++) {
					if (nodelist.item(jj).getNodeType() == Node.ELEMENT_NODE) {
						((JComboBox<String>)component).addItem(((Element)nodelist.item(jj)).getAttribute("value"));
					}
				}

				((JComboBox)component).setSelectedItem(optionElements[i].getAttribute("value"));
				((JComboBox)component).addItemListener(this);
				break;
			}

			int width = 1;

			if ((label != null) && (!label.equals(""))){
				if (label.length() > 12) width = 2; //for oversize labels

				if (column+width >= 7) {
					column = 1;
					row++;
				}
				lblcomponent = new JLabel(label);
				valueComponents.add(lblcomponent);

				this.add(lblcomponent, "cell " + column + " " + row + " " + width + " 1");
				column += width;
				width = 1;
				lblcomponent.setVisible(true);
			}

			if (column >= 7) {
				column = 1;
				row++;
			}
			this.add(component, "cell " + column + " " + row + " " + width + " 1");
			column++;

			valueComponents.add(component);
			component.setVisible(true);
		}
	}

	@Override
	public void close() {
		// TODO Auto-generated method stub

	}

	@Override
	public void set(HashMap hm, String inkey, String invalue) {
		if (hm==null) {
			System.out.println("\thashmap null in CombinedPanel. Set canceled for " + inkey + " value " + invalue);
			return;
		}
		if (!hm.containsKey(hashKey)) return;
		final String val = (String)hm.get(hashKey);

		//Now loop through valueComponents and see if the value string can be matched
		if (val == null) return;

		final String[] tokens = val.split("\\s+");
		int counter = 0;

		for (final JComponent comp : valueComponents) {
			boolean digested = false;
			if (tokens.length <= counter) break;

			if (comp.getClass() == JCheckBox.class) {
				final String cbString = ((JCheckBox) comp).getText();
				if (tokens[counter].equals(cbString)) {
					((JCheckBox) comp).setSelected(true);
					digested = true;
				}
			} else if (comp.getClass() == JComboBox.class) {
				JComboBox box = (JComboBox) comp;

				for (int ii=0; ii<box.getItemCount(); ii++) {
					final String itemString = (String) box.getItemAt(ii);
					if (itemString.equals(tokens[counter])) {
						box.setSelectedIndex(ii);
						digested = true;
					}
				}
			} else if (comp.getClass() == JTextField.class) {
				((JTextField) comp).setText(tokens[counter]);
				digested = true;
			} else if (comp.getClass() == JFormattedTextField.class) {
				((JFormattedTextField) comp).setText(tokens[counter]);
				digested = true;
			}

			if ((comp.getClass() != JLabel.class) && (digested))
				counter++;
		}
	}

	@Override
	public String getValue() {
		String value = "";

		for (final JComponent comp : valueComponents) {
			if (comp.getClass() == JTextField.class) {
				value += ((JTextField) comp).getText();
			} else if (comp.getClass() == JFormattedTextField.class) {
				value += ((JFormattedTextField) comp).getText();
			} else if (comp.getClass() == JCheckBox.class) {
				if (((JCheckBox) comp).isSelected())
					value += ((JCheckBox) comp).getText();
			} else if (comp.getClass() == JComboBox.class) {
				value += ((String)((JComboBox)comp).getSelectedItem());
			}

			if (value.length() > 0)
				value += " ";
		}

		if (value.equals("")) return null;
		else	return value;
	}

	@Override
	public void itemStateChanged(ItemEvent arg0) {
		String selectedOption = null;

		if (arg0.getItemSelectable().getClass() == JCheckBox.class) {
			selectedOption = ((JCheckBox) arg0.getItemSelectable()).getText();
			super.handleItemEvent(arg0, selectedOption);
		} else if (arg0.getItemSelectable().getClass() == JComboBox.class) {
			selectedOption = ((String)((JComboBox)arg0.getItemSelectable()).getSelectedItem());
			//HACK: should call super...
		}
	}
}
