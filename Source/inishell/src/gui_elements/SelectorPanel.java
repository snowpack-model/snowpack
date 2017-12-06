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

import javax.swing.JButton;
import javax.swing.JOptionPane;


import main.GUIBuilder;
import main.XMLHelper;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * This class introduces more flexibility into the inishell interface:
 * The SelectorPanel offers the user to add a specific parameter to
 * the interface.
 *
 * @author Thomas Egger
 */
public class SelectorPanel extends ControlledPanel implements ActionListener {

	private static final long serialVersionUID = 8610877538909126920L;

	String[] parameters;

	JButton addButton;
	JButton removeButton;
	Element template;
	String  templateKey;
	HashMap<String, ControlledPanel> mapOfParameters;
	ControlledPanel cp = null;

	public SelectorPanel(Element element, ControlledPanel parent) {
		super(element, parent);
		templateKey = ""; //will store the pattern of the key

		addButton = new JButton("+");
		this.add(addButton, "cell 1 0, split 2, grow");
		addButton.addActionListener(this);

		removeButton = new JButton("-");
		this.add(removeButton, "cell 1 0, grow, wrap");
		removeButton.addActionListener(this);

		//The following lines try to extract the template key (e.g. %::filters#, %::resample, COPY::%)
		final NodeList nodelist = element.getChildNodes();
		for (int ii = 0; ii < nodelist.getLength(); ii++) {
			if (nodelist.item(ii).getNodeType() == Node.ELEMENT_NODE) {
				final Element currentElement = (Element)nodelist.item(ii);
				if (currentElement.getAttribute("template").equals("true")){
					template = (Element)currentElement.cloneNode(true);
					templateKey = template.getAttribute("key");
					template.setAttribute("template", "false");
				}
			}
		}

		final Element[] options = XMLHelper.getElementArray(element, "option");
		parameters = new String[options.length + 1];
		for (int ii=0; ii < options.length; ii++) {
			parameters[ii] = options[ii].getAttribute("value");
		}
		parameters[options.length] = "Other...";

		key = templateKey; //IMPORTANT: this changes the key only, not the label!
		hashKey = this.getSection().toUpperCase() + "::" + this.getKey().toUpperCase();
		mapOfParameters  = new HashMap<String, ControlledPanel>();
	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		if (arg0.getSource() == addButton) {
			add();
		} else if (arg0.getSource() == removeButton) {
			remove();
		}
	}

	@Override
	public void set(HashMap hm, String key, String value) {
		String myKey = this.section.toUpperCase() + "::" + templateKey.toUpperCase().replace("%","([a-zA-Z0-9_]+)");
		myKey = myKey.replace("#","[1-9]+[0-9]*"); //all integer numbers from 1 to infinity
		final String myKey2 = key.replaceAll(myKey, "$1");
		final String myKey3 = key.replaceAll(myKey, "$0");

		boolean toAdd = true;
		if (mapOfParameters.get(myKey2) != null) toAdd = false; //Parameter already exists

		if (toAdd) {
			ControlledPanel child = add(myKey2, true, hm, key, value);
			cp = child;
			if (child != null) child.set(hm, key, value);
		} else {
			if (cp != null) cp.set(hm, key, value);
		}
	}

	public void add() {
		String param = "";
		String choice = "Other...";

		if (parameters.length > 1)
			choice = (String)JOptionPane.showInputDialog(this, "Add a parameter", "Add", JOptionPane.QUESTION_MESSAGE, null, parameters, null);

		if (choice == null) return;

		if (choice.equals("Other...")) {
			param = JOptionPane.showInputDialog(this, "Enter the name of a meteo parameter\n(e.g. VW_AVG, P2");
		} else {
			param = choice;
		}

		add(param, false, null, null, null);
	}

	public ControlledPanel add(String param, boolean doSet, HashMap hm, String key, String value) {
		if ((param == null) || (param.length() == 0)) return null;
		param = param.toUpperCase();

		if (mapOfParameters.get(param) != null) return null; //Parameter already exists

		Element child = (Element)template.cloneNode(true);
		recursiveKeyReplace(child, param);
		element.appendChild(child);

		ControlledPanel childPanel = null;
		try {
			childPanel = ControlledPanel.createSingleParameterPanel(child, this);
			childPanel.setKey(childPanel.getKey());
			GUIBuilder.gui.addToTab(childPanel, section, this);
			childPanel.hold();
			GUIBuilder.control(childPanel);
			GUIBuilder.recursiveBuild(child, this, doSet, hm, key, value);

			mapOfParameters.put(param, childPanel);
		} catch (final GUIBuildException e1) {
			e1.printStackTrace();
		}

		this.invalidate();

		return childPanel;
	}

	@Override
	public void close() {

	}

	@Override
	public String getValue() {
		return null;
	}

	public void remove() {
		if (mapOfParameters.size() == 0) return;

		final String choice = (String)JOptionPane.showInputDialog(this, "Remove a parameter", "Remove", JOptionPane.QUESTION_MESSAGE, null, mapOfParameters.keySet().toArray(), null);

		if (choice == null) return; //nothing to do

		final String childKey = templateKey.replaceAll("%", choice);
		final Element child = XMLHelper.getCertainChildElement(element, "parameter", "key", childKey);

		if (child == null) return;

		try { //this makes sure that all child nodes are destroyed
			GUIBuilder.recursiveDestruct(child, this);
			element.removeChild(child);
			final ControlledPanel childPanel = GUIBuilder.panelControl.get(section, childKey);

			childPanel.release();
			childPanel.close();

			if (!childPanel.isNeeded())
				GUIBuilder.gui.removeFromTab(childPanel, section, this);

			mapOfParameters.remove(choice);
		} catch (final GUIBuildException e) {
			e.printStackTrace();
		}

		this.invalidate();
	}

	private void recursiveKeyReplace(Element element, String string) {
		if (element.hasAttribute("key"))
			element.setAttribute("key", element.getAttribute("key").replaceAll("%", string));

		final NodeList nodelist = element.getChildNodes();
		for (int i = 0; i < nodelist.getLength(); i++) {
			if (nodelist.item(i).getNodeType() == Node.ELEMENT_NODE) {
				recursiveKeyReplace(((Element) nodelist.item(i)), string);
			}
		}
	}
}
