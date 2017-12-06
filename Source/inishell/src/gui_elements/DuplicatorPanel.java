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
import java.util.*;

import javax.swing.JButton;

import main.GUIBuilder;
import main.XMLHelper;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class DuplicatorPanel extends ControlledPanel implements ActionListener {
	private static final long serialVersionUID = 6801877538909126920L;
	JButton addButton;
	Element clone;

	JButton removeButton;

	public DuplicatorPanel(Element element, ControlledPanel parent) {
		super(element, parent);

		addButton = new JButton("+");
		this.add(addButton, "cell 1 0, split 2, grow");
		addButton.addActionListener(this);

		removeButton = new JButton("-");
		this.add(removeButton, "cell 1 0, grow, wrap");
		removeButton.addActionListener(this);

		clone = (Element) element.cloneNode(true);
		clone.removeAttribute("counter");
		clone.setAttribute("counted", "true");
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
		if (key.indexOf('#') == -1) {
			ControlledPanel child = add(true, hm, key, value);
			child.set(hm, key, value);
		}
	}

	public ControlledPanel add() {
		return add(false, null, null, null);
	}

	public ControlledPanel add(boolean doSet, HashMap hm, String key, String value) {
		final Element child = (Element) this.clone.cloneNode(true);
		final int counter = Integer.parseInt(element.getAttribute("counter"));
		element.setAttribute("counter", (counter + 1) + "");

		recursiveKeyReplace(child, counter + "");
		element.appendChild(child);

		ControlledPanel childPanel = null;
		try {
			childPanel = ControlledPanel.createSingleParameterPanel(child, this);
			childPanel.setKey(childPanel.getKey());
			GUIBuilder.gui.addToTab(childPanel, section, this);
			childPanel.hold();
			GUIBuilder.control(childPanel);
			GUIBuilder.recursiveBuild(child, this, doSet, hm, key, value);

		} catch (final GUIBuildException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		this.invalidate();

		return childPanel;
	}

	@Override
	public void close() {
		boolean doRemove = true;

		do {
			doRemove = remove();
		} while (doRemove);
	}

	@Override
	public String getValue() {
		return null;
	}

	public boolean remove() {
		final int counter = Integer.parseInt(element.getAttribute("counter")) - 1;
		final String childKey = element.getAttribute("key").replace("#", counter + "");
		final Element child = XMLHelper.getCertainChildElement(element, "parameter", "key", childKey);

		if (child == null) return false;

		element.removeChild(child);
		element.setAttribute("counter", (counter) + "");

		final ControlledPanel childPanel = GUIBuilder.panelControl.get(section, childKey);

		if (childPanel != null) {
			childPanel.release();
			childPanel.close();

			try {
				if (!childPanel.isNeeded()) {
					GUIBuilder.gui.removeFromTab(childPanel, section, this);
				}
			} catch (final GUIBuildException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		this.invalidate();

		return true;
	}

	private void recursiveKeyReplace(Element element, String string) {
		if (element.hasAttribute("key")) {
			element.setAttribute("key", element.getAttribute("key").replaceAll("#", string));
		}
		final NodeList nodelist = element.getChildNodes();
		for (int i = 0; i < nodelist.getLength(); i++) {
			if (nodelist.item(i).getNodeType() == Node.ELEMENT_NODE) {
				recursiveKeyReplace(((Element) nodelist.item(i)), string);
			}
		}
	}

}
