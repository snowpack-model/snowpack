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
import java.awt.event.ItemListener;
import java.util.LinkedList;
import java.util.HashMap;
import javax.swing.event.EventListenerList;

import main.GUIBuilder;
import main.XMLHelper;

import org.w3c.dom.Element;
import org.w3c.dom.Node;

/**
 * A ControlledPanel that maintains a group of choices and invokes panel
 * creation according to the current choices.
 *
 * @author korhammer et egger
 *
 */
public abstract class OptionPanel extends ControlledPanel implements ItemListener {

	private static final long serialVersionUID = 88L;

	protected LinkedList<String> defaultTrues;
	protected EventListenerList listenerList;

	protected String[] values;
	protected String[] labels;
	protected String[] helptexts;

	/**
	 * The constructor for an OptionPanel.
	 *
	 * @param element
	 *            : The xml element which the OptionPanel is created for.
	 * @throws GUIBuildException
	 */
	public OptionPanel(Element element, ControlledPanel parent) throws GUIBuildException {
		this(element, parent, true);
	}

	public OptionPanel(Element element, ControlledPanel parent, boolean doLayout) throws GUIBuildException {
		super(element, parent, doLayout);

		String elementdefault = "";

		if (element.hasAttribute("default"))
			elementdefault = element.getAttribute("default");

		listenerList = new EventListenerList();

		final Element[] options = XMLHelper.getElementArray(element, "option");
		defaultTrues = new LinkedList<String>();

		values = new String[options.length];
		labels = new String[options.length];
		helptexts = new String[options.length];

		for (int i=0;i<options.length;i++) {
			values[i] = options[i].getAttribute("value");
			labels[i] = options[i].getAttribute("label");

			if (labels[i] == "")
				labels[i] = values[i];

			if (elementdefault.equals(values[i]))
				defaultTrues.add(values[i]);

			if (options[i].hasAttribute("default") && options[i].getAttribute("default").equals("true"))
				defaultTrues.add(values[i]);

			helptexts[i] = XMLHelper.getChildElementContent(options[i],"help");
		}
	}

	@Override
	public void set(HashMap hm, String key, String value) {
		if (hm==null) {
			System.out.println("\thashmap null in OptionPanel. Set canceled for " + key + " value " + value);
			return;
		}
		//System.out.println("Trying to change OptionPanel");
		//ControlledPanel child = add("TA2");
		//child.set(element);
	}


	/**
	 * Invokes the appropriate build or destruct methods when changes are made
	 * to choices.
	 *
	 * @param arg0
	 * @param selectedOption
	 */
	public void handleItemEvent(ItemEvent arg0, String selectedOption) {
		handleItemEvent(arg0, selectedOption, false, null, null, null);
	}

	public void handleItemEvent(ItemEvent arg0, String selectedOption, boolean doSet, HashMap hm, String setKey, String setValue) {
		Element optionElement = XMLHelper.getCertainChildElement(this.element, "option", "label", selectedOption);

		if (optionElement == null)
			optionElement = XMLHelper.getCertainChildElement(this.element, "option", "value", selectedOption);

		if (arg0.getStateChange() == ItemEvent.SELECTED) {
			try {
				GUIBuilder.recursiveBuild(optionElement, this, doSet, hm, setKey, setValue);
			} catch (final GUIBuildException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else {
			try {
				GUIBuilder.recursiveDestruct(optionElement, this);
			} catch (final GUIBuildException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}

}
