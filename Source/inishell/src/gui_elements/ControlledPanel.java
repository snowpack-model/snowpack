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

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;

import javax.swing.JEditorPane;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.*;
import java.util.HashMap;
import javax.swing.border.*;

import main.GUI;
import main.GUIBuilder;
import main.XMLHelper;
import net.miginfocom.swing.MigLayout;

import org.w3c.dom.Element;
import org.w3c.dom.Node;

/**
 * Abstract class for JPanels that are maintained with their key and designated
 * section and can return their current value.
 */
public abstract class ControlledPanel extends JPanel {

	private static final long serialVersionUID = 1L;

	private int dependencies;

	protected JLabel jlabel;

	protected String section;
	protected String key;
	protected String label;
	protected String hashKey;

	protected Element element;

	protected boolean optional;
	protected String replace;

	protected int hierarchy;
	protected String helptext;
	protected JEditorPane help;

	/**
	 * Constructor for a controlled panel.
	 *
	 * @param element
	 */
	public ControlledPanel(Element element, ControlledPanel parent) {
		this(element, parent, true);
	}

	public ControlledPanel(Element element, ControlledPanel parent, boolean doLayout) {
		super();
		initialize(element, parent);

		if (doLayout) {
			final int column1 = 200 - getHierarchy() * 5; //Indentation for the label, 5px per hierarchy level

			setLayout(new MigLayout("gapx 6, wrap 8, top, ins n 5 n n", "["+ column1 + "!][100!][100!][100!][100!][100!][300:400:400]", ""));
			add(jlabel, "growx, ay top, gaptop 5, gapleft " + getHierarchy() * 5);
			add(help, "cell 6 0, wrap");
		}
	}

	private void initialize(Element element, ControlledPanel parent) {
		dependencies = 0;

		this.element = element;
		this.key = element.getAttribute("key");
		this.label = this.key;

		this.section = element.getAttribute("section");
		this.hashKey = this.section.toUpperCase() + "::" + this.key.toUpperCase();
		this.optional = element.getAttribute("optional").equals("true");

		jlabel = new JLabel(this.label + ":      ");
		jlabel.setForeground(GUI.LABEL_COLOR);

		this.setBackground(Color.white);
		this.setVisible(true);

		/* Add help to the right side of the panel. */
		helptext = XMLHelper.getChildElementContent(element, "help");
		help = ControlledPanel.createHelpPane("<html>"+helptext+"</html>");
		if (helptext != null) help.setText(helptext);

		if (parent != null) {
			setHierarchy(parent.getHierarchy() + 1);
		} else {
			setHierarchy(0);
		}
	}

	/**
	 * Creates a ControlledPanel for a parameter that matches its type. (i.e.
	 * radio buttons for an option parameter, ...
	 *
	 * @param element
	 * @return the Controlled Panel
	 * @throws GUIBuildException
	 */
	public static ControlledPanel createSingleParameterPanel(Element element, ControlledPanel parent) throws GUIBuildException {
		final String elem_type = element.getAttribute("type");

		if (element.hasAttribute("counter")) {
			return new DuplicatorPanel(element, parent);
		} else if (element.getTagName().equals("frame")) {
			return new FramePanel(element, parent);
		} else if (elem_type.equals("decimal")) {
			return new DecimalPanel(element, parent);
		} else if (elem_type.equals("integer") || elem_type.equals("integer+")){
			return new IntegerPanel(element, parent);
		} else if (elem_type.equals("choice")) {
			return new CheckBoxPanel(element, parent);
		} else if (elem_type.equals("combination")) {
			return new CombinedPanel(element, parent);
		} else if (elem_type.equals("selector")) {
			return new SelectorPanel(element, parent);
		} else if (elem_type.equals("alternative")) {
			return new AlternativePanel(element, parent);
		} else if (elem_type.equals("string")) {
			return new TextfieldPanel(element, parent);
		} else if (elem_type.equals("path") || elem_type.equals("file")) {
			return new PathPanel(element, parent);
		}

		throw new GUIBuildException("Cannot create a GUI panel for a parameter of type " + elem_type);
	}

	public static JEditorPane createHelpPane(String helptext){
		JEditorPane help = new JEditorPane();
		help.putClientProperty(JEditorPane.HONOR_DISPLAY_PROPERTIES, Boolean.TRUE);
		help.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 12));
		help.setEditable(false);
		help.setContentType("text/html");
		help.setForeground(Color.GRAY);
		help.addHyperlinkListener(GUIBuilder.gui);

		if (helptext!=null) {
			help.setText("<html>"+helptext+"</html>");
		} else {
			help.setText("<html></html>");
		}

		return help;
	}

	public void setHierarchy(int val) {
		hierarchy = val;
	}

	public int getHierarchy() {
		return hierarchy;
	}

	/**
	 * This method should be called when the panel is removed.
	 */
	public abstract void close();

	/**
	 * @return the GUI component.
	 */
	public Component getComponent() {
		return this;
	}

	/**
	 * @return the number of dependencies on the panel.
	 */
	public int getDependencies() {
		return dependencies;
	}

	/**
	 * @return the key of the corresponding key / value set
	 */
	public String getKey() {
		return key;
	}

	/**
	 * @return the current value of the key / value set null if the key value
	 *         set should not be printed
	 */
	public abstract String getValue();

	/**
	 * This method should be called when another panel depends needs this panel
	 * to exist.
	 */
	public synchronized void hold() {
		this.dependencies++;
	}

	/**
	 * @return true if dependencies on this panel exist
	 * @throws GUIBuildException
	 */
	public synchronized boolean isNeeded() throws GUIBuildException {
		if (dependencies == 0) {
			return false;
		} else if (dependencies > 0) {
			return true;
		} else {
			throw new GUIBuildException("Inplausible number of dependencies on panel");
		}
	}

	/**
	 * @return true if a value for the underlying parameter is optional
	 */
	public boolean isOptional() {
		return optional;
	}

	/**
	 * This method should be called when a panel no longer needs this panel to
	 * exist. The panel will be removed if no dependency exists.
	 */
	public synchronized void release() {
		this.dependencies--;
	}

	public void set(HashMap hm, String key, String value) {
		setKey(key);
		if (element!=null) this.section = element.getAttribute("section");
		this.hashKey = this.section.toUpperCase() + "::" + this.key.toUpperCase();
		
		System.out.println("\tSet not implemented for ControlledPanel #" + serialVersionUID);
	}
	
	
	public String getHashKey() {
		return hashKey;
	}

	/**
	 * This method sets a new key for a panel and changes its label.
	 *
	 * @param key
	 *            the new key
	 */
	public void setKey(String key) {
		this.key = key;
		element.setAttribute("key", key);
		jlabel.setText(this.label + ":      ");
	}

	/**
	 * Getter method for the section.
	 * @return the section
	 */
	public String getSection() {
		return section;
	}
}
