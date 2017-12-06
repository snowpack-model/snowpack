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
import java.util.HashMap;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JPanel;

import org.w3c.dom.Element;

public class PathPanel extends TextfieldPanel implements ActionListener {

	private static final long serialVersionUID = 1486664773406118295L;

	JButton select_button;

	/**
	 * Creates a panel with a text field and a file / path chooser button.
	 *
	 * @param element
	 * @throws GUIBuildException
	 */
	public PathPanel(Element element, ControlledPanel parent) throws GUIBuildException {
		super(element, parent);

		if (element.getAttribute("type").equals("path")) {
			select_button = new JButton("Select Path");
		} else if (element.getAttribute("type").equals("file")) {
			select_button = new JButton("Select File");
		} else {
			throw new GUIBuildException("PathPanel is not applicable for parameter types other than 'file' and 'path'.");
		}

		this.remove(textfield);
		this.add(textfield, "cell 1 0, growx, span2");

		select_button.addActionListener(this);
		select_button.setActionCommand("open chooser");
		select_button.setVisible(true);
		this.add(select_button, "cell 3 0, growx, wrap");

	}

	@Override
	public void actionPerformed(ActionEvent event) {
		if (event.getActionCommand().equals("open chooser") && event.getSource() == this.select_button) {
			try {
				chooseFile();
			} catch (final GUIBuildException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	@Override
	public void set(HashMap hm, String key, String value) {
		if (hm==null) {
			System.out.println("\thashmap null in PathPanel. Set canceled for " + key + " value " + value);
			return;
		}
		final String val = (String)hm.get(hashKey);
		if (val != null) textfield.setText(val);
		this.checkBgColor();
	}

	/**
	 * Opens the file chooser dialog and writes the path / filepath into the
	 * texfield.
	 *
	 * @throws GUIBuildException
	 */
	public void chooseFile() throws GUIBuildException {
		final JFileChooser filechooser = new JFileChooser(System.getProperty("user.dir"));

		if (element.getAttribute("type").equals("path")) {
			filechooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		} else if (element.getAttribute("type").equals("file")) {
			filechooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
		} else {
			throw new GUIBuildException("PathPanel is not applicable for parameter types other than 'file' and 'path'.");
		}

		final int returnVal = filechooser.showOpenDialog(new JPanel());

		if (returnVal == JFileChooser.APPROVE_OPTION) {
			if (element.getAttribute("type").equals("path")) {
				textfield.setText(filechooser.getSelectedFile().getAbsolutePath());
			} else {
				textfield.setText(filechooser.getSelectedFile().getName());
			}

			checkBgColor();
		}
	}

}
