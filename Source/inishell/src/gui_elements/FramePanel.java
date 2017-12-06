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

import main.GUIBuilder;
import main.XMLHelper;

import java.util.HashMap;
import javax.swing.JPanel;
import javax.swing.border.*;
import javax.swing.BorderFactory;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import net.miginfocom.swing.MigLayout;

public class FramePanel extends ControlledPanel {

	private static final long serialVersionUID = 1987664773406118295L;
	protected String title = "";

	/**
	 *
	 *
	 * @param element
	 * @throws GUIBuildException
	 */
	public FramePanel(Element element, ControlledPanel parent) throws GUIBuildException {
		super(element, parent, false);
		this.setLayout(new MigLayout("wrap 8, top, ins n 5 n n", "[200!][100!][100!][100!][100!][100!][300:400:400]", ""));

		if (parent == null) hierarchy = 1;

		title = element.getAttribute("label");

		final Border mytitledborder = BorderFactory.createTitledBorder(title);
		//Border loweredetched = BorderFactory.createEtchedBorder(EtchedBorder.LOWERED);
		this.setBorder(mytitledborder);

		this.invalidate();
	}

	@Override
	public String getValue() {
		return null;
	}

	@Override
	public void close() {

	}
}
