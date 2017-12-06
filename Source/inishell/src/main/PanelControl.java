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

package main;

import gui_elements.ControlledPanel;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Set;
import java.util.TreeMap;
import java.util.Iterator;

/**
 * Class that maintains all ControlledPanels so that values can be retrieved
 * from panels, ini-files can be written. And a panel is not multiply created.
 *
 * @author korhammer
 *
 */
public class PanelControl {

	LinkedList<String> order;
	HashMap<String, TreeMap<String, ControlledPanel>> panels;
	HashMap<String, String> comments = null;
	HashMap<String, String> added_comments = null;
	HashMap<String, HashMap<String, String>> unusedKeys = null;
	PanelNode rootNode = null;

	/**
	 * Constructor for the panel control.
	 */
	public PanelControl() {
		panels = new HashMap<String, TreeMap<String, ControlledPanel>>();
		order = new LinkedList<String>();

		unusedKeys = new HashMap<String, HashMap<String, String>>();
	}

	public PanelControl(PanelControl p) {
		panels = new HashMap<String, TreeMap<String, ControlledPanel>>(p.panels);
		order = new LinkedList<String>(p.order);
		if (added_comments != null)
			added_comments = new HashMap<String, String>(p.added_comments);

		if (comments != null)
			comments = new HashMap<String, String>(p.comments);
	}

	public void setUnusedKeys(String key, String value) {
		//Parse key
		final int offset = key.indexOf("::");
		final String section = key.substring(0, offset);
		final String sectionkey = key.substring(offset + 2);

		final HashMap<String, String> hm = unusedKeys.get(section);
		if (hm == null) unusedKeys.put(section, new HashMap<String, String>());

		unusedKeys.get(section).put(sectionkey, value);
	}

	public void setComments(HashMap<String, String> inaddedcomments, HashMap<String, String> incomments) {
		this.added_comments = new HashMap<String, String>(inaddedcomments);
		this.comments = new HashMap<String, String>(incomments);
	}

	/**
	 * @param section
	 * @param key
	 * @return true if a panel with such key is maintained in the specified
	 *         section
	 */
	public boolean contains(String section, String key) {
		if (!panels.containsKey(section))
			return false;

		return panels.get(section).containsKey(key);
	}

	/**
	 *
	 * @param section
	 * @param key
	 * @return the panel with the specified section / key pair
	 */
	public ControlledPanel get(String section, String key) {
		return panels.get(section).get(key);
	}

	/**
	 * Prints the ini-file at the specified path.
	 *
	 * @param path
	 */
	public void printToFile(String path) {
		final File file = new File(path);

		try {
			final FileWriter writer = new FileWriter(file);
			writer.write(this.toString());
			writer.close();
		} catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * Adds the panel to the specified section.
	 *
	 * @param section
	 * @param panel
	 */
	public void put(String section, ControlledPanel panel) {
		final String key = panel.getKey();
		if (!panels.containsKey(section)) {
			order.add(section);
			panels.put(section, new TreeMap<String, ControlledPanel>());
		}
		panels.get(section).put(key, panel);
	}

	/**
	 * Removes the panel with the specified section and key.
	 *
	 * @param section
	 * @param key
	 * @return the panel
	 */
	public ControlledPanel remove(String section, String key) {
		return panels.get(section).remove(key);
	}

	/**
	 * Prints the maintained panels to a string.
	 *
	 * @return an ini-file-like structure
	 */
	public void setRootNode(PanelNode inrootNode) {
		rootNode = inrootNode;
	}

	public String extractKey(String instring) {
		final int offset = instring.indexOf("::");
		if (offset == -1) return null;

		return instring.substring(offset+2);
	}

	@Override
	public String toString() {
		String out = "";

		for (final String section : order) {
			/* print section */
			LinkedList<String> keys = new LinkedList<String>();
			final PanelNode sectionNode = rootNode.get(section.toUpperCase());
			if (sectionNode != null) sectionNode.getKeyList(keys);

			//final Set<String> keys = panels.get(section).keySet();
			String keyvaluestring = "";
			for (final String hashKey : keys) {
				final String key = extractKey(hashKey);
				if (key == null) continue;

				final boolean hasKey = panels.get(section).containsKey(key);
				if (hasKey==false) {
					System.out.println("The key '"+key+"' could not be found in section '"+section+"', this is not normal...");
					continue;
				}
				final String value = panels.get(section).get(key).getValue();
				if (value != null) {
					if ((comments != null) && comments.containsKey(hashKey.toUpperCase())) {
						final String precomments =  comments.get(hashKey.toUpperCase());
						if (!precomments.equals(""))
						 	keyvaluestring += precomments;
					}

					keyvaluestring += key + "\t=\t" + value;

					if ((added_comments != null) && added_comments.containsKey(hashKey.toUpperCase()))
						keyvaluestring += "\t" + added_comments.get(hashKey.toUpperCase());

					keyvaluestring += "\n";
				}
			}

			//Add unused keys
			final HashMap<String, String> extra = unusedKeys.get(section.toUpperCase());
			if (extra != null) {
				if (!keyvaluestring.equals("")) keyvaluestring += "\n";

				for (final String extrakey : extra.keySet()) {
					final String hashKey = section.toUpperCase() + "::" + extrakey.toUpperCase();
					if ((comments != null) && comments.containsKey(hashKey)) {
						final String precomments =  comments.get(hashKey);
						if (!precomments.equals(""))
						 	keyvaluestring += precomments;
					}

					keyvaluestring += extrakey + "\t=\t" + extra.get(extrakey);

					if ((added_comments != null) && added_comments.containsKey(hashKey))
					    keyvaluestring += "\t" + added_comments.get(hashKey);

					keyvaluestring += "\n";
				}
			}

			//check if there aren't any comments in that section
			final String hashKey = section.toUpperCase() + "::";
			if ((comments != null) && comments.containsKey(hashKey)) {
				String precomments =  comments.get(hashKey);
				if (!precomments.equals("")) {
					if (precomments.indexOf("\n") == -1) precomments += "\n";
					keyvaluestring = precomments + keyvaluestring + "\n";
				}
			}

			/* add section only if keys / values are printed */
			final String sectionstring = "[" + section.toUpperCase() + "]\n";
			if (!keyvaluestring.equals("")) {
				out += sectionstring + keyvaluestring + "\n";
			}
		}

		return out;
	}
}
