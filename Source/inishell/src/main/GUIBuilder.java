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
import gui_elements.GUIBuildException;

import java.io.IOException;

import javax.swing.ProgressMonitor;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.filechooser.FileNameExtensionFilter;

import java.util.*;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * The main class to build the .ini file creator. Upon construction it reads an
 * XML file and creates the corresponding GUI. Maintains an ElementControl
 * object.
 *
 * @author korhammer et egger
 *
 */
public class GUIBuilder {
	public static ProgressMonitor monitor;

	public static final String xsd = GUIBuilder.class.getClassLoader().getResource("resources/config_schema.xsd").toString();
	public static String res_filepath = GUIBuilder.class.getClassLoader().getResource("resources/default_config.xml").toString(); //the default config to start with
	public static Document doc;

	public static GUI gui;
	public static PanelControl panelControl;

	public static String application;
	public static String currentConfigFile;

	public static String saveFilePath = System.getProperty("user.dir");


	private static String getFileExtension(final String file_and_path) {
		final String file_sep = System.getProperty("file.separator");

		// get filename without path
		final int fileSepIndex = file_and_path.lastIndexOf(file_sep);
		if (fileSepIndex == -1) {
			return "";
		}
		final String filename = file_and_path.substring(fileSepIndex+1); //remove path

		// get the extension from the filename
		final int extIndex = filename.lastIndexOf(".");
		if (extIndex == -1)
			return "";

		return filename.substring(extIndex+1);
	}

	/**
	 * The main method. Opens the GUI Builder with a menu bar on the side.
	 * If an XML filepath is handed over as an argument,
	 * the XML file is opened.
	 *
	 * @param args
	 * @throws GUIBuildException
	 */
	public static void main(String args[]) throws GUIBuildException {
		new GUIBuilder();

		if ((args.length != 0) && (args[0].length() > 0))
			res_filepath = args[0]; //User may provice xml file as first argument

		buildGUI(res_filepath);
	}

	/**
	 * Builds a GUI from the XML-file specified with validation through
	 * the specified schema.
	 *
	 * @param filepath
	 * @param schemapath
	 * @throws GUIBuildException
	 */
	public GUIBuilder()	throws GUIBuildException {
		if (System.getProperty("mrj.version") != null) {  //detect if running on mac
			System.setProperty("com.apple.mrj.application.apple.menu.about.name", "Inishell");
			System.setProperty("apple.awt.brushMetalLook", "true");
		}
		gui = new GUI();
	}

	public static void setComments(HashMap<String, String> added_comments, HashMap<String, String> comments) {
		panelControl.setComments(added_comments, comments);
	}

	/**
	 * Builds the gui for a filepath that can either be handed over as an argument
	 * on application start or chosen by using the open button.
	 *
	 * @param filepath
	 * @throws GUIBuildException
	 */
	public static void buildGUI(String filepath) throws GUIBuildException{
		currentConfigFile = filepath;

		doc = XMLHelper.readXML(filepath, xsd);
		doc = XMLHelper.resolveReferences(doc);

		application = doc.getDocumentElement().getAttribute("application");
		gui.setApplicationForTitle(application);

		panelControl = new PanelControl();

		final Element root = doc.getDocumentElement();
		recursiveBuild(root, null, false);
	}


	public static void setValues(final HashMap<String, String> hm, final HashMap<String, String> added_comments, final HashMap<String, String> comments) {
		if (hm == null) return;

		monitor = new ProgressMonitor(gui, "Loading INI file", "0 out of " + panelControl.order.size() + " sections loaded", 0, panelControl.order.size());
		monitor.setMillisToPopup(0);

		Thread myrunner = new Thread() {
				public void run() {
					int progress_counter = 0; //Effectively counts the sections already parsed into the GUI
					try {
						final int returnval = GUIBuilder.closeFile();
						if (returnval == JOptionPane.CANCEL_OPTION) return;
						monitor.setProgress(progress_counter);

						GUIBuilder.buildGUI(currentConfigFile);
					} catch (GUIBuildException e) {
						e.printStackTrace();
					}

					gui.setEnabled(false); //Disable user interaction while the ProgressMonitor is on

					HashSet<String> usedKeys = new HashSet<String>(); //save all keys that were used from the HashMap
					HashSet<String> unusedKeys = new HashSet<String>(); //save all keys that were used from the HashMap

					Iterator<Map.Entry<String, TreeMap<String, ControlledPanel>>> sectionIterator = panelControl.panels.entrySet().iterator();
					while (sectionIterator.hasNext()) {
						progress_counter++;
						monitor.setProgress(progress_counter);
						monitor.setNote(progress_counter + " out of " + panelControl.order.size() + " sections loaded");

						if (monitor.isCanceled()) { //user cancelled loading
							monitor.close();
							gui.setEnabled(true);
							break; //HACK, should reset the whole environment before breaking
						}

						final Map.Entry<String, TreeMap<String, ControlledPanel>> entry = sectionIterator.next();
						final TreeMap<String, ControlledPanel> value = entry.getValue();

						final List<String> sectionKeys = new ArrayList<String>(hm.keySet()); //all keys as defined in the existing ini
						final String sectionKey = entry.getKey().toUpperCase() + "::"; //The section string
						filterList(sectionKey, sectionKeys); //filter keys for current section only (optimization)

						LinkedList<ControlledPanel> todos = new LinkedList<ControlledPanel>();
						LinkedList<String> keys = new LinkedList<String>();
						LinkedList<String> values = new LinkedList<String>();

						do {
							todos.clear();
							keys.clear();
							values.clear();

							final Iterator<Map.Entry<String, ControlledPanel>> keyIterator = value.entrySet().iterator();
							while(keyIterator.hasNext()) { //go through all keys of current section
								final ControlledPanel mypanel = keyIterator.next().getValue();

								String pattern = mypanel.getHashKey();
								pattern = pattern.replace("%","[a-zA-Z0-9_]+");
								pattern = pattern.replace("#","[1-9]+[0-9]*"); //all integer numbers from 1 to infinity

								//Loop through all defined keys of this section and compare them with pattern
								final Iterator<String> keyit = sectionKeys.iterator();
								while (keyit.hasNext()) {
									final String currkey = keyit.next();

									if ((!usedKeys.contains(currkey)) && currkey.matches(pattern)) {
										//Add mypanel to list of todos
										todos.add(mypanel);
										keys.add(currkey);
										values.add(hm.get(currkey));

										usedKeys.add(currkey);
									}
								}
							}

							for (int jj = 0; jj < todos.size(); jj++) {
								todos.get(jj).set(hm, keys.get(jj), values.get(jj));
							}
						} while (todos.size() != 0);

						if (progress_counter == panelControl.order.size()) gui.setEnabled(true);
					}

					//Now find out all keys that were not used in the hashmap
					for (final String key : hm.keySet()) {
						if (!usedKeys.contains(key)) {
							//System.out.println("Unknown key = " + key);
							unusedKeys.add(key);
							panelControl.setUnusedKeys(key, hm.get(key));
						}
					}

					setComments(added_comments, comments);
				}
			};//thread
		myrunner.start();
	}

	private static void filterList(String filter, List<String> list) {
		for (Iterator<String> it=list.iterator(); it.hasNext();) {
			final String next = it.next();
			if (!next.startsWith(filter)) {
				it.remove();
			}
		}
	}

	public static int closeFile(){
		final String msg = "You will lose all changes made to the current ini-file. " + "Do you want this?";
		final String title = "Inishell is already open for " + application + "!";
		final int returnval = JOptionPane.showConfirmDialog(null, msg, title, JOptionPane.OK_CANCEL_OPTION, JOptionPane.WARNING_MESSAGE);
		
		if (returnval == JOptionPane.OK_OPTION) {
			gui.closeAllTabs();
			gui.rootNode = new PanelNode("");
			application = null;
			doc = null;
			panelControl = null;
		}

		return returnval;
	}


	/**
	 * Puts the specified panel into the panel control.
	 * @param panel
	 */
	public static void control(ControlledPanel panel) {
		final String section = panel.getSection();
		panelControl.put(section, panel);
	}


	/**
	 * Recursively traverses all children of the specified element
	 * and builds the corresponding panels.
	 *
	 * @param parentElement
	 * @param parentPanel
	 * @throws GUIBuildException
	 */
	public static void recursiveBuild(Element parentElement, ControlledPanel parentPanel, boolean doSet) throws GUIBuildException {
		recursiveBuild(parentElement, parentPanel, doSet, null, null, null);
	}


	public static void recursiveBuild(Element parentElement, ControlledPanel parentPanel, boolean doSet, HashMap hm, String key, String value) throws GUIBuildException {
		/* get all child elements for the root */
		final String[] tags = {"parameter", "frame"};
		final Element[] parameterElements = XMLHelper.getElementArray(parentElement, tags);

		for (final Element element : parameterElements) {
			final String section = element.getAttribute("section");
			ControlledPanel parameterPanel;

			if (!panelControl.contains(section, element.getAttribute("key"))) {
				parameterPanel = ControlledPanel.createSingleParameterPanel(element, parentPanel);
				parameterPanel.setKey(parameterPanel.getKey());

				if (!element.getAttribute("template").equals("true")) {
					if (doSet) parameterPanel.set(hm, parameterPanel.getKey(), value);

					gui.addToTab(parameterPanel, section, parentPanel);
					control(parameterPanel);
				}
			} else {
				parameterPanel = panelControl.get(section, element.getAttribute("key"));
			}

			parameterPanel.hold();
			recursiveBuild(element, parameterPanel, doSet);
		}

		gui.validate();
	}

	/**
	 * Recursively destructs child elements of the specified root.
	 *
	 * @param rootElement
	 * @param rootPanel
	 * @throws GUIBuildException
	 */
	public static void recursiveDestruct(Element rootElement, ControlledPanel rootPanel)
		throws GUIBuildException {
		
		/* get all child elements for the root */
		final Element[] parameterElements = XMLHelper.getElementArray(rootElement, "parameter");

		for (final Element element : parameterElements) {
			final String section = element.getAttribute("section");
			final ControlledPanel parameterPanel = panelControl.get(section, element.getAttribute("key"));
			recursiveDestruct(element, parameterPanel);

			parameterPanel.release();
			parameterPanel.close();

			if (!parameterPanel.isNeeded()) {
				gui.removeFromTab(parameterPanel, section, rootPanel);
			}
		}

		gui.validate();
	}


	/**
	 * Opens a prompt for the path to save the ini file to. Then prints the ini
	 * to the path specified.
	 *
	 * @throws IOException
	 * @throws GUIBuildException
	 */
	public static void printIOFile() throws IOException, GUIBuildException {
		final FileNameExtensionFilter inifilter = new FileNameExtensionFilter(".ini files", "ini");

		final JFileChooser filechooser = new JFileChooser(saveFilePath);
		filechooser.setFileFilter(inifilter);

		final int returnval = filechooser.showSaveDialog(new JPanel());

		if (returnval == JFileChooser.APPROVE_OPTION) {
			String file_and_path = filechooser.getSelectedFile().getAbsolutePath();
			final String ext = getFileExtension(file_and_path).toLowerCase();
			if(!ext.equals("ini")) {
				file_and_path += ".ini";
			}

			if (filechooser.getSelectedFile().isDirectory()) return;
			saveFilePath = filechooser.getSelectedFile().getParent();

			panelControl.printToFile(file_and_path);
			gui.setIniFileName(file_and_path);
			gui.setApplicationForTitle(application);
		}
	}
}
