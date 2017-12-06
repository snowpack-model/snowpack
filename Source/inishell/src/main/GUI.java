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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Desktop;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Arrays;

import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JEditorPane;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JToolBar;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.io.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

import net.miginfocom.swing.MigLayout;

/**
 * The GUI created by the author.
 *
 * @author korhammer et egger
 *
 */
public class GUI extends JFrame implements WindowListener, ActionListener, HyperlinkListener {

	/**
	 *
	 */
	private static final long serialVersionUID = -2428732881181670225L;
	public static PanelNode rootNode = null;

	private static String title = GUI.class.getPackage().getImplementationTitle() + " " + GUI.class.getPackage().getImplementationVersion();

	public static final Color DFLT_VAL_COLOR = Color.decode("#AFADAD");
	public static final Color VAL_COLOR = Color.decode("#000000");
	public static final Color LABEL_COLOR = Color.decode("#000080");
	public static Color NORMAL_BACKGROUND_COLOR = Color.decode("#EEEEEE");
	public static Color WARNING_BACKGROUND_COLOR = Color.decode("#FF0011");

	public static ImageIcon icon16 = new ImageIcon(GUI.class.getClassLoader().getResource("resources/icons/inishell16.png"));
	public static ImageIcon icon32 = new ImageIcon(GUI.class.getClassLoader().getResource("resources/icons/inishell32.png"));
	public static ImageIcon icon48 = new ImageIcon(GUI.class.getClassLoader().getResource("resources/icons/inishell48.png"));
	public static ImageIcon icon256 = new ImageIcon(GUI.class.getClassLoader().getResource("resources/icons/inishell256.png"));

	public static ImageIcon openicon = new ImageIcon(GUI.class.getClassLoader().getResource("resources/icons/xml_open.png"));
	public static ImageIcon openiniicon = new ImageIcon(GUI.class.getClassLoader().getResource("resources/icons/fileopen.png"));
	public static ImageIcon exporticon = new ImageIcon(GUI.class.getClassLoader().getResource("resources/icons/filesave.png"));
	public static ImageIcon previewicon = new ImageIcon(GUI.class.getClassLoader().getResource("resources/icons/kghostview.png"));
	public static ImageIcon reseticon = new ImageIcon(GUI.class.getClassLoader().getResource("resources/icons/reset.png"));

	private final JButton createbutton;
	private final JButton previewButton;
	private final JButton openbutton;
	private final JButton openinibutton;
	private final JButton resetbutton;

	private final JTabbedPane tabpane;
	private final JToolBar toolbar;
	private HashMap<String,JPanel> tabpanels;

	private String configFilePath, iniFilePath, iniFileName;

	/*
	 * Creates a window with the static components of the GUI, i.e. components
	 * that do not depend on the configuration file.
	 */
	public GUI() throws GUIBuildException {
		super(title);
		//System.setProperty("com.apple.mrj.application.apple.menu.about.name", "AppName");
		//System.setProperty("apple.laf.useScreenMenuBar", "true");
		final List<Image> app_icons = Arrays.asList( icon16.getImage(),
		                                             icon32.getImage(),
		                                             icon48.getImage(),
		                                             icon256.getImage());
		this.setIconImages(app_icons);

		configFilePath = System.getProperty("user.dir");
		iniFilePath = System.getProperty("user.dir");
		iniFileName = "";

		tabpanels = new HashMap<String,JPanel>();

		Dimension screenDim = new Dimension(1280, 800);
		Dimension maxDim = null;

		if (this.getMaximizedBounds() == null) {
			Toolkit tk = Toolkit.getDefaultToolkit();
			maxDim = tk.getScreenSize();
			this.setMaximizedBounds(new Rectangle(0, 0, (int)maxDim.getWidth(), (int)maxDim.getHeight()));
		} else {
			maxDim = this.getMaximizedBounds().getSize();
		}

		if (maxDim.height < screenDim.height) {
			screenDim.height = maxDim.height;
			if (screenDim.height > 100)
				screenDim.height -= 25;
		}

		if (maxDim.width < screenDim.width) {
			screenDim.width = maxDim.width;
			if (screenDim.width > 100)
				screenDim.width -= 25;
		}

		this.setSize(screenDim);
		//this.setExtendedState(MAXIMIZED_BOTH);
		this.setVisible(true);
		this.setLayout(new BorderLayout());

		this.toolbar = new JToolBar();
		toolbar.setLayout(new BoxLayout(toolbar, 1));
		toolbar.setFloatable(false);
		this.setVisible(true);
		this.add(toolbar, BorderLayout.LINE_START);

		this.tabpane = new JTabbedPane();
		this.add(tabpane);
		tabpane.setVisible(true);

		this.openbutton = new JButton(openicon);
		this.openbutton.setActionCommand("open config");
		this.openbutton.setToolTipText("Open XML-configuration");
		openbutton.addActionListener(this);
		toolbar.add(this.openbutton);

		this.openinibutton = new JButton(openiniicon);
		this.openinibutton.setActionCommand("open ini");
		this.openinibutton.setToolTipText("Open existing INI file");
		openinibutton.addActionListener(this);
		toolbar.add(this.openinibutton);

		this.resetbutton = new JButton(reseticon);
		this.resetbutton.setActionCommand("reset");
		this.resetbutton.setToolTipText("Reset Interface");
		resetbutton.addActionListener(this);
		toolbar.add(this.resetbutton);

		this.createbutton = new JButton(exporticon);
		this.createbutton.setActionCommand("write config");
		this.createbutton.setToolTipText("Write to INI-file");
		createbutton.addActionListener(this);
		toolbar.add(this.createbutton);

		this.previewButton = new JButton(previewicon);
		this.previewButton.setActionCommand("preview");
		this.previewButton.setToolTipText("Preview INI-file");
		previewButton.addActionListener(this);
		toolbar.add(previewButton);

		this.addWindowListener(this);

		rootNode = new PanelNode("");
	}

	/**
	 * Writes the application name into the GUI window title bar.
	 *
	 * @param application
	 */
	public void setApplicationForTitle(String application){
		String full_title = title;
		if (!application.equals("")) {
			full_title += " for " + application;
		}
		if (!iniFileName.equals("")) {
			full_title += " - " + iniFileName;
		}
		
		this.setTitle(full_title);
	}
	
	/**
	 * Set the cached ini file name, so we know the last 
	 * loaded or saved ini file
	 */
	public void setIniFileName(String ini_file_name){
		iniFileName = ini_file_name;
	}

	public void hideOpenButton(){
		openbutton.setVisible(false);
	}

	/**
	 * Closes all tabs.
	 */
	public void closeAllTabs(){
		tabpane.removeAll();
	}

	/**
	 * Adds a component to the tab specified and creates the tab if it does not
	 * already exist.
	 *
	 * @param comp
	 * @param tabName
	 * @throws GUIBuildException
	 */
	public void addToTab(ControlledPanel comp, String tabName, ControlledPanel parentPanel) throws GUIBuildException {
		tabName = Character.toUpperCase(tabName.charAt(0)) + tabName.substring(1).toLowerCase();

		final String section = ((ControlledPanel)comp).getSection().toUpperCase();
		final String key = section + "::" + ((ControlledPanel)comp).getKey();//((ControlledPanel)comp).getHashKey();

		PanelNode sectionNode = rootNode.getChild(section);
		if (sectionNode == null) {
			final PanelNode tmp = new PanelNode(section);
			sectionNode = rootNode.add(tmp);
		}

		if (parentPanel != null) {
			sectionNode = sectionNode.get(parentPanel.getSection().toUpperCase() + "::" + parentPanel.getKey());
		}

		if (sectionNode != null) {
			PanelNode keyNode = sectionNode.get(key);
			if (keyNode == null) {
				final PanelNode tmp = new PanelNode(key);
				keyNode = sectionNode.add(tmp);
			}
		}

		final JPanel tab;
		if (tabpane.indexOfTab(tabName) == -1) {
			tab = new JPanel();
			tab.setBackground(Color.white);
			tab.setVisible(true);
			tab.setLayout(new MigLayout("wrap 1"));

			final JScrollPane scrollpane = new JScrollPane(tab);
			tabpane.add(tabName, scrollpane);
			tabpanels.put(tabName, tab);
		} else {
			tab = tabpanels.get(tabName);
		}

		if (parentPanel != null && containsPanel(tab, parentPanel)) {
			parentPanel.add(comp, "span");
		} else {
			tab.add(comp);
		}

		this.validate();
	}

	public void removeFromTab(ControlledPanel parameterPanel, String section, ControlledPanel parentPanel) {
		final String tabTitle = section.substring(0, 1).toUpperCase() + section.substring(1).toLowerCase();
		final int index = tabpane.indexOfTab(tabTitle);

		final Container container = (Container) tabpane.getComponentAt(index);

		final PanelNode tmp = rootNode.get( parameterPanel.getSection().toUpperCase() + "::" + parameterPanel.getKey() );
		if (tmp==null) //we have a bug here, I don't understand what is going on with OUTPUT::AGGREGATE_PRF
			System.out.println("Could not remove widget " + parameterPanel.getSection().toUpperCase() + "::" + parameterPanel.getKey());
		else 
			tmp.getParent().remove(tmp);

		if (parentPanel != null && containsPanel(container, parentPanel)) {
			parentPanel.remove(parameterPanel);
		} else {
			container.remove(parameterPanel);
		}

		GUIBuilder.panelControl.remove(section, parameterPanel.getKey());

		container.repaint();
		this.validate();
	}

	public boolean containsPanel(Container container, JPanel panel) {

		for (final Component component : container.getComponents()) {
			if (component == panel)
				return true;
			if (containsPanel((Container) component, panel))
				return true;
		}
		return false;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("write config") && GUIBuilder.application!=null ) {
			GUIBuilder.panelControl.setRootNode(rootNode);

			try {
				GUIBuilder.printIOFile();
			} catch (final IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			} catch (final GUIBuildException e2) {
				// TODO Auto-generated catch block
				e2.printStackTrace();
			}
		} else if (e.getActionCommand().equals("preview") && GUIBuilder.application!=null) {
			GUIBuilder.panelControl.setRootNode(rootNode);
			previewINIFile();
		} else if (e.getActionCommand().equals("reset") && GUIBuilder.application!=null) {
			resetGUI();
		} else if( (e.getActionCommand().equals("open config"))) {
			openFile();
		} else if( (e.getActionCommand().equals("open ini"))) {
			openINIFile();
		}
	}

	public void resetGUI() {
		try {
			final int returnval = GUIBuilder.closeFile();
			if (returnval == JOptionPane.CANCEL_OPTION) return;
			GUIBuilder.buildGUI(GUIBuilder.currentConfigFile);
		} catch (GUIBuildException e) {
			e.printStackTrace();
		}
	}

	public void previewINIFile(){
		JFrame frame = new JFrame();
		frame.setSize(600,800);
		frame.setVisible(true);

		JEditorPane pane = new JEditorPane();

		pane.setText(GUIBuilder.panelControl.toString());
		pane.setVisible(true);
		pane.setEditable(false);

		JScrollPane scrollPane = new JScrollPane(pane);
		scrollPane.setVisible(true);
		frame.add(scrollPane);

		frame.validate();
		frame.repaint();
	}

	public void openINIFile(){
		final FileNameExtensionFilter inifilter = new FileNameExtensionFilter(".ini files", "ini");
		final JFileChooser filechooser = new JFileChooser(iniFilePath);
		filechooser.setFileFilter(inifilter);
		final int returnval = filechooser.showOpenDialog(new JPanel());

		HashMap<String, String> hm = new HashMap<String, String>();
		HashMap<String, String> comments = new HashMap<String, String>();
		HashMap<String, String> added_comments = new HashMap<String, String>();

		if (returnval == JFileChooser.APPROVE_OPTION) {
			final String path = filechooser.getSelectedFile().toString();
			if (filechooser.getSelectedFile().isDirectory()) return;
			iniFilePath = filechooser.getSelectedFile().getParent();
			System.out.println("Opening INI file: " + path);
			iniFileName = path;

			try {
				BufferedReader br = new BufferedReader(new FileReader(path));
				String line;

				String section = "GENERAL::";
				String clines = "";

				while((line = br.readLine()) != null) {
					int offset = line.indexOf(";");
					final int offset2 = line.indexOf("#");

					if ((offset2 != -1) && (offset != -1)) {
						if (offset2 < offset) offset = offset2;
					} else if ((offset2 != -1) && (offset == -1)) {
						offset = offset2;
					}

					String comment = null;
					if (offset != -1) {
						comment = line.substring(offset);
					}

					if (-1 != offset) line = line.substring(0, offset);

					offset = line.indexOf("#");
					if (-1 != offset) line = line.substring(0, offset);

					line = line.trim(); //take away ws

					if (line.length()>1) {

						if (line.charAt(0) == '[') { //section
							offset = line.indexOf("]");
							if ((offset == -1) || (offset <= 1)) continue;

							if (!clines.equals("")) {
								int lbr = clines.indexOf("\n");
								clines = clines.substring(0, lbr);
								comments.put(section, clines);
								clines = "";
							}

							section = line.substring(1, offset).toUpperCase() + "::";

						} else { //key-value pair
							final String[] tokens = line.split("=");
							if (tokens.length == 2) {
								final String key = tokens[0].trim();
								final String value = tokens[1].trim();
								if (key.length() > 0) {
									hm.put(section + key.toUpperCase(), value);
									if (comment != null) added_comments.put(section + key.toUpperCase(), comment);
									if (!clines.equals("")) {
										comments.put(section + key.toUpperCase(), clines);
										clines = "";
									}
								}
							}
						}
					} else {
						//it was a comment line only, associate it with the next key
						if (comment != null) clines += comment;
						clines += "\n";
					}
				}
				br.close();

				GUIBuilder.setValues(hm, added_comments, comments);
			} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
			}
		}
	}

	public void openFile(){
		/* If a configuration is already open, ask for close. */
		final FileNameExtensionFilter xmlfilter = new FileNameExtensionFilter(".xml files", "xml");
		final JFileChooser filechooser = new JFileChooser(configFilePath);
		filechooser.setFileFilter(xmlfilter);
		final int returnval = filechooser.showOpenDialog(new JPanel());
		if (returnval == JFileChooser.APPROVE_OPTION) {
			if(GUIBuilder.application != null){
				final int returnval2 = GUIBuilder.closeFile();
				if (returnval2 == JOptionPane.CANCEL_OPTION) return;
			}

			final String path = filechooser.getSelectedFile().toString();
			if (filechooser.getSelectedFile().isDirectory()) return;

			configFilePath = filechooser.getSelectedFile().getParent();

			try {
				GUIBuilder.buildGUI(path);
			} catch (GUIBuildException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	@Override
	public void windowActivated(WindowEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void windowClosed(WindowEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void windowClosing(WindowEvent e) {
		System.exit(NORMAL);

	}

	@Override
	public void windowDeactivated(WindowEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void windowDeiconified(WindowEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void windowIconified(WindowEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void windowOpened(WindowEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void hyperlinkUpdate(HyperlinkEvent e) {

		if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
		      try {
				Desktop.getDesktop().browse(e.getURL().toURI());
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			} catch (URISyntaxException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}
	}
}
