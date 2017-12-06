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

import gui_elements.GUIBuildException;

import java.io.IOException;
import java.util.LinkedList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

//for the XML printing
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.*;
import javax.xml.transform.dom.*;
import java.io.StringWriter;

public class XMLHelper {

	/**
	 * Find the parameter values for which the element should be cloned
	 * e.g. for section tags.
	 *
	 * @param element
	 * @param tagName
	 * @param attributeName
	 * @param parentArgument
	 * @return the
	 */
	public static String[] findApplicableValues(Element element,
			String tagName, String attributeName, String parentArgument) {

		final String[] argumentsFromElement = getAttributesFromChildren(
				element, tagName, attributeName);

		/*
		 * If no such arguments are specified for the current element but there
		 * are for the parent element, take that of the parent.
		 */
		if (argumentsFromElement.length == 0 && parentArgument != null)
			return new String[] { parentArgument };
		else if ((argumentsFromElement.length > 0 && parentArgument != null)) {
			if (getCertainChildElement(element, tagName, attributeName,
					parentArgument) != null)
				return new String[] { parentArgument };
			else
				/*
				 * The empty intersection is an empty array so that his case can
				 * be handled differently from the null case below.
				 */
				return new String[0];
		}

		/*
		 * If no arguments are specified for the parent use all element
		 * arguments.
		 */
		else if (argumentsFromElement.length > 0 && parentArgument == null)
			return argumentsFromElement;
		else if (argumentsFromElement.length == 0 && parentArgument == null)
			return null;
		return null;
	}

	/**
	 * Get all children of an element that have a certain element tag and
	 * extract the values for a certain attribute.
	 *
	 * @param parent
	 *            the parent element
	 * @param tag
	 *            the element tag
	 * @param attribute
	 *            the attribute name
	 *
	 * @return a String array with the attribute values
	 */
	public static String[] getAttributesFromChildren(Element parent, String tag, String attribute) {
		final Element[] elements = getElementArray(parent, tag);

		final String[] attributes = new String[elements.length];
		for (int i = 0; i < elements.length; i++) {
			attributes[i] = elements[i].getAttribute(attribute);
		}
		return attributes;
	}

	/**
	 * Returns the first child element that has the specified element tag and
	 * the value for the specified attribute.
	 *
	 * @param parent
	 * @param tag
	 * @param attribute
	 * @param value
	 * @return the first child element that meets the criteria
	 */
	public static Element getCertainChildElement(Element parent, String tag,
			String attribute, String value) {
		final Element[] allElements = getElementArray(parent, tag);
		for (final Element element : allElements) {
			if (element.hasAttribute(attribute)) {
				if (element.getAttribute(attribute).equals(value))
					return element;
			}
		}
		return null;
	}

	/**
	 * For the specified element gets the text that is contained in the first
	 * child that has a certain element name.
	 *
	 *
	 * @param parent
	 *            the parent element
	 * @param tag
	 *            the element tag for the child
	 * @return the text content
	 */
	public static String getChildElementContent(Element parent, String tag) {
		final Element[] elements = getElementArray(parent, tag);
		if (elements.length == 0)
			return null;

		return elements[0].getTextContent();
	}



	/**
	 * @param parent
	 * @param tag
	 * @return all child elements of the parent with the specified tag
	 */
	public static Element[] getElementArray(Element parent, String tag) {
		if (!parent.hasChildNodes())
			return new Element[0];

		final NodeList nodeList = parent.getChildNodes();
		final LinkedList<Element> elementlist = new LinkedList<Element>();
		for (int i = 0; i < nodeList.getLength(); i++) {

			final Node element = nodeList.item(i);
			if (element.getNodeType() == Node.ELEMENT_NODE && ((Element) element).getTagName().equals(tag)) {
				elementlist.add((Element) element);
			}
		}
		final Element[] elements = new Element[elementlist.size()];

		return elementlist.toArray(elements);
	}

	public static Element[] getElementArray(Element parent, String[] tags) {
		if (!parent.hasChildNodes())
			return new Element[0];

		final NodeList nodeList = parent.getChildNodes();
		final LinkedList<Element> elementlist = new LinkedList<Element>();
		for (int ii = 0; ii < nodeList.getLength(); ii++) {

			final Node element = nodeList.item(ii);
			if (element.getNodeType() == Node.ELEMENT_NODE) {
				//Loop through all possible tags
				for (int jj=0; jj<tags.length; jj++) {
					if (((Element) element).getTagName().equals(tags[jj])) {
						elementlist.add((Element) element);
						break;
					}
				}
			}
		}

		final Element[] elements = new Element[elementlist.size()];
		return elementlist.toArray(elements);
	}

	/**
	 * Prints the Document object to the console. as XML.
	 *
	 * @param doc
	 */
	public static String printToScreen(Document doc) {

		final javax.xml.transform.TransformerFactory tfactory = TransformerFactory
				.newInstance();

		javax.xml.transform.Transformer xform;
		try {
			xform = tfactory.newTransformer();

			final javax.xml.transform.Source src = new javax.xml.transform.dom.DOMSource(
					doc);

			final java.io.StringWriter writer = new java.io.StringWriter();
			final javax.xml.transform.Result result = new javax.xml.transform.stream.StreamResult(
					writer);

			xform.transform(src, result);

			return writer.toString();

		} catch (final TransformerConfigurationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (final TransformerException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}



	/**
	 * Reads an XML file and writes it into a Document structure.
	 *
	 * @param filepath
	 * @param xsdpath
	 *            the path of the W3 XML schema that the file is based on.
	 * @throws GUIBuildException
	 */
	public static Document readXML(String filepath, String xsdpath)
			throws GUIBuildException {

		/* read XML file and write it into Document structure */
		final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		dbf.setIgnoringComments(true);
		dbf.setValidating(true);
		dbf.setIgnoringElementContentWhitespace(true);

		dbf.setAttribute(
				"http://java.sun.com/xml/jaxp/properties/schemaLanguage",
				"http://www.w3.org/2001/XMLSchema");
		dbf.setAttribute(
				"http://java.sun.com/xml/jaxp/properties/schemaSource", xsdpath);


		DocumentBuilder db;
		Document doc;

		try {
			db = dbf.newDocumentBuilder();
			doc = db.parse(filepath);
		} catch (final ParserConfigurationException e) {
			throw new GUIBuildException(e.getMessage());
		} catch (final SAXException e) {
			throw new GUIBuildException(e.getMessage());
		} catch (final IOException e) {
			throw new GUIBuildException(e.getMessage());
		}

		doc.normalize();
		doc.normalizeDocument();

		resolveIncludes(doc, xsdpath);
		resolveReferences(doc);
		recursivelyResolveSections(doc.getDocumentElement(), null);
		recursivelyResolveReplaces(doc.getDocumentElement(), null);

		return doc;
	}

	/**
	 * This method recursively traverses the document tree and replaces all
	 * parameters that contain replace elements with the corresponding number of
	 * cloned elements.
	 *
	 * @param current
	 *            the current element
	 * @param parentReplace
	 *
	 */
	public static void recursivelyResolveReplaces(Element current,
			String parentReplace) {

		final String[] replacements = XMLHelper.findApplicableValues(
				current, "replace", "name", parentReplace);

		/* delete replace children */
		final Element[] replaceChildren = XMLHelper.getElementArray(current,
				"replace");
		for (final Element replace : replaceChildren) {
			current.removeChild(replace);
		}

		if (replacements == null) {
			final NodeList childNodes = current.getChildNodes();
			final Node[] nodes = new Node[childNodes.getLength()];
			for (int i = 0; i < childNodes.getLength(); i++) {
				nodes[i] = childNodes.item(i);
			}

			for (final Node node : nodes) {
				if (node.getNodeType() == Node.ELEMENT_NODE) {
					XMLHelper.recursivelyResolveReplaces((Element) node,
							parentReplace);
				}
			}
		} else if (replacements.length == 0) {
			current.getParentNode().removeChild(current);
		} else if (replacements != null) {
			/* clone element */
			for (final String replace : replacements) {
				final Element clone = (Element) current.cloneNode(true);
				clone.setAttribute("key",
						clone.getAttribute("key").replaceAll("%", replace));
				clone.setAttribute("replaced", replace);
				current.getParentNode().appendChild(clone);

				/* Recursion: Do not do recursion on NodeList as it is updated */
				final NodeList childNodes = clone.getChildNodes();
				final Node[] nodes = new Node[childNodes.getLength()];
				for (int i = 0; i < childNodes.getLength(); i++) {
					nodes[i] = childNodes.item(i);
				}

				for (final Node node : nodes) {
					if (node.getNodeType() == Node.ELEMENT_NODE) {
						XMLHelper.recursivelyResolveReplaces((Element) node,
								replace);
					}
				}

			}
			current.getParentNode().removeChild(current);
		}
	}

	/**
	 * This method recursively traverses the document tree and replaces all
	 * parameters that contain section elements with the corresponding number of
	 * cloned elements.
	 *
	 * @param current
	 */
	public static void recursivelyResolveSections(Element current,
			String parentSection) {

		final String[] sections = XMLHelper.findApplicableValues(current,
				"section", "name", parentSection);


		final Element[] sectionChildren = XMLHelper.getElementArray(current,
				"section");
		for (final Element section : sectionChildren) {
			current.removeChild(section);
		}

		if (sections == null) {
			final NodeList childNodes = current.getChildNodes();
			final Node[] nodes = new Node[childNodes.getLength()];
			for (int i = 0; i < childNodes.getLength(); i++) {
				nodes[i] = childNodes.item(i);
			}

			for (final Node node : nodes) {
				if (node.getNodeType() == Node.ELEMENT_NODE) {
					XMLHelper.recursivelyResolveSections((Element) node,
							parentSection);
				}
			}
		} else if (sections != null) {
			/* clone element */
			for (final String section : sections) {
				final Element clone = (Element) current.cloneNode(true);
				clone.setAttribute("section", section);

				current.getParentNode().appendChild(clone);

				/* Recursion: Do not do recursion on NodeList as it is updated */
				final NodeList childNodes = clone.getChildNodes();
				final Node[] nodes = new Node[childNodes.getLength()];
				for (int i = 0; i < childNodes.getLength(); i++) {
					nodes[i] = childNodes.item(i);
				}

				for (final Node node : nodes) {
					if (node.getNodeType() == Node.ELEMENT_NODE) {
						XMLHelper.recursivelyResolveSections((Element) node,
								section);
					}
				}

			}
			current.getParentNode().removeChild(current);
		}
	}

	/**
	 * Parses an XML file as specified in the path of an "include" element and
	 * replaces the include element with the content of the XML file.
	 *
	 * @param doc
	 */
	public static Document resolveIncludes(Document doc, String xsdpath) {
		final NodeList include = doc.getElementsByTagName("include");

		for (int i = 0; i < include.getLength(); i++) {
			final String includePath = ((Element) include.item(i))
					.getAttribute("path");

			try {
				final Document partialDocument = readXML(includePath, xsdpath);
				final NodeList children = partialDocument.getChildNodes()
						.item(0).getChildNodes();
				final Node parent = include.item(i).getParentNode();

				for (int j = 0; j < children.getLength(); j++) {
					final Node imported = doc
							.importNode(children.item(j), true);
					parent.appendChild(imported);
				}

				// parent.replaceChild(fragment, include.item(i));
			} catch (final GUIBuildException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}

		return doc;
	}

	/**
	 * Manipulates the document structure, so that "reference" elements are
	 * replaced by the corresponding "parametergroup".
	 *
	 * @param doc
	 * @return
	 */
	public static Document resolveReferences(Document doc) {
		final NodeList references = doc.getElementsByTagName("reference");
		final NodeList parametergroups = doc.getElementsByTagName("parametergroup");

		while (references.getLength() != 0) {
			final Element referenceElement = (Element) references.item(0);
			final Element parent = (Element) referenceElement.getParentNode();

			parent.removeChild(referenceElement);

			for (int i = 0; i < parametergroups.getLength(); i++) {
				final Node parametergroupElement = parametergroups.item(i);

				if (((Element) parametergroupElement).getAttribute("name")
						.equals(referenceElement.getAttribute("name"))) {
					final NodeList children = parametergroupElement
							.getChildNodes();

					for (int j = 0; j < children.getLength(); j++) {
						final Node clone = children.item(j).cloneNode(true);
						parent.appendChild(clone);
					}
				}
			}
		}

		while (parametergroups.getLength() != 0) {
			parametergroups.item(0).getParentNode()
					.removeChild(parametergroups.item(0));
		}

		return doc;
	}


	/**
	 * @brief  The funciton expects a Node element as argument and transform the XML
	 *         (sub)tree starting with the given Node into a string
	 * @return A string representing the (sub)tree starting with node
	 * @author Thomas Egger
	 */
	public static String getString(Node node){
		String str = "";

		try {
			Transformer transformer = TransformerFactory.newInstance().newTransformer();
			final StringWriter buffer = new StringWriter();
			transformer.transform(new DOMSource(node), new StreamResult(buffer));
			str = buffer.toString();
		} catch (Exception ex) {
			throw new RuntimeException("Error converting to String", ex);
		}

		return str;
	}
}
