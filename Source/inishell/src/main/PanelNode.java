package main;

import java.util.List;
import java.util.LinkedList;

public class PanelNode {
	private String data;
	private PanelNode parent;
	private LinkedList<PanelNode> children = new LinkedList<PanelNode>();

	public PanelNode(String indata) {
		data = indata;
		parent = null;
	}

	public PanelNode(PanelNode p) {
		data = p.data;
		parent = p.parent;
		children = new LinkedList<PanelNode>(p.children);
	}

	public String getData() {
		return data;
	}

	public PanelNode getChild(String data) {
		for (final PanelNode child : children) {
			if (child.data.equals(data))
				return child;
		}

		return null;
	}

	public PanelNode add(PanelNode c) {
		c.parent = this;
		children.add(c);

		return children.getLast();
	}

	public PanelNode get(String searchstring) {
		if (data.equals(searchstring)) return this;

		for (final PanelNode child : children) {
			final PanelNode tmp = child.get(searchstring);
			if (tmp != null) return tmp;
		}
		return null; //Nothing found in subtree
	}

	public boolean remove(PanelNode tmp) {
		return children.remove(tmp);
	}

	public PanelNode getParent() {
		return parent;
	}

	public void getKeyList(LinkedList<String> keyList) {
		if ((data.indexOf('#') == -1) && (data.indexOf('%') == -1))
			keyList.add(data);

		for (PanelNode child : children) {
			child.getKeyList(keyList);
		}
	}
}
