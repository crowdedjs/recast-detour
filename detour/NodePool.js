/*
Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
Recast4J Copyright (c) 2015 Piotr Piastucki piotr@jtilia.org

This software is provided 'as-is', without any express or implied
warranty.  In no event will the authors be held liable for any damages
arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:
1. The origin of this software must not be misrepresented; you must not
 claim that you wrote the original software. If you use this software
 in a product, an acknowledgment in the product documentation would be
 appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
 misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

import Node from "./Node.js"


class NodePool {

	m_map = [];
	m_nodes = [];

	constructor() {

	}

	clear() {
		this.m_nodes = [];
		this.m_map = [];
	}

	findNodes(id) {
		let nodes = this.m_map[id];
		if (nodes == null) {
			nodes = [];
		}
		return nodes;
	}

	findNode(id) {
		let nodes = this.m_map[id];
		if (nodes != null && !nodes.length == 0) {
			return nodes[0];
		}
		return null;
	}

	findNode(id, state) {
		let nodes = this.m_map[id];
		if (nodes != null) {
			for (let node of nodes) {
				if (node.state == state) {
					return node;
				}
			}
		}
		return null;
	}

	getNode(id, state = 0) {
		let nodes = this.m_map[id];
		if (nodes != null) {
			for (let node of nodes) {
				if (node.state == state) {
					return node;
				}
			}
		}
		return this.create(id, state);
	}

create(id, state) {
	let node = new Node(this.m_nodes.length + 1);
	node.id = id;
	node.state = state;
	this.m_nodes.push(node);
	let nodes = this.m_map[id];
	if (nodes == null) {
		nodes = [];
		this.m_map[id] = nodes;
	}
	nodes.push(node);
	return node;
}

getNodeIdx(node) {
	return node != null ? node.index : 0;
}

getNodeAtIdx(idx) {
	return idx != 0 ? this.m_nodes[idx - 1] : null;
}

getNodeCount() {
	return this.m_nodes.length;
}

// getNode(ref) {
// return this.getNode(ref, 0);
// }

	/*
	
	inline let getMaxNodes() const { return m_maxNodes; }
	inline dtNodeIndex getFirst(bucket) const { return m_first[bucket]; }
	inline dtNodeIndex getNext(i) const { return m_next[i]; }
	*/
}

export default NodePool;