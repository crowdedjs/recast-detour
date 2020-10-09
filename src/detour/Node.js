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

class Node {

	static DT_NODE_OPEN = 0x01;
	static DT_NODE_CLOSED = 0x02;
	/** parent of the node is not adjacent. Found using raycast. */
	static DT_NODE_PARENT_DETACHED = 0x04;

	index = 0;

	/** Position of the node. */
	pos = new Array(3);
	/** Cost from previous node to current node. */
	cost = 0;
	/** Cost up to the node. */
	total = 0;
	/** Index to parent node. */
	pidx = 0;
	/** extra state information. A polyRef can have multiple nodes with different extra info. see DT_MAX_STATES_PER_NODE */
	state = 0;
	/** Node flags. A combination of dtNodeFlags. */
	flags = 0;
	/** Polygon ref the node corresponds to. */
	id = 0;

	constructor(index) {
		this.index = index;
	}

	toString() {
		return "Node [id=" + id + "]";
	}

}

export default Node;
