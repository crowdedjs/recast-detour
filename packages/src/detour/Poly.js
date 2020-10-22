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

import NavMesh from "./NavMesh.js"

/** Defines a polyogn within a dtPoly object. */
class Poly {

	index;
	/** The polygon is a standard convex polygon that is part of the surface of the mesh. */
	static DT_POLYTYPE_GROUND = 0;
	/** The polygon is an off-mesh connection consisting of two vertices. */
	static DT_POLYTYPE_OFFMESH_CONNECTION = 1;
	/** Index to first link in linked list. (Or #DT_NULL_LINK if there is no link.) */
	firstLink = 0;
	/** The indices of the polygon's vertices. The actual vertices are located in MeshTile::verts. */
	verts = [];
	/** Packed data representing neighbor polygons references and flags for each edge. */
	neis = [];
	/** The user defined polygon flags. */
	flags = 0;
	/** The number of vertices in the polygon. */
	vertCount = 0;
	/**
	 * The bit packed area id and polygon type.
	 * 
	 * @note Use the structure's set and get methods to access this value.
	 */
	areaAndtype;

	constructor(index, maxVertsPerPoly) {
		this.index = index;
		this.firstLink = NavMesh.DT_NULL_LINK;
		this.verts = new Array(maxVertsPerPoly);
		this.neis = new Array(maxVertsPerPoly);
		for(let i = 0; i < this.verts.length; i++){
			this.verts[i] = 0;
		}
		for(let i = 0; i < this.neis.length; i++){
			this.neis[i] = 0;
		}
	}
	and(v1, v2) {
		var hi = 0x80000000;
		var low = 0x7fffffff;
		var hi1 = ~~(v1 / hi);
		var hi2 = ~~(v2 / hi);
		var low1 = v1 & low;
		var low2 = v2 & low;
		var h = hi1 & hi2;
		var l = low1 & low2;
		return h * hi + l;
	}
	or(v1, v2) {
		var hi = 0x80000000;
		var low = 0x7fffffff;
		var hi1 = ~~(v1 / hi);
		var hi2 = ~~(v2 / hi);
		var low1 = v1 & low;
		var low2 = v2 & low;
		var h = hi1 | hi2;
		var l = low1 | low2;
		return h * hi + l;
	}


	/** Sets the user defined area id. [Limit: < #DT_MAX_AREAS] */
	setArea(a) {
		this.areaAndtype = this.or(this.and(this.areaAndtype , 0xc0) , this.and(a , 0x3f));
	}

	/** Sets the polygon type. (See: #dtPolyTypes.) */
	setType(t) {
		this.areaAndtype = this.or(this.and(this.areaAndtype , 0x3f) | (t << 6));
	}

	/** Gets the user defined area id. */
	getArea() {
		return this.areaAndtype & 0x3;
	}

	/** Gets the polygon type. (See: #dtPolyTypes) */
	getType() {
		return this.areaAndtype >> 6;
	}

}

export default Poly;
