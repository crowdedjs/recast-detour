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

import DetourCommon from "./DetourCommon.js"
import VectorPtr from "./VectorPtr.js"
import MeshHeader from "./MeshHeader.js"
import Poly from "./Poly.js"
import PolyDetail from "./PolyDetail.js"
import BVNode from "./BVNode.js"
import MeshData from "./MeshData.js"
import OffMeshConnection from "./OffMeshConnection.js"
import NavMesh from "./NavMesh.js"

function arraycopy(one, oneStart, two, twoStart, len) {
	for (let i = 0; i < len; i++) {
		two[twoStart + i] = one[oneStart + i];
	}
}

class BVItem {
	bmin = new Array(3).fill(0);
	bmax = new Array(3).fill(0);
	i = 0;
};

function compareX(a, b) {
	if (a.bmin[0] < b.bmin[0])
		return -1;
	if (a.bmin[0] > b.bmin[0])
		return 1;
	return 0;
}

function compareY(a, b) {
	if (a.bmin[1] < b.bmin[1])
		return -1;
	if (a.bmin[1] > b.bmin[1])
		return 1;
	return 0;
}

function compareZ(a, b) {
	if (a.bmin[2] < b.bmin[2])
		return -1;
	if (a.bmin[2] > b.bmin[2])
		return 1;
	return 0;
}

class CompareItemX {

	compare(a, b) {
		if (a.bmin[0] < b.bmin[0])
			return -1;
		if (a.bmin[0] > b.bmin[0])
			return 1;
		return 0;
	}

}



class CompareItemZ {

	compare(a, b) {
		if (a.bmin[2] < b.bmin[2])
			return -1;
		if (a.bmin[2] > b.bmin[2])
			return 1;
		return 0;
	}

}

class CompareItemY {

	compare(a, b) {
		if (a.bmin[1] < b.bmin[1])
			return -1;
		if (a.bmin[1] > b.bmin[1])
			return 1;
		return 0;
	}

}


class NavMeshBuilder {

	static MESH_NULL_IDX = 0xffff;




	static calcExtends(items, nitems, imin, imax) {
		let bmin = new Array(3);
		let bmax = new Array(3);
		bmin[0] = items[imin].bmin[0];
		bmin[1] = items[imin].bmin[1];
		bmin[2] = items[imin].bmin[2];

		bmax[0] = items[imin].bmax[0];
		bmax[1] = items[imin].bmax[1];
		bmax[2] = items[imin].bmax[2];

		for (let i = imin + 1; i < imax; ++i) {
			let it = items[i];
			if (it.bmin[0] < bmin[0])
				bmin[0] = it.bmin[0];
			if (it.bmin[1] < bmin[1])
				bmin[1] = it.bmin[1];
			if (it.bmin[2] < bmin[2])
				bmin[2] = it.bmin[2];

			if (it.bmax[0] > bmax[0])
				bmax[0] = it.bmax[0];
			if (it.bmax[1] > bmax[1])
				bmax[1] = it.bmax[1];
			if (it.bmax[2] > bmax[2])
				bmax[2] = it.bmax[2];
		}
		return [bmin, bmax];
	}

	static longestAxis(x, y, z) {
		let axis = 0;
		let maxVal = x;
		if (y > maxVal) {
			axis = 1;
			maxVal = y;
		}
		if (z > maxVal) {
			axis = 2;
			maxVal = z;
		}
		return axis;
	}

	static subdivide(items, nitems, imin, imax, curNode, nodes) {
		let inum = imax - imin;
		let icur = curNode;

		let node = new BVNode();
		nodes[curNode++] = node;

		if (inum == 1) {
			// Leaf
			node.bmin[0] = items[imin].bmin[0];
			node.bmin[1] = items[imin].bmin[1];
			node.bmin[2] = items[imin].bmin[2];

			node.bmax[0] = items[imin].bmax[0];
			node.bmax[1] = items[imin].bmax[1];
			node.bmax[2] = items[imin].bmax[2];

			node.i = items[imin].i;
		} else {
			// Split
			let minmax = NavMeshBuilder.calcExtends(items, nitems, imin, imax);
			node.bmin = minmax[0];
			node.bmax = minmax[1];

			let axis = NavMeshBuilder.longestAxis(node.bmax[0] - node.bmin[0], node.bmax[1] - node.bmin[1],
				node.bmax[2] - node.bmin[2]);

			if (axis == 0) {
				// Sort aPoly x-axis
				//Arrays.sort(items, imin, imin + inum, new CompareItemX());
				//TODO: is this right?
				let shallowCopy = items.slice(imin, imin + inum);
				shallowCopy = shallowCopy.sort(compareX);
				for(let i = imin; i < imin+inum; i++){
					items[i] = shallowCopy[i-imin];
					if(!items[i])
						console.log("uh-oh")

				}
			} else if (axis == 1) {
				// Sort aPoly y-axis
				//Arrays.sort(items, imin, imin + inum, new CompareItemY());
				let shallowCopy = items.slice(imin, imin + inum);
				shallowCopy = shallowCopy.sort(compareY);
				for(let i = imin; i < imin+inum; i++){
					items[i] = shallowCopy[i-imin];
					if(!items[i])
						console.log("uh-oh")
				}
			} else {
				// Sort aPoly z-axis
				//Arrays.sort(items, imin, imin + inum, new CompareItemZ());
				let shallowCopy = items.slice(imin, imin + inum);
				shallowCopy = shallowCopy.sort(compareZ);
				for(let i = imin; i < imin+inum; i++){
					items[i] = shallowCopy[i-imin];
					if(!items[i])
						console.log("uh-oh")
				}
			}

			let isplit = Math.floor(imin + inum / 2);

			// Left
			curNode = NavMeshBuilder.subdivide(items, nitems, imin, isplit, curNode, nodes);
			// Right
			curNode = NavMeshBuilder.subdivide(items, nitems, isplit, imax, curNode, nodes);

			let iescape = curNode - icur;
			// Negative index means escape.
			node.i = -iescape;
		}
		return curNode;
	}

	static createBVTree(params, nodes) {
		// Build tree
		let quantFactor = 1 / params.cs;
		let items = new Array(params.polyCount).fill(null);
		for (let i = 0; i < params.polyCount; i++) {
			let it = new BVItem();
			items[i] = it;
			it.i = i;
			// Calc polygon bounds. Use detail meshes if available.
			if (params.detailMeshes != null) {
				let vb = params.detailMeshes[i * 4 + 0];
				let ndv = params.detailMeshes[i * 4 + 1];
				let bmin = new Array(3);
				let bmax = new Array(3);
				let dv = vb * 3;
				DetourCommon.vCopy(bmin, params.detailVerts, dv);
				DetourCommon.vCopy(bmax, params.detailVerts, dv);
				for (let j = 1; j < ndv; j++) {
					DetourCommon.vMin(bmin, params.detailVerts, dv + j * 3);
					DetourCommon.vMax(bmax, params.detailVerts, dv + j * 3);
				}

				// BV-tree uses cs for all dimensions
				it.bmin[0] = DetourCommon.clamp(Math.round((bmin[0] - params.bmin[0]) * quantFactor), 0, 0xffff);
				it.bmin[1] = DetourCommon.clamp(Math.round((bmin[1] - params.bmin[1]) * quantFactor), 0, 0xffff);
				it.bmin[2] = DetourCommon.clamp(Math.round((bmin[2] - params.bmin[2]) * quantFactor), 0, 0xffff);

				it.bmax[0] = DetourCommon.clamp(Math.round((bmax[0] - params.bmin[0]) * quantFactor), 0, 0xffff);
				it.bmax[1] = DetourCommon.clamp(Math.round((bmax[1] - params.bmin[1]) * quantFactor), 0, 0xffff);
				it.bmax[2] = DetourCommon.clamp(Math.round((bmax[2] - params.bmin[2]) * quantFactor), 0, 0xffff);
			} else {
				let p = i * params.nvp * 2;
				it.bmin[0] = it.bmax[0] = params.verts[params.polys[p] * 3 + 0];
				it.bmin[1] = it.bmax[1] = params.verts[params.polys[p] * 3 + 1];
				it.bmin[2] = it.bmax[2] = params.verts[params.polys[p] * 3 + 2];

				for (let j = 1; j < params.nvp; ++j) {
					if (params.polys[p + j] == NavMeshBuilder.MESH_NULL_IDX)
						break;
					let x = params.verts[params.polys[p + j] * 3 + 0];
					let y = params.verts[params.polys[p + j] * 3 + 1];
					let z = params.verts[params.polys[p + j] * 3 + 2];

					if (x < it.bmin[0])
						it.bmin[0] = x;
					if (y < it.bmin[1])
						it.bmin[1] = y;
					if (z < it.bmin[2])
						it.bmin[2] = z;

					if (x > it.bmax[0])
						it.bmax[0] = x;
					if (y > it.bmax[1])
						it.bmax[1] = y;
					if (z > it.bmax[2])
						it.bmax[2] = z;
				}
				// Remap y
				it.bmin[1] = Math.floor(it.bmin[1] * params.ch / params.cs);
				it.bmax[1] = Math.ceil(it.bmax[1] * params.ch / params.cs);
			}
		}

		return NavMeshBuilder.subdivide(items, params.polyCount, 0, params.polyCount, 0, nodes);
	}

	static XP = 1 << 0;
	static ZP = 1 << 1;
	static XM = 1 << 2;
	static ZM = 1 << 3;

	static classifyOffMeshPoint(pt, bmin, bmax) {

		let outcode = 0;
		outcode |= (pt[0] >= bmax[0]) ? NavMeshBuilder.XP : 0;
		outcode |= (pt[2] >= bmax[2]) ? ZP : 0;
		outcode |= (pt[0] < bmin[0]) ? XM : 0;
		outcode |= (pt[2] < bmin[2]) ? ZM : 0;

		switch (outcode) {
			case NavMeshBuilder.XP:
				return 0;
			case NavMeshBuilder.XP | NavMeshBuilder.ZP:
				return 1;
			case NavMeshBuilder.ZP:
				return 2;
			case NavMeshBuilder.XM | NavMeshBuilder.ZP:
				return 3;
			case NavMeshBuilder.XM:
				return 4;
			case NavMeshBuilder.XM | NavMeshBuilder.ZM:
				return 5;
			case NavMeshBuilder.ZM:
				return 6;
			case NavMeshBuilder.XP | NavMeshBuilder.ZM:
				return 7;
		}

		return 0xff;
	}

	/**
	 * Builds navigation mesh tile data from the provided tile creation data.
	 * 
	 * @param params
	 *            Tile creation data.
	 * 
	 * @return created tile data
	 */
	static createNavMeshData(params) {
		if (params.vertCount >= 0xffff)
			return null;
		if (params.vertCount == 0 || params.verts == null)
			return null;
		if (params.polyCount == 0 || params.polys == null)
			return null;

		let nvp = params.nvp;

		// Classify off-mesh connection points. We store only the connections
		// whose start poPoly is inside the tile.
		let offMeshConClass = null;
		let storedOffMeshConCount = 0;
		let offMeshConLinkCount = 0;

		if (params.offMeshConCount > 0) {
			offMeshConClass = new Array(params.offMeshConCount * 2);

			// Find tight heigh bounds, used for culling out off-mesh start
			// locations.
			let hmin = Number.MAX_VALUE;
			let hmax = -Number.MAX_VALUE;

			if (params.detailVerts != null && params.detailVertsCount != 0) {
				for (let i = 0; i < params.detailVertsCount; ++i) {
					let h = params.detailVerts[i * 3 + 1];
					hmin = Math.min(hmin, h);
					hmax = Math.max(hmax, h);
				}
			} else {
				for (let i = 0; i < params.vertCount; ++i) {
					let iv = i * 3;
					h = params.bmin[1] + params.verts[iv + 1] * params.ch;
					hmin = Math.min(hmin, h);
					hmax = Math.max(hmax, h);
				}
			}
			hmin -= params.walkableClimb;
			hmax += params.walkableClimb;
			let bmin = new Array(3);
			let bmax = new Array(3);
			DetourCommon.vCopy(bmin, params.bmin);
			DetourCommon.vCopy(bmax, params.bmax);
			bmin[1] = hmin;
			bmax[1] = hmax;

			for (let i = 0; i < params.offMeshConCount; ++i) {
				let p0 = new VectorPtr(params.offMeshConVerts, (i * 2 + 0) * 3);
				let p1 = new VectorPtr(params.offMeshConVerts, (i * 2 + 1) * 3);

				offMeshConClass[i * 2 + 0] = NavMeshBuilder.classifyOffMeshPoint(p0, bmin, bmax);
				offMeshConClass[i * 2 + 1] = NavMeshBuilder.classifyOffMeshPoint(p1, bmin, bmax);

				// Zero out off-mesh start positions which are not even
				// potentially touching the mesh.
				if (offMeshConClass[i * 2 + 0] == 0xff) {
					if (p0[1] < bmin[1] || p0[1] > bmax[1])
						offMeshConClass[i * 2 + 0] = 0;
				}

				// Count how many links should be allocated for off-mesh
				// connections.
				if (offMeshConClass[i * 2 + 0] == 0xff)
					offMeshConLinkCount++;
				if (offMeshConClass[i * 2 + 1] == 0xff)
					offMeshConLinkCount++;

				if (offMeshConClass[i * 2 + 0] == 0xff)
					storedOffMeshConCount++;
			}
		}

		// Off-mesh connectionss are stored as polygons, adjust values.
		let totPolyCount = params.polyCount + storedOffMeshConCount;
		let totVertCount = params.vertCount + storedOffMeshConCount * 2;

		// Find portal edges which are at tile borders.
		let edgeCount = 0;
		let portalCount = 0;
		for (let i = 0; i < params.polyCount; ++i) {
			let p = i * 2 * nvp;
			for (let j = 0; j < nvp; ++j) {
				if (params.polys[p + j] == NavMeshBuilder.MESH_NULL_IDX)
					break;
				edgeCount++;

				if ((params.polys[p + nvp + j] & 0x8000) != 0) {
					let dir = params.polys[p + nvp + j] & 0xf;
					if (dir != 0xf)
						portalCount++;
				}
			}
		}

		let maxLinkCount = edgeCount + portalCount * 2 + offMeshConLinkCount * 2;

		// Find unique detail vertices.
		let uniqueDetailVertCount = 0;
		let detailTriCount = 0;
		if (params.detailMeshes != null) {
			// Has detail mesh, count unique detail vertex count and use input
			// detail tri count.
			detailTriCount = params.detailTriCount;
			for (let i = 0; i < params.polyCount; ++i) {
				let p = i * nvp * 2;
				let ndv = params.detailMeshes[i * 4 + 1];
				let nv = 0;
				for (let j = 0; j < nvp; ++j) {
					if (params.polys[p + j] == NavMeshBuilder.MESH_NULL_IDX)
						break;
					nv++;
				}
				ndv -= nv;
				uniqueDetailVertCount += ndv;
			}
		} else {
			// No input detail mesh, build detail mesh from nav polys.
			uniqueDetailVertCount = 0; // No extra detail verts.
			detailTriCount = 0;
			for (let i = 0; i < params.polyCount; ++i) {
				let p = i * nvp * 2;
				let nv = 0;
				for (let j = 0; j < nvp; ++j) {
					if (params.polys[p + j] == NavMeshBuilder.MESH_NULL_IDX)
						break;
					nv++;
				}
				detailTriCount += nv - 2;
			}
		}

		let bvTreeSize = params.buildBvTree ? params.polyCount * 2 : 0;
		let header = new MeshHeader();
		let navVerts = new Array(3 * totVertCount);
		let navPolys = new Array(totPolyCount);
		let navDMeshes = new Array(params.polyCount);
		let navDVerts = new Array(3 * uniqueDetailVertCount);
		let navDTris = new Array(4 * detailTriCount);
		let navBvtree = new Array(bvTreeSize);
		let offMeshCons = new Array(storedOffMeshConCount);

		// Store header
		header.magic = MeshHeader.DT_NAVMESH_MAGIC;
		header.version = MeshHeader.DT_NAVMESH_VERSION;
		header.x = params.tileX;
		header.y = params.tileY;
		header.layer = params.tileLayer;
		header.userId = params.userId;
		header.polyCount = totPolyCount;
		header.vertCount = totVertCount;
		header.maxLinkCount = maxLinkCount;
		DetourCommon.vCopy(header.bmin, params.bmin);
		DetourCommon.vCopy(header.bmax, params.bmax);
		header.detailMeshCount = params.polyCount;
		header.detailVertCount = uniqueDetailVertCount;
		header.detailTriCount = detailTriCount;
		header.bvQuantFactor = 1.0 / params.cs;
		header.offMeshBase = params.polyCount;
		header.walkableHeight = params.walkableHeight;
		header.walkableRadius = params.walkableRadius;
		header.walkableClimb = params.walkableClimb;
		header.offMeshConCount = storedOffMeshConCount;
		header.bvNodeCount = bvTreeSize;

		let offMeshVertsBase = params.vertCount;
		let offMeshPolyBase = params.polyCount;

		// Store vertices
		// Mesh vertices
		for (let i = 0; i < params.vertCount; ++i) {
			let iv = i * 3;
			let v = i * 3;
			navVerts[v] = params.bmin[0] + params.verts[iv] * params.cs;
			navVerts[v + 1] = params.bmin[1] + params.verts[iv + 1] * params.ch;
			navVerts[v + 2] = params.bmin[2] + params.verts[iv + 2] * params.cs;
			// console.log(navVerts)
		}
		// Off-mesh link vertices.
		let n = 0;
		for (let i = 0; i < params.offMeshConCount; ++i) {
			// Only store connections which start from this tile.
			if (offMeshConClass[i * 2 + 0] == 0xff) {
				let linkv = i * 2 * 3;
				let v = (offMeshVertsBase + n * 2) * 3;
				arraycopy(params.offMeshConVerts, linkv, navVerts, v, 6);
				n++;
			}
		}

		// Store polygons
		// Mesh polys
		let src = 0;
		for (let i = 0; i < params.polyCount; ++i) {
			let p = new Poly(i, nvp);
			navPolys[i] = p;
			p.vertCount = 0;
			p.flags = params.polyFlags[i];
			p.setArea(params.polyAreas[i]);
			p.setType(Poly.DT_POLYTYPE_GROUND);
			for (let j = 0; j < nvp; ++j) {
				if (params.polys[src + j] == NavMeshBuilder.MESH_NULL_IDX)
					break;
				p.verts[j] = params.polys[src + j];
				if ((params.polys[src + nvp + j] & 0x8000) != 0) {
					// Border or portal edge.
					let dir = params.polys[src + nvp + j] & 0xf;
					if (dir == 0xf) // Border
						p.neis[j] = 0;
					else if (dir == 0) // Portal x-
						p.neis[j] = NavMesh.DT_EXT_LINK | 4;
					else if (dir == 1) // Portal z+
						p.neis[j] = NavMesh.DT_EXT_LINK | 2;
					else if (dir == 2) // Portal x+
						p.neis[j] = NavMesh.DT_EXT_LINK | 0;
					else if (dir == 3) // Portal z-
						p.neis[j] = NavMesh.DT_EXT_LINK | 6;
				} else {
					// Normal connection
					p.neis[j] = params.polys[src + nvp + j] + 1;
				}

				p.vertCount++;
			}
			src += nvp * 2;
		}
		// Off-mesh connection vertices.
		n = 0;
		for (let i = 0; i < params.offMeshConCount; ++i) {
			// Only store connections which start from this tile.
			if (offMeshConClass[i * 2 + 0] == 0xff) {
				let p = new Poly(offMeshPolyBase + n, nvp);
				navPolys[offMeshPolyBase + n] = p;
				p.vertCount = 2;
				p.verts[0] = offMeshVertsBase + n * 2;
				p.verts[1] = offMeshVertsBase + n * 2 + 1;
				p.flags = params.offMeshConFlags[i];
				p.setArea(params.offMeshConAreas[i]);
				p.setType(Poly.DT_POLYTYPE_OFFMESH_CONNECTION);
				n++;
			}
		}

		// Store detail meshes and vertices.
		// The nav polygon vertices are stored as the first vertices on each
		// mesh.
		// We compress the mesh data by skipping them and using the navmesh
		// coordinates.
		if (params.detailMeshes != null) {
			let vbase = 0;
			for (let i = 0; i < params.polyCount; ++i) {
				let dtl = new PolyDetail();
				navDMeshes[i] = dtl;
				let vb = params.detailMeshes[i * 4 + 0];
				let ndv = params.detailMeshes[i * 4 + 1];
				let nv = navPolys[i].vertCount;
				dtl.vertBase = vbase;
				dtl.vertCount = (ndv - nv);
				dtl.triBase = params.detailMeshes[i * 4 + 2];
				dtl.triCount = params.detailMeshes[i * 4 + 3];
				// Copy vertices except the first 'nv' verts which are equal to
				// nav let verts.
				if (ndv - nv != 0) {
					arraycopy(params.detailVerts, (vb + nv) * 3, navDVerts, vbase * 3, 3 * (ndv - nv));
					vbase += ndv - nv;
				}
			}
			// Store triangles.
			arraycopy(params.detailTris, 0, navDTris, 0, 4 * params.detailTriCount);
		} else {
			// Create dummy detail mesh by triangulating polys.
			let tbase = 0;
			for (let i = 0; i < params.polyCount; ++i) {
				let dtl = new PolyDetail();
				navDMeshes[i] = dtl;
				let nv = navPolys[i].vertCount;
				dtl.vertBase = 0;
				dtl.vertCount = 0;
				dtl.triBase = tbase;
				dtl.triCount = (nv - 2);
				// Triangulate polygon (local indices).
				for (let j = 2; j < nv; ++j) {
					let t = tbase * 4;
					navDTris[t + 0] = 0;
					navDTris[t + 1] = (j - 1);
					navDTris[t + 2] = j;
					// Bit for each edge that belongs to let boundary.
					navDTris[t + 3] = (1 << 2);
					if (j == 2)
						navDTris[t + 3] |= (1 << 0);
					if (j == nv - 1)
						navDTris[t + 3] |= (1 << 4);
					tbase++;
				}
			}
		}

		// Store and create BVtree.
		// TODO: take detail mesh into account! use byte per bbox extent?
		if (params.buildBvTree) {
			// Do not set header.bvNodeCount set to make it work look exactly the same as in original Detour  
			header.bvNodeCount = NavMeshBuilder.createBVTree(params, navBvtree);
		}

		// Store Off-Mesh connections.
		n = 0;
		for (let i = 0; i < params.offMeshConCount; ++i) {
			// Only store connections which start from this tile.
			if (offMeshConClass[i * 2 + 0] == 0xff) {
				let con = new OffMeshConnection();
				offMeshCons[n] = con;
				con.poly = (offMeshPolyBase + n);
				// Copy connection end-points.
				let endPts = i * 2 * 3;
				arraycopy(params.offMeshConVerts, endPts, con.pos, 0, 6);
				con.rad = params.offMeshConRad[i];
				con.flags = params.offMeshConDir[i] != 0 ? NavMesh.DT_OFFMESH_CON_BIDIR : 0;
				con.side = offMeshConClass[i * 2 + 1];
				if (params.offMeshConUserID != null)
					con.userId = params.offMeshConUserID[i];
				n++;
			}
		}

		let nmd = new MeshData();
		nmd.header = header;
		nmd.verts = navVerts;
		nmd.polys = navPolys;
		nmd.detailMeshes = navDMeshes;
		nmd.detailVerts = navDVerts;
		nmd.detailTris = navDTris;
		nmd.bvTree = navBvtree;
		nmd.offMeshCons = offMeshCons;
		return nmd;
	}

}

export default NavMeshBuilder;
