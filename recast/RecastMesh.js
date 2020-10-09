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

import PolyMesh from "./PolyMesh.js"
import RecastVectors from "./RecastVectors.js"
import RecastConstants from "./RecastConstants.js"


function arraycopy(one, oneStart, two, twoStart, len) {
	for (let i = 0; i < len; i++) {
		two[twoStart + i] = one[oneStart + i];
	}
}

class Edge {
	vert = new Array(2);
	polyEdge = new Array(2);
	poly = new Array(2);

}

class RecastMesh {

	static VERTEX_BUCKET_COUNT = (1 << 12);



	static buildMeshAdjacency(polys, npolys, nverts, vertsPerPoly) {
		// Based on code by Eric Lengyel from:
		// http://www.terathon.com/code/edges.php

		let maxEdgeCount = npolys * vertsPerPoly;
		let firstEdge = new Array(nverts + maxEdgeCount);
		let nextEdge = nverts;
		let edgeCount = 0;

		let edges = new Array(maxEdgeCount);

		for (let i = 0; i < nverts; i++)
			firstEdge[i] = RecastConstants.RC_MESH_NULL_IDX;

		for (let i = 0; i < npolys; ++i) {
			let t = i * vertsPerPoly * 2;
			for (let j = 0; j < vertsPerPoly; ++j) {
				if (polys[t + j] == RecastConstants.RC_MESH_NULL_IDX)
					break;
				let v0 = polys[t + j];
				let v1 = (j + 1 >= vertsPerPoly || polys[t + j + 1] == RecastConstants.RC_MESH_NULL_IDX) ? polys[t + 0]
					: polys[t + j + 1];
				if (v0 < v1) {
					let edge = new Edge();
					edges[edgeCount] = edge;
					edge.vert[0] = v0;
					edge.vert[1] = v1;
					edge.poly[0] = i;
					edge.polyEdge[0] = j;
					edge.poly[1] = i;
					edge.polyEdge[1] = 0;
					// Insert edge
					firstEdge[nextEdge + edgeCount] = firstEdge[v0];
					firstEdge[v0] = edgeCount;
					edgeCount++;
				}
			}
		}

		for (let i = 0; i < npolys; ++i) {
			let t = i * vertsPerPoly * 2;
			for (let j = 0; j < vertsPerPoly; ++j) {
				if (polys[t + j] == RecastConstants.RC_MESH_NULL_IDX)
					break;
				let v0 = polys[t + j];
				let v1 = (j + 1 >= vertsPerPoly || polys[t + j + 1] == RecastConstants.RC_MESH_NULL_IDX) ? polys[t + 0]
					: polys[t + j + 1];
				if (v0 > v1) {
					for (let e = firstEdge[v1]; e != RecastConstants.RC_MESH_NULL_IDX; e = firstEdge[nextEdge + e]) {
						let edge = edges[e];
						if (edge.vert[1] == v0 && edge.poly[0] == edge.poly[1]) {
							edge.poly[1] = i;
							edge.polyEdge[1] = j;
							break;
						}
					}
				}
			}
		}

		// Store adjacency
		for (let i = 0; i < edgeCount; ++i) {
			let e = edges[i];
			if (e.poly[0] != e.poly[1]) {
				let p0 = e.poly[0] * vertsPerPoly * 2;
				let p1 = e.poly[1] * vertsPerPoly * 2;
				polys[p0 + vertsPerPoly + e.polyEdge[0]] = e.poly[1];
				polys[p1 + vertsPerPoly + e.polyEdge[1]] = e.poly[0];
			}
		}

	}

	static computeVertexHash(x, y, z) {
		let h1 = 0x8da6b343; // Large multiplicative constants;
		let h2 = 0xd8163841; // here arbitrarily chosen primes
		let h3 = 0xcb1ab31;
		let n = h1 * x + h2 * y + h3 * z;
		return n & (RecastMesh.VERTEX_BUCKET_COUNT - 1);
	}

	static addVertex(x, y, z, verts, firstVert, nextVert, nv) {
		let bucket = RecastMesh.computeVertexHash(x, 0, z);
		let i = firstVert[bucket];

		while (i != -1) {
			let v = i * 3;
			if (verts[v + 0] == x && (Math.abs(verts[v + 1] - y) <= 2) && verts[v + 2] == z)
				return [i, nv];
			i = nextVert[i]; // next
		}

		// Could not find, create new.
		i = nv;
		nv++;
		let v = i * 3;
		verts[v + 0] = x;
		verts[v + 1] = y;
		verts[v + 2] = z;
		nextVert[i] = firstVert[bucket];
		firstVert[bucket] = i;

		return [i, nv];
	}

	static prev(i, n) {
		return i - 1 >= 0 ? i - 1 : n - 1;
	}

	static next(i, n) {
		return i + 1 < n ? i + 1 : 0;
	}

	static area2(verts, a, b, c) {
		return (verts[b + 0] - verts[a + 0]) * (verts[c + 2] - verts[a + 2])
			- (verts[c + 0] - verts[a + 0]) * (verts[b + 2] - verts[a + 2]);
	}

	// Returns true iff c is strictly to the left of the directed
	// line through a to b.
	static left(verts, a, b, c) {
		return RecastMesh.area2(verts, a, b, c) < 0;
	}

	static leftOn(verts, a, b, c) {
		return RecastMesh.area2(verts, a, b, c) <= 0;
	}

	static collinear(verts, a, b, c) {
		return RecastMesh.area2(verts, a, b, c) == 0;
	}

	// Returns true iff ab properly intersects cd: they share
	// a poPoly interior to both segments. The properness of the
	// intersection is ensured by using strict leftness.
	static intersectProp(verts, a, b, c, d) {
		// Eliminate improper cases.
		if (RecastMesh.collinear(verts, a, b, c) || RecastMesh.collinear(verts, a, b, d) || RecastMesh.collinear(verts, c, d, a)
			|| RecastMesh.collinear(verts, c, d, b))
			return false;

		return (RecastMesh.left(verts, a, b, c) ^ RecastMesh.left(verts, a, b, d)) && (RecastMesh.left(verts, c, d, a) ^ RecastMesh.left(verts, c, d, b));
	}

	// Returns T iff (a,b,c) are collinear and poPoly c lies
	// on the closed segement ab.
	static between(verts, a, b, c) {
		if (!RecastMesh.collinear(verts, a, b, c))
			return false;
		// If ab not vertical, check betweenness on x; else on y.
		if (verts[a + 0] != verts[b + 0])
			return ((verts[a + 0] <= verts[c + 0]) && (verts[c + 0] <= verts[b + 0]))
				|| ((verts[a + 0] >= verts[c + 0]) && (verts[c + 0] >= verts[b + 0]));
		else
			return ((verts[a + 2] <= verts[c + 2]) && (verts[c + 2] <= verts[b + 2]))
				|| ((verts[a + 2] >= verts[c + 2]) && (verts[c + 2] >= verts[b + 2]));
	}

	// Returns true iff segments ab and cd intersect, properly or improperly.
	static intersect(verts, a, b, c, d) {
		if (RecastMesh.intersectProp(verts, a, b, c, d))
			return true;
		else if (RecastMesh.between(verts, a, b, c) || RecastMesh.between(verts, a, b, d) || RecastMesh.between(verts, c, d, a)
			|| RecastMesh.between(verts, c, d, b))
			return true;
		else
			return false;
	}

	static vequal(verts, a, b) {
		return verts[a + 0] == verts[b + 0] && verts[a + 2] == verts[b + 2];
	}

	// Returns T iff (v_i, v_j) is a proper internal *or* external
	// diagonal of P, *ignoring edges incident to v_i and v_j*.
	static diagonalie(i, j, n, verts, indices) {
		let d0 = (indices[i] & 0x0fffffff) * 4;
		let d1 = (indices[j] & 0x0fffffff) * 4;

		// For each edge (k,k+1) of P
		for (let k = 0; k < n; k++) {
			let k1 = RecastMesh.next(k, n);
			// Skip edges incident to i or j
			if (!((k == i) || (k1 == i) || (k == j) || (k1 == j))) {
				let p0 = (indices[k] & 0x0fffffff) * 4;
				let p1 = (indices[k1] & 0x0fffffff) * 4;

				if (RecastMesh.vequal(verts, d0, p0) || RecastMesh.vequal(verts, d1, p0) || RecastMesh.vequal(verts, d0, p1) || RecastMesh.vequal(verts, d1, p1))
					continue;

				if (RecastMesh.intersect(verts, d0, d1, p0, p1))
					return false;
			}
		}
		return true;
	}

	// Returns true iff the diagonal (i,j) is strictly internal to the
	// polygon P in the neighborhood of the i endpoint.
	static inCone(i, j, n, verts, indices) {
		let pi = (indices[i] & 0x0fffffff) * 4;
		let pj = (indices[j] & 0x0fffffff) * 4;
		let pi1 = (indices[RecastMesh.next(i, n)] & 0x0fffffff) * 4;
		let pin1 = (indices[RecastMesh.prev(i, n)] & 0x0fffffff) * 4;
		// If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].
		if (RecastMesh.leftOn(verts, pin1, pi, pi1)) {
			return RecastMesh.left(verts, pi, pj, pin1) && RecastMesh.left(verts, pj, pi, pi1);
		}
		// Assume (i-1,i,i+1) not collinear.
		// else P[i] is reflex.
		return !(RecastMesh.leftOn(verts, pi, pj, pi1) && RecastMesh.leftOn(verts, pj, pi, pin1));
	}

	// Returns T iff (v_i, v_j) is a proper internal
	// diagonal of P.
	static diagonal(i, j, n, verts, indices) {
		return RecastMesh.inCone(i, j, n, verts, indices) && RecastMesh.diagonalie(i, j, n, verts, indices);
	}

	static diagonalieLoose(i, j, n, verts, indices) {
		let d0 = (indices[i] & 0x0fffffff) * 4;
		let d1 = (indices[j] & 0x0fffffff) * 4;

		// For each edge (k,k+1) of P
		for (let k = 0; k < n; k++) {
			let k1 = RecastMesh.next(k, n);
			// Skip edges incident to i or j
			if (!((k == i) || (k1 == i) || (k == j) || (k1 == j))) {
				let p0 = (indices[k] & 0x0fffffff) * 4;
				let p1 = (indices[k1] & 0x0fffffff) * 4;

				if (RecastMesh.vequal(verts, d0, p0) || RecastMesh.vequal(verts, d1, p0) || RecastMesh.vequal(verts, d0, p1) || RecastMesh.vequal(verts, d1, p1))
					continue;

				if (RecastMesh.intersectProp(verts, d0, d1, p0, p1))
					return false;
			}
		}
		return true;
	}

	static inConeLoose(i, j, n, verts, indices) {
		let pi = (indices[i] & 0x0fffffff) * 4;
		let pj = (indices[j] & 0x0fffffff) * 4;
		let pi1 = (indices[RecastMesh.next(i, n)] & 0x0fffffff) * 4;
		let pin1 = (indices[RecastMesh.prev(i, n)] & 0x0fffffff) * 4;

		// If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].
		if (RecastMesh.leftOn(verts, pin1, pi, pi1))
			return RecastMesh.leftOn(verts, pi, pj, pin1) && RecastMesh.leftOn(verts, pj, pi, pi1);
		// Assume (i-1,i,i+1) not collinear.
		// else P[i] is reflex.
		return !(RecastMesh.leftOn(verts, pi, pj, pi1) && RecastMesh.leftOn(verts, pj, pi, pin1));
	}

	static diagonalLoose(i, j, n, verts, indices) {
		return RecastMesh.inConeLoose(i, j, n, verts, indices) && RecastMesh.diagonalieLoose(i, j, n, verts, indices);
	}

	static triangulate(n, verts, indices, tris) {
		let ntris = 0;

		// The last bit of the index is used to indicate if the vertex can be removed.
		for (let i = 0; i < n; i++) {
			let i1 = RecastMesh.next(i, n);
			let i2 = RecastMesh.next(i1, n);
			if (RecastMesh.diagonal(i, i2, n, verts, indices)) {
				indices[i1] |= 0x80000000;
			}
		}

		while (n > 3) {
			let minLen = -1;
			let mini = -1;
			for (let i = 0; i < n; i++) {
				let i1 = RecastMesh.next(i, n);
				if ((indices[i1] & 0x80000000) != 0) {
					let p0 = (indices[i] & 0x0fffffff) * 4;
					let p2 = (indices[RecastMesh.next(i1, n)] & 0x0fffffff) * 4;

					let dx = verts[p2 + 0] - verts[p0 + 0];
					let dy = verts[p2 + 2] - verts[p0 + 2];
					let len = dx * dx + dy * dy;

					if (minLen < 0 || len < minLen) {
						minLen = len;
						mini = i;
					}
				}
			}

			if (mini == -1) {
				// We might get here because the contour has overlapping segments, like this:
				//
				// A o-o=====o---o B
				// / |C D| \
				// o o o o
				// : : : :
				// We'll try to recover by loosing up the inCone test a bit so that a diagonal
				// like A-B or C-D can be found and we can continue.
				minLen = -1;
				mini = -1;
				for (let i = 0; i < n; i++) {
					let i1 = RecastMesh.next(i, n);
					let i2 = RecastMesh.next(i1, n);
					if (RecastMesh.diagonalLoose(i, i2, n, verts, indices)) {
						let p0 = (indices[i] & 0x0fffffff) * 4;
						let p2 = (indices[RecastMesh.next(i2, n)] & 0x0fffffff) * 4;
						let dx = verts[p2 + 0] - verts[p0 + 0];
						let dy = verts[p2 + 2] - verts[p0 + 2];
						let len = dx * dx + dy * dy;

						if (minLen < 0 || len < minLen) {
							minLen = len;
							mini = i;
						}
					}
				}
				if (mini == -1) {
					// The contour is messed up. This sometimes happens
					// if the contour simplification is too aggressive.
					return -ntris;
				}
			}

			let i = mini;
			let i1 = RecastMesh.next(i, n);
			let i2 = RecastMesh.next(i1, n);

			tris[ntris * 3] = indices[i] & 0x0ffffff;
			tris[ntris * 3 + 1] = indices[i1] & 0x0ffffff;
			tris[ntris * 3 + 2] = indices[i2] & 0x0ffffff;
			ntris++;

			// Removes P[i1] by copying P[i+1]...P[n-1] left one index.
			n--;
			for (let k = i1; k < n; k++)
				indices[k] = indices[k + 1];

			if (i1 >= n)
				i1 = 0;
			i = RecastMesh.prev(i1, n);
			// Update diagonal flags.
			if (RecastMesh.diagonal(RecastMesh.prev(i, n), i1, n, verts, indices))
				indices[i] |= 0x80000000;
			else
				indices[i] &= 0x0ffffff;

			if (RecastMesh.diagonal(i, RecastMesh.next(i1, n), n, verts, indices))
				indices[i1] |= 0x80000000;
			else
				indices[i1] &= 0x0ffffff;
		}

		// Append the remaining triangle.
		tris[ntris * 3] = indices[0] & 0x0ffffff;
		tris[ntris * 3 + 1] = indices[1] & 0x0ffffff;
		tris[ntris * 3 + 2] = indices[2] & 0x0ffffff;
		ntris++;

		return ntris;
	}

	static countPolyVerts(p, j, nvp) {
		for (let i = 0; i < nvp; ++i)
			if (p[i + j] == RecastConstants.RC_MESH_NULL_IDX)
				return i;
		return nvp;
	}

	static uleft(verts, a, b, c) {
		return (verts[b + 0] - verts[a + 0]) * (verts[c + 2] - verts[a + 2])
			- (verts[c + 0] - verts[a + 0]) * (verts[b + 2] - verts[a + 2]) < 0;
	}

	static getPolyMergeValue(polys, pa, pb, verts, nvp) {
		let ea = -1;
		let eb = -1;
		let na = RecastMesh.countPolyVerts(polys, pa, nvp);
		let nb = RecastMesh.countPolyVerts(polys, pb, nvp);

		// If the merged polygon would be too big, do not merge.
		if (na + nb - 2 > nvp)
			return [-1, ea, eb];

		// Check if the polygons share an edge.

		for (let i = 0; i < na; ++i) {
			let va0 = polys[pa + i];
			let va1 = polys[pa + (i + 1) % na];
			if (va0 > va1) {
				let temp = va0;
				va0 = va1;
				va1 = temp;
			}
			for (let j = 0; j < nb; ++j) {
				let vb0 = polys[pb + j];
				let vb1 = polys[pb + (j + 1) % nb];
				if (vb0 > vb1) {
					let temp = vb0;
					vb0 = vb1;
					vb1 = temp;
				}
				if (va0 == vb0 && va1 == vb1) {
					ea = i;
					eb = j;
					break;
				}
			}
		}

		// No common edge, cannot merge.
		if (ea == -1 || eb == -1)
			return [-1, ea, eb];

		// Check to see if the merged polygon would be convex.
		let va, vb, vc;

		va = polys[pa + (ea + na - 1) % na];
		vb = polys[pa + ea];
		vc = polys[pb + (eb + 2) % nb];
		if (!RecastMesh.uleft(verts, va * 3, vb * 3, vc * 3))
			return [-1, ea, eb];

		va = polys[pb + (eb + nb - 1) % nb];
		vb = polys[pb + eb];
		vc = polys[pa + (ea + 2) % na];
		if (!RecastMesh.uleft(verts, va * 3, vb * 3, vc * 3))
			return [-1, ea, eb];

		va = polys[pa + ea];
		vb = polys[pa + (ea + 1) % na];

		let dx = verts[va * 3 + 0] - verts[vb * 3 + 0];
		let dy = verts[va * 3 + 2] - verts[vb * 3 + 2];

		return [dx * dx + dy * dy, ea, eb];
	}

	static mergePolyVerts(polys, pa, pb, ea, eb, tmp, nvp) {
		let na = RecastMesh.countPolyVerts(polys, pa, nvp);
		let nb = RecastMesh.countPolyVerts(polys, pb, nvp);

		// Merge polygons.
		polys.fill(RecastConstants.RC_MESH_NULL_IDX, tmp, tmp + nvp, );
		let n = 0;
		// Add pa
		for (let i = 0; i < na - 1; ++i) {
			polys[tmp + n] = polys[pa + (ea + 1 + i) % na];
			n++;
		}
		// Add pb
		for (let i = 0; i < nb - 1; ++i) {
			polys[tmp + n] = polys[pb + (eb + 1 + i) % nb];
			n++;
		}
		//arraycopy(polys, tmp, polys, pa, nvp);
		arraycopy(polys, tmp, polys, pa, nvp)
		// for (let i = 0; i < nvp; i++) {
		// 	polys[pa + i] = polys[tmp + i];
		// }
	}

	static pushFront(v, arr, an) {
		an++;
		for (let i = an - 1; i > 0; --i)
			arr[i] = arr[i - 1];
		arr[0] = v;
		return an;
	}

	static pushBack(v, arr, an) {
		arr[an] = v;
		an++;
		return an;
	}

	static canRemoveVertex(ctx, mesh, rem) {
		let nvp = mesh.nvp;

		// Count number of polygons to remove.
		let numTouchedVerts = 0;
		let numRemainingEdges = 0;
		for (let i = 0; i < mesh.npolys; ++i) {
			let p = i * nvp * 2;
			let nv = RecastMesh.countPolyVerts(mesh.polys, p, nvp);
			let numRemoved = 0;
			let numVerts = 0;
			for (let j = 0; j < nv; ++j) {
				if (mesh.polys[p + j] == rem) {
					numTouchedVerts++;
					numRemoved++;
				}
				numVerts++;
			}
			if (numRemoved != 0) {
				numRemainingEdges += numVerts - (numRemoved + 1);
			}
		}
		// There would be too few edges remaining to create a polygon.
		// This can happen for example when a tip of a triangle is marked
		// as dePolyion, but there are no other polys that share the vertex.
		// In this case, the vertex should not be removed.
		if (numRemainingEdges <= 2)
			return false;

		// Find edges which share the removed vertex.
		let maxEdges = numTouchedVerts * 2;
		let nedges = 0;
		let edges = new Array(maxEdges * 3);

		for (let i = 0; i < mesh.npolys; ++i) {
			let p = i * nvp * 2;
			let nv = RecastMesh.countPolyVerts(mesh.polys, p, nvp);

			// Collect edges which touches the removed vertex.
			for (let j = 0, k = nv - 1; j < nv; k = j++) {
				if (mesh.polys[p + j] == rem || mesh.polys[p + k] == rem) {
					// Arrange edge so that a=rem.
					let a = mesh.polys[p + j], b = mesh.polys[p + k];
					if (b == rem) {
						let temp = a;
						a = b;
						b = temp;
					}
					// Check if the edge exists
					let exists = false;
					for (let m = 0; m < nedges; ++m) {
						let e = m * 3;
						if (edges[e + 1] == b) {
							// Exists, increment vertex share count.
							edges[e + 2]++;
							exists = true;
						}
					}
					// Add new edge.
					if (!exists) {
						let e = nedges * 3;
						edges[e + 0] = a;
						edges[e + 1] = b;
						edges[e + 2] = 1;
						nedges++;
					}
				}
			}
		}

		// There should be no more than 2 open edges.
		// This catches the case that two non-adjacent polygons
		// share the removed vertex. In that case, do not remove the vertex.
		let numOpenEdges = 0;
		for (let i = 0; i < nedges; ++i) {
			if (edges[i * 3 + 2] < 2)
				numOpenEdges++;
		}
		if (numOpenEdges > 2)
			return false;

		return true;
	}

	static removeVertex(ctx, mesh, rem, maxTris) {
		let nvp = mesh.nvp;

		// Count number of polygons to remove.
		let numRemovedVerts = 0;
		for (let i = 0; i < mesh.npolys; ++i) {
			let p = i * nvp * 2;
			let nv = RecastMesh.countPolyVerts(mesh.polys, p, nvp);
			for (let j = 0; j < nv; ++j) {
				if (mesh.polys[p + j] == rem)
					numRemovedVerts++;
			}
		}

		let nedges = 0;
		let edges = new Array(numRemovedVerts * nvp * 4);

		let nhole = 0;
		let hole = new Array(numRemovedVerts * nvp);

		let nhreg = 0;
		let hreg = new Array(numRemovedVerts * nvp);

		let nharea = 0;
		let harea = new Array(numRemovedVerts * nvp);

		for (let i = 0; i < mesh.npolys; ++i) {
			let p = i * nvp * 2;
			let nv = RecastMesh.countPolyVerts(mesh.polys, p, nvp);
			let hasRem = false;
			for (let j = 0; j < nv; ++j)
				if (mesh.polys[p + j] == rem)
					hasRem = true;
			if (hasRem) {
				// Collect edges which does not touch the removed vertex.
				for (let j = 0, k = nv - 1; j < nv; k = j++) {
					if (mesh.polys[p + j] != rem && mesh.polys[p + k] != rem) {
						let e = nedges * 4;
						edges[e + 0] = mesh.polys[p + k];
						edges[e + 1] = mesh.polys[p + j];
						edges[e + 2] = mesh.regs[i];
						edges[e + 3] = mesh.areas[i];
						nedges++;
					}
				}
				// Remove the polygon.
				let p2 = (mesh.npolys - 1) * nvp * 2;
				if (p != p2) {
					arraycopy(mesh.polys, p2, mesh.polys, p, nvp);
				}
				mesh.polys.fill( RecastConstants.RC_MESH_NULL_IDX, p + nvp, p + nvp + nvp);
				mesh.regs[i] = mesh.regs[mesh.npolys - 1];
				mesh.areas[i] = mesh.areas[mesh.npolys - 1];
				mesh.npolys--;
				--i;
			}
		}

		// Remove vertex.
		for (let i = rem; i < mesh.nverts - 1; ++i) {
			mesh.verts[i * 3 + 0] = mesh.verts[(i + 1) * 3 + 0];
			mesh.verts[i * 3 + 1] = mesh.verts[(i + 1) * 3 + 1];
			mesh.verts[i * 3 + 2] = mesh.verts[(i + 1) * 3 + 2];
		}
		mesh.nverts--;

		// Adjust indices to match the removed vertex layout.
		for (let i = 0; i < mesh.npolys; ++i) {
			let p = i * nvp * 2;
			let nv = RecastMesh.countPolyVerts(mesh.polys, p, nvp);
			for (let j = 0; j < nv; ++j)
				if (mesh.polys[p + j] > rem)
					mesh.polys[p + j]--;
		}
		for (let i = 0; i < nedges; ++i) {
			if (edges[i * 4 + 0] > rem)
				edges[i * 4 + 0]--;
			if (edges[i * 4 + 1] > rem)
				edges[i * 4 + 1]--;
		}

		if (nedges == 0)
			return;

		// Start with one vertex, keep appending connected
		// segments to the start and end of the hole.
		RecastMesh.pushBack(edges[0], hole, nhole);
		RecastMesh.pushBack(edges[2], hreg, nhreg);
		RecastMesh.pushBack(edges[3], harea, nharea);

		while (nedges != 0) {
			let match = false;

			for (let i = 0; i < nedges; ++i) {
				let ea = edges[i * 4 + 0];
				let eb = edges[i * 4 + 1];
				let r = edges[i * 4 + 2];
				let a = edges[i * 4 + 3];
				let add = false;
				if (hole[0] == eb) {
					// The segment matches the beginning of the hole boundary.
					RecastMesh.pushFront(ea, hole, nhole);
					RecastMesh.pushFront(r, hreg, nhreg);
					RecastMesh.pushFront(a, harea, nharea);
					add = true;
				} else if (hole[nhole - 1] == ea) {
					// The segment matches the end of the hole boundary.
					nhole = RecastMesh.pushBack(eb, hole, nhole);
					nhreg = RecastMesh.pushBack(r, hreg, nhreg);
					nharea = RecastMesh.pushBack(a, harea, nharea);
					add = true;
				}
				if (add) {
					// The edge segment was added, remove it.
					edges[i * 4 + 0] = edges[(nedges - 1) * 4 + 0];
					edges[i * 4 + 1] = edges[(nedges - 1) * 4 + 1];
					edges[i * 4 + 2] = edges[(nedges - 1) * 4 + 2];
					edges[i * 4 + 3] = edges[(nedges - 1) * 4 + 3];
					--nedges;
					match = true;
					--i;
				}
			}

			if (!match)
				break;
		}

		let tris = new Array(nhole * 3);

		let tverts = new Array(nhole * 4);

		let thole = new Array(nhole);

		// Generate temp vertex array for triangulation.
		for (let i = 0; i < nhole; ++i) {
			let pi = hole[i];
			tverts[i * 4 + 0] = mesh.verts[pi * 3 + 0];
			tverts[i * 4 + 1] = mesh.verts[pi * 3 + 1];
			tverts[i * 4 + 2] = mesh.verts[pi * 3 + 2];
			tverts[i * 4 + 3] = 0;
			thole[i] = i;
		}

		// Triangulate the hole.
		let ntris = RecastMesh.triangulate(nhole, tverts, thole, tris);
		if (ntris < 0) {
			ntris = -ntris;
			ctx.warn("removeVertex: triangulate() returned bad results.");
		}

		// Merge the hole triangles back to polygons.
		let polys = new Array((ntris + 1) * nvp);
		let pregs = new Array(ntris);
		let pareas = new Array(ntris);

		let tpmPoly = ntris * nvp;

		// Build initial polygons.
		let npolys = 0;
		polys.fill(RecastConstants.RC_MESH_NULL_IDX, 0, ntris * nvp);
		for (let j = 0; j < ntris; ++j) {
			let t = j * 3;
			if (tris[t + 0] != tris[t + 1] && tris[t + 0] != tris[t + 2] && tris[t + 1] != tris[t + 2]) {
				polys[npolys * nvp + 0] = hole[tris[t + 0]];
				polys[npolys * nvp + 1] = hole[tris[t + 1]];
				polys[npolys * nvp + 2] = hole[tris[t + 2]];

				// If this polygon covers multiple region types then
				// mark it as such
				if (hreg[tris[t + 0]] != hreg[tris[t + 1]] || hreg[tris[t + 1]] != hreg[tris[t + 2]])
					pregs[npolys] = RecastConstants.RC_MULTIPLE_REGS;
				else
					pregs[npolys] = hreg[tris[t + 0]];

				pareas[npolys] = harea[tris[t + 0]];
				npolys++;
			}
		}
		if (npolys == 0)
			return;

		// Merge polygons.
		if (nvp > 3) {
			for (; ;) {
				// Find best polygons to merge.
				let bestMergeVal = 0;
				let bestPa = 0, bestPb = 0, bestEa = 0, bestEb = 0;

				for (let j = 0; j < npolys - 1; ++j) {
					let pj = j * nvp;
					for (let k = j + 1; k < npolys; ++k) {
						let pk = k * nvp;
						let veaeb = RecastMesh.getPolyMergeValue(polys, pj, pk, mesh.verts, nvp);
						let v = veaeb[0];
						let ea = veaeb[1];
						let eb = veaeb[2];
						if (v > bestMergeVal) {
							bestMergeVal = v;
							bestPa = j;
							bestPb = k;
							bestEa = ea;
							bestEb = eb;
						}
					}
				}

				if (bestMergeVal > 0) {
					// Found best, merge.
					let pa = bestPa * nvp;
					let pb = bestPb * nvp;
					RecastMesh.mergePolyVerts(polys, pa, pb, bestEa, bestEb, tmpPoly, nvp);
					if (pregs[bestPa] != pregs[bestPb])
						pregs[bestPa] = RecastConstants.RC_MULTIPLE_REGS;
					let last = (npolys - 1) * nvp;
					if (pb != last) {
						arraycopy(polys, last, polys, pb, nvp);
					}
					pregs[bestPb] = pregs[npolys - 1];
					pareas[bestPb] = pareas[npolys - 1];
					npolys--;
				} else {
					// Could not merge any polygons, stop.
					break;
				}
			}
		}

		// Store polygons.
		for (let i = 0; i < npolys; ++i) {
			if (mesh.npolys >= maxTris)
				break;
			let p = mesh.npolys * nvp * 2;
			mesh.polys.fill(RecastConstants.RC_MESH_NULL_IDX, p, p + nvp * 2 );
			for (let j = 0; j < nvp; ++j)
				mesh.polys[p + j] = polys[i * nvp + j];
			mesh.regs[mesh.npolys] = pregs[i];
			mesh.areas[mesh.npolys] = pareas[i];
			mesh.npolys++;
			if (mesh.npolys > maxTris) {
				throw new RuntimeException("removeVertex: Too many polygons " + mesh.npolys + " (max:" + maxTris + ".");
			}
		}

	}

	/// @par
	///
	/// @note If the mesh data is to be used to construct a Detour navigation mesh, then the upper
	/// limit must be retricted to <= #DT_VERTS_PER_POLYGON.
	///
	/// @see rcAllocPolyMesh, rcContourSet, rcPolyMesh, rcConfig
	static buildPolyMesh(ctx, cset, nvp) {
		ctx.startTimer("BUILD_POLYMESH");
		let mesh = new PolyMesh();
		RecastVectors.copy3(mesh.bmin, cset.bmin, 0);
		RecastVectors.copy3(mesh.bmax, cset.bmax, 0);
		mesh.cs = cset.cs;
		mesh.ch = cset.ch;
		mesh.borderSize = cset.borderSize;
		mesh.maxEdgeError = cset.maxError;

		let maxVertices = 0;
		let maxTris = 0;
		let maxVertsPerCont = 0;
		for (let i = 0; i < cset.conts.length; ++i) {
			// Skip null contours.
			if (cset.conts[i].nverts < 3)
				continue;
			maxVertices += cset.conts[i].nverts;
			maxTris += cset.conts[i].nverts - 2;
			maxVertsPerCont = Math.max(maxVertsPerCont, cset.conts[i].nverts);
		}
		if (maxVertices >= 0xfffe) {
			throw new RuntimeException("rcBuildPolyMesh: Too many vertices " + maxVertices);
		}
		let vflags = new Array(maxVertices).fill(0);
		vflags.fill(0);

		mesh.verts = new Array(maxVertices * 3).fill(0);
		mesh.polys = new Array(maxTris * nvp * 2).fill(RecastConstants.RC_MESH_NULL_IDX);

		// Arrays.fill(mesh.polys, RecastConstants.RC_MESH_NULL_IDX);
		mesh.regs = new Array(maxTris).fill(0);
		mesh.areas = new Array(maxTris).fill(0);

		mesh.nverts = 0;
		mesh.npolys = 0;
		mesh.nvp = nvp;
		mesh.maxpolys = maxTris;

		let nextVert = new Array(maxVertices);

		let firstVert = new Array(RecastMesh.VERTEX_BUCKET_COUNT);
		for (let i = 0; i < RecastMesh.VERTEX_BUCKET_COUNT; ++i)
			firstVert[i] = -1;

		let indices = new Array(maxVertsPerCont);
		let tris = new Array(maxVertsPerCont * 3);
		let polys = new Array((maxVertsPerCont + 1) * nvp);

		let tmpPoly = maxVertsPerCont * nvp;

		for (let i = 0; i < cset.conts.length; ++i) {
			let cont = cset.conts[i];

			// Skip null contours.
			if (cont.nverts < 3)
				continue;

			// Triangulate contour
			for (let j = 0; j < cont.nverts; ++j)
				indices[j] = j;
			let ntris = RecastMesh.triangulate(cont.nverts, cont.verts, indices, tris);
			if (ntris <= 0) {
				// Bad triangulation, should not happen.
				ctx.warn("buildPolyMesh: Bad triangulation Contour " + i + ".");
				ntris = -ntris;
			}

			// Add and merge vertices.
			for (let j = 0; j < cont.nverts; ++j) {
				let v = j * 4;
				let inv = RecastMesh.addVertex(cont.verts[v + 0], cont.verts[v + 1], cont.verts[v + 2], mesh.verts, firstVert,
					nextVert, mesh.nverts);
				indices[j] = inv[0];
				mesh.nverts = inv[1];
				if ((cont.verts[v + 3] & RecastConstants.RC_BORDER_VERTEX) != 0) {
					// This vertex should be removed.
					vflags[indices[j]] = 1;
				}
			}

			// Build initial polygons.
			let npolys = 0;
			// Arrays.fill(polys, RecastConstants.RC_MESH_NULL_IDX);
			polys.fill(RecastConstants.RC_MESH_NULL_IDX);
			for (let j = 0; j < ntris; ++j) {
				let t = j * 3;
				if (tris[t + 0] != tris[t + 1] && tris[t + 0] != tris[t + 2] && tris[t + 1] != tris[t + 2]) {
					polys[npolys * nvp + 0] = indices[tris[t + 0]];
					polys[npolys * nvp + 1] = indices[tris[t + 1]];
					polys[npolys * nvp + 2] = indices[tris[t + 2]];
					npolys++;
				}
			}
			if (npolys == 0)
				continue;

			// Merge polygons.
			if (nvp > 3) {
				for (; ;) {
					// Find best polygons to merge.
					let bestMergeVal = 0;
					let bestPa = 0, bestPb = 0, bestEa = 0, bestEb = 0;

					for (let j = 0; j < npolys - 1; ++j) {
						let pj = j * nvp;
						for (let k = j + 1; k < npolys; ++k) {
							let pk = k * nvp;
							let veaeb = RecastMesh.getPolyMergeValue(polys, pj, pk, mesh.verts, nvp);
							let v = veaeb[0];
							let ea = veaeb[1];
							let eb = veaeb[2];
							if (v > bestMergeVal) {
								bestMergeVal = v;
								bestPa = j;
								bestPb = k;
								bestEa = ea;
								bestEb = eb;
							}
						}
					}

					if (bestMergeVal > 0) {
						// Found best, merge.
						let pa = bestPa * nvp;
						let pb = bestPb * nvp;
						RecastMesh.mergePolyVerts(polys, pa, pb, bestEa, bestEb, tmpPoly, nvp);
						let lastPoly = (npolys - 1) * nvp;
						if (pb != lastPoly) {
							// arraycopy(polys, lastPoly, polys, pb, nvp);
							for (let i = 0; i < nvp; i++) {
								polys[pb + i] = polys[lastPoly + i];
								if (polys[28] == 80)
									console.log("break")
							}
						}
						npolys--;
					} else {
						// Could not merge any polygons, stop.
						break;
					}
				}
			}

			//Here they are the same
			// Store polygons.
			for (let j = 0; j < npolys; ++j) {
				let p = mesh.npolys * nvp * 2;
				let q = j * nvp;
				for (let k = 0; k < nvp; ++k)
					mesh.polys[p + k] = polys[q + k];
				mesh.regs[mesh.npolys] = cont.reg;
				mesh.areas[mesh.npolys] = cont.area;
				mesh.npolys++;
				if (mesh.npolys > maxTris) {
					throw new RuntimeException(
						"rcBuildPolyMesh: Too many polygons " + mesh.npolys + " (max:" + maxTris + ").");
				}
			}
		}
		// fs.writeFileSync("./pmeshjs.txt", JSON.stringify(mesh.polys))


		//Now we are the same here
		// Remove edge vertices.
		for (let i = 0; i < mesh.nverts; ++i) {
			if (vflags[i] != 0) {
				if (!RecastMesh.canRemoveVertex(ctx, mesh, i))
					continue;
				RecastMesh.removeVertex(ctx, mesh, i, maxTris);
				// Remove vertex
				// Note: mesh.nverts is already decremented inside removeVertex()!
				// Fixup vertex flags
				for (let j = i; j < mesh.nverts; ++j)
					vflags[j] = vflags[j + 1];
				--i;
			}
		}

		// Calculate adjacency.
		this.buildMeshAdjacency(mesh.polys, mesh.npolys, mesh.nverts, nvp);

		// Find portal edges
		if (mesh.borderSize > 0) {
			let w = cset.width;
			let h = cset.height;
			for (let i = 0; i < mesh.npolys; ++i) {
				let p = i * 2 * nvp;
				for (let j = 0; j < nvp; ++j) {
					if (mesh.polys[p + j] == RecastConstants.RC_MESH_NULL_IDX)
						break;
					// Skip connected edges.
					if (mesh.polys[p + nvp + j] != RecastConstants.RC_MESH_NULL_IDX)
						continue;
					let nj = j + 1;
					if (nj >= nvp || mesh.polys[p + nj] == RecastConstants.RC_MESH_NULL_IDX)
						nj = 0;
					let va = mesh.polys[p + j] * 3;
					let vb = mesh.polys[p + nj] * 3;

					if (mesh.verts[va + 0] == 0 && mesh.verts[vb + 0] == 0)
						mesh.polys[p + nvp + j] = 0x8000 | 0;
					else if (mesh.verts[va + 2] == h && mesh.verts[vb + 2] == h)
						mesh.polys[p + nvp + j] = 0x8000 | 1;
					else if (mesh.verts[va + 0] == w && mesh.verts[vb + 0] == w)
						mesh.polys[p + nvp + j] = 0x8000 | 2;
					else if (mesh.verts[va + 2] == 0 && mesh.verts[vb + 2] == 0)
						mesh.polys[p + nvp + j] = 0x8000 | 3;
				}
			}
		}

		// Just allocate the mesh flags array. The user is resposible to fill it.
		mesh.flags = new Array(mesh.npolys);
		mesh.flags.fill(0);

		if (mesh.nverts > 0xffff) {
			throw new RuntimeException("rcBuildPolyMesh: The resulting mesh has too many vertices " + mesh.nverts
				+ " (max " + 0xffff + "). Data can be corrupted.");
		}
		if (mesh.npolys > 0xffff) {
			throw new RuntimeException("rcBuildPolyMesh: The resulting mesh has too many polygons " + mesh.npolys
				+ " (max " + 0xffff + "). Data can be corrupted.");
		}

		ctx.stopTimer("BUILD_POLYMESH");
		return mesh;

	}

	/// @see rcAllocPolyMesh, rcPolyMesh
	static mergePolyMeshes(ctx, meshes, nmeshes) {

		if (nmeshes == 0 || meshes == null)
			return null;

		ctx.startTimer("MERGE_POLYMESH");
		let mesh = new PolyMesh();
		mesh.nvp = meshes[0].nvp;
		mesh.cs = meshes[0].cs;
		mesh.ch = meshes[0].ch;
		RecastVectors.copy(mesh.bmin, meshes[0].bmin, 0);
		RecastVectors.copy(mesh.bmax, meshes[0].bmax, 0);

		let maxVerts = 0;
		let maxPolys = 0;
		let maxVertsPerMesh = 0;
		for (let i = 0; i < nmeshes; ++i) {
			RecastVectors.min(mesh.bmin, meshes[i].bmin, 0);
			RecastVectors.max(mesh.bmax, meshes[i].bmax, 0);
			maxVertsPerMesh = Math.max(maxVertsPerMesh, meshes[i].nverts);
			maxVerts += meshes[i].nverts;
			maxPolys += meshes[i].npolys;
		}

		mesh.nverts = 0;
		mesh.verts = new Array(maxVerts * 3);

		mesh.npolys = 0;
		mesh.polys = new Array(maxPolys * 2 * mesh.nvp);
		mesh.polys.fill(RecastConstants.RC_MESH_NULL_IDX, 0, mesh.polys.length);
		mesh.regs = new Array(maxPolys);
		mesh.areas = new Array(maxPolys);
		mesh.flags = new Array(maxPolys);

		let nextVert = new Array(maxVerts);

		let firstVert = new Array(RecastMesh.VERTEX_BUCKET_COUNT);
		for (let i = 0; i < RecastMesh.VERTEX_BUCKET_COUNT; ++i)
			firstVert[i] = -1;

		let vremap = new Array(maxVertsPerMesh);

		for (let i = 0; i < nmeshes; ++i) {
			let pmesh = meshes[i];

			let ox = Math.floor((pmesh.bmin[0] - mesh.bmin[0]) / mesh.cs + 0.5);
			let oz = Math.floor((pmesh.bmin[2] - mesh.bmin[2]) / mesh.cs + 0.5);

			let isMinX = (ox == 0);
			let isMinZ = (oz == 0);
			let isMaxX = (Math.floor((mesh.bmax[0] - pmesh.bmax[0]) / mesh.cs + 0.5)) == 0;
			let isMaxZ = (Math.floor((mesh.bmax[2] - pmesh.bmax[2]) / mesh.cs + 0.5)) == 0;
			let isOnBorder = (isMinX || isMinZ || isMaxX || isMaxZ);

			for (let j = 0; j < pmesh.nverts; ++j) {
				let v = j * 3;
				let inv = addVertex(pmesh.verts[v + 0] + ox, pmesh.verts[v + 1], pmesh.verts[v + 2] + oz, mesh.verts,
					firstVert, nextVert, mesh.nverts);

				vremap[j] = inv[0];
				mesh.nverts = inv[1];
			}

			for (let j = 0; j < pmesh.npolys; ++j) {
				let tgt = mesh.npolys * 2 * mesh.nvp;
				let src = j * 2 * mesh.nvp;
				mesh.regs[mesh.npolys] = pmesh.regs[j];
				mesh.areas[mesh.npolys] = pmesh.areas[j];
				mesh.flags[mesh.npolys] = pmesh.flags[j];
				mesh.npolys++;
				for (let k = 0; k < mesh.nvp; ++k) {
					if (pmesh.polys[src + k] == RecastConstants.RC_MESH_NULL_IDX)
						break;
					mesh.polys[tgt + k] = vremap[pmesh.polys[src + k]];
				}

				if (isOnBorder) {
					for (let k = mesh.nvp; k < mesh.nvp * 2; ++k) {
						if ((pmesh.polys[src + k] & 0x8000) != 0 && pmesh.polys[src + k] != 0xffff) {
							let dir = pmesh.polys[src + k] & 0xf;
							switch (dir) {
								case 0: // Portal x-
									if (isMinX)
										mesh.polys[tgt + k] = pmesh.polys[src + k];
									break;
								case 1: // Portal z+
									if (isMaxZ)
										mesh.polys[tgt + k] = pmesh.polys[src + k];
									break;
								case 2: // Portal x+
									if (isMaxX)
										mesh.polys[tgt + k] = pmesh.polys[src + k];
									break;
								case 3: // Portal z-
									if (isMinZ)
										mesh.polys[tgt + k] = pmesh.polys[src + k];
									break;
							}
						}
					}
				}
			}
		}

		// Calculate adjacency.
		buildMeshAdjacency(mesh.polys, mesh.npolys, mesh.nverts, mesh.nvp);
		if (mesh.nverts > 0xffff) {
			throw new RuntimeException("rcBuildPolyMesh: The resulting mesh has too many vertices " + mesh.nverts
				+ " (max " + 0xffff + "). Data can be corrupted.");
		}
		if (mesh.npolys > 0xffff) {
			throw new RuntimeException("rcBuildPolyMesh: The resulting mesh has too many polygons " + mesh.npolys
				+ " (max " + 0xffff + "). Data can be corrupted.");
		}

		ctx.stopTimer("MERGE_POLYMESH");

		return mesh;
	}

	static copyPolyMesh(ctx, src) {
		let dst = new PolyMesh();

		dst.nverts = src.nverts;
		dst.npolys = src.npolys;
		dst.maxpolys = src.npolys;
		dst.nvp = src.nvp;
		RecastVectors.copy(dst.bmin, src.bmin, 0);
		RecastVectors.copy(dst.bmax, src.bmax, 0);
		dst.cs = src.cs;
		dst.ch = src.ch;
		dst.borderSize = src.borderSize;
		dst.maxEdgeError = src.maxEdgeError;

		dst.verts = new Array(src.nverts * 3);
		arraycopy(src.verts, 0, dst.verts, 0, dst.verts.length);
		dst.polys = new Array(src.npolys * 2 * src.nvp);
		arraycopy(src.polys, 0, dst.polys, 0, dst.polys.length);
		dst.regs = new Array(src.npolys);
		arraycopy(src.regs, 0, dst.regs, 0, dst.regs.length);
		dst.areas = new Array(src.npolys);
		arraycopy(src.areas, 0, dst.areas, 0, dst.areas.length);
		dst.flags = new Array(src.npolys);
		arraycopy(src.flags, 0, dst.flags, 0, dst.flags.length);
		return dst;
	}
}

export default RecastMesh;
