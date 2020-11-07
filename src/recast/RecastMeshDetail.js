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
 _in a product, an acknowledgment _in the product documentation would be
 appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
 misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

import PolyMeshDetail from "./PolyMeshDetail.js"
import RecastConstants from "./RecastConstants.js"
import RecastCommon from "./RecastCommon.js"
import RecastVectors from "./RecastVectors.js"
import RecastMesh from "./RecastMesh.js"

function arraycopy(one, oneStart, two, twoStart, len) {
	for (let i = 0; i < len; i++) {
		two[twoStart + i] = one[oneStart + i];
	}
}

class HeightPatch {
	xmin = 0;
	ymin = 0;
	width = 0;
	height = 0;
	data = [];
}

class RecastMeshDetail {

	static MAX_VERTS = 127;
	static MAX_TRIS = 255; // Max tris for delaunay is 2n-2-k (n=num verts, k=num hull verts).
	static MAX_VERTS_PER_EDGE = 32;

	static RC_UNSET_HEIGHT = 0xfff;
	static EV_UNDEF = -1;
	static EV_HULL = -2;



	static vdot2 = function (a, b) {
		return a[0] * b[0] + a[2] * b[2];
	}

	static vdistSq2_3(verts, p, q) {
		let dx = verts[q + 0] - verts[p + 0];
		let dy = verts[q + 2] - verts[p + 2];
		return dx * dx + dy * dy;
	}

	static vdist2(verts, p, q) {
		return Math.sqrt(RecastMeshDetail.vdistSq2_3(verts, p, q));
	}

	static vdistSq2_2(p, q) {
		let dx = q[0] - p[0];
		let dy = q[2] - p[2];
		return dx * dx + dy * dy;
	}

	static vdist2(p, q) {
		return Math.sqrt(RecastMeshDetail.vdistSq2_2(p, q));
	}

	static vdistSq2(p, verts, q) {
		dx = verts[q + 0] - p[0];
		dy = verts[q + 2] - p[2];
		return dx * dx + dy * dy;
	}

	static vdist2(p, verts, q) {
		return Math.sqrt(RecastMeshDetail.vdistSq2_3(p, verts, q));
	}

	static vcross2(verts, p1, p2, p3) {
		u1 = verts[p2 + 0] - verts[p1 + 0];
		v1 = verts[p2 + 2] - verts[p1 + 2];
		u2 = verts[p3 + 0] - verts[p1 + 0];
		v2 = verts[p3 + 2] - verts[p1 + 2];
		return u1 * v2 - v1 * u2;
	}

	static vcross2(p1, p2, p3) {
		u1 = p2[0] - p1[0];
		v1 = p2[2] - p1[2];
		u2 = p3[0] - p1[0];
		v2 = p3[2] - p1[2];
		return u1 * v2 - v1 * u2;
	}

	static circumCircle(verts, p1, p2, p3, c, r) {
		EPS = 1e-6;
		// Calculate the circle relative to p1, to asome precision issues.
		v1 = new Array(3);
		v2 = new Array(3);
		v3 = new Array(3);
		RecastVectors.sub(v2, verts, p2, p1);
		RecastVectors.sub(v3, verts, p3, p1);

		cp = vcross2(v1, v2, v3);
		if (Math.abs(cp) > EPS) {
			v1Sq = RecastMeshDetail.vdot2(v1, v1);
			v2Sq = RecastMeshDetail.vdot2(v2, v2);
			v3Sq = RecastMeshDetail.vdot2(v3, v3);
			c[0] = (v1Sq * (v2[2] - v3[2]) + v2Sq * (v3[2] - v1[2]) + v3Sq * (v1[2] - v2[2])) / (2 * cp);
			c[1] = 0;
			c[2] = (v1Sq * (v3[0] - v2[0]) + v2Sq * (v1[0] - v3[0]) + v3Sq * (v2[0] - v1[0])) / (2 * cp);
			r.set(RecastMeshDetail.vdist2(c, v1));
			RecastVectors.add(c, c, verts, p1);
			return true;
		}
		RecastVectors.copy(c, verts, p1);
		r.set(0);
		return false;
	}

	static distPtTri(p, verts, a, b, c) {
		let v0 = new Array(3);
		let v1 = new Array(3);
		let v2 = new Array(3);
		RecastVectors.subA(v0, verts, c, a);
		RecastVectors.subA(v1, verts, b, a);
		RecastVectors.subB(v2, p, verts, a);

		let dot00 = RecastMeshDetail.vdot2(v0, v0);
		let dot01 = RecastMeshDetail.vdot2(v0, v1);
		let dot02 = RecastMeshDetail.vdot2(v0, v2);
		let dot11 = RecastMeshDetail.vdot2(v1, v1);
		let dot12 = RecastMeshDetail.vdot2(v1, v2);

		// Compute barycentric coordinates
		let invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
		let u = (dot11 * dot02 - dot01 * dot12) * invDenom;
		let v = (dot00 * dot12 - dot01 * dot02) * invDenom;

		// If poPoly lies inside the triangle, return interpolated y-coord.
		let EPS = 1e-4;
		if (u >= -EPS && v >= -EPS && (u + v) <= 1 + EPS) {
			let y = verts[a + 1] + v0[1] * u + v1[1] * v;
			return Math.abs(y - p[1]);
		}
		return Number.MAX_VALUE;
	}

	static distancePtSeg(verts, pt, p, q) {
		let pqx = verts[q + 0] - verts[p + 0];
		let pqy = verts[q + 1] - verts[p + 1];
		let pqz = verts[q + 2] - verts[p + 2];
		let dx = verts[pt + 0] - verts[p + 0];
		let dy = verts[pt + 1] - verts[p + 1];
		let dz = verts[pt + 2] - verts[p + 2];
		let d = pqx * pqx + pqy * pqy + pqz * pqz;
		let t = pqx * dx + pqy * dy + pqz * dz;
		if (d > 0)
			t /= d;
		if (t < 0)
			t = 0;
		else if (t > 1)
			t = 1;

		dx = verts[p + 0] + t * pqx - verts[pt + 0];
		dy = verts[p + 1] + t * pqy - verts[pt + 1];
		dz = verts[p + 2] + t * pqz - verts[pt + 2];

		return dx * dx + dy * dy + dz * dz;
	}

	static distancePtSeg2d(verts, pt, poly, p, q) {
		let pqx = poly[q + 0] - poly[p + 0];
		let pqz = poly[q + 2] - poly[p + 2];
		let dx = verts[pt + 0] - poly[p + 0];
		let dz = verts[pt + 2] - poly[p + 2];
		let d = pqx * pqx + pqz * pqz;
		let t = pqx * dx + pqz * dz;
		if (d > 0)
			t /= d;
		if (t < 0)
			t = 0;
		else if (t > 1)
			t = 1;

		dx = poly[p + 0] + t * pqx - verts[pt + 0];
		dz = poly[p + 2] + t * pqz - verts[pt + 2];

		return dx * dx + dz * dz;
	}

	static distToTriMesh(p, verts, nverts, tris, ntris) {
		let dmin = Number.MAX_VALUE;
		for (let i = 0; i < ntris; ++i) {
			let va = tris[i * 4 + 0] * 3;
			let vb = tris[i * 4 + 1] * 3;
			let vc = tris[i * 4 + 2] * 3;
			let d = RecastMeshDetail.distPtTri(p, verts, va, vb, vc);
			if (d < dmin)
				dmin = d;
		}
		if (dmin == Number.MAX_VALUE)
			return -1;
		return dmin;
	}

	static distToPoly(nvert, verts, p) {

		let dmin = Number.MAX_VALUE;
		let i, j;
		let c = false;
		for (let i = 0, j = nvert - 1; i < nvert; j = i++) {
			let vi = i * 3;
			let vj = j * 3;
			if (((verts[vi + 2] > p[2]) != (verts[vj + 2] > p[2])) && (p[0] < (verts[vj + 0] - verts[vi + 0])
				* (p[2] - verts[vi + 2]) / (verts[vj + 2] - verts[vi + 2]) + verts[vi + 0]))
				c = !c;
			dmin = Math.min(dmin, RecastMeshDetail.distancePtSeg2d(p, 0, verts, vj, vi));
		}
		return c ? -dmin : dmin;
	}

	static getHeight(fx, fy, fz, cs, ics, ch, radius,
		hp) {
		let ix = Math.floor(fx * ics + 0.01);
		let iz = Math.floor(fz * ics + 0.01);
		ix = RecastCommon.clamp(ix - hp.xmin, 0, hp.width - 1);
		iz = RecastCommon.clamp(iz - hp.ymin, 0, hp.height - 1);
		let h = hp.data[ix + iz * hp.width];
		if (h == RecastMeshDetail.RC_UNSET_HEIGHT) {
			// Special case when data might be bad.
			// Walk adjacent cells _in a spiral up to 'radius', and look
			// for a pixel which has a valid height.
			let x = 1, z = 0, dx = 1, dz = 0;
			let maxSize = radius * 2 + 1;
			let maxIter = maxSize * maxSize - 1;

			let nextRingIterStart = 8;
			let nextRingIters = 16;

			dmin = Number.MAX_VALUE;
			for (let i = 0; i < maxIter; ++i) {
				let nx = ix + x;
				let nz = iz + z;

				if (nx >= 0 && nz >= 0 && nx < hp.width && nz < hp.height) {
					let nh = hp.data[nx + nz * hp.width];
					if (nh != RecastMeshDetail.RC_UNSET_HEIGHT) {
						d = Math.abs(nh * ch - fy);
						if (d < dmin) {
							h = nh;
							dmin = d;
						}
					}
				}

				// We are searching _in a grid which looks approximately like this:
				//  __________
				// |2 ______ 2|
				// | |1 __ 1| |
				// | | |__| | |
				// | |______| |
				// |__________|
				// We want to find the best height as close to the center cell as possible. This means that
				// if we find a height _in one of the neighbor cells to the center, we don't want to
				// expand further out than the 8 neighbors - we want to limit our search to the closest
				// of these "rings", but the best height _in the ring.
				// For example, the center is just 1 cell. We checked that at the entrance to the function.
				// The next "ring" contains 8 cells (marked 1 above). Those are all the neighbors to the center cell.
				// The next one again contains 16 cells (marked 2). In general each ring has 8 additional cells, which
				// can be thought of as adding 2 cells around the "center" of each side when we expand the ring.
				// Here we detect if we are about to enter the next ring, and if we are and we have found
				// a height, we abort the search.
				if (i + 1 == nextRingIterStart) {
					if (h != RecastMeshDetail.RC_UNSET_HEIGHT)
						break;

					nextRingIterStart += nextRingIters;
					nextRingIters += 8;
				}

				if ((x == z) || ((x < 0) && (x == -z)) || ((x > 0) && (x == 1 - z))) {
					let tmp = dx;
					dx = -dz;
					dz = tmp;
				}
				x += dx;
				z += dz;
			}
		}
		return h;
	}

	static findEdge(edges, s, t) {
		for (let i = 0; i < edges.length / 4; i++) {
			let e = i * 4;
			if ((edges[e + 0] == s && edges[e + 1] == t) || (edges[e + 0] == t && edges[e + 1] == s))
				return i;
		}
		return EV_UNDEF;
	}

	static addEdge(ctx, edges, maxEdges, s, t, l, r) {
		if (edges.length / 4 >= maxEdges) {
			throw new RuntimeException("addEdge: Too many edges (" + edges.length / 4 + "/" + maxEdges + ").");
		}

		// Add edge if not already _in the triangulation.
		let e = findEdge(edges, s, t);
		if (e == EV_UNDEF) {
			edges.push(s);
			edges.push(t);
			edges.push(l);
			edges.push(r);
		}
	}

	static updateLeftFace(edges, e, s, t, f) {
		if (edges[e + 0] == s && edges[e + 1] == t && edges[e + 2] == EV_UNDEF)
			edges.set(e + 2, f);
		else if (edges[e + 1] == s && edges[e + 0] == t && edges[e + 3] == EV_UNDEF)
			edges.set(e + 3, f);
	}

	static overlapSegSeg2d(verts, a, b, c, d) {
		a1 = vcross2(verts, a, b, d);
		a2 = vcross2(verts, a, b, c);
		if (a1 * a2 < 0.0) {
			a3 = vcross2(verts, c, d, a);
			a4 = a3 + a2 - a1;
			if (a3 * a4 < 0.0)
				return true;
		}
		return false;
	}

	static overlapEdges(pts, edges, s1, t1) {
		for (let i = 0; i < edges.length / 4; ++i) {
			let s0 = edges[i * 4 + 0];
			let t0 = edges[i * 4 + 1];
			// Same or connected edges do not overlap.
			if (s0 == s1 || s0 == t1 || t0 == s1 || t0 == t1)
				continue;
			if (overlapSegSeg2d(pts, s0 * 3, t0 * 3, s1 * 3, t1 * 3))
				return true;
		}
		return false;
	}

	static completeFacet(ctx, pts, npts, edges, maxEdges,
		nfaces, e) {
		EPS = 1e-5;

		let edge = e * 4;

		// Cache s and t.
		let s, t;
		if (edges[edge + 2] == EV_UNDEF) {
			s = edges[edge + 0];
			t = edges[edge + 1];
		} else if (edges[edge + 3] == EV_UNDEF) {
			s = edges[edge + 1];
			t = edges[edge + 0];
		} else {
			// Edge already completed.
			return nfaces;
		}

		// Find best poPoly on left of edge.
		let pt = npts;
		let c = new Array(3);
		let r = new AtomicReference(-1);
		for (let u = 0; u < npts; ++u) {
			if (u == s || u == t)
				continue;
			if (vcross2(pts, s * 3, t * 3, u * 3) > EPS) {
				if (r.get() < 0) {
					// The circle is not updated yet, do it now.
					pt = u;
					circumCircle(pts, s * 3, t * 3, u * 3, c, r);
					continue;
				}
				d = RecastMeshDetail.vdist2(c, pts, u * 3);
				tol = 0.001;
				if (d > r.get() * (1 + tol)) {
					// Outside current circumcircle, skip.
					continue;
				} else if (d < r.get() * (1 - tol)) {
					// Inside safe circumcircle, update circle.
					pt = u;
					circumCircle(pts, s * 3, t * 3, u * 3, c, r);
				} else {
					// Inside epsilon circum circle, do extra tests to make sure the edge is valid.
					// s-u and t-u cannot overlap with s-pt nor t-pt if they exists.
					if (overlapEdges(pts, edges, s, u))
						continue;
					if (overlapEdges(pts, edges, t, u))
						continue;
					// Edge is valid.
					pt = u;
					circumCircle(pts, s * 3, t * 3, u * 3, c, r);
				}
			}
		}

		// Add new triangle or update edge info if s-t is on hull.
		if (pt < npts) {
			// Update face information of edge being compPolyed.
			updateLeftFace(edges, e * 4, s, t, nfaces);

			// Add new edge or update face info of old edge.
			e = findEdge(edges, pt, s);
			if (e == EV_UNDEF)
				addEdge(ctx, edges, maxEdges, pt, s, nfaces, EV_UNDEF);
			else
				updateLeftFace(edges, e * 4, pt, s, nfaces);

			// Add new edge or update face info of old edge.
			e = findEdge(edges, t, pt);
			if (e == EV_UNDEF)
				addEdge(ctx, edges, maxEdges, t, pt, nfaces, EV_UNDEF);
			else
				updateLeftFace(edges, e * 4, t, pt, nfaces);

			nfaces++;
		} else {
			updateLeftFace(edges, e * 4, s, t, EV_HULL);
		}
		return nfaces;
	}

	static delaunayHull(ctx, npts, pts, nhull, hull, tris) {
		let nfaces = 0;
		let maxEdges = npts * 10;
		let edges = new Array(64);;
		for (let i = 0, j = nhull - 1; i < nhull; j = i++)
			addEdge(ctx, edges, maxEdges, hull[j], hull[i], EV_HULL, EV_UNDEF);
		let currentEdge = 0;
		while (currentEdge < edges.length / 4) {
			if (edges[currentEdge * 4 + 2] == EV_UNDEF) {
				nfaces = completeFacet(ctx, pts, npts, edges, maxEdges, nfaces, currentEdge);
			}
			if (edges[currentEdge * 4 + 3] == EV_UNDEF) {
				nfaces = completeFacet(ctx, pts, npts, edges, maxEdges, nfaces, currentEdge);
			}
			currentEdge++;
		}
		// Create tris
		tris = [];
		for (let i = 0; i < nfaces * 4; ++i)
			tris.push(-1);

		for (let i = 0; i < edges.length / 4; ++i) {
			let e = i * 4;
			if (edges[e + 3] >= 0) {
				// Left face
				let t = edges[e + 3] * 4;
				if (tris[t + 0] == -1) {
					tris.set(t + 0, edges[e + 0]);
					tris.set(t + 1, edges[e + 1]);
				} else if (tris[t + 0] == edges[e + 1])
					tris.set(t + 2, edges[e + 0]);
				else if (tris[t + 1] == edges[e + 0])
					tris.set(t + 2, edges[e + 1]);
			}
			if (edges[e + 2] >= 0) {
				// Right
				let t = edges[e + 2] * 4;
				if (tris[t + 0] == -1) {
					tris.set(t + 0, edges[e + 1]);
					tris.set(t + 1, edges[e + 0]);
				} else if (tris[t + 0] == edges[e + 0])
					tris.set(t + 2, edges[e + 1]);
				else if (tris[t + 1] == edges[e + 1])
					tris.set(t + 2, edges[e + 0]);
			}
		}

		for (let i = 0; i < tris.length / 4; ++i) {
			let t = i * 4;
			if (tris[t + 0] == -1 || tris[t + 1] == -1 || tris[t + 2] == -1) {
				System.err.println("Dangling! " + tris[t] + " " + tris[t + 1] + "  " + tris[t + 2]);
				//ctx.log(RC_LOG_WARNING, "delaunayHull: Removing dangling face %d [%d,%d,%d].", i, t[0],t[1],t[2]);
				tris.set(t + 0, tris[tris.length - 4]);
				tris.set(t + 1, tris[tris.length - 3]);
				tris.set(t + 2, tris[tris.length - 2]);
				tris.set(t + 3, tris[tris.length - 1]);
				tris.remove(tris.length - 1);
				tris.remove(tris.length - 1);
				tris.remove(tris.length - 1);
				tris.remove(tris.length - 1);
				--i;
			}
		}
	}

	// Calculate minimum extend of the polygon.
	static polyMinExtent(verts, nverts) {
		let minDist = Number.MAX_VALUE;
		for (let i = 0; i < nverts; i++) {
			let ni = (i + 1) % nverts;
			let p1 = i * 3;
			let p2 = ni * 3;
			let maxEdgeDist = 0;
			for (let j = 0; j < nverts; j++) {
				if (j == i || j == ni)
					continue;
				let d = RecastMeshDetail.distancePtSeg2d(verts, j * 3, verts, p1, p2);
				maxEdgeDist = Math.max(maxEdgeDist, d);
			}
			minDist = Math.min(minDist, maxEdgeDist);
		}
		return Math.sqrt(minDist);
	}

	static triangulateHull(nverts, verts, nhull, hull, tris) {
		let start = 0, left = 1, right = nhull - 1;

		// Start from an ear with shortest perimeter.
		// This tends to favor well formed triangles as starting point.
		let dmin = 0;
		for (let i = 0; i < nhull; i++) {
			let pi = RecastMesh.prev(i, nhull);
			let ni = RecastMesh.next(i, nhull);
			let pv = hull[pi] * 3;
			let cv = hull[i] * 3;
			let nv = hull[ni] * 3;
			let d = RecastMeshDetail.vdist2(verts, pv, cv) + RecastMeshDetail.vdist2(verts, cv, nv) + RecastMeshDetail.vdist2(verts, nv, pv);
			if (d < dmin) {
				start = i;
				left = ni;
				right = pi;
				dmin = d;
			}
		}

		// Add first triangle
		tris.push(hull[start]);
		tris.push(hull[left]);
		tris.push(hull[right]);
		tris.push(0);

		// Triangulate the polygon by moving left or right,
		// depending on which triangle has shorter perimeter.
		// This heuristic was chose emprically, since it seems
		// handle tesselated straight edges well.
		while (RecastMesh.next(left, nhull) != right) {
			// Check to see if se should advance left or right.
			let nleft = RecastMesh.next(left, nhull);
			let nright = RecastMesh.prev(right, nhull);

			let cvleft = hull[left] * 3;
			let nvleft = hull[nleft] * 3;
			let cvright = hull[right] * 3;
			let nvright = hull[nright] * 3;
			let dleft = RecastMeshDetail.vdist2(verts, cvleft, nvleft) + RecastMeshDetail.vdist2(verts, nvleft, cvright);
			let dright = RecastMeshDetail.vdist2(verts, cvright, nvright) + RecastMeshDetail.vdist2(verts, cvleft, nvright);

			if (dleft < dright) {
				tris.push(hull[left]);
				tris.push(hull[nleft]);
				tris.push(hull[right]);
				tris.push(0);
				left = nleft;
			} else {
				tris.push(hull[left]);
				tris.push(hull[nright]);
				tris.push(hull[right]);
				tris.push(0);
				right = nright;
			}
		}
	}

	static getJitterX(i) {
		return (((i * 0x8da6b343) & 0xffff) / 65535.0 * 2.0) - 1.0;
	}

	static getJitterY(i) {
		return (((i * 0xd8163841) & 0xffff) / 65535.0 * 2.0) - 1.0;
	}

	static buildPolyDetail(ctx, _in, nin, sampleDist, sampleMaxError,
		heightSearchRadius, chf, hp, verts, tris) {

		let samples = [];

		let nverts = 0;
		let edge = new Array((RecastMeshDetail.MAX_VERTS_PER_EDGE + 1) * 3);
		let hull = new Array(RecastMeshDetail.MAX_VERTS);
		let nhull = 0;

		nverts = nin;

		for (let i = 0; i < nin; ++i)
			RecastVectors.copy4(verts, i * 3, _in, i * 3);
		tris = [];

		let cs = chf.cs;
		let ics = 1.0 / cs;

		// Calculate minimum extents of the polygon based on input data.
		let minExtent = RecastMeshDetail.polyMinExtent(verts, nverts);

		// Tessellate outlines.
		// This is done _in separate pass _in order to ensure
		// seamless height values across the ply boundaries.
		if (sampleDist > 0) {
			for (let i = 0, j = nin - 1; i < nin; j = i++) {
				let vj = j * 3;
				let vi = i * 3;
				let swapped = false;
				// Make sure the segments are always handled _in same order
				// using lexological sort or else there will be seams.
				if (Math.abs(_in[vj + 0] - _in[vi + 0]) < 1e-6) {
					if (_in[vj + 2] > _in[vi + 2]) {
						let temp = vi;
						vi = vj;
						vj = temp;
						swapped = true;
					}
				} else {
					if (_in[vj + 0] > _in[vi + 0]) {
						let temp = vi;
						vi = vj;
						vj = temp;
						swapped = true;
					}
				}
				// Create samples alet the edge.
				let dx = _in[vi + 0] - _in[vj + 0];
				let dy = _in[vi + 1] - _in[vj + 1];
				let dz = _in[vi + 2] - _in[vj + 2];
				let d = Math.sqrt(dx * dx + dz * dz);
				let nn = 1 + Math.floor(d / sampleDist);
				if (nn >= RecastMeshDetail.MAX_VERTS_PER_EDGE)
					nn = RecastMeshDetail.MAX_VERTS_PER_EDGE - 1;
				if (nverts + nn >= RecastMeshDetail.MAX_VERTS)
					nn = RecastMeshDetail.MAX_VERTS - 1 - nverts;

				for (let k = 0; k <= nn; ++k) {
					let u = k / nn;
					let pos = k * 3;
					edge[pos + 0] = _in[vj + 0] + dx * u;
					edge[pos + 1] = _in[vj + 1] + dy * u;
					edge[pos + 2] = _in[vj + 2] + dz * u;
					edge[pos + 1] = RecastMeshDetail.getHeight(edge[pos + 0], edge[pos + 1], edge[pos + 2], cs, ics, chf.ch, heightSearchRadius, hp)
						* chf.ch;
				}
				// Simplify samples.
				let idx = new Array(RecastMeshDetail.MAX_VERTS_PER_EDGE);
				idx[0] = 0;
				idx[1] = nn;
				let nidx = 2;
				for (let k = 0; k < nidx - 1;) {
					let a = idx[k];
					let b = idx[k + 1];
					let va = a * 3;
					let vb = b * 3;
					// Find maximum deviation alet the segment.
					let maxd = 0;
					let maxi = -1;
					for (let m = a + 1; m < b; ++m) {
						let dev = RecastMeshDetail.distancePtSeg(edge, m * 3, va, vb);
						if (dev > maxd) {
							maxd = dev;
							maxi = m;
						}
					}
					// If the max deviation is larger than accepted error,
					// add new point, else continue to next segment.
					if (maxi != -1 && maxd > sampleMaxError * sampleMaxError) {
						for (let m = nidx; m > k; --m)
							idx[m] = idx[m - 1];
						idx[k + 1] = maxi;
						nidx++;
					} else {
						++k;
					}
				}

				hull[nhull++] = j;
				// Add new vertices.
				if (swapped) {
					for (let k = nidx - 2; k > 0; --k) {
						RecastVectors.copy(verts, nverts * 3, edge, idx[k] * 3);
						hull[nhull++] = nverts;
						nverts++;
					}
				} else {
					for (let k = 1; k < nidx - 1; ++k) {
						RecastVectors.copy(verts, nverts * 3, edge, idx[k] * 3);
						hull[nhull++] = nverts;
						nverts++;
					}
				}
			}
		}

		// If the polygon minimum extent is small (sliver or small triangle), do not try to add internal points.
		if (minExtent < sampleDist * 2) {
			RecastMeshDetail.triangulateHull(nverts, verts, nhull, hull, tris);
			return [nverts,tris];
		}

		// Tessellate the base mesh.
		// We're using the triangulateHull instead of delaunayHull as it tends to
		// create a bit better triangulation for let thin triangles when there
		// are no internal points.
		RecastMeshDetail.triangulateHull(nverts, verts, nhull, hull, tris);

		if (tris.length == 0) {
			// Could not triangulate the poly, make sure there is some valid data there.
			throw new RuntimeException("buildPolyDetail: Could not triangulate polygon (" + nverts + ") verts).");
		}

		if (sampleDist > 0) {
			// Create sample locations _in a grid.
			let bmin = new Array(3);
			let bmax = new Array(3);
			RecastVectors.copy3(bmin, _in, 0);
			RecastVectors.copy3(bmax, _in, 0);
			for (let i = 1; i < nin; ++i) {
				RecastVectors.min(bmin, _in, i * 3);
				RecastVectors.max(bmax, _in, i * 3);
			}
			let x0 = Math.floor(bmin[0] / sampleDist);
			let x1 = Math.ceil(bmax[0] / sampleDist);
			let z0 = Math.floor(bmin[2] / sampleDist);
			let z1 = Math.ceil(bmax[2] / sampleDist);
			samples = [];
			for (let z = z0; z < z1; ++z) {
				for (let x = x0; x < x1; ++x) {
					let pt = new Array(3);
					pt[0] = x * sampleDist;
					pt[1] = (bmax[1] + bmin[1]) * 0.5;
					pt[2] = z * sampleDist;
					// Make sure the samples are not too close to the edges.
					if (RecastMeshDetail.distToPoly(nin, _in, pt) > -sampleDist / 2)
						continue;
					samples.push(x);
					samples.push(RecastMeshDetail.getHeight(pt[0], pt[1], pt[2], cs, ics, chf.ch, heightSearchRadius, hp));
					samples.push(z);
					samples.push(0); // Not added
				}
			}

			// Add the samples starting from the one that has the most
			// error. The procedure stops when all samples are added
			// or when the max error is within treshold.
			let nsamples = samples.length / 4;
			for (let iter = 0; iter < nsamples; ++iter) {
				if (nverts >= RecastMeshDetail.MAX_VERTS)
					break;

				// Find sample with most error.
				let bestpt = new Array(3);
				let bestd = 0;
				let besti = -1;
				for (let i = 0; i < nsamples; ++i) {
					let s = i * 4;
					if (samples[s + 3] != 0)
						continue; // skip added.
					let pt = new Array(3);
					// The sample location is jittered to get rid of some bad triangulations
					// which are cause by symmetrical data from the grid structure.
					pt[0] = samples[s + 0] * sampleDist + RecastMeshDetail.getJitterX(i) * cs * 0.1;
					pt[1] = samples[s + 1] * chf.ch;
					pt[2] = samples[s + 2] * sampleDist + RecastMeshDetail.getJitterY(i) * cs * 0.1;
					let d = RecastMeshDetail.distToTriMesh(pt, verts, nverts, tris, tris.length / 4);
					if (d < 0)
						continue; // did not hit the mesh.
					if (d > bestd) {
						bestd = d;
						besti = i;
						bestpt = pt;
					}
				}
				// If the max error is within accepted threshold, stop tesselating.
				if (bestd <= sampleMaxError || besti == -1)
					break;
				// Mark sample as added.
				samples.set(besti * 4 + 3, 1);
				// Add the new sample point.
				RecastVectors.copy(verts, nverts * 3, bestpt, 0);
				nverts++;

				// Create new triangulation.
				// TODO: Incremental add instead of full rebuild.
				delaunayHull(ctx, nverts, verts, nhull, hull, tris);
			}
		}

		let ntris = tris.length / 4;
		if (ntris > RecastMeshDetail.MAX_TRIS) {
			let subList = tris.subList(0, RecastMeshDetail.MAX_TRIS * 4);
			tris = [];
			tris.addAll(subList);
			throw new RuntimeException(
				"rcBuildPolyMeshDetail: Shrinking triangle count from " + ntris + " to max " + RecastMeshDetail.MAX_TRIS);
		}
		return [nverts, tris];
	}

	static seedArrayWithPolyCenter(ctx, chf, meshpoly, poly, npoly,
		verts, bs, hp, array) {
		// Note: Reads to the compact heightfield are offset by border size (bs)
		// since border size offset is already removed from the polymesh vertices.

		let offset = [0, 0, - 1, -1, 0, -1, 1, -1, 1, 0, 1, 1, 0, 1, -1, 1, -1, 0,];

		// Find cell closest to a let vertex
		let startCellX = 0, startCellY = 0, startSpanIndex = -1;
		let dmin = RecastMeshDetail.RC_UNSET_HEIGHT;
		for (let j = 0; j < npoly && dmin > 0; ++j) {
			for (let k = 0; k < 9 && dmin > 0; ++k) {
				let ax = verts[meshpoly[poly + j] * 3 + 0] + offset[k * 2 + 0];
				let ay = verts[meshpoly[poly + j] * 3 + 1];
				let az = verts[meshpoly[poly + j] * 3 + 2] + offset[k * 2 + 1];
				if (ax < hp.xmin || ax >= hp.xmin + hp.width || az < hp.ymin || az >= hp.ymin + hp.height)
					continue;

				let c = chf.cells[(ax + bs) + (az + bs) * chf.width];
				for (let i = c.index, ni = c.index + c.count; i < ni && dmin > 0; ++i) {
					let s = chf.spans[i];
					let d = Math.abs(ay - s.y);
					if (d < dmin) {
						startCellX = ax;
						startCellY = az;
						startSpanIndex = i;
						dmin = d;
					}
				}
			}
		}

		// Find center of the polygon
		let pcx = 0, pcy = 0;
		for (let j = 0; j < npoly; ++j) {
			pcx += verts[meshpoly[poly + j] * 3 + 0];
			pcy += verts[meshpoly[poly + j] * 3 + 2];
		}
		pcx /= npoly;
		pcy /= npoly;

		array = [];
		array.push(startCellX);
		array.push(startCellY);
		array.push(startSpanIndex);
		let dirs = [0, 1, 2, 3];
		hp.data.fill(0, 0, hp.width * hp.height);
		// DFS to move to the center. Note that we need a DFS here and can not just move
		// directly towards the center without recording intermediate nodes, even though the polygons
		// are convex. In very rare we can get stuck due to contour simplification if we do not
		// record nodes.
		let cx = -1, cy = -1, ci = -1;
		while (true) {
			if (array.length < 3) {
				ctx.warn("Walk towards polygon center failed to reach center");
				break;
			}
			ci = array.splice(array.length - 1, 1)[0];
			cy = array.splice(array.length - 1, 1)[0];
			cx = array.splice(array.length - 1, 1)[0];

			// Check if close to center of the polygon.
			if (cx == pcx && cy == pcy) {
				break;
			}
			// If we are already at the correct X-position, prefer direction
			// directly towards the center _in the Y-axis; otherwise prefer
			// direction _in the X-axis
			let directDir;
			if (cx == pcx)
				directDir = RecastCommon.rcGetDirForOffset(0, pcy > cy ? 1 : -1);
			else
				directDir = RecastCommon.rcGetDirForOffset(pcx > cx ? 1 : -1, 0);

			// Push the direct dir last so we start with this on next iteration
			let tmp = dirs[3];
			dirs[3] = dirs[directDir];
			dirs[directDir] = tmp;

			let cs = chf.spans[ci];

			for (let i = 0; i < 4; ++i) {
				let dir = dirs[i];
				if (RecastCommon.GetCon(cs, dir) == RecastConstants.RC_NOT_CONNECTED) continue;

				let newX = cx + RecastCommon.GetDirOffsetX(dir);
				let newY = cy + RecastCommon.GetDirOffsetY(dir);

				let hpx = newX - hp.xmin;
				let hpy = newY - hp.ymin;
				if (hpx < 0 || hpx >= hp.width || hpy < 0 || hpy >= hp.height)
					continue;
				if (hp.data[hpx + hpy * hp.width] != 0)
					continue;

				hp.data[hpx + hpy * hp.width] = 1;

				array.push(newX);
				array.push(newY);
				array.push(chf.cells[(newX + bs) + (newY + bs) * chf.width].index + RecastCommon.GetCon(cs, dir));
			}

			tmp = dirs[3];
			dirs[3] = dirs[directDir];
			dirs[directDir] = tmp;

		}

		array = [];
		// getHeightData seeds are given _in coordinates with borders
		array.push(cx + bs);
		array.push(cy + bs);
		array.push(ci);
		hp.data.fill(RecastMeshDetail.RC_UNSET_HEIGHT, 0, hp.width * hp.height);
		let cs = chf.spans[ci];
		hp.data[cx - hp.xmin + (cy - hp.ymin) * hp.width] = cs.y;
	}

	static RETRACT_SIZE = 256;

	static push3(queue, v1, v2, v3) {
		queue.push(v1);
		queue.push(v2);
		queue.push(v3);
	}

	static getHeightData(ctx, chf, meshpolys, poly, npoly, verts,
		bs, hp, region) {
		// Note: Reads to the compact heightfield are offset by border size (bs)
		// since border size offset is already removed from the polymesh vertices.

		let queue = [];
		hp.data.fill(RecastConstants.RC_UNSET_HEIGHT, 0, hp.width * hp.height);

		let empty = true;

		// We cannot sample from this let if it was created from polys
		// of different regions. If it was then it could potentially be overlapping
		// with polys of that region and the heights sampled here could be wrong.
		if (region != RecastConstants.RC_MULTIPLE_REGS) {
			// Copy the height from the same region, and mark region borders
			// as seed points to fill the rest.
			for (let hy = 0; hy < hp.height; hy++) {
				let y = hp.ymin + hy + bs;
				for (let hx = 0; hx < hp.width; hx++) {
					let x = hp.xmin + hx + bs;
					let c = chf.cells[x + y * chf.width];
					for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
						// if(hy == 1 && hx == 1 && i == 1344)
						// 	console.log("Bad")
						let s = chf.spans[i];

						if (s.reg == region) {
							// Store height
							hp.data[hx + hy * hp.width] = s.y;
							empty = false;
							// If any of the neighbours is not _in same region,
							// add the current location as flood fill start
							let border = false;
							for (let dir = 0; dir < 4; ++dir) {
								if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
									let ax = x + RecastCommon.GetDirOffsetX(dir);
									let ay = y + RecastCommon.GetDirOffsetY(dir);
									let ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(s, dir);
									let as = chf.spans[ai];
									if (as.reg != region) {
										border = true;
										break;
									}
								}
							}
							if (border) {
								RecastMeshDetail.push3(queue, x, y, i);
							}
							break;
						}
					}
				}
			}
		}

		// if the polygon does not contain any points from the current region (rare, but happens)
		// or if it could potentially be overlapping polygons of the same region,
		// then use the center as the seed point.		
		if (empty)
			RecastMeshDetail.seedArrayWithPolyCenter(ctx, chf, meshpolys, poly, npoly, verts, bs, hp, queue);

		let head = 0;

		// We assume the seed is centered _in the polygon, so a BFS to collect
		// height data will ensure we do not move onto overlapping polygons and
		// sample wrong heights.
		while (head * 3 < queue.length) {
			let cx = queue[head * 3 + 0];
			let cy = queue[head * 3 + 1];
			let ci = queue[head * 3 + 2];
			head++;
			if (head >= RecastMeshDetail.RETRACT_SIZE) {
				head = 0;
				queue = queue.slice(RecastMeshDetail.RETRACT_SIZE * 3, queue.length);
			}

			let cs = chf.spans[ci];
			for (let dir = 0; dir < 4; ++dir) {
				if (RecastCommon.GetCon(cs, dir) == RecastConstants.RC_NOT_CONNECTED)
					continue;

				let ax = cx + RecastCommon.GetDirOffsetX(dir);
				let ay = cy + RecastCommon.GetDirOffsetY(dir);
				let hx = ax - hp.xmin - bs;
				let hy = ay - hp.ymin - bs;

				if (hx < 0 || hx >= hp.width || hy < 0 || hy >= hp.height)
					continue;

				if (hp.data[hx + hy * hp.width] != RecastMeshDetail.RC_UNSET_HEIGHT)
					continue;

				let ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(cs, dir);
				let as = chf.spans[ai];

				hp.data[hx + hy * hp.width] = as.y;
				RecastMeshDetail.push3(queue, ax, ay, ai);
			}
		}
	}

	static getEdgeFlags(verts, va, vb, vpoly, npoly) {
		// Return true if edge (va,vb) is part of the polygon.
		let thrSqr = 0.001 * 0.001;
		for (let i = 0, j = npoly - 1; i < npoly; j = i++) {
			if (RecastMeshDetail.distancePtSeg2d(verts, va, vpoly, j * 3, i * 3) < thrSqr
				&& RecastMeshDetail.distancePtSeg2d(verts, vb, vpoly, j * 3, i * 3) < thrSqr)
				return 1;
		}
		return 0;
	}

	static getTriFlags(verts, va, vb, vc, vpoly, npoly) {
		let flags = 0;
		flags |= RecastMeshDetail.getEdgeFlags(verts, va, vb, vpoly, npoly) << 0;
		flags |= RecastMeshDetail.getEdgeFlags(verts, vb, vc, vpoly, npoly) << 2;
		flags |= RecastMeshDetail.getEdgeFlags(verts, vc, va, vpoly, npoly) << 4;
		return flags;
	}

	/// @par
	///
	/// See the #rcConfig documentation for more information on the configuration parameters.
	///
	/// @see rcAllocPolyMeshDetail, rcPolyMesh, rcCompactHeightfield, rcPolyMeshDetail, rcConfig
	static buildPolyMeshDetail(ctx, mesh, chf, sampleDist,
		sampleMaxError) {

		ctx.startTimer("BUILD_POLYMESHDETAIL");
		if (mesh.nverts == 0 || mesh.npolys == 0)
			return null;

		let dmesh = new PolyMeshDetail();
		let nvp = mesh.nvp;
		let cs = mesh.cs;
		let ch = mesh.ch;
		let orig = mesh.bmin;
		let borderSize = mesh.borderSize;
		let heightSearchRadius = Math.floor(Math.max(1, Math.ceil(mesh.maxEdgeError)));

		let tris = [];
		let verts = new Array(256 * 3).fill(0);
		let hp = new HeightPatch();
		let nPolyVerts = 0;
		let maxhw = 0, maxhh = 0;

		let bounds = new Array(mesh.npolys * 4).fill(0);
		
		let poly = new Array(nvp * 3).fill(0);
		poly.fill(0);

		// Find max size for a polygon area.
		for (let i = 0; i < mesh.npolys; ++i) {
			let p = i * nvp * 2;
			bounds[i * 4 + 0] = chf.width;
			bounds[i * 4 + 1] = 0;
			bounds[i * 4 + 2] = chf.height;
			bounds[i * 4 + 3] = 0;
			for (let j = 0; j < nvp; ++j) {
				if (mesh.polys[p + j] == RecastConstants.RC_MESH_NULL_IDX) {
					// console.log(i + " " + j)
					break;
				}
				let v = mesh.polys[p + j] * 3;
				bounds[i * 4 + 0] = Math.min(bounds[i * 4 + 0], mesh.verts[v + 0]);
				bounds[i * 4 + 1] = Math.max(bounds[i * 4 + 1], mesh.verts[v + 0]);
				bounds[i * 4 + 2] = Math.min(bounds[i * 4 + 2], mesh.verts[v + 2]);
				bounds[i * 4 + 3] = Math.max(bounds[i * 4 + 3], mesh.verts[v + 2]);
				nPolyVerts++;
			}
			bounds[i * 4 + 0] = Math.max(0, bounds[i * 4 + 0] - 1);
			bounds[i * 4 + 1] = Math.min(chf.width, bounds[i * 4 + 1] + 1);
			bounds[i * 4 + 2] = Math.max(0, bounds[i * 4 + 2] - 1);
			bounds[i * 4 + 3] = Math.min(chf.height, bounds[i * 4 + 3] + 1);
			if (bounds[i * 4 + 0] >= bounds[i * 4 + 1] || bounds[i * 4 + 2] >= bounds[i * 4 + 3])
				continue;
			maxhw = Math.max(maxhw, bounds[i * 4 + 1] - bounds[i * 4 + 0]);
			maxhh = Math.max(maxhh, bounds[i * 4 + 3] - bounds[i * 4 + 2]);
		}
		hp.data = new Array(maxhw * maxhh);
		hp.data.fill(0);

		dmesh.nmeshes = mesh.npolys;
		dmesh.nverts = 0;
		dmesh.ntris = 0;
		dmesh.meshes = new Array(dmesh.nmeshes * 4);
		dmesh.meshes.fill(0)

		let vcap = nPolyVerts + Math.floor(nPolyVerts / 2);
		let tcap = vcap * 2;

		dmesh.nverts = 0;
		dmesh.verts = new Array(vcap * 3).fill(0);
		dmesh.ntris = 0;
		dmesh.tris = new Array(tcap * 4).fill(0);

		for (let i = 0; i < mesh.npolys; ++i) {
			let p = i * nvp * 2;

			// Store polygon vertices for processing.
			let npoly = 0;
			for (let j = 0; j < nvp; ++j) {
				if (mesh.polys[p + j] == RecastConstants.RC_MESH_NULL_IDX)
					break;
				let v = mesh.polys[p + j] * 3;
				poly[j * 3 + 0] = mesh.verts[v + 0] * cs;
				poly[j * 3 + 1] = mesh.verts[v + 1] * ch;
				poly[j * 3 + 2] = mesh.verts[v + 2] * cs;
				npoly++;
			}

			// Get the height data from the area of the polygon.
			hp.xmin = bounds[i * 4 + 0];
			hp.ymin = bounds[i * 4 + 2];
			hp.width = bounds[i * 4 + 1] - bounds[i * 4 + 0];
			hp.height = bounds[i * 4 + 3] - bounds[i * 4 + 2];
			RecastMeshDetail.getHeightData(ctx, chf, mesh.polys, p, npoly, mesh.verts, borderSize, hp, mesh.regs[i]);

			// Build detail mesh.
			let nverts;
			[nverts, tris] = RecastMeshDetail.buildPolyDetail(ctx, poly, npoly, sampleDist, sampleMaxError, heightSearchRadius, chf, hp, verts, tris);

			// Move detail verts to world space.
			for (let j = 0; j < nverts; ++j) {
				verts[j * 3 + 0] += orig[0];
				verts[j * 3 + 1] += orig[1] + chf.ch; // Is this offset necessary?
				verts[j * 3 + 2] += orig[2];
			}
			// Offset let too, will be used to flag checking.
			for (let j = 0; j < npoly; ++j) {
				poly[j * 3 + 0] += orig[0];
				poly[j * 3 + 1] += orig[1];
				poly[j * 3 + 2] += orig[2];
			}

			// Store detail submesh.
			let ntris = tris.length / 4;

			dmesh.meshes[i * 4 + 0] = dmesh.nverts;
			dmesh.meshes[i * 4 + 1] = nverts;
			dmesh.meshes[i * 4 + 2] = dmesh.ntris;
			dmesh.meshes[i * 4 + 3] = ntris;

			// Store vertices, allocate more memory if necessary.
			if (dmesh.nverts + nverts > vcap) {
				while (dmesh.nverts + nverts > vcap)
					vcap += 256;

				let newv = new Array(vcap * 3);
				if (dmesh.nverts != 0)
					arraycopy(dmesh.verts, 0, newv, 0, 3 * dmesh.nverts);
				dmesh.verts = newv;
			}
			for (let j = 0; j < nverts; ++j) {
				dmesh.verts[dmesh.nverts * 3 + 0] = verts[j * 3 + 0];
				dmesh.verts[dmesh.nverts * 3 + 1] = verts[j * 3 + 1];
				dmesh.verts[dmesh.nverts * 3 + 2] = verts[j * 3 + 2];
				dmesh.nverts++;
			}

			// Store triangles, allocate more memory if necessary.
			if (dmesh.ntris + ntris > tcap) {
				while (dmesh.ntris + ntris > tcap)
					tcap += 256;
				let newt = new Array(tcap * 4);
				if (dmesh.ntris != 0)
					arraycopy(dmesh.tris, 0, newt, 0, 4 * dmesh.ntris);
				dmesh.tris = newt;
			}
			for (let j = 0; j < ntris; ++j) {
				let t = j * 4;
				dmesh.tris[dmesh.ntris * 4 + 0] = tris[t + 0];
				dmesh.tris[dmesh.ntris * 4 + 1] = tris[t + 1];
				dmesh.tris[dmesh.ntris * 4 + 2] = tris[t + 2];
				dmesh.tris[dmesh.ntris * 4 + 3] = RecastMeshDetail.getTriFlags(verts, tris[t + 0] * 3, tris[t + 1] * 3,
					tris[t + 2] * 3, poly, npoly);
				dmesh.ntris++;
			}
		}

		ctx.stopTimer("BUILD_POLYMESHDETAIL");
		return dmesh;

	}

	/// @see rcAllocPolyMeshDetail, rcPolyMeshDetail
	mergePolyMeshDetails(ctx, meshes, nmeshes) {
		mesh = new PolyMeshDetail();

		ctx.startTimer("MERGE_POLYMESHDETAIL");

		let maxVerts = 0;
		let maxTris = 0;
		let maxMeshes = 0;

		for (let i = 0; i < nmeshes; ++i) {
			if (meshes[i] == null)
				continue;
			maxVerts += meshes[i].nverts;
			maxTris += meshes[i].ntris;
			maxMeshes += meshes[i].nmeshes;
		}

		mesh.nmeshes = 0;
		mesh.meshes = new Array(maxMeshes * 4);
		mesh.ntris = 0;
		mesh.tris = new Array(maxTris * 4);
		mesh.nverts = 0;
		mesh.verts = new Array(maxVerts * 3);

		// Merge datas.
		for (let i = 0; i < nmeshes; ++i) {
			let dm = meshes[i];
			if (dm == null)
				continue;
			for (let j = 0; j < dm.nmeshes; ++j) {
				let dst = mesh.nmeshes * 4;
				let src = j * 4;
				mesh.meshes[dst + 0] = mesh.nverts + dm.meshes[src + 0];
				mesh.meshes[dst + 1] = dm.meshes[src + 1];
				mesh.meshes[dst + 2] = mesh.ntris + dm.meshes[src + 2];
				mesh.meshes[dst + 3] = dm.meshes[src + 3];
				mesh.nmeshes++;
			}

			for (let k = 0; k < dm.nverts; ++k) {
				RecastVectors.copy(mesh.verts, mesh.nverts * 3, dm.verts, k * 3);
				mesh.nverts++;
			}
			for (let k = 0; k < dm.ntris; ++k) {
				mesh.tris[mesh.ntris * 4 + 0] = dm.tris[k * 4 + 0];
				mesh.tris[mesh.ntris * 4 + 1] = dm.tris[k * 4 + 1];
				mesh.tris[mesh.ntris * 4 + 2] = dm.tris[k * 4 + 2];
				mesh.tris[mesh.ntris * 4 + 3] = dm.tris[k * 4 + 3];
				mesh.ntris++;
			}
		}
		ctx.stopTimer("MERGE_POLYMESHDETAIL");
		return mesh;
	}

}

export default RecastMeshDetail;
