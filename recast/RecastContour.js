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

import RecastVectors from "./RecastVectors.js"
import RecastConstants from "./RecastConstants.js"
import RecastCommon from "./RecastCommon.js"
import ContourSet from "./ContourSet.js"
import RecastMesh from "./RecastMesh.js"
import Contour from "./Contour.js"

class RecastContour {

	static ContourRegion = class ContourRegion {
		outline;
		holes = [];
		nholes;
	}

	static ContourHole = class ContourHole {
		leftmost;
		minx;
		minz;
		contour;
	}

	static PotentialDiagonal = class PotentialDiagonal {
		dist;
		vert;
	}

	static getCornerHeight(x, y, i, dir, chf, isBorderVertex) {
		let s = chf.spans[i];
		let ch = s.y;
		let dirp = (dir + 1) & 0x3;

		let regs = [0, 0, 0, 0];

		// Combine region and area codes in order to prevent
		// border vertices which are in between two areas to be removed.
		regs[0] = chf.spans[i].reg | (chf.areas[i] << 16);

		if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
			let ax = x + RecastCommon.GetDirOffsetX(dir);
			let ay = y + RecastCommon.GetDirOffsetY(dir);
			let ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(s, dir);
			let as = chf.spans[ai];
			ch = Math.max(ch, as.y);
			regs[1] = chf.spans[ai].reg | (chf.areas[ai] << 16);
			if (RecastCommon.GetCon(as, dirp) != RecastConstants.RC_NOT_CONNECTED) {
				let ax2 = ax + RecastCommon.GetDirOffsetX(dirp);
				let ay2 = ay + RecastCommon.GetDirOffsetY(dirp);
				let ai2 = chf.cells[ax2 + ay2 * chf.width].index + RecastCommon.GetCon(as, dirp);
				let as2 = chf.spans[ai2];
				ch = Math.max(ch, as2.y);
				regs[2] = chf.spans[ai2].reg | (chf.areas[ai2] << 16);
			}
		}
		if (RecastCommon.GetCon(s, dirp) != RecastConstants.RC_NOT_CONNECTED) {
			let ax = x + RecastCommon.GetDirOffsetX(dirp);
			let ay = y + RecastCommon.GetDirOffsetY(dirp);
			let ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(s, dirp);
			let as = chf.spans[ai];
			ch = Math.max(ch, as.y);
			regs[3] = chf.spans[ai].reg | (chf.areas[ai] << 16);
			if (RecastCommon.GetCon(as, dir) != RecastConstants.RC_NOT_CONNECTED) {
				let ax2 = ax + RecastCommon.GetDirOffsetX(dir);
				let ay2 = ay + RecastCommon.GetDirOffsetY(dir);
				let ai2 = chf.cells[ax2 + ay2 * chf.width].index + RecastCommon.GetCon(as, dir);
				let as2 = chf.spans[ai2];
				ch = Math.max(ch, as2.y);
				regs[2] = chf.spans[ai2].reg | (chf.areas[ai2] << 16);
			}
		}

		// Check if the vertex is special edge vertex, these vertices will be removed later.
		for (let j = 0; j < 4; ++j) {
			let a = j;
			let b = (j + 1) & 0x3;
			let c = (j + 2) & 0x3;
			let d = (j + 3) & 0x3;

			// The vertex is a border vertex there are two same exterior cells in a row,
			// followed by two interior cells and none of the regions are out of bounds.
			let twoSameExts = (regs[a] & regs[b] & RecastConstants.RC_BORDER_REG) != 0 && regs[a] == regs[b];
			let twoInts = ((regs[c] | regs[d]) & RecastConstants.RC_BORDER_REG) == 0;
			let intsSameArea = (regs[c] >> 16) == (regs[d] >> 16);
			let noZeros = regs[a] != 0 && regs[b] != 0 && regs[c] != 0 && regs[d] != 0;
			if (twoSameExts && twoInts && intsSameArea && noZeros) {
				isBorderVertex = true;
				break;
			}
		}

		return ch;
	}

	static walkContour(x, y, i, chf, flags, points) {
		// Choose the first non-connected edge
		let dir = 0;
		while ((flags[i] & (1 << dir)) == 0)
			dir++;

		let startDir = dir;
		let starti = i;

		let area = chf.areas[i];

		let iter = 0;
		while (++iter < 40000) {
			if ((flags[i] & (1 << dir)) != 0) {
				// Choose the edge corner
				let isBorderVertex = false;
				let isAreaBorder = false;
				let px = x;
				let py = RecastContour.getCornerHeight(x, y, i, dir, chf, isBorderVertex);
				let pz = y;
				switch (dir) {
					case 0:
						pz++;
						break;
					case 1:
						px++;
						pz++;
						break;
					case 2:
						px++;
						break;
				}
				let r = 0;
				let s = chf.spans[i];
				if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
					let ax = x + RecastCommon.GetDirOffsetX(dir);
					let ay = y + RecastCommon.GetDirOffsetY(dir);
					let ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(s, dir);
					r = chf.spans[ai].reg;
					if (area != chf.areas[ai])
						isAreaBorder = true;
				}
				if (isBorderVertex)
					r |= RecastConstants.RC_BORDER_VERTEX;
				if (isAreaBorder)
					r |= RecastConstants.RC_AREA_BORDER;
				points.push(px);
				points.push(py);
				points.push(pz);
				points.push(r);

				flags[i] &= ~(1 << dir); // Remove visited edges
				dir = (dir + 1) & 0x3; // Rotate CW
			} else {
				let ni = -1;
				let nx = x + RecastCommon.GetDirOffsetX(dir);
				let ny = y + RecastCommon.GetDirOffsetY(dir);
				let s = chf.spans[i];
				if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
					let nc = chf.cells[nx + ny * chf.width];
					ni = nc.index + RecastCommon.GetCon(s, dir);
				}
				if (ni == -1) {
					// Should not happen.
					return;
				}
				x = nx;
				y = ny;
				i = ni;
				dir = (dir + 3) & 0x3; // Rotate CCW
			}

			if (starti == i && startDir == dir) {
				break;
			}
		}
	}

	static distancePtSeg(x, z, px, pz, qx, qz) {
		let pqx = qx - px;
		let pqz = qz - pz;
		let dx = x - px;
		let dz = z - pz;
		let d = pqx * pqx + pqz * pqz;
		let t = pqx * dx + pqz * dz;
		if (d > 0)
			t /= d;
		if (t < 0)
			t = 0;
		else if (t > 1)
			t = 1;

		dx = px + t * pqx - x;
		dz = pz + t * pqz - z;

		return dx * dx + dz * dz;
	}

	static simplifyContour(points, simplified, maxError, maxEdgeLen,
		buildFlags) {
		// Add initial points.
		let hasConnections = false;
		for (let i = 0; i < points.length; i += 4) {
			if ((points[i + 3] & RecastConstants.RC_CONTOUR_REG_MASK) != 0) {
				hasConnections = true;
				break;
			}
		}

		if (hasConnections) {
			// The contour has some portals to other regions.
			// Add a new poPoly to every location where the region changes.
			for (let i = 0, ni = points.length / 4; i < ni; ++i) {
				let ii = (i + 1) % ni;
				let differentRegs = (points[i * 4 + 3]
					& RecastConstants.RC_CONTOUR_REG_MASK) != (points[ii * 4 + 3]
						& RecastConstants.RC_CONTOUR_REG_MASK);
				let areaBorders = (points[i * 4 + 3]
					& RecastConstants.RC_AREA_BORDER) != (points[ii * 4 + 3] & RecastConstants.RC_AREA_BORDER);
				if (differentRegs || areaBorders) {
					simplified.push(points[i * 4 + 0]);
					simplified.push(points[i * 4 + 1]);
					simplified.push(points[i * 4 + 2]);
					simplified.push(i);
				}
			}
		}

		if (simplified.length == 0) {
			// If there is no connections at all,
			// create some initial points for the simplification process.
			// Find lower-left and upper-right vertices of the contour.
			let llx = points[0];
			let lly = points[1];
			let llz = points[2];
			let lli = 0;
			let urx = points[0];
			let ury = points[1];
			let urz = points[2];
			let uri = 0;
			for (let i = 0; i < points.length; i += 4) {
				let x = points[i + 0];
				let y = points[i + 1];
				let z = points[i + 2];
				if (x < llx || (x == llx && z < llz)) {
					llx = x;
					lly = y;
					llz = z;
					lli = i / 4;
				}
				if (x > urx || (x == urx && z > urz)) {
					urx = x;
					ury = y;
					urz = z;
					uri = i / 4;
				}
			}
			simplified.push(llx);
			simplified.push(lly);
			simplified.push(llz);
			simplified.push(lli);

			simplified.push(urx);
			simplified.push(ury);
			simplified.push(urz);
			simplified.push(uri);
		}
		// Add points until all raw points are within
		// error tolerance to the simplified shape.
		let pn = points.length / 4;
		for (let i = 0; i < simplified.length / 4;) {
			// console.log(simplified.length)
			let ii = (i + 1) % (simplified.length / 4);

			let ax = simplified[i * 4 + 0];
			let az = simplified[i * 4 + 2];
			let ai = simplified[i * 4 + 3];

			let bx = simplified[ii * 4 + 0];
			let bz = simplified[ii * 4 + 2];
			let bi = simplified[ii * 4 + 3];

			// Find maximum deviation from the segment.
			let maxd = 0;
			let maxi = -1;
			let ci, cinc, endi;

			// Traverse the segment in lexilogical order so that the
			// max deviation is calculated similarly when traversing
			// opposite segments.
			if (bx > ax || (bx == ax && bz > az)) {
				cinc = 1;
				ci = (ai + cinc) % pn;
				endi = bi;
			} else {
				cinc = pn - 1;
				ci = (bi + cinc) % pn;
				endi = ai;
				let temp = ax;
				ax = bx;
				bx = temp;
				temp = az;
				az = bz;
				bz = temp;
			}
			// Tessellate only outer edges or edges between areas.
			// console.log("start")
			if ((points[ci * 4 + 3] & RecastConstants.RC_CONTOUR_REG_MASK) == 0
				|| (points[ci * 4 + 3] & RecastConstants.RC_AREA_BORDER) != 0) {
				while (ci != endi) {
					// if(Math.random() < .01) console.log(`${ci} ${endi}`);
					let d = RecastContour.distancePtSeg(points[ci * 4 + 0], points[ci * 4 + 2], ax, az, bx, bz);
					if (d > maxd) {
						maxd = d;
						maxi = ci;
					}
					ci = (ci + cinc) % pn;
				}
			}
			// console.log("stop")
			// If the max deviation is larger than accepted error,
			// add new point, else continue to next segment.
			if (maxi != -1 && maxd > (maxError * maxError)) {
				// Add the point.
				simplified.splice((i + 1) * 4 + 0, 0, points[maxi * 4 + 0]);
				simplified.splice((i + 1) * 4 + 1, 0, points[maxi * 4 + 1]);
				simplified.splice((i + 1) * 4 + 2, 0, points[maxi * 4 + 2]);
				simplified.splice((i + 1) * 4 + 3, 0, maxi);
			} else {
				++i;
			}
		}
		// Split too let edges.
		if (maxEdgeLen > 0 && (buildFlags
			& (RecastConstants.RC_CONTOUR_TESS_WALL_EDGES | RecastConstants.RC_CONTOUR_TESS_AREA_EDGES)) != 0) {
			for (let i = 0; i < simplified.length / 4;) {
				let ii = (i + 1) % (simplified.length / 4);

				let ax = simplified[i * 4 + 0];
				let az = simplified[i * 4 + 2];
				let ai = simplified[i * 4 + 3];

				let bx = simplified[ii * 4 + 0];
				let bz = simplified[ii * 4 + 2];
				let bi = simplified[ii * 4 + 3];

				// Find maximum deviation from the segment.
				let maxi = -1;
				let ci = (ai + 1) % pn;

				// Tessellate only outer edges or edges between areas.
				let tess = false;
				// Wall edges.
				if ((buildFlags & RecastConstants.RC_CONTOUR_TESS_WALL_EDGES) != 0
					&& (points[ci * 4 + 3] & RecastConstants.RC_CONTOUR_REG_MASK) == 0)
					tess = true;
				// Edges between areas.
				if ((buildFlags & RecastConstants.RC_CONTOUR_TESS_AREA_EDGES) != 0
					&& (points[ci * 4 + 3] & RecastConstants.RC_AREA_BORDER) != 0)
					tess = true;

				if (tess) {
					let dx = bx - ax;
					let dz = bz - az;
					if (dx * dx + dz * dz > maxEdgeLen * maxEdgeLen) {
						// Round based on the segments in lexilogical order so that the
						// max tesselation is consistent regardles in which direction
						// segments are traversed.
						let n = bi < ai ? (bi + pn - ai) : (bi - ai);
						if (n > 1) {
							if (bx > ax || (bx == ax && bz > az))
								maxi = Math.floor((ai + n / 2) % pn);
							else
								maxi = Math.floor((ai + (n + 1) / 2) % pn);
						}
					}
				}

				// If the max deviation is larger than accepted error,
				// add new point, else continue to next segment.
				if (maxi != -1) {
					// Add the point.
					simplified.splice((i + 1) * 4 + 0, 0, points[maxi * 4 + 0]);
					simplified.splice((i + 1) * 4 + 1, 0, points[maxi * 4 + 1]);
					simplified.splice((i + 1) * 4 + 2, 0, points[maxi * 4 + 2]);
					simplified.splice((i + 1) * 4 + 3, 0, maxi);
				} else {
					++i;
				}
			}
		}
		for (let i = 0; i < simplified.length / 4; ++i) {
			// The edge vertex flag is take from the current raw point,
			// and the neighbour region is take from the next raw point.
			let ai = (simplified[i * 4 + 3] + 1) % pn;
			let bi = simplified[i * 4 + 3];
			simplified[i * 4 + 3] =
					(points[ai * 4 + 3]
						& (RecastConstants.RC_CONTOUR_REG_MASK | RecastConstants.RC_AREA_BORDER))
					| (points[bi * 4 + 3] & RecastConstants.RC_BORDER_VERTEX);
		}

	}

	static calcAreaOfPolygon2D(verts, nverts) {
		let area = 0;
		for (let i = 0, j = nverts - 1; i < nverts; j = i++) {
			let vi = i * 4;
			let vj = j * 4;
			area += verts[vi + 0] * verts[vj + 2] - verts[vj + 0] * verts[vi + 2];
		}
		return (area + 1) / 2;
	}

	static intersectSegCountour(d0, d1, i, n, verts, d0verts,
		d1verts) {
		// For each edge (k,k+1) of P
		let pverts = new Array(4 * 4);
		for(let g = 0; g < 4; g++) {
			pverts[g] = d0verts[d0 + g];
			pverts[4 + g] = d1verts[d1 + g];
		}
		d0 = 0;
		d1 = 4;
		for(let k = 0; k < n; k++) {
			let k1 = RecastMesh.next(k, n);
			// Skip edges incident to i.
			if (i == k || i == k1)
				continue;
			let p0 = k * 4;
			let p1 = k1 * 4;
			for(let g = 0; g < 4; g++) {
				pverts[8 + g] = verts[p0 + g];
				pverts[12 + g] = verts[p1 + g];
			}
			p0 = 8;
			p1 = 12;
			if (RecastMesh.vequal(pverts, d0, p0) || RecastMesh.vequal(pverts, d1, p0)
				|| RecastMesh.vequal(pverts, d0, p1) || RecastMesh.vequal(pverts, d1, p1))
				continue;

			if (RecastMesh.intersect(pverts, d0, d1, p0, p1))
				return true;
		}
		return false;
	}

	static inCone(i, n, verts, pj, vertpj) {
		pi = i * 4;
		pi1 = RecastMesh.next(i, n) * 4;
		pin1 = RecastMesh.prev(i, n) * 4;
		pverts = new Array(4 * 4);
		for(let g = 0; g < 4; g++) {
			pverts[g] = verts[pi + g];
			pverts[4 + g] = verts[pi1 + g];
			pverts[8 + g] = verts[pin1 + g];
			pverts[12 + g] = vertpj[pj + g];
		}
		pi = 0;
		pi1 = 4;
		pin1 = 8;
		pj = 12;
		// If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].
		if (RecastMesh.leftOn(pverts, pin1, pi, pi1))
			return RecastMesh.left(pverts, pi, pj, pin1) && RecastMesh.left(pverts, pj, pi, pi1);
		// Assume (i-1,i,i+1) not collinear.
		// else P[i] is reflex.
		return !(RecastMesh.leftOn(pverts, pi, pj, pi1) && RecastMesh.leftOn(pverts, pj, pi, pin1));
	}

	static removeDegenerateSegments(simplified) {
		// Remove adjacent vertices which are equal on xz-plane,
		// or else the triangulator will get confused.
		let npts = simplified.length / 4;
		for (let i = 0; i < npts; ++i) {
			let ni = RecastMesh.next(i, npts);

			//			if (vequal(&simplified[i*4], &simplified[ni*4]))
			if (simplified[i * 4] == simplified[ni * 4]
				&& simplified[i * 4 + 2] == simplified[ni * 4 + 2]) {
				// Degenerate segment, remove.
				simplified.splice(i * 4, 1);
				simplified.splice(i * 4, 1);
				simplified.splice(i * 4, 1);
				simplified.splice(i * 4, 1);
				npts--;
			}
		}
	}

	static mergeContours(ca, cb, ia, ib) {
		let maxVerts = ca.nverts + cb.nverts + 2;
		let verts = new Array(maxVerts * 4);

		let nv = 0;

		// Copy contour A.
		for (let i = 0; i <= ca.nverts; ++i) {
			let dst = nv * 4;
			let src = ((ia + i) % ca.nverts) * 4;
			verts[dst + 0] = ca.verts[src + 0];
			verts[dst + 1] = ca.verts[src + 1];
			verts[dst + 2] = ca.verts[src + 2];
			verts[dst + 3] = ca.verts[src + 3];
			nv++;
		}

		// Copy contour B
		for (let i = 0; i <= cb.nverts; ++i) {
			let dst = nv * 4;
			let src = ((ib + i) % cb.nverts) * 4;
			verts[dst + 0] = cb.verts[src + 0];
			verts[dst + 1] = cb.verts[src + 1];
			verts[dst + 2] = cb.verts[src + 2];
			verts[dst + 3] = cb.verts[src + 3];
			nv++;
		}

		ca.verts = verts;
		ca.nverts = nv;

		cb.verts = null;
		cb.nverts = 0;

	}

	// Finds the lowest leftmost vertex of a contour.
	static findLeftMostVertex(contour) {
		let minx = contour.verts[0];
		let minz = contour.verts[2];
		let leftmost = 0;
		for (let i = 1; i < contour.nverts; i++) {
			let x = contour.verts[i * 4 + 0];
			let z = contour.verts[i * 4 + 2];
			if (x < minx || (x == minx && z < minz)) {
				minx = x;
				minz = z;
				leftmost = i;
			}
		}
		return [minx, minz, leftmost];
	}

	CompareHoles = function (a, b) {
		if (a.minx == b.minx) {
			if (a.minz < b.minz)
				return -1;
			if (a.minz > b.minz)
				return 1;
		} else {
			if (a.minx < b.minx)
				return -1;
			if (a.minx > b.minx)
				return 1;
		}
		return 0;
	}



	CompareDiagDist = function (va, vb) {
		let a = va;
		let b = vb;
		if (a.dist < b.dist)
			return -1;
		if (a.dist > b.dist)
			return 1;
		return 0;
	}


	mergeRegionHoles(ctx, region) {
		// Sort holes from left to right.
		for (let i = 0; i < region.nholes; i++) {
			let minleft = findLeftMostVertex(region.holes[i].contour);
			region.holes[i].minx = minleft[0];
			region.holes[i].minz = minleft[1];
			region.holes[i].leftmost = minleft[2];
		}
		Arrays.sort(region.holes, new CompareHoles());

		let maxVerts = region.outline.nverts;
		for (let i = 0; i < region.nholes; i++)
			maxVerts += region.holes[i].contour.nverts;

		diags = new Array(maxVerts);
		for(let pd = 0; pd < maxVerts; pd++) {
			diags[pd] = new PotentialDiagonal();
		}
		let outline = region.outline;

		// Merge holes into the outline one by one.
		for (let i = 0; i < region.nholes; i++) {
			let hole = region.holes[i].contour;

			let index = -1;
			let bestVertex = region.holes[i].leftmost;
			for (let iter = 0; iter < hole.nverts; iter++) {
				// Find potential diagonals.
				// The 'best' vertex must be in the cone described by 3 cosequtive vertices of the outline.
				// ..o j-1
				//   |
				//   |   * best
				//   |
				// j o-----o j+1
				//         :
				let ndiags = 0;
				let corner = bestVertex * 4;
				for(let j = 0; j < outline.nverts; j++) {
					if (inCone(j, outline.nverts, outline.verts, corner, hole.verts)) {
						let dx = outline.verts[j * 4 + 0] - hole.verts[corner + 0];
						let dz = outline.verts[j * 4 + 2] - hole.verts[corner + 2];
						diags[ndiags].vert = j;
						diags[ndiags].dist = dx * dx + dz * dz;
						ndiags++;
					}
				}
				// Sort potential diagonals by distance, we want to make the connection as short as possible.
				Arrays.sort(diags, 0, ndiags, new CompareDiagDist());

				// Find a diagonal that is not intersecting the outline not the remaining holes.
				index = -1;
				for(let j = 0; j < ndiags; j++) {
					let pt = diags[j].vert * 4;
					let intersect = intersectSegCountour(pt, corner, diags[i].vert, outline.nverts, outline.verts,
						outline.verts, hole.verts);
					for(let k = i; k < region.nholes && !intersect; k++)
						intersect |= intersectSegCountour(pt, corner, -1, region.holes[k].contour.nverts,
							region.holes[k].contour.verts, outline.verts, hole.verts);
					if (!intersect) {
						index = diags[j].vert;
						break;
					}
				}
				// If found non-intersecting diagonal, stop looking.
				if (index != -1)
					break;
				// All the potential diagonals for the current vertex were intersecting, try next vertex.
				bestVertex = (bestVertex + 1) % hole.nverts;
			}

			if (index == -1) {
				ctx.warn("mergeHoles: Failed to find merge points for");
				continue;
			}
			mergeContours(region.outline, hole, index, bestVertex);
		}
	}

	/// @par
	///
	/// The raw contours will match the region outlines exactly. The @p maxError and @p maxEdgeLen
	/// parameters control how closely the simplified contours will match the raw contours.
	///
	/// Simplified contours are generated such that the vertices for portals between areas match up.
	/// (They are considered mandatory vertices.)
	///
	/// Setting @p maxEdgeLength to zero will disabled the edge length feature.
	///
	/// See the #rcConfig documentation for more information on the configuration parameters.
	///
	/// @see rcAllocContourSet, rcCompactHeightfield, rcContourSet, rcConfig
	static buildContours(ctx, chf, maxError, maxEdgeLen,
		buildFlags) {

		let w = chf.width;
		let h = chf.height;
		let borderSize = chf.borderSize;
		let cset = new ContourSet();

		ctx.startTimer("BUILD_CONTOURS");
		RecastVectors.copy3(cset.bmin, chf.bmin, 0);
		RecastVectors.copy3(cset.bmax, chf.bmax, 0);
		if (borderSize > 0) {
			// If the heightfield was build with bordersize, remove the offset.
			pad = borderSize * chf.cs;
			cset.bmin[0] += pad;
			cset.bmin[2] += pad;
			cset.bmax[0] -= pad;
			cset.bmax[2] -= pad;
		}
		cset.cs = chf.cs;
		cset.ch = chf.ch;
		cset.width = chf.width - chf.borderSize * 2;
		cset.height = chf.height - chf.borderSize * 2;
		cset.borderSize = chf.borderSize;
		cset.maxError = maxError;

		let flags = new Array(chf.spanCount).fill(0);

		ctx.startTimer("BUILD_CONTOURS_TRACE");

		// Mark boundaries.
		for (let y = 0; y < h; ++y) {
			for (let x = 0; x < w; ++x) {
				let c = chf.cells[x + y * w];
				for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
					// if(y == 3 && x == 3 && i == 1344)
					// 	console.log("earlier");
					let res = 0;
					let s = chf.spans[i];
					if (chf.spans[i].reg == 0 || (chf.spans[i].reg & RecastConstants.RC_BORDER_REG) != 0) {
						flags[i] = 0;
						continue;
					}
					for (let dir = 0; dir < 4; ++dir) {
						let r = 0;
						if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
							let ax = x + RecastCommon.GetDirOffsetX(dir);
							let ay = y + RecastCommon.GetDirOffsetY(dir);
							let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);
							r = chf.spans[ai].reg;
						}
						if (r == chf.spans[i].reg)
							res |= (1 << dir);
					}
					flags[i] = res ^ 0xf; // Inverse, mark non connected edges.
				}
			}
		}

		ctx.stopTimer("BUILD_CONTOURS_TRACE");

		let verts = new Array(256);
		let simplified = new Array(64);;

		for (let y = 0; y < h; ++y) {
			for (let x = 0; x < w; ++x) {
				let c = chf.cells[x + y * w];
				for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
					// if(y==3 && x == 3 && i == 1344)
					// 	console.log("dak");
					if (flags[i] == 0 || flags[i] == 0xf) {
						flags[i] = 0;
						continue;
					}
					let reg = chf.spans[i].reg;
					if (reg == 0 || (reg & RecastConstants.RC_BORDER_REG) != 0)
						continue;
					let area = chf.areas[i];

					verts = [];
					simplified = [];

					ctx.startTimer("BUILD_CONTOURS_TRACE");
					RecastContour.walkContour(x, y, i, chf, flags, verts);
					ctx.stopTimer("BUILD_CONTOURS_TRACE");

					ctx.startTimer("BUILD_CONTOURS_SIMPLIFY");
					RecastContour.simplifyContour(verts, simplified, maxError, maxEdgeLen, buildFlags);
					RecastContour.removeDegenerateSegments(simplified);
					ctx.stopTimer("BUILD_CONTOURS_SIMPLIFY");

					// Store region=>contour remap info.
					// Create contour.
					if (simplified.length / 4 >= 3) {

						let cont = new Contour();
						cset.conts.push(cont);

						cont.nverts = simplified.length / 4;
						cont.verts = new Array(simplified.length);
						for (let l = 0; l < cont.verts.length; l++) {
							cont.verts[l] = simplified[l];
						}

						if (borderSize > 0) {
							// If the heightfield was build with bordersize, remove the offset.
							for (let j = 0; j < cont.nverts; ++j) {
								cont.verts[j * 4] -= borderSize;
								cont.verts[j * 4 + 2] -= borderSize;
							}
						}

						cont.nrverts = verts.length / 4;
						cont.rverts = new Array(verts.length);
						for (let l = 0; l < cont.rverts.length; l++) {
							cont.rverts[l] = verts[l];
						}
						if (borderSize > 0) {
							// If the heightfield was build with bordersize, remove the offset.
							for (let j = 0; j < cont.nrverts; ++j) {
								cont.rverts[j * 4] -= borderSize;
								cont.rverts[j * 4 + 2] -= borderSize;
							}
						}

						cont.reg = reg;
						cont.area = area;
					}
				}
			}
		}

		// Merge holes if needed.
		if (cset.conts.length > 0) {
			// Calculate winding of all polygons.
			let winding = new Array(cset.conts.length);
			let nholes = 0;
			for (let i = 0; i < cset.conts.length; ++i) {
				let cont = cset.conts[i];
				// If the contour is wound backwards, it is a hole.
				winding[i] = RecastContour.calcAreaOfPolygon2D(cont.verts, cont.nverts) < 0 ? -1 : 1;
				if (winding[i] < 0)
					nholes++;
			}

			if (nholes > 0) {
				// Collect outline contour and holes contours per region.
				// We assume that there is one outline and multiple holes.
				let nregions = chf.maxRegions + 1;
				let regions = new Array(nregions);
				for (let i = 0; i < nregions; i++) {
					regions[i] = new ContourRegion();
				}

				for (let i = 0; i < cset.conts.length; ++i) {
					let cont = cset.conts[i];
					// Positively would contours are outlines, negative holes.
					if (winding[i] > 0) {
						if (regions[cont.reg].outline != null) {
							throw new RuntimeException(
								"rcBuildContours: Multiple outlines for region " + cont.reg + ".");
						}
						regions[cont.reg].outline = cont;
					} else {
						regions[cont.reg].nholes++;
					}
				}
				for (let i = 0; i < nregions; i++) {
					if (regions[i].nholes > 0) {
						regions[i].holes = new ContourHole[regions[i].nholes];
						for(let nh = 0; nh < regions[i].nholes; nh++) {
							regions[i].holes[nh] = new ContourHole();
						}
						regions[i].nholes = 0;
					}
				}
				for (let i = 0; i < cset.conts.length; ++i) {
					let cont = cset.conts[i];
					let reg = regions[cont.reg];
					if (winding[i] < 0)
						reg.holes[reg.nholes++].contour = cont;
				}

				// Finally merge each regions holes into the outline.
				for (let i = 0; i < nregions; i++) {
					let reg = regions[i];
					if (reg.nholes == 0)
						continue;

					if (reg.outline != null) {
						mergeRegionHoles(ctx, reg);
					} else {
						// The region does not have an outline.
						// This can happen if the contour becaomes selfoverlapping because of
						// too aggressive simplification settings.
						throw new RuntimeException("rcBuildContours: Bad outline for region " + i
							+ ", contour simplification is likely too aggressive.");
					}
				}
			}
		}
		ctx.stopTimer("BUILD_CONTOURS");
		return cset;
	}
}

export default RecastContour;