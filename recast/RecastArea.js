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

import RecastConstants from "./RecastConstants.js"
import RecastCommon from "./RecastCommon.js"



class RecastArea {

	/// @par
	///
	/// Basically, any spans that are closer to a boundary or obstruction than the specified radius
	/// are marked as unwalkable.
	///
	/// This method is usually called immediately after the heightfield has been built.
	///
	/// @see rcCompactHeightfield, rcBuildCompactHeightfield, rcConfig::walkableRadius
static erodeWalkableArea(ctx,  radius, chf) {
		let w = chf.width;
		let h = chf.height;
		ctx.startTimer("ERODE_AREA");

		let dist = new Array(chf.spanCount).fill(255);
		// Arrays.fill(dist, 255);
		// Mark boundary cells.
		for (let y = 0; y < h; ++y) {
			for (let x = 0; x < w; ++x) {
				let c = chf.cells[x + y * w];
				for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
					if (chf.areas[i] == RecastConstants.RC_NULL_AREA) {
						dist[i] = 0;
					} else {
						let s = chf.spans[i];
						let nc = 0;
						for (let dir = 0; dir < 4; ++dir) {
							if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
								let nx = x + RecastCommon.GetDirOffsetX(dir);
								let ny = y + RecastCommon.GetDirOffsetY(dir);
								let nidx = chf.cells[nx + ny * w].index + RecastCommon.GetCon(s, dir);
								if (chf.areas[nidx] != RecastConstants.RC_NULL_AREA) {
									nc++;
								}
							}
						}
						// At least one missing neighbour.
						if (nc != 4)
							dist[i] = 0;
					}
				}
			}
		}

		let nd;
		// console.log(chf.spans[450240]);

		// Pass 1
		for (let y = 0; y < h; ++y) {
			for (let x = 0; x < w; ++x) {
				let c = chf.cells[x + y * w];
				for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
					let s = chf.spans[i];

					if (RecastCommon.GetCon(s, 0) != RecastConstants.RC_NOT_CONNECTED) {
						// (-1,0)
						let ax = x + RecastCommon.GetDirOffsetX(0);
						let ay = y + RecastCommon.GetDirOffsetY(0);
						let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 0);
						let as = chf.spans[ai];
						nd = Math.min(dist[ai] + 2, 255);
						if (nd < dist[i])
							dist[i] = nd;

						// (-1,-1)
						if (RecastCommon.GetCon(as, 3) != RecastConstants.RC_NOT_CONNECTED) {
							let aax = ax + RecastCommon.GetDirOffsetX(3);
							let aay = ay + RecastCommon.GetDirOffsetY(3);
							let aai = chf.cells[aax + aay * w].index + RecastCommon.GetCon(as, 3);
							nd = Math.min(dist[aai] + 3, 255);
							if (nd < dist[i])
								dist[i] = nd;
						}
					}
					if (RecastCommon.GetCon(s, 3) != RecastConstants.RC_NOT_CONNECTED) {
						// (0,-1)
						let ax = x + RecastCommon.GetDirOffsetX(3);
						let ay = y + RecastCommon.GetDirOffsetY(3);
						let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 3);
						let as = chf.spans[ai];
						nd = Math.min(dist[ai] + 2, 255);
						if (nd < dist[i])
							dist[i] = nd;

						// (1,-1)
						if (RecastCommon.GetCon(as, 2) != RecastConstants.RC_NOT_CONNECTED) {
							let aax = ax + RecastCommon.GetDirOffsetX(2);
							let aay = ay + RecastCommon.GetDirOffsetY(2);
							let aai = chf.cells[aax + aay * w].index + RecastCommon.GetCon(as, 2);
							nd = Math.min(dist[aai] + 3, 255);
							if (nd < dist[i])
								dist[i] = nd;
						}
					}
				}
			}
		}

		// Pass 2
		for (let y = h - 1; y >= 0; --y) {
			for (let x = w - 1; x >= 0; --x) {
				// if(y == 671&& x == 671)
				// 	console.log("There")
				let c = chf.cells[x + y * w];
				for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
					let s = chf.spans[i];

					if (RecastCommon.GetCon(s, 2) != RecastConstants.RC_NOT_CONNECTED) {
						// (1,0)
						let ax = x + RecastCommon.GetDirOffsetX(2);
						let ay = y + RecastCommon.GetDirOffsetY(2);
						let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 2);
						if(ai == 450241)
							console.log("Here")
						let as = chf.spans[ai];
						nd = Math.min(dist[ai] + 2, 255);
						if (nd < dist[i])
							dist[i] = nd;

						// (1,1)
						if (RecastCommon.GetCon(as, 1) != RecastConstants.RC_NOT_CONNECTED) {
							let aax = ax + RecastCommon.GetDirOffsetX(1);
							let aay = ay + RecastCommon.GetDirOffsetY(1);
							let aai = chf.cells[aax + aay * w].index + RecastCommon.GetCon(as, 1);
							nd = Math.min(dist[aai] + 3, 255);
							if (nd < dist[i])
								dist[i] = nd;
						}
					}
					if (RecastCommon.GetCon(s, 1) != RecastConstants.RC_NOT_CONNECTED) {
						// (0,1)
						let ax = x + RecastCommon.GetDirOffsetX(1);
						let ay = y + RecastCommon.GetDirOffsetY(1);
						let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 1);
						let as = chf.spans[ai];
						nd = Math.min(dist[ai] + 2, 255);
						if (nd < dist[i])
							dist[i] = nd;

						// (-1,1)
						if (RecastCommon.GetCon(as, 0) != RecastConstants.RC_NOT_CONNECTED) {
							let aax = ax + RecastCommon.GetDirOffsetX(0);
							let aay = ay + RecastCommon.GetDirOffsetY(0);
							let aai = chf.cells[aax + aay * w].index + RecastCommon.GetCon(as, 0);
							nd = Math.min(dist[aai] + 3, 255);
							if (nd < dist[i])
								dist[i] = nd;
						}
					}
				}
			}
		}

		let thr = radius * 2;
		for (let i = 0; i < chf.spanCount; ++i)
			if (dist[i] < thr)
				chf.areas[i] = RecastConstants.RC_NULL_AREA;

		ctx.stopTimer("ERODE_AREA");
	}

	/// @par
	///
	/// This filter is usually applied after applying area id's using functions
	/// such as #rcMarkBoxArea, #rcMarkConvexPolyArea, and #rcMarkCylinderArea.
	///
	/// @see rcCompactHeightfield
medianFilterWalkableArea(ctx, chf) {

		let w = chf.width;
		let h = chf.height;

		ctx.startTimer("MEDIAN_AREA");

		let areas = new Array(chf.spanCount);

		for(let y = 0; y < h; ++y) {
			for(let x = 0; x < w; ++x) {
				let c = chf.cells[x + y * w];
				for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
					s = chf.spans[i];
					if (chf.areas[i] == RecastConstants.RC_NULL_AREA) {
						areas[i] = chf.areas[i];
						continue;
					}

					nei = new Array(9);
					for(let j = 0; j < 9; ++j)
						nei[j] = chf.areas[i];

					for(let dir = 0; dir < 4; ++dir) {
						if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
							let ax = x + RecastCommon.GetDirOffsetX(dir);
							let ay = y + RecastCommon.GetDirOffsetY(dir);
							let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);
							if (chf.areas[ai] != RecastConstants.RC_NULL_AREA)
								nei[dir * 2 + 0] = chf.areas[ai];

							let as = chf.spans[ai];
							let dir2 = (dir + 1) & 0x3;
							if (RecastCommon.GetCon(as, dir2) != RecastConstants.RC_NOT_CONNECTED) {
								let ax2 = ax + RecastCommon.GetDirOffsetX(dir2);
								let ay2 = ay + RecastCommon.GetDirOffsetY(dir2);
								let ai2 = chf.cells[ax2 + ay2 * w].index + RecastCommon.GetCon(as, dir2);
								if (chf.areas[ai2] != RecastConstants.RC_NULL_AREA)
									nei[dir * 2 + 1] = chf.areas[ai2];
							}
						}
					}
					Arrays.sort(nei);
					areas[i] = nei[4];
				}
			}
		}
		chf.areas = areas;

		ctx.stopTimer("MEDIAN_AREA");

		return true;
	}

	/// @par
	///
	/// The value of spacial parameters are in world units.
	///
	/// @see rcCompactHeightfield, rcMedianFilterWalkableArea
markBoxArea(ctx, bmin, bmax, areaMod, chf) {
		ctx.startTimer("MARK_BOX_AREA");

		let minx =  Math.floor((bmin[0] - chf.bmin[0]) / chf.cs);
		let miny =  Math.floor((bmin[1] - chf.bmin[1]) / chf.ch);
		let minz =  Math.floor((bmin[2] - chf.bmin[2]) / chf.cs);
		let maxx =  Math.floor((bmax[0] - chf.bmin[0]) / chf.cs);
		let maxy =  Math.floor((bmax[1] - chf.bmin[1]) / chf.ch);
		let maxz =  Math.floor((bmax[2] - chf.bmin[2]) / chf.cs);

		if (maxx < 0)
			return;
		if (minx >= chf.width)
			return;
		if (maxz < 0)
			return;
		if (minz >= chf.height)
			return;

		if (minx < 0)
			minx = 0;
		if (maxx >= chf.width)
			maxx = chf.width - 1;
		if (minz < 0)
			minz = 0;
		if (maxz >= chf.height)
			maxz = chf.height - 1;

		for (let z = minz; z <= maxz; ++z) {
			for (let x = minx; x <= maxx; ++x) {
				let c = chf.cells[x + z * chf.width];
				for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
					let  s = chf.spans[i];
					if (s.y >= miny && s.y <= maxy) {
						if (chf.areas[i] != RecastConstants.RC_NULL_AREA)
							chf.areas[i] = areaMod.apply(chf.areas[i]);
					}
				}
			}
		}

		ctx.stopTimer("MARK_BOX_AREA");

	}

	static pointInPoly( verts, p) {
		let c = false;
		let i, j;
		for (let i = 0, j = verts.length - 1; i < verts.length; j = i++) {
			let vi = i * 3;
			let vj = j * 3;
			if (((verts[vi + 2] > p[2]) != (verts[vj + 2] > p[2]))
					&& (p[0] < (verts[vj] - verts[vi]) * (p[2] - verts[vi + 2]) / (verts[vj + 2] - verts[vi + 2])
							+ verts[vi]))
				c = !c;
		}
		return c;
	}

	/// @par
	///
	/// The value of spacial parameters are in world units.
	///
	/// The y-values of the polygon vertices are ignored. So the polygon is effectively
	/// projected onto the xz-plane at @p hmin, then extruded to @p hmax.
	///
	/// @see rcCompactHeightfield, rcMedianFilterWalkableArea
static markConvexPolyArea(ctx, verts,  hmin,  hmax, areaMod,
			chf) {
		ctx.startTimer("MARK_CONVEXPOLY_AREA");

		 bmin = new Array(3);
		 max = new Array(3);
		RecastVectors.copy(bmin, verts, 0);
		RecastVectors.copy(bmax, verts, 0);
		for (let i = 1; i < verts.length; ++i) {
			RecastVectors.min(bmin, verts, i * 3);
			RecastVectors.max(bmax, verts, i * 3);
		}
		bmin[1] = hmin;
		bmax[1] = hmax;

		let minx =  Math.floor((bmin[0] - chf.bmin[0]) / chf.cs);
		let miny =  Math.floor((bmin[1] - chf.bmin[1]) / chf.ch);
		let minz =  Math.floor((bmin[2] - chf.bmin[2]) / chf.cs);
		let maxx =  Math.floor((bmax[0] - chf.bmin[0]) / chf.cs);
		let maxy =  Math.floor((bmax[1] - chf.bmin[1]) / chf.ch);
		let maxz =  Math.floor((bmax[2] - chf.bmin[2]) / chf.cs);

		if (maxx < 0)
			return;
		if (minx >= chf.width)
			return;
		if (maxz < 0)
			return;
		if (minz >= chf.height)
			return;

		if (minx < 0)
			minx = 0;
		if (maxx >= chf.width)
			maxx = chf.width - 1;
		if (minz < 0)
			minz = 0;
		if (maxz >= chf.height)
			maxz = chf.height - 1;

		// TODO: Optimize.
		for (let z = minz; z <= maxz; ++z) {
			for (let x = minx; x <= maxx; ++x) {
				let c = chf.cells[x + z * chf.width];
				for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
					let s = chf.spans[i];
					if (chf.areas[i] == RecastConstants.RC_NULL_AREA)
						continue;
					if (s.y >= miny && s.y <= maxy) {
						 p = new Array(3);
						p[0] = chf.bmin[0] + (x +0.5) * chf.cs;
						p[1] = 0;
						p[2] = chf.bmin[2] + (z +0.5) * chf.cs;

						if (pointInPoly(verts, p)) {
							chf.areas[i] = areaMod.apply(chf.areas[i]);
						}
					}
				}
			}
		}

		ctx.stopTimer("MARK_CONVEXPOLY_AREA");
	}

	offsetPoly( verts,  nverts,  offset, outVerts,  maxOutVerts) {
		 MITER_LIMIT = 1.20;

		let n = 0;

		for (let i = 0; i < nverts; i++) {
			let a = (i + nverts - 1) % nverts;
			let b = i;
			let c = (i + 1) % nverts;
			let va = a * 3;
			let vb = b * 3;
			let vc = c * 3;
			 dx0 = verts[vb] - verts[va];
			 dy0 = verts[vb + 2] - verts[va + 2];
			 d0 = dx0 * dx0 + dy0 * dy0;
			if (d0 > 1e-6) {
				d0 =  (1.0 / Math.sqrt(d0));
				dx0 *= d0;
				dy0 *= d0;
			}
			 dx1 = verts[vc] - verts[vb];
			 dy1 = verts[vc + 2] - verts[vb + 2];
			 d1 = dx1 * dx1 + dy1 * dy1;
			if (d1 > 1e-6) {
				d1 =  (1.0 / Math.sqrt(d1));
				dx1 *= d1;
				dy1 *= d1;
			}
			 dlx0 = -dy0;
			 dly0 = dx0;
			 dlx1 = -dy1;
			 dly1 = dx1;
			 cross = dx1 * dy0 - dx0 * dy1;
			 dmx = (dlx0 + dlx1) * 0.5;
			 dmy = (dly0 + dly1) * 0.5;
			 dmr2 = dmx * dmx + dmy * dmy;
			let bevel = dmr2 * MITER_LIMIT * MITER_LIMIT < 1.0;
			if (dmr2 > 1e-6) {
				 scale = 1.0 / dmr2;
				dmx *= scale;
				dmy *= scale;
			}

			if (bevel && cross < 0.0) {
				if (n + 2 >= maxOutVerts)
					return 0;
				 d = (1.0 - (dx0 * dx1 + dy0 * dy1)) * 0.5;
				outVerts[n * 3 + 0] = verts[vb] + (-dlx0 + dx0 * d) * offset;
				outVerts[n * 3 + 1] = verts[vb + 1];
				outVerts[n * 3 + 2] = verts[vb + 2] + (-dly0 + dy0 * d) * offset;
				n++;
				outVerts[n * 3 + 0] = verts[vb] + (-dlx1 - dx1 * d) * offset;
				outVerts[n * 3 + 1] = verts[vb + 1];
				outVerts[n * 3 + 2] = verts[vb + 2] + (-dly1 - dy1 * d) * offset;
				n++;
			} else {
				if (n + 1 >= maxOutVerts)
					return 0;
				outVerts[n * 3 + 0] = verts[vb] - dmx * offset;
				outVerts[n * 3 + 1] = verts[vb + 1];
				outVerts[n * 3 + 2] = verts[vb + 2] - dmy * offset;
				n++;
			}
		}

		return n;
	}

	/// @par
	///
	/// The value of spacial parameters are in world units.
	///
	/// @see rcCompactHeightfield, rcMedianFilterWalkableArea
markCylinderArea(ctx, pos,  r,  h, areaMod,
			chf) {

		ctx.startTimer("MARK_CYLINDER_AREA");

		 bmin = new Array(3);
		 bmax= new Array(3);
		bmin[0] = pos[0] - r;
		bmin[1] = pos[1];
		bmin[2] = pos[2] - r;
		bmax[0] = pos[0] + r;
		bmax[1] = pos[1] + h;
		bmax[2] = pos[2] + r;
		 r2 = r * r;

		let minx =  Math.floor((bmin[0] - chf.bmin[0]) / chf.cs);
		let miny =  Math.floor((bmin[1] - chf.bmin[1]) / chf.ch);
		let minz =  Math.floor((bmin[2] - chf.bmin[2]) / chf.cs);
		let maxx =  Math.floor((bmax[0] - chf.bmin[0]) / chf.cs);
		let maxy =  Math.floor((bmax[1] - chf.bmin[1]) / chf.ch);
		let maxz =  Math.floor((bmax[2] - chf.bmin[2]) / chf.cs);

		if (maxx < 0)
			return;
		if (minx >= chf.width)
			return;
		if (maxz < 0)
			return;
		if (minz >= chf.height)
			return;

		if (minx < 0)
			minx = 0;
		if (maxx >= chf.width)
			maxx = chf.width - 1;
		if (minz < 0)
			minz = 0;
		if (maxz >= chf.height)
			maxz = chf.height - 1;

		for (let z = minz; z <= maxz; ++z) {
			for (let x = minx; x <= maxx; ++x) {
				let c = chf.cells[x + z * chf.width];
				for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
					let s = chf.spans[i];

					if (chf.areas[i] == RecastConstants.RC_NULL_AREA)
						continue;

					if (s.y >= miny && s.y <= maxy) {
						 sx = chf.bmin[0] + (x +0.5) * chf.cs;
						 sz = chf.bmin[2] + (z +0.5) * chf.cs;
						 dx = sx - pos[0];
						 dz = sz - pos[2];

						if (dx * dx + dz * dz < r2) {
							chf.areas[i] = areaMod.apply(chf.areas[i]);
						}
					}
				}
			}
		}
		ctx.stopTimer("MARK_CYLINDER_AREA");
	}

}

export default RecastArea;
