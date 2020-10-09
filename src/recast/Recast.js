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
import CompactHeightfield from "./CompactHeightfield.js"
import RecastConstants from "./RecastConstants.js"
import CompactCell from "./CompactCell.js"
import CompactSpan from "./CompactSpan.js"
import RecastCommon from "./RecastCommon.js"

class Recast {

	calcBounds( verts,  nv, bmin, bmax) {
		for (let i = 0; i < 3; i++) {
			bmin[i] = verts[i];
			bmax[i] = verts[i];
		}
		for (let i = 1; i < nv; ++i) {
			for(let j = 0; j < 3; j++) {
				bmin[j] = Math.min(bmin[j], verts[i * 3 + j]);
				bmax[j] = Math.max(bmax[j], verts[i * 3 + j]);
			}
		}
		// Calculate bounding box.
	}

static calcGridSize( bmin, bmax,  cs) {
		return [  Math.floor((bmax[0] - bmin[0]) / cs +0.5),  Math.floor((bmax[2] - bmin[2]) / cs +0.5) ];
	}

static calcTileCount( bmin, bmax,  cs,  tileSize) {
		let gwh = Recast.calcGridSize(bmin, bmax, cs);
		let gw = gwh[0];
		let gh = gwh[1];
		let ts = tileSize;
		let tw = (gw + ts - 1) / ts;
		let th = (gh + ts - 1) / ts;
		return [ tw, th ];
	}

	/// @par
	///
	/// Modifies the area id of all triangles with a slope below the specified value.
	///
	/// See the #rcConfig documentation for more information on the configuration parameters.
	///
	/// @see rcHeightfield, rcClearUnwalkableTriangles, rcRasterizeTriangles
static markWalkableTriangles(ctx,  walkableSlopeAngle, verts, tris,  nt,
			areaMod) {
		let areas = new  Array(nt).fill(0);
		 let walkableThr = Math.cos(walkableSlopeAngle / 180.0 * Math.PI);
		 let norm = new Array(3).fill(0);
		for (let i = 0; i < nt; ++i) {
			let tri = i * 3;
			Recast.calcTriNormal(verts, tris[tri], tris[tri + 1], tris[tri + 2], norm);
			// Check if the face is walkable.
			if (norm[1] > walkableThr)
				areas[i] = areaMod.apply(areas[i]);
		}
		return areas;
	}

	static calcTriNormal( verts,  v0,  v1,  v2, norm) {
		 let e0 = new Array(3);
		 let e1 = new Array(3);
		RecastVectors.subA(e0, verts, v1 * 3, v0 * 3);
		RecastVectors.subA(e1, verts, v2 * 3, v0 * 3);
		RecastVectors.cross(norm, e0, e1);
		RecastVectors.normalize(norm);
	}

	/// @par
	///
	/// Only sets the area id's for the unwalkable triangles. Does not alter the
	/// area id's for walkable triangles.
	///
	/// See the #rcConfig documentation for more information on the configuration parameters.
	///
	/// @see rcHeightfield, rcClearUnwalkableTriangles, rcRasterizeTriangles
	clearUnwalkableTriangles(ctx,  walkableSlopeAngle, verts,  nv, tris,  nt,
			areas) {
		 walkableThr = Math.cos(walkableSlopeAngle / 180.0 * Math.PI);

		 norm = new Array(3);

		for (let i = 0; i < nt; ++i) {
			let tri = i * 3;
			calcTriNormal(verts, tris[tri], tris[tri + 1], tris[tri + 2], norm);
			// Check if the face is walkable.
			if (norm[1] <= walkableThr)
				areas[i] = RecastConstants.RC_NULL_AREA;
		}
	}

	static getHeightFieldSpanCount(ctx, hf) {
		let w = hf.width;
		let h = hf.height;
		let spanCount = 0;
		for (let y = 0; y < h; ++y) {
			for (let x = 0; x < w; ++x) {
				for (let s = hf.spans[x + y * w]; s != null; s = s.next) {
					if (s.area != RecastConstants.RC_NULL_AREA)
						spanCount++;
				}
			}
		}
		return spanCount;
	}

	/// @par
	///
	/// This is just the beginning of the process of fully building a compact heightfield.
	/// Various filters may be applied, then the distance field and regions built.
	/// E.g: #rcBuildDistanceField and #rcBuildRegions
	///
	/// See the #rcConfig documentation for more information on the configuration parameters.
	///
	/// @see rcAllocCompactHeightfield, rcHeightfield, rcCompactHeightfield, rcConfig

static buildCompactHeightfield(ctx,  walkableHeight,  walkableClimb,
			 hf) {

		ctx.startTimer("BUILD_COMPACTHEIGHTFIELD");

		let chf = new CompactHeightfield();
		let w = hf.width;
		let h = hf.height;
		let spanCount = this.getHeightFieldSpanCount(ctx, hf);

		// Fill in header.
		chf.width = w;
		chf.height = h;
		chf.spanCount = spanCount;
		chf.walkableHeight = walkableHeight;
		chf.walkableClimb = walkableClimb;
		chf.maxRegions = 0;
		RecastVectors.copy2(chf.bmin, hf.bmin);
		RecastVectors.copy2(chf.bmax, hf.bmax);
		chf.bmax[1] += walkableHeight * hf.ch;
		chf.cs = hf.cs;
		chf.ch = hf.ch;
		// let bigSize = w*h;
		chf.cells = new Array(Math.floor(w)*Math.floor(h));
		chf.spans = new Array(spanCount);
		chf.areas = new Array(spanCount);
		let MAX_HEIGHT = 0xffff;
		for (let i = 0; i < chf.cells.length; i++) {
			chf.cells[i] = new CompactCell();
		}
		for (let i = 0; i < chf.spans.length; i++) {
			chf.spans[i] = new CompactSpan();
		}
		// Fill in cells and spans.
		let idx = 0;
		for (let y = 0; y < h; ++y) {
			for (let x = 0; x < w; ++x) {
				let s = hf.spans[x + y * w];
				// If there are no spans at this cell, just leave the data to index=0, count=0.
				if (s == null)
					continue;
				let c = chf.cells[x + y * w];
				c.index = idx;
				c.count = 0;
				while (s != null) {
					if (s.area != RecastConstants.RC_NULL_AREA) {
						let bot = s.smax;
						let top = s.next != null ?  Math.floor(s.next.smin) : MAX_HEIGHT;
						chf.spans[idx].y = RecastCommon.clamp(bot, 0, 0xffff);
						chf.spans[idx].h = RecastCommon.clamp(top - bot, 0, 0xff);
						chf.areas[idx] = s.area;
						idx++;
						// if(idx == 450240)
						// 	console.log("here")
						c.count++;
					}
					s = s.next;
				}
			}
		}

		// Find neighbour connections.
		let MAX_LAYERS = RecastConstants.RC_NOT_CONNECTED - 1;
		let tooHighNeighbour = 0;
		for (var y = 0; y < h; ++y) {
			for (var x = 0; x < w; ++x) {
				let c = chf.cells[x + y * w];
				for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
					// if(i == 450240)
					// 	console.log("stop")
					let s = chf.spans[i];

					for (var dir = 0; dir < 4; ++dir) {
						RecastCommon.SetCon(s, dir, RecastConstants.RC_NOT_CONNECTED);
						let nx = x + RecastCommon.GetDirOffsetX(dir);
						let ny = y + RecastCommon.GetDirOffsetY(dir);
						// First check that the neighbour cell is in bounds.
						if (nx < 0 || ny < 0 || nx >= w || ny >= h)
							continue;

						// Iterate over all neighbour spans and check if any of the is
						// accessible from current cell.
						let nc = chf.cells[nx + ny * w];
						for (let k = nc.index, nk = nc.index + nc.count; k < nk; ++k) {
							let ns = chf.spans[k];
							let bot = Math.max(s.y, ns.y);
							let top = Math.min(s.y + s.h, ns.y + ns.h);

							// Check that the gap between the spans is walkable,
							// and that the climb height between the gaps is not too high.
							if ((top - bot) >= walkableHeight && Math.abs(ns.y - s.y) <= walkableClimb) {
								// Mark direction as walkable.
								let lidx = k - nc.index;
								if (lidx < 0 || lidx > MAX_LAYERS) {
									tooHighNeighbour = Math.max(tooHighNeighbour, lidx);
									continue;
								}
								RecastCommon.SetCon(s, dir, lidx);
								break;
							}
						}

					}
				}
			}
		}

		if (tooHighNeighbour > MAX_LAYERS) {
			throw new RuntimeException("rcBuildCompactHeightfield: Heightfield has too many layers " + tooHighNeighbour
					+ " (max: " + MAX_LAYERS + ")");
		}
		ctx.stopTimer("BUILD_COMPACTHEIGHTFIELD");
		return chf;
	}
}

export default Recast;
