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

import RecastVectors from "./RecastVectors.js"
import RecastCommon from "./RecastCommon.js"
import RecastConstants from "./RecastConstants.js"
import Span from "./Span.js"

class RecastRasterization {

	static overlapBounds(amin, amax, bmin, bmax) {
		let overlap = true;
		overlap = (amin[0] > bmax[0] || amax[0] < bmin[0]) ? false : overlap;
		overlap = (amin[1] > bmax[1] || amax[1] < bmin[1]) ? false : overlap;
		overlap = (amin[2] > bmax[2] || amax[2] < bmin[2]) ? false : overlap;
		return overlap;
	}

	/**
	 * The span addition can be set to favor flags. If the span is merged to another span and the new 'smax' is
	 * within 'flagMergeThr' units from the existing span, the span flags are merged.
	 * 
	 * @see Heightfield, Span.
	 */
	static addSpan(hf, x, y, smin, smax, area, flagMergeThr) {

		let idx = x + y * hf.width;

		let s = new Span();
		s.smin = smin;
		s.smax = smax;
		s.area = area;
		s.next = null;

		// Empty cell, add the first span.
		if (hf.spans[idx] == null) {
			hf.spans[idx] = s;
			return;
		}
		let prev = null;
		let cur = hf.spans[idx];

		// Insert and merge spans.
		while (cur != null) {
			if (cur.smin > s.smax) {
				// Current span is further than the new span, break.
				break;
			} else if (cur.smax < s.smin) {
				// Current span is before the new span advance.
				prev = cur;
				cur = cur.next;
			} else {
				// Merge spans.
				if (cur.smin < s.smin)
					s.smin = cur.smin;
				if (cur.smax > s.smax)
					s.smax = cur.smax;

				// Merge flags.
				if (Math.abs(s.smax - cur.smax) <= flagMergeThr)
					s.area = Math.max(s.area, cur.area);

				// Remove current span.
				let next = cur.next;
				if (prev != null)
					prev.next = next;
				else
					hf.spans[idx] = next;
				cur = next;
			}
		}

		// Insert new span.
		if (prev != null) {
			s.next = prev.next;
			prev.next = s;
		} else {
			s.next = hf.spans[idx];
			hf.spans[idx] = s;
		}
	}

	//divides a convex polygons into two convex polygons on both sides of a line
	static dividePoly(buf, _in, nin, out1, out2, x, axis) {
		let d = new Array(12);
		for (let i = 0; i < nin; ++i)
			d[i] = x - buf[_in + i * 3 + axis];

		let m = 0, n = 0;
		for (let i = 0, j = nin - 1; i < nin; j = i, ++i) {
			let ina = d[j] >= 0;
			let inb = d[i] >= 0;
			if (ina != inb) {
				let s = d[j] / (d[j] - d[i]);
				buf[out1 + m * 3 + 0] = buf[_in + j * 3 + 0] + (buf[_in + i * 3 + 0] - buf[_in + j * 3 + 0]) * s;
				buf[out1 + m * 3 + 1] = buf[_in + j * 3 + 1] + (buf[_in + i * 3 + 1] - buf[_in + j * 3 + 1]) * s;
				buf[out1 + m * 3 + 2] = buf[_in + j * 3 + 2] + (buf[_in + i * 3 + 2] - buf[_in + j * 3 + 2]) * s;
				RecastVectors.copy4(buf, out2 + n * 3, buf, out1 + m * 3);
				m++;
				n++;
				// add the i'th poPoly to the right polygon. Do NOT add points that are on the dividing line
				// since these were already added above
				if (d[i] > 0) {
					RecastVectors.copy4(buf, out1 + m * 3, buf, _in + i * 3);
					m++;
				} else if (d[i] < 0) {
					RecastVectors.copy4(buf, out2 + n * 3, buf, _in + i * 3);
					n++;
				}
			} else // same side
			{
				// add the i'th poPoly to the right polygon. Addition is done even for points on the dividing line
				if (d[i] >= 0) {
					RecastVectors.copy4(buf, out1 + m * 3, buf, _in + i * 3);
					m++;
					if (d[i] != 0)
						continue;
				}
				RecastVectors.copy4(buf, out2 + n * 3, buf, _in + i * 3);
				n++;
			}
		}
		return [m, n];
	}

	static rasterizeTri(verts, v0, v1, v2, area, hf, bmin,
		bmax, cs, ics, ich, flagMergeThr) {
		let w = hf.width;
		let h = hf.height;
		let tmin = new Array(3);
		let tmax = new Array(3);
		let by = bmax[1] - bmin[1];

		// Calculate the bounding box of the triangle.
		RecastVectors.copy3(tmin, verts, v0 * 3);
		RecastVectors.copy3(tmax, verts, v0 * 3);
		RecastVectors.min(tmin, verts, v1 * 3);
		RecastVectors.min(tmin, verts, v2 * 3);
		RecastVectors.max(tmax, verts, v1 * 3);
		RecastVectors.max(tmax, verts, v2 * 3);

		// If the triangle does not touch the bbox of the heightfield, skip the triagle.
		if (!RecastRasterization.overlapBounds(bmin, bmax, tmin, tmax))
			return;

		// Calculate the footprPoly of the triangle on the grid's y-axis
		let y0 = Math.floor((tmin[2] - bmin[2]) * ics);
		let y1 = Math.floor((tmax[2] - bmin[2]) * ics);
		y0 = RecastCommon.clamp(y0, 0, h - 1);
		y1 = RecastCommon.clamp(y1, 0, h - 1);

		// Clip the triangle into all grid cells it touches.
		let buf = new Array(7 * 3 * 4).fill(0);
		let _in = 0;
		let inrow = 7 * 3;
		let p1 = inrow + 7 * 3;
		let p2 = p1 + 7 * 3;

		RecastVectors.copy4(buf, 0, verts, v0 * 3);
		RecastVectors.copy4(buf, 3, verts, v1 * 3);
		RecastVectors.copy4(buf, 6, verts, v2 * 3);
		let nvrow, nvIn = 3;

		for (let y = y0; y <= y1; ++y) {
			// Clip polygon to row. Store the remaining polygon as well
			let cz = bmin[2] + y * cs;
			let nvrowin = RecastRasterization.dividePoly(buf, _in, nvIn, inrow, p1, cz + cs, 2);
			nvrow = nvrowin[0];
			nvIn = nvrowin[1];
			{
				let temp = _in;
				_in = p1;
				p1 = temp;
			}
			if (nvrow < 3)
				continue;

			// find the horizontal bounds _in the row
			let minX = buf[inrow];
			let maxX = buf[inrow];
			for (let i = 1; i < nvrow; ++i) {
				if (minX > buf[inrow + i * 3])
					minX = buf[inrow + i * 3];
				if (maxX < buf[inrow + i * 3])
					maxX = buf[inrow + i * 3];
			}
			let x0 = Math.floor((minX - bmin[0]) * ics);
			let x1 = Math.floor((maxX - bmin[0]) * ics);
			x0 = RecastCommon.clamp(x0, 0, w - 1);
			x1 = RecastCommon.clamp(x1, 0, w - 1);

			let nv, nv2 = nvrow;
			for (let x = x0; x <= x1; ++x) {
				// Clip polygon to column. store the remaining polygon as well
				let cx = bmin[0] + x * cs;
				let nvnv2 = RecastRasterization.dividePoly(buf, inrow, nv2, p1, p2, cx + cs, 0);
				nv = nvnv2[0];
				nv2 = nvnv2[1];
				{
					let temp = inrow;
					inrow = p2;
					p2 = temp;
				}
				if (nv < 3)
					continue;

				// Calculate min and max of the span.
				let smin = buf[p1 + 1];
				let smax = buf[p1 + 1];
				for (let i = 1; i < nv; ++i) {
					smin = Math.min(smin, buf[p1 + i * 3 + 1]);
					smax = Math.max(smax, buf[p1 + i * 3 + 1]);
				}
				smin -= bmin[1];
				smax -= bmin[1];
				// Skip the span if it is outside the heightfield bbox
				if (smax < 0.0)
					continue;
				if (smin > by)
					continue;
				// Clamp the span to the heightfield bbox.
				if (smin < 0.0)
					smin = 0;
				if (smax > by)
					smax = by;

				// Snap the span to the heightfield height grid.
				let ismin = RecastCommon.clamp(Math.floor(smin * ich), 0, RecastConstants.RC_SPAN_MAX_HEIGHT);
				let ismax = RecastCommon.clamp(Math.ceil(smax * ich), ismin + 1,
					RecastConstants.RC_SPAN_MAX_HEIGHT);

				RecastRasterization.addSpan(hf, x, y, ismin, ismax, area, flagMergeThr);
			}
		}
	}

	/**
	 * No spans will be added if the triangle does not overlap the heightfield grid.
	 * 
	 * @see Heightfield
	 */
	static rasterizeTriangle(ctx, verts, v0, v1, v2, area,
		solid, flagMergeThr) {

		ctx.startTimer("RASTERIZE_TRIANGLES");

		ics = 1.0 / solid.cs;
		ich = 1.0 / solid.ch;
		rasterizeTri(verts, v0, v1, v2, area, solid, solid.bmin, solid.bmax, solid.cs, ics, ich, flagMergeThr);

		ctx.stopTimer("RASTERIZE_TRIANGLES");
	}

	/**
	 * Spans will only be added for triangles that overlap the heightfield grid.
	 * 
	 * @see Heightfield
	 */
	static rasterizeTrianglesA( ctx, verts, tris, areas, nt,
		solid, flagMergeThr) {

		ctx.startTimer("RASTERIZE_TRIANGLES");

		let ics = 1.0 / solid.cs;
		let ich = 1.0 / solid.ch;
		// Rasterize triangles.
		for (let i = 0; i < nt; ++i) {
			let v0 = tris[i * 3 + 0];
			let v1 = tris[i * 3 + 1];
			let v2 = tris[i * 3 + 2];
			// Rasterize.
			RecastRasterization.rasterizeTri(verts, v0, v1, v2, areas[i], solid, solid.bmin, solid.bmax, solid.cs, ics, ich, flagMergeThr);
		}

		ctx.stopTimer("RASTERIZE_TRIANGLES");
	}

	/**
	 * Spans will only be added for triangles that overlap the heightfield grid.
	 * 
	 * @see Heightfield
	 */
	static rasterizeTrianglesB(ctx, verts, areas, nt, solid,
		flagMergeThr) {
		ctx.startTimer("RASTERIZE_TRIANGLES");

		let ics = 1.0 / solid.cs;
		let ich = 1.0 / solid.ch;
		// Rasterize triangles.
		for (let i = 0; i < nt; ++i) {
			let v0 = (i * 3 + 0);
			let v1 = (i * 3 + 1);
			let v2 = (i * 3 + 2);
			// Rasterize.
			RecastRasterization.rasterizeTri(verts, v0, v1, v2, areas[i], solid, solid.bmin, solid.bmax, solid.cs, ics, ich, flagMergeThr);
		}
		ctx.stopTimer("RASTERIZE_TRIANGLES");
	}
}

export default RecastRasterization;
