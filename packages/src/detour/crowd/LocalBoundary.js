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

import DetourCommon from "../DetourCommon.js"

function arraycopy(one, oneStart, two, twoStart, len) {
	for (let i = 0; i < len; i++) {
		two[twoStart + i] = one[oneStart + i];
	}
}

class LocalBoundary {

	static MAX_LOCAL_SEGS = 8;

	static Segment = class Segment {
		/** Segment start/end */
		s = new Array(6).fill(0);
		/** Distance for pruning. */
		d;
	}

	m_center = new Array(3);
	m_segs = [];
	m_polys = [];

	constructor() {
		this.m_center[0] = this.m_center[1] = this.m_center[2] = Number.MAX_VALUE;
	}

	reset() {
		this.m_center[0] = this.m_center[1] = this.m_center[2] = Number.MAX_VALUE;
		this.m_polys = [];
		this.m_segs = [];
	}

	addSegment(dist, s) {
		// Insert neighbour based on the distance.
		let seg = new LocalBoundary.Segment();
		arraycopy(s, 0, seg.s, 0, 6);
		seg.d = dist;
		if (this.m_segs.length == 0) {
			this.m_segs.push(seg);
		} else if (dist >= this.m_segs[this.m_segs.length - 1].d) {
			if (this.m_segs.length >= LocalBoundary.MAX_LOCAL_SEGS) {
				return;
			}
			this.m_segs.push(seg);
		} else {
			// Insert inbetween.
			let i;
			for (i = 0; i < this.m_segs.length; ++i)
				if (dist <= this.m_segs[i].d)
					break;
			this.m_segs.splice(i,0, seg);
		}
		while (this.m_segs.length > LocalBoundary.MAX_LOCAL_SEGS) {
			this.m_segs.splice(this.m_segs.length - 1,1);
		}
	}

	update(ref, pos, collisionQueryRange, navquery, filter) {
		if (ref == 0) {
			reset();
			return;
		}
		DetourCommon.vCopy(this.m_center, pos);
		// First query non-overlapping polygons.
		let res = navquery.findLocalNeighbourhood(ref, pos, collisionQueryRange, filter);
		this.m_polys = res.getRefs();
		this.m_segs = [];
		// Secondly, store all polygon edges.
		for(let j = 0; j < this.m_polys.length; ++j) {
			let gpws = navquery.getPolyWallSegments(this.m_polys[j], false, filter);
			for(let k = 0; k < gpws.getSegmentRefs().length; ++k) {
				let s = gpws.getSegmentVerts()[k];
				// Skip too distant segments.
				let distseg = DetourCommon.distancePtSegSqr2D4(pos, s, 0, 3);
				if (distseg[0] > DetourCommon.sqr(collisionQueryRange))
					continue;
				this.addSegment(distseg[0], s);
			}
		}
	}

	isValid(navquery, filter) {
		if (this.m_polys.length == 0)
			return false;

		// Check that all polygons still pass query filter.
		for (let ref of this.m_polys) {
			if (!navquery.isValidPolyRef(ref, filter))
				return false;
		}

		return true;
	}

	getCenter() {
		return this.m_center;
	}

	getSegment(j) {
		return this.m_segs[j].s;
	}

	getSegmentCount() {
		return this.m_segs.length;
	}
}

export default LocalBoundary;
