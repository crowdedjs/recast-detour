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

class RecastFilter {

	/// @par
	///
	/// Allows the formation of walkable regions that will flow over low lying 
	/// objects such as curbs, and up structures such as stairways. 
	/// 
	/// Two neighboring spans are walkable if: <tt>rcAbs(currentSpan.smax - neighborSpan.smax) < waklableClimb</tt>
	/// 
	/// @warning Will override the effect of #rcFilterLedgeSpans.  So if both filters are used, call
	/// #rcFilterLedgeSpans after calling this filter. 
	///
	/// @see rcHeightfield, rcConfig
static filterLowHangingWalkableObstacles(ctx,  walkableClimb, solid) {

		ctx.startTimer("FILTER_LOW_OBSTACLES");

		let w = solid.width;
		let h = solid.height;

		for (let y = 0; y < h; ++y) {
			for (let x = 0; x < w; ++x) {
				let ps = null;
				let previousWalkable = false;
				let previousArea = RecastConstants.RC_NULL_AREA;

				for (let s = solid.spans[x + y * w]; s != null; ps = s, s = s.next) {
					let walkable = s.area != RecastConstants.RC_NULL_AREA;
					// If current span is not walkable, but there is walkable
					// span just below it, mark the span above it walkable too.
					if (!walkable && previousWalkable) {
						if (Math.abs(s.smax - ps.smax) <= walkableClimb)
							s.area = previousArea;
					}
					// Copy walkable flag so that it cannot propagate
					// past multiple non-walkable objects.
					previousWalkable = walkable;
					previousArea = s.area;
				}
			}
		}

		ctx.stopTimer("FILTER_LOW_OBSTACLES");
	}

	/// @par
	///
	/// A ledge is a span with one or more neighbors whose maximum is further away than @p walkableClimb
	/// from the current span's maximum.
	/// This method removes the impact of the overestimation of conservative voxelization 
	/// so the resulting mesh will not have regions hanging in the air over ledges.
	/// 
	/// A span is a ledge if: <tt>rcAbs(currentSpan.smax - neighborSpan.smax) > walkableClimb</tt>
	/// 
	/// @see rcHeightfield, rcConfig
static filterLedgeSpans(ctx,  walkableHeight,  walkableClimb, solid) {
		ctx.startTimer("FILTER_BORDER");

		let w = solid.width;
		let h = solid.height;
		let MAX_HEIGHT = 0xfff;

		// Mark border spans.
		for (let y = 0; y < h; ++y) {
			for (let x = 0; x < w; ++x) {
				for (let s = solid.spans[x + y * w]; s != null; s = s.next) {
					// Skip non walkable spans.
					if (s.area == RecastConstants.RC_NULL_AREA)
						continue;

					let bot = (s.smax);
					let top = s.next != null ? s.next.smin : MAX_HEIGHT;

					// Find neighbours minimum height.
					let minh = MAX_HEIGHT;

					// Min and max height of accessible neighbours.
					let asmin = s.smax;
					let asmax = s.smax;

					for (let dir = 0; dir < 4; ++dir) {
						let dx = x + RecastCommon.GetDirOffsetX(dir);
						let dy = y + RecastCommon.GetDirOffsetY(dir);
						// Skip neighbours which are out of bounds.
						if (dx < 0 || dy < 0 || dx >= w || dy >= h) {
							minh = Math.min(minh, -walkableClimb - bot);
							continue;
						}

						// From minus infinity to the first span.
						let ns = solid.spans[dx + dy * w];
						let nbot = -walkableClimb;
						let ntop = ns != null ? ns.smin : MAX_HEIGHT;
						// Skip neightbour if the gap between the spans is too small.
						if (Math.min(top, ntop) - Math.max(bot, nbot) > walkableHeight)
							minh = Math.min(minh, nbot - bot);

						// Rest of the spans.
						for(let ns = solid.spans[dx + dy * w]; ns != null; ns = ns.next) {
							nbot = ns.smax;
							ntop = ns.next != null ? ns.next.smin : MAX_HEIGHT;
							// Skip neightbour if the gap between the spans is too small.
							if (Math.min(top, ntop) - Math.max(bot, nbot) > walkableHeight) {
								minh = Math.min(minh, nbot - bot);

								// Find min/max accessible neighbour height. 
								if (Math.abs(nbot - bot) <= walkableClimb) {
									if (nbot < asmin)
										asmin = nbot;
									if (nbot > asmax)
										asmax = nbot;
								}

							}
						}
					}

					// The current span is close to a ledge if the drop to any
					// neighbour span is less than the walkableClimb.
					if (minh < -walkableClimb)
						s.area = RecastConstants.RC_NULL_AREA;

					// If the difference between all neighbours is too large,
					// we are at steep slope, mark the span as ledge.
					if ((asmax - asmin) > walkableClimb) {
						s.area = RecastConstants.RC_NULL_AREA;
					}
				}
			}
		}

		ctx.stopTimer("FILTER_BORDER");
	}

	/// @par
	///
	/// For this filter, the clearance above the span is the distance from the span's 
	/// maximum to the next higher span's minimum. (Same grid column.)
	/// 
	/// @see rcHeightfield, rcConfig
static filterWalkableLowHeightSpans(ctx,  walkableHeight, solid) {
		ctx.startTimer("FILTER_WALKABLE");

		let w = solid.width;
		let h = solid.height;
		let MAX_HEIGHT = 0xfff;

		// Remove walkable flag from spans which do not have enough
		// space above them for the agent to stand there.
		for (let y = 0; y < h; ++y) {
			for (let x = 0; x < w; ++x) {
				for (let s = solid.spans[x + y * w]; s != null; s = s.next) {
					let bot = (s.smax);
					let top = s.next != null ? s.next.smin : MAX_HEIGHT;
					if ((top - bot) <= walkableHeight)
						s.area = RecastConstants.RC_NULL_AREA;
				}
			}
		}
		ctx.stopTimer("FILTER_WALKABLE");
	}
}

export default RecastFilter;
