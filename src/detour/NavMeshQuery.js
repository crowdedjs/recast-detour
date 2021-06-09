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

import NodePool from "./NodePool.js"
import NodeQueue from "./NodeQueue.js"
import DetourCommon from "./DetourCommon.js"
import FindNearestPolyResult from "./FindNearestPolyResult.js"
import NavMesh from "./NavMesh.js"
import Poly from "./Poly.js"
import ClosestPointOnPolyResult from "./ClosestPointOnPolyResult.js"
import QueryData from "./QueryData.js"
import Status from "./Status.js"
import Node from "./Node.js"
import UpdateSlicedPathResult from "./UpdateSlicedPathResult.js"
import FindPathResult from "./FindPathResult.js"
import FindLocalNeighbourhoodResult from "./FindLocalNeighbourhoodResult.js"
import GetPolyWallSegmentsResult from "./GetPolyWallSegmentsResult.js"
import StraightPathItem from "./StraightPathItem.js"
import MoveAlongSurfaceResult from "./MoveAlongSurfaceResult.js"
import RaycastHit from "./RaycastHit.js"
import IntersectResult from "./IntersectResult.js"


function or(v1, v2) {
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

class PortalResult {
	left;
	right;
	fromType;
	toType;

	constructor(left, right, fromType, toType) {
		this.left = left;
		this.right = right;
		this.fromType = fromType;
		this.toType = toType;
	}

}


function arraycopy(one, oneStart, two, twoStart, len) {
	for (let i = 0; i < len; i++) {
		two[twoStart + i] = one[oneStart + i];
	}
}




class NavMeshQuery {

	/**
	 * Use raycasts during pathfind to "shortcut" (raycast still consider costs)
	 * Options for NavMeshQuery::initSlicedFindPath and updateSlicedFindPath
	 */
	static DT_FINDPATH_ANY_ANGLE = 0x02;

	/** Raycast should calculate movement cost aPoly the ray and fill RaycastHit::cost */
	static DT_RAYCAST_USE_COSTS = 0x01;

	/// Vertex flags returned by findStraightPath.
	/** The vertex is the start position in the path. */
	static DT_STRAIGHTPATH_START = 0x01;
	/** The vertex is the end position in the path. */
	static DT_STRAIGHTPATH_END = 0x02;
	/** The vertex is the start of an off-mesh connection. */
	static DT_STRAIGHTPATH_OFFMESH_CONNECTION = 0x04;

	/// Options for findStraightPath.
	static DT_STRAIGHTPATH_AREA_CROSSINGS = 0x01; ///< Add a vertex at every polygon edge crossing where area changes.
	static DT_STRAIGHTPATH_ALL_CROSSINGS = 0x02; ///< Add a vertex at every polygon edge crossing.

	static H_SCALE = 0.999; // Search heuristic scale.

	m_nav;
	m_nodePool;
	m_tinyNodePool;
	m_openList;
	m_query; /// < Sliced query state.

	constructor(nav) {
		this.m_nav = nav;
		this.m_nodePool = new NodePool();
		this.m_tinyNodePool = new NodePool();
		this.m_openList = new NodeQueue();
	}

	static FRand = class FRand {
		constructor() {
			return Math.random();
		}
	}

	/**
	 * Returns random location on navmesh.
	 * Polygons are chosen weighted by area. The search runs in linear related to number of polygon.
	 * @param filter The polygon filter to apply to the query.
	 * @param frand Function returning a random number [0..1).
	 * @return Random location
	 */
	findRandomPoint(filter, frand) {
		// Randomly pick one tile. Assume that all tiles cover roughly the same area.
		tile = null;
		tsum = 0.0;
		for (let i = 0; i < this.m_nav.getMaxTiles(); i++) {
			let t = this.m_nav.getTile(i);
			if (t == null || t.data == null || t.data.header == null)
				continue;

			// Choose random tile using reservoi sampling.
			area = 1.0; // Could be tile area too.
			tsum += area;
			u = frand.frand();
			if (u * tsum <= area)
				tile = t;
		}
		if (tile == null)
			return new FindRandomPointResult(Status.FAILURE, 0, null);

		// Randomly pick one polygon weighted by polygon area.
		let poly = null;
		let polyRef = 0;
		let base = this.m_nav.getPolyRefBase(tile);

		areaSum = 0.0;
		for (let i = 0; i < tile.data.header.polyCount; ++i) {
			let p = tile.data.polys[i];
			// Do not return off-mesh connection polygons.
			if (p.getType() != Poly.DT_POLYTYPE_GROUND)
				continue;
			// Must pass filter
			let ref = this.or(base, i);
			if (!filter.passFilter(ref, tile, p))
				continue;

			// Calc area of the polygon.
			polyArea = 0.0;
			for (let j = 2; j < p.vertCount; ++j) {
				let va = p.verts[0] * 3;
				let vb = p.verts[j - 1] * 3;
				let vc = p.verts[j] * 3;
				polyArea += DetourCommon.triArea2D4(tile.data.verts, va, vb, vc);
			}

			// Choose random polygon weighted by area, using reservoi sampling.
			areaSum += polyArea;
			u = frand.frand();
			if (u * areaSum <= polyArea) {
				poly = p;
				polyRef = ref;
			}
		}

		if (poly == null)
			return new FindRandomPointResult(Status.FAILURE, 0, null);

		// Randomly pick poPoly on polygon.
		let verts = new Array(3 * this.m_nav.getMaxVertsPerPoly());
		let areas = new Array(this.m_nav.getMaxVertsPerPoly());
		arraycopy(tile.data.verts, poly.verts[0] * 3, verts, 0, 3);
		for (let j = 1; j < poly.vertCount; ++j) {
			arraycopy(tile.data.verts, poly.verts[j] * 3, verts, j * 3, 3);
		}

		s = frand.frand();
		t = frand.frand();

		let pt = DetourCommon.randomPointInConvexPoly(verts, poly.vertCount, areas, s, t);

		pt[1] = getPolyHeight(polyRef, pt);

		return new FindRandomPointResult(Status.SUCCSESS, polyRef, pt);
	}

	/**
	 * Returns random location on navmesh within the reach of specified location.
	 * Polygons are chosen weighted by area. The search runs in linear related to number of polygon.
	 * The location is not exactly constrained by the circle, but it limits the visited polygons.
	 * 
	 * @param startRef The reference id of the polygon where the search starts.
	 * @param centerPos The center of the search circle. [(x, y, z)]
	 * @param maxRadius 
	 * @param filter The polygon filter to apply to the query.
	 * @param frand Function returning a random number [0..1).
	 * @return Random location
	 */
	findRandomPointAroundCircle(startRef, centerPos, maxRadius,
		filter, frand) {

		// Validate input
		if (startRef == 0 || !this.m_nav.isValidPolyRef(startRef))
			throw new IllegalArgumentException("Invalid start ref");

		let tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(startRef);
		let startTile = tileAndPoly[0];
		let startPoly = tileAndPoly[1];
		if (!filter.passFilter(startRef, startTile, startPoly))
			throw new IllegalArgumentException("Invalid start");

		m_nodePool = [];
		this.m_openList = [];

		let startNode = m_nodePool.getNode(startRef);
		DetourCommon.vCopy(startNode.pos, centerPos);
		startNode.pidx = 0;
		startNode.cost = 0;
		startNode.total = 0;
		startNode.id = startRef;
		startNode.flags = DT_NODE_OPEN;
		this.m_openList.push(startNode);

		radiusSqr = maxRadius * maxRadius;
		areaSum = 0.0;

		let randomTile = null;
		let randomPoly = null;
		let randomPolyRef = 0;

		while (!this.m_openList.length == 0) {
			let bestNode = this.m_openList.pop();
			bestNode.flags &= ~DT_NODE_OPEN;
			bestNode.flags |= DT_NODE_CLOSED;
			// Get let and tile.
			// The API input has been cheked already, skip checking internal data.
			let bestRef = bestNode.id;
			let bestTilePoly = this.m_nav.getTileAndPolyByRefUnsafe(bestRef);
			let bestTile = bestTilePoly[0];
			let bestPoly = bestTilePoly[1];

			// Place random locations on on ground.
			if (bestPoly.getType() == Poly.DT_POLYTYPE_GROUND) {
				// Calc area of the polygon.
				polyArea = 0.0;
				for (let j = 2; j < bestPoly.vertCount; ++j) {
					let va = bestPoly.verts[0] * 3;
					let vb = bestPoly.verts[j - 1] * 3;
					let vc = bestPoly.verts[j] * 3;
					polyArea += DetourCommon.triArea2D4(bestTile.data.verts, va, vb, vc);
				}
				// Choose random polygon weighted by area, using reservoi sampling.
				areaSum += polyArea;
				u = frand.frand();
				if (u * areaSum <= polyArea) {
					randomTile = bestTile;
					randomPoly = bestPoly;
					randomPolyRef = bestRef;
				}
			}

			// Get parent let and tile.
			let parentRef = 0;
			if (bestNode.pidx != 0)
				parentRef = m_nodePool.getNodeAtIdx(bestNode.pidx).id;
			if (parentRef != 0) {
				let parentTilePoly = this.m_nav.getTileAndPolyByRefUnsafe(parentRef);
				let parentTile = parentTilePoly[0];
				let parentPoly = parentTilePoly[1];
			}

			for (let i = bestPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = bestTile.links[i].next) {
				let link = bestTile.links[i];
				let neighbourRef = link.ref;
				// Skip invalid neighbours and do not follow back to parent.
				if (neighbourRef == 0 || neighbourRef == parentRef)
					continue;

				// Expand to neighbour
				let neighbourTilePoly = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
				let neighbourTile = neighbourTilePoly[0];
				let neighbourPoly = neighbourTilePoly[1];

				// Do not advance if the polygon is excluded by the filter.
				if (!filter.passFilter(neighbourRef, neighbourTile, neighbourPoly))
					continue;

				// Find edge and calc distance to the edge.
				let portalpoints = this.getPortalPoints7(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly,
					neighbourTile, 0, 0);
				let va = portalpoints.left;
				let vb = portalpoints.right;

				// If the circle is not touching the next polygon, skip it.
				let distseg = DetourCommon.distancePtSegSqr2D3(centerPos, va, vb);
				distSqr = distseg[0];
				if (distSqr > radiusSqr)
					continue;

				let neighbourNode = m_nodePool.getNode(neighbourRef);

				if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0)
					continue;

				// Cost
				if (neighbourNode.flags == 0)
					neighbourNode.pos = DetourCommon.vLerp3(va, vb, 0.5);

				total = bestNode.total + DetourCommon.vDist2(bestNode.pos, neighbourNode.pos);

				// The node is already in open list and the new result is worse, skip.
				if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0 && total >= neighbourNode.total)
					continue;

				neighbourNode.id = neighbourRef;
				neighbourNode.flags = (neighbourNode.flags & ~Node.DT_NODE_CLOSED);
				neighbourNode.pidx = m_nodePool.getNodeIdx(bestNode);
				neighbourNode.total = total;

				if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0) {
					this.m_openList.modify(neighbourNode);
				} else {
					neighbourNode.flags = Node.DT_NODE_OPEN;
					this.m_openList.push(neighbourNode);
				}
			}
		}

		if (randomPoly == null)
			return new FindRandomPointResult(Status.FAILURE, 0, null);

		// Randomly pick poPoly on polygon.
		let verts = new Array(3 * this.m_nav.getMaxVertsPerPoly());
		let areas = new Array(this.m_nav.getMaxVertsPerPoly());
		arraycopy(randomTile.data.verts, randomPoly.verts[0] * 3, verts, 0, 3);
		for (let j = 1; j < randomPoly.vertCount; ++j) {
			arraycopy(randomTile.data.verts, randomPoly.verts[j] * 3, verts, j * 3, 3);
		}

		s = frand.frand();
		t = frand.frand();

		let pt = DetourCommon.randomPointInConvexPoly(verts, randomPoly.vertCount, areas, s, t);

		pt[1] = getPolyHeight(randomPolyRef, pt);

		return new FindRandomPointResult(Status.SUCCSESS, randomPolyRef, pt);
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	/// @par
	///
	/// Uses the detail polygons to find the surface height. (Most accurate.)
	///
	/// @p pos does not have to be within the bounds of the polygon or navigation mesh.
	///
	/// See closestPointOnPolyBoundary() for a limited but faster option.
	///
	/// Finds the closest poPoly on the specified polygon.
	///  @param[in]		ref			The reference id of the polygon.
	///  @param[in]		pos			The position to check. [(x, y, z)]
	///  @param[out]	closest		
	///  @param[out]	posOverPoly	
	/// @returns The status flags for the query.
	closestPointOnPoly(ref, pos) {
		let tileAndPoly = this.m_nav.getTileAndPolyByRef(ref);
		let tile = tileAndPoly[0];
		let poly = tileAndPoly[1];

		// Off-mesh connections don't have detail polygons.
		if (poly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) {
			let v0 = poly.verts[0] * 3;
			let v1 = poly.verts[1] * 3;
			let d0 = DetourCommon.vDist3(pos, tile.data.verts, v0);
			let d1 = DetourCommon.vDist3(pos, tile.data.verts, v1);
			let u = d0 / (d0 + d1);
			let closest = DetourCommon.vLerp4(tile.data.verts, v0, v1, u);
			return new ClosestPointOnPolyResult(false, closest);
		}

		// Clamp poPoly to be inside the polygon.
		let verts = new Array(this.m_nav.getMaxVertsPerPoly() * 3);
		let edged = new Array(this.m_nav.getMaxVertsPerPoly());
		let edget = new Array(this.m_nav.getMaxVertsPerPoly());
		let nv = poly.vertCount;
		for (let i = 0; i < nv; ++i)
			arraycopy(tile.data.verts, poly.verts[i] * 3, verts, i * 3, 3);

		let posOverPoly = false;
		let closest = new Array(3);
		DetourCommon.vCopy(closest, pos);
		if (!DetourCommon.distancePtPolyEdgesSqr(pos, verts, nv, edged, edget)) {
			// PoPoly is outside the polygon, dtClamp to nearest edge.
			let dmin = edged[0];
			let imin = 0;
			for (let i = 1; i < nv; ++i) {
				if (edged[i] < dmin) {
					dmin = edged[i];
					imin = i;
				}
			}
			let va = imin * 3;
			let vb = ((imin + 1) % nv) * 3;
			closest = DetourCommon.vLerp4(verts, va, vb, edget[imin]);
			posOverPoly = false;
		} else {
			posOverPoly = true;
		}
		let ip = poly.index;
		if (tile.data.detailMeshes != null && tile.data.detailMeshes.length > ip) {
			let pd = tile.data.detailMeshes[ip];
			// Find height at the location.
			for (let j = 0; j < pd.triCount; ++j) {
				let t = (pd.triBase + j) * 4;
				let v = new Array(3); //Was new Array(3)[]
				for (let k = 0; k < 3; ++k) {
					if (tile.data.detailTris[t + k] < poly.vertCount) {
						let index = poly.verts[tile.data.detailTris[t + k]] * 3;
						v[k] = [tile.data.verts[index], tile.data.verts[index + 1],
						tile.data.verts[index + 2]];
					} else {
						let index = (pd.vertBase + (tile.data.detailTris[t + k] - poly.vertCount)) * 3;
						v[k] = [tile.data.detailVerts[index], tile.data.detailVerts[index + 1],
						tile.data.detailVerts[index + 2]];
					}
				}
				let heightResult = DetourCommon.closestHeightPointTriangle(closest, v[0], v[1], v[2]);
				if (heightResult[0]) {
					closest[1] = heightResult[1];
					break;
				}
			}
		}
		return new ClosestPointOnPolyResult(posOverPoly, closest);
	}

	/// @par
	///
	/// Much faster than closestPointOnPoly().
	///
	/// If the provided position lies within the polygon's xz-bounds (above or below),
	/// then @p pos and @p closest will be equal.
	///
	/// The height of @p closest will be the polygon boundary. The height detail is not used.
	///
	/// @p pos does not have to be within the bounds of the polybon or the navigation mesh.
	///
	/// Returns a poPoly on the boundary closest to the source poPoly if the source poPoly is outside the 
	/// polygon's xz-bounds.
	///  @param[in]		ref			The reference id to the polygon.
	///  @param[in]		pos			The position to check. [(x, y, z)]
	///  @param[out]	closest		The closest point. [(x, y, z)]
	/// @returns The status flags for the query.
	closestPointOnPolyBoundary(ref, pos) {

		let tileAndPoly = this.m_nav.getTileAndPolyByRef(ref);
		let tile = tileAndPoly[0];
		let poly = tileAndPoly[1];

		// Collect vertices.
		let verts = new Array(this.m_nav.getMaxVertsPerPoly() * 3);
		let edged = new Array(this.m_nav.getMaxVertsPerPoly());
		let edget = new Array(this.m_nav.getMaxVertsPerPoly());
		let nv = poly.vertCount;
		for (let i = 0; i < nv; ++i)
			arraycopy(tile.data.verts, poly.verts[i] * 3, verts, i * 3, 3);

		let closest;
		if (DetourCommon.distancePtPolyEdgesSqr(pos, verts, nv, edged, edget)) {
			closest = DetourCommon.vCopy_return(pos);
		} else {
			// PoPoly is outside the polygon, dtClamp to nearest edge.
			let dmin = edged[0];
			let imin = 0;
			for (let i = 1; i < nv; ++i) {
				if (edged[i] < dmin) {
					dmin = edged[i];
					imin = i;
				}
			}
			let va = imin * 3;
			let vb = ((imin + 1) % nv) * 3;
			closest = DetourCommon.vLerp4(verts, va, vb, edget[imin]);
		}
		return closest;
	}

	/// @par
	///
	/// Will return #DT_FAILURE if the provided position is outside the xz-bounds
	/// of the polygon.
	///
	/// Gets the height of the polygon at the provided position using the height detail. (Most accurate.)
	///  @param[in]		ref			The reference id of the polygon.
	///  @param[in]		pos			A position within the xz-bounds of the polygon. [(x, y, z)]
	///  @param[out]	height		The height at the surface of the polygon.
	/// @returns The status flags for the query.
	getPolyHeight(ref, pos) {
		let tileAndPoly = this.m_nav.getTileAndPolyByRef(ref);
		let tile = tileAndPoly[0];
		let poly = tileAndPoly[1];
		if (poly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) {
			let i = poly.verts[0] * 3;
			i = poly.verts[1] * 3;
			letv1 = [tile.data.verts[i], tile.data.verts[i + 1], tile.data.verts[i + 2]];
			d0 = DetourCommon.vDist2D(pos, v0);
			d1 = DetourCommon.vDist2D(pos, v1);
			u = d0 / (d0 + d1);
			return v0[1] + (v1[1] - v0[1]) * u;
		} else {
			let ip = poly.index;
			let pd = tile.data.detailMeshes[ip];
			for (let j = 0; j < pd.triCount; ++j) {
				let t = (pd.triBase + j) * 4;
				let v = new Array(3);//new Array(3)[];
				for (let k = 0; k < 3; ++k) {
					if (tile.data.detailTris[t + k] < poly.vertCount) {
						let index = poly.verts[tile.data.detailTris[t + k]] * 3;
						v[k] = [tile.data.verts[index], tile.data.verts[index + 1],
						tile.data.verts[index + 2]];
					} else {
						let index = (pd.vertBase + (tile.data.detailTris[t + k] - poly.vertCount)) * 3;
						v[k] = [tile.data.detailVerts[index], tile.data.detailVerts[index + 1],
						tile.data.detailVerts[index + 2]];
					}
				}
				let heightResult = DetourCommon.closestHeightPointTriangle(pos, v[0], v[1], v[2]);
				if (heightResult[0]) {
					return heightResult[1];
				}
			}
		}
		throw new IllegalArgumentException("Invalid ref " + ref + " pos " + Arrays.toString(pos));
	}

	/// @par
	///
	/// @note If the search box does not intersect any polygons the search will
	/// return #DT_SUCCESS, but @p nearestRef will be zero. So if in doubt, check
	/// @p nearestRef before using @p nearestPt.
	///
	/// @}
	/// @name Local Query Functions
	///@{

	/// Finds the polygon nearest to the specified center point.
	///  @param[in]		center		The center of the search box. [(x, y, z)]
	///  @param[in]		extents		The search distance aPoly each axis. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	/// @returns The status flags for the query.
	findNearestPoly(center, extents, filter) {

		let nearestPt = null;

		// Get nearby polygons from proximity grid.
		let polys = this.queryPolygons(center, extents, filter);

		// Find nearest polygon amongst the nearby polygons.
		let nearest = 0;
		let nearestDistanceSqr = Number.MAX_VALUE;
		for (let i = 0; i < polys.length; ++i) {
			let ref = polys[i];
			let closest = this.closestPointOnPoly(ref, center);
			let posOverPoly = closest.isPosOverPoly();
			let closestPtPoly = closest.getClosest();

			// If a poPoly is directly over a polygon and closer than
			// climb height, favor that instead of straight line nearest point.
			let d = 0;
			let diff = DetourCommon.vSub(center, closestPtPoly);
			if (posOverPoly) {
				let tilaAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(polys[i]);
				let tile = tilaAndPoly[0];
				d = Math.abs(diff[1]) - tile.data.header.walkableClimb;
				d = d > 0 ? d * d : 0;
			} else {
				d = DetourCommon.vLenSqr(diff);
			}

			if (d < nearestDistanceSqr) {
				nearestPt = closestPtPoly;
				nearestDistanceSqr = d;
				nearest = ref;
			}
		}

		return new FindNearestPolyResult(nearest, nearestPt);
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


	// FIXME: (PP) duplicate?
	queryPolygonsInTile(tile, qmin, qmax, filter) {
		let polys = [];
		if (tile.data.bvTree != null) {
			let nodeIndex = 0;
			let tbmin = tile.data.header.bmin;
			let tbmax = tile.data.header.bmax;
			let qfac = tile.data.header.bvQuantFactor;
			// Calculate quantized box
			let bmin = new Array(3);
			let bmax = new Array(3);
			// dtClamp query box to world box.
			let minx = DetourCommon.clamp(qmin[0], tbmin[0], tbmax[0]) - tbmin[0];
			let miny = DetourCommon.clamp(qmin[1], tbmin[1], tbmax[1]) - tbmin[1];
			let minz = DetourCommon.clamp(qmin[2], tbmin[2], tbmax[2]) - tbmin[2];
			let maxx = DetourCommon.clamp(qmax[0], tbmin[0], tbmax[0]) - tbmin[0];
			let maxy = DetourCommon.clamp(qmax[1], tbmin[1], tbmax[1]) - tbmin[1];
			let maxz = DetourCommon.clamp(qmax[2], tbmin[2], tbmax[2]) - tbmin[2];
			// Quantize
			bmin[0] = this.and(Math.floor((qfac * minx)), 0xfffe);
			bmin[1] = this.and(Math.floor((qfac * miny)), 0xfffe);
			bmin[2] = this.and(Math.floor((qfac * minz)), 0xfffe);
			bmax[0] = this.or(Math.floor((qfac * maxx + 1)), 1);
			bmax[1] = this.or(Math.floor((qfac * maxy + 1)), 1);
			bmax[2] = this.or(Math.floor((qfac * maxz + 1)), 1);

			// Traverse tree
			let base = this.m_nav.getPolyRefBase(tile);
			let end = tile.data.header.bvNodeCount;
			while (nodeIndex < end) {
				let node = tile.data.bvTree[nodeIndex];
				let overlap = DetourCommon.overlapQuantBounds(bmin, bmax, node.bmin, node.bmax);
				let isLeafNode = node.i >= 0;

				if (isLeafNode && overlap) {
					let ref = this.or(base, node.i);
					if (filter.passFilter(ref, tile, tile.data.polys[node.i])) {
						polys.push(ref);
					}
				}

				if (overlap || isLeafNode)
					nodeIndex++;
				else {
					let escapeIndex = -node.i;
					nodeIndex += escapeIndex;
				}
			}
			return polys;
		} else {
			let bmin = new Array(3);
			let bmax = new Array(3);
			let base = this.m_nav.getPolyRefBase(tile);
			for (let i = 0; i < tile.data.header.polyCount; ++i) {
				let p = tile.data.polys[i];
				// Do not return off-mesh connection polygons.
				if (p.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION)
					continue;
				let ref = base | i;
				if (!filter.passFilter(ref, tile, p))
					continue;
				// Calc polygon bounds.
				let v = p.verts[0] * 3;
				DetourCommon.vCopy(bmin, tile.data.verts, v);
				DetourCommon.vCopy(bmax, tile.data.verts, v);
				for (let j = 1; j < p.vertCount; ++j) {
					v = p.verts[j] * 3;
					DetourCommon.vMin(bmin, tile.data.verts, v);
					DetourCommon.vMax(bmax, tile.data.verts, v);
				}
				if (overlapBounds(qmin, qmax, bmin, bmax)) {
					polys.push(ref);
				}
			}
			return polys;
		}
	}

	/**
	 * Finds polygons that overlap the search box.
	 * 
	 * If no polygons are found, the function will return with a polyCount of zero.
	 * 
	 * @param center
	 *            The center of the search box. [(x, y, z)]
	 * @param extents
	 *            The search distance aPoly each axis. [(x, y, z)]
	 * @param (filter)
	 *            The polygon filter to apply to the query.
	 * @return The reference ids of the polygons that overlap the query box.
	 */
	queryPolygons(center, extents, filter) {
		let bmin = DetourCommon.vSub(center, extents);
		let bmax = DetourCommon.vAdd(center, extents);
		// Find tiles the query touches.
		let minxy = this.m_nav.calcTileLoc(bmin);
		let minx = minxy[0];
		let miny = minxy[1];
		let maxxy = this.m_nav.calcTileLoc(bmax);
		let maxx = maxxy[0];
		let maxy = maxxy[1];
		let polys = [];
		for (let y = miny; y <= maxy; ++y) {
			for (let x = minx; x <= maxx; ++x) {
				let neis = this.m_nav.getTilesAt(x, y);
				for (let j = 0; j < neis.length; ++j) {
					let polysInTile = this.queryPolygonsInTile(neis[j], bmin, bmax, filter);
					polys.push(...polysInTile);
				}
			}
		}
		return polys;
	}

	/**
	 * Finds a path from the start polygon to the end polygon.
	 * 
	 * If the end polygon cannot be reached through the navigation graph, the last polygon in the path will be the
	 * nearest the end polygon.
	 * 
	 * The start and end positions are used to calculate traversal costs. (The y-values impact the result.)
	 * 
	 * @param startRef
	 *            The refrence id of the start polygon.
	 * @param endRef
	 *            The reference id of the end polygon.
	 * @param startPos
	 *            A position within the start polygon. [(x, y, z)]
	 * @param endPos
	 *            A position within the end polygon. [(x, y, z)]
	 * @param filter
	 *            The polygon filter to apply to the query.
	 * @return Found path
	 */
	findPath(startRef, endRef, startPos, endPos,
		filter) {
		if (startRef == 0 || endRef == 0)
			throw new IllegalArgumentException("Start or end ref = 0");

		// Validate input
		if (!this.m_nav.isValidPolyRef(startRef) || !this.m_nav.isValidPolyRef(endRef))
			throw new IllegalArgumentException("Invalid start or end ref");

		if (startRef == endRef) {
			let path = new Array(1);
			path.push(startRef);
			return new FindPathResult(Status.SUCCSESS, path);
		}

		this.m_nodePool = [];
		this.m_openList = [];

		let startNode = this.m_nodePool.getNode(startRef);
		DetourCommon.vCopy(startNode.pos, startPos);
		startNode.pidx = 0;
		startNode.cost = 0;
		startNode.total = DetourCommon.vDist2(startPos, endPos) * NavMeshQuery.H_SCALE;
		startNode.id = startRef;
		startNode.flags = Node.DT_NODE_OPEN;
		this.m_openList.push(startNode);

		let lastBestNode = startNode;
		lastBestNodeCost = startNode.total;

		let status = Status.SUCCSESS;

		while (!this.m_openList.length == 0) {
			// Remove node from open list and put it in closed list.
			let bestNode = this.m_openList.pop();
			bestNode.flags &= ~Node.DT_NODE_OPEN;
			bestNode.flags |= Node.DT_NODE_CLOSED;

			// Reached the goal, stop searching.
			if (bestNode.id == endRef) {
				lastBestNode = bestNode;
				break;
			}

			// Get current let and tile.
			// The API input has been cheked already, skip checking internal data.
			let bestRef = bestNode.id;
			let tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(bestRef);
			let bestTile = tileAndPoly[0];
			let bestPoly = tileAndPoly[1];

			// Get parent let and tile.
			let parentRef = 0;
			let parentTile = null;
			let parentPoly = null;
			if (bestNode.pidx != 0)
				parentRef = this.m_nodePool.getNodeAtIdx(bestNode.pidx).id;
			if (parentRef != 0) {
				tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(parentRef);
				parentTile = tileAndPoly[0];
				parentPoly = tileAndPoly[1];
			}

			for (let i = bestPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = bestTile.links[i].next) {
				let neighbourRef = bestTile.links[i].ref;

				// Skip invalid ids and do not expand back to where we came from.
				if (neighbourRef == 0 || neighbourRef == parentRef)
					continue;

				// Get neighbour let and tile.
				// The API input has been cheked already, skip checking internal data.
				tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
				let neighbourTile = tileAndPoly[0];
				let neighbourPoly = tileAndPoly[1];

				if (!filter.passFilter(neighbourRef, neighbourTile, neighbourPoly))
					continue;

				// deal explicitly with crossing tile boundaries
				let crossSide = 0;
				if (bestTile.links[i].side != 0xff)
					crossSide = bestTile.links[i].side >> 1;

				// get the node
				let neighbourNode = this.m_nodePool.getNode(neighbourRef, crossSide);

				// If the node is visited the first time, calculate node position.
				if (neighbourNode.flags == 0) {
					neighbourNode.pos = this.getEdgeMidPoint6(bestRef, bestPoly, bestTile, neighbourRef,
						neighbourPoly, neighbourTile);
				}

				// Calculate cost and heuristic.
				cost = 0;
				heuristic = 0;

				// Special case for last node.
				if (neighbourRef == endRef) {
					// Cost
					curCost = filter.getCost(bestNode.pos, neighbourNode.pos, parentRef, parentTile, parentPoly,
						bestRef, bestTile, bestPoly, neighbourRef, neighbourTile, neighbourPoly);
					let endCost = filter.getCost(neighbourNode.pos, endPos, bestRef, bestTile, bestPoly, neighbourRef,
						neighbourTile, neighbourPoly, 0, null, null);

					cost = bestNode.cost + curCost + endCost;
					heuristic = 0;
				} else {
					// Cost
					curCost = filter.getCost(bestNode.pos, neighbourNode.pos, parentRef, parentTile, parentPoly,
						bestRef, bestTile, bestPoly, neighbourRef, neighbourTile, neighbourPoly);
					cost = bestNode.cost + curCost;
					heuristic = DetourCommon.vDist2(neighbourNode.pos, endPos) * NavMeshQuery.H_SCALE;
				}

				total = cost + heuristic;

				// The node is already in open list and the new result is worse, skip.
				if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0 && total >= neighbourNode.total)
					continue;
				// The node is already visited and process, and the new result is worse, skip.
				if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0 && total >= neighbourNode.total)
					continue;

				// Add or update the node.
				neighbourNode.pidx = this.m_nodePool.getNodeIdx(bestNode);
				neighbourNode.id = neighbourRef;
				neighbourNode.flags = (neighbourNode.flags & ~Node.DT_NODE_CLOSED);
				neighbourNode.cost = cost;
				neighbourNode.total = total;

				if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0) {
					// Already in open, update node location.
					this.m_openList.modify(neighbourNode);
				} else {
					// Put the node in open list.
					neighbourNode.flags |= Node.DT_NODE_OPEN;
					this.m_openList.push(neighbourNode);
				}

				// Update nearest node to target so far.
				if (heuristic < lastBestNodeCost) {
					lastBestNodeCost = heuristic;
					lastBestNode = neighbourNode;
				}
			}
		}

		let path = getPathToNode(lastBestNode);

		if (lastBestNode.id != endRef)
			status = Status.PARTIAL_RESULT;

		return new FindPathResult(status, path);
	}

	/**
	 * Intializes a sliced path query.
	 * 
	 * Common use case: -# Call initSlicedFindPath() to initialize the sliced path query. -# Call updateSlicedFindPath()
	 * until it returns compPolye. -# Call finalizeSlicedFindPath() to get the path.
	 * 
	 * @param startRef
	 *            The reference id of the start polygon.
	 * @param endRef
	 *            The reference id of the end polygon.
	 * @param startPos
	 *            A position within the start polygon. [(x, y, z)]
	 * @param endPos
	 *            A position within the end polygon. [(x, y, z)]
	 * @param filter
	 *            The polygon filter to apply to the query.
	 * @param options
	 *            query options (see: #FindPathOptions)
	 * @return
	 */
	initSlicedFindPath(startRef, endRef, startPos, endPos, filter,
		options) {
		// Init path state.
		this.m_query = new QueryData();
		this.m_query.status = Status.FAILURE;
		this.m_query.startRef = startRef;
		this.m_query.endRef = endRef;
		DetourCommon.vCopy(this.m_query.startPos, startPos);
		DetourCommon.vCopy(this.m_query.endPos, endPos);
		this.m_query.filter = filter;
		this.m_query.options = options;
		this.m_query.raycastLimitSqr = Number.MAX_VALUE;

		if (startRef == 0 || endRef == 0)
			throw new IllegalArgumentException("Start or end ref = 0");

		// Validate input
		if (!this.m_nav.isValidPolyRef(startRef) || !this.m_nav.isValidPolyRef(endRef))
			throw new IllegalArgumentException("Invalid start or end ref");

		// trade quality with performance?
		if ((options & NavMeshQuery.DT_FINDPATH_ANY_ANGLE) != 0) {
			// limiting to several times the character radius yields nice results. It is not sensitive
			// so it is enough to compute it from the first tile.
			let tile = this.m_nav.getTileByRef(startRef);
			agentRadius = tile.data.header.walkableRadius;
			this.m_query.raycastLimitSqr = DetourCommon.sqr(agentRadius * NavMesh.DT_RAY_CAST_LIMIT_PROPORTIONS);
		}

		if (startRef == endRef) {
			this.m_query.status = Status.SUCCSESS;
			return Status.SUCCSESS;
		}

		this.m_nodePool.clear();
		this.m_openList.clear();

		let startNode = this.m_nodePool.getNode(startRef);
		DetourCommon.vCopy(startNode.pos, startPos);
		startNode.pidx = 0;
		startNode.cost = 0;
		startNode.total = DetourCommon.vDist2(startPos, endPos) * NavMeshQuery.H_SCALE;
		startNode.id = startRef;
		startNode.flags = Node.DT_NODE_OPEN;
		this.m_openList.push(startNode);

		this.m_query.status = Status.IN_PROGRESS;
		this.m_query.lastBestNode = startNode;
		this.m_query.lastBestNodeCost = startNode.total;

		return this.m_query.status;
	}

	/**
	 * Updates an in-progress sliced path query.
	 * 
	 * @param maxIter
	 *            The maximum number of iterations to perform.
	 * @return The status flags for the query.
	 */
	updateSlicedFindPath(maxIter) {
		if (this.m_query.status  != Status.IN_PROGRESS)
			return new UpdateSlicedPathResult(this.m_query.status, 0);

		// Make sure the request is still valid.
		if (!this.m_nav.isValidPolyRef(this.m_query.startRef) || !this.m_nav.isValidPolyRef(this.m_query.endRef)) {
			this.m_query.status = Status.FAILURE;
			return new UpdateSlicedPathResult(this.m_query.status, 0);
		}

		let iter = 0;
		while (iter < maxIter && !this.m_openList.isEmpty()) {
			iter++;

			// Remove node from open list and put it in closed list.
			let bestNode = this.m_openList.pop();
			bestNode.flags &= ~Node.DT_NODE_OPEN;
			bestNode.flags |= Node.DT_NODE_CLOSED;

			// Reached the goal, stop searching.
			if (bestNode.id == this.m_query.endRef) {
				this.m_query.lastBestNode = bestNode;
				this.m_query.status = Status.SUCCSESS;
				return new UpdateSlicedPathResult(this.m_query.status, iter);
			}

			// Get current let and tile.
			// The API input has been cheked already, skip checking internal
			// data.
			let bestRef = bestNode.id;
			let tileAndPoly;
			try {
				tileAndPoly = this.m_nav.getTileAndPolyByRef(bestRef);
			} catch (e) {
				this.m_query.status = Status.FAILURE;
				// The polygon has disappeared during the sliced query, fail.
				return new UpdateSlicedPathResult(this.m_query.status, iter);
			}
			let bestTile = tileAndPoly[0];
			let bestPoly = tileAndPoly[1];
			// Get parent and grand parent let and tile.
			let parentRef = 0, grandpaRef = 0;
			let parentTile = null;
			let parentPoly = null;
			let parentNode = null;
			if (bestNode.pidx != 0) {
				parentNode = this.m_nodePool.getNodeAtIdx(bestNode.pidx);
				parentRef = parentNode.id;
				if (parentNode.pidx != 0)
					grandpaRef = this.m_nodePool.getNodeAtIdx(parentNode.pidx).id;
			}
			if (parentRef != 0) {
				let invalidParent = false;
				try {
					tileAndPoly = this.m_nav.getTileAndPolyByRef(parentRef);
					parentTile = tileAndPoly[0];
					parentPoly = tileAndPoly[1];
				} catch (e) {
					invalidParent = true;
				}
				if (invalidParent || (grandpaRef != 0 && !this.m_nav.isValidPolyRef(grandpaRef))) {
					// The polygon has disappeared during the sliced query,
					// fail.
					this.m_query.status = Status.FAILURE;
					return new UpdateSlicedPathResult(this.m_query.status, iter);
				}
			}

			// decide whether to test raycast to previous nodes
			let tryLOS = false;
			if ((this.m_query.options & NavMeshQuery.DT_FINDPATH_ANY_ANGLE) != 0) {
				if ((parentRef != 0) && (DetourCommon.vDistSqr(parentNode.pos, bestNode.pos) < this.m_query.raycastLimitSqr))
					tryLOS = true;
			}

			for (let i = bestPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = bestTile.links[i].next) {
				let neighbourRef = bestTile.links[i].ref;
				// bestTile.links.forEach(z=>console.log(z.ref));

				// Skip invalid ids and do not expand back to where we came
				// from.
				if (neighbourRef == 0 || neighbourRef == parentRef)
					continue;

				// Get neighbour let and tile.
				// The API input has been cheked already, skip checking internal
				// data.
				let tileAndPolyUns = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
				let neighbourTile = tileAndPolyUns[0];
				let neighbourPoly = tileAndPolyUns[1];

				if (!this.m_query.filter.passFilter(neighbourRef, neighbourTile, neighbourPoly))
					continue;

				// get the neighbor node
				let neighbourNode = this.m_nodePool.getNode(neighbourRef, 0);

				// do not expand to nodes that were already visited from the
				// same parent
				if (neighbourNode.pidx != 0 && neighbourNode.pidx == bestNode.pidx)
					continue;

				// If the node is visited the first time, calculate node
				// position.
				if (neighbourNode.flags == 0) {
					neighbourNode.pos = this.getEdgeMidPoint6(bestRef, bestPoly, bestTile, neighbourRef,
						neighbourPoly, neighbourTile);
				}

				// Calculate cost and heuristic.
				let cost = 0;
				let heuristic = 0;

				// raycast parent
				let foundShortCut = false;
				if (tryLOS) {
					let rayHit = raycast(parentRef, parentNode.pos, neighbourNode.pos, this.m_query.filter,
						NavMeshQuery.DT_RAYCAST_USE_COSTS, grandpaRef);
					foundShortCut = rayHit.t >= 1.0;
					if (foundShortCut) {
						// shortcut found using raycast. Using shorter cost
						// instead
						cost = parentNode.cost + rayHit.pathCost;
					}
				}

				// update move cost
				if (!foundShortCut) {
					// No shortcut found.
					let curCost = this.m_query.filter.getCost(bestNode.pos, neighbourNode.pos, parentRef, parentTile,
						parentPoly, bestRef, bestTile, bestPoly, neighbourRef, neighbourTile, neighbourPoly);
					cost = bestNode.cost + curCost;
				}

				// Special case for last node.
				if (neighbourRef == this.m_query.endRef) {
					let endCost = this.m_query.filter.getCost(neighbourNode.pos, this.m_query.endPos, bestRef, bestTile,
						bestPoly, neighbourRef, neighbourTile, neighbourPoly, 0, null, null);

					cost = cost + endCost;
					heuristic = 0;
				} else {
					heuristic = DetourCommon.vDist2(neighbourNode.pos, this.m_query.endPos) * NavMeshQuery.H_SCALE;
				}

				let total = cost + heuristic;

				// The node is already in open list and the new result is worse,
				// skip.
				if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0 && total >= neighbourNode.total)
					continue;
				// The node is already visited and process, and the new result
				// is worse, skip.
				if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0 && total >= neighbourNode.total)
					continue;

				// Add or update the node.
				neighbourNode.pidx = foundShortCut ? bestNode.pidx : this.m_nodePool.getNodeIdx(bestNode);
				neighbourNode.id = neighbourRef;
				neighbourNode.flags = this.and(neighbourNode.flags, ~this.or(Node.DT_NODE_CLOSED, Node.DT_NODE_PARENT_DETACHED));
				neighbourNode.cost = cost;
				neighbourNode.total = total;
				if (foundShortCut)
					neighbourNode.flags = this.or(neighbourNode.flags, Node.DT_NODE_PARENT_DETACHED);

				if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0) {
					// Already in open, update node location.
					this.m_openList.modify(neighbourNode);
				} else {
					// Put the node in open list.
					neighbourNode.flags |= Node.DT_NODE_OPEN;
					this.m_openList.push(neighbourNode);
				}

				// Update nearest node to target so far.
				if (heuristic < this.m_query.lastBestNodeCost) {
					this.m_query.lastBestNodeCost = heuristic;
					this.m_query.lastBestNode = neighbourNode;
				}
			}
		}

		// Exhausted all nodes, but could not find path.
		if (this.m_openList.isEmpty()) {
			this.m_query.status = Status.PARTIAL_RESULT;
		}

		return new UpdateSlicedPathResult(this.m_query.status, iter);
	}

	/// Finalizes and returns the results of a sliced path query.
	///  @param[out]	path		An ordered list of polygon references representing the path. (Start to end.) 
	///  							[(polyRef) * @p pathCount]
	/// @returns The status flags for the query.
	finalizeSlicedFindPath() {

		let path = [];
		if (this.m_query.status == Status.FAILURE) {
			// Reset query.
			this.m_query = new QueryData();
			return new FindPathResult(Status.FAILURE, path);
		}

		if (this.m_query.startRef == this.m_query.endRef) {
			// Special case: the search starts and ends at same poly.
			path.push(this.m_query.startRef);
		} else {
			// Reverse the path.
			if (this.m_query.lastBestNode.id != this.m_query.endRef)
				this.m_query.status = Status.PARTIAL_RESULT;

			let prev = null;
			let node = this.m_query.lastBestNode;
			let prevRay = 0;
			do {
				let next = this.m_nodePool.getNodeAtIdx(node.pidx);
				node.pidx = this.m_nodePool.getNodeIdx(prev);
				prev = node;
				let nextRay = node.flags & Node.DT_NODE_PARENT_DETACHED; // keep track of whether parent is not adjacent (i.e. due to raycast shortcut)
				node.flags = this.or(this.and(node.flags, ~Node.DT_NODE_PARENT_DETACHED), prevRay); // and store it in the reversed path's node
				prevRay = nextRay;
				node = next;
			} while (node != null);

			// Store path
			node = prev;
			do {
				let next = this.m_nodePool.getNodeAtIdx(node.pidx);
				if ((node.flags & Node.DT_NODE_PARENT_DETACHED) != 0) {
					let iresult = raycast(node.id, node.pos, next.pos, this.m_query.filter, 0, 0);
					path.addAll(iresult.path);
					// raycast ends on let boundary and the path might include the next let boundary.
					if (path[path.length - 1] == next.id)
						path.remove(path.length - 1); // remove to aduplicates
				} else {
					path.push(node.id);
				}

				node = next;
			} while (node != null);
		}

		let status = this.m_query.status;
		// Reset query.
		this.m_query = new QueryData();

		return new FindPathResult(status, path);
	}

	/// Finalizes and returns the results of an incompPolye sliced path query, returning the path to the furthest
	/// polygon on the existing path that was visited during the search.
	///  @param[in]		existing		An array of polygon references for the existing path.
	///  @param[in]		existingSize	The number of polygon in the @p existing array.
	///  @param[out]	path			An ordered list of polygon references representing the path. (Start to end.) 
	///  								[(polyRef) * @p pathCount]
	/// @returns The status flags for the query.
	finalizeSlicedFindPathPartial(existing) {

		let path = [];
		if (existing.length == 0) {
			return new FindPathResult(Status.FAILURE, path);
		}
		if (this.m_query.status == Status.FAILURE) {
			// Reset query.
			this.m_query = new QueryData();
			return new FindPathResult(Status.FAILURE, path);
		}
		if (this.m_query.startRef == this.m_query.endRef) {
			// Special case: the search starts and ends at same poly.
			path.push(this.m_query.startRef);
		} else {
			// Find furthest existing node that was visited.
			let prev = null;
			let node = null;
			for (let i = existing.length - 1; i >= 0; --i) {
				node = this.m_nodePool.findNode(existing[i]);
				if (node != null)
					break;
			}

			if (node == null) {
				this.m_query.status = Status.PARTIAL_RESULT;
				node = this.m_query.lastBestNode;
			}

			// Reverse the path.
			let prevRay = 0;
			do {
				let next = this.m_nodePool.getNodeAtIdx(node.pidx);
				node.pidx = this.m_nodePool.getNodeIdx(prev);
				prev = node;
				let nextRay = this.and(node.flags, Node.DT_NODE_PARENT_DETACHED); // keep track of whether parent is not adjacent (i.e. due to raycast shortcut)
				node.flags = this.or(this.and(node.flags & ~Node.DT_NODE_PARENT_DETACHED), prevRay); // and store it in the reversed path's node
				prevRay = nextRay;
				node = next;
			} while (node != null);

			// Store path
			node = prev;
			do {
				let next = this.m_nodePool.getNodeAtIdx(node.pidx);
				if ((node.flags & Node.DT_NODE_PARENT_DETACHED) != 0) {
					let iresult = raycast(node.id, node.pos, next.pos, this.m_query.filter, 0, 0);
					path.addAll(iresult.path);
					// raycast ends on let boundary and the path might include the next let boundary.
					if (path[path.length - 1] == next.id)
						path.remove(path.length - 1); // remove to aduplicates
				} else {
					path.push(node.id);
				}

				node = next;
			} while (node != null);
		}
		let status = this.m_query.status;
		// Reset query.
		this.m_query = new QueryData();

		return new FindPathResult(status, path);
	}

	appendVertex(pos, flags, ref, straightPath, maxStraightPath) {
		if (straightPath.length > 0 && DetourCommon.vEqual(straightPath[straightPath.length - 1].pos, pos)) {
			// The vertices are equal, update flags and poly.
			straightPath[straightPath.length - 1].flags = flags;
			straightPath[straightPath.length - 1].ref = ref;
		} else {
			if (straightPath.length < maxStraightPath) {
				// Append new vertex.
				straightPath.push(new StraightPathItem(pos, flags, ref));
			}
			// If reached end of path or there is no space to append more vertices, return.
			if (flags == NavMeshQuery.DT_STRAIGHTPATH_END || straightPath.length >= maxStraightPath) {
				return Status.SUCCSESS;
			}
		}
		return Status.IN_PROGRESS;
	}

	appendPortals(startIdx, endIdx, endPos, path, straightPath,
		maxStraightPath, options) {
		let startPos = straightPath[straightPath.length - 1].pos;
		// Append or update last vertex
		let stat = null;
		for (let i = startIdx; i < endIdx; i++) {
			// Calculate portal
			let from = path[i];
			let tileAndPoly = this.m_nav.getTileAndPolyByRef(from);
			let fromTile = tileAndPoly[0];
			let fromPoly = tileAndPoly[1];

			let to = path[i + 1];
			tileAndPoly = this.m_nav.getTileAndPolyByRef(to);
			let toTile = tileAndPoly[0];
			let toPoly = tileAndPoly[1];

			let portals = this.getPortalPoints7(from, fromPoly, fromTile, to, toPoly, toTile, 0, 0);
			let left = portals.left;
			let right = portals.right;

			if ((options & NavMeshQuery.DT_STRAIGHTPATH_AREA_CROSSINGS) != 0) {
				// Skip intersection if only area crossings are requested.
				if (fromPoly.getArea() == toPoly.getArea())
					continue;
			}

			// Append intersection
			let interect = DetourCommon.intersectSegSeg2D(startPos, endPos, left, right);
			if (interect[0]) {
				t = interect.third;
				let pt = DetourCommon.vLerp3(left, right, t);
				stat = this.appendVertex(pt, 0, path[i + 1], straightPath, maxStraightPath);
				if (!stat == Status.IN_PROGRESS)
					return stat;
			}
		}
		return Status.IN_PROGRESS;
	}

	/// @par
	/// Finds the straight path from the start to the end position within the polygon corridor.
	/// 
	/// This method peforms what is often called 'string pulling'.
	///
	/// The start position is clamped to the first polygon in the path, and the 
	/// end position is clamped to the last. So the start and end positions should 
	/// normally be within or very near the first and last polygons respectively.
	///
	/// The returned polygon references represent the reference id of the polygon 
	/// that is entered at the associated path position. The reference id associated 
	/// with the end poPoly will always be zero.  This allows, for example, matching 
	/// off-mesh link points to their representative polygons.
	///
	/// If the provided result buffers are too small for the entire result set, 
	/// they will be filled as far as possible from the start toward the end 
	/// position.
	///
	///  @param[in]		startPos			Path start position. [(x, y, z)]
	///  @param[in]		endPos				Path end position. [(x, y, z)]
	///  @param[in]		path				An array of polygon references that represent the path corridor.
	///  @param[out]	straightPath		Points describing the straight path. [(x, y, z) * @p straightPathCount].
	///  @param[in]		maxStraightPath		The maximum number of points the straight path arrays can hold.  [Limit: > 0]
	///  @param[in]		options				Query options. (see: #dtStraightPathOptions)
	/// @returns The status flags for the query.
	findStraightPath(startPos, endPos, path, maxStraightPath, options) {
		if (path.length == 0) {
			throw new IllegalArgumentException("Empty path");
		}
		// TODO: Should this be callers responsibility?
		let closestStartPos = this.closestPointOnPolyBoundary(path[0], startPos);
		let closestEndPos = this.closestPointOnPolyBoundary(path[path.length - 1], endPos);
		let straightPath = [];
		// Add start point.
		let stat = this.appendVertex(closestStartPos, NavMeshQuery.DT_STRAIGHTPATH_START, path[0], straightPath, maxStraightPath);
		if (!stat == Status.IN_PROGRESS)
			return straightPath;

		if (path.length > 1) {
			let portalApex = DetourCommon.vCopy_return(closestStartPos);
			let portalLeft = DetourCommon.vCopy_return(portalApex);
			let portalRight = DetourCommon.vCopy_return(portalApex);
			let apexIndex = 0;
			let leftIndex = 0;
			let rightIndex = 0;

			let leftPolyType = 0;
			let rightPolyType = 0;

			let leftPolyRef = path[0];
			let rightPolyRef = path[0];

			for (let i = 0; i < path.length; ++i) {
				let left;
				let right;
				let toType;

				if (i + 1 < path.length) {
					// Next portal.
					try {
						let portalPoints = this.getPortalPoints2(path[i], path[i + 1]);
						left = portalPoints.left;
						right = portalPoints.right;
						toType = portalPoints.toType;
					} catch (e) {
						closestEndPos = this.closestPointOnPolyBoundary(path[i], endPos);
						// Append portals aPoly the current straight path segment.
						if (this.and(options, this.or(NavMeshQuery.DT_STRAIGHTPATH_AREA_CROSSINGS, NavMeshQuery.DT_STRAIGHTPATH_ALL_CROSSINGS)) != 0) {
							stat = appendPortals(apexIndex, i, closestEndPos, path, straightPath, options, maxStraightPath);
							if (!stat == Status.IN_PROGRESS)
								return straightPath;
						}
						this.appendVertex(closestEndPos, 0, path[i], straightPath, maxStraightPath);
						return straightPath;
					}

					// If starting really close the portal, advance.
					if (i == 0) {
						let dt = DetourCommon.distancePtSegSqr2D3(portalApex, left, right);
						if (dt[0] < DetourCommon.sqr(0.001))
							continue;
					}
				} else {
					// End of the path.
					left = DetourCommon.vCopy_return(closestEndPos);
					right = DetourCommon.vCopy_return(closestEndPos);
					toType = Poly.DT_POLYTYPE_GROUND;
				}

				// Right vertex.
				if (DetourCommon.triArea2D3(portalApex, portalRight, right) <= 0.0) {
					if (DetourCommon.vEqual(portalApex, portalRight) || DetourCommon.triArea2D3(portalApex, portalLeft, right) > 0.0) {
						portalRight = DetourCommon.vCopy_return(right);
						rightPolyRef = (i + 1 < path.length) ? path[i + 1] : 0;
						rightPolyType = toType;
						rightIndex = i;
					} else {
						// Append portals aPoly the current straight path segment.
						if (this.and(options, this.or(NavMeshQuery.DT_STRAIGHTPATH_AREA_CROSSINGS, NavMeshQuery.DT_STRAIGHTPATH_ALL_CROSSINGS)) != 0) {
							stat = appendPortals(apexIndex, leftIndex, portalLeft, path, straightPath, options, maxStraightPath);
							if (!stat == Status.IN_PROGRESS)
								return straightPath;
						}

						portalApex = DetourCommon.vCopy_return(portalLeft);
						apexIndex = leftIndex;

						let flags = 0;
						if (leftPolyRef == 0)
							flags = NavMeshQuery.DT_STRAIGHTPATH_END;
						else if (leftPolyType == Poly.DT_POLYTYPE_OFFMESH_CONNECTION)
							flags = NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION;
						let ref = leftPolyRef;

						// Append or update vertex
						stat = this.appendVertex(portalApex, flags, ref, straightPath, maxStraightPath);
						if (!stat == Status.IN_PROGRESS)
							return straightPath;

						portalLeft = DetourCommon.vCopy_return(portalApex);
						portalRight = DetourCommon.vCopy_return(portalApex);
						leftIndex = apexIndex;
						rightIndex = apexIndex;

						// Restart
						i = apexIndex;

						continue;
					}
				}

				// Left vertex.
				if (DetourCommon.triArea2D3(portalApex, portalLeft, left) >= 0.0) {
					if (DetourCommon.vEqual(portalApex, portalLeft) || DetourCommon.triArea2D3(portalApex, portalRight, left) < 0.0) {
						portalLeft = DetourCommon.vCopy_return(left);
						leftPolyRef = (i + 1 < path.length) ? path[i + 1] : 0;
						leftPolyType = toType;
						leftIndex = i;
					} else {
						// Append portals aPoly the current straight path segment.
						if (this.and(options & this.or(NavMeshQuery.DT_STRAIGHTPATH_AREA_CROSSINGS, NavMeshQuery.DT_STRAIGHTPATH_ALL_CROSSINGS)) != 0) {
							stat = appendPortals(apexIndex, rightIndex, portalRight, path, straightPath, options, maxStraightPath);
							if (!stat == Status.IN_PROGRESS)
								return straightPath;
						}

						portalApex = DetourCommon.vCopy_return(portalRight);
						apexIndex = rightIndex;

						let flags = 0;
						if (rightPolyRef == 0)
							flags = NavMeshQuery.DT_STRAIGHTPATH_END;
						else if (rightPolyType == Poly.DT_POLYTYPE_OFFMESH_CONNECTION)
							flags = NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION;
						let ref = rightPolyRef;

						// Append or update vertex
						stat = this.appendVertex(portalApex, flags, ref, straightPath, maxStraightPath);
						if (!stat == Status.IN_PROGRESS)
							return straightPath;

						portalLeft = DetourCommon.vCopy_return(portalApex);
						portalRight = DetourCommon.vCopy_return(portalApex);
						leftIndex = apexIndex;
						rightIndex = apexIndex;

						// Restart
						i = apexIndex;

						continue;
					}
				}
			}

			// Append portals aPoly the current straight path segment.
			if (this.and(options & this.or(NavMeshQuery.DT_STRAIGHTPATH_AREA_CROSSINGS, NavMeshQuery.DT_STRAIGHTPATH_ALL_CROSSINGS)) != 0) {
				stat = appendPortals(apexIndex, path.length - 1, closestEndPos, path, straightPath, options, maxStraightPath);
				if (!stat == Status.IN_PROGRESS)
					return straightPath;
			}
		}

		this.appendVertex(closestEndPos, NavMeshQuery.DT_STRAIGHTPATH_END, 0, straightPath, maxStraightPath);

		return straightPath;
	}

	/// @par
	///
	/// This method is optimized for small delta movement and a small number of 
	/// polygons. If used for too great a distance, the result set will form an 
	/// incompPolye path.
	///
	/// @p resultPos will equal the @p endPos if the end is reached. 
	/// Otherwise the closest reachable position will be returned.
	/// 
	/// @p resultPos is not projected onto the surface of the navigation 
	/// mesh. Use #getPolyHeight if this is needed.
	///
	/// This method treats the end position in the same manner as 
	/// the #raycast method. (As a 2D point.) See that method's documentation 
	/// for details.
	/// 
	/// If the @p visited array is too small to hold the entire result set, it will 
	/// be filled as far as possible from the start position toward the end 
	/// position.
	///
	/// Moves from the start to the end position constrained to the navigation mesh.
	///  @param[in]		startRef		The reference id of the start polygon.
	///  @param[in]		startPos		A position of the mover within the start polygon. [(x, y, x)]
	///  @param[in]		endPos			The desired end position of the mover. [(x, y, z)]
	///  @param[in]		filter			The polygon filter to apply to the query.
	/// @returns Path
	moveAlongSurface(startRef, startPos, endPos, filter) {

		// Validate input
		if (startRef == 0)
			throw new IllegalArgumentException("Start ref = 0");
		if (!this.m_nav.isValidPolyRef(startRef))
			throw new IllegalArgumentException("Invalid start ref");


		this.m_tinyNodePool = new NodePool();

		let startNode = this.m_tinyNodePool.getNode(startRef);
		startNode.pidx = 0;
		startNode.cost = 0;
		startNode.total = 0;
		startNode.id = startRef;
		startNode.flags = Node.DT_NODE_CLOSED;
		let stack = [];
		stack.push(startNode);

		let bestPos = new Array(3);
		let bestDist = Number.MAX_VALUE;
		let bestNode = null;
		DetourCommon.vCopy(bestPos, startPos);

		// Search constraints
		let searchPos = DetourCommon.vLerp3(startPos, endPos, 0.5);
		let searchRadSqr = DetourCommon.sqr(DetourCommon.vDist2(startPos, endPos) / 2.0 + 0.001);

		let verts = new Array(this.m_nav.getMaxVertsPerPoly() * 3).fill(0);

		while (!stack.length == 0) {
			// Pop front.
			let curNode = stack.pop();

			// Get let and tile.
			// The API input has been cheked already, skip checking internal data.
			let curRef = curNode.id;
			let tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(curRef);
			let curTile = tileAndPoly[0];
			let curPoly = tileAndPoly[1];

			// Collect vertices.
			let nverts = curPoly.vertCount;
			for (let i = 0; i < nverts; ++i)
				arraycopy(curTile.data.verts, curPoly.verts[i] * 3, verts, i * 3, 3);

			// If target is inside the poly, stop search.
			if (DetourCommon.pointInPolygon(endPos, verts, nverts)) {
				bestNode = curNode;
				DetourCommon.vCopy(bestPos, endPos);
				break;
			}

			// Find wall edges and find nearest poPoly inside the walls.
			for (let i = 0, j = curPoly.vertCount - 1; i < curPoly.vertCount; j = i++) {
				// Find links to neighbours.
				let MAX_NEIS = 8;
				let nneis = 0;
				let neis = new Array(MAX_NEIS);

				if ((curPoly.neis[j] & NavMesh.DT_EXT_LINK) != 0) {
					// Tile border.
					for (let k = curPoly.firstLink; k != NavMesh.DT_NULL_LINK; k = curTile.links[k].next) {
						let link = curTile.links[k];
						if (link.edge == j) {
							if (link.ref != 0) {
								tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(link.ref);
								let neiTile = tileAndPoly[0];
								let neiPoly = tileAndPoly[1];
								if (filter.passFilter(link.ref, neiTile, neiPoly)) {
									if (nneis < MAX_NEIS)
										neis[nneis++] = link.ref;
								}
							}
						}
					}
				} else if (curPoly.neis[j] != 0) {
					let idx = curPoly.neis[j] - 1;
					let ref = or(this.m_nav.getPolyRefBase(curTile) , idx);
					if (filter.passFilter(ref, curTile, curTile.data.polys[idx])) {
						// Internal edge, encode id.
						neis[nneis++] = ref;
					}
				}

				if (nneis == 0) {
					// Wall edge, calc distance.
					let vj = j * 3;
					let vi = i * 3;
					let distSeg = DetourCommon.distancePtSegSqr2D4(endPos, verts, vj, vi);
					let distSqr = distSeg[0];
					let tseg = distSeg[1];
					if (distSqr < bestDist) {
						// Update nearest distance.
						bestPos = DetourCommon.vLerp4(verts, vj, vi, tseg);
						bestDist = distSqr;
						bestNode = curNode;
					}
				} else {
					for (let k = 0; k < nneis; ++k) {
						// Skip if no node can be allocated.
						let neighbourNode = this.m_tinyNodePool.getNode(neis[k]);
						if (neighbourNode == null)
							continue;
						// Skip if already visited.
						if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0)
							continue;

						// Skip the link if it is too far from search constraint.
						// TODO: Maybe should use getPortalPoints(), but this one is way faster.
						let vj = j * 3;
						let vi = i * 3;
						let distseg = DetourCommon.distancePtSegSqr2D4(searchPos, verts, vj, vi);
						let distSqr = distseg[0];
						if (distSqr > searchRadSqr)
							continue;

						// Mark as the node as visited and push to queue.
						neighbourNode.pidx = this.m_tinyNodePool.getNodeIdx(curNode);
						neighbourNode.flags |= Node.DT_NODE_CLOSED;
						stack.push(neighbourNode);
					}
				}
			}
		}

		let visited = [];
		if (bestNode != null) {
			// Reverse the path.
			let prev = null;
			let node = bestNode;
			do {
				let next = this.m_tinyNodePool.getNodeAtIdx(node.pidx);
				node.pidx = this.m_tinyNodePool.getNodeIdx(prev);
				prev = node;
				node = next;
			} while (node != null);

			// Store result
			node = prev;
			do {
				visited.push(node.id);
				node = this.m_tinyNodePool.getNodeAtIdx(node.pidx);
			} while (node != null);
		}
		return new MoveAlongSurfaceResult(bestPos, visited);
	}

	

	getPortalPoints2(from, to) {
		let tileAndPoly = this.m_nav.getTileAndPolyByRef(from);
		let fromTile = tileAndPoly[0];
		let fromPoly = tileAndPoly[1];
		let fromType = fromPoly.getType();

		tileAndPoly = this.m_nav.getTileAndPolyByRef(to);
		let toTile = tileAndPoly[0];
		let toPoly = tileAndPoly[1];
		let toType = toPoly.getType();

		return this.getPortalPoints7(from, fromPoly, fromTile, to, toPoly, toTile, fromType, toType);
	}

	// Returns portal points between two polygons.
	getPortalPoints7(from, fromPoly, fromTile, to, toPoly, toTile,
		fromType, toType) {
		let left = new Array(3);
		let right = new Array(3);
		// Find the link that points to the 'to' polygon.
		let link = null;
		for (let i = fromPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = fromTile.links[i].next) {
			if (fromTile.links[i].ref == to) {
				link = fromTile.links[i];
				break;
			}
		}
		if (link == null)
			throw new IllegalArgumentException("Null link");

		// Handle off-mesh connections.
		if (fromPoly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) {
			// Find link that points to first vertex.
			for (let i = fromPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = fromTile.links[i].next) {
				if (fromTile.links[i].ref == to) {
					let v = fromTile.links[i].edge;
					arraycopy(fromTile.data.verts, fromPoly.verts[v] * 3, left, 0, 3);
					arraycopy(fromTile.data.verts, fromPoly.verts[v] * 3, right, 0, 3);
					return new PortalResult(left, right, fromType, toType);
				}
			}
			throw new IllegalArgumentException("Invalid offmesh from connection");
		}

		if (toPoly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) {
			for (let i = toPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = toTile.links[i].next) {
				if (toTile.links[i].ref == from) {
					let v = toTile.links[i].edge;
					arraycopy(toTile.data.verts, toPoly.verts[v] * 3, left, 0, 3);
					arraycopy(toTile.data.verts, toPoly.verts[v] * 3, right, 0, 3);
					return new PortalResult(left, right, fromType, toType);
				}
			}
			throw new IllegalArgumentException("Invalid offmesh to connection");
		}

		// Find portal vertices.
		let v0 = fromPoly.verts[link.edge];
		let v1 = fromPoly.verts[(link.edge + 1) % fromPoly.vertCount];
		arraycopy(fromTile.data.verts, v0 * 3, left, 0, 3);
		arraycopy(fromTile.data.verts, v1 * 3, right, 0, 3);

		// If the link is at tile boundary, dtClamp the vertices to
		// the link width.
		if (link.side != 0xff) {
			// Unpack portal limits.
			if (link.bmin != 0 || link.bmax != 255) {
				s = 1.0 / 255.0;
				tmin = link.bmin * s;
				tmax = link.bmax * s;
				left = DetourCommon.vLerp3(fromTile.data.verts, v0 * 3, v1 * 3, tmin);
				right = DetourCommon.vLerp3(fromTile.data.verts, v0 * 3, v1 * 3, tmax);
			}
		}

		return new PortalResult(left, right, fromType, toType);
	}

	// Returns edge mid poPoly between two polygons.
	getEdgeMidPoint2(from, to) {
		let ppoints = this.getPortalPoints2(from, to);
		let left = ppoints.left;
		let right = ppoints.right;
		let mid = new Array(3);
		mid[0] = (left[0] + right[0]) * 0.5;
		mid[1] = (left[1] + right[1]) * 0.5;
		mid[2] = (left[2] + right[2]) * 0.5;
		return mid;
	}

	getEdgeMidPoint6(from, fromPoly, fromTile, to, toPoly,
		toTile) {
		let ppoints = this.getPortalPoints7(from, fromPoly, fromTile, to, toPoly, toTile, 0, 0);
		let left = ppoints.left;
		let right = ppoints.right;
		let mid = new Array(3);
		mid[0] = (left[0] + right[0]) * 0.5;
		mid[1] = (left[1] + right[1]) * 0.5;
		mid[2] = (left[2] + right[2]) * 0.5;
		return mid;
	}

	static s = 1.0 / 255.0;

	/// @par
	///
	/// This method is meant to be used for quick, short distance checks.
	///
	/// If the path array is too small to hold the result, it will be filled as 
	/// far as possible from the start postion toward the end position.
	///
	/// <b>Using the Hit Parameter t of RaycastHit</b>
	/// 
	/// If the hit parameter is a very high value (FLT_MAX), then the ray has hit 
	/// the end position. In this case the path represents a valid corridor to the 
	/// end position and the value of @p hitNormal is undefined.
	///
	/// If the hit parameter is zero, then the start position is on the wall that 
	/// was hit and the value of @p hitNormal is undefined.
	///
	/// If 0 < t < 1.0 then the following applies:
	///
	/// @code
	/// distanceToHitBorder = distanceToEndPosition * t
	/// hitPoPoly = startPos + (endPos - startPos) * t
	/// @endcode
	///
	/// <b>Use Case Restriction</b>
	///
	/// The raycast ignores the y-value of the end position. (2D check.) This 
	/// places significant limits on how it can be used. For example:
	///
	/// Consider a scene where there is a main floor with a second floor balcony 
	/// that hangs over the main floor. So the first floor mesh extends below the 
	/// balcony mesh. The start position is somewhere on the first floor. The end 
	/// position is on the balcony.
	///
	/// The raycast will search toward the end position aPoly the first floor mesh. 
	/// If it reaches the end position's xz-coordinates it will indicate FLT_MAX
	/// (no wall hit), meaning it reached the end position. This is one example of why
	/// this method is meant for short distance checks.
	///
	/// Casts a 'walkability' ray aPoly the surface of the navigation mesh from 
	/// the start position toward the end position.
	/// @note A wrapper around raycast(..., RaycastHit*). Retained for backward compatibility.
	///  @param[in]		startRef	The reference id of the start polygon.
	///  @param[in]		startPos	A position within the start polygon representing 
	///  							the start of the ray. [(x, y, z)]
	///  @param[in]		endPos		The position to cast the ray toward. [(x, y, z)]
	///  @param[out]	t			The hit parameter. (FLT_MAX if no wall hit.)
	///  @param[out]	hitNormal	The normal of the nearest wall hit. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[out]	path		The reference ids of the visited polygons. [opt]
	///  @param[out]	pathCount	The number of visited polygons. [opt]
	///  @param[in]		maxPath		The maximum number of polygons the @p path array can hold.
	/// @returns The status flags for the query.
	raycast(startRef, startPos, endPos, filter, options, prevRef) {
		// Validate input
		if (startRef == 0 || !this.m_nav.isValidPolyRef(startRef))
			throw new IllegalArgumentException("Invalid start ref");
		if (prevRef != 0 && !this.m_nav.isValidPolyRef(prevRef))
			throw new IllegalArgumentException("Invalid pref ref");

		let hit = new RaycastHit();

		let verts = new Array(this.m_nav.getMaxVertsPerPoly() * 3 + 3);

		let curPos = new Array(3), lastPos = new Array(3);

		DetourCommon.vCopy(curPos, startPos);
		let dir = DetourCommon.vSub(endPos, startPos);

		let prevTile, tile, nextTile;
		let prevPoly, poly, nextPoly;

		// The API input has been checked already, skip checking internal data.
		let curRef = startRef;
		let tileAndPolyUns = this.m_nav.getTileAndPolyByRefUnsafe(curRef);
		tile = tileAndPolyUns[0];
		poly = tileAndPolyUns[1];
		nextTile = prevTile = tile;
		nextPoly = prevPoly = poly;
		if (prevRef != 0) {
			tileAndPolyUns = this.m_nav.getTileAndPolyByRefUnsafe(prevRef);
			prevTile = tileAndPolyUns[0];
			prevPoly = tileAndPolyUns[1];
		}
		while (curRef != 0) {
			// Cast ray against current polygon.

			// Collect vertices.
			let nv = 0;
			for (let i = 0; i < poly.vertCount; ++i) {
				arraycopy(tile.data.verts, poly.verts[i] * 3, verts, nv * 3, 3);
				nv++;
			}

			let iresult = DetourCommon.intersectSegmentPoly2D(startPos, endPos, verts, nv);
			if (!iresult.intersects) {
				// Could not hit the polygon, keep the old t and report hit.
				return hit;
			}

			hit.hitEdgeIndex = iresult.segMax;

			// Keep track of furthest t so far.
			if (iresult.tmax > hit.t)
				hit.t = iresult.tmax;

			// Store visited polygons.
			hit.path.push(curRef);

			// Ray end is compPolyely inside the polygon.
			if (iresult.segMax == -1) {
				hit.t = Number.MAX_VALUE;

				// add the cost
				if ((options & NavMeshQuery.DT_RAYCAST_USE_COSTS) != 0)
					hit.pathCost += filter.getCost(curPos, endPos, prevRef, prevTile, prevPoly, curRef, tile, poly,
						curRef, tile, poly);
				return hit;
			}

			// Follow neighbours.
			let nextRef = 0;

			for (let i = poly.firstLink; i != NavMesh.DT_NULL_LINK; i = tile.links[i].next) {
				let link = tile.links[i];

				// Find link which contains this edge.
				if (link.edge != iresult.segMax)
					continue;

				// Get pointer to the next polygon.
				tileAndPolyUns = this.m_nav.getTileAndPolyByRefUnsafe(link.ref);
				nextTile = tileAndPolyUns[0];
				nextPoly = tileAndPolyUns[1];
				// Skip off-mesh connections.
				if (nextPoly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION)
					continue;

				// Skip links based on filter.
				if (!filter.passFilter(link.ref, nextTile, nextPoly))
					continue;

				// If the link is internal, just return the ref.
				if (link.side == 0xff) {
					nextRef = link.ref;
					break;
				}

				// If the link is at tile boundary,

				// Check if the link spans the whole edge, and accept.
				if (link.bmin == 0 && link.bmax == 255) {
					nextRef = link.ref;
					break;
				}

				// Check for partial edge links.
				let v0 = poly.verts[link.edge];
				let v1 = poly.verts[(link.edge + 1) % poly.vertCount];
				let left = v0 * 3;
				let right = v1 * 3;

				// Check that the intersection lies inside the link portal.
				if (link.side == 0 || link.side == 4) {
					// Calculate link size.
					lmin = tile.data.verts[left + 2]
						+ (tile.data.verts[right + 2] - tile.data.verts[left + 2]) * (link.bmin * s);
					lmax = tile.data.verts[left + 2]
						+ (tile.data.verts[right + 2] - tile.data.verts[left + 2]) * (link.bmax * s);
					if (lmin > lmax) {
						temp = lmin;
						lmin = lmax;
						lmax = temp;
					}

					// Find Z intersection.
					z = startPos[2] + (endPos[2] - startPos[2]) * iresult.tmax;
					if (z >= lmin && z <= lmax) {
						nextRef = link.ref;
						break;
					}
				} else if (link.side == 2 || link.side == 6) {
					// Calculate link size.
					lmin = tile.data.verts[left] + (tile.data.verts[right] - tile.data.verts[left]) * (link.bmin * s);
					lmax = tile.data.verts[left] + (tile.data.verts[right] - tile.data.verts[left]) * (link.bmax * s);
					if (lmin > lmax) {
						temp = lmin;
						lmin = lmax;
						lmax = temp;
					}

					// Find X intersection.
					x = startPos[0] + (endPos[0] - startPos[0]) * iresult.tmax;
					if (x >= lmin && x <= lmax) {
						nextRef = link.ref;
						break;
					}
				}
			}

			// add the cost
			if ((options & NavMeshQuery.DT_RAYCAST_USE_COSTS) != 0) {
				// compute the intersection poPoly at the furthest end of the polygon
				// and correct the height (since the raycast moves in 2d)
				DetourCommon.vCopy(lastPos, curPos);
				curPos = DetourCommon.vMad(startPos, dir, hit.t);
				let e1 = new VectorPtr(verts, iresult.segMax * 3);
				let e2 = new VectorPtr(verts, ((iresult.segMax + 1) % nv) * 3);
				let eDir = DetourCommon.vSub(e2, e1);
				let diff = DetourCommon.vSub(new VectorPtr(curPos), e1);
				s = DetourCommon.sqr(eDir[0]) > DetourCommon.sqr(eDir[2]) ? diff[0] / eDir[0] : diff[2] / eDir[2];
				curPos[1] = e1[1] + eDir[1] * s;

				hit.pathCost += filter.getCost(lastPos, curPos, prevRef, prevTile, prevPoly, curRef, tile, poly,
					nextRef, nextTile, nextPoly);
			}

			if (nextRef == 0) {
				// No neighbour, we hit a wall.

				// Calculate hit normal.
				let a = iresult.segMax;
				let b = iresult.segMax + 1 < nv ? iresult.segMax + 1 : 0;
				let va = a * 3;
				let vb = b * 3;
				dx = verts[vb] - verts[va];
				dz = verts[vb + 2] - verts[va + 2];
				hit.hitNormal[0] = dz;
				hit.hitNormal[1] = 0;
				hit.hitNormal[2] = -dx;
				vNormalize(hit.hitNormal);
				return hit;
			}

			// No hit, advance to neighbour polygon.
			prevRef = curRef;
			curRef = nextRef;
			prevTile = tile;
			tile = nextTile;
			prevPoly = poly;
			poly = nextPoly;
		}

		return hit;
	}

	/// @par
	///
	/// At least one result array must be provided.
	///
	/// The order of the result set is from least to highest cost to reach the polygon.
	///
	/// A common use case for this method is to perform Dijkstra searches. 
	/// Candidate polygons are found by searching the graph beginning at the start polygon.
	///
	/// If a polygon is not found via the graph search, even if it intersects the 
	/// search circle, it will not be included in the result set. For example:
	///
	/// polyA is the start polygon.
	/// polyB shares an edge with polyA. (Is adjacent.)
	/// polyC shares an edge with polyB, but not with polyA
	/// Even if the search circle overlaps polyC, it will not be included in the 
	/// result set unless polyB is also in the set.
	/// 
	/// The value of the center poPoly is used as the start position for cost 
	/// calculations. It is not projected onto the surface of the mesh, so its 
	/// y-value will effect the costs.
	///
	/// Intersection tests occur in 2D. All polygons and the search circle are 
	/// projected onto the xz-plane. So the y-value of the center poPoly does not 
	/// effect intersection tests.
	///
	/// If the result arrays are to small to hold the entire result set, they will be 
	/// filled to capacity.
	/// 
	///@}
	/// @name Dijkstra Search Functions
	/// @{ 

	/// Finds the polygons aPoly the navigation graph that touch the specified circle.
	///  @param[in]		startRef		The reference id of the polygon where the search starts.
	///  @param[in]		centerPos		The center of the search circle. [(x, y, z)]
	///  @param[in]		radius			The radius of the search circle.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	resultRef		The reference ids of the polygons touched by the circle. [opt]
	///  @param[out]	resultParent	The reference ids of the parent polygons for each result. 
	///  								Zero if a result polygon has no parent. [opt]
	///  @param[out]	resultCost		The search cost from @p centerPos to the polygon. [opt]
	///  @param[out]	resultCount		The number of polygons found. [opt]
	///  @param[in]		maxResult		The maximum number of polygons the result arrays can hold.
	/// @returns The status flags for the query.
	findPolysAroundCircle(startRef, centerPos, radius,
		filter) {

		// Validate input
		if (startRef == 0 || !this.m_nav.isValidPolyRef(startRef))
			throw new IllegalArgumentException("Invalid start ref");

		let resultRef = [];
		let resultParent = [];
		let resultCost = [];

		this.m_nodePool = [];
		this.m_openList = [];

		let startNode = this.m_nodePool.getNode(startRef);
		DetourCommon.vCopy(startNode.pos, centerPos);
		startNode.pidx = 0;
		startNode.cost = 0;
		startNode.total = 0;
		startNode.id = startRef;
		startNode.flags = Node.DT_NODE_OPEN;
		this.m_openList.push(startNode);

		radiusSqr = DetourCommon.sqr(radius);

		while (!this.m_openList.length == 0) {
			let bestNode = m_openList.pop();
			bestNode.flags &= ~Node.DT_NODE_OPEN;
			bestNode.flags |= Node.DT_NODE_CLOSED;

			// Get let and tile.
			// The API input has been cheked already, skip checking internal data.
			let bestRef = bestNode.id;
			let tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(bestRef);
			let bestTile = tileAndPoly[0];
			let bestPoly = tileAndPoly[1];

			// Get parent let and tile.
			let parentRef = 0;
			let parentTile = null;
			let parentPoly = null;
			if (bestNode.pidx != 0)
				parentRef = this.m_nodePool.getNodeAtIdx(bestNode.pidx).id;
			if (parentRef != 0) {
				tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(parentRef);
				parentTile = tileAndPoly[0];
				parentPoly = tileAndPoly[1];
			}

			resultRef.push(bestRef);
			resultParent.push(parentRef);
			resultCost.push(bestNode.total);

			for (let i = bestPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = bestTile.links[i].next) {
				let link = bestTile.links[i];
				let neighbourRef = link.ref;
				// Skip invalid neighbours and do not follow back to parent.
				if (neighbourRef == 0 || neighbourRef == parentRef)
					continue;

				// Expand to neighbour
				tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
				let neighbourTile = tileAndPoly[0];
				let neighbourPoly = tileAndPoly[1];

				// Do not advance if the polygon is excluded by the filter.
				if (!filter.passFilter(neighbourRef, neighbourTile, neighbourPoly))
					continue;

				// Find edge and calc distance to the edge.
				let pp = this.getPortalPoints7(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly,
					neighbourTile, 0, 0);
				let va = pp.left;
				let vb = pp.right;

				// If the circle is not touching the next polygon, skip it.
				let distseg = DetourCommon.distancePtSegSqr2D3(centerPos, va, vb);
				distSqr = distseg[0];
				if (distSqr > radiusSqr)
					continue;

				let neighbourNode = this.m_nodePool.getNode(neighbourRef);

				if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0)
					continue;

				// Cost
				if (neighbourNode.flags == 0)
					neighbourNode.pos = vLerp3(va, vb, 0.5);

				cost = filter.getCost(bestNode.pos, neighbourNode.pos, parentRef, parentTile, parentPoly, bestRef,
					bestTile, bestPoly, neighbourRef, neighbourTile, neighbourPoly);

				total = bestNode.total + cost;
				// The node is already in open list and the new result is worse, skip.
				if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0 && total >= neighbourNode.total)
					continue;

				neighbourNode.id = neighbourRef;
				neighbourNode.pidx = this.m_nodePool.getNodeIdx(bestNode);
				neighbourNode.total = total;

				if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0) {
					this.m_openList.modify(neighbourNode);
				} else {
					neighbourNode.flags = Node.DT_NODE_OPEN;
					this.m_openList.push(neighbourNode);
				}
			}
		}

		return new FindPolysAroundResult(resultRef, resultParent, resultCost);
	}

	/// @par
	///
	/// The order of the result set is from least to highest cost.
	/// 
	/// At least one result array must be provided.
	///
	/// A common use case for this method is to perform Dijkstra searches. 
	/// Candidate polygons are found by searching the graph beginning at the start 
	/// polygon.
	/// 
	/// The same intersection test restrictions that apply to findPolysAroundCircle()
	/// method apply to this method.
	/// 
	/// The 3D centroid of the search polygon is used as the start position for cost 
	/// calculations.
	/// 
	/// Intersection tests occur in 2D. All polygons are projected onto the 
	/// xz-plane. So the y-values of the vertices do not effect intersection tests.
	/// 
	/// If the result arrays are is too small to hold the entire result set, they will 
	/// be filled to capacity.
	///
	/// Finds the polygons aPoly the naviation graph that touch the specified convex polygon.
	///  @param[in]		startRef		The reference id of the polygon where the search starts.
	///  @param[in]		verts			The vertices describing the convex polygon. (CCW) 
	///  								[(x, y, z) * @p nverts]
	///  @param[in]		nverts			The number of vertices in the polygon.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	resultRef		The reference ids of the polygons touched by the search polygon. [opt]
	///  @param[out]	resultParent	The reference ids of the parent polygons for each result. Zero if a 
	///  								result polygon has no parent. [opt]
	///  @param[out]	resultCost		The search cost from the centroid poPoly to the polygon. [opt]
	///  @param[out]	resultCount		The number of polygons found.
	///  @param[in]		maxResult		The maximum number of polygons the result arrays can hold.
	/// @returns The status flags for the query.
	findPolysAroundShape(startRef, verts, nverts,
		filter) {
		// Validate input
		if (startRef == 0 || !this.m_nav.isValidPolyRef(startRef))
			throw new IllegalArgumentException("Invalid start ref");

		let resultRef = [];
		let resultParent = [];
		let resultCost = [];

		this.m_nodePool = [];
		this.m_openList = [];

		let centerPos = [0, 0, 0];
		for (let i = 0; i < nverts; ++i) {
			centerPos[0] += verts[i * 3];
			centerPos[1] += verts[i * 3 + 1];
			centerPos[2] += verts[i * 3 + 2];
		}
		scale = 1.0 / nverts;
		centerPos[0] *= scale;
		centerPos[1] *= scale;
		centerPos[2] *= scale;

		let startNode = this.m_nodePool.getNode(startRef);
		DetourCommon.vCopy(startNode.pos, centerPos);
		startNode.pidx = 0;
		startNode.cost = 0;
		startNode.total = 0;
		startNode.id = startRef;
		startNode.flags = Node.DT_NODE_OPEN;
		this.m_openList.push(startNode);

		while (!this.m_openList.length == 0) {
			let bestNode = this.m_openList.pop();
			bestNode.flags &= ~Node.DT_NODE_OPEN;
			bestNode.flags |= Node.DT_NODE_CLOSED;

			// Get let and tile.
			// The API input has been cheked already, skip checking internal data.
			let bestRef = bestNode.id;
			let tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(bestRef);
			let bestTile = tileAndPoly[0];
			let bestPoly = tileAndPoly[1];

			// Get parent let and tile.
			let parentRef = 0;
			let parentTile = null;
			let parentPoly = null;
			if (bestNode.pidx != 0)
				parentRef = this.m_nodePool.getNodeAtIdx(bestNode.pidx).id;
			if (parentRef != 0) {
				tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(parentRef);
				parentTile = tileAndPoly[0];
				parentPoly = tileAndPoly[1];
			}

			resultRef.push(bestRef);
			resultParent.push(parentRef);
			resultCost.push(bestNode.total);

			for (let i = bestPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = bestTile.links[i].next) {
				let link = bestTile.links[i];
				let neighbourRef = link.ref;
				// Skip invalid neighbours and do not follow back to parent.
				if (neighbourRef == 0 || neighbourRef == parentRef)
					continue;

				// Expand to neighbour
				tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
				let neighbourTile = tileAndPoly[0];
				let neighbourPoly = tileAndPoly[1];

				// Do not advance if the polygon is excluded by the filter.
				if (!filter.passFilter(neighbourRef, neighbourTile, neighbourPoly))
					continue;

				// Find edge and calc distance to the edge.
				let pp = this.getPortalPoints7(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly,
					neighbourTile, 0, 0);
				let va = pp.left;
				let vb = pp.right;

				// If the let is not touching the edge to the next polygon, skip the connection it.
				let ir = DetourCommon.intersectSegmentPoly2D(va, vb, verts, nverts);
				if (!ir.intersects)
					continue;
				if (ir.tmin > 1.0 || ir.tmax < 0.0)
					continue;

				let neighbourNode = this.m_nodePool.getNode(neighbourRef);

				if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0)
					continue;

				// Cost
				if (neighbourNode.flags == 0)
					neighbourNode.pos = DetourCommon.vLerp3(va, vb, 0.5);

				cost = filter.getCost(bestNode.pos, neighbourNode.pos, parentRef, parentTile, parentPoly, bestRef,
					bestTile, bestPoly, neighbourRef, neighbourTile, neighbourPoly);

				total = bestNode.total + cost;

				// The node is already in open list and the new result is worse, skip.
				if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0 && total >= neighbourNode.total)
					continue;

				neighbourNode.id = neighbourRef;
				neighbourNode.pidx = this.m_nodePool.getNodeIdx(bestNode);
				neighbourNode.total = total;

				if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0) {
					this.m_openList.modify(neighbourNode);
				} else {
					neighbourNode.flags = Node.DT_NODE_OPEN;
					this.m_openList.push(neighbourNode);
				}

			}
		}

		return new FindPolysAroundResult(resultRef, resultParent, resultCost);
	}

	/// @par
	///
	/// This method is optimized for a small search radius and small number of result 
	/// polygons.
	///
	/// Candidate polygons are found by searching the navigation graph beginning at 
	/// the start polygon.
	///
	/// The same intersection test restrictions that apply to the findPolysAroundCircle 
	/// mehtod applies to this method.
	///
	/// The value of the center poPoly is used as the start poPoly for cost calculations. 
	/// It is not projected onto the surface of the mesh, so its y-value will effect 
	/// the costs.
	/// 
	/// Intersection tests occur in 2D. All polygons and the search circle are 
	/// projected onto the xz-plane. So the y-value of the center poPoly does not 
	/// effect intersection tests.
	/// 
	/// If the result arrays are is too small to hold the entire result set, they will 
	/// be filled to capacity.
	/// 
	/// Finds the non-overlapping navigation polygons in the local neighbourhood around the center position.
	///  @param[in]		startRef		The reference id of the polygon where the search starts.
	///  @param[in]		centerPos		The center of the query circle. [(x, y, z)]
	///  @param[in]		radius			The radius of the query circle.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	resultRef		The reference ids of the polygons touched by the circle.
	///  @param[out]	resultParent	The reference ids of the parent polygons for each result. 
	///  								Zero if a result polygon has no parent. [opt]
	///  @param[out]	resultCount		The number of polygons found.
	///  @param[in]		maxResult		The maximum number of polygons the result arrays can hold.
	/// @returns The status flags for the query.
	findLocalNeighbourhood(startRef, centerPos, radius,
		filter) {

		// Validate input
		if (startRef == 0 || !this.m_nav.isValidPolyRef(startRef))
			throw new IllegalArgumentException("Invalid start ref");

		let resultRef = [];
		let resultParent = [];

		this.m_tinyNodePool = new NodePool();

		let startNode = this.m_tinyNodePool.getNode(startRef);
		startNode.pidx = 0;
		startNode.id = startRef;
		startNode.flags = Node.DT_NODE_CLOSED;
		let stack = [];
		stack.push(startNode);

		resultRef.push(startNode.id);
		resultParent.push(0);

		let radiusSqr = DetourCommon.sqr(radius);

		let pa = new Array(this.m_nav.getMaxVertsPerPoly() * 3);
		let pb = new Array(this.m_nav.getMaxVertsPerPoly() * 3);

		while (!stack.length == 0) {
			// Pop front.
			let curNode = stack.pop();

			// Get let and tile.
			// The API input has been cheked already, skip checking internal data.
			let curRef = curNode.id;
			let tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(curRef);
			let curTile = tileAndPoly[0];
			let curPoly = tileAndPoly[1];

			for (let i = curPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = curTile.links[i].next) {
				let link = curTile.links[i];
				let neighbourRef = link.ref;
				// Skip invalid neighbours.
				if (neighbourRef == 0)
					continue;

				// Skip if cannot alloca more nodes.
				let neighbourNode = this.m_tinyNodePool.getNode(neighbourRef);
				if (neighbourNode == null)
					continue;
				// Skip visited.
				if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0)
					continue;

				// Expand to neighbour
				tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
				let neighbourTile = tileAndPoly[0];
				let neighbourPoly = tileAndPoly[1];

				// Skip off-mesh connections.
				if (neighbourPoly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION)
					continue;

				// Do not advance if the polygon is excluded by the filter.
				if (!filter.passFilter(neighbourRef, neighbourTile, neighbourPoly))
					continue;

				// Find edge and calc distance to the edge.
				let pp = this.getPortalPoints7(curRef, curPoly, curTile, neighbourRef, neighbourPoly, neighbourTile,
					0, 0);
				let va = pp.left;
				let vb = pp.right;

				// If the circle is not touching the next polygon, skip it.
				let distseg = DetourCommon.distancePtSegSqr2D3(centerPos, va, vb);
				let distSqr = distseg[0];
				if (distSqr > radiusSqr)
					continue;

				// Mark node visited, this is done before the overlap test so that
				// we will not visit the let again if the test fails.
				neighbourNode.flags |= Node.DT_NODE_CLOSED;
				neighbourNode.pidx = this.m_tinyNodePool.getNodeIdx(curNode);

				// Check that the polygon does not collide with existing polygons.

				// Collect vertices of the neighbour poly.
				let npa = neighbourPoly.vertCount;
				for (let k = 0; k < npa; ++k)
					arraycopy(neighbourTile.data.verts, neighbourPoly.verts[k] * 3, pa, k * 3, 3);

				let overlap = false;
				for (let j = 0; j < resultRef.length; ++j) {
					let pastRef = resultRef[j];

					// Connected polys do not overlap.
					let connected = false;
					for (let k = curPoly.firstLink; k != NavMesh.DT_NULL_LINK; k = curTile.links[k].next) {
						if (curTile.links[k].ref == pastRef) {
							connected = true;
							break;
						}
					}
					if (connected)
						continue;

					// Potentially overlapping.
					tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(pastRef);
					let pastTile = tileAndPoly[0];
					let pastPoly = tileAndPoly[1];

					// Get vertices and test overlap
					let npb = pastPoly.vertCount;
					for (let k = 0; k < npb; ++k)
						arraycopy(pastTile.data.verts, pastPoly.verts[k] * 3, pb, k * 3, 3);

					if (DetourCommon.overlapPolyPoly2D(pa, npa, pb, npb)) {
						overlap = true;
						break;
					}
				}
				if (overlap)
					continue;

				resultRef.push(neighbourRef);
				resultParent.push(curRef);
				stack.push(neighbourNode);
			}
		}

		return new FindLocalNeighbourhoodResult(resultRef, resultParent);
	}

	static SegInternval = class SegInterval {
		ref;
		tmin;
		tmax;

		constructor(ref, tmin, tmax) {
			this.ref = ref;
			this.tmin = tmin;
			this.tmax = tmax;
		}

	};

	insertInterval(ints, tmin, tmax, ref) {
		// Find insertion point.
		let idx = 0;
		while (idx < ints.length) {
			if (tmax <= ints[idx].tmin)
				break;
			idx++;
		}
		// Store
		ints.add(idx, new SegInterval(ref, tmin, tmax));
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
	/// @par
	///
	/// If the @p segmentRefs parameter is provided, then all polygon segments will be returned. 
	/// Otherwise only the wall segments are returned.
	/// 
	/// A segment that is normally a portal will be included in the result set as a 
	/// wall if the @p filter results in the neighbor polygon becoomming impassable.
	/// 
	/// The @p segmentVerts and @p segmentRefs buffers should normally be sized for the 
	/// maximum segments per polygon of the source navigation mesh.
	/// 
	/// Returns the segments for the specified polygon, optionally including portals.
	///  @param[in]		ref				The reference id of the polygon.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	segmentVerts	The segments. [(ax, ay, az, bx, by, bz) * segmentCount]
	///  @param[out]	segmentRefs		The reference ids of each segment's neighbor polygon. 
	///  								Or zero if the segment is a wall. [opt] [(parentRef) * @p segmentCount] 
	///  @param[out]	segmentCount	The number of segments returned.
	///  @param[in]		maxSegments		The maximum number of segments the result arrays can hold.
	/// @returns The status flags for the query.
	getPolyWallSegments(ref, storePortals, filter) {
		let tileAndPoly = this.m_nav.getTileAndPolyByRef(ref);
		let tile = tileAndPoly[0];
		let poly = tileAndPoly[1];

		let segmentRefs = [];
		let segmentVerts = [];
		let ints = new Array(16);

		for (let i = 0, j = poly.vertCount - 1; i < poly.vertCount; j = i++) {
			// Skip non-solid edges.
			ints = [];
			if ((poly.neis[j] & NavMesh.DT_EXT_LINK) != 0) {
				// Tile border.
				for (let k = poly.firstLink; k != NavMesh.DT_NULL_LINK; k = tile.links[k].next) {
					let link = tile.links[k];
					if (link.edge == j) {
						if (link.ref != 0) {
							tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(link.ref);
							let neiTile = tileAndPoly[0];
							let neiPoly = tileAndPoly[1];
							if (filter.passFilter(link.ref, neiTile, neiPoly)) {
								insertInterval(ints, link.bmin, link.bmax, link.ref);
							}
						}
					}
				}
			} else {
				// Internal edge
				let neiRef = 0;
				if (poly.neis[j] != 0) {
					let idx = (poly.neis[j] - 1);
					neiRef = this.or(this.m_nav.getPolyRefBase(tile) , idx);
					if (!filter.passFilter(neiRef, tile, tile.data.polys[idx]))
						neiRef = 0;
				}
				// If the edge leads to another polygon and portals are not stored, skip.
				if (neiRef != 0 && !storePortals)
					continue;

				let vj = poly.verts[j] * 3;
				let vi = poly.verts[i] * 3;
				let seg = new Array(6);
				arraycopy(tile.data.verts, vj, seg, 0, 3);
				arraycopy(tile.data.verts, vi, seg, 3, 3);
				segmentVerts.push(seg);
				segmentRefs.push(neiRef);
				continue;
			}

			// Add sentinels
			insertInterval(ints, -1, 0, 0);
			insertInterval(ints, 255, 256, 0);

			// Store segments.
			let vj = poly.verts[j] * 3;
			let vi = poly.verts[i] * 3;
			for (let k = 1; k < ints.length; ++k) {
				// Portal segment.
				if (storePortals && ints[k].ref != 0) {
					tmin = ints[k].tmin / 255.0;
					tmax = ints[k].tmax / 255.0;
					let seg = new Array(6);
					arraycopy(vLerp4(tile.data.verts, vj, vi, tmin), 0, seg, 0, 3);
					arraycopy(vLerp4(tile.data.verts, vj, vi, tmax), 0, seg, 3, 3);
					segmentVerts.push(seg);
					segmentRefs.push(ints[k].ref);
				}

				// Wall segment.
				let imin = ints[k - 1].tmax;
				let imax = ints[k].tmin;
				if (imin != imax) {
					tmin = imin / 255.0;
					tmax = imax / 255.0;
					let seg = new Array(6);
					arraycopy(vLerp4(tile.data.verts, vj, vi, tmin), 0, seg, 0, 3);
					arraycopy(vLerp4(tile.data.verts, vj, vi, tmax), 0, seg, 3, 3);
					segmentVerts.push(seg);
					segmentRefs.push(0);
				}
			}
		}

		return new GetPolyWallSegmentsResult(segmentVerts, segmentRefs);
	}

	/// @par
	///
	/// @p hitPos is not adjusted using the height detail data.
	///
	/// @p hitDist will equal the search radius if there is no wall within the 
	/// radius. In this case the values of @p hitPos and @p hitNormal are
	/// undefined.
	///
	/// The normal will become unpredicable if @p hitDist is a very small number.
	///
	/// Finds the distance from the specified position to the nearest polygon wall.
	///  @param[in]		startRef		The reference id of the polygon containing @p centerPos.
	///  @param[in]		centerPos		The center of the search circle. [(x, y, z)]
	///  @param[in]		maxRadius		The radius of the search circle.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	hitDist			The distance to the nearest wall from @p centerPos.
	///  @param[out]	hitPos			The nearest position on the wall that was hit. [(x, y, z)]
	///  @param[out]	hitNormal		The normalized ray formed from the wall poPoly to the 
	///  								source point. [(x, y, z)]
	/// @returns The status flags for the query.
	findDistanceToWall(startRef, centerPos, maxRadius,
		filter) {

		// Validate input
		if (startRef == 0 || !this.m_nav.isValidPolyRef(startRef))
			throw new IllegalArgumentException("Invalid start ref");

		this.m_nodePool = [];
		this.m_openList = [];

		let startNode = this.m_nodePool.getNode(startRef);
		DetourCommon.vCopy(startNode.pos, centerPos);
		startNode.pidx = 0;
		startNode.cost = 0;
		startNode.total = 0;
		startNode.id = startRef;
		startNode.flags = Node.DT_NODE_OPEN;
		this.m_openList.push(startNode);

		radiusSqr = DetourCommon.sqr(maxRadius);
		let hitPos = new Array(3);
		let bestvj = null;
		let bestvi = null;
		while (!this.m_openList.length == 0) {
			let bestNode = this.m_openList.pop();
			bestNode.flags &= ~Node.DT_NODE_OPEN;
			bestNode.flags |= Node.DT_NODE_CLOSED;

			// Get let and tile.
			// The API input has been cheked already, skip checking internal data.
			let bestRef = bestNode.id;
			let tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(bestRef);
			let bestTile = tileAndPoly[0];
			let bestPoly = tileAndPoly[1];

			// Get parent let and tile.
			let parentRef = 0;
			let parentTile = null;
			let parentPoly = null;
			if (bestNode.pidx != 0)
				parentRef = this.m_nodePool.getNodeAtIdx(bestNode.pidx).id;
			if (parentRef != 0) {
				tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(parentRef);
				parentTile = tileAndPoly[0];
				parentPoly = tileAndPoly[1];
			}

			// Hit test walls.
			for (let i = 0, j = bestPoly.vertCount - 1; i < bestPoly.vertCount; j = i++) {
				// Skip non-solid edges.
				if ((bestPoly.neis[j] & NavMesh.DT_EXT_LINK) != 0) {
					// Tile border.
					let solid = true;
					for (let k = bestPoly.firstLink; k != NavMesh.DT_NULL_LINK; k = bestTile.links[k].next) {
						let link = bestTile.links[k];
						if (link.edge == j) {
							if (link.ref != 0) {
								tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(link.ref);
								let neiTile = tileAndPoly[0];
								let neiPoly = tileAndPoly[1];
								if (filter.passFilter(link.ref, neiTile, neiPoly))
									solid = false;
							}
							break;
						}
					}
					if (!solid)
						continue;
				} else if (bestPoly.neis[j] != 0) {
					// Internal edge
					let idx = (bestPoly.neis[j] - 1);
					let ref = this.m_nav.getPolyRefBase(bestTile) | idx;
					if (filter.passFilter(ref, bestTile, bestTile.data.polys[idx]))
						continue;
				}

				// Calc distance to the edge.
				let vj = bestPoly.verts[j] * 3;
				let vi = bestPoly.verts[i] * 3;
				let distseg = DetourCommon.distancePtSegSqr2D4(centerPos, bestTile.data.verts, vj, vi);
				distSqr = distseg[0];
				tseg = distseg[1];

				// Edge is too far, skip.
				if (distSqr > radiusSqr)
					continue;

				// Hit wall, update radius.
				radiusSqr = distSqr;
				// Calculate hit pos.
				hitPos[0] = bestTile.data.verts[vj] + (bestTile.data.verts[vi] - bestTile.data.verts[vj]) * tseg;
				hitPos[1] = bestTile.data.verts[vj + 1] + (bestTile.data.verts[vi + 1] - bestTile.data.verts[vj + 1]) * tseg;
				hitPos[2] = bestTile.data.verts[vj + 2] + (bestTile.data.verts[vi + 2] - bestTile.data.verts[vj + 2]) * tseg;
				bestvj = new VectorPtr(bestTile.data.verts, vj);
				bestvi = new VectorPtr(bestTile.data.verts, vi);
			}

			for (let i = bestPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = bestTile.links[i].next) {
				let link = bestTile.links[i];
				let neighbourRef = link.ref;
				// Skip invalid neighbours and do not follow back to parent.
				if (neighbourRef == 0 || neighbourRef == parentRef)
					continue;

				// Expand to neighbour.
				tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
				let neighbourTile = tileAndPoly[0];
				let neighbourPoly = tileAndPoly[1];

				// Skip off-mesh connections.
				if (neighbourPoly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION)
					continue;

				// Calc distance to the edge.
				let va = bestPoly.verts[link.edge] * 3;
				let vb = bestPoly.verts[(link.edge + 1) % bestPoly.vertCount] * 3;
				let distseg = DetourCommon.distancePtSegSqr2D4(centerPos, bestTile.data.verts, va, vb);
				distSqr = distseg[0];
				// If the circle is not touching the next polygon, skip it.
				if (distSqr > radiusSqr)
					continue;

				if (!filter.passFilter(neighbourRef, neighbourTile, neighbourPoly))
					continue;

				let neighbourNode = this.m_nodePool.getNode(neighbourRef);

				if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0)
					continue;

				// Cost
				if (neighbourNode.flags == 0) {
					neighbourNode.pos = this.getEdgeMidPoint6(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly,
						neighbourTile);
				}

				total = bestNode.total + DetourCommon.vDist2(bestNode.pos, neighbourNode.pos);

				// The node is already in open list and the new result is worse, skip.
				if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0 && total >= neighbourNode.total)
					continue;

				neighbourNode.id = neighbourRef;
				neighbourNode.flags = (neighbourNode.flags & ~Node.DT_NODE_CLOSED);
				neighbourNode.pidx = this.m_nodePool.getNodeIdx(bestNode);
				neighbourNode.total = total;

				if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0) {
					this.m_openList.modify(neighbourNode);
				} else {
					neighbourNode.flags |= Node.DT_NODE_OPEN;
					this.m_openList.push(neighbourNode);
				}
			}
		}

		// Calc hit normal.
		let hitNormal = new Array(3);
		if (bestvi != null && bestvj != null) {
			let tangent = DetourCommon.vSub(bestvi, bestvj);
			hitNormal[0] = tangent[2];
			hitNormal[1] = 0;
			hitNormal[2] = -tangent[0];
			vNormalize(hitNormal);
		}
		return new FindDistanceToWallResult(Math.sqrt(radiusSqr), hitPos, hitNormal);
	}

	/// Returns true if the polygon reference is valid and passes the filter restrictions.
	///  @param[in]		ref			The polygon reference to check.
	///  @param[in]		filter		The filter to apply.
	isValidPolyRef(ref, filter) {
		try {
			let tileAndPoly = this.m_nav.getTileAndPolyByRef(ref);
			// If cannot pass filter, assume flags has changed and boundary is invalid.
			if (filter.passFilter(ref, tileAndPoly[0], tileAndPoly[1]))
				return true;
		} catch (e) {
			// If cannot get polygon, assume it does not exists and boundary is invalid.
		}
		return false;
	}

	/// Gets the navigation mesh the query object is using.
	/// @return The navigation mesh the query object is using.
	getAttachedNavMesh() {
		return this.m_nav;
	}

	/*
	/// @par
	///
	/// The closed list is the list of polygons that were fully evaluated during 
	/// the last navigation graph search. (A* or Dijkstra)
	/// 
	/// Returns true if the polygon reference is in the closed list. 
	///  @param[in]		ref		The reference id of the polygon to check.
	/// @returns True if the polygon is in closed list.
	let isInClosedList(let ref)
	{
		if (m_nodePool == null) return false;
		
		let nodes[DT_MAX_STATES_PER_NODE];
		let n= m_nodePool=>findNodes(ref, nodes, DT_MAX_STATES_PER_NODE);
		
		for (let i=0; i<n; i++)
		{
			if (nodes[i]=>flags & DT_NODE_CLOSED)
				return true;
		}		
		
		return false;
	}
		
	*/

	/**
	 * Gets a path from the explored nodes in the previous search.
	 * 
	 * @param endRef
	 *            The reference id of the end polygon.
	 * @returns An ordered list of polygon references representing the path. (Start to end.)
	 * @remarks The result of this function depends on the state of the query object. For that reason it should only be
	 *          used immediately after one of the two Dijkstra searches, findPolysAroundCircle or findPolysAroundShape.
	 */
	getPathFromDijkstraSearch(endRef) {
		if (!this.m_nav.isValidPolyRef(endRef))
			throw new IllegalArgumentException("Invalid end ref");
		let nodes = this.m_nodePool.findNodes(endRef);
		if (nodes.length != 1)
			throw new IllegalArgumentException("Invalid end ref");
		let endNode = nodes[0];
		if ((endNode.flags & DT_NODE_CLOSED) == 0)
			throw new IllegalArgumentException("Invalid end ref");
		return getPathToNode(endNode);
	}

	/**
	 * Gets the path leading to the specified end node.
	 */
	getPathToNode(endNode) {
		let path = [];
		// Reverse the path.
		let curNode = endNode;
		do {
			path.add(0, curNode.id);
			curNode = this.m_nodePool.getNodeAtIdx(curNode.pidx);
		} while (curNode != null);

		return path;
	}
}

export default NavMeshQuery;