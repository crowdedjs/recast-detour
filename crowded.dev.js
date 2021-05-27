(function (global, factory) {
	typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
	typeof define === 'function' && define.amd ? define(['exports'], factory) :
	(global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.crowded = {}));
}(this, (function (exports) { 'use strict';

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

	/// Configuration parameters for a crowd agent.
	/// @ingroup crowd
	class CrowdAgentParams {
		radius = 0;						///<  radius. [Limit: >= 0]
		height = 0;						///<  height. [Limit: > 0]
		maxAcceleration = 0;			///< Maximum allowed acceleration. [Limit: >= 0]
		maxSpeed = 0;						///< Maximum allowed speed. [Limit: >= 0]

		/// Defines how close a collision element must be before it is considered for steering behaviors. [Limits: > 0]
		collisionQueryRange = 0;

		pathOptimizationRange = 0;		///< The path visibility optimization range. [Limit: > 0]

		/// How aggresive the agent manager should be at avoiding collisions with this agent. [Limit: >= 0]
		separationWeight = 0;

		///  agent update flags.
		static DT_CROWD_ANTICIPATE_TURNS = 1;
		static DT_CROWD_OBSTACLE_AVOIDANCE = 2;
		static DT_CROWD_SEPARATION = 4;
		static DT_CROWD_OPTIMIZE_VIS = 8;			///< Use #dtPathCorridor::optimizePathVisibility() to optimize the agent path.
		static DT_CROWD_OPTIMIZE_TOPO = 16;		///< Use dtPathCorridor::optimizePathTopology() to optimize the agent path.

		/// Flags that impact steering behavior. (See: #UpdateFlags)
		updateFlags = 0;

		/// The index of the avoidance configuration to use for the agent. 
		/// [Limits: 0 <= value <= #DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS]
		obstacleAvoidanceType = 0;

		/// The index of the query filter used by this agent.
		queryFilterType = 0;

		/// User defined data attached to the agent.
		userData = null;
	}

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

	/**
	 *  Configuration parameters used to define multi-tile navigation meshes.
	 *  The values are used to allocate space during the initialization of a navigation mesh.
	 *  @see NavMesh
	 */
	class NavMeshParams {
		/** The world space origin of the navigation mesh's tile space. [(x, y, z)] */
	orig = new Array(3); 
		/** The width of each tile. (APoly the x-axis.) */
	 tileWidth = 0;
		/** The height of each tile. (APoly the z-axis.) */
	 tileHeight = 0;
		/** The maximum number of tiles the navigation mesh can contain. */
	maxTiles = 0;
		/** The maximum number of polygons each tile can contain. */
	maxPolys = 0; 
	}

	class IntersectResult {
	  intersects = false;
	  tmin = 0;;
	  tmax = 1;
	  segMin = -1;
	  segMax = -1;
	}

	/*
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

	/**
	 * Wrapper for 3-element pieces (3D vectors) of a bigger  array.
	 *
	 */
	class VectorPtr$1 {

	array;
	 index;

	// constructor( array) {
	// 		this(array, 0);
	// 	}

	constructor( array,  index=0) {
			this.array = array;
			this.index = index;
		}

	 get(offset) {
			return array[index + offset];
		}

	}

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

	class DetourCommon {

		static EPS = 1e-4;

		/// Performs a scaled vector addition. (@p v1 + (@p v2 * @p s))
		/// @param[out] dest The result vector. [(x, y, z)]
		/// @param[_in] v1 The base vector. [(x, y, z)]
		/// @param[_in] v2 The vector to scale and add to @p v1. [(x, y, z)]
		/// @param[_in] s The amount to scale @p v2 by before adding to @p v1.
		static vMad(v1, v2, s) {
			let dest = new Array(3);
			dest[0] = v1[0] + v2[0] * s;
			dest[1] = v1[1] + v2[1] * s;
			dest[2] = v1[2] + v2[2] * s;
			return dest;
		}

		/// Performs a linear interpolation between two vectors. (@p v1 toward @p
		/// v2)
		/// @param[out] dest The result vector. [(x, y, x)]
		/// @param[_in] v1 The starting vector.
		/// @param[_in] v2 The destination vector.
		/// @param[_in] t The interpolation factor. [Limits: 0 <= value <= 1.0]
		static vLerp4(verts, v1, v2, t) {
			let dest = new Array(3);
			dest[0] = verts[v1 + 0] + (verts[v2 + 0] - verts[v1 + 0]) * t;
			dest[1] = verts[v1 + 1] + (verts[v2 + 1] - verts[v1 + 1]) * t;
			dest[2] = verts[v1 + 2] + (verts[v2 + 2] - verts[v1 + 2]) * t;
			return dest;
		}

		static vLerp3(v1, v2, t) {
			let dest = new Array(3);
			dest[0] = v1[0] + (v2[0] - v1[0]) * t;
			dest[1] = v1[1] + (v2[1] - v1[1]) * t;
			dest[2] = v1[2] + (v2[2] - v1[2]) * t;
			return dest;
		}

		static vSub(v1, v2) {
			let dest = new Array(3);
			dest[0] = v1[0] - v2[0];
			dest[1] = v1[1] - v2[1];
			dest[2] = v1[2] - v2[2];
			return dest;
		}

		// static vSub(v1, v2) {
		// 	let dest = new Array(3);
		// 	dest[0] = v1[0] - v2[0];
		// 	dest[1] = v1[1] - v2[1];
		// 	dest[2] = v1[2] - v2[2];
		// 	return dest;
		// }


		static vAdd(v1, v2) {
			let dest = new Array(3);
			dest[0] = v1[0] + v2[0];
			dest[1] = v1[1] + v2[1];
			dest[2] = v1[2] + v2[2];
			return dest;
		}

		static vCopy_return(i) {
			let out = new Array(3);
			out[0] = i[0];
			out[1] = i[1];
			out[2] = i[2];
			return out;
		}

		static vSet(out, a, b, c) {
			out[0] = a;
			out[1] = b;
			out[2] = c;
		}

		// static vCopy(out, i) {
		// 	out[0] = i[0];
		// 	out[1] = i[1];
		// 	out[2] = i[2];
		// }

		//See vCopy_return for the 1 parameter version
		static vCopy(out, _in, i = 0) {
			out[0] = _in[i];
			out[1] = _in[i + 1];
			out[2] = _in[i + 2];
		}

		static vMin(out, _in, i) {
			out[0] = Math.min(out[0], _in[i]);
			out[1] = Math.min(out[1], _in[i + 1]);
			out[2] = Math.min(out[2], _in[i + 2]);
		}

		static vMax(out, _in, i) {
			out[0] = Math.max(out[0], _in[i]);
			out[1] = Math.max(out[1], _in[i + 1]);
			out[2] = Math.max(out[2], _in[i + 2]);
		}

		/// Returns the distance between two points.
		/// @param[_in] v1 A point. [(x, y, z)]
		/// @param[_in] v2 A point. [(x, y, z)]
		/// @return The distance between the two points.
		static vDist2(v1, v2) {
			let dx = v2[0] - v1[0];
			let dy = v2[1] - v1[1];
			let dz = v2[2] - v1[2];
			return Math.sqrt(dx * dx + dy * dy + dz * dz);
		}

		/// Returns the distance between two points.
		/// @param[_in] v1 A point. [(x, y, z)]
		/// @param[_in] v2 A point. [(x, y, z)]
		/// @return The distance between the two points.
		static vDistSqr(v1, v2) {
			let dx = v2[0] - v1[0];
			let dy = v2[1] - v1[1];
			let dz = v2[2] - v1[2];
			return dx * dx + dy * dy + dz * dz;
		}

		static sqr(a) {
			return a * a;
		}

		/// Derives the square of the scalar length of the vector. (len * len)
		/// @param[_in] v The vector. [(x, y, z)]
		/// @return The square of the scalar length of the vector.
		static vLenSqr(v) {
			return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
		}

		static vLen(v) {
			return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
		}

		static vDist3(v1, verts, i) {
			let dx = verts[i] - v1[0];
			let dy = verts[i + 1] - v1[1];
			let dz = verts[i + 2] - v1[2];
			return Math.sqrt(dx * dx + dy * dy + dz * dz);
		}

		static clamp(v, min, max) {
			return Math.max(Math.min(v, max), min);
		}

		static clamp(v, min, max) {
			return Math.max(Math.min(v, max), min);
		}

		/// Derives the distance between the specified points on the xz-plane.
		/// @param[_in] v1 A point. [(x, y, z)]
		/// @param[_in] v2 A point. [(x, y, z)]
		/// @return The distance between the poPoly on the xz-plane.
		///
		/// The vectors are projected onto the xz-plane, so the y-values are
		/// ignored.
		static vDist2D(v1, v2) {
			let dx = v2[0] - v1[0];
			let dz = v2[2] - v1[2];
			return Math.sqrt(dx * dx + dz * dz);
		}

		static vDist2DSqr(v1, v2) {
			let dx = v2[0] - v1[0];
			let dz = v2[2] - v1[2];
			return dx * dx + dz * dz;
		}

		/// Normalizes the vector.
		/// @param[_in,out] v The vector to normalize. [(x, y, z)]
		static vNormalize(v) {
			let d = (1.0 / Math.sqrt(DetourCommon.sqr(v[0]) + DetourCommon.sqr(v[1]) + DetourCommon.sqr(v[2])));
			if (d != 0) {
				v[0] *= d;
				v[1] *= d;
				v[2] *= d;
			}
		}

		static thr = DetourCommon.sqr(1.0 / 16384.0);

		/// Performs a 'sloppy' colocation check of the specified points.
		/// @param[_in] p0 A point. [(x, y, z)]
		/// @param[_in] p1 A point. [(x, y, z)]
		/// @return True if the points are considered to be at the same location.
		///
		/// Basically, this function will return true if the specified points are
		/// close enough to eachother to be considered colocated.
		static vEqual(p0, p1) {
			let d = DetourCommon.vDistSqr(p0, p1);
			return d < DetourCommon. thr;
		}

		/// Derives the dot product of two vectors on the xz-plane. (@p u . @p v)
		/// @param[_in] u A vector [(x, y, z)]
		/// @param[_in] v A vector [(x, y, z)]
		/// @return The dot product on the xz-plane.
		///
		/// The vectors are projected onto the xz-plane, so the y-values are
		/// ignored.
		// static DetourCommon.vDot2D(u, v) {
		// 	return u[0] * v[0] + u[2] * v[2];
		// }

		static vDot2D(u, v, vi=0) {
			return u[0] * v[vi] + u[2] * v[vi + 2];
		}

		/// Derives the xz-plane 2D perp product of the two vectors. (uz*vx - ux*vz)
		/// @param[_in] u The LHV vector [(x, y, z)]
		/// @param[_in] v The RHV vector [(x, y, z)]
		/// @return The dot product on the xz-plane.
		///
		/// The vectors are projected onto the xz-plane, so the y-values are
		/// ignored.
		static vPerp2D(u, v) {
			return u[2] * v[0] - u[0] * v[2];
		}

		/// @}
		/// @name Computational geometry helper functions.
		/// @{

		/// Derives the signed xz-plane area of the triangle ABC, or the
		/// relationship of line AB to poPoly C.
		/// @param[_in] a Vertex A. [(x, y, z)]
		/// @param[_in] b Vertex B. [(x, y, z)]
		/// @param[_in] c Vertex C. [(x, y, z)]
		/// @return The signed xz-plane area of the triangle.
		static triArea2D4(verts, a, b, c) {
			abx = verts[b] - verts[a];
			abz = verts[b + 2] - verts[a + 2];
			acx = verts[c] - verts[a];
			acz = verts[c + 2] - verts[a + 2];
			return acx * abz - abx * acz;
		}

		static triArea2D3(a, b, c) {
			let abx = b[0] - a[0];
			let abz = b[2] - a[2];
			let acx = c[0] - a[0];
			let acz = c[2] - a[2];
			return acx * abz - abx * acz;
		}

		/// Determines if two axis-aligned bounding boxes overlap.
		/// @param[_in] amin Minimum bounds of box A. [(x, y, z)]
		/// @param[_in] amax Maximum bounds of box A. [(x, y, z)]
		/// @param[_in] bmin Minimum bounds of box B. [(x, y, z)]
		/// @param[_in] bmax Maximum bounds of box B. [(x, y, z)]
		/// @return True if the two AABB's overlap.
		/// @see dtOverlapBounds
		static overlapQuantBounds(amin, amax, bmin, bmax) {
			let overlap = true;
			overlap = (amin[0] > bmax[0] || amax[0] < bmin[0]) ? false : overlap;
			overlap = (amin[1] > bmax[1] || amax[1] < bmin[1]) ? false : overlap;
			overlap = (amin[2] > bmax[2] || amax[2] < bmin[2]) ? false : overlap;
			return overlap;
		}

		/// Determines if two axis-aligned bounding boxes overlap.
		/// @param[_in] amin Minimum bounds of box A. [(x, y, z)]
		/// @param[_in] amax Maximum bounds of box A. [(x, y, z)]
		/// @param[_in] bmin Minimum bounds of box B. [(x, y, z)]
		/// @param[_in] bmax Maximum bounds of box B. [(x, y, z)]
		/// @return True if the two AABB's overlap.
		/// @see dtOverlapQuantBounds
		// static overlapBounds(amin, amax, bmin, bmax) {
		// 	let overlap = true;
		// 	overlap = (amin[0] > bmax[0] || amax[0] < bmin[0]) ? false : overlap;
		// 	overlap = (amin[1] > bmax[1] || amax[1] < bmin[1]) ? false : overlap;
		// 	overlap = (amin[2] > bmax[2] || amax[2] < bmin[2]) ? false : overlap;
		// 	return overlap;
		// }

		static distancePtSegSqr2D3(pt, p, q) {
			let pqx = q[0] - p[0];
			let pqz = q[2] - p[2];
			let dx = pt[0] - p[0];
			let dz = pt[2] - p[2];
			let d = pqx * pqx + pqz * pqz;
			let t = pqx * dx + pqz * dz;
			if (d > 0)
				t /= d;
			if (t < 0)
				t = 0;
			else if (t > 1)
				t = 1;
			dx = p[0] + t * pqx - pt[0];
			dz = p[2] + t * pqz - pt[2];
			return [dx * dx + dz * dz, t];
		}

		static closestHeightPointTriangle(p, a, b, c) {
			let v0 = DetourCommon.vSub(c, a);
			let v1 = DetourCommon.vSub(b, a);
			let v2 = DetourCommon.vSub(p, a);

			let dot00 = DetourCommon.vDot2D(v0, v0);
			let dot01 = DetourCommon.vDot2D(v0, v1);
			let dot02 = DetourCommon.vDot2D(v0, v2);
			let dot11 = DetourCommon.vDot2D(v1, v1);
			let dot12 = DetourCommon.vDot2D(v1, v2);

			// Compute barycentric coordinates
			let invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
			let u = (dot11 * dot02 - dot01 * dot12) * invDenom;
			let v = (dot00 * dot12 - dot01 * dot02) * invDenom;

			// The (sloppy) epsilon is needed to allow to get height of points which
			// are interpolated aPoly the edges of the triangles.

			// If point lies inside the triangle, return interpolated ycoord.
			if (u >= -DetourCommon.EPS && v >= -DetourCommon.EPS && (u + v) <= 1 + DetourCommon.EPS) {
				let h = a[1] + v0[1] * u + v1[1] * v;
				return [true, h];
			}

			return [false, null];
		}

		/// @par
		///
		/// All points are projected onto the xz-plane, so the y-values are ignored.
		static pointInPolygon(pt, verts, nverts) {
			let c = false;
			for (let i = 0, j = nverts - 1; i < nverts; j = i++) {
				let vi = i * 3;
				let vj = j * 3;
				if (((verts[vi + 2] > pt[2]) != (verts[vj + 2] > pt[2])) && (pt[0] < (verts[vj + 0] - verts[vi + 0])
					* (pt[2] - verts[vi + 2]) / (verts[vj + 2] - verts[vi + 2]) + verts[vi + 0]))
					c = !c;
			}
			return c;
		}

		static distancePtPolyEdgesSqr(pt, verts, nverts, ed, et) {
			let c = false;
			for (let i = 0, j = nverts - 1; i < nverts; j = i++) {
				let vi = i * 3;
				let vj = j * 3;
				if (((verts[vi + 2] > pt[2]) != (verts[vj + 2] > pt[2])) && (pt[0] < (verts[vj + 0] - verts[vi + 0])
					* (pt[2] - verts[vi + 2]) / (verts[vj + 2] - verts[vi + 2]) + verts[vi + 0]))
					c = !c;
				let edet = DetourCommon.distancePtSegSqr2D4(pt, verts, vj, vi);
				ed[j] = edet[0];
				et[j] = edet[1];
			}
			return c;
		}

		static projectPoly(axis, poly, npoly) {
			let rmin;
			let rmax;
			rmin = rmax = DetourCommon.vDot2D(axis, poly, 0);
			for (let i = 1; i < npoly; ++i) {
				let d = DetourCommon.vDot2D(axis, poly, i * 3);
				rmin = Math.min(rmin, d);
				rmax = Math.max(rmax, d);
			}
			return [rmin, rmax];
		}

		static overlapRange(amin, amax, bmin, bmax, eps) {
			return ((amin + eps) > bmax || (amax - eps) < bmin) ? false : true;
		}

		static eps = 1e-4;

		/// @par
		///
		/// All vertices are projected onto the xz-plane, so the y-values are ignored.
		static overlapPolyPoly2D(polya, npolya, polyb, npolyb) {

			for (let i = 0, j = npolya - 1; i < npolya; j = i++) {
				let va = j * 3;
				let vb = i * 3;

				let n = [polya[vb + 2] - polya[va + 2], 0, -(polya[vb + 0] - polya[va + 0])];

				let aminmax = DetourCommon.projectPoly(n, polya, npolya);
				let bminmax = DetourCommon.projectPoly(n, polyb, npolyb);
				if (!DetourCommon.overlapRange(aminmax[0], aminmax[1], bminmax[0], bminmax[1], DetourCommon.eps)) {
					// Found separating axis
					return false;
				}
			}
			for (let i = 0, j = npolyb - 1; i < npolyb; j = i++) {
				let va = j * 3;
				let vb = i * 3;

				let n = [polyb[vb + 2] - polyb[va + 2], 0, -(polyb[vb + 0] - polyb[va + 0])];

				let aminmax = DetourCommon.projectPoly(n, polya, npolya);
				let bminmax = DetourCommon.projectPoly(n, polyb, npolyb);
				if (!DetourCommon.overlapRange(aminmax[0], aminmax[1], bminmax[0], bminmax[1], DetourCommon.eps)) {
					// Found separating axis
					return false;
				}
			}
			return true;
		}

		// Returns a random poPoly _in a convex polygon.
		// Adapted from Graphics Gems article.
		static randomPointInConvexPoly(pts, npts, areas, s, t) {
			// Calc triangle araes
			areasum = 0.0;
			for (let i = 2; i < npts; i++) {
				areas[i] = DetourCommon.triArea2D4(pts, 0, (i - 1) * 3, i * 3);
				areasum += Math.max(0.001, areas[i]);
			}
			// Find sub triangle weighted by area.
			thr = s * areasum;
			acc = 0.0;
			u = 1.0;
			let tri = npts - 1;
			for (let i = 2; i < npts; i++) {
				dacc = areas[i];
				if (thr >= acc && thr < (acc + dacc)) {
					u = (thr - acc) / dacc;
					tri = i;
					break;
				}
				acc += dacc;
			}

			v = Math.sqrt(t);

			a = 1 - v;
			b = (1 - u) * v;
			c = u * v;
			let pa = 0;
			let pb = (tri - 1) * 3;
			let pc = tri * 3;

			return [a * pts[pa] + b * pts[pb] + c * pts[pc],
			a * pts[pa + 1] + b * pts[pb + 1] + c * pts[pc + 1],
			a * pts[pa + 2] + b * pts[pb + 2] + c * pts[pc + 2]];
		}

		static nextPow2(v) {
			v--;
			v |= v >> 1;
			v |= v >> 2;
			v |= v >> 4;
			v |= v >> 8;
			v |= v >> 16;
			v++;
			return v;
		}

		static ilog2(v) {
			let r;
			let shift;
			r = (v > 0xffff ? 1 : 0) << 4;
			v >>= r;
			shift = (v > 0xff ? 1 : 0) << 3;
			v >>= shift;
			r |= shift;
			shift = (v > 0xf ? 1 : 0) << 2;
			v >>= shift;
			r |= shift;
			shift = (v > 0x3 ? 1 : 0) << 1;
			v >>= shift;
			r |= shift;
			r |= (v >> 1);
			return r;
		}

		

		static intersectSegmentPoly2D(p0, p1, verts, nverts) {

			let result = new IntersectResult();
			let EPS = 0.00000001;
			let dir = DetourCommon.vSub(p1, p0);

			let p0v = new VectorPtr$1(p0);
			for (let i = 0, j = nverts - 1; i < nverts; j = i++) {
				let vpj = new VectorPtr$1(verts, j * 3);
				let edge = DetourCommon.vSub(new VectorPtr$1(verts, i * 3), vpj);
				let diff = DetourCommon.vSub(p0v, vpj);
				let n = DetourCommon.vPerp2D(edge, diff);
				let d = DetourCommon.vPerp2D(dir, edge);
				if (Math.abs(d) < EPS) {
					// S is nearly parallel to this edge
					if (n < 0)
						return result;
					else
						continue;
				}
				let t = n / d;
				if (d < 0) {
					// segment S is entering across this edge
					if (t > result.tmin) {
						result.tmin = t;
						result.segMin = j;
						// S enters after leaving polygon
						if (result.tmin > result.tmax)
							return result;
					}
				} else {
					// segment S is leaving across this edge
					if (t < result.tmax) {
						result.tmax = t;
						result.segMax = j;
						// S leaves before entering polygon
						if (result.tmax < result.tmin)
							return result;
					}
				}
			}
			result.intersects = true;
			return result;
		}

		static distancePtSegSqr2D4(pt, verts, p, q) {
			let pqx = verts[q + 0] - verts[p + 0];
			let pqz = verts[q + 2] - verts[p + 2];
			let dx = pt[0] - verts[p + 0];
			let dz = pt[2] - verts[p + 2];
			let d = pqx * pqx + pqz * pqz;
			let t = pqx * dx + pqz * dz;
			if (d > 0)
				t /= d;
			if (t < 0)
				t = 0;
			else if (t > 1)
				t = 1;
			dx = verts[p + 0] + t * pqx - pt[0];
			dz = verts[p + 2] + t * pqz - pt[2];
			return [dx * dx + dz * dz, t];
		}

		static oppositeTile(side) {
			return (side + 4) & 0x7;
		}

		static vperpXZ(a, b) {
			return a[0] * b[2] - a[2] * b[0];
		}

		static intersectSegSeg2D(ap, aq, bp, bq) {
			let u = DetourCommon.vSub(aq, ap);
			let v = DetourCommon.vSub(bq, bp);
			let w = DetourCommon.vSub(ap, bp);
			d = DetourCommon.vperpXZ(u, v);
			if (Math.abs(d) < 1e-6)
				return [false, 0, 0];
			s = DetourCommon.vperpXZ(v, w) / d;
			t = DetourCommon.vperpXZ(u, w) / d;
			return [true, s, t];
		}

		static vScale(_in, scale) {
			let out = new Array(3);
			out[0] = _in[0] * scale;
			out[1] = _in[1] * scale;
			out[2] = _in[2] * scale;
			return out;
		}

	}

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

	/**
	 * Defines a navigation mesh tile.
	 */
	class MeshTile {
		index = 0;
		/** Counter describing modifications to the tile. */
		salt = 0;
		/** The tile data. */
		data = null;
		/** The tile links. */
		links = [];
		/** Index to the next free link. */
		linksFreeList = NavMesh.DT_NULL_LINK;
		/** Tile flags. (See: #dtTileFlags) */
		flags = 0;
		/** The next free tile, or the next tile in the spatial grid. */
		next = null;
		constructor(index) {
			this.index = index;
		}

	}

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

	/** Defines a polyogn within a dtPoly object. */
	class Poly {

		index;
		/** The polygon is a standard convex polygon that is part of the surface of the mesh. */
		static DT_POLYTYPE_GROUND = 0;
		/** The polygon is an off-mesh connection consisting of two vertices. */
		static DT_POLYTYPE_OFFMESH_CONNECTION = 1;
		/** Index to first link in linked list. (Or #DT_NULL_LINK if there is no link.) */
		firstLink = 0;
		/** The indices of the polygon's vertices. The actual vertices are located in MeshTile::verts. */
		verts = [];
		/** Packed data representing neighbor polygons references and flags for each edge. */
		neis = [];
		/** The user defined polygon flags. */
		flags = 0;
		/** The number of vertices in the polygon. */
		vertCount = 0;
		/**
		 * The bit packed area id and polygon type.
		 * 
		 * @note Use the structure's set and get methods to access this value.
		 */
		areaAndtype;

		constructor(index, maxVertsPerPoly) {
			this.index = index;
			this.firstLink = NavMesh.DT_NULL_LINK;
			this.verts = new Array(maxVertsPerPoly);
			this.neis = new Array(maxVertsPerPoly);
			for(let i = 0; i < this.verts.length; i++){
				this.verts[i] = 0;
			}
			for(let i = 0; i < this.neis.length; i++){
				this.neis[i] = 0;
			}
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


		/** Sets the user defined area id. [Limit: < #DT_MAX_AREAS] */
		setArea(a) {
			this.areaAndtype = this.or(this.and(this.areaAndtype , 0xc0) , this.and(a , 0x3f));
		}

		/** Sets the polygon type. (See: #dtPolyTypes.) */
		setType(t) {
			this.areaAndtype = this.or(this.and(this.areaAndtype , 0x3f) | (t << 6));
		}

		/** Gets the user defined area id. */
		getArea() {
			return this.areaAndtype & 0x3;
		}

		/** Gets the polygon type. (See: #dtPolyTypes) */
		getType() {
			return this.areaAndtype >> 6;
		}

	}

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

	/**
	 * Defines a link between polygons.
	 * 
	 * @note This structure is rarely if ever used by the end user.
	 * @see MeshTile
	 */
	class Link {
		/** Neighbour reference. (The neighbor that is linked to.) */
		 ref = 0;
		/** Index of the next link. */
		next = 0;
		/** Index of the polygon edge that owns this link. */
		 edge = 0;
		/** If a boundary link, defines on which side the link is. */
		 side = 0;
		/** If a boundary link, defines the minimum sub-edge area. */
		bmin = 0;
		/** If a boundary link, defines the maximum sub-edge area. */
		bmax = 0;

	}

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

	class ClosestPointOnPolyResult {

		posOverPoly = false;
		closest = [];

		constructor(posOverPoly, closest) {
			this.posOverPoly = posOverPoly;
			this.closest = closest;
		}

		/** Returns true if the position is over the polygon. */
		isPosOverPoly() {
			return this.posOverPoly;
		}

		/** Returns the closest poPoly on the polygon. [(x, y, z)] */
		getClosest() {
			return this.closest;
		}

	}

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

	class FindNearestPolyResult {
		nearestRef = 0;
		nearestPos = [];

		constructor(nearestRef, nearestPos) {
			this.nearestRef = nearestRef;
			this.nearestPos = nearestPos;
		}

		/** Returns the reference id of the nearest polygon. */
		getNearestRef() {
			return this.nearestRef;
		}

		/** Returns the nearest poPoly on the polygon. [opt] [(x, y, z)] */
		getNearestPos() {
			return this.nearestPos;
		}

	}

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

	function arraycopy(one, oneStart, two, twoStart, len) {
		for (let i = 0; i < len; i++) {
			two[twoStart + i] = one[oneStart + i];
		}
	}

	class NavMesh {

		static DT_SALT_BITS = 16;
		static DT_TILE_BITS = 28;
		static DT_POLY_BITS = 20;

		/// A flag that indicates that an entity links to an external entity.
		/// (E.g. A polygon edge is a portal that links to another polygon.)
		static DT_EXT_LINK = 0x8000;

		/// A value that indicates the entity does not link to anything.
		// static DT_NULL_LINK = 0xffffffff;
		static DT_NULL_LINK = -1;

		/// A flag that indicates that an off-mesh connection can be traversed in
		/// both directions. (Is bidirectional.)
		static DT_OFFMESH_CON_BIDIR = 1;

		/// The maximum number of user defined area ids.
		static DT_MAX_AREAS = 64;

		/// Limit raycasting during any angle pahfinding
		/// The limit is given as a multiple of the character radius
		static DT_RAY_CAST_LIMIT_PROPORTIONS = 50.0;

		m_params = null/// < Current initialization params. TODO: do not store this info twice.
		m_orig = []; /// < Origin of the tile (0,0)
		//  m_orig[3]; ///< Origin of the tile (0,0)
		m_tileWidth = 0;
		m_tileHeight = 0; /// < Dimensions of each tile.
		m_maxTiles = 0; /// < Max number of tiles.
		m_tileLutSize = 0; /// < Tile hash lookup size (must be pot).
		m_tileLutMask = 0; /// < Tile hash lookup mask.
		m_posLookup = []; /// < Tile hash lookup.
		m_nextFree = null; /// < Freelist of tiles.
		m_tiles = []; /// < List of tiles.
		/** The maximum number of vertices per navigation polygon. */
		m_maxVertPerPoly = 0;
		m_tileCount = 0;

		/**
		 * The maximum number of tiles supported by the navigation mesh.
		 * 
		 * @return The maximum number of tiles supported by the navigation mesh.
		 */
		getMaxTiles() {
			return this.m_maxTiles;
		}

		/**
		 * Returns tile in the tile array.
		 */
		getTile(i) {
			return this.m_tiles[i];
		}

		/**
		 * Gets the polygon reference for the tile's base polygon.
		 * 
		 * @param tile
		 *            The tile.
		 * @return The polygon reference for the base polygon in the specified tile.
		 */
		getPolyRefBase(tile) {
			if (tile == null)
				return 0;
			let it = tile.index;
			return NavMesh.encodePolyId(tile.salt, it, 0);
		}

		/**
		 * Derives a standard polygon reference.
		 * 
		 * @note This function is generally meant for internal use only.
		 * @param salt
		 *            The tile's salt value.
		 * @param it
		 *            The index of the tile.
		 * @param ip
		 *            The index of the polygon within the tile.
		 * @return encoded polygon reference
		 */
		//https://stackoverflow.com/a/337572/10047920
		static lshift(num, bits) {
			return num * Math.pow(2, bits);
		}
		static rshift(num, bits) {
			return Math.floor(num / Math.pow(2, bits));
		}
		//https://stackoverflow.com/a/43666199/10047920
		static and(v1, v2) {
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
		static or(v1, v2) {
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
		static encodePolyId(salt, it, ip) {
			let a = (NavMesh.lshift(salt, NavMesh.DT_POLY_BITS + NavMesh.DT_TILE_BITS));
			let b = (NavMesh.lshift(it, NavMesh.DT_POLY_BITS));
			return NavMesh.or(NavMesh.or(a, b), ip);
		}

		/// Decodes a standard polygon reference.
		/// @note This function is generally meant for internal use only.
		/// @param[in] ref The polygon reference to decode.
		/// @param[out] salt The tile's salt value.
		/// @param[out] it The index of the tile.
		/// @param[out] ip The index of the polygon within the tile.
		/// @see #encodePolyId
		static decodePolyId(ref) {
			let salt;
			let it;
			let ip;
			let saltMask = NavMesh.lshift(1, NavMesh.DT_SALT_BITS) - 1;
			let tileMask = NavMesh.lshift(1, NavMesh.DT_TILE_BITS) - 1;
			let polyMask = NavMesh.lshift(1, NavMesh.DT_POLY_BITS) - 1;
			salt = Math.floor(NavMesh.and(NavMesh.rshift(ref, (NavMesh.DT_POLY_BITS + NavMesh.DT_TILE_BITS)), saltMask));
			it = Math.floor(NavMesh.and(NavMesh.rshift(ref, NavMesh.DT_POLY_BITS), tileMask));
			ip = Math.floor(NavMesh.and(ref, polyMask));
			return [salt, it, ip];
		}

		/// Extracts a tile's salt value from the specified polygon reference.
		/// @note This function is generally meant for internal use only.
		/// @param[in] ref The polygon reference.
		/// @see #encodePolyId
		static decodePolyIdSalt(ref) {
			let saltMask = (1 << NavMesh.DT_SALT_BITS) - 1;
			return Math.floor((ref >> (NavMesh.DT_POLY_BITS + NavMesh.DT_TILE_BITS)) & saltMask);
		}

		/// Extracts the tile's index from the specified polygon reference.
		/// @note This function is generally meant for internal use only.
		/// @param[in] ref The polygon reference.
		/// @see #encodePolyId
		static decodePolyIdTile(ref) {
			let tileMask = (1 << NavMesh.DT_TILE_BITS) - 1;
			return Math.floor((ref >> NavMesh.DT_POLY_BITS) & tileMask);
		}

		/// Extracts the polygon's index (within its tile) from the specified
		/// polygon reference.
		/// @note This function is generally meant for internal use only.
		/// @param[in] ref The polygon reference.
		/// @see #encodePolyId
		static decodePolyIdPoly(ref) {
			let polyMask = (1 << NavMesh.DT_POLY_BITS) - 1;
			return Math.floor(ref & polyMask);
		}

		allocLink(tile) {
			if (tile.linksFreeList == NavMesh.DT_NULL_LINK) {
				let link = new Link();
				link.next = NavMesh.DT_NULL_LINK;
				tile.links.push(link);
				return tile.links.length - 1;
			}
			let link = tile.linksFreeList;
			tile.linksFreeList = tile.links[link].next;
			return link;
		}

		freeLink(tile, link) {
			tile.links[link].next = tile.linksFreeList;
			tile.linksFreeList = link;
		}

		/**
		 * Calculates the tile grid location for the specified world position.
		 * 
		 * @param pos
		 *            The world position for the query. [(x, y, z)]
		 * @return 2-element let array with (tx,ty) tile location
		 */
		calcTileLoc(pos) {
			let tx = Math.floor((pos[0] - this.m_orig[0]) / this.m_tileWidth);
			let ty = Math.floor((pos[2] - this.m_orig[2]) / this.m_tileHeight);
			return [tx, ty];
		}

		getTileAndPolyByRef(ref) {
			if (ref == 0) {
				throw new IllegalArgumentException("ref = 0");
			}
			let saltitip = NavMesh.decodePolyId(ref);
			let salt = saltitip[0];
			let it = saltitip[1];
			let ip = saltitip[2];
			if (it >= this.m_maxTiles)
				throw new IllegalArgumentException("tile > m_maxTiles");
			if (this.m_tiles[it].salt != salt || this.m_tiles[it].data.header == null)
				throw new IllegalArgumentException("Invalid salt or header");
			if (ip >= this.m_tiles[it].data.header.polyCount)
				throw new IllegalArgumentException("poly > polyCount");
			return [this.m_tiles[it], this.m_tiles[it].data.polys[ip]];
		}

		/// @par
		///
		/// @warning Only use this function if it is known that the provided polygon
		/// reference is valid. This function is faster than #getTileAndPolyByRef,
		/// but
		/// it does not validate the reference.
		getTileAndPolyByRefUnsafe(ref) {
			let saltitip = NavMesh.decodePolyId(ref);
			let it = saltitip[1];
			let ip = saltitip[2];
			return [this.m_tiles[it], this.m_tiles[it].data.polys[ip]];
		}

		isValidPolyRef(ref) {
			if (ref == 0)
				return false;
			let saltitip = NavMesh.decodePolyId(ref);
			let salt = saltitip[0];
			let it = saltitip[1];
			let ip = saltitip[2];
			if (it >= this.m_maxTiles)
				return false;
			if (this.m_tiles[it].salt != salt || this.m_tiles[it].data == null)
				return false;
			if (ip >= this.m_tiles[it].data.header.polyCount)
				return false;
			return true;
		}

		getParams() {
			return this.m_params;
		}

		constructor(one, maxVertsPerPoly, flags) {

			if (flags || flags == 0) {
				this._constructor(NavMesh.getNavMeshParams(one), maxVertsPerPoly);
				this.addTile(one, flags, 0);
			}
			else {
				this._constructor(one, maxVertsPerPoly);
			}
		}

		_constructor(params, maxVertsPerPoly) {
			this.m_params = params;
			this.m_orig = params.orig;
			this.m_tileWidth = params.tileWidth;
			this.m_tileHeight = params.tileHeight;
			// Init tiles
			this.m_maxTiles = params.maxTiles;
			this.m_maxVertPerPoly = maxVertsPerPoly;
			let lutsize = DetourCommon.nextPow2(params.maxTiles / 4);
			if (lutsize == 0)
				lutsize = 1;
			this.m_tileLutSize = lutsize;
			this.m_tileLutMask = this.m_tileLutSize - 1;
			this.m_tiles = new Array(this.m_maxTiles);
			this.m_posLookup = new Array(this.m_tileLutSize);
			this.m_nextFree = null;
			for (let i = this.m_maxTiles - 1; i >= 0; --i) {
				this.m_tiles[i] = new MeshTile(i);
				this.m_tiles[i].salt = 1;
				this.m_tiles[i].next = this.m_nextFree;
				this.m_nextFree = this.m_tiles[i];
			}

		}

		static getNavMeshParams(data) {
			let params = new NavMeshParams();
			DetourCommon.vCopy(params.orig, data.header.bmin);
			params.tileWidth = data.header.bmax[0] - data.header.bmin[0];
			params.tileHeight = data.header.bmax[2] - data.header.bmin[2];
			params.maxTiles = 1;
			params.maxPolys = data.header.polyCount;
			return params;
		}

		// TODO: These methods are duplicates from dtNavMeshQuery, but are needed
		// for off-mesh connection finding.

		queryPolygonsInTile(tile, qmin, qmax) {
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
				bmin[0] = NavMesh.and(Math.floor((qfac * minx)), 0xfffe);
				bmin[1] = NavMesh.and(Math.floor((qfac * miny)), 0xfffe);
				bmin[2] = NavMesh.and(Math.floor((qfac * minz)), 0xfffe);
				bmax[0] = NavMesh.or(Math.floor((qfac * maxx + 1)), 1);
				bmax[1] = NavMesh.or(Math.floor((qfac * maxy + 1)), 1);
				bmax[2] = NavMesh.or(Math.floor((qfac * maxz + 1)), 1);

				// Traverse tree
				let base = this.getPolyRefBase(tile);
				let end = tile.data.header.bvNodeCount;
				while (nodeIndex < end) {
					let node = tile.data.bvTree[nodeIndex];
					let overlap = DetourCommon.overlapQuantBounds(bmin, bmax, node.bmin, node.bmax);
					let isLeafNode = node.i >= 0;

					if (isLeafNode && overlap) {
						polys.push(NavMesh.or(base, node.i));
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
				bmin = [null, null, null];
				bmax = [null, null, null];
				let base = this.getPolyRefBase(tile);
				for (let i = 0; i < tile.data.header.polyCount; ++i) {
					let p = tile.data.polys[i];
					// Do not return off-mesh connection polygons.
					if (p.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION)
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
						polys.push(NavMesh.or(base, i));
					}
				}
				return polys;
			}
		}

		/// Adds a tile to the navigation mesh.
		/// @param[in] data Data for the new tile mesh. (See: #dtCreateNavMeshData)
		/// @param[in] dataSize Data size of the new tile mesh.
		/// @param[in] flags Tile flags. (See: #dtTileFlags)
		/// @param[in] lastRef The desired reference for the tile. (When reloading a
		/// tile.) [opt] [Default: 0]
		/// @param[out] result The tile reference. (If the tile was succesfully
		/// added.) [opt]
		/// @return The status flags for the operation.
		/// @par
		///
		/// The add operation will fail if the data is in the wrong format, the
		/// allocated tile
		/// space is full, or there is a tile already at the specified reference.
		///
		/// The lastRef parameter is used to restore a tile with the same tile
		/// reference it had previously used. In this case the #dtPolyRef's for the
		/// tile will be restored to the same values they were before the tile was
		/// removed.
		///
		/// The nav mesh assumes exclusive access to the data passed and will make
		/// changes to the dynamic portion of the data. For that reason the data
		/// should not be reused in other nav meshes until the tile has been successfully
		/// removed from this nav mesh.
		///
		/// @see dtCreateNavMeshData, #removeTile
		addTile(data, flags, lastRef) {
			// Make sure the data is in right format.
			let header = data.header;

			// Make sure the location is free.
			if (this.getTileAt(header.x, header.y, header.layer) != null)
				throw new RuntimeException("Tile already exists");

			// Allocate a tile.
			let tile = null;
			if (lastRef == 0) {
				if (this.m_nextFree != null) {
					tile = this.m_nextFree;
					this.m_nextFree = tile.next;
					tile.next = null;
					this.m_tileCount++;
				}
			} else {
				// Try to relocate the tile to specific index with same salt.
				let tileIndex = decodePolyIdTile(lastRef);
				if (tileIndex >= this.m_maxTiles)
					throw new RuntimeException("Tile index too high");
				// Try to find the specific tile id from the free list.
				let target = this.m_tiles[tileIndex];
				let prev = null;
				tile = this.m_nextFree;
				while (tile != null && tile != target) {
					prev = tile;
					tile = tile.next;
				}
				// Could not find the correct location.
				if (tile != target)
					throw new RuntimeException("Could not find tile");
				// Remove from freelist
				if (prev == null)
					this.m_nextFree = tile.next;
				else
					prev.next = tile.next;

				// Restore salt.
				tile.salt = decodePolyIdSalt(lastRef);
			}

			// Make sure we could allocate a tile.
			if (tile == null)
				throw new RuntimeException("Could not allocate a tile");

			tile.data = data;
			tile.flags = flags;
			tile.links = [];

			// Insert tile into the position lut.
			let h = NavMesh.computeTileHash(header.x, header.y, this.m_tileLutMask);
			tile.next = this.m_posLookup[h];
			this.m_posLookup[h] = tile;

			// Patch header pointers.

			// If there are no items in the bvtree, reset the tree pointer.
			if (tile.data.bvTree != null && tile.data.bvTree.length == 0)
				tile.data.bvTree = null;

			// Init tile.

			this.connectIntLinks(tile);
			// Base off-mesh connections to their starting polygons and connect connections inside the tile.
			this.baseOffMeshLinks(tile);
			this.connectExtOffMeshLinks(tile, tile, -1);

			// Connect with layers in current tile.
			let neis = this.getTilesAt(header.x, header.y);
			for (let j = 0; j < neis.length; ++j) {
				if (neis[j] == tile) {
					continue;
				}
				connectExtLinks(tile, neis[j], -1);
				connectExtLinks(neis[j], tile, -1);
				this.connectExtOffMeshLinks(tile, neis[j], -1);
				this.connectExtOffMeshLinks(neis[j], tile, -1);
			}

			// Connect with neighbour tiles.
			for (let i = 0; i < 8; ++i) {
				neis = this.getNeighbourTilesAt(header.x, header.y, i);
				for (let j = 0; j < neis.length; ++j) {
					connectExtLinks(tile, neis[j], i);
					connectExtLinks(neis[j], tile, DetourCommon.oppositeTile(i));
					this.connectExtOffMeshLinks(tile, neis[j], i);
					this.connectExtOffMeshLinks(neis[j], tile, DetourCommon.oppositeTile(i));
				}
			}

			return this.getTileRef(tile);
		}

		/// Removes the specified tile from the navigation mesh.
		/// @param[in] ref The reference of the tile to remove.
		/// @param[out] data Data associated with deleted tile.
		/// @param[out] dataSize Size of the data associated with deleted tile.
		/// @return The status flags for the operation.
		// dtStatus removeTile(dtTileRef ref, char** data, int* dataSize);
		/// @par
		///
		/// This function returns the data for the tile so that, if desired,
		/// it can be added back to the navigation mesh at a later point.
		///
		/// @see #addTile
		removeTile(ref) {
			if (ref == 0) {
				return null;
			}
			let tileIndex = decodePolyIdTile(ref);
			let tileSalt = decodePolyIdSalt(ref);
			if (tileIndex >= this.m_maxTiles)
				throw new RuntimeException("Invalid tile index");
			let tile = this.m_tiles[tileIndex];
			if (tile.salt != tileSalt)
				throw new RuntimeException("Invalid tile salt");

			// Remove tile from hash lookup.
			let h = NavMesh.computeTileHash(tile.data.header.x, tile.data.header.y, this.m_tileLutMask);
			let prev = null;
			let cur = this.m_posLookup[h];
			while (cur != null) {
				if (cur == tile) {
					if (prev != null)
						prev.next = cur.next;
					else
						this.m_posLookup[h] = cur.next;
					break;
				}
				prev = cur;
				cur = cur.next;
			}

			// Remove connections to neighbour tiles.
			// Create connections with neighbour tiles.

			// Disconnect from other layers in current tile.
			nneis = this.getTilesAt(tile.data.header.x, tile.data.header.y);
			for (let j of nneis) {
				if (j == tile)
					continue;
				unconnectLinks(j, tile);
			}

			// Disconnect from neighbour tiles.
			for (let i = 0; i < 8; ++i) {
				nneis = this.getNeighbourTilesAt(tile.data.header.x, tile.data.header.y, i);
				for (let j of nneis)
					unconnectLinks(j, tile);
			}
			let data = tile.data;
			// Reset tile.
			tile.data = null;

			tile.flags = 0;
			tile.links = [];

			// Update salt, salt should never be zero.
			tile.salt = (tile.salt + 1) & ((1 << NavMesh.DT_SALT_BITS) - 1);
			if (tile.salt == 0)
				tile.salt++;

			// Add to free list.
			tile.next = this.m_nextFree;
			this.m_nextFree = tile;
			this.m_tileCount--;
			return data;
		}

		/// Builds internal polygons links for a tile.
		connectIntLinks(tile) {
			if (tile == null)
				return;

			let base = this.getPolyRefBase(tile);

			for (let i = 0; i < tile.data.header.polyCount; ++i) {
				let poly = tile.data.polys[i];
				poly.firstLink = NavMesh.DT_NULL_LINK;

				if (poly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION)
					continue;

				// Build edge links backwards so that the links will be
				// in the linked list from lowest index to highest.
				for (let j = poly.vertCount - 1; j >= 0; --j) {
					// Skip hard and non-internal edges.
					if (poly.neis[j] == 0 || (poly.neis[j] & NavMesh.DT_EXT_LINK) != 0)
						continue;

					let idx = this.allocLink(tile);
					let link = tile.links[idx];
					link.ref = NavMesh.or(base, (poly.neis[j] - 1));
					link.edge = j;
					link.side = 0xff;
					link.bmin = link.bmax = 0;
					// Add to linked list.
					link.next = poly.firstLink;
					poly.firstLink = idx;
				}
			}
		}
		unconnectLinks(tile, target) {
			if (tile == null || target == null)
				return;

			let targetNum = decodePolyIdTile(this.getTileRef(target));

			for (let i = 0; i < tile.data.header.polyCount; ++i) {
				let poly = tile.data.polys[i];
				let j = poly.firstLink;
				let pj = NavMesh.DT_NULL_LINK;
				while (j != NavMesh.DT_NULL_LINK) {
					if (decodePolyIdTile(tile.links[j].ref) == targetNum) {
						// Remove link.
						let nj = tile.links[j].next;
						if (pj == NavMesh.DT_NULL_LINK)
							poly.firstLink = nj;
						else
							tile.links[pj].next = nj;
						freeLink(tile, j);
						j = nj;
					} else {
						// Advance
						pj = j;
						j = tile.links[j].next;
					}
				}
			}
		}

		connectExtLinks(tile, target, side) {
			if (tile == null)
				return;

			// Connect border links.
			for (let i = 0; i < tile.data.header.polyCount; ++i) {
				let poly = tile.data.polys[i];

				// Create new links.
				// short m = NavMesh.DT_EXT_LINK | (short)side;

				let nv = poly.vertCount;
				for (let j = 0; j < nv; ++j) {
					// Skip non-portal edges.
					if ((poly.neis[j] & NavMesh.DT_EXT_LINK) == 0)
						continue;

					let dir = poly.neis[j] & 0xff;
					if (side != -1 && dir != side)
						continue;

					// Create new links
					let va = poly.verts[j] * 3;
					let vb = poly.verts[(j + 1) % nv] * 3;
					connectedPolys = findConnectingPolys(tile.data.verts, va, vb, target,
						DetourCommon.oppositeTile(dir), 4);
					nei = connectedPolys[0];
					neia = connectedPolys[1];
					let nnei = connectedPolys.third;
					for (let k = 0; k < nnei; ++k) {
						let idx = this.allocLink(tile);
						let link = tile.links[idx];
						link.ref = nei[k];
						link.edge = j;
						link.side = dir;

						link.next = poly.firstLink;
						poly.firstLink = idx;

						// Compress portal limits to a byte value.
						if (dir == 0 || dir == 4) {
							tmin = (neia[k * 2 + 0] - tile.data.verts[va + 2])
								/ (tile.data.verts[vb + 2] - tile.data.verts[va + 2]);
							tmax = (neia[k * 2 + 1] - tile.data.verts[va + 2])
								/ (tile.data.verts[vb + 2] - tile.data.verts[va + 2]);
							if (tmin > tmax) {
								temp = tmin;
								tmin = tmax;
								tmax = temp;
							}
							link.bmin = Math.floor(DetourCommon.clamp(tmin, 0.0, 1.0) * 255.0);
							link.bmax = Math.floor(DetourCommon.clamp(tmax, 0.0, 1.0) * 255.0);
						} else if (dir == 2 || dir == 6) {
							tmin = (neia[k * 2 + 0] - tile.data.verts[va])
								/ (tile.data.verts[vb] - tile.data.verts[va]);
							tmax = (neia[k * 2 + 1] - tile.data.verts[va])
								/ (tile.data.verts[vb] - tile.data.verts[va]);
							if (tmin > tmax) {
								temp = tmin;
								tmin = tmax;
								tmax = temp;
							}
							link.bmin = Math.floor(DetourCommon.clamp(tmin, 0.0, 1.0) * 255.0);
							link.bmax = Math.floor(DetourCommon.clamp(tmax, 0.0, 1.0) * 255.0);
						}
					}
				}
			}
		}

		connectExtOffMeshLinks(tile, target, side) {
			if (tile == null)
				return;

			// Connect off-mesh links.
			// We are interested on links which land from target tile to this tile.
			let oppositeSide = (side == -1) ? 0xff : DetourCommon.oppositeTile(side);

			for (let i = 0; i < target.data.header.offMeshConCount; ++i) {
				let targetCon = target.data.offMeshCons[i];
				if (targetCon.side != oppositeSide)
					continue;

				let targetPoly = target.data.polys[targetCon.poly];
				// Skip off-mesh connections which start location could not be
				// connected at all.
				if (targetPoly.firstLink == NavMesh.DT_NULL_LINK)
					continue;

				let ext = [targetCon.rad, target.data.header.walkableClimb, targetCon.rad];

				// Find polygon to connect to.
				let p = new Array(3);
				p[0] = targetCon.pos[3];
				p[1] = targetCon.pos[4];
				p[2] = targetCon.pos[5];
				let nearest = this.findNearestPolyInTile(tile, p, ext);
				let ref = nearest.getNearestRef();
				if (ref == 0)
					continue;
				let nearestPt = nearest.getNearestPos();
				// findNearestPoly may return too optimistic results, further check
				// to make sure.

				if (DetourCommon.sqr(nearestPt[0] - p[0]) + DetourCommon.sqr(nearestPt[2] - p[2]) > DetourCommon.sqr(targetCon.rad))
					continue;
				// Make sure the location is on curren mesh.
				target.data.verts[targetPoly.verts[1] * 3] = nearestPt[0];
				target.data.verts[targetPoly.verts[1] * 3 + 1] = nearestPt[1];
				target.data.verts[targetPoly.verts[1] * 3 + 2] = nearestPt[2];

				// let off-mesh connection to target poly.
				let idx = this.allocLink(target);
				let link = target.links[idx];
				link.ref = ref;
				link.edge = 1;
				link.side = oppositeSide;
				link.bmin = link.bmax = 0;
				// Add to linked list.
				link.next = targetPoly.firstLink;
				targetPoly.firstLink = idx;

				// let target poly to off-mesh connection.
				if ((targetCon.flags & NavMesh.DT_OFFMESH_CON_BIDIR) != 0) {
					let tidx = this.allocLink(tile);
					let landPolyIdx = NavMesh.decodePolyIdPoly(ref);
					let landPoly = tile.data.polys[landPolyIdx];
					link = tile.links[tidx];
					link.ref = NavMesh.or(this.getPolyRefBase(target), (targetCon.poly));
					link.edge = 0xff;
					link.side = (side == -1 ? 0xff : side);
					link.bmin = link.bmax = 0;
					// Add to linked list.
					link.next = landPoly.firstLink;
					landPoly.firstLink = tidx;
				}
			}
		}

		findConnectingPolys(verts, va, vb, tile, side, maxcon) {
			if (tile == null)
				return [null, null, 0];
			let con = new Array(maxcon);
			let conarea = new Array(maxcon * 2);
			let amin = new Array(2);
			let amax = new Array(2);
			calcSlabEndPoints(verts, va, vb, amin, amax, side);
			apos = getSlabCoord(verts, va, side);

			// Remove links pointing to 'side' and compact the links array.
			let bmin = new Array(2);
			let bmax = new Array(2);
			let m = NavMesh.or(NavMesh.DT_EXT_LINK, side);
			let n = 0;
			let base = this.getPolyRefBase(tile);

			for (let i = 0; i < tile.data.header.polyCount; ++i) {
				let poly = tile.data.polys[i];
				let nv = poly.vertCount;
				for (let j = 0; j < nv; ++j) {
					// Skip edges which do not poPoly to the right side.
					if (poly.neis[j] != m)
						continue;
					let vc = poly.verts[j] * 3;
					let vd = poly.verts[(j + 1) % nv] * 3;
					bpos = getSlabCoord(tile.data.verts, vc, side);
					// Segments are not close enough.
					if (Math.abs(apos - bpos) > 0.01)
						continue;

					// Check if the segments touch.
					calcSlabEndPoints(tile.data.verts, vc, vd, bmin, bmax, side);

					if (!overlapSlabs(amin, amax, bmin, bmax, 0.01, tile.data.header.walkableClimb))
						continue;

					// Add return value.
					if (n < maxcon) {
						conarea[n * 2 + 0] = Math.max(amin[0], bmin[0]);
						conarea[n * 2 + 1] = Math.min(amax[0], bmax[0]);
						con[n] = NavMesh.or(base, i);
						n++;
					}
					break;
				}
			}
			return [con, conarea, n];
		}

		static getSlabCoord(verts, va, side) {
			if (side == 0 || side == 4)
				return verts[va];
			else if (side == 2 || side == 6)
				return verts[va + 2];
			return 0;
		}

		static calcSlabEndPoints(verts, va, vb, bmin, bmax, side) {
			if (side == 0 || side == 4) {
				if (verts[va + 2] < verts[vb + 2]) {
					bmin[0] = verts[va + 2];
					bmin[1] = verts[va + 1];
					bmax[0] = verts[vb + 2];
					bmax[1] = verts[vb + 1];
				} else {
					bmin[0] = verts[vb + 2];
					bmin[1] = verts[vb + 1];
					bmax[0] = verts[va + 2];
					bmax[1] = verts[va + 1];
				}
			} else if (side == 2 || side == 6) {
				if (verts[va + 0] < verts[vb + 0]) {
					bmin[0] = verts[va + 0];
					bmin[1] = verts[va + 1];
					bmax[0] = verts[vb + 0];
					bmax[1] = verts[vb + 1];
				} else {
					bmin[0] = verts[vb + 0];
					bmin[1] = verts[vb + 1];
					bmax[0] = verts[va + 0];
					bmax[1] = verts[va + 1];
				}
			}
		}

		overlapSlabs(amin, amax, bmin, bmax, px, py) {
			// Check for horizontal overlap.
			// The segment is shrunken a little so that slabs which touch
			// at end points are not connected.
			minx = Math.max(amin[0] + px, bmin[0] + px);
			maxx = Math.min(amax[0] - px, bmax[0] - px);
			if (minx > maxx)
				return false;

			// Check vertical overlap.
			ad = (amax[1] - amin[1]) / (amax[0] - amin[0]);
			ak = amin[1] - ad * amin[0];
			bd = (bmax[1] - bmin[1]) / (bmax[0] - bmin[0]);
			bk = bmin[1] - bd * bmin[0];
			aminy = ad * minx + ak;
			amaxy = ad * maxx + ak;
			bminy = bd * minx + bk;
			bmaxy = bd * maxx + bk;
			dmin = bminy - aminy;
			dmax = bmaxy - amaxy;

			// Crossing segments always overlap.
			if (dmin * dmax < 0)
				return true;

			// Check for overlap at endpoints.
			thr = (py * 2) * (py * 2);
			if (dmin * dmin <= thr || dmax * dmax <= thr)
				return true;

			return false;
		}

		/**
		 * Builds internal polygons links for a tile.
		 * 
		 * @param tile
		 */
		baseOffMeshLinks(tile) {
			if (tile == null)
				return;

			let base = this.getPolyRefBase(tile);

			// Base off-mesh connection start points.
			for (let i = 0; i < tile.data.header.offMeshConCount; ++i) {
				let con = tile.data.offMeshCons[i];
				let poly = tile.data.polys[con.poly];

				let ext = [con.rad, tile.data.header.walkableClimb, con.rad];

				// Find polygon to connect to.
				let nearestPoly = this.findNearestPolyInTile(tile, con.pos, ext);
				let ref = nearestPoly.getNearestRef();
				if (ref == 0)
					continue;
				let p = con.pos; // First vertex
				let nearestPt = nearestPoly.getNearestPos();
				// findNearestPoly may return too optimistic results, further check
				// to make sure.
				if (DetourCommon.sqr(nearestPt[0] - p[0]) + DetourCommon.sqr(nearestPt[2] - p[2]) > DetourCommon.sqr(con.rad))
					continue;
				// Make sure the location is on current mesh.
				tile.data.verts[poly.verts[0] * 3] = nearestPt[0];
				tile.data.verts[poly.verts[0] * 3 + 1] = nearestPt[1];
				tile.data.verts[poly.verts[0] * 3 + 2] = nearestPt[2];

				// let off-mesh connection to target poly.
				let idx = this.allocLink(tile);
				let link = tile.links[idx];
				link.ref = ref;
				link.edge = 0;
				link.side = 0xff;
				link.bmin = link.bmax = 0;
				// Add to linked list.
				link.next = poly.firstLink;
				poly.firstLink = idx;

				// Start end-poPoly is always connect back to off-mesh connection.
				let tidx = this.allocLink(tile);
				let landPolyIdx = NavMesh.decodePolyIdPoly(ref);
				let landPoly = tile.data.polys[landPolyIdx];
				link = tile.links[tidx];
				link.ref = NavMesh.or(base, (con.poly));
				link.edge = 0xff;
				link.side = 0xff;
				link.bmin = link.bmax = 0;
				// Add to linked list.
				link.next = landPoly.firstLink;
				landPoly.firstLink = tidx;
			}
		}

		/**
		 * Returns closest poPoly on polygon.
		 * 
		 * @param ref
		 * @param pos
		 * @return
		 */
		closestPointOnPoly(ref, pos) {
			let tileAndPoly = this.getTileAndPolyByRefUnsafe(ref);
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
			let verts = new Array(this.m_maxVertPerPoly * 3);
			let edged = new Array(this.m_maxVertPerPoly);
			let edget = new Array(this.m_maxVertPerPoly);
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

			// Find height at the location.
			let ip = poly.index;
			if (tile.data.detailMeshes != null && tile.data.detailMeshes.length > ip) {
				let pd = tile.data.detailMeshes[ip];
				for (let j = 0; j < pd.triCount; ++j) {
					let t = (pd.triBase + j) * 4;
					let v = [];//Was [3][]
					for (let k = 0; k < 3; ++k) {
						if (tile.data.detailTris[t + k] < poly.vertCount) {
							let index = poly.verts[tile.data.detailTris[t + k]] * 3;
							v[k] = [
								tile.data.verts[index], tile.data.verts[index + 1],
								tile.data.verts[index + 2]
							];
						} else {
							let index = (pd.vertBase + (tile.data.detailTris[t + k] - poly.vertCount)) * 3;
							v[k] = [
								tile.data.detailVerts[index], tile.data.detailVerts[index + 1],
								tile.data.detailVerts[index + 2]
							];
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

		findNearestPolyInTile(tile, center, extents) {
			let nearestPt = null;
			let bmin = DetourCommon.vSub(center, extents);
			let bmax = DetourCommon.vAdd(center, extents);

			// Get nearby polygons from proximity grid.
			let polys = this.queryPolygonsInTile(tile, bmin, bmax);

			// Find nearest polygon amongst the nearby polygons.
			let nearest = 0;
			let nearestDistanceSqr = Number.MAX_VALUE;
			for (let i = 0; i < polys.length; ++i) {
				let ref = polys[i];
				let d = 0;
				let cpp = this.closestPointOnPoly(ref, center);
				let posOverPoly = cpp.isPosOverPoly();
				let closestPtPoly = cpp.getClosest();

				// If a poPoly is directly over a polygon and closer than
				// climb height, favor that instead of straight line nearest point.
				let diff = DetourCommon.vSub(center, closestPtPoly);
				if (posOverPoly) {
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

		getTileAt(x, y, layer) {
			// Find tile based on hash.
			let h = NavMesh.computeTileHash(x, y, this.m_tileLutMask);
			let tile = this.m_posLookup[h];
			while (tile != null) {
				if (tile.data.header != null && tile.data.header.x == x && tile.data.header.y == y
					&& tile.data.header.layer == layer) {
					return tile;
				}
				tile = tile.next;
			}
			return null;
		}

		getNeighbourTilesAt(x, y, side) {
			let nx = x, ny = y;
			switch (side) {
				case 0:
					nx++;
					break;
				case 1:
					nx++;
					ny++;
					break;
				case 2:
					ny++;
					break;
				case 3:
					nx--;
					ny++;
					break;
				case 4:
					nx--;
					break;
				case 5:
					nx--;
					ny--;
					break;
				case 6:
					ny--;
					break;
				case 7:
					nx++;
					ny--;
					break;
			}
			return this.getTilesAt(nx, ny);
		}

		getTilesAt(x, y) {
			let tiles = [];
			// Find tile based on hash.
			let h = NavMesh.computeTileHash(x, y, this.m_tileLutMask);
			let tile = this.m_posLookup[h];
			while (tile != null) {
				if (tile.data.header != null && tile.data.header.x == x && tile.data.header.y == y) {
					tiles.push(tile);
				}
				tile = tile.next;
			}
			return tiles;
		}

		getTileRefAt(x, y, layer) {
			// Find tile based on hash.
			let h = NavMesh.computeTileHash(x, y, this.m_tileLutMask);
			let tile = this.m_posLookup[h];
			while (tile != null) {
				if (tile.data.header != null && tile.data.header.x == x && tile.data.header.y == y
					&& tile.data.header.layer == layer) {
					return this.getTileRef(tile);
				}
				tile = tile.next;
			}
			return 0;
		}

		getTileByRef(ref) {
			if (ref == 0)
				return null;
			let tileIndex = decodePolyIdTile(ref);
			let tileSalt = decodePolyIdSalt(ref);
			if (tileIndex >= this.m_maxTiles)
				return null;
			let tile = this.m_tiles[tileIndex];
			if (tile.salt != tileSalt)
				return null;
			return tile;
		}

		getTileRef(tile) {
			if (tile == null)
				return 0;
			let it = tile.index;
			return NavMesh.encodePolyId(tile.salt, it, 0);
		}

		static computeTileHash(x, y, mask) {
			let h1 = 0x8da6b343; // Large multiplicative constants;
			let h2 = 0xd8163841; // here arbitrarily chosen primes
			let n = h1 * x + h2 * y;
			return n & mask;
		}

		/// @par
		///
		/// Off-mesh connections are stored in the navigation mesh as special
		/// 2-vertex
		/// polygons with a single edge. At least one of the vertices is expected to
		/// be
		/// inside a normal polygon. So an off-mesh connection is "entered" from a
		/// normal polygon at one of its endpoints. This is the polygon identified
		/// by
		/// the prevRef parameter.
		getOffMeshConnectionPolyEndPoints(prevRef, polyRef) {
			if (polyRef == 0)
				throw new IllegalArgumentException("polyRef = 0");

			// Get current polygon
			let saltitip = NavMesh.decodePolyId(polyRef);
			let salt = saltitip[0];
			let it = saltitip[1];
			let ip = saltitip[2];
			if (it >= this.m_maxTiles) {
				throw new IllegalArgumentException("Invalid tile ID > max tiles");
			}
			if (this.m_tiles[it].salt != salt || this.m_tiles[it].data.header == null) {
				throw new IllegalArgumentException("Invalid salt or missing tile header");
			}
			let tile = this.m_tiles[it];
			if (ip >= tile.data.header.polyCount) {
				throw new IllegalArgumentException("Invalid poly ID > poly count");
			}
			let poly = tile.data.polys[ip];

			// Make sure that the current poly is indeed off-mesh link.
			if (poly.getType() != Poly.DT_POLYTYPE_OFFMESH_CONNECTION)
				throw new IllegalArgumentException("Invalid poly type");

			// Figure out which way to hand out the vertices.
			let idx0 = 0, idx1 = 1;

			// Find link that points to first vertex.
			for (let i = poly.firstLink; i != NavMesh.DT_NULL_LINK; i = tile.links[i].next) {
				if (tile.links[i].edge == 0) {
					if (tile.links[i].ref != prevRef) {
						idx0 = 1;
						idx1 = 0;
					}
					break;
				}
			}
			let startPos = new Array(3);
			let endPos = new Array(3);
			DetourCommon.vCopy(startPos, tile.data.verts, poly.verts[idx0] * 3);
			DetourCommon.vCopy(endPos, tile.data.verts, poly.verts[idx1] * 3);
			return [startPos, endPos];

		}

		getMaxVertsPerPoly() {
			return this.m_maxVertPerPoly;
		}

		getTileCount() {
			return this.m_tileCount;
		}
	}

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

	class RecastVectors$1 {

	static min( a, b,  i) {
			a[0] = Math.min(a[0], b[i + 0]);
			a[1] = Math.min(a[1], b[i + 1]);
			a[2] = Math.min(a[2], b[i + 2]);
		}

	static max( a, b,  i) {
			a[0] = Math.max(a[0], b[i + 0]);
			a[1] = Math.max(a[1], b[i + 1]);
			a[2] = Math.max(a[2], b[i + 2]);
		}

	static copy3( out, _in,  i) {
			RecastVectors$1.copy4(out, 0, _in, i);
		}

	static copy2( out, _in) {
			RecastVectors$1.copy4(out, 0, _in, 0);
		}

	static copy4( out,  n, _in,  m) {
			out[n] = _in[m];
			out[n + 1] = _in[m + 1];
			out[n + 2] = _in[m + 2];
		}

	static add( e0, a, verts,  i) {
			e0[0] = a[0] + verts[i];
			e0[1] = a[1] + verts[i + 1];
			e0[2] = a[2] + verts[i + 2];
		}

	static subA( e0, verts,  i,  j) {
			e0[0] = verts[i] - verts[j];
			e0[1] = verts[i + 1] - verts[j + 1];
			e0[2] = verts[i + 2] - verts[j + 2];
		}

	static subB( e0, i, verts,  j) {
			e0[0] = i[0] - verts[j];
			e0[1] = i[1] - verts[j + 1];
			e0[2] = i[2] - verts[j + 2];
		}

	static cross( dest, v1, v2) {
			dest[0] = v1[1] * v2[2] - v1[2] * v2[1];
			dest[1] = v1[2] * v2[0] - v1[0] * v2[2];
			dest[2] = v1[0] * v2[1] - v1[1] * v2[0];
		}

	static normalize( v) {
			let  d = (1.0 / Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
			v[0] *= d;
			v[1] *= d;
			v[2] *= d;
		}

	}

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



	class ChunkyTriMesh {

	    static BoundsItem = class BoundsItem {
	        bmin = new Array(2);
	        bmax = new Array(2);
	        i;
	    }

	    static ChunkyTriMeshNode = class ChunkyTriMeshNode {
	        bmin = new Array(2);
	        bmax = new Array(2);
	        i;
	        tris = [];
	    }

	    CompareItemX = function compare(a, b) {
	        if (a.bmin[0] < b.bmin[0]) {
	            return -1;
	        }
	        if (a.bmin[0] > b.bmin[0]) {
	            return 1;
	        }
	        return 0;
	    }


	    CompareItemY = function compare(a, b) {
	        if (a.bmin[1] < b.bmin[1]) {
	            return -1;
	        }
	        if (a.bmin[1] > b.bmin[1]) {
	            return 1;
	        }
	        return 0;
	    }




	    nodes = [];
	    ntris;
	    maxTrisPerChunk;

	    calcExtends(items, imin, imax, bmin, bmax) {
	        bmin[0] = items[imin].bmin[0];
	        bmin[1] = items[imin].bmin[1];

	        bmax[0] = items[imin].bmax[0];
	        bmax[1] = items[imin].bmax[1];

	        for (let i = imin + 1; i < imax; ++i) {
	            let it = items[i];
	            if (it.bmin[0] < bmin[0]) {
	                bmin[0] = it.bmin[0];
	            }
	            if (it.bmin[1] < bmin[1]) {
	                bmin[1] = it.bmin[1];
	            }

	            if (it.bmax[0] > bmax[0]) {
	                bmax[0] = it.bmax[0];
	            }
	            if (it.bmax[1] > bmax[1]) {
	                bmax[1] = it.bmax[1];
	            }
	        }
	    }

	    longestAxis(x, y) {
	        return y > x ? 1 : 0;
	    }


	    // https://stackoverflow.com/a/45245772/10047920
	    partialSort(arr, start, end, sortFx) {
	        let preSorted = arr.slice(0, start), postSorted = arr.slice(end);
	        let sorted = arr.slice(start, end).sort(sortFx);
	        arr.length = 0;
	        arr.push.apply(arr, preSorted.concat(sorted).concat(postSorted));
	        return arr;
	    }


	    subdivide(items, imin, imax, trisPerChunk, nodes, inTris) {
	        let inum = imax - imin;

	        let node = new ChunkyTriMesh.ChunkyTriMeshNode();
	        nodes.push(node);

	        if (inum <= trisPerChunk) {
	            // Leaf
	            this.calcExtends(items, imin, imax, node.bmin, node.bmax);

	            // Copy triangles.
	            node.i = nodes.length;
	            node.tris = new Array(inum * 3);

	            let dst = 0;
	            for (let i = imin; i < imax; ++i) {
	                let src = items[i].i * 3;
	                node.tris[dst++] = inTris[src];
	                node.tris[dst++] = inTris[src + 1];
	                node.tris[dst++] = inTris[src + 2];
	            }
	        } else {
	            // Split
	            this.calcExtends(items, imin, imax, node.bmin, node.bmax);

	            let axis = this.longestAxis(node.bmax[0] - node.bmin[0], node.bmax[1] - node.bmin[1]);

	            if (axis == 0) {
	                //Arrays.sort(items, imin, imax, new CompareItemX());
	                this.partialSort(items, imin, imax, (a,b)=>(this.CompareItemX(a,b)));
	                // let subarray = items.slice(imin, imax);
	                // subarray.sort((a, b) => (this.CompareItemX(a, b)))
	                // // let newitems = items.slice(0, imin);
	                // // newitems.push(...subarray)
	                // // newitems.push(...items.slice(imax, items.length));
	                // for (let i = imin; i < imax; i++) {
	                //     items[i].bmin[0] = subarray[i].bmin[0];
	                //     items[i].bmin[1] = subarray[i].bmin[1];
	                //     items[i].bmax[0] = subarray[i].bmax[0];
	                //     items[i].bmax[1] = subarray[i].bmax[1];
	                // }
	                //items = newitems;
	                // Sort aPoly x-axis
	            } else if (axis == 1) {
	                // Arrays.sort(items, imin, imax, new CompareItemY());
	                // Sort aPoly y-axis
	                this.partialSort(items, imin, imax, (a,b)=>(this.CompareItemY(a,b)));

	                // let subarray = items.slice(imin, imax);
	                // subarray.sort((a, b) => (this.CompareItemY(a, b)))
	                // // let newitems = items.slice(0, imin);
	                // // newitems.push(...subarray)
	                // // newitems.push(...items.slice(imax, items.length));
	                // for (let i = imin; i < imax; i++) {
	                //     items[i].bmin[0] = subarray[i].bmin[0];
	                //     items[i].bmin[1] = subarray[i].bmin[1];
	                //     items[i].bmax[0] = subarray[i].bmax[0];
	                //     items[i].bmax[1] = subarray[i].bmax[1];
	                // }
	                //items = newitems;
	            }

	            let isplit = Math.floor(imin + Math.floor(inum / 2));

	            let s = "";
	            for (let i = 0; i < items.length; i++) {
	                let item = items[i];
	                s += "" + i + " " + item.bmin[0] + " " + item.bmin[1] + "\n";
	                s += "" + i + " " + item.bmax[0] + " " + item.bmax[1] + "\n";
	            }
	            // fs.writeFileSync("items_" + imin + "_" + imax + ".txt", s);
	            // console.log("done")

	            // Left
	            // console.log("Before left " + imin + " " + isplit);
	            // console.log(items[279].bmin[1])
	            this.subdivide(items, imin, isplit, trisPerChunk, nodes, inTris);
	            // console.log("Done left " + imin + " " + isplit);
	            // console.log(items[279].bmin[1])
	            // Right
	            // console.log("Before right " + isplit + " " + imax);
	            // console.log(items[279].bmin[1])
	            this.subdivide(items, isplit, imax, trisPerChunk, nodes, inTris);
	            // console.log("Done right " + isplit + " " + imax);
	            // console.log(items[279].bmin[1])

	            // Negative index means escape.
	            node.i = -nodes.length;
	        }
	    }

	    constructor(verts, tris, ntris, trisPerChunk) {

	        this.nodes = [];
	        this.ntris = ntris;

	        // Build tree
	        let items = [];//new Array(ntris);

	        for (let i = 0; i < ntris; i++) {
	            let t = i * 3;
	            let it = items[i] = new ChunkyTriMesh.BoundsItem();
	            it.i = i;
	            // Calc triangle XZ bounds.
	            it.bmin[0] = it.bmax[0] = verts[tris[t] * 3 + 0];
	            it.bmin[1] = it.bmax[1] = verts[tris[t] * 3 + 2];
	            for (let j = 1; j < 3; ++j) {
	                let v = tris[t + j] * 3;
	                if (verts[v] < it.bmin[0]) {
	                    it.bmin[0] = verts[v];
	                }
	                if (verts[v + 2] < it.bmin[1]) {
	                    it.bmin[1] = verts[v + 2];
	                }

	                if (verts[v] > it.bmax[0]) {
	                    it.bmax[0] = verts[v];
	                }
	                if (verts[v + 2] > it.bmax[1]) {
	                    it.bmax[1] = verts[v + 2];
	                }
	            }
	        }

	        this.subdivide(items, 0, ntris, trisPerChunk, this.nodes, tris);

	        // Calc max tris per node.
	        this.maxTrisPerChunk = 0;
	        for (let node of this.nodes) {
	            let isLeaf = node.i >= 0;
	            if (!isLeaf) {
	                continue;
	            }
	            if (node.tris.length / 3 > this.maxTrisPerChunk) {
	               this. maxTrisPerChunk = node.tris.length / 3;
	            }
	        }

	    }

	    checkOverlapRect(amin, amax, bmin, bmax) {
	        let overlap = true;
	        overlap = (amin[0] > bmax[0] || amax[0] < bmin[0]) ? false : overlap;
	        overlap = (amin[1] > bmax[1] || amax[1] < bmin[1]) ? false : overlap;
	        return overlap;
	    }

	    getChunksOverlappingRect(bmin, bmax) {
	        // Traverse tree
	        ids = [];
	        let i = 0;
	        while (i < nodes.length) {
	            node = nodes[i];
	            let overlap = checkOverlapRect(bmin, bmax, node.bmin, node.bmax);
	            let isLeafNode = node.i >= 0;

	            if (isLeafNode && overlap) {
	                ids.push(node);
	            }

	            if (overlap || isLeafNode) {
	                i++;
	            } else {
	                i = -node.i;
	            }
	        }
	        return ids;
	    }

	}

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

	class TriMesh {

	vertices;
	 faces;
	 chunkyTriMesh;

	constructor( vertices, faces) {
	        this.vertices = vertices;
	        this.faces = faces;
	        this.chunkyTriMesh = new ChunkyTriMesh(vertices, faces, faces.length / 3, 256);
	    }

	getTris() {
	        return this.faces;
	    }

	 getVerts() {
	        return this.vertices;
	    }

	 getChunksOverlappingRect( bmin, bmax) {
	        return chunkyTriMesh.getChunksOverlappingRect(bmin, bmax);
	    }

	}

	class SimpleInputGeomProvider /*extends InputGeomProvider*/ {

	    vertices;
	    faces;
	    bmin;
	    bmax;
	    volumes = [];

	    static fromIndeces = function(vertexPositions, meshFaces) {
	       return new SimpleInputGeomProvider(this.mapVertices(vertexPositions), this.mapFaces(meshFaces));
	    }

	    static mapFaces(meshFaces) {
	        let faces = new Array(meshFaces.length);
	        for (let i = 0; i < faces.length; i++) {
	            faces[i] = meshFaces[i];
	        }
	        return faces;
	    }

	    static mapVertices(vertexPositions) {
	        let vertices = new Array(vertexPositions.length);
	        for (let i = 0; i < vertices.length; i++) {
	            vertices[i] = vertexPositions[i];
	        }
	        return vertices;
	    }

	    constructor(vertices, faces) {
	        this.vertices = vertices;
	        this.faces = faces;
	        this.bmin = new Array(3);
	        this.bmax = new Array(3);
	        RecastVectors$1.copy3(this.bmin, vertices, 0);
	        RecastVectors$1.copy3(this.bmax, vertices, 0);
	        for (let i = 1; i < vertices.length / 3; i++) {
	            RecastVectors$1.min(this.bmin, vertices, i * 3);
	            RecastVectors$1.max(this.bmax, vertices, i * 3);
	        }
	    }

	    getMeshBoundsMin() {
	        return this.bmin;
	    }

	    getMeshBoundsMax() {
	        return this.bmax;
	    }

	    getConvexVolumes() {
	        return [];
	    }

	    addConvexVolume(verts, minh, maxh, areaMod) {
	        let vol = new ConvexVolume();
	        vol.hmin = minh;
	        vol.hmax = maxh;
	        vol.verts = verts;
	        vol.areaMod = areaMod;
	        this.volumes.push(vol);
	    }

	    meshes() {
	        // return Collections.singPolyonList(new TriMesh(vertices, faces));
	       return [new TriMesh(this.vertices, this.faces)];
	    }
	}

	/*
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

	class ObjImporterContext {
	    vertexPositions = [];
	    meshFaces = [];
	}


	class ObjImporter {

	    // OBJImporterContext = 

	    load(is) {
	        let context = new ObjImporterContext();
	        try {
	            let slurp = is;
	            let lines = slurp.split(/\r?\n/);
	            for (let i = 0; i < lines.length; i++) {
	                let line = lines[i];
	                this.readLine(line, context);
	            }
	            // reader = new BufferedReader(new InputStreamReader(is));
	            // let line;
	            // while ((line = reader.readLine()) != null) {
	            //     line = line.trim();
	            //     readLine(line, context);
	            // }
	        } catch (e) {
	            throw e;
	        } finally {
	        }
	        return SimpleInputGeomProvider.fromIndeces(context.vertexPositions, context.meshFaces);

	    }

	    readLine(line, context) {
	        if (line.startsWith("v")) {
	            this.readVertex(line, context);
	        } else if (line.startsWith("f")) {
	            this.readFace(line, context);
	        }
	    }

	    readVertex(line, context) {
	        if (line.startsWith("v ")) {
	            let vert = this.readVector3f(line);
	            for (let vp of vert) {
	                context.vertexPositions.push(vp);
	            }
	        }
	    }

	    readVector3f(line) {
	        let v = line.split(/\s+/);
	        if (v.length < 4) {
	            throw new RuntimeException("Invalid vector, expected 3 coordinates, found " + (v.length - 1));
	        }
	        return [parseFloat(v[1]), parseFloat(v[2]), parseFloat(v[3])];
	    }

	    readFace(line, context) {
	        let v = line.split(/\s+/);
	        if (v.length < 4) {
	            throw new RuntimeException("Invalid number of face vertices: 3 coordinates expected, found "
	                + v.length);
	        }
	        for (let j = 0; j < v.length - 3; j++) {
	            context.meshFaces.push(this.readFaceVertex(v[1], context));
	            for (let i = 0; i < 2; i++) {
	                context.meshFaces.push(this.readFaceVertex(v[2 + j + i], context));
	            }
	        }
	    }

	    readFaceVertex(face, context) {
	        let v = face.split(/\//);
	        return this.getIndex(parseInt(v[0]), context.vertexPositions.length);
	    }

	    getIndex(posi, size) {
	        if (posi > 0) {
	            posi--;
	        } else if (posi < 0) {
	            posi = size + posi;
	        } else {
	            throw new RuntimeException("0 vertex index");
	        }
	        return posi;
	    }

	}

	/*
	Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
	Recast4J Copyright (c) 2015-2018 Piotr Piastucki piotr@jtilia.org

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

	class RecastConstants {

	    static RC_NULL_AREA = 0;
	    static RC_NOT_CONNECTED = 0x3f;
	    static RC_WALKABLE_AREA = 63;
	    /// Defines the number of bits allocated to rcSpan::smin and rcSpan::smax.
	    static RC_SPAN_HEIGHT_BITS = 13;
	    /// Defines the maximum value for rcSpan::smin and rcSpan::smax.
	    static RC_SPAN_MAX_HEIGHT = (1 << this.RC_SPAN_HEIGHT_BITS) - 1;
	    /// Heighfield border flag.
	    /// If a heightfield region ID has this bit set, then the region is a border
	    /// region and its spans are considered unwalkable.
	    /// (Used during the region and contour build process.)
	    /// @see rcCompactSpan::reg
	    static RC_BORDER_REG = 0x8000;
	    /// Polygon touches multiple regions.
	    /// If a polygon has this region ID it was merged with or created
	    /// from polygons of different regions during the polymesh
	    /// build step that removes redundant border vertices.
	    /// (Used during the polymesh and detail polymesh build processes)
	    /// @see rcPolyMesh::regs
	    static RC_MULTIPLE_REGS = 0;
	    // Border vertex flag.
	    /// If a region ID has this bit set, then the associated element lies on
	    /// a tile border. If a contour vertex's region ID has this bit set, the
	    /// vertex will later be removed in order to match the segments and vertices
	    /// at tile boundaries.
	    /// (Used during the build process.)
	    /// @see rcCompactSpan::reg, #rcContour::verts, #rcContour::rverts
	    static RC_BORDER_VERTEX = 0x10000;
	    /// Area border flag.
	    /// If a region ID has this bit set, then the associated element lies on
	    /// the border of an area.
	    /// (Used during the region and contour build process.)
	    /// @see rcCompactSpan::reg, #rcContour::verts, #rcContour::rverts
	    static RC_AREA_BORDER = 0x20000;
	    /// Applied to the region id field of contour vertices in order to extract the region id.
	    /// The region id field of a vertex may have several flags applied to it. So the
	    /// fields value can't be used directly.
	    /// @see rcContour::verts, rcContour::rverts
	    static RC_CONTOUR_REG_MASK = 0xffff;
	    /// A value which indicates an invalid index within a mesh.
	    /// @note This does not necessarily indicate an error.
	    /// @see rcPolyMesh::polys
	    static RC_MESH_NULL_IDX = 0xffff;

	    static RC_CONTOUR_TESS_WALL_EDGES = 0x01; /// < Tessellate solid (impassable) edges during contour simplification.
	    static RC_CONTOUR_TESS_AREA_EDGES = 0x02; /// < Tessellate edges between areas during contour simplification.


	    static WATERSHED = 10;
	    static MONOTONE = 20;
	    static LAYERS = 30;


	    static RC_LOG_WARNING = 1;

	}

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


	class RecastConfig {
		partitionType = null;

		/** The width/height size of tile's on the xz-plane. [Limit: >= 0] [Units: vx] **/
		tileSize;

		/** The xz-plane cell size to use for fields. [Limit: > 0] [Units: wu] **/
		cs;

		/** The y-axis cell size to use for fields. [Limit: > 0] [Units: wu] **/
		ch;

		/** The maximum slope that is considered walkable. [Limits: 0 <= value < 90] [Units: Degrees] **/
		walkableSlopeAngle;

		/**
		 * Minimum floor to 'ceiling' height that will still allow the floor area to be considered walkable. [Limit: >= 3]
		 * [Units: vx]
		 **/
		walkableHeight;

		/** Maximum ledge height that is considered to still be traversable. [Limit: >=0] [Units: vx] **/
		walkableClimb;

		/**
		 * The distance to erode/shrink the walkable area of the heightfield away from obstructions. [Limit: >=0] [Units:
		 * vx]
		 **/
		walkableRadius;

		/** The maximum allowed length for contour edges aPoly the border of the mesh. [Limit: >=0] [Units: vx] **/
		maxEdgeLen;

		/**
		 * The maximum distance a simplfied contour's border edges should deviate the original raw contour. [Limit: >=0]
		 * [Units: vx]
		 **/
		maxSimplificationError;

		/** The minimum number of cells allowed to form isolated island areas. [Limit: >=0] [Units: vx] **/
		minRegionArea;

		/**
		 * Any regions with a span count smaller than this value will, if possible, be merged with larger regions. [Limit:
		 * >=0] [Units: vx]
		 **/
		mergeRegionArea;

		/**
		 * The maximum number of vertices allowed for polygons generated during the contour to polygon conversion process.
		 * [Limit: >= 3]
		 **/
		maxVertsPerPoly = 0;

		/**
		 * Sets the sampling distance to use when generating the detail mesh. (For height detail only.) [Limits: 0 or >=
		 * 0.9] [Units: wu]
		 **/
		detailSampleDist;

		/**
		 * The maximum distance the detail mesh surface should deviate from heightfield data. (For height detail only.)
		 * [Limit: >=0] [Units: wu]
		 **/
		detailSampleMaxError;

		walkableAreaMod;

		constructor(partitionType, cellSize, cellHeight, agentHeight,
			agentRadius, agentMaxClimb, agentMaxSlope, regionMinSize, regionMergeSize,
			edgeMaxLen, edgeMaxError, vertsPerPoly, detailSampleDist, detailSampleMaxError,
			tileSize, walkableAreaMod) {
			this.partitionType = partitionType;
			this.cs = cellSize;
			this.ch = cellHeight;
			this.walkableSlopeAngle = agentMaxSlope;
			this.walkableHeight = Math.ceil(agentHeight / this.ch);
			this.walkableClimb = Math.floor(agentMaxClimb / this.ch);
			this.walkableRadius = Math.ceil(agentRadius / this.cs);
			this.maxEdgeLen = Math.floor(edgeMaxLen / cellSize);
			this.maxSimplificationError = edgeMaxError;
			this.minRegionArea = regionMinSize * regionMinSize; // Note: area = size*size
			this.mergeRegionArea = regionMergeSize * regionMergeSize; // Note: area = size*size
			this.maxVertsPerPoly = vertsPerPoly;
			this.detailSampleDist = detailSampleDist < 0.9 ? 0 : cellSize * detailSampleDist;
			this.detailSampleMaxError = cellHeight * detailSampleMaxError;
			this.tileSize = tileSize;
			this.walkableAreaMod = walkableAreaMod;
		}

	}

	class AreaModification {

	static  RC_AREA_FLAGS_MASK = 0x3F;
	value;
	mask;

		/**
		 * Mask is set to all available bits, which means value is fully applied
		 * 
		 * @param value
		 *            The area id to apply. [Limit: <= #RC_AREA_FLAGS_MASK]
		 */
	// constructor(value) {
	// 		this.value = value;
	// 		this.mask = AreaModification.RC_AREA_FLAGS_MASK;
	// 	}

		/**
		 * 
		 * @param value
		 *            The area id to apply. [Limit: <= #RC_AREA_FLAGS_MASK]
		 * @param mask
		 *            Bitwise mask used when applying value. [Limit: <= #RC_AREA_FLAGS_MASK]
		 */
	constructor(value,  mask = AreaModification.RC_AREA_FLAGS_MASK) {
			this.value = value;
			this.mask = mask;
		}

	static clone = function(other) {
		return new AreaModification(other.value, other.mask);
			// this.value = other.value;
			// this.mask = other.mask;
		}

	getMaskedValue() {
			return this.value & this.mask;
		}

	apply(area) {
			return ((this.value & this.mask) | (area & ~this.mask));
		}
	}

	class SampleAreaModifications {

	static SAMPLE_POLYAREA_TYPE_MASK = 0x07;
	static SAMPLE_POLYAREA_TYPE_GROUND = 0x1;
	static SAMPLE_POLYAREA_TYPE_WATER = 0x2;
	static SAMPLE_POLYAREA_TYPE_ROAD = 0x3;
	static SAMPLE_POLYAREA_TYPE_DOOR = 0x4;
	static SAMPLE_POLYAREA_TYPE_GRASS = 0x5;
	static SAMPLE_POLYAREA_TYPE_JUMP = 0x6;

	static SAMPLE_AREAMOD_GROUND = new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_GROUND,
	  SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK);
	static SAMPLE_AREAMOD_WATER = new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_WATER,
	  SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK);
	static SAMPLE_AREAMOD_ROAD = new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_ROAD,
	  SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK);
	static SAMPLE_AREAMOD_GRASS = new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_GRASS,
	  SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK);
	static SAMPLE_AREAMOD_DOOR = new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_DOOR,
	  SampleAreaModifications.SAMPLE_POLYAREA_TYPE_DOOR);
	static SAMPLE_AREAMOD_JUMP = new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_JUMP,
	  SampleAreaModifications.SAMPLE_POLYAREA_TYPE_JUMP);

	static  SAMPLE_POLYFLAGS_WALK = 0x01;	// Ability to walk (ground, grass, road)
	static  SAMPLE_POLYFLAGS_SWIM = 0x02;   // Ability to swim (water).
	static  SAMPLE_POLYFLAGS_DOOR = 0x04;   // Ability to move through doors.
	static  SAMPLE_POLYFLAGS_JUMP = 0x08;   // Ability to jump.
	static  SAMPLE_POLYFLAGS_DISABLED = 0x10; // Disabled polygon
	static  SAMPLE_POLYFLAGS_ALL = 0xfff; // All abilities.
	}

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

	/** A compact, static heightfield representing unobstructed space. */
	class CompactHeightfield {

		/** The width of the heightfield. (APoly the x-axis in cell units.) */
	 width = 0;
		/** The height of the heightfield. (APoly the z-axis in cell units.) */
	 height = 0;
		/** The number of spans in the heightfield. */
	 spanCount = 0;
		/** The walkable height used during the build of the field.  (See: RecastConfig::walkableHeight) */
	 walkableHeight = 0;
		/** The walkable climb used during the build of the field. (See: RecastConfig::walkableClimb) */
	 walkableClimb = 0;
		/** The AABB border size used during the build of the field. (See: RecastConfig::borderSize) */
	 borderSize = 0 ;
		/** The maximum distance value of any span within the field. */
	 maxDistance = 0;
		/** The maximum region id of any span within the field. */
	 maxRegions;
		/** The minimum bounds in world space. [(x, y, z)] */
	 bmin = new Array(3);
		/** The maximum bounds in world space. [(x, y, z)] */
	 bmax = new Array(3);
		/** The size of each cell. (On the xz-plane.) */
	 cs = 0;
		/** The height of each cell. (The minimum increment aPoly the y-axis.) */
	 ch = 0;
		/** Array of cells. [Size: #width*#height] */
	cells = [];
		/** Array of spans. [Size: #spanCount] */
	spans = [];
		/** Array containing border distance data. [Size: #spanCount] */
	 dist = [];
		/** Array containing area id data. [Size: #spanCount] */
	 areas = [];

	}

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

	/** Provides information on the content of a cell column in a compact heightfield. */
	class CompactCell {

		/** Index to the first span in the column. */
		index = 0;
		/** Number of spans in the column. */
		count = 0;

	}

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

	/**  Represents a span of unobstructed space within a compact heightfield. */
	class CompactSpan {

		/** The lower extent of the span. (Measured from the heightfield's base.) */
	y=0;
		/** The id of the region the span belongs to. (Or zero if not in a region.) */
	reg=0;
		/** Packed neighbor connection data. */
	con=0;
		/** The height of the span.  (Measured from #y.) */
	h=0;

	}

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

	class RecastCommon {

		/// Gets neighbor connection data for the specified direction.
		///  @param[in]		s		The span to check.
		///  @param[in]		dir		The direction to check. [Limits: 0 <= value < 4]
		///  @return The neighbor connection data for the specified direction,
		///  	or #RC_NOT_CONNECTED if there is no connection.	
		static GetCon = function (s, dir) {
			let shift = dir * 6;
			return (s.con >> shift) & 0x3f;
		}

		/// Gets the standard width (x-axis) offset for the specified direction.
		///  @param[in]		dir		The direction. [Limits: 0 <= value < 4]
		///  @return The width offset to apply to the current cell position to move
		///  	in the direction.
		static GetDirOffsetX = function (dir) {
			let offset = [- 1, 0, 1, 0,];
			return offset[dir & 0x03];
		}

		/// Gets the standard height (z-axis) offset for the specified direction.
		///  @param[in]		dir		The direction. [Limits: 0 <= value < 4]
		///  @return The height offset to apply to the current cell position to move
		///  	in the direction.
		static GetDirOffsetY = function (dir) {
			let offset = [0, 1, 0, - 1];

			return offset[dir & 0x03];
		}

		/// Gets the direction for the specified offset. One of x and y should be 0.
		/// @param[in] x The x offset. [Limits: -1 <= value <= 1]
		/// @param[in] y The y offset. [Limits: -1 <= value <= 1]
		/// @return The direction that represents the offset.
		static rcGetDirForOffset = function (x, y) {
			let dirs = [3, 0, - 1, 2, 1];

			return dirs[((y + 1) << 1) + x];
		}

		/// Sets the neighbor connection data for the specified direction.
		///  @param[in]		s		The span to update.
		///  @param[in]		dir		The direction to set. [Limits: 0 <= value < 4]
		///  @param[in]		i		The index of the neighbor span.
		static SetCon = function (s, dir, i) {
			let shift = dir * 6;
			let con = s.con;
			s.con = (con & ~(0x3f << shift)) | ((i & 0x3f) << shift);
		}

		static clamp = function (v, min, max) {
			return Math.max(Math.min(max, v), min);
		}

	}

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
			RecastVectors$1.subA(e0, verts, v1 * 3, v0 * 3);
			RecastVectors$1.subA(e1, verts, v2 * 3, v0 * 3);
			RecastVectors$1.cross(norm, e0, e1);
			RecastVectors$1.normalize(norm);
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
			RecastVectors$1.copy2(chf.bmin, hf.bmin);
			RecastVectors$1.copy2(chf.bmax, hf.bmax);
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

	class RecastBuilderConfig$1 {

	 cfg = null;

	    /** The width of the field aPoly the x-axis. [Limit: >= 0] [Units: vx] **/
	 width = 0;

	    /** The height of the field aPoly the z-axis. [Limit: >= 0] [Units: vx] **/
	 height = 0;

	    /** The minimum bounds of the field's AABB. [(x, y, z)] [Units: wu] **/
	 bmin = new Array(3);

	    /** The maximum bounds of the field's AABB. [(x, y, z)] [Units: wu] **/
	bmax = new Array(3);

	    /** The size of the non-navigable border around the heightfield. [Limit: >=0] [Units: vx] **/
	 borderSize = 0;

	    /** Set to true for tiled build **/
	 tiled = 0;

	    /** Set to false to disable building detailed mesh **/
	 buildMeshDetail = false;

	// RecastBuilderConfig(RecastConfig cfg, bmin, bmax) {
	//         this(cfg, bmin, bmax, 0, 0, false);
	//     }

	// RecastBuilderConfig(RecastConfig cfg, bmin, bmax,  tx,  ty,  tiled) {
	//         this(cfg, bmin, bmax, tx, ty, tiled, true);
	//     }

	constructor(cfg, bmin, bmax,  tx = 0,  ty = 0,  tiled = false,  buildMeshDetail = true) {
	        this.cfg = cfg;
	        this.tiled = tiled;
	        this.buildMeshDetail = buildMeshDetail;
	        RecastVectors$1.copy2(this.bmin, bmin);
	        RecastVectors$1.copy2(this.bmax, bmax);
	        if (tiled) {
	             ts = cfg.tileSize * cfg.cs;
	            this.bmin[0] += tx * ts;
	            this.bmin[2] += ty * ts;
	            this.bmax[0] = this.bmin[0] + ts;
	            this.bmax[2] = this.bmin[2] + ts;
	            // Expand the heighfield bounding box by border size to find the extents of geometry we need to build this
	            // tile.
	            //
	            // This is done in order to make sure that the navmesh tiles connect correctly at the borders,
	            // and the obstacles close to the border work correctly with the dilation process.
	            // No polygons (or contours) will be created on the border area.
	            //
	            // IMPORTANT!
	            //
	            // :''''''''':
	            // : +-----+ :
	            // : | | :
	            // : | |<--- tile to build
	            // : | | :
	            // : +-----+ :<-- geometry needed
	            // :.........:
	            //
	            // You should use this bounding box to query your input geometry.
	            //
	            // For example if you build a navmesh for terrain, and want the navmesh tiles to match the terrain tile size
	            // you will need to pass in data from neighbour terrain tiles too! In a simple case, just pass in all the 8
	            // neighbours,
	            // or use the bounding box below to only pass in a sliver of each of the 8 neighbours.
	            this.borderSize = cfg.walkableRadius + 3; // Reserve enough padding.
	            this.bmin[0] -= this.borderSize * cfg.cs;
	            this.bmin[2] -= this.borderSize * cfg.cs;
	            this.bmax[0] += this.borderSize * cfg.cs;
	            this.bmax[2] += this.borderSize * cfg.cs;
	            this.width = cfg.tileSize + this.borderSize * 2;
	            this.height = cfg.tileSize + this.borderSize * 2;
	            
	        } else {
	            let wh = Recast.calcGridSize(this.bmin, this.bmax, cfg.cs);
	            this.width = wh[0];
	            this.height = wh[1];
	            this.borderSize = 0;
	        }
	        //ADDED:
	        this.width = Math.floor(this.width);
	        this.height = Math.floor(this.height);
	    }

	}

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

	class Context {

		startTimer(string) {
			// TODO Auto-generated method stub
		}

		stopTimer(string) {
			// TODO Auto-generated method stub
		}

		warn(string) {
			console.log(string);
		}

	}

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

	/** Represents a heightfield layer within a layer set. */
	class Heightfield {

		/** The width of the heightfield. (APoly the x-axis in cell units.) */
	  width = 0;
		/** The height of the heightfield. (APoly the z-axis in cell units.) */
	  height = 0;
		/** The minimum bounds in world space. [(x, y, z)] */
	bmin;
		/** The maximum bounds in world space. [(x, y, z)] */
	 bmax;
		/** The size of each cell. (On the xz-plane.) */
	  cs;
		/** The height of each cell. (The minimum increment aPoly the y-axis.) */
	  ch;
		/** Heightfield of spans (width*height). */
	 spans = [];

	constructor(width,  height, bmin, bmax,  cs,  ch) {
			this.width = width;
			this.height = height;
			this.bmin = bmin;
			this.bmax = bmax;
			this.cs = cs;
			this.ch = ch;
			this.spans = new Array(width * height);

		}
	}

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

	/** Represents a span in a heightfield. */
	class Span {

		/** The lower limit of the span. [Limit: < #smax] */
	smin = 0;
		/** The upper limit of the span. [Limit: <= #RC_SPAN_MAX_HEIGHT] */
	smax = 0;
		/** The area id assigned to the span. */
	area = 0;
		/** The next span higher up in column. */
	next = null;

	}

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
					RecastVectors$1.copy4(buf, out2 + n * 3, buf, out1 + m * 3);
					m++;
					n++;
					// add the i'th poPoly to the right polygon. Do NOT add points that are on the dividing line
					// since these were already added above
					if (d[i] > 0) {
						RecastVectors$1.copy4(buf, out1 + m * 3, buf, _in + i * 3);
						m++;
					} else if (d[i] < 0) {
						RecastVectors$1.copy4(buf, out2 + n * 3, buf, _in + i * 3);
						n++;
					}
				} else // same side
				{
					// add the i'th poPoly to the right polygon. Addition is done even for points on the dividing line
					if (d[i] >= 0) {
						RecastVectors$1.copy4(buf, out1 + m * 3, buf, _in + i * 3);
						m++;
						if (d[i] != 0)
							continue;
					}
					RecastVectors$1.copy4(buf, out2 + n * 3, buf, _in + i * 3);
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
			RecastVectors$1.copy3(tmin, verts, v0 * 3);
			RecastVectors$1.copy3(tmax, verts, v0 * 3);
			RecastVectors$1.min(tmin, verts, v1 * 3);
			RecastVectors$1.min(tmin, verts, v2 * 3);
			RecastVectors$1.max(tmax, verts, v1 * 3);
			RecastVectors$1.max(tmax, verts, v2 * 3);

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

			RecastVectors$1.copy4(buf, 0, verts, v0 * 3);
			RecastVectors$1.copy4(buf, 3, verts, v1 * 3);
			RecastVectors$1.copy4(buf, 6, verts, v2 * 3);
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
								console.log("Here");
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

	/*
	Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
	Recast4J Copyright (c) 2015-2018 Piotr Piastucki piotr@jtilia.org

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

	class RecastRegion {

	    static RC_NULL_NEI = 0xfff;

	    static SweepSpan = class SweepSpan {
	        rid; // row id
	        id; // region id
	        ns; // number samples
	        nei; // neighbour id
	    }

	    static calculateDistanceField(chf, src) {
	        let maxDist;
	        let w = chf.width;
	        let h = chf.height;

	        // Init distance and points.
	        for (let i = 0; i < chf.spanCount; ++i) {
	            src[i] = 0xffff;
	        }

	        // Mark boundary cells.
	        for (let y = 0; y < h; ++y) {
	            for (let x = 0; x < w; ++x) {
	                let c = chf.cells[x + y * w];
	                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
	                    let s = chf.spans[i];
	                    let area = chf.areas[i];

	                    let nc = 0;
	                    for (let dir = 0; dir < 4; ++dir) {
	                        if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
	                            let ax = x + RecastCommon.GetDirOffsetX(dir);
	                            let ay = y + RecastCommon.GetDirOffsetY(dir);
	                            let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);
	                            if (area == chf.areas[ai]) {
	                                nc++;
	                            }
	                        }
	                    }
	                    if (nc != 4) {
	                        src[i] = 0;
	                    }
	                }
	            }
	        }

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
	                        if (src[ai] + 2 < src[i]) {
	                            src[i] = src[ai] + 2;
	                        }

	                        // (-1,-1)
	                        if (RecastCommon.GetCon(as, 3) != RecastConstants.RC_NOT_CONNECTED) {
	                            let aax = ax + RecastCommon.GetDirOffsetX(3);
	                            let aay = ay + RecastCommon.GetDirOffsetY(3);
	                            let aai = chf.cells[aax + aay * w].index + RecastCommon.GetCon(as, 3);
	                            if (src[aai] + 3 < src[i]) {
	                                src[i] = src[aai] + 3;
	                            }
	                        }
	                    }
	                    if (RecastCommon.GetCon(s, 3) != RecastConstants.RC_NOT_CONNECTED) {
	                        // (0,-1)
	                        let ax = x + RecastCommon.GetDirOffsetX(3);
	                        let ay = y + RecastCommon.GetDirOffsetY(3);
	                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 3);
	                        let as = chf.spans[ai];
	                        if (src[ai] + 2 < src[i]) {
	                            src[i] = src[ai] + 2;
	                        }

	                        // (1,-1)
	                        if (RecastCommon.GetCon(as, 2) != RecastConstants.RC_NOT_CONNECTED) {
	                            let aax = ax + RecastCommon.GetDirOffsetX(2);
	                            let aay = ay + RecastCommon.GetDirOffsetY(2);
	                            let aai = chf.cells[aax + aay * w].index + RecastCommon.GetCon(as, 2);
	                            if (src[aai] + 3 < src[i]) {
	                                src[i] = src[aai] + 3;
	                            }
	                        }
	                    }
	                }
	            }
	        }

	        // Pass 2
	        for (let y = h - 1; y >= 0; --y) {
	            for (let x = w - 1; x >= 0; --x) {
	                let c = chf.cells[x + y * w];
	                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
	                    let s = chf.spans[i];

	                    if (RecastCommon.GetCon(s, 2) != RecastConstants.RC_NOT_CONNECTED) {
	                        // (1,0)
	                        let ax = x + RecastCommon.GetDirOffsetX(2);
	                        let ay = y + RecastCommon.GetDirOffsetY(2);
	                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 2);
	                        let as = chf.spans[ai];
	                        if (src[ai] + 2 < src[i]) {
	                            src[i] = src[ai] + 2;
	                        }

	                        // (1,1)
	                        if (RecastCommon.GetCon(as, 1) != RecastConstants.RC_NOT_CONNECTED) {
	                            let aax = ax + RecastCommon.GetDirOffsetX(1);
	                            let aay = ay + RecastCommon.GetDirOffsetY(1);
	                            let aai = chf.cells[aax + aay * w].index + RecastCommon.GetCon(as, 1);
	                            if (src[aai] + 3 < src[i]) {
	                                src[i] = src[aai] + 3;
	                            }
	                        }
	                    }
	                    if (RecastCommon.GetCon(s, 1) != RecastConstants.RC_NOT_CONNECTED) {
	                        // (0,1)
	                        let ax = x + RecastCommon.GetDirOffsetX(1);
	                        let ay = y + RecastCommon.GetDirOffsetY(1);
	                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 1);
	                        let as = chf.spans[ai];
	                        if (src[ai] + 2 < src[i]) {
	                            src[i] = src[ai] + 2;
	                        }

	                        // (-1,1)
	                        if (RecastCommon.GetCon(as, 0) != RecastConstants.RC_NOT_CONNECTED) {
	                            let aax = ax + RecastCommon.GetDirOffsetX(0);
	                            let aay = ay + RecastCommon.GetDirOffsetY(0);
	                            let aai = chf.cells[aax + aay * w].index + RecastCommon.GetCon(as, 0);
	                            if (src[aai] + 3 < src[i]) {
	                                src[i] = src[aai] + 3;
	                            }
	                        }
	                    }
	                }
	            }
	        }

	        maxDist = 0;
	        for (let i = 0; i < chf.spanCount; ++i) {
	            maxDist = Math.max(src[i], maxDist);
	        }

	        return maxDist;
	    }

	    static boxBlur(chf, thr, src) {
	        let w = chf.width;
	        let h = chf.height;
	        let dst = new Array(chf.spanCount);

	        thr *= 2;

	        for (let y = 0; y < h; ++y) {
	            for (let x = 0; x < w; ++x) {
	                let c = chf.cells[x + y * w];
	                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
	                    let s = chf.spans[i];
	                    let cd = src[i];
	                    if (cd <= thr) {
	                        dst[i] = cd;
	                        continue;
	                    }

	                    let d = cd;
	                    for (let dir = 0; dir < 4; ++dir) {
	                        if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
	                            let ax = x + RecastCommon.GetDirOffsetX(dir);
	                            let ay = y + RecastCommon.GetDirOffsetY(dir);
	                            let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);
	                            d += src[ai];

	                            let as = chf.spans[ai];
	                            let dir2 = (dir + 1) & 0x3;
	                            if (RecastCommon.GetCon(as, dir2) != RecastConstants.RC_NOT_CONNECTED) {
	                                let ax2 = ax + RecastCommon.GetDirOffsetX(dir2);
	                                let ay2 = ay + RecastCommon.GetDirOffsetY(dir2);
	                                let ai2 = chf.cells[ax2 + ay2 * w].index + RecastCommon.GetCon(as, dir2);
	                                d += src[ai2];
	                            } else {
	                                d += cd;
	                            }
	                        } else {
	                            d += cd * 2;
	                        }
	                    }
	                    dst[i] = Math.floor((d + 5) / 9);
	                }
	            }
	        }
	        return dst;
	    }

	    static floodRegion(x, y, i, level, r, chf, srcReg,
	        srcDist, stack) {
	        let w = chf.width;

	        let area = chf.areas[i];

	        // Flood fill mark region.
	        stack = [];
	        stack.push(x);
	        stack.push(y);
	        stack.push(i);
	        srcReg[i] = r;
	        srcDist[i] = 0;

	        let lev = level >= 2 ? level - 2 : 0;
	        let count = 0;

	        while (stack.length > 0) {
	            let ci = stack.splice(stack.length - 1,1)[0];
	            let cy = stack.splice(stack.length - 1,1)[0];
	            let cx = stack.splice(stack.length - 1,1)[0];

	            let cs = chf.spans[ci];

	            // Check if any of the neighbours already have a valid region set.
	            let ar = 0;
	            for (let dir = 0; dir < 4; ++dir) {
	                // 8 connected
	                if (RecastCommon.GetCon(cs, dir) != RecastConstants.RC_NOT_CONNECTED) {
	                    let ax = cx + RecastCommon.GetDirOffsetX(dir);
	                    let ay = cy + RecastCommon.GetDirOffsetY(dir);
	                    let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(cs, dir);
	                    if (chf.areas[ai] != area) {
	                        continue;
	                    }
	                    let nr = srcReg[ai];
	                    if ((nr & RecastConstants.RC_BORDER_REG) != 0) {
	                        continue;
	                    }
	                    if (nr != 0 && nr != r) {
	                        ar = nr;
	                        break;
	                    }

	                    let as = chf.spans[ai];

	                    let dir2 = (dir + 1) & 0x3;
	                    if (RecastCommon.GetCon(as, dir2) != RecastConstants.RC_NOT_CONNECTED) {
	                        let ax2 = ax + RecastCommon.GetDirOffsetX(dir2);
	                        let ay2 = ay + RecastCommon.GetDirOffsetY(dir2);
	                        let ai2 = chf.cells[ax2 + ay2 * w].index + RecastCommon.GetCon(as, dir2);
	                        if (chf.areas[ai2] != area) {
	                            continue;
	                        }
	                        let nr2 = srcReg[ai2];
	                        if (nr2 != 0 && nr2 != r) {
	                            ar = nr2;
	                            break;
	                        }
	                    }
	                }
	            }
	            if (ar != 0) {
	                srcReg[ci] = 0;
	                continue;
	            }

	            count++;

	            // Expand neighbours.
	            for (let dir = 0; dir < 4; ++dir) {
	                if (RecastCommon.GetCon(cs, dir) != RecastConstants.RC_NOT_CONNECTED) {
	                    let ax = cx + RecastCommon.GetDirOffsetX(dir);
	                    let ay = cy + RecastCommon.GetDirOffsetY(dir);
	                    let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(cs, dir);
	                    if (chf.areas[ai] != area) {
	                        continue;
	                    }
	                    if (chf.dist[ai] >= lev && srcReg[ai] == 0) {
	                        srcReg[ai] = r;
	                        srcDist[ai] = 0;
	                        stack.push(ax);
	                        stack.push(ay);
	                        stack.push(ai);
	                    }
	                }
	            }
	        }

	        return count > 0;
	    }

	    static expandRegions(maxIter, level, chf, srcReg, srcDist,
	        stack, fillStack) {
	        let w = chf.width;
	        let h = chf.height;

	        if (fillStack) {
	            // Find cells revealed by the raised level.
	            //stack = [];
	            stack = [];
	            for (let y = 0; y < h; ++y) {
	                for (let x = 0; x < w; ++x) {
	                    let c = chf.cells[x + y * w];
	                    for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
	                        if (chf.dist[i] >= level && srcReg[i] == 0 && chf.areas[i] != RecastConstants.RC_NULL_AREA) {
	                            stack.push(x);
	                            stack.push(y);
	                            stack.push(i);
	                        }
	                    }
	                }
	            }
	        } else // use cells in the input stack
	        {
	            // mark all cells which already have a region
	            for (let j = 0; j < stack.length; j += 3) {
	                let i = stack[j + 2];
	                if (srcReg[i] != 0) {
	                    stack[j + 2] = -1;
	                }
	            }
	        }

	        let dirtyEntries = [];
	        let iter = 0;
	        while (stack.length > 0) {
	            let failed = 0;
	            // dirtyEntries = [];
	            dirtyEntries = [];

	            for (let j = 0; j < stack.length; j += 3) {
	                let x = stack[j + 0];
	                let y = stack[j + 1];
	                let i = stack[j + 2];
	                if (i < 0) {
	                    failed++;
	                    continue;
	                }

	                let r = srcReg[i];
	                let d2 = 0xfff;
	                let area = chf.areas[i];
	                let s = chf.spans[i];
	                for (let dir = 0; dir < 4; ++dir) {
	                    if (RecastCommon.GetCon(s, dir) == RecastConstants.RC_NOT_CONNECTED) {
	                        continue;
	                    }
	                    let ax = x + RecastCommon.GetDirOffsetX(dir);
	                    let ay = y + RecastCommon.GetDirOffsetY(dir);
	                    let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);
	                    if (chf.areas[ai] != area) {
	                        continue;
	                    }
	                    if (srcReg[ai] > 0 && (srcReg[ai] & RecastConstants.RC_BORDER_REG) == 0) {
	                        if (srcDist[ai] + 2 < d2) {
	                            r = srcReg[ai];
	                            d2 = srcDist[ai] + 2;
	                        }
	                    }
	                }
	                if (r != 0) {
	                    stack[j + 2]= -1; // mark as used
	                    dirtyEntries.push(i);
	                    dirtyEntries.push(r);
	                    dirtyEntries.push(d2);
	                } else {
	                    failed++;
	                }
	            }

	            // Copy entries that differ between src and dst to keep them in sync.
	            for (let i = 0; i < dirtyEntries.length; i += 3) {
	                let idx = dirtyEntries[i];
	                // if (idx == 1344)
	                //     console.log("dirty")
	                srcReg[idx] = dirtyEntries[i + 1];
	                srcDist[idx] = dirtyEntries[i + 2];
	            }

	            if (failed * 3 == stack.length) {
	                break;
	            }

	            if (level > 0) {
	                ++iter;
	                if (iter >= maxIter) {
	                    break;
	                }
	            }
	        }

	        return srcReg;
	    }

	    static sortCellsByLevel(startLevel, chf, srcReg, nbStacks,
	        stacks, loglevelsPerStack) // the levels per stack (2 in our case) as a bit shift
	    {
	        let w = chf.width;
	        let h = chf.height;
	        startLevel = startLevel >> loglevelsPerStack;

	        for (let j = 0; j < nbStacks; ++j) {
	            // stacks[j] = new Array(1024);
	            // stacks[j] = [];
	            stacks[j] = [];
	        }

	        // put all cells in the level range into the appropriate stacks
	        for (let y = 0; y < h; ++y) {
	            for (let x = 0; x < w; ++x) {
	                let c = chf.cells[x + y * w];
	                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
	                    if (chf.areas[i] == RecastConstants.RC_NULL_AREA || srcReg[i] != 0) {
	                        continue;
	                    }

	                    let level = chf.dist[i] >> loglevelsPerStack;
	                    let sId = startLevel - level;
	                    if (sId >= nbStacks) {
	                        continue;
	                    }
	                    if (sId < 0) {
	                        sId = 0;
	                    }

	                    stacks[sId].push(x);
	                    stacks[sId].push(y);
	                    stacks[sId].push(i);
	                }
	            }
	        }
	    }

	    static appendStacks(srcStack, dstStack, srcReg) {
	        for (let j = 0; j < srcStack.length; j += 3) {
	            let i = srcStack[j + 2];
	            if ((i < 0) || (srcReg[i] != 0)) {
	                continue;
	            }
	            dstStack.push(srcStack[j]);
	            dstStack.push(srcStack[j + 1]);
	            dstStack.push(srcStack[j + 2]);
	        }
	    }

	    static Region = class Region {
	        spanCount = 0; // Number of spans belonging to this region
	        id = 0; // ID of the region
	        areaType = 0; // Are type.
	        remap = false;
	        visited = false
	        overlap = false;
	        connectsToBorder = false;
	        ymin = 0;
	        ymax = 0;
	        connections = [];
	        floors = [];

	        constructor(i) {
	            this.id = i;
	            this.ymin = 0xFFFF;
	            this.connections = [];
	            this.floors = [];
	        }

	    };

	    static removeAdjacentNeighbours(reg) {
	        // Remove adjacent duplicates.
	        for (let i = 0; i < reg.connections.length && reg.connections.length > 1;) {
	            let ni = (i + 1) % reg.connections.length;
	            if (reg.connections[i] == reg.connections[ni] && reg.connections[i] >=-128 && reg.connections[i] <= 127) {
	                reg.connections.splice(i,1);
	            } else {
	                ++i;
	            }
	        }
	    }

	    static replaceNeighbour(reg, oldId, newId) {
	        let neiChanged = false;
	        for (let i = 0; i < reg.connections.length; ++i) {
	            if (reg.connections[i] == oldId) {
	                reg.connections[i] = newId;
	                neiChanged = true;
	            }
	        }
	        for (let i = 0; i < reg.floors.length; ++i) {
	            if (reg.floors[i] == oldId) {
	                reg.floors.set(i, newId);
	            }
	        }
	        if (neiChanged) {
	            RecastRegion.removeAdjacentNeighbours(reg);
	        }
	    }

	    static canMergeWithRegion(rega, regb) {
	        if (rega.areaType != regb.areaType) {
	            return false;
	        }
	        let n = 0;
	        for (let i = 0; i < rega.connections.length; ++i) {
	            if (rega.connections[i] == regb.id) {
	                n++;
	            }
	        }
	        if (n > 1) {
	            return false;
	        }
	        for (let i = 0; i < rega.floors.length; ++i) {
	            if (rega.floors[i] == regb.id) {
	                return false;
	            }
	        }
	        return true;
	    }

	    static addUniqueFloorRegion(reg, n) {
	        if (!reg.floors.includes(n)) {
	            reg.floors.push(n);
	        }
	    }

	    static mergeRegions(rega, regb) {
	        let aid = rega.id;
	        let bid = regb.id;

	        // Duplicate current neighbourhood.
	        let acon = rega.connections;
	        let bcon = regb.connections;

	        // Find insertion poPoly on A.
	        let insa = -1;
	        for (let i = 0; i < acon.length; ++i) {
	            if (acon[i] == bid) {
	                insa = i;
	                break;
	            }
	        }
	        if (insa == -1) {
	            return false;
	        }

	        // Find insertion poPoly on B.
	        let insb = -1;
	        for (let i = 0; i < bcon.length; ++i) {
	            if (bcon[i] == aid) {
	                insb = i;
	                break;
	            }
	        }
	        if (insb == -1) {
	            return false;
	        }

	        // Merge neighbours.
	        rega.connections = [];
	        for (let i = 0, ni = acon.length; i < ni - 1; ++i) {
	            rega.connections.push(acon[(insa + 1 + i) % ni]);
	        }

	        for (let i = 0, ni = bcon.length; i < ni - 1; ++i) {
	            rega.connections.push(bcon[(insb + 1 + i) % ni]);
	        }

	        RecastRegion.removeAdjacentNeighbours(rega);

	        for (let j = 0; j < regb.floors.length; ++j) {
	            RecastRegion.addUniqueFloorRegion(rega, regb.floors[j]);
	        }
	        rega.spanCount += regb.spanCount;
	        regb.spanCount = 0;
	        regb.connections = [];

	        return true;
	    }

	    static isRegionConnectedToBorder(reg) {
	        // Region is connected to border if
	        // one of the neighbours is null id.
	        return reg.connections.includes(0);
	    }

	    static isSolidEdge(chf, srcReg, x, y, i, dir) {
	        let s = chf.spans[i];
	        let r = 0;
	        if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
	            let ax = x + RecastCommon.GetDirOffsetX(dir);
	            let ay = y + RecastCommon.GetDirOffsetY(dir);
	            let ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(s, dir);
	            r = srcReg[ai];
	        }
	        if (r == srcReg[i]) {
	            return false;
	        }
	        return true;
	    }

	    static walkContour(x, y, i, dir, chf, srcReg,
	        cont) {
	        let startDir = dir;
	        let starti = i;

	        let ss = chf.spans[i];
	        let curReg = 0;
	        if (RecastCommon.GetCon(ss, dir) != RecastConstants.RC_NOT_CONNECTED) {
	            let ax = x + RecastCommon.GetDirOffsetX(dir);
	            let ay = y + RecastCommon.GetDirOffsetY(dir);
	            let ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(ss, dir);
	            curReg = srcReg[ai];
	        }
	        cont.push(curReg);

	        let iter = 0;
	        while (++iter < 40000) {
	            let s = chf.spans[i];

	            if (RecastRegion.isSolidEdge(chf, srcReg, x, y, i, dir)) {
	                // Choose the edge corner
	                let r = 0;
	                if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
	                    let ax = x + RecastCommon.GetDirOffsetX(dir);
	                    let ay = y + RecastCommon.GetDirOffsetY(dir);
	                    let ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(s, dir);
	                    r = srcReg[ai];
	                }
	                if (r != curReg) {
	                    curReg = r;
	                    cont.push(curReg);
	                }

	                dir = (dir + 1) & 0x3; // Rotate CW
	            } else {
	                let ni = -1;
	                let nx = x + RecastCommon.GetDirOffsetX(dir);
	                let ny = y + RecastCommon.GetDirOffsetY(dir);
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

	        // Remove adjacent duplicates.
	        if (cont.length > 1) {
	            for (let j = 0; j < cont.length;) {
	                let nj = (j + 1) % cont.length;
	                if (cont[j] == cont[nj] && cont[j]>=-128 && cont[j] <=127) {
	                    cont.splice(j,1);
	                } else {
	                    ++j;
	                }
	            }
	        }
	    }

	    static mergeAndFilterRegions(ctx, minRegionArea, mergeRegionSize, maxRegionId,
	        chf, srcReg, overlaps) {
	        let w = chf.width;
	        let h = chf.height;

	        let nreg = maxRegionId + 1;
	        let regions = new Array(nreg);

	        // Construct regions
	        for (let i = 0; i < nreg; ++i) {
	            regions[i] = new RecastRegion.Region(i);
	        }

	        // Find edge of a region and find connections around the contour.
	        for (let y = 0; y < h; ++y) {
	            for (let x = 0; x < w; ++x) {
	                let c = chf.cells[x + y * w];
	                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
	                    let r = srcReg[i];
	                    if (r == 0 || r >= nreg) {
	                        continue;
	                    }

	                    let reg = regions[r];
	                    reg.spanCount++;

	                    // Update floors.
	                    for (let j = c.index; j < ni; ++j) {
	                        if (i == j) {
	                            continue;
	                        }
	                        let floorId = srcReg[j];
	                        if (floorId == 0 || floorId >= nreg) {
	                            continue;
	                        }
	                        if (floorId == r) {
	                            reg.overlap = true;
	                        }
	                        RecastRegion.addUniqueFloorRegion(reg, floorId);
	                    }

	                    // Have found contour
	                    if (reg.connections.length > 0) {
	                        continue;
	                    }

	                    reg.areaType = chf.areas[i];

	                    // Check if this cell is next to a border.
	                    let ndir = -1;
	                    for (let dir = 0; dir < 4; ++dir) {
	                        if (RecastRegion.isSolidEdge(chf, srcReg, x, y, i, dir)) {
	                            ndir = dir;
	                            break;
	                        }
	                    }

	                    if (ndir != -1) {
	                        // The cell is at border.
	                        // Walk around the contour to find all the neighbours.
	                        RecastRegion.walkContour(x, y, i, ndir, chf, srcReg, reg.connections);
	                    }
	                }
	            }
	        }

	        // Remove too small regions.
	        let stack = new Array(32);
	        let trace = new Array(32);
	        for (let i = 0; i < nreg; ++i) {
	            let reg = regions[i];
	            if (reg.id == 0 || (reg.id & RecastConstants.RC_BORDER_REG) != 0) {
	                continue;
	            }
	            if (reg.spanCount == 0) {
	                continue;
	            }
	            if (reg.visited) {
	                continue;
	            }

	            // Count the total size of all the connected regions.
	            // Also keep track of the regions connects to a tile border.
	            let connectsToBorder = false;
	            let spanCount = 0;
	            //stack = [];
	            //trace = [];
	            stack = [];
	            trace = [];

	            reg.visited = true;
	            stack.push(i);

	            while (stack.length > 0) {
	                // Pop
	                let ri = stack.splice(stack.length - 1, 1);

	                let creg = regions[ri];

	                spanCount += creg.spanCount;
	                trace.push(ri);

	                for (let j = 0; j < creg.connections.length; ++j) {
	                    if ((creg.connections[j] & RecastConstants.RC_BORDER_REG) != 0) {
	                        connectsToBorder = true;
	                        continue;
	                    }
	                    let neireg = regions[creg.connections[j]];
	                    if (neireg.visited) {
	                        continue;
	                    }
	                    if (neireg.id == 0 || (neireg.id & RecastConstants.RC_BORDER_REG) != 0) {
	                        continue;
	                    }
	                    // Visit
	                    stack.push(neireg.id);
	                    neireg.visited = true;
	                }
	            }

	            // If the accumulated regions size is too small, remove it.
	            // Do not remove areas which connect to tile borders
	            // as their size cannot be estimated correctly and removing them
	            // can potentially remove necessary areas.
	            if (spanCount < minRegionArea && !connectsToBorder) {
	                // Kill all visited regions.
	                for (let j = 0; j < trace.length; ++j) {
	                    regions[trace[j]].spanCount = 0;
	                    regions[trace[j]].id = 0;
	                }
	            }
	        }

	        // Merge too small regions to neighbour regions.
	        let mergeCount = 0;
	        do {
	            mergeCount = 0;
	            for (let i = 0; i < nreg; ++i) {
	                let reg = regions[i];
	                if (reg.id == 0 || (reg.id & RecastConstants.RC_BORDER_REG) != 0) {
	                    continue;
	                }
	                if (reg.overlap) {
	                    continue;
	                }
	                if (reg.spanCount == 0) {
	                    continue;
	                }

	                // Check to see if the region should be merged.
	                if (reg.spanCount > mergeRegionSize && RecastRegion.isRegionConnectedToBorder(reg)) {
	                    continue;
	                }

	                // Small region with more than 1 connection.
	                // Or region which is not connected to a border at all.
	                // Find smallest neighbour region that connects to this one.
	                let smallest = 0xfffffff;
	                let mergeId = reg.id;
	                for (let j = 0; j < reg.connections.length; ++j) {
	                    if ((reg.connections[j] & RecastConstants.RC_BORDER_REG) != 0) {
	                        continue;
	                    }
	                    let mreg = regions[reg.connections[j]];
	                    if (mreg.id == 0 || (mreg.id & RecastConstants.RC_BORDER_REG) != 0 || mreg.overlap) {
	                        continue;
	                    }
	                    if (mreg.spanCount < smallest && RecastRegion.canMergeWithRegion(reg, mreg) && RecastRegion.canMergeWithRegion(mreg, reg)) {
	                        smallest = mreg.spanCount;
	                        mergeId = mreg.id;
	                    }
	                }
	                // Found new id.
	                if (mergeId != reg.id) {
	                    let oldId = reg.id;
	                    let target = regions[mergeId];

	                    // Merge neighbours.
	                    if (RecastRegion.mergeRegions(target, reg)) {
	                        // Fixup regions pointing to current region.
	                        for (let j = 0; j < nreg; ++j) {
	                            if (regions[j].id == 0 || (regions[j].id & RecastConstants.RC_BORDER_REG) != 0) {
	                                continue;
	                            }
	                            // If another region was already merged into current region
	                            // change the nid of the previous region too.
	                            if (regions[j].id == oldId) {
	                                regions[j].id = mergeId;
	                            }
	                            // Replace the current region with the new one if the
	                            // current regions is neighbour.
	                            RecastRegion.replaceNeighbour(regions[j], oldId, mergeId);
	                        }
	                        mergeCount++;
	                    }
	                }
	            }
	        } while (mergeCount > 0);

	        // Compress region Ids.
	        for (let i = 0; i < nreg; ++i) {
	            regions[i].remap = false;
	            if (regions[i].id == 0) {
	                continue; // Skip nil regions.
	            }
	            if ((regions[i].id & RecastConstants.RC_BORDER_REG) != 0) {
	                continue; // Skip external regions.
	            }
	            regions[i].remap = true;
	        }

	        let regIdGen = 0;
	        for (let i = 0; i < nreg; ++i) {
	            if (!regions[i].remap) {
	                continue;
	            }
	            let oldId = regions[i].id;
	            let newId = ++regIdGen;
	            for (let j = i; j < nreg; ++j) {
	                if (regions[j].id == oldId) {
	                    regions[j].id = newId;
	                    regions[j].remap = false;
	                }
	            }
	        }
	        maxRegionId = regIdGen;

	        // Remap regions.
	        for (let i = 0; i < chf.spanCount; ++i) {
	            if ((srcReg[i] & RecastConstants.RC_BORDER_REG) == 0) {
	                srcReg[i] = regions[srcReg[i]].id;
	            }
	        }

	        // Return regions that we found to be overlapping.
	        for (let i = 0; i < nreg; ++i) {
	            if (regions[i].overlap) {
	                overlaps.push(regions[i].id);
	            }
	        }

	        return maxRegionId;
	    }

	    static addUniqueConnection(reg, n) {
	        if (!reg.connections.contains(n)) {
	            reg.connections.push(n);
	        }
	    }

	    static mergeAndFilterLayerRegions(ctx, minRegionArea, maxRegionId,
	        chf, srcReg, overlaps) {
	        let w = chf.width;
	        let h = chf.height;

	        let nreg = maxRegionId + 1;
	        let regions = new Array(nreg);

	        // Construct regions
	        for (let i = 0; i < nreg; ++i) {
	            regions[i] = new Region(i);
	        }

	        // Find region neighbours and overlapping regions.
	        let lregs = new Array(32);
	        for (let y = 0; y < h; ++y) {
	            for (let x = 0; x < w; ++x) {
	                let c = chf.cells[x + y * w];

	                lregs = [];

	                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
	                    let s = chf.spans[i];
	                    let ri = srcReg[i];
	                    if (ri == 0 || ri >= nreg) {
	                        continue;
	                    }
	                    let reg = regions[ri];

	                    reg.spanCount++;

	                    reg.ymin = Math.min(reg.ymin, s.y);
	                    reg.ymax = Math.max(reg.ymax, s.y);
	                    // Collect all region layers.
	                    lregs.push(ri);

	                    // Update neighbours
	                    for (let dir = 0; dir < 4; ++dir) {
	                        if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
	                            let ax = x + RecastCommon.GetDirOffsetX(dir);
	                            let ay = y + RecastCommon.GetDirOffsetY(dir);
	                            let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);
	                            let rai = srcReg[ai];
	                            if (rai > 0 && rai < nreg && rai != ri) {
	                                addUniqueConnection(reg, rai);
	                            }
	                            if ((rai & RecastConstants.RC_BORDER_REG) != 0) {
	                                reg.connectsToBorder = true;
	                            }
	                        }
	                    }

	                }

	                // Update overlapping regions.
	                for (let i = 0; i < lregs.length - 1; ++i) {
	                    for (let j = i + 1; j < lregs.length; ++j) {
	                        if (lregs[i] != lregs[j]) {
	                            let ri = regions[lregs[i]];
	                            let rj = regions[lregs[j]];
	                            RecastRegion.addUniqueFloorRegion(ri, lregs[j]);
	                            RecastRegion.addUniqueFloorRegion(rj, lregs[i]);
	                        }
	                    }
	                }

	            }
	        }

	        // Create 2D layers from regions.
	        let layerId = 1;

	        for (let i = 0; i < nreg; ++i) {
	            regions[i].id = 0;
	        }

	        // Merge montone regions to create non-overlapping areas.
	        let stack = new Array(32);
	        for (let i = 1; i < nreg; ++i) {
	            let root = regions[i];
	            // Skip already visited.
	            if (root.id != 0) {
	                continue;
	            }

	            // Start search.
	            root.id = layerId;

	            stack = [];
	            stack.push(i);

	            while (stack.length > 0) {
	                // Pop front
	                let reg = regions[stack.remove(0)];

	                let ncons = reg.connections.length;
	                for (let j = 0; j < ncons; ++j) {
	                    let nei = reg.connections[j];
	                    let regn = regions[nei];
	                    // Skip already visited.
	                    if (regn.id != 0) {
	                        continue;
	                    }
	                    // Skip if the neighbour is overlapping root region.
	                    let overlap = false;
	                    for(let k = 0; k < root.floors.length; k++) {
	                        if (root.floors[k] == nei) {
	                            overlap = true;
	                            break;
	                        }
	                    }
	                    if (overlap) {
	                        continue;
	                    }

	                    // Deepen
	                    stack.push(nei);

	                    // Mark layer id
	                    regn.id = layerId;
	                    // Merge current layers to root.
	                    for(let k = 0; k < regn.floors.length; ++k) {
	                        RecastRegion.addUniqueFloorRegion(root, regn.floors[k]);
	                    }
	                    root.ymin = Math.min(root.ymin, regn.ymin);
	                    root.ymax = Math.max(root.ymax, regn.ymax);
	                    root.spanCount += regn.spanCount;
	                    regn.spanCount = 0;
	                    root.connectsToBorder = root.connectsToBorder || regn.connectsToBorder;
	                }
	            }

	            layerId++;
	        }

	        // Remove small regions
	        for (let i = 0; i < nreg; ++i) {
	            if (regions[i].spanCount > 0 && regions[i].spanCount < minRegionArea && !regions[i].connectsToBorder) {
	                let reg = regions[i].id;
	                for (let j = 0; j < nreg; ++j) {
	                    if (regions[j].id == reg) {
	                        regions[j].id = 0;
	                    }
	                }
	            }
	        }

	        // Compress region Ids.
	        for (let i = 0; i < nreg; ++i) {
	            regions[i].remap = false;
	            if (regions[i].id == 0) {
	                continue; // Skip nil regions.
	            }
	            if ((regions[i].id & RecastConstants.RC_BORDER_REG) != 0) {
	                continue; // Skip external regions.
	            }
	            regions[i].remap = true;
	        }

	        let regIdGen = 0;
	        for (let i = 0; i < nreg; ++i) {
	            if (!regions[i].remap) {
	                continue;
	            }
	            let oldId = regions[i].id;
	            let newId = ++regIdGen;
	            for (let j = i; j < nreg; ++j) {
	                if (regions[j].id == oldId) {
	                    regions[j].id = newId;
	                    regions[j].remap = false;
	                }
	            }
	        }
	        maxRegionId = regIdGen;

	        // Remap regions.
	        for (let i = 0; i < chf.spanCount; ++i) {
	            if ((srcReg[i] & RecastConstants.RC_BORDER_REG) == 0) {
	                srcReg[i] = regions[srcReg[i]].id;
	            }
	        }

	        return maxRegionId;
	    }

	    /// @par
	    ///
	    /// This is usually the second to the last step in creating a fully built
	    /// compact heightfield. This step is required before regions are built
	    /// using #rcBuildRegions or #rcBuildRegionsMonotone.
	    ///
	    /// After this step, the distance data is available via the rcCompactHeightfield::maxDistance
	    /// and rcCompactHeightfield::dist fields.
	    ///
	    /// @see rcCompactHeightfield, rcBuildRegions, rcBuildRegionsMonotone
	    static buildDistanceField(ctx, chf) {

	        ctx.startTimer("BUILD_DISTANCEFIELD");
	        let src = new Array(chf.spanCount);
	        ctx.startTimer("DISTANCEFIELD_DIST");

	        let maxDist = this.calculateDistanceField(chf, src);
	        chf.maxDistance = maxDist;

	        ctx.stopTimer("DISTANCEFIELD_DIST");

	        ctx.startTimer("DISTANCEFIELD_BLUR");

	        // Blur
	        src = this.boxBlur(chf, 1, src);

	        // Store distance.
	        chf.dist = src;

	        ctx.stopTimer("DISTANCEFIELD_BLUR");

	        ctx.stopTimer("BUILD_DISTANCEFIELD");

	    }

	    static paintRectRegion(minx, maxx, miny, maxy, regId, chf,
	        srcReg) {
	        let w = chf.width;
	        for (let y = miny; y < maxy; ++y) {
	            for (let y = minx; x < maxx; ++x) {
	                let c = chf.cells[x + y * w];
	                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
	                    if (chf.areas[i] != RecastConstants.RC_NULL_AREA) {
	                        srcReg[i] = regId;
	                    }
	                }
	            }
	        }
	    }

	    /// @par
	    ///
	    /// Non-null regions will consist of connected, non-overlapping walkable spans that form a single contour.
	    /// Contours will form simple polygons.
	    ///
	    /// If multiple regions form an area that is smaller than @p minRegionArea, then all spans will be
	    /// re-assigned to the zero (null) region.
	    ///
	    /// Partitioning can result in smaller than necessary regions. @p mergeRegionArea helps
	    /// reduce unecessarily small regions.
	    ///
	    /// See the #rcConfig documentation for more information on the configuration parameters.
	    ///
	    /// The region data will be available via the rcCompactHeightfield::maxRegions
	    /// and rcCompactSpan::reg fields.
	    ///
	    /// @warning The distance field must be created using #rcBuildDistanceField before attempting to build regions.
	    ///
	    /// @see rcCompactHeightfield, rcCompactSpan, rcBuildDistanceField, rcBuildRegionsMonotone, rcConfig
	    static buildRegionsMonotone(ctx, chf, borderSize, minRegionArea,
	        mergeRegionArea) {
	        ctx.startTimer("BUILD_REGIONS");

	        let w = chf.width;
	        let h = chf.height;
	        let id = 1;

	        let srcReg = new Array(chf.spanCount);

	        let nsweeps = Math.max(chf.width, chf.height);
	        sweeps = new Array(nsweeps);
	        for (let i = 0; i < sweeps.length; i++) {
	            sweeps[i] = new SweepSpan();
	        }

	        // Mark border regions.
	        if (borderSize > 0) {
	            // Make sure border will not overflow.
	            let bw = Math.min(w, borderSize);
	            let bh = Math.min(h, borderSize);
	            // PaPoly regions
	            paintRectRegion(0, bw, 0, h, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
	            id++;
	            paintRectRegion(w - bw, w, 0, h, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
	            id++;
	            paintRectRegion(0, w, 0, bh, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
	            id++;
	            paintRectRegion(0, w, h - bh, h, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
	            id++;

	        }

	        chf.borderSize = borderSize;

	        let prev = new Array(256);

	        // Sweep one line at a time.
	        for (let y = borderSize; y < h - borderSize; ++y) {
	            // Collect spans from this row.
	            prev.fill(0, 0, id);
	            let rid = 1;

	            for (let y = borderSize; x < w - borderSize; ++x) {
	                let c = chf.cells[x + y * w];

	                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
	                    let s = chf.spans[i];
	                    if (chf.areas[i] == RecastConstants.RC_NULL_AREA) {
	                        continue;
	                    }

	                    // -x
	                    let previd = 0;
	                    if (RecastCommon.GetCon(s, 0) != RecastConstants.RC_NOT_CONNECTED) {
	                        let ax = x + RecastCommon.GetDirOffsetX(0);
	                        let ay = y + RecastCommon.GetDirOffsetY(0);
	                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 0);
	                        if ((srcReg[ai] & RecastConstants.RC_BORDER_REG) == 0 && chf.areas[i] == chf.areas[ai]) {
	                            previd = srcReg[ai];
	                        }
	                    }

	                    if (previd == 0) {
	                        previd = rid++;
	                        sweeps[previd].rid = previd;
	                        sweeps[previd].ns = 0;
	                        sweeps[previd].nei = 0;
	                    }

	                    // -y
	                    if (RecastCommon.GetCon(s, 3) != RecastConstants.RC_NOT_CONNECTED) {
	                        let ax = x + RecastCommon.GetDirOffsetX(3);
	                        let ay = y + RecastCommon.GetDirOffsetY(3);
	                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 3);
	                        if (srcReg[ai] != 0 && (srcReg[ai] & RecastConstants.RC_BORDER_REG) == 0 && chf.areas[i] == chf.areas[ai]) {
	                            let nr = srcReg[ai];
	                            if (sweeps[previd].nei == 0 || sweeps[previd].nei == nr) {
	                                sweeps[previd].nei = nr;
	                                sweeps[previd].ns++;
	                                prev[nr]++;
	                            } else {
	                                sweeps[previd].nei = RC_NULL_NEI;
	                            }
	                        }
	                    }

	                    srcReg[i] = previd;
	                }
	            }

	            // Create unique ID.
	            for (let i = 1; i < rid; ++i) {
	                if (sweeps[i].nei != RC_NULL_NEI && sweeps[i].nei != 0 && prev[sweeps[i].nei] == sweeps[i].ns) {
	                    sweeps[i].id = sweeps[i].nei;
	                } else {
	                    sweeps[i].id = id++;
	                }
	            }

	            // Remap IDs
	            for (let y = borderSize; x < w - borderSize; ++x) {
	                let c = chf.cells[x + y * w];

	                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
	                    if (srcReg[i] > 0 && srcReg[i] < rid) {
	                        srcReg[i] = sweeps[srcReg[i]].id;
	                    }
	                }
	            }
	        }

	        ctx.startTimer("BUILD_REGIONS_FILTER");

	        // Merge regions and filter out small regions.
	        let overlaps = [];
	        chf.maxRegions = mergeAndFilterRegions(ctx, minRegionArea, mergeRegionArea, id, chf, srcReg, overlaps);

	        // Monotone partitioning does not generate overlapping regions.

	        ctx.stopTimer("BUILD_REGIONS_FILTER");

	        // Store the result out.
	        for (let i = 0; i < chf.spanCount; ++i) {
	            // if (i == 1344)
	            //     console.log("4431")
	            chf.spans[i].reg = srcReg[i];
	        }

	        ctx.stopTimer("BUILD_REGIONS");

	    }

	    /// @par
	    ///
	    /// Non-null regions will consist of connected, non-overlapping walkable spans that form a single contour.
	    /// Contours will form simple polygons.
	    ///
	    /// If multiple regions form an area that is smaller than @p minRegionArea, then all spans will be
	    /// re-assigned to the zero (null) region.
	    ///
	    /// Watershed partitioning can result in smaller than necessary regions, especially in diagonal corridors.
	    /// @p mergeRegionArea helps reduce unecessarily small regions.
	    ///
	    /// See the #rcConfig documentation for more information on the configuration parameters.
	    ///
	    /// The region data will be available via the rcCompactHeightfield::maxRegions
	    /// and rcCompactSpan::reg fields.
	    ///
	    /// @warning The distance field must be created using #rcBuildDistanceField before attempting to build regions.
	    ///
	    /// @see rcCompactHeightfield, rcCompactSpan, rcBuildDistanceField, rcBuildRegionsMonotone, rcConfig
	    static buildRegions(ctx, chf, borderSize, minRegionArea,
	        mergeRegionArea) {
	        ctx.startTimer("BUILD_REGIONS");

	        let w = chf.width;
	        let h = chf.height;

	        ctx.startTimer("REGIONS_WATERSHED");

	        let LOG_NB_STACKS = 3;
	        let NB_STACKS = 1 << LOG_NB_STACKS;
	        let lvlStacks = [];
	        for (let i = 0; i < NB_STACKS; ++i) {
	            lvlStacks.push([]);
	        }

	        let stack = new Array(1024);

	        let srcReg = new Array(chf.spanCount).fill(0);
	        let srcDist = new Array(chf.spanCount).fill(0);

	        let regionId = 1;
	        let level = (chf.maxDistance + 1) & ~1;

	        // TODO: Figure better formula, expandIters defines how much the
	        // watershed "overflows" and simplifies the regions. Tying it to
	        // agent radius was usually good indication how greedy it could be.
	        // const let expandIters = 4 + walkableRadius * 2;
	        let expandIters = 8;

	        if (borderSize > 0) {
	            // Make sure border will not overflow.
	            let bw = Math.min(w, borderSize);
	            let bh = Math.min(h, borderSize);
	            // PaPoly regions
	            paintRectRegion(0, bw, 0, h, regionId | RecastConstants.RC_BORDER_REG, chf, srcReg);
	            regionId++;
	            paintRectRegion(w - bw, w, 0, h, regionId | RecastConstants.RC_BORDER_REG, chf, srcReg);
	            regionId++;
	            paintRectRegion(0, w, 0, bh, regionId | RecastConstants.RC_BORDER_REG, chf, srcReg);
	            regionId++;
	            paintRectRegion(0, w, h - bh, h, regionId | RecastConstants.RC_BORDER_REG, chf, srcReg);
	            regionId++;

	        }

	        chf.borderSize = borderSize;

	        let sId = -1;
	        while (level > 0) {
	            level = level >= 2 ? level - 2 : 0;
	            sId = (sId + 1) & (NB_STACKS - 1);

	            // ctx=>startTimer(RC_TIMER_DIVIDE_TO_LEVELS);

	            if (sId == 0) {
	                RecastRegion.sortCellsByLevel(level, chf, srcReg, NB_STACKS, lvlStacks, 1);
	            } else {
	                RecastRegion.appendStacks(lvlStacks[sId - 1], lvlStacks[sId], srcReg); // copy left overs from last level
	            }

	            // ctx=>stopTimer(RC_TIMER_DIVIDE_TO_LEVELS);

	            ctx.startTimer("BUILD_REGIONS_EXPAND");

	            // Expand current regions until no empty connected cells found.
	            RecastRegion.expandRegions(expandIters, level, chf, srcReg, srcDist, lvlStacks[sId], false);

	            ctx.stopTimer("BUILD_REGIONS_EXPAND");

	            ctx.startTimer("BUILD_REGIONS_FLOOD");

	            // Mark new regions with IDs.
	            for (let j = 0; j < lvlStacks[sId].length; j += 3) {
	                let x = lvlStacks[sId][j];
	                let y = lvlStacks[sId][j + 1];
	                let i = lvlStacks[sId][j + 2];
	                if (i >= 0 && srcReg[i] == 0) {
	                    if (RecastRegion.floodRegion(x, y, i, level, regionId, chf, srcReg, srcDist, stack)) {
	                        regionId++;
	                    }
	                }
	            }

	            ctx.stopTimer("BUILD_REGIONS_FLOOD");
	        }

	        // Expand current regions until no empty connected cells found.
	        RecastRegion.expandRegions(expandIters * 8, 0, chf, srcReg, srcDist, stack, true);

	        ctx.stopTimer("BUILD_REGIONS_WATERSHED");

	        ctx.startTimer("BUILD_REGIONS_FILTER");

	        // Merge regions and filter out smalle regions.
	        let overlaps = [];
	        chf.maxRegions = this.mergeAndFilterRegions(ctx, minRegionArea, mergeRegionArea, regionId, chf, srcReg, overlaps);

	        // If overlapping regions were found during merging, split those regions.
	        if (overlaps.length > 0) {
	            throw new RuntimeException("rcBuildRegions: " + overlaps.length + " overlapping regions.");
	        }

	        ctx.stopTimer("BUILD_REGIONS_FILTER");

	        // Write the result out.
	        for (let i = 0; i < chf.spanCount; ++i) {
	            // if (i == 1344)
	            //     console.log("4431")
	            chf.spans[i].reg = srcReg[i];
	        }

	        ctx.stopTimer("BUILD_REGIONS");

	    }

	    static buildLayerRegions(ctx, chf, borderSize, minRegionArea) {

	        ctx.startTimer("BUILD_REGIONS");

	        let w = chf.width;
	        let h = chf.height;
	        let id = 1;

	        let srcReg = new Array(chf.spanCount);
	        let nsweeps = Math.max(chf.width, chf.height);
	        let sweeps = new ArrayList(nsweeps);
	        for (let i = 0; i < sweeps.length; i++) {
	            sweeps[i] = new SweepSpan();
	        }

	        // Mark border regions.
	        if (borderSize > 0) {
	            // Make sure border will not overflow.
	            let bw = Math.min(w, borderSize);
	            let bh = Math.min(h, borderSize);
	            // PaPoly regions
	            paintRectRegion(0, bw, 0, h, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
	            id++;
	            paintRectRegion(w - bw, w, 0, h, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
	            id++;
	            paintRectRegion(0, w, 0, bh, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
	            id++;
	            paintRectRegion(0, w, h - bh, h, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
	            id++;

	        }

	        chf.borderSize = borderSize;

	        let prev = new Array(256);

	        // Sweep one line at a time.
	        for (let y = borderSize; y < h - borderSize; ++y) {
	            // Collect spans from this row.
	            prev.fill(0, 0, id);
	            let rid = 1;

	            for (let y = borderSize; x < w - borderSize; ++x) {
	                let c = chf.cells[x + y * w];

	                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
	                    let s = chf.spans[i];
	                    if (chf.areas[i] == RecastConstants.RC_NULL_AREA) {
	                        continue;
	                    }

	                    // -x
	                    let previd = 0;
	                    if (RecastCommon.GetCon(s, 0) != RecastConstants.RC_NOT_CONNECTED) {
	                        let ax = x + RecastCommon.GetDirOffsetX(0);
	                        let ay = y + RecastCommon.GetDirOffsetY(0);
	                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 0);
	                        if ((srcReg[ai] & RecastConstants.RC_BORDER_REG) == 0 && chf.areas[i] == chf.areas[ai]) {
	                            previd = srcReg[ai];
	                        }
	                    }

	                    if (previd == 0) {
	                        previd = rid++;
	                        sweeps[previd].rid = previd;
	                        sweeps[previd].ns = 0;
	                        sweeps[previd].nei = 0;
	                    }

	                    // -y
	                    if (RecastCommon.GetCon(s, 3) != RecastConstants.RC_NOT_CONNECTED) {
	                        let ax = x + RecastCommon.GetDirOffsetX(3);
	                        let ay = y + RecastCommon.GetDirOffsetY(3);
	                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 3);
	                        if (srcReg[ai] != 0 && (srcReg[ai] & RecastConstants.RC_BORDER_REG) == 0 && chf.areas[i] == chf.areas[ai]) {
	                            let nr = srcReg[ai];
	                            if (sweeps[previd].nei == 0 || sweeps[previd].nei == nr) {
	                                sweeps[previd].nei = nr;
	                                sweeps[previd].ns++;
	                                prev[nr]++;
	                            } else {
	                                sweeps[previd].nei = RC_NULL_NEI;
	                            }
	                        }
	                    }

	                    srcReg[i] = previd;
	                }
	            }

	            // Create unique ID.
	            for (let i = 1; i < rid; ++i) {
	                if (sweeps[i].nei != RC_NULL_NEI && sweeps[i].nei != 0 && prev[sweeps[i].nei] == sweeps[i].ns) {
	                    sweeps[i].id = sweeps[i].nei;
	                } else {
	                    sweeps[i].id = id++;
	                }
	            }

	            // Remap IDs
	            for (let y = borderSize; x < w - borderSize; ++x) {
	                let c = chf.cells[x + y * w];

	                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
	                    if (srcReg[i] > 0 && srcReg[i] < rid) {
	                        srcReg[i] = sweeps[srcReg[i]].id;
	                    }
	                }
	            }
	        }

	        ctx.startTimer("BUILD_REGIONS_FILTER");

	        // Merge monotone regions to layers and remove small regions.
	        let overlaps = [];
	        chf.maxRegions = mergeAndFilterLayerRegions(ctx, minRegionArea, id, chf, srcReg, overlaps);

	        ctx.stopTimer("BUILD_REGIONS_FILTER");

	        // Store the result out.
	        for (let i = 0; i < chf.spanCount; ++i) {
	            chf.spans[i].reg = srcReg[i];
	        }

	        ctx.stopTimer("BUILD_REGIONS");

	    }
	}

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


	/** Represents a group of related contours. */
	class ContourSet {

		/** A list of the contours in the set. */
	  conts = [];
		/** The minimum bounds in world space. [(x, y, z)] */
	 bmin = new Array(3);
		/** The maximum bounds in world space. [(x, y, z)] */
	 bmax = new Array(3);
		/** The size of each cell. (On the xz-plane.) */
	 cs = 0;
		/** The height of each cell. (The minimum increment aPoly the y-axis.) */
	 ch = 0;
		/** The width of the set. (APoly the x-axis in cell units.) */
	 width = 0;
		/** The height of the set. (APoly the z-axis in cell units.) */
	 height = 0;
		/** The AABB border size used to generate the source data from which the contours were derived. */
	 borderSize = 0;
		/** The max edge error that this contour set was simplified with. */
	 maxError = 0;
	}

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

	/** Represents a polygon mesh suitable for use in building a navigation mesh. */
	class PolyMesh {

		/** The mesh vertices. [Form: (x, y, z) coordinates * #nverts] */
		verts = [];
		/** Polygon and neighbor data. [Length: #maxpolys * 2 * #nvp] */
		polys = [];
		/** The region id assigned to each polygon. [Length: #maxpolys] */
		regs = [];
		/** The area id assigned to each polygon. [Length: #maxpolys] */
		areas = [];
		/** The number of vertices. */
		nverts = 0;
		/** The number of polygons. */
		npolys = 0;
		/** The maximum number of vertices per polygon. */
		nvp = 0;
		/** The number of allocated polygons. */
		maxpolys = 0;
		/** The user defined flags for each polygon. [Length: #maxpolys] */
		flags = [];
		/** The minimum bounds in world space. [(x, y, z)] */
		bmin = new Array(3);
		/** The maximum bounds in world space. [(x, y, z)] */
		bmax = new Array(3);
		/** The size of each cell. (On the xz-plane.) */
		cs = 0;
		/** The height of each cell. (The minimum increment aPoly the y-axis.) */
		ch = 0;
		/** The AABB border size used to generate the source data from which the mesh was derived. */
		borderSize = 0;
		/** The max error of the polygon edges in the mesh. */
		maxEdgeError = 0;
	}

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


	function arraycopy$1(one, oneStart, two, twoStart, len) {
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
			arraycopy$1(polys, tmp, polys, pa, nvp);
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
						arraycopy$1(mesh.polys, p2, mesh.polys, p, nvp);
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
							arraycopy$1(polys, last, polys, pb, nvp);
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
			RecastVectors$1.copy3(mesh.bmin, cset.bmin, 0);
			RecastVectors$1.copy3(mesh.bmax, cset.bmax, 0);
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
										console.log("break");
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
			RecastVectors$1.copy(mesh.bmin, meshes[0].bmin, 0);
			RecastVectors$1.copy(mesh.bmax, meshes[0].bmax, 0);

			let maxVerts = 0;
			let maxPolys = 0;
			let maxVertsPerMesh = 0;
			for (let i = 0; i < nmeshes; ++i) {
				RecastVectors$1.min(mesh.bmin, meshes[i].bmin, 0);
				RecastVectors$1.max(mesh.bmax, meshes[i].bmax, 0);
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
			RecastVectors$1.copy(dst.bmin, src.bmin, 0);
			RecastVectors$1.copy(dst.bmax, src.bmax, 0);
			dst.cs = src.cs;
			dst.ch = src.ch;
			dst.borderSize = src.borderSize;
			dst.maxEdgeError = src.maxEdgeError;

			dst.verts = new Array(src.nverts * 3);
			arraycopy$1(src.verts, 0, dst.verts, 0, dst.verts.length);
			dst.polys = new Array(src.npolys * 2 * src.nvp);
			arraycopy$1(src.polys, 0, dst.polys, 0, dst.polys.length);
			dst.regs = new Array(src.npolys);
			arraycopy$1(src.regs, 0, dst.regs, 0, dst.regs.length);
			dst.areas = new Array(src.npolys);
			arraycopy$1(src.areas, 0, dst.areas, 0, dst.areas.length);
			dst.flags = new Array(src.npolys);
			arraycopy$1(src.flags, 0, dst.flags, 0, dst.flags.length);
			return dst;
		}
	}

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

	/** Represents a simple, non-overlapping contour in field space. */
	class Contour {

		/** Simplified contour vertex and connection data. [Size: 4 * #nverts] */
	verts = [];
		/** The number of vertices in the simplified contour. */
	nverts = 0;
		/** Raw contour vertex and connection data. [Size: 4 * #nrverts] */
	rverts = [];
		/** The number of vertices in the raw contour.  */
	nrverts = 0;
		/** The region id of the contour. */
	area = 0;
		/** The area id of the contour. */
	reg = 0;

	}

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
			RecastVectors$1.copy3(cset.bmin, chf.bmin, 0);
			RecastVectors$1.copy3(cset.bmax, chf.bmax, 0);
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
			let simplified = new Array(64);
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

	/** Contains triangle meshes that represent detailed height data associated with the polygons in its associated polygon mesh object. */
	class PolyMeshDetail {

		/** The sub-mesh data. [Size: 4*#nmeshes] */
	 meshes = [];
		/** The mesh vertices. [Size: 3*#nverts] */
	 verts = [];
		/** The mesh triangles. [Size: 4*#ntris] */
	 tris = [];
		/** The number of sub-meshes defined by #meshes. */
	 nmeshes =0;
		/** The number of vertices in #verts. */
	 nverts =0;
		/** The number of triangles in #tris. */
	 ntris=0;

	}

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

	function arraycopy$2(one, oneStart, two, twoStart, len) {
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
			RecastVectors$1.sub(v2, verts, p2, p1);
			RecastVectors$1.sub(v3, verts, p3, p1);

			cp = vcross2(v1, v2, v3);
			if (Math.abs(cp) > EPS) {
				v1Sq = RecastMeshDetail.vdot2(v1, v1);
				v2Sq = RecastMeshDetail.vdot2(v2, v2);
				v3Sq = RecastMeshDetail.vdot2(v3, v3);
				c[0] = (v1Sq * (v2[2] - v3[2]) + v2Sq * (v3[2] - v1[2]) + v3Sq * (v1[2] - v2[2])) / (2 * cp);
				c[1] = 0;
				c[2] = (v1Sq * (v3[0] - v2[0]) + v2Sq * (v1[0] - v3[0]) + v3Sq * (v2[0] - v1[0])) / (2 * cp);
				r.set(RecastMeshDetail.vdist2(c, v1));
				RecastVectors$1.add(c, c, verts, p1);
				return true;
			}
			RecastVectors$1.copy(c, verts, p1);
			r.set(0);
			return false;
		}

		static distPtTri(p, verts, a, b, c) {
			let v0 = new Array(3);
			let v1 = new Array(3);
			let v2 = new Array(3);
			RecastVectors$1.subA(v0, verts, c, a);
			RecastVectors$1.subA(v1, verts, b, a);
			RecastVectors$1.subB(v2, p, verts, a);

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
			let edges = new Array(64);		for (let i = 0, j = nhull - 1; i < nhull; j = i++)
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
				RecastVectors$1.copy4(verts, i * 3, _in, i * 3);
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
							RecastVectors$1.copy(verts, nverts * 3, edge, idx[k] * 3);
							hull[nhull++] = nverts;
							nverts++;
						}
					} else {
						for (let k = 1; k < nidx - 1; ++k) {
							RecastVectors$1.copy(verts, nverts * 3, edge, idx[k] * 3);
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
				RecastVectors$1.copy3(bmin, _in, 0);
				RecastVectors$1.copy3(bmax, _in, 0);
				for (let i = 1; i < nin; ++i) {
					RecastVectors$1.min(bmin, _in, i * 3);
					RecastVectors$1.max(bmax, _in, i * 3);
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
					RecastVectors$1.copy(verts, nverts * 3, bestpt, 0);
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
			dmesh.meshes.fill(0);

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
						arraycopy$2(dmesh.verts, 0, newv, 0, 3 * dmesh.nverts);
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
						arraycopy$2(dmesh.tris, 0, newt, 0, 4 * dmesh.ntris);
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
					RecastVectors$1.copy(mesh.verts, mesh.nverts * 3, dm.verts, k * 3);
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

	class RecastBuilderResult {
	    pmesh;
	    dmesh;

	    constructor(pmesh, dmesh) {
	        this.pmesh = pmesh;
	        this.dmesh = dmesh;
	    }

	    getMesh() {
	        return this.pmesh;
	    }

	    getMeshDetail() {
	        return this.dmesh;
	    }
	}


	class RecastBuilder {

	    // class RecastBuilderProgressListener {
	    //     onProgress(completed, total);
	    // }

	    progressListener;

	    // RecastBuilder() {
	    //     progressListener = null;
	    // }

	    constructor(progressListener = null) {
	        this.progressListener = progressListener;
	    }



	    buildTiles(geom, cfg, threads) {
	        let bmin = geom.getMeshBoundsMin();
	        let bmax = geom.getMeshBoundsMax();
	        let twh = Recast.calcTileCount(bmin, bmax, cfg.cs, cfg.tileSize);
	        let tw = twh[0];
	        let th = twh[1];
	        result = null;
	        if (threads == 1) {
	            result = buildSingleThread(geom, cfg, bmin, bmax, tw, th);
	        } else {
	            result = buildMultiThread(geom, cfg, bmin, bmax, tw, th, threads);
	        }
	        return result;
	    }

	    buildSingleThread(geom, cfg, bmin, bmax, tw,
	        th) {
	        result = new RecastBuilderResult[tw][th];
	        counter = new AtomicInteger();
	        for(let x = 0; x < tw; ++x) {
	            for(let y = 0; y < th; ++y) {
	                result[x][y] = buildTile(geom, cfg, bmin, bmax, x, y, counter, tw * th);
	            }
	        }
	        return result;
	    }

	    buildMultiThread(geom, cfg, bmin, bmax, tw, th,
	        threads) {
	        ec = Executors.newFixedThreadPool(threads);
	        result = new RecastBuilderResult[tw][th];
	        counter = new AtomicInteger();
	        for(let x = 0; x < tw; ++x) {
	            for(let y = 0; y < th; ++y) {
	                let tx = x;
	                let ty = y;
	                ec.submit(() => {
	                    result[tx][ty] = buildTile(geom, cfg, bmin, bmax, tx, ty, counter, tw * th);
	                });
	            }
	        }
	        ec.shutdown();
	        try {
	            ec.awaitTermination(1000, TimeUnit.HOURS);
	        } catch (e) {
	        }
	        return result;
	    }

	    buildTile(geom, cfg, bmin, bmax, tx, ty,
	        counter, total) {
	        result = build(geom, new RecastBuilderConfig(cfg, bmin, bmax, tx, ty, true));
	        if (progressListener != null) {
	            progressListener.onProgress(counter.incrementAndGet(), total);
	        }
	        return result;
	    }

	    build(geom, builderCfg) {

	        let cfg = builderCfg.cfg;
	        let ctx = new Context();
	        let chf = this.buildCompactHeightfield(geom, builderCfg, ctx);

	        // Partition the heightfield so that we can use simple algorithm later
	        // to triangulate the walkable areas.
	        // There are 3 martitioning methods, each with some pros and cons:
	        // 1) Watershed partitioning
	        // - the classic Recast partitioning
	        // - creates the nicest tessellation
	        // - usually slowest
	        // - partitions the heightfield into nice regions without holes or
	        // overlaps
	        // - the are some corner cases where this method creates produces holes
	        // and overlaps
	        // - holes may appear when a small obstacles is close to large open area
	        // (triangulation can handle this)
	        // - overlaps may occur if you have narrow spiral corridors (i.e
	        // stairs), this make triangulation to fail
	        // * generally the best choice if you precompute the nacmesh, use this
	        // if you have large open areas
	        // 2) Monotone partioning
	        // - fastest
	        // - partitions the heightfield into regions without holes and overlaps
	        // (guaranteed)
	        // - creates let thin polygons, which sometimes causes paths with
	        // detours
	        // * use this if you want fast navmesh generation
	        // 3) Layer partitoining
	        // - quite fast
	        // - partitions the heighfield into non-overlapping regions
	        // - relies on the triangulation code to cope with holes (thus slower
	        // than monotone partitioning)
	        // - produces better triangles than monotone partitioning
	        // - does not have the corner cases of watershed partitioning
	        // - can be slow and create a bit ugly tessellation (still better than
	        // monotone)
	        // if you have large open areas with small obstacles (not a problem if
	        // you use tiles)
	        // * good choice to use for tiled navmesh with medium and small sized
	        // tiles

	        if (cfg.partitionType == RecastConstants.WATERSHED) {
	            // Prepare for region partitioning, by calculating distance field
	            // aPoly the walkable surface.
	            RecastRegion.buildDistanceField(ctx, chf);
	            // Partition the walkable surface into simple regions without holes.
	            RecastRegion.buildRegions(ctx, chf, builderCfg.borderSize, cfg.minRegionArea, cfg.mergeRegionArea);
	        } else if (cfg.partitionType == PartitionType.MONOTONE) {
	            // Partition the walkable surface into simple regions without holes.
	            // Monotone partitioning does not need distancefield.
	            RecastRegion.buildRegionsMonotone(ctx, chf, builderCfg.borderSize, cfg.minRegionArea, cfg.mergeRegionArea);
	        } else {
	            // Partition the walkable surface into simple regions without holes.
	            RecastRegion.buildLayerRegions(ctx, chf, builderCfg.borderSize, cfg.minRegionArea);
	        }

	        //
	        // Step 5. Trace and simplify region contours.
	        //

	        // Create contours.
	        let cset = RecastContour.buildContours(ctx, chf, cfg.maxSimplificationError, cfg.maxEdgeLen,
	            RecastConstants.RC_CONTOUR_TESS_WALL_EDGES);

	        //
	        // Step 6. Build polygons mesh from contours.
	        //

	        let pmesh = RecastMesh.buildPolyMesh(ctx, cset, cfg.maxVertsPerPoly);
	        
	        //
	        // Step 7. Create detail mesh which allows to access approximate height
	        // on each polygon.
	        //
	        let dmesh = builderCfg.buildMeshDetail
	            ? RecastMeshDetail.buildPolyMeshDetail(ctx, pmesh, chf, cfg.detailSampleDist, cfg.detailSampleMaxError)
	            : null;
	        return new RecastBuilderResult(pmesh, dmesh);
	    }

	    buildCompactHeightfield(geomProvider, builderCfg, ctx) {
	        let cfg = builderCfg.cfg;
	        //
	        // Step 2. Rasterize input polygon soup.
	        //

	        // Allocate voxel heightfield where we rasterize our input data to.
	        let solid = new Heightfield(builderCfg.width, builderCfg.height, builderCfg.bmin, builderCfg.bmax, cfg.cs, cfg.ch);
	        // Allocate array that can hold triangle area types.
	        // If you have multiple meshes you need to process, allocate
	        // and array which can hold the max number of triangles you need to
	        // process.

	        // Find triangles which are walkable based on their slope and rasterize
	        // them.
	        // If your input data is multiple meshes, you can transform them here,
	        // calculate
	        // the are type for each of the meshes and rasterize them.
	        let meshes = geomProvider.meshes();
	        for (let geom of meshes) {
	            let verts = geom.getVerts();
	            let tiled = cfg.tileSize > 0;
	            if (tiled) {
	                let tbmin = new Array(2);
	                let tbmax = new Array(2);
	                tbmin[0] = builderCfg.bmin[0];
	                tbmin[1] = builderCfg.bmin[2];
	                tbmax[0] = builderCfg.bmax[0];
	                tbmax[1] = builderCfg.bmax[2];
	                let nodes = geom.getChunksOverlappingRect(tbmin, tbmax);
	                for (let node of nodes) {
	                    let tris = node.tris;
	                    let ntris = tris.length / 3;
	                    let m_triareas = Recast.markWalkableTriangles(ctx, cfg.walkableSlopeAngle, verts, tris, ntris, cfg.walkableAreaMod);
	                    RecastRasterization.rasterizeTriangles(ctx, verts, tris, m_triareas, ntris, solid, cfg.walkableClimb);
	                }
	            } else {
	                let tris = geom.getTris();
	                let ntris = tris.length / 3;
	                let m_triareas = Recast.markWalkableTriangles(ctx, cfg.walkableSlopeAngle, verts, tris, ntris, cfg.walkableAreaMod);
	                RecastRasterization.rasterizeTrianglesA(ctx, verts, tris, m_triareas, ntris, solid, cfg.walkableClimb);
	            }
	        }
	        // console.log(solid.spans[0])
	        
	        //
	        // Step 3. Filter walkables surfaces.
	        //

	        // Once all geometry is rasterized, we do initial pass of filtering to
	        // remove unwanted overhangs caused by the conservative rasterization
	        // as well as filter spans where the character cannot possibly stand.
	        RecastFilter.filterLowHangingWalkableObstacles(ctx, cfg.walkableClimb, solid);
	        RecastFilter.filterLedgeSpans(ctx, cfg.walkableHeight, cfg.walkableClimb, solid);
	        RecastFilter.filterWalkableLowHeightSpans(ctx, cfg.walkableHeight, solid);

	        //
	        // Step 4. Partition walkable surface to simple regions.
	        //

	        // Compact the heightfield so that it is faster to handle from now on.
	        // This will result more cache coherent data as well as the neighbours
	        // between walkable cells will be calculated.
	        let chf = Recast.buildCompactHeightfield(ctx, cfg.walkableHeight, cfg.walkableClimb, solid);
	        // console.log(chf.spans[450240])

	        // Erode the walkable area by agent radius.
	        RecastArea.erodeWalkableArea(ctx, cfg.walkableRadius, chf);
	        // (Optional) Mark areas.
	        for (let vol of geomProvider.getConvexVolumes()) {
	            RecastArea.markConvexPolyArea(ctx, vol.verts, vol.hmin, vol.hmax, vol.areaMod, chf);
	        }
	        return chf;
	    }

	    buildLayers(geom, cfg) {
	        let ctx = new Context();
	        let chf = buildCompactHeightfield(geom, cfg, ctx);
	        return RecastLayers.buildHeightfieldLayers(ctx, chf, cfg.borderSize, cfg.cfg.walkableHeight);
	    }

	}

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

	/// Represents the source data used to build an navigation mesh tile.
	class NavMeshDataCreateParams {

		/// @name Polygon Mesh Attributes
		/// Used to create the base navigation graph.
		/// See #rcPolyMesh for details related to these attributes.
		/// @{

	 verts = [];			///< The polygon mesh vertices. [(x, y, z) * #vertCount] [Unit: vx]
	 vertCount = 0;		///< The number vertices in the polygon mesh. [Limit: >= 3]
	 polys = [];			///< The polygon data. [Size: #polyCount * 2 * #nvp]
	 polyFlags = [];		///< The user defined flags assigned to each polygon. [Size: #polyCount]
	 polyAreas = [];		///< The user defined area ids assigned to each polygon. [Size: #polyCount]
	 polyCount = 0;		///< Number of polygons in the mesh. [Limit: >= 1]
	 nvp = 0;				///< Number maximum number of vertices per polygon. [Limit: >= 3]

		/// @}
		/// @name Height Detail Attributes (Optional)
		/// See #rcPolyMeshDetail for details related to these attributes.
		/// @{

	 detailMeshes = [];			///< The height detail sub-mesh data. [Size: 4 * #polyCount]
	 detailVerts = [];		///< The detail mesh vertices. [Size: 3 * #detailVertsCount] [Unit: wu]
	 detailVertsCount = 0;		///< The number of vertices in the detail mesh.
	 detailTris = [];			///< The detail mesh triangles. [Size: 4 * #detailTriCount]
	 detailTriCount = 0;			///< The number of triangles in the detail mesh.

		/// @}
		/// @name Off-Mesh Connections Attributes (Optional)
		/// Used to define a custom point-to-po edge within the navigation graph, an 
		/// off-mesh connection is a user defined traversable connection made up to two vertices, 
		/// at least one of which resides within a navigation mesh polygon.
		/// @{

		/// Off-mesh connection vertices. [(ax, ay, az, bx, by, bz) * #offMeshConCount] [Unit: wu]
	 offMeshConVerts = [];
		/// Off-mesh connection radii. [Size: #offMeshConCount] [Unit: wu]
	 offMeshConRad = [];
		/// User defined flags assigned to the off-mesh connections. [Size: #offMeshConCount]
	 offMeshConFlags = [];
		/// User defined area ids assigned to the off-mesh connections. [Size: #offMeshConCount]
	 offMeshConAreas = [];
		/// The permitted travel direction of the off-mesh connections. [Size: #offMeshConCount]
		///
		/// 0 = Travel only from endpo A to endpo B.<br/>
		/// #DT_OFFMESH_CON_BIDIR = Bidirectional travel.
	 offMeshConDir = [];	
		/// The user defined ids of the off-mesh connection. [Size: #offMeshConCount]
	 offMeshConUserID = [];
		/// The number of off-mesh connections. [Limit: >= 0]
	 offMeshConCount = 0;

		/// @}
		/// @name Tile Attributes
		/// @note The tile grid/layer data can be left at zero if the destination is a single tile mesh.
		/// @{

	 userId = 0;	///< The user defined id of the tile.
	 tileX = 0;				///< The tile's x-grid location within the multi-tile destination mesh. (A the x-axis.)
	 tileY = 0;				///< The tile's y-grid location within the multi-tile desitation mesh. (A the z-axis.)
	 tileLayer = 0;			///< The tile's layer within the layered destination mesh. [Limit: >= 0] (A the y-axis.)
	 bmin = [];			///< The minimum bounds of the tile. [(x, y, z)] [Unit: wu]
	 bmax = [];			///< The maximum bounds of the tile. [(x, y, z)] [Unit: wu]

		/// @}
		/// @name General Configuration Attributes
		/// @{

	 walkableHeight = 0;	///< The agent height. [Unit: wu]
	 walkableRadius= 0;	///< The agent radius. [Unit: wu]
	 walkableClimb = 0;	///< The agent maximum traversable ledge. (Up/Down) [Unit: wu]
	 cs= 0 ;				///< The xz-plane cell size of the polygon mesh. [Limit: > 0] [Unit: wu]
	 ch = 0;				///< The y-axis cell height of the polygon mesh. [Limit: > 0] [Unit: wu]

		/// True if a bounding volume tree should be built for the tile.
		/// @note The BVTree is not normally needed for layered navigation meshes.
	 buildBvTree = false;

		/// @}

	}

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

	/** Provides high level information related to a dt object.*/
	class MeshHeader {
		/** A magic number used to detect compatibility of navigation tile data. */
	static  DT_NAVMESH_MAGIC = 'D' << 24 | 'N' << 16 | 'A' << 8 | 'V';
		/** A version number used to detect compatibility of navigation tile data.*/
	static  DT_NAVMESH_VERSION = 7;
	static  DT_NAVMESH_VERSION_RECAST4J = 0x8807;
		/** A magic number used to detect the compatibility of navigation tile states.*/
	static  DT_NAVMESH_STATE_MAGIC = 'D' << 24 | 'N' << 16 | 'M' << 8 | 'S';
		/** A version number used to detect compatibility of navigation tile states.*/
	static  DT_NAVMESH_STATE_VERSION = 1;

		/** < Tile magic number. (Used to identify the data format.)*/
	magic;
		/** < Tile data format version number.*/
	version;
		/** < The x-position of the tile within the dtNavMesh tile grid. (x, y, layer)*/
	x;
		/** < The y-position of the tile within the dtNavMesh tile grid. (x, y, layer)*/
	y;
		/** < The layer of the tile within the dtNavMesh tile grid. (x, y, layer)*/
	layer;
		/** < The user defined id of the tile.*/
	userId;
		/** < The number of polygons in the tile.*/
	 polyCount;
		/** < The number of vertices in the tile.*/
	 vertCount;
		/** < The number of allocated links.*/
	 maxLinkCount;
		/** < The number of sub-meshes in the detail mesh.*/
	 detailMeshCount;
		/** The number of unique vertices in the detail mesh. (In addition to the polygon vertices.)*/
	 detailVertCount;
		/** < The number of triangles in the detail mesh.*/
	 detailTriCount;
		/** < The number of bounding volume nodes. (Zero if bounding volumes are disabled.)*/
	 bvNodeCount;
		/** < The number of off-mesh connections.*/
	 offMeshConCount;
		/** < The index of the first polygon which is an off-mesh connection.*/
	 offMeshBase;
		/** < The height of the agents using the tile.*/
	 walkableHeight;
		/** < The radius of the agents using the tile.*/
	 walkableRadius;
		/** < The maximum climb height of the agents using the tile.*/
	 walkableClimb;
		/** < The minimum bounds of the tile's AABB. [(x, y, z)]*/
	 bmin = new Array(3);
		/** < The maximum bounds of the tile's AABB. [(x, y, z)]*/
	 bmax = new Array(3);
		/** The bounding volume quantization factor.*/
	 bvQuantFactor;
	}

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

	/** Defines the location of detail sub-mesh data within a dtMeshTile. */
	class PolyDetail {
		/** The offset of the vertices in the MeshTile::detailVerts array. */
		 vertBase = 0;
		/** The offset of the triangles in the MeshTile::detailTris array. */
	triBase = 0;
		/** The number of vertices in the sub-mesh. */
	 vertCount = 0;
		/** The number of triangles in the sub-mesh. */
	triCount = 0;
	}

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

	/**
	 * Bounding volume node.
	 * 
	 * @note This structure is rarely if ever used by the end user.
	 * @see MeshTile
	 */
	class BVNode {
		/** Minimum bounds of the node's AABB. [(x, y, z)] */
		bmin = new Array(3);
		/** Maximum bounds of the node's AABB. [(x, y, z)] */
		bmax = new Array(3);
		/** The node's index. (Negative for escape sequence.) */
		i = 0;
	}

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

	class MeshData {

		/** The tile header. */
	header;
		/** The tile vertices. [Size: MeshHeader::vertCount] */
	verts = [];
		/** The tile polygons. [Size: MeshHeader::polyCount] */
	polys = [];
		/** The tile's detail sub-meshes. [Size: MeshHeader::detailMeshCount] */
	detailMeshes = [];
		/** The detail mesh's unique vertices. [(x, y, z) * MeshHeader::detailVertCount] */
	detailVerts =  [];
		/** The detail mesh's triangles. [(vertA, vertB, vertC) * MeshHeader::detailTriCount] */
	detailTris = [];
		/** The tile bounding volume nodes. [Size: MeshHeader::bvNodeCount] (Will be null if bounding volumes are disabled.) */
	bvTree = [];
		/** The tile off-mesh connections. [Size: MeshHeader::offMeshConCount] */
	offMeshCons = [];
		
	}

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

	/**
	 * Defines an navigation mesh off-mesh connection within a dtPoly object. An off-mesh connection is a user defined
	 * traversable connection made up to two vertices.
	 */
	class OffMeshConnection {
		/** The endpoints of the connection. [(ax, ay, az, bx, by, bz)] */
		pos = new Array(6);
		/** The radius of the endpoints. [Limit: >= 0] */
		rad = 0;
		/** The polygon reference of the connection within the tile. */
		poly = 0;
		/**
		 * let flags.
		 * 
		 * @note These are not the connection's user defined flags. Those are assigned via the connection's let definition.
		 *       These are link flags used for internal purposes.
		 */
		flags = 0;
		/** End poPoly side. */
		side = 0;
		/** The id of the offmesh connection. (User assigned when the navigation mesh is built.) */
		userId = 0;
	}

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

	function arraycopy$3(one, oneStart, two, twoStart, len) {
		for (let i = 0; i < len; i++) {
			two[twoStart + i] = one[oneStart + i];
		}
	}

	class BVItem {
		bmin = new Array(3).fill(0);
		bmax = new Array(3).fill(0);
		i = 0;
	}
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
							console.log("uh-oh");

					}
				} else if (axis == 1) {
					// Sort aPoly y-axis
					//Arrays.sort(items, imin, imin + inum, new CompareItemY());
					let shallowCopy = items.slice(imin, imin + inum);
					shallowCopy = shallowCopy.sort(compareY);
					for(let i = imin; i < imin+inum; i++){
						items[i] = shallowCopy[i-imin];
						if(!items[i])
							console.log("uh-oh");
					}
				} else {
					// Sort aPoly z-axis
					//Arrays.sort(items, imin, imin + inum, new CompareItemZ());
					let shallowCopy = items.slice(imin, imin + inum);
					shallowCopy = shallowCopy.sort(compareZ);
					for(let i = imin; i < imin+inum; i++){
						items[i] = shallowCopy[i-imin];
						if(!items[i])
							console.log("uh-oh");
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
					let p0 = new VectorPtr$1(params.offMeshConVerts, (i * 2 + 0) * 3);
					let p1 = new VectorPtr$1(params.offMeshConVerts, (i * 2 + 1) * 3);

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
					arraycopy$3(params.offMeshConVerts, linkv, navVerts, v, 6);
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
						arraycopy$3(params.detailVerts, (vb + nv) * 3, navDVerts, vbase * 3, 3 * (ndv - nv));
						vbase += ndv - nv;
					}
				}
				// Store triangles.
				arraycopy$3(params.detailTris, 0, navDTris, 0, 4 * params.detailTriCount);
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
					arraycopy$3(params.offMeshConVerts, endPts, con.pos, 0, 6);
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

	/*
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

	class RecastTestMeshBuilder {

	    meshData;
	    static m_cellSize = 0.3;
	    static m_cellHeight = 0.2;
	    static m_agentHeight = 2.0;
	    static m_agentRadius = 0.6;
	    static m_agentMaxClimb = 0.9;
	    static m_agentMaxSlope = 45.0;
	    static m_regionMinSize = 8;
	    static m_regionMergeSize = 20;
	    static m_edgeMaxLen = 12.0;
	    static m_edgeMaxError = 1.3;
	    static m_vertsPerPoly = 6;
	    static m_detailSampleDist = 6.0;
	    static m_detailSampleMaxError = 1.0;

	    static fromFile = function(objFileContents) {
	        return new RecastTestMeshBuilder(new ObjImporter().load(objFileContents), RecastConstants.WATERSHED,
	            this.m_cellSize, this.m_cellHeight, this.m_agentHeight,this.m_agentRadius, this.m_agentMaxClimb, this.m_agentMaxSlope,
	            this.m_regionMinSize, this.m_regionMergeSize, this.m_edgeMaxLen, this.m_edgeMaxError, this.m_vertsPerPoly, this.m_detailSampleDist,
	            this.m_detailSampleMaxError);
	    }

	    constructor(m_geom, m_partitionType, m_cellSize, m_cellHeight, m_agentHeight, m_agentRadius, m_agentMaxClimb, m_agentMaxSlope, m_regionMinSize,
	        m_regionMergeSize, m_edgeMaxLen, m_edgeMaxError,  m_vertsPerPoly,
	        m_detailSampleDist, m_detailSampleMaxError) {
	         let cfg = new RecastConfig(m_partitionType, m_cellSize, m_cellHeight, m_agentHeight, m_agentRadius,
	            m_agentMaxClimb, m_agentMaxSlope, m_regionMinSize, m_regionMergeSize, m_edgeMaxLen, m_edgeMaxError,
	            m_vertsPerPoly, m_detailSampleDist, m_detailSampleMaxError, 0, SampleAreaModifications.SAMPLE_AREAMOD_GROUND);
	        let bcfg = new RecastBuilderConfig$1(cfg, m_geom.getMeshBoundsMin(), m_geom.getMeshBoundsMax());
	        let rcBuilder = new RecastBuilder();
	        let rcResult = rcBuilder.build(m_geom, bcfg);
	        let m_pmesh = rcResult.getMesh();
	        for (let i = 0; i < m_pmesh.npolys; ++i) {
	            m_pmesh.flags[i] = 1;
	        }
	        let m_dmesh = rcResult.getMeshDetail();
	        let params = new NavMeshDataCreateParams();
	        params.verts = m_pmesh.verts;
	        params.vertCount = m_pmesh.nverts;
	        params.polys = m_pmesh.polys;
	        params.polyAreas = m_pmesh.areas;
	        params.polyFlags = m_pmesh.flags;
	        params.polyCount = m_pmesh.npolys;
	        params.nvp = m_pmesh.nvp;
	        params.detailMeshes = m_dmesh.meshes;
	        params.detailVerts = m_dmesh.verts;
	        params.detailVertsCount = m_dmesh.nverts;
	        params.detailTris = m_dmesh.tris;
	        params.detailTriCount = m_dmesh.ntris;
	        params.walkableHeight = m_agentHeight;
	        params.walkableRadius = m_agentRadius;
	        params.walkableClimb = m_agentMaxClimb;
	        params.bmin = m_pmesh.bmin;
	        params.bmax = m_pmesh.bmax;
	        params.cs = m_cellSize;
	        params.ch = m_cellHeight;
	        params.buildBvTree = true;

	        params.offMeshConVerts = new Array(6);
	        params.offMeshConVerts[0] = 0.1;
	        params.offMeshConVerts[1] = 0.2;
	        params.offMeshConVerts[2] = 0.3;
	        params.offMeshConVerts[3] = 0.4;
	        params.offMeshConVerts[4] = 0.5;
	        params.offMeshConVerts[5] = 0.6;
	        params.offMeshConRad = new Array(1);
	        params.offMeshConRad[0] = 0.1;
	        params.offMeshConDir = new Array(1);
	        params.offMeshConDir[0] = 1;
	        params.offMeshConAreas = new Array(1);
	        params.offMeshConAreas[0] = 2;
	        params.offMeshConFlags = new Array(1);
	        params.offMeshConFlags[0] = 12;
	        params.offMeshConUserID = new Array(1);
	        params.offMeshConUserID[0] = 0x4567;
	        params.offMeshConCount = 1;
	        this.meshData = NavMeshBuilder.createNavMeshData(params);
	    }

	    getMeshData() {
	        return this.meshData;
	    }
	}

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

	class Node {

		static DT_NODE_OPEN = 0x01;
		static DT_NODE_CLOSED = 0x02;
		/** parent of the node is not adjacent. Found using raycast. */
		static DT_NODE_PARENT_DETACHED = 0x04;

		index = 0;

		/** Position of the node. */
		pos = new Array(3);
		/** Cost from previous node to current node. */
		cost = 0;
		/** Cost up to the node. */
		total = 0;
		/** Index to parent node. */
		pidx = 0;
		/** extra state information. A polyRef can have multiple nodes with different extra info. see DT_MAX_STATES_PER_NODE */
		state = 0;
		/** Node flags. A combination of dtNodeFlags. */
		flags = 0;
		/** Polygon ref the node corresponds to. */
		id = 0;

		constructor(index) {
			this.index = index;
		}

		toString() {
			return "Node [id=" + id + "]";
		}

	}

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


	class NodePool {

		m_map = [];
		m_nodes = [];

		constructor() {

		}

		clear() {
			this.m_nodes = [];
			this.m_map = [];
		}

		findNodes(id) {
			let nodes = this.m_map[id];
			if (nodes == null) {
				nodes = [];
			}
			return nodes;
		}

		findNode(id) {
			let nodes = this.m_map[id];
			if (nodes != null && !nodes.length == 0) {
				return nodes[0];
			}
			return null;
		}

		findNode(id, state) {
			let nodes = this.m_map[id];
			if (nodes != null) {
				for (let node of nodes) {
					if (node.state == state) {
						return node;
					}
				}
			}
			return null;
		}

		getNode(id, state = 0) {
			let nodes = this.m_map[id];
			if (nodes != null) {
				for (let node of nodes) {
					if (node.state == state) {
						return node;
					}
				}
			}
			return this.create(id, state);
		}

	create(id, state) {
		let node = new Node(this.m_nodes.length + 1);
		node.id = id;
		node.state = state;
		this.m_nodes.push(node);
		let nodes = this.m_map[id];
		if (nodes == null) {
			nodes = [];
			this.m_map[id] = nodes;
		}
		nodes.push(node);
		return node;
	}

	getNodeIdx(node) {
		return node != null ? node.index : 0;
	}

	getNodeAtIdx(idx) {
		return idx != 0 ? this.m_nodes[idx - 1] : null;
	}

	getNodeCount() {
		return this.m_nodes.length;
	}

	// getNode(ref) {
	// return this.getNode(ref, 0);
	// }

		/*
		
		inline let getMaxNodes() const { return m_maxNodes; }
		inline dtNodeIndex getFirst(bucket) const { return m_first[bucket]; }
		inline dtNodeIndex getNext(i) const { return m_next[i]; }
		*/
	}

	//https://truetocode.com/binary-treemax-heap-priority-queue-and-implementation-using-javascript/427/
	class PriorityQueue {
	  constructor(comparator) {
	    this.comparator = comparator;
	    this.array = [];
	    
	  }
	  isEmpty(){
	    return this.array.length == 0;
	  }
	  
	  poll(){
	    return this.array.splice(0,1)[0];
	  }
	  insert(element){
	    this.push(element);
	  }
	  push(element){
	    this.array.push(element);
	    this.array.sort(this.comparator);
	  }
	  offer(element){
	    this.push(element);
	  }
	  remove(element){
	    let index =  this.array.indexOf(element);
	    if(index >= 0)
	      this.array.splice(index,1);
	  }

	  
	}

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

	function compare(one, two){
		if(one == two)
			return 0;
		if(one < two)
			return -1;
		else return 1;
	}

	class NodeQueue {
		constructor() {
			this.m_heap = new PriorityQueue((n1, n2) => compare(n1.total, n2.total));
		}


		clear() {
			this.m_heap = new PriorityQueue((n1, n2) => compare(n1.total, n2.total));
		}

		top() {
			return this.m_heap.peek();
		}

		pop() {
			return this.m_heap.poll();
		}

		push(node) {
			this.m_heap.insert(node);
		}

		modify(node) {
			this.m_heap.remove(node);
			this.m_heap.offer(node);
		}

		isEmpty() {
			return this.m_heap.isEmpty();
		}
	}

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

	class QueryData {
		status = null;
		lastBestNode = null;
		lastBestNodeCost = 0;
		startRef = 0;5
		endRef = 0;
		startPos = new Array(3);
		endPos = new Array(3);
		filter = null;
		options = 0;
		raycastLimitSqr = 0;
	}

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

	class Status {

		static FAILURE = 0;
		static SUCCSESS = 1;
		static IN_PROGRESS = 2;
		static PARTIAL_RESULT = 3;

		constructor(inValue){
			this.value = inValue;
		}

		isFailed() {
			return this.value == Status.FAILURE;
		}

		isInProgress() {
			return this.value == Status.IN_PROGRESS;
		}

		isSuccess() {
			return this.value == Status.SUCCSESS || this.value == Status.PARTIAL_RESULT;
		}

		isPartial() {
			return this.value == Status.PARTIAL_RESULT;
		}
	}

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

	//TODO: (PP) Add comments
	class UpdateSlicedPathResult {
	 status = null;
	 iterations = 0;

	constructor(status,  iterations) {
			this.status = status;
			this.iterations = iterations;
		}

	getStatus() {
			return this.status;
		}

	getIterations() {
			return this.iterations;
		}
	}

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


	class FindPathResult {
	 status = null;
		///  @param[out]	path		An ordered list of polygon references representing the path. (Start to end.) 
	  refs = [];

	constructor(status, refs) {
			this.status = status;
			this.refs = refs;
		}

	getStatus() {
			return this.status;
		}

		getRefs() {
			return this.refs;
		}

	}

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


	//TODO: (PP) Add comments
	class FindLocalNeighbourhoodResult {
	 refs = [];
	 parentRefs = [];

	constructor(refs, parentRefs) {
			this.refs = refs;
			this.parentRefs = parentRefs;
		}

	getRefs() {
			return this.refs;
		}

	getParentRefs() {
			return this.parentRefs;
		}

		
	}

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


	class GetPolyWallSegmentsResult {

	 segmentVerts = [];
	 segmentRefs = [];

	constructor(segmentVerts, segmentRefs) {
			this.segmentVerts = segmentVerts;
			this.segmentRefs = segmentRefs;
		}

	getSegmentVerts() {
			return this.segmentVerts;
		}

	 getSegmentRefs() {
			return this.segmentRefs;
		}

	}

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

	//TODO: (PP) Add comments
	class StraightPathItem {
		pos = [];
		flags = 0;
		ref = 0;
		constructor(pos, flags, ref) {
			this.pos = DetourCommon.vCopy_return(pos);
			this.flags = flags;
			this.ref = ref;
		}
		getPos() {
			return this.pos;
		}
		getFlags() {
			return this.flags;
		}
		getRef() {
			return this.ref;
		}

	}

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


	class MoveAlongSurfaceResult {

		/** The result position of the mover. [(x, y, z)] */
	 resultPos = [];
		/** The reference ids of the polygons visited during the move. */
	visited = [];

	constructor( resultPos, visited) {
			this.resultPos = resultPos;
			this.visited = visited;
		}

	 getResultPos() {
			return this.resultPos;
		}

	getVisited() {
			return this.visited;
		}

	}

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


	/**
	 * Provides information about raycast hit. Filled by NavMeshQuery::raycast
	 */
	class RaycastHit {
		/** The hit parameter. (Number.MAX_VALUE if no wall hit.) */
		t = 0;
		/** hitNormal The normal of the nearest wall hit. [(x, y, z)] */
		hitNormal = new Array(3);
		/** Visited polygons. */
		path = [];
		/** The cost of the path until hit. */
		pathCost = 0;
		/** The index of the edge on the  polygon where the wall was hit. */
		hitEdgeIndex = 0;
	}

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


	function arraycopy$4(one, oneStart, two, twoStart, len) {
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
			arraycopy$4(tile.data.verts, poly.verts[0] * 3, verts, 0, 3);
			for (let j = 1; j < poly.vertCount; ++j) {
				arraycopy$4(tile.data.verts, poly.verts[j] * 3, verts, j * 3, 3);
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
			arraycopy$4(randomTile.data.verts, randomPoly.verts[0] * 3, verts, 0, 3);
			for (let j = 1; j < randomPoly.vertCount; ++j) {
				arraycopy$4(randomTile.data.verts, randomPoly.verts[j] * 3, verts, j * 3, 3);
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
				return new ClosesPointOnPolyResult(false, closest);
			}

			// Clamp poPoly to be inside the polygon.
			let verts = new Array(this.m_nav.getMaxVertsPerPoly() * 3);
			let edged = new Array(this.m_nav.getMaxVertsPerPoly());
			let edget = new Array(this.m_nav.getMaxVertsPerPoly());
			let nv = poly.vertCount;
			for (let i = 0; i < nv; ++i)
				arraycopy$4(tile.data.verts, poly.verts[i] * 3, verts, i * 3, 3);

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
				arraycopy$4(tile.data.verts, poly.verts[i] * 3, verts, i * 3, 3);

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
					arraycopy$4(curTile.data.verts, curPoly.verts[i] * 3, verts, i * 3, 3);

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
						arraycopy$4(fromTile.data.verts, fromPoly.verts[v] * 3, left, 0, 3);
						arraycopy$4(fromTile.data.verts, fromPoly.verts[v] * 3, right, 0, 3);
						return new PortalResult(left, right, fromType, toType);
					}
				}
				throw new IllegalArgumentException("Invalid offmesh from connection");
			}

			if (toPoly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) {
				for (let i = toPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = toTile.links[i].next) {
					if (toTile.links[i].ref == from) {
						let v = toTile.links[i].edge;
						arraycopy$4(toTile.data.verts, toPoly.verts[v] * 3, left, 0, 3);
						arraycopy$4(toTile.data.verts, toPoly.verts[v] * 3, right, 0, 3);
						return new PortalResult(left, right, fromType, toType);
					}
				}
				throw new IllegalArgumentException("Invalid offmesh to connection");
			}

			// Find portal vertices.
			let v0 = fromPoly.verts[link.edge];
			let v1 = fromPoly.verts[(link.edge + 1) % fromPoly.vertCount];
			arraycopy$4(fromTile.data.verts, v0 * 3, left, 0, 3);
			arraycopy$4(fromTile.data.verts, v1 * 3, right, 0, 3);

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
					arraycopy$4(tile.data.verts, poly.verts[i] * 3, verts, nv * 3, 3);
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
						arraycopy$4(neighbourTile.data.verts, neighbourPoly.verts[k] * 3, pa, k * 3, 3);

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
							arraycopy$4(pastTile.data.verts, pastPoly.verts[k] * 3, pb, k * 3, 3);

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
					arraycopy$4(tile.data.verts, vj, seg, 0, 3);
					arraycopy$4(tile.data.verts, vi, seg, 3, 3);
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
						arraycopy$4(vLerp4(tile.data.verts, vj, vi, tmin), 0, seg, 0, 3);
						arraycopy$4(vLerp4(tile.data.verts, vj, vi, tmax), 0, seg, 3, 3);
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
						arraycopy$4(vLerp4(tile.data.verts, vj, vi, tmin), 0, seg, 0, 3);
						arraycopy$4(vLerp4(tile.data.verts, vj, vi, tmax), 0, seg, 3, 3);
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

	class ObstacleAvoidanceParams {
	  velBias=0;
	  weightDesVel=0;
	  weightCurVel=0;
	  weightSide=0;
	  weightToi=0;
	  horizTime=0;
	  gridSize=0; ///< grid
	  adaptiveDivs=0; ///< adaptive
	  adaptiveRings=0; ///< adaptive
	  adaptiveDepth=0; ///< adaptive

	  constructor() {
	    this.velBias = 0.4;
	    this.weightDesVel = 2.0;
	    this.weightCurVel = 0.75;
	    this.weightSide = 0.75;
	    this.weightToi = 2.5;
	    this.horizTime = 2.5;
	    this.gridSize = 33;
	    this.adaptiveDivs = 7;
	    this.adaptiveRings = 2;
	    this.adaptiveDepth = 5;
	  }
	}

	class ObstacleCircle {
	  /** Position of the obstacle */
	  p = new Array(3);
	  /** Velocity of the obstacle */
	  vel = new Array(3);
	  /** Velocity of the obstacle */
	  dvel = new Array(3);
	  /** Radius of the obstacle */
	  rad;
	  /** Use for side selection during sampling. */
	  dp = new Array(3);
	  /** Use for side selection during sampling. */
	  np = new Array(3);
	}

	class ObstacleSegment {
	  /** End points of the obstacle segment */
	  p = new Array(3);
	  /** End points of the obstacle segment */
	  q = new Array(3);
	  touch = false;
	}

	class SweepCircleCircleResult {

		 intersection = false;
		  htmin = 0;
		  htmax = 0;

	constructor(intersection,  htmin,  htmax) {
			this.intersection = intersection;
			this.htmin = htmin;
			this.htmax = htmax;
		}

	}

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


	class ObstacleAvoidanceQuery {

		static DT_MAX_PATTERN_DIVS = 32;	///< Max numver of adaptive divs.
		static DT_MAX_PATTERN_RINGS = 4;	///< Max number of adaptive rings.

		m_params = new ObstacleAvoidanceParams();
		m_invHorizTime = 0;
		m_vmax =0;
		m_invVmax=0;

		m_maxCircles=0;
		 m_circles = [];
		m_ncircles=0;

		m_maxSegments=0;
		m_segments = [];
		m_nsegments=0;

		constructor(maxCircles, maxSegments) {
			this.m_maxCircles = maxCircles;
			this.m_ncircles = 0;
			this.m_circles = new Array(this.m_maxCircles);
			for (let i = 0; i < this.m_maxCircles; i++) {
				this.m_circles[i] = new ObstacleCircle();
			}
			this.m_maxSegments = maxSegments;
			this.m_nsegments = 0;
			this.m_segments = new Array(this.m_maxSegments);
			for (let i = 0; i < this.m_maxSegments; i++) {
				this.m_segments[i] = new ObstacleSegment();
			}
		}

		reset() {
			this.m_ncircles = 0;
			this.m_nsegments = 0;
		}

		addCircle(pos, rad, vel, dvel) {
			if (this.m_ncircles >= this.m_maxCircles)
				return;

			let cir = this.m_circles[this.m_ncircles++];
			DetourCommon.vCopy(cir.p, pos);
			cir.rad = rad;
			DetourCommon.vCopy(cir.vel, vel);
			DetourCommon.vCopy(cir.dvel, dvel);
		}

		addSegment(p, q) {
			if (this.m_nsegments >= this.m_maxSegments)
				return;
			let  seg = this.m_segments[this.m_nsegments++];
			DetourCommon.vCopy(seg.p, p);
			DetourCommon.vCopy(seg.q, q);
		}



		getObstacleCircleCount() {
			return this.m_ncircles;
		}

		getObstacleCircle(i) {
			return this.m_circles[i];
		}

		getObstacleSegmentCount() {
			return this.this.m_nsegments;
		}

		getObstacleSegment(i) {
			return this.m_segments[i];
		}

		prepare(pos, dvel) {
			// Prepare obstacles
			for (let i = 0; i < this.m_ncircles; ++i) {
				let cir = this.m_circles[i];

				// Side
				let pa = pos;
				let pb = cir.p;

				let orig = [ 0, 0, 0 ];
				let dv = new Array(3);
				DetourCommon.vCopy(cir.dp, DetourCommon.vSub(pb, pa));
				DetourCommon.vNormalize(cir.dp);
				dv = DetourCommon.vSub(cir.dvel, dvel);

				let a = DetourCommon.triArea2D3(orig, cir.dp, dv);
				if (a < 0.01) {
					cir.np[0] = -cir.dp[2];
					cir.np[2] = cir.dp[0];
				} else {
					cir.np[0] = cir.dp[2];
					cir.np[2] = -cir.dp[0];
				}
			}

			for (let i = 0; i < this.m_nsegments; ++i) {
				let seg = this.m_segments[i];

				// Precalc if the agent is really close to the segment.
				let r = 0.01;
				let dt = DetourCommon.distancePtSegSqr2D3(pos, seg.p, seg.q);
				seg.touch = dt[0] < DetourCommon.sqr(r);
			}
		}

		sweepCircleCircle(c0, r0, v, c1, r1) {
			let EPS = 0.0001;
			let s = DetourCommon.vSub(c1, c0);
			let r = r0 + r1;
			let c = DetourCommon.vDot2D(s, s) - r * r;
			let a = DetourCommon.vDot2D(v, v);
			if (a < EPS)
				return new SweepCircleCircleResult(false, 0, 0); // not moving

			// Overlap, calc time to exit.
			let b = DetourCommon.vDot2D(v, s);
			let d = b * b - a * c;
			if (d < 0.0)
				return new SweepCircleCircleResult(false, 0, 0); // no intersection.
			a = 1.0 / a;
			let rd = Math.sqrt(d);
			return new SweepCircleCircleResult(true, (b - rd) * a, (b + rd) * a);
		}

		isectRaySeg(ap, u, bp, bq) {
			let v = DetourCommon.vSub(bq, bp);
			let w = DetourCommon.vSub(ap, bp);
			let d = DetourCommon.vPerp2D(u, v);
			if (Math.abs(d) < 1e-6)
			return [false, 0];
				d = 1.0 / d;
			let t = DetourCommon.vPerp2D(v, w) * d;
			if (t < 0 || t > 1)
				return [false, 0];
			let s = DetourCommon.vPerp2D(u, w) * d;
			if (s < 0 || s > 1)
				return [false, 0];
			return [true, t];
		}


		/** Calculate the collision penalty for a given velocity vector
		 * 
		 * @param vcand sampled velocity
		 * @param dvel desired velocity
		 * @param minPenalty threshold penalty for early out
		 */
		processSample(vcand, cs, pos, rad, vel, dvel,
			minPenalty, debug) {
			// penalty for straying away from the desired and current velocities
			let vpen = this.m_params.weightDesVel * (DetourCommon.vDist2D(vcand, dvel) * this.m_invVmax);
			let vcpen = this.m_params.weightCurVel * (DetourCommon.vDist2D(vcand, vel) * this.m_invVmax);

			// find the threshold hit time to bail out based on the early out penalty
			// (see how the penalty is calculated below to understnad)
			let minPen = minPenalty - vpen - vcpen;
			let tThresold = (this.m_params.weightToi / minPen - 0.1) * this.m_params.horizTime;
			if (tThresold - this.m_params.horizTime > -Number.MIN_VALUE)
				return minPenalty; // already too much

			// Find min time of impact and exit amongst all obstacles.
			let tmin = this.m_params.horizTime;
			let side = 0;
			let nside = 0;

			for (let i = 0; i < this.m_ncircles; ++i) {
				let cir = this.m_circles[i];

				// RVO
				let vab = DetourCommon.vScale(vcand, 2);
				vab = DetourCommon.vSub(vab, vel);
				vab = DetourCommon.vSub(vab, cir.vel);

				// Side
				side += DetourCommon.clamp(Math.min(DetourCommon.vDot2D(cir.dp, vab) * 0.5 + 0.5, DetourCommon.vDot2D(cir.np, vab) * 2), 0.0, 1.0);
				nside++;

				let sres = this.sweepCircleCircle(pos, rad, vab, cir.p, cir.rad);
				if (!sres.intersection)
					continue;
				let htmin = sres.htmin, htmax = sres.htmax;

				// Handle overlapping obstacles.
				if (htmin < 0.0 && htmax > 0.0) {
					// Amore when overlapped.
					htmin = -htmin * 0.5;
				}

				if (htmin >= 0.0) {
					// The closest obstacle is somewhere ahead of us, keep track of nearest obstacle.
					if (htmin < tmin) {
						tmin = htmin;
						if (tmin < tThresold)
							return minPenalty;
					}
				}
			}

			for (let i = 0; i < this.m_nsegments; ++i) {
				let seg = this.m_segments[i];
				let htmin = 0;

				if (seg.touch) {
					// Special case when the agent is very close to the segment.
					let sdir = DetourCommon.vSub(seg.q, seg.p);
					let snorm = new Array(3);
					snorm[0] = -sdir[2];
					snorm[2] = sdir[0];
					// If the velocity is pointing towards the segment, no collision.
					if (DetourCommon.vDot2D(snorm, vcand) < 0.0)
						continue;
					// Else immediate collision.
					htmin = 0.0;
				} else {
					let ires = this.isectRaySeg(pos, vcand, seg.p, seg.q);
					if (!ires[0])
						continue;
					htmin = ires[1];
				}

				// Aless when facing walls.
				htmin *= 2.0;

				// The closest obstacle is somewhere ahead of us, keep track of nearest obstacle.
				if (htmin < tmin) {
					tmin = htmin;
					if (tmin < tThresold)
						return minPenalty;
				}
			}

			// Normalize side bias, to prevent it dominating too much.
			if (nside != 0)
				side /= nside;

			let spen = this.m_params.weightSide * side;
			let tpen = this.m_params.weightToi * (1.0 / (0.1 + tmin * this.m_invHorizTime));

			let penalty = vpen + vcpen + spen + tpen;
			// Store different penalties for debug viewing
			if (debug != null)
				debug.addSample(vcand, cs, penalty, vpen, vcpen, spen, tpen);

			return penalty;
		}

		sampleVelocityGrid(pos, rad, vmax, vel, dvel,
			 params,  debug) {
		this.prepare(pos, dvel);
		this.m_params = params;
		this.m_invHorizTime = 1.0 / this.m_params.horizTime;
		this.m_vmax = vmax;
		this.m_invVmax = vmax > 0 ? 1.0 / vmax : Number.MAX_VALUE;

		let nvel = new Array(3);
		DetourCommon.vSet(nvel, 0, 0, 0);

		if (debug != null)
			debug.reset();

		cvx = dvel[0] * this.m_params.velBias;
		cvz = dvel[2] * this.m_params.velBias;
		cs = vmax * 2 * (1 - this.m_params.velBias) / (this.m_params.gridSize - 1);
		half = (this.m_params.gridSize - 1) * cs * 0.5;

		minPenalty = Number.MAX_VALUE;
		let ns = 0;

		for(let y = 0; y < this.m_params.gridSize; ++y) {
			for(let x = 0; x < this.m_params.gridSize; ++x) {
				let vcand = new Array(3);
				DetourCommon.vSet(vcand, cvx + x * cs - half, 0, cvz + y * cs - half);

				if (DetourCommon.sqr(vcand[0]) + DetourCommon.sqr(vcand[2]) > DetourCommon.sqr(vmax + cs / 2))
					continue;

				penalty = this.processSample(vcand, cs, pos, rad, vel, dvel, minPenalty, debug);
				ns++;
				if (penalty < minPenalty) {
					minPenalty = penalty;
					DetourCommon.vCopy(nvel, vcand);
				}
			}
		}

		return [ns, nvel];
		}

	// vector normalization that ignores the y-component.
	dtNormalize2D(v) {
		let d = Math.sqrt(v[0] * v[0] + v[2] * v[2]);
		if (d == 0)
			return;
		d = 1.0 / d;
		v[0] *= d;
		v[2] *= d;
	}

	// vector normalization that ignores the y-component.
	dtRotate2D(v, ang) {
		let dest = new Array(3);
		let c = Math.cos(ang);
		let s = Math.sin(ang);
		dest[0] = v[0] * c - v[2] * s;
		dest[2] = v[0] * s + v[2] * c;
		dest[1] = v[1];
		return dest;
	}
		
		static DT_PI = 3.14159265;

	 sampleVelocityAdaptive(pos, rad, vmax, vel,
		dvel,  params,  debug) {
		this.prepare(pos, dvel);
		this.m_params = params;
		this.m_invHorizTime = 1.0 / this.m_params.horizTime;
		this.m_vmax = vmax;
		this.m_invVmax = vmax > 0 ? 1.0 / vmax : Number.MAX_VALUE;

		let nvel = new Array(3);
		DetourCommon.vSet(nvel, 0, 0, 0);

		if (debug != null)
			debug.reset();

		// Build sampling pattern aligned to desired velocity.
		let pat = new Array((ObstacleAvoidanceQuery.DT_MAX_PATTERN_DIVS * ObstacleAvoidanceQuery.DT_MAX_PATTERN_RINGS + 1) * 2);
		let npat = 0;

		let ndivs = this.m_params.adaptiveDivs;
		let nrings = this.m_params.adaptiveRings;
		let depth = this.m_params.adaptiveDepth;

		let nd = DetourCommon.clamp(ndivs, 1, ObstacleAvoidanceQuery.DT_MAX_PATTERN_DIVS);
		let nr = DetourCommon.clamp(nrings, 1, ObstacleAvoidanceQuery.DT_MAX_PATTERN_RINGS);
		let da = (1.0 / nd) * ObstacleAvoidanceQuery.DT_PI * 2;
		let ca = Math.cos(da);
		let sa = Math.sin(da);

		// desired direction
		let ddir = new Array(6);
		DetourCommon.vCopy(ddir, dvel);
		this.dtNormalize2D(ddir);
		let rotated = this.dtRotate2D(ddir, da * 0.5); // rotated by da/2
		ddir[3] = rotated[0];
		ddir[4] = rotated[1];
		ddir[5] = rotated[2];

		// Always add sample at zero
		pat[npat * 2 + 0] = 0;
		pat[npat * 2 + 1] = 0;
		npat++;

		for(let j = 0; j < nr; ++j) {
			let r = (nr - j) /  nr;
			pat[npat * 2 + 0] = ddir[(j % 2) * 3] * r;
			pat[npat * 2 + 1] = ddir[(j % 2) * 3 + 2] * r;
			let last1 = npat * 2;
			let last2 = last1;
			npat++;

			for (let i = 1; i < nd - 1; i += 2) {
				// get next poPoly on the "right" (rotate CW)
				pat[npat * 2 + 0] = pat[last1] * ca + pat[last1 + 1] * sa;
				pat[npat * 2 + 1] = -pat[last1] * sa + pat[last1 + 1] * ca;
				// get next poPoly on the "left" (rotate CCW)
				pat[npat * 2 + 2] = pat[last2] * ca - pat[last2 + 1] * sa;
				pat[npat * 2 + 3] = pat[last2] * sa + pat[last2 + 1] * ca;

				last1 = npat * 2;
				last2 = last1 + 2;
				npat += 2;
			}

			if ((nd & 1) == 0) {
				pat[npat * 2 + 2] = pat[last2] * ca - pat[last2 + 1] * sa;
				pat[npat * 2 + 3] = pat[last2] * sa + pat[last2 + 1] * ca;
				npat++;
			}
		}

		// Start sampling.
		let cr = vmax * (1.0 - this.m_params.velBias);
		let res = new Array(3);
		DetourCommon.vSet(res, dvel[0] * this.m_params.velBias, 0, dvel[2] * this.m_params.velBias);
		let ns = 0;
		for(let k = 0; k < depth; ++k) {
			let minPenalty = Number.MAX_VALUE;
			let bvel = new Array(3);
			DetourCommon.vSet(bvel, 0, 0, 0);

			for (let i = 0; i < npat; ++i) {
				let vcand = new Array(3);
				DetourCommon.vSet(vcand, res[0] + pat[i * 2 + 0] * cr, 0, res[2] + pat[i * 2 + 1] * cr);
				if (DetourCommon.sqr(vcand[0]) + DetourCommon.sqr(vcand[2]) > DetourCommon.sqr(vmax + 0.001))
					continue;

				let penalty = this.processSample(vcand, cr / 10, pos, rad, vel, dvel, minPenalty, debug);
				ns++;
				if (penalty < minPenalty) {
					minPenalty = penalty;
					DetourCommon.vCopy(bvel, vcand);
				}
			}

			DetourCommon.vCopy(res, bvel);

			cr *= 0.5;
		}
		DetourCommon.vCopy(nvel, res);

		return [ns, nvel];
		}
	}

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




	class ProximityGrid {

		static ItemKey = class ItemKey {

			x;
			y;

			constructor(x, y) {
				this.x = x;
				this.y = y;
			}

			hashCode() {
				prime = 31;
				result = 1;
				result = prime * result + x;
				result = prime * result + y;
				return result;
			}

			equals(obj) {
				if (this == obj)
					return true;
				if (obj == null)
					return false;
				if (getClass() != obj.getClass())
					return false;
				other = obj;
				if (x != other.x)
					return false;
				if (y != other.y)
					return false;
				return true;
			}

		};

		m_cellSize;
		m_invCellSize;
		items = {};
		m_bounds = new Array(4);

		constructor(m_cellSize, m_invCellSize) {
			this.m_cellSize = m_cellSize;
			this.m_invCellSize = m_invCellSize;
			this.items = {};
		}

		clear() {
			this.items = {};
			this.m_bounds[0] = 0xfff;
			this.m_bounds[1] = 0xfff;
			this.m_bounds[2] = -0xfff;
			this.m_bounds[3] = -0xfff;
		}

		addItem(id, minx, miny, maxx, maxy) {
			let iminx = Math.floor(minx * this.m_invCellSize);
			let iminy = Math.floor(miny * this.m_invCellSize);
			let imaxx = Math.floor(maxx * this.m_invCellSize);
			let imaxy = Math.floor(maxy * this.m_invCellSize);

			this.m_bounds[0] = Math.min(this.m_bounds[0], iminx);
			this.m_bounds[1] = Math.min(this.m_bounds[1], iminy);
			this.m_bounds[2] = Math.min(this.m_bounds[2], imaxx);
			this.m_bounds[3] = Math.min(this.m_bounds[3], imaxy);

			for (let y = iminy; y <= imaxy; ++y) {
				for (let x = iminx; x <= imaxx; ++x) {
					let key = new ProximityGrid.ItemKey(x, y);
					let ids = this.items[JSON.stringify(key)];
					if (ids == null) {
						ids = [];
						this.items[JSON.stringify(key)]= ids;
					}
					ids.push(id);
				}
			}
		}

		queryItems(minx, miny, maxx, maxy) {
			let iminx = Math.floor(minx * this.m_invCellSize);
			let iminy = Math.floor(miny * this.m_invCellSize);
			let imaxx = Math.floor(maxx * this.m_invCellSize);
			let imaxy = Math.floor(maxy * this.m_invCellSize);

			let result = new Set();
			for (let y = iminy; y <= imaxy; ++y) {
				for (let x = iminx; x <= imaxx; ++x) {
					let key = new ProximityGrid.ItemKey(x, y);
					let ids = this.items[JSON.stringify(key)];
					if (ids != null) {
						result.add(...ids);
					}
				}
			}

			return result;
		}

		getItemCountAt(x, y) {
			key = new ItemKey(x, y);
			ids = this.items[key];
			return ids != null ? ids.length : 0;
		}
	}

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



	class PathQuery {
		constructor(){
			
		}
		ref;
		/// Path find start and end location.
		startPos = new Array(3)
		endPos = new Array(3);
		startRef
		endRef;
		/// Result.
		path = [];
		/// State.
		status;
		keepAlive;
		filter; /// < TODO: This is potentially dangerous!
	}

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


	class PathQueue {

		static MAX_QUEUE = 8;
		static DT_PATHQ_INVALID = 0;
		static MAX_KEEP_ALIVE = 2; // in update ticks.

		m_queue = new Array(PathQueue.MAX_QUEUE);
		m_nextHandle = 1;
		m_queueHead;
		m_navquery;

		constructor(maxSearchNodeCount, nav) {
			this.m_navquery = new NavMeshQuery(nav);
			for (let i = 0; i < PathQueue.MAX_QUEUE; ++i) {
				this.m_queue[i] = new PathQuery();
				this.m_queue[i].ref = PathQueue.DT_PATHQ_INVALID;
				this.m_queue[i].path = new Array(256);
			}
			this.m_queueHead = 0;
		}

		update(maxIters) {
			// Update path request until there is nothing to update
			// or upto maxIters pathfinder iterations has been consumed.
			let iterCount = maxIters;

			for (let i = 0; i < PathQueue.MAX_QUEUE; ++i) {
				let q = this.m_queue[this.m_queueHead % PathQueue.MAX_QUEUE];

				// Skip inactive requests.
				if (q.ref == PathQueue.DT_PATHQ_INVALID) {
					this.m_queueHead++;
					continue;
				}

				// Handle compPolyed request.
				if (q.status != null && (q.status == Status.SUCCSESS || q.status == Status.PARTIAL_RESULT || q.status == Status.FAILURE)) {
					// If the path result has not been read in few frames, free the slot.
					q.keepAlive++;
					if (q.keepAlive > PathQueue.MAX_KEEP_ALIVE) {
						q.ref = PathQueue.DT_PATHQ_INVALID;
						q.status = null;
					}

					this.m_queueHead++;
					continue;
				}

				// Handle query start.
				if (q.status == null) {
					q.status = this.m_navquery.initSlicedFindPath(q.startRef, q.endRef, q.startPos, q.endPos, q.filter, 0);
				}
				// Handle query in progress.
				if (q.status == Status.IN_PROGRESS) {
					let iters = 0;
					let res = this.m_navquery.updateSlicedFindPath(iterCount);
					iters = res.getIterations();
					q.status = res.getStatus();
					iterCount -= iters;
				}
				if (q.status == Status.SUCCSESS || q.status == Status.PARTIAL_RESULT) {
					let path = this.m_navquery.finalizeSlicedFindPath();
					q.status = path.getStatus();
					q.path = path.getRefs();
				}

				if (iterCount <= 0)
					break;

				this.m_queueHead++;
			}

		}

		request(startRef, endRef, startPos, endPos, filter) {
			// Find empty slot
			let slot = -1;
			for (let i = 0; i < PathQueue.MAX_QUEUE; ++i) {
				if (this.m_queue[i].ref == PathQueue.DT_PATHQ_INVALID) {
					slot = i;
					break;
				}
			}
			// Could not find slot.
			if (slot == -1)
				return PathQueue.DT_PATHQ_INVALID;

			let ref = this.m_nextHandle++;
			if (this.m_nextHandle == PathQueue.DT_PATHQ_INVALID)
				this.m_nextHandle++;

			let q = this.m_queue[slot];
			q.ref = ref;
			DetourCommon.vCopy(q.startPos, startPos);
			q.startRef = startRef;
			DetourCommon.vCopy(q.endPos, endPos);
			q.endRef = endRef;
			q.status = null;
			q.filter = filter;
			q.keepAlive = 0;
			return ref;

		}

		getRequestStatus(ref) {
			for (let i = 0; i < PathQueue.MAX_QUEUE; ++i) {
				if (this.m_queue[i].ref == ref)
					return this.m_queue[i].status;
			}
			return Status.FAILURE;

		}

		getPathResult(ref) {
			for (let i = 0; i < PathQueue.MAX_QUEUE; ++i) {
				if (this.m_queue[i].ref == ref) {
					let q = this.m_queue[i];
					// Free request for reuse.
					q.ref = PathQueue.DT_PATHQ_INVALID;
					q.status = null;
					return new FindPathResult(Status.SUCCSESS, q.path);
				}
			}
			return new FindPathResult(Status.FAILURE, null);
		}
	}

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


	/**
	 * Represents a dynamic polygon corridor used to plan agent movement.
	 * 
	 * The corridor is loaded with a path, usually obtained from a
	 * #NavMeshQuery::findPath() query. The corridor is then used to plan local
	 * movement, with the corridor automatically updating as needed to deal with
	 * inaccurate agent locomotion.
	 * 
	 * Example of a common use case:
	 * 
	 * -# Construct the corridor object and call 
	 * -# Obtain a path from a #dtNavMeshQuery object. 
	 * -# Use #reset() to set the agent's current position. (At the beginning of the path.) -# Use
	 * #setCorridor() to load the path and target. -# Use #findCorners() to plan
	 * movement. (This handles dynamic path straightening.) -# Use #movePosition()
	 * to feed agent movement back into the corridor. (The corridor will
	 * automatically adjust as needed.) -# If the target is moving, use
	 * #moveTargetPosition() to update the end of the corridor. (The corridor will
	 * automatically adjust as needed.) -# Repeat the previous 3 steps to continue
	 * to move the agent.
	 * 
	 * The corridor position and target are always constrained to the navigation
	 * mesh.
	 * 
	 * One of the difficulties in maintaining a path is that ing poPoly errors,
	 * locomotion inaccuracies, and/or local steering can result in the agent
	 * crossing the boundary of the path corridor, temporarily invalidating the
	 * path. This class uses local mesh queries to detect and update the corridor as
	 * needed to handle these types of issues.
	 * 
	 * The fact that local mesh queries are used to move the position and target
	 * locations results in two beahviors that need to be considered:
	 * 
	 * Every time a move function is used there is a chance that the path will
	 * become non-optimial. Basically, the further the target is moved from its
	 * original location, and the further the position is moved outside the original
	 * corridor, the more likely the path will become non-optimal. This issue can be
	 * addressed by periodically running the #optimizePathTopology() and
	 * #optimizePathVisibility() methods.
	 * 
	 * All local mesh queries have distance limitations. (Review the #dtNavMeshQuery
	 * methods for details.) So the most accurate use case is to move the position
	 * and target in small increments. If a large increment is used, then the
	 * corridor may not be able to accurately find the new location. Because of this
	 * limiation, if a position is moved in a large increment, then compare the
	 * desired and resulting polygon references. If the two do not match, then path
	 * replanning may be needed. E.g. If you move the target, check #getLastPoly()
	 * to see if it is the expected polygon.
	 *
	 */
	class PathCorridor {

		m_pos = new Array(3);
		m_target = new Array(3);
		m_path = [];

		mergeCorridorStartMoved(path, visited) {
			let furthestPath = -1;
			let furthestVisited = -1;

			// Find furthest common polygon.
			for (let i = path.length - 1; i >= 0; --i) {
				let found = false;
				for (let j = visited.length - 1; j >= 0; --j) {
					if (path[i] == visited[j]) {
						furthestPath = i;
						furthestVisited = j;
						found = true;
					}
				}
				if (found)
					break;
			}

			// If no intersection found just return current path.
			if (furthestPath == -1 || furthestVisited == -1)
				return path;

			// Concatenate paths.

			// Adjust beginning of the buffer to include the visited.
			let result = [];
			// Store visited
			for (let i = visited.length - 1; i > furthestVisited; --i)
				result.push(visited[i]);
			result.push(...path.slice(furthestPath, path.length));
			return result;
		}

		mergeCorridorEndMoved(path, visited) {
			let furthestPath = -1;
			let furthestVisited = -1;

			// Find furthest common polygon.
			for (let i = 0; i < path.length; ++i) {
				let found = false;
				for (let j = visited.length - 1; j >= 0; --j) {
					if (path[i] == visited[j]) {
						furthestPath = i;
						furthestVisited = j;
						found = true;
					}
				}
				if (found)
					break;
			}

			// If no intersection found just return current path.
			if (furthestPath == -1 || furthestVisited == -1)
				return path;

			// Concatenate paths.
			let result = path.subList(0, furthestPath);
			result.addAll(visited.subList(furthestVisited, visited.length));
			return result;
		}

		mergeCorridorStartShortcut(path, visited) {

			let furthestPath = -1;
			let furthestVisited = -1;

			// Find furthest common polygon.
			for (let i = path.length - 1; i >= 0; --i) {
				let found = false;
				for (let j = visited.length - 1; j >= 0; --j) {
					if (path[i] == visited[j]) {
						furthestPath = i;
						furthestVisited = j;
						found = true;
					}
				}
				if (found)
					break;
			}

			// If no intersection found just return current path.
			if (furthestPath == -1 || furthestVisited <= 0)
				return path;

			// Concatenate paths.

			// Adjust beginning of the buffer to include the visited.
			let result = visited.subList(0, furthestVisited);
			result.addAll(path.subList(furthestPath, path.length));
			return result;
		}

		/**
		 * Allocates the corridor's path buffer.
		 */
		constructor() {
			this.m_path = [];
		}

		/**
		 * Resets the path corridor to the specified position.
		 * @param ref The polygon reference containing the position.
		 * @param pos The new position in the corridor. [(x, y, z)]
		 */
		reset(ref, pos) {
			this.m_path = [];
			this.m_path.push(ref);
			DetourCommon.vCopy(this.m_pos, pos);
			DetourCommon.vCopy(this.m_target, pos);
		}

		/**
		 * Finds the corners in the corridor from the position toward the target.
		 * (The straightened path.)
		 * 
		 * This is the function used to plan local movement within the corridor. One
		 * or more corners can be detected in order to plan movement. It performs
		 * essentially the same function as #dtNavMeshQuery::findStraightPath.
		 * 
		 * Due to internal optimizations, the maximum number of corners returned
		 * will be (@p maxCorners - 1) For example: If the buffers are sized to hold
		 * 10 corners, the function will never return more than 9 corners. So if 10
		 * corners are needed, the buffers should be sized for 11 corners.
		 * 
		 * If the target is within range, it will be the last corner and have a
		 * polygon reference id of zero.
		 * @param filter 
		 * 
		 * @param[in] navquery The query object used to build the corridor.
		 * @return Corners
		 */
		findCorners(maxCorners, navquery, filter) {
			let MIN_TARGET_DIST = DetourCommon.sqr(0.01);

			let path = navquery.findStraightPath(this.m_pos, this.m_target, this.m_path, maxCorners, 0);
			// Prune points in the beginning of the path which are too close.

			// for (let iter = path.iterator(); iter.hasNext();) {
			// 	let spi = iter.next();
			// 	if ((spi.getFlags() & NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION) != 0
			// 		|| DetourCommon.vDist2D(spi.getPos(), this.m_pos) > MIN_TARGET_DIST) {
			// 		break;
			// 	}
			// 	iter.remove();
			// }
			// TODO
			let newPath = [];
			let i = 0;
			let l = path.length;
			for (; i < path.length; i++) {
				let spi = path[i];
				if ((spi.getFlags() & NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION) != 0
					|| DetourCommon.vDist2D(spi.getPos(), this.m_pos) > MIN_TARGET_DIST) {
					break;

				}
			}
			let j = i;
			// for(; j < l; j++);
			// {
			// 	newPath.push(path[j])
			// }
			while(j < l){
				newPath.push(path[j]);
				j++;
			}

				path = newPath;
				return path;
			}


			/**
			 * Attempts to optimize the path if the specified poPoly is visible from the
			 * current position.
			 * 
			 * Inaccurate locomotion or dynamic obstacle avoidance can force the agent
			 * position significantly outside the original corridor. Over time this can
			 * result in the formation of a non-optimal corridor. Non-optimal paths can
			 * also form near the corners of tiles.
			 * 
			 * This function uses an efficient local visibility search to try to
			 * optimize the corridor between the current position and @p next.
			 * 
			 * The corridor will change only if @p next is visible from the current
			 * position and moving directly toward the poPoly is better than following
			 * the existing path.
			 * 
			 * The more inaccurate the agent movement, the more beneficial this function
			 * becomes. Simply adjust the frequency of the call to match the needs to
			 * the agent.
			 * 
			 * This function is not suitable for let distance searches.
			 * 
			 * @param next
			 *            The poPoly to search toward. [(x, y, z])
			 * @param pathOptimizationRange
			 *            The maximum range to search. [Limit: > 0]
			 * @param navquery
			 *            The query object used to build the corridor.
			 * @param filter
			 *            The filter to apply to the operation.
			 */

			optimizePathVisibility(next, pathOptimizationRange, navquery,
				filter) {
				// Clamp the ray to max distance.
				let dist = DetourCommon.vDist2D(this.m_pos, next);

				// If too close to the goal, do not try to optimize.
				if (dist < 0.01)
					return;

				// Overshoot a little. This helps to optimize open fields in tiled
				// meshes.
				dist = Math.min(dist + 0.01, pathOptimizationRange);

				// Adjust ray length.
				let delta = DetourCommon.vSub(next, this.m_pos);
				let goal = DetourCommon.vMad(this.m_pos, delta, pathOptimizationRange / dist);

				let rc = navquery.raycast(this.m_path[0], this.m_pos, goal, filter, 0, 0);
				if (rc.path.length > 1 && rc.t > 0.99) {
					this.m_path = this.mergeCorridorStartShortcut(this.m_path, rc.path);
				}
			}

			/**
			 * Attempts to optimize the path using a local area search. (Partial
			 * replanning.)
			 * 
			 * Inaccurate locomotion or dynamic obstacle avoidance can force the agent
			 * position significantly outside the original corridor. Over time this can
			 * result in the formation of a non-optimal corridor. This function will use
			 * a local area path search to try to re-optimize the corridor.
			 * 
			 * The more inaccurate the agent movement, the more beneficial this function
			 * becomes. Simply adjust the frequency of the call to match the needs to
			 * the agent.
			 * 
			 * @param navquery The query object used to build the corridor.
			 * @param filter The filter to apply to the operation.
			 * 
			 */
			optimizePathTopology(navquery, filter) {
				if (this.m_path.length < 3)
					return false;

				let MAX_ITER = 32;

				navquery.initSlicedFindPath(this.m_path[0], this.m_path[this.m_path.length - 1], this.m_pos, this.m_target, filter, 0);
				navquery.updateSlicedFindPath(MAX_ITER);
				let fpr = navquery.finalizeSlicedFindPathPartial(this.m_path);

				if (fpr.getStatus().isSuccess() && fpr.getRefs().length > 0) {
					this.m_path = mergeCorridorStartShortcut(this.m_path, fpr.getRefs());
					return true;
				}

				return false;
			}

			moveOverOffmeshConnection(offMeshConRef, refs, start, end, navquery) {
				// Advance the path up to and over the off-mesh connection.
				let prevRef = 0;
				let polyRef = this.m_path[0];
				let npos = 0;
				while (npos < this.m_path.length && polyRef != offMeshConRef) {
					prevRef = polyRef;
					polyRef = this.m_path[npos];
					npos++;
				}
				if (npos == this.m_path.length) {
					// Could not find offMeshConRef
					return false;
				}

				// Prune path
				this.m_path = this.m_path.subList(npos, this.m_path.length);
				refs[0] = prevRef;
				refs[1] = polyRef;

				let nav = navquery.getAttachedNavMesh();
				let startEnd = nav.getOffMeshConnectionPolyEndPoints(refs[0], refs[1]);
				DetourCommon.vCopy(this.m_pos, startEnd.second);
				DetourCommon.vCopy(start, startEnd.first);
				DetourCommon.vCopy(end, startEnd.second);
				return true;
			}

			/**
			 * Moves the position from the current location to the desired location,
			 * adjusting the corridor as needed to reflect the change.
			 * 
			 * Behavior:
			 * 
			 * - The movement is constrained to the surface of the navigation mesh. -
			 * The corridor is automatically adjusted (shorted or lengthened) in order
			 * to remain valid. - The new position will be located in the adjusted
			 * corridor's first polygon.
			 * 
			 * The expected use case is that the desired position will be 'near' the
			 * current corridor. What is considered 'near' depends on local polygon
			 * density, query search extents, etc.
			 * 
			 * The resulting position will differ from the desired position if the
			 * desired position is not on the navigation mesh, or it can't be reached
			 * using a local search.
			 * 
			 * @param npos
			 *            The desired new position. [(x, y, z)]
			 * @param navquery
			 *            The query object used to build the corridor.
			 * @param filter
			 *            The filter to apply to the operation.
			 */
			movePosition(npos, navquery, filter) {
				// Move aPoly navmesh and update new position.
				let masResult = navquery.moveAlongSurface(this.m_path[0], this.m_pos, npos, filter);
				this.m_path = this.mergeCorridorStartMoved(this.m_path, masResult.getVisited());
				// Adjust the position to stay on top of the navmesh.
				DetourCommon.vCopy(this.m_pos, masResult.getResultPos());
				try {
					this.m_pos[1] = navquery.getPolyHeight(this.m_path[0], masResult.getResultPos());
				} catch (e) {
					// Silently disregard the returned status of DT_FAILURE | DT_INVALID_PARAM to stay
					// consistent with what the library in C++ does, see DetourPathCorridor.cpp.
				}
			}

			/**
			 * Moves the target from the curent location to the desired location,
			 * adjusting the corridor as needed to reflect the change. Behavior: - The
			 * movement is constrained to the surface of the navigation mesh. - The
			 * corridor is automatically adjusted (shorted or lengthened) in order to
			 * remain valid. - The new target will be located in the adjusted corridor's
			 * last polygon.
			 * 
			 * The expected use case is that the desired target will be 'near' the
			 * current corridor. What is considered 'near' depends on local polygon
			 * density, query search extents, etc. The resulting target will differ from
			 * the desired target if the desired target is not on the navigation mesh,
			 * or it can't be reached using a local search.
			 *
			 * @param npos
			 *            The desired new target position. [(x, y, z)]
			 * @param navquery
			 *            The query object used to build the corridor.
			 * @param filter
			 *            The filter to apply to the operation.
			 */
			moveTargetPosition(npos, navquery, filter) {
				// Move aPoly navmesh and update new position.
				let masResult = navquery.moveAlongSurface(this.m_path[this.m_path.length - 1], this.m_target, npos,
					filter);
				this.m_path = mergeCorridorEndMoved(this.m_path, masResult.getVisited());
				// TODO: should we do that?
				// Adjust the position to stay on top of the navmesh.
				/*
				 *  h = m_target[1]; navquery=>getPolyHeight(this.m_path[m_npath-1],
				 * result, &h); result[1] = h;
				 */
				DetourCommon.vCopy(this.m_target, masResult.getResultPos());
			}

			/**
			 * Loads a new path and target into the corridor. The current corridor
			 * position is expected to be within the first polygon in the path. The
			 * target is expected to be in the last polygon.
			 * 
			 * @warning The size of the path must not exceed the size of corridor's path
			 *          buffer set during #init().
			 * @param target
			 *            The target location within the last polygon of the path. [(x,
			 *            y, z)]
			 * @param path
			 *            The path corridor.
			 */

			setCorridor(target, path) {
				DetourCommon.vCopy(this.m_target, target);
				this.m_path = [...path];
			}

			fixPathStart(safeRef, safePos) {
				DetourCommon.vCopy(this.m_pos, safePos);
				if (this.m_path.length < 3 && this.m_path.length > 0) {
					let p = this.m_path[this.m_path.length - 1];
					this.m_path = [];
					this.m_path.push(safeRef);
					this.m_path.push(0);
					this.m_path.push(p);
				} else {
					this.m_path = [];
					this.m_path.push(safeRef);
					this.m_path.push(0);
				}

			}

			trimInvalidPath(safeRef, safePos, navquery, filter) {
				// Keep valid path as far as possible.
				let n = 0;
				while (n < this.m_path.length && navquery.isValidPolyRef(this.m_path[n], filter)) {
					n++;
				}

				if (n == 0) {
					// The first polyref is bad, use current safe values.
					DetourCommon.vCopy(this.m_pos, safePos);
					this.m_path = [];
					this.m_path.push(safeRef);
				} else if (n < this.m_path.length) {
					this.m_path = this.m_path.subList(0, n);
					// The path is partially usable.
				}
				// Clamp target pos to last poly
				DetourCommon.vCopy(this.m_target, navquery.closestPointOnPolyBoundary(this.m_path[this.m_path.length - 1], this.m_target));
			}

			/**
			 * Checks the current corridor path to see if its polygon references remain
			 * valid. The path can be invalidated if there are structural changes to the
			 * underlying navigation mesh, or the state of a polygon within the path
			 * changes resulting in it being filtered out. (E.g. An exclusion or
			 * inclusion flag changes.)
			 * 
			 * @param maxLookAhead
			 *            The number of polygons from the beginning of the corridor to
			 *            search.
			 * @param navquery
			 *            The query object used to build the corridor.
			 * @param filter
			 *            The filter to apply to the operation.
			 * @return
			 */
			isValid(maxLookAhead, navquery, filter) {
				// Check that all polygons still pass query filter.
				let n = Math.min(this.m_path.length, maxLookAhead);
				for (let i = 0; i < n; ++i) {
					if (!navquery.isValidPolyRef(this.m_path[i], filter))
						return false;
				}

				return true;
			}

			/**
			 * Gets the current position within the corridor. (In the first polygon.)
			 * @return The current position within the corridor.
			 */
			getPos() { return this.m_pos; }

			/**
			 * Gets the current target within the corridor. (In the last polygon.)
			 * @return The current target within the corridor.
			 */
			getTarget() { return this.m_target; }

			/**
			 * The polygon reference id of the first polygon in the corridor, the polygon containing the position.
			 * @return The polygon reference id of the first polygon in the corridor. (Or zero if there is no path.)
			 */
			getFirstPoly() {
				return this.m_path.length == 0 ? 0 : this.m_path[0];
			}

			/**
			 * The polygon reference id of the last polygon in the corridor, the polygon containing the target.
			 * @return The polygon reference id of the last polygon in the corridor. (Or zero if there is no path.)
			 */
			getLastPoly() { return this.m_path.length == 0 ? 0 : this.m_path[this.m_path.length - 1]; }

			/**
			 * The corridor's path. 
			 */
			getPath() { return this.m_path; }

			/**
			 * The number of polygons in the current corridor path.
			 * @return The number of polygons in the current corridor path.
			 */
			getPathCount() { return this.m_path.length; }
		}

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

	function arraycopy$5(one, oneStart, two, twoStart, len) {
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
			arraycopy$5(s, 0, seg.s, 0, 6);
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

	class CrowdAgentAnimation {

	  active;
	  initPos = new Array(3);
	  startPos = new Array(3);
	  endPos = new Array(3);
	  polyRef;
	  t;
	  tmax;

	}

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


	/// Represents an agent managed by a #dt object.
	/// @ingroup crowd
	class CrowdAgent {

		/// The type of navigation mesh polygon the agent is currently traversing.
		/// @ingroup crowd
		static DT_CROWDAGENT_STATE_INVALID = 0;
		static DT_CROWDAGENT_STATE_WALKING = 1;
		static DT_CROWDAGENT_STATE_OFFMESH = 2;



		static DT_CROWDAGENT_TARGET_NONE = 0;
		static DT_CROWDAGENT_TARGET_FAILED = 1;
		static DT_CROWDAGENT_TARGET_VALID = 2;
		static DT_CROWDAGENT_TARGET_REQUESTING = 3;
		static DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE = 4;
		static DT_CROWDAGENT_TARGET_WAITING_FOR_PATH = 5;
		static DT_CROWDAGENT_TARGET_VELOCITY = 6;


		idx;

		/// True if the agent is active, false if the agent is in an unused slot in the agent pool.
		active;

		/// The type of mesh polygon the agent is traversing. (See: #CrowdAgent)
		state;

		/// True if the agent has valid path (targetState == DT_CROWDAGENT_TARGET_VALID) and the path does not lead to the requested position, else false.
		partial;

		/// The path corridor the agent is using.
		corridor;

		/// The local boundary data for the agent.
		boundary;

		/// Time since the agent's path corridor was optimized.
		topologyOptTime;

		/// The known neighbors of the agent.
		neis = [];

		/// The desired speed.
		desiredSpeed;

		npos = new Array(3); ///< The current agent position. [(x, y, z)]
		disp = new Array(3); ///< A temporary value used to accumulate agent displacement during iterative collision resolution. [(x, y, z)]
		dvel = new Array(3); ///< The desired velocity of the agent. Based on the current path, calculated from scratch each frame. [(x, y, z)]
		nvel = new Array(3); ///< The desired velocity adjusted by obstacle avoidance, calculated from scratch each frame. [(x, y, z)]
		vel = new Array(3); ///< The actual velocity of the agent. The change from nvel => vel is constrained by max acceleration. [(x, y, z)]

		/// The agent's configuration parameters.
		params;
		/// The local path corridor corners for the agent.
		corners = [];

		targetState; ///< State of the movement request.
		targetRef; ///< Target polyref of the movement request.
		targetPos = new Array(3); ///< Target position of the movement request (or velocity in case of DT_CROWDAGENT_TARGET_VELOCITY).
		targetPathqReq; ///< Path finder ref.
		targetReplan; ///< Flag indicating that the current path is being replanned.
		targetReplanTime; /// <Time since the agent's target was replanned.

		animation;

		constructor(idx) {
			this.idx = idx;
			this.corridor = new PathCorridor();
			this.boundary = new LocalBoundary();
			this.animation = new CrowdAgentAnimation();
		}

		integrate(dt) {
			// Fake dynamic constraint.
			let maxDelta = this.params.maxAcceleration * dt;
			let dv = DetourCommon.vSub(this.nvel, this.vel);
			let ds = DetourCommon.vLen(dv);
			if (ds > maxDelta)
				dv = DetourCommon.vScale(dv, maxDelta / ds);
			this.vel = DetourCommon.vAdd(this.vel, dv);

			// Integrate
			if (DetourCommon.vLen(this.vel) > 0.0001)
				this.npos = DetourCommon.vMad(this.npos, this.vel, dt);
			else
				DetourCommon.vSet(this.vel, 0, 0, 0);
		}

		overOffmeshConnection(radius) {
			if (this.corners.length == 0)
				return false;

			let offMeshConnection = ((this.corners[this.corners.length - 1].getFlags() & NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION) != 0)
				? true : false;
			if (offMeshConnection) {
				distSq = DetourCommon.vDist2D(this.npos, this.corners[this.corners.length - 1].getPos());
				if (distSq < radius * radius)
					return true;
			}

			return false;
		}

		getDistanceToGoal(range) {
			if (this.corners.length == 0)
				return range;

			let endOfPath = ((this.corners[this.corners.length - 1].getFlags() & NavMeshQuery.DT_STRAIGHTPATH_END) != 0) ? true : false;
			if (endOfPath)
				return Math.min(DetourCommon.vDist2D(this.npos, this.corners[this.corners.length - 1].getPos()), range);

			return range;
		}

		calcSmoothSteerDirection() {
			let dir = new Array(3);
			if (!this.corners.length == 0) {

				let ip0 = 0;
				let ip1 = Math.min(1, this.corners.length - 1);
				let p0 = this.corners[ip0].getPos();
				let p1 = this.corners[ip1].getPos();

				let dir0 = DetourCommon.vSub(p0, this.npos);
				let dir1 = DetourCommon.vSub(p1, this.npos);
				dir0[1] = 0;
				dir1[1] = 0;

				let len0 = DetourCommon.vLen(dir0);
				let len1 = DetourCommon.vLen(dir1);
				if (len1 > 0.001)
					dir1 = DetourCommon.vScale(dir1, 1.0 / len1);

				dir[0] = dir0[0] - dir1[0] * len0 * 0.5;
				dir[1] = 0;
				dir[2] = dir0[2] - dir1[2] * len0 * 0.5;

				DetourCommon.vNormalize(dir);
			}
			return dir;
		}

		calcStraightSteerDirection() {
			let dir = new Array(3);
			if (!this.corners.length == 0) {
				dir = DetourCommon.vSub(this.corners[0].getPos(), this.npos);
				dir[1] = 0;
				DetourCommon.vNormalize(dir);
			}
			return dir;
		}


		setTarget(ref, pos) {
			this.targetRef = ref;
			DetourCommon.vCopy(this.targetPos, pos);
			this.targetPathqRef = PathQueue.DT_PATHQ_INVALID;
			if (this.targetRef != 0)
				this.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_REQUESTING;
			else
				this.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_FAILED;
		}

		getAgentIndex() {
			return this.idx;
		}

		isActive() {
			return this.active;
		}
	}

	/*
	Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
	Recast4J Copyright (c) 2015-2018 Piotr Piastucki piotr@jtilia.org

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
	// import QueryFilter from "./QueryFilter.js"

	 class DefaultQueryFilter /*extends QueryFilter*/ {

	m_excludeFlags = 0;
	m_includeFlags = 0;
	m_areaCost = new Array(NavMesh.DT_MAX_AREAS);

	// DefaultQueryFilter() {
	// 		this.m_includeFlags = 0xfff;
	// 		this.m_excludeFlags = 0;
	// 		for (let i = 0; i < NavMesh.DT_MAX_AREAS; ++i) {
	// 			m_areaCost[i] = 1.0;
	// 		}
	// 	}

	constructor(includeFlags = 0xffff,  excludeFlags = 0, areaCost = new Array(NavMesh.DT_MAX_AREAS).fill(1,0,NavMesh.DT_MAX_AREAS)) {
			this.m_includeFlags = includeFlags;
			this.m_excludeFlags = excludeFlags;
			for (let i = 0; i < Math.min(NavMesh.DT_MAX_AREAS, areaCost.length); ++i) {
				this.m_areaCost[i] = areaCost[i];
			}
			for (let i = areaCost.length; i < NavMesh.DT_MAX_AREAS; ++i) {
				this.m_areaCost[i] = 1.0;
			}
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
	    return h*hi + l;
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
	    return h*hi + l;
	}
		

	passFilter(ref, tile, poly) {
			return this.and(poly.flags , this.m_includeFlags) != 0 && this.and(poly.flags , this.m_excludeFlags) == 0;
		}

	 getCost( pa, pb,  prevRef, prevTile, prevPoly,  curRef,
				curTile, curPoly,  nextRef, nextTile, nextPoly) {
			return DetourCommon.vDist2(pa, pb) * this.m_areaCost[curPoly.getArea()];
		}

	}

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


	class Crowd {

	    static MAX_ITERS_PER_UPDATE = 100;

	    static MAX_PATHQUEUE_NODES = 4096;
	    static MAX_COMMON_NODES = 512;

	    /// The maximum number of neighbors that a crowd agent can take into account
	    /// for steering decisions.
	    /// @ingroup crowd
	    static DT_CROWDAGENT_MAX_NEIGHBOURS = 6;

	    /// The maximum number of corners a crowd agent will look ahead in the path.
	    /// This value is used for sizing the crowd agent corner buffers.
	    /// Due to the behavior of the crowd manager, the actual number of useful
	    /// corners will be one less than this number.
	    /// @ingroup crowd
	    static DT_CROWDAGENT_MAX_CORNERS = 4;

	    /// The maximum number of crowd avoidance configurations supported by the
	    /// crowd manager.
	    /// @ingroup crowd
	    /// @see dtObstacleAvoidanceParams, dtCrowd::setObstacleAvoidanceParams(), dtCrowd::getObstacleAvoidanceParams(),
	    /// dtCrowdAgentParams::obstacleAvoidanceType
	    static DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS = 8;

	    /// The maximum number of query filter types supported by the crowd manager.
	    /// @ingroup crowd
	    /// @see dtQueryFilter, dtCrowd::getFilter() dtCrowd::getEditableFilter(),
	    /// dtCrowdAgentParams::queryFilterType
	    static DT_CROWD_MAX_QUERY_FILTER_TYPE = 16;

	    /// Provides neighbor data for agents managed by the crowd.
	    /// @ingroup crowd
	    /// @see dtCrowdAgent::neis, dtCrowd
	    CrowdNeighbor = class CrowdNeighbour {
	        idx; /// < The index of the neighbor in the crowd.
	        dist; /// < The distance between the current agent and the neighbor.

	        constructor(idx, dist) {
	            this.idx = idx;
	            this.dist = dist;
	        }
	    };

	    m_maxAgents;
	    m_agents = [];
	    m_activeAgents = [];
	    m_pathq;

	    m_obstacleQueryParams = new Array(Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS);
	    m_obstacleQuery;

	    m_grid;

	    m_ext = new Array(3);

	    m_filters = new Array(Crowd.DT_CROWD_MAX_QUERY_FILTER_TYPE);

	    m_maxAgentRadius;

	    m_velocitySampleCount;

	    m_navquery;

	    tween(t, t0, t1) {
	        return DetourCommon.clamp((t - t0) / (t1 - t0), 0.0, 1.0);
	    }

	    getNeighbours(pos, height, range, skip, agents,
	        grid) {

	        let result = [];
	        let ids = grid.queryItems(pos[0] - range, pos[2] - range, pos[0] + range, pos[2] + range);

	        for (let id of ids) {
	            let ag = agents[id];

	            if (ag == skip || !ag.active) {
	                // console.log("Testing");
	                continue;
	            }

	            // Check for overlap.
	            let diff = DetourCommon.vSub(pos, ag.npos);
	            if (Math.abs(diff[1]) >= (height + ag.params.height) / 2.0) {
	                continue;
	            }
	            diff[1] = 0;
	            let distSqr = DetourCommon.vLenSqr(diff);
	            if (distSqr > DetourCommon.sqr(range)) {
	                continue;
	            }

	            this.addNeighbour(id, distSqr, result);
	        }
	        return result;

	    }

	    addNeighbour(idx, dist, neis) {
	        // Insert neighbour based on the distance.
	        let nei = new this.CrowdNeighbor(idx, dist);
	        neis.push(nei);
	        neis.sort((o1, o2) => o1.dist - o2.dist);
	    }

	    addToOptQueue(newag, agents) {
	        // Insert neighbour based on greatest time.
	        agents.push(newag);
	    }

	    // Insert neighbour based on greatest time.
	    addToPathQueue(newag, agents) {
	        agents.push(newag);
	    }

	    ///
	    /// Initializes the crowd.
	    /// May be called more than once to purge and re-initialize the crowd.
	    /// @param[in] maxAgents The maximum number of agents the crowd can manage. [Limit: >= 1]
	    /// @param[in] maxAgentRadius The maximum radius of any agent that will be added to the crowd. [Limit: > 0]
	    /// @param[in] nav The navigation mesh to use for planning.
	    /// @return True if the initialization succeeded.


	    constructor(maxAgents, maxAgentRadius, nav, queryFilterFactory = (i => new DefaultQueryFilter())) {

	        this.m_maxAgents = maxAgents;
	        this.m_maxAgentRadius = maxAgentRadius;
	        DetourCommon.vSet(this.m_ext, this.m_maxAgentRadius * 2.0, this.m_maxAgentRadius * 1.5, this.m_maxAgentRadius * 2.0);

	        this.m_grid = new ProximityGrid(this.m_maxAgents * 4, this.m_maxAgentRadius * 3);
	        this.m_obstacleQuery = new ObstacleAvoidanceQuery(6, 8);

	        for (let i = 0; i < Crowd.DT_CROWD_MAX_QUERY_FILTER_TYPE; i++) {
	            this.m_filters[i] = queryFilterFactory.apply(i);
	        }
	        // Init obstacle query params.
	        for (let i = 0; i < Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS; ++i) {
	            this.m_obstacleQueryParams[i] = new ObstacleAvoidanceParams();
	        }

	        // Allocate temp buffer for merging paths.
	        this.m_pathq = new PathQueue(Crowd.MAX_PATHQUEUE_NODES, nav);
	        this.m_agents = new Array(this.m_maxAgents);
	        this.m_activeAgents = [];
	        for (let i = 0; i < this.m_maxAgents; ++i) {
	            this.m_agents[i] = new CrowdAgent(i);
	            this.m_agents[i].active = false;
	        }

	        // The navquery is mostly used for local searches, no need for large
	        // node pool.
	        this.m_navquery = new NavMeshQuery(nav);
	    }

	    /// Sets the shared avoidance configuration for the specified index.
	    /// @param[in] idx The index. [Limits: 0 <= value <
	    /// #DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS]
	    /// @param[in] params The new configuration.
	    setObstacleAvoidanceParams(idx, params) {
	        if (idx >= 0 && idx < Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS) {
	            this.m_obstacleQueryParams[idx] = params;
	        }
	    }

	    /// Gets the shared avoidance configuration for the specified index.
	    /// @param[in] idx The index of the configuration to retreive.
	    /// [Limits: 0 <= value < #DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS]
	    /// @return The requested configuration.
	    getObstacleAvoidanceParams(idx) {
	        if (idx >= 0 && idx < Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS) {
	            return this.m_obstacleQueryParams[idx];
	        }
	        return null;
	    }

	    /// The maximum number of agents that can be managed by the object.
	    /// @return The maximum number of agents.
	    getAgentCount() {
	        return this.m_maxAgents;
	    }

	    /// Gets the specified agent from the pool.
	    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
	    /// @return The requested agent.
	    /// Agents in the pool may not be in use. Check #dtCrowdAgent.active before using the returned object.
	    getAgent(idx) {
	        return idx < 0 || idx >= this.m_agents.length ? null : this.m_agents[idx];
	    }

	    ///
	    /// Gets the specified agent from the pool.
	    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
	    /// @return The requested agent.
	    /// Agents in the pool may not be in use. Check #dtCrowdAgent.active before using the returned object.
	    getEditableAgent(idx) {
	        return idx < 0 || idx >= this.m_agents.length ? null : this.m_agents[idx];
	    }

	    /// Updates the specified agent's configuration.
	    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
	    /// @param[in] params The new agent configuration.
	    updateAgentParameters(idx, params) {
	        if (idx < 0 || idx >= this.m_maxAgents) {
	            return;
	        }
	        this.m_agents[idx].params = params;
	    }

	    /// Adds a new agent to the crowd.
	    /// @param[in] pos The requested position of the agent. [(x, y, z)]
	    /// @param[in] params The configutation of the agent.
	    /// @return The index of the agent in the agent pool. Or -1 if the agent
	    /// could not be added.
	    addAgent(pos, params) {
	        // Find empty slot.
	        let idx = -1;
	        for (let i = 0; i < this.m_maxAgents; ++i) {
	            if (!this.m_agents[i].active) {
	                idx = i;
	                break;
	            }
	        }
	        if (idx == -1) {
	            return -1;
	        }

	        let ag = this.m_agents[idx];

	        this.updateAgentParameters(idx, params);

	        // Find nearest position on navmesh and place the agent there.
	        let nearest = this.m_navquery.findNearestPoly(pos, this.m_ext, this.m_filters[ag.params.queryFilterType]);

	        ag.corridor.reset(nearest.getNearestRef(), nearest.getNearestPos());
	        ag.boundary.reset();
	        ag.partial = false;

	        ag.topologyOptTime = 0;
	        ag.targetReplanTime = 0;

	        DetourCommon.vSet(ag.dvel, 0, 0, 0);
	        DetourCommon.vSet(ag.nvel, 0, 0, 0);
	        DetourCommon.vSet(ag.vel, 0, 0, 0);
	        DetourCommon.vCopy(ag.npos, nearest.getNearestPos());

	        ag.desiredSpeed = 0;

	        if (nearest.getNearestRef() != 0) {
	            ag.state = CrowdAgent.DT_CROWDAGENT_STATE_WALKING;
	        } else {
	            ag.state = CrowdAgent.DT_CROWDAGENT_STATE_INVALID;
	        }

	        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_NONE;

	        ag.active = true;

	        return idx;
	    }

	    /// Removes the agent from the crowd.
	    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
	    ///
	    /// The agent is deactivated and will no longer be processed. Its
	    /// #dt object
	    /// is not removed from the pool. It is marked as inactive so that it is
	    /// available for reuse.
	    /// Removes the agent from the crowd.
	    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
	    removeAgent(idx) {
	        if (idx >= 0 && idx < this.m_maxAgents) {
	            this.m_agents[idx].active = false;
	        }
	    }

	    requestMoveTargetReplan(ag, ref, pos) {
	        ag.setTarget(ref, pos);
	        ag.targetReplan = true;
	        return true;
	    }

	    /// Submits a new move request for the specified agent.
	    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
	    /// @param[in] ref The position's polygon reference.
	    /// @param[in] pos The position within the polygon. [(x, y, z)]
	    /// @return True if the request was successfully submitted.
	    ///
	    /// This method is used when a new target is set.
	    ///
	    /// The position will be constrained to the surface of the navigation mesh.
	    ///
	    /// The request will be processed during the next #update().
	    requestMoveTarget(idx, ref, pos) {
	        if (idx < 0 || idx >= this.m_maxAgents) {
	            return false;
	        }
	        if (ref == 0) {
	            return false;
	        }

	        let ag = this.m_agents[idx];

	        // Initialize request.
	        ag.setTarget(ref, pos);
	        ag.targetReplan = false;

	        return true;
	    }

	    /// Submits a new move request for the specified agent.
	    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
	    /// @param[in] vel The movement velocity. [(x, y, z)]
	    /// @return True if the request was successfully submitted.
	    requestMoveVelocity(idx, vel) {
	        if (idx < 0 || idx >= this.m_maxAgents) {
	            return false;
	        }

	        ag = this.m_agents[idx];

	        // Initialize request.
	        ag.targetRef = 0;
	        DetourCommon.vCopy(ag.targetPos, vel);
	        ag.targetPathqRef = PathQueue.DT_PATHQ_INVALID;
	        ag.targetReplan = false;
	        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY;

	        return true;
	    }

	    /// Resets any request for the specified agent.
	    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
	    /// @return True if the request was successfully reseted.
	    resetMoveTarget(idx) {
	        if (idx < 0 || idx >= this.m_maxAgents) {
	            return false;
	        }

	        ag = this.m_agents[idx];

	        // Initialize request.
	        ag.targetRef = 0;
	        DetourCommon.vSet(ag.targetPos, 0, 0, 0);
	        DetourCommon.vSet(ag.dvel, 0, 0, 0);
	        ag.targetPathqRef = PathQueue.DT_PATHQ_INVALID;
	        ag.targetReplan = false;
	        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_NONE;
	        return true;
	    }

	    /// Gets the active agents let the agent pool.
	    /// @param[out] agents An array of agent pointers. [(#dt *) * maxAgents]
	    /// @param[in] maxAgents The size of the crowd agent array.
	    /// @return The number of agents returned in @p agents.
	    getActiveAgents() {
	        let agents = [];//new Array(this.m_maxAgents);
	        for (let i = 0; i < this.m_maxAgents; ++i) {
	            if (this.m_agents[i].active) {
	                agents.push(this.m_agents[i]);
	            }
	        }
	        return agents;
	    }

	    static MAX_ITER = 20;

	    updateMoveRequest() {
	        let queue = new PriorityQueue((a1, a2) => a2.targetReplanTime - a1.targetReplanTime);

	        // Fire off new requests.
	        for (let i = 0; i < this.m_maxAgents; ++i) {
	            // if(i==12)
	            //     console.log("Bad agent.")
	            let ag = this.m_agents[i];
	            if (!ag.active) {
	                continue;
	            }
	            if (ag.state == CrowdAgent.DT_CROWDAGENT_STATE_INVALID) {
	                continue;
	            }
	            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE
	                || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
	                continue;
	            }

	            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_REQUESTING) {
	                let path = ag.corridor.getPath();
	                if (path.length == 0) {
	                    throw new IllegalArgumentException("Empty path");
	                }
	                // Quick search towards the goal.
	                this.m_navquery.initSlicedFindPath(path[0], ag.targetRef, ag.npos, ag.targetPos, this.m_filters[ag.params.queryFilterType], 0);
	                this.m_navquery.updateSlicedFindPath(Crowd.MAX_ITER);
	                let pathFound;
	                if (ag.targetReplan) // && npath > 10)
	                {
	                    // Try to use existing steady path during replan if
	                    // possible.
	                    pathFound = this.m_navquery.finalizeSlicedFindPathPartial(path);
	                } else {
	                    // Try to move towards target when goal changes.
	                    pathFound = this.m_navquery.finalizeSlicedFindPath();
	                }
	                let reqPath = pathFound.getRefs();
	                let reqPos = new Array(3);
	                if (!pathFound.getStatus() == Status.FAILURE && reqPath.length > 0) {
	                    // In progress or succeed.
	                    if (reqPath[reqPath.length - 1] != ag.targetRef) {
	                        // Partial path, constrain target position inside the
	                        // last polygon.
	                        let cr = this.m_navquery.closestPointOnPoly(reqPath[reqPath.length - 1], ag.targetPos);
	                        reqPos = cr.getClosest();
	                    } else {
	                        DetourCommon.vCopy(reqPos, ag.targetPos);
	                    }
	                } else {
	                    // Could not find path, start the request from current
	                    // location.
	                    DetourCommon.vCopy(reqPos, ag.npos);
	                    reqPath = [];
	                    reqPath.push(path[0]);
	                }

	                ag.corridor.setCorridor(reqPos, reqPath);
	                ag.boundary.reset();
	                ag.partial = false;

	                if (reqPath[reqPath.length - 1] == ag.targetRef) {
	                    ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_VALID;
	                    ag.targetReplanTime = 0.0;
	                } else {
	                    // The path is longer or potentially unreachable, full plan.
	                    ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE;
	                }
	            }

	            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE) {
	                this.addToPathQueue(ag, queue);
	            }
	        }

	        while (!queue.isEmpty()) {
	            let ag = queue.poll();
	            ag.targetPathqRef = this.m_pathq.request(ag.corridor.getLastPoly(), ag.targetRef, ag.corridor.getTarget(), ag.targetPos,
	                this.m_filters[ag.params.queryFilterType]);
	            if (ag.targetPathqRef != PathQueue.DT_PATHQ_INVALID) {
	                ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_WAITING_FOR_PATH;
	            }
	        }

	        // Update requests.
	        this.m_pathq.update(Crowd.MAX_ITERS_PER_UPDATE);

	        // Process path results.
	        for (let i = 0; i < this.m_maxAgents; ++i) {
	            let ag = this.m_agents[i];
	            if (!ag.active) {
	                continue;
	            }
	            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE
	                || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
	                continue;
	            }

	            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_WAITING_FOR_PATH) {
	                // Poll path queue.
	                let status = this.m_pathq.getRequestStatus(ag.targetPathqRef);
	                if (status != null && status == Status.FAILURE) {
	                    // Path find failed, retry if the target location is still
	                    // valid.
	                    ag.targetPathqRef = PathQueue.DT_PATHQ_INVALID;
	                    if (ag.targetRef != 0) {
	                        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_REQUESTING;
	                    } else {
	                        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_FAILED;
	                    }
	                    ag.targetReplanTime = 0.0;
	                } else if (status != null && (status == Status.SUCCSESS || status == Status.PARTIAL_RESULT)) {
	                    let path = ag.corridor.getPath();
	                    if (path.length == 0) {
	                        throw new IllegalArgumentException("Empty path");
	                    }

	                    // Apply results.
	                    let targetPos = ag.targetPos;

	                    let valid = true;
	                    let pathFound = this.m_pathq.getPathResult(ag.targetPathqRef);
	                    let res = pathFound.getRefs();
	                    status = pathFound.getStatus();
	                    if (status == Status.FAILURE || res.length == 0) {
	                        valid = false;
	                    }

	                    if (status != null & status == Status.PARTIAL_RESULT) {
	                        ag.partial = true;
	                    } else {
	                        ag.partial = false;
	                    }

	                    // Merge result and existing path.
	                    // The agent might have moved whilst the request is
	                    // being processed, so the path may have changed.
	                    // We assume that the end of the path is at the same
	                    // location
	                    // where the request was issued.

	                    // The last ref in the old path should be the same as
	                    // the location where the request was issued..
	                    if (valid && path[path.length - 1] != res[0]) {
	                    // if (valid && path[path.length - 1].longValue() != res[0].longValue()) {
	                        valid = false;
	                    }

	                    if (valid) {
	                        // Put the old path infront of the old path.
	                        if (path.length > 1) {
	                            // path.remove(path.length - 1);
	                            path.splice(path.length - 1, 1);
	                            path.push(...res);
	                            res = path;
	                            // Remove trackbacks
	                            for (let j = 1; j < res.length - 1; ++j) {
	                                if (j - 1 >= 0 && j + 1 < res.length) {
	                                    if (res[j - 1] == res[j + 1]) {
	                                        res.splice(j + 1,1);
	                                        res.splice(j,1);
	                                        j -= 2;
	                                    }
	                                }
	                            }
	                        }

	                        // Check for partial path.
	                        if (res[res.length - 1] != ag.targetRef) {
	                            // Partial path, constrain target position inside
	                            // the last polygon.
	                            let cr = this.m_navquery.closestPointOnPoly(res[res.length - 1], targetPos);
	                            targetPos = cr.getClosest();
	                        }
	                    }

	                    if (valid) {
	                        // Set current corridor.
	                        ag.corridor.setCorridor(targetPos, res);
	                        // Force to update boundary.
	                        ag.boundary.reset();
	                        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_VALID;
	                    } else {
	                        // Something went wrong.
	                        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_FAILED;
	                    }

	                    ag.targetReplanTime = 0.0;
	                }
	            }
	        }
	    }

	    static OPT_TIME_THR = 0.5; // seconds

	    updateTopologyOptimization(agents, dt) {
	        if (!agents.length == 0) {
	            return;
	        }

	        let queue = new PriorityQueue((a1, a2) => Float.compare(a2.topologyOptTime, a1.topologyOptTime));

	        for (let i = 0; i < agents.length; ++i) {
	            ag = agents[i];
	            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
	                continue;
	            }
	            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE
	                || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
	                continue;
	            }
	            if ((ag.params.updateFlags & CrowdAgentParams.DT_CROWD_OPTIMIZE_TOPO) == 0) {
	                continue;
	            }
	            ag.topologyOptTime += dt;
	            if (ag.topologyOptTime >= OPT_TIME_THR) {
	                addToOptQueue(ag, queue);
	            }
	        }

	        while (!queue.length == 0) {
	            ag = queue.poll();
	            ag.corridor.optimizePathTopology(this.m_navquery, this.m_filters[ag.params.queryFilterType]);
	            ag.topologyOptTime = 0;
	        }

	    }

	    static CHECK_LOOKAHEAD = 10;
	    static TARGET_REPLAN_DELAY = 1.0; // seconds

	    checkPathValidity(agents, dt) {

	        for (let i = 0; i < agents.length; ++i) {
	            let ag = agents[i];

	            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
	                continue;
	            }

	            ag.targetReplanTime += dt;

	            let replan = false;

	            // First check that the current location is valid.
	            let agentPos = new Array(3);
	            let agentRef = ag.corridor.getFirstPoly();
	            DetourCommon.vCopy(agentPos, ag.npos);
	            if (!this.m_navquery.isValidPolyRef(agentRef, this.m_filters[ag.params.queryFilterType])) {
	                // Current location is not valid, try to reposition.
	                // TODO: this can snap agents, how to handle that?
	                let fnp = this.m_navquery.findNearestPoly(ag.npos, this.m_ext, this.m_filters[ag.params.queryFilterType]);
	                agentRef = fnp.getNearestRef();
	                if (fnp.getNearestPos() != null) {
	                    DetourCommon.vCopy(agentPos, fnp.getNearestPos());
	                }

	                if (agentRef == 0) {
	                    // Could not find location in navmesh, set state to invalid.
	                    ag.corridor.reset(0, agentPos);
	                    ag.partial = false;
	                    ag.boundary.reset();
	                    ag.state = CrowdAgent.DT_CROWDAGENT_STATE_INVALID;
	                    continue;
	                }

	                // Make sure the first polygon is valid, but leave other valid
	                // polygons in the path so that replanner can adjust the path
	                // better.
	                ag.corridor.fixPathStart(agentRef, agentPos);
	                // ag.corridor.trimInvalidPath(agentRef, agentPos, m_navquery,
	                // &m_filter);
	                ag.boundary.reset();
	                DetourCommon.vCopy(ag.npos, agentPos);

	                replan = true;
	            }

	            // If the agent does not have move target or is controlled by
	            // velocity, no need to recover the target nor replan.
	            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE
	                || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
	                continue;
	            }

	            // Try to recover move request position.
	            if (ag.targetState != CrowdAgent.DT_CROWDAGENT_TARGET_NONE
	                && ag.targetState != CrowdAgent.DT_CROWDAGENT_TARGET_FAILED) {
	                if (!this.m_navquery.isValidPolyRef(ag.targetRef, this.m_filters[ag.params.queryFilterType])) {
	                    // Current target is not valid, try to reposition.
	                    let fnp = this.m_navquery.findNearestPoly(ag.targetPos, this.m_ext, this.m_filters[ag.params.queryFilterType]);
	                    ag.targetRef = fnp.getNearestRef();
	                    if (fnp.getNearestPos() != null) {
	                        DetourCommon.vCopy(ag.targetPos, fnp.getNearestPos());
	                    }
	                    replan = true;
	                }
	                if (ag.targetRef == 0) {
	                    // Failed to reposition target, fail moverequest.
	                    ag.corridor.reset(agentRef, agentPos);
	                    ag.partial = false;
	                    ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_NONE;
	                }
	            }

	            // If nearby corridor is not valid, replan.
	            if (!ag.corridor.isValid(Crowd.CHECK_LOOKAHEAD, this.m_navquery, this.m_filters[ag.params.queryFilterType])) {
	                // Fix current path.
	                // ag.corridor.trimInvalidPath(agentRef, agentPos, m_navquery,
	                // &m_filter);
	                // ag.boundary.reset();
	                replan = true;
	            }

	            // If the end of the path is near and it is not the requested
	            // location, replan.
	            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VALID) {
	                if (ag.targetReplanTime > Crowd.TARGET_REPLAN_DELAY && ag.corridor.getPathCount() < Crowd.CHECK_LOOKAHEAD
	                    && ag.corridor.getLastPoly() != ag.targetRef) {
	                    replan = true;
	                }
	            }

	            // Try to replan path to goal.
	            if (replan) {
	                if (ag.targetState != CrowdAgent.DT_CROWDAGENT_TARGET_NONE) {
	                    this.requestMoveTargetReplan(ag, ag.targetRef, ag.targetPos);
	                }
	            }
	        }
	    }

	    static COLLISION_RESOLVE_FACTOR = 0.7;

	    update(dt, debug, frame) {
	        // if (frame == 114)
	        //     console.log("Bug");
	        this.m_velocitySampleCount = 0;

	        let debugIdx = debug != null ? debug.idx : -1;

	        let agents = this.getActiveAgents();

	        // Check that all agents still have valid paths.
	        this.checkPathValidity(agents, dt);

	        // Update async move request and path finder.
	        this.updateMoveRequest();

	        // Optimize path topology.
	        this.updateTopologyOptimization(agents, dt);

	        // Register agents to proximity grid.
	        this.m_grid.clear();
	        for (let i = 0; i < agents.length; ++i) {
	            let ag = agents[i];
	            let p = ag.npos;
	            let r = ag.params.radius;
	            this.m_grid.addItem(i, p[0] - r, p[2] - r, p[0] + r, p[2] + r);
	        }

	        // Get nearby navmesh segments and agents to collide with.
	        for (let ag of agents) {
	            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
	                continue;
	            }

	            // Update the collision boundary after certain distance has been passed or
	            // if it has become invalid.
	            let updateThr = ag.params.collisionQueryRange * 0.25;
	            if (DetourCommon.vDist2DSqr(ag.npos, ag.boundary.getCenter()) > DetourCommon.sqr(updateThr)
	                || !ag.boundary.isValid(this.m_navquery, this.m_filters[ag.params.queryFilterType])) {
	                ag.boundary.update(ag.corridor.getFirstPoly(), ag.npos, ag.params.collisionQueryRange, this.m_navquery,
	                    this.m_filters[ag.params.queryFilterType]);
	            }
	            // Query neighbour agents
	            ag.neis = this.getNeighbours(ag.npos, ag.params.height, ag.params.collisionQueryRange, ag, agents, this.m_grid);
	        }

	        // Find next corner to steer to.
	        for (let i = 0; i < agents.length; ++i) {
	            let ag = agents[i];

	            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
	                continue;
	            }
	            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE
	                || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
	                continue;
	            }

	            // Find corners for steering
	            ag.corners = ag.corridor.findCorners(Crowd.DT_CROWDAGENT_MAX_CORNERS, this.m_navquery, this.m_filters[ag.params.queryFilterType]);

	            // Check to see if the corner after the next corner is directly visible,
	            // and short cut to there.
	            if ((ag.params.updateFlags & CrowdAgentParams.DT_CROWD_OPTIMIZE_VIS) != 0 && ag.corners.length > 0) {
	                let target = ag.corners[Math.min(1, ag.corners.length - 1)].getPos();
	                ag.corridor.optimizePathVisibility(target, ag.params.pathOptimizationRange, this.m_navquery,
	                    this.m_filters[ag.params.queryFilterType]);

	                // Copy data for debug purposes.
	                if (debugIdx == i) {
	                    DetourCommon.vCopy(debug.optStart, ag.corridor.getPos());
	                    DetourCommon.vCopy(debug.optEnd, target);
	                }
	            } else {
	                // Copy data for debug purposes.
	                if (debugIdx == i) {
	                    DetourCommon.vSet(debug.optStart, 0, 0, 0);
	                    DetourCommon.vSet(debug.optEnd, 0, 0, 0);
	                }
	            }
	        }

	        // Trigger off-mesh connections (depends on corners).
	        for (let ag of agents) {

	            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
	                continue;
	            }
	            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE
	                || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
	                continue;
	            }

	            // Check
	            let triggerRadius = ag.params.radius * 2.25;
	            if (ag.overOffmeshConnection(triggerRadius)) {
	                // Prepare to off-mesh connection.
	                anim = ag.animation;

	                // Adjust the path over the off-mesh connection.
	                let refs = new long[2];
	                if (ag.corridor.moveOverOffmeshConnection(ag.corners[ag.corners.length - 1].getRef(), refs, anim.startPos, anim.endPos,
	                    this.m_navquery)) {
	                    DetourCommon.vCopy(anim.initPos, ag.npos);
	                    anim.polyRef = refs[1];
	                    anim.active = true;
	                    anim.t = 0.0;
	                    anim.tmax = (DetourCommon.vDist2D(anim.startPos, anim.endPos) / ag.params.maxSpeed) * 0.5;

	                    ag.state = CrowdAgent.DT_CROWDAGENT_STATE_OFFMESH;
	                    ag.corners = [];
	                    ag.neis = [];
	                    continue;
	                }
	            }
	        }

	        // Calculate steering.
	        for (let ag of agents) {

	            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
	                continue;
	            }
	            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE) {
	                continue;
	            }

	            let dvel = new Array(3);

	            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
	                DetourCommon.vCopy(dvel, ag.targetPos);
	                ag.desiredSpeed = DetourCommon.vLen(ag.targetPos);
	            } else {
	                // Calculate steering direction.
	                if ((ag.params.updateFlags & CrowdAgentParams.DT_CROWD_ANTICIPATE_TURNS) != 0) {
	                    dvel = ag.calcSmoothSteerDirection();
	                } else {
	                    dvel = ag.calcStraightSteerDirection();
	                }
	                // Calculate speed scale, which tells the agent to slowdown at the end of the path.
	                let slowDownRadius = ag.params.radius * 2; // TODO: make less hacky.
	                let speedScale = ag.getDistanceToGoal(slowDownRadius) / slowDownRadius;

	                ag.desiredSpeed = ag.params.maxSpeed;
	                dvel = DetourCommon.vScale(dvel, ag.desiredSpeed * speedScale);
	            }

	            // Separation
	            if ((ag.params.updateFlags & CrowdAgentParams.DT_CROWD_SEPARATION) != 0) {
	                separationDist = ag.params.collisionQueryRange;
	                invSeparationDist = 1.0 / separationDist;
	                separationWeight = ag.params.separationWeight;

	                w = 0;
	                let disp = new Array(3);

	                for (let j = 0; j < ag.neis.length; ++j) {
	                    nei = agents.get(ag.neis[j].idx);

	                    let diff = DetourCommon.vSub(ag.npos, nei.npos);
	                    diff[1] = 0;

	                    distSqr = DetourCommon.vLenSqr(diff);
	                    if (distSqr < 0.00001) {
	                        continue;
	                    }
	                    if (distSqr > DetourCommon.sqr(separationDist)) {
	                        continue;
	                    }
	                    dist = Math.sqrt(distSqr);
	                    weight = separationWeight * (1.0 - DetourCommon.sqr(dist * invSeparationDist));

	                    disp = DetourCommon.vMad(disp, diff, weight / dist);
	                    w += 1.0;
	                }

	                if (w > 0.0001) {
	                    // Adjust desired velocity.
	                    dvel = DetourCommon.vMad(dvel, disp, 1.0 / w);
	                    // Clamp desired velocity to desired speed.
	                    speedSqr = DetourCommon.vLenSqr(dvel);
	                    desiredSqr = DetourCommon.sqr(ag.desiredSpeed);
	                    if (speedSqr > desiredSqr) {
	                        dvel = DetourCommon.vScale(dvel, desiredSqr / speedSqr);
	                    }
	                }
	            }

	            // Set the desired velocity.
	            DetourCommon.vCopy(ag.dvel, dvel);
	        }

	        // Velocity planning.
	        for (let i = 0; i < agents.length; ++i) {
	            let ag = agents[i];

	            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
	                continue;
	            }

	            if ((ag.params.updateFlags & CrowdAgentParams.DT_CROWD_OBSTACLE_AVOIDANCE) != 0) {
	                this.m_obstacleQuery.reset();

	                // Add neighbours as obstacles.
	                for (let j = 0; j < ag.neis.length; ++j) {
	                    let nei = agents[ag.neis[j].idx];
	                    this.m_obstacleQuery.addCircle(nei.npos, nei.params.radius, nei.vel, nei.dvel);
	                }

	                // if (ag.neis.length > 0 && i == 0) {
	                //     console.log("Frame " + frame)
	                // }

	                // Append neighbour segments as obstacles.
	                for (let j = 0; j < ag.boundary.getSegmentCount(); ++j) {
	                    let s = ag.boundary.getSegment(j);
	                    //let s3 = Arrays.copyOfRange(s, 3, 6);
	                    //let s3 = Arrays.copyOfRange(s, 3, 6);
	                    let s3 = s.slice(3,6);
	                    if (DetourCommon.triArea2D3(ag.npos, s, s3) < 0.0) {
	                        continue;
	                    }
	                    this.m_obstacleQuery.addSegment(s, s3);
	                }

	                let vod = null;
	                if (debugIdx == i) {
	                    vod = debug.vod;
	                }
	                let ns = 0;

	                let params = this.m_obstacleQueryParams[ag.params.obstacleAvoidanceType];

	                {
	                    let nsnvel = this.m_obstacleQuery.sampleVelocityAdaptive(ag.npos, ag.params.radius, ag.desiredSpeed,
	                        ag.vel, ag.dvel, params, vod);
	                    ns = nsnvel[0];
	                    ag.nvel = nsnvel[1];
	                }
	                this.m_velocitySampleCount += ns;
	            } else {
	                // If not using velocity planning, new velocity is directly the desired velocity.
	                DetourCommon.vCopy(ag.nvel, ag.dvel);
	            }
	        }

	        // Integrate.
	        for (let i = 0; i < agents.length; ++i) {
	            let ag = agents[i];
	            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
	                continue;
	            }
	            ag.integrate(dt);
	        }

	        // Handle collisions.

	        for (let iter = 0; iter < 4; ++iter) {
	            for (let i = 0; i < agents.length; ++i) {
	                let ag = agents[i];
	                let idx0 = ag.getAgentIndex();
	                if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
	                    continue;
	                }

	                DetourCommon.vSet(ag.disp, 0, 0, 0);

	                let w = 0;

	                for (let j = 0; j < ag.neis.length; ++j) {
	                    let nei = agents[ag.neis[j].idx];
	                    let idx1 = nei.getAgentIndex();
	                    let diff = DetourCommon.vSub(ag.npos, nei.npos);
	                    diff[1] = 0;

	                    let dist = DetourCommon.vLenSqr(diff);
	                    if (dist > DetourCommon.sqr(ag.params.radius + nei.params.radius)) {
	                        continue;
	                    }
	                    dist = Math.sqrt(dist);
	                    let pen = (ag.params.radius + nei.params.radius) - dist;
	                    if (dist < 0.0001) {
	                        // Agents on top of each other, try to choose diverging separation directions.
	                        if (idx0 > idx1) {
	                            DetourCommon.vSet(diff, -ag.dvel[2], 0, ag.dvel[0]);
	                        } else {
	                            DetourCommon.vSet(diff, ag.dvel[2], 0, -ag.dvel[0]);
	                        }
	                        pen = 0.01;
	                    } else {
	                        pen = (1.0 / dist) * (pen * 0.5) * Crowd.COLLISION_RESOLVE_FACTOR;
	                    }

	                    ag.disp = DetourCommon.vMad(ag.disp, diff, pen);

	                    w += 1.0;
	                }

	                if (w > 0.0001) {
	                    let iw = 1.0 / w;
	                    ag.disp = DetourCommon.vScale(ag.disp, iw);
	                }
	            }

	            for (let i = 0; i < agents.length; ++i) {
	                let ag = agents[i];
	                if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
	                    continue;
	                }

	                ag.npos = DetourCommon.vAdd(ag.npos, ag.disp);
	            }
	        }

	        for (let i = 0; i < agents.length; ++i) {
	            // if (frame == 492)
	            //     console.log("Bad agent")
	            let ag = agents[i];
	            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
	                continue;
	            }

	            // Move along navmesh.
	            ag.corridor.movePosition(ag.npos, this.m_navquery, this.m_filters[ag.params.queryFilterType]);
	            // Get valid constrained position back.
	            DetourCommon.vCopy(ag.npos, ag.corridor.getPos());

	            // If not using path, truncate the corridor to just one poly.
	            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE
	                || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
	                ag.corridor.reset(ag.corridor.getFirstPoly(), ag.npos);
	                ag.partial = false;
	            }

	        }

	        // Update agents using off-mesh connection.
	        for (let i = 0; i < this.m_maxAgents; ++i) {
	            let anim = this.m_agents[i].animation;
	            if (!anim.active) {
	                continue;
	            }
	            ag = this.m_agents[i];

	            anim.t += dt;
	            if (anim.t > anim.tmax) {
	                // Reset animation
	                anim.active = false;
	                // Prepare agent for walking.
	                ag.state = CrowdAgent.DT_CROWDAGENT_STATE_WALKING;
	                continue;
	            }

	            // Update position
	            ta = anim.tmax * 0.15;
	            tb = anim.tmax;
	            if (anim.t < ta) {
	                u = tween(anim.t, 0.0, ta);
	                ag.npos = DetourCommon.vLerp3(anim.initPos, anim.startPos, u);
	            } else {
	                u = tween(anim.t, ta, tb);
	                ag.npos = DetourCommon.vLerp3(anim.startPos, anim.endPos, u);
	            }

	            // Update velocity.
	            DetourCommon.vSet(ag.vel, 0, 0, 0);
	            DetourCommon.vSet(ag.dvel, 0, 0, 0);
	        }
	    }

	    getQueryExtents() {
	        return this.m_ext;
	    }

	    getFilter(i) {
	        return i >= 0 && i < Crowd.DT_CROWD_MAX_QUERY_FILTER_TYPE ? this.m_filters[i] : null;
	    }

	}

	exports.Crowd = Crowd;
	exports.CrowdAgentParams = CrowdAgentParams;
	exports.NavMesh = NavMesh;
	exports.NavMeshQuery = NavMeshQuery;
	exports.ObstacleAvoidanceParams = ObstacleAvoidanceParams;
	exports.RecastTestMeshBuilder = RecastTestMeshBuilder;

	Object.defineProperty(exports, '__esModule', { value: true });

})));
