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

import IntersectResult from "./IntersectResult.js"
import VectorPtr from "./VectorPtr.js"

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
		// TODO: Replace pnPoly with triArea2D tests?
		let i, j;
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
		// TODO: Replace pnPoly with triArea2D tests?
		let i, j;
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

		let p0v = new VectorPtr(p0);
		for (let i = 0, j = nverts - 1; i < nverts; j = i++) {
			let vpj = new VectorPtr(verts, j * 3);
			let edge = DetourCommon.vSub(new VectorPtr(verts, i * 3), vpj);
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

export default DetourCommon;
