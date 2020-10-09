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

import ObstacleAvoidanceParams from "./ObstacleAvoidanceParams.js"
import ObstacleCircle from "./ObstacleCircle.js"
import ObstacleSegment from "./ObstacleSegment.js"
import DetourCommon from "../DetourCommon.js"
import SweepCircleCircleResult from "./SweepCircleCircleResult.js"


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
		let d = DetourCommon.vDot2D(u, v);
		if (Math.abs(d) < 1e-6)
		return [false, 0];
			d = 1.0 / d;
		let t = DetourCommon.vDot2D(v, w) * d;
		if (t < 0 || t > 1)
			return [false, 0];
		let s = DetourCommon.vDot2D(u, w) * d;
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

export default ObstacleAvoidanceQuery;