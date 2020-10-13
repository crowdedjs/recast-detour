"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _ObstacleAvoidanceParams = _interopRequireDefault(require("./ObstacleAvoidanceParams.js"));

var _ObstacleCircle = _interopRequireDefault(require("./ObstacleCircle.js"));

var _ObstacleSegment = _interopRequireDefault(require("./ObstacleSegment.js"));

var _DetourCommon = _interopRequireDefault(require("../DetourCommon.js"));

var _SweepCircleCircleResult = _interopRequireDefault(require("./SweepCircleCircleResult.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var ObstacleAvoidanceQuery = /*#__PURE__*/function () {
  ///< Max numver of adaptive divs.
  ///< Max number of adaptive rings.
  function ObstacleAvoidanceQuery(maxCircles, maxSegments) {
    _classCallCheck(this, ObstacleAvoidanceQuery);

    _defineProperty(this, "m_params", new _ObstacleAvoidanceParams["default"]());

    _defineProperty(this, "m_invHorizTime", 0);

    _defineProperty(this, "m_vmax", 0);

    _defineProperty(this, "m_invVmax", 0);

    _defineProperty(this, "m_maxCircles", 0);

    _defineProperty(this, "m_circles", []);

    _defineProperty(this, "m_ncircles", 0);

    _defineProperty(this, "m_maxSegments", 0);

    _defineProperty(this, "m_segments", []);

    _defineProperty(this, "m_nsegments", 0);

    this.m_maxCircles = maxCircles;
    this.m_ncircles = 0;
    this.m_circles = new Array(this.m_maxCircles);

    for (var i = 0; i < this.m_maxCircles; i++) {
      this.m_circles[i] = new _ObstacleCircle["default"]();
    }

    this.m_maxSegments = maxSegments;
    this.m_nsegments = 0;
    this.m_segments = new Array(this.m_maxSegments);

    for (var _i = 0; _i < this.m_maxSegments; _i++) {
      this.m_segments[_i] = new _ObstacleSegment["default"]();
    }
  }

  _createClass(ObstacleAvoidanceQuery, [{
    key: "reset",
    value: function reset() {
      this.m_ncircles = 0;
      this.m_nsegments = 0;
    }
  }, {
    key: "addCircle",
    value: function addCircle(pos, rad, vel, dvel) {
      if (this.m_ncircles >= this.m_maxCircles) return;
      var cir = this.m_circles[this.m_ncircles++];

      _DetourCommon["default"].vCopy(cir.p, pos);

      cir.rad = rad;

      _DetourCommon["default"].vCopy(cir.vel, vel);

      _DetourCommon["default"].vCopy(cir.dvel, dvel);
    }
  }, {
    key: "addSegment",
    value: function addSegment(p, q) {
      if (this.m_nsegments >= this.m_maxSegments) return;
      var seg = this.m_segments[this.m_nsegments++];

      _DetourCommon["default"].vCopy(seg.p, p);

      _DetourCommon["default"].vCopy(seg.q, q);
    }
  }, {
    key: "getObstacleCircleCount",
    value: function getObstacleCircleCount() {
      return this.m_ncircles;
    }
  }, {
    key: "getObstacleCircle",
    value: function getObstacleCircle(i) {
      return this.m_circles[i];
    }
  }, {
    key: "getObstacleSegmentCount",
    value: function getObstacleSegmentCount() {
      return this["this"].m_nsegments;
    }
  }, {
    key: "getObstacleSegment",
    value: function getObstacleSegment(i) {
      return this.m_segments[i];
    }
  }, {
    key: "prepare",
    value: function prepare(pos, dvel) {
      // Prepare obstacles
      for (var i = 0; i < this.m_ncircles; ++i) {
        var cir = this.m_circles[i]; // Side

        var pa = pos;
        var pb = cir.p;
        var orig = [0, 0, 0];
        var dv = new Array(3);

        _DetourCommon["default"].vCopy(cir.dp, _DetourCommon["default"].vSub(pb, pa));

        _DetourCommon["default"].vNormalize(cir.dp);

        dv = _DetourCommon["default"].vSub(cir.dvel, dvel);

        var a = _DetourCommon["default"].triArea2D3(orig, cir.dp, dv);

        if (a < 0.01) {
          cir.np[0] = -cir.dp[2];
          cir.np[2] = cir.dp[0];
        } else {
          cir.np[0] = cir.dp[2];
          cir.np[2] = -cir.dp[0];
        }
      }

      for (var _i2 = 0; _i2 < this.m_nsegments; ++_i2) {
        var seg = this.m_segments[_i2]; // Precalc if the agent is really close to the segment.

        var r = 0.01;

        var dt = _DetourCommon["default"].distancePtSegSqr2D3(pos, seg.p, seg.q);

        seg.touch = dt[0] < _DetourCommon["default"].sqr(r);
      }
    }
  }, {
    key: "sweepCircleCircle",
    value: function sweepCircleCircle(c0, r0, v, c1, r1) {
      var EPS = 0.0001;

      var s = _DetourCommon["default"].vSub(c1, c0);

      var r = r0 + r1;
      var c = _DetourCommon["default"].vDot2D(s, s) - r * r;

      var a = _DetourCommon["default"].vDot2D(v, v);

      if (a < EPS) return new _SweepCircleCircleResult["default"](false, 0, 0); // not moving
      // Overlap, calc time to exit.

      var b = _DetourCommon["default"].vDot2D(v, s);

      var d = b * b - a * c;
      if (d < 0.0) return new _SweepCircleCircleResult["default"](false, 0, 0); // no intersection.

      a = 1.0 / a;
      var rd = Math.sqrt(d);
      return new _SweepCircleCircleResult["default"](true, (b - rd) * a, (b + rd) * a);
    }
  }, {
    key: "isectRaySeg",
    value: function isectRaySeg(ap, u, bp, bq) {
      var v = _DetourCommon["default"].vSub(bq, bp);

      var w = _DetourCommon["default"].vSub(ap, bp);

      var d = _DetourCommon["default"].vPerp2D(u, v);

      if (Math.abs(d) < 1e-6) return [false, 0];
      d = 1.0 / d;
      var t = _DetourCommon["default"].vPerp2D(v, w) * d;
      if (t < 0 || t > 1) return [false, 0];
      var s = _DetourCommon["default"].vPerp2D(u, w) * d;
      if (s < 0 || s > 1) return [false, 0];
      return [true, t];
    }
    /** Calculate the collision penalty for a given velocity vector
     * 
     * @param vcand sampled velocity
     * @param dvel desired velocity
     * @param minPenalty threshold penalty for early out
     */

  }, {
    key: "processSample",
    value: function processSample(vcand, cs, pos, rad, vel, dvel, minPenalty, debug) {
      // penalty for straying away from the desired and current velocities
      var vpen = this.m_params.weightDesVel * (_DetourCommon["default"].vDist2D(vcand, dvel) * this.m_invVmax);
      var vcpen = this.m_params.weightCurVel * (_DetourCommon["default"].vDist2D(vcand, vel) * this.m_invVmax); // find the threshold hit time to bail out based on the early out penalty
      // (see how the penalty is calculated below to understnad)

      var minPen = minPenalty - vpen - vcpen;
      var tThresold = (this.m_params.weightToi / minPen - 0.1) * this.m_params.horizTime;
      if (tThresold - this.m_params.horizTime > -Number.MIN_VALUE) return minPenalty; // already too much
      // Find min time of impact and exit amongst all obstacles.

      var tmin = this.m_params.horizTime;
      var side = 0;
      var nside = 0;

      for (var i = 0; i < this.m_ncircles; ++i) {
        var cir = this.m_circles[i]; // RVO

        var vab = _DetourCommon["default"].vScale(vcand, 2);

        vab = _DetourCommon["default"].vSub(vab, vel);
        vab = _DetourCommon["default"].vSub(vab, cir.vel); // Side

        side += _DetourCommon["default"].clamp(Math.min(_DetourCommon["default"].vDot2D(cir.dp, vab) * 0.5 + 0.5, _DetourCommon["default"].vDot2D(cir.np, vab) * 2), 0.0, 1.0);
        nside++;
        var sres = this.sweepCircleCircle(pos, rad, vab, cir.p, cir.rad);
        if (!sres.intersection) continue;
        var htmin = sres.htmin,
            htmax = sres.htmax; // Handle overlapping obstacles.

        if (htmin < 0.0 && htmax > 0.0) {
          // Amore when overlapped.
          htmin = -htmin * 0.5;
        }

        if (htmin >= 0.0) {
          // The closest obstacle is somewhere ahead of us, keep track of nearest obstacle.
          if (htmin < tmin) {
            tmin = htmin;
            if (tmin < tThresold) return minPenalty;
          }
        }
      }

      for (var _i3 = 0; _i3 < this.m_nsegments; ++_i3) {
        var seg = this.m_segments[_i3];
        var _htmin = 0;

        if (seg.touch) {
          // Special case when the agent is very close to the segment.
          var sdir = _DetourCommon["default"].vSub(seg.q, seg.p);

          var snorm = new Array(3);
          snorm[0] = -sdir[2];
          snorm[2] = sdir[0]; // If the velocity is pointing towards the segment, no collision.

          if (_DetourCommon["default"].vDot2D(snorm, vcand) < 0.0) continue; // Else immediate collision.

          _htmin = 0.0;
        } else {
          var ires = this.isectRaySeg(pos, vcand, seg.p, seg.q);
          if (!ires[0]) continue;
          _htmin = ires[1];
        } // Aless when facing walls.


        _htmin *= 2.0; // The closest obstacle is somewhere ahead of us, keep track of nearest obstacle.

        if (_htmin < tmin) {
          tmin = _htmin;
          if (tmin < tThresold) return minPenalty;
        }
      } // Normalize side bias, to prevent it dominating too much.


      if (nside != 0) side /= nside;
      var spen = this.m_params.weightSide * side;
      var tpen = this.m_params.weightToi * (1.0 / (0.1 + tmin * this.m_invHorizTime));
      var penalty = vpen + vcpen + spen + tpen; // Store different penalties for debug viewing

      if (debug != null) debug.addSample(vcand, cs, penalty, vpen, vcpen, spen, tpen);
      return penalty;
    }
  }, {
    key: "sampleVelocityGrid",
    value: function sampleVelocityGrid(pos, rad, vmax, vel, dvel, params, debug) {
      this.prepare(pos, dvel);
      this.m_params = params;
      this.m_invHorizTime = 1.0 / this.m_params.horizTime;
      this.m_vmax = vmax;
      this.m_invVmax = vmax > 0 ? 1.0 / vmax : Number.MAX_VALUE;
      var nvel = new Array(3);

      _DetourCommon["default"].vSet(nvel, 0, 0, 0);

      if (debug != null) debug.reset();
      cvx = dvel[0] * this.m_params.velBias;
      cvz = dvel[2] * this.m_params.velBias;
      cs = vmax * 2 * (1 - this.m_params.velBias) / (this.m_params.gridSize - 1);
      half = (this.m_params.gridSize - 1) * cs * 0.5;
      minPenalty = Number.MAX_VALUE;
      var ns = 0;

      for (var y = 0; y < this.m_params.gridSize; ++y) {
        for (var x = 0; x < this.m_params.gridSize; ++x) {
          var vcand = new Array(3);

          _DetourCommon["default"].vSet(vcand, cvx + x * cs - half, 0, cvz + y * cs - half);

          if (_DetourCommon["default"].sqr(vcand[0]) + _DetourCommon["default"].sqr(vcand[2]) > _DetourCommon["default"].sqr(vmax + cs / 2)) continue;
          penalty = this.processSample(vcand, cs, pos, rad, vel, dvel, minPenalty, debug);
          ns++;

          if (penalty < minPenalty) {
            minPenalty = penalty;

            _DetourCommon["default"].vCopy(nvel, vcand);
          }
        }
      }

      return [ns, nvel];
    } // vector normalization that ignores the y-component.

  }, {
    key: "dtNormalize2D",
    value: function dtNormalize2D(v) {
      var d = Math.sqrt(v[0] * v[0] + v[2] * v[2]);
      if (d == 0) return;
      d = 1.0 / d;
      v[0] *= d;
      v[2] *= d;
    } // vector normalization that ignores the y-component.

  }, {
    key: "dtRotate2D",
    value: function dtRotate2D(v, ang) {
      var dest = new Array(3);
      var c = Math.cos(ang);
      var s = Math.sin(ang);
      dest[0] = v[0] * c - v[2] * s;
      dest[2] = v[0] * s + v[2] * c;
      dest[1] = v[1];
      return dest;
    }
  }, {
    key: "sampleVelocityAdaptive",
    value: function sampleVelocityAdaptive(pos, rad, vmax, vel, dvel, params, debug) {
      this.prepare(pos, dvel);
      this.m_params = params;
      this.m_invHorizTime = 1.0 / this.m_params.horizTime;
      this.m_vmax = vmax;
      this.m_invVmax = vmax > 0 ? 1.0 / vmax : Number.MAX_VALUE;
      var nvel = new Array(3);

      _DetourCommon["default"].vSet(nvel, 0, 0, 0);

      if (debug != null) debug.reset(); // Build sampling pattern aligned to desired velocity.

      var pat = new Array((ObstacleAvoidanceQuery.DT_MAX_PATTERN_DIVS * ObstacleAvoidanceQuery.DT_MAX_PATTERN_RINGS + 1) * 2);
      var npat = 0;
      var ndivs = this.m_params.adaptiveDivs;
      var nrings = this.m_params.adaptiveRings;
      var depth = this.m_params.adaptiveDepth;

      var nd = _DetourCommon["default"].clamp(ndivs, 1, ObstacleAvoidanceQuery.DT_MAX_PATTERN_DIVS);

      var nr = _DetourCommon["default"].clamp(nrings, 1, ObstacleAvoidanceQuery.DT_MAX_PATTERN_RINGS);

      var da = 1.0 / nd * ObstacleAvoidanceQuery.DT_PI * 2;
      var ca = Math.cos(da);
      var sa = Math.sin(da); // desired direction

      var ddir = new Array(6);

      _DetourCommon["default"].vCopy(ddir, dvel);

      this.dtNormalize2D(ddir);
      var rotated = this.dtRotate2D(ddir, da * 0.5); // rotated by da/2

      ddir[3] = rotated[0];
      ddir[4] = rotated[1];
      ddir[5] = rotated[2]; // Always add sample at zero

      pat[npat * 2 + 0] = 0;
      pat[npat * 2 + 1] = 0;
      npat++;

      for (var j = 0; j < nr; ++j) {
        var r = (nr - j) / nr;
        pat[npat * 2 + 0] = ddir[j % 2 * 3] * r;
        pat[npat * 2 + 1] = ddir[j % 2 * 3 + 2] * r;
        var last1 = npat * 2;
        var last2 = last1;
        npat++;

        for (var i = 1; i < nd - 1; i += 2) {
          // get next poPoly on the "right" (rotate CW)
          pat[npat * 2 + 0] = pat[last1] * ca + pat[last1 + 1] * sa;
          pat[npat * 2 + 1] = -pat[last1] * sa + pat[last1 + 1] * ca; // get next poPoly on the "left" (rotate CCW)

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
      } // Start sampling.


      var cr = vmax * (1.0 - this.m_params.velBias);
      var res = new Array(3);

      _DetourCommon["default"].vSet(res, dvel[0] * this.m_params.velBias, 0, dvel[2] * this.m_params.velBias);

      var ns = 0;

      for (var k = 0; k < depth; ++k) {
        var _minPenalty = Number.MAX_VALUE;
        var bvel = new Array(3);

        _DetourCommon["default"].vSet(bvel, 0, 0, 0);

        for (var _i4 = 0; _i4 < npat; ++_i4) {
          var vcand = new Array(3);

          _DetourCommon["default"].vSet(vcand, res[0] + pat[_i4 * 2 + 0] * cr, 0, res[2] + pat[_i4 * 2 + 1] * cr);

          if (_DetourCommon["default"].sqr(vcand[0]) + _DetourCommon["default"].sqr(vcand[2]) > _DetourCommon["default"].sqr(vmax + 0.001)) continue;

          var _penalty = this.processSample(vcand, cr / 10, pos, rad, vel, dvel, _minPenalty, debug);

          ns++;

          if (_penalty < _minPenalty) {
            _minPenalty = _penalty;

            _DetourCommon["default"].vCopy(bvel, vcand);
          }
        }

        _DetourCommon["default"].vCopy(res, bvel);

        cr *= 0.5;
      }

      _DetourCommon["default"].vCopy(nvel, res);

      return [ns, nvel];
    }
  }]);

  return ObstacleAvoidanceQuery;
}();

_defineProperty(ObstacleAvoidanceQuery, "DT_MAX_PATTERN_DIVS", 32);

_defineProperty(ObstacleAvoidanceQuery, "DT_MAX_PATTERN_RINGS", 4);

_defineProperty(ObstacleAvoidanceQuery, "DT_PI", 3.14159265);

var _default = ObstacleAvoidanceQuery;
exports["default"] = _default;