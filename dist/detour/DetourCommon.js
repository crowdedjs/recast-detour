"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _IntersectResult = _interopRequireDefault(require("./IntersectResult.js"));

var _VectorPtr = _interopRequireDefault(require("./VectorPtr.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var DetourCommon = /*#__PURE__*/function () {
  function DetourCommon() {
    _classCallCheck(this, DetourCommon);
  }

  _createClass(DetourCommon, null, [{
    key: "vMad",
    /// Performs a scaled vector addition. (@p v1 + (@p v2 * @p s))
    /// @param[out] dest The result vector. [(x, y, z)]
    /// @param[_in] v1 The base vector. [(x, y, z)]
    /// @param[_in] v2 The vector to scale and add to @p v1. [(x, y, z)]
    /// @param[_in] s The amount to scale @p v2 by before adding to @p v1.
    value: function vMad(v1, v2, s) {
      var dest = new Array(3);
      dest[0] = v1[0] + v2[0] * s;
      dest[1] = v1[1] + v2[1] * s;
      dest[2] = v1[2] + v2[2] * s;
      return dest;
    } /// Performs a linear interpolation between two vectors. (@p v1 toward @p
    /// v2)
    /// @param[out] dest The result vector. [(x, y, x)]
    /// @param[_in] v1 The starting vector.
    /// @param[_in] v2 The destination vector.
    /// @param[_in] t The interpolation factor. [Limits: 0 <= value <= 1.0]

  }, {
    key: "vLerp4",
    value: function vLerp4(verts, v1, v2, t) {
      var dest = new Array(3);
      dest[0] = verts[v1 + 0] + (verts[v2 + 0] - verts[v1 + 0]) * t;
      dest[1] = verts[v1 + 1] + (verts[v2 + 1] - verts[v1 + 1]) * t;
      dest[2] = verts[v1 + 2] + (verts[v2 + 2] - verts[v1 + 2]) * t;
      return dest;
    }
  }, {
    key: "vLerp3",
    value: function vLerp3(v1, v2, t) {
      var dest = new Array(3);
      dest[0] = v1[0] + (v2[0] - v1[0]) * t;
      dest[1] = v1[1] + (v2[1] - v1[1]) * t;
      dest[2] = v1[2] + (v2[2] - v1[2]) * t;
      return dest;
    }
  }, {
    key: "vSub",
    value: function vSub(v1, v2) {
      var dest = new Array(3);
      dest[0] = v1[0] - v2[0];
      dest[1] = v1[1] - v2[1];
      dest[2] = v1[2] - v2[2];
      return dest;
    } // static vSub(v1, v2) {
    // 	let dest = new Array(3);
    // 	dest[0] = v1[0] - v2[0];
    // 	dest[1] = v1[1] - v2[1];
    // 	dest[2] = v1[2] - v2[2];
    // 	return dest;
    // }

  }, {
    key: "vAdd",
    value: function vAdd(v1, v2) {
      var dest = new Array(3);
      dest[0] = v1[0] + v2[0];
      dest[1] = v1[1] + v2[1];
      dest[2] = v1[2] + v2[2];
      return dest;
    }
  }, {
    key: "vCopy_return",
    value: function vCopy_return(i) {
      var out = new Array(3);
      out[0] = i[0];
      out[1] = i[1];
      out[2] = i[2];
      return out;
    }
  }, {
    key: "vSet",
    value: function vSet(out, a, b, c) {
      out[0] = a;
      out[1] = b;
      out[2] = c;
    } // static vCopy(out, i) {
    // 	out[0] = i[0];
    // 	out[1] = i[1];
    // 	out[2] = i[2];
    // }
    //See vCopy_return for the 1 parameter version

  }, {
    key: "vCopy",
    value: function vCopy(out, _in) {
      var i = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : 0;
      out[0] = _in[i];
      out[1] = _in[i + 1];
      out[2] = _in[i + 2];
    }
  }, {
    key: "vMin",
    value: function vMin(out, _in, i) {
      out[0] = Math.min(out[0], _in[i]);
      out[1] = Math.min(out[1], _in[i + 1]);
      out[2] = Math.min(out[2], _in[i + 2]);
    }
  }, {
    key: "vMax",
    value: function vMax(out, _in, i) {
      out[0] = Math.max(out[0], _in[i]);
      out[1] = Math.max(out[1], _in[i + 1]);
      out[2] = Math.max(out[2], _in[i + 2]);
    } /// Returns the distance between two points.
    /// @param[_in] v1 A point. [(x, y, z)]
    /// @param[_in] v2 A point. [(x, y, z)]
    /// @return The distance between the two points.

  }, {
    key: "vDist2",
    value: function vDist2(v1, v2) {
      var dx = v2[0] - v1[0];
      var dy = v2[1] - v1[1];
      var dz = v2[2] - v1[2];
      return Math.sqrt(dx * dx + dy * dy + dz * dz);
    } /// Returns the distance between two points.
    /// @param[_in] v1 A point. [(x, y, z)]
    /// @param[_in] v2 A point. [(x, y, z)]
    /// @return The distance between the two points.

  }, {
    key: "vDistSqr",
    value: function vDistSqr(v1, v2) {
      var dx = v2[0] - v1[0];
      var dy = v2[1] - v1[1];
      var dz = v2[2] - v1[2];
      return dx * dx + dy * dy + dz * dz;
    }
  }, {
    key: "sqr",
    value: function sqr(a) {
      return a * a;
    } /// Derives the square of the scalar length of the vector. (len * len)
    /// @param[_in] v The vector. [(x, y, z)]
    /// @return The square of the scalar length of the vector.

  }, {
    key: "vLenSqr",
    value: function vLenSqr(v) {
      return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    }
  }, {
    key: "vLen",
    value: function vLen(v) {
      return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
  }, {
    key: "vDist3",
    value: function vDist3(v1, verts, i) {
      var dx = verts[i] - v1[0];
      var dy = verts[i + 1] - v1[1];
      var dz = verts[i + 2] - v1[2];
      return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }
  }, {
    key: "clamp",
    value: function clamp(v, min, max) {
      return Math.max(Math.min(v, max), min);
    }
  }, {
    key: "clamp",
    value: function clamp(v, min, max) {
      return Math.max(Math.min(v, max), min);
    } /// Derives the distance between the specified points on the xz-plane.
    /// @param[_in] v1 A point. [(x, y, z)]
    /// @param[_in] v2 A point. [(x, y, z)]
    /// @return The distance between the poPoly on the xz-plane.
    ///
    /// The vectors are projected onto the xz-plane, so the y-values are
    /// ignored.

  }, {
    key: "vDist2D",
    value: function vDist2D(v1, v2) {
      var dx = v2[0] - v1[0];
      var dz = v2[2] - v1[2];
      return Math.sqrt(dx * dx + dz * dz);
    }
  }, {
    key: "vDist2DSqr",
    value: function vDist2DSqr(v1, v2) {
      var dx = v2[0] - v1[0];
      var dz = v2[2] - v1[2];
      return dx * dx + dz * dz;
    } /// Normalizes the vector.
    /// @param[_in,out] v The vector to normalize. [(x, y, z)]

  }, {
    key: "vNormalize",
    value: function vNormalize(v) {
      var d = 1.0 / Math.sqrt(DetourCommon.sqr(v[0]) + DetourCommon.sqr(v[1]) + DetourCommon.sqr(v[2]));

      if (d != 0) {
        v[0] *= d;
        v[1] *= d;
        v[2] *= d;
      }
    }
  }, {
    key: "vEqual",
    /// Performs a 'sloppy' colocation check of the specified points.
    /// @param[_in] p0 A point. [(x, y, z)]
    /// @param[_in] p1 A point. [(x, y, z)]
    /// @return True if the points are considered to be at the same location.
    ///
    /// Basically, this function will return true if the specified points are
    /// close enough to eachother to be considered colocated.
    value: function vEqual(p0, p1) {
      var d = DetourCommon.vDistSqr(p0, p1);
      return d < DetourCommon.thr;
    } /// Derives the dot product of two vectors on the xz-plane. (@p u . @p v)
    /// @param[_in] u A vector [(x, y, z)]
    /// @param[_in] v A vector [(x, y, z)]
    /// @return The dot product on the xz-plane.
    ///
    /// The vectors are projected onto the xz-plane, so the y-values are
    /// ignored.
    // static DetourCommon.vDot2D(u, v) {
    // 	return u[0] * v[0] + u[2] * v[2];
    // }

  }, {
    key: "vDot2D",
    value: function vDot2D(u, v) {
      var vi = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : 0;
      return u[0] * v[vi] + u[2] * v[vi + 2];
    } /// Derives the xz-plane 2D perp product of the two vectors. (uz*vx - ux*vz)
    /// @param[_in] u The LHV vector [(x, y, z)]
    /// @param[_in] v The RHV vector [(x, y, z)]
    /// @return The dot product on the xz-plane.
    ///
    /// The vectors are projected onto the xz-plane, so the y-values are
    /// ignored.

  }, {
    key: "vPerp2D",
    value: function vPerp2D(u, v) {
      return u[2] * v[0] - u[0] * v[2];
    } /// @}
    /// @name Computational geometry helper functions.
    /// @{
    /// Derives the signed xz-plane area of the triangle ABC, or the
    /// relationship of line AB to poPoly C.
    /// @param[_in] a Vertex A. [(x, y, z)]
    /// @param[_in] b Vertex B. [(x, y, z)]
    /// @param[_in] c Vertex C. [(x, y, z)]
    /// @return The signed xz-plane area of the triangle.

  }, {
    key: "triArea2D4",
    value: function triArea2D4(verts, a, b, c) {
      abx = verts[b] - verts[a];
      abz = verts[b + 2] - verts[a + 2];
      acx = verts[c] - verts[a];
      acz = verts[c + 2] - verts[a + 2];
      return acx * abz - abx * acz;
    }
  }, {
    key: "triArea2D3",
    value: function triArea2D3(a, b, c) {
      var abx = b[0] - a[0];
      var abz = b[2] - a[2];
      var acx = c[0] - a[0];
      var acz = c[2] - a[2];
      return acx * abz - abx * acz;
    } /// Determines if two axis-aligned bounding boxes overlap.
    /// @param[_in] amin Minimum bounds of box A. [(x, y, z)]
    /// @param[_in] amax Maximum bounds of box A. [(x, y, z)]
    /// @param[_in] bmin Minimum bounds of box B. [(x, y, z)]
    /// @param[_in] bmax Maximum bounds of box B. [(x, y, z)]
    /// @return True if the two AABB's overlap.
    /// @see dtOverlapBounds

  }, {
    key: "overlapQuantBounds",
    value: function overlapQuantBounds(amin, amax, bmin, bmax) {
      var overlap = true;
      overlap = amin[0] > bmax[0] || amax[0] < bmin[0] ? false : overlap;
      overlap = amin[1] > bmax[1] || amax[1] < bmin[1] ? false : overlap;
      overlap = amin[2] > bmax[2] || amax[2] < bmin[2] ? false : overlap;
      return overlap;
    } /// Determines if two axis-aligned bounding boxes overlap.
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

  }, {
    key: "distancePtSegSqr2D3",
    value: function distancePtSegSqr2D3(pt, p, q) {
      var pqx = q[0] - p[0];
      var pqz = q[2] - p[2];
      var dx = pt[0] - p[0];
      var dz = pt[2] - p[2];
      var d = pqx * pqx + pqz * pqz;
      var t = pqx * dx + pqz * dz;
      if (d > 0) t /= d;
      if (t < 0) t = 0;else if (t > 1) t = 1;
      dx = p[0] + t * pqx - pt[0];
      dz = p[2] + t * pqz - pt[2];
      return [dx * dx + dz * dz, t];
    }
  }, {
    key: "closestHeightPointTriangle",
    value: function closestHeightPointTriangle(p, a, b, c) {
      var v0 = DetourCommon.vSub(c, a);
      var v1 = DetourCommon.vSub(b, a);
      var v2 = DetourCommon.vSub(p, a);
      var dot00 = DetourCommon.vDot2D(v0, v0);
      var dot01 = DetourCommon.vDot2D(v0, v1);
      var dot02 = DetourCommon.vDot2D(v0, v2);
      var dot11 = DetourCommon.vDot2D(v1, v1);
      var dot12 = DetourCommon.vDot2D(v1, v2); // Compute barycentric coordinates

      var invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
      var u = (dot11 * dot02 - dot01 * dot12) * invDenom;
      var v = (dot00 * dot12 - dot01 * dot02) * invDenom; // The (sloppy) epsilon is needed to allow to get height of points which
      // are interpolated aPoly the edges of the triangles.
      // If point lies inside the triangle, return interpolated ycoord.

      if (u >= -DetourCommon.EPS && v >= -DetourCommon.EPS && u + v <= 1 + DetourCommon.EPS) {
        var h = a[1] + v0[1] * u + v1[1] * v;
        return [true, h];
      }

      return [false, null];
    } /// @par
    ///
    /// All points are projected onto the xz-plane, so the y-values are ignored.

  }, {
    key: "pointInPolygon",
    value: function pointInPolygon(pt, verts, nverts) {
      // TODO: Replace pnPoly with triArea2D tests?
      var i, j;
      var c = false;

      for (var _i = 0, _j = nverts - 1; _i < nverts; _j = _i++) {
        var vi = _i * 3;
        var vj = _j * 3;
        if (verts[vi + 2] > pt[2] != verts[vj + 2] > pt[2] && pt[0] < (verts[vj + 0] - verts[vi + 0]) * (pt[2] - verts[vi + 2]) / (verts[vj + 2] - verts[vi + 2]) + verts[vi + 0]) c = !c;
      }

      return c;
    }
  }, {
    key: "distancePtPolyEdgesSqr",
    value: function distancePtPolyEdgesSqr(pt, verts, nverts, ed, et) {
      // TODO: Replace pnPoly with triArea2D tests?
      var i, j;
      var c = false;

      for (var _i2 = 0, _j2 = nverts - 1; _i2 < nverts; _j2 = _i2++) {
        var vi = _i2 * 3;
        var vj = _j2 * 3;
        if (verts[vi + 2] > pt[2] != verts[vj + 2] > pt[2] && pt[0] < (verts[vj + 0] - verts[vi + 0]) * (pt[2] - verts[vi + 2]) / (verts[vj + 2] - verts[vi + 2]) + verts[vi + 0]) c = !c;
        var edet = DetourCommon.distancePtSegSqr2D4(pt, verts, vj, vi);
        ed[_j2] = edet[0];
        et[_j2] = edet[1];
      }

      return c;
    }
  }, {
    key: "projectPoly",
    value: function projectPoly(axis, poly, npoly) {
      var rmin;
      var rmax;
      rmin = rmax = DetourCommon.vDot2D(axis, poly, 0);

      for (var i = 1; i < npoly; ++i) {
        var _d = DetourCommon.vDot2D(axis, poly, i * 3);

        rmin = Math.min(rmin, _d);
        rmax = Math.max(rmax, _d);
      }

      return [rmin, rmax];
    }
  }, {
    key: "overlapRange",
    value: function overlapRange(amin, amax, bmin, bmax, eps) {
      return amin + eps > bmax || amax - eps < bmin ? false : true;
    }
  }, {
    key: "overlapPolyPoly2D",
    /// @par
    ///
    /// All vertices are projected onto the xz-plane, so the y-values are ignored.
    value: function overlapPolyPoly2D(polya, npolya, polyb, npolyb) {
      for (var i = 0, j = npolya - 1; i < npolya; j = i++) {
        var va = j * 3;
        var vb = i * 3;
        var n = [polya[vb + 2] - polya[va + 2], 0, -(polya[vb + 0] - polya[va + 0])];
        var aminmax = DetourCommon.projectPoly(n, polya, npolya);
        var bminmax = DetourCommon.projectPoly(n, polyb, npolyb);

        if (!DetourCommon.overlapRange(aminmax[0], aminmax[1], bminmax[0], bminmax[1], DetourCommon.eps)) {
          // Found separating axis
          return false;
        }
      }

      for (var _i3 = 0, _j3 = npolyb - 1; _i3 < npolyb; _j3 = _i3++) {
        var _va = _j3 * 3;

        var _vb = _i3 * 3;

        var _n = [polyb[_vb + 2] - polyb[_va + 2], 0, -(polyb[_vb + 0] - polyb[_va + 0])];

        var _aminmax = DetourCommon.projectPoly(_n, polya, npolya);

        var _bminmax = DetourCommon.projectPoly(_n, polyb, npolyb);

        if (!DetourCommon.overlapRange(_aminmax[0], _aminmax[1], _bminmax[0], _bminmax[1], DetourCommon.eps)) {
          // Found separating axis
          return false;
        }
      }

      return true;
    } // Returns a random poPoly _in a convex polygon.
    // Adapted from Graphics Gems article.

  }, {
    key: "randomPointInConvexPoly",
    value: function randomPointInConvexPoly(pts, npts, areas, s, t) {
      // Calc triangle araes
      areasum = 0.0;

      for (var i = 2; i < npts; i++) {
        areas[i] = DetourCommon.triArea2D4(pts, 0, (i - 1) * 3, i * 3);
        areasum += Math.max(0.001, areas[i]);
      } // Find sub triangle weighted by area.


      thr = s * areasum;
      acc = 0.0;
      u = 1.0;
      var tri = npts - 1;

      for (var _i4 = 2; _i4 < npts; _i4++) {
        dacc = areas[_i4];

        if (thr >= acc && thr < acc + dacc) {
          u = (thr - acc) / dacc;
          tri = _i4;
          break;
        }

        acc += dacc;
      }

      v = Math.sqrt(t);
      a = 1 - v;
      b = (1 - u) * v;
      c = u * v;
      var pa = 0;
      var pb = (tri - 1) * 3;
      var pc = tri * 3;
      return [a * pts[pa] + b * pts[pb] + c * pts[pc], a * pts[pa + 1] + b * pts[pb + 1] + c * pts[pc + 1], a * pts[pa + 2] + b * pts[pb + 2] + c * pts[pc + 2]];
    }
  }, {
    key: "nextPow2",
    value: function nextPow2(v) {
      v--;
      v |= v >> 1;
      v |= v >> 2;
      v |= v >> 4;
      v |= v >> 8;
      v |= v >> 16;
      v++;
      return v;
    }
  }, {
    key: "ilog2",
    value: function ilog2(v) {
      var r;
      var shift;
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
      r |= v >> 1;
      return r;
    }
  }, {
    key: "intersectSegmentPoly2D",
    value: function intersectSegmentPoly2D(p0, p1, verts, nverts) {
      var result = new _IntersectResult["default"]();
      var EPS = 0.00000001;
      var dir = DetourCommon.vSub(p1, p0);
      var p0v = new _VectorPtr["default"](p0);

      for (var i = 0, j = nverts - 1; i < nverts; j = i++) {
        var vpj = new _VectorPtr["default"](verts, j * 3);
        var edge = DetourCommon.vSub(new _VectorPtr["default"](verts, i * 3), vpj);
        var diff = DetourCommon.vSub(p0v, vpj);
        var n = DetourCommon.vDot2D(edge, diff);

        var _d2 = DetourCommon.vDot2D(dir, edge);

        if (Math.abs(_d2) < EPS) {
          // S is nearly parallel to this edge
          if (n < 0) return result;else continue;
        }

        var _t = n / _d2;

        if (_d2 < 0) {
          // segment S is entering across this edge
          if (_t > result.tmin) {
            result.tmin = _t;
            result.segMin = j; // S enters after leaving polygon

            if (result.tmin > result.tmax) return result;
          }
        } else {
          // segment S is leaving across this edge
          if (_t < result.tmax) {
            result.tmax = _t;
            result.segMax = j; // S leaves before entering polygon

            if (result.tmax < result.tmin) return result;
          }
        }
      }

      result.intersects = true;
      return result;
    }
  }, {
    key: "distancePtSegSqr2D4",
    value: function distancePtSegSqr2D4(pt, verts, p, q) {
      var pqx = verts[q + 0] - verts[p + 0];
      var pqz = verts[q + 2] - verts[p + 2];
      var dx = pt[0] - verts[p + 0];
      var dz = pt[2] - verts[p + 2];
      var d = pqx * pqx + pqz * pqz;
      var t = pqx * dx + pqz * dz;
      if (d > 0) t /= d;
      if (t < 0) t = 0;else if (t > 1) t = 1;
      dx = verts[p + 0] + t * pqx - pt[0];
      dz = verts[p + 2] + t * pqz - pt[2];
      return [dx * dx + dz * dz, t];
    }
  }, {
    key: "oppositeTile",
    value: function oppositeTile(side) {
      return side + 4 & 0x7;
    }
  }, {
    key: "vperpXZ",
    value: function vperpXZ(a, b) {
      return a[0] * b[2] - a[2] * b[0];
    }
  }, {
    key: "intersectSegSeg2D",
    value: function intersectSegSeg2D(ap, aq, bp, bq) {
      var u = DetourCommon.vSub(aq, ap);
      var v = DetourCommon.vSub(bq, bp);
      var w = DetourCommon.vSub(ap, bp);
      d = DetourCommon.vperpXZ(u, v);
      if (Math.abs(d) < 1e-6) return [false, 0, 0];
      s = DetourCommon.vperpXZ(v, w) / d;
      t = DetourCommon.vperpXZ(u, w) / d;
      return [true, s, t];
    }
  }, {
    key: "vScale",
    value: function vScale(_in, scale) {
      var out = new Array(3);
      out[0] = _in[0] * scale;
      out[1] = _in[1] * scale;
      out[2] = _in[2] * scale;
      return out;
    }
  }]);

  return DetourCommon;
}();

_defineProperty(DetourCommon, "EPS", 1e-4);

_defineProperty(DetourCommon, "thr", DetourCommon.sqr(1.0 / 16384.0));

_defineProperty(DetourCommon, "eps", 1e-4);

var _default = DetourCommon;
exports["default"] = _default;