(function (global, factory) {
  typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
  typeof define === 'function' && define.amd ? define(['exports'], factory) :
  (global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.crowded = {}));
}(this, (function (exports) { 'use strict';

  function _classCallCheck(instance, Constructor) {
    if (!(instance instanceof Constructor)) {
      throw new TypeError("Cannot call a class as a function");
    }
  }

  function _defineProperties(target, props) {
    for (var i = 0; i < props.length; i++) {
      var descriptor = props[i];
      descriptor.enumerable = descriptor.enumerable || false;
      descriptor.configurable = true;
      if ("value" in descriptor) descriptor.writable = true;
      Object.defineProperty(target, descriptor.key, descriptor);
    }
  }

  function _createClass(Constructor, protoProps, staticProps) {
    if (protoProps) _defineProperties(Constructor.prototype, protoProps);
    if (staticProps) _defineProperties(Constructor, staticProps);
    return Constructor;
  }

  function _defineProperty(obj, key, value) {
    if (key in obj) {
      Object.defineProperty(obj, key, {
        value: value,
        enumerable: true,
        configurable: true,
        writable: true
      });
    } else {
      obj[key] = value;
    }

    return obj;
  }

  function _slicedToArray(arr, i) {
    return _arrayWithHoles(arr) || _iterableToArrayLimit(arr, i) || _unsupportedIterableToArray(arr, i) || _nonIterableRest();
  }

  function _toConsumableArray(arr) {
    return _arrayWithoutHoles(arr) || _iterableToArray(arr) || _unsupportedIterableToArray(arr) || _nonIterableSpread();
  }

  function _arrayWithoutHoles(arr) {
    if (Array.isArray(arr)) return _arrayLikeToArray(arr);
  }

  function _arrayWithHoles(arr) {
    if (Array.isArray(arr)) return arr;
  }

  function _iterableToArray(iter) {
    if (typeof Symbol !== "undefined" && Symbol.iterator in Object(iter)) return Array.from(iter);
  }

  function _iterableToArrayLimit(arr, i) {
    if (typeof Symbol === "undefined" || !(Symbol.iterator in Object(arr))) return;
    var _arr = [];
    var _n = true;
    var _d = false;
    var _e = undefined;

    try {
      for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) {
        _arr.push(_s.value);

        if (i && _arr.length === i) break;
      }
    } catch (err) {
      _d = true;
      _e = err;
    } finally {
      try {
        if (!_n && _i["return"] != null) _i["return"]();
      } finally {
        if (_d) throw _e;
      }
    }

    return _arr;
  }

  function _unsupportedIterableToArray(o, minLen) {
    if (!o) return;
    if (typeof o === "string") return _arrayLikeToArray(o, minLen);
    var n = Object.prototype.toString.call(o).slice(8, -1);
    if (n === "Object" && o.constructor) n = o.constructor.name;
    if (n === "Map" || n === "Set") return Array.from(o);
    if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen);
  }

  function _arrayLikeToArray(arr, len) {
    if (len == null || len > arr.length) len = arr.length;

    for (var i = 0, arr2 = new Array(len); i < len; i++) arr2[i] = arr[i];

    return arr2;
  }

  function _nonIterableSpread() {
    throw new TypeError("Invalid attempt to spread non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.");
  }

  function _nonIterableRest() {
    throw new TypeError("Invalid attempt to destructure non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.");
  }

  function _createForOfIteratorHelper(o, allowArrayLike) {
    var it;

    if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) {
      if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") {
        if (it) o = it;
        var i = 0;

        var F = function () {};

        return {
          s: F,
          n: function () {
            if (i >= o.length) return {
              done: true
            };
            return {
              done: false,
              value: o[i++]
            };
          },
          e: function (e) {
            throw e;
          },
          f: F
        };
      }

      throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.");
    }

    var normalCompletion = true,
        didErr = false,
        err;
    return {
      s: function () {
        it = o[Symbol.iterator]();
      },
      n: function () {
        var step = it.next();
        normalCompletion = step.done;
        return step;
      },
      e: function (e) {
        didErr = true;
        err = e;
      },
      f: function () {
        try {
          if (!normalCompletion && it.return != null) it.return();
        } finally {
          if (didErr) throw err;
        }
      }
    };
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
  /// Configuration parameters for a crowd agent.
  /// @ingroup crowd
  var CrowdAgentParams = function CrowdAgentParams() {
    _classCallCheck(this, CrowdAgentParams);

    _defineProperty(this, "radius", 0);

    _defineProperty(this, "height", 0);

    _defineProperty(this, "maxAcceleration", 0);

    _defineProperty(this, "maxSpeed", 0);

    _defineProperty(this, "collisionQueryRange", 0);

    _defineProperty(this, "pathOptimizationRange", 0);

    _defineProperty(this, "separationWeight", 0);

    _defineProperty(this, "updateFlags", 0);

    _defineProperty(this, "obstacleAvoidanceType", 0);

    _defineProperty(this, "queryFilterType", 0);

    _defineProperty(this, "userData", null);
  };

  _defineProperty(CrowdAgentParams, "DT_CROWD_ANTICIPATE_TURNS", 1);

  _defineProperty(CrowdAgentParams, "DT_CROWD_OBSTACLE_AVOIDANCE", 2);

  _defineProperty(CrowdAgentParams, "DT_CROWD_SEPARATION", 4);

  _defineProperty(CrowdAgentParams, "DT_CROWD_OPTIMIZE_VIS", 8);

  _defineProperty(CrowdAgentParams, "DT_CROWD_OPTIMIZE_TOPO", 16);

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
  var NavMeshParams = function NavMeshParams() {
    _classCallCheck(this, NavMeshParams);

    _defineProperty(this, "orig", new Array(3));

    _defineProperty(this, "tileWidth", 0);

    _defineProperty(this, "tileHeight", 0);

    _defineProperty(this, "maxTiles", 0);

    _defineProperty(this, "maxPolys", 0);
  };

  var IntersectResult = function IntersectResult() {
    _classCallCheck(this, IntersectResult);

    _defineProperty(this, "intersects", false);

    _defineProperty(this, "tmin", 0);

    _defineProperty(this, "tmax", 1);

    _defineProperty(this, "segMin", -1);

    _defineProperty(this, "segMax", -1);
  };

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
  var VectorPtr$1 = /*#__PURE__*/function () {
    // constructor( array) {
    // 		this(array, 0);
    // 	}
    function VectorPtr(array) {
      var index = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;

      _classCallCheck(this, VectorPtr);

      _defineProperty(this, "array", void 0);

      _defineProperty(this, "index", void 0);

      this.array = array;
      this.index = index;
    }

    _createClass(VectorPtr, [{
      key: "get",
      value: function get(offset) {
        return array[index + offset];
      }
    }]);

    return VectorPtr;
  }();

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
        var result = new IntersectResult();
        var EPS = 0.00000001;
        var dir = DetourCommon.vSub(p1, p0);
        var p0v = new VectorPtr$1(p0);

        for (var i = 0, j = nverts - 1; i < nverts; j = i++) {
          var vpj = new VectorPtr$1(verts, j * 3);
          var edge = DetourCommon.vSub(new VectorPtr$1(verts, i * 3), vpj);
          var diff = DetourCommon.vSub(p0v, vpj);
          var n = DetourCommon.vPerp2D(edge, diff);

          var _d2 = DetourCommon.vPerp2D(dir, edge);

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

  /**
   * Defines a navigation mesh tile.
   */

  var MeshTile =
  /** Counter describing modifications to the tile. */

  /** The tile data. */

  /** The tile links. */

  /** Index to the next free link. */

  /** Tile flags. (See: #dtTileFlags) */

  /** The next free tile, or the next tile in the spatial grid. */
  function MeshTile(index) {
    _classCallCheck(this, MeshTile);

    _defineProperty(this, "index", 0);

    _defineProperty(this, "salt", 0);

    _defineProperty(this, "data", null);

    _defineProperty(this, "links", []);

    _defineProperty(this, "linksFreeList", NavMesh.DT_NULL_LINK);

    _defineProperty(this, "flags", 0);

    _defineProperty(this, "next", null);

    this.index = index;
  };

  /** Defines a polyogn within a dtPoly object. */

  var Poly = /*#__PURE__*/function () {
    /** The polygon is a standard convex polygon that is part of the surface of the mesh. */

    /** The polygon is an off-mesh connection consisting of two vertices. */

    /** Index to first link in linked list. (Or #DT_NULL_LINK if there is no link.) */

    /** The indices of the polygon's vertices. The actual vertices are located in MeshTile::verts. */

    /** Packed data representing neighbor polygons references and flags for each edge. */

    /** The user defined polygon flags. */

    /** The number of vertices in the polygon. */

    /**
     * The bit packed area id and polygon type.
     * 
     * @note Use the structure's set and get methods to access this value.
     */
    function Poly(index, maxVertsPerPoly) {
      _classCallCheck(this, Poly);

      _defineProperty(this, "index", void 0);

      _defineProperty(this, "firstLink", 0);

      _defineProperty(this, "verts", []);

      _defineProperty(this, "neis", []);

      _defineProperty(this, "flags", 0);

      _defineProperty(this, "vertCount", 0);

      _defineProperty(this, "areaAndtype", void 0);

      this.index = index;
      this.firstLink = NavMesh.DT_NULL_LINK;
      this.verts = new Array(maxVertsPerPoly);
      this.neis = new Array(maxVertsPerPoly);

      for (var i = 0; i < this.verts.length; i++) {
        this.verts[i] = 0;
      }

      for (var _i = 0; _i < this.neis.length; _i++) {
        this.neis[_i] = 0;
      }
    }

    _createClass(Poly, [{
      key: "and",
      value: function and(v1, v2) {
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
    }, {
      key: "or",
      value: function or(v1, v2) {
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

    }, {
      key: "setArea",
      value: function setArea(a) {
        this.areaAndtype = this.or(this.and(this.areaAndtype, 0xc0), this.and(a, 0x3f));
      }
      /** Sets the polygon type. (See: #dtPolyTypes.) */

    }, {
      key: "setType",
      value: function setType(t) {
        this.areaAndtype = this.or(this.and(this.areaAndtype, 0x3f) | t << 6);
      }
      /** Gets the user defined area id. */

    }, {
      key: "getArea",
      value: function getArea() {
        return this.areaAndtype & 0x3;
      }
      /** Gets the polygon type. (See: #dtPolyTypes) */

    }, {
      key: "getType",
      value: function getType() {
        return this.areaAndtype >> 6;
      }
    }]);

    return Poly;
  }();

  _defineProperty(Poly, "DT_POLYTYPE_GROUND", 0);

  _defineProperty(Poly, "DT_POLYTYPE_OFFMESH_CONNECTION", 1);

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
  var Link = function Link() {
    _classCallCheck(this, Link);

    _defineProperty(this, "ref", 0);

    _defineProperty(this, "next", 0);

    _defineProperty(this, "edge", 0);

    _defineProperty(this, "side", 0);

    _defineProperty(this, "bmin", 0);

    _defineProperty(this, "bmax", 0);
  };

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
  var ClosestPointOnPolyResult = /*#__PURE__*/function () {
    function ClosestPointOnPolyResult(posOverPoly, closest) {
      _classCallCheck(this, ClosestPointOnPolyResult);

      _defineProperty(this, "posOverPoly", false);

      _defineProperty(this, "closest", []);

      this.posOverPoly = posOverPoly;
      this.closest = closest;
    }
    /** Returns true if the position is over the polygon. */


    _createClass(ClosestPointOnPolyResult, [{
      key: "isPosOverPoly",
      value: function isPosOverPoly() {
        return this.posOverPoly;
      }
      /** Returns the closest poPoly on the polygon. [(x, y, z)] */

    }, {
      key: "getClosest",
      value: function getClosest() {
        return this.closest;
      }
    }]);

    return ClosestPointOnPolyResult;
  }();

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
  var FindNearestPolyResult = /*#__PURE__*/function () {
    function FindNearestPolyResult(nearestRef, nearestPos) {
      _classCallCheck(this, FindNearestPolyResult);

      _defineProperty(this, "nearestRef", 0);

      _defineProperty(this, "nearestPos", []);

      this.nearestRef = nearestRef;
      this.nearestPos = nearestPos;
    }
    /** Returns the reference id of the nearest polygon. */


    _createClass(FindNearestPolyResult, [{
      key: "getNearestRef",
      value: function getNearestRef() {
        return this.nearestRef;
      }
      /** Returns the nearest poPoly on the polygon. [opt] [(x, y, z)] */

    }, {
      key: "getNearestPos",
      value: function getNearestPos() {
        return this.nearestPos;
      }
    }]);

    return FindNearestPolyResult;
  }();

  function arraycopy(one, oneStart, two, twoStart, len) {
    for (var i = 0; i < len; i++) {
      two[twoStart + i] = one[oneStart + i];
    }
  }

  var NavMesh = /*#__PURE__*/function () {
    _createClass(NavMesh, [{
      key: "getMaxTiles",
      /// A flag that indicates that an entity links to an external entity.
      /// (E.g. A polygon edge is a portal that links to another polygon.)
      /// A value that indicates the entity does not link to anything.
      // static DT_NULL_LINK = 0xffffffff;
      /// A flag that indicates that an off-mesh connection can be traversed in
      /// both directions. (Is bidirectional.)
      /// The maximum number of user defined area ids.
      /// Limit raycasting during any angle pahfinding
      /// The limit is given as a multiple of the character radius
      /// < Current initialization params. TODO: do not store this info twice.
      /// < Origin of the tile (0,0)
      //  m_orig[3]; ///< Origin of the tile (0,0)
      /// < Dimensions of each tile.
      /// < Max number of tiles.
      /// < Tile hash lookup size (must be pot).
      /// < Tile hash lookup mask.
      /// < Tile hash lookup.
      /// < Freelist of tiles.
      /// < List of tiles.

      /** The maximum number of vertices per navigation polygon. */

      /**
       * The maximum number of tiles supported by the navigation mesh.
       * 
       * @return The maximum number of tiles supported by the navigation mesh.
       */
      value: function getMaxTiles() {
        return this.m_maxTiles;
      }
      /**
       * Returns tile in the tile array.
       */

    }, {
      key: "getTile",
      value: function getTile(i) {
        return this.m_tiles[i];
      }
      /**
       * Gets the polygon reference for the tile's base polygon.
       * 
       * @param tile
       *            The tile.
       * @return The polygon reference for the base polygon in the specified tile.
       */

    }, {
      key: "getPolyRefBase",
      value: function getPolyRefBase(tile) {
        if (tile == null) return 0;
        var it = tile.index;
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

    }, {
      key: "allocLink",
      value: function allocLink(tile) {
        if (tile.linksFreeList == NavMesh.DT_NULL_LINK) {
          var _link = new Link();

          _link.next = NavMesh.DT_NULL_LINK;
          tile.links.push(_link);
          return tile.links.length - 1;
        }

        var link = tile.linksFreeList;
        tile.linksFreeList = tile.links[link].next;
        return link;
      }
    }, {
      key: "freeLink",
      value: function freeLink(tile, link) {
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

    }, {
      key: "calcTileLoc",
      value: function calcTileLoc(pos) {
        var tx = Math.floor((pos[0] - this.m_orig[0]) / this.m_tileWidth);
        var ty = Math.floor((pos[2] - this.m_orig[2]) / this.m_tileHeight);
        return [tx, ty];
      }
    }, {
      key: "getTileAndPolyByRef",
      value: function getTileAndPolyByRef(ref) {
        if (ref == 0) {
          throw new IllegalArgumentException("ref = 0");
        }

        var saltitip = NavMesh.decodePolyId(ref);
        var salt = saltitip[0];
        var it = saltitip[1];
        var ip = saltitip[2];
        if (it >= this.m_maxTiles) throw new IllegalArgumentException("tile > m_maxTiles");
        if (this.m_tiles[it].salt != salt || this.m_tiles[it].data.header == null) throw new IllegalArgumentException("Invalid salt or header");
        if (ip >= this.m_tiles[it].data.header.polyCount) throw new IllegalArgumentException("poly > polyCount");
        return [this.m_tiles[it], this.m_tiles[it].data.polys[ip]];
      } /// @par
      ///
      /// @warning Only use this function if it is known that the provided polygon
      /// reference is valid. This function is faster than #getTileAndPolyByRef,
      /// but
      /// it does not validate the reference.

    }, {
      key: "getTileAndPolyByRefUnsafe",
      value: function getTileAndPolyByRefUnsafe(ref) {
        var saltitip = NavMesh.decodePolyId(ref);
        var it = saltitip[1];
        var ip = saltitip[2];
        return [this.m_tiles[it], this.m_tiles[it].data.polys[ip]];
      }
    }, {
      key: "isValidPolyRef",
      value: function isValidPolyRef(ref) {
        if (ref == 0) return false;
        var saltitip = NavMesh.decodePolyId(ref);
        var salt = saltitip[0];
        var it = saltitip[1];
        var ip = saltitip[2];
        if (it >= this.m_maxTiles) return false;
        if (this.m_tiles[it].salt != salt || this.m_tiles[it].data == null) return false;
        if (ip >= this.m_tiles[it].data.header.polyCount) return false;
        return true;
      }
    }, {
      key: "getParams",
      value: function getParams() {
        return this.m_params;
      }
    }], [{
      key: "lshift",
      value: function lshift(num, bits) {
        return num * Math.pow(2, bits);
      }
    }, {
      key: "rshift",
      value: function rshift(num, bits) {
        return Math.floor(num / Math.pow(2, bits));
      } //https://stackoverflow.com/a/43666199/10047920

    }, {
      key: "and",
      value: function and(v1, v2) {
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
    }, {
      key: "or",
      value: function or(v1, v2) {
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
    }, {
      key: "encodePolyId",
      value: function encodePolyId(salt, it, ip) {
        var a = NavMesh.lshift(salt, NavMesh.DT_POLY_BITS + NavMesh.DT_TILE_BITS);
        var b = NavMesh.lshift(it, NavMesh.DT_POLY_BITS);
        return NavMesh.or(NavMesh.or(a, b), ip);
      } /// Decodes a standard polygon reference.
      /// @note This function is generally meant for internal use only.
      /// @param[in] ref The polygon reference to decode.
      /// @param[out] salt The tile's salt value.
      /// @param[out] it The index of the tile.
      /// @param[out] ip The index of the polygon within the tile.
      /// @see #encodePolyId

    }, {
      key: "decodePolyId",
      value: function decodePolyId(ref) {
        var salt;
        var it;
        var ip;
        var saltMask = NavMesh.lshift(1, NavMesh.DT_SALT_BITS) - 1;
        var tileMask = NavMesh.lshift(1, NavMesh.DT_TILE_BITS) - 1;
        var polyMask = NavMesh.lshift(1, NavMesh.DT_POLY_BITS) - 1;
        salt = Math.floor(NavMesh.and(NavMesh.rshift(ref, NavMesh.DT_POLY_BITS + NavMesh.DT_TILE_BITS), saltMask));
        it = Math.floor(NavMesh.and(NavMesh.rshift(ref, NavMesh.DT_POLY_BITS), tileMask));
        ip = Math.floor(NavMesh.and(ref, polyMask));
        return [salt, it, ip];
      } /// Extracts a tile's salt value from the specified polygon reference.
      /// @note This function is generally meant for internal use only.
      /// @param[in] ref The polygon reference.
      /// @see #encodePolyId

    }, {
      key: "decodePolyIdSalt",
      value: function decodePolyIdSalt(ref) {
        var saltMask = (1 << NavMesh.DT_SALT_BITS) - 1;
        return Math.floor(ref >> NavMesh.DT_POLY_BITS + NavMesh.DT_TILE_BITS & saltMask);
      } /// Extracts the tile's index from the specified polygon reference.
      /// @note This function is generally meant for internal use only.
      /// @param[in] ref The polygon reference.
      /// @see #encodePolyId

    }, {
      key: "decodePolyIdTile",
      value: function decodePolyIdTile(ref) {
        var tileMask = (1 << NavMesh.DT_TILE_BITS) - 1;
        return Math.floor(ref >> NavMesh.DT_POLY_BITS & tileMask);
      } /// Extracts the polygon's index (within its tile) from the specified
      /// polygon reference.
      /// @note This function is generally meant for internal use only.
      /// @param[in] ref The polygon reference.
      /// @see #encodePolyId

    }, {
      key: "decodePolyIdPoly",
      value: function decodePolyIdPoly(ref) {
        var polyMask = (1 << NavMesh.DT_POLY_BITS) - 1;
        return Math.floor(ref & polyMask);
      }
    }]);

    function NavMesh(one, maxVertsPerPoly, flags) {
      _classCallCheck(this, NavMesh);

      _defineProperty(this, "m_params", null);

      _defineProperty(this, "m_orig", []);

      _defineProperty(this, "m_tileWidth", 0);

      _defineProperty(this, "m_tileHeight", 0);

      _defineProperty(this, "m_maxTiles", 0);

      _defineProperty(this, "m_tileLutSize", 0);

      _defineProperty(this, "m_tileLutMask", 0);

      _defineProperty(this, "m_posLookup", []);

      _defineProperty(this, "m_nextFree", null);

      _defineProperty(this, "m_tiles", []);

      _defineProperty(this, "m_maxVertPerPoly", 0);

      _defineProperty(this, "m_tileCount", 0);

      if (flags || flags == 0) {
        this._constructor(NavMesh.getNavMeshParams(one), maxVertsPerPoly);

        this.addTile(one, flags, 0);
      } else {
        this._constructor(one, maxVertsPerPoly);
      }
    }

    _createClass(NavMesh, [{
      key: "_constructor",
      value: function _constructor(params, maxVertsPerPoly) {
        this.m_params = params;
        this.m_orig = params.orig;
        this.m_tileWidth = params.tileWidth;
        this.m_tileHeight = params.tileHeight; // Init tiles

        this.m_maxTiles = params.maxTiles;
        this.m_maxVertPerPoly = maxVertsPerPoly;
        var lutsize = DetourCommon.nextPow2(params.maxTiles / 4);
        if (lutsize == 0) lutsize = 1;
        this.m_tileLutSize = lutsize;
        this.m_tileLutMask = this.m_tileLutSize - 1;
        this.m_tiles = new Array(this.m_maxTiles);
        this.m_posLookup = new Array(this.m_tileLutSize);
        this.m_nextFree = null;

        for (var i = this.m_maxTiles - 1; i >= 0; --i) {
          this.m_tiles[i] = new MeshTile(i);
          this.m_tiles[i].salt = 1;
          this.m_tiles[i].next = this.m_nextFree;
          this.m_nextFree = this.m_tiles[i];
        }
      }
    }, {
      key: "queryPolygonsInTile",
      // TODO: These methods are duplicates from dtNavMeshQuery, but are needed
      // for off-mesh connection finding.
      value: function queryPolygonsInTile(tile, qmin, qmax) {
        var polys = [];

        if (tile.data.bvTree != null) {
          var nodeIndex = 0;
          var tbmin = tile.data.header.bmin;
          var tbmax = tile.data.header.bmax;
          var qfac = tile.data.header.bvQuantFactor; // Calculate quantized box

          var _bmin = new Array(3);

          var _bmax = new Array(3); // dtClamp query box to world box.


          var _minx = DetourCommon.clamp(qmin[0], tbmin[0], tbmax[0]) - tbmin[0];

          var miny = DetourCommon.clamp(qmin[1], tbmin[1], tbmax[1]) - tbmin[1];
          var minz = DetourCommon.clamp(qmin[2], tbmin[2], tbmax[2]) - tbmin[2];

          var _maxx = DetourCommon.clamp(qmax[0], tbmin[0], tbmax[0]) - tbmin[0];

          var maxy = DetourCommon.clamp(qmax[1], tbmin[1], tbmax[1]) - tbmin[1];
          var maxz = DetourCommon.clamp(qmax[2], tbmin[2], tbmax[2]) - tbmin[2]; // Quantize

          _bmin[0] = NavMesh.and(Math.floor(qfac * _minx), 0xfffe);
          _bmin[1] = NavMesh.and(Math.floor(qfac * miny), 0xfffe);
          _bmin[2] = NavMesh.and(Math.floor(qfac * minz), 0xfffe);
          _bmax[0] = NavMesh.or(Math.floor(qfac * _maxx + 1), 1);
          _bmax[1] = NavMesh.or(Math.floor(qfac * maxy + 1), 1);
          _bmax[2] = NavMesh.or(Math.floor(qfac * maxz + 1), 1); // Traverse tree

          var base = this.getPolyRefBase(tile);
          var end = tile.data.header.bvNodeCount;

          while (nodeIndex < end) {
            var node = tile.data.bvTree[nodeIndex];
            var overlap = DetourCommon.overlapQuantBounds(_bmin, _bmax, node.bmin, node.bmax);
            var isLeafNode = node.i >= 0;

            if (isLeafNode && overlap) {
              polys.push(NavMesh.or(base, node.i));
            }

            if (overlap || isLeafNode) nodeIndex++;else {
              var escapeIndex = -node.i;
              nodeIndex += escapeIndex;
            }
          }

          return polys;
        } else {
          bmin = [null, null, null];
          bmax = [null, null, null];

          var _base = this.getPolyRefBase(tile);

          for (var i = 0; i < tile.data.header.polyCount; ++i) {
            var p = tile.data.polys[i]; // Do not return off-mesh connection polygons.

            if (p.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) continue; // Calc polygon bounds.

            var v = p.verts[0] * 3;
            DetourCommon.vCopy(bmin, tile.data.verts, v);
            DetourCommon.vCopy(bmax, tile.data.verts, v);

            for (var j = 1; j < p.vertCount; ++j) {
              v = p.verts[j] * 3;
              DetourCommon.vMin(bmin, tile.data.verts, v);
              DetourCommon.vMax(bmax, tile.data.verts, v);
            }

            if (overlapBounds(qmin, qmax, bmin, bmax)) {
              polys.push(NavMesh.or(_base, i));
            }
          }

          return polys;
        }
      } /// Adds a tile to the navigation mesh.
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

    }, {
      key: "addTile",
      value: function addTile(data, flags, lastRef) {
        // Make sure the data is in right format.
        var header = data.header; // Make sure the location is free.

        if (this.getTileAt(header.x, header.y, header.layer) != null) throw new RuntimeException("Tile already exists"); // Allocate a tile.

        var tile = null;

        if (lastRef == 0) {
          if (this.m_nextFree != null) {
            tile = this.m_nextFree;
            this.m_nextFree = tile.next;
            tile.next = null;
            this.m_tileCount++;
          }
        } else {
          // Try to relocate the tile to specific index with same salt.
          var tileIndex = decodePolyIdTile(lastRef);
          if (tileIndex >= this.m_maxTiles) throw new RuntimeException("Tile index too high"); // Try to find the specific tile id from the free list.

          var target = this.m_tiles[tileIndex];
          var prev = null;
          tile = this.m_nextFree;

          while (tile != null && tile != target) {
            prev = tile;
            tile = tile.next;
          } // Could not find the correct location.


          if (tile != target) throw new RuntimeException("Could not find tile"); // Remove from freelist

          if (prev == null) this.m_nextFree = tile.next;else prev.next = tile.next; // Restore salt.

          tile.salt = decodePolyIdSalt(lastRef);
        } // Make sure we could allocate a tile.


        if (tile == null) throw new RuntimeException("Could not allocate a tile");
        tile.data = data;
        tile.flags = flags;
        tile.links = []; // Insert tile into the position lut.

        var h = NavMesh.computeTileHash(header.x, header.y, this.m_tileLutMask);
        tile.next = this.m_posLookup[h];
        this.m_posLookup[h] = tile; // Patch header pointers.
        // If there are no items in the bvtree, reset the tree pointer.

        if (tile.data.bvTree != null && tile.data.bvTree.length == 0) tile.data.bvTree = null; // Init tile.

        this.connectIntLinks(tile); // Base off-mesh connections to their starting polygons and connect connections inside the tile.

        this.baseOffMeshLinks(tile);
        this.connectExtOffMeshLinks(tile, tile, -1); // Connect with layers in current tile.

        var neis = this.getTilesAt(header.x, header.y);

        for (var j = 0; j < neis.length; ++j) {
          if (neis[j] == tile) {
            continue;
          }

          connectExtLinks(tile, neis[j], -1);
          connectExtLinks(neis[j], tile, -1);
          this.connectExtOffMeshLinks(tile, neis[j], -1);
          this.connectExtOffMeshLinks(neis[j], tile, -1);
        } // Connect with neighbour tiles.


        for (var i = 0; i < 8; ++i) {
          neis = this.getNeighbourTilesAt(header.x, header.y, i);

          for (var _j = 0; _j < neis.length; ++_j) {
            connectExtLinks(tile, neis[_j], i);
            connectExtLinks(neis[_j], tile, DetourCommon.oppositeTile(i));
            this.connectExtOffMeshLinks(tile, neis[_j], i);
            this.connectExtOffMeshLinks(neis[_j], tile, DetourCommon.oppositeTile(i));
          }
        }

        return this.getTileRef(tile);
      } /// Removes the specified tile from the navigation mesh.
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

    }, {
      key: "removeTile",
      value: function removeTile(ref) {
        if (ref == 0) {
          return null;
        }

        var tileIndex = decodePolyIdTile(ref);
        var tileSalt = decodePolyIdSalt(ref);
        if (tileIndex >= this.m_maxTiles) throw new RuntimeException("Invalid tile index");
        var tile = this.m_tiles[tileIndex];
        if (tile.salt != tileSalt) throw new RuntimeException("Invalid tile salt"); // Remove tile from hash lookup.

        var h = NavMesh.computeTileHash(tile.data.header.x, tile.data.header.y, this.m_tileLutMask);
        var prev = null;
        var cur = this.m_posLookup[h];

        while (cur != null) {
          if (cur == tile) {
            if (prev != null) prev.next = cur.next;else this.m_posLookup[h] = cur.next;
            break;
          }

          prev = cur;
          cur = cur.next;
        } // Remove connections to neighbour tiles.
        // Create connections with neighbour tiles.
        // Disconnect from other layers in current tile.


        nneis = this.getTilesAt(tile.data.header.x, tile.data.header.y);

        var _iterator = _createForOfIteratorHelper(nneis),
            _step;

        try {
          for (_iterator.s(); !(_step = _iterator.n()).done;) {
            var _j2 = _step.value;
            if (_j2 == tile) continue;
            unconnectLinks(_j2, tile);
          } // Disconnect from neighbour tiles.

        } catch (err) {
          _iterator.e(err);
        } finally {
          _iterator.f();
        }

        for (var i = 0; i < 8; ++i) {
          nneis = this.getNeighbourTilesAt(tile.data.header.x, tile.data.header.y, i);

          var _iterator2 = _createForOfIteratorHelper(nneis),
              _step2;

          try {
            for (_iterator2.s(); !(_step2 = _iterator2.n()).done;) {
              var j = _step2.value;
              unconnectLinks(j, tile);
            }
          } catch (err) {
            _iterator2.e(err);
          } finally {
            _iterator2.f();
          }
        }

        var data = tile.data; // Reset tile.

        tile.data = null;
        tile.flags = 0;
        tile.links = []; // Update salt, salt should never be zero.

        tile.salt = tile.salt + 1 & (1 << NavMesh.DT_SALT_BITS) - 1;
        if (tile.salt == 0) tile.salt++; // Add to free list.

        tile.next = this.m_nextFree;
        this.m_nextFree = tile;
        this.m_tileCount--;
        return data;
      } /// Builds internal polygons links for a tile.

    }, {
      key: "connectIntLinks",
      value: function connectIntLinks(tile) {
        if (tile == null) return;
        var base = this.getPolyRefBase(tile);

        for (var i = 0; i < tile.data.header.polyCount; ++i) {
          var poly = tile.data.polys[i];
          poly.firstLink = NavMesh.DT_NULL_LINK;
          if (poly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) continue; // Build edge links backwards so that the links will be
          // in the linked list from lowest index to highest.

          for (var j = poly.vertCount - 1; j >= 0; --j) {
            // Skip hard and non-internal edges.
            if (poly.neis[j] == 0 || (poly.neis[j] & NavMesh.DT_EXT_LINK) != 0) continue;
            var idx = this.allocLink(tile);
            var link = tile.links[idx];
            link.ref = NavMesh.or(base, poly.neis[j] - 1);
            link.edge = j;
            link.side = 0xff;
            link.bmin = link.bmax = 0; // Add to linked list.

            link.next = poly.firstLink;
            poly.firstLink = idx;
          }
        }
      }
    }, {
      key: "unconnectLinks",
      value: function unconnectLinks(tile, target) {
        if (tile == null || target == null) return;
        var targetNum = decodePolyIdTile(this.getTileRef(target));

        for (var i = 0; i < tile.data.header.polyCount; ++i) {
          var poly = tile.data.polys[i];
          var j = poly.firstLink;
          var pj = NavMesh.DT_NULL_LINK;

          while (j != NavMesh.DT_NULL_LINK) {
            if (decodePolyIdTile(tile.links[j].ref) == targetNum) {
              // Remove link.
              var nj = tile.links[j].next;
              if (pj == NavMesh.DT_NULL_LINK) poly.firstLink = nj;else tile.links[pj].next = nj;
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
    }, {
      key: "connectExtLinks",
      value: function connectExtLinks(tile, target, side) {
        if (tile == null) return; // Connect border links.

        for (var i = 0; i < tile.data.header.polyCount; ++i) {
          var poly = tile.data.polys[i]; // Create new links.
          // short m = NavMesh.DT_EXT_LINK | (short)side;

          var nv = poly.vertCount;

          for (var j = 0; j < nv; ++j) {
            // Skip non-portal edges.
            if ((poly.neis[j] & NavMesh.DT_EXT_LINK) == 0) continue;
            var dir = poly.neis[j] & 0xff;
            if (side != -1 && dir != side) continue; // Create new links

            var va = poly.verts[j] * 3;
            var vb = poly.verts[(j + 1) % nv] * 3;
            connectedPolys = findConnectingPolys(tile.data.verts, va, vb, target, DetourCommon.oppositeTile(dir), 4);
            nei = connectedPolys[0];
            neia = connectedPolys[1];
            var nnei = connectedPolys.third;

            for (var k = 0; k < nnei; ++k) {
              var idx = this.allocLink(tile);
              var link = tile.links[idx];
              link.ref = nei[k];
              link.edge = j;
              link.side = dir;
              link.next = poly.firstLink;
              poly.firstLink = idx; // Compress portal limits to a byte value.

              if (dir == 0 || dir == 4) {
                tmin = (neia[k * 2 + 0] - tile.data.verts[va + 2]) / (tile.data.verts[vb + 2] - tile.data.verts[va + 2]);
                tmax = (neia[k * 2 + 1] - tile.data.verts[va + 2]) / (tile.data.verts[vb + 2] - tile.data.verts[va + 2]);

                if (tmin > tmax) {
                  temp = tmin;
                  tmin = tmax;
                  tmax = temp;
                }

                link.bmin = Math.floor(DetourCommon.clamp(tmin, 0.0, 1.0) * 255.0);
                link.bmax = Math.floor(DetourCommon.clamp(tmax, 0.0, 1.0) * 255.0);
              } else if (dir == 2 || dir == 6) {
                tmin = (neia[k * 2 + 0] - tile.data.verts[va]) / (tile.data.verts[vb] - tile.data.verts[va]);
                tmax = (neia[k * 2 + 1] - tile.data.verts[va]) / (tile.data.verts[vb] - tile.data.verts[va]);

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
    }, {
      key: "connectExtOffMeshLinks",
      value: function connectExtOffMeshLinks(tile, target, side) {
        if (tile == null) return; // Connect off-mesh links.
        // We are interested on links which land from target tile to this tile.

        var oppositeSide = side == -1 ? 0xff : DetourCommon.oppositeTile(side);

        for (var i = 0; i < target.data.header.offMeshConCount; ++i) {
          var targetCon = target.data.offMeshCons[i];
          if (targetCon.side != oppositeSide) continue;
          var targetPoly = target.data.polys[targetCon.poly]; // Skip off-mesh connections which start location could not be
          // connected at all.

          if (targetPoly.firstLink == NavMesh.DT_NULL_LINK) continue;
          var ext = [targetCon.rad, target.data.header.walkableClimb, targetCon.rad]; // Find polygon to connect to.

          var p = new Array(3);
          p[0] = targetCon.pos[3];
          p[1] = targetCon.pos[4];
          p[2] = targetCon.pos[5];
          var nearest = this.findNearestPolyInTile(tile, p, ext);
          var ref = nearest.getNearestRef();
          if (ref == 0) continue;
          var nearestPt = nearest.getNearestPos(); // findNearestPoly may return too optimistic results, further check
          // to make sure.

          if (DetourCommon.sqr(nearestPt[0] - p[0]) + DetourCommon.sqr(nearestPt[2] - p[2]) > DetourCommon.sqr(targetCon.rad)) continue; // Make sure the location is on curren mesh.

          target.data.verts[targetPoly.verts[1] * 3] = nearestPt[0];
          target.data.verts[targetPoly.verts[1] * 3 + 1] = nearestPt[1];
          target.data.verts[targetPoly.verts[1] * 3 + 2] = nearestPt[2]; // let off-mesh connection to target poly.

          var idx = this.allocLink(target);
          var link = target.links[idx];
          link.ref = ref;
          link.edge = 1;
          link.side = oppositeSide;
          link.bmin = link.bmax = 0; // Add to linked list.

          link.next = targetPoly.firstLink;
          targetPoly.firstLink = idx; // let target poly to off-mesh connection.

          if ((targetCon.flags & NavMesh.DT_OFFMESH_CON_BIDIR) != 0) {
            var tidx = this.allocLink(tile);
            var landPolyIdx = NavMesh.decodePolyIdPoly(ref);
            var landPoly = tile.data.polys[landPolyIdx];
            link = tile.links[tidx];
            link.ref = NavMesh.or(this.getPolyRefBase(target), targetCon.poly);
            link.edge = 0xff;
            link.side = side == -1 ? 0xff : side;
            link.bmin = link.bmax = 0; // Add to linked list.

            link.next = landPoly.firstLink;
            landPoly.firstLink = tidx;
          }
        }
      }
    }, {
      key: "findConnectingPolys",
      value: function findConnectingPolys(verts, va, vb, tile, side, maxcon) {
        if (tile == null) return [null, null, 0];
        var con = new Array(maxcon);
        var conarea = new Array(maxcon * 2);
        var amin = new Array(2);
        var amax = new Array(2);
        calcSlabEndPoints(verts, va, vb, amin, amax, side);
        apos = getSlabCoord(verts, va, side); // Remove links pointing to 'side' and compact the links array.

        var bmin = new Array(2);
        var bmax = new Array(2);
        var m = NavMesh.or(NavMesh.DT_EXT_LINK, side);
        var n = 0;
        var base = this.getPolyRefBase(tile);

        for (var i = 0; i < tile.data.header.polyCount; ++i) {
          var poly = tile.data.polys[i];
          var nv = poly.vertCount;

          for (var j = 0; j < nv; ++j) {
            // Skip edges which do not poPoly to the right side.
            if (poly.neis[j] != m) continue;
            var vc = poly.verts[j] * 3;
            var vd = poly.verts[(j + 1) % nv] * 3;
            bpos = getSlabCoord(tile.data.verts, vc, side); // Segments are not close enough.

            if (Math.abs(apos - bpos) > 0.01) continue; // Check if the segments touch.

            calcSlabEndPoints(tile.data.verts, vc, vd, bmin, bmax, side);
            if (!overlapSlabs(amin, amax, bmin, bmax, 0.01, tile.data.header.walkableClimb)) continue; // Add return value.

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
    }, {
      key: "overlapSlabs",
      value: function overlapSlabs(amin, amax, bmin, bmax, px, py) {
        // Check for horizontal overlap.
        // The segment is shrunken a little so that slabs which touch
        // at end points are not connected.
        minx = Math.max(amin[0] + px, bmin[0] + px);
        maxx = Math.min(amax[0] - px, bmax[0] - px);
        if (minx > maxx) return false; // Check vertical overlap.

        ad = (amax[1] - amin[1]) / (amax[0] - amin[0]);
        ak = amin[1] - ad * amin[0];
        bd = (bmax[1] - bmin[1]) / (bmax[0] - bmin[0]);
        bk = bmin[1] - bd * bmin[0];
        aminy = ad * minx + ak;
        amaxy = ad * maxx + ak;
        bminy = bd * minx + bk;
        bmaxy = bd * maxx + bk;
        dmin = bminy - aminy;
        dmax = bmaxy - amaxy; // Crossing segments always overlap.

        if (dmin * dmax < 0) return true; // Check for overlap at endpoints.

        thr = py * 2 * (py * 2);
        if (dmin * dmin <= thr || dmax * dmax <= thr) return true;
        return false;
      }
      /**
       * Builds internal polygons links for a tile.
       * 
       * @param tile
       */

    }, {
      key: "baseOffMeshLinks",
      value: function baseOffMeshLinks(tile) {
        if (tile == null) return;
        var base = this.getPolyRefBase(tile); // Base off-mesh connection start points.

        for (var i = 0; i < tile.data.header.offMeshConCount; ++i) {
          var con = tile.data.offMeshCons[i];
          var poly = tile.data.polys[con.poly];
          var ext = [con.rad, tile.data.header.walkableClimb, con.rad]; // Find polygon to connect to.

          var nearestPoly = this.findNearestPolyInTile(tile, con.pos, ext);
          var ref = nearestPoly.getNearestRef();
          if (ref == 0) continue;
          var p = con.pos; // First vertex

          var nearestPt = nearestPoly.getNearestPos(); // findNearestPoly may return too optimistic results, further check
          // to make sure.

          if (DetourCommon.sqr(nearestPt[0] - p[0]) + DetourCommon.sqr(nearestPt[2] - p[2]) > DetourCommon.sqr(con.rad)) continue; // Make sure the location is on current mesh.

          tile.data.verts[poly.verts[0] * 3] = nearestPt[0];
          tile.data.verts[poly.verts[0] * 3 + 1] = nearestPt[1];
          tile.data.verts[poly.verts[0] * 3 + 2] = nearestPt[2]; // let off-mesh connection to target poly.

          var idx = this.allocLink(tile);
          var link = tile.links[idx];
          link.ref = ref;
          link.edge = 0;
          link.side = 0xff;
          link.bmin = link.bmax = 0; // Add to linked list.

          link.next = poly.firstLink;
          poly.firstLink = idx; // Start end-poPoly is always connect back to off-mesh connection.

          var tidx = this.allocLink(tile);
          var landPolyIdx = NavMesh.decodePolyIdPoly(ref);
          var landPoly = tile.data.polys[landPolyIdx];
          link = tile.links[tidx];
          link.ref = NavMesh.or(base, con.poly);
          link.edge = 0xff;
          link.side = 0xff;
          link.bmin = link.bmax = 0; // Add to linked list.

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

    }, {
      key: "closestPointOnPoly",
      value: function closestPointOnPoly(ref, pos) {
        var tileAndPoly = this.getTileAndPolyByRefUnsafe(ref);
        var tile = tileAndPoly[0];
        var poly = tileAndPoly[1]; // Off-mesh connections don't have detail polygons.

        if (poly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) {
          var v0 = poly.verts[0] * 3;
          var v1 = poly.verts[1] * 3;
          var d0 = DetourCommon.vDist3(pos, tile.data.verts, v0);
          var d1 = DetourCommon.vDist3(pos, tile.data.verts, v1);
          var u = d0 / (d0 + d1);

          var _closest = DetourCommon.vLerp4(tile.data.verts, v0, v1, u);

          return new ClosestPointOnPolyResult(false, _closest);
        } // Clamp poPoly to be inside the polygon.


        var verts = new Array(this.m_maxVertPerPoly * 3);
        var edged = new Array(this.m_maxVertPerPoly);
        var edget = new Array(this.m_maxVertPerPoly);
        var nv = poly.vertCount;

        for (var i = 0; i < nv; ++i) {
          arraycopy(tile.data.verts, poly.verts[i] * 3, verts, i * 3, 3);
        }

        var posOverPoly = false;
        var closest = new Array(3);
        DetourCommon.vCopy(closest, pos);

        if (!DetourCommon.distancePtPolyEdgesSqr(pos, verts, nv, edged, edget)) {
          // PoPoly is outside the polygon, dtClamp to nearest edge.
          var _dmin = edged[0];
          var imin = 0;

          for (var _i = 1; _i < nv; ++_i) {
            if (edged[_i] < _dmin) {
              _dmin = edged[_i];
              imin = _i;
            }
          }

          var va = imin * 3;
          var vb = (imin + 1) % nv * 3;
          closest = DetourCommon.vLerp4(verts, va, vb, edget[imin]);
          posOverPoly = false;
        } else {
          posOverPoly = true;
        } // Find height at the location.


        var ip = poly.index;

        if (tile.data.detailMeshes != null && tile.data.detailMeshes.length > ip) {
          var pd = tile.data.detailMeshes[ip];

          for (var j = 0; j < pd.triCount; ++j) {
            var t = (pd.triBase + j) * 4;
            var v = []; //Was [3][]

            for (var k = 0; k < 3; ++k) {
              if (tile.data.detailTris[t + k] < poly.vertCount) {
                var index = poly.verts[tile.data.detailTris[t + k]] * 3;
                v[k] = [tile.data.verts[index], tile.data.verts[index + 1], tile.data.verts[index + 2]];
              } else {
                var _index = (pd.vertBase + (tile.data.detailTris[t + k] - poly.vertCount)) * 3;

                v[k] = [tile.data.detailVerts[_index], tile.data.detailVerts[_index + 1], tile.data.detailVerts[_index + 2]];
              }
            }

            var heightResult = DetourCommon.closestHeightPointTriangle(closest, v[0], v[1], v[2]);

            if (heightResult[0]) {
              closest[1] = heightResult[1];
              break;
            }
          }
        }

        return new ClosestPointOnPolyResult(posOverPoly, closest);
      }
    }, {
      key: "findNearestPolyInTile",
      value: function findNearestPolyInTile(tile, center, extents) {
        var nearestPt = null;
        var bmin = DetourCommon.vSub(center, extents);
        var bmax = DetourCommon.vAdd(center, extents); // Get nearby polygons from proximity grid.

        var polys = this.queryPolygonsInTile(tile, bmin, bmax); // Find nearest polygon amongst the nearby polygons.

        var nearest = 0;
        var nearestDistanceSqr = Number.MAX_VALUE;

        for (var i = 0; i < polys.length; ++i) {
          var ref = polys[i];
          var d = 0;
          var cpp = this.closestPointOnPoly(ref, center);
          var posOverPoly = cpp.isPosOverPoly();
          var closestPtPoly = cpp.getClosest(); // If a poPoly is directly over a polygon and closer than
          // climb height, favor that instead of straight line nearest point.

          var diff = DetourCommon.vSub(center, closestPtPoly);

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
    }, {
      key: "getTileAt",
      value: function getTileAt(x, y, layer) {
        // Find tile based on hash.
        var h = NavMesh.computeTileHash(x, y, this.m_tileLutMask);
        var tile = this.m_posLookup[h];

        while (tile != null) {
          if (tile.data.header != null && tile.data.header.x == x && tile.data.header.y == y && tile.data.header.layer == layer) {
            return tile;
          }

          tile = tile.next;
        }

        return null;
      }
    }, {
      key: "getNeighbourTilesAt",
      value: function getNeighbourTilesAt(x, y, side) {
        var nx = x,
            ny = y;

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
    }, {
      key: "getTilesAt",
      value: function getTilesAt(x, y) {
        var tiles = []; // Find tile based on hash.

        var h = NavMesh.computeTileHash(x, y, this.m_tileLutMask);
        var tile = this.m_posLookup[h];

        while (tile != null) {
          if (tile.data.header != null && tile.data.header.x == x && tile.data.header.y == y) {
            tiles.push(tile);
          }

          tile = tile.next;
        }

        return tiles;
      }
    }, {
      key: "getTileRefAt",
      value: function getTileRefAt(x, y, layer) {
        // Find tile based on hash.
        var h = NavMesh.computeTileHash(x, y, this.m_tileLutMask);
        var tile = this.m_posLookup[h];

        while (tile != null) {
          if (tile.data.header != null && tile.data.header.x == x && tile.data.header.y == y && tile.data.header.layer == layer) {
            return this.getTileRef(tile);
          }

          tile = tile.next;
        }

        return 0;
      }
    }, {
      key: "getTileByRef",
      value: function getTileByRef(ref) {
        if (ref == 0) return null;
        var tileIndex = decodePolyIdTile(ref);
        var tileSalt = decodePolyIdSalt(ref);
        if (tileIndex >= this.m_maxTiles) return null;
        var tile = this.m_tiles[tileIndex];
        if (tile.salt != tileSalt) return null;
        return tile;
      }
    }, {
      key: "getTileRef",
      value: function getTileRef(tile) {
        if (tile == null) return 0;
        var it = tile.index;
        return NavMesh.encodePolyId(tile.salt, it, 0);
      }
    }, {
      key: "getOffMeshConnectionPolyEndPoints",
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
      value: function getOffMeshConnectionPolyEndPoints(prevRef, polyRef) {
        if (polyRef == 0) throw new IllegalArgumentException("polyRef = 0"); // Get current polygon

        var saltitip = NavMesh.decodePolyId(polyRef);
        var salt = saltitip[0];
        var it = saltitip[1];
        var ip = saltitip[2];

        if (it >= this.m_maxTiles) {
          throw new IllegalArgumentException("Invalid tile ID > max tiles");
        }

        if (this.m_tiles[it].salt != salt || this.m_tiles[it].data.header == null) {
          throw new IllegalArgumentException("Invalid salt or missing tile header");
        }

        var tile = this.m_tiles[it];

        if (ip >= tile.data.header.polyCount) {
          throw new IllegalArgumentException("Invalid poly ID > poly count");
        }

        var poly = tile.data.polys[ip]; // Make sure that the current poly is indeed off-mesh link.

        if (poly.getType() != Poly.DT_POLYTYPE_OFFMESH_CONNECTION) throw new IllegalArgumentException("Invalid poly type"); // Figure out which way to hand out the vertices.

        var idx0 = 0,
            idx1 = 1; // Find link that points to first vertex.

        for (var i = poly.firstLink; i != NavMesh.DT_NULL_LINK; i = tile.links[i].next) {
          if (tile.links[i].edge == 0) {
            if (tile.links[i].ref != prevRef) {
              idx0 = 1;
              idx1 = 0;
            }

            break;
          }
        }

        var startPos = new Array(3);
        var endPos = new Array(3);
        DetourCommon.vCopy(startPos, tile.data.verts, poly.verts[idx0] * 3);
        DetourCommon.vCopy(endPos, tile.data.verts, poly.verts[idx1] * 3);
        return [startPos, endPos];
      }
    }, {
      key: "getMaxVertsPerPoly",
      value: function getMaxVertsPerPoly() {
        return this.m_maxVertPerPoly;
      }
    }, {
      key: "getTileCount",
      value: function getTileCount() {
        return this.m_tileCount;
      }
    }], [{
      key: "getNavMeshParams",
      value: function getNavMeshParams(data) {
        var params = new NavMeshParams();
        DetourCommon.vCopy(params.orig, data.header.bmin);
        params.tileWidth = data.header.bmax[0] - data.header.bmin[0];
        params.tileHeight = data.header.bmax[2] - data.header.bmin[2];
        params.maxTiles = 1;
        params.maxPolys = data.header.polyCount;
        return params;
      }
    }, {
      key: "getSlabCoord",
      value: function getSlabCoord(verts, va, side) {
        if (side == 0 || side == 4) return verts[va];else if (side == 2 || side == 6) return verts[va + 2];
        return 0;
      }
    }, {
      key: "calcSlabEndPoints",
      value: function calcSlabEndPoints(verts, va, vb, bmin, bmax, side) {
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
    }, {
      key: "computeTileHash",
      value: function computeTileHash(x, y, mask) {
        var h1 = 0x8da6b343; // Large multiplicative constants;

        var h2 = 0xd8163841; // here arbitrarily chosen primes

        var n = h1 * x + h2 * y;
        return n & mask;
      }
    }]);

    return NavMesh;
  }();

  _defineProperty(NavMesh, "DT_SALT_BITS", 16);

  _defineProperty(NavMesh, "DT_TILE_BITS", 28);

  _defineProperty(NavMesh, "DT_POLY_BITS", 20);

  _defineProperty(NavMesh, "DT_EXT_LINK", 0x8000);

  _defineProperty(NavMesh, "DT_NULL_LINK", -1);

  _defineProperty(NavMesh, "DT_OFFMESH_CON_BIDIR", 1);

  _defineProperty(NavMesh, "DT_MAX_AREAS", 64);

  _defineProperty(NavMesh, "DT_RAY_CAST_LIMIT_PROPORTIONS", 50.0);

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
  var RecastVectors$1 = /*#__PURE__*/function () {
    function RecastVectors() {
      _classCallCheck(this, RecastVectors);
    }

    _createClass(RecastVectors, null, [{
      key: "min",
      value: function min(a, b, i) {
        a[0] = Math.min(a[0], b[i + 0]);
        a[1] = Math.min(a[1], b[i + 1]);
        a[2] = Math.min(a[2], b[i + 2]);
      }
    }, {
      key: "max",
      value: function max(a, b, i) {
        a[0] = Math.max(a[0], b[i + 0]);
        a[1] = Math.max(a[1], b[i + 1]);
        a[2] = Math.max(a[2], b[i + 2]);
      }
    }, {
      key: "copy3",
      value: function copy3(out, _in, i) {
        RecastVectors.copy4(out, 0, _in, i);
      }
    }, {
      key: "copy2",
      value: function copy2(out, _in) {
        RecastVectors.copy4(out, 0, _in, 0);
      }
    }, {
      key: "copy4",
      value: function copy4(out, n, _in, m) {
        out[n] = _in[m];
        out[n + 1] = _in[m + 1];
        out[n + 2] = _in[m + 2];
      }
    }, {
      key: "add",
      value: function add(e0, a, verts, i) {
        e0[0] = a[0] + verts[i];
        e0[1] = a[1] + verts[i + 1];
        e0[2] = a[2] + verts[i + 2];
      }
    }, {
      key: "subA",
      value: function subA(e0, verts, i, j) {
        e0[0] = verts[i] - verts[j];
        e0[1] = verts[i + 1] - verts[j + 1];
        e0[2] = verts[i + 2] - verts[j + 2];
      }
    }, {
      key: "subB",
      value: function subB(e0, i, verts, j) {
        e0[0] = i[0] - verts[j];
        e0[1] = i[1] - verts[j + 1];
        e0[2] = i[2] - verts[j + 2];
      }
    }, {
      key: "cross",
      value: function cross(dest, v1, v2) {
        dest[0] = v1[1] * v2[2] - v1[2] * v2[1];
        dest[1] = v1[2] * v2[0] - v1[0] * v2[2];
        dest[2] = v1[0] * v2[1] - v1[1] * v2[0];
      }
    }, {
      key: "normalize",
      value: function normalize(v) {
        var d = 1.0 / Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        v[0] *= d;
        v[1] *= d;
        v[2] *= d;
      }
    }]);

    return RecastVectors;
  }();

  var _temp, _temp2;

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
  var ChunkyTriMesh = /*#__PURE__*/function () {
    _createClass(ChunkyTriMesh, [{
      key: "calcExtends",
      value: function calcExtends(items, imin, imax, bmin, bmax) {
        bmin[0] = items[imin].bmin[0];
        bmin[1] = items[imin].bmin[1];
        bmax[0] = items[imin].bmax[0];
        bmax[1] = items[imin].bmax[1];

        for (var i = imin + 1; i < imax; ++i) {
          var it = items[i];

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
    }, {
      key: "longestAxis",
      value: function longestAxis(x, y) {
        return y > x ? 1 : 0;
      } // https://stackoverflow.com/a/45245772/10047920

    }, {
      key: "partialSort",
      value: function partialSort(arr, start, end, sortFx) {
        var preSorted = arr.slice(0, start),
            postSorted = arr.slice(end);
        var sorted = arr.slice(start, end).sort(sortFx);
        arr.length = 0;
        arr.push.apply(arr, preSorted.concat(sorted).concat(postSorted));
        return arr;
      }
    }, {
      key: "subdivide",
      value: function subdivide(items, imin, imax, trisPerChunk, nodes, inTris) {
        var _this = this;

        var inum = imax - imin;
        var node = new ChunkyTriMesh.ChunkyTriMeshNode();
        nodes.push(node);

        if (inum <= trisPerChunk) {
          // Leaf
          this.calcExtends(items, imin, imax, node.bmin, node.bmax); // Copy triangles.

          node.i = nodes.length;
          node.tris = new Array(inum * 3);
          var dst = 0;

          for (var i = imin; i < imax; ++i) {
            var src = items[i].i * 3;
            node.tris[dst++] = inTris[src];
            node.tris[dst++] = inTris[src + 1];
            node.tris[dst++] = inTris[src + 2];
          }
        } else {
          // Split
          this.calcExtends(items, imin, imax, node.bmin, node.bmax);
          var axis = this.longestAxis(node.bmax[0] - node.bmin[0], node.bmax[1] - node.bmin[1]);

          if (axis == 0) {
            //Arrays.sort(items, imin, imax, new CompareItemX());
            this.partialSort(items, imin, imax, function (a, b) {
              return _this.CompareItemX(a, b);
            }); // let subarray = items.slice(imin, imax);
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
            this.partialSort(items, imin, imax, function (a, b) {
              return _this.CompareItemY(a, b);
            }); // let subarray = items.slice(imin, imax);
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

          var isplit = Math.floor(imin + Math.floor(inum / 2));
          var s = "";

          for (var _i = 0; _i < items.length; _i++) {
            var item = items[_i];
            s += "" + _i + " " + item.bmin[0] + " " + item.bmin[1] + "\n";
            s += "" + _i + " " + item.bmax[0] + " " + item.bmax[1] + "\n";
          } // fs.writeFileSync("items_" + imin + "_" + imax + ".txt", s);
          // console.log("done")
          // Left
          // console.log("Before left " + imin + " " + isplit);
          // console.log(items[279].bmin[1])


          this.subdivide(items, imin, isplit, trisPerChunk, nodes, inTris); // console.log("Done left " + imin + " " + isplit);
          // console.log(items[279].bmin[1])
          // Right
          // console.log("Before right " + isplit + " " + imax);
          // console.log(items[279].bmin[1])

          this.subdivide(items, isplit, imax, trisPerChunk, nodes, inTris); // console.log("Done right " + isplit + " " + imax);
          // console.log(items[279].bmin[1])
          // Negative index means escape.

          node.i = -nodes.length;
        }
      }
    }]);

    function ChunkyTriMesh(verts, tris, ntris, trisPerChunk) {
      _classCallCheck(this, ChunkyTriMesh);

      _defineProperty(this, "CompareItemX", function compare(a, b) {
        if (a.bmin[0] < b.bmin[0]) {
          return -1;
        }

        if (a.bmin[0] > b.bmin[0]) {
          return 1;
        }

        return 0;
      });

      _defineProperty(this, "CompareItemY", function compare(a, b) {
        if (a.bmin[1] < b.bmin[1]) {
          return -1;
        }

        if (a.bmin[1] > b.bmin[1]) {
          return 1;
        }

        return 0;
      });

      _defineProperty(this, "nodes", []);

      _defineProperty(this, "ntris", void 0);

      _defineProperty(this, "maxTrisPerChunk", void 0);
      this.nodes = [];
      this.ntris = ntris; // Build tree

      var items = []; //new Array(ntris);

      for (var i = 0; i < ntris; i++) {
        var t = i * 3;
        var it = items[i] = new ChunkyTriMesh.BoundsItem();
        it.i = i; // Calc triangle XZ bounds.

        it.bmin[0] = it.bmax[0] = verts[tris[t] * 3 + 0];
        it.bmin[1] = it.bmax[1] = verts[tris[t] * 3 + 2];

        for (var j = 1; j < 3; ++j) {
          var v = tris[t + j] * 3;

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

      this.subdivide(items, 0, ntris, trisPerChunk, this.nodes, tris); // Calc max tris per node.

      this.maxTrisPerChunk = 0;

      var _iterator = _createForOfIteratorHelper(this.nodes),
          _step;

      try {
        for (_iterator.s(); !(_step = _iterator.n()).done;) {
          var _node = _step.value;
          var isLeaf = _node.i >= 0;

          if (!isLeaf) {
            continue;
          }

          if (_node.tris.length / 3 > this.maxTrisPerChunk) {
            this.maxTrisPerChunk = _node.tris.length / 3;
          }
        }
      } catch (err) {
        _iterator.e(err);
      } finally {
        _iterator.f();
      }
    }

    _createClass(ChunkyTriMesh, [{
      key: "checkOverlapRect",
      value: function checkOverlapRect(amin, amax, bmin, bmax) {
        var overlap = true;
        overlap = amin[0] > bmax[0] || amax[0] < bmin[0] ? false : overlap;
        overlap = amin[1] > bmax[1] || amax[1] < bmin[1] ? false : overlap;
        return overlap;
      }
    }, {
      key: "getChunksOverlappingRect",
      value: function getChunksOverlappingRect(bmin, bmax) {
        // Traverse tree
        ids = [];
        var i = 0;

        while (i < nodes.length) {
          node = nodes[i];
          var overlap = checkOverlapRect(bmin, bmax, node.bmin, node.bmax);
          var isLeafNode = node.i >= 0;

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
    }]);

    return ChunkyTriMesh;
  }();

  _defineProperty(ChunkyTriMesh, "BoundsItem", (_temp = function BoundsItem() {
    _classCallCheck(this, BoundsItem);

    _defineProperty(this, "bmin", new Array(2));

    _defineProperty(this, "bmax", new Array(2));

    _defineProperty(this, "i", void 0);
  }, _temp));

  _defineProperty(ChunkyTriMesh, "ChunkyTriMeshNode", (_temp2 = function ChunkyTriMeshNode() {
    _classCallCheck(this, ChunkyTriMeshNode);

    _defineProperty(this, "bmin", new Array(2));

    _defineProperty(this, "bmax", new Array(2));

    _defineProperty(this, "i", void 0);

    _defineProperty(this, "tris", []);
  }, _temp2));

  var TriMesh = /*#__PURE__*/function () {
    function TriMesh(vertices, faces) {
      _classCallCheck(this, TriMesh);

      _defineProperty(this, "vertices", void 0);

      _defineProperty(this, "faces", void 0);

      _defineProperty(this, "chunkyTriMesh", void 0);

      this.vertices = vertices;
      this.faces = faces;
      this.chunkyTriMesh = new ChunkyTriMesh(vertices, faces, faces.length / 3, 256);
    }

    _createClass(TriMesh, [{
      key: "getTris",
      value: function getTris() {
        return this.faces;
      }
    }, {
      key: "getVerts",
      value: function getVerts() {
        return this.vertices;
      }
    }, {
      key: "getChunksOverlappingRect",
      value: function getChunksOverlappingRect(bmin, bmax) {
        return chunkyTriMesh.getChunksOverlappingRect(bmin, bmax);
      }
    }]);

    return TriMesh;
  }();

  var SimpleInputGeomProvider
  /*extends InputGeomProvider*/
  = /*#__PURE__*/function () {
    _createClass(SimpleInputGeomProvider, null, [{
      key: "mapFaces",
      value: function mapFaces(meshFaces) {
        var faces = new Array(meshFaces.length);

        for (var i = 0; i < faces.length; i++) {
          faces[i] = meshFaces[i];
        }

        return faces;
      }
    }, {
      key: "mapVertices",
      value: function mapVertices(vertexPositions) {
        var vertices = new Array(vertexPositions.length);

        for (var i = 0; i < vertices.length; i++) {
          vertices[i] = vertexPositions[i];
        }

        return vertices;
      }
    }]);

    function SimpleInputGeomProvider(vertices, faces) {
      _classCallCheck(this, SimpleInputGeomProvider);

      _defineProperty(this, "vertices", void 0);

      _defineProperty(this, "faces", void 0);

      _defineProperty(this, "bmin", void 0);

      _defineProperty(this, "bmax", void 0);

      _defineProperty(this, "volumes", []);

      this.vertices = vertices;
      this.faces = faces;
      this.bmin = new Array(3);
      this.bmax = new Array(3);
      RecastVectors$1.copy3(this.bmin, vertices, 0);
      RecastVectors$1.copy3(this.bmax, vertices, 0);

      for (var i = 1; i < vertices.length / 3; i++) {
        RecastVectors$1.min(this.bmin, vertices, i * 3);
        RecastVectors$1.max(this.bmax, vertices, i * 3);
      }
    }

    _createClass(SimpleInputGeomProvider, [{
      key: "getMeshBoundsMin",
      value: function getMeshBoundsMin() {
        return this.bmin;
      }
    }, {
      key: "getMeshBoundsMax",
      value: function getMeshBoundsMax() {
        return this.bmax;
      }
    }, {
      key: "getConvexVolumes",
      value: function getConvexVolumes() {
        return [];
      }
    }, {
      key: "addConvexVolume",
      value: function addConvexVolume(verts, minh, maxh, areaMod) {
        var vol = new ConvexVolume();
        vol.hmin = minh;
        vol.hmax = maxh;
        vol.verts = verts;
        vol.areaMod = areaMod;
        this.volumes.push(vol);
      }
    }, {
      key: "meshes",
      value: function meshes() {
        // return Collections.singPolyonList(new TriMesh(vertices, faces));
        return [new TriMesh(this.vertices, this.faces)];
      }
    }]);

    return SimpleInputGeomProvider;
  }();

  _defineProperty(SimpleInputGeomProvider, "fromIndeces", function (vertexPositions, meshFaces) {
    return new SimpleInputGeomProvider(this.mapVertices(vertexPositions), this.mapFaces(meshFaces));
  });

  var ObjImporterContext = function ObjImporterContext() {
    _classCallCheck(this, ObjImporterContext);

    _defineProperty(this, "vertexPositions", []);

    _defineProperty(this, "meshFaces", []);
  };

  var ObjImporter = /*#__PURE__*/function () {
    function ObjImporter() {
      _classCallCheck(this, ObjImporter);
    }

    _createClass(ObjImporter, [{
      key: "load",
      // OBJImporterContext = 
      value: function load(is) {
        var context = new ObjImporterContext();

        try {
          var slurp = is;
          var lines = slurp.split(/\r?\n/);

          for (var i = 0; i < lines.length; i++) {
            var line = lines[i];
            this.readLine(line, context);
          } // reader = new BufferedReader(new InputStreamReader(is));
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
    }, {
      key: "readLine",
      value: function readLine(line, context) {
        if (line.startsWith("v")) {
          this.readVertex(line, context);
        } else if (line.startsWith("f")) {
          this.readFace(line, context);
        }
      }
    }, {
      key: "readVertex",
      value: function readVertex(line, context) {
        if (line.startsWith("v ")) {
          var vert = this.readVector3f(line);

          var _iterator = _createForOfIteratorHelper(vert),
              _step;

          try {
            for (_iterator.s(); !(_step = _iterator.n()).done;) {
              var vp = _step.value;
              context.vertexPositions.push(vp);
            }
          } catch (err) {
            _iterator.e(err);
          } finally {
            _iterator.f();
          }
        }
      }
    }, {
      key: "readVector3f",
      value: function readVector3f(line) {
        var v = line.split(/\s+/);

        if (v.length < 4) {
          throw new RuntimeException("Invalid vector, expected 3 coordinates, found " + (v.length - 1));
        }

        return [parseFloat(v[1]), parseFloat(v[2]), parseFloat(v[3])];
      }
    }, {
      key: "readFace",
      value: function readFace(line, context) {
        var v = line.split(/\s+/);

        if (v.length < 4) {
          throw new RuntimeException("Invalid number of face vertices: 3 coordinates expected, found " + v.length);
        }

        for (var j = 0; j < v.length - 3; j++) {
          context.meshFaces.push(this.readFaceVertex(v[1], context));

          for (var i = 0; i < 2; i++) {
            context.meshFaces.push(this.readFaceVertex(v[2 + j + i], context));
          }
        }
      }
    }, {
      key: "readFaceVertex",
      value: function readFaceVertex(face, context) {
        var v = face.split(/\//);
        return this.getIndex(parseInt(v[0]), context.vertexPositions.length);
      }
    }, {
      key: "getIndex",
      value: function getIndex(posi, size) {
        if (posi > 0) {
          posi--;
        } else if (posi < 0) {
          posi = size + posi;
        } else {
          throw new RuntimeException("0 vertex index");
        }

        return posi;
      }
    }]);

    return ObjImporter;
  }();

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
  var RecastConstants = function RecastConstants() {
    _classCallCheck(this, RecastConstants);
  };

  _defineProperty(RecastConstants, "RC_NULL_AREA", 0);

  _defineProperty(RecastConstants, "RC_NOT_CONNECTED", 0x3f);

  _defineProperty(RecastConstants, "RC_WALKABLE_AREA", 63);

  _defineProperty(RecastConstants, "RC_SPAN_HEIGHT_BITS", 13);

  _defineProperty(RecastConstants, "RC_SPAN_MAX_HEIGHT", (1 << RecastConstants.RC_SPAN_HEIGHT_BITS) - 1);

  _defineProperty(RecastConstants, "RC_BORDER_REG", 0x8000);

  _defineProperty(RecastConstants, "RC_MULTIPLE_REGS", 0);

  _defineProperty(RecastConstants, "RC_BORDER_VERTEX", 0x10000);

  _defineProperty(RecastConstants, "RC_AREA_BORDER", 0x20000);

  _defineProperty(RecastConstants, "RC_CONTOUR_REG_MASK", 0xffff);

  _defineProperty(RecastConstants, "RC_MESH_NULL_IDX", 0xffff);

  _defineProperty(RecastConstants, "RC_CONTOUR_TESS_WALL_EDGES", 0x01);

  _defineProperty(RecastConstants, "RC_CONTOUR_TESS_AREA_EDGES", 0x02);

  _defineProperty(RecastConstants, "WATERSHED", 10);

  _defineProperty(RecastConstants, "MONOTONE", 20);

  _defineProperty(RecastConstants, "LAYERS", 30);

  _defineProperty(RecastConstants, "RC_LOG_WARNING", 1);

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
  var RecastConfig =
  /** The width/height size of tile's on the xz-plane. [Limit: >= 0] [Units: vx] **/

  /** The xz-plane cell size to use for fields. [Limit: > 0] [Units: wu] **/

  /** The y-axis cell size to use for fields. [Limit: > 0] [Units: wu] **/

  /** The maximum slope that is considered walkable. [Limits: 0 <= value < 90] [Units: Degrees] **/

  /**
   * Minimum floor to 'ceiling' height that will still allow the floor area to be considered walkable. [Limit: >= 3]
   * [Units: vx]
   **/

  /** Maximum ledge height that is considered to still be traversable. [Limit: >=0] [Units: vx] **/

  /**
   * The distance to erode/shrink the walkable area of the heightfield away from obstructions. [Limit: >=0] [Units:
   * vx]
   **/

  /** The maximum allowed length for contour edges aPoly the border of the mesh. [Limit: >=0] [Units: vx] **/

  /**
   * The maximum distance a simplfied contour's border edges should deviate the original raw contour. [Limit: >=0]
   * [Units: vx]
   **/

  /** The minimum number of cells allowed to form isolated island areas. [Limit: >=0] [Units: vx] **/

  /**
   * Any regions with a span count smaller than this value will, if possible, be merged with larger regions. [Limit:
   * >=0] [Units: vx]
   **/

  /**
   * The maximum number of vertices allowed for polygons generated during the contour to polygon conversion process.
   * [Limit: >= 3]
   **/

  /**
   * Sets the sampling distance to use when generating the detail mesh. (For height detail only.) [Limits: 0 or >=
   * 0.9] [Units: wu]
   **/

  /**
   * The maximum distance the detail mesh surface should deviate from heightfield data. (For height detail only.)
   * [Limit: >=0] [Units: wu]
   **/
  function RecastConfig(partitionType, cellSize, cellHeight, agentHeight, agentRadius, agentMaxClimb, agentMaxSlope, regionMinSize, regionMergeSize, edgeMaxLen, edgeMaxError, vertsPerPoly, detailSampleDist, detailSampleMaxError, tileSize, walkableAreaMod) {
    _classCallCheck(this, RecastConfig);

    _defineProperty(this, "partitionType", null);

    _defineProperty(this, "tileSize", void 0);

    _defineProperty(this, "cs", void 0);

    _defineProperty(this, "ch", void 0);

    _defineProperty(this, "walkableSlopeAngle", void 0);

    _defineProperty(this, "walkableHeight", void 0);

    _defineProperty(this, "walkableClimb", void 0);

    _defineProperty(this, "walkableRadius", void 0);

    _defineProperty(this, "maxEdgeLen", void 0);

    _defineProperty(this, "maxSimplificationError", void 0);

    _defineProperty(this, "minRegionArea", void 0);

    _defineProperty(this, "mergeRegionArea", void 0);

    _defineProperty(this, "maxVertsPerPoly", 0);

    _defineProperty(this, "detailSampleDist", void 0);

    _defineProperty(this, "detailSampleMaxError", void 0);

    _defineProperty(this, "walkableAreaMod", void 0);

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
  };

  var AreaModification = /*#__PURE__*/function () {
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
    function AreaModification(value) {
      var mask = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : AreaModification.RC_AREA_FLAGS_MASK;

      _classCallCheck(this, AreaModification);

      _defineProperty(this, "value", void 0);

      _defineProperty(this, "mask", void 0);

      this.value = value;
      this.mask = mask;
    }

    _createClass(AreaModification, [{
      key: "getMaskedValue",
      value: function getMaskedValue() {
        return this.value & this.mask;
      }
    }, {
      key: "apply",
      value: function apply(area) {
        return this.value & this.mask | area & ~this.mask;
      }
    }]);

    return AreaModification;
  }();

  _defineProperty(AreaModification, "RC_AREA_FLAGS_MASK", 0x3F);

  _defineProperty(AreaModification, "clone", function (other) {
    return new AreaModification(other.value, other.mask); // this.value = other.value;
    // this.mask = other.mask;
  });

  var SampleAreaModifications = function SampleAreaModifications() {
    _classCallCheck(this, SampleAreaModifications);
  };

  _defineProperty(SampleAreaModifications, "SAMPLE_POLYAREA_TYPE_MASK", 0x07);

  _defineProperty(SampleAreaModifications, "SAMPLE_POLYAREA_TYPE_GROUND", 0x1);

  _defineProperty(SampleAreaModifications, "SAMPLE_POLYAREA_TYPE_WATER", 0x2);

  _defineProperty(SampleAreaModifications, "SAMPLE_POLYAREA_TYPE_ROAD", 0x3);

  _defineProperty(SampleAreaModifications, "SAMPLE_POLYAREA_TYPE_DOOR", 0x4);

  _defineProperty(SampleAreaModifications, "SAMPLE_POLYAREA_TYPE_GRASS", 0x5);

  _defineProperty(SampleAreaModifications, "SAMPLE_POLYAREA_TYPE_JUMP", 0x6);

  _defineProperty(SampleAreaModifications, "SAMPLE_AREAMOD_GROUND", new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_GROUND, SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK));

  _defineProperty(SampleAreaModifications, "SAMPLE_AREAMOD_WATER", new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_WATER, SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK));

  _defineProperty(SampleAreaModifications, "SAMPLE_AREAMOD_ROAD", new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_ROAD, SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK));

  _defineProperty(SampleAreaModifications, "SAMPLE_AREAMOD_GRASS", new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_GRASS, SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK));

  _defineProperty(SampleAreaModifications, "SAMPLE_AREAMOD_DOOR", new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_DOOR, SampleAreaModifications.SAMPLE_POLYAREA_TYPE_DOOR));

  _defineProperty(SampleAreaModifications, "SAMPLE_AREAMOD_JUMP", new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_JUMP, SampleAreaModifications.SAMPLE_POLYAREA_TYPE_JUMP));

  _defineProperty(SampleAreaModifications, "SAMPLE_POLYFLAGS_WALK", 0x01);

  _defineProperty(SampleAreaModifications, "SAMPLE_POLYFLAGS_SWIM", 0x02);

  _defineProperty(SampleAreaModifications, "SAMPLE_POLYFLAGS_DOOR", 0x04);

  _defineProperty(SampleAreaModifications, "SAMPLE_POLYFLAGS_JUMP", 0x08);

  _defineProperty(SampleAreaModifications, "SAMPLE_POLYFLAGS_DISABLED", 0x10);

  _defineProperty(SampleAreaModifications, "SAMPLE_POLYFLAGS_ALL", 0xfff);

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
  var CompactHeightfield = function CompactHeightfield() {
    _classCallCheck(this, CompactHeightfield);

    _defineProperty(this, "width", 0);

    _defineProperty(this, "height", 0);

    _defineProperty(this, "spanCount", 0);

    _defineProperty(this, "walkableHeight", 0);

    _defineProperty(this, "walkableClimb", 0);

    _defineProperty(this, "borderSize", 0);

    _defineProperty(this, "maxDistance", 0);

    _defineProperty(this, "maxRegions", void 0);

    _defineProperty(this, "bmin", new Array(3));

    _defineProperty(this, "bmax", new Array(3));

    _defineProperty(this, "cs", 0);

    _defineProperty(this, "ch", 0);

    _defineProperty(this, "cells", []);

    _defineProperty(this, "spans", []);

    _defineProperty(this, "dist", []);

    _defineProperty(this, "areas", []);
  };

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
  var CompactCell = function CompactCell() {
    _classCallCheck(this, CompactCell);

    _defineProperty(this, "index", 0);

    _defineProperty(this, "count", 0);
  };

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
  var CompactSpan = function CompactSpan() {
    _classCallCheck(this, CompactSpan);

    _defineProperty(this, "y", 0);

    _defineProperty(this, "reg", 0);

    _defineProperty(this, "con", 0);

    _defineProperty(this, "h", 0);
  };

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
  var RecastCommon = function RecastCommon() {
    _classCallCheck(this, RecastCommon);
  };

  _defineProperty(RecastCommon, "GetCon", function (s, dir) {
    var shift = dir * 6;
    return s.con >> shift & 0x3f;
  });

  _defineProperty(RecastCommon, "GetDirOffsetX", function (dir) {
    var offset = [-1, 0, 1, 0];
    return offset[dir & 0x03];
  });

  _defineProperty(RecastCommon, "GetDirOffsetY", function (dir) {
    var offset = [0, 1, 0, -1];
    return offset[dir & 0x03];
  });

  _defineProperty(RecastCommon, "rcGetDirForOffset", function (x, y) {
    var dirs = [3, 0, -1, 2, 1];
    return dirs[(y + 1 << 1) + x];
  });

  _defineProperty(RecastCommon, "SetCon", function (s, dir, i) {
    var shift = dir * 6;
    var con = s.con;
    s.con = con & ~(0x3f << shift) | (i & 0x3f) << shift;
  });

  _defineProperty(RecastCommon, "clamp", function (v, min, max) {
    return Math.max(Math.min(max, v), min);
  });

  var Recast = /*#__PURE__*/function () {
    function Recast() {
      _classCallCheck(this, Recast);
    }

    _createClass(Recast, [{
      key: "calcBounds",
      value: function calcBounds(verts, nv, bmin, bmax) {
        for (var i = 0; i < 3; i++) {
          bmin[i] = verts[i];
          bmax[i] = verts[i];
        }

        for (var _i = 1; _i < nv; ++_i) {
          for (var j = 0; j < 3; j++) {
            bmin[j] = Math.min(bmin[j], verts[_i * 3 + j]);
            bmax[j] = Math.max(bmax[j], verts[_i * 3 + j]);
          }
        } // Calculate bounding box.

      }
    }, {
      key: "clearUnwalkableTriangles",
      /// @par
      ///
      /// Only sets the area id's for the unwalkable triangles. Does not alter the
      /// area id's for walkable triangles.
      ///
      /// See the #rcConfig documentation for more information on the configuration parameters.
      ///
      /// @see rcHeightfield, rcClearUnwalkableTriangles, rcRasterizeTriangles
      value: function clearUnwalkableTriangles(ctx, walkableSlopeAngle, verts, nv, tris, nt, areas) {
        walkableThr = Math.cos(walkableSlopeAngle / 180.0 * Math.PI);
        norm = new Array(3);

        for (var i = 0; i < nt; ++i) {
          var tri = i * 3;
          calcTriNormal(verts, tris[tri], tris[tri + 1], tris[tri + 2], norm); // Check if the face is walkable.

          if (norm[1] <= walkableThr) areas[i] = RecastConstants.RC_NULL_AREA;
        }
      }
    }], [{
      key: "calcGridSize",
      value: function calcGridSize(bmin, bmax, cs) {
        return [Math.floor((bmax[0] - bmin[0]) / cs + 0.5), Math.floor((bmax[2] - bmin[2]) / cs + 0.5)];
      }
    }, {
      key: "calcTileCount",
      value: function calcTileCount(bmin, bmax, cs, tileSize) {
        var gwh = Recast.calcGridSize(bmin, bmax, cs);
        var gw = gwh[0];
        var gh = gwh[1];
        var ts = tileSize;
        var tw = (gw + ts - 1) / ts;
        var th = (gh + ts - 1) / ts;
        return [tw, th];
      } /// @par
      ///
      /// Modifies the area id of all triangles with a slope below the specified value.
      ///
      /// See the #rcConfig documentation for more information on the configuration parameters.
      ///
      /// @see rcHeightfield, rcClearUnwalkableTriangles, rcRasterizeTriangles

    }, {
      key: "markWalkableTriangles",
      value: function markWalkableTriangles(ctx, walkableSlopeAngle, verts, tris, nt, areaMod) {
        var areas = new Array(nt).fill(0);
        var walkableThr = Math.cos(walkableSlopeAngle / 180.0 * Math.PI);
        var norm = new Array(3).fill(0);

        for (var i = 0; i < nt; ++i) {
          var tri = i * 3;
          Recast.calcTriNormal(verts, tris[tri], tris[tri + 1], tris[tri + 2], norm); // Check if the face is walkable.

          if (norm[1] > walkableThr) areas[i] = areaMod.apply(areas[i]);
        }

        return areas;
      }
    }, {
      key: "calcTriNormal",
      value: function calcTriNormal(verts, v0, v1, v2, norm) {
        var e0 = new Array(3);
        var e1 = new Array(3);
        RecastVectors$1.subA(e0, verts, v1 * 3, v0 * 3);
        RecastVectors$1.subA(e1, verts, v2 * 3, v0 * 3);
        RecastVectors$1.cross(norm, e0, e1);
        RecastVectors$1.normalize(norm);
      }
    }, {
      key: "getHeightFieldSpanCount",
      value: function getHeightFieldSpanCount(ctx, hf) {
        var w = hf.width;
        var h = hf.height;
        var spanCount = 0;

        for (var y = 0; y < h; ++y) {
          for (var x = 0; x < w; ++x) {
            for (var s = hf.spans[x + y * w]; s != null; s = s.next) {
              if (s.area != RecastConstants.RC_NULL_AREA) spanCount++;
            }
          }
        }

        return spanCount;
      } /// @par
      ///
      /// This is just the beginning of the process of fully building a compact heightfield.
      /// Various filters may be applied, then the distance field and regions built.
      /// E.g: #rcBuildDistanceField and #rcBuildRegions
      ///
      /// See the #rcConfig documentation for more information on the configuration parameters.
      ///
      /// @see rcAllocCompactHeightfield, rcHeightfield, rcCompactHeightfield, rcConfig

    }, {
      key: "buildCompactHeightfield",
      value: function buildCompactHeightfield(ctx, walkableHeight, walkableClimb, hf) {
        ctx.startTimer("BUILD_COMPACTHEIGHTFIELD");
        var chf = new CompactHeightfield();
        var w = hf.width;
        var h = hf.height;
        var spanCount = this.getHeightFieldSpanCount(ctx, hf); // Fill in header.

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
        chf.ch = hf.ch; // let bigSize = w*h;

        chf.cells = new Array(Math.floor(w) * Math.floor(h));
        chf.spans = new Array(spanCount);
        chf.areas = new Array(spanCount);
        var MAX_HEIGHT = 0xffff;

        for (var i = 0; i < chf.cells.length; i++) {
          chf.cells[i] = new CompactCell();
        }

        for (var _i2 = 0; _i2 < chf.spans.length; _i2++) {
          chf.spans[_i2] = new CompactSpan();
        } // Fill in cells and spans.


        var idx = 0;

        for (var _y = 0; _y < h; ++_y) {
          for (var _x = 0; _x < w; ++_x) {
            var s = hf.spans[_x + _y * w]; // If there are no spans at this cell, just leave the data to index=0, count=0.

            if (s == null) continue;
            var c = chf.cells[_x + _y * w];
            c.index = idx;
            c.count = 0;

            while (s != null) {
              if (s.area != RecastConstants.RC_NULL_AREA) {
                var bot = s.smax;
                var top = s.next != null ? Math.floor(s.next.smin) : MAX_HEIGHT;
                chf.spans[idx].y = RecastCommon.clamp(bot, 0, 0xffff);
                chf.spans[idx].h = RecastCommon.clamp(top - bot, 0, 0xff);
                chf.areas[idx] = s.area;
                idx++; // if(idx == 450240)
                // 	console.log("here")

                c.count++;
              }

              s = s.next;
            }
          }
        } // Find neighbour connections.


        var MAX_LAYERS = RecastConstants.RC_NOT_CONNECTED - 1;
        var tooHighNeighbour = 0;

        for (var y = 0; y < h; ++y) {
          for (var x = 0; x < w; ++x) {
            var _c = chf.cells[x + y * w];

            for (var _i3 = _c.index, ni = _c.index + _c.count; _i3 < ni; ++_i3) {
              // if(i == 450240)
              // 	console.log("stop")
              var _s = chf.spans[_i3];

              for (var dir = 0; dir < 4; ++dir) {
                RecastCommon.SetCon(_s, dir, RecastConstants.RC_NOT_CONNECTED);
                var nx = x + RecastCommon.GetDirOffsetX(dir);
                var ny = y + RecastCommon.GetDirOffsetY(dir); // First check that the neighbour cell is in bounds.

                if (nx < 0 || ny < 0 || nx >= w || ny >= h) continue; // Iterate over all neighbour spans and check if any of the is
                // accessible from current cell.

                var nc = chf.cells[nx + ny * w];

                for (var k = nc.index, nk = nc.index + nc.count; k < nk; ++k) {
                  var ns = chf.spans[k];

                  var _bot = Math.max(_s.y, ns.y);

                  var _top = Math.min(_s.y + _s.h, ns.y + ns.h); // Check that the gap between the spans is walkable,
                  // and that the climb height between the gaps is not too high.


                  if (_top - _bot >= walkableHeight && Math.abs(ns.y - _s.y) <= walkableClimb) {
                    // Mark direction as walkable.
                    var lidx = k - nc.index;

                    if (lidx < 0 || lidx > MAX_LAYERS) {
                      tooHighNeighbour = Math.max(tooHighNeighbour, lidx);
                      continue;
                    }

                    RecastCommon.SetCon(_s, dir, lidx);
                    break;
                  }
                }
              }
            }
          }
        }

        if (tooHighNeighbour > MAX_LAYERS) {
          throw new RuntimeException("rcBuildCompactHeightfield: Heightfield has too many layers " + tooHighNeighbour + " (max: " + MAX_LAYERS + ")");
        }

        ctx.stopTimer("BUILD_COMPACTHEIGHTFIELD");
        return chf;
      }
    }]);

    return Recast;
  }();

  var RecastBuilderConfig$1 =
  /** The width of the field aPoly the x-axis. [Limit: >= 0] [Units: vx] **/

  /** The height of the field aPoly the z-axis. [Limit: >= 0] [Units: vx] **/

  /** The minimum bounds of the field's AABB. [(x, y, z)] [Units: wu] **/

  /** The maximum bounds of the field's AABB. [(x, y, z)] [Units: wu] **/

  /** The size of the non-navigable border around the heightfield. [Limit: >=0] [Units: vx] **/

  /** Set to true for tiled build **/

  /** Set to false to disable building detailed mesh **/
  // RecastBuilderConfig(RecastConfig cfg, bmin, bmax) {
  //         this(cfg, bmin, bmax, 0, 0, false);
  //     }
  // RecastBuilderConfig(RecastConfig cfg, bmin, bmax,  tx,  ty,  tiled) {
  //         this(cfg, bmin, bmax, tx, ty, tiled, true);
  //     }
  function RecastBuilderConfig(cfg, bmin, bmax) {
    var tx = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : 0;
    var ty = arguments.length > 4 && arguments[4] !== undefined ? arguments[4] : 0;
    var tiled = arguments.length > 5 && arguments[5] !== undefined ? arguments[5] : false;
    var buildMeshDetail = arguments.length > 6 && arguments[6] !== undefined ? arguments[6] : true;

    _classCallCheck(this, RecastBuilderConfig);

    _defineProperty(this, "cfg", null);

    _defineProperty(this, "width", 0);

    _defineProperty(this, "height", 0);

    _defineProperty(this, "bmin", new Array(3));

    _defineProperty(this, "bmax", new Array(3));

    _defineProperty(this, "borderSize", 0);

    _defineProperty(this, "tiled", 0);

    _defineProperty(this, "buildMeshDetail", false);

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
      this.bmax[2] = this.bmin[2] + ts; // Expand the heighfield bounding box by border size to find the extents of geometry we need to build this
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
      var wh = Recast.calcGridSize(this.bmin, this.bmax, cfg.cs);
      this.width = wh[0];
      this.height = wh[1];
      this.borderSize = 0;
    } //ADDED:


    this.width = Math.floor(this.width);
    this.height = Math.floor(this.height);
  };

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
  var Context = /*#__PURE__*/function () {
    function Context() {
      _classCallCheck(this, Context);
    }

    _createClass(Context, [{
      key: "startTimer",
      value: function startTimer(string) {// TODO Auto-generated method stub
      }
    }, {
      key: "stopTimer",
      value: function stopTimer(string) {// TODO Auto-generated method stub
      }
    }, {
      key: "warn",
      value: function warn(string) {
        console.log(string);
      }
    }]);

    return Context;
  }();

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
  var Heightfield =
  /** The width of the heightfield. (APoly the x-axis in cell units.) */

  /** The height of the heightfield. (APoly the z-axis in cell units.) */

  /** The minimum bounds in world space. [(x, y, z)] */

  /** The maximum bounds in world space. [(x, y, z)] */

  /** The size of each cell. (On the xz-plane.) */

  /** The height of each cell. (The minimum increment aPoly the y-axis.) */

  /** Heightfield of spans (width*height). */
  function Heightfield(width, height, bmin, bmax, cs, ch) {
    _classCallCheck(this, Heightfield);

    _defineProperty(this, "width", 0);

    _defineProperty(this, "height", 0);

    _defineProperty(this, "bmin", void 0);

    _defineProperty(this, "bmax", void 0);

    _defineProperty(this, "cs", void 0);

    _defineProperty(this, "ch", void 0);

    _defineProperty(this, "spans", []);

    this.width = width;
    this.height = height;
    this.bmin = bmin;
    this.bmax = bmax;
    this.cs = cs;
    this.ch = ch;
    this.spans = new Array(width * height);
  };

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
  var Span = function Span() {
    _classCallCheck(this, Span);

    _defineProperty(this, "smin", 0);

    _defineProperty(this, "smax", 0);

    _defineProperty(this, "area", 0);

    _defineProperty(this, "next", null);
  };

  var RecastRasterization = /*#__PURE__*/function () {
    function RecastRasterization() {
      _classCallCheck(this, RecastRasterization);
    }

    _createClass(RecastRasterization, null, [{
      key: "overlapBounds",
      value: function overlapBounds(amin, amax, bmin, bmax) {
        var overlap = true;
        overlap = amin[0] > bmax[0] || amax[0] < bmin[0] ? false : overlap;
        overlap = amin[1] > bmax[1] || amax[1] < bmin[1] ? false : overlap;
        overlap = amin[2] > bmax[2] || amax[2] < bmin[2] ? false : overlap;
        return overlap;
      }
      /**
       * The span addition can be set to favor flags. If the span is merged to another span and the new 'smax' is
       * within 'flagMergeThr' units from the existing span, the span flags are merged.
       * 
       * @see Heightfield, Span.
       */

    }, {
      key: "addSpan",
      value: function addSpan(hf, x, y, smin, smax, area, flagMergeThr) {
        var idx = x + y * hf.width;
        var s = new Span();
        s.smin = smin;
        s.smax = smax;
        s.area = area;
        s.next = null; // Empty cell, add the first span.

        if (hf.spans[idx] == null) {
          hf.spans[idx] = s;
          return;
        }

        var prev = null;
        var cur = hf.spans[idx]; // Insert and merge spans.

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
            if (cur.smin < s.smin) s.smin = cur.smin;
            if (cur.smax > s.smax) s.smax = cur.smax; // Merge flags.

            if (Math.abs(s.smax - cur.smax) <= flagMergeThr) s.area = Math.max(s.area, cur.area); // Remove current span.

            var next = cur.next;
            if (prev != null) prev.next = next;else hf.spans[idx] = next;
            cur = next;
          }
        } // Insert new span.


        if (prev != null) {
          s.next = prev.next;
          prev.next = s;
        } else {
          s.next = hf.spans[idx];
          hf.spans[idx] = s;
        }
      } //divides a convex polygons into two convex polygons on both sides of a line

    }, {
      key: "dividePoly",
      value: function dividePoly(buf, _in, nin, out1, out2, x, axis) {
        var d = new Array(12);

        for (var i = 0; i < nin; ++i) {
          d[i] = x - buf[_in + i * 3 + axis];
        }

        var m = 0,
            n = 0;

        for (var _i = 0, j = nin - 1; _i < nin; j = _i, ++_i) {
          var ina = d[j] >= 0;
          var inb = d[_i] >= 0;

          if (ina != inb) {
            var s = d[j] / (d[j] - d[_i]);
            buf[out1 + m * 3 + 0] = buf[_in + j * 3 + 0] + (buf[_in + _i * 3 + 0] - buf[_in + j * 3 + 0]) * s;
            buf[out1 + m * 3 + 1] = buf[_in + j * 3 + 1] + (buf[_in + _i * 3 + 1] - buf[_in + j * 3 + 1]) * s;
            buf[out1 + m * 3 + 2] = buf[_in + j * 3 + 2] + (buf[_in + _i * 3 + 2] - buf[_in + j * 3 + 2]) * s;
            RecastVectors$1.copy4(buf, out2 + n * 3, buf, out1 + m * 3);
            m++;
            n++; // add the i'th poPoly to the right polygon. Do NOT add points that are on the dividing line
            // since these were already added above

            if (d[_i] > 0) {
              RecastVectors$1.copy4(buf, out1 + m * 3, buf, _in + _i * 3);
              m++;
            } else if (d[_i] < 0) {
              RecastVectors$1.copy4(buf, out2 + n * 3, buf, _in + _i * 3);
              n++;
            }
          } else // same side
            {
              // add the i'th poPoly to the right polygon. Addition is done even for points on the dividing line
              if (d[_i] >= 0) {
                RecastVectors$1.copy4(buf, out1 + m * 3, buf, _in + _i * 3);
                m++;
                if (d[_i] != 0) continue;
              }

              RecastVectors$1.copy4(buf, out2 + n * 3, buf, _in + _i * 3);
              n++;
            }
        }

        return [m, n];
      }
    }, {
      key: "rasterizeTri",
      value: function rasterizeTri(verts, v0, v1, v2, area, hf, bmin, bmax, cs, ics, ich, flagMergeThr) {
        var w = hf.width;
        var h = hf.height;
        var tmin = new Array(3);
        var tmax = new Array(3);
        var by = bmax[1] - bmin[1]; // Calculate the bounding box of the triangle.

        RecastVectors$1.copy3(tmin, verts, v0 * 3);
        RecastVectors$1.copy3(tmax, verts, v0 * 3);
        RecastVectors$1.min(tmin, verts, v1 * 3);
        RecastVectors$1.min(tmin, verts, v2 * 3);
        RecastVectors$1.max(tmax, verts, v1 * 3);
        RecastVectors$1.max(tmax, verts, v2 * 3); // If the triangle does not touch the bbox of the heightfield, skip the triagle.

        if (!RecastRasterization.overlapBounds(bmin, bmax, tmin, tmax)) return; // Calculate the footprPoly of the triangle on the grid's y-axis

        var y0 = Math.floor((tmin[2] - bmin[2]) * ics);
        var y1 = Math.floor((tmax[2] - bmin[2]) * ics);
        y0 = RecastCommon.clamp(y0, 0, h - 1);
        y1 = RecastCommon.clamp(y1, 0, h - 1); // Clip the triangle into all grid cells it touches.

        var buf = new Array(7 * 3 * 4).fill(0);
        var _in = 0;
        var inrow = 7 * 3;
        var p1 = inrow + 7 * 3;
        var p2 = p1 + 7 * 3;
        RecastVectors$1.copy4(buf, 0, verts, v0 * 3);
        RecastVectors$1.copy4(buf, 3, verts, v1 * 3);
        RecastVectors$1.copy4(buf, 6, verts, v2 * 3);
        var nvrow,
            nvIn = 3;

        for (var y = y0; y <= y1; ++y) {
          // Clip polygon to row. Store the remaining polygon as well
          var cz = bmin[2] + y * cs;
          var nvrowin = RecastRasterization.dividePoly(buf, _in, nvIn, inrow, p1, cz + cs, 2);
          nvrow = nvrowin[0];
          nvIn = nvrowin[1];
          {
            var temp = _in;
            _in = p1;
            p1 = temp;
          }
          if (nvrow < 3) continue; // find the horizontal bounds _in the row

          var minX = buf[inrow];
          var maxX = buf[inrow];

          for (var i = 1; i < nvrow; ++i) {
            if (minX > buf[inrow + i * 3]) minX = buf[inrow + i * 3];
            if (maxX < buf[inrow + i * 3]) maxX = buf[inrow + i * 3];
          }

          var x0 = Math.floor((minX - bmin[0]) * ics);
          var x1 = Math.floor((maxX - bmin[0]) * ics);
          x0 = RecastCommon.clamp(x0, 0, w - 1);
          x1 = RecastCommon.clamp(x1, 0, w - 1);
          var nv = void 0,
              nv2 = nvrow;

          for (var x = x0; x <= x1; ++x) {
            // Clip polygon to column. store the remaining polygon as well
            var cx = bmin[0] + x * cs;
            var nvnv2 = RecastRasterization.dividePoly(buf, inrow, nv2, p1, p2, cx + cs, 0);
            nv = nvnv2[0];
            nv2 = nvnv2[1];
            {
              var _temp = inrow;
              inrow = p2;
              p2 = _temp;
            }
            if (nv < 3) continue; // Calculate min and max of the span.

            var smin = buf[p1 + 1];
            var smax = buf[p1 + 1];

            for (var _i2 = 1; _i2 < nv; ++_i2) {
              smin = Math.min(smin, buf[p1 + _i2 * 3 + 1]);
              smax = Math.max(smax, buf[p1 + _i2 * 3 + 1]);
            }

            smin -= bmin[1];
            smax -= bmin[1]; // Skip the span if it is outside the heightfield bbox

            if (smax < 0.0) continue;
            if (smin > by) continue; // Clamp the span to the heightfield bbox.

            if (smin < 0.0) smin = 0;
            if (smax > by) smax = by; // Snap the span to the heightfield height grid.

            var ismin = RecastCommon.clamp(Math.floor(smin * ich), 0, RecastConstants.RC_SPAN_MAX_HEIGHT);
            var ismax = RecastCommon.clamp(Math.ceil(smax * ich), ismin + 1, RecastConstants.RC_SPAN_MAX_HEIGHT);
            RecastRasterization.addSpan(hf, x, y, ismin, ismax, area, flagMergeThr);
          }
        }
      }
      /**
       * No spans will be added if the triangle does not overlap the heightfield grid.
       * 
       * @see Heightfield
       */

    }, {
      key: "rasterizeTriangle",
      value: function rasterizeTriangle(ctx, verts, v0, v1, v2, area, solid, flagMergeThr) {
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

    }, {
      key: "rasterizeTrianglesA",
      value: function rasterizeTrianglesA(ctx, verts, tris, areas, nt, solid, flagMergeThr) {
        ctx.startTimer("RASTERIZE_TRIANGLES");
        var ics = 1.0 / solid.cs;
        var ich = 1.0 / solid.ch; // Rasterize triangles.

        for (var i = 0; i < nt; ++i) {
          var v0 = tris[i * 3 + 0];
          var v1 = tris[i * 3 + 1];
          var v2 = tris[i * 3 + 2]; // Rasterize.

          RecastRasterization.rasterizeTri(verts, v0, v1, v2, areas[i], solid, solid.bmin, solid.bmax, solid.cs, ics, ich, flagMergeThr);
        }

        ctx.stopTimer("RASTERIZE_TRIANGLES");
      }
      /**
       * Spans will only be added for triangles that overlap the heightfield grid.
       * 
       * @see Heightfield
       */

    }, {
      key: "rasterizeTrianglesB",
      value: function rasterizeTrianglesB(ctx, verts, areas, nt, solid, flagMergeThr) {
        ctx.startTimer("RASTERIZE_TRIANGLES");
        var ics = 1.0 / solid.cs;
        var ich = 1.0 / solid.ch; // Rasterize triangles.

        for (var i = 0; i < nt; ++i) {
          var v0 = i * 3 + 0;
          var v1 = i * 3 + 1;
          var v2 = i * 3 + 2; // Rasterize.

          RecastRasterization.rasterizeTri(verts, v0, v1, v2, areas[i], solid, solid.bmin, solid.bmax, solid.cs, ics, ich, flagMergeThr);
        }

        ctx.stopTimer("RASTERIZE_TRIANGLES");
      }
    }]);

    return RecastRasterization;
  }();

  var RecastFilter = /*#__PURE__*/function () {
    function RecastFilter() {
      _classCallCheck(this, RecastFilter);
    }

    _createClass(RecastFilter, null, [{
      key: "filterLowHangingWalkableObstacles",
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
      value: function filterLowHangingWalkableObstacles(ctx, walkableClimb, solid) {
        ctx.startTimer("FILTER_LOW_OBSTACLES");
        var w = solid.width;
        var h = solid.height;

        for (var y = 0; y < h; ++y) {
          for (var x = 0; x < w; ++x) {
            var ps = null;
            var previousWalkable = false;
            var previousArea = RecastConstants.RC_NULL_AREA;

            for (var s = solid.spans[x + y * w]; s != null; ps = s, s = s.next) {
              var walkable = s.area != RecastConstants.RC_NULL_AREA; // If current span is not walkable, but there is walkable
              // span just below it, mark the span above it walkable too.

              if (!walkable && previousWalkable) {
                if (Math.abs(s.smax - ps.smax) <= walkableClimb) s.area = previousArea;
              } // Copy walkable flag so that it cannot propagate
              // past multiple non-walkable objects.


              previousWalkable = walkable;
              previousArea = s.area;
            }
          }
        }

        ctx.stopTimer("FILTER_LOW_OBSTACLES");
      } /// @par
      ///
      /// A ledge is a span with one or more neighbors whose maximum is further away than @p walkableClimb
      /// from the current span's maximum.
      /// This method removes the impact of the overestimation of conservative voxelization 
      /// so the resulting mesh will not have regions hanging in the air over ledges.
      /// 
      /// A span is a ledge if: <tt>rcAbs(currentSpan.smax - neighborSpan.smax) > walkableClimb</tt>
      /// 
      /// @see rcHeightfield, rcConfig

    }, {
      key: "filterLedgeSpans",
      value: function filterLedgeSpans(ctx, walkableHeight, walkableClimb, solid) {
        ctx.startTimer("FILTER_BORDER");
        var w = solid.width;
        var h = solid.height;
        var MAX_HEIGHT = 0xfff; // Mark border spans.

        for (var y = 0; y < h; ++y) {
          for (var x = 0; x < w; ++x) {
            for (var s = solid.spans[x + y * w]; s != null; s = s.next) {
              // Skip non walkable spans.
              if (s.area == RecastConstants.RC_NULL_AREA) continue;
              var bot = s.smax;
              var top = s.next != null ? s.next.smin : MAX_HEIGHT; // Find neighbours minimum height.

              var minh = MAX_HEIGHT; // Min and max height of accessible neighbours.

              var asmin = s.smax;
              var asmax = s.smax;

              for (var dir = 0; dir < 4; ++dir) {
                var dx = x + RecastCommon.GetDirOffsetX(dir);
                var dy = y + RecastCommon.GetDirOffsetY(dir); // Skip neighbours which are out of bounds.

                if (dx < 0 || dy < 0 || dx >= w || dy >= h) {
                  minh = Math.min(minh, -walkableClimb - bot);
                  continue;
                } // From minus infinity to the first span.


                var ns = solid.spans[dx + dy * w];
                var nbot = -walkableClimb;
                var ntop = ns != null ? ns.smin : MAX_HEIGHT; // Skip neightbour if the gap between the spans is too small.

                if (Math.min(top, ntop) - Math.max(bot, nbot) > walkableHeight) minh = Math.min(minh, nbot - bot); // Rest of the spans.

                for (var _ns = solid.spans[dx + dy * w]; _ns != null; _ns = _ns.next) {
                  nbot = _ns.smax;
                  ntop = _ns.next != null ? _ns.next.smin : MAX_HEIGHT; // Skip neightbour if the gap between the spans is too small.

                  if (Math.min(top, ntop) - Math.max(bot, nbot) > walkableHeight) {
                    minh = Math.min(minh, nbot - bot); // Find min/max accessible neighbour height. 

                    if (Math.abs(nbot - bot) <= walkableClimb) {
                      if (nbot < asmin) asmin = nbot;
                      if (nbot > asmax) asmax = nbot;
                    }
                  }
                }
              } // The current span is close to a ledge if the drop to any
              // neighbour span is less than the walkableClimb.


              if (minh < -walkableClimb) s.area = RecastConstants.RC_NULL_AREA; // If the difference between all neighbours is too large,
              // we are at steep slope, mark the span as ledge.

              if (asmax - asmin > walkableClimb) {
                s.area = RecastConstants.RC_NULL_AREA;
              }
            }
          }
        }

        ctx.stopTimer("FILTER_BORDER");
      } /// @par
      ///
      /// For this filter, the clearance above the span is the distance from the span's 
      /// maximum to the next higher span's minimum. (Same grid column.)
      /// 
      /// @see rcHeightfield, rcConfig

    }, {
      key: "filterWalkableLowHeightSpans",
      value: function filterWalkableLowHeightSpans(ctx, walkableHeight, solid) {
        ctx.startTimer("FILTER_WALKABLE");
        var w = solid.width;
        var h = solid.height;
        var MAX_HEIGHT = 0xfff; // Remove walkable flag from spans which do not have enough
        // space above them for the agent to stand there.

        for (var y = 0; y < h; ++y) {
          for (var x = 0; x < w; ++x) {
            for (var s = solid.spans[x + y * w]; s != null; s = s.next) {
              var bot = s.smax;
              var top = s.next != null ? s.next.smin : MAX_HEIGHT;
              if (top - bot <= walkableHeight) s.area = RecastConstants.RC_NULL_AREA;
            }
          }
        }

        ctx.stopTimer("FILTER_WALKABLE");
      }
    }]);

    return RecastFilter;
  }();

  var RecastArea = /*#__PURE__*/function () {
    function RecastArea() {
      _classCallCheck(this, RecastArea);
    }

    _createClass(RecastArea, [{
      key: "medianFilterWalkableArea",
      /// @par
      ///
      /// This filter is usually applied after applying area id's using functions
      /// such as #rcMarkBoxArea, #rcMarkConvexPolyArea, and #rcMarkCylinderArea.
      ///
      /// @see rcCompactHeightfield
      value: function medianFilterWalkableArea(ctx, chf) {
        var w = chf.width;
        var h = chf.height;
        ctx.startTimer("MEDIAN_AREA");
        var areas = new Array(chf.spanCount);

        for (var y = 0; y < h; ++y) {
          for (var x = 0; x < w; ++x) {
            var c = chf.cells[x + y * w];

            for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
              s = chf.spans[i];

              if (chf.areas[i] == RecastConstants.RC_NULL_AREA) {
                areas[i] = chf.areas[i];
                continue;
              }

              nei = new Array(9);

              for (var j = 0; j < 9; ++j) {
                nei[j] = chf.areas[i];
              }

              for (var dir = 0; dir < 4; ++dir) {
                if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
                  var ax = x + RecastCommon.GetDirOffsetX(dir);
                  var ay = y + RecastCommon.GetDirOffsetY(dir);
                  var ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);
                  if (chf.areas[ai] != RecastConstants.RC_NULL_AREA) nei[dir * 2 + 0] = chf.areas[ai];
                  var as = chf.spans[ai];
                  var dir2 = dir + 1 & 0x3;

                  if (RecastCommon.GetCon(as, dir2) != RecastConstants.RC_NOT_CONNECTED) {
                    var ax2 = ax + RecastCommon.GetDirOffsetX(dir2);
                    var ay2 = ay + RecastCommon.GetDirOffsetY(dir2);
                    var ai2 = chf.cells[ax2 + ay2 * w].index + RecastCommon.GetCon(as, dir2);
                    if (chf.areas[ai2] != RecastConstants.RC_NULL_AREA) nei[dir * 2 + 1] = chf.areas[ai2];
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
      } /// @par
      ///
      /// The value of spacial parameters are in world units.
      ///
      /// @see rcCompactHeightfield, rcMedianFilterWalkableArea

    }, {
      key: "markBoxArea",
      value: function markBoxArea(ctx, bmin, bmax, areaMod, chf) {
        ctx.startTimer("MARK_BOX_AREA");
        var minx = Math.floor((bmin[0] - chf.bmin[0]) / chf.cs);
        var miny = Math.floor((bmin[1] - chf.bmin[1]) / chf.ch);
        var minz = Math.floor((bmin[2] - chf.bmin[2]) / chf.cs);
        var maxx = Math.floor((bmax[0] - chf.bmin[0]) / chf.cs);
        var maxy = Math.floor((bmax[1] - chf.bmin[1]) / chf.ch);
        var maxz = Math.floor((bmax[2] - chf.bmin[2]) / chf.cs);
        if (maxx < 0) return;
        if (minx >= chf.width) return;
        if (maxz < 0) return;
        if (minz >= chf.height) return;
        if (minx < 0) minx = 0;
        if (maxx >= chf.width) maxx = chf.width - 1;
        if (minz < 0) minz = 0;
        if (maxz >= chf.height) maxz = chf.height - 1;

        for (var z = minz; z <= maxz; ++z) {
          for (var x = minx; x <= maxx; ++x) {
            var c = chf.cells[x + z * chf.width];

            for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
              var _s = chf.spans[i];

              if (_s.y >= miny && _s.y <= maxy) {
                if (chf.areas[i] != RecastConstants.RC_NULL_AREA) chf.areas[i] = areaMod.apply(chf.areas[i]);
              }
            }
          }
        }

        ctx.stopTimer("MARK_BOX_AREA");
      }
    }, {
      key: "offsetPoly",
      value: function offsetPoly(verts, nverts, offset, outVerts, maxOutVerts) {
        MITER_LIMIT = 1.20;
        var n = 0;

        for (var i = 0; i < nverts; i++) {
          var a = (i + nverts - 1) % nverts;
          var b = i;
          var c = (i + 1) % nverts;
          var va = a * 3;
          var vb = b * 3;
          var vc = c * 3;
          dx0 = verts[vb] - verts[va];
          dy0 = verts[vb + 2] - verts[va + 2];
          d0 = dx0 * dx0 + dy0 * dy0;

          if (d0 > 1e-6) {
            d0 = 1.0 / Math.sqrt(d0);
            dx0 *= d0;
            dy0 *= d0;
          }

          dx1 = verts[vc] - verts[vb];
          dy1 = verts[vc + 2] - verts[vb + 2];
          d1 = dx1 * dx1 + dy1 * dy1;

          if (d1 > 1e-6) {
            d1 = 1.0 / Math.sqrt(d1);
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
          var bevel = dmr2 * MITER_LIMIT * MITER_LIMIT < 1.0;

          if (dmr2 > 1e-6) {
            scale = 1.0 / dmr2;
            dmx *= scale;
            dmy *= scale;
          }

          if (bevel && cross < 0.0) {
            if (n + 2 >= maxOutVerts) return 0;
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
            if (n + 1 >= maxOutVerts) return 0;
            outVerts[n * 3 + 0] = verts[vb] - dmx * offset;
            outVerts[n * 3 + 1] = verts[vb + 1];
            outVerts[n * 3 + 2] = verts[vb + 2] - dmy * offset;
            n++;
          }
        }

        return n;
      } /// @par
      ///
      /// The value of spacial parameters are in world units.
      ///
      /// @see rcCompactHeightfield, rcMedianFilterWalkableArea

    }, {
      key: "markCylinderArea",
      value: function markCylinderArea(ctx, pos, r, h, areaMod, chf) {
        ctx.startTimer("MARK_CYLINDER_AREA");
        bmin = new Array(3);
        bmax = new Array(3);
        bmin[0] = pos[0] - r;
        bmin[1] = pos[1];
        bmin[2] = pos[2] - r;
        bmax[0] = pos[0] + r;
        bmax[1] = pos[1] + h;
        bmax[2] = pos[2] + r;
        r2 = r * r;
        var minx = Math.floor((bmin[0] - chf.bmin[0]) / chf.cs);
        var miny = Math.floor((bmin[1] - chf.bmin[1]) / chf.ch);
        var minz = Math.floor((bmin[2] - chf.bmin[2]) / chf.cs);
        var maxx = Math.floor((bmax[0] - chf.bmin[0]) / chf.cs);
        var maxy = Math.floor((bmax[1] - chf.bmin[1]) / chf.ch);
        var maxz = Math.floor((bmax[2] - chf.bmin[2]) / chf.cs);
        if (maxx < 0) return;
        if (minx >= chf.width) return;
        if (maxz < 0) return;
        if (minz >= chf.height) return;
        if (minx < 0) minx = 0;
        if (maxx >= chf.width) maxx = chf.width - 1;
        if (minz < 0) minz = 0;
        if (maxz >= chf.height) maxz = chf.height - 1;

        for (var z = minz; z <= maxz; ++z) {
          for (var x = minx; x <= maxx; ++x) {
            var c = chf.cells[x + z * chf.width];

            for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
              var _s2 = chf.spans[i];
              if (chf.areas[i] == RecastConstants.RC_NULL_AREA) continue;

              if (_s2.y >= miny && _s2.y <= maxy) {
                sx = chf.bmin[0] + (x + 0.5) * chf.cs;
                sz = chf.bmin[2] + (z + 0.5) * chf.cs;
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
    }], [{
      key: "erodeWalkableArea",
      /// @par
      ///
      /// Basically, any spans that are closer to a boundary or obstruction than the specified radius
      /// are marked as unwalkable.
      ///
      /// This method is usually called immediately after the heightfield has been built.
      ///
      /// @see rcCompactHeightfield, rcBuildCompactHeightfield, rcConfig::walkableRadius
      value: function erodeWalkableArea(ctx, radius, chf) {
        var w = chf.width;
        var h = chf.height;
        ctx.startTimer("ERODE_AREA");
        var dist = new Array(chf.spanCount).fill(255); // Arrays.fill(dist, 255);
        // Mark boundary cells.

        for (var y = 0; y < h; ++y) {
          for (var x = 0; x < w; ++x) {
            var c = chf.cells[x + y * w];

            for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
              if (chf.areas[i] == RecastConstants.RC_NULL_AREA) {
                dist[i] = 0;
              } else {
                var _s3 = chf.spans[i];
                var nc = 0;

                for (var dir = 0; dir < 4; ++dir) {
                  if (RecastCommon.GetCon(_s3, dir) != RecastConstants.RC_NOT_CONNECTED) {
                    var nx = x + RecastCommon.GetDirOffsetX(dir);
                    var ny = y + RecastCommon.GetDirOffsetY(dir);
                    var nidx = chf.cells[nx + ny * w].index + RecastCommon.GetCon(_s3, dir);

                    if (chf.areas[nidx] != RecastConstants.RC_NULL_AREA) {
                      nc++;
                    }
                  }
                } // At least one missing neighbour.


                if (nc != 4) dist[i] = 0;
              }
            }
          }
        }

        var nd; // console.log(chf.spans[450240]);
        // Pass 1

        for (var _y = 0; _y < h; ++_y) {
          for (var _x = 0; _x < w; ++_x) {
            var _c = chf.cells[_x + _y * w];

            for (var _i = _c.index, _ni = _c.index + _c.count; _i < _ni; ++_i) {
              var _s4 = chf.spans[_i];

              if (RecastCommon.GetCon(_s4, 0) != RecastConstants.RC_NOT_CONNECTED) {
                // (-1,0)
                var ax = _x + RecastCommon.GetDirOffsetX(0);

                var ay = _y + RecastCommon.GetDirOffsetY(0);

                var ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(_s4, 0);
                var as = chf.spans[ai];
                nd = Math.min(dist[ai] + 2, 255);
                if (nd < dist[_i]) dist[_i] = nd; // (-1,-1)

                if (RecastCommon.GetCon(as, 3) != RecastConstants.RC_NOT_CONNECTED) {
                  var aax = ax + RecastCommon.GetDirOffsetX(3);
                  var aay = ay + RecastCommon.GetDirOffsetY(3);
                  var aai = chf.cells[aax + aay * w].index + RecastCommon.GetCon(as, 3);
                  nd = Math.min(dist[aai] + 3, 255);
                  if (nd < dist[_i]) dist[_i] = nd;
                }
              }

              if (RecastCommon.GetCon(_s4, 3) != RecastConstants.RC_NOT_CONNECTED) {
                // (0,-1)
                var _ax = _x + RecastCommon.GetDirOffsetX(3);

                var _ay = _y + RecastCommon.GetDirOffsetY(3);

                var _ai = chf.cells[_ax + _ay * w].index + RecastCommon.GetCon(_s4, 3);

                var _as = chf.spans[_ai];
                nd = Math.min(dist[_ai] + 2, 255);
                if (nd < dist[_i]) dist[_i] = nd; // (1,-1)

                if (RecastCommon.GetCon(_as, 2) != RecastConstants.RC_NOT_CONNECTED) {
                  var _aax = _ax + RecastCommon.GetDirOffsetX(2);

                  var _aay = _ay + RecastCommon.GetDirOffsetY(2);

                  var _aai = chf.cells[_aax + _aay * w].index + RecastCommon.GetCon(_as, 2);

                  nd = Math.min(dist[_aai] + 3, 255);
                  if (nd < dist[_i]) dist[_i] = nd;
                }
              }
            }
          }
        } // Pass 2


        for (var _y2 = h - 1; _y2 >= 0; --_y2) {
          for (var _x2 = w - 1; _x2 >= 0; --_x2) {
            // if(y == 671&& x == 671)
            // 	console.log("There")
            var _c2 = chf.cells[_x2 + _y2 * w];

            for (var _i2 = _c2.index, _ni2 = _c2.index + _c2.count; _i2 < _ni2; ++_i2) {
              var _s5 = chf.spans[_i2];

              if (RecastCommon.GetCon(_s5, 2) != RecastConstants.RC_NOT_CONNECTED) {
                // (1,0)
                var _ax2 = _x2 + RecastCommon.GetDirOffsetX(2);

                var _ay2 = _y2 + RecastCommon.GetDirOffsetY(2);

                var _ai2 = chf.cells[_ax2 + _ay2 * w].index + RecastCommon.GetCon(_s5, 2);

                if (_ai2 == 450241) console.log("Here");
                var _as2 = chf.spans[_ai2];
                nd = Math.min(dist[_ai2] + 2, 255);
                if (nd < dist[_i2]) dist[_i2] = nd; // (1,1)

                if (RecastCommon.GetCon(_as2, 1) != RecastConstants.RC_NOT_CONNECTED) {
                  var _aax2 = _ax2 + RecastCommon.GetDirOffsetX(1);

                  var _aay2 = _ay2 + RecastCommon.GetDirOffsetY(1);

                  var _aai2 = chf.cells[_aax2 + _aay2 * w].index + RecastCommon.GetCon(_as2, 1);

                  nd = Math.min(dist[_aai2] + 3, 255);
                  if (nd < dist[_i2]) dist[_i2] = nd;
                }
              }

              if (RecastCommon.GetCon(_s5, 1) != RecastConstants.RC_NOT_CONNECTED) {
                // (0,1)
                var _ax3 = _x2 + RecastCommon.GetDirOffsetX(1);

                var _ay3 = _y2 + RecastCommon.GetDirOffsetY(1);

                var _ai3 = chf.cells[_ax3 + _ay3 * w].index + RecastCommon.GetCon(_s5, 1);

                var _as3 = chf.spans[_ai3];
                nd = Math.min(dist[_ai3] + 2, 255);
                if (nd < dist[_i2]) dist[_i2] = nd; // (-1,1)

                if (RecastCommon.GetCon(_as3, 0) != RecastConstants.RC_NOT_CONNECTED) {
                  var _aax3 = _ax3 + RecastCommon.GetDirOffsetX(0);

                  var _aay3 = _ay3 + RecastCommon.GetDirOffsetY(0);

                  var _aai3 = chf.cells[_aax3 + _aay3 * w].index + RecastCommon.GetCon(_as3, 0);

                  nd = Math.min(dist[_aai3] + 3, 255);
                  if (nd < dist[_i2]) dist[_i2] = nd;
                }
              }
            }
          }
        }

        var thr = radius * 2;

        for (var _i3 = 0; _i3 < chf.spanCount; ++_i3) {
          if (dist[_i3] < thr) chf.areas[_i3] = RecastConstants.RC_NULL_AREA;
        }

        ctx.stopTimer("ERODE_AREA");
      }
    }, {
      key: "pointInPoly",
      value: function pointInPoly(verts, p) {
        var c = false;

        for (var _i4 = 0, _j = verts.length - 1; _i4 < verts.length; _j = _i4++) {
          var vi = _i4 * 3;
          var vj = _j * 3;
          if (verts[vi + 2] > p[2] != verts[vj + 2] > p[2] && p[0] < (verts[vj] - verts[vi]) * (p[2] - verts[vi + 2]) / (verts[vj + 2] - verts[vi + 2]) + verts[vi]) c = !c;
        }

        return c;
      } /// @par
      ///
      /// The value of spacial parameters are in world units.
      ///
      /// The y-values of the polygon vertices are ignored. So the polygon is effectively
      /// projected onto the xz-plane at @p hmin, then extruded to @p hmax.
      ///
      /// @see rcCompactHeightfield, rcMedianFilterWalkableArea

    }, {
      key: "markConvexPolyArea",
      value: function markConvexPolyArea(ctx, verts, hmin, hmax, areaMod, chf) {
        ctx.startTimer("MARK_CONVEXPOLY_AREA");
        bmin = new Array(3);
        max = new Array(3);
        RecastVectors.copy(bmin, verts, 0);
        RecastVectors.copy(bmax, verts, 0);

        for (var i = 1; i < verts.length; ++i) {
          RecastVectors.min(bmin, verts, i * 3);
          RecastVectors.max(bmax, verts, i * 3);
        }

        bmin[1] = hmin;
        bmax[1] = hmax;
        var minx = Math.floor((bmin[0] - chf.bmin[0]) / chf.cs);
        var miny = Math.floor((bmin[1] - chf.bmin[1]) / chf.ch);
        var minz = Math.floor((bmin[2] - chf.bmin[2]) / chf.cs);
        var maxx = Math.floor((bmax[0] - chf.bmin[0]) / chf.cs);
        var maxy = Math.floor((bmax[1] - chf.bmin[1]) / chf.ch);
        var maxz = Math.floor((bmax[2] - chf.bmin[2]) / chf.cs);
        if (maxx < 0) return;
        if (minx >= chf.width) return;
        if (maxz < 0) return;
        if (minz >= chf.height) return;
        if (minx < 0) minx = 0;
        if (maxx >= chf.width) maxx = chf.width - 1;
        if (minz < 0) minz = 0;
        if (maxz >= chf.height) maxz = chf.height - 1; // TODO: Optimize.

        for (var z = minz; z <= maxz; ++z) {
          for (var x = minx; x <= maxx; ++x) {
            var c = chf.cells[x + z * chf.width];

            for (var _i5 = c.index, ni = c.index + c.count; _i5 < ni; ++_i5) {
              var _s6 = chf.spans[_i5];
              if (chf.areas[_i5] == RecastConstants.RC_NULL_AREA) continue;

              if (_s6.y >= miny && _s6.y <= maxy) {
                p = new Array(3);
                p[0] = chf.bmin[0] + (x + 0.5) * chf.cs;
                p[1] = 0;
                p[2] = chf.bmin[2] + (z + 0.5) * chf.cs;

                if (pointInPoly(verts, p)) {
                  chf.areas[_i5] = areaMod.apply(chf.areas[_i5]);
                }
              }
            }
          }
        }

        ctx.stopTimer("MARK_CONVEXPOLY_AREA");
      }
    }]);

    return RecastArea;
  }();

  var _temp$1, _temp2$1;

  var RecastRegion = /*#__PURE__*/function () {
    function RecastRegion() {
      _classCallCheck(this, RecastRegion);
    }

    _createClass(RecastRegion, null, [{
      key: "calculateDistanceField",
      value: function calculateDistanceField(chf, src) {
        var maxDist;
        var w = chf.width;
        var h = chf.height; // Init distance and points.

        for (var i = 0; i < chf.spanCount; ++i) {
          src[i] = 0xffff;
        } // Mark boundary cells.


        for (var y = 0; y < h; ++y) {
          for (var _x = 0; _x < w; ++_x) {
            var c = chf.cells[_x + y * w];

            for (var _i = c.index, ni = c.index + c.count; _i < ni; ++_i) {
              var s = chf.spans[_i];
              var area = chf.areas[_i];
              var nc = 0;

              for (var dir = 0; dir < 4; ++dir) {
                if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
                  var ax = _x + RecastCommon.GetDirOffsetX(dir);

                  var ay = y + RecastCommon.GetDirOffsetY(dir);
                  var ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);

                  if (area == chf.areas[ai]) {
                    nc++;
                  }
                }
              }

              if (nc != 4) {
                src[_i] = 0;
              }
            }
          }
        } // Pass 1


        for (var _y = 0; _y < h; ++_y) {
          for (var _x2 = 0; _x2 < w; ++_x2) {
            var _c = chf.cells[_x2 + _y * w];

            for (var _i2 = _c.index, _ni = _c.index + _c.count; _i2 < _ni; ++_i2) {
              var _s = chf.spans[_i2];

              if (RecastCommon.GetCon(_s, 0) != RecastConstants.RC_NOT_CONNECTED) {
                // (-1,0)
                var _ax = _x2 + RecastCommon.GetDirOffsetX(0);

                var _ay = _y + RecastCommon.GetDirOffsetY(0);

                var _ai = chf.cells[_ax + _ay * w].index + RecastCommon.GetCon(_s, 0);

                var as = chf.spans[_ai];

                if (src[_ai] + 2 < src[_i2]) {
                  src[_i2] = src[_ai] + 2;
                } // (-1,-1)


                if (RecastCommon.GetCon(as, 3) != RecastConstants.RC_NOT_CONNECTED) {
                  var aax = _ax + RecastCommon.GetDirOffsetX(3);

                  var aay = _ay + RecastCommon.GetDirOffsetY(3);

                  var aai = chf.cells[aax + aay * w].index + RecastCommon.GetCon(as, 3);

                  if (src[aai] + 3 < src[_i2]) {
                    src[_i2] = src[aai] + 3;
                  }
                }
              }

              if (RecastCommon.GetCon(_s, 3) != RecastConstants.RC_NOT_CONNECTED) {
                // (0,-1)
                var _ax2 = _x2 + RecastCommon.GetDirOffsetX(3);

                var _ay2 = _y + RecastCommon.GetDirOffsetY(3);

                var _ai2 = chf.cells[_ax2 + _ay2 * w].index + RecastCommon.GetCon(_s, 3);

                var _as = chf.spans[_ai2];

                if (src[_ai2] + 2 < src[_i2]) {
                  src[_i2] = src[_ai2] + 2;
                } // (1,-1)


                if (RecastCommon.GetCon(_as, 2) != RecastConstants.RC_NOT_CONNECTED) {
                  var _aax = _ax2 + RecastCommon.GetDirOffsetX(2);

                  var _aay = _ay2 + RecastCommon.GetDirOffsetY(2);

                  var _aai = chf.cells[_aax + _aay * w].index + RecastCommon.GetCon(_as, 2);

                  if (src[_aai] + 3 < src[_i2]) {
                    src[_i2] = src[_aai] + 3;
                  }
                }
              }
            }
          }
        } // Pass 2


        for (var _y2 = h - 1; _y2 >= 0; --_y2) {
          for (var _x3 = w - 1; _x3 >= 0; --_x3) {
            var _c2 = chf.cells[_x3 + _y2 * w];

            for (var _i3 = _c2.index, _ni2 = _c2.index + _c2.count; _i3 < _ni2; ++_i3) {
              var _s2 = chf.spans[_i3];

              if (RecastCommon.GetCon(_s2, 2) != RecastConstants.RC_NOT_CONNECTED) {
                // (1,0)
                var _ax3 = _x3 + RecastCommon.GetDirOffsetX(2);

                var _ay3 = _y2 + RecastCommon.GetDirOffsetY(2);

                var _ai3 = chf.cells[_ax3 + _ay3 * w].index + RecastCommon.GetCon(_s2, 2);

                var _as2 = chf.spans[_ai3];

                if (src[_ai3] + 2 < src[_i3]) {
                  src[_i3] = src[_ai3] + 2;
                } // (1,1)


                if (RecastCommon.GetCon(_as2, 1) != RecastConstants.RC_NOT_CONNECTED) {
                  var _aax2 = _ax3 + RecastCommon.GetDirOffsetX(1);

                  var _aay2 = _ay3 + RecastCommon.GetDirOffsetY(1);

                  var _aai2 = chf.cells[_aax2 + _aay2 * w].index + RecastCommon.GetCon(_as2, 1);

                  if (src[_aai2] + 3 < src[_i3]) {
                    src[_i3] = src[_aai2] + 3;
                  }
                }
              }

              if (RecastCommon.GetCon(_s2, 1) != RecastConstants.RC_NOT_CONNECTED) {
                // (0,1)
                var _ax4 = _x3 + RecastCommon.GetDirOffsetX(1);

                var _ay4 = _y2 + RecastCommon.GetDirOffsetY(1);

                var _ai4 = chf.cells[_ax4 + _ay4 * w].index + RecastCommon.GetCon(_s2, 1);

                var _as3 = chf.spans[_ai4];

                if (src[_ai4] + 2 < src[_i3]) {
                  src[_i3] = src[_ai4] + 2;
                } // (-1,1)


                if (RecastCommon.GetCon(_as3, 0) != RecastConstants.RC_NOT_CONNECTED) {
                  var _aax3 = _ax4 + RecastCommon.GetDirOffsetX(0);

                  var _aay3 = _ay4 + RecastCommon.GetDirOffsetY(0);

                  var _aai3 = chf.cells[_aax3 + _aay3 * w].index + RecastCommon.GetCon(_as3, 0);

                  if (src[_aai3] + 3 < src[_i3]) {
                    src[_i3] = src[_aai3] + 3;
                  }
                }
              }
            }
          }
        }

        maxDist = 0;

        for (var _i4 = 0; _i4 < chf.spanCount; ++_i4) {
          maxDist = Math.max(src[_i4], maxDist);
        }

        return maxDist;
      }
    }, {
      key: "boxBlur",
      value: function boxBlur(chf, thr, src) {
        var w = chf.width;
        var h = chf.height;
        var dst = new Array(chf.spanCount);
        thr *= 2;

        for (var y = 0; y < h; ++y) {
          for (var _x4 = 0; _x4 < w; ++_x4) {
            var c = chf.cells[_x4 + y * w];

            for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
              var s = chf.spans[i];
              var cd = src[i];

              if (cd <= thr) {
                dst[i] = cd;
                continue;
              }

              var d = cd;

              for (var dir = 0; dir < 4; ++dir) {
                if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
                  var ax = _x4 + RecastCommon.GetDirOffsetX(dir);

                  var ay = y + RecastCommon.GetDirOffsetY(dir);
                  var ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);
                  d += src[ai];
                  var as = chf.spans[ai];
                  var dir2 = dir + 1 & 0x3;

                  if (RecastCommon.GetCon(as, dir2) != RecastConstants.RC_NOT_CONNECTED) {
                    var ax2 = ax + RecastCommon.GetDirOffsetX(dir2);
                    var ay2 = ay + RecastCommon.GetDirOffsetY(dir2);
                    var ai2 = chf.cells[ax2 + ay2 * w].index + RecastCommon.GetCon(as, dir2);
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
    }, {
      key: "floodRegion",
      value: function floodRegion(x, y, i, level, r, chf, srcReg, srcDist, stack) {
        var w = chf.width;
        var area = chf.areas[i]; // Flood fill mark region.

        stack = [];
        stack.push(x);
        stack.push(y);
        stack.push(i);
        srcReg[i] = r;
        srcDist[i] = 0;
        var lev = level >= 2 ? level - 2 : 0;
        var count = 0;

        while (stack.length > 0) {
          var ci = stack.splice(stack.length - 1, 1)[0];
          var cy = stack.splice(stack.length - 1, 1)[0];
          var cx = stack.splice(stack.length - 1, 1)[0];
          var cs = chf.spans[ci]; // Check if any of the neighbours already have a valid region set.

          var ar = 0;

          for (var dir = 0; dir < 4; ++dir) {
            // 8 connected
            if (RecastCommon.GetCon(cs, dir) != RecastConstants.RC_NOT_CONNECTED) {
              var ax = cx + RecastCommon.GetDirOffsetX(dir);
              var ay = cy + RecastCommon.GetDirOffsetY(dir);
              var ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(cs, dir);

              if (chf.areas[ai] != area) {
                continue;
              }

              var nr = srcReg[ai];

              if ((nr & RecastConstants.RC_BORDER_REG) != 0) {
                continue;
              }

              if (nr != 0 && nr != r) {
                ar = nr;
                break;
              }

              var as = chf.spans[ai];
              var dir2 = dir + 1 & 0x3;

              if (RecastCommon.GetCon(as, dir2) != RecastConstants.RC_NOT_CONNECTED) {
                var ax2 = ax + RecastCommon.GetDirOffsetX(dir2);
                var ay2 = ay + RecastCommon.GetDirOffsetY(dir2);
                var ai2 = chf.cells[ax2 + ay2 * w].index + RecastCommon.GetCon(as, dir2);

                if (chf.areas[ai2] != area) {
                  continue;
                }

                var nr2 = srcReg[ai2];

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

          count++; // Expand neighbours.

          for (var _dir = 0; _dir < 4; ++_dir) {
            if (RecastCommon.GetCon(cs, _dir) != RecastConstants.RC_NOT_CONNECTED) {
              var _ax5 = cx + RecastCommon.GetDirOffsetX(_dir);

              var _ay5 = cy + RecastCommon.GetDirOffsetY(_dir);

              var _ai5 = chf.cells[_ax5 + _ay5 * w].index + RecastCommon.GetCon(cs, _dir);

              if (chf.areas[_ai5] != area) {
                continue;
              }

              if (chf.dist[_ai5] >= lev && srcReg[_ai5] == 0) {
                srcReg[_ai5] = r;
                srcDist[_ai5] = 0;
                stack.push(_ax5);
                stack.push(_ay5);
                stack.push(_ai5);
              }
            }
          }
        }

        return count > 0;
      }
    }, {
      key: "expandRegions",
      value: function expandRegions(maxIter, level, chf, srcReg, srcDist, stack, fillStack) {
        var w = chf.width;
        var h = chf.height;

        if (fillStack) {
          // Find cells revealed by the raised level.
          //stack = [];
          stack = [];

          for (var y = 0; y < h; ++y) {
            for (var _x5 = 0; _x5 < w; ++_x5) {
              var c = chf.cells[_x5 + y * w];

              for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
                if (chf.dist[i] >= level && srcReg[i] == 0 && chf.areas[i] != RecastConstants.RC_NULL_AREA) {
                  stack.push(_x5);
                  stack.push(y);
                  stack.push(i);
                }
              }
            }
          }
        } else // use cells in the input stack
          {
            // mark all cells which already have a region
            for (var j = 0; j < stack.length; j += 3) {
              var _i5 = stack[j + 2];

              if (srcReg[_i5] != 0) {
                stack[j + 2] = -1;
              }
            }
          }

        var dirtyEntries = [];
        var iter = 0;

        while (stack.length > 0) {
          var failed = 0; // dirtyEntries = [];

          dirtyEntries = [];

          for (var _j = 0; _j < stack.length; _j += 3) {
            var _x6 = stack[_j + 0];
            var _y3 = stack[_j + 1];
            var _i6 = stack[_j + 2];

            if (_i6 < 0) {
              failed++;
              continue;
            }

            var r = srcReg[_i6];
            var d2 = 0xfff;
            var area = chf.areas[_i6];
            var s = chf.spans[_i6];

            for (var dir = 0; dir < 4; ++dir) {
              if (RecastCommon.GetCon(s, dir) == RecastConstants.RC_NOT_CONNECTED) {
                continue;
              }

              var ax = _x6 + RecastCommon.GetDirOffsetX(dir);

              var ay = _y3 + RecastCommon.GetDirOffsetY(dir);

              var ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);

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
              stack[_j + 2] = -1; // mark as used

              dirtyEntries.push(_i6);
              dirtyEntries.push(r);
              dirtyEntries.push(d2);
            } else {
              failed++;
            }
          } // Copy entries that differ between src and dst to keep them in sync.


          for (var _i7 = 0; _i7 < dirtyEntries.length; _i7 += 3) {
            var idx = dirtyEntries[_i7]; // if (idx == 1344)
            //     console.log("dirty")

            srcReg[idx] = dirtyEntries[_i7 + 1];
            srcDist[idx] = dirtyEntries[_i7 + 2];
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
    }, {
      key: "sortCellsByLevel",
      value: function sortCellsByLevel(startLevel, chf, srcReg, nbStacks, stacks, loglevelsPerStack) // the levels per stack (2 in our case) as a bit shift
      {
        var w = chf.width;
        var h = chf.height;
        startLevel = startLevel >> loglevelsPerStack;

        for (var j = 0; j < nbStacks; ++j) {
          // stacks[j] = new Array(1024);
          // stacks[j] = [];
          stacks[j] = [];
        }

        for (var y = 0; y < h; ++y) {
          for (var _x7 = 0; _x7 < w; ++_x7) {
            var c = chf.cells[_x7 + y * w];

            for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
              if (chf.areas[i] == RecastConstants.RC_NULL_AREA || srcReg[i] != 0) {
                continue;
              }

              var level = chf.dist[i] >> loglevelsPerStack;
              var sId = startLevel - level;

              if (sId >= nbStacks) {
                continue;
              }

              if (sId < 0) {
                sId = 0;
              }

              stacks[sId].push(_x7);
              stacks[sId].push(y);
              stacks[sId].push(i);
            }
          }
        }
      }
    }, {
      key: "appendStacks",
      value: function appendStacks(srcStack, dstStack, srcReg) {
        for (var j = 0; j < srcStack.length; j += 3) {
          var i = srcStack[j + 2];

          if (i < 0 || srcReg[i] != 0) {
            continue;
          }

          dstStack.push(srcStack[j]);
          dstStack.push(srcStack[j + 1]);
          dstStack.push(srcStack[j + 2]);
        }
      }
    }, {
      key: "removeAdjacentNeighbours",
      value: function removeAdjacentNeighbours(reg) {
        // Remove adjacent duplicates.
        for (var i = 0; i < reg.connections.length && reg.connections.length > 1;) {
          var ni = (i + 1) % reg.connections.length;

          if (reg.connections[i] == reg.connections[ni] && reg.connections[i] >= -128 && reg.connections[i] <= 127) {
            reg.connections.splice(i, 1);
          } else {
            ++i;
          }
        }
      }
    }, {
      key: "replaceNeighbour",
      value: function replaceNeighbour(reg, oldId, newId) {
        var neiChanged = false;

        for (var i = 0; i < reg.connections.length; ++i) {
          if (reg.connections[i] == oldId) {
            reg.connections[i] = newId;
            neiChanged = true;
          }
        }

        for (var _i8 = 0; _i8 < reg.floors.length; ++_i8) {
          if (reg.floors[_i8] == oldId) {
            reg.floors.set(_i8, newId);
          }
        }

        if (neiChanged) {
          RecastRegion.removeAdjacentNeighbours(reg);
        }
      }
    }, {
      key: "canMergeWithRegion",
      value: function canMergeWithRegion(rega, regb) {
        if (rega.areaType != regb.areaType) {
          return false;
        }

        var n = 0;

        for (var i = 0; i < rega.connections.length; ++i) {
          if (rega.connections[i] == regb.id) {
            n++;
          }
        }

        if (n > 1) {
          return false;
        }

        for (var _i9 = 0; _i9 < rega.floors.length; ++_i9) {
          if (rega.floors[_i9] == regb.id) {
            return false;
          }
        }

        return true;
      }
    }, {
      key: "addUniqueFloorRegion",
      value: function addUniqueFloorRegion(reg, n) {
        if (!reg.floors.includes(n)) {
          reg.floors.push(n);
        }
      }
    }, {
      key: "mergeRegions",
      value: function mergeRegions(rega, regb) {
        var aid = rega.id;
        var bid = regb.id; // Duplicate current neighbourhood.

        var acon = rega.connections;
        var bcon = regb.connections; // Find insertion poPoly on A.

        var insa = -1;

        for (var i = 0; i < acon.length; ++i) {
          if (acon[i] == bid) {
            insa = i;
            break;
          }
        }

        if (insa == -1) {
          return false;
        } // Find insertion poPoly on B.


        var insb = -1;

        for (var _i10 = 0; _i10 < bcon.length; ++_i10) {
          if (bcon[_i10] == aid) {
            insb = _i10;
            break;
          }
        }

        if (insb == -1) {
          return false;
        } // Merge neighbours.


        rega.connections = [];

        for (var _i11 = 0, ni = acon.length; _i11 < ni - 1; ++_i11) {
          rega.connections.push(acon[(insa + 1 + _i11) % ni]);
        }

        for (var _i12 = 0, _ni3 = bcon.length; _i12 < _ni3 - 1; ++_i12) {
          rega.connections.push(bcon[(insb + 1 + _i12) % _ni3]);
        }

        RecastRegion.removeAdjacentNeighbours(rega);

        for (var j = 0; j < regb.floors.length; ++j) {
          RecastRegion.addUniqueFloorRegion(rega, regb.floors[j]);
        }

        rega.spanCount += regb.spanCount;
        regb.spanCount = 0;
        regb.connections = [];
        return true;
      }
    }, {
      key: "isRegionConnectedToBorder",
      value: function isRegionConnectedToBorder(reg) {
        // Region is connected to border if
        // one of the neighbours is null id.
        return reg.connections.includes(0);
      }
    }, {
      key: "isSolidEdge",
      value: function isSolidEdge(chf, srcReg, x, y, i, dir) {
        var s = chf.spans[i];
        var r = 0;

        if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
          var ax = x + RecastCommon.GetDirOffsetX(dir);
          var ay = y + RecastCommon.GetDirOffsetY(dir);
          var ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(s, dir);
          r = srcReg[ai];
        }

        if (r == srcReg[i]) {
          return false;
        }

        return true;
      }
    }, {
      key: "walkContour",
      value: function walkContour(x, y, i, dir, chf, srcReg, cont) {
        var startDir = dir;
        var starti = i;
        var ss = chf.spans[i];
        var curReg = 0;

        if (RecastCommon.GetCon(ss, dir) != RecastConstants.RC_NOT_CONNECTED) {
          var ax = x + RecastCommon.GetDirOffsetX(dir);
          var ay = y + RecastCommon.GetDirOffsetY(dir);
          var ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(ss, dir);
          curReg = srcReg[ai];
        }

        cont.push(curReg);
        var iter = 0;

        while (++iter < 40000) {
          var s = chf.spans[i];

          if (RecastRegion.isSolidEdge(chf, srcReg, x, y, i, dir)) {
            // Choose the edge corner
            var r = 0;

            if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
              var _ax6 = x + RecastCommon.GetDirOffsetX(dir);

              var _ay6 = y + RecastCommon.GetDirOffsetY(dir);

              var _ai6 = chf.cells[_ax6 + _ay6 * chf.width].index + RecastCommon.GetCon(s, dir);

              r = srcReg[_ai6];
            }

            if (r != curReg) {
              curReg = r;
              cont.push(curReg);
            }

            dir = dir + 1 & 0x3; // Rotate CW
          } else {
            var ni = -1;
            var nx = x + RecastCommon.GetDirOffsetX(dir);
            var ny = y + RecastCommon.GetDirOffsetY(dir);

            if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
              var nc = chf.cells[nx + ny * chf.width];
              ni = nc.index + RecastCommon.GetCon(s, dir);
            }

            if (ni == -1) {
              // Should not happen.
              return;
            }

            x = nx;
            y = ny;
            i = ni;
            dir = dir + 3 & 0x3; // Rotate CCW
          }

          if (starti == i && startDir == dir) {
            break;
          }
        } // Remove adjacent duplicates.


        if (cont.length > 1) {
          for (var j = 0; j < cont.length;) {
            var nj = (j + 1) % cont.length;

            if (cont[j] == cont[nj] && cont[j] >= -128 && cont[j] <= 127) {
              cont.splice(j, 1);
            } else {
              ++j;
            }
          }
        }
      }
    }, {
      key: "mergeAndFilterRegions",
      value: function mergeAndFilterRegions(ctx, minRegionArea, mergeRegionSize, maxRegionId, chf, srcReg, overlaps) {
        var w = chf.width;
        var h = chf.height;
        var nreg = maxRegionId + 1;
        var regions = new Array(nreg); // Construct regions

        for (var i = 0; i < nreg; ++i) {
          regions[i] = new RecastRegion.Region(i);
        } // Find edge of a region and find connections around the contour.


        for (var y = 0; y < h; ++y) {
          for (var _x8 = 0; _x8 < w; ++_x8) {
            var c = chf.cells[_x8 + y * w];

            for (var _i13 = c.index, ni = c.index + c.count; _i13 < ni; ++_i13) {
              var r = srcReg[_i13];

              if (r == 0 || r >= nreg) {
                continue;
              }

              var reg = regions[r];
              reg.spanCount++; // Update floors.

              for (var j = c.index; j < ni; ++j) {
                if (_i13 == j) {
                  continue;
                }

                var floorId = srcReg[j];

                if (floorId == 0 || floorId >= nreg) {
                  continue;
                }

                if (floorId == r) {
                  reg.overlap = true;
                }

                RecastRegion.addUniqueFloorRegion(reg, floorId);
              } // Have found contour


              if (reg.connections.length > 0) {
                continue;
              }

              reg.areaType = chf.areas[_i13]; // Check if this cell is next to a border.

              var ndir = -1;

              for (var dir = 0; dir < 4; ++dir) {
                if (RecastRegion.isSolidEdge(chf, srcReg, _x8, y, _i13, dir)) {
                  ndir = dir;
                  break;
                }
              }

              if (ndir != -1) {
                // The cell is at border.
                // Walk around the contour to find all the neighbours.
                RecastRegion.walkContour(_x8, y, _i13, ndir, chf, srcReg, reg.connections);
              }
            }
          }
        } // Remove too small regions.


        var stack = new Array(32);
        var trace = new Array(32);

        for (var _i14 = 0; _i14 < nreg; ++_i14) {
          var _reg = regions[_i14];

          if (_reg.id == 0 || (_reg.id & RecastConstants.RC_BORDER_REG) != 0) {
            continue;
          }

          if (_reg.spanCount == 0) {
            continue;
          }

          if (_reg.visited) {
            continue;
          } // Count the total size of all the connected regions.
          // Also keep track of the regions connects to a tile border.


          var connectsToBorder = false;
          var spanCount = 0; //stack = [];
          //trace = [];

          stack = [];
          trace = [];
          _reg.visited = true;
          stack.push(_i14);

          while (stack.length > 0) {
            // Pop
            var ri = stack.splice(stack.length - 1, 1);
            var creg = regions[ri];
            spanCount += creg.spanCount;
            trace.push(ri);

            for (var _j2 = 0; _j2 < creg.connections.length; ++_j2) {
              if ((creg.connections[_j2] & RecastConstants.RC_BORDER_REG) != 0) {
                connectsToBorder = true;
                continue;
              }

              var neireg = regions[creg.connections[_j2]];

              if (neireg.visited) {
                continue;
              }

              if (neireg.id == 0 || (neireg.id & RecastConstants.RC_BORDER_REG) != 0) {
                continue;
              } // Visit


              stack.push(neireg.id);
              neireg.visited = true;
            }
          } // If the accumulated regions size is too small, remove it.
          // Do not remove areas which connect to tile borders
          // as their size cannot be estimated correctly and removing them
          // can potentially remove necessary areas.


          if (spanCount < minRegionArea && !connectsToBorder) {
            // Kill all visited regions.
            for (var _j3 = 0; _j3 < trace.length; ++_j3) {
              regions[trace[_j3]].spanCount = 0;
              regions[trace[_j3]].id = 0;
            }
          }
        } // Merge too small regions to neighbour regions.


        var mergeCount = 0;

        do {
          mergeCount = 0;

          for (var _i15 = 0; _i15 < nreg; ++_i15) {
            var _reg2 = regions[_i15];

            if (_reg2.id == 0 || (_reg2.id & RecastConstants.RC_BORDER_REG) != 0) {
              continue;
            }

            if (_reg2.overlap) {
              continue;
            }

            if (_reg2.spanCount == 0) {
              continue;
            } // Check to see if the region should be merged.


            if (_reg2.spanCount > mergeRegionSize && RecastRegion.isRegionConnectedToBorder(_reg2)) {
              continue;
            } // Small region with more than 1 connection.
            // Or region which is not connected to a border at all.
            // Find smallest neighbour region that connects to this one.


            var smallest = 0xfffffff;
            var mergeId = _reg2.id;

            for (var _j4 = 0; _j4 < _reg2.connections.length; ++_j4) {
              if ((_reg2.connections[_j4] & RecastConstants.RC_BORDER_REG) != 0) {
                continue;
              }

              var mreg = regions[_reg2.connections[_j4]];

              if (mreg.id == 0 || (mreg.id & RecastConstants.RC_BORDER_REG) != 0 || mreg.overlap) {
                continue;
              }

              if (mreg.spanCount < smallest && RecastRegion.canMergeWithRegion(_reg2, mreg) && RecastRegion.canMergeWithRegion(mreg, _reg2)) {
                smallest = mreg.spanCount;
                mergeId = mreg.id;
              }
            } // Found new id.


            if (mergeId != _reg2.id) {
              var oldId = _reg2.id;
              var target = regions[mergeId]; // Merge neighbours.

              if (RecastRegion.mergeRegions(target, _reg2)) {
                // Fixup regions pointing to current region.
                for (var _j5 = 0; _j5 < nreg; ++_j5) {
                  if (regions[_j5].id == 0 || (regions[_j5].id & RecastConstants.RC_BORDER_REG) != 0) {
                    continue;
                  } // If another region was already merged into current region
                  // change the nid of the previous region too.


                  if (regions[_j5].id == oldId) {
                    regions[_j5].id = mergeId;
                  } // Replace the current region with the new one if the
                  // current regions is neighbour.


                  RecastRegion.replaceNeighbour(regions[_j5], oldId, mergeId);
                }

                mergeCount++;
              }
            }
          }
        } while (mergeCount > 0); // Compress region Ids.


        for (var _i16 = 0; _i16 < nreg; ++_i16) {
          regions[_i16].remap = false;

          if (regions[_i16].id == 0) {
            continue; // Skip nil regions.
          }

          if ((regions[_i16].id & RecastConstants.RC_BORDER_REG) != 0) {
            continue; // Skip external regions.
          }

          regions[_i16].remap = true;
        }

        var regIdGen = 0;

        for (var _i17 = 0; _i17 < nreg; ++_i17) {
          if (!regions[_i17].remap) {
            continue;
          }

          var _oldId = regions[_i17].id;
          var newId = ++regIdGen;

          for (var _j6 = _i17; _j6 < nreg; ++_j6) {
            if (regions[_j6].id == _oldId) {
              regions[_j6].id = newId;
              regions[_j6].remap = false;
            }
          }
        }

        maxRegionId = regIdGen; // Remap regions.

        for (var _i18 = 0; _i18 < chf.spanCount; ++_i18) {
          if ((srcReg[_i18] & RecastConstants.RC_BORDER_REG) == 0) {
            srcReg[_i18] = regions[srcReg[_i18]].id;
          }
        } // Return regions that we found to be overlapping.


        for (var _i19 = 0; _i19 < nreg; ++_i19) {
          if (regions[_i19].overlap) {
            overlaps.push(regions[_i19].id);
          }
        }

        return maxRegionId;
      }
    }, {
      key: "addUniqueConnection",
      value: function addUniqueConnection(reg, n) {
        if (!reg.connections.contains(n)) {
          reg.connections.push(n);
        }
      }
    }, {
      key: "mergeAndFilterLayerRegions",
      value: function mergeAndFilterLayerRegions(ctx, minRegionArea, maxRegionId, chf, srcReg, overlaps) {
        var w = chf.width;
        var h = chf.height;
        var nreg = maxRegionId + 1;
        var regions = new Array(nreg); // Construct regions

        for (var i = 0; i < nreg; ++i) {
          regions[i] = new Region(i);
        } // Find region neighbours and overlapping regions.


        var lregs = new Array(32);

        for (var y = 0; y < h; ++y) {
          for (var _x9 = 0; _x9 < w; ++_x9) {
            var c = chf.cells[_x9 + y * w];
            lregs = [];

            for (var _i20 = c.index, ni = c.index + c.count; _i20 < ni; ++_i20) {
              var s = chf.spans[_i20];
              var ri = srcReg[_i20];

              if (ri == 0 || ri >= nreg) {
                continue;
              }

              var reg = regions[ri];
              reg.spanCount++;
              reg.ymin = Math.min(reg.ymin, s.y);
              reg.ymax = Math.max(reg.ymax, s.y); // Collect all region layers.

              lregs.push(ri); // Update neighbours

              for (var dir = 0; dir < 4; ++dir) {
                if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
                  var ax = _x9 + RecastCommon.GetDirOffsetX(dir);

                  var ay = y + RecastCommon.GetDirOffsetY(dir);
                  var ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);
                  var rai = srcReg[ai];

                  if (rai > 0 && rai < nreg && rai != ri) {
                    addUniqueConnection(reg, rai);
                  }

                  if ((rai & RecastConstants.RC_BORDER_REG) != 0) {
                    reg.connectsToBorder = true;
                  }
                }
              }
            } // Update overlapping regions.


            for (var _i21 = 0; _i21 < lregs.length - 1; ++_i21) {
              for (var j = _i21 + 1; j < lregs.length; ++j) {
                if (lregs[_i21] != lregs[j]) {
                  var _ri = regions[lregs[_i21]];
                  var rj = regions[lregs[j]];
                  RecastRegion.addUniqueFloorRegion(_ri, lregs[j]);
                  RecastRegion.addUniqueFloorRegion(rj, lregs[_i21]);
                }
              }
            }
          }
        } // Create 2D layers from regions.


        var layerId = 1;

        for (var _i22 = 0; _i22 < nreg; ++_i22) {
          regions[_i22].id = 0;
        } // Merge montone regions to create non-overlapping areas.


        var stack = new Array(32);

        for (var _i23 = 1; _i23 < nreg; ++_i23) {
          var root = regions[_i23]; // Skip already visited.

          if (root.id != 0) {
            continue;
          } // Start search.


          root.id = layerId;
          stack = [];
          stack.push(_i23);

          while (stack.length > 0) {
            // Pop front
            var _reg3 = regions[stack.remove(0)];
            var ncons = _reg3.connections.length;

            for (var _j7 = 0; _j7 < ncons; ++_j7) {
              var nei = _reg3.connections[_j7];
              var regn = regions[nei]; // Skip already visited.

              if (regn.id != 0) {
                continue;
              } // Skip if the neighbour is overlapping root region.


              var overlap = false;

              for (var k = 0; k < root.floors.length; k++) {
                if (root.floors[k] == nei) {
                  overlap = true;
                  break;
                }
              }

              if (overlap) {
                continue;
              } // Deepen


              stack.push(nei); // Mark layer id

              regn.id = layerId; // Merge current layers to root.

              for (var _k = 0; _k < regn.floors.length; ++_k) {
                RecastRegion.addUniqueFloorRegion(root, regn.floors[_k]);
              }

              root.ymin = Math.min(root.ymin, regn.ymin);
              root.ymax = Math.max(root.ymax, regn.ymax);
              root.spanCount += regn.spanCount;
              regn.spanCount = 0;
              root.connectsToBorder = root.connectsToBorder || regn.connectsToBorder;
            }
          }

          layerId++;
        } // Remove small regions


        for (var _i24 = 0; _i24 < nreg; ++_i24) {
          if (regions[_i24].spanCount > 0 && regions[_i24].spanCount < minRegionArea && !regions[_i24].connectsToBorder) {
            var _reg4 = regions[_i24].id;

            for (var _j8 = 0; _j8 < nreg; ++_j8) {
              if (regions[_j8].id == _reg4) {
                regions[_j8].id = 0;
              }
            }
          }
        } // Compress region Ids.


        for (var _i25 = 0; _i25 < nreg; ++_i25) {
          regions[_i25].remap = false;

          if (regions[_i25].id == 0) {
            continue; // Skip nil regions.
          }

          if ((regions[_i25].id & RecastConstants.RC_BORDER_REG) != 0) {
            continue; // Skip external regions.
          }

          regions[_i25].remap = true;
        }

        var regIdGen = 0;

        for (var _i26 = 0; _i26 < nreg; ++_i26) {
          if (!regions[_i26].remap) {
            continue;
          }

          var oldId = regions[_i26].id;
          var newId = ++regIdGen;

          for (var _j9 = _i26; _j9 < nreg; ++_j9) {
            if (regions[_j9].id == oldId) {
              regions[_j9].id = newId;
              regions[_j9].remap = false;
            }
          }
        }

        maxRegionId = regIdGen; // Remap regions.

        for (var _i27 = 0; _i27 < chf.spanCount; ++_i27) {
          if ((srcReg[_i27] & RecastConstants.RC_BORDER_REG) == 0) {
            srcReg[_i27] = regions[srcReg[_i27]].id;
          }
        }

        return maxRegionId;
      } /// @par
      ///
      /// This is usually the second to the last step in creating a fully built
      /// compact heightfield. This step is required before regions are built
      /// using #rcBuildRegions or #rcBuildRegionsMonotone.
      ///
      /// After this step, the distance data is available via the rcCompactHeightfield::maxDistance
      /// and rcCompactHeightfield::dist fields.
      ///
      /// @see rcCompactHeightfield, rcBuildRegions, rcBuildRegionsMonotone

    }, {
      key: "buildDistanceField",
      value: function buildDistanceField(ctx, chf) {
        ctx.startTimer("BUILD_DISTANCEFIELD");
        var src = new Array(chf.spanCount);
        ctx.startTimer("DISTANCEFIELD_DIST");
        var maxDist = this.calculateDistanceField(chf, src);
        chf.maxDistance = maxDist;
        ctx.stopTimer("DISTANCEFIELD_DIST");
        ctx.startTimer("DISTANCEFIELD_BLUR"); // Blur

        src = this.boxBlur(chf, 1, src); // Store distance.

        chf.dist = src;
        ctx.stopTimer("DISTANCEFIELD_BLUR");
        ctx.stopTimer("BUILD_DISTANCEFIELD");
      }
    }, {
      key: "paintRectRegion",
      value: function paintRectRegion(minx, maxx, miny, maxy, regId, chf, srcReg) {
        var w = chf.width;

        for (var y = miny; y < maxy; ++y) {
          for (var _y4 = minx; x < maxx; ++x) {
            var c = chf.cells[x + _y4 * w];

            for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
              if (chf.areas[i] != RecastConstants.RC_NULL_AREA) {
                srcReg[i] = regId;
              }
            }
          }
        }
      } /// @par
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

    }, {
      key: "buildRegionsMonotone",
      value: function buildRegionsMonotone(ctx, chf, borderSize, minRegionArea, mergeRegionArea) {
        ctx.startTimer("BUILD_REGIONS");
        var w = chf.width;
        var h = chf.height;
        var id = 1;
        var srcReg = new Array(chf.spanCount);
        var nsweeps = Math.max(chf.width, chf.height);
        sweeps = new Array(nsweeps);

        for (var i = 0; i < sweeps.length; i++) {
          sweeps[i] = new SweepSpan();
        } // Mark border regions.


        if (borderSize > 0) {
          // Make sure border will not overflow.
          var bw = Math.min(w, borderSize);
          var bh = Math.min(h, borderSize); // PaPoly regions

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
        var prev = new Array(256); // Sweep one line at a time.

        for (var y = borderSize; y < h - borderSize; ++y) {
          // Collect spans from this row.
          prev.fill(0, 0, id);
          var rid = 1;

          for (var _y5 = borderSize; x < w - borderSize; ++x) {
            var c = chf.cells[x + _y5 * w];

            for (var _i28 = c.index, ni = c.index + c.count; _i28 < ni; ++_i28) {
              var s = chf.spans[_i28];

              if (chf.areas[_i28] == RecastConstants.RC_NULL_AREA) {
                continue;
              } // -x


              var previd = 0;

              if (RecastCommon.GetCon(s, 0) != RecastConstants.RC_NOT_CONNECTED) {
                var ax = x + RecastCommon.GetDirOffsetX(0);

                var ay = _y5 + RecastCommon.GetDirOffsetY(0);

                var ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 0);

                if ((srcReg[ai] & RecastConstants.RC_BORDER_REG) == 0 && chf.areas[_i28] == chf.areas[ai]) {
                  previd = srcReg[ai];
                }
              }

              if (previd == 0) {
                previd = rid++;
                sweeps[previd].rid = previd;
                sweeps[previd].ns = 0;
                sweeps[previd].nei = 0;
              } // -y


              if (RecastCommon.GetCon(s, 3) != RecastConstants.RC_NOT_CONNECTED) {
                var _ax7 = x + RecastCommon.GetDirOffsetX(3);

                var _ay7 = _y5 + RecastCommon.GetDirOffsetY(3);

                var _ai7 = chf.cells[_ax7 + _ay7 * w].index + RecastCommon.GetCon(s, 3);

                if (srcReg[_ai7] != 0 && (srcReg[_ai7] & RecastConstants.RC_BORDER_REG) == 0 && chf.areas[_i28] == chf.areas[_ai7]) {
                  var nr = srcReg[_ai7];

                  if (sweeps[previd].nei == 0 || sweeps[previd].nei == nr) {
                    sweeps[previd].nei = nr;
                    sweeps[previd].ns++;
                    prev[nr]++;
                  } else {
                    sweeps[previd].nei = RC_NULL_NEI;
                  }
                }
              }

              srcReg[_i28] = previd;
            }
          } // Create unique ID.


          for (var _i29 = 1; _i29 < rid; ++_i29) {
            if (sweeps[_i29].nei != RC_NULL_NEI && sweeps[_i29].nei != 0 && prev[sweeps[_i29].nei] == sweeps[_i29].ns) {
              sweeps[_i29].id = sweeps[_i29].nei;
            } else {
              sweeps[_i29].id = id++;
            }
          } // Remap IDs


          for (var _y6 = borderSize; x < w - borderSize; ++x) {
            var _c3 = chf.cells[x + _y6 * w];

            for (var _i30 = _c3.index, _ni4 = _c3.index + _c3.count; _i30 < _ni4; ++_i30) {
              if (srcReg[_i30] > 0 && srcReg[_i30] < rid) {
                srcReg[_i30] = sweeps[srcReg[_i30]].id;
              }
            }
          }
        }

        ctx.startTimer("BUILD_REGIONS_FILTER"); // Merge regions and filter out small regions.

        var overlaps = [];
        chf.maxRegions = mergeAndFilterRegions(ctx, minRegionArea, mergeRegionArea, id, chf, srcReg, overlaps); // Monotone partitioning does not generate overlapping regions.

        ctx.stopTimer("BUILD_REGIONS_FILTER"); // Store the result out.

        for (var _i31 = 0; _i31 < chf.spanCount; ++_i31) {
          // if (i == 1344)
          //     console.log("4431")
          chf.spans[_i31].reg = srcReg[_i31];
        }

        ctx.stopTimer("BUILD_REGIONS");
      } /// @par
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

    }, {
      key: "buildRegions",
      value: function buildRegions(ctx, chf, borderSize, minRegionArea, mergeRegionArea) {
        ctx.startTimer("BUILD_REGIONS");
        var w = chf.width;
        var h = chf.height;
        ctx.startTimer("REGIONS_WATERSHED");
        var LOG_NB_STACKS = 3;
        var NB_STACKS = 1 << LOG_NB_STACKS;
        var lvlStacks = [];

        for (var i = 0; i < NB_STACKS; ++i) {
          lvlStacks.push([]);
        }

        var stack = new Array(1024);
        var srcReg = new Array(chf.spanCount).fill(0);
        var srcDist = new Array(chf.spanCount).fill(0);
        var regionId = 1;
        var level = chf.maxDistance + 1 & ~1; // TODO: Figure better formula, expandIters defines how much the
        // watershed "overflows" and simplifies the regions. Tying it to
        // agent radius was usually good indication how greedy it could be.
        // const let expandIters = 4 + walkableRadius * 2;

        var expandIters = 8;

        if (borderSize > 0) {
          // Make sure border will not overflow.
          var bw = Math.min(w, borderSize);
          var bh = Math.min(h, borderSize); // PaPoly regions

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
        var sId = -1;

        while (level > 0) {
          level = level >= 2 ? level - 2 : 0;
          sId = sId + 1 & NB_STACKS - 1; // ctx=>startTimer(RC_TIMER_DIVIDE_TO_LEVELS);

          if (sId == 0) {
            RecastRegion.sortCellsByLevel(level, chf, srcReg, NB_STACKS, lvlStacks, 1);
          } else {
            RecastRegion.appendStacks(lvlStacks[sId - 1], lvlStacks[sId], srcReg); // copy left overs from last level
          } // ctx=>stopTimer(RC_TIMER_DIVIDE_TO_LEVELS);


          ctx.startTimer("BUILD_REGIONS_EXPAND"); // Expand current regions until no empty connected cells found.

          RecastRegion.expandRegions(expandIters, level, chf, srcReg, srcDist, lvlStacks[sId], false);
          ctx.stopTimer("BUILD_REGIONS_EXPAND");
          ctx.startTimer("BUILD_REGIONS_FLOOD"); // Mark new regions with IDs.

          for (var j = 0; j < lvlStacks[sId].length; j += 3) {
            var _x10 = lvlStacks[sId][j];
            var y = lvlStacks[sId][j + 1];
            var _i32 = lvlStacks[sId][j + 2];

            if (_i32 >= 0 && srcReg[_i32] == 0) {
              if (RecastRegion.floodRegion(_x10, y, _i32, level, regionId, chf, srcReg, srcDist, stack)) {
                regionId++;
              }
            }
          }

          ctx.stopTimer("BUILD_REGIONS_FLOOD");
        } // Expand current regions until no empty connected cells found.


        RecastRegion.expandRegions(expandIters * 8, 0, chf, srcReg, srcDist, stack, true);
        ctx.stopTimer("BUILD_REGIONS_WATERSHED");
        ctx.startTimer("BUILD_REGIONS_FILTER"); // Merge regions and filter out smalle regions.

        var overlaps = [];
        chf.maxRegions = this.mergeAndFilterRegions(ctx, minRegionArea, mergeRegionArea, regionId, chf, srcReg, overlaps); // If overlapping regions were found during merging, split those regions.

        if (overlaps.length > 0) {
          throw new RuntimeException("rcBuildRegions: " + overlaps.length + " overlapping regions.");
        }

        ctx.stopTimer("BUILD_REGIONS_FILTER"); // Write the result out.

        for (var _i33 = 0; _i33 < chf.spanCount; ++_i33) {
          // if (i == 1344)
          //     console.log("4431")
          chf.spans[_i33].reg = srcReg[_i33];
        }

        ctx.stopTimer("BUILD_REGIONS");
      }
    }, {
      key: "buildLayerRegions",
      value: function buildLayerRegions(ctx, chf, borderSize, minRegionArea) {
        ctx.startTimer("BUILD_REGIONS");
        var w = chf.width;
        var h = chf.height;
        var id = 1;
        var srcReg = new Array(chf.spanCount);
        var nsweeps = Math.max(chf.width, chf.height);
        var sweeps = new ArrayList(nsweeps);

        for (var i = 0; i < sweeps.length; i++) {
          sweeps[i] = new SweepSpan();
        } // Mark border regions.


        if (borderSize > 0) {
          // Make sure border will not overflow.
          var bw = Math.min(w, borderSize);
          var bh = Math.min(h, borderSize); // PaPoly regions

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
        var prev = new Array(256); // Sweep one line at a time.

        for (var y = borderSize; y < h - borderSize; ++y) {
          // Collect spans from this row.
          prev.fill(0, 0, id);
          var rid = 1;

          for (var _y7 = borderSize; x < w - borderSize; ++x) {
            var c = chf.cells[x + _y7 * w];

            for (var _i34 = c.index, ni = c.index + c.count; _i34 < ni; ++_i34) {
              var s = chf.spans[_i34];

              if (chf.areas[_i34] == RecastConstants.RC_NULL_AREA) {
                continue;
              } // -x


              var previd = 0;

              if (RecastCommon.GetCon(s, 0) != RecastConstants.RC_NOT_CONNECTED) {
                var ax = x + RecastCommon.GetDirOffsetX(0);

                var ay = _y7 + RecastCommon.GetDirOffsetY(0);

                var ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 0);

                if ((srcReg[ai] & RecastConstants.RC_BORDER_REG) == 0 && chf.areas[_i34] == chf.areas[ai]) {
                  previd = srcReg[ai];
                }
              }

              if (previd == 0) {
                previd = rid++;
                sweeps[previd].rid = previd;
                sweeps[previd].ns = 0;
                sweeps[previd].nei = 0;
              } // -y


              if (RecastCommon.GetCon(s, 3) != RecastConstants.RC_NOT_CONNECTED) {
                var _ax8 = x + RecastCommon.GetDirOffsetX(3);

                var _ay8 = _y7 + RecastCommon.GetDirOffsetY(3);

                var _ai8 = chf.cells[_ax8 + _ay8 * w].index + RecastCommon.GetCon(s, 3);

                if (srcReg[_ai8] != 0 && (srcReg[_ai8] & RecastConstants.RC_BORDER_REG) == 0 && chf.areas[_i34] == chf.areas[_ai8]) {
                  var nr = srcReg[_ai8];

                  if (sweeps[previd].nei == 0 || sweeps[previd].nei == nr) {
                    sweeps[previd].nei = nr;
                    sweeps[previd].ns++;
                    prev[nr]++;
                  } else {
                    sweeps[previd].nei = RC_NULL_NEI;
                  }
                }
              }

              srcReg[_i34] = previd;
            }
          } // Create unique ID.


          for (var _i35 = 1; _i35 < rid; ++_i35) {
            if (sweeps[_i35].nei != RC_NULL_NEI && sweeps[_i35].nei != 0 && prev[sweeps[_i35].nei] == sweeps[_i35].ns) {
              sweeps[_i35].id = sweeps[_i35].nei;
            } else {
              sweeps[_i35].id = id++;
            }
          } // Remap IDs


          for (var _y8 = borderSize; x < w - borderSize; ++x) {
            var _c4 = chf.cells[x + _y8 * w];

            for (var _i36 = _c4.index, _ni5 = _c4.index + _c4.count; _i36 < _ni5; ++_i36) {
              if (srcReg[_i36] > 0 && srcReg[_i36] < rid) {
                srcReg[_i36] = sweeps[srcReg[_i36]].id;
              }
            }
          }
        }

        ctx.startTimer("BUILD_REGIONS_FILTER"); // Merge monotone regions to layers and remove small regions.

        var overlaps = [];
        chf.maxRegions = mergeAndFilterLayerRegions(ctx, minRegionArea, id, chf, srcReg, overlaps);
        ctx.stopTimer("BUILD_REGIONS_FILTER"); // Store the result out.

        for (var _i37 = 0; _i37 < chf.spanCount; ++_i37) {
          chf.spans[_i37].reg = srcReg[_i37];
        }

        ctx.stopTimer("BUILD_REGIONS");
      }
    }]);

    return RecastRegion;
  }();

  _defineProperty(RecastRegion, "RC_NULL_NEI", 0xfff);

  _defineProperty(RecastRegion, "SweepSpan", (_temp$1 = function SweepSpan() {
    _classCallCheck(this, SweepSpan);

    _defineProperty(this, "rid", void 0);

    _defineProperty(this, "id", void 0);

    _defineProperty(this, "ns", void 0);

    _defineProperty(this, "nei", void 0);
  } // neighbour id
  , _temp$1));

  _defineProperty(RecastRegion, "Region", (_temp2$1 = // Number of spans belonging to this region
  // ID of the region
  // Are type.
  function Region(i) {
    _classCallCheck(this, Region);

    _defineProperty(this, "spanCount", 0);

    _defineProperty(this, "id", 0);

    _defineProperty(this, "areaType", 0);

    _defineProperty(this, "remap", false);

    _defineProperty(this, "visited", false);

    _defineProperty(this, "overlap", false);

    _defineProperty(this, "connectsToBorder", false);

    _defineProperty(this, "ymin", 0);

    _defineProperty(this, "ymax", 0);

    _defineProperty(this, "connections", []);

    _defineProperty(this, "floors", []);

    this.id = i;
    this.ymin = 0xFFFF;
    this.connections = [];
    this.floors = [];
  }, _temp2$1));

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
  var ContourSet = function ContourSet() {
    _classCallCheck(this, ContourSet);

    _defineProperty(this, "conts", []);

    _defineProperty(this, "bmin", new Array(3));

    _defineProperty(this, "bmax", new Array(3));

    _defineProperty(this, "cs", 0);

    _defineProperty(this, "ch", 0);

    _defineProperty(this, "width", 0);

    _defineProperty(this, "height", 0);

    _defineProperty(this, "borderSize", 0);

    _defineProperty(this, "maxError", 0);
  };

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
  var PolyMesh = function PolyMesh() {
    _classCallCheck(this, PolyMesh);

    _defineProperty(this, "verts", []);

    _defineProperty(this, "polys", []);

    _defineProperty(this, "regs", []);

    _defineProperty(this, "areas", []);

    _defineProperty(this, "nverts", 0);

    _defineProperty(this, "npolys", 0);

    _defineProperty(this, "nvp", 0);

    _defineProperty(this, "maxpolys", 0);

    _defineProperty(this, "flags", []);

    _defineProperty(this, "bmin", new Array(3));

    _defineProperty(this, "bmax", new Array(3));

    _defineProperty(this, "cs", 0);

    _defineProperty(this, "ch", 0);

    _defineProperty(this, "borderSize", 0);

    _defineProperty(this, "maxEdgeError", 0);
  };

  function arraycopy$1(one, oneStart, two, twoStart, len) {
    for (var i = 0; i < len; i++) {
      two[twoStart + i] = one[oneStart + i];
    }
  }

  var Edge = function Edge() {
    _classCallCheck(this, Edge);

    _defineProperty(this, "vert", new Array(2));

    _defineProperty(this, "polyEdge", new Array(2));

    _defineProperty(this, "poly", new Array(2));
  };

  var RecastMesh = /*#__PURE__*/function () {
    function RecastMesh() {
      _classCallCheck(this, RecastMesh);
    }

    _createClass(RecastMesh, null, [{
      key: "buildMeshAdjacency",
      value: function buildMeshAdjacency(polys, npolys, nverts, vertsPerPoly) {
        // Based on code by Eric Lengyel from:
        // http://www.terathon.com/code/edges.php
        var maxEdgeCount = npolys * vertsPerPoly;
        var firstEdge = new Array(nverts + maxEdgeCount);
        var nextEdge = nverts;
        var edgeCount = 0;
        var edges = new Array(maxEdgeCount);

        for (var i = 0; i < nverts; i++) {
          firstEdge[i] = RecastConstants.RC_MESH_NULL_IDX;
        }

        for (var _i = 0; _i < npolys; ++_i) {
          var t = _i * vertsPerPoly * 2;

          for (var j = 0; j < vertsPerPoly; ++j) {
            if (polys[t + j] == RecastConstants.RC_MESH_NULL_IDX) break;
            var v0 = polys[t + j];
            var v1 = j + 1 >= vertsPerPoly || polys[t + j + 1] == RecastConstants.RC_MESH_NULL_IDX ? polys[t + 0] : polys[t + j + 1];

            if (v0 < v1) {
              var edge = new Edge();
              edges[edgeCount] = edge;
              edge.vert[0] = v0;
              edge.vert[1] = v1;
              edge.poly[0] = _i;
              edge.polyEdge[0] = j;
              edge.poly[1] = _i;
              edge.polyEdge[1] = 0; // Insert edge

              firstEdge[nextEdge + edgeCount] = firstEdge[v0];
              firstEdge[v0] = edgeCount;
              edgeCount++;
            }
          }
        }

        for (var _i2 = 0; _i2 < npolys; ++_i2) {
          var _t = _i2 * vertsPerPoly * 2;

          for (var _j = 0; _j < vertsPerPoly; ++_j) {
            if (polys[_t + _j] == RecastConstants.RC_MESH_NULL_IDX) break;
            var _v = polys[_t + _j];

            var _v2 = _j + 1 >= vertsPerPoly || polys[_t + _j + 1] == RecastConstants.RC_MESH_NULL_IDX ? polys[_t + 0] : polys[_t + _j + 1];

            if (_v > _v2) {
              for (var e = firstEdge[_v2]; e != RecastConstants.RC_MESH_NULL_IDX; e = firstEdge[nextEdge + e]) {
                var _edge = edges[e];

                if (_edge.vert[1] == _v && _edge.poly[0] == _edge.poly[1]) {
                  _edge.poly[1] = _i2;
                  _edge.polyEdge[1] = _j;
                  break;
                }
              }
            }
          }
        } // Store adjacency


        for (var _i3 = 0; _i3 < edgeCount; ++_i3) {
          var _e = edges[_i3];

          if (_e.poly[0] != _e.poly[1]) {
            var p0 = _e.poly[0] * vertsPerPoly * 2;
            var p1 = _e.poly[1] * vertsPerPoly * 2;
            polys[p0 + vertsPerPoly + _e.polyEdge[0]] = _e.poly[1];
            polys[p1 + vertsPerPoly + _e.polyEdge[1]] = _e.poly[0];
          }
        }
      }
    }, {
      key: "computeVertexHash",
      value: function computeVertexHash(x, y, z) {
        var h1 = 0x8da6b343; // Large multiplicative constants;

        var h2 = 0xd8163841; // here arbitrarily chosen primes

        var h3 = 0xcb1ab31;
        var n = h1 * x + h2 * y + h3 * z;
        return n & RecastMesh.VERTEX_BUCKET_COUNT - 1;
      }
    }, {
      key: "addVertex",
      value: function addVertex(x, y, z, verts, firstVert, nextVert, nv) {
        var bucket = RecastMesh.computeVertexHash(x, 0, z);
        var i = firstVert[bucket];

        while (i != -1) {
          var _v3 = i * 3;

          if (verts[_v3 + 0] == x && Math.abs(verts[_v3 + 1] - y) <= 2 && verts[_v3 + 2] == z) return [i, nv];
          i = nextVert[i]; // next
        } // Could not find, create new.


        i = nv;
        nv++;
        var v = i * 3;
        verts[v + 0] = x;
        verts[v + 1] = y;
        verts[v + 2] = z;
        nextVert[i] = firstVert[bucket];
        firstVert[bucket] = i;
        return [i, nv];
      }
    }, {
      key: "prev",
      value: function prev(i, n) {
        return i - 1 >= 0 ? i - 1 : n - 1;
      }
    }, {
      key: "next",
      value: function next(i, n) {
        return i + 1 < n ? i + 1 : 0;
      }
    }, {
      key: "area2",
      value: function area2(verts, a, b, c) {
        return (verts[b + 0] - verts[a + 0]) * (verts[c + 2] - verts[a + 2]) - (verts[c + 0] - verts[a + 0]) * (verts[b + 2] - verts[a + 2]);
      } // Returns true iff c is strictly to the left of the directed
      // line through a to b.

    }, {
      key: "left",
      value: function left(verts, a, b, c) {
        return RecastMesh.area2(verts, a, b, c) < 0;
      }
    }, {
      key: "leftOn",
      value: function leftOn(verts, a, b, c) {
        return RecastMesh.area2(verts, a, b, c) <= 0;
      }
    }, {
      key: "collinear",
      value: function collinear(verts, a, b, c) {
        return RecastMesh.area2(verts, a, b, c) == 0;
      } // Returns true iff ab properly intersects cd: they share
      // a poPoly interior to both segments. The properness of the
      // intersection is ensured by using strict leftness.

    }, {
      key: "intersectProp",
      value: function intersectProp(verts, a, b, c, d) {
        // Eliminate improper cases.
        if (RecastMesh.collinear(verts, a, b, c) || RecastMesh.collinear(verts, a, b, d) || RecastMesh.collinear(verts, c, d, a) || RecastMesh.collinear(verts, c, d, b)) return false;
        return RecastMesh.left(verts, a, b, c) ^ RecastMesh.left(verts, a, b, d) && RecastMesh.left(verts, c, d, a) ^ RecastMesh.left(verts, c, d, b);
      } // Returns T iff (a,b,c) are collinear and poPoly c lies
      // on the closed segement ab.

    }, {
      key: "between",
      value: function between(verts, a, b, c) {
        if (!RecastMesh.collinear(verts, a, b, c)) return false; // If ab not vertical, check betweenness on x; else on y.

        if (verts[a + 0] != verts[b + 0]) return verts[a + 0] <= verts[c + 0] && verts[c + 0] <= verts[b + 0] || verts[a + 0] >= verts[c + 0] && verts[c + 0] >= verts[b + 0];else return verts[a + 2] <= verts[c + 2] && verts[c + 2] <= verts[b + 2] || verts[a + 2] >= verts[c + 2] && verts[c + 2] >= verts[b + 2];
      } // Returns true iff segments ab and cd intersect, properly or improperly.

    }, {
      key: "intersect",
      value: function intersect(verts, a, b, c, d) {
        if (RecastMesh.intersectProp(verts, a, b, c, d)) return true;else if (RecastMesh.between(verts, a, b, c) || RecastMesh.between(verts, a, b, d) || RecastMesh.between(verts, c, d, a) || RecastMesh.between(verts, c, d, b)) return true;else return false;
      }
    }, {
      key: "vequal",
      value: function vequal(verts, a, b) {
        return verts[a + 0] == verts[b + 0] && verts[a + 2] == verts[b + 2];
      } // Returns T iff (v_i, v_j) is a proper internal *or* external
      // diagonal of P, *ignoring edges incident to v_i and v_j*.

    }, {
      key: "diagonalie",
      value: function diagonalie(i, j, n, verts, indices) {
        var d0 = (indices[i] & 0x0fffffff) * 4;
        var d1 = (indices[j] & 0x0fffffff) * 4; // For each edge (k,k+1) of P

        for (var k = 0; k < n; k++) {
          var k1 = RecastMesh.next(k, n); // Skip edges incident to i or j

          if (!(k == i || k1 == i || k == j || k1 == j)) {
            var p0 = (indices[k] & 0x0fffffff) * 4;
            var p1 = (indices[k1] & 0x0fffffff) * 4;
            if (RecastMesh.vequal(verts, d0, p0) || RecastMesh.vequal(verts, d1, p0) || RecastMesh.vequal(verts, d0, p1) || RecastMesh.vequal(verts, d1, p1)) continue;
            if (RecastMesh.intersect(verts, d0, d1, p0, p1)) return false;
          }
        }

        return true;
      } // Returns true iff the diagonal (i,j) is strictly internal to the
      // polygon P in the neighborhood of the i endpoint.

    }, {
      key: "inCone",
      value: function inCone(i, j, n, verts, indices) {
        var pi = (indices[i] & 0x0fffffff) * 4;
        var pj = (indices[j] & 0x0fffffff) * 4;
        var pi1 = (indices[RecastMesh.next(i, n)] & 0x0fffffff) * 4;
        var pin1 = (indices[RecastMesh.prev(i, n)] & 0x0fffffff) * 4; // If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].

        if (RecastMesh.leftOn(verts, pin1, pi, pi1)) {
          return RecastMesh.left(verts, pi, pj, pin1) && RecastMesh.left(verts, pj, pi, pi1);
        } // Assume (i-1,i,i+1) not collinear.
        // else P[i] is reflex.


        return !(RecastMesh.leftOn(verts, pi, pj, pi1) && RecastMesh.leftOn(verts, pj, pi, pin1));
      } // Returns T iff (v_i, v_j) is a proper internal
      // diagonal of P.

    }, {
      key: "diagonal",
      value: function diagonal(i, j, n, verts, indices) {
        return RecastMesh.inCone(i, j, n, verts, indices) && RecastMesh.diagonalie(i, j, n, verts, indices);
      }
    }, {
      key: "diagonalieLoose",
      value: function diagonalieLoose(i, j, n, verts, indices) {
        var d0 = (indices[i] & 0x0fffffff) * 4;
        var d1 = (indices[j] & 0x0fffffff) * 4; // For each edge (k,k+1) of P

        for (var k = 0; k < n; k++) {
          var k1 = RecastMesh.next(k, n); // Skip edges incident to i or j

          if (!(k == i || k1 == i || k == j || k1 == j)) {
            var p0 = (indices[k] & 0x0fffffff) * 4;
            var p1 = (indices[k1] & 0x0fffffff) * 4;
            if (RecastMesh.vequal(verts, d0, p0) || RecastMesh.vequal(verts, d1, p0) || RecastMesh.vequal(verts, d0, p1) || RecastMesh.vequal(verts, d1, p1)) continue;
            if (RecastMesh.intersectProp(verts, d0, d1, p0, p1)) return false;
          }
        }

        return true;
      }
    }, {
      key: "inConeLoose",
      value: function inConeLoose(i, j, n, verts, indices) {
        var pi = (indices[i] & 0x0fffffff) * 4;
        var pj = (indices[j] & 0x0fffffff) * 4;
        var pi1 = (indices[RecastMesh.next(i, n)] & 0x0fffffff) * 4;
        var pin1 = (indices[RecastMesh.prev(i, n)] & 0x0fffffff) * 4; // If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].

        if (RecastMesh.leftOn(verts, pin1, pi, pi1)) return RecastMesh.leftOn(verts, pi, pj, pin1) && RecastMesh.leftOn(verts, pj, pi, pi1); // Assume (i-1,i,i+1) not collinear.
        // else P[i] is reflex.

        return !(RecastMesh.leftOn(verts, pi, pj, pi1) && RecastMesh.leftOn(verts, pj, pi, pin1));
      }
    }, {
      key: "diagonalLoose",
      value: function diagonalLoose(i, j, n, verts, indices) {
        return RecastMesh.inConeLoose(i, j, n, verts, indices) && RecastMesh.diagonalieLoose(i, j, n, verts, indices);
      }
    }, {
      key: "triangulate",
      value: function triangulate(n, verts, indices, tris) {
        var ntris = 0; // The last bit of the index is used to indicate if the vertex can be removed.

        for (var i = 0; i < n; i++) {
          var i1 = RecastMesh.next(i, n);
          var i2 = RecastMesh.next(i1, n);

          if (RecastMesh.diagonal(i, i2, n, verts, indices)) {
            indices[i1] |= 0x80000000;
          }
        }

        while (n > 3) {
          var minLen = -1;
          var mini = -1;

          for (var _i7 = 0; _i7 < n; _i7++) {
            var _i8 = RecastMesh.next(_i7, n);

            if ((indices[_i8] & 0x80000000) != 0) {
              var p0 = (indices[_i7] & 0x0fffffff) * 4;
              var p2 = (indices[RecastMesh.next(_i8, n)] & 0x0fffffff) * 4;
              var dx = verts[p2 + 0] - verts[p0 + 0];
              var dy = verts[p2 + 2] - verts[p0 + 2];
              var len = dx * dx + dy * dy;

              if (minLen < 0 || len < minLen) {
                minLen = len;
                mini = _i7;
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

            for (var _i9 = 0; _i9 < n; _i9++) {
              var _i10 = RecastMesh.next(_i9, n);

              var _i11 = RecastMesh.next(_i10, n);

              if (RecastMesh.diagonalLoose(_i9, _i11, n, verts, indices)) {
                var _p = (indices[_i9] & 0x0fffffff) * 4;

                var _p2 = (indices[RecastMesh.next(_i11, n)] & 0x0fffffff) * 4;

                var _dx = verts[_p2 + 0] - verts[_p + 0];

                var _dy = verts[_p2 + 2] - verts[_p + 2];

                var _len = _dx * _dx + _dy * _dy;

                if (minLen < 0 || _len < minLen) {
                  minLen = _len;
                  mini = _i9;
                }
              }
            }

            if (mini == -1) {
              // The contour is messed up. This sometimes happens
              // if the contour simplification is too aggressive.
              return -ntris;
            }
          }

          var _i4 = mini;

          var _i5 = RecastMesh.next(_i4, n);

          var _i6 = RecastMesh.next(_i5, n);

          tris[ntris * 3] = indices[_i4] & 0x0ffffff;
          tris[ntris * 3 + 1] = indices[_i5] & 0x0ffffff;
          tris[ntris * 3 + 2] = indices[_i6] & 0x0ffffff;
          ntris++; // Removes P[i1] by copying P[i+1]...P[n-1] left one index.

          n--;

          for (var k = _i5; k < n; k++) {
            indices[k] = indices[k + 1];
          }

          if (_i5 >= n) _i5 = 0;
          _i4 = RecastMesh.prev(_i5, n); // Update diagonal flags.

          if (RecastMesh.diagonal(RecastMesh.prev(_i4, n), _i5, n, verts, indices)) indices[_i4] |= 0x80000000;else indices[_i4] &= 0x0ffffff;
          if (RecastMesh.diagonal(_i4, RecastMesh.next(_i5, n), n, verts, indices)) indices[_i5] |= 0x80000000;else indices[_i5] &= 0x0ffffff;
        } // Append the remaining triangle.


        tris[ntris * 3] = indices[0] & 0x0ffffff;
        tris[ntris * 3 + 1] = indices[1] & 0x0ffffff;
        tris[ntris * 3 + 2] = indices[2] & 0x0ffffff;
        ntris++;
        return ntris;
      }
    }, {
      key: "countPolyVerts",
      value: function countPolyVerts(p, j, nvp) {
        for (var i = 0; i < nvp; ++i) {
          if (p[i + j] == RecastConstants.RC_MESH_NULL_IDX) return i;
        }

        return nvp;
      }
    }, {
      key: "uleft",
      value: function uleft(verts, a, b, c) {
        return (verts[b + 0] - verts[a + 0]) * (verts[c + 2] - verts[a + 2]) - (verts[c + 0] - verts[a + 0]) * (verts[b + 2] - verts[a + 2]) < 0;
      }
    }, {
      key: "getPolyMergeValue",
      value: function getPolyMergeValue(polys, pa, pb, verts, nvp) {
        var ea = -1;
        var eb = -1;
        var na = RecastMesh.countPolyVerts(polys, pa, nvp);
        var nb = RecastMesh.countPolyVerts(polys, pb, nvp); // If the merged polygon would be too big, do not merge.

        if (na + nb - 2 > nvp) return [-1, ea, eb]; // Check if the polygons share an edge.

        for (var i = 0; i < na; ++i) {
          var va0 = polys[pa + i];
          var va1 = polys[pa + (i + 1) % na];

          if (va0 > va1) {
            var temp = va0;
            va0 = va1;
            va1 = temp;
          }

          for (var j = 0; j < nb; ++j) {
            var vb0 = polys[pb + j];
            var vb1 = polys[pb + (j + 1) % nb];

            if (vb0 > vb1) {
              var _temp = vb0;
              vb0 = vb1;
              vb1 = _temp;
            }

            if (va0 == vb0 && va1 == vb1) {
              ea = i;
              eb = j;
              break;
            }
          }
        } // No common edge, cannot merge.


        if (ea == -1 || eb == -1) return [-1, ea, eb]; // Check to see if the merged polygon would be convex.

        var va, vb, vc;
        va = polys[pa + (ea + na - 1) % na];
        vb = polys[pa + ea];
        vc = polys[pb + (eb + 2) % nb];
        if (!RecastMesh.uleft(verts, va * 3, vb * 3, vc * 3)) return [-1, ea, eb];
        va = polys[pb + (eb + nb - 1) % nb];
        vb = polys[pb + eb];
        vc = polys[pa + (ea + 2) % na];
        if (!RecastMesh.uleft(verts, va * 3, vb * 3, vc * 3)) return [-1, ea, eb];
        va = polys[pa + ea];
        vb = polys[pa + (ea + 1) % na];
        var dx = verts[va * 3 + 0] - verts[vb * 3 + 0];
        var dy = verts[va * 3 + 2] - verts[vb * 3 + 2];
        return [dx * dx + dy * dy, ea, eb];
      }
    }, {
      key: "mergePolyVerts",
      value: function mergePolyVerts(polys, pa, pb, ea, eb, tmp, nvp) {
        var na = RecastMesh.countPolyVerts(polys, pa, nvp);
        var nb = RecastMesh.countPolyVerts(polys, pb, nvp); // Merge polygons.

        polys.fill(RecastConstants.RC_MESH_NULL_IDX, tmp, tmp + nvp);
        var n = 0; // Add pa

        for (var i = 0; i < na - 1; ++i) {
          polys[tmp + n] = polys[pa + (ea + 1 + i) % na];
          n++;
        } // Add pb


        for (var _i12 = 0; _i12 < nb - 1; ++_i12) {
          polys[tmp + n] = polys[pb + (eb + 1 + _i12) % nb];
          n++;
        } //arraycopy(polys, tmp, polys, pa, nvp);


        arraycopy$1(polys, tmp, polys, pa, nvp); // for (let i = 0; i < nvp; i++) {
        // 	polys[pa + i] = polys[tmp + i];
        // }
      }
    }, {
      key: "pushFront",
      value: function pushFront(v, arr, an) {
        an++;

        for (var i = an - 1; i > 0; --i) {
          arr[i] = arr[i - 1];
        }

        arr[0] = v;
        return an;
      }
    }, {
      key: "pushBack",
      value: function pushBack(v, arr, an) {
        arr[an] = v;
        an++;
        return an;
      }
    }, {
      key: "canRemoveVertex",
      value: function canRemoveVertex(ctx, mesh, rem) {
        var nvp = mesh.nvp; // Count number of polygons to remove.

        var numTouchedVerts = 0;
        var numRemainingEdges = 0;

        for (var i = 0; i < mesh.npolys; ++i) {
          var p = i * nvp * 2;
          var nv = RecastMesh.countPolyVerts(mesh.polys, p, nvp);
          var numRemoved = 0;
          var numVerts = 0;

          for (var j = 0; j < nv; ++j) {
            if (mesh.polys[p + j] == rem) {
              numTouchedVerts++;
              numRemoved++;
            }

            numVerts++;
          }

          if (numRemoved != 0) {
            numRemainingEdges += numVerts - (numRemoved + 1);
          }
        } // There would be too few edges remaining to create a polygon.
        // This can happen for example when a tip of a triangle is marked
        // as dePolyion, but there are no other polys that share the vertex.
        // In this case, the vertex should not be removed.


        if (numRemainingEdges <= 2) return false; // Find edges which share the removed vertex.

        var maxEdges = numTouchedVerts * 2;
        var nedges = 0;
        var edges = new Array(maxEdges * 3);

        for (var _i13 = 0; _i13 < mesh.npolys; ++_i13) {
          var _p3 = _i13 * nvp * 2;

          var _nv = RecastMesh.countPolyVerts(mesh.polys, _p3, nvp); // Collect edges which touches the removed vertex.


          for (var _j2 = 0, k = _nv - 1; _j2 < _nv; k = _j2++) {
            if (mesh.polys[_p3 + _j2] == rem || mesh.polys[_p3 + k] == rem) {
              // Arrange edge so that a=rem.
              var a = mesh.polys[_p3 + _j2],
                  b = mesh.polys[_p3 + k];

              if (b == rem) {
                var temp = a;
                a = b;
                b = temp;
              } // Check if the edge exists


              var exists = false;

              for (var m = 0; m < nedges; ++m) {
                var e = m * 3;

                if (edges[e + 1] == b) {
                  // Exists, increment vertex share count.
                  edges[e + 2]++;
                  exists = true;
                }
              } // Add new edge.


              if (!exists) {
                var _e2 = nedges * 3;

                edges[_e2 + 0] = a;
                edges[_e2 + 1] = b;
                edges[_e2 + 2] = 1;
                nedges++;
              }
            }
          }
        } // There should be no more than 2 open edges.
        // This catches the case that two non-adjacent polygons
        // share the removed vertex. In that case, do not remove the vertex.


        var numOpenEdges = 0;

        for (var _i14 = 0; _i14 < nedges; ++_i14) {
          if (edges[_i14 * 3 + 2] < 2) numOpenEdges++;
        }

        if (numOpenEdges > 2) return false;
        return true;
      }
    }, {
      key: "removeVertex",
      value: function removeVertex(ctx, mesh, rem, maxTris) {
        var nvp = mesh.nvp; // Count number of polygons to remove.

        var numRemovedVerts = 0;

        for (var i = 0; i < mesh.npolys; ++i) {
          var p = i * nvp * 2;
          var nv = RecastMesh.countPolyVerts(mesh.polys, p, nvp);

          for (var j = 0; j < nv; ++j) {
            if (mesh.polys[p + j] == rem) numRemovedVerts++;
          }
        }

        var nedges = 0;
        var edges = new Array(numRemovedVerts * nvp * 4);
        var nhole = 0;
        var hole = new Array(numRemovedVerts * nvp);
        var nhreg = 0;
        var hreg = new Array(numRemovedVerts * nvp);
        var nharea = 0;
        var harea = new Array(numRemovedVerts * nvp);

        for (var _i15 = 0; _i15 < mesh.npolys; ++_i15) {
          var _p4 = _i15 * nvp * 2;

          var _nv2 = RecastMesh.countPolyVerts(mesh.polys, _p4, nvp);

          var hasRem = false;

          for (var _j3 = 0; _j3 < _nv2; ++_j3) {
            if (mesh.polys[_p4 + _j3] == rem) hasRem = true;
          }

          if (hasRem) {
            // Collect edges which does not touch the removed vertex.
            for (var _j4 = 0, k = _nv2 - 1; _j4 < _nv2; k = _j4++) {
              if (mesh.polys[_p4 + _j4] != rem && mesh.polys[_p4 + k] != rem) {
                var e = nedges * 4;
                edges[e + 0] = mesh.polys[_p4 + k];
                edges[e + 1] = mesh.polys[_p4 + _j4];
                edges[e + 2] = mesh.regs[_i15];
                edges[e + 3] = mesh.areas[_i15];
                nedges++;
              }
            } // Remove the polygon.


            var p2 = (mesh.npolys - 1) * nvp * 2;

            if (_p4 != p2) {
              arraycopy$1(mesh.polys, p2, mesh.polys, _p4, nvp);
            }

            mesh.polys.fill(RecastConstants.RC_MESH_NULL_IDX, _p4 + nvp, _p4 + nvp + nvp);
            mesh.regs[_i15] = mesh.regs[mesh.npolys - 1];
            mesh.areas[_i15] = mesh.areas[mesh.npolys - 1];
            mesh.npolys--;
            --_i15;
          }
        } // Remove vertex.


        for (var _i16 = rem; _i16 < mesh.nverts - 1; ++_i16) {
          mesh.verts[_i16 * 3 + 0] = mesh.verts[(_i16 + 1) * 3 + 0];
          mesh.verts[_i16 * 3 + 1] = mesh.verts[(_i16 + 1) * 3 + 1];
          mesh.verts[_i16 * 3 + 2] = mesh.verts[(_i16 + 1) * 3 + 2];
        }

        mesh.nverts--; // Adjust indices to match the removed vertex layout.

        for (var _i17 = 0; _i17 < mesh.npolys; ++_i17) {
          var _p5 = _i17 * nvp * 2;

          var _nv3 = RecastMesh.countPolyVerts(mesh.polys, _p5, nvp);

          for (var _j5 = 0; _j5 < _nv3; ++_j5) {
            if (mesh.polys[_p5 + _j5] > rem) mesh.polys[_p5 + _j5]--;
          }
        }

        for (var _i18 = 0; _i18 < nedges; ++_i18) {
          if (edges[_i18 * 4 + 0] > rem) edges[_i18 * 4 + 0]--;
          if (edges[_i18 * 4 + 1] > rem) edges[_i18 * 4 + 1]--;
        }

        if (nedges == 0) return; // Start with one vertex, keep appending connected
        // segments to the start and end of the hole.

        RecastMesh.pushBack(edges[0], hole, nhole);
        RecastMesh.pushBack(edges[2], hreg, nhreg);
        RecastMesh.pushBack(edges[3], harea, nharea);

        while (nedges != 0) {
          var match = false;

          for (var _i19 = 0; _i19 < nedges; ++_i19) {
            var ea = edges[_i19 * 4 + 0];
            var eb = edges[_i19 * 4 + 1];
            var r = edges[_i19 * 4 + 2];
            var a = edges[_i19 * 4 + 3];
            var add = false;

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
              edges[_i19 * 4 + 0] = edges[(nedges - 1) * 4 + 0];
              edges[_i19 * 4 + 1] = edges[(nedges - 1) * 4 + 1];
              edges[_i19 * 4 + 2] = edges[(nedges - 1) * 4 + 2];
              edges[_i19 * 4 + 3] = edges[(nedges - 1) * 4 + 3];
              --nedges;
              match = true;
              --_i19;
            }
          }

          if (!match) break;
        }

        var tris = new Array(nhole * 3);
        var tverts = new Array(nhole * 4);
        var thole = new Array(nhole); // Generate temp vertex array for triangulation.

        for (var _i20 = 0; _i20 < nhole; ++_i20) {
          var pi = hole[_i20];
          tverts[_i20 * 4 + 0] = mesh.verts[pi * 3 + 0];
          tverts[_i20 * 4 + 1] = mesh.verts[pi * 3 + 1];
          tverts[_i20 * 4 + 2] = mesh.verts[pi * 3 + 2];
          tverts[_i20 * 4 + 3] = 0;
          thole[_i20] = _i20;
        } // Triangulate the hole.


        var ntris = RecastMesh.triangulate(nhole, tverts, thole, tris);

        if (ntris < 0) {
          ntris = -ntris;
          ctx.warn("removeVertex: triangulate() returned bad results.");
        } // Merge the hole triangles back to polygons.


        var polys = new Array((ntris + 1) * nvp);
        var pregs = new Array(ntris);
        var pareas = new Array(ntris);

        var npolys = 0;
        polys.fill(RecastConstants.RC_MESH_NULL_IDX, 0, ntris * nvp);

        for (var _j6 = 0; _j6 < ntris; ++_j6) {
          var t = _j6 * 3;

          if (tris[t + 0] != tris[t + 1] && tris[t + 0] != tris[t + 2] && tris[t + 1] != tris[t + 2]) {
            polys[npolys * nvp + 0] = hole[tris[t + 0]];
            polys[npolys * nvp + 1] = hole[tris[t + 1]];
            polys[npolys * nvp + 2] = hole[tris[t + 2]]; // If this polygon covers multiple region types then
            // mark it as such

            if (hreg[tris[t + 0]] != hreg[tris[t + 1]] || hreg[tris[t + 1]] != hreg[tris[t + 2]]) pregs[npolys] = RecastConstants.RC_MULTIPLE_REGS;else pregs[npolys] = hreg[tris[t + 0]];
            pareas[npolys] = harea[tris[t + 0]];
            npolys++;
          }
        }

        if (npolys == 0) return; // Merge polygons.

        if (nvp > 3) {
          for (;;) {
            // Find best polygons to merge.
            var bestMergeVal = 0;
            var bestPa = 0,
                bestPb = 0,
                bestEa = 0,
                bestEb = 0;

            for (var _j7 = 0; _j7 < npolys - 1; ++_j7) {
              var pj = _j7 * nvp;

              for (var _k = _j7 + 1; _k < npolys; ++_k) {
                var pk = _k * nvp;
                var veaeb = RecastMesh.getPolyMergeValue(polys, pj, pk, mesh.verts, nvp);
                var v = veaeb[0];
                var _ea = veaeb[1];
                var _eb = veaeb[2];

                if (v > bestMergeVal) {
                  bestMergeVal = v;
                  bestPa = _j7;
                  bestPb = _k;
                  bestEa = _ea;
                  bestEb = _eb;
                }
              }
            }

            if (bestMergeVal > 0) {
              // Found best, merge.
              var pa = bestPa * nvp;
              var pb = bestPb * nvp;
              RecastMesh.mergePolyVerts(polys, pa, pb, bestEa, bestEb, tmpPoly, nvp);
              if (pregs[bestPa] != pregs[bestPb]) pregs[bestPa] = RecastConstants.RC_MULTIPLE_REGS;
              var last = (npolys - 1) * nvp;

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
        } // Store polygons.


        for (var _i21 = 0; _i21 < npolys; ++_i21) {
          if (mesh.npolys >= maxTris) break;

          var _p6 = mesh.npolys * nvp * 2;

          mesh.polys.fill(RecastConstants.RC_MESH_NULL_IDX, _p6, _p6 + nvp * 2);

          for (var _j8 = 0; _j8 < nvp; ++_j8) {
            mesh.polys[_p6 + _j8] = polys[_i21 * nvp + _j8];
          }

          mesh.regs[mesh.npolys] = pregs[_i21];
          mesh.areas[mesh.npolys] = pareas[_i21];
          mesh.npolys++;

          if (mesh.npolys > maxTris) {
            throw new RuntimeException("removeVertex: Too many polygons " + mesh.npolys + " (max:" + maxTris + ".");
          }
        }
      } /// @par
      ///
      /// @note If the mesh data is to be used to construct a Detour navigation mesh, then the upper
      /// limit must be retricted to <= #DT_VERTS_PER_POLYGON.
      ///
      /// @see rcAllocPolyMesh, rcContourSet, rcPolyMesh, rcConfig

    }, {
      key: "buildPolyMesh",
      value: function buildPolyMesh(ctx, cset, nvp) {
        ctx.startTimer("BUILD_POLYMESH");
        var mesh = new PolyMesh();
        RecastVectors$1.copy3(mesh.bmin, cset.bmin, 0);
        RecastVectors$1.copy3(mesh.bmax, cset.bmax, 0);
        mesh.cs = cset.cs;
        mesh.ch = cset.ch;
        mesh.borderSize = cset.borderSize;
        mesh.maxEdgeError = cset.maxError;
        var maxVertices = 0;
        var maxTris = 0;
        var maxVertsPerCont = 0;

        for (var i = 0; i < cset.conts.length; ++i) {
          // Skip null contours.
          if (cset.conts[i].nverts < 3) continue;
          maxVertices += cset.conts[i].nverts;
          maxTris += cset.conts[i].nverts - 2;
          maxVertsPerCont = Math.max(maxVertsPerCont, cset.conts[i].nverts);
        }

        if (maxVertices >= 0xfffe) {
          throw new RuntimeException("rcBuildPolyMesh: Too many vertices " + maxVertices);
        }

        var vflags = new Array(maxVertices).fill(0);
        vflags.fill(0);
        mesh.verts = new Array(maxVertices * 3).fill(0);
        mesh.polys = new Array(maxTris * nvp * 2).fill(RecastConstants.RC_MESH_NULL_IDX); // Arrays.fill(mesh.polys, RecastConstants.RC_MESH_NULL_IDX);

        mesh.regs = new Array(maxTris).fill(0);
        mesh.areas = new Array(maxTris).fill(0);
        mesh.nverts = 0;
        mesh.npolys = 0;
        mesh.nvp = nvp;
        mesh.maxpolys = maxTris;
        var nextVert = new Array(maxVertices);
        var firstVert = new Array(RecastMesh.VERTEX_BUCKET_COUNT);

        for (var _i22 = 0; _i22 < RecastMesh.VERTEX_BUCKET_COUNT; ++_i22) {
          firstVert[_i22] = -1;
        }

        var indices = new Array(maxVertsPerCont);
        var tris = new Array(maxVertsPerCont * 3);
        var polys = new Array((maxVertsPerCont + 1) * nvp);
        var tmpPoly = maxVertsPerCont * nvp;

        for (var _i23 = 0; _i23 < cset.conts.length; ++_i23) {
          var cont = cset.conts[_i23]; // Skip null contours.

          if (cont.nverts < 3) continue; // Triangulate contour

          for (var j = 0; j < cont.nverts; ++j) {
            indices[j] = j;
          }

          var ntris = RecastMesh.triangulate(cont.nverts, cont.verts, indices, tris);

          if (ntris <= 0) {
            // Bad triangulation, should not happen.
            ctx.warn("buildPolyMesh: Bad triangulation Contour " + _i23 + ".");
            ntris = -ntris;
          } // Add and merge vertices.


          for (var _j9 = 0; _j9 < cont.nverts; ++_j9) {
            var v = _j9 * 4;
            var inv = RecastMesh.addVertex(cont.verts[v + 0], cont.verts[v + 1], cont.verts[v + 2], mesh.verts, firstVert, nextVert, mesh.nverts);
            indices[_j9] = inv[0];
            mesh.nverts = inv[1];

            if ((cont.verts[v + 3] & RecastConstants.RC_BORDER_VERTEX) != 0) {
              // This vertex should be removed.
              vflags[indices[_j9]] = 1;
            }
          } // Build initial polygons.


          var npolys = 0; // Arrays.fill(polys, RecastConstants.RC_MESH_NULL_IDX);

          polys.fill(RecastConstants.RC_MESH_NULL_IDX);

          for (var _j10 = 0; _j10 < ntris; ++_j10) {
            var t = _j10 * 3;

            if (tris[t + 0] != tris[t + 1] && tris[t + 0] != tris[t + 2] && tris[t + 1] != tris[t + 2]) {
              polys[npolys * nvp + 0] = indices[tris[t + 0]];
              polys[npolys * nvp + 1] = indices[tris[t + 1]];
              polys[npolys * nvp + 2] = indices[tris[t + 2]];
              npolys++;
            }
          }

          if (npolys == 0) continue; // Merge polygons.

          if (nvp > 3) {
            for (;;) {
              // Find best polygons to merge.
              var bestMergeVal = 0;
              var bestPa = 0,
                  bestPb = 0,
                  bestEa = 0,
                  bestEb = 0;

              for (var _j11 = 0; _j11 < npolys - 1; ++_j11) {
                var pj = _j11 * nvp;

                for (var k = _j11 + 1; k < npolys; ++k) {
                  var pk = k * nvp;
                  var veaeb = RecastMesh.getPolyMergeValue(polys, pj, pk, mesh.verts, nvp);
                  var _v4 = veaeb[0];
                  var ea = veaeb[1];
                  var eb = veaeb[2];

                  if (_v4 > bestMergeVal) {
                    bestMergeVal = _v4;
                    bestPa = _j11;
                    bestPb = k;
                    bestEa = ea;
                    bestEb = eb;
                  }
                }
              }

              if (bestMergeVal > 0) {
                // Found best, merge.
                var pa = bestPa * nvp;
                var pb = bestPb * nvp;
                RecastMesh.mergePolyVerts(polys, pa, pb, bestEa, bestEb, tmpPoly, nvp);
                var lastPoly = (npolys - 1) * nvp;

                if (pb != lastPoly) {
                  // arraycopy(polys, lastPoly, polys, pb, nvp);
                  for (var _i24 = 0; _i24 < nvp; _i24++) {
                    polys[pb + _i24] = polys[lastPoly + _i24];
                    if (polys[28] == 80) console.log("break");
                  }
                }

                npolys--;
              } else {
                // Could not merge any polygons, stop.
                break;
              }
            }
          } //Here they are the same
          // Store polygons.


          for (var _j12 = 0; _j12 < npolys; ++_j12) {
            var p = mesh.npolys * nvp * 2;
            var q = _j12 * nvp;

            for (var _k2 = 0; _k2 < nvp; ++_k2) {
              mesh.polys[p + _k2] = polys[q + _k2];
            }

            mesh.regs[mesh.npolys] = cont.reg;
            mesh.areas[mesh.npolys] = cont.area;
            mesh.npolys++;

            if (mesh.npolys > maxTris) {
              throw new RuntimeException("rcBuildPolyMesh: Too many polygons " + mesh.npolys + " (max:" + maxTris + ").");
            }
          }
        } // fs.writeFileSync("./pmeshjs.txt", JSON.stringify(mesh.polys))
        //Now we are the same here
        // Remove edge vertices.


        for (var _i25 = 0; _i25 < mesh.nverts; ++_i25) {
          if (vflags[_i25] != 0) {
            if (!RecastMesh.canRemoveVertex(ctx, mesh, _i25)) continue;
            RecastMesh.removeVertex(ctx, mesh, _i25, maxTris); // Remove vertex
            // Note: mesh.nverts is already decremented inside removeVertex()!
            // Fixup vertex flags

            for (var _j13 = _i25; _j13 < mesh.nverts; ++_j13) {
              vflags[_j13] = vflags[_j13 + 1];
            }

            --_i25;
          }
        } // Calculate adjacency.


        this.buildMeshAdjacency(mesh.polys, mesh.npolys, mesh.nverts, nvp); // Find portal edges

        if (mesh.borderSize > 0) {
          var w = cset.width;
          var h = cset.height;

          for (var _i26 = 0; _i26 < mesh.npolys; ++_i26) {
            var _p7 = _i26 * 2 * nvp;

            for (var _j14 = 0; _j14 < nvp; ++_j14) {
              if (mesh.polys[_p7 + _j14] == RecastConstants.RC_MESH_NULL_IDX) break; // Skip connected edges.

              if (mesh.polys[_p7 + nvp + _j14] != RecastConstants.RC_MESH_NULL_IDX) continue;
              var nj = _j14 + 1;
              if (nj >= nvp || mesh.polys[_p7 + nj] == RecastConstants.RC_MESH_NULL_IDX) nj = 0;
              var va = mesh.polys[_p7 + _j14] * 3;
              var vb = mesh.polys[_p7 + nj] * 3;
              if (mesh.verts[va + 0] == 0 && mesh.verts[vb + 0] == 0) mesh.polys[_p7 + nvp + _j14] = 0x8000 | 0;else if (mesh.verts[va + 2] == h && mesh.verts[vb + 2] == h) mesh.polys[_p7 + nvp + _j14] = 0x8000 | 1;else if (mesh.verts[va + 0] == w && mesh.verts[vb + 0] == w) mesh.polys[_p7 + nvp + _j14] = 0x8000 | 2;else if (mesh.verts[va + 2] == 0 && mesh.verts[vb + 2] == 0) mesh.polys[_p7 + nvp + _j14] = 0x8000 | 3;
            }
          }
        } // Just allocate the mesh flags array. The user is resposible to fill it.


        mesh.flags = new Array(mesh.npolys);
        mesh.flags.fill(0);

        if (mesh.nverts > 0xffff) {
          throw new RuntimeException("rcBuildPolyMesh: The resulting mesh has too many vertices " + mesh.nverts + " (max " + 0xffff + "). Data can be corrupted.");
        }

        if (mesh.npolys > 0xffff) {
          throw new RuntimeException("rcBuildPolyMesh: The resulting mesh has too many polygons " + mesh.npolys + " (max " + 0xffff + "). Data can be corrupted.");
        }

        ctx.stopTimer("BUILD_POLYMESH");
        return mesh;
      } /// @see rcAllocPolyMesh, rcPolyMesh

    }, {
      key: "mergePolyMeshes",
      value: function mergePolyMeshes(ctx, meshes, nmeshes) {
        if (nmeshes == 0 || meshes == null) return null;
        ctx.startTimer("MERGE_POLYMESH");
        var mesh = new PolyMesh();
        mesh.nvp = meshes[0].nvp;
        mesh.cs = meshes[0].cs;
        mesh.ch = meshes[0].ch;
        RecastVectors$1.copy(mesh.bmin, meshes[0].bmin, 0);
        RecastVectors$1.copy(mesh.bmax, meshes[0].bmax, 0);
        var maxVerts = 0;
        var maxPolys = 0;
        var maxVertsPerMesh = 0;

        for (var i = 0; i < nmeshes; ++i) {
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
        var nextVert = new Array(maxVerts);
        var firstVert = new Array(RecastMesh.VERTEX_BUCKET_COUNT);

        for (var _i27 = 0; _i27 < RecastMesh.VERTEX_BUCKET_COUNT; ++_i27) {
          firstVert[_i27] = -1;
        }

        var vremap = new Array(maxVertsPerMesh);

        for (var _i28 = 0; _i28 < nmeshes; ++_i28) {
          var pmesh = meshes[_i28];
          var ox = Math.floor((pmesh.bmin[0] - mesh.bmin[0]) / mesh.cs + 0.5);
          var oz = Math.floor((pmesh.bmin[2] - mesh.bmin[2]) / mesh.cs + 0.5);
          var isMinX = ox == 0;
          var isMinZ = oz == 0;
          var isMaxX = Math.floor((mesh.bmax[0] - pmesh.bmax[0]) / mesh.cs + 0.5) == 0;
          var isMaxZ = Math.floor((mesh.bmax[2] - pmesh.bmax[2]) / mesh.cs + 0.5) == 0;
          var isOnBorder = isMinX || isMinZ || isMaxX || isMaxZ;

          for (var j = 0; j < pmesh.nverts; ++j) {
            var v = j * 3;
            var inv = addVertex(pmesh.verts[v + 0] + ox, pmesh.verts[v + 1], pmesh.verts[v + 2] + oz, mesh.verts, firstVert, nextVert, mesh.nverts);
            vremap[j] = inv[0];
            mesh.nverts = inv[1];
          }

          for (var _j15 = 0; _j15 < pmesh.npolys; ++_j15) {
            var tgt = mesh.npolys * 2 * mesh.nvp;
            var src = _j15 * 2 * mesh.nvp;
            mesh.regs[mesh.npolys] = pmesh.regs[_j15];
            mesh.areas[mesh.npolys] = pmesh.areas[_j15];
            mesh.flags[mesh.npolys] = pmesh.flags[_j15];
            mesh.npolys++;

            for (var k = 0; k < mesh.nvp; ++k) {
              if (pmesh.polys[src + k] == RecastConstants.RC_MESH_NULL_IDX) break;
              mesh.polys[tgt + k] = vremap[pmesh.polys[src + k]];
            }

            if (isOnBorder) {
              for (var _k3 = mesh.nvp; _k3 < mesh.nvp * 2; ++_k3) {
                if ((pmesh.polys[src + _k3] & 0x8000) != 0 && pmesh.polys[src + _k3] != 0xffff) {
                  var dir = pmesh.polys[src + _k3] & 0xf;

                  switch (dir) {
                    case 0:
                      // Portal x-
                      if (isMinX) mesh.polys[tgt + _k3] = pmesh.polys[src + _k3];
                      break;

                    case 1:
                      // Portal z+
                      if (isMaxZ) mesh.polys[tgt + _k3] = pmesh.polys[src + _k3];
                      break;

                    case 2:
                      // Portal x+
                      if (isMaxX) mesh.polys[tgt + _k3] = pmesh.polys[src + _k3];
                      break;

                    case 3:
                      // Portal z-
                      if (isMinZ) mesh.polys[tgt + _k3] = pmesh.polys[src + _k3];
                      break;
                  }
                }
              }
            }
          }
        } // Calculate adjacency.


        buildMeshAdjacency(mesh.polys, mesh.npolys, mesh.nverts, mesh.nvp);

        if (mesh.nverts > 0xffff) {
          throw new RuntimeException("rcBuildPolyMesh: The resulting mesh has too many vertices " + mesh.nverts + " (max " + 0xffff + "). Data can be corrupted.");
        }

        if (mesh.npolys > 0xffff) {
          throw new RuntimeException("rcBuildPolyMesh: The resulting mesh has too many polygons " + mesh.npolys + " (max " + 0xffff + "). Data can be corrupted.");
        }

        ctx.stopTimer("MERGE_POLYMESH");
        return mesh;
      }
    }, {
      key: "copyPolyMesh",
      value: function copyPolyMesh(ctx, src) {
        var dst = new PolyMesh();
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
    }]);

    return RecastMesh;
  }();

  _defineProperty(RecastMesh, "VERTEX_BUCKET_COUNT", 1 << 12);

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
  var Contour = function Contour() {
    _classCallCheck(this, Contour);

    _defineProperty(this, "verts", []);

    _defineProperty(this, "nverts", 0);

    _defineProperty(this, "rverts", []);

    _defineProperty(this, "nrverts", 0);

    _defineProperty(this, "area", 0);

    _defineProperty(this, "reg", 0);
  };

  var _temp$2, _temp2$2, _temp3;

  var RecastContour = /*#__PURE__*/function () {
    function RecastContour() {
      _classCallCheck(this, RecastContour);

      _defineProperty(this, "CompareHoles", function (a, b) {
        if (a.minx == b.minx) {
          if (a.minz < b.minz) return -1;
          if (a.minz > b.minz) return 1;
        } else {
          if (a.minx < b.minx) return -1;
          if (a.minx > b.minx) return 1;
        }

        return 0;
      });

      _defineProperty(this, "CompareDiagDist", function (va, vb) {
        var a = va;
        var b = vb;
        if (a.dist < b.dist) return -1;
        if (a.dist > b.dist) return 1;
        return 0;
      });
    }

    _createClass(RecastContour, [{
      key: "mergeRegionHoles",
      value: function mergeRegionHoles(ctx, region) {
        // Sort holes from left to right.
        for (var i = 0; i < region.nholes; i++) {
          var minleft = findLeftMostVertex(region.holes[i].contour);
          region.holes[i].minx = minleft[0];
          region.holes[i].minz = minleft[1];
          region.holes[i].leftmost = minleft[2];
        }

        Arrays.sort(region.holes, new CompareHoles());
        var maxVerts = region.outline.nverts;

        for (var _i = 0; _i < region.nholes; _i++) {
          maxVerts += region.holes[_i].contour.nverts;
        }

        diags = new Array(maxVerts);

        for (var pd = 0; pd < maxVerts; pd++) {
          diags[pd] = new PotentialDiagonal();
        }

        var outline = region.outline; // Merge holes into the outline one by one.

        for (var _i2 = 0; _i2 < region.nholes; _i2++) {
          var hole = region.holes[_i2].contour;
          var index = -1;
          var bestVertex = region.holes[_i2].leftmost;

          for (var iter = 0; iter < hole.nverts; iter++) {
            // Find potential diagonals.
            // The 'best' vertex must be in the cone described by 3 cosequtive vertices of the outline.
            // ..o j-1
            //   |
            //   |   * best
            //   |
            // j o-----o j+1
            //         :
            var ndiags = 0;
            var corner = bestVertex * 4;

            for (var j = 0; j < outline.nverts; j++) {
              if (inCone(j, outline.nverts, outline.verts, corner, hole.verts)) {
                var dx = outline.verts[j * 4 + 0] - hole.verts[corner + 0];
                var dz = outline.verts[j * 4 + 2] - hole.verts[corner + 2];
                diags[ndiags].vert = j;
                diags[ndiags].dist = dx * dx + dz * dz;
                ndiags++;
              }
            } // Sort potential diagonals by distance, we want to make the connection as short as possible.


            Arrays.sort(diags, 0, ndiags, new CompareDiagDist()); // Find a diagonal that is not intersecting the outline not the remaining holes.

            index = -1;

            for (var _j = 0; _j < ndiags; _j++) {
              var pt = diags[_j].vert * 4;
              var intersect = intersectSegCountour(pt, corner, diags[_i2].vert, outline.nverts, outline.verts, outline.verts, hole.verts);

              for (var k = _i2; k < region.nholes && !intersect; k++) {
                intersect |= intersectSegCountour(pt, corner, -1, region.holes[k].contour.nverts, region.holes[k].contour.verts, outline.verts, hole.verts);
              }

              if (!intersect) {
                index = diags[_j].vert;
                break;
              }
            } // If found non-intersecting diagonal, stop looking.


            if (index != -1) break; // All the potential diagonals for the current vertex were intersecting, try next vertex.

            bestVertex = (bestVertex + 1) % hole.nverts;
          }

          if (index == -1) {
            ctx.warn("mergeHoles: Failed to find merge points for");
            continue;
          }

          mergeContours(region.outline, hole, index, bestVertex);
        }
      } /// @par
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

    }], [{
      key: "getCornerHeight",
      value: function getCornerHeight(x, y, i, dir, chf, isBorderVertex) {
        var s = chf.spans[i];
        var ch = s.y;
        var dirp = dir + 1 & 0x3;
        var regs = [0, 0, 0, 0]; // Combine region and area codes in order to prevent
        // border vertices which are in between two areas to be removed.

        regs[0] = chf.spans[i].reg | chf.areas[i] << 16;

        if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
          var ax = x + RecastCommon.GetDirOffsetX(dir);
          var ay = y + RecastCommon.GetDirOffsetY(dir);
          var ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(s, dir);
          var as = chf.spans[ai];
          ch = Math.max(ch, as.y);
          regs[1] = chf.spans[ai].reg | chf.areas[ai] << 16;

          if (RecastCommon.GetCon(as, dirp) != RecastConstants.RC_NOT_CONNECTED) {
            var ax2 = ax + RecastCommon.GetDirOffsetX(dirp);
            var ay2 = ay + RecastCommon.GetDirOffsetY(dirp);
            var ai2 = chf.cells[ax2 + ay2 * chf.width].index + RecastCommon.GetCon(as, dirp);
            var as2 = chf.spans[ai2];
            ch = Math.max(ch, as2.y);
            regs[2] = chf.spans[ai2].reg | chf.areas[ai2] << 16;
          }
        }

        if (RecastCommon.GetCon(s, dirp) != RecastConstants.RC_NOT_CONNECTED) {
          var _ax = x + RecastCommon.GetDirOffsetX(dirp);

          var _ay = y + RecastCommon.GetDirOffsetY(dirp);

          var _ai = chf.cells[_ax + _ay * chf.width].index + RecastCommon.GetCon(s, dirp);

          var _as = chf.spans[_ai];
          ch = Math.max(ch, _as.y);
          regs[3] = chf.spans[_ai].reg | chf.areas[_ai] << 16;

          if (RecastCommon.GetCon(_as, dir) != RecastConstants.RC_NOT_CONNECTED) {
            var _ax2 = _ax + RecastCommon.GetDirOffsetX(dir);

            var _ay2 = _ay + RecastCommon.GetDirOffsetY(dir);

            var _ai2 = chf.cells[_ax2 + _ay2 * chf.width].index + RecastCommon.GetCon(_as, dir);

            var _as2 = chf.spans[_ai2];
            ch = Math.max(ch, _as2.y);
            regs[2] = chf.spans[_ai2].reg | chf.areas[_ai2] << 16;
          }
        } // Check if the vertex is special edge vertex, these vertices will be removed later.

        return ch;
      }
    }, {
      key: "walkContour",
      value: function walkContour(x, y, i, chf, flags, points) {
        // Choose the first non-connected edge
        var dir = 0;

        while ((flags[i] & 1 << dir) == 0) {
          dir++;
        }

        var startDir = dir;
        var starti = i;
        var area = chf.areas[i];
        var iter = 0;

        while (++iter < 40000) {
          if ((flags[i] & 1 << dir) != 0) {
            // Choose the edge corner
            var isBorderVertex = false;
            var isAreaBorder = false;
            var px = x;
            var py = RecastContour.getCornerHeight(x, y, i, dir, chf, isBorderVertex);
            var pz = y;

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

            var r = 0;
            var s = chf.spans[i];

            if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
              var ax = x + RecastCommon.GetDirOffsetX(dir);
              var ay = y + RecastCommon.GetDirOffsetY(dir);
              var ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(s, dir);
              r = chf.spans[ai].reg;
              if (area != chf.areas[ai]) isAreaBorder = true;
            }

            if (isBorderVertex) r |= RecastConstants.RC_BORDER_VERTEX;
            if (isAreaBorder) r |= RecastConstants.RC_AREA_BORDER;
            points.push(px);
            points.push(py);
            points.push(pz);
            points.push(r);
            flags[i] &= ~(1 << dir); // Remove visited edges

            dir = dir + 1 & 0x3; // Rotate CW
          } else {
            var ni = -1;
            var nx = x + RecastCommon.GetDirOffsetX(dir);
            var ny = y + RecastCommon.GetDirOffsetY(dir);
            var _s = chf.spans[i];

            if (RecastCommon.GetCon(_s, dir) != RecastConstants.RC_NOT_CONNECTED) {
              var nc = chf.cells[nx + ny * chf.width];
              ni = nc.index + RecastCommon.GetCon(_s, dir);
            }

            if (ni == -1) {
              // Should not happen.
              return;
            }

            x = nx;
            y = ny;
            i = ni;
            dir = dir + 3 & 0x3; // Rotate CCW
          }

          if (starti == i && startDir == dir) {
            break;
          }
        }
      }
    }, {
      key: "distancePtSeg",
      value: function distancePtSeg(x, z, px, pz, qx, qz) {
        var pqx = qx - px;
        var pqz = qz - pz;
        var dx = x - px;
        var dz = z - pz;
        var d = pqx * pqx + pqz * pqz;
        var t = pqx * dx + pqz * dz;
        if (d > 0) t /= d;
        if (t < 0) t = 0;else if (t > 1) t = 1;
        dx = px + t * pqx - x;
        dz = pz + t * pqz - z;
        return dx * dx + dz * dz;
      }
    }, {
      key: "simplifyContour",
      value: function simplifyContour(points, simplified, maxError, maxEdgeLen, buildFlags) {
        // Add initial points.
        var hasConnections = false;

        for (var i = 0; i < points.length; i += 4) {
          if ((points[i + 3] & RecastConstants.RC_CONTOUR_REG_MASK) != 0) {
            hasConnections = true;
            break;
          }
        }

        if (hasConnections) {
          // The contour has some portals to other regions.
          // Add a new poPoly to every location where the region changes.
          for (var _i3 = 0, ni = points.length / 4; _i3 < ni; ++_i3) {
            var ii = (_i3 + 1) % ni;
            var differentRegs = (points[_i3 * 4 + 3] & RecastConstants.RC_CONTOUR_REG_MASK) != (points[ii * 4 + 3] & RecastConstants.RC_CONTOUR_REG_MASK);
            var areaBorders = (points[_i3 * 4 + 3] & RecastConstants.RC_AREA_BORDER) != (points[ii * 4 + 3] & RecastConstants.RC_AREA_BORDER);

            if (differentRegs || areaBorders) {
              simplified.push(points[_i3 * 4 + 0]);
              simplified.push(points[_i3 * 4 + 1]);
              simplified.push(points[_i3 * 4 + 2]);
              simplified.push(_i3);
            }
          }
        }

        if (simplified.length == 0) {
          // If there is no connections at all,
          // create some initial points for the simplification process.
          // Find lower-left and upper-right vertices of the contour.
          var llx = points[0];
          var lly = points[1];
          var llz = points[2];
          var lli = 0;
          var urx = points[0];
          var ury = points[1];
          var urz = points[2];
          var uri = 0;

          for (var _i4 = 0; _i4 < points.length; _i4 += 4) {
            var x = points[_i4 + 0];
            var y = points[_i4 + 1];
            var z = points[_i4 + 2];

            if (x < llx || x == llx && z < llz) {
              llx = x;
              lly = y;
              llz = z;
              lli = _i4 / 4;
            }

            if (x > urx || x == urx && z > urz) {
              urx = x;
              ury = y;
              urz = z;
              uri = _i4 / 4;
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
        } // Add points until all raw points are within
        // error tolerance to the simplified shape.


        var pn = points.length / 4;

        for (var _i5 = 0; _i5 < simplified.length / 4;) {
          // console.log(simplified.length)
          var _ii = (_i5 + 1) % (simplified.length / 4);

          var ax = simplified[_i5 * 4 + 0];
          var az = simplified[_i5 * 4 + 2];
          var ai = simplified[_i5 * 4 + 3];
          var bx = simplified[_ii * 4 + 0];
          var bz = simplified[_ii * 4 + 2];
          var bi = simplified[_ii * 4 + 3]; // Find maximum deviation from the segment.

          var maxd = 0;
          var maxi = -1;
          var ci = void 0,
              cinc = void 0,
              endi = void 0; // Traverse the segment in lexilogical order so that the
          // max deviation is calculated similarly when traversing
          // opposite segments.

          if (bx > ax || bx == ax && bz > az) {
            cinc = 1;
            ci = (ai + cinc) % pn;
            endi = bi;
          } else {
            cinc = pn - 1;
            ci = (bi + cinc) % pn;
            endi = ai;
            var temp = ax;
            ax = bx;
            bx = temp;
            temp = az;
            az = bz;
            bz = temp;
          } // Tessellate only outer edges or edges between areas.
          // console.log("start")


          if ((points[ci * 4 + 3] & RecastConstants.RC_CONTOUR_REG_MASK) == 0 || (points[ci * 4 + 3] & RecastConstants.RC_AREA_BORDER) != 0) {
            while (ci != endi) {
              // if(Math.random() < .01) console.log(`${ci} ${endi}`);
              var d = RecastContour.distancePtSeg(points[ci * 4 + 0], points[ci * 4 + 2], ax, az, bx, bz);

              if (d > maxd) {
                maxd = d;
                maxi = ci;
              }

              ci = (ci + cinc) % pn;
            }
          } // console.log("stop")
          // If the max deviation is larger than accepted error,
          // add new point, else continue to next segment.


          if (maxi != -1 && maxd > maxError * maxError) {
            // Add the point.
            simplified.splice((_i5 + 1) * 4 + 0, 0, points[maxi * 4 + 0]);
            simplified.splice((_i5 + 1) * 4 + 1, 0, points[maxi * 4 + 1]);
            simplified.splice((_i5 + 1) * 4 + 2, 0, points[maxi * 4 + 2]);
            simplified.splice((_i5 + 1) * 4 + 3, 0, maxi);
          } else {
            ++_i5;
          }
        } // Split too let edges.


        if (maxEdgeLen > 0 && (buildFlags & (RecastConstants.RC_CONTOUR_TESS_WALL_EDGES | RecastConstants.RC_CONTOUR_TESS_AREA_EDGES)) != 0) {
          for (var _i6 = 0; _i6 < simplified.length / 4;) {
            var _ii2 = (_i6 + 1) % (simplified.length / 4);

            var _ax3 = simplified[_i6 * 4 + 0];
            var _az = simplified[_i6 * 4 + 2];
            var _ai3 = simplified[_i6 * 4 + 3];
            var _bx = simplified[_ii2 * 4 + 0];
            var _bz = simplified[_ii2 * 4 + 2];
            var _bi = simplified[_ii2 * 4 + 3]; // Find maximum deviation from the segment.

            var _maxi = -1;

            var _ci = (_ai3 + 1) % pn; // Tessellate only outer edges or edges between areas.


            var tess = false; // Wall edges.

            if ((buildFlags & RecastConstants.RC_CONTOUR_TESS_WALL_EDGES) != 0 && (points[_ci * 4 + 3] & RecastConstants.RC_CONTOUR_REG_MASK) == 0) tess = true; // Edges between areas.

            if ((buildFlags & RecastConstants.RC_CONTOUR_TESS_AREA_EDGES) != 0 && (points[_ci * 4 + 3] & RecastConstants.RC_AREA_BORDER) != 0) tess = true;

            if (tess) {
              var dx = _bx - _ax3;
              var dz = _bz - _az;

              if (dx * dx + dz * dz > maxEdgeLen * maxEdgeLen) {
                // Round based on the segments in lexilogical order so that the
                // max tesselation is consistent regardles in which direction
                // segments are traversed.
                var n = _bi < _ai3 ? _bi + pn - _ai3 : _bi - _ai3;

                if (n > 1) {
                  if (_bx > _ax3 || _bx == _ax3 && _bz > _az) _maxi = Math.floor((_ai3 + n / 2) % pn);else _maxi = Math.floor((_ai3 + (n + 1) / 2) % pn);
                }
              }
            } // If the max deviation is larger than accepted error,
            // add new point, else continue to next segment.


            if (_maxi != -1) {
              // Add the point.
              simplified.splice((_i6 + 1) * 4 + 0, 0, points[_maxi * 4 + 0]);
              simplified.splice((_i6 + 1) * 4 + 1, 0, points[_maxi * 4 + 1]);
              simplified.splice((_i6 + 1) * 4 + 2, 0, points[_maxi * 4 + 2]);
              simplified.splice((_i6 + 1) * 4 + 3, 0, _maxi);
            } else {
              ++_i6;
            }
          }
        }

        for (var _i7 = 0; _i7 < simplified.length / 4; ++_i7) {
          // The edge vertex flag is take from the current raw point,
          // and the neighbour region is take from the next raw point.
          var _ai4 = (simplified[_i7 * 4 + 3] + 1) % pn;

          var _bi2 = simplified[_i7 * 4 + 3];
          simplified[_i7 * 4 + 3] = points[_ai4 * 4 + 3] & (RecastConstants.RC_CONTOUR_REG_MASK | RecastConstants.RC_AREA_BORDER) | points[_bi2 * 4 + 3] & RecastConstants.RC_BORDER_VERTEX;
        }
      }
    }, {
      key: "calcAreaOfPolygon2D",
      value: function calcAreaOfPolygon2D(verts, nverts) {
        var area = 0;

        for (var i = 0, j = nverts - 1; i < nverts; j = i++) {
          var vi = i * 4;
          var vj = j * 4;
          area += verts[vi + 0] * verts[vj + 2] - verts[vj + 0] * verts[vi + 2];
        }

        return (area + 1) / 2;
      }
    }, {
      key: "intersectSegCountour",
      value: function intersectSegCountour(d0, d1, i, n, verts, d0verts, d1verts) {
        // For each edge (k,k+1) of P
        var pverts = new Array(4 * 4);

        for (var g = 0; g < 4; g++) {
          pverts[g] = d0verts[d0 + g];
          pverts[4 + g] = d1verts[d1 + g];
        }

        d0 = 0;
        d1 = 4;

        for (var k = 0; k < n; k++) {
          var k1 = RecastMesh.next(k, n); // Skip edges incident to i.

          if (i == k || i == k1) continue;
          var p0 = k * 4;
          var p1 = k1 * 4;

          for (var _g = 0; _g < 4; _g++) {
            pverts[8 + _g] = verts[p0 + _g];
            pverts[12 + _g] = verts[p1 + _g];
          }

          p0 = 8;
          p1 = 12;
          if (RecastMesh.vequal(pverts, d0, p0) || RecastMesh.vequal(pverts, d1, p0) || RecastMesh.vequal(pverts, d0, p1) || RecastMesh.vequal(pverts, d1, p1)) continue;
          if (RecastMesh.intersect(pverts, d0, d1, p0, p1)) return true;
        }

        return false;
      }
    }, {
      key: "inCone",
      value: function inCone(i, n, verts, pj, vertpj) {
        pi = i * 4;
        pi1 = RecastMesh.next(i, n) * 4;
        pin1 = RecastMesh.prev(i, n) * 4;
        pverts = new Array(4 * 4);

        for (var g = 0; g < 4; g++) {
          pverts[g] = verts[pi + g];
          pverts[4 + g] = verts[pi1 + g];
          pverts[8 + g] = verts[pin1 + g];
          pverts[12 + g] = vertpj[pj + g];
        }

        pi = 0;
        pi1 = 4;
        pin1 = 8;
        pj = 12; // If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].

        if (RecastMesh.leftOn(pverts, pin1, pi, pi1)) return RecastMesh.left(pverts, pi, pj, pin1) && RecastMesh.left(pverts, pj, pi, pi1); // Assume (i-1,i,i+1) not collinear.
        // else P[i] is reflex.

        return !(RecastMesh.leftOn(pverts, pi, pj, pi1) && RecastMesh.leftOn(pverts, pj, pi, pin1));
      }
    }, {
      key: "removeDegenerateSegments",
      value: function removeDegenerateSegments(simplified) {
        // Remove adjacent vertices which are equal on xz-plane,
        // or else the triangulator will get confused.
        var npts = simplified.length / 4;

        for (var i = 0; i < npts; ++i) {
          var ni = RecastMesh.next(i, npts); //			if (vequal(&simplified[i*4], &simplified[ni*4]))

          if (simplified[i * 4] == simplified[ni * 4] && simplified[i * 4 + 2] == simplified[ni * 4 + 2]) {
            // Degenerate segment, remove.
            simplified.splice(i * 4, 1);
            simplified.splice(i * 4, 1);
            simplified.splice(i * 4, 1);
            simplified.splice(i * 4, 1);
            npts--;
          }
        }
      }
    }, {
      key: "mergeContours",
      value: function mergeContours(ca, cb, ia, ib) {
        var maxVerts = ca.nverts + cb.nverts + 2;
        var verts = new Array(maxVerts * 4);
        var nv = 0; // Copy contour A.

        for (var i = 0; i <= ca.nverts; ++i) {
          var dst = nv * 4;
          var src = (ia + i) % ca.nverts * 4;
          verts[dst + 0] = ca.verts[src + 0];
          verts[dst + 1] = ca.verts[src + 1];
          verts[dst + 2] = ca.verts[src + 2];
          verts[dst + 3] = ca.verts[src + 3];
          nv++;
        } // Copy contour B


        for (var _i8 = 0; _i8 <= cb.nverts; ++_i8) {
          var _dst = nv * 4;

          var _src = (ib + _i8) % cb.nverts * 4;

          verts[_dst + 0] = cb.verts[_src + 0];
          verts[_dst + 1] = cb.verts[_src + 1];
          verts[_dst + 2] = cb.verts[_src + 2];
          verts[_dst + 3] = cb.verts[_src + 3];
          nv++;
        }

        ca.verts = verts;
        ca.nverts = nv;
        cb.verts = null;
        cb.nverts = 0;
      } // Finds the lowest leftmost vertex of a contour.

    }, {
      key: "findLeftMostVertex",
      value: function findLeftMostVertex(contour) {
        var minx = contour.verts[0];
        var minz = contour.verts[2];
        var leftmost = 0;

        for (var i = 1; i < contour.nverts; i++) {
          var x = contour.verts[i * 4 + 0];
          var z = contour.verts[i * 4 + 2];

          if (x < minx || x == minx && z < minz) {
            minx = x;
            minz = z;
            leftmost = i;
          }
        }

        return [minx, minz, leftmost];
      }
    }, {
      key: "buildContours",
      value: function buildContours(ctx, chf, maxError, maxEdgeLen, buildFlags) {
        var w = chf.width;
        var h = chf.height;
        var borderSize = chf.borderSize;
        var cset = new ContourSet();
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
        var flags = new Array(chf.spanCount).fill(0);
        ctx.startTimer("BUILD_CONTOURS_TRACE"); // Mark boundaries.

        for (var y = 0; y < h; ++y) {
          for (var x = 0; x < w; ++x) {
            var c = chf.cells[x + y * w];

            for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
              // if(y == 3 && x == 3 && i == 1344)
              // 	console.log("earlier");
              var res = 0;
              var s = chf.spans[i];

              if (chf.spans[i].reg == 0 || (chf.spans[i].reg & RecastConstants.RC_BORDER_REG) != 0) {
                flags[i] = 0;
                continue;
              }

              for (var dir = 0; dir < 4; ++dir) {
                var r = 0;

                if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
                  var ax = x + RecastCommon.GetDirOffsetX(dir);
                  var ay = y + RecastCommon.GetDirOffsetY(dir);
                  var ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);
                  r = chf.spans[ai].reg;
                }

                if (r == chf.spans[i].reg) res |= 1 << dir;
              }

              flags[i] = res ^ 0xf; // Inverse, mark non connected edges.
            }
          }
        }

        ctx.stopTimer("BUILD_CONTOURS_TRACE");
        var verts = new Array(256);
        var simplified = new Array(64);

        for (var _y = 0; _y < h; ++_y) {
          for (var _x = 0; _x < w; ++_x) {
            var _c = chf.cells[_x + _y * w];

            for (var _i9 = _c.index, _ni = _c.index + _c.count; _i9 < _ni; ++_i9) {
              // if(y==3 && x == 3 && i == 1344)
              // 	console.log("dak");
              if (flags[_i9] == 0 || flags[_i9] == 0xf) {
                flags[_i9] = 0;
                continue;
              }

              var reg = chf.spans[_i9].reg;
              if (reg == 0 || (reg & RecastConstants.RC_BORDER_REG) != 0) continue;
              var area = chf.areas[_i9];
              verts = [];
              simplified = [];
              ctx.startTimer("BUILD_CONTOURS_TRACE");
              RecastContour.walkContour(_x, _y, _i9, chf, flags, verts);
              ctx.stopTimer("BUILD_CONTOURS_TRACE");
              ctx.startTimer("BUILD_CONTOURS_SIMPLIFY");
              RecastContour.simplifyContour(verts, simplified, maxError, maxEdgeLen, buildFlags);
              RecastContour.removeDegenerateSegments(simplified);
              ctx.stopTimer("BUILD_CONTOURS_SIMPLIFY"); // Store region=>contour remap info.
              // Create contour.

              if (simplified.length / 4 >= 3) {
                var cont = new Contour();
                cset.conts.push(cont);
                cont.nverts = simplified.length / 4;
                cont.verts = new Array(simplified.length);

                for (var l = 0; l < cont.verts.length; l++) {
                  cont.verts[l] = simplified[l];
                }

                if (borderSize > 0) {
                  // If the heightfield was build with bordersize, remove the offset.
                  for (var j = 0; j < cont.nverts; ++j) {
                    cont.verts[j * 4] -= borderSize;
                    cont.verts[j * 4 + 2] -= borderSize;
                  }
                }

                cont.nrverts = verts.length / 4;
                cont.rverts = new Array(verts.length);

                for (var _l = 0; _l < cont.rverts.length; _l++) {
                  cont.rverts[_l] = verts[_l];
                }

                if (borderSize > 0) {
                  // If the heightfield was build with bordersize, remove the offset.
                  for (var _j2 = 0; _j2 < cont.nrverts; ++_j2) {
                    cont.rverts[_j2 * 4] -= borderSize;
                    cont.rverts[_j2 * 4 + 2] -= borderSize;
                  }
                }

                cont.reg = reg;
                cont.area = area;
              }
            }
          }
        } // Merge holes if needed.


        if (cset.conts.length > 0) {
          // Calculate winding of all polygons.
          var winding = new Array(cset.conts.length);
          var nholes = 0;

          for (var _i10 = 0; _i10 < cset.conts.length; ++_i10) {
            var _cont = cset.conts[_i10]; // If the contour is wound backwards, it is a hole.

            winding[_i10] = RecastContour.calcAreaOfPolygon2D(_cont.verts, _cont.nverts) < 0 ? -1 : 1;
            if (winding[_i10] < 0) nholes++;
          }

          if (nholes > 0) {
            // Collect outline contour and holes contours per region.
            // We assume that there is one outline and multiple holes.
            var nregions = chf.maxRegions + 1;
            var regions = new Array(nregions);

            for (var _i11 = 0; _i11 < nregions; _i11++) {
              regions[_i11] = new ContourRegion();
            }

            for (var _i12 = 0; _i12 < cset.conts.length; ++_i12) {
              var _cont2 = cset.conts[_i12]; // Positively would contours are outlines, negative holes.

              if (winding[_i12] > 0) {
                if (regions[_cont2.reg].outline != null) {
                  throw new RuntimeException("rcBuildContours: Multiple outlines for region " + _cont2.reg + ".");
                }

                regions[_cont2.reg].outline = _cont2;
              } else {
                regions[_cont2.reg].nholes++;
              }
            }

            for (var _i13 = 0; _i13 < nregions; _i13++) {
              if (regions[_i13].nholes > 0) {
                regions[_i13].holes = new ContourHole[regions[_i13].nholes]();

                for (var nh = 0; nh < regions[_i13].nholes; nh++) {
                  regions[_i13].holes[nh] = new ContourHole();
                }

                regions[_i13].nholes = 0;
              }
            }

            for (var _i14 = 0; _i14 < cset.conts.length; ++_i14) {
              var _cont3 = cset.conts[_i14];
              var _reg = regions[_cont3.reg];
              if (winding[_i14] < 0) _reg.holes[_reg.nholes++].contour = _cont3;
            } // Finally merge each regions holes into the outline.


            for (var _i15 = 0; _i15 < nregions; _i15++) {
              var _reg2 = regions[_i15];
              if (_reg2.nholes == 0) continue;

              if (_reg2.outline != null) {
                mergeRegionHoles(ctx, _reg2);
              } else {
                // The region does not have an outline.
                // This can happen if the contour becaomes selfoverlapping because of
                // too aggressive simplification settings.
                throw new RuntimeException("rcBuildContours: Bad outline for region " + _i15 + ", contour simplification is likely too aggressive.");
              }
            }
          }
        }

        ctx.stopTimer("BUILD_CONTOURS");
        return cset;
      }
    }]);

    return RecastContour;
  }();

  _defineProperty(RecastContour, "ContourRegion", (_temp$2 = function ContourRegion() {
    _classCallCheck(this, ContourRegion);

    _defineProperty(this, "outline", void 0);

    _defineProperty(this, "holes", []);

    _defineProperty(this, "nholes", void 0);
  }, _temp$2));

  _defineProperty(RecastContour, "ContourHole", (_temp2$2 = function ContourHole() {
    _classCallCheck(this, ContourHole);

    _defineProperty(this, "leftmost", void 0);

    _defineProperty(this, "minx", void 0);

    _defineProperty(this, "minz", void 0);

    _defineProperty(this, "contour", void 0);
  }, _temp2$2));

  _defineProperty(RecastContour, "PotentialDiagonal", (_temp3 = function PotentialDiagonal() {
    _classCallCheck(this, PotentialDiagonal);

    _defineProperty(this, "dist", void 0);

    _defineProperty(this, "vert", void 0);
  }, _temp3));

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
  var PolyMeshDetail = function PolyMeshDetail() {
    _classCallCheck(this, PolyMeshDetail);

    _defineProperty(this, "meshes", []);

    _defineProperty(this, "verts", []);

    _defineProperty(this, "tris", []);

    _defineProperty(this, "nmeshes", 0);

    _defineProperty(this, "nverts", 0);

    _defineProperty(this, "ntris", 0);
  };

  function arraycopy$2(one, oneStart, two, twoStart, len) {
    for (var i = 0; i < len; i++) {
      two[twoStart + i] = one[oneStart + i];
    }
  }

  var HeightPatch = function HeightPatch() {
    _classCallCheck(this, HeightPatch);

    _defineProperty(this, "xmin", 0);

    _defineProperty(this, "ymin", 0);

    _defineProperty(this, "width", 0);

    _defineProperty(this, "height", 0);

    _defineProperty(this, "data", []);
  };

  var RecastMeshDetail = /*#__PURE__*/function () {
    function RecastMeshDetail() {
      _classCallCheck(this, RecastMeshDetail);
    }

    _createClass(RecastMeshDetail, [{
      key: "mergePolyMeshDetails",
      /// @see rcAllocPolyMeshDetail, rcPolyMeshDetail
      value: function mergePolyMeshDetails(ctx, meshes, nmeshes) {
        mesh = new PolyMeshDetail();
        ctx.startTimer("MERGE_POLYMESHDETAIL");
        var maxVerts = 0;
        var maxTris = 0;
        var maxMeshes = 0;

        for (var i = 0; i < nmeshes; ++i) {
          if (meshes[i] == null) continue;
          maxVerts += meshes[i].nverts;
          maxTris += meshes[i].ntris;
          maxMeshes += meshes[i].nmeshes;
        }

        mesh.nmeshes = 0;
        mesh.meshes = new Array(maxMeshes * 4);
        mesh.ntris = 0;
        mesh.tris = new Array(maxTris * 4);
        mesh.nverts = 0;
        mesh.verts = new Array(maxVerts * 3); // Merge datas.

        for (var _i = 0; _i < nmeshes; ++_i) {
          var dm = meshes[_i];
          if (dm == null) continue;

          for (var j = 0; j < dm.nmeshes; ++j) {
            var dst = mesh.nmeshes * 4;
            var src = j * 4;
            mesh.meshes[dst + 0] = mesh.nverts + dm.meshes[src + 0];
            mesh.meshes[dst + 1] = dm.meshes[src + 1];
            mesh.meshes[dst + 2] = mesh.ntris + dm.meshes[src + 2];
            mesh.meshes[dst + 3] = dm.meshes[src + 3];
            mesh.nmeshes++;
          }

          for (var k = 0; k < dm.nverts; ++k) {
            RecastVectors$1.copy(mesh.verts, mesh.nverts * 3, dm.verts, k * 3);
            mesh.nverts++;
          }

          for (var _k = 0; _k < dm.ntris; ++_k) {
            mesh.tris[mesh.ntris * 4 + 0] = dm.tris[_k * 4 + 0];
            mesh.tris[mesh.ntris * 4 + 1] = dm.tris[_k * 4 + 1];
            mesh.tris[mesh.ntris * 4 + 2] = dm.tris[_k * 4 + 2];
            mesh.tris[mesh.ntris * 4 + 3] = dm.tris[_k * 4 + 3];
            mesh.ntris++;
          }
        }

        ctx.stopTimer("MERGE_POLYMESHDETAIL");
        return mesh;
      }
    }], [{
      key: "vdistSq2_3",
      // Max tris for delaunay is 2n-2-k (n=num verts, k=num hull verts).
      value: function vdistSq2_3(verts, p, q) {
        var dx = verts[q + 0] - verts[p + 0];
        var dy = verts[q + 2] - verts[p + 2];
        return dx * dx + dy * dy;
      }
    }, {
      key: "vdist2",
      value: function vdist2(verts, p, q) {
        return Math.sqrt(RecastMeshDetail.vdistSq2_3(verts, p, q));
      }
    }, {
      key: "vdistSq2_2",
      value: function vdistSq2_2(p, q) {
        var dx = q[0] - p[0];
        var dy = q[2] - p[2];
        return dx * dx + dy * dy;
      }
    }, {
      key: "vdist2",
      value: function vdist2(p, q) {
        return Math.sqrt(RecastMeshDetail.vdistSq2_2(p, q));
      }
    }, {
      key: "vdistSq2",
      value: function vdistSq2(p, verts, q) {
        dx = verts[q + 0] - p[0];
        dy = verts[q + 2] - p[2];
        return dx * dx + dy * dy;
      }
    }, {
      key: "vdist2",
      value: function vdist2(p, verts, q) {
        return Math.sqrt(RecastMeshDetail.vdistSq2_3(p, verts, q));
      }
    }, {
      key: "vcross2",
      value: function vcross2(verts, p1, p2, p3) {
        u1 = verts[p2 + 0] - verts[p1 + 0];
        v1 = verts[p2 + 2] - verts[p1 + 2];
        u2 = verts[p3 + 0] - verts[p1 + 0];
        v2 = verts[p3 + 2] - verts[p1 + 2];
        return u1 * v2 - v1 * u2;
      }
    }, {
      key: "vcross2",
      value: function vcross2(p1, p2, p3) {
        u1 = p2[0] - p1[0];
        v1 = p2[2] - p1[2];
        u2 = p3[0] - p1[0];
        v2 = p3[2] - p1[2];
        return u1 * v2 - v1 * u2;
      }
    }, {
      key: "circumCircle",
      value: function circumCircle(verts, p1, p2, p3, c, r) {
        EPS = 1e-6; // Calculate the circle relative to p1, to asome precision issues.

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
    }, {
      key: "distPtTri",
      value: function distPtTri(p, verts, a, b, c) {
        var v0 = new Array(3);
        var v1 = new Array(3);
        var v2 = new Array(3);
        RecastVectors$1.subA(v0, verts, c, a);
        RecastVectors$1.subA(v1, verts, b, a);
        RecastVectors$1.subB(v2, p, verts, a);
        var dot00 = RecastMeshDetail.vdot2(v0, v0);
        var dot01 = RecastMeshDetail.vdot2(v0, v1);
        var dot02 = RecastMeshDetail.vdot2(v0, v2);
        var dot11 = RecastMeshDetail.vdot2(v1, v1);
        var dot12 = RecastMeshDetail.vdot2(v1, v2); // Compute barycentric coordinates

        var invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
        var u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        var v = (dot00 * dot12 - dot01 * dot02) * invDenom; // If poPoly lies inside the triangle, return interpolated y-coord.

        var EPS = 1e-4;

        if (u >= -EPS && v >= -EPS && u + v <= 1 + EPS) {
          var y = verts[a + 1] + v0[1] * u + v1[1] * v;
          return Math.abs(y - p[1]);
        }

        return Number.MAX_VALUE;
      }
    }, {
      key: "distancePtSeg",
      value: function distancePtSeg(verts, pt, p, q) {
        var pqx = verts[q + 0] - verts[p + 0];
        var pqy = verts[q + 1] - verts[p + 1];
        var pqz = verts[q + 2] - verts[p + 2];
        var dx = verts[pt + 0] - verts[p + 0];
        var dy = verts[pt + 1] - verts[p + 1];
        var dz = verts[pt + 2] - verts[p + 2];
        var d = pqx * pqx + pqy * pqy + pqz * pqz;
        var t = pqx * dx + pqy * dy + pqz * dz;
        if (d > 0) t /= d;
        if (t < 0) t = 0;else if (t > 1) t = 1;
        dx = verts[p + 0] + t * pqx - verts[pt + 0];
        dy = verts[p + 1] + t * pqy - verts[pt + 1];
        dz = verts[p + 2] + t * pqz - verts[pt + 2];
        return dx * dx + dy * dy + dz * dz;
      }
    }, {
      key: "distancePtSeg2d",
      value: function distancePtSeg2d(verts, pt, poly, p, q) {
        var pqx = poly[q + 0] - poly[p + 0];
        var pqz = poly[q + 2] - poly[p + 2];
        var dx = verts[pt + 0] - poly[p + 0];
        var dz = verts[pt + 2] - poly[p + 2];
        var d = pqx * pqx + pqz * pqz;
        var t = pqx * dx + pqz * dz;
        if (d > 0) t /= d;
        if (t < 0) t = 0;else if (t > 1) t = 1;
        dx = poly[p + 0] + t * pqx - verts[pt + 0];
        dz = poly[p + 2] + t * pqz - verts[pt + 2];
        return dx * dx + dz * dz;
      }
    }, {
      key: "distToTriMesh",
      value: function distToTriMesh(p, verts, nverts, tris, ntris) {
        var dmin = Number.MAX_VALUE;

        for (var i = 0; i < ntris; ++i) {
          var va = tris[i * 4 + 0] * 3;
          var vb = tris[i * 4 + 1] * 3;
          var vc = tris[i * 4 + 2] * 3;

          var _d = RecastMeshDetail.distPtTri(p, verts, va, vb, vc);

          if (_d < dmin) dmin = _d;
        }

        if (dmin == Number.MAX_VALUE) return -1;
        return dmin;
      }
    }, {
      key: "distToPoly",
      value: function distToPoly(nvert, verts, p) {
        var dmin = Number.MAX_VALUE;
        var c = false;

        for (var _i2 = 0, _j = nvert - 1; _i2 < nvert; _j = _i2++) {
          var vi = _i2 * 3;
          var vj = _j * 3;
          if (verts[vi + 2] > p[2] != verts[vj + 2] > p[2] && p[0] < (verts[vj + 0] - verts[vi + 0]) * (p[2] - verts[vi + 2]) / (verts[vj + 2] - verts[vi + 2]) + verts[vi + 0]) c = !c;
          dmin = Math.min(dmin, RecastMeshDetail.distancePtSeg2d(p, 0, verts, vj, vi));
        }

        return c ? -dmin : dmin;
      }
    }, {
      key: "getHeight",
      value: function getHeight(fx, fy, fz, cs, ics, ch, radius, hp) {
        var ix = Math.floor(fx * ics + 0.01);
        var iz = Math.floor(fz * ics + 0.01);
        ix = RecastCommon.clamp(ix - hp.xmin, 0, hp.width - 1);
        iz = RecastCommon.clamp(iz - hp.ymin, 0, hp.height - 1);
        var h = hp.data[ix + iz * hp.width];

        if (h == RecastMeshDetail.RC_UNSET_HEIGHT) {
          // Special case when data might be bad.
          // Walk adjacent cells _in a spiral up to 'radius', and look
          // for a pixel which has a valid height.
          var x = 1,
              z = 0,
              _dx = 1,
              dz = 0;
          var maxSize = radius * 2 + 1;
          var maxIter = maxSize * maxSize - 1;
          var nextRingIterStart = 8;
          var nextRingIters = 16;
          dmin = Number.MAX_VALUE;

          for (var i = 0; i < maxIter; ++i) {
            var nx = ix + x;
            var nz = iz + z;

            if (nx >= 0 && nz >= 0 && nx < hp.width && nz < hp.height) {
              var nh = hp.data[nx + nz * hp.width];

              if (nh != RecastMeshDetail.RC_UNSET_HEIGHT) {
                d = Math.abs(nh * ch - fy);

                if (d < dmin) {
                  h = nh;
                  dmin = d;
                }
              }
            } // We are searching _in a grid which looks approximately like this:
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
              if (h != RecastMeshDetail.RC_UNSET_HEIGHT) break;
              nextRingIterStart += nextRingIters;
              nextRingIters += 8;
            }

            if (x == z || x < 0 && x == -z || x > 0 && x == 1 - z) {
              var tmp = _dx;
              _dx = -dz;
              dz = tmp;
            }

            x += _dx;
            z += dz;
          }
        }

        return h;
      }
    }, {
      key: "findEdge",
      value: function findEdge(edges, s, t) {
        for (var i = 0; i < edges.length / 4; i++) {
          var e = i * 4;
          if (edges[e + 0] == s && edges[e + 1] == t || edges[e + 0] == t && edges[e + 1] == s) return i;
        }

        return EV_UNDEF;
      }
    }, {
      key: "addEdge",
      value: function addEdge(ctx, edges, maxEdges, s, t, l, r) {
        if (edges.length / 4 >= maxEdges) {
          throw new RuntimeException("addEdge: Too many edges (" + edges.length / 4 + "/" + maxEdges + ").");
        } // Add edge if not already _in the triangulation.


        var e = findEdge(edges, s, t);

        if (e == EV_UNDEF) {
          edges.push(s);
          edges.push(t);
          edges.push(l);
          edges.push(r);
        }
      }
    }, {
      key: "updateLeftFace",
      value: function updateLeftFace(edges, e, s, t, f) {
        if (edges[e + 0] == s && edges[e + 1] == t && edges[e + 2] == EV_UNDEF) edges.set(e + 2, f);else if (edges[e + 1] == s && edges[e + 0] == t && edges[e + 3] == EV_UNDEF) edges.set(e + 3, f);
      }
    }, {
      key: "overlapSegSeg2d",
      value: function overlapSegSeg2d(verts, a, b, c, d) {
        a1 = vcross2(verts, a, b, d);
        a2 = vcross2(verts, a, b, c);

        if (a1 * a2 < 0.0) {
          a3 = vcross2(verts, c, d, a);
          a4 = a3 + a2 - a1;
          if (a3 * a4 < 0.0) return true;
        }

        return false;
      }
    }, {
      key: "overlapEdges",
      value: function overlapEdges(pts, edges, s1, t1) {
        for (var i = 0; i < edges.length / 4; ++i) {
          var s0 = edges[i * 4 + 0];
          var t0 = edges[i * 4 + 1]; // Same or connected edges do not overlap.

          if (s0 == s1 || s0 == t1 || t0 == s1 || t0 == t1) continue;
          if (overlapSegSeg2d(pts, s0 * 3, t0 * 3, s1 * 3, t1 * 3)) return true;
        }

        return false;
      }
    }, {
      key: "completeFacet",
      value: function completeFacet(ctx, pts, npts, edges, maxEdges, nfaces, e) {
        EPS = 1e-5;
        var edge = e * 4; // Cache s and t.

        var s, t;

        if (edges[edge + 2] == EV_UNDEF) {
          s = edges[edge + 0];
          t = edges[edge + 1];
        } else if (edges[edge + 3] == EV_UNDEF) {
          s = edges[edge + 1];
          t = edges[edge + 0];
        } else {
          // Edge already completed.
          return nfaces;
        } // Find best poPoly on left of edge.


        var pt = npts;
        var c = new Array(3);
        var r = new AtomicReference(-1);

        for (var u = 0; u < npts; ++u) {
          if (u == s || u == t) continue;

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
              if (overlapEdges(pts, edges, s, u)) continue;
              if (overlapEdges(pts, edges, t, u)) continue; // Edge is valid.

              pt = u;
              circumCircle(pts, s * 3, t * 3, u * 3, c, r);
            }
          }
        } // Add new triangle or update edge info if s-t is on hull.


        if (pt < npts) {
          // Update face information of edge being compPolyed.
          updateLeftFace(edges, e * 4, s, t, nfaces); // Add new edge or update face info of old edge.

          e = findEdge(edges, pt, s);
          if (e == EV_UNDEF) addEdge(ctx, edges, maxEdges, pt, s, nfaces, EV_UNDEF);else updateLeftFace(edges, e * 4, pt, s, nfaces); // Add new edge or update face info of old edge.

          e = findEdge(edges, t, pt);
          if (e == EV_UNDEF) addEdge(ctx, edges, maxEdges, t, pt, nfaces, EV_UNDEF);else updateLeftFace(edges, e * 4, t, pt, nfaces);
          nfaces++;
        } else {
          updateLeftFace(edges, e * 4, s, t, EV_HULL);
        }

        return nfaces;
      }
    }, {
      key: "delaunayHull",
      value: function delaunayHull(ctx, npts, pts, nhull, hull, tris) {
        var nfaces = 0;
        var maxEdges = npts * 10;
        var edges = new Array(64);

        for (var i = 0, j = nhull - 1; i < nhull; j = i++) {
          addEdge(ctx, edges, maxEdges, hull[j], hull[i], EV_HULL, EV_UNDEF);
        }

        var currentEdge = 0;

        while (currentEdge < edges.length / 4) {
          if (edges[currentEdge * 4 + 2] == EV_UNDEF) {
            nfaces = completeFacet(ctx, pts, npts, edges, maxEdges, nfaces, currentEdge);
          }

          if (edges[currentEdge * 4 + 3] == EV_UNDEF) {
            nfaces = completeFacet(ctx, pts, npts, edges, maxEdges, nfaces, currentEdge);
          }

          currentEdge++;
        } // Create tris


        tris = [];

        for (var _i3 = 0; _i3 < nfaces * 4; ++_i3) {
          tris.push(-1);
        }

        for (var _i4 = 0; _i4 < edges.length / 4; ++_i4) {
          var e = _i4 * 4;

          if (edges[e + 3] >= 0) {
            // Left face
            var t = edges[e + 3] * 4;

            if (tris[t + 0] == -1) {
              tris.set(t + 0, edges[e + 0]);
              tris.set(t + 1, edges[e + 1]);
            } else if (tris[t + 0] == edges[e + 1]) tris.set(t + 2, edges[e + 0]);else if (tris[t + 1] == edges[e + 0]) tris.set(t + 2, edges[e + 1]);
          }

          if (edges[e + 2] >= 0) {
            // Right
            var _t = edges[e + 2] * 4;

            if (tris[_t + 0] == -1) {
              tris.set(_t + 0, edges[e + 1]);
              tris.set(_t + 1, edges[e + 0]);
            } else if (tris[_t + 0] == edges[e + 0]) tris.set(_t + 2, edges[e + 1]);else if (tris[_t + 1] == edges[e + 1]) tris.set(_t + 2, edges[e + 0]);
          }
        }

        for (var _i5 = 0; _i5 < tris.length / 4; ++_i5) {
          var _t2 = _i5 * 4;

          if (tris[_t2 + 0] == -1 || tris[_t2 + 1] == -1 || tris[_t2 + 2] == -1) {
            System.err.println("Dangling! " + tris[_t2] + " " + tris[_t2 + 1] + "  " + tris[_t2 + 2]); //ctx.log(RC_LOG_WARNING, "delaunayHull: Removing dangling face %d [%d,%d,%d].", i, t[0],t[1],t[2]);

            tris.set(_t2 + 0, tris[tris.length - 4]);
            tris.set(_t2 + 1, tris[tris.length - 3]);
            tris.set(_t2 + 2, tris[tris.length - 2]);
            tris.set(_t2 + 3, tris[tris.length - 1]);
            tris.remove(tris.length - 1);
            tris.remove(tris.length - 1);
            tris.remove(tris.length - 1);
            tris.remove(tris.length - 1);
            --_i5;
          }
        }
      } // Calculate minimum extend of the polygon.

    }, {
      key: "polyMinExtent",
      value: function polyMinExtent(verts, nverts) {
        var minDist = Number.MAX_VALUE;

        for (var i = 0; i < nverts; i++) {
          var ni = (i + 1) % nverts;
          var p1 = i * 3;
          var p2 = ni * 3;
          var maxEdgeDist = 0;

          for (var j = 0; j < nverts; j++) {
            if (j == i || j == ni) continue;

            var _d2 = RecastMeshDetail.distancePtSeg2d(verts, j * 3, verts, p1, p2);

            maxEdgeDist = Math.max(maxEdgeDist, _d2);
          }

          minDist = Math.min(minDist, maxEdgeDist);
        }

        return Math.sqrt(minDist);
      }
    }, {
      key: "triangulateHull",
      value: function triangulateHull(nverts, verts, nhull, hull, tris) {
        var start = 0,
            left = 1,
            right = nhull - 1; // Start from an ear with shortest perimeter.
        // This tends to favor well formed triangles as starting point.

        var dmin = 0;

        for (var i = 0; i < nhull; i++) {
          var pi = RecastMesh.prev(i, nhull);
          var ni = RecastMesh.next(i, nhull);
          var pv = hull[pi] * 3;
          var cv = hull[i] * 3;
          var nv = hull[ni] * 3;

          var _d3 = RecastMeshDetail.vdist2(verts, pv, cv) + RecastMeshDetail.vdist2(verts, cv, nv) + RecastMeshDetail.vdist2(verts, nv, pv);

          if (_d3 < dmin) {
            start = i;
            left = ni;
            right = pi;
            dmin = _d3;
          }
        } // Add first triangle


        tris.push(hull[start]);
        tris.push(hull[left]);
        tris.push(hull[right]);
        tris.push(0); // Triangulate the polygon by moving left or right,
        // depending on which triangle has shorter perimeter.
        // This heuristic was chose emprically, since it seems
        // handle tesselated straight edges well.

        while (RecastMesh.next(left, nhull) != right) {
          // Check to see if se should advance left or right.
          var nleft = RecastMesh.next(left, nhull);
          var nright = RecastMesh.prev(right, nhull);
          var cvleft = hull[left] * 3;
          var nvleft = hull[nleft] * 3;
          var cvright = hull[right] * 3;
          var nvright = hull[nright] * 3;
          var dleft = RecastMeshDetail.vdist2(verts, cvleft, nvleft) + RecastMeshDetail.vdist2(verts, nvleft, cvright);
          var dright = RecastMeshDetail.vdist2(verts, cvright, nvright) + RecastMeshDetail.vdist2(verts, cvleft, nvright);

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
    }, {
      key: "getJitterX",
      value: function getJitterX(i) {
        return (i * 0x8da6b343 & 0xffff) / 65535.0 * 2.0 - 1.0;
      }
    }, {
      key: "getJitterY",
      value: function getJitterY(i) {
        return (i * 0xd8163841 & 0xffff) / 65535.0 * 2.0 - 1.0;
      }
    }, {
      key: "buildPolyDetail",
      value: function buildPolyDetail(ctx, _in, nin, sampleDist, sampleMaxError, heightSearchRadius, chf, hp, verts, tris) {
        var samples = [];
        var nverts = 0;
        var edge = new Array((RecastMeshDetail.MAX_VERTS_PER_EDGE + 1) * 3);
        var hull = new Array(RecastMeshDetail.MAX_VERTS);
        var nhull = 0;
        nverts = nin;

        for (var i = 0; i < nin; ++i) {
          RecastVectors$1.copy4(verts, i * 3, _in, i * 3);
        }

        tris = [];
        var cs = chf.cs;
        var ics = 1.0 / cs; // Calculate minimum extents of the polygon based on input data.

        var minExtent = RecastMeshDetail.polyMinExtent(verts, nverts); // Tessellate outlines.
        // This is done _in separate pass _in order to ensure
        // seamless height values across the ply boundaries.

        if (sampleDist > 0) {
          for (var _i6 = 0, j = nin - 1; _i6 < nin; j = _i6++) {
            var vj = j * 3;
            var vi = _i6 * 3;
            var swapped = false; // Make sure the segments are always handled _in same order
            // using lexological sort or else there will be seams.

            if (Math.abs(_in[vj + 0] - _in[vi + 0]) < 1e-6) {
              if (_in[vj + 2] > _in[vi + 2]) {
                var temp = vi;
                vi = vj;
                vj = temp;
                swapped = true;
              }
            } else {
              if (_in[vj + 0] > _in[vi + 0]) {
                var _temp = vi;
                vi = vj;
                vj = _temp;
                swapped = true;
              }
            } // Create samples alet the edge.


            var _dx2 = _in[vi + 0] - _in[vj + 0];

            var _dy = _in[vi + 1] - _in[vj + 1];

            var dz = _in[vi + 2] - _in[vj + 2];

            var _d4 = Math.sqrt(_dx2 * _dx2 + dz * dz);

            var nn = 1 + Math.floor(_d4 / sampleDist);
            if (nn >= RecastMeshDetail.MAX_VERTS_PER_EDGE) nn = RecastMeshDetail.MAX_VERTS_PER_EDGE - 1;
            if (nverts + nn >= RecastMeshDetail.MAX_VERTS) nn = RecastMeshDetail.MAX_VERTS - 1 - nverts;

            for (var k = 0; k <= nn; ++k) {
              var u = k / nn;
              var pos = k * 3;
              edge[pos + 0] = _in[vj + 0] + _dx2 * u;
              edge[pos + 1] = _in[vj + 1] + _dy * u;
              edge[pos + 2] = _in[vj + 2] + dz * u;
              edge[pos + 1] = RecastMeshDetail.getHeight(edge[pos + 0], edge[pos + 1], edge[pos + 2], cs, ics, chf.ch, heightSearchRadius, hp) * chf.ch;
            } // Simplify samples.


            var idx = new Array(RecastMeshDetail.MAX_VERTS_PER_EDGE);
            idx[0] = 0;
            idx[1] = nn;
            var nidx = 2;

            for (var _k2 = 0; _k2 < nidx - 1;) {
              var a = idx[_k2];
              var b = idx[_k2 + 1];
              var va = a * 3;
              var vb = b * 3; // Find maximum deviation alet the segment.

              var maxd = 0;
              var maxi = -1;

              for (var m = a + 1; m < b; ++m) {
                var dev = RecastMeshDetail.distancePtSeg(edge, m * 3, va, vb);

                if (dev > maxd) {
                  maxd = dev;
                  maxi = m;
                }
              } // If the max deviation is larger than accepted error,
              // add new point, else continue to next segment.


              if (maxi != -1 && maxd > sampleMaxError * sampleMaxError) {
                for (var _m = nidx; _m > _k2; --_m) {
                  idx[_m] = idx[_m - 1];
                }

                idx[_k2 + 1] = maxi;
                nidx++;
              } else {
                ++_k2;
              }
            }

            hull[nhull++] = j; // Add new vertices.

            if (swapped) {
              for (var _k3 = nidx - 2; _k3 > 0; --_k3) {
                RecastVectors$1.copy(verts, nverts * 3, edge, idx[_k3] * 3);
                hull[nhull++] = nverts;
                nverts++;
              }
            } else {
              for (var _k4 = 1; _k4 < nidx - 1; ++_k4) {
                RecastVectors$1.copy(verts, nverts * 3, edge, idx[_k4] * 3);
                hull[nhull++] = nverts;
                nverts++;
              }
            }
          }
        } // If the polygon minimum extent is small (sliver or small triangle), do not try to add internal points.


        if (minExtent < sampleDist * 2) {
          RecastMeshDetail.triangulateHull(nverts, verts, nhull, hull, tris);
          return [nverts, tris];
        } // Tessellate the base mesh.
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
          var bmin = new Array(3);
          var bmax = new Array(3);
          RecastVectors$1.copy3(bmin, _in, 0);
          RecastVectors$1.copy3(bmax, _in, 0);

          for (var _i7 = 1; _i7 < nin; ++_i7) {
            RecastVectors$1.min(bmin, _in, _i7 * 3);
            RecastVectors$1.max(bmax, _in, _i7 * 3);
          }

          var x0 = Math.floor(bmin[0] / sampleDist);
          var x1 = Math.ceil(bmax[0] / sampleDist);
          var z0 = Math.floor(bmin[2] / sampleDist);
          var z1 = Math.ceil(bmax[2] / sampleDist);
          samples = [];

          for (var z = z0; z < z1; ++z) {
            for (var x = x0; x < x1; ++x) {
              var pt = new Array(3);
              pt[0] = x * sampleDist;
              pt[1] = (bmax[1] + bmin[1]) * 0.5;
              pt[2] = z * sampleDist; // Make sure the samples are not too close to the edges.

              if (RecastMeshDetail.distToPoly(nin, _in, pt) > -sampleDist / 2) continue;
              samples.push(x);
              samples.push(RecastMeshDetail.getHeight(pt[0], pt[1], pt[2], cs, ics, chf.ch, heightSearchRadius, hp));
              samples.push(z);
              samples.push(0); // Not added
            }
          } // Add the samples starting from the one that has the most
          // error. The procedure stops when all samples are added
          // or when the max error is within treshold.


          var nsamples = samples.length / 4;

          for (var iter = 0; iter < nsamples; ++iter) {
            if (nverts >= RecastMeshDetail.MAX_VERTS) break; // Find sample with most error.

            var bestpt = new Array(3);
            var bestd = 0;
            var besti = -1;

            for (var _i8 = 0; _i8 < nsamples; ++_i8) {
              var s = _i8 * 4;
              if (samples[s + 3] != 0) continue; // skip added.

              var _pt = new Array(3); // The sample location is jittered to get rid of some bad triangulations
              // which are cause by symmetrical data from the grid structure.


              _pt[0] = samples[s + 0] * sampleDist + RecastMeshDetail.getJitterX(_i8) * cs * 0.1;
              _pt[1] = samples[s + 1] * chf.ch;
              _pt[2] = samples[s + 2] * sampleDist + RecastMeshDetail.getJitterY(_i8) * cs * 0.1;

              var _d5 = RecastMeshDetail.distToTriMesh(_pt, verts, nverts, tris, tris.length / 4);

              if (_d5 < 0) continue; // did not hit the mesh.

              if (_d5 > bestd) {
                bestd = _d5;
                besti = _i8;
                bestpt = _pt;
              }
            } // If the max error is within accepted threshold, stop tesselating.


            if (bestd <= sampleMaxError || besti == -1) break; // Mark sample as added.

            samples.set(besti * 4 + 3, 1); // Add the new sample point.

            RecastVectors$1.copy(verts, nverts * 3, bestpt, 0);
            nverts++; // Create new triangulation.
            // TODO: Incremental add instead of full rebuild.

            delaunayHull(ctx, nverts, verts, nhull, hull, tris);
          }
        }

        var ntris = tris.length / 4;

        if (ntris > RecastMeshDetail.MAX_TRIS) {
          var subList = tris.subList(0, RecastMeshDetail.MAX_TRIS * 4);
          tris = [];
          tris.addAll(subList);
          throw new RuntimeException("rcBuildPolyMeshDetail: Shrinking triangle count from " + ntris + " to max " + RecastMeshDetail.MAX_TRIS);
        }

        return [nverts, tris];
      }
    }, {
      key: "seedArrayWithPolyCenter",
      value: function seedArrayWithPolyCenter(ctx, chf, meshpoly, poly, npoly, verts, bs, hp, array) {
        // Note: Reads to the compact heightfield are offset by border size (bs)
        // since border size offset is already removed from the polymesh vertices.
        var offset = [0, 0, -1, -1, 0, -1, 1, -1, 1, 0, 1, 1, 0, 1, -1, 1, -1, 0]; // Find cell closest to a let vertex

        var startCellX = 0,
            startCellY = 0,
            startSpanIndex = -1;
        var dmin = RecastMeshDetail.RC_UNSET_HEIGHT;

        for (var j = 0; j < npoly && dmin > 0; ++j) {
          for (var k = 0; k < 9 && dmin > 0; ++k) {
            var ax = verts[meshpoly[poly + j] * 3 + 0] + offset[k * 2 + 0];
            var ay = verts[meshpoly[poly + j] * 3 + 1];
            var az = verts[meshpoly[poly + j] * 3 + 2] + offset[k * 2 + 1];
            if (ax < hp.xmin || ax >= hp.xmin + hp.width || az < hp.ymin || az >= hp.ymin + hp.height) continue;
            var c = chf.cells[ax + bs + (az + bs) * chf.width];

            for (var i = c.index, ni = c.index + c.count; i < ni && dmin > 0; ++i) {
              var s = chf.spans[i];

              var _d6 = Math.abs(ay - s.y);

              if (_d6 < dmin) {
                startCellX = ax;
                startCellY = az;
                startSpanIndex = i;
                dmin = _d6;
              }
            }
          }
        } // Find center of the polygon


        var pcx = 0,
            pcy = 0;

        for (var _j2 = 0; _j2 < npoly; ++_j2) {
          pcx += verts[meshpoly[poly + _j2] * 3 + 0];
          pcy += verts[meshpoly[poly + _j2] * 3 + 2];
        }

        pcx /= npoly;
        pcy /= npoly;
        array = [];
        array.push(startCellX);
        array.push(startCellY);
        array.push(startSpanIndex);
        var dirs = [0, 1, 2, 3];
        hp.data.fill(0, 0, hp.width * hp.height); // DFS to move to the center. Note that we need a DFS here and can not just move
        // directly towards the center without recording intermediate nodes, even though the polygons
        // are convex. In very rare we can get stuck due to contour simplification if we do not
        // record nodes.

        var cx = -1,
            cy = -1,
            ci = -1;

        while (true) {
          if (array.length < 3) {
            ctx.warn("Walk towards polygon center failed to reach center");
            break;
          }

          ci = array.splice(array.length - 1, 1)[0];
          cy = array.splice(array.length - 1, 1)[0];
          cx = array.splice(array.length - 1, 1)[0]; // Check if close to center of the polygon.

          if (cx == pcx && cy == pcy) {
            break;
          } // If we are already at the correct X-position, prefer direction
          // directly towards the center _in the Y-axis; otherwise prefer
          // direction _in the X-axis


          var directDir = void 0;
          if (cx == pcx) directDir = RecastCommon.rcGetDirForOffset(0, pcy > cy ? 1 : -1);else directDir = RecastCommon.rcGetDirForOffset(pcx > cx ? 1 : -1, 0); // Push the direct dir last so we start with this on next iteration

          var tmp = dirs[3];
          dirs[3] = dirs[directDir];
          dirs[directDir] = tmp;
          var _cs = chf.spans[ci];

          for (var _i9 = 0; _i9 < 4; ++_i9) {
            var dir = dirs[_i9];
            if (RecastCommon.GetCon(_cs, dir) == RecastConstants.RC_NOT_CONNECTED) continue;
            var newX = cx + RecastCommon.GetDirOffsetX(dir);
            var newY = cy + RecastCommon.GetDirOffsetY(dir);
            var hpx = newX - hp.xmin;
            var hpy = newY - hp.ymin;
            if (hpx < 0 || hpx >= hp.width || hpy < 0 || hpy >= hp.height) continue;
            if (hp.data[hpx + hpy * hp.width] != 0) continue;
            hp.data[hpx + hpy * hp.width] = 1;
            array.push(newX);
            array.push(newY);
            array.push(chf.cells[newX + bs + (newY + bs) * chf.width].index + RecastCommon.GetCon(_cs, dir));
          }

          tmp = dirs[3];
          dirs[3] = dirs[directDir];
          dirs[directDir] = tmp;
        }

        array = []; // getHeightData seeds are given _in coordinates with borders

        array.push(cx + bs);
        array.push(cy + bs);
        array.push(ci);
        hp.data.fill(RecastMeshDetail.RC_UNSET_HEIGHT, 0, hp.width * hp.height);
        var cs = chf.spans[ci];
        hp.data[cx - hp.xmin + (cy - hp.ymin) * hp.width] = cs.y;
      }
    }, {
      key: "push3",
      value: function push3(queue, v1, v2, v3) {
        queue.push(v1);
        queue.push(v2);
        queue.push(v3);
      }
    }, {
      key: "getHeightData",
      value: function getHeightData(ctx, chf, meshpolys, poly, npoly, verts, bs, hp, region) {
        // Note: Reads to the compact heightfield are offset by border size (bs)
        // since border size offset is already removed from the polymesh vertices.
        var queue = [];
        hp.data.fill(RecastConstants.RC_UNSET_HEIGHT, 0, hp.width * hp.height);
        var empty = true; // We cannot sample from this let if it was created from polys
        // of different regions. If it was then it could potentially be overlapping
        // with polys of that region and the heights sampled here could be wrong.

        if (region != RecastConstants.RC_MULTIPLE_REGS) {
          // Copy the height from the same region, and mark region borders
          // as seed points to fill the rest.
          for (var hy = 0; hy < hp.height; hy++) {
            var y = hp.ymin + hy + bs;

            for (var hx = 0; hx < hp.width; hx++) {
              var x = hp.xmin + hx + bs;
              var c = chf.cells[x + y * chf.width];

              for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
                // if(hy == 1 && hx == 1 && i == 1344)
                // 	console.log("Bad")
                var s = chf.spans[i];

                if (s.reg == region) {
                  // Store height
                  hp.data[hx + hy * hp.width] = s.y;
                  empty = false; // If any of the neighbours is not _in same region,
                  // add the current location as flood fill start

                  var border = false;

                  for (var dir = 0; dir < 4; ++dir) {
                    if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
                      var ax = x + RecastCommon.GetDirOffsetX(dir);
                      var ay = y + RecastCommon.GetDirOffsetY(dir);
                      var ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(s, dir);
                      var as = chf.spans[ai];

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
        } // if the polygon does not contain any points from the current region (rare, but happens)
        // or if it could potentially be overlapping polygons of the same region,
        // then use the center as the seed point.		


        if (empty) RecastMeshDetail.seedArrayWithPolyCenter(ctx, chf, meshpolys, poly, npoly, verts, bs, hp, queue);
        var head = 0; // We assume the seed is centered _in the polygon, so a BFS to collect
        // height data will ensure we do not move onto overlapping polygons and
        // sample wrong heights.

        while (head * 3 < queue.length) {
          var cx = queue[head * 3 + 0];
          var cy = queue[head * 3 + 1];
          var ci = queue[head * 3 + 2];
          head++;

          if (head >= RecastMeshDetail.RETRACT_SIZE) {
            head = 0;
            queue = queue.slice(RecastMeshDetail.RETRACT_SIZE * 3, queue.length);
          }

          var cs = chf.spans[ci];

          for (var _dir = 0; _dir < 4; ++_dir) {
            if (RecastCommon.GetCon(cs, _dir) == RecastConstants.RC_NOT_CONNECTED) continue;

            var _ax = cx + RecastCommon.GetDirOffsetX(_dir);

            var _ay = cy + RecastCommon.GetDirOffsetY(_dir);

            var _hx = _ax - hp.xmin - bs;

            var _hy = _ay - hp.ymin - bs;

            if (_hx < 0 || _hx >= hp.width || _hy < 0 || _hy >= hp.height) continue;
            if (hp.data[_hx + _hy * hp.width] != RecastMeshDetail.RC_UNSET_HEIGHT) continue;

            var _ai = chf.cells[_ax + _ay * chf.width].index + RecastCommon.GetCon(cs, _dir);

            var _as = chf.spans[_ai];
            hp.data[_hx + _hy * hp.width] = _as.y;
            RecastMeshDetail.push3(queue, _ax, _ay, _ai);
          }
        }
      }
    }, {
      key: "getEdgeFlags",
      value: function getEdgeFlags(verts, va, vb, vpoly, npoly) {
        // Return true if edge (va,vb) is part of the polygon.
        var thrSqr = 0.001 * 0.001;

        for (var i = 0, j = npoly - 1; i < npoly; j = i++) {
          if (RecastMeshDetail.distancePtSeg2d(verts, va, vpoly, j * 3, i * 3) < thrSqr && RecastMeshDetail.distancePtSeg2d(verts, vb, vpoly, j * 3, i * 3) < thrSqr) return 1;
        }

        return 0;
      }
    }, {
      key: "getTriFlags",
      value: function getTriFlags(verts, va, vb, vc, vpoly, npoly) {
        var flags = 0;
        flags |= RecastMeshDetail.getEdgeFlags(verts, va, vb, vpoly, npoly) << 0;
        flags |= RecastMeshDetail.getEdgeFlags(verts, vb, vc, vpoly, npoly) << 2;
        flags |= RecastMeshDetail.getEdgeFlags(verts, vc, va, vpoly, npoly) << 4;
        return flags;
      } /// @par
      ///
      /// See the #rcConfig documentation for more information on the configuration parameters.
      ///
      /// @see rcAllocPolyMeshDetail, rcPolyMesh, rcCompactHeightfield, rcPolyMeshDetail, rcConfig

    }, {
      key: "buildPolyMeshDetail",
      value: function buildPolyMeshDetail(ctx, mesh, chf, sampleDist, sampleMaxError) {
        ctx.startTimer("BUILD_POLYMESHDETAIL");
        if (mesh.nverts == 0 || mesh.npolys == 0) return null;
        var dmesh = new PolyMeshDetail();
        var nvp = mesh.nvp;
        var cs = mesh.cs;
        var ch = mesh.ch;
        var orig = mesh.bmin;
        var borderSize = mesh.borderSize;
        var heightSearchRadius = Math.floor(Math.max(1, Math.ceil(mesh.maxEdgeError)));
        var tris = [];
        var verts = new Array(256 * 3).fill(0);
        var hp = new HeightPatch();
        var nPolyVerts = 0;
        var maxhw = 0,
            maxhh = 0;
        var bounds = new Array(mesh.npolys * 4).fill(0);
        var poly = new Array(nvp * 3).fill(0);
        poly.fill(0); // Find max size for a polygon area.

        for (var i = 0; i < mesh.npolys; ++i) {
          var p = i * nvp * 2;
          bounds[i * 4 + 0] = chf.width;
          bounds[i * 4 + 1] = 0;
          bounds[i * 4 + 2] = chf.height;
          bounds[i * 4 + 3] = 0;

          for (var j = 0; j < nvp; ++j) {
            if (mesh.polys[p + j] == RecastConstants.RC_MESH_NULL_IDX) {
              // console.log(i + " " + j)
              break;
            }

            var v = mesh.polys[p + j] * 3;
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
          if (bounds[i * 4 + 0] >= bounds[i * 4 + 1] || bounds[i * 4 + 2] >= bounds[i * 4 + 3]) continue;
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
        var vcap = nPolyVerts + Math.floor(nPolyVerts / 2);
        var tcap = vcap * 2;
        dmesh.nverts = 0;
        dmesh.verts = new Array(vcap * 3).fill(0);
        dmesh.ntris = 0;
        dmesh.tris = new Array(tcap * 4).fill(0);

        for (var _i10 = 0; _i10 < mesh.npolys; ++_i10) {
          var _p = _i10 * nvp * 2; // Store polygon vertices for processing.


          var npoly = 0;

          for (var _j3 = 0; _j3 < nvp; ++_j3) {
            if (mesh.polys[_p + _j3] == RecastConstants.RC_MESH_NULL_IDX) break;

            var _v = mesh.polys[_p + _j3] * 3;

            poly[_j3 * 3 + 0] = mesh.verts[_v + 0] * cs;
            poly[_j3 * 3 + 1] = mesh.verts[_v + 1] * ch;
            poly[_j3 * 3 + 2] = mesh.verts[_v + 2] * cs;
            npoly++;
          } // Get the height data from the area of the polygon.


          hp.xmin = bounds[_i10 * 4 + 0];
          hp.ymin = bounds[_i10 * 4 + 2];
          hp.width = bounds[_i10 * 4 + 1] - bounds[_i10 * 4 + 0];
          hp.height = bounds[_i10 * 4 + 3] - bounds[_i10 * 4 + 2];
          RecastMeshDetail.getHeightData(ctx, chf, mesh.polys, _p, npoly, mesh.verts, borderSize, hp, mesh.regs[_i10]); // Build detail mesh.

          var nverts = void 0;

          var _RecastMeshDetail$bui = RecastMeshDetail.buildPolyDetail(ctx, poly, npoly, sampleDist, sampleMaxError, heightSearchRadius, chf, hp, verts, tris);

          var _RecastMeshDetail$bui2 = _slicedToArray(_RecastMeshDetail$bui, 2);

          nverts = _RecastMeshDetail$bui2[0];
          tris = _RecastMeshDetail$bui2[1];

          // Move detail verts to world space.
          for (var _j4 = 0; _j4 < nverts; ++_j4) {
            verts[_j4 * 3 + 0] += orig[0];
            verts[_j4 * 3 + 1] += orig[1] + chf.ch; // Is this offset necessary?

            verts[_j4 * 3 + 2] += orig[2];
          } // Offset let too, will be used to flag checking.


          for (var _j5 = 0; _j5 < npoly; ++_j5) {
            poly[_j5 * 3 + 0] += orig[0];
            poly[_j5 * 3 + 1] += orig[1];
            poly[_j5 * 3 + 2] += orig[2];
          } // Store detail submesh.


          var ntris = tris.length / 4;
          dmesh.meshes[_i10 * 4 + 0] = dmesh.nverts;
          dmesh.meshes[_i10 * 4 + 1] = nverts;
          dmesh.meshes[_i10 * 4 + 2] = dmesh.ntris;
          dmesh.meshes[_i10 * 4 + 3] = ntris; // Store vertices, allocate more memory if necessary.

          if (dmesh.nverts + nverts > vcap) {
            while (dmesh.nverts + nverts > vcap) {
              vcap += 256;
            }

            var newv = new Array(vcap * 3);
            if (dmesh.nverts != 0) arraycopy$2(dmesh.verts, 0, newv, 0, 3 * dmesh.nverts);
            dmesh.verts = newv;
          }

          for (var _j6 = 0; _j6 < nverts; ++_j6) {
            dmesh.verts[dmesh.nverts * 3 + 0] = verts[_j6 * 3 + 0];
            dmesh.verts[dmesh.nverts * 3 + 1] = verts[_j6 * 3 + 1];
            dmesh.verts[dmesh.nverts * 3 + 2] = verts[_j6 * 3 + 2];
            dmesh.nverts++;
          } // Store triangles, allocate more memory if necessary.


          if (dmesh.ntris + ntris > tcap) {
            while (dmesh.ntris + ntris > tcap) {
              tcap += 256;
            }

            var newt = new Array(tcap * 4);
            if (dmesh.ntris != 0) arraycopy$2(dmesh.tris, 0, newt, 0, 4 * dmesh.ntris);
            dmesh.tris = newt;
          }

          for (var _j7 = 0; _j7 < ntris; ++_j7) {
            var t = _j7 * 4;
            dmesh.tris[dmesh.ntris * 4 + 0] = tris[t + 0];
            dmesh.tris[dmesh.ntris * 4 + 1] = tris[t + 1];
            dmesh.tris[dmesh.ntris * 4 + 2] = tris[t + 2];
            dmesh.tris[dmesh.ntris * 4 + 3] = RecastMeshDetail.getTriFlags(verts, tris[t + 0] * 3, tris[t + 1] * 3, tris[t + 2] * 3, poly, npoly);
            dmesh.ntris++;
          }
        }

        ctx.stopTimer("BUILD_POLYMESHDETAIL");
        return dmesh;
      }
    }]);

    return RecastMeshDetail;
  }();

  _defineProperty(RecastMeshDetail, "MAX_VERTS", 127);

  _defineProperty(RecastMeshDetail, "MAX_TRIS", 255);

  _defineProperty(RecastMeshDetail, "MAX_VERTS_PER_EDGE", 32);

  _defineProperty(RecastMeshDetail, "RC_UNSET_HEIGHT", 0xfff);

  _defineProperty(RecastMeshDetail, "EV_UNDEF", -1);

  _defineProperty(RecastMeshDetail, "EV_HULL", -2);

  _defineProperty(RecastMeshDetail, "vdot2", function (a, b) {
    return a[0] * b[0] + a[2] * b[2];
  });

  _defineProperty(RecastMeshDetail, "RETRACT_SIZE", 256);

  var RecastBuilderResult = /*#__PURE__*/function () {
    function RecastBuilderResult(pmesh, dmesh) {
      _classCallCheck(this, RecastBuilderResult);

      _defineProperty(this, "pmesh", void 0);

      _defineProperty(this, "dmesh", void 0);

      this.pmesh = pmesh;
      this.dmesh = dmesh;
    }

    _createClass(RecastBuilderResult, [{
      key: "getMesh",
      value: function getMesh() {
        return this.pmesh;
      }
    }, {
      key: "getMeshDetail",
      value: function getMeshDetail() {
        return this.dmesh;
      }
    }]);

    return RecastBuilderResult;
  }();

  var RecastBuilder = /*#__PURE__*/function () {
    // class RecastBuilderProgressListener {
    //     onProgress(completed, total);
    // }
    // RecastBuilder() {
    //     progressListener = null;
    // }
    function RecastBuilder() {
      var progressListener = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : null;

      _classCallCheck(this, RecastBuilder);

      _defineProperty(this, "progressListener", void 0);

      this.progressListener = progressListener;
    }

    _createClass(RecastBuilder, [{
      key: "buildTiles",
      value: function buildTiles(geom, cfg, threads) {
        var bmin = geom.getMeshBoundsMin();
        var bmax = geom.getMeshBoundsMax();
        var twh = Recast.calcTileCount(bmin, bmax, cfg.cs, cfg.tileSize);
        var tw = twh[0];
        var th = twh[1];
        result = null;

        if (threads == 1) {
          result = buildSingleThread(geom, cfg, bmin, bmax, tw, th);
        } else {
          result = buildMultiThread(geom, cfg, bmin, bmax, tw, th, threads);
        }

        return result;
      }
    }, {
      key: "buildSingleThread",
      value: function buildSingleThread(geom, cfg, bmin, bmax, tw, th) {
        result = new RecastBuilderResult[tw][th]();
        counter = new AtomicInteger();

        for (var x = 0; x < tw; ++x) {
          for (var y = 0; y < th; ++y) {
            result[x][y] = buildTile(geom, cfg, bmin, bmax, x, y, counter, tw * th);
          }
        }

        return result;
      }
    }, {
      key: "buildMultiThread",
      value: function buildMultiThread(geom, cfg, bmin, bmax, tw, th, threads) {
        ec = Executors.newFixedThreadPool(threads);
        result = new RecastBuilderResult[tw][th]();
        counter = new AtomicInteger();

        for (var x = 0; x < tw; ++x) {
          var _loop = function _loop(y) {
            var tx = x;
            var ty = y;
            ec.submit(function () {
              result[tx][ty] = buildTile(geom, cfg, bmin, bmax, tx, ty, counter, tw * th);
            });
          };

          for (var y = 0; y < th; ++y) {
            _loop(y);
          }
        }

        ec.shutdown();

        try {
          ec.awaitTermination(1000, TimeUnit.HOURS);
        } catch (e) {}

        return result;
      }
    }, {
      key: "buildTile",
      value: function buildTile(geom, cfg, bmin, bmax, tx, ty, counter, total) {
        result = build(geom, new RecastBuilderConfig(cfg, bmin, bmax, tx, ty, true));

        if (progressListener != null) {
          progressListener.onProgress(counter.incrementAndGet(), total);
        }

        return result;
      }
    }, {
      key: "build",
      value: function build(geom, builderCfg) {
        var cfg = builderCfg.cfg;
        var ctx = new Context();
        var chf = this.buildCompactHeightfield(geom, builderCfg, ctx); // Partition the heightfield so that we can use simple algorithm later
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
          RecastRegion.buildDistanceField(ctx, chf); // Partition the walkable surface into simple regions without holes.

          RecastRegion.buildRegions(ctx, chf, builderCfg.borderSize, cfg.minRegionArea, cfg.mergeRegionArea);
        } else if (cfg.partitionType == PartitionType.MONOTONE) {
          // Partition the walkable surface into simple regions without holes.
          // Monotone partitioning does not need distancefield.
          RecastRegion.buildRegionsMonotone(ctx, chf, builderCfg.borderSize, cfg.minRegionArea, cfg.mergeRegionArea);
        } else {
          // Partition the walkable surface into simple regions without holes.
          RecastRegion.buildLayerRegions(ctx, chf, builderCfg.borderSize, cfg.minRegionArea);
        } //
        // Step 5. Trace and simplify region contours.
        //
        // Create contours.


        var cset = RecastContour.buildContours(ctx, chf, cfg.maxSimplificationError, cfg.maxEdgeLen, RecastConstants.RC_CONTOUR_TESS_WALL_EDGES); //
        // Step 6. Build polygons mesh from contours.
        //

        var pmesh = RecastMesh.buildPolyMesh(ctx, cset, cfg.maxVertsPerPoly); //
        // Step 7. Create detail mesh which allows to access approximate height
        // on each polygon.
        //

        var dmesh = builderCfg.buildMeshDetail ? RecastMeshDetail.buildPolyMeshDetail(ctx, pmesh, chf, cfg.detailSampleDist, cfg.detailSampleMaxError) : null;
        return new RecastBuilderResult(pmesh, dmesh);
      }
    }, {
      key: "buildCompactHeightfield",
      value: function buildCompactHeightfield(geomProvider, builderCfg, ctx) {
        var cfg = builderCfg.cfg; //
        // Step 2. Rasterize input polygon soup.
        //
        // Allocate voxel heightfield where we rasterize our input data to.

        var solid = new Heightfield(builderCfg.width, builderCfg.height, builderCfg.bmin, builderCfg.bmax, cfg.cs, cfg.ch); // Allocate array that can hold triangle area types.
        // If you have multiple meshes you need to process, allocate
        // and array which can hold the max number of triangles you need to
        // process.
        // Find triangles which are walkable based on their slope and rasterize
        // them.
        // If your input data is multiple meshes, you can transform them here,
        // calculate
        // the are type for each of the meshes and rasterize them.

        var meshes = geomProvider.meshes();

        var _iterator = _createForOfIteratorHelper(meshes),
            _step;

        try {
          for (_iterator.s(); !(_step = _iterator.n()).done;) {
            var geom = _step.value;
            var verts = geom.getVerts();
            var tiled = cfg.tileSize > 0;
            var totaltris = 0;

            if (tiled) {
              var tbmin = new Array(2);
              var tbmax = new Array(2);
              tbmin[0] = builderCfg.bmin[0];
              tbmin[1] = builderCfg.bmin[2];
              tbmax[0] = builderCfg.bmax[0];
              tbmax[1] = builderCfg.bmax[2];
              var nodes = geom.getChunksOverlappingRect(tbmin, tbmax);

              var _iterator3 = _createForOfIteratorHelper(nodes),
                  _step3;

              try {
                for (_iterator3.s(); !(_step3 = _iterator3.n()).done;) {
                  var node = _step3.value;
                  var tris = node.tris;
                  var ntris = tris.length / 3;
                  totaltris += ntris;
                  var m_triareas = Recast.markWalkableTriangles(ctx, cfg.walkableSlopeAngle, verts, tris, ntris, cfg.walkableAreaMod);
                  RecastRasterization.rasterizeTriangles(ctx, verts, tris, m_triareas, ntris, solid, cfg.walkableClimb);
                }
              } catch (err) {
                _iterator3.e(err);
              } finally {
                _iterator3.f();
              }
            } else {
              var _tris = geom.getTris();

              var _ntris = _tris.length / 3;

              var _m_triareas = Recast.markWalkableTriangles(ctx, cfg.walkableSlopeAngle, verts, _tris, _ntris, cfg.walkableAreaMod);

              totaltris = _ntris;
              RecastRasterization.rasterizeTrianglesA(ctx, verts, _tris, _m_triareas, _ntris, solid, cfg.walkableClimb);
            }
          } // console.log(solid.spans[0])
          //
          // Step 3. Filter walkables surfaces.
          //
          // Once all geometry is rasterized, we do initial pass of filtering to
          // remove unwanted overhangs caused by the conservative rasterization
          // as well as filter spans where the character cannot possibly stand.

        } catch (err) {
          _iterator.e(err);
        } finally {
          _iterator.f();
        }

        RecastFilter.filterLowHangingWalkableObstacles(ctx, cfg.walkableClimb, solid);
        RecastFilter.filterLedgeSpans(ctx, cfg.walkableHeight, cfg.walkableClimb, solid);
        RecastFilter.filterWalkableLowHeightSpans(ctx, cfg.walkableHeight, solid); //
        // Step 4. Partition walkable surface to simple regions.
        //
        // Compact the heightfield so that it is faster to handle from now on.
        // This will result more cache coherent data as well as the neighbours
        // between walkable cells will be calculated.

        var chf = Recast.buildCompactHeightfield(ctx, cfg.walkableHeight, cfg.walkableClimb, solid); // console.log(chf.spans[450240])
        // Erode the walkable area by agent radius.

        RecastArea.erodeWalkableArea(ctx, cfg.walkableRadius, chf); // (Optional) Mark areas.

        var _iterator2 = _createForOfIteratorHelper(geomProvider.getConvexVolumes()),
            _step2;

        try {
          for (_iterator2.s(); !(_step2 = _iterator2.n()).done;) {
            var vol = _step2.value;
            RecastArea.markConvexPolyArea(ctx, vol.verts, vol.hmin, vol.hmax, vol.areaMod, chf);
          }
        } catch (err) {
          _iterator2.e(err);
        } finally {
          _iterator2.f();
        }

        return chf;
      }
    }, {
      key: "buildLayers",
      value: function buildLayers(geom, cfg) {
        var ctx = new Context();
        var chf = buildCompactHeightfield(geom, cfg, ctx);
        return RecastLayers.buildHeightfieldLayers(ctx, chf, cfg.borderSize, cfg.cfg.walkableHeight);
      }
    }]);

    return RecastBuilder;
  }();

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
  var NavMeshDataCreateParams = function NavMeshDataCreateParams() {
    _classCallCheck(this, NavMeshDataCreateParams);

    _defineProperty(this, "verts", []);

    _defineProperty(this, "vertCount", 0);

    _defineProperty(this, "polys", []);

    _defineProperty(this, "polyFlags", []);

    _defineProperty(this, "polyAreas", []);

    _defineProperty(this, "polyCount", 0);

    _defineProperty(this, "nvp", 0);

    _defineProperty(this, "detailMeshes", []);

    _defineProperty(this, "detailVerts", []);

    _defineProperty(this, "detailVertsCount", 0);

    _defineProperty(this, "detailTris", []);

    _defineProperty(this, "detailTriCount", 0);

    _defineProperty(this, "offMeshConVerts", []);

    _defineProperty(this, "offMeshConRad", []);

    _defineProperty(this, "offMeshConFlags", []);

    _defineProperty(this, "offMeshConAreas", []);

    _defineProperty(this, "offMeshConDir", []);

    _defineProperty(this, "offMeshConUserID", []);

    _defineProperty(this, "offMeshConCount", 0);

    _defineProperty(this, "userId", 0);

    _defineProperty(this, "tileX", 0);

    _defineProperty(this, "tileY", 0);

    _defineProperty(this, "tileLayer", 0);

    _defineProperty(this, "bmin", []);

    _defineProperty(this, "bmax", []);

    _defineProperty(this, "walkableHeight", 0);

    _defineProperty(this, "walkableRadius", 0);

    _defineProperty(this, "walkableClimb", 0);

    _defineProperty(this, "cs", 0);

    _defineProperty(this, "ch", 0);

    _defineProperty(this, "buildBvTree", false);
  } /// @}
  ;

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
  var MeshHeader = function MeshHeader() {
    _classCallCheck(this, MeshHeader);

    _defineProperty(this, "magic", void 0);

    _defineProperty(this, "version", void 0);

    _defineProperty(this, "x", void 0);

    _defineProperty(this, "y", void 0);

    _defineProperty(this, "layer", void 0);

    _defineProperty(this, "userId", void 0);

    _defineProperty(this, "polyCount", void 0);

    _defineProperty(this, "vertCount", void 0);

    _defineProperty(this, "maxLinkCount", void 0);

    _defineProperty(this, "detailMeshCount", void 0);

    _defineProperty(this, "detailVertCount", void 0);

    _defineProperty(this, "detailTriCount", void 0);

    _defineProperty(this, "bvNodeCount", void 0);

    _defineProperty(this, "offMeshConCount", void 0);

    _defineProperty(this, "offMeshBase", void 0);

    _defineProperty(this, "walkableHeight", void 0);

    _defineProperty(this, "walkableRadius", void 0);

    _defineProperty(this, "walkableClimb", void 0);

    _defineProperty(this, "bmin", new Array(3));

    _defineProperty(this, "bmax", new Array(3));

    _defineProperty(this, "bvQuantFactor", void 0);
  };

  _defineProperty(MeshHeader, "DT_NAVMESH_MAGIC", 'D' << 24 | 'N' << 16 | 'A' << 8 | 'V');

  _defineProperty(MeshHeader, "DT_NAVMESH_VERSION", 7);

  _defineProperty(MeshHeader, "DT_NAVMESH_VERSION_RECAST4J", 0x8807);

  _defineProperty(MeshHeader, "DT_NAVMESH_STATE_MAGIC", 'D' << 24 | 'N' << 16 | 'M' << 8 | 'S');

  _defineProperty(MeshHeader, "DT_NAVMESH_STATE_VERSION", 1);

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
  var PolyDetail = function PolyDetail() {
    _classCallCheck(this, PolyDetail);

    _defineProperty(this, "vertBase", 0);

    _defineProperty(this, "triBase", 0);

    _defineProperty(this, "vertCount", 0);

    _defineProperty(this, "triCount", 0);
  };

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
  var BVNode = function BVNode() {
    _classCallCheck(this, BVNode);

    _defineProperty(this, "bmin", new Array(3));

    _defineProperty(this, "bmax", new Array(3));

    _defineProperty(this, "i", 0);
  };

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
  var MeshData = function MeshData() {
    _classCallCheck(this, MeshData);

    _defineProperty(this, "header", void 0);

    _defineProperty(this, "verts", []);

    _defineProperty(this, "polys", []);

    _defineProperty(this, "detailMeshes", []);

    _defineProperty(this, "detailVerts", []);

    _defineProperty(this, "detailTris", []);

    _defineProperty(this, "bvTree", []);

    _defineProperty(this, "offMeshCons", []);
  };

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
  var OffMeshConnection = function OffMeshConnection() {
    _classCallCheck(this, OffMeshConnection);

    _defineProperty(this, "pos", new Array(6));

    _defineProperty(this, "rad", 0);

    _defineProperty(this, "poly", 0);

    _defineProperty(this, "flags", 0);

    _defineProperty(this, "side", 0);

    _defineProperty(this, "userId", 0);
  };

  function arraycopy$3(one, oneStart, two, twoStart, len) {
    for (var i = 0; i < len; i++) {
      two[twoStart + i] = one[oneStart + i];
    }
  }

  var BVItem = function BVItem() {
    _classCallCheck(this, BVItem);

    _defineProperty(this, "bmin", new Array(3).fill(0));

    _defineProperty(this, "bmax", new Array(3).fill(0));

    _defineProperty(this, "i", 0);
  };

  function compareX(a, b) {
    if (a.bmin[0] < b.bmin[0]) return -1;
    if (a.bmin[0] > b.bmin[0]) return 1;
    return 0;
  }

  function compareY(a, b) {
    if (a.bmin[1] < b.bmin[1]) return -1;
    if (a.bmin[1] > b.bmin[1]) return 1;
    return 0;
  }

  function compareZ(a, b) {
    if (a.bmin[2] < b.bmin[2]) return -1;
    if (a.bmin[2] > b.bmin[2]) return 1;
    return 0;
  }

  var NavMeshBuilder = /*#__PURE__*/function () {
    function NavMeshBuilder() {
      _classCallCheck(this, NavMeshBuilder);
    }

    _createClass(NavMeshBuilder, null, [{
      key: "calcExtends",
      value: function calcExtends(items, nitems, imin, imax) {
        var bmin = new Array(3);
        var bmax = new Array(3);
        bmin[0] = items[imin].bmin[0];
        bmin[1] = items[imin].bmin[1];
        bmin[2] = items[imin].bmin[2];
        bmax[0] = items[imin].bmax[0];
        bmax[1] = items[imin].bmax[1];
        bmax[2] = items[imin].bmax[2];

        for (var i = imin + 1; i < imax; ++i) {
          var it = items[i];
          if (it.bmin[0] < bmin[0]) bmin[0] = it.bmin[0];
          if (it.bmin[1] < bmin[1]) bmin[1] = it.bmin[1];
          if (it.bmin[2] < bmin[2]) bmin[2] = it.bmin[2];
          if (it.bmax[0] > bmax[0]) bmax[0] = it.bmax[0];
          if (it.bmax[1] > bmax[1]) bmax[1] = it.bmax[1];
          if (it.bmax[2] > bmax[2]) bmax[2] = it.bmax[2];
        }

        return [bmin, bmax];
      }
    }, {
      key: "longestAxis",
      value: function longestAxis(x, y, z) {
        var axis = 0;
        var maxVal = x;

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
    }, {
      key: "subdivide",
      value: function subdivide(items, nitems, imin, imax, curNode, nodes) {
        var inum = imax - imin;
        var icur = curNode;
        var node = new BVNode();
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
          var minmax = NavMeshBuilder.calcExtends(items, nitems, imin, imax);
          node.bmin = minmax[0];
          node.bmax = minmax[1];
          var axis = NavMeshBuilder.longestAxis(node.bmax[0] - node.bmin[0], node.bmax[1] - node.bmin[1], node.bmax[2] - node.bmin[2]);

          if (axis == 0) {
            // Sort aPoly x-axis
            //Arrays.sort(items, imin, imin + inum, new CompareItemX());
            //TODO: is this right?
            var shallowCopy = items.slice(imin, imin + inum);
            shallowCopy = shallowCopy.sort(compareX);

            for (var i = imin; i < imin + inum; i++) {
              items[i] = shallowCopy[i - imin];
              if (!items[i]) console.log("uh-oh");
            }
          } else if (axis == 1) {
            // Sort aPoly y-axis
            //Arrays.sort(items, imin, imin + inum, new CompareItemY());
            var _shallowCopy = items.slice(imin, imin + inum);

            _shallowCopy = _shallowCopy.sort(compareY);

            for (var _i = imin; _i < imin + inum; _i++) {
              items[_i] = _shallowCopy[_i - imin];
              if (!items[_i]) console.log("uh-oh");
            }
          } else {
            // Sort aPoly z-axis
            //Arrays.sort(items, imin, imin + inum, new CompareItemZ());
            var _shallowCopy2 = items.slice(imin, imin + inum);

            _shallowCopy2 = _shallowCopy2.sort(compareZ);

            for (var _i2 = imin; _i2 < imin + inum; _i2++) {
              items[_i2] = _shallowCopy2[_i2 - imin];
              if (!items[_i2]) console.log("uh-oh");
            }
          }

          var isplit = Math.floor(imin + inum / 2); // Left

          curNode = NavMeshBuilder.subdivide(items, nitems, imin, isplit, curNode, nodes); // Right

          curNode = NavMeshBuilder.subdivide(items, nitems, isplit, imax, curNode, nodes);
          var iescape = curNode - icur; // Negative index means escape.

          node.i = -iescape;
        }

        return curNode;
      }
    }, {
      key: "createBVTree",
      value: function createBVTree(params, nodes) {
        // Build tree
        var quantFactor = 1 / params.cs;
        var items = new Array(params.polyCount).fill(null);

        for (var i = 0; i < params.polyCount; i++) {
          var it = new BVItem();
          items[i] = it;
          it.i = i; // Calc polygon bounds. Use detail meshes if available.

          if (params.detailMeshes != null) {
            var vb = params.detailMeshes[i * 4 + 0];
            var ndv = params.detailMeshes[i * 4 + 1];
            var bmin = new Array(3);
            var bmax = new Array(3);
            var dv = vb * 3;
            DetourCommon.vCopy(bmin, params.detailVerts, dv);
            DetourCommon.vCopy(bmax, params.detailVerts, dv);

            for (var j = 1; j < ndv; j++) {
              DetourCommon.vMin(bmin, params.detailVerts, dv + j * 3);
              DetourCommon.vMax(bmax, params.detailVerts, dv + j * 3);
            } // BV-tree uses cs for all dimensions


            it.bmin[0] = DetourCommon.clamp(Math.round((bmin[0] - params.bmin[0]) * quantFactor), 0, 0xffff);
            it.bmin[1] = DetourCommon.clamp(Math.round((bmin[1] - params.bmin[1]) * quantFactor), 0, 0xffff);
            it.bmin[2] = DetourCommon.clamp(Math.round((bmin[2] - params.bmin[2]) * quantFactor), 0, 0xffff);
            it.bmax[0] = DetourCommon.clamp(Math.round((bmax[0] - params.bmin[0]) * quantFactor), 0, 0xffff);
            it.bmax[1] = DetourCommon.clamp(Math.round((bmax[1] - params.bmin[1]) * quantFactor), 0, 0xffff);
            it.bmax[2] = DetourCommon.clamp(Math.round((bmax[2] - params.bmin[2]) * quantFactor), 0, 0xffff);
          } else {
            var p = i * params.nvp * 2;
            it.bmin[0] = it.bmax[0] = params.verts[params.polys[p] * 3 + 0];
            it.bmin[1] = it.bmax[1] = params.verts[params.polys[p] * 3 + 1];
            it.bmin[2] = it.bmax[2] = params.verts[params.polys[p] * 3 + 2];

            for (var _j = 1; _j < params.nvp; ++_j) {
              if (params.polys[p + _j] == NavMeshBuilder.MESH_NULL_IDX) break;
              var x = params.verts[params.polys[p + _j] * 3 + 0];
              var y = params.verts[params.polys[p + _j] * 3 + 1];
              var z = params.verts[params.polys[p + _j] * 3 + 2];
              if (x < it.bmin[0]) it.bmin[0] = x;
              if (y < it.bmin[1]) it.bmin[1] = y;
              if (z < it.bmin[2]) it.bmin[2] = z;
              if (x > it.bmax[0]) it.bmax[0] = x;
              if (y > it.bmax[1]) it.bmax[1] = y;
              if (z > it.bmax[2]) it.bmax[2] = z;
            } // Remap y


            it.bmin[1] = Math.floor(it.bmin[1] * params.ch / params.cs);
            it.bmax[1] = Math.ceil(it.bmax[1] * params.ch / params.cs);
          }
        }

        return NavMeshBuilder.subdivide(items, params.polyCount, 0, params.polyCount, 0, nodes);
      }
    }, {
      key: "classifyOffMeshPoint",
      value: function classifyOffMeshPoint(pt, bmin, bmax) {
        var outcode = 0;
        outcode |= pt[0] >= bmax[0] ? NavMeshBuilder.XP : 0;
        outcode |= pt[2] >= bmax[2] ? ZP : 0;
        outcode |= pt[0] < bmin[0] ? XM : 0;
        outcode |= pt[2] < bmin[2] ? ZM : 0;

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

    }, {
      key: "createNavMeshData",
      value: function createNavMeshData(params) {
        if (params.vertCount >= 0xffff) return null;
        if (params.vertCount == 0 || params.verts == null) return null;
        if (params.polyCount == 0 || params.polys == null) return null;
        var nvp = params.nvp; // Classify off-mesh connection points. We store only the connections
        // whose start poPoly is inside the tile.

        var offMeshConClass = null;
        var storedOffMeshConCount = 0;
        var offMeshConLinkCount = 0;

        if (params.offMeshConCount > 0) {
          offMeshConClass = new Array(params.offMeshConCount * 2); // Find tight heigh bounds, used for culling out off-mesh start
          // locations.

          var hmin = Number.MAX_VALUE;
          var hmax = -Number.MAX_VALUE;

          if (params.detailVerts != null && params.detailVertsCount != 0) {
            for (var i = 0; i < params.detailVertsCount; ++i) {
              var _h = params.detailVerts[i * 3 + 1];
              hmin = Math.min(hmin, _h);
              hmax = Math.max(hmax, _h);
            }
          } else {
            for (var _i3 = 0; _i3 < params.vertCount; ++_i3) {
              var iv = _i3 * 3;
              h = params.bmin[1] + params.verts[iv + 1] * params.ch;
              hmin = Math.min(hmin, h);
              hmax = Math.max(hmax, h);
            }
          }

          hmin -= params.walkableClimb;
          hmax += params.walkableClimb;
          var bmin = new Array(3);
          var bmax = new Array(3);
          DetourCommon.vCopy(bmin, params.bmin);
          DetourCommon.vCopy(bmax, params.bmax);
          bmin[1] = hmin;
          bmax[1] = hmax;

          for (var _i4 = 0; _i4 < params.offMeshConCount; ++_i4) {
            var p0 = new VectorPtr$1(params.offMeshConVerts, (_i4 * 2 + 0) * 3);
            var p1 = new VectorPtr$1(params.offMeshConVerts, (_i4 * 2 + 1) * 3);
            offMeshConClass[_i4 * 2 + 0] = NavMeshBuilder.classifyOffMeshPoint(p0, bmin, bmax);
            offMeshConClass[_i4 * 2 + 1] = NavMeshBuilder.classifyOffMeshPoint(p1, bmin, bmax); // Zero out off-mesh start positions which are not even
            // potentially touching the mesh.

            if (offMeshConClass[_i4 * 2 + 0] == 0xff) {
              if (p0[1] < bmin[1] || p0[1] > bmax[1]) offMeshConClass[_i4 * 2 + 0] = 0;
            } // Count how many links should be allocated for off-mesh
            // connections.


            if (offMeshConClass[_i4 * 2 + 0] == 0xff) offMeshConLinkCount++;
            if (offMeshConClass[_i4 * 2 + 1] == 0xff) offMeshConLinkCount++;
            if (offMeshConClass[_i4 * 2 + 0] == 0xff) storedOffMeshConCount++;
          }
        } // Off-mesh connectionss are stored as polygons, adjust values.


        var totPolyCount = params.polyCount + storedOffMeshConCount;
        var totVertCount = params.vertCount + storedOffMeshConCount * 2; // Find portal edges which are at tile borders.

        var edgeCount = 0;
        var portalCount = 0;

        for (var _i5 = 0; _i5 < params.polyCount; ++_i5) {
          var p = _i5 * 2 * nvp;

          for (var j = 0; j < nvp; ++j) {
            if (params.polys[p + j] == NavMeshBuilder.MESH_NULL_IDX) break;
            edgeCount++;

            if ((params.polys[p + nvp + j] & 0x8000) != 0) {
              var dir = params.polys[p + nvp + j] & 0xf;
              if (dir != 0xf) portalCount++;
            }
          }
        }

        var maxLinkCount = edgeCount + portalCount * 2 + offMeshConLinkCount * 2; // Find unique detail vertices.

        var uniqueDetailVertCount = 0;
        var detailTriCount = 0;

        if (params.detailMeshes != null) {
          // Has detail mesh, count unique detail vertex count and use input
          // detail tri count.
          detailTriCount = params.detailTriCount;

          for (var _i6 = 0; _i6 < params.polyCount; ++_i6) {
            var _p = _i6 * nvp * 2;

            var ndv = params.detailMeshes[_i6 * 4 + 1];
            var nv = 0;

            for (var _j2 = 0; _j2 < nvp; ++_j2) {
              if (params.polys[_p + _j2] == NavMeshBuilder.MESH_NULL_IDX) break;
              nv++;
            }

            ndv -= nv;
            uniqueDetailVertCount += ndv;
          }
        } else {
          // No input detail mesh, build detail mesh from nav polys.
          uniqueDetailVertCount = 0; // No extra detail verts.

          detailTriCount = 0;

          for (var _i7 = 0; _i7 < params.polyCount; ++_i7) {
            var _p2 = _i7 * nvp * 2;

            var _nv = 0;

            for (var _j3 = 0; _j3 < nvp; ++_j3) {
              if (params.polys[_p2 + _j3] == NavMeshBuilder.MESH_NULL_IDX) break;
              _nv++;
            }

            detailTriCount += _nv - 2;
          }
        }

        var bvTreeSize = params.buildBvTree ? params.polyCount * 2 : 0;
        var header = new MeshHeader();
        var navVerts = new Array(3 * totVertCount);
        var navPolys = new Array(totPolyCount);
        var navDMeshes = new Array(params.polyCount);
        var navDVerts = new Array(3 * uniqueDetailVertCount);
        var navDTris = new Array(4 * detailTriCount);
        var navBvtree = new Array(bvTreeSize);
        var offMeshCons = new Array(storedOffMeshConCount); // Store header

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
        var offMeshVertsBase = params.vertCount;
        var offMeshPolyBase = params.polyCount; // Store vertices
        // Mesh vertices

        for (var _i8 = 0; _i8 < params.vertCount; ++_i8) {
          var _iv = _i8 * 3;

          var v = _i8 * 3;
          navVerts[v] = params.bmin[0] + params.verts[_iv] * params.cs;
          navVerts[v + 1] = params.bmin[1] + params.verts[_iv + 1] * params.ch;
          navVerts[v + 2] = params.bmin[2] + params.verts[_iv + 2] * params.cs; // console.log(navVerts)
        } // Off-mesh link vertices.


        var n = 0;

        for (var _i9 = 0; _i9 < params.offMeshConCount; ++_i9) {
          // Only store connections which start from this tile.
          if (offMeshConClass[_i9 * 2 + 0] == 0xff) {
            var linkv = _i9 * 2 * 3;

            var _v = (offMeshVertsBase + n * 2) * 3;

            arraycopy$3(params.offMeshConVerts, linkv, navVerts, _v, 6);
            n++;
          }
        } // Store polygons
        // Mesh polys


        var src = 0;

        for (var _i10 = 0; _i10 < params.polyCount; ++_i10) {
          var _p3 = new Poly(_i10, nvp);

          navPolys[_i10] = _p3;
          _p3.vertCount = 0;
          _p3.flags = params.polyFlags[_i10];

          _p3.setArea(params.polyAreas[_i10]);

          _p3.setType(Poly.DT_POLYTYPE_GROUND);

          for (var _j4 = 0; _j4 < nvp; ++_j4) {
            if (params.polys[src + _j4] == NavMeshBuilder.MESH_NULL_IDX) break;
            _p3.verts[_j4] = params.polys[src + _j4];

            if ((params.polys[src + nvp + _j4] & 0x8000) != 0) {
              // Border or portal edge.
              var _dir = params.polys[src + nvp + _j4] & 0xf;

              if (_dir == 0xf) // Border
                _p3.neis[_j4] = 0;else if (_dir == 0) // Portal x-
                _p3.neis[_j4] = NavMesh.DT_EXT_LINK | 4;else if (_dir == 1) // Portal z+
                _p3.neis[_j4] = NavMesh.DT_EXT_LINK | 2;else if (_dir == 2) // Portal x+
                _p3.neis[_j4] = NavMesh.DT_EXT_LINK | 0;else if (_dir == 3) // Portal z-
                _p3.neis[_j4] = NavMesh.DT_EXT_LINK | 6;
            } else {
              // Normal connection
              _p3.neis[_j4] = params.polys[src + nvp + _j4] + 1;
            }

            _p3.vertCount++;
          }

          src += nvp * 2;
        } // Off-mesh connection vertices.


        n = 0;

        for (var _i11 = 0; _i11 < params.offMeshConCount; ++_i11) {
          // Only store connections which start from this tile.
          if (offMeshConClass[_i11 * 2 + 0] == 0xff) {
            var _p4 = new Poly(offMeshPolyBase + n, nvp);

            navPolys[offMeshPolyBase + n] = _p4;
            _p4.vertCount = 2;
            _p4.verts[0] = offMeshVertsBase + n * 2;
            _p4.verts[1] = offMeshVertsBase + n * 2 + 1;
            _p4.flags = params.offMeshConFlags[_i11];

            _p4.setArea(params.offMeshConAreas[_i11]);

            _p4.setType(Poly.DT_POLYTYPE_OFFMESH_CONNECTION);

            n++;
          }
        } // Store detail meshes and vertices.
        // The nav polygon vertices are stored as the first vertices on each
        // mesh.
        // We compress the mesh data by skipping them and using the navmesh
        // coordinates.


        if (params.detailMeshes != null) {
          var vbase = 0;

          for (var _i12 = 0; _i12 < params.polyCount; ++_i12) {
            var dtl = new PolyDetail();
            navDMeshes[_i12] = dtl;
            var vb = params.detailMeshes[_i12 * 4 + 0];
            var _ndv = params.detailMeshes[_i12 * 4 + 1];
            var _nv2 = navPolys[_i12].vertCount;
            dtl.vertBase = vbase;
            dtl.vertCount = _ndv - _nv2;
            dtl.triBase = params.detailMeshes[_i12 * 4 + 2];
            dtl.triCount = params.detailMeshes[_i12 * 4 + 3]; // Copy vertices except the first 'nv' verts which are equal to
            // nav let verts.

            if (_ndv - _nv2 != 0) {
              arraycopy$3(params.detailVerts, (vb + _nv2) * 3, navDVerts, vbase * 3, 3 * (_ndv - _nv2));
              vbase += _ndv - _nv2;
            }
          } // Store triangles.


          arraycopy$3(params.detailTris, 0, navDTris, 0, 4 * params.detailTriCount);
        } else {
          // Create dummy detail mesh by triangulating polys.
          var tbase = 0;

          for (var _i13 = 0; _i13 < params.polyCount; ++_i13) {
            var _dtl = new PolyDetail();

            navDMeshes[_i13] = _dtl;
            var _nv3 = navPolys[_i13].vertCount;
            _dtl.vertBase = 0;
            _dtl.vertCount = 0;
            _dtl.triBase = tbase;
            _dtl.triCount = _nv3 - 2; // Triangulate polygon (local indices).

            for (var _j5 = 2; _j5 < _nv3; ++_j5) {
              var t = tbase * 4;
              navDTris[t + 0] = 0;
              navDTris[t + 1] = _j5 - 1;
              navDTris[t + 2] = _j5; // Bit for each edge that belongs to let boundary.

              navDTris[t + 3] = 1 << 2;
              if (_j5 == 2) navDTris[t + 3] |= 1 << 0;
              if (_j5 == _nv3 - 1) navDTris[t + 3] |= 1 << 4;
              tbase++;
            }
          }
        } // Store and create BVtree.
        // TODO: take detail mesh into account! use byte per bbox extent?


        if (params.buildBvTree) {
          // Do not set header.bvNodeCount set to make it work look exactly the same as in original Detour  
          header.bvNodeCount = NavMeshBuilder.createBVTree(params, navBvtree);
        } // Store Off-Mesh connections.


        n = 0;

        for (var _i14 = 0; _i14 < params.offMeshConCount; ++_i14) {
          // Only store connections which start from this tile.
          if (offMeshConClass[_i14 * 2 + 0] == 0xff) {
            var con = new OffMeshConnection();
            offMeshCons[n] = con;
            con.poly = offMeshPolyBase + n; // Copy connection end-points.

            var endPts = _i14 * 2 * 3;
            arraycopy$3(params.offMeshConVerts, endPts, con.pos, 0, 6);
            con.rad = params.offMeshConRad[_i14];
            con.flags = params.offMeshConDir[_i14] != 0 ? NavMesh.DT_OFFMESH_CON_BIDIR : 0;
            con.side = offMeshConClass[_i14 * 2 + 1];
            if (params.offMeshConUserID != null) con.userId = params.offMeshConUserID[_i14];
            n++;
          }
        }

        var nmd = new MeshData();
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
    }]);

    return NavMeshBuilder;
  }();

  _defineProperty(NavMeshBuilder, "MESH_NULL_IDX", 0xffff);

  _defineProperty(NavMeshBuilder, "XP", 1 << 0);

  _defineProperty(NavMeshBuilder, "ZP", 1 << 1);

  _defineProperty(NavMeshBuilder, "XM", 1 << 2);

  _defineProperty(NavMeshBuilder, "ZM", 1 << 3);

  var RecastTestMeshBuilder = /*#__PURE__*/function () {
    function RecastTestMeshBuilder(m_geom, m_partitionType, m_cellSize, m_cellHeight, m_agentHeight, m_agentRadius, m_agentMaxClimb, m_agentMaxSlope, m_regionMinSize, m_regionMergeSize, m_edgeMaxLen, m_edgeMaxError, m_vertsPerPoly, m_detailSampleDist, m_detailSampleMaxError) {
      _classCallCheck(this, RecastTestMeshBuilder);

      _defineProperty(this, "meshData", void 0);

      var cfg = new RecastConfig(m_partitionType, m_cellSize, m_cellHeight, m_agentHeight, m_agentRadius, m_agentMaxClimb, m_agentMaxSlope, m_regionMinSize, m_regionMergeSize, m_edgeMaxLen, m_edgeMaxError, m_vertsPerPoly, m_detailSampleDist, m_detailSampleMaxError, 0, SampleAreaModifications.SAMPLE_AREAMOD_GROUND);
      var bcfg = new RecastBuilderConfig$1(cfg, m_geom.getMeshBoundsMin(), m_geom.getMeshBoundsMax());
      var rcBuilder = new RecastBuilder();
      var rcResult = rcBuilder.build(m_geom, bcfg);
      var m_pmesh = rcResult.getMesh();

      for (var i = 0; i < m_pmesh.npolys; ++i) {
        m_pmesh.flags[i] = 1;
      }

      var m_dmesh = rcResult.getMeshDetail();
      var params = new NavMeshDataCreateParams();
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

    _createClass(RecastTestMeshBuilder, [{
      key: "getMeshData",
      value: function getMeshData() {
        return this.meshData;
      }
    }]);

    return RecastTestMeshBuilder;
  }();

  _defineProperty(RecastTestMeshBuilder, "m_cellSize", 0.3);

  _defineProperty(RecastTestMeshBuilder, "m_cellHeight", 0.2);

  _defineProperty(RecastTestMeshBuilder, "m_agentHeight", 2.0);

  _defineProperty(RecastTestMeshBuilder, "m_agentRadius", 0.6);

  _defineProperty(RecastTestMeshBuilder, "m_agentMaxClimb", 0.9);

  _defineProperty(RecastTestMeshBuilder, "m_agentMaxSlope", 45.0);

  _defineProperty(RecastTestMeshBuilder, "m_regionMinSize", 8);

  _defineProperty(RecastTestMeshBuilder, "m_regionMergeSize", 20);

  _defineProperty(RecastTestMeshBuilder, "m_edgeMaxLen", 12.0);

  _defineProperty(RecastTestMeshBuilder, "m_edgeMaxError", 1.3);

  _defineProperty(RecastTestMeshBuilder, "m_vertsPerPoly", 6);

  _defineProperty(RecastTestMeshBuilder, "m_detailSampleDist", 6.0);

  _defineProperty(RecastTestMeshBuilder, "m_detailSampleMaxError", 1.0);

  _defineProperty(RecastTestMeshBuilder, "fromFile", function (objFileContents) {
    return new RecastTestMeshBuilder(new ObjImporter().load(objFileContents), RecastConstants.WATERSHED, this.m_cellSize, this.m_cellHeight, this.m_agentHeight, this.m_agentRadius, this.m_agentMaxClimb, this.m_agentMaxSlope, this.m_regionMinSize, this.m_regionMergeSize, this.m_edgeMaxLen, this.m_edgeMaxError, this.m_vertsPerPoly, this.m_detailSampleDist, this.m_detailSampleMaxError);
  });

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
  var Node = /*#__PURE__*/function () {
    /** parent of the node is not adjacent. Found using raycast. */

    /** Position of the node. */

    /** Cost from previous node to current node. */

    /** Cost up to the node. */

    /** Index to parent node. */

    /** extra state information. A polyRef can have multiple nodes with different extra info. see DT_MAX_STATES_PER_NODE */

    /** Node flags. A combination of dtNodeFlags. */

    /** Polygon ref the node corresponds to. */
    function Node(index) {
      _classCallCheck(this, Node);

      _defineProperty(this, "index", 0);

      _defineProperty(this, "pos", new Array(3));

      _defineProperty(this, "cost", 0);

      _defineProperty(this, "total", 0);

      _defineProperty(this, "pidx", 0);

      _defineProperty(this, "state", 0);

      _defineProperty(this, "flags", 0);

      _defineProperty(this, "id", 0);

      this.index = index;
    }

    _createClass(Node, [{
      key: "toString",
      value: function toString() {
        return "Node [id=" + id + "]";
      }
    }]);

    return Node;
  }();

  _defineProperty(Node, "DT_NODE_OPEN", 0x01);

  _defineProperty(Node, "DT_NODE_CLOSED", 0x02);

  _defineProperty(Node, "DT_NODE_PARENT_DETACHED", 0x04);

  var NodePool = /*#__PURE__*/function () {
    function NodePool() {
      _classCallCheck(this, NodePool);

      _defineProperty(this, "m_map", []);

      _defineProperty(this, "m_nodes", []);
    }

    _createClass(NodePool, [{
      key: "clear",
      value: function clear() {
        this.m_nodes = [];
        this.m_map = [];
      }
    }, {
      key: "findNodes",
      value: function findNodes(id) {
        var nodes = this.m_map[id];

        if (nodes == null) {
          nodes = [];
        }

        return nodes;
      }
    }, {
      key: "findNode",
      value: function findNode(id) {
        var nodes = this.m_map[id];

        if (nodes != null && !nodes.length == 0) {
          return nodes[0];
        }

        return null;
      }
    }, {
      key: "findNode",
      value: function findNode(id, state) {
        var nodes = this.m_map[id];

        if (nodes != null) {
          var _iterator = _createForOfIteratorHelper(nodes),
              _step;

          try {
            for (_iterator.s(); !(_step = _iterator.n()).done;) {
              var node = _step.value;

              if (node.state == state) {
                return node;
              }
            }
          } catch (err) {
            _iterator.e(err);
          } finally {
            _iterator.f();
          }
        }

        return null;
      }
    }, {
      key: "getNode",
      value: function getNode(id) {
        var state = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
        var nodes = this.m_map[id];

        if (nodes != null) {
          var _iterator2 = _createForOfIteratorHelper(nodes),
              _step2;

          try {
            for (_iterator2.s(); !(_step2 = _iterator2.n()).done;) {
              var node = _step2.value;

              if (node.state == state) {
                return node;
              }
            }
          } catch (err) {
            _iterator2.e(err);
          } finally {
            _iterator2.f();
          }
        }

        return this.create(id, state);
      }
    }, {
      key: "create",
      value: function create(id, state) {
        var node = new Node(this.m_nodes.length + 1);
        node.id = id;
        node.state = state;
        this.m_nodes.push(node);
        var nodes = this.m_map[id];

        if (nodes == null) {
          nodes = [];
          this.m_map[id] = nodes;
        }

        nodes.push(node);
        return node;
      }
    }, {
      key: "getNodeIdx",
      value: function getNodeIdx(node) {
        return node != null ? node.index : 0;
      }
    }, {
      key: "getNodeAtIdx",
      value: function getNodeAtIdx(idx) {
        return idx != 0 ? this.m_nodes[idx - 1] : null;
      }
    }, {
      key: "getNodeCount",
      value: function getNodeCount() {
        return this.m_nodes.length;
      } // getNode(ref) {
      // return this.getNode(ref, 0);
      // }

      /*
      
      inline let getMaxNodes() const { return m_maxNodes; }
      inline dtNodeIndex getFirst(bucket) const { return m_first[bucket]; }
      inline dtNodeIndex getNext(i) const { return m_next[i]; }
      */

    }]);

    return NodePool;
  }();

  //https://truetocode.com/binary-treemax-heap-priority-queue-and-implementation-using-javascript/427/
  var PriorityQueue = /*#__PURE__*/function () {
    function PriorityQueue(comparator) {
      _classCallCheck(this, PriorityQueue);

      this.comparator = comparator;
      this.array = [];
    }

    _createClass(PriorityQueue, [{
      key: "isEmpty",
      value: function isEmpty() {
        return this.array.length == 0;
      }
    }, {
      key: "poll",
      value: function poll() {
        return this.array.splice(0, 1)[0];
      }
    }, {
      key: "insert",
      value: function insert(element) {
        this.push(element);
      }
    }, {
      key: "push",
      value: function push(element) {
        this.array.push(element);
        this.array.sort(this.comparator);
      }
    }, {
      key: "offer",
      value: function offer(element) {
        this.push(element);
      }
    }, {
      key: "remove",
      value: function remove(element) {
        var index = this.array.indexOf(element);
        if (index >= 0) this.array.splice(index, 1);
      }
    }]);

    return PriorityQueue;
  }();

  function compare(one, two) {
    if (one == two) return 0;
    if (one < two) return -1;else return 1;
  }

  var NodeQueue = /*#__PURE__*/function () {
    function NodeQueue() {
      _classCallCheck(this, NodeQueue);

      this.m_heap = new PriorityQueue(function (n1, n2) {
        return compare(n1.total, n2.total);
      });
    }

    _createClass(NodeQueue, [{
      key: "clear",
      value: function clear() {
        this.m_heap = new PriorityQueue(function (n1, n2) {
          return compare(n1.total, n2.total);
        });
      }
    }, {
      key: "top",
      value: function top() {
        return this.m_heap.peek();
      }
    }, {
      key: "pop",
      value: function pop() {
        return this.m_heap.poll();
      }
    }, {
      key: "push",
      value: function push(node) {
        this.m_heap.insert(node);
      }
    }, {
      key: "modify",
      value: function modify(node) {
        this.m_heap.remove(node);
        this.m_heap.offer(node);
      }
    }, {
      key: "isEmpty",
      value: function isEmpty() {
        return this.m_heap.isEmpty();
      }
    }]);

    return NodeQueue;
  }();

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
  var QueryData = function QueryData() {
    _classCallCheck(this, QueryData);

    _defineProperty(this, "status", null);

    _defineProperty(this, "lastBestNode", null);

    _defineProperty(this, "lastBestNodeCost", 0);

    _defineProperty(this, "startRef", 0);

    _defineProperty(this, 5, void 0);

    _defineProperty(this, "endRef", 0);

    _defineProperty(this, "startPos", new Array(3));

    _defineProperty(this, "endPos", new Array(3));

    _defineProperty(this, "filter", null);

    _defineProperty(this, "options", 0);

    _defineProperty(this, "raycastLimitSqr", 0);
  };

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
  var Status = /*#__PURE__*/function () {
    function Status(inValue) {
      _classCallCheck(this, Status);

      this.value = inValue;
    }

    _createClass(Status, [{
      key: "isFailed",
      value: function isFailed() {
        return this.value == Status.FAILURE;
      }
    }, {
      key: "isInProgress",
      value: function isInProgress() {
        return this.value == Status.IN_PROGRESS;
      }
    }, {
      key: "isSuccess",
      value: function isSuccess() {
        return this.value == Status.SUCCSESS || this.value == Status.PARTIAL_RESULT;
      }
    }, {
      key: "isPartial",
      value: function isPartial() {
        return this.value == Status.PARTIAL_RESULT;
      }
    }]);

    return Status;
  }();

  _defineProperty(Status, "FAILURE", 0);

  _defineProperty(Status, "SUCCSESS", 1);

  _defineProperty(Status, "IN_PROGRESS", 2);

  _defineProperty(Status, "PARTIAL_RESULT", 3);

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
  var UpdateSlicedPathResult = /*#__PURE__*/function () {
    function UpdateSlicedPathResult(status, iterations) {
      _classCallCheck(this, UpdateSlicedPathResult);

      _defineProperty(this, "status", null);

      _defineProperty(this, "iterations", 0);

      this.status = status;
      this.iterations = iterations;
    }

    _createClass(UpdateSlicedPathResult, [{
      key: "getStatus",
      value: function getStatus() {
        return this.status;
      }
    }, {
      key: "getIterations",
      value: function getIterations() {
        return this.iterations;
      }
    }]);

    return UpdateSlicedPathResult;
  }();

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
  var FindPathResult = /*#__PURE__*/function () {
    ///  @param[out]	path		An ordered list of polygon references representing the path. (Start to end.) 
    function FindPathResult(status, refs) {
      _classCallCheck(this, FindPathResult);

      _defineProperty(this, "status", null);

      _defineProperty(this, "refs", []);

      this.status = status;
      this.refs = refs;
    }

    _createClass(FindPathResult, [{
      key: "getStatus",
      value: function getStatus() {
        return this.status;
      }
    }, {
      key: "getRefs",
      value: function getRefs() {
        return this.refs;
      }
    }]);

    return FindPathResult;
  }();

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
  var FindLocalNeighbourhoodResult = /*#__PURE__*/function () {
    function FindLocalNeighbourhoodResult(refs, parentRefs) {
      _classCallCheck(this, FindLocalNeighbourhoodResult);

      _defineProperty(this, "refs", []);

      _defineProperty(this, "parentRefs", []);

      this.refs = refs;
      this.parentRefs = parentRefs;
    }

    _createClass(FindLocalNeighbourhoodResult, [{
      key: "getRefs",
      value: function getRefs() {
        return this.refs;
      }
    }, {
      key: "getParentRefs",
      value: function getParentRefs() {
        return this.parentRefs;
      }
    }]);

    return FindLocalNeighbourhoodResult;
  }();

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
  var GetPolyWallSegmentsResult = /*#__PURE__*/function () {
    function GetPolyWallSegmentsResult(segmentVerts, segmentRefs) {
      _classCallCheck(this, GetPolyWallSegmentsResult);

      _defineProperty(this, "segmentVerts", []);

      _defineProperty(this, "segmentRefs", []);

      this.segmentVerts = segmentVerts;
      this.segmentRefs = segmentRefs;
    }

    _createClass(GetPolyWallSegmentsResult, [{
      key: "getSegmentVerts",
      value: function getSegmentVerts() {
        return this.segmentVerts;
      }
    }, {
      key: "getSegmentRefs",
      value: function getSegmentRefs() {
        return this.segmentRefs;
      }
    }]);

    return GetPolyWallSegmentsResult;
  }();

  var StraightPathItem = /*#__PURE__*/function () {
    function StraightPathItem(pos, flags, ref) {
      _classCallCheck(this, StraightPathItem);

      _defineProperty(this, "pos", []);

      _defineProperty(this, "flags", 0);

      _defineProperty(this, "ref", 0);

      this.pos = DetourCommon.vCopy_return(pos);
      this.flags = flags;
      this.ref = ref;
    }

    _createClass(StraightPathItem, [{
      key: "getPos",
      value: function getPos() {
        return this.pos;
      }
    }, {
      key: "getFlags",
      value: function getFlags() {
        return this.flags;
      }
    }, {
      key: "getRef",
      value: function getRef() {
        return this.ref;
      }
    }]);

    return StraightPathItem;
  }();

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
  var MoveAlongSurfaceResult = /*#__PURE__*/function () {
    /** The result position of the mover. [(x, y, z)] */

    /** The reference ids of the polygons visited during the move. */
    function MoveAlongSurfaceResult(resultPos, visited) {
      _classCallCheck(this, MoveAlongSurfaceResult);

      _defineProperty(this, "resultPos", []);

      _defineProperty(this, "visited", []);

      this.resultPos = resultPos;
      this.visited = visited;
    }

    _createClass(MoveAlongSurfaceResult, [{
      key: "getResultPos",
      value: function getResultPos() {
        return this.resultPos;
      }
    }, {
      key: "getVisited",
      value: function getVisited() {
        return this.visited;
      }
    }]);

    return MoveAlongSurfaceResult;
  }();

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
  var RaycastHit = function RaycastHit() {
    _classCallCheck(this, RaycastHit);

    _defineProperty(this, "t", 0);

    _defineProperty(this, "hitNormal", new Array(3));

    _defineProperty(this, "path", []);

    _defineProperty(this, "pathCost", 0);

    _defineProperty(this, "hitEdgeIndex", 0);
  };

  var _temp$3;

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

  var PortalResult = function PortalResult(left, right, fromType, toType) {
    _classCallCheck(this, PortalResult);

    _defineProperty(this, "left", void 0);

    _defineProperty(this, "right", void 0);

    _defineProperty(this, "fromType", void 0);

    _defineProperty(this, "toType", void 0);

    this.left = left;
    this.right = right;
    this.fromType = fromType;
    this.toType = toType;
  };

  function arraycopy$4(one, oneStart, two, twoStart, len) {
    for (var i = 0; i < len; i++) {
      two[twoStart + i] = one[oneStart + i];
    }
  }

  var NavMeshQuery = /*#__PURE__*/function () {
    /**
     * Use raycasts during pathfind to "shortcut" (raycast still consider costs)
     * Options for NavMeshQuery::initSlicedFindPath and updateSlicedFindPath
     */

    /** Raycast should calculate movement cost aPoly the ray and fill RaycastHit::cost */
    /// Vertex flags returned by findStraightPath.

    /** The vertex is the start position in the path. */

    /** The vertex is the end position in the path. */

    /** The vertex is the start of an off-mesh connection. */
    /// Options for findStraightPath.
    ///< Add a vertex at every polygon edge crossing where area changes.
    ///< Add a vertex at every polygon edge crossing.
    // Search heuristic scale.
    /// < Sliced query state.
    function NavMeshQuery(nav) {
      _classCallCheck(this, NavMeshQuery);

      _defineProperty(this, "m_nav", void 0);

      _defineProperty(this, "m_nodePool", void 0);

      _defineProperty(this, "m_tinyNodePool", void 0);

      _defineProperty(this, "m_openList", void 0);

      _defineProperty(this, "m_query", void 0);

      this.m_nav = nav;
      this.m_nodePool = new NodePool();
      this.m_tinyNodePool = new NodePool();
      this.m_openList = new NodeQueue();
    }

    _createClass(NavMeshQuery, [{
      key: "findRandomPoint",

      /**
       * Returns random location on navmesh.
       * Polygons are chosen weighted by area. The search runs in linear related to number of polygon.
       * @param filter The polygon filter to apply to the query.
       * @param frand Function returning a random number [0..1).
       * @return Random location
       */
      value: function findRandomPoint(filter, frand) {
        // Randomly pick one tile. Assume that all tiles cover roughly the same area.
        tile = null;
        tsum = 0.0;

        for (var i = 0; i < this.m_nav.getMaxTiles(); i++) {
          var _t = this.m_nav.getTile(i);

          if (_t == null || _t.data == null || _t.data.header == null) continue; // Choose random tile using reservoi sampling.

          area = 1.0; // Could be tile area too.

          tsum += area;
          u = frand.frand();
          if (u * tsum <= area) tile = _t;
        }

        if (tile == null) return new FindRandomPointResult(Status.FAILURE, 0, null); // Randomly pick one polygon weighted by polygon area.

        var poly = null;
        var polyRef = 0;
        var base = this.m_nav.getPolyRefBase(tile);
        areaSum = 0.0;

        for (var _i = 0; _i < tile.data.header.polyCount; ++_i) {
          var p = tile.data.polys[_i]; // Do not return off-mesh connection polygons.

          if (p.getType() != Poly.DT_POLYTYPE_GROUND) continue; // Must pass filter

          var ref = this.or(base, _i);
          if (!filter.passFilter(ref, tile, p)) continue; // Calc area of the polygon.

          polyArea = 0.0;

          for (var j = 2; j < p.vertCount; ++j) {
            var va = p.verts[0] * 3;
            var vb = p.verts[j - 1] * 3;
            var vc = p.verts[j] * 3;
            polyArea += DetourCommon.triArea2D4(tile.data.verts, va, vb, vc);
          } // Choose random polygon weighted by area, using reservoi sampling.


          areaSum += polyArea;
          u = frand.frand();

          if (u * areaSum <= polyArea) {
            poly = p;
            polyRef = ref;
          }
        }

        if (poly == null) return new FindRandomPointResult(Status.FAILURE, 0, null); // Randomly pick poPoly on polygon.

        var verts = new Array(3 * this.m_nav.getMaxVertsPerPoly());
        var areas = new Array(this.m_nav.getMaxVertsPerPoly());
        arraycopy$4(tile.data.verts, poly.verts[0] * 3, verts, 0, 3);

        for (var _j = 1; _j < poly.vertCount; ++_j) {
          arraycopy$4(tile.data.verts, poly.verts[_j] * 3, verts, _j * 3, 3);
        }

        s = frand.frand();
        t = frand.frand();
        var pt = DetourCommon.randomPointInConvexPoly(verts, poly.vertCount, areas, s, t);
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

    }, {
      key: "findRandomPointAroundCircle",
      value: function findRandomPointAroundCircle(startRef, centerPos, maxRadius, filter, frand) {
        // Validate input
        if (startRef == 0 || !this.m_nav.isValidPolyRef(startRef)) throw new IllegalArgumentException("Invalid start ref");
        var tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(startRef);
        var startTile = tileAndPoly[0];
        var startPoly = tileAndPoly[1];
        if (!filter.passFilter(startRef, startTile, startPoly)) throw new IllegalArgumentException("Invalid start");
        m_nodePool = [];
        this.m_openList = [];
        var startNode = m_nodePool.getNode(startRef);
        DetourCommon.vCopy(startNode.pos, centerPos);
        startNode.pidx = 0;
        startNode.cost = 0;
        startNode.total = 0;
        startNode.id = startRef;
        startNode.flags = DT_NODE_OPEN;
        this.m_openList.push(startNode);
        radiusSqr = maxRadius * maxRadius;
        areaSum = 0.0;
        var randomTile = null;
        var randomPoly = null;
        var randomPolyRef = 0;

        while (!this.m_openList.length == 0) {
          var bestNode = this.m_openList.pop();
          bestNode.flags &= ~DT_NODE_OPEN;
          bestNode.flags |= DT_NODE_CLOSED; // Get let and tile.
          // The API input has been cheked already, skip checking internal data.

          var bestRef = bestNode.id;
          var bestTilePoly = this.m_nav.getTileAndPolyByRefUnsafe(bestRef);
          var bestTile = bestTilePoly[0];
          var bestPoly = bestTilePoly[1]; // Place random locations on on ground.

          if (bestPoly.getType() == Poly.DT_POLYTYPE_GROUND) {
            // Calc area of the polygon.
            polyArea = 0.0;

            for (var j = 2; j < bestPoly.vertCount; ++j) {
              var va = bestPoly.verts[0] * 3;
              var vb = bestPoly.verts[j - 1] * 3;
              var vc = bestPoly.verts[j] * 3;
              polyArea += DetourCommon.triArea2D4(bestTile.data.verts, va, vb, vc);
            } // Choose random polygon weighted by area, using reservoi sampling.


            areaSum += polyArea;
            u = frand.frand();

            if (u * areaSum <= polyArea) {
              randomTile = bestTile;
              randomPoly = bestPoly;
              randomPolyRef = bestRef;
            }
          } // Get parent let and tile.


          var parentRef = 0;
          if (bestNode.pidx != 0) parentRef = m_nodePool.getNodeAtIdx(bestNode.pidx).id;

          if (parentRef != 0) {
            var parentTilePoly = this.m_nav.getTileAndPolyByRefUnsafe(parentRef);
            var parentTile = parentTilePoly[0];
            var parentPoly = parentTilePoly[1];
          }

          for (var i = bestPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = bestTile.links[i].next) {
            var link = bestTile.links[i];
            var neighbourRef = link.ref; // Skip invalid neighbours and do not follow back to parent.

            if (neighbourRef == 0 || neighbourRef == parentRef) continue; // Expand to neighbour

            var neighbourTilePoly = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
            var neighbourTile = neighbourTilePoly[0];
            var neighbourPoly = neighbourTilePoly[1]; // Do not advance if the polygon is excluded by the filter.

            if (!filter.passFilter(neighbourRef, neighbourTile, neighbourPoly)) continue; // Find edge and calc distance to the edge.

            var portalpoints = this.getPortalPoints7(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile, 0, 0);
            var _va = portalpoints.left;
            var _vb = portalpoints.right; // If the circle is not touching the next polygon, skip it.

            var distseg = DetourCommon.distancePtSegSqr2D3(centerPos, _va, _vb);
            distSqr = distseg[0];
            if (distSqr > radiusSqr) continue;
            var neighbourNode = m_nodePool.getNode(neighbourRef);
            if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0) continue; // Cost

            if (neighbourNode.flags == 0) neighbourNode.pos = DetourCommon.vLerp3(_va, _vb, 0.5);
            total = bestNode.total + DetourCommon.vDist2(bestNode.pos, neighbourNode.pos); // The node is already in open list and the new result is worse, skip.

            if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0 && total >= neighbourNode.total) continue;
            neighbourNode.id = neighbourRef;
            neighbourNode.flags = neighbourNode.flags & ~Node.DT_NODE_CLOSED;
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

        if (randomPoly == null) return new FindRandomPointResult(Status.FAILURE, 0, null); // Randomly pick poPoly on polygon.

        var verts = new Array(3 * this.m_nav.getMaxVertsPerPoly());
        var areas = new Array(this.m_nav.getMaxVertsPerPoly());
        arraycopy$4(randomTile.data.verts, randomPoly.verts[0] * 3, verts, 0, 3);

        for (var _j2 = 1; _j2 < randomPoly.vertCount; ++_j2) {
          arraycopy$4(randomTile.data.verts, randomPoly.verts[_j2] * 3, verts, _j2 * 3, 3);
        }

        s = frand.frand();
        t = frand.frand();
        var pt = DetourCommon.randomPointInConvexPoly(verts, randomPoly.vertCount, areas, s, t);
        pt[1] = getPolyHeight(randomPolyRef, pt);
        return new FindRandomPointResult(Status.SUCCSESS, randomPolyRef, pt);
      } //////////////////////////////////////////////////////////////////////////////////////////
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

    }, {
      key: "closestPointOnPoly",
      value: function closestPointOnPoly(ref, pos) {
        var tileAndPoly = this.m_nav.getTileAndPolyByRef(ref);
        var tile = tileAndPoly[0];
        var poly = tileAndPoly[1]; // Off-mesh connections don't have detail polygons.

        if (poly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) {
          var _v = poly.verts[0] * 3;

          var _v2 = poly.verts[1] * 3;

          var _d = DetourCommon.vDist3(pos, tile.data.verts, _v);

          var _d2 = DetourCommon.vDist3(pos, tile.data.verts, _v2);

          var _u = _d / (_d + _d2);

          var _closest = DetourCommon.vLerp4(tile.data.verts, _v, _v2, _u);

          return new ClosesPointOnPolyResult(false, _closest);
        } // Clamp poPoly to be inside the polygon.


        var verts = new Array(this.m_nav.getMaxVertsPerPoly() * 3);
        var edged = new Array(this.m_nav.getMaxVertsPerPoly());
        var edget = new Array(this.m_nav.getMaxVertsPerPoly());
        var nv = poly.vertCount;

        for (var i = 0; i < nv; ++i) {
          arraycopy$4(tile.data.verts, poly.verts[i] * 3, verts, i * 3, 3);
        }

        var posOverPoly = false;
        var closest = new Array(3);
        DetourCommon.vCopy(closest, pos);

        if (!DetourCommon.distancePtPolyEdgesSqr(pos, verts, nv, edged, edget)) {
          // PoPoly is outside the polygon, dtClamp to nearest edge.
          var dmin = edged[0];
          var imin = 0;

          for (var _i2 = 1; _i2 < nv; ++_i2) {
            if (edged[_i2] < dmin) {
              dmin = edged[_i2];
              imin = _i2;
            }
          }

          var va = imin * 3;
          var vb = (imin + 1) % nv * 3;
          closest = DetourCommon.vLerp4(verts, va, vb, edget[imin]);
          posOverPoly = false;
        } else {
          posOverPoly = true;
        }

        var ip = poly.index;

        if (tile.data.detailMeshes != null && tile.data.detailMeshes.length > ip) {
          var pd = tile.data.detailMeshes[ip]; // Find height at the location.

          for (var j = 0; j < pd.triCount; ++j) {
            var _t2 = (pd.triBase + j) * 4;

            var v = new Array(3); //Was new Array(3)[]

            for (var k = 0; k < 3; ++k) {
              if (tile.data.detailTris[_t2 + k] < poly.vertCount) {
                var index = poly.verts[tile.data.detailTris[_t2 + k]] * 3;
                v[k] = [tile.data.verts[index], tile.data.verts[index + 1], tile.data.verts[index + 2]];
              } else {
                var _index = (pd.vertBase + (tile.data.detailTris[_t2 + k] - poly.vertCount)) * 3;

                v[k] = [tile.data.detailVerts[_index], tile.data.detailVerts[_index + 1], tile.data.detailVerts[_index + 2]];
              }
            }

            var heightResult = DetourCommon.closestHeightPointTriangle(closest, v[0], v[1], v[2]);

            if (heightResult[0]) {
              closest[1] = heightResult[1];
              break;
            }
          }
        }

        return new ClosestPointOnPolyResult(posOverPoly, closest);
      } /// @par
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

    }, {
      key: "closestPointOnPolyBoundary",
      value: function closestPointOnPolyBoundary(ref, pos) {
        var tileAndPoly = this.m_nav.getTileAndPolyByRef(ref);
        var tile = tileAndPoly[0];
        var poly = tileAndPoly[1]; // Collect vertices.

        var verts = new Array(this.m_nav.getMaxVertsPerPoly() * 3);
        var edged = new Array(this.m_nav.getMaxVertsPerPoly());
        var edget = new Array(this.m_nav.getMaxVertsPerPoly());
        var nv = poly.vertCount;

        for (var i = 0; i < nv; ++i) {
          arraycopy$4(tile.data.verts, poly.verts[i] * 3, verts, i * 3, 3);
        }

        var closest;

        if (DetourCommon.distancePtPolyEdgesSqr(pos, verts, nv, edged, edget)) {
          closest = DetourCommon.vCopy_return(pos);
        } else {
          // PoPoly is outside the polygon, dtClamp to nearest edge.
          var dmin = edged[0];
          var imin = 0;

          for (var _i3 = 1; _i3 < nv; ++_i3) {
            if (edged[_i3] < dmin) {
              dmin = edged[_i3];
              imin = _i3;
            }
          }

          var va = imin * 3;
          var vb = (imin + 1) % nv * 3;
          closest = DetourCommon.vLerp4(verts, va, vb, edget[imin]);
        }

        return closest;
      } /// @par
      ///
      /// Will return #DT_FAILURE if the provided position is outside the xz-bounds
      /// of the polygon.
      ///
      /// Gets the height of the polygon at the provided position using the height detail. (Most accurate.)
      ///  @param[in]		ref			The reference id of the polygon.
      ///  @param[in]		pos			A position within the xz-bounds of the polygon. [(x, y, z)]
      ///  @param[out]	height		The height at the surface of the polygon.
      /// @returns The status flags for the query.

    }, {
      key: "getPolyHeight",
      value: function getPolyHeight(ref, pos) {
        var tileAndPoly = this.m_nav.getTileAndPolyByRef(ref);
        var tile = tileAndPoly[0];
        var poly = tileAndPoly[1];

        if (poly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) {
          var i = poly.verts[0] * 3;
          i = poly.verts[1] * 3;
          letv1 = [tile.data.verts[i], tile.data.verts[i + 1], tile.data.verts[i + 2]];
          d0 = DetourCommon.vDist2D(pos, v0);
          d1 = DetourCommon.vDist2D(pos, v1);
          u = d0 / (d0 + d1);
          return v0[1] + (v1[1] - v0[1]) * u;
        } else {
          var ip = poly.index;
          var pd = tile.data.detailMeshes[ip];

          for (var j = 0; j < pd.triCount; ++j) {
            var _t3 = (pd.triBase + j) * 4;

            var v = new Array(3); //new Array(3)[];

            for (var k = 0; k < 3; ++k) {
              if (tile.data.detailTris[_t3 + k] < poly.vertCount) {
                var index = poly.verts[tile.data.detailTris[_t3 + k]] * 3;
                v[k] = [tile.data.verts[index], tile.data.verts[index + 1], tile.data.verts[index + 2]];
              } else {
                var _index2 = (pd.vertBase + (tile.data.detailTris[_t3 + k] - poly.vertCount)) * 3;

                v[k] = [tile.data.detailVerts[_index2], tile.data.detailVerts[_index2 + 1], tile.data.detailVerts[_index2 + 2]];
              }
            }

            var heightResult = DetourCommon.closestHeightPointTriangle(pos, v[0], v[1], v[2]);

            if (heightResult[0]) {
              return heightResult[1];
            }
          }
        }

        throw new IllegalArgumentException("Invalid ref " + ref + " pos " + Arrays.toString(pos));
      } /// @par
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

    }, {
      key: "findNearestPoly",
      value: function findNearestPoly(center, extents, filter) {
        var nearestPt = null; // Get nearby polygons from proximity grid.

        var polys = this.queryPolygons(center, extents, filter); // Find nearest polygon amongst the nearby polygons.

        var nearest = 0;
        var nearestDistanceSqr = Number.MAX_VALUE;

        for (var i = 0; i < polys.length; ++i) {
          var ref = polys[i];
          var closest = this.closestPointOnPoly(ref, center);
          var posOverPoly = closest.isPosOverPoly();
          var closestPtPoly = closest.getClosest(); // If a poPoly is directly over a polygon and closer than
          // climb height, favor that instead of straight line nearest point.

          var d = 0;
          var diff = DetourCommon.vSub(center, closestPtPoly);

          if (posOverPoly) {
            var tilaAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(polys[i]);
            var _tile = tilaAndPoly[0];
            d = Math.abs(diff[1]) - _tile.data.header.walkableClimb;
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
    }, {
      key: "and",
      value: function and(v1, v2) {
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
    }, {
      key: "or",
      value: function or(v1, v2) {
        var hi = 0x80000000;
        var low = 0x7fffffff;
        var hi1 = ~~(v1 / hi);
        var hi2 = ~~(v2 / hi);
        var low1 = v1 & low;
        var low2 = v2 & low;
        var h = hi1 | hi2;
        var l = low1 | low2;
        return h * hi + l;
      } // FIXME: (PP) duplicate?

    }, {
      key: "queryPolygonsInTile",
      value: function queryPolygonsInTile(tile, qmin, qmax, filter) {
        var polys = [];

        if (tile.data.bvTree != null) {
          var nodeIndex = 0;
          var tbmin = tile.data.header.bmin;
          var tbmax = tile.data.header.bmax;
          var qfac = tile.data.header.bvQuantFactor; // Calculate quantized box

          var bmin = new Array(3);
          var bmax = new Array(3); // dtClamp query box to world box.

          var minx = DetourCommon.clamp(qmin[0], tbmin[0], tbmax[0]) - tbmin[0];
          var miny = DetourCommon.clamp(qmin[1], tbmin[1], tbmax[1]) - tbmin[1];
          var minz = DetourCommon.clamp(qmin[2], tbmin[2], tbmax[2]) - tbmin[2];
          var maxx = DetourCommon.clamp(qmax[0], tbmin[0], tbmax[0]) - tbmin[0];
          var maxy = DetourCommon.clamp(qmax[1], tbmin[1], tbmax[1]) - tbmin[1];
          var maxz = DetourCommon.clamp(qmax[2], tbmin[2], tbmax[2]) - tbmin[2]; // Quantize

          bmin[0] = this.and(Math.floor(qfac * minx), 0xfffe);
          bmin[1] = this.and(Math.floor(qfac * miny), 0xfffe);
          bmin[2] = this.and(Math.floor(qfac * minz), 0xfffe);
          bmax[0] = this.or(Math.floor(qfac * maxx + 1), 1);
          bmax[1] = this.or(Math.floor(qfac * maxy + 1), 1);
          bmax[2] = this.or(Math.floor(qfac * maxz + 1), 1); // Traverse tree

          var base = this.m_nav.getPolyRefBase(tile);
          var end = tile.data.header.bvNodeCount;

          while (nodeIndex < end) {
            var node = tile.data.bvTree[nodeIndex];
            var overlap = DetourCommon.overlapQuantBounds(bmin, bmax, node.bmin, node.bmax);
            var isLeafNode = node.i >= 0;

            if (isLeafNode && overlap) {
              var ref = this.or(base, node.i);

              if (filter.passFilter(ref, tile, tile.data.polys[node.i])) {
                polys.push(ref);
              }
            }

            if (overlap || isLeafNode) nodeIndex++;else {
              var escapeIndex = -node.i;
              nodeIndex += escapeIndex;
            }
          }

          return polys;
        } else {
          var _bmin = new Array(3);

          var _bmax = new Array(3);

          var _base = this.m_nav.getPolyRefBase(tile);

          for (var i = 0; i < tile.data.header.polyCount; ++i) {
            var p = tile.data.polys[i]; // Do not return off-mesh connection polygons.

            if (p.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) continue;

            var _ref = _base | i;

            if (!filter.passFilter(_ref, tile, p)) continue; // Calc polygon bounds.

            var v = p.verts[0] * 3;
            DetourCommon.vCopy(_bmin, tile.data.verts, v);
            DetourCommon.vCopy(_bmax, tile.data.verts, v);

            for (var j = 1; j < p.vertCount; ++j) {
              v = p.verts[j] * 3;
              DetourCommon.vMin(_bmin, tile.data.verts, v);
              DetourCommon.vMax(_bmax, tile.data.verts, v);
            }

            if (overlapBounds(qmin, qmax, _bmin, _bmax)) {
              polys.push(_ref);
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

    }, {
      key: "queryPolygons",
      value: function queryPolygons(center, extents, filter) {
        var bmin = DetourCommon.vSub(center, extents);
        var bmax = DetourCommon.vAdd(center, extents); // Find tiles the query touches.

        var minxy = this.m_nav.calcTileLoc(bmin);
        var minx = minxy[0];
        var miny = minxy[1];
        var maxxy = this.m_nav.calcTileLoc(bmax);
        var maxx = maxxy[0];
        var maxy = maxxy[1];
        var polys = [];

        for (var y = miny; y <= maxy; ++y) {
          for (var _x = minx; _x <= maxx; ++_x) {
            var neis = this.m_nav.getTilesAt(_x, y);

            for (var j = 0; j < neis.length; ++j) {
              var polysInTile = this.queryPolygonsInTile(neis[j], bmin, bmax, filter);
              polys.push.apply(polys, _toConsumableArray(polysInTile));
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

    }, {
      key: "findPath",
      value: function findPath(startRef, endRef, startPos, endPos, filter) {
        if (startRef == 0 || endRef == 0) throw new IllegalArgumentException("Start or end ref = 0"); // Validate input

        if (!this.m_nav.isValidPolyRef(startRef) || !this.m_nav.isValidPolyRef(endRef)) throw new IllegalArgumentException("Invalid start or end ref");

        if (startRef == endRef) {
          var _path = new Array(1);

          _path.push(startRef);

          return new FindPathResult(Status.SUCCSESS, _path);
        }

        this.m_nodePool = [];
        this.m_openList = [];
        var startNode = this.m_nodePool.getNode(startRef);
        DetourCommon.vCopy(startNode.pos, startPos);
        startNode.pidx = 0;
        startNode.cost = 0;
        startNode.total = DetourCommon.vDist2(startPos, endPos) * NavMeshQuery.H_SCALE;
        startNode.id = startRef;
        startNode.flags = Node.DT_NODE_OPEN;
        this.m_openList.push(startNode);
        var lastBestNode = startNode;
        lastBestNodeCost = startNode.total;
        var status = Status.SUCCSESS;

        while (!this.m_openList.length == 0) {
          // Remove node from open list and put it in closed list.
          var bestNode = this.m_openList.pop();
          bestNode.flags &= ~Node.DT_NODE_OPEN;
          bestNode.flags |= Node.DT_NODE_CLOSED; // Reached the goal, stop searching.

          if (bestNode.id == endRef) {
            lastBestNode = bestNode;
            break;
          } // Get current let and tile.
          // The API input has been cheked already, skip checking internal data.


          var bestRef = bestNode.id;
          var tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(bestRef);
          var bestTile = tileAndPoly[0];
          var bestPoly = tileAndPoly[1]; // Get parent let and tile.

          var parentRef = 0;
          var parentTile = null;
          var parentPoly = null;
          if (bestNode.pidx != 0) parentRef = this.m_nodePool.getNodeAtIdx(bestNode.pidx).id;

          if (parentRef != 0) {
            tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(parentRef);
            parentTile = tileAndPoly[0];
            parentPoly = tileAndPoly[1];
          }

          for (var i = bestPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = bestTile.links[i].next) {
            var neighbourRef = bestTile.links[i].ref; // Skip invalid ids and do not expand back to where we came from.

            if (neighbourRef == 0 || neighbourRef == parentRef) continue; // Get neighbour let and tile.
            // The API input has been cheked already, skip checking internal data.

            tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
            var neighbourTile = tileAndPoly[0];
            var neighbourPoly = tileAndPoly[1];
            if (!filter.passFilter(neighbourRef, neighbourTile, neighbourPoly)) continue; // deal explicitly with crossing tile boundaries

            var crossSide = 0;
            if (bestTile.links[i].side != 0xff) crossSide = bestTile.links[i].side >> 1; // get the node

            var neighbourNode = this.m_nodePool.getNode(neighbourRef, crossSide); // If the node is visited the first time, calculate node position.

            if (neighbourNode.flags == 0) {
              neighbourNode.pos = this.getEdgeMidPoint6(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile);
            } // Calculate cost and heuristic.


            cost = 0;
            heuristic = 0; // Special case for last node.

            if (neighbourRef == endRef) {
              // Cost
              curCost = filter.getCost(bestNode.pos, neighbourNode.pos, parentRef, parentTile, parentPoly, bestRef, bestTile, bestPoly, neighbourRef, neighbourTile, neighbourPoly);
              var endCost = filter.getCost(neighbourNode.pos, endPos, bestRef, bestTile, bestPoly, neighbourRef, neighbourTile, neighbourPoly, 0, null, null);
              cost = bestNode.cost + curCost + endCost;
              heuristic = 0;
            } else {
              // Cost
              curCost = filter.getCost(bestNode.pos, neighbourNode.pos, parentRef, parentTile, parentPoly, bestRef, bestTile, bestPoly, neighbourRef, neighbourTile, neighbourPoly);
              cost = bestNode.cost + curCost;
              heuristic = DetourCommon.vDist2(neighbourNode.pos, endPos) * NavMeshQuery.H_SCALE;
            }

            total = cost + heuristic; // The node is already in open list and the new result is worse, skip.

            if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0 && total >= neighbourNode.total) continue; // The node is already visited and process, and the new result is worse, skip.

            if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0 && total >= neighbourNode.total) continue; // Add or update the node.

            neighbourNode.pidx = this.m_nodePool.getNodeIdx(bestNode);
            neighbourNode.id = neighbourRef;
            neighbourNode.flags = neighbourNode.flags & ~Node.DT_NODE_CLOSED;
            neighbourNode.cost = cost;
            neighbourNode.total = total;

            if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0) {
              // Already in open, update node location.
              this.m_openList.modify(neighbourNode);
            } else {
              // Put the node in open list.
              neighbourNode.flags |= Node.DT_NODE_OPEN;
              this.m_openList.push(neighbourNode);
            } // Update nearest node to target so far.


            if (heuristic < lastBestNodeCost) {
              lastBestNodeCost = heuristic;
              lastBestNode = neighbourNode;
            }
          }
        }

        var path = getPathToNode(lastBestNode);
        if (lastBestNode.id != endRef) status = Status.PARTIAL_RESULT;
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

    }, {
      key: "initSlicedFindPath",
      value: function initSlicedFindPath(startRef, endRef, startPos, endPos, filter, options) {
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
        if (startRef == 0 || endRef == 0) throw new IllegalArgumentException("Start or end ref = 0"); // Validate input

        if (!this.m_nav.isValidPolyRef(startRef) || !this.m_nav.isValidPolyRef(endRef)) throw new IllegalArgumentException("Invalid start or end ref"); // trade quality with performance?

        if ((options & NavMeshQuery.DT_FINDPATH_ANY_ANGLE) != 0) {
          // limiting to several times the character radius yields nice results. It is not sensitive
          // so it is enough to compute it from the first tile.
          var _tile2 = this.m_nav.getTileByRef(startRef);

          agentRadius = _tile2.data.header.walkableRadius;
          this.m_query.raycastLimitSqr = DetourCommon.sqr(agentRadius * NavMesh.DT_RAY_CAST_LIMIT_PROPORTIONS);
        }

        if (startRef == endRef) {
          this.m_query.status = Status.SUCCSESS;
          return Status.SUCCSESS;
        }

        this.m_nodePool.clear();
        this.m_openList.clear();
        var startNode = this.m_nodePool.getNode(startRef);
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

    }, {
      key: "updateSlicedFindPath",
      value: function updateSlicedFindPath(maxIter) {
        if (this.m_query.status != Status.IN_PROGRESS) return new UpdateSlicedPathResult(this.m_query.status, 0); // Make sure the request is still valid.

        if (!this.m_nav.isValidPolyRef(this.m_query.startRef) || !this.m_nav.isValidPolyRef(this.m_query.endRef)) {
          this.m_query.status = Status.FAILURE;
          return new UpdateSlicedPathResult(this.m_query.status, 0);
        }

        var iter = 0;

        while (iter < maxIter && !this.m_openList.isEmpty()) {
          iter++; // Remove node from open list and put it in closed list.

          var bestNode = this.m_openList.pop();
          bestNode.flags &= ~Node.DT_NODE_OPEN;
          bestNode.flags |= Node.DT_NODE_CLOSED; // Reached the goal, stop searching.

          if (bestNode.id == this.m_query.endRef) {
            this.m_query.lastBestNode = bestNode;
            this.m_query.status = Status.SUCCSESS;
            return new UpdateSlicedPathResult(this.m_query.status, iter);
          } // Get current let and tile.
          // The API input has been cheked already, skip checking internal
          // data.


          var bestRef = bestNode.id;
          var tileAndPoly = void 0;

          try {
            tileAndPoly = this.m_nav.getTileAndPolyByRef(bestRef);
          } catch (e) {
            this.m_query.status = Status.FAILURE; // The polygon has disappeared during the sliced query, fail.

            return new UpdateSlicedPathResult(this.m_query.status, iter);
          }

          var bestTile = tileAndPoly[0];
          var bestPoly = tileAndPoly[1]; // Get parent and grand parent let and tile.

          var parentRef = 0,
              grandpaRef = 0;
          var parentTile = null;
          var parentPoly = null;
          var parentNode = null;

          if (bestNode.pidx != 0) {
            parentNode = this.m_nodePool.getNodeAtIdx(bestNode.pidx);
            parentRef = parentNode.id;
            if (parentNode.pidx != 0) grandpaRef = this.m_nodePool.getNodeAtIdx(parentNode.pidx).id;
          }

          if (parentRef != 0) {
            var invalidParent = false;

            try {
              tileAndPoly = this.m_nav.getTileAndPolyByRef(parentRef);
              parentTile = tileAndPoly[0];
              parentPoly = tileAndPoly[1];
            } catch (e) {
              invalidParent = true;
            }

            if (invalidParent || grandpaRef != 0 && !this.m_nav.isValidPolyRef(grandpaRef)) {
              // The polygon has disappeared during the sliced query,
              // fail.
              this.m_query.status = Status.FAILURE;
              return new UpdateSlicedPathResult(this.m_query.status, iter);
            }
          } // decide whether to test raycast to previous nodes


          var tryLOS = false;

          if ((this.m_query.options & NavMeshQuery.DT_FINDPATH_ANY_ANGLE) != 0) {
            if (parentRef != 0 && DetourCommon.vDistSqr(parentNode.pos, bestNode.pos) < this.m_query.raycastLimitSqr) tryLOS = true;
          }

          for (var i = bestPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = bestTile.links[i].next) {
            var neighbourRef = bestTile.links[i].ref; // bestTile.links.forEach(z=>console.log(z.ref));
            // Skip invalid ids and do not expand back to where we came
            // from.

            if (neighbourRef == 0 || neighbourRef == parentRef) continue; // Get neighbour let and tile.
            // The API input has been cheked already, skip checking internal
            // data.

            var tileAndPolyUns = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
            var neighbourTile = tileAndPolyUns[0];
            var neighbourPoly = tileAndPolyUns[1];
            if (!this.m_query.filter.passFilter(neighbourRef, neighbourTile, neighbourPoly)) continue; // get the neighbor node

            var neighbourNode = this.m_nodePool.getNode(neighbourRef, 0); // do not expand to nodes that were already visited from the
            // same parent

            if (neighbourNode.pidx != 0 && neighbourNode.pidx == bestNode.pidx) continue; // If the node is visited the first time, calculate node
            // position.

            if (neighbourNode.flags == 0) {
              neighbourNode.pos = this.getEdgeMidPoint6(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile);
            } // Calculate cost and heuristic.


            var _cost = 0;
            var _heuristic = 0; // raycast parent

            var foundShortCut = false;

            if (tryLOS) {
              var rayHit = raycast(parentRef, parentNode.pos, neighbourNode.pos, this.m_query.filter, NavMeshQuery.DT_RAYCAST_USE_COSTS, grandpaRef);
              foundShortCut = rayHit.t >= 1.0;

              if (foundShortCut) {
                // shortcut found using raycast. Using shorter cost
                // instead
                _cost = parentNode.cost + rayHit.pathCost;
              }
            } // update move cost


            if (!foundShortCut) {
              // No shortcut found.
              var _curCost = this.m_query.filter.getCost(bestNode.pos, neighbourNode.pos, parentRef, parentTile, parentPoly, bestRef, bestTile, bestPoly, neighbourRef, neighbourTile, neighbourPoly);

              _cost = bestNode.cost + _curCost;
            } // Special case for last node.


            if (neighbourRef == this.m_query.endRef) {
              var endCost = this.m_query.filter.getCost(neighbourNode.pos, this.m_query.endPos, bestRef, bestTile, bestPoly, neighbourRef, neighbourTile, neighbourPoly, 0, null, null);
              _cost = _cost + endCost;
              _heuristic = 0;
            } else {
              _heuristic = DetourCommon.vDist2(neighbourNode.pos, this.m_query.endPos) * NavMeshQuery.H_SCALE;
            }

            var _total = _cost + _heuristic; // The node is already in open list and the new result is worse,
            // skip.


            if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0 && _total >= neighbourNode.total) continue; // The node is already visited and process, and the new result
            // is worse, skip.

            if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0 && _total >= neighbourNode.total) continue; // Add or update the node.

            neighbourNode.pidx = foundShortCut ? bestNode.pidx : this.m_nodePool.getNodeIdx(bestNode);
            neighbourNode.id = neighbourRef;
            neighbourNode.flags = this.and(neighbourNode.flags, ~this.or(Node.DT_NODE_CLOSED, Node.DT_NODE_PARENT_DETACHED));
            neighbourNode.cost = _cost;
            neighbourNode.total = _total;
            if (foundShortCut) neighbourNode.flags = this.or(neighbourNode.flags, Node.DT_NODE_PARENT_DETACHED);

            if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0) {
              // Already in open, update node location.
              this.m_openList.modify(neighbourNode);
            } else {
              // Put the node in open list.
              neighbourNode.flags |= Node.DT_NODE_OPEN;
              this.m_openList.push(neighbourNode);
            } // Update nearest node to target so far.


            if (_heuristic < this.m_query.lastBestNodeCost) {
              this.m_query.lastBestNodeCost = _heuristic;
              this.m_query.lastBestNode = neighbourNode;
            }
          }
        } // Exhausted all nodes, but could not find path.


        if (this.m_openList.isEmpty()) {
          this.m_query.status = Status.PARTIAL_RESULT;
        }

        return new UpdateSlicedPathResult(this.m_query.status, iter);
      } /// Finalizes and returns the results of a sliced path query.
      ///  @param[out]	path		An ordered list of polygon references representing the path. (Start to end.) 
      ///  							[(polyRef) * @p pathCount]
      /// @returns The status flags for the query.

    }, {
      key: "finalizeSlicedFindPath",
      value: function finalizeSlicedFindPath() {
        var path = [];

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
          if (this.m_query.lastBestNode.id != this.m_query.endRef) this.m_query.status = Status.PARTIAL_RESULT;
          var prev = null;
          var node = this.m_query.lastBestNode;
          var prevRay = 0;

          do {
            var next = this.m_nodePool.getNodeAtIdx(node.pidx);
            node.pidx = this.m_nodePool.getNodeIdx(prev);
            prev = node;
            var nextRay = node.flags & Node.DT_NODE_PARENT_DETACHED; // keep track of whether parent is not adjacent (i.e. due to raycast shortcut)

            node.flags = this.or(this.and(node.flags, ~Node.DT_NODE_PARENT_DETACHED), prevRay); // and store it in the reversed path's node

            prevRay = nextRay;
            node = next;
          } while (node != null); // Store path


          node = prev;

          do {
            var _next = this.m_nodePool.getNodeAtIdx(node.pidx);

            if ((node.flags & Node.DT_NODE_PARENT_DETACHED) != 0) {
              var iresult = raycast(node.id, node.pos, _next.pos, this.m_query.filter, 0, 0);
              path.addAll(iresult.path); // raycast ends on let boundary and the path might include the next let boundary.

              if (path[path.length - 1] == _next.id) path.remove(path.length - 1); // remove to aduplicates
            } else {
              path.push(node.id);
            }

            node = _next;
          } while (node != null);
        }

        var status = this.m_query.status; // Reset query.

        this.m_query = new QueryData();
        return new FindPathResult(status, path);
      } /// Finalizes and returns the results of an incompPolye sliced path query, returning the path to the furthest
      /// polygon on the existing path that was visited during the search.
      ///  @param[in]		existing		An array of polygon references for the existing path.
      ///  @param[in]		existingSize	The number of polygon in the @p existing array.
      ///  @param[out]	path			An ordered list of polygon references representing the path. (Start to end.) 
      ///  								[(polyRef) * @p pathCount]
      /// @returns The status flags for the query.

    }, {
      key: "finalizeSlicedFindPathPartial",
      value: function finalizeSlicedFindPathPartial(existing) {
        var path = [];

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
          var prev = null;
          var node = null;

          for (var i = existing.length - 1; i >= 0; --i) {
            node = this.m_nodePool.findNode(existing[i]);
            if (node != null) break;
          }

          if (node == null) {
            this.m_query.status = Status.PARTIAL_RESULT;
            node = this.m_query.lastBestNode;
          } // Reverse the path.


          var prevRay = 0;

          do {
            var next = this.m_nodePool.getNodeAtIdx(node.pidx);
            node.pidx = this.m_nodePool.getNodeIdx(prev);
            prev = node;
            var nextRay = this.and(node.flags, Node.DT_NODE_PARENT_DETACHED); // keep track of whether parent is not adjacent (i.e. due to raycast shortcut)

            node.flags = this.or(this.and(node.flags & ~Node.DT_NODE_PARENT_DETACHED), prevRay); // and store it in the reversed path's node

            prevRay = nextRay;
            node = next;
          } while (node != null); // Store path


          node = prev;

          do {
            var _next2 = this.m_nodePool.getNodeAtIdx(node.pidx);

            if ((node.flags & Node.DT_NODE_PARENT_DETACHED) != 0) {
              var iresult = raycast(node.id, node.pos, _next2.pos, this.m_query.filter, 0, 0);
              path.addAll(iresult.path); // raycast ends on let boundary and the path might include the next let boundary.

              if (path[path.length - 1] == _next2.id) path.remove(path.length - 1); // remove to aduplicates
            } else {
              path.push(node.id);
            }

            node = _next2;
          } while (node != null);
        }

        var status = this.m_query.status; // Reset query.

        this.m_query = new QueryData();
        return new FindPathResult(status, path);
      }
    }, {
      key: "appendVertex",
      value: function appendVertex(pos, flags, ref, straightPath, maxStraightPath) {
        if (straightPath.length > 0 && DetourCommon.vEqual(straightPath[straightPath.length - 1].pos, pos)) {
          // The vertices are equal, update flags and poly.
          straightPath[straightPath.length - 1].flags = flags;
          straightPath[straightPath.length - 1].ref = ref;
        } else {
          if (straightPath.length < maxStraightPath) {
            // Append new vertex.
            straightPath.push(new StraightPathItem(pos, flags, ref));
          } // If reached end of path or there is no space to append more vertices, return.


          if (flags == NavMeshQuery.DT_STRAIGHTPATH_END || straightPath.length >= maxStraightPath) {
            return Status.SUCCSESS;
          }
        }

        return Status.IN_PROGRESS;
      }
    }, {
      key: "appendPortals",
      value: function appendPortals(startIdx, endIdx, endPos, path, straightPath, maxStraightPath, options) {
        var startPos = straightPath[straightPath.length - 1].pos; // Append or update last vertex

        var stat = null;

        for (var i = startIdx; i < endIdx; i++) {
          // Calculate portal
          var from = path[i];
          var tileAndPoly = this.m_nav.getTileAndPolyByRef(from);
          var fromTile = tileAndPoly[0];
          var fromPoly = tileAndPoly[1];
          var to = path[i + 1];
          tileAndPoly = this.m_nav.getTileAndPolyByRef(to);
          var toTile = tileAndPoly[0];
          var toPoly = tileAndPoly[1];
          var portals = this.getPortalPoints7(from, fromPoly, fromTile, to, toPoly, toTile, 0, 0);
          var left = portals.left;
          var right = portals.right;

          if ((options & NavMeshQuery.DT_STRAIGHTPATH_AREA_CROSSINGS) != 0) {
            // Skip intersection if only area crossings are requested.
            if (fromPoly.getArea() == toPoly.getArea()) continue;
          } // Append intersection


          var interect = DetourCommon.intersectSegSeg2D(startPos, endPos, left, right);

          if (interect[0]) {
            t = interect.third;
            var pt = DetourCommon.vLerp3(left, right, t);
            stat = this.appendVertex(pt, 0, path[i + 1], straightPath, maxStraightPath);
            if (!stat == Status.IN_PROGRESS) return stat;
          }
        }

        return Status.IN_PROGRESS;
      } /// @par
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

    }, {
      key: "findStraightPath",
      value: function findStraightPath(startPos, endPos, path, maxStraightPath, options) {
        if (path.length == 0) {
          throw new IllegalArgumentException("Empty path");
        } // TODO: Should this be callers responsibility?


        var closestStartPos = this.closestPointOnPolyBoundary(path[0], startPos);
        var closestEndPos = this.closestPointOnPolyBoundary(path[path.length - 1], endPos);
        var straightPath = []; // Add start point.

        var stat = this.appendVertex(closestStartPos, NavMeshQuery.DT_STRAIGHTPATH_START, path[0], straightPath, maxStraightPath);
        if (!stat == Status.IN_PROGRESS) return straightPath;

        if (path.length > 1) {
          var portalApex = DetourCommon.vCopy_return(closestStartPos);
          var portalLeft = DetourCommon.vCopy_return(portalApex);
          var portalRight = DetourCommon.vCopy_return(portalApex);
          var apexIndex = 0;
          var leftIndex = 0;
          var rightIndex = 0;
          var leftPolyType = 0;
          var rightPolyType = 0;
          var leftPolyRef = path[0];
          var rightPolyRef = path[0];

          for (var i = 0; i < path.length; ++i) {
            var left = void 0;
            var right = void 0;
            var toType = void 0;

            if (i + 1 < path.length) {
              // Next portal.
              try {
                var portalPoints = this.getPortalPoints2(path[i], path[i + 1]);
                left = portalPoints.left;
                right = portalPoints.right;
                toType = portalPoints.toType;
              } catch (e) {
                closestEndPos = this.closestPointOnPolyBoundary(path[i], endPos); // Append portals aPoly the current straight path segment.

                if (this.and(options, this.or(NavMeshQuery.DT_STRAIGHTPATH_AREA_CROSSINGS, NavMeshQuery.DT_STRAIGHTPATH_ALL_CROSSINGS)) != 0) {
                  stat = appendPortals(apexIndex, i, closestEndPos, path, straightPath, options, maxStraightPath);
                  if (!stat == Status.IN_PROGRESS) return straightPath;
                }

                this.appendVertex(closestEndPos, 0, path[i], straightPath, maxStraightPath);
                return straightPath;
              } // If starting really close the portal, advance.


              if (i == 0) {
                var dt = DetourCommon.distancePtSegSqr2D3(portalApex, left, right);
                if (dt[0] < DetourCommon.sqr(0.001)) continue;
              }
            } else {
              // End of the path.
              left = DetourCommon.vCopy_return(closestEndPos);
              right = DetourCommon.vCopy_return(closestEndPos);
              toType = Poly.DT_POLYTYPE_GROUND;
            } // Right vertex.


            if (DetourCommon.triArea2D3(portalApex, portalRight, right) <= 0.0) {
              if (DetourCommon.vEqual(portalApex, portalRight) || DetourCommon.triArea2D3(portalApex, portalLeft, right) > 0.0) {
                portalRight = DetourCommon.vCopy_return(right);
                rightPolyRef = i + 1 < path.length ? path[i + 1] : 0;
                rightPolyType = toType;
                rightIndex = i;
              } else {
                // Append portals aPoly the current straight path segment.
                if (this.and(options, this.or(NavMeshQuery.DT_STRAIGHTPATH_AREA_CROSSINGS, NavMeshQuery.DT_STRAIGHTPATH_ALL_CROSSINGS)) != 0) {
                  stat = appendPortals(apexIndex, leftIndex, portalLeft, path, straightPath, options, maxStraightPath);
                  if (!stat == Status.IN_PROGRESS) return straightPath;
                }

                portalApex = DetourCommon.vCopy_return(portalLeft);
                apexIndex = leftIndex;
                var flags = 0;
                if (leftPolyRef == 0) flags = NavMeshQuery.DT_STRAIGHTPATH_END;else if (leftPolyType == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) flags = NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION;
                var ref = leftPolyRef; // Append or update vertex

                stat = this.appendVertex(portalApex, flags, ref, straightPath, maxStraightPath);
                if (!stat == Status.IN_PROGRESS) return straightPath;
                portalLeft = DetourCommon.vCopy_return(portalApex);
                portalRight = DetourCommon.vCopy_return(portalApex);
                leftIndex = apexIndex;
                rightIndex = apexIndex; // Restart

                i = apexIndex;
                continue;
              }
            } // Left vertex.


            if (DetourCommon.triArea2D3(portalApex, portalLeft, left) >= 0.0) {
              if (DetourCommon.vEqual(portalApex, portalLeft) || DetourCommon.triArea2D3(portalApex, portalRight, left) < 0.0) {
                portalLeft = DetourCommon.vCopy_return(left);
                leftPolyRef = i + 1 < path.length ? path[i + 1] : 0;
                leftPolyType = toType;
                leftIndex = i;
              } else {
                // Append portals aPoly the current straight path segment.
                if (this.and(options & this.or(NavMeshQuery.DT_STRAIGHTPATH_AREA_CROSSINGS, NavMeshQuery.DT_STRAIGHTPATH_ALL_CROSSINGS)) != 0) {
                  stat = appendPortals(apexIndex, rightIndex, portalRight, path, straightPath, options, maxStraightPath);
                  if (!stat == Status.IN_PROGRESS) return straightPath;
                }

                portalApex = DetourCommon.vCopy_return(portalRight);
                apexIndex = rightIndex;
                var _flags = 0;
                if (rightPolyRef == 0) _flags = NavMeshQuery.DT_STRAIGHTPATH_END;else if (rightPolyType == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) _flags = NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION;
                var _ref2 = rightPolyRef; // Append or update vertex

                stat = this.appendVertex(portalApex, _flags, _ref2, straightPath, maxStraightPath);
                if (!stat == Status.IN_PROGRESS) return straightPath;
                portalLeft = DetourCommon.vCopy_return(portalApex);
                portalRight = DetourCommon.vCopy_return(portalApex);
                leftIndex = apexIndex;
                rightIndex = apexIndex; // Restart

                i = apexIndex;
                continue;
              }
            }
          } // Append portals aPoly the current straight path segment.


          if (this.and(options & this.or(NavMeshQuery.DT_STRAIGHTPATH_AREA_CROSSINGS, NavMeshQuery.DT_STRAIGHTPATH_ALL_CROSSINGS)) != 0) {
            stat = appendPortals(apexIndex, path.length - 1, closestEndPos, path, straightPath, options, maxStraightPath);
            if (!stat == Status.IN_PROGRESS) return straightPath;
          }
        }

        this.appendVertex(closestEndPos, NavMeshQuery.DT_STRAIGHTPATH_END, 0, straightPath, maxStraightPath);
        return straightPath;
      } /// @par
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

    }, {
      key: "moveAlongSurface",
      value: function moveAlongSurface(startRef, startPos, endPos, filter) {
        // Validate input
        if (startRef == 0) throw new IllegalArgumentException("Start ref = 0");
        if (!this.m_nav.isValidPolyRef(startRef)) throw new IllegalArgumentException("Invalid start ref");
        this.m_tinyNodePool = new NodePool();
        var startNode = this.m_tinyNodePool.getNode(startRef);
        startNode.pidx = 0;
        startNode.cost = 0;
        startNode.total = 0;
        startNode.id = startRef;
        startNode.flags = Node.DT_NODE_CLOSED;
        var stack = [];
        stack.push(startNode);
        var bestPos = new Array(3);
        var bestDist = Number.MAX_VALUE;
        var bestNode = null;
        DetourCommon.vCopy(bestPos, startPos); // Search constraints

        var searchPos = DetourCommon.vLerp3(startPos, endPos, 0.5);
        var searchRadSqr = DetourCommon.sqr(DetourCommon.vDist2(startPos, endPos) / 2.0 + 0.001);
        var verts = new Array(this.m_nav.getMaxVertsPerPoly() * 3).fill(0);

        while (!stack.length == 0) {
          // Pop front.
          var curNode = stack.pop(); // Get let and tile.
          // The API input has been cheked already, skip checking internal data.

          var curRef = curNode.id;
          var tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(curRef);
          var curTile = tileAndPoly[0];
          var curPoly = tileAndPoly[1]; // Collect vertices.

          var nverts = curPoly.vertCount;

          for (var i = 0; i < nverts; ++i) {
            arraycopy$4(curTile.data.verts, curPoly.verts[i] * 3, verts, i * 3, 3);
          } // If target is inside the poly, stop search.


          if (DetourCommon.pointInPolygon(endPos, verts, nverts)) {
            bestNode = curNode;
            DetourCommon.vCopy(bestPos, endPos);
            break;
          } // Find wall edges and find nearest poPoly inside the walls.


          for (var _i4 = 0, j = curPoly.vertCount - 1; _i4 < curPoly.vertCount; j = _i4++) {
            // Find links to neighbours.
            var MAX_NEIS = 8;
            var nneis = 0;
            var neis = new Array(MAX_NEIS);

            if ((curPoly.neis[j] & NavMesh.DT_EXT_LINK) != 0) {
              // Tile border.
              for (var k = curPoly.firstLink; k != NavMesh.DT_NULL_LINK; k = curTile.links[k].next) {
                var link = curTile.links[k];

                if (link.edge == j) {
                  if (link.ref != 0) {
                    tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(link.ref);
                    var neiTile = tileAndPoly[0];
                    var neiPoly = tileAndPoly[1];

                    if (filter.passFilter(link.ref, neiTile, neiPoly)) {
                      if (nneis < MAX_NEIS) neis[nneis++] = link.ref;
                    }
                  }
                }
              }
            } else if (curPoly.neis[j] != 0) {
              var idx = curPoly.neis[j] - 1;
              var ref = or(this.m_nav.getPolyRefBase(curTile), idx);

              if (filter.passFilter(ref, curTile, curTile.data.polys[idx])) {
                // Internal edge, encode id.
                neis[nneis++] = ref;
              }
            }

            if (nneis == 0) {
              // Wall edge, calc distance.
              var vj = j * 3;
              var vi = _i4 * 3;
              var distSeg = DetourCommon.distancePtSegSqr2D4(endPos, verts, vj, vi);
              var _distSqr = distSeg[0];
              var _tseg = distSeg[1];

              if (_distSqr < bestDist) {
                // Update nearest distance.
                bestPos = DetourCommon.vLerp4(verts, vj, vi, _tseg);
                bestDist = _distSqr;
                bestNode = curNode;
              }
            } else {
              for (var _k = 0; _k < nneis; ++_k) {
                // Skip if no node can be allocated.
                var neighbourNode = this.m_tinyNodePool.getNode(neis[_k]);
                if (neighbourNode == null) continue; // Skip if already visited.

                if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0) continue; // Skip the link if it is too far from search constraint.
                // TODO: Maybe should use getPortalPoints(), but this one is way faster.

                var _vj = j * 3;

                var _vi = _i4 * 3;

                var distseg = DetourCommon.distancePtSegSqr2D4(searchPos, verts, _vj, _vi);
                var _distSqr2 = distseg[0];
                if (_distSqr2 > searchRadSqr) continue; // Mark as the node as visited and push to queue.

                neighbourNode.pidx = this.m_tinyNodePool.getNodeIdx(curNode);
                neighbourNode.flags |= Node.DT_NODE_CLOSED;
                stack.push(neighbourNode);
              }
            }
          }
        }

        var visited = [];

        if (bestNode != null) {
          // Reverse the path.
          var prev = null;
          var node = bestNode;

          do {
            var next = this.m_tinyNodePool.getNodeAtIdx(node.pidx);
            node.pidx = this.m_tinyNodePool.getNodeIdx(prev);
            prev = node;
            node = next;
          } while (node != null); // Store result


          node = prev;

          do {
            visited.push(node.id);
            node = this.m_tinyNodePool.getNodeAtIdx(node.pidx);
          } while (node != null);
        }

        return new MoveAlongSurfaceResult(bestPos, visited);
      }
    }, {
      key: "getPortalPoints2",
      value: function getPortalPoints2(from, to) {
        var tileAndPoly = this.m_nav.getTileAndPolyByRef(from);
        var fromTile = tileAndPoly[0];
        var fromPoly = tileAndPoly[1];
        var fromType = fromPoly.getType();
        tileAndPoly = this.m_nav.getTileAndPolyByRef(to);
        var toTile = tileAndPoly[0];
        var toPoly = tileAndPoly[1];
        var toType = toPoly.getType();
        return this.getPortalPoints7(from, fromPoly, fromTile, to, toPoly, toTile, fromType, toType);
      } // Returns portal points between two polygons.

    }, {
      key: "getPortalPoints7",
      value: function getPortalPoints7(from, fromPoly, fromTile, to, toPoly, toTile, fromType, toType) {
        var left = new Array(3);
        var right = new Array(3); // Find the link that points to the 'to' polygon.

        var link = null;

        for (var i = fromPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = fromTile.links[i].next) {
          if (fromTile.links[i].ref == to) {
            link = fromTile.links[i];
            break;
          }
        }

        if (link == null) throw new IllegalArgumentException("Null link"); // Handle off-mesh connections.

        if (fromPoly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) {
          // Find link that points to first vertex.
          for (var _i5 = fromPoly.firstLink; _i5 != NavMesh.DT_NULL_LINK; _i5 = fromTile.links[_i5].next) {
            if (fromTile.links[_i5].ref == to) {
              var v = fromTile.links[_i5].edge;
              arraycopy$4(fromTile.data.verts, fromPoly.verts[v] * 3, left, 0, 3);
              arraycopy$4(fromTile.data.verts, fromPoly.verts[v] * 3, right, 0, 3);
              return new PortalResult(left, right, fromType, toType);
            }
          }

          throw new IllegalArgumentException("Invalid offmesh from connection");
        }

        if (toPoly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) {
          for (var _i6 = toPoly.firstLink; _i6 != NavMesh.DT_NULL_LINK; _i6 = toTile.links[_i6].next) {
            if (toTile.links[_i6].ref == from) {
              var _v3 = toTile.links[_i6].edge;
              arraycopy$4(toTile.data.verts, toPoly.verts[_v3] * 3, left, 0, 3);
              arraycopy$4(toTile.data.verts, toPoly.verts[_v3] * 3, right, 0, 3);
              return new PortalResult(left, right, fromType, toType);
            }
          }

          throw new IllegalArgumentException("Invalid offmesh to connection");
        } // Find portal vertices.


        var v0 = fromPoly.verts[link.edge];
        var v1 = fromPoly.verts[(link.edge + 1) % fromPoly.vertCount];
        arraycopy$4(fromTile.data.verts, v0 * 3, left, 0, 3);
        arraycopy$4(fromTile.data.verts, v1 * 3, right, 0, 3); // If the link is at tile boundary, dtClamp the vertices to
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
      } // Returns edge mid poPoly between two polygons.

    }, {
      key: "getEdgeMidPoint2",
      value: function getEdgeMidPoint2(from, to) {
        var ppoints = this.getPortalPoints2(from, to);
        var left = ppoints.left;
        var right = ppoints.right;
        var mid = new Array(3);
        mid[0] = (left[0] + right[0]) * 0.5;
        mid[1] = (left[1] + right[1]) * 0.5;
        mid[2] = (left[2] + right[2]) * 0.5;
        return mid;
      }
    }, {
      key: "getEdgeMidPoint6",
      value: function getEdgeMidPoint6(from, fromPoly, fromTile, to, toPoly, toTile) {
        var ppoints = this.getPortalPoints7(from, fromPoly, fromTile, to, toPoly, toTile, 0, 0);
        var left = ppoints.left;
        var right = ppoints.right;
        var mid = new Array(3);
        mid[0] = (left[0] + right[0]) * 0.5;
        mid[1] = (left[1] + right[1]) * 0.5;
        mid[2] = (left[2] + right[2]) * 0.5;
        return mid;
      }
    }, {
      key: "raycast",
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
      value: function raycast(startRef, startPos, endPos, filter, options, prevRef) {
        // Validate input
        if (startRef == 0 || !this.m_nav.isValidPolyRef(startRef)) throw new IllegalArgumentException("Invalid start ref");
        if (prevRef != 0 && !this.m_nav.isValidPolyRef(prevRef)) throw new IllegalArgumentException("Invalid pref ref");
        var hit = new RaycastHit();
        var verts = new Array(this.m_nav.getMaxVertsPerPoly() * 3 + 3);
        var curPos = new Array(3),
            lastPos = new Array(3);
        DetourCommon.vCopy(curPos, startPos);
        var dir = DetourCommon.vSub(endPos, startPos);
        var prevTile, tile, nextTile;
        var prevPoly, poly, nextPoly; // The API input has been checked already, skip checking internal data.

        var curRef = startRef;
        var tileAndPolyUns = this.m_nav.getTileAndPolyByRefUnsafe(curRef);
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
          var nv = 0;

          for (var i = 0; i < poly.vertCount; ++i) {
            arraycopy$4(tile.data.verts, poly.verts[i] * 3, verts, nv * 3, 3);
            nv++;
          }

          var iresult = DetourCommon.intersectSegmentPoly2D(startPos, endPos, verts, nv);

          if (!iresult.intersects) {
            // Could not hit the polygon, keep the old t and report hit.
            return hit;
          }

          hit.hitEdgeIndex = iresult.segMax; // Keep track of furthest t so far.

          if (iresult.tmax > hit.t) hit.t = iresult.tmax; // Store visited polygons.

          hit.path.push(curRef); // Ray end is compPolyely inside the polygon.

          if (iresult.segMax == -1) {
            hit.t = Number.MAX_VALUE; // add the cost

            if ((options & NavMeshQuery.DT_RAYCAST_USE_COSTS) != 0) hit.pathCost += filter.getCost(curPos, endPos, prevRef, prevTile, prevPoly, curRef, tile, poly, curRef, tile, poly);
            return hit;
          } // Follow neighbours.


          var nextRef = 0;

          for (var _i7 = poly.firstLink; _i7 != NavMesh.DT_NULL_LINK; _i7 = tile.links[_i7].next) {
            var link = tile.links[_i7]; // Find link which contains this edge.

            if (link.edge != iresult.segMax) continue; // Get pointer to the next polygon.

            tileAndPolyUns = this.m_nav.getTileAndPolyByRefUnsafe(link.ref);
            nextTile = tileAndPolyUns[0];
            nextPoly = tileAndPolyUns[1]; // Skip off-mesh connections.

            if (nextPoly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) continue; // Skip links based on filter.

            if (!filter.passFilter(link.ref, nextTile, nextPoly)) continue; // If the link is internal, just return the ref.

            if (link.side == 0xff) {
              nextRef = link.ref;
              break;
            } // If the link is at tile boundary,
            // Check if the link spans the whole edge, and accept.


            if (link.bmin == 0 && link.bmax == 255) {
              nextRef = link.ref;
              break;
            } // Check for partial edge links.


            var _v4 = poly.verts[link.edge];
            var _v5 = poly.verts[(link.edge + 1) % poly.vertCount];
            var left = _v4 * 3;
            var right = _v5 * 3; // Check that the intersection lies inside the link portal.

            if (link.side == 0 || link.side == 4) {
              // Calculate link size.
              lmin = tile.data.verts[left + 2] + (tile.data.verts[right + 2] - tile.data.verts[left + 2]) * (link.bmin * s);
              lmax = tile.data.verts[left + 2] + (tile.data.verts[right + 2] - tile.data.verts[left + 2]) * (link.bmax * s);

              if (lmin > lmax) {
                temp = lmin;
                lmin = lmax;
                lmax = temp;
              } // Find Z intersection.


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
              } // Find X intersection.


              x = startPos[0] + (endPos[0] - startPos[0]) * iresult.tmax;

              if (x >= lmin && x <= lmax) {
                nextRef = link.ref;
                break;
              }
            }
          } // add the cost


          if ((options & NavMeshQuery.DT_RAYCAST_USE_COSTS) != 0) {
            // compute the intersection poPoly at the furthest end of the polygon
            // and correct the height (since the raycast moves in 2d)
            DetourCommon.vCopy(lastPos, curPos);
            curPos = DetourCommon.vMad(startPos, dir, hit.t);
            var e1 = new VectorPtr(verts, iresult.segMax * 3);
            var e2 = new VectorPtr(verts, (iresult.segMax + 1) % nv * 3);
            var eDir = DetourCommon.vSub(e2, e1);
            var diff = DetourCommon.vSub(new VectorPtr(curPos), e1);
            s = DetourCommon.sqr(eDir[0]) > DetourCommon.sqr(eDir[2]) ? diff[0] / eDir[0] : diff[2] / eDir[2];
            curPos[1] = e1[1] + eDir[1] * s;
            hit.pathCost += filter.getCost(lastPos, curPos, prevRef, prevTile, prevPoly, curRef, tile, poly, nextRef, nextTile, nextPoly);
          }

          if (nextRef == 0) {
            // No neighbour, we hit a wall.
            // Calculate hit normal.
            var a = iresult.segMax;
            var b = iresult.segMax + 1 < nv ? iresult.segMax + 1 : 0;
            var va = a * 3;
            var vb = b * 3;
            dx = verts[vb] - verts[va];
            dz = verts[vb + 2] - verts[va + 2];
            hit.hitNormal[0] = dz;
            hit.hitNormal[1] = 0;
            hit.hitNormal[2] = -dx;
            vNormalize(hit.hitNormal);
            return hit;
          } // No hit, advance to neighbour polygon.


          prevRef = curRef;
          curRef = nextRef;
          prevTile = tile;
          tile = nextTile;
          prevPoly = poly;
          poly = nextPoly;
        }

        return hit;
      } /// @par
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

    }, {
      key: "findPolysAroundCircle",
      value: function findPolysAroundCircle(startRef, centerPos, radius, filter) {
        // Validate input
        if (startRef == 0 || !this.m_nav.isValidPolyRef(startRef)) throw new IllegalArgumentException("Invalid start ref");
        var resultRef = [];
        var resultParent = [];
        var resultCost = [];
        this.m_nodePool = [];
        this.m_openList = [];
        var startNode = this.m_nodePool.getNode(startRef);
        DetourCommon.vCopy(startNode.pos, centerPos);
        startNode.pidx = 0;
        startNode.cost = 0;
        startNode.total = 0;
        startNode.id = startRef;
        startNode.flags = Node.DT_NODE_OPEN;
        this.m_openList.push(startNode);
        radiusSqr = DetourCommon.sqr(radius);

        while (!this.m_openList.length == 0) {
          var bestNode = m_openList.pop();
          bestNode.flags &= ~Node.DT_NODE_OPEN;
          bestNode.flags |= Node.DT_NODE_CLOSED; // Get let and tile.
          // The API input has been cheked already, skip checking internal data.

          var bestRef = bestNode.id;
          var tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(bestRef);
          var bestTile = tileAndPoly[0];
          var bestPoly = tileAndPoly[1]; // Get parent let and tile.

          var parentRef = 0;
          var parentTile = null;
          var parentPoly = null;
          if (bestNode.pidx != 0) parentRef = this.m_nodePool.getNodeAtIdx(bestNode.pidx).id;

          if (parentRef != 0) {
            tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(parentRef);
            parentTile = tileAndPoly[0];
            parentPoly = tileAndPoly[1];
          }

          resultRef.push(bestRef);
          resultParent.push(parentRef);
          resultCost.push(bestNode.total);

          for (var i = bestPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = bestTile.links[i].next) {
            var link = bestTile.links[i];
            var neighbourRef = link.ref; // Skip invalid neighbours and do not follow back to parent.

            if (neighbourRef == 0 || neighbourRef == parentRef) continue; // Expand to neighbour

            tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
            var neighbourTile = tileAndPoly[0];
            var neighbourPoly = tileAndPoly[1]; // Do not advance if the polygon is excluded by the filter.

            if (!filter.passFilter(neighbourRef, neighbourTile, neighbourPoly)) continue; // Find edge and calc distance to the edge.

            var pp = this.getPortalPoints7(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile, 0, 0);
            var va = pp.left;
            var vb = pp.right; // If the circle is not touching the next polygon, skip it.

            var distseg = DetourCommon.distancePtSegSqr2D3(centerPos, va, vb);
            distSqr = distseg[0];
            if (distSqr > radiusSqr) continue;
            var neighbourNode = this.m_nodePool.getNode(neighbourRef);
            if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0) continue; // Cost

            if (neighbourNode.flags == 0) neighbourNode.pos = vLerp3(va, vb, 0.5);
            cost = filter.getCost(bestNode.pos, neighbourNode.pos, parentRef, parentTile, parentPoly, bestRef, bestTile, bestPoly, neighbourRef, neighbourTile, neighbourPoly);
            total = bestNode.total + cost; // The node is already in open list and the new result is worse, skip.

            if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0 && total >= neighbourNode.total) continue;
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
      } /// @par
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

    }, {
      key: "findPolysAroundShape",
      value: function findPolysAroundShape(startRef, verts, nverts, filter) {
        // Validate input
        if (startRef == 0 || !this.m_nav.isValidPolyRef(startRef)) throw new IllegalArgumentException("Invalid start ref");
        var resultRef = [];
        var resultParent = [];
        var resultCost = [];
        this.m_nodePool = [];
        this.m_openList = [];
        var centerPos = [0, 0, 0];

        for (var i = 0; i < nverts; ++i) {
          centerPos[0] += verts[i * 3];
          centerPos[1] += verts[i * 3 + 1];
          centerPos[2] += verts[i * 3 + 2];
        }

        scale = 1.0 / nverts;
        centerPos[0] *= scale;
        centerPos[1] *= scale;
        centerPos[2] *= scale;
        var startNode = this.m_nodePool.getNode(startRef);
        DetourCommon.vCopy(startNode.pos, centerPos);
        startNode.pidx = 0;
        startNode.cost = 0;
        startNode.total = 0;
        startNode.id = startRef;
        startNode.flags = Node.DT_NODE_OPEN;
        this.m_openList.push(startNode);

        while (!this.m_openList.length == 0) {
          var bestNode = this.m_openList.pop();
          bestNode.flags &= ~Node.DT_NODE_OPEN;
          bestNode.flags |= Node.DT_NODE_CLOSED; // Get let and tile.
          // The API input has been cheked already, skip checking internal data.

          var bestRef = bestNode.id;
          var tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(bestRef);
          var bestTile = tileAndPoly[0];
          var bestPoly = tileAndPoly[1]; // Get parent let and tile.

          var parentRef = 0;
          var parentTile = null;
          var parentPoly = null;
          if (bestNode.pidx != 0) parentRef = this.m_nodePool.getNodeAtIdx(bestNode.pidx).id;

          if (parentRef != 0) {
            tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(parentRef);
            parentTile = tileAndPoly[0];
            parentPoly = tileAndPoly[1];
          }

          resultRef.push(bestRef);
          resultParent.push(parentRef);
          resultCost.push(bestNode.total);

          for (var _i8 = bestPoly.firstLink; _i8 != NavMesh.DT_NULL_LINK; _i8 = bestTile.links[_i8].next) {
            var link = bestTile.links[_i8];
            var neighbourRef = link.ref; // Skip invalid neighbours and do not follow back to parent.

            if (neighbourRef == 0 || neighbourRef == parentRef) continue; // Expand to neighbour

            tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
            var neighbourTile = tileAndPoly[0];
            var neighbourPoly = tileAndPoly[1]; // Do not advance if the polygon is excluded by the filter.

            if (!filter.passFilter(neighbourRef, neighbourTile, neighbourPoly)) continue; // Find edge and calc distance to the edge.

            var pp = this.getPortalPoints7(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile, 0, 0);
            var va = pp.left;
            var vb = pp.right; // If the let is not touching the edge to the next polygon, skip the connection it.

            var ir = DetourCommon.intersectSegmentPoly2D(va, vb, verts, nverts);
            if (!ir.intersects) continue;
            if (ir.tmin > 1.0 || ir.tmax < 0.0) continue;
            var neighbourNode = this.m_nodePool.getNode(neighbourRef);
            if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0) continue; // Cost

            if (neighbourNode.flags == 0) neighbourNode.pos = DetourCommon.vLerp3(va, vb, 0.5);
            cost = filter.getCost(bestNode.pos, neighbourNode.pos, parentRef, parentTile, parentPoly, bestRef, bestTile, bestPoly, neighbourRef, neighbourTile, neighbourPoly);
            total = bestNode.total + cost; // The node is already in open list and the new result is worse, skip.

            if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0 && total >= neighbourNode.total) continue;
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
      } /// @par
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

    }, {
      key: "findLocalNeighbourhood",
      value: function findLocalNeighbourhood(startRef, centerPos, radius, filter) {
        // Validate input
        if (startRef == 0 || !this.m_nav.isValidPolyRef(startRef)) throw new IllegalArgumentException("Invalid start ref");
        var resultRef = [];
        var resultParent = [];
        this.m_tinyNodePool = new NodePool();
        var startNode = this.m_tinyNodePool.getNode(startRef);
        startNode.pidx = 0;
        startNode.id = startRef;
        startNode.flags = Node.DT_NODE_CLOSED;
        var stack = [];
        stack.push(startNode);
        resultRef.push(startNode.id);
        resultParent.push(0);
        var radiusSqr = DetourCommon.sqr(radius);
        var pa = new Array(this.m_nav.getMaxVertsPerPoly() * 3);
        var pb = new Array(this.m_nav.getMaxVertsPerPoly() * 3);

        while (!stack.length == 0) {
          // Pop front.
          var curNode = stack.pop(); // Get let and tile.
          // The API input has been cheked already, skip checking internal data.

          var curRef = curNode.id;
          var tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(curRef);
          var curTile = tileAndPoly[0];
          var curPoly = tileAndPoly[1];

          for (var i = curPoly.firstLink; i != NavMesh.DT_NULL_LINK; i = curTile.links[i].next) {
            var link = curTile.links[i];
            var neighbourRef = link.ref; // Skip invalid neighbours.

            if (neighbourRef == 0) continue; // Skip if cannot alloca more nodes.

            var neighbourNode = this.m_tinyNodePool.getNode(neighbourRef);
            if (neighbourNode == null) continue; // Skip visited.

            if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0) continue; // Expand to neighbour

            tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
            var neighbourTile = tileAndPoly[0];
            var neighbourPoly = tileAndPoly[1]; // Skip off-mesh connections.

            if (neighbourPoly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) continue; // Do not advance if the polygon is excluded by the filter.

            if (!filter.passFilter(neighbourRef, neighbourTile, neighbourPoly)) continue; // Find edge and calc distance to the edge.

            var pp = this.getPortalPoints7(curRef, curPoly, curTile, neighbourRef, neighbourPoly, neighbourTile, 0, 0);
            var va = pp.left;
            var vb = pp.right; // If the circle is not touching the next polygon, skip it.

            var distseg = DetourCommon.distancePtSegSqr2D3(centerPos, va, vb);
            var _distSqr3 = distseg[0];
            if (_distSqr3 > radiusSqr) continue; // Mark node visited, this is done before the overlap test so that
            // we will not visit the let again if the test fails.

            neighbourNode.flags |= Node.DT_NODE_CLOSED;
            neighbourNode.pidx = this.m_tinyNodePool.getNodeIdx(curNode); // Check that the polygon does not collide with existing polygons.
            // Collect vertices of the neighbour poly.

            var npa = neighbourPoly.vertCount;

            for (var k = 0; k < npa; ++k) {
              arraycopy$4(neighbourTile.data.verts, neighbourPoly.verts[k] * 3, pa, k * 3, 3);
            }

            var overlap = false;

            for (var j = 0; j < resultRef.length; ++j) {
              var pastRef = resultRef[j]; // Connected polys do not overlap.

              var connected = false;

              for (var _k2 = curPoly.firstLink; _k2 != NavMesh.DT_NULL_LINK; _k2 = curTile.links[_k2].next) {
                if (curTile.links[_k2].ref == pastRef) {
                  connected = true;
                  break;
                }
              }

              if (connected) continue; // Potentially overlapping.

              tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(pastRef);
              var pastTile = tileAndPoly[0];
              var pastPoly = tileAndPoly[1]; // Get vertices and test overlap

              var npb = pastPoly.vertCount;

              for (var _k3 = 0; _k3 < npb; ++_k3) {
                arraycopy$4(pastTile.data.verts, pastPoly.verts[_k3] * 3, pb, _k3 * 3, 3);
              }

              if (DetourCommon.overlapPolyPoly2D(pa, npa, pb, npb)) {
                overlap = true;
                break;
              }
            }

            if (overlap) continue;
            resultRef.push(neighbourRef);
            resultParent.push(curRef);
            stack.push(neighbourNode);
          }
        }

        return new FindLocalNeighbourhoodResult(resultRef, resultParent);
      }
    }, {
      key: "insertInterval",
      value: function insertInterval(ints, tmin, tmax, ref) {
        // Find insertion point.
        var idx = 0;

        while (idx < ints.length) {
          if (tmax <= ints[idx].tmin) break;
          idx++;
        } // Store


        ints.add(idx, new SegInterval(ref, tmin, tmax));
      }
    }, {
      key: "or",
      value: function or(v1, v2) {
        var hi = 0x80000000;
        var low = 0x7fffffff;
        var hi1 = ~~(v1 / hi);
        var hi2 = ~~(v2 / hi);
        var low1 = v1 & low;
        var low2 = v2 & low;
        var h = hi1 | hi2;
        var l = low1 | low2;
        return h * hi + l;
      } /// @par
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

    }, {
      key: "getPolyWallSegments",
      value: function getPolyWallSegments(ref, storePortals, filter) {
        var tileAndPoly = this.m_nav.getTileAndPolyByRef(ref);
        var tile = tileAndPoly[0];
        var poly = tileAndPoly[1];
        var segmentRefs = [];
        var segmentVerts = [];
        var ints = new Array(16);

        for (var i = 0, j = poly.vertCount - 1; i < poly.vertCount; j = i++) {
          // Skip non-solid edges.
          ints = [];

          if ((poly.neis[j] & NavMesh.DT_EXT_LINK) != 0) {
            // Tile border.
            for (var k = poly.firstLink; k != NavMesh.DT_NULL_LINK; k = tile.links[k].next) {
              var link = tile.links[k];

              if (link.edge == j) {
                if (link.ref != 0) {
                  tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(link.ref);
                  var neiTile = tileAndPoly[0];
                  var neiPoly = tileAndPoly[1];

                  if (filter.passFilter(link.ref, neiTile, neiPoly)) {
                    insertInterval(ints, link.bmin, link.bmax, link.ref);
                  }
                }
              }
            }
          } else {
            // Internal edge
            var neiRef = 0;

            if (poly.neis[j] != 0) {
              var idx = poly.neis[j] - 1;
              neiRef = this.or(this.m_nav.getPolyRefBase(tile), idx);
              if (!filter.passFilter(neiRef, tile, tile.data.polys[idx])) neiRef = 0;
            } // If the edge leads to another polygon and portals are not stored, skip.


            if (neiRef != 0 && !storePortals) continue;

            var _vj2 = poly.verts[j] * 3;

            var _vi2 = poly.verts[i] * 3;

            var seg = new Array(6);
            arraycopy$4(tile.data.verts, _vj2, seg, 0, 3);
            arraycopy$4(tile.data.verts, _vi2, seg, 3, 3);
            segmentVerts.push(seg);
            segmentRefs.push(neiRef);
            continue;
          } // Add sentinels


          insertInterval(ints, -1, 0, 0);
          insertInterval(ints, 255, 256, 0); // Store segments.

          var vj = poly.verts[j] * 3;
          var vi = poly.verts[i] * 3;

          for (var _k4 = 1; _k4 < ints.length; ++_k4) {
            // Portal segment.
            if (storePortals && ints[_k4].ref != 0) {
              tmin = ints[_k4].tmin / 255.0;
              tmax = ints[_k4].tmax / 255.0;

              var _seg = new Array(6);

              arraycopy$4(vLerp4(tile.data.verts, vj, vi, tmin), 0, _seg, 0, 3);
              arraycopy$4(vLerp4(tile.data.verts, vj, vi, tmax), 0, _seg, 3, 3);
              segmentVerts.push(_seg);
              segmentRefs.push(ints[_k4].ref);
            } // Wall segment.


            var imin = ints[_k4 - 1].tmax;
            var imax = ints[_k4].tmin;

            if (imin != imax) {
              tmin = imin / 255.0;
              tmax = imax / 255.0;

              var _seg2 = new Array(6);

              arraycopy$4(vLerp4(tile.data.verts, vj, vi, tmin), 0, _seg2, 0, 3);
              arraycopy$4(vLerp4(tile.data.verts, vj, vi, tmax), 0, _seg2, 3, 3);
              segmentVerts.push(_seg2);
              segmentRefs.push(0);
            }
          }
        }

        return new GetPolyWallSegmentsResult(segmentVerts, segmentRefs);
      } /// @par
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

    }, {
      key: "findDistanceToWall",
      value: function findDistanceToWall(startRef, centerPos, maxRadius, filter) {
        // Validate input
        if (startRef == 0 || !this.m_nav.isValidPolyRef(startRef)) throw new IllegalArgumentException("Invalid start ref");
        this.m_nodePool = [];
        this.m_openList = [];
        var startNode = this.m_nodePool.getNode(startRef);
        DetourCommon.vCopy(startNode.pos, centerPos);
        startNode.pidx = 0;
        startNode.cost = 0;
        startNode.total = 0;
        startNode.id = startRef;
        startNode.flags = Node.DT_NODE_OPEN;
        this.m_openList.push(startNode);
        radiusSqr = DetourCommon.sqr(maxRadius);
        var hitPos = new Array(3);
        var bestvj = null;
        var bestvi = null;

        while (!this.m_openList.length == 0) {
          var bestNode = this.m_openList.pop();
          bestNode.flags &= ~Node.DT_NODE_OPEN;
          bestNode.flags |= Node.DT_NODE_CLOSED; // Get let and tile.
          // The API input has been cheked already, skip checking internal data.

          var bestRef = bestNode.id;
          var tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(bestRef);
          var bestTile = tileAndPoly[0];
          var bestPoly = tileAndPoly[1]; // Get parent let and tile.

          var parentRef = 0;
          var parentTile = null;
          var parentPoly = null;
          if (bestNode.pidx != 0) parentRef = this.m_nodePool.getNodeAtIdx(bestNode.pidx).id;

          if (parentRef != 0) {
            tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(parentRef);
            parentTile = tileAndPoly[0];
            parentPoly = tileAndPoly[1];
          } // Hit test walls.


          for (var i = 0, j = bestPoly.vertCount - 1; i < bestPoly.vertCount; j = i++) {
            // Skip non-solid edges.
            if ((bestPoly.neis[j] & NavMesh.DT_EXT_LINK) != 0) {
              // Tile border.
              var solid = true;

              for (var k = bestPoly.firstLink; k != NavMesh.DT_NULL_LINK; k = bestTile.links[k].next) {
                var link = bestTile.links[k];

                if (link.edge == j) {
                  if (link.ref != 0) {
                    tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(link.ref);
                    var neiTile = tileAndPoly[0];
                    var neiPoly = tileAndPoly[1];
                    if (filter.passFilter(link.ref, neiTile, neiPoly)) solid = false;
                  }

                  break;
                }
              }

              if (!solid) continue;
            } else if (bestPoly.neis[j] != 0) {
              // Internal edge
              var idx = bestPoly.neis[j] - 1;
              var ref = this.m_nav.getPolyRefBase(bestTile) | idx;
              if (filter.passFilter(ref, bestTile, bestTile.data.polys[idx])) continue;
            } // Calc distance to the edge.


            var vj = bestPoly.verts[j] * 3;
            var vi = bestPoly.verts[i] * 3;
            var distseg = DetourCommon.distancePtSegSqr2D4(centerPos, bestTile.data.verts, vj, vi);
            distSqr = distseg[0];
            tseg = distseg[1]; // Edge is too far, skip.

            if (distSqr > radiusSqr) continue; // Hit wall, update radius.

            radiusSqr = distSqr; // Calculate hit pos.

            hitPos[0] = bestTile.data.verts[vj] + (bestTile.data.verts[vi] - bestTile.data.verts[vj]) * tseg;
            hitPos[1] = bestTile.data.verts[vj + 1] + (bestTile.data.verts[vi + 1] - bestTile.data.verts[vj + 1]) * tseg;
            hitPos[2] = bestTile.data.verts[vj + 2] + (bestTile.data.verts[vi + 2] - bestTile.data.verts[vj + 2]) * tseg;
            bestvj = new VectorPtr(bestTile.data.verts, vj);
            bestvi = new VectorPtr(bestTile.data.verts, vi);
          }

          for (var _i9 = bestPoly.firstLink; _i9 != NavMesh.DT_NULL_LINK; _i9 = bestTile.links[_i9].next) {
            var _link = bestTile.links[_i9];
            var neighbourRef = _link.ref; // Skip invalid neighbours and do not follow back to parent.

            if (neighbourRef == 0 || neighbourRef == parentRef) continue; // Expand to neighbour.

            tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
            var neighbourTile = tileAndPoly[0];
            var neighbourPoly = tileAndPoly[1]; // Skip off-mesh connections.

            if (neighbourPoly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) continue; // Calc distance to the edge.

            var va = bestPoly.verts[_link.edge] * 3;
            var vb = bestPoly.verts[(_link.edge + 1) % bestPoly.vertCount] * 3;

            var _distseg = DetourCommon.distancePtSegSqr2D4(centerPos, bestTile.data.verts, va, vb);

            distSqr = _distseg[0]; // If the circle is not touching the next polygon, skip it.

            if (distSqr > radiusSqr) continue;
            if (!filter.passFilter(neighbourRef, neighbourTile, neighbourPoly)) continue;
            var neighbourNode = this.m_nodePool.getNode(neighbourRef);
            if ((neighbourNode.flags & Node.DT_NODE_CLOSED) != 0) continue; // Cost

            if (neighbourNode.flags == 0) {
              neighbourNode.pos = this.getEdgeMidPoint6(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile);
            }

            total = bestNode.total + DetourCommon.vDist2(bestNode.pos, neighbourNode.pos); // The node is already in open list and the new result is worse, skip.

            if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0 && total >= neighbourNode.total) continue;
            neighbourNode.id = neighbourRef;
            neighbourNode.flags = neighbourNode.flags & ~Node.DT_NODE_CLOSED;
            neighbourNode.pidx = this.m_nodePool.getNodeIdx(bestNode);
            neighbourNode.total = total;

            if ((neighbourNode.flags & Node.DT_NODE_OPEN) != 0) {
              this.m_openList.modify(neighbourNode);
            } else {
              neighbourNode.flags |= Node.DT_NODE_OPEN;
              this.m_openList.push(neighbourNode);
            }
          }
        } // Calc hit normal.


        var hitNormal = new Array(3);

        if (bestvi != null && bestvj != null) {
          var tangent = DetourCommon.vSub(bestvi, bestvj);
          hitNormal[0] = tangent[2];
          hitNormal[1] = 0;
          hitNormal[2] = -tangent[0];
          vNormalize(hitNormal);
        }

        return new FindDistanceToWallResult(Math.sqrt(radiusSqr), hitPos, hitNormal);
      } /// Returns true if the polygon reference is valid and passes the filter restrictions.
      ///  @param[in]		ref			The polygon reference to check.
      ///  @param[in]		filter		The filter to apply.

    }, {
      key: "isValidPolyRef",
      value: function isValidPolyRef(ref, filter) {
        try {
          var tileAndPoly = this.m_nav.getTileAndPolyByRef(ref); // If cannot pass filter, assume flags has changed and boundary is invalid.

          if (filter.passFilter(ref, tileAndPoly[0], tileAndPoly[1])) return true;
        } catch (e) {// If cannot get polygon, assume it does not exists and boundary is invalid.
        }

        return false;
      } /// Gets the navigation mesh the query object is using.
      /// @return The navigation mesh the query object is using.

    }, {
      key: "getAttachedNavMesh",
      value: function getAttachedNavMesh() {
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

    }, {
      key: "getPathFromDijkstraSearch",
      value: function getPathFromDijkstraSearch(endRef) {
        if (!this.m_nav.isValidPolyRef(endRef)) throw new IllegalArgumentException("Invalid end ref");
        var nodes = this.m_nodePool.findNodes(endRef);
        if (nodes.length != 1) throw new IllegalArgumentException("Invalid end ref");
        var endNode = nodes[0];
        if ((endNode.flags & DT_NODE_CLOSED) == 0) throw new IllegalArgumentException("Invalid end ref");
        return getPathToNode(endNode);
      }
      /**
       * Gets the path leading to the specified end node.
       */

    }, {
      key: "getPathToNode",
      value: function getPathToNode(endNode) {
        var path = []; // Reverse the path.

        var curNode = endNode;

        do {
          path.add(0, curNode.id);
          curNode = this.m_nodePool.getNodeAtIdx(curNode.pidx);
        } while (curNode != null);

        return path;
      }
    }]);

    return NavMeshQuery;
  }();

  _defineProperty(NavMeshQuery, "DT_FINDPATH_ANY_ANGLE", 0x02);

  _defineProperty(NavMeshQuery, "DT_RAYCAST_USE_COSTS", 0x01);

  _defineProperty(NavMeshQuery, "DT_STRAIGHTPATH_START", 0x01);

  _defineProperty(NavMeshQuery, "DT_STRAIGHTPATH_END", 0x02);

  _defineProperty(NavMeshQuery, "DT_STRAIGHTPATH_OFFMESH_CONNECTION", 0x04);

  _defineProperty(NavMeshQuery, "DT_STRAIGHTPATH_AREA_CROSSINGS", 0x01);

  _defineProperty(NavMeshQuery, "DT_STRAIGHTPATH_ALL_CROSSINGS", 0x02);

  _defineProperty(NavMeshQuery, "H_SCALE", 0.999);

  _defineProperty(NavMeshQuery, "FRand", function FRand() {
    _classCallCheck(this, FRand);

    return Math.random();
  });

  _defineProperty(NavMeshQuery, "s", 1.0 / 255.0);

  _defineProperty(NavMeshQuery, "SegInternval", (_temp$3 = function SegInterval(ref, tmin, tmax) {
    _classCallCheck(this, SegInterval);

    _defineProperty(this, "ref", void 0);

    _defineProperty(this, "tmin", void 0);

    _defineProperty(this, "tmax", void 0);

    this.ref = ref;
    this.tmin = tmin;
    this.tmax = tmax;
  }, _temp$3));

  var ObstacleAvoidanceParams = ///< grid
  ///< adaptive
  ///< adaptive
  ///< adaptive
  function ObstacleAvoidanceParams() {
    _classCallCheck(this, ObstacleAvoidanceParams);

    _defineProperty(this, "velBias", 0);

    _defineProperty(this, "weightDesVel", 0);

    _defineProperty(this, "weightCurVel", 0);

    _defineProperty(this, "weightSide", 0);

    _defineProperty(this, "weightToi", 0);

    _defineProperty(this, "horizTime", 0);

    _defineProperty(this, "gridSize", 0);

    _defineProperty(this, "adaptiveDivs", 0);

    _defineProperty(this, "adaptiveRings", 0);

    _defineProperty(this, "adaptiveDepth", 0);

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
  };

  var ObstacleCircle = function ObstacleCircle() {
    _classCallCheck(this, ObstacleCircle);

    _defineProperty(this, "p", new Array(3));

    _defineProperty(this, "vel", new Array(3));

    _defineProperty(this, "dvel", new Array(3));

    _defineProperty(this, "rad", void 0);

    _defineProperty(this, "dp", new Array(3));

    _defineProperty(this, "np", new Array(3));
  };

  var ObstacleSegment = function ObstacleSegment() {
    _classCallCheck(this, ObstacleSegment);

    _defineProperty(this, "p", new Array(3));

    _defineProperty(this, "q", new Array(3));

    _defineProperty(this, "touch", false);
  };

  var SweepCircleCircleResult = function SweepCircleCircleResult(intersection, htmin, htmax) {
    _classCallCheck(this, SweepCircleCircleResult);

    _defineProperty(this, "intersection", false);

    _defineProperty(this, "htmin", 0);

    _defineProperty(this, "htmax", 0);

    this.intersection = intersection;
    this.htmin = htmin;
    this.htmax = htmax;
  };

  var ObstacleAvoidanceQuery = /*#__PURE__*/function () {
    ///< Max numver of adaptive divs.
    ///< Max number of adaptive rings.
    function ObstacleAvoidanceQuery(maxCircles, maxSegments) {
      _classCallCheck(this, ObstacleAvoidanceQuery);

      _defineProperty(this, "m_params", new ObstacleAvoidanceParams());

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
        this.m_circles[i] = new ObstacleCircle();
      }

      this.m_maxSegments = maxSegments;
      this.m_nsegments = 0;
      this.m_segments = new Array(this.m_maxSegments);

      for (var _i = 0; _i < this.m_maxSegments; _i++) {
        this.m_segments[_i] = new ObstacleSegment();
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
        DetourCommon.vCopy(cir.p, pos);
        cir.rad = rad;
        DetourCommon.vCopy(cir.vel, vel);
        DetourCommon.vCopy(cir.dvel, dvel);
      }
    }, {
      key: "addSegment",
      value: function addSegment(p, q) {
        if (this.m_nsegments >= this.m_maxSegments) return;
        var seg = this.m_segments[this.m_nsegments++];
        DetourCommon.vCopy(seg.p, p);
        DetourCommon.vCopy(seg.q, q);
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
          DetourCommon.vCopy(cir.dp, DetourCommon.vSub(pb, pa));
          DetourCommon.vNormalize(cir.dp);
          dv = DetourCommon.vSub(cir.dvel, dvel);
          var a = DetourCommon.triArea2D3(orig, cir.dp, dv);

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
          var dt = DetourCommon.distancePtSegSqr2D3(pos, seg.p, seg.q);
          seg.touch = dt[0] < DetourCommon.sqr(r);
        }
      }
    }, {
      key: "sweepCircleCircle",
      value: function sweepCircleCircle(c0, r0, v, c1, r1) {
        var EPS = 0.0001;
        var s = DetourCommon.vSub(c1, c0);
        var r = r0 + r1;
        var c = DetourCommon.vDot2D(s, s) - r * r;
        var a = DetourCommon.vDot2D(v, v);
        if (a < EPS) return new SweepCircleCircleResult(false, 0, 0); // not moving
        // Overlap, calc time to exit.

        var b = DetourCommon.vDot2D(v, s);
        var d = b * b - a * c;
        if (d < 0.0) return new SweepCircleCircleResult(false, 0, 0); // no intersection.

        a = 1.0 / a;
        var rd = Math.sqrt(d);
        return new SweepCircleCircleResult(true, (b - rd) * a, (b + rd) * a);
      }
    }, {
      key: "isectRaySeg",
      value: function isectRaySeg(ap, u, bp, bq) {
        var v = DetourCommon.vSub(bq, bp);
        var w = DetourCommon.vSub(ap, bp);
        var d = DetourCommon.vPerp2D(u, v);
        if (Math.abs(d) < 1e-6) return [false, 0];
        d = 1.0 / d;
        var t = DetourCommon.vPerp2D(v, w) * d;
        if (t < 0 || t > 1) return [false, 0];
        var s = DetourCommon.vPerp2D(u, w) * d;
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
        var vpen = this.m_params.weightDesVel * (DetourCommon.vDist2D(vcand, dvel) * this.m_invVmax);
        var vcpen = this.m_params.weightCurVel * (DetourCommon.vDist2D(vcand, vel) * this.m_invVmax); // find the threshold hit time to bail out based on the early out penalty
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

          var vab = DetourCommon.vScale(vcand, 2);
          vab = DetourCommon.vSub(vab, vel);
          vab = DetourCommon.vSub(vab, cir.vel); // Side

          side += DetourCommon.clamp(Math.min(DetourCommon.vDot2D(cir.dp, vab) * 0.5 + 0.5, DetourCommon.vDot2D(cir.np, vab) * 2), 0.0, 1.0);
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
            var sdir = DetourCommon.vSub(seg.q, seg.p);
            var snorm = new Array(3);
            snorm[0] = -sdir[2];
            snorm[2] = sdir[0]; // If the velocity is pointing towards the segment, no collision.

            if (DetourCommon.vDot2D(snorm, vcand) < 0.0) continue; // Else immediate collision.

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
        DetourCommon.vSet(nvel, 0, 0, 0);
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
            DetourCommon.vSet(vcand, cvx + x * cs - half, 0, cvz + y * cs - half);
            if (DetourCommon.sqr(vcand[0]) + DetourCommon.sqr(vcand[2]) > DetourCommon.sqr(vmax + cs / 2)) continue;
            penalty = this.processSample(vcand, cs, pos, rad, vel, dvel, minPenalty, debug);
            ns++;

            if (penalty < minPenalty) {
              minPenalty = penalty;
              DetourCommon.vCopy(nvel, vcand);
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
        DetourCommon.vSet(nvel, 0, 0, 0);
        if (debug != null) debug.reset(); // Build sampling pattern aligned to desired velocity.

        var pat = new Array((ObstacleAvoidanceQuery.DT_MAX_PATTERN_DIVS * ObstacleAvoidanceQuery.DT_MAX_PATTERN_RINGS + 1) * 2);
        var npat = 0;
        var ndivs = this.m_params.adaptiveDivs;
        var nrings = this.m_params.adaptiveRings;
        var depth = this.m_params.adaptiveDepth;
        var nd = DetourCommon.clamp(ndivs, 1, ObstacleAvoidanceQuery.DT_MAX_PATTERN_DIVS);
        var nr = DetourCommon.clamp(nrings, 1, ObstacleAvoidanceQuery.DT_MAX_PATTERN_RINGS);
        var da = 1.0 / nd * ObstacleAvoidanceQuery.DT_PI * 2;
        var ca = Math.cos(da);
        var sa = Math.sin(da); // desired direction

        var ddir = new Array(6);
        DetourCommon.vCopy(ddir, dvel);
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
        DetourCommon.vSet(res, dvel[0] * this.m_params.velBias, 0, dvel[2] * this.m_params.velBias);
        var ns = 0;

        for (var k = 0; k < depth; ++k) {
          var _minPenalty = Number.MAX_VALUE;
          var bvel = new Array(3);
          DetourCommon.vSet(bvel, 0, 0, 0);

          for (var _i4 = 0; _i4 < npat; ++_i4) {
            var vcand = new Array(3);
            DetourCommon.vSet(vcand, res[0] + pat[_i4 * 2 + 0] * cr, 0, res[2] + pat[_i4 * 2 + 1] * cr);
            if (DetourCommon.sqr(vcand[0]) + DetourCommon.sqr(vcand[2]) > DetourCommon.sqr(vmax + 0.001)) continue;

            var _penalty = this.processSample(vcand, cr / 10, pos, rad, vel, dvel, _minPenalty, debug);

            ns++;

            if (_penalty < _minPenalty) {
              _minPenalty = _penalty;
              DetourCommon.vCopy(bvel, vcand);
            }
          }

          DetourCommon.vCopy(res, bvel);
          cr *= 0.5;
        }

        DetourCommon.vCopy(nvel, res);
        return [ns, nvel];
      }
    }]);

    return ObstacleAvoidanceQuery;
  }();

  _defineProperty(ObstacleAvoidanceQuery, "DT_MAX_PATTERN_DIVS", 32);

  _defineProperty(ObstacleAvoidanceQuery, "DT_MAX_PATTERN_RINGS", 4);

  _defineProperty(ObstacleAvoidanceQuery, "DT_PI", 3.14159265);

  var _temp$4;

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
  var ProximityGrid = /*#__PURE__*/function () {
    // Optimization
    //items = {};
    function ProximityGrid(m_cellSize, m_invCellSize) {
      _classCallCheck(this, ProximityGrid);

      _defineProperty(this, "m_cellSize", void 0);

      _defineProperty(this, "m_invCellSize", void 0);

      _defineProperty(this, "items", new Map());

      _defineProperty(this, "m_bounds", new Array(4));

      this.m_cellSize = m_cellSize;
      this.m_invCellSize = m_invCellSize; //this.items = {};
      //TODO: I think this line can be removed, it's redundant

      this.items = new Map();
    }

    _createClass(ProximityGrid, [{
      key: "clear",
      value: function clear() {
        //this.items = {};
        this.items = new Map();
        this.m_bounds[0] = 0xfff;
        this.m_bounds[1] = 0xfff;
        this.m_bounds[2] = -0xfff;
        this.m_bounds[3] = -0xfff;
      }
    }, {
      key: "addItem",
      value: function addItem(id, minx, miny, maxx, maxy) {
        var iminx = Math.floor(minx * this.m_invCellSize);
        var iminy = Math.floor(miny * this.m_invCellSize);
        var imaxx = Math.floor(maxx * this.m_invCellSize);
        var imaxy = Math.floor(maxy * this.m_invCellSize);
        this.m_bounds[0] = Math.min(this.m_bounds[0], iminx);
        this.m_bounds[1] = Math.min(this.m_bounds[1], iminy);
        this.m_bounds[2] = Math.min(this.m_bounds[2], imaxx);
        this.m_bounds[3] = Math.min(this.m_bounds[3], imaxy);

        for (var _y = iminy; _y <= imaxy; ++_y) {
          for (var _x = iminx; _x <= imaxx; ++_x) {
            var _key = new ProximityGrid.ItemKey(_x, _y); //let ids = this.items[JSON.stringify(key)];


            var _ids = this.items.get(_key);

            if (_ids == null) {
              _ids = []; //this.items[JSON.stringify(key)]= ids;

              this.items.set(_key, _ids);
            }

            _ids.push(id);
          }
        }
      }
    }, {
      key: "queryItems",
      value: function queryItems(minx, miny, maxx, maxy) {
        var iminx = Math.floor(minx * this.m_invCellSize);
        var iminy = Math.floor(miny * this.m_invCellSize);
        var imaxx = Math.floor(maxx * this.m_invCellSize);
        var imaxy = Math.floor(maxy * this.m_invCellSize);
        var result = new Set();

        for (var _y2 = iminy; _y2 <= imaxy; ++_y2) {
          for (var _x2 = iminx; _x2 <= imaxx; ++_x2) {
            var _key2 = new ProximityGrid.ItemKey(_x2, _y2); //let ids = this.items[JSON.stringify(key)];


            var _ids2 = this.items.get(_key2);

            if (_ids2 != null) {
              result.add.apply(result, _toConsumableArray(_ids2));
            }
          }
        }

        return result;
      }
    }, {
      key: "getItemCountAt",
      value: function getItemCountAt(x, y) {
        key = new ItemKey(x, y);
        ids = this.items.get(key); //ids = this.items[key];

        return ids != null ? ids.length : 0;
      }
    }]);

    return ProximityGrid;
  }();

  _defineProperty(ProximityGrid, "ItemKey", (_temp$4 = /*#__PURE__*/function () {
    function ItemKey(x, y) {
      _classCallCheck(this, ItemKey);

      _defineProperty(this, "x", void 0);

      _defineProperty(this, "y", void 0);

      this.x = x;
      this.y = y;
    }

    _createClass(ItemKey, [{
      key: "hashCode",
      value: function hashCode() {
        prime = 31;
        result = 1;
        result = prime * result + x;
        result = prime * result + y;
        return result;
      }
    }, {
      key: "equals",
      value: function equals(obj) {
        if (this == obj) return true;
        if (obj == null) return false;
        if (getClass() != obj.getClass()) return false;
        other = obj;
        if (x != other.x) return false;
        if (y != other.y) return false;
        return true;
      }
    }]);

    return ItemKey;
  }(), _temp$4));

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
  var PathQuery = function PathQuery() {
    _classCallCheck(this, PathQuery);

    _defineProperty(this, "ref", void 0);

    _defineProperty(this, "startPos", new Array(3));

    _defineProperty(this, "endPos", new Array(3));

    _defineProperty(this, "startRef", void 0);

    _defineProperty(this, "endRef", void 0);

    _defineProperty(this, "path", []);

    _defineProperty(this, "status", void 0);

    _defineProperty(this, "keepAlive", void 0);

    _defineProperty(this, "filter", void 0);
  } /// < TODO: This is potentially dangerous!
  ;

  var PathQueue = /*#__PURE__*/function () {
    // in update ticks.
    function PathQueue(maxSearchNodeCount, nav) {
      _classCallCheck(this, PathQueue);

      _defineProperty(this, "m_queue", new Array(PathQueue.MAX_QUEUE));

      _defineProperty(this, "m_nextHandle", 1);

      _defineProperty(this, "m_queueHead", void 0);

      _defineProperty(this, "m_navquery", void 0);

      this.m_navquery = new NavMeshQuery(nav);

      for (var i = 0; i < PathQueue.MAX_QUEUE; ++i) {
        this.m_queue[i] = new PathQuery();
        this.m_queue[i].ref = PathQueue.DT_PATHQ_INVALID;
        this.m_queue[i].path = new Array(256);
      }

      this.m_queueHead = 0;
    }

    _createClass(PathQueue, [{
      key: "update",
      value: function update(maxIters) {
        // Update path request until there is nothing to update
        // or upto maxIters pathfinder iterations has been consumed.
        var iterCount = maxIters;

        for (var i = 0; i < PathQueue.MAX_QUEUE; ++i) {
          var q = this.m_queue[this.m_queueHead % PathQueue.MAX_QUEUE]; // Skip inactive requests.

          if (q.ref == PathQueue.DT_PATHQ_INVALID) {
            this.m_queueHead++;
            continue;
          } // Handle compPolyed request.


          if (q.status != null && (q.status == Status.SUCCSESS || q.status == Status.PARTIAL_RESULT || q.status == Status.FAILURE)) {
            // If the path result has not been read in few frames, free the slot.
            q.keepAlive++;

            if (q.keepAlive > PathQueue.MAX_KEEP_ALIVE) {
              q.ref = PathQueue.DT_PATHQ_INVALID;
              q.status = null;
            }

            this.m_queueHead++;
            continue;
          } // Handle query start.


          if (q.status == null) {
            q.status = this.m_navquery.initSlicedFindPath(q.startRef, q.endRef, q.startPos, q.endPos, q.filter, 0);
          } // Handle query in progress.


          if (q.status == Status.IN_PROGRESS) {
            var iters = 0;
            var res = this.m_navquery.updateSlicedFindPath(iterCount);
            iters = res.getIterations();
            q.status = res.getStatus();
            iterCount -= iters;
          }

          if (q.status == Status.SUCCSESS || q.status == Status.PARTIAL_RESULT) {
            var path = this.m_navquery.finalizeSlicedFindPath();
            q.status = path.getStatus();
            q.path = path.getRefs();
          }

          if (iterCount <= 0) break;
          this.m_queueHead++;
        }
      }
    }, {
      key: "request",
      value: function request(startRef, endRef, startPos, endPos, filter) {
        // Find empty slot
        var slot = -1;

        for (var i = 0; i < PathQueue.MAX_QUEUE; ++i) {
          if (this.m_queue[i].ref == PathQueue.DT_PATHQ_INVALID) {
            slot = i;
            break;
          }
        } // Could not find slot.


        if (slot == -1) return PathQueue.DT_PATHQ_INVALID;
        var ref = this.m_nextHandle++;
        if (this.m_nextHandle == PathQueue.DT_PATHQ_INVALID) this.m_nextHandle++;
        var q = this.m_queue[slot];
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
    }, {
      key: "getRequestStatus",
      value: function getRequestStatus(ref) {
        for (var i = 0; i < PathQueue.MAX_QUEUE; ++i) {
          if (this.m_queue[i].ref == ref) return this.m_queue[i].status;
        }

        return Status.FAILURE;
      }
    }, {
      key: "getPathResult",
      value: function getPathResult(ref) {
        for (var i = 0; i < PathQueue.MAX_QUEUE; ++i) {
          if (this.m_queue[i].ref == ref) {
            var q = this.m_queue[i]; // Free request for reuse.

            q.ref = PathQueue.DT_PATHQ_INVALID;
            q.status = null;
            return new FindPathResult(Status.SUCCSESS, q.path);
          }
        }

        return new FindPathResult(Status.FAILURE, null);
      }
    }]);

    return PathQueue;
  }();

  _defineProperty(PathQueue, "MAX_QUEUE", 8);

  _defineProperty(PathQueue, "DT_PATHQ_INVALID", 0);

  _defineProperty(PathQueue, "MAX_KEEP_ALIVE", 2);

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

  var PathCorridor = /*#__PURE__*/function () {
    _createClass(PathCorridor, [{
      key: "mergeCorridorStartMoved",
      value: function mergeCorridorStartMoved(path, visited) {
        var furthestPath = -1;
        var furthestVisited = -1; // Find furthest common polygon.

        for (var i = path.length - 1; i >= 0; --i) {
          var found = false;

          for (var j = visited.length - 1; j >= 0; --j) {
            if (path[i] == visited[j]) {
              furthestPath = i;
              furthestVisited = j;
              found = true;
            }
          }

          if (found) break;
        } // If no intersection found just return current path.


        if (furthestPath == -1 || furthestVisited == -1) return path; // Concatenate paths.
        // Adjust beginning of the buffer to include the visited.

        var result = []; // Store visited

        for (var _i = visited.length - 1; _i > furthestVisited; --_i) {
          result.push(visited[_i]);
        }

        result.push.apply(result, _toConsumableArray(path.slice(furthestPath, path.length)));
        return result;
      }
    }, {
      key: "mergeCorridorEndMoved",
      value: function mergeCorridorEndMoved(path, visited) {
        var furthestPath = -1;
        var furthestVisited = -1; // Find furthest common polygon.

        for (var i = 0; i < path.length; ++i) {
          var found = false;

          for (var j = visited.length - 1; j >= 0; --j) {
            if (path[i] == visited[j]) {
              furthestPath = i;
              furthestVisited = j;
              found = true;
            }
          }

          if (found) break;
        } // If no intersection found just return current path.


        if (furthestPath == -1 || furthestVisited == -1) return path; // Concatenate paths.

        var result = path.subList(0, furthestPath);
        result.addAll(visited.subList(furthestVisited, visited.length));
        return result;
      }
    }, {
      key: "mergeCorridorStartShortcut",
      value: function mergeCorridorStartShortcut(path, visited) {
        var furthestPath = -1;
        var furthestVisited = -1; // Find furthest common polygon.

        for (var i = path.length - 1; i >= 0; --i) {
          var found = false;

          for (var j = visited.length - 1; j >= 0; --j) {
            if (path[i] == visited[j]) {
              furthestPath = i;
              furthestVisited = j;
              found = true;
            }
          }

          if (found) break;
        } // If no intersection found just return current path.


        if (furthestPath == -1 || furthestVisited <= 0) return path; // Concatenate paths.
        // Adjust beginning of the buffer to include the visited.

        var result = visited.subList(0, furthestVisited);
        result.addAll(path.subList(furthestPath, path.length));
        return result;
      }
      /**
       * Allocates the corridor's path buffer.
       */

    }]);

    function PathCorridor() {
      _classCallCheck(this, PathCorridor);

      _defineProperty(this, "m_pos", new Array(3));

      _defineProperty(this, "m_target", new Array(3));

      _defineProperty(this, "m_path", []);

      this.m_path = [];
    }
    /**
     * Resets the path corridor to the specified position.
     * @param ref The polygon reference containing the position.
     * @param pos The new position in the corridor. [(x, y, z)]
     */


    _createClass(PathCorridor, [{
      key: "reset",
      value: function reset(ref, pos) {
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

    }, {
      key: "findCorners",
      value: function findCorners(maxCorners, navquery, filter) {
        var MIN_TARGET_DIST = DetourCommon.sqr(0.01);
        var path = navquery.findStraightPath(this.m_pos, this.m_target, this.m_path, maxCorners, 0); // Prune points in the beginning of the path which are too close.
        // for (let iter = path.iterator(); iter.hasNext();) {
        // 	let spi = iter.next();
        // 	if ((spi.getFlags() & NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION) != 0
        // 		|| DetourCommon.vDist2D(spi.getPos(), this.m_pos) > MIN_TARGET_DIST) {
        // 		break;
        // 	}
        // 	iter.remove();
        // }
        // TODO

        var newPath = [];
        var i = 0;
        var l = path.length;

        for (; i < path.length; i++) {
          var spi = path[i];

          if ((spi.getFlags() & NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION) != 0 || DetourCommon.vDist2D(spi.getPos(), this.m_pos) > MIN_TARGET_DIST) {
            break;
          }
        }

        var j = i; // for(; j < l; j++);
        // {
        // 	newPath.push(path[j])
        // }

        while (j < l) {
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

    }, {
      key: "optimizePathVisibility",
      value: function optimizePathVisibility(next, pathOptimizationRange, navquery, filter) {
        // Clamp the ray to max distance.
        var dist = DetourCommon.vDist2D(this.m_pos, next); // If too close to the goal, do not try to optimize.

        if (dist < 0.01) return; // Overshoot a little. This helps to optimize open fields in tiled
        // meshes.

        dist = Math.min(dist + 0.01, pathOptimizationRange); // Adjust ray length.

        var delta = DetourCommon.vSub(next, this.m_pos);
        var goal = DetourCommon.vMad(this.m_pos, delta, pathOptimizationRange / dist);
        var rc = navquery.raycast(this.m_path[0], this.m_pos, goal, filter, 0, 0);

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

    }, {
      key: "optimizePathTopology",
      value: function optimizePathTopology(navquery, filter) {
        if (this.m_path.length < 3) return false;
        var MAX_ITER = 32;
        navquery.initSlicedFindPath(this.m_path[0], this.m_path[this.m_path.length - 1], this.m_pos, this.m_target, filter, 0);
        navquery.updateSlicedFindPath(MAX_ITER);
        var fpr = navquery.finalizeSlicedFindPathPartial(this.m_path);

        if (fpr.getStatus().isSuccess() && fpr.getRefs().length > 0) {
          this.m_path = mergeCorridorStartShortcut(this.m_path, fpr.getRefs());
          return true;
        }

        return false;
      }
    }, {
      key: "moveOverOffmeshConnection",
      value: function moveOverOffmeshConnection(offMeshConRef, refs, start, end, navquery) {
        // Advance the path up to and over the off-mesh connection.
        var prevRef = 0;
        var polyRef = this.m_path[0];
        var npos = 0;

        while (npos < this.m_path.length && polyRef != offMeshConRef) {
          prevRef = polyRef;
          polyRef = this.m_path[npos];
          npos++;
        }

        if (npos == this.m_path.length) {
          // Could not find offMeshConRef
          return false;
        } // Prune path


        this.m_path = this.m_path.subList(npos, this.m_path.length);
        refs[0] = prevRef;
        refs[1] = polyRef;
        var nav = navquery.getAttachedNavMesh();
        var startEnd = nav.getOffMeshConnectionPolyEndPoints(refs[0], refs[1]);
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

    }, {
      key: "movePosition",
      value: function movePosition(npos, navquery, filter) {
        // Move aPoly navmesh and update new position.
        var masResult = navquery.moveAlongSurface(this.m_path[0], this.m_pos, npos, filter);
        this.m_path = this.mergeCorridorStartMoved(this.m_path, masResult.getVisited()); // Adjust the position to stay on top of the navmesh.

        DetourCommon.vCopy(this.m_pos, masResult.getResultPos());

        try {
          this.m_pos[1] = navquery.getPolyHeight(this.m_path[0], masResult.getResultPos());
        } catch (e) {// Silently disregard the returned status of DT_FAILURE | DT_INVALID_PARAM to stay
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

    }, {
      key: "moveTargetPosition",
      value: function moveTargetPosition(npos, navquery, filter) {
        // Move aPoly navmesh and update new position.
        var masResult = navquery.moveAlongSurface(this.m_path[this.m_path.length - 1], this.m_target, npos, filter);
        this.m_path = mergeCorridorEndMoved(this.m_path, masResult.getVisited()); // TODO: should we do that?
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

    }, {
      key: "setCorridor",
      value: function setCorridor(target, path) {
        DetourCommon.vCopy(this.m_target, target);
        this.m_path = _toConsumableArray(path);
      }
    }, {
      key: "fixPathStart",
      value: function fixPathStart(safeRef, safePos) {
        DetourCommon.vCopy(this.m_pos, safePos);

        if (this.m_path.length < 3 && this.m_path.length > 0) {
          var p = this.m_path[this.m_path.length - 1];
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
    }, {
      key: "trimInvalidPath",
      value: function trimInvalidPath(safeRef, safePos, navquery, filter) {
        // Keep valid path as far as possible.
        var n = 0;

        while (n < this.m_path.length && navquery.isValidPolyRef(this.m_path[n], filter)) {
          n++;
        }

        if (n == 0) {
          // The first polyref is bad, use current safe values.
          DetourCommon.vCopy(this.m_pos, safePos);
          this.m_path = [];
          this.m_path.push(safeRef);
        } else if (n < this.m_path.length) {
          this.m_path = this.m_path.subList(0, n); // The path is partially usable.
        } // Clamp target pos to last poly


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

    }, {
      key: "isValid",
      value: function isValid(maxLookAhead, navquery, filter) {
        // Check that all polygons still pass query filter.
        var n = Math.min(this.m_path.length, maxLookAhead);

        for (var i = 0; i < n; ++i) {
          if (!navquery.isValidPolyRef(this.m_path[i], filter)) return false;
        }

        return true;
      }
      /**
       * Gets the current position within the corridor. (In the first polygon.)
       * @return The current position within the corridor.
       */

    }, {
      key: "getPos",
      value: function getPos() {
        return this.m_pos;
      }
      /**
       * Gets the current target within the corridor. (In the last polygon.)
       * @return The current target within the corridor.
       */

    }, {
      key: "getTarget",
      value: function getTarget() {
        return this.m_target;
      }
      /**
       * The polygon reference id of the first polygon in the corridor, the polygon containing the position.
       * @return The polygon reference id of the first polygon in the corridor. (Or zero if there is no path.)
       */

    }, {
      key: "getFirstPoly",
      value: function getFirstPoly() {
        return this.m_path.length == 0 ? 0 : this.m_path[0];
      }
      /**
       * The polygon reference id of the last polygon in the corridor, the polygon containing the target.
       * @return The polygon reference id of the last polygon in the corridor. (Or zero if there is no path.)
       */

    }, {
      key: "getLastPoly",
      value: function getLastPoly() {
        return this.m_path.length == 0 ? 0 : this.m_path[this.m_path.length - 1];
      }
      /**
       * The corridor's path. 
       */

    }, {
      key: "getPath",
      value: function getPath() {
        return this.m_path;
      }
      /**
       * The number of polygons in the current corridor path.
       * @return The number of polygons in the current corridor path.
       */

    }, {
      key: "getPathCount",
      value: function getPathCount() {
        return this.m_path.length;
      }
    }]);

    return PathCorridor;
  }();

  var _temp$5;

  function arraycopy$5(one, oneStart, two, twoStart, len) {
    for (var i = 0; i < len; i++) {
      two[twoStart + i] = one[oneStart + i];
    }
  }

  var LocalBoundary = /*#__PURE__*/function () {
    function LocalBoundary() {
      _classCallCheck(this, LocalBoundary);

      _defineProperty(this, "m_center", new Array(3));

      _defineProperty(this, "m_segs", []);

      _defineProperty(this, "m_polys", []);

      this.m_center[0] = this.m_center[1] = this.m_center[2] = Number.MAX_VALUE;
    }

    _createClass(LocalBoundary, [{
      key: "reset",
      value: function reset() {
        this.m_center[0] = this.m_center[1] = this.m_center[2] = Number.MAX_VALUE;
        this.m_polys = [];
        this.m_segs = [];
      }
    }, {
      key: "addSegment",
      value: function addSegment(dist, s) {
        // Insert neighbour based on the distance.
        var seg = new LocalBoundary.Segment();
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
          var i;

          for (i = 0; i < this.m_segs.length; ++i) {
            if (dist <= this.m_segs[i].d) break;
          }

          this.m_segs.splice(i, 0, seg);
        }

        while (this.m_segs.length > LocalBoundary.MAX_LOCAL_SEGS) {
          this.m_segs.splice(this.m_segs.length - 1, 1);
        }
      }
    }, {
      key: "update",
      value: function update(ref, pos, collisionQueryRange, navquery, filter) {
        if (ref == 0) {
          reset();
          return;
        }

        DetourCommon.vCopy(this.m_center, pos); // First query non-overlapping polygons.

        var res = navquery.findLocalNeighbourhood(ref, pos, collisionQueryRange, filter);
        this.m_polys = res.getRefs();
        this.m_segs = []; // Secondly, store all polygon edges.

        for (var j = 0; j < this.m_polys.length; ++j) {
          var gpws = navquery.getPolyWallSegments(this.m_polys[j], false, filter);

          for (var k = 0; k < gpws.getSegmentRefs().length; ++k) {
            var s = gpws.getSegmentVerts()[k]; // Skip too distant segments.

            var distseg = DetourCommon.distancePtSegSqr2D4(pos, s, 0, 3);
            if (distseg[0] > DetourCommon.sqr(collisionQueryRange)) continue;
            this.addSegment(distseg[0], s);
          }
        }
      }
    }, {
      key: "isValid",
      value: function isValid(navquery, filter) {
        if (this.m_polys.length == 0) return false; // Check that all polygons still pass query filter.

        var _iterator = _createForOfIteratorHelper(this.m_polys),
            _step;

        try {
          for (_iterator.s(); !(_step = _iterator.n()).done;) {
            var ref = _step.value;
            if (!navquery.isValidPolyRef(ref, filter)) return false;
          }
        } catch (err) {
          _iterator.e(err);
        } finally {
          _iterator.f();
        }

        return true;
      }
    }, {
      key: "getCenter",
      value: function getCenter() {
        return this.m_center;
      }
    }, {
      key: "getSegment",
      value: function getSegment(j) {
        return this.m_segs[j].s;
      }
    }, {
      key: "getSegmentCount",
      value: function getSegmentCount() {
        return this.m_segs.length;
      }
    }]);

    return LocalBoundary;
  }();

  _defineProperty(LocalBoundary, "MAX_LOCAL_SEGS", 8);

  _defineProperty(LocalBoundary, "Segment", (_temp$5 = function Segment() {
    _classCallCheck(this, Segment);

    _defineProperty(this, "s", new Array(6).fill(0));

    _defineProperty(this, "d", void 0);
  }, _temp$5));

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
  var CrowdAgentAnimation = function CrowdAgentAnimation() {
    _classCallCheck(this, CrowdAgentAnimation);

    _defineProperty(this, "active", void 0);

    _defineProperty(this, "initPos", new Array(3));

    _defineProperty(this, "startPos", new Array(3));

    _defineProperty(this, "endPos", new Array(3));

    _defineProperty(this, "polyRef", void 0);

    _defineProperty(this, "t", void 0);

    _defineProperty(this, "tmax", void 0);
  };

  /// @ingroup crowd

  var CrowdAgent = /*#__PURE__*/function () {
    /// The type of navigation mesh polygon the agent is currently traversing.
    /// @ingroup crowd
    /// True if the agent is active, false if the agent is in an unused slot in the agent pool.
    /// The type of mesh polygon the agent is traversing. (See: #CrowdAgent)
    /// True if the agent has valid path (targetState == DT_CROWDAGENT_TARGET_VALID) and the path does not lead to the requested position, else false.
    /// The path corridor the agent is using.
    /// The local boundary data for the agent.
    /// Time since the agent's path corridor was optimized.
    /// The known neighbors of the agent.
    /// The desired speed.
    ///< The current agent position. [(x, y, z)]
    ///< A temporary value used to accumulate agent displacement during iterative collision resolution. [(x, y, z)]
    ///< The desired velocity of the agent. Based on the current path, calculated from scratch each frame. [(x, y, z)]
    ///< The desired velocity adjusted by obstacle avoidance, calculated from scratch each frame. [(x, y, z)]
    ///< The actual velocity of the agent. The change from nvel => vel is constrained by max acceleration. [(x, y, z)]
    /// The agent's configuration parameters.
    /// The local path corridor corners for the agent.
    ///< State of the movement request.
    ///< Target polyref of the movement request.
    ///< Target position of the movement request (or velocity in case of DT_CROWDAGENT_TARGET_VELOCITY).
    ///< Path finder ref.
    ///< Flag indicating that the current path is being replanned.
    /// <Time since the agent's target was replanned.
    function CrowdAgent(idx) {
      _classCallCheck(this, CrowdAgent);

      _defineProperty(this, "idx", void 0);

      _defineProperty(this, "active", void 0);

      _defineProperty(this, "state", void 0);

      _defineProperty(this, "partial", void 0);

      _defineProperty(this, "corridor", void 0);

      _defineProperty(this, "boundary", void 0);

      _defineProperty(this, "topologyOptTime", void 0);

      _defineProperty(this, "neis", []);

      _defineProperty(this, "desiredSpeed", void 0);

      _defineProperty(this, "npos", new Array(3));

      _defineProperty(this, "disp", new Array(3));

      _defineProperty(this, "dvel", new Array(3));

      _defineProperty(this, "nvel", new Array(3));

      _defineProperty(this, "vel", new Array(3));

      _defineProperty(this, "params", void 0);

      _defineProperty(this, "corners", []);

      _defineProperty(this, "targetState", void 0);

      _defineProperty(this, "targetRef", void 0);

      _defineProperty(this, "targetPos", new Array(3));

      _defineProperty(this, "targetPathqReq", void 0);

      _defineProperty(this, "targetReplan", void 0);

      _defineProperty(this, "targetReplanTime", void 0);

      _defineProperty(this, "animation", void 0);

      this.idx = idx;
      this.corridor = new PathCorridor();
      this.boundary = new LocalBoundary();
      this.animation = new CrowdAgentAnimation();
    }

    _createClass(CrowdAgent, [{
      key: "integrate",
      value: function integrate(dt) {
        // Fake dynamic constraint.
        var maxDelta = this.params.maxAcceleration * dt;
        var dv = DetourCommon.vSub(this.nvel, this.vel);
        var ds = DetourCommon.vLen(dv);
        if (ds > maxDelta) dv = DetourCommon.vScale(dv, maxDelta / ds);
        this.vel = DetourCommon.vAdd(this.vel, dv); // Integrate

        if (DetourCommon.vLen(this.vel) > 0.0001) this.npos = DetourCommon.vMad(this.npos, this.vel, dt);else DetourCommon.vSet(this.vel, 0, 0, 0);
      }
    }, {
      key: "overOffmeshConnection",
      value: function overOffmeshConnection(radius) {
        if (this.corners.length == 0) return false;
        var offMeshConnection = (this.corners[this.corners.length - 1].getFlags() & NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION) != 0 ? true : false;

        if (offMeshConnection) {
          distSq = DetourCommon.vDist2D(this.npos, this.corners[this.corners.length - 1].getPos());
          if (distSq < radius * radius) return true;
        }

        return false;
      }
    }, {
      key: "getDistanceToGoal",
      value: function getDistanceToGoal(range) {
        if (this.corners.length == 0) return range;
        var endOfPath = (this.corners[this.corners.length - 1].getFlags() & NavMeshQuery.DT_STRAIGHTPATH_END) != 0 ? true : false;
        if (endOfPath) return Math.min(DetourCommon.vDist2D(this.npos, this.corners[this.corners.length - 1].getPos()), range);
        return range;
      }
    }, {
      key: "calcSmoothSteerDirection",
      value: function calcSmoothSteerDirection() {
        var dir = new Array(3);

        if (!this.corners.length == 0) {
          var ip0 = 0;
          var ip1 = Math.min(1, this.corners.length - 1);
          var p0 = this.corners[ip0].getPos();
          var p1 = this.corners[ip1].getPos();
          var dir0 = DetourCommon.vSub(p0, this.npos);
          var dir1 = DetourCommon.vSub(p1, this.npos);
          dir0[1] = 0;
          dir1[1] = 0;
          var len0 = DetourCommon.vLen(dir0);
          var len1 = DetourCommon.vLen(dir1);
          if (len1 > 0.001) dir1 = DetourCommon.vScale(dir1, 1.0 / len1);
          dir[0] = dir0[0] - dir1[0] * len0 * 0.5;
          dir[1] = 0;
          dir[2] = dir0[2] - dir1[2] * len0 * 0.5;
          DetourCommon.vNormalize(dir);
        }

        return dir;
      }
    }, {
      key: "calcStraightSteerDirection",
      value: function calcStraightSteerDirection() {
        var dir = new Array(3);

        if (!this.corners.length == 0) {
          dir = DetourCommon.vSub(this.corners[0].getPos(), this.npos);
          dir[1] = 0;
          DetourCommon.vNormalize(dir);
        }

        return dir;
      }
    }, {
      key: "setTarget",
      value: function setTarget(ref, pos) {
        this.targetRef = ref;
        DetourCommon.vCopy(this.targetPos, pos);
        this.targetPathqRef = PathQueue.DT_PATHQ_INVALID;
        if (this.targetRef != 0) this.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_REQUESTING;else this.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_FAILED;
      }
    }, {
      key: "getAgentIndex",
      value: function getAgentIndex() {
        return this.idx;
      }
    }, {
      key: "isActive",
      value: function isActive() {
        return this.active;
      }
    }]);

    return CrowdAgent;
  }();

  _defineProperty(CrowdAgent, "DT_CROWDAGENT_STATE_INVALID", 0);

  _defineProperty(CrowdAgent, "DT_CROWDAGENT_STATE_WALKING", 1);

  _defineProperty(CrowdAgent, "DT_CROWDAGENT_STATE_OFFMESH", 2);

  _defineProperty(CrowdAgent, "DT_CROWDAGENT_TARGET_NONE", 0);

  _defineProperty(CrowdAgent, "DT_CROWDAGENT_TARGET_FAILED", 1);

  _defineProperty(CrowdAgent, "DT_CROWDAGENT_TARGET_VALID", 2);

  _defineProperty(CrowdAgent, "DT_CROWDAGENT_TARGET_REQUESTING", 3);

  _defineProperty(CrowdAgent, "DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE", 4);

  _defineProperty(CrowdAgent, "DT_CROWDAGENT_TARGET_WAITING_FOR_PATH", 5);

  _defineProperty(CrowdAgent, "DT_CROWDAGENT_TARGET_VELOCITY", 6);

  var DefaultQueryFilter
  /*extends QueryFilter*/
  = /*#__PURE__*/function () {
    // DefaultQueryFilter() {
    // 		this.m_includeFlags = 0xfff;
    // 		this.m_excludeFlags = 0;
    // 		for (let i = 0; i < NavMesh.DT_MAX_AREAS; ++i) {
    // 			m_areaCost[i] = 1.0;
    // 		}
    // 	}
    function DefaultQueryFilter() {
      var includeFlags = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : 0xffff;
      var excludeFlags = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
      var areaCost = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : new Array(NavMesh.DT_MAX_AREAS).fill(1, 0, NavMesh.DT_MAX_AREAS);

      _classCallCheck(this, DefaultQueryFilter);

      _defineProperty(this, "m_excludeFlags", 0);

      _defineProperty(this, "m_includeFlags", 0);

      _defineProperty(this, "m_areaCost", new Array(NavMesh.DT_MAX_AREAS));

      this.m_includeFlags = includeFlags;
      this.m_excludeFlags = excludeFlags;

      for (var i = 0; i < Math.min(NavMesh.DT_MAX_AREAS, areaCost.length); ++i) {
        this.m_areaCost[i] = areaCost[i];
      }

      for (var _i = areaCost.length; _i < NavMesh.DT_MAX_AREAS; ++_i) {
        this.m_areaCost[_i] = 1.0;
      }
    }

    _createClass(DefaultQueryFilter, [{
      key: "and",
      value: function and(v1, v2) {
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
    }, {
      key: "or",
      value: function or(v1, v2) {
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
    }, {
      key: "passFilter",
      value: function passFilter(ref, tile, poly) {
        return this.and(poly.flags, this.m_includeFlags) != 0 && this.and(poly.flags, this.m_excludeFlags) == 0;
      }
    }, {
      key: "getCost",
      value: function getCost(pa, pb, prevRef, prevTile, prevPoly, curRef, curTile, curPoly, nextRef, nextTile, nextPoly) {
        return DetourCommon.vDist2(pa, pb) * this.m_areaCost[curPoly.getArea()];
      }
    }]);

    return DefaultQueryFilter;
  }();

  var Crowd = /*#__PURE__*/function () {
    _createClass(Crowd, [{
      key: "tween",
      /// The maximum number of neighbors that a crowd agent can take into account
      /// for steering decisions.
      /// @ingroup crowd
      /// The maximum number of corners a crowd agent will look ahead in the path.
      /// This value is used for sizing the crowd agent corner buffers.
      /// Due to the behavior of the crowd manager, the actual number of useful
      /// corners will be one less than this number.
      /// @ingroup crowd
      /// The maximum number of crowd avoidance configurations supported by the
      /// crowd manager.
      /// @ingroup crowd
      /// @see dtObstacleAvoidanceParams, dtCrowd::setObstacleAvoidanceParams(), dtCrowd::getObstacleAvoidanceParams(),
      /// dtCrowdAgentParams::obstacleAvoidanceType
      /// The maximum number of query filter types supported by the crowd manager.
      /// @ingroup crowd
      /// @see dtQueryFilter, dtCrowd::getFilter() dtCrowd::getEditableFilter(),
      /// dtCrowdAgentParams::queryFilterType
      /// Provides neighbor data for agents managed by the crowd.
      /// @ingroup crowd
      /// @see dtCrowdAgent::neis, dtCrowd
      value: function tween(t, t0, t1) {
        return DetourCommon.clamp((t - t0) / (t1 - t0), 0.0, 1.0);
      }
    }, {
      key: "getNeighbours",
      value: function getNeighbours(pos, height, range, skip, agents, grid) {
        var result = [];
        var ids = grid.queryItems(pos[0] - range, pos[2] - range, pos[0] + range, pos[2] + range);

        var _iterator = _createForOfIteratorHelper(ids),
            _step;

        try {
          for (_iterator.s(); !(_step = _iterator.n()).done;) {
            var id = _step.value;
            var _ag = agents[id];

            if (_ag == skip || !_ag.active) {
              // console.log("Testing");
              continue;
            } // Check for overlap.


            var diff = DetourCommon.vSub(pos, _ag.npos);

            if (Math.abs(diff[1]) >= (height + _ag.params.height) / 2.0) {
              continue;
            }

            diff[1] = 0;

            var _distSqr = DetourCommon.vLenSqr(diff);

            if (_distSqr > DetourCommon.sqr(range)) {
              continue;
            }

            this.addNeighbour(id, _distSqr, result);
          }
        } catch (err) {
          _iterator.e(err);
        } finally {
          _iterator.f();
        }

        return result;
      }
    }, {
      key: "addNeighbour",
      value: function addNeighbour(idx, dist, neis) {
        // Insert neighbour based on the distance.
        var nei = new this.CrowdNeighbor(idx, dist);
        neis.push(nei);
        neis.sort(function (o1, o2) {
          return o1.dist - o2.dist;
        });
      }
    }, {
      key: "addToOptQueue",
      value: function addToOptQueue(newag, agents) {
        // Insert neighbour based on greatest time.
        agents.push(newag);
      } // Insert neighbour based on greatest time.

    }, {
      key: "addToPathQueue",
      value: function addToPathQueue(newag, agents) {
        agents.push(newag);
      } ///
      /// Initializes the crowd.
      /// May be called more than once to purge and re-initialize the crowd.
      /// @param[in] maxAgents The maximum number of agents the crowd can manage. [Limit: >= 1]
      /// @param[in] maxAgentRadius The maximum radius of any agent that will be added to the crowd. [Limit: > 0]
      /// @param[in] nav The navigation mesh to use for planning.
      /// @return True if the initialization succeeded.

    }]);

    function Crowd(maxAgents, maxAgentRadius, nav) {
      var _temp;

      var queryFilterFactory = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : function (i) {
        return new DefaultQueryFilter();
      };

      _classCallCheck(this, Crowd);

      _defineProperty(this, "CrowdNeighbor", (_temp = /// < The index of the neighbor in the crowd.
      /// < The distance between the current agent and the neighbor.
      function CrowdNeighbour(idx, dist) {
        _classCallCheck(this, CrowdNeighbour);

        _defineProperty(this, "idx", void 0);

        _defineProperty(this, "dist", void 0);

        this.idx = idx;
        this.dist = dist;
      }, _temp));

      _defineProperty(this, "m_maxAgents", void 0);

      _defineProperty(this, "m_agents", []);

      _defineProperty(this, "m_activeAgents", []);

      _defineProperty(this, "m_pathq", void 0);

      _defineProperty(this, "m_obstacleQueryParams", new Array(Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS));

      _defineProperty(this, "m_obstacleQuery", void 0);

      _defineProperty(this, "m_grid", void 0);

      _defineProperty(this, "m_ext", new Array(3));

      _defineProperty(this, "m_filters", new Array(Crowd.DT_CROWD_MAX_QUERY_FILTER_TYPE));

      _defineProperty(this, "m_maxAgentRadius", void 0);

      _defineProperty(this, "m_velocitySampleCount", void 0);

      _defineProperty(this, "m_navquery", void 0);

      this.m_maxAgents = maxAgents;
      this.m_maxAgentRadius = maxAgentRadius;
      DetourCommon.vSet(this.m_ext, this.m_maxAgentRadius * 2.0, this.m_maxAgentRadius * 1.5, this.m_maxAgentRadius * 2.0);
      this.m_grid = new ProximityGrid(this.m_maxAgents * 4, this.m_maxAgentRadius * 3);
      this.m_obstacleQuery = new ObstacleAvoidanceQuery(6, 8);

      for (var i = 0; i < Crowd.DT_CROWD_MAX_QUERY_FILTER_TYPE; i++) {
        this.m_filters[i] = queryFilterFactory.apply(i);
      } // Init obstacle query params.


      for (var _i = 0; _i < Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS; ++_i) {
        this.m_obstacleQueryParams[_i] = new ObstacleAvoidanceParams();
      } // Allocate temp buffer for merging paths.


      this.m_pathq = new PathQueue(Crowd.MAX_PATHQUEUE_NODES, nav);
      this.m_agents = new Array(this.m_maxAgents);
      this.m_activeAgents = [];

      for (var _i2 = 0; _i2 < this.m_maxAgents; ++_i2) {
        this.m_agents[_i2] = new CrowdAgent(_i2);
        this.m_agents[_i2].active = false;
      } // The navquery is mostly used for local searches, no need for large
      // node pool.


      this.m_navquery = new NavMeshQuery(nav);
    } /// Sets the shared avoidance configuration for the specified index.
    /// @param[in] idx The index. [Limits: 0 <= value <
    /// #DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS]
    /// @param[in] params The new configuration.


    _createClass(Crowd, [{
      key: "setObstacleAvoidanceParams",
      value: function setObstacleAvoidanceParams(idx, params) {
        if (idx >= 0 && idx < Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS) {
          this.m_obstacleQueryParams[idx] = params;
        }
      } /// Gets the shared avoidance configuration for the specified index.
      /// @param[in] idx The index of the configuration to retreive.
      /// [Limits: 0 <= value < #DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS]
      /// @return The requested configuration.

    }, {
      key: "getObstacleAvoidanceParams",
      value: function getObstacleAvoidanceParams(idx) {
        if (idx >= 0 && idx < Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS) {
          return this.m_obstacleQueryParams[idx];
        }

        return null;
      } /// The maximum number of agents that can be managed by the object.
      /// @return The maximum number of agents.

    }, {
      key: "getAgentCount",
      value: function getAgentCount() {
        return this.m_maxAgents;
      } /// Gets the specified agent from the pool.
      /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
      /// @return The requested agent.
      /// Agents in the pool may not be in use. Check #dtCrowdAgent.active before using the returned object.

    }, {
      key: "getAgent",
      value: function getAgent(idx) {
        return idx < 0 || idx >= this.m_agents.length ? null : this.m_agents[idx];
      } ///
      /// Gets the specified agent from the pool.
      /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
      /// @return The requested agent.
      /// Agents in the pool may not be in use. Check #dtCrowdAgent.active before using the returned object.

    }, {
      key: "getEditableAgent",
      value: function getEditableAgent(idx) {
        return idx < 0 || idx >= this.m_agents.length ? null : this.m_agents[idx];
      } /// Updates the specified agent's configuration.
      /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
      /// @param[in] params The new agent configuration.

    }, {
      key: "updateAgentParameters",
      value: function updateAgentParameters(idx, params) {
        if (idx < 0 || idx >= this.m_maxAgents) {
          return;
        }

        this.m_agents[idx].params = params;
      } /// Adds a new agent to the crowd.
      /// @param[in] pos The requested position of the agent. [(x, y, z)]
      /// @param[in] params The configutation of the agent.
      /// @return The index of the agent in the agent pool. Or -1 if the agent
      /// could not be added.

    }, {
      key: "addAgent",
      value: function addAgent(pos, params) {
        // Find empty slot.
        var idx = -1;

        for (var i = 0; i < this.m_maxAgents; ++i) {
          if (!this.m_agents[i].active) {
            idx = i;
            break;
          }
        }

        if (idx == -1) {
          return -1;
        }

        var ag = this.m_agents[idx];
        this.updateAgentParameters(idx, params); // Find nearest position on navmesh and place the agent there.

        var nearest = this.m_navquery.findNearestPoly(pos, this.m_ext, this.m_filters[ag.params.queryFilterType]);
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
      } /// Removes the agent from the crowd.
      /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
      ///
      /// The agent is deactivated and will no longer be processed. Its
      /// #dt object
      /// is not removed from the pool. It is marked as inactive so that it is
      /// available for reuse.
      /// Removes the agent from the crowd.
      /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]

    }, {
      key: "removeAgent",
      value: function removeAgent(idx) {
        if (idx >= 0 && idx < this.m_maxAgents) {
          this.m_agents[idx].active = false;
        }
      }
    }, {
      key: "requestMoveTargetReplan",
      value: function requestMoveTargetReplan(ag, ref, pos) {
        ag.setTarget(ref, pos);
        ag.targetReplan = true;
        return true;
      } /// Submits a new move request for the specified agent.
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

    }, {
      key: "requestMoveTarget",
      value: function requestMoveTarget(idx, ref, pos) {
        if (idx < 0 || idx >= this.m_maxAgents) {
          return false;
        }

        if (ref == 0) {
          return false;
        }

        var ag = this.m_agents[idx]; // Initialize request.

        ag.setTarget(ref, pos);
        ag.targetReplan = false;
        return true;
      } /// Submits a new move request for the specified agent.
      /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
      /// @param[in] vel The movement velocity. [(x, y, z)]
      /// @return True if the request was successfully submitted.

    }, {
      key: "requestMoveVelocity",
      value: function requestMoveVelocity(idx, vel) {
        if (idx < 0 || idx >= this.m_maxAgents) {
          return false;
        }

        ag = this.m_agents[idx]; // Initialize request.

        ag.targetRef = 0;
        DetourCommon.vCopy(ag.targetPos, vel);
        ag.targetPathqRef = PathQueue.DT_PATHQ_INVALID;
        ag.targetReplan = false;
        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY;
        return true;
      } /// Resets any request for the specified agent.
      /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
      /// @return True if the request was successfully reseted.

    }, {
      key: "resetMoveTarget",
      value: function resetMoveTarget(idx) {
        if (idx < 0 || idx >= this.m_maxAgents) {
          return false;
        }

        ag = this.m_agents[idx]; // Initialize request.

        ag.targetRef = 0;
        DetourCommon.vSet(ag.targetPos, 0, 0, 0);
        DetourCommon.vSet(ag.dvel, 0, 0, 0);
        ag.targetPathqRef = PathQueue.DT_PATHQ_INVALID;
        ag.targetReplan = false;
        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_NONE;
        return true;
      } /// Gets the active agents let the agent pool.
      /// @param[out] agents An array of agent pointers. [(#dt *) * maxAgents]
      /// @param[in] maxAgents The size of the crowd agent array.
      /// @return The number of agents returned in @p agents.

    }, {
      key: "getActiveAgents",
      value: function getActiveAgents() {
        var agents = []; //new Array(this.m_maxAgents);

        for (var i = 0; i < this.m_maxAgents; ++i) {
          if (this.m_agents[i].active) {
            agents.push(this.m_agents[i]);
          }
        }

        return agents;
      }
    }, {
      key: "updateMoveRequest",
      value: function updateMoveRequest() {
        var queue = new PriorityQueue(function (a1, a2) {
          return a2.targetReplanTime - a1.targetReplanTime;
        }); // Fire off new requests.

        for (var i = 0; i < this.m_maxAgents; ++i) {
          // if(i==12)
          //     console.log("Bad agent.")
          var _ag2 = this.m_agents[i];

          if (!_ag2.active) {
            continue;
          }

          if (_ag2.state == CrowdAgent.DT_CROWDAGENT_STATE_INVALID) {
            continue;
          }

          if (_ag2.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE || _ag2.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
            continue;
          }

          if (_ag2.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_REQUESTING) {
            var path = _ag2.corridor.getPath();

            if (path.length == 0) {
              throw new IllegalArgumentException("Empty path");
            } // Quick search towards the goal.


            this.m_navquery.initSlicedFindPath(path[0], _ag2.targetRef, _ag2.npos, _ag2.targetPos, this.m_filters[_ag2.params.queryFilterType], 0);
            this.m_navquery.updateSlicedFindPath(Crowd.MAX_ITER);
            var pathFound = void 0;

            if (_ag2.targetReplan) // && npath > 10)
              {
                // Try to use existing steady path during replan if
                // possible.
                pathFound = this.m_navquery.finalizeSlicedFindPathPartial(path);
              } else {
              // Try to move towards target when goal changes.
              pathFound = this.m_navquery.finalizeSlicedFindPath();
            }

            var reqPath = pathFound.getRefs();
            var reqPos = new Array(3);

            if (!pathFound.getStatus() == Status.FAILURE && reqPath.length > 0) {
              // In progress or succeed.
              if (reqPath[reqPath.length - 1] != _ag2.targetRef) {
                // Partial path, constrain target position inside the
                // last polygon.
                var cr = this.m_navquery.closestPointOnPoly(reqPath[reqPath.length - 1], _ag2.targetPos);
                reqPos = cr.getClosest();
              } else {
                DetourCommon.vCopy(reqPos, _ag2.targetPos);
              }
            } else {
              // Could not find path, start the request from current
              // location.
              DetourCommon.vCopy(reqPos, _ag2.npos);
              reqPath = [];
              reqPath.push(path[0]);
            }

            _ag2.corridor.setCorridor(reqPos, reqPath);

            _ag2.boundary.reset();

            _ag2.partial = false;

            if (reqPath[reqPath.length - 1] == _ag2.targetRef) {
              _ag2.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_VALID;
              _ag2.targetReplanTime = 0.0;
            } else {
              // The path is longer or potentially unreachable, full plan.
              _ag2.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE;
            }
          }

          if (_ag2.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE) {
            this.addToPathQueue(_ag2, queue);
          }
        }

        while (!queue.isEmpty()) {
          var _ag3 = queue.poll();

          _ag3.targetPathqRef = this.m_pathq.request(_ag3.corridor.getLastPoly(), _ag3.targetRef, _ag3.corridor.getTarget(), _ag3.targetPos, this.m_filters[_ag3.params.queryFilterType]);

          if (_ag3.targetPathqRef != PathQueue.DT_PATHQ_INVALID) {
            _ag3.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_WAITING_FOR_PATH;
          }
        } // Update requests.


        this.m_pathq.update(Crowd.MAX_ITERS_PER_UPDATE); // Process path results.

        for (var _i3 = 0; _i3 < this.m_maxAgents; ++_i3) {
          var _ag4 = this.m_agents[_i3];

          if (!_ag4.active) {
            continue;
          }

          if (_ag4.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE || _ag4.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
            continue;
          }

          if (_ag4.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_WAITING_FOR_PATH) {
            // Poll path queue.
            var status = this.m_pathq.getRequestStatus(_ag4.targetPathqRef);

            if (status != null && status == Status.FAILURE) {
              // Path find failed, retry if the target location is still
              // valid.
              _ag4.targetPathqRef = PathQueue.DT_PATHQ_INVALID;

              if (_ag4.targetRef != 0) {
                _ag4.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_REQUESTING;
              } else {
                _ag4.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_FAILED;
              }

              _ag4.targetReplanTime = 0.0;
            } else if (status != null && (status == Status.SUCCSESS || status == Status.PARTIAL_RESULT)) {
              var _path = _ag4.corridor.getPath();

              if (_path.length == 0) {
                throw new IllegalArgumentException("Empty path");
              } // Apply results.


              var targetPos = _ag4.targetPos;
              var valid = true;

              var _pathFound = this.m_pathq.getPathResult(_ag4.targetPathqRef);

              var res = _pathFound.getRefs();

              status = _pathFound.getStatus();

              if (status == Status.FAILURE || res.length == 0) {
                valid = false;
              }

              if (status != null & status == Status.PARTIAL_RESULT) {
                _ag4.partial = true;
              } else {
                _ag4.partial = false;
              } // Merge result and existing path.
              // The agent might have moved whilst the request is
              // being processed, so the path may have changed.
              // We assume that the end of the path is at the same
              // location
              // where the request was issued.
              // The last ref in the old path should be the same as
              // the location where the request was issued..


              if (valid && _path[_path.length - 1] != res[0]) {
                // if (valid && path[path.length - 1].longValue() != res[0].longValue()) {
                valid = false;
              }

              if (valid) {
                // Put the old path infront of the old path.
                if (_path.length > 1) {
                  // path.remove(path.length - 1);
                  _path.splice(_path.length - 1, 1);

                  _path.push.apply(_path, _toConsumableArray(res));

                  res = _path; // Remove trackbacks

                  for (var j = 1; j < res.length - 1; ++j) {
                    if (j - 1 >= 0 && j + 1 < res.length) {
                      if (res[j - 1] == res[j + 1]) {
                        res.splice(j + 1, 1);
                        res.splice(j, 1);
                        j -= 2;
                      }
                    }
                  }
                } // Check for partial path.


                if (res[res.length - 1] != _ag4.targetRef) {
                  // Partial path, constrain target position inside
                  // the last polygon.
                  var _cr = this.m_navquery.closestPointOnPoly(res[res.length - 1], targetPos);

                  targetPos = _cr.getClosest();
                }
              }

              if (valid) {
                // Set current corridor.
                _ag4.corridor.setCorridor(targetPos, res); // Force to update boundary.


                _ag4.boundary.reset();

                _ag4.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_VALID;
              } else {
                // Something went wrong.
                _ag4.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_FAILED;
              }

              _ag4.targetReplanTime = 0.0;
            }
          }
        }
      }
    }, {
      key: "updateTopologyOptimization",
      // seconds
      value: function updateTopologyOptimization(agents, dt) {
        if (!agents.length == 0) {
          return;
        }

        var queue = new PriorityQueue(function (a1, a2) {
          return Float.compare(a2.topologyOptTime, a1.topologyOptTime);
        });

        for (var i = 0; i < agents.length; ++i) {
          ag = agents[i];

          if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
            continue;
          }

          if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
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
    }, {
      key: "checkPathValidity",
      // seconds
      value: function checkPathValidity(agents, dt) {
        for (var i = 0; i < agents.length; ++i) {
          var _ag5 = agents[i];

          if (_ag5.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
            continue;
          }

          _ag5.targetReplanTime += dt;
          var replan = false; // First check that the current location is valid.

          var agentPos = new Array(3);

          var agentRef = _ag5.corridor.getFirstPoly();

          DetourCommon.vCopy(agentPos, _ag5.npos);

          if (!this.m_navquery.isValidPolyRef(agentRef, this.m_filters[_ag5.params.queryFilterType])) {
            // Current location is not valid, try to reposition.
            // TODO: this can snap agents, how to handle that?
            var fnp = this.m_navquery.findNearestPoly(_ag5.npos, this.m_ext, this.m_filters[_ag5.params.queryFilterType]);
            agentRef = fnp.getNearestRef();

            if (fnp.getNearestPos() != null) {
              DetourCommon.vCopy(agentPos, fnp.getNearestPos());
            }

            if (agentRef == 0) {
              // Could not find location in navmesh, set state to invalid.
              _ag5.corridor.reset(0, agentPos);

              _ag5.partial = false;

              _ag5.boundary.reset();

              _ag5.state = CrowdAgent.DT_CROWDAGENT_STATE_INVALID;
              continue;
            } // Make sure the first polygon is valid, but leave other valid
            // polygons in the path so that replanner can adjust the path
            // better.


            _ag5.corridor.fixPathStart(agentRef, agentPos); // ag.corridor.trimInvalidPath(agentRef, agentPos, m_navquery,
            // &m_filter);


            _ag5.boundary.reset();

            DetourCommon.vCopy(_ag5.npos, agentPos);
            replan = true;
          } // If the agent does not have move target or is controlled by
          // velocity, no need to recover the target nor replan.


          if (_ag5.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE || _ag5.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
            continue;
          } // Try to recover move request position.


          if (_ag5.targetState != CrowdAgent.DT_CROWDAGENT_TARGET_NONE && _ag5.targetState != CrowdAgent.DT_CROWDAGENT_TARGET_FAILED) {
            if (!this.m_navquery.isValidPolyRef(_ag5.targetRef, this.m_filters[_ag5.params.queryFilterType])) {
              // Current target is not valid, try to reposition.
              var _fnp = this.m_navquery.findNearestPoly(_ag5.targetPos, this.m_ext, this.m_filters[_ag5.params.queryFilterType]);

              _ag5.targetRef = _fnp.getNearestRef();

              if (_fnp.getNearestPos() != null) {
                DetourCommon.vCopy(_ag5.targetPos, _fnp.getNearestPos());
              }

              replan = true;
            }

            if (_ag5.targetRef == 0) {
              // Failed to reposition target, fail moverequest.
              _ag5.corridor.reset(agentRef, agentPos);

              _ag5.partial = false;
              _ag5.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_NONE;
            }
          } // If nearby corridor is not valid, replan.


          if (!_ag5.corridor.isValid(Crowd.CHECK_LOOKAHEAD, this.m_navquery, this.m_filters[_ag5.params.queryFilterType])) {
            // Fix current path.
            // ag.corridor.trimInvalidPath(agentRef, agentPos, m_navquery,
            // &m_filter);
            // ag.boundary.reset();
            replan = true;
          } // If the end of the path is near and it is not the requested
          // location, replan.


          if (_ag5.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VALID) {
            if (_ag5.targetReplanTime > Crowd.TARGET_REPLAN_DELAY && _ag5.corridor.getPathCount() < Crowd.CHECK_LOOKAHEAD && _ag5.corridor.getLastPoly() != _ag5.targetRef) {
              replan = true;
            }
          } // Try to replan path to goal.


          if (replan) {
            if (_ag5.targetState != CrowdAgent.DT_CROWDAGENT_TARGET_NONE) {
              this.requestMoveTargetReplan(_ag5, _ag5.targetRef, _ag5.targetPos);
            }
          }
        }
      }
    }, {
      key: "update",
      value: function update(dt, debug, frame) {
        // if (frame == 114)
        //     console.log("Bug");
        this.m_velocitySampleCount = 0;
        var debugIdx = debug != null ? debug.idx : -1;
        var agents = this.getActiveAgents(); // Check that all agents still have valid paths.

        this.checkPathValidity(agents, dt); // Update async move request and path finder.

        this.updateMoveRequest(); // Optimize path topology.

        this.updateTopologyOptimization(agents, dt); // Register agents to proximity grid.

        this.m_grid.clear();

        for (var i = 0; i < agents.length; ++i) {
          var _ag6 = agents[i];
          var p = _ag6.npos;
          var r = _ag6.params.radius;
          this.m_grid.addItem(i, p[0] - r, p[2] - r, p[0] + r, p[2] + r);
        } // Get nearby navmesh segments and agents to collide with.


        var _iterator2 = _createForOfIteratorHelper(agents),
            _step2;

        try {
          for (_iterator2.s(); !(_step2 = _iterator2.n()).done;) {
            var _ag13 = _step2.value;

            if (_ag13.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
              continue;
            } // Update the collision boundary after certain distance has been passed or
            // if it has become invalid.


            var updateThr = _ag13.params.collisionQueryRange * 0.25;

            if (DetourCommon.vDist2DSqr(_ag13.npos, _ag13.boundary.getCenter()) > DetourCommon.sqr(updateThr) || !_ag13.boundary.isValid(this.m_navquery, this.m_filters[_ag13.params.queryFilterType])) {
              _ag13.boundary.update(_ag13.corridor.getFirstPoly(), _ag13.npos, _ag13.params.collisionQueryRange, this.m_navquery, this.m_filters[_ag13.params.queryFilterType]);
            } // Query neighbour agents


            _ag13.neis = this.getNeighbours(_ag13.npos, _ag13.params.height, _ag13.params.collisionQueryRange, _ag13, agents, this.m_grid);
          } // Find next corner to steer to.

        } catch (err) {
          _iterator2.e(err);
        } finally {
          _iterator2.f();
        }

        for (var _i4 = 0; _i4 < agents.length; ++_i4) {
          var _ag7 = agents[_i4];

          if (_ag7.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
            continue;
          }

          if (_ag7.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE || _ag7.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
            continue;
          } // Find corners for steering


          _ag7.corners = _ag7.corridor.findCorners(Crowd.DT_CROWDAGENT_MAX_CORNERS, this.m_navquery, this.m_filters[_ag7.params.queryFilterType]); // Check to see if the corner after the next corner is directly visible,
          // and short cut to there.

          if ((_ag7.params.updateFlags & CrowdAgentParams.DT_CROWD_OPTIMIZE_VIS) != 0 && _ag7.corners.length > 0) {
            var target = _ag7.corners[Math.min(1, _ag7.corners.length - 1)].getPos();

            _ag7.corridor.optimizePathVisibility(target, _ag7.params.pathOptimizationRange, this.m_navquery, this.m_filters[_ag7.params.queryFilterType]); // Copy data for debug purposes.


            if (debugIdx == _i4) {
              DetourCommon.vCopy(debug.optStart, _ag7.corridor.getPos());
              DetourCommon.vCopy(debug.optEnd, target);
            }
          } else {
            // Copy data for debug purposes.
            if (debugIdx == _i4) {
              DetourCommon.vSet(debug.optStart, 0, 0, 0);
              DetourCommon.vSet(debug.optEnd, 0, 0, 0);
            }
          }
        } // Trigger off-mesh connections (depends on corners).


        var _iterator3 = _createForOfIteratorHelper(agents),
            _step3;

        try {
          for (_iterator3.s(); !(_step3 = _iterator3.n()).done;) {
            var _ag14 = _step3.value;

            if (_ag14.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
              continue;
            }

            if (_ag14.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE || _ag14.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
              continue;
            } // Check


            var triggerRadius = _ag14.params.radius * 2.25;

            if (_ag14.overOffmeshConnection(triggerRadius)) {
              // Prepare to off-mesh connection.
              anim = _ag14.animation; // Adjust the path over the off-mesh connection.

              var refs = new long[2]();

              if (_ag14.corridor.moveOverOffmeshConnection(_ag14.corners[_ag14.corners.length - 1].getRef(), refs, anim.startPos, anim.endPos, this.m_navquery)) {
                DetourCommon.vCopy(anim.initPos, _ag14.npos);
                anim.polyRef = refs[1];
                anim.active = true;
                anim.t = 0.0;
                anim.tmax = DetourCommon.vDist2D(anim.startPos, anim.endPos) / _ag14.params.maxSpeed * 0.5;
                _ag14.state = CrowdAgent.DT_CROWDAGENT_STATE_OFFMESH;
                _ag14.corners = [];
                _ag14.neis = [];
                continue;
              } else {// Path validity check will ensure that bad/blocked connections will be replanned.
              }
            }
          } // Calculate steering.

        } catch (err) {
          _iterator3.e(err);
        } finally {
          _iterator3.f();
        }

        var _iterator4 = _createForOfIteratorHelper(agents),
            _step4;

        try {
          for (_iterator4.s(); !(_step4 = _iterator4.n()).done;) {
            var _ag15 = _step4.value;

            if (_ag15.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
              continue;
            }

            if (_ag15.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE) {
              continue;
            }

            var dvel = new Array(3);

            if (_ag15.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
              DetourCommon.vCopy(dvel, _ag15.targetPos);
              _ag15.desiredSpeed = DetourCommon.vLen(_ag15.targetPos);
            } else {
              // Calculate steering direction.
              if ((_ag15.params.updateFlags & CrowdAgentParams.DT_CROWD_ANTICIPATE_TURNS) != 0) {
                dvel = _ag15.calcSmoothSteerDirection();
              } else {
                dvel = _ag15.calcStraightSteerDirection();
              } // Calculate speed scale, which tells the agent to slowdown at the end of the path.


              var slowDownRadius = _ag15.params.radius * 2; // TODO: make less hacky.

              var speedScale = _ag15.getDistanceToGoal(slowDownRadius) / slowDownRadius;
              _ag15.desiredSpeed = _ag15.params.maxSpeed;
              dvel = DetourCommon.vScale(dvel, _ag15.desiredSpeed * speedScale);
            } // Separation


            if ((_ag15.params.updateFlags & CrowdAgentParams.DT_CROWD_SEPARATION) != 0) {
              separationDist = _ag15.params.collisionQueryRange;
              invSeparationDist = 1.0 / separationDist;
              separationWeight = _ag15.params.separationWeight;
              w = 0;
              var disp = new Array(3);

              for (var _j3 = 0; _j3 < _ag15.neis.length; ++_j3) {
                nei = agents.get(_ag15.neis[_j3].idx);

                var _diff = DetourCommon.vSub(_ag15.npos, nei.npos);

                _diff[1] = 0;
                distSqr = DetourCommon.vLenSqr(_diff);

                if (distSqr < 0.00001) {
                  continue;
                }

                if (distSqr > DetourCommon.sqr(separationDist)) {
                  continue;
                }

                dist = Math.sqrt(distSqr);
                weight = separationWeight * (1.0 - DetourCommon.sqr(dist * invSeparationDist));
                disp = DetourCommon.vMad(disp, _diff, weight / dist);
                w += 1.0;
              }

              if (w > 0.0001) {
                // Adjust desired velocity.
                dvel = DetourCommon.vMad(dvel, disp, 1.0 / w); // Clamp desired velocity to desired speed.

                speedSqr = DetourCommon.vLenSqr(dvel);
                desiredSqr = DetourCommon.sqr(_ag15.desiredSpeed);

                if (speedSqr > desiredSqr) {
                  dvel = DetourCommon.vScale(dvel, desiredSqr / speedSqr);
                }
              }
            } // Set the desired velocity.


            DetourCommon.vCopy(_ag15.dvel, dvel);
          } // Velocity planning.

        } catch (err) {
          _iterator4.e(err);
        } finally {
          _iterator4.f();
        }

        for (var _i5 = 0; _i5 < agents.length; ++_i5) {
          var _ag8 = agents[_i5];

          if (_ag8.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
            continue;
          }

          if ((_ag8.params.updateFlags & CrowdAgentParams.DT_CROWD_OBSTACLE_AVOIDANCE) != 0) {
            this.m_obstacleQuery.reset(); // Add neighbours as obstacles.

            for (var j = 0; j < _ag8.neis.length; ++j) {
              var _nei = agents[_ag8.neis[j].idx];
              this.m_obstacleQuery.addCircle(_nei.npos, _nei.params.radius, _nei.vel, _nei.dvel);
            } // if (ag.neis.length > 0 && i == 0) {
            //     console.log("Frame " + frame)
            // }
            // Append neighbour segments as obstacles.


            for (var _j = 0; _j < _ag8.boundary.getSegmentCount(); ++_j) {
              var s = _ag8.boundary.getSegment(_j); //let s3 = Arrays.copyOfRange(s, 3, 6);
              //let s3 = Arrays.copyOfRange(s, 3, 6);


              var s3 = s.slice(3, 6);

              if (DetourCommon.triArea2D3(_ag8.npos, s, s3) < 0.0) {
                continue;
              }

              this.m_obstacleQuery.addSegment(s, s3);
            }

            var vod = null;

            if (debugIdx == _i5) {
              vod = debug.vod;
            } // Sample new safe velocity.


            var adaptive = true;
            var ns = 0;
            var params = this.m_obstacleQueryParams[_ag8.params.obstacleAvoidanceType];

            if (adaptive) {
              var nsnvel = this.m_obstacleQuery.sampleVelocityAdaptive(_ag8.npos, _ag8.params.radius, _ag8.desiredSpeed, _ag8.vel, _ag8.dvel, params, vod);
              ns = nsnvel[0];
              _ag8.nvel = nsnvel[1];
            } else {
              var _nsnvel = this.m_obstacleQuery.sampleVelocityGrid(_ag8.npos, _ag8.params.radius, _ag8.desiredSpeed, _ag8.vel, _ag8.dvel, params, vod);

              ns = _nsnvel[0];
              _ag8.nvel = _nsnvel[1];
            }

            this.m_velocitySampleCount += ns;
          } else {
            // If not using velocity planning, new velocity is directly the desired velocity.
            DetourCommon.vCopy(_ag8.nvel, _ag8.dvel);
          }
        } // Integrate.


        for (var _i6 = 0; _i6 < agents.length; ++_i6) {
          var _ag9 = agents[_i6];

          if (_ag9.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
            continue;
          }

          _ag9.integrate(dt);
        } // Handle collisions.


        for (var iter = 0; iter < 4; ++iter) {
          for (var _i7 = 0; _i7 < agents.length; ++_i7) {
            var _ag10 = agents[_i7];

            var idx0 = _ag10.getAgentIndex();

            if (_ag10.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
              continue;
            }

            DetourCommon.vSet(_ag10.disp, 0, 0, 0);
            var _w = 0;

            for (var _j2 = 0; _j2 < _ag10.neis.length; ++_j2) {
              var _nei2 = agents[_ag10.neis[_j2].idx];

              var idx1 = _nei2.getAgentIndex();

              var diff = DetourCommon.vSub(_ag10.npos, _nei2.npos);
              diff[1] = 0;

              var _dist = DetourCommon.vLenSqr(diff);

              if (_dist > DetourCommon.sqr(_ag10.params.radius + _nei2.params.radius)) {
                continue;
              }

              _dist = Math.sqrt(_dist);
              var pen = _ag10.params.radius + _nei2.params.radius - _dist;

              if (_dist < 0.0001) {
                // Agents on top of each other, try to choose diverging separation directions.
                if (idx0 > idx1) {
                  DetourCommon.vSet(diff, -_ag10.dvel[2], 0, _ag10.dvel[0]);
                } else {
                  DetourCommon.vSet(diff, _ag10.dvel[2], 0, -_ag10.dvel[0]);
                }

                pen = 0.01;
              } else {
                pen = 1.0 / _dist * (pen * 0.5) * Crowd.COLLISION_RESOLVE_FACTOR;
              }

              _ag10.disp = DetourCommon.vMad(_ag10.disp, diff, pen);
              _w += 1.0;
            }

            if (_w > 0.0001) {
              var iw = 1.0 / _w;
              _ag10.disp = DetourCommon.vScale(_ag10.disp, iw);
            }
          }

          for (var _i8 = 0; _i8 < agents.length; ++_i8) {
            var _ag11 = agents[_i8];

            if (_ag11.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
              continue;
            }

            _ag11.npos = DetourCommon.vAdd(_ag11.npos, _ag11.disp);
          }
        }

        for (var _i9 = 0; _i9 < agents.length; ++_i9) {
          // if (frame == 492)
          //     console.log("Bad agent")
          var _ag12 = agents[_i9];

          if (_ag12.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
            continue;
          } // Move along navmesh.


          _ag12.corridor.movePosition(_ag12.npos, this.m_navquery, this.m_filters[_ag12.params.queryFilterType]); // Get valid constrained position back.


          DetourCommon.vCopy(_ag12.npos, _ag12.corridor.getPos()); // If not using path, truncate the corridor to just one poly.

          if (_ag12.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE || _ag12.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
            _ag12.corridor.reset(_ag12.corridor.getFirstPoly(), _ag12.npos);

            _ag12.partial = false;
          }
        } // Update agents using off-mesh connection.


        for (var _i10 = 0; _i10 < this.m_maxAgents; ++_i10) {
          var _anim = this.m_agents[_i10].animation;

          if (!_anim.active) {
            continue;
          }

          ag = this.m_agents[_i10];
          _anim.t += dt;

          if (_anim.t > _anim.tmax) {
            // Reset animation
            _anim.active = false; // Prepare agent for walking.

            ag.state = CrowdAgent.DT_CROWDAGENT_STATE_WALKING;
            continue;
          } // Update position


          ta = _anim.tmax * 0.15;
          tb = _anim.tmax;

          if (_anim.t < ta) {
            u = tween(_anim.t, 0.0, ta);
            ag.npos = DetourCommon.vLerp3(_anim.initPos, _anim.startPos, u);
          } else {
            u = tween(_anim.t, ta, tb);
            ag.npos = DetourCommon.vLerp3(_anim.startPos, _anim.endPos, u);
          } // Update velocity.


          DetourCommon.vSet(ag.vel, 0, 0, 0);
          DetourCommon.vSet(ag.dvel, 0, 0, 0);
        }
      }
    }, {
      key: "getQueryExtents",
      value: function getQueryExtents() {
        return this.m_ext;
      }
    }, {
      key: "getFilter",
      value: function getFilter(i) {
        return i >= 0 && i < Crowd.DT_CROWD_MAX_QUERY_FILTER_TYPE ? this.m_filters[i] : null;
      }
    }]);

    return Crowd;
  }();

  _defineProperty(Crowd, "MAX_ITERS_PER_UPDATE", 100);

  _defineProperty(Crowd, "MAX_PATHQUEUE_NODES", 4096);

  _defineProperty(Crowd, "MAX_COMMON_NODES", 512);

  _defineProperty(Crowd, "DT_CROWDAGENT_MAX_NEIGHBOURS", 6);

  _defineProperty(Crowd, "DT_CROWDAGENT_MAX_CORNERS", 4);

  _defineProperty(Crowd, "DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS", 8);

  _defineProperty(Crowd, "DT_CROWD_MAX_QUERY_FILTER_TYPE", 16);

  _defineProperty(Crowd, "MAX_ITER", 20);

  _defineProperty(Crowd, "OPT_TIME_THR", 0.5);

  _defineProperty(Crowd, "CHECK_LOOKAHEAD", 10);

  _defineProperty(Crowd, "TARGET_REPLAN_DELAY", 1.0);

  _defineProperty(Crowd, "COLLISION_RESOLVE_FACTOR", 0.7);

  exports.Crowd = Crowd;
  exports.CrowdAgentParams = CrowdAgentParams;
  exports.NavMesh = NavMesh;
  exports.NavMeshQuery = NavMeshQuery;
  exports.ObstacleAvoidanceParams = ObstacleAvoidanceParams;
  exports.RecastTestMeshBuilder = RecastTestMeshBuilder;

  Object.defineProperty(exports, '__esModule', { value: true });

})));
