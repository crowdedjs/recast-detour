"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _PolyMeshDetail = _interopRequireDefault(require("./PolyMeshDetail.js"));

var _RecastConstants = _interopRequireDefault(require("./RecastConstants.js"));

var _RecastCommon = _interopRequireDefault(require("./RecastCommon.js"));

var _RecastVectors = _interopRequireDefault(require("./RecastVectors.js"));

var _RecastMesh = _interopRequireDefault(require("./RecastMesh.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _slicedToArray(arr, i) { return _arrayWithHoles(arr) || _iterableToArrayLimit(arr, i) || _unsupportedIterableToArray(arr, i) || _nonIterableRest(); }

function _nonIterableRest() { throw new TypeError("Invalid attempt to destructure non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function _iterableToArrayLimit(arr, i) { if (typeof Symbol === "undefined" || !(Symbol.iterator in Object(arr))) return; var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"] != null) _i["return"](); } finally { if (_d) throw _e; } } return _arr; }

function _arrayWithHoles(arr) { if (Array.isArray(arr)) return arr; }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

function arraycopy(one, oneStart, two, twoStart, len) {
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
      mesh = new _PolyMeshDetail["default"]();
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
          _RecastVectors["default"].copy(mesh.verts, mesh.nverts * 3, dm.verts, k * 3);

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

      _RecastVectors["default"].sub(v2, verts, p2, p1);

      _RecastVectors["default"].sub(v3, verts, p3, p1);

      cp = vcross2(v1, v2, v3);

      if (Math.abs(cp) > EPS) {
        v1Sq = RecastMeshDetail.vdot2(v1, v1);
        v2Sq = RecastMeshDetail.vdot2(v2, v2);
        v3Sq = RecastMeshDetail.vdot2(v3, v3);
        c[0] = (v1Sq * (v2[2] - v3[2]) + v2Sq * (v3[2] - v1[2]) + v3Sq * (v1[2] - v2[2])) / (2 * cp);
        c[1] = 0;
        c[2] = (v1Sq * (v3[0] - v2[0]) + v2Sq * (v1[0] - v3[0]) + v3Sq * (v2[0] - v1[0])) / (2 * cp);
        r.set(RecastMeshDetail.vdist2(c, v1));

        _RecastVectors["default"].add(c, c, verts, p1);

        return true;
      }

      _RecastVectors["default"].copy(c, verts, p1);

      r.set(0);
      return false;
    }
  }, {
    key: "distPtTri",
    value: function distPtTri(p, verts, a, b, c) {
      var v0 = new Array(3);
      var v1 = new Array(3);
      var v2 = new Array(3);

      _RecastVectors["default"].subA(v0, verts, c, a);

      _RecastVectors["default"].subA(v1, verts, b, a);

      _RecastVectors["default"].subB(v2, p, verts, a);

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
      var i, j;
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
      ix = _RecastCommon["default"].clamp(ix - hp.xmin, 0, hp.width - 1);
      iz = _RecastCommon["default"].clamp(iz - hp.ymin, 0, hp.height - 1);
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
      ;

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
        var pi = _RecastMesh["default"].prev(i, nhull);

        var ni = _RecastMesh["default"].next(i, nhull);

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

      while (_RecastMesh["default"].next(left, nhull) != right) {
        // Check to see if se should advance left or right.
        var nleft = _RecastMesh["default"].next(left, nhull);

        var nright = _RecastMesh["default"].prev(right, nhull);

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
        _RecastVectors["default"].copy4(verts, i * 3, _in, i * 3);
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
              _RecastVectors["default"].copy(verts, nverts * 3, edge, idx[_k3] * 3);

              hull[nhull++] = nverts;
              nverts++;
            }
          } else {
            for (var _k4 = 1; _k4 < nidx - 1; ++_k4) {
              _RecastVectors["default"].copy(verts, nverts * 3, edge, idx[_k4] * 3);

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

        _RecastVectors["default"].copy3(bmin, _in, 0);

        _RecastVectors["default"].copy3(bmax, _in, 0);

        for (var _i7 = 1; _i7 < nin; ++_i7) {
          _RecastVectors["default"].min(bmin, _in, _i7 * 3);

          _RecastVectors["default"].max(bmax, _in, _i7 * 3);
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

          _RecastVectors["default"].copy(verts, nverts * 3, bestpt, 0);

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
        if (cx == pcx) directDir = _RecastCommon["default"].rcGetDirForOffset(0, pcy > cy ? 1 : -1);else directDir = _RecastCommon["default"].rcGetDirForOffset(pcx > cx ? 1 : -1, 0); // Push the direct dir last so we start with this on next iteration

        var tmp = dirs[3];
        dirs[3] = dirs[directDir];
        dirs[directDir] = tmp;
        var _cs = chf.spans[ci];

        for (var _i9 = 0; _i9 < 4; ++_i9) {
          var dir = dirs[_i9];
          if (_RecastCommon["default"].GetCon(_cs, dir) == _RecastConstants["default"].RC_NOT_CONNECTED) continue;

          var newX = cx + _RecastCommon["default"].GetDirOffsetX(dir);

          var newY = cy + _RecastCommon["default"].GetDirOffsetY(dir);

          var hpx = newX - hp.xmin;
          var hpy = newY - hp.ymin;
          if (hpx < 0 || hpx >= hp.width || hpy < 0 || hpy >= hp.height) continue;
          if (hp.data[hpx + hpy * hp.width] != 0) continue;
          hp.data[hpx + hpy * hp.width] = 1;
          array.push(newX);
          array.push(newY);
          array.push(chf.cells[newX + bs + (newY + bs) * chf.width].index + _RecastCommon["default"].GetCon(_cs, dir));
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
      hp.data.fill(_RecastConstants["default"].RC_UNSET_HEIGHT, 0, hp.width * hp.height);
      var empty = true; // We cannot sample from this let if it was created from polys
      // of different regions. If it was then it could potentially be overlapping
      // with polys of that region and the heights sampled here could be wrong.

      if (region != _RecastConstants["default"].RC_MULTIPLE_REGS) {
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
                  if (_RecastCommon["default"].GetCon(s, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                    var ax = x + _RecastCommon["default"].GetDirOffsetX(dir);

                    var ay = y + _RecastCommon["default"].GetDirOffsetY(dir);

                    var ai = chf.cells[ax + ay * chf.width].index + _RecastCommon["default"].GetCon(s, dir);

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
          if (_RecastCommon["default"].GetCon(cs, _dir) == _RecastConstants["default"].RC_NOT_CONNECTED) continue;

          var _ax = cx + _RecastCommon["default"].GetDirOffsetX(_dir);

          var _ay = cy + _RecastCommon["default"].GetDirOffsetY(_dir);

          var _hx = _ax - hp.xmin - bs;

          var _hy = _ay - hp.ymin - bs;

          if (_hx < 0 || _hx >= hp.width || _hy < 0 || _hy >= hp.height) continue;
          if (hp.data[_hx + _hy * hp.width] != RecastMeshDetail.RC_UNSET_HEIGHT) continue;

          var _ai = chf.cells[_ax + _ay * chf.width].index + _RecastCommon["default"].GetCon(cs, _dir);

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
      var dmesh = new _PolyMeshDetail["default"]();
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
          if (mesh.polys[p + j] == _RecastConstants["default"].RC_MESH_NULL_IDX) {
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
          if (mesh.polys[_p + _j3] == _RecastConstants["default"].RC_MESH_NULL_IDX) break;

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
          if (dmesh.nverts != 0) arraycopy(dmesh.verts, 0, newv, 0, 3 * dmesh.nverts);
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
          if (dmesh.ntris != 0) arraycopy(dmesh.tris, 0, newt, 0, 4 * dmesh.ntris);
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

var _default = RecastMeshDetail;
exports["default"] = _default;