"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _PolyMesh = _interopRequireDefault(require("./PolyMesh.js"));

var _RecastVectors = _interopRequireDefault(require("./RecastVectors.js"));

var _RecastConstants = _interopRequireDefault(require("./RecastConstants.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

function arraycopy(one, oneStart, two, twoStart, len) {
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
        firstEdge[i] = _RecastConstants["default"].RC_MESH_NULL_IDX;
      }

      for (var _i = 0; _i < npolys; ++_i) {
        var t = _i * vertsPerPoly * 2;

        for (var j = 0; j < vertsPerPoly; ++j) {
          if (polys[t + j] == _RecastConstants["default"].RC_MESH_NULL_IDX) break;
          var v0 = polys[t + j];
          var v1 = j + 1 >= vertsPerPoly || polys[t + j + 1] == _RecastConstants["default"].RC_MESH_NULL_IDX ? polys[t + 0] : polys[t + j + 1];

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
          if (polys[_t + _j] == _RecastConstants["default"].RC_MESH_NULL_IDX) break;
          var _v = polys[_t + _j];

          var _v2 = _j + 1 >= vertsPerPoly || polys[_t + _j + 1] == _RecastConstants["default"].RC_MESH_NULL_IDX ? polys[_t + 0] : polys[_t + _j + 1];

          if (_v > _v2) {
            for (var e = firstEdge[_v2]; e != _RecastConstants["default"].RC_MESH_NULL_IDX; e = firstEdge[nextEdge + e]) {
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
        if (p[i + j] == _RecastConstants["default"].RC_MESH_NULL_IDX) return i;
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

      polys.fill(_RecastConstants["default"].RC_MESH_NULL_IDX, tmp, tmp + nvp);
      var n = 0; // Add pa

      for (var i = 0; i < na - 1; ++i) {
        polys[tmp + n] = polys[pa + (ea + 1 + i) % na];
        n++;
      } // Add pb


      for (var _i12 = 0; _i12 < nb - 1; ++_i12) {
        polys[tmp + n] = polys[pb + (eb + 1 + _i12) % nb];
        n++;
      } //arraycopy(polys, tmp, polys, pa, nvp);


      arraycopy(polys, tmp, polys, pa, nvp); // for (let i = 0; i < nvp; i++) {
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
            arraycopy(mesh.polys, p2, mesh.polys, _p4, nvp);
          }

          mesh.polys.fill(_RecastConstants["default"].RC_MESH_NULL_IDX, _p4 + nvp, _p4 + nvp + nvp);
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
      var tpmPoly = ntris * nvp; // Build initial polygons.

      var npolys = 0;
      polys.fill(_RecastConstants["default"].RC_MESH_NULL_IDX, 0, ntris * nvp);

      for (var _j6 = 0; _j6 < ntris; ++_j6) {
        var t = _j6 * 3;

        if (tris[t + 0] != tris[t + 1] && tris[t + 0] != tris[t + 2] && tris[t + 1] != tris[t + 2]) {
          polys[npolys * nvp + 0] = hole[tris[t + 0]];
          polys[npolys * nvp + 1] = hole[tris[t + 1]];
          polys[npolys * nvp + 2] = hole[tris[t + 2]]; // If this polygon covers multiple region types then
          // mark it as such

          if (hreg[tris[t + 0]] != hreg[tris[t + 1]] || hreg[tris[t + 1]] != hreg[tris[t + 2]]) pregs[npolys] = _RecastConstants["default"].RC_MULTIPLE_REGS;else pregs[npolys] = hreg[tris[t + 0]];
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
            if (pregs[bestPa] != pregs[bestPb]) pregs[bestPa] = _RecastConstants["default"].RC_MULTIPLE_REGS;
            var last = (npolys - 1) * nvp;

            if (pb != last) {
              arraycopy(polys, last, polys, pb, nvp);
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

        mesh.polys.fill(_RecastConstants["default"].RC_MESH_NULL_IDX, _p6, _p6 + nvp * 2);

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
      var mesh = new _PolyMesh["default"]();

      _RecastVectors["default"].copy3(mesh.bmin, cset.bmin, 0);

      _RecastVectors["default"].copy3(mesh.bmax, cset.bmax, 0);

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
      mesh.polys = new Array(maxTris * nvp * 2).fill(_RecastConstants["default"].RC_MESH_NULL_IDX); // Arrays.fill(mesh.polys, RecastConstants.RC_MESH_NULL_IDX);

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

          if ((cont.verts[v + 3] & _RecastConstants["default"].RC_BORDER_VERTEX) != 0) {
            // This vertex should be removed.
            vflags[indices[_j9]] = 1;
          }
        } // Build initial polygons.


        var npolys = 0; // Arrays.fill(polys, RecastConstants.RC_MESH_NULL_IDX);

        polys.fill(_RecastConstants["default"].RC_MESH_NULL_IDX);

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
            if (mesh.polys[_p7 + _j14] == _RecastConstants["default"].RC_MESH_NULL_IDX) break; // Skip connected edges.

            if (mesh.polys[_p7 + nvp + _j14] != _RecastConstants["default"].RC_MESH_NULL_IDX) continue;
            var nj = _j14 + 1;
            if (nj >= nvp || mesh.polys[_p7 + nj] == _RecastConstants["default"].RC_MESH_NULL_IDX) nj = 0;
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
      var mesh = new _PolyMesh["default"]();
      mesh.nvp = meshes[0].nvp;
      mesh.cs = meshes[0].cs;
      mesh.ch = meshes[0].ch;

      _RecastVectors["default"].copy(mesh.bmin, meshes[0].bmin, 0);

      _RecastVectors["default"].copy(mesh.bmax, meshes[0].bmax, 0);

      var maxVerts = 0;
      var maxPolys = 0;
      var maxVertsPerMesh = 0;

      for (var i = 0; i < nmeshes; ++i) {
        _RecastVectors["default"].min(mesh.bmin, meshes[i].bmin, 0);

        _RecastVectors["default"].max(mesh.bmax, meshes[i].bmax, 0);

        maxVertsPerMesh = Math.max(maxVertsPerMesh, meshes[i].nverts);
        maxVerts += meshes[i].nverts;
        maxPolys += meshes[i].npolys;
      }

      mesh.nverts = 0;
      mesh.verts = new Array(maxVerts * 3);
      mesh.npolys = 0;
      mesh.polys = new Array(maxPolys * 2 * mesh.nvp);
      mesh.polys.fill(_RecastConstants["default"].RC_MESH_NULL_IDX, 0, mesh.polys.length);
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
            if (pmesh.polys[src + k] == _RecastConstants["default"].RC_MESH_NULL_IDX) break;
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
      var dst = new _PolyMesh["default"]();
      dst.nverts = src.nverts;
      dst.npolys = src.npolys;
      dst.maxpolys = src.npolys;
      dst.nvp = src.nvp;

      _RecastVectors["default"].copy(dst.bmin, src.bmin, 0);

      _RecastVectors["default"].copy(dst.bmax, src.bmax, 0);

      dst.cs = src.cs;
      dst.ch = src.ch;
      dst.borderSize = src.borderSize;
      dst.maxEdgeError = src.maxEdgeError;
      dst.verts = new Array(src.nverts * 3);
      arraycopy(src.verts, 0, dst.verts, 0, dst.verts.length);
      dst.polys = new Array(src.npolys * 2 * src.nvp);
      arraycopy(src.polys, 0, dst.polys, 0, dst.polys.length);
      dst.regs = new Array(src.npolys);
      arraycopy(src.regs, 0, dst.regs, 0, dst.regs.length);
      dst.areas = new Array(src.npolys);
      arraycopy(src.areas, 0, dst.areas, 0, dst.areas.length);
      dst.flags = new Array(src.npolys);
      arraycopy(src.flags, 0, dst.flags, 0, dst.flags.length);
      return dst;
    }
  }]);

  return RecastMesh;
}();

_defineProperty(RecastMesh, "VERTEX_BUCKET_COUNT", 1 << 12);

var _default = RecastMesh;
exports["default"] = _default;