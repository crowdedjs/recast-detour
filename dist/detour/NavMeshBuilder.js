"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _DetourCommon = _interopRequireDefault(require("./DetourCommon.js"));

var _VectorPtr = _interopRequireDefault(require("./VectorPtr.js"));

var _MeshHeader = _interopRequireDefault(require("./MeshHeader.js"));

var _Poly = _interopRequireDefault(require("./Poly.js"));

var _PolyDetail = _interopRequireDefault(require("./PolyDetail.js"));

var _BVNode = _interopRequireDefault(require("./BVNode.js"));

var _MeshData = _interopRequireDefault(require("./MeshData.js"));

var _OffMeshConnection = _interopRequireDefault(require("./OffMeshConnection.js"));

var _NavMesh = _interopRequireDefault(require("./NavMesh.js"));

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

var BVItem = function BVItem() {
  _classCallCheck(this, BVItem);

  _defineProperty(this, "bmin", new Array(3).fill(0));

  _defineProperty(this, "bmax", new Array(3).fill(0));

  _defineProperty(this, "i", 0);
};

;

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

var CompareItemX = /*#__PURE__*/function () {
  function CompareItemX() {
    _classCallCheck(this, CompareItemX);
  }

  _createClass(CompareItemX, [{
    key: "compare",
    value: function compare(a, b) {
      if (a.bmin[0] < b.bmin[0]) return -1;
      if (a.bmin[0] > b.bmin[0]) return 1;
      return 0;
    }
  }]);

  return CompareItemX;
}();

var CompareItemZ = /*#__PURE__*/function () {
  function CompareItemZ() {
    _classCallCheck(this, CompareItemZ);
  }

  _createClass(CompareItemZ, [{
    key: "compare",
    value: function compare(a, b) {
      if (a.bmin[2] < b.bmin[2]) return -1;
      if (a.bmin[2] > b.bmin[2]) return 1;
      return 0;
    }
  }]);

  return CompareItemZ;
}();

var CompareItemY = /*#__PURE__*/function () {
  function CompareItemY() {
    _classCallCheck(this, CompareItemY);
  }

  _createClass(CompareItemY, [{
    key: "compare",
    value: function compare(a, b) {
      if (a.bmin[1] < b.bmin[1]) return -1;
      if (a.bmin[1] > b.bmin[1]) return 1;
      return 0;
    }
  }]);

  return CompareItemY;
}();

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
      var node = new _BVNode["default"]();
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

          _DetourCommon["default"].vCopy(bmin, params.detailVerts, dv);

          _DetourCommon["default"].vCopy(bmax, params.detailVerts, dv);

          for (var j = 1; j < ndv; j++) {
            _DetourCommon["default"].vMin(bmin, params.detailVerts, dv + j * 3);

            _DetourCommon["default"].vMax(bmax, params.detailVerts, dv + j * 3);
          } // BV-tree uses cs for all dimensions


          it.bmin[0] = _DetourCommon["default"].clamp(Math.round((bmin[0] - params.bmin[0]) * quantFactor), 0, 0xffff);
          it.bmin[1] = _DetourCommon["default"].clamp(Math.round((bmin[1] - params.bmin[1]) * quantFactor), 0, 0xffff);
          it.bmin[2] = _DetourCommon["default"].clamp(Math.round((bmin[2] - params.bmin[2]) * quantFactor), 0, 0xffff);
          it.bmax[0] = _DetourCommon["default"].clamp(Math.round((bmax[0] - params.bmin[0]) * quantFactor), 0, 0xffff);
          it.bmax[1] = _DetourCommon["default"].clamp(Math.round((bmax[1] - params.bmin[1]) * quantFactor), 0, 0xffff);
          it.bmax[2] = _DetourCommon["default"].clamp(Math.round((bmax[2] - params.bmin[2]) * quantFactor), 0, 0xffff);
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

        _DetourCommon["default"].vCopy(bmin, params.bmin);

        _DetourCommon["default"].vCopy(bmax, params.bmax);

        bmin[1] = hmin;
        bmax[1] = hmax;

        for (var _i4 = 0; _i4 < params.offMeshConCount; ++_i4) {
          var p0 = new _VectorPtr["default"](params.offMeshConVerts, (_i4 * 2 + 0) * 3);
          var p1 = new _VectorPtr["default"](params.offMeshConVerts, (_i4 * 2 + 1) * 3);
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
      var header = new _MeshHeader["default"]();
      var navVerts = new Array(3 * totVertCount);
      var navPolys = new Array(totPolyCount);
      var navDMeshes = new Array(params.polyCount);
      var navDVerts = new Array(3 * uniqueDetailVertCount);
      var navDTris = new Array(4 * detailTriCount);
      var navBvtree = new Array(bvTreeSize);
      var offMeshCons = new Array(storedOffMeshConCount); // Store header

      header.magic = _MeshHeader["default"].DT_NAVMESH_MAGIC;
      header.version = _MeshHeader["default"].DT_NAVMESH_VERSION;
      header.x = params.tileX;
      header.y = params.tileY;
      header.layer = params.tileLayer;
      header.userId = params.userId;
      header.polyCount = totPolyCount;
      header.vertCount = totVertCount;
      header.maxLinkCount = maxLinkCount;

      _DetourCommon["default"].vCopy(header.bmin, params.bmin);

      _DetourCommon["default"].vCopy(header.bmax, params.bmax);

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

          arraycopy(params.offMeshConVerts, linkv, navVerts, _v, 6);
          n++;
        }
      } // Store polygons
      // Mesh polys


      var src = 0;

      for (var _i10 = 0; _i10 < params.polyCount; ++_i10) {
        var _p3 = new _Poly["default"](_i10, nvp);

        navPolys[_i10] = _p3;
        _p3.vertCount = 0;
        _p3.flags = params.polyFlags[_i10];

        _p3.setArea(params.polyAreas[_i10]);

        _p3.setType(_Poly["default"].DT_POLYTYPE_GROUND);

        for (var _j4 = 0; _j4 < nvp; ++_j4) {
          if (params.polys[src + _j4] == NavMeshBuilder.MESH_NULL_IDX) break;
          _p3.verts[_j4] = params.polys[src + _j4];

          if ((params.polys[src + nvp + _j4] & 0x8000) != 0) {
            // Border or portal edge.
            var _dir = params.polys[src + nvp + _j4] & 0xf;

            if (_dir == 0xf) // Border
              _p3.neis[_j4] = 0;else if (_dir == 0) // Portal x-
              _p3.neis[_j4] = _NavMesh["default"].DT_EXT_LINK | 4;else if (_dir == 1) // Portal z+
              _p3.neis[_j4] = _NavMesh["default"].DT_EXT_LINK | 2;else if (_dir == 2) // Portal x+
              _p3.neis[_j4] = _NavMesh["default"].DT_EXT_LINK | 0;else if (_dir == 3) // Portal z-
              _p3.neis[_j4] = _NavMesh["default"].DT_EXT_LINK | 6;
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
          var _p4 = new _Poly["default"](offMeshPolyBase + n, nvp);

          navPolys[offMeshPolyBase + n] = _p4;
          _p4.vertCount = 2;
          _p4.verts[0] = offMeshVertsBase + n * 2;
          _p4.verts[1] = offMeshVertsBase + n * 2 + 1;
          _p4.flags = params.offMeshConFlags[_i11];

          _p4.setArea(params.offMeshConAreas[_i11]);

          _p4.setType(_Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION);

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
          var dtl = new _PolyDetail["default"]();
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
            arraycopy(params.detailVerts, (vb + _nv2) * 3, navDVerts, vbase * 3, 3 * (_ndv - _nv2));
            vbase += _ndv - _nv2;
          }
        } // Store triangles.


        arraycopy(params.detailTris, 0, navDTris, 0, 4 * params.detailTriCount);
      } else {
        // Create dummy detail mesh by triangulating polys.
        var tbase = 0;

        for (var _i13 = 0; _i13 < params.polyCount; ++_i13) {
          var _dtl = new _PolyDetail["default"]();

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
          var con = new _OffMeshConnection["default"]();
          offMeshCons[n] = con;
          con.poly = offMeshPolyBase + n; // Copy connection end-points.

          var endPts = _i14 * 2 * 3;
          arraycopy(params.offMeshConVerts, endPts, con.pos, 0, 6);
          con.rad = params.offMeshConRad[_i14];
          con.flags = params.offMeshConDir[_i14] != 0 ? _NavMesh["default"].DT_OFFMESH_CON_BIDIR : 0;
          con.side = offMeshConClass[_i14 * 2 + 1];
          if (params.offMeshConUserID != null) con.userId = params.offMeshConUserID[_i14];
          n++;
        }
      }

      var nmd = new _MeshData["default"]();
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

var _default = NavMeshBuilder;
exports["default"] = _default;