"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _NodePool = _interopRequireDefault(require("./NodePool.js"));

var _NodeQueue = _interopRequireDefault(require("./NodeQueue.js"));

var _DetourCommon = _interopRequireDefault(require("./DetourCommon.js"));

var _FindNearestPolyResult = _interopRequireDefault(require("./FindNearestPolyResult.js"));

var _NavMesh = _interopRequireDefault(require("./NavMesh.js"));

var _Poly = _interopRequireDefault(require("./Poly.js"));

var _ClosestPointOnPolyResult = _interopRequireDefault(require("./ClosestPointOnPolyResult.js"));

var _QueryData = _interopRequireDefault(require("./QueryData.js"));

var _Status = _interopRequireDefault(require("./Status.js"));

var _Node = _interopRequireDefault(require("./Node.js"));

var _UpdateSlicedPathResult = _interopRequireDefault(require("./UpdateSlicedPathResult.js"));

var _FindPathResult = _interopRequireDefault(require("./FindPathResult.js"));

var _FindLocalNeighbourhoodResult = _interopRequireDefault(require("./FindLocalNeighbourhoodResult.js"));

var _GetPolyWallSegmentsResult = _interopRequireDefault(require("./GetPolyWallSegmentsResult.js"));

var _StraightPathItem = _interopRequireDefault(require("./StraightPathItem.js"));

var _MoveAlongSurfaceResult = _interopRequireDefault(require("./MoveAlongSurfaceResult.js"));

var _RaycastHit = _interopRequireDefault(require("./RaycastHit.js"));

var _IntersectResult = _interopRequireDefault(require("./IntersectResult.js"));

var _temp;

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _toConsumableArray(arr) { return _arrayWithoutHoles(arr) || _iterableToArray(arr) || _unsupportedIterableToArray(arr) || _nonIterableSpread(); }

function _nonIterableSpread() { throw new TypeError("Invalid attempt to spread non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _iterableToArray(iter) { if (typeof Symbol !== "undefined" && Symbol.iterator in Object(iter)) return Array.from(iter); }

function _arrayWithoutHoles(arr) { if (Array.isArray(arr)) return _arrayLikeToArray(arr); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

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

function arraycopy(one, oneStart, two, twoStart, len) {
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
    this.m_nodePool = new _NodePool["default"]();
    this.m_tinyNodePool = new _NodePool["default"]();
    this.m_openList = new _NodeQueue["default"]();
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

      if (tile == null) return new FindRandomPointResult(_Status["default"].FAILURE, 0, null); // Randomly pick one polygon weighted by polygon area.

      var poly = null;
      var polyRef = 0;
      var base = this.m_nav.getPolyRefBase(tile);
      areaSum = 0.0;

      for (var _i = 0; _i < tile.data.header.polyCount; ++_i) {
        var p = tile.data.polys[_i]; // Do not return off-mesh connection polygons.

        if (p.getType() != _Poly["default"].DT_POLYTYPE_GROUND) continue; // Must pass filter

        var ref = this.or(base, _i);
        if (!filter.passFilter(ref, tile, p)) continue; // Calc area of the polygon.

        polyArea = 0.0;

        for (var j = 2; j < p.vertCount; ++j) {
          var va = p.verts[0] * 3;
          var vb = p.verts[j - 1] * 3;
          var vc = p.verts[j] * 3;
          polyArea += _DetourCommon["default"].triArea2D4(tile.data.verts, va, vb, vc);
        } // Choose random polygon weighted by area, using reservoi sampling.


        areaSum += polyArea;
        u = frand.frand();

        if (u * areaSum <= polyArea) {
          poly = p;
          polyRef = ref;
        }
      }

      if (poly == null) return new FindRandomPointResult(_Status["default"].FAILURE, 0, null); // Randomly pick poPoly on polygon.

      var verts = new Array(3 * this.m_nav.getMaxVertsPerPoly());
      var areas = new Array(this.m_nav.getMaxVertsPerPoly());
      arraycopy(tile.data.verts, poly.verts[0] * 3, verts, 0, 3);

      for (var _j = 1; _j < poly.vertCount; ++_j) {
        arraycopy(tile.data.verts, poly.verts[_j] * 3, verts, _j * 3, 3);
      }

      s = frand.frand();
      t = frand.frand();

      var pt = _DetourCommon["default"].randomPointInConvexPoly(verts, poly.vertCount, areas, s, t);

      pt[1] = getPolyHeight(polyRef, pt);
      return new FindRandomPointResult(_Status["default"].SUCCSESS, polyRef, pt);
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

      _DetourCommon["default"].vCopy(startNode.pos, centerPos);

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

        if (bestPoly.getType() == _Poly["default"].DT_POLYTYPE_GROUND) {
          // Calc area of the polygon.
          polyArea = 0.0;

          for (var j = 2; j < bestPoly.vertCount; ++j) {
            var va = bestPoly.verts[0] * 3;
            var vb = bestPoly.verts[j - 1] * 3;
            var vc = bestPoly.verts[j] * 3;
            polyArea += _DetourCommon["default"].triArea2D4(bestTile.data.verts, va, vb, vc);
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

        for (var i = bestPoly.firstLink; i != _NavMesh["default"].DT_NULL_LINK; i = bestTile.links[i].next) {
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

          var distseg = _DetourCommon["default"].distancePtSegSqr2D3(centerPos, _va, _vb);

          distSqr = distseg[0];
          if (distSqr > radiusSqr) continue;
          var neighbourNode = m_nodePool.getNode(neighbourRef);
          if ((neighbourNode.flags & _Node["default"].DT_NODE_CLOSED) != 0) continue; // Cost

          if (neighbourNode.flags == 0) neighbourNode.pos = _DetourCommon["default"].vLerp3(_va, _vb, 0.5);
          total = bestNode.total + _DetourCommon["default"].vDist2(bestNode.pos, neighbourNode.pos); // The node is already in open list and the new result is worse, skip.

          if ((neighbourNode.flags & _Node["default"].DT_NODE_OPEN) != 0 && total >= neighbourNode.total) continue;
          neighbourNode.id = neighbourRef;
          neighbourNode.flags = neighbourNode.flags & ~_Node["default"].DT_NODE_CLOSED;
          neighbourNode.pidx = m_nodePool.getNodeIdx(bestNode);
          neighbourNode.total = total;

          if ((neighbourNode.flags & _Node["default"].DT_NODE_OPEN) != 0) {
            this.m_openList.modify(neighbourNode);
          } else {
            neighbourNode.flags = _Node["default"].DT_NODE_OPEN;
            this.m_openList.push(neighbourNode);
          }
        }
      }

      if (randomPoly == null) return new FindRandomPointResult(_Status["default"].FAILURE, 0, null); // Randomly pick poPoly on polygon.

      var verts = new Array(3 * this.m_nav.getMaxVertsPerPoly());
      var areas = new Array(this.m_nav.getMaxVertsPerPoly());
      arraycopy(randomTile.data.verts, randomPoly.verts[0] * 3, verts, 0, 3);

      for (var _j2 = 1; _j2 < randomPoly.vertCount; ++_j2) {
        arraycopy(randomTile.data.verts, randomPoly.verts[_j2] * 3, verts, _j2 * 3, 3);
      }

      s = frand.frand();
      t = frand.frand();

      var pt = _DetourCommon["default"].randomPointInConvexPoly(verts, randomPoly.vertCount, areas, s, t);

      pt[1] = getPolyHeight(randomPolyRef, pt);
      return new FindRandomPointResult(_Status["default"].SUCCSESS, randomPolyRef, pt);
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

      if (poly.getType() == _Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION) {
        var _v = poly.verts[0] * 3;

        var _v2 = poly.verts[1] * 3;

        var _d = _DetourCommon["default"].vDist3(pos, tile.data.verts, _v);

        var _d2 = _DetourCommon["default"].vDist3(pos, tile.data.verts, _v2);

        var _u = _d / (_d + _d2);

        var _closest = _DetourCommon["default"].vLerp4(tile.data.verts, _v, _v2, _u);

        return new ClosesPointOnPolyResult(false, _closest);
      } // Clamp poPoly to be inside the polygon.


      var verts = new Array(this.m_nav.getMaxVertsPerPoly() * 3);
      var edged = new Array(this.m_nav.getMaxVertsPerPoly());
      var edget = new Array(this.m_nav.getMaxVertsPerPoly());
      var nv = poly.vertCount;

      for (var i = 0; i < nv; ++i) {
        arraycopy(tile.data.verts, poly.verts[i] * 3, verts, i * 3, 3);
      }

      var posOverPoly = false;
      var closest = new Array(3);

      _DetourCommon["default"].vCopy(closest, pos);

      if (!_DetourCommon["default"].distancePtPolyEdgesSqr(pos, verts, nv, edged, edget)) {
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
        closest = _DetourCommon["default"].vLerp4(verts, va, vb, edget[imin]);
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

          var heightResult = _DetourCommon["default"].closestHeightPointTriangle(closest, v[0], v[1], v[2]);

          if (heightResult[0]) {
            closest[1] = heightResult[1];
            break;
          }
        }
      }

      return new _ClosestPointOnPolyResult["default"](posOverPoly, closest);
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
        arraycopy(tile.data.verts, poly.verts[i] * 3, verts, i * 3, 3);
      }

      var closest;

      if (_DetourCommon["default"].distancePtPolyEdgesSqr(pos, verts, nv, edged, edget)) {
        closest = _DetourCommon["default"].vCopy_return(pos);
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
        closest = _DetourCommon["default"].vLerp4(verts, va, vb, edget[imin]);
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

      if (poly.getType() == _Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION) {
        var i = poly.verts[0] * 3;
        i = poly.verts[1] * 3;
        letv1 = [tile.data.verts[i], tile.data.verts[i + 1], tile.data.verts[i + 2]];
        d0 = _DetourCommon["default"].vDist2D(pos, v0);
        d1 = _DetourCommon["default"].vDist2D(pos, v1);
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

          var heightResult = _DetourCommon["default"].closestHeightPointTriangle(pos, v[0], v[1], v[2]);

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

        var diff = _DetourCommon["default"].vSub(center, closestPtPoly);

        if (posOverPoly) {
          var tilaAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(polys[i]);
          var _tile = tilaAndPoly[0];
          d = Math.abs(diff[1]) - _tile.data.header.walkableClimb;
          d = d > 0 ? d * d : 0;
        } else {
          d = _DetourCommon["default"].vLenSqr(diff);
        }

        if (d < nearestDistanceSqr) {
          nearestPt = closestPtPoly;
          nearestDistanceSqr = d;
          nearest = ref;
        }
      }

      return new _FindNearestPolyResult["default"](nearest, nearestPt);
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

        var minx = _DetourCommon["default"].clamp(qmin[0], tbmin[0], tbmax[0]) - tbmin[0];
        var miny = _DetourCommon["default"].clamp(qmin[1], tbmin[1], tbmax[1]) - tbmin[1];
        var minz = _DetourCommon["default"].clamp(qmin[2], tbmin[2], tbmax[2]) - tbmin[2];
        var maxx = _DetourCommon["default"].clamp(qmax[0], tbmin[0], tbmax[0]) - tbmin[0];
        var maxy = _DetourCommon["default"].clamp(qmax[1], tbmin[1], tbmax[1]) - tbmin[1];
        var maxz = _DetourCommon["default"].clamp(qmax[2], tbmin[2], tbmax[2]) - tbmin[2]; // Quantize

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

          var overlap = _DetourCommon["default"].overlapQuantBounds(bmin, bmax, node.bmin, node.bmax);

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

          if (p.getType() == _Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION) continue;

          var _ref = _base | i;

          if (!filter.passFilter(_ref, tile, p)) continue; // Calc polygon bounds.

          var v = p.verts[0] * 3;

          _DetourCommon["default"].vCopy(_bmin, tile.data.verts, v);

          _DetourCommon["default"].vCopy(_bmax, tile.data.verts, v);

          for (var j = 1; j < p.vertCount; ++j) {
            v = p.verts[j] * 3;

            _DetourCommon["default"].vMin(_bmin, tile.data.verts, v);

            _DetourCommon["default"].vMax(_bmax, tile.data.verts, v);
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
      var bmin = _DetourCommon["default"].vSub(center, extents);

      var bmax = _DetourCommon["default"].vAdd(center, extents); // Find tiles the query touches.


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

        return new _FindPathResult["default"](_Status["default"].SUCCSESS, _path);
      }

      this.m_nodePool = [];
      this.m_openList = [];
      var startNode = this.m_nodePool.getNode(startRef);

      _DetourCommon["default"].vCopy(startNode.pos, startPos);

      startNode.pidx = 0;
      startNode.cost = 0;
      startNode.total = _DetourCommon["default"].vDist2(startPos, endPos) * NavMeshQuery.H_SCALE;
      startNode.id = startRef;
      startNode.flags = _Node["default"].DT_NODE_OPEN;
      this.m_openList.push(startNode);
      var lastBestNode = startNode;
      lastBestNodeCost = startNode.total;
      var status = _Status["default"].SUCCSESS;

      while (!this.m_openList.length == 0) {
        // Remove node from open list and put it in closed list.
        var bestNode = this.m_openList.pop();
        bestNode.flags &= ~_Node["default"].DT_NODE_OPEN;
        bestNode.flags |= _Node["default"].DT_NODE_CLOSED; // Reached the goal, stop searching.

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

        for (var i = bestPoly.firstLink; i != _NavMesh["default"].DT_NULL_LINK; i = bestTile.links[i].next) {
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
            heuristic = _DetourCommon["default"].vDist2(neighbourNode.pos, endPos) * NavMeshQuery.H_SCALE;
          }

          total = cost + heuristic; // The node is already in open list and the new result is worse, skip.

          if ((neighbourNode.flags & _Node["default"].DT_NODE_OPEN) != 0 && total >= neighbourNode.total) continue; // The node is already visited and process, and the new result is worse, skip.

          if ((neighbourNode.flags & _Node["default"].DT_NODE_CLOSED) != 0 && total >= neighbourNode.total) continue; // Add or update the node.

          neighbourNode.pidx = this.m_nodePool.getNodeIdx(bestNode);
          neighbourNode.id = neighbourRef;
          neighbourNode.flags = neighbourNode.flags & ~_Node["default"].DT_NODE_CLOSED;
          neighbourNode.cost = cost;
          neighbourNode.total = total;

          if ((neighbourNode.flags & _Node["default"].DT_NODE_OPEN) != 0) {
            // Already in open, update node location.
            this.m_openList.modify(neighbourNode);
          } else {
            // Put the node in open list.
            neighbourNode.flags |= _Node["default"].DT_NODE_OPEN;
            this.m_openList.push(neighbourNode);
          } // Update nearest node to target so far.


          if (heuristic < lastBestNodeCost) {
            lastBestNodeCost = heuristic;
            lastBestNode = neighbourNode;
          }
        }
      }

      var path = getPathToNode(lastBestNode);
      if (lastBestNode.id != endRef) status = _Status["default"].PARTIAL_RESULT;
      return new _FindPathResult["default"](status, path);
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
      this.m_query = new _QueryData["default"]();
      this.m_query.status = _Status["default"].FAILURE;
      this.m_query.startRef = startRef;
      this.m_query.endRef = endRef;

      _DetourCommon["default"].vCopy(this.m_query.startPos, startPos);

      _DetourCommon["default"].vCopy(this.m_query.endPos, endPos);

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
        this.m_query.raycastLimitSqr = _DetourCommon["default"].sqr(agentRadius * _NavMesh["default"].DT_RAY_CAST_LIMIT_PROPORTIONS);
      }

      if (startRef == endRef) {
        this.m_query.status = _Status["default"].SUCCSESS;
        return _Status["default"].SUCCSESS;
      }

      this.m_nodePool.clear();
      this.m_openList.clear();
      var startNode = this.m_nodePool.getNode(startRef);

      _DetourCommon["default"].vCopy(startNode.pos, startPos);

      startNode.pidx = 0;
      startNode.cost = 0;
      startNode.total = _DetourCommon["default"].vDist2(startPos, endPos) * NavMeshQuery.H_SCALE;
      startNode.id = startRef;
      startNode.flags = _Node["default"].DT_NODE_OPEN;
      this.m_openList.push(startNode);
      this.m_query.status = _Status["default"].IN_PROGRESS;
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
      if (this.m_query.status != _Status["default"].IN_PROGRESS) return new _UpdateSlicedPathResult["default"](this.m_query.status, 0); // Make sure the request is still valid.

      if (!this.m_nav.isValidPolyRef(this.m_query.startRef) || !this.m_nav.isValidPolyRef(this.m_query.endRef)) {
        this.m_query.status = _Status["default"].FAILURE;
        return new _UpdateSlicedPathResult["default"](this.m_query.status, 0);
      }

      var iter = 0;

      while (iter < maxIter && !this.m_openList.isEmpty()) {
        iter++; // Remove node from open list and put it in closed list.

        var bestNode = this.m_openList.pop();
        bestNode.flags &= ~_Node["default"].DT_NODE_OPEN;
        bestNode.flags |= _Node["default"].DT_NODE_CLOSED; // Reached the goal, stop searching.

        if (bestNode.id == this.m_query.endRef) {
          this.m_query.lastBestNode = bestNode;
          this.m_query.status = _Status["default"].SUCCSESS;
          return new _UpdateSlicedPathResult["default"](this.m_query.status, iter);
        } // Get current let and tile.
        // The API input has been cheked already, skip checking internal
        // data.


        var bestRef = bestNode.id;
        var tileAndPoly = void 0;

        try {
          tileAndPoly = this.m_nav.getTileAndPolyByRef(bestRef);
        } catch (e) {
          this.m_query.status = _Status["default"].FAILURE; // The polygon has disappeared during the sliced query, fail.

          return new _UpdateSlicedPathResult["default"](this.m_query.status, iter);
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
            this.m_query.status = _Status["default"].FAILURE;
            return new _UpdateSlicedPathResult["default"](this.m_query.status, iter);
          }
        } // decide whether to test raycast to previous nodes


        var tryLOS = false;

        if ((this.m_query.options & NavMeshQuery.DT_FINDPATH_ANY_ANGLE) != 0) {
          if (parentRef != 0 && _DetourCommon["default"].vDistSqr(parentNode.pos, bestNode.pos) < this.m_query.raycastLimitSqr) tryLOS = true;
        }

        for (var i = bestPoly.firstLink; i != _NavMesh["default"].DT_NULL_LINK; i = bestTile.links[i].next) {
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
            _heuristic = _DetourCommon["default"].vDist2(neighbourNode.pos, this.m_query.endPos) * NavMeshQuery.H_SCALE;
          }

          var _total = _cost + _heuristic; // The node is already in open list and the new result is worse,
          // skip.


          if ((neighbourNode.flags & _Node["default"].DT_NODE_OPEN) != 0 && _total >= neighbourNode.total) continue; // The node is already visited and process, and the new result
          // is worse, skip.

          if ((neighbourNode.flags & _Node["default"].DT_NODE_CLOSED) != 0 && _total >= neighbourNode.total) continue; // Add or update the node.

          neighbourNode.pidx = foundShortCut ? bestNode.pidx : this.m_nodePool.getNodeIdx(bestNode);
          neighbourNode.id = neighbourRef;
          neighbourNode.flags = this.and(neighbourNode.flags, ~this.or(_Node["default"].DT_NODE_CLOSED, _Node["default"].DT_NODE_PARENT_DETACHED));
          neighbourNode.cost = _cost;
          neighbourNode.total = _total;
          if (foundShortCut) neighbourNode.flags = this.or(neighbourNode.flags, _Node["default"].DT_NODE_PARENT_DETACHED);

          if ((neighbourNode.flags & _Node["default"].DT_NODE_OPEN) != 0) {
            // Already in open, update node location.
            this.m_openList.modify(neighbourNode);
          } else {
            // Put the node in open list.
            neighbourNode.flags |= _Node["default"].DT_NODE_OPEN;
            this.m_openList.push(neighbourNode);
          } // Update nearest node to target so far.


          if (_heuristic < this.m_query.lastBestNodeCost) {
            this.m_query.lastBestNodeCost = _heuristic;
            this.m_query.lastBestNode = neighbourNode;
          }
        }
      } // Exhausted all nodes, but could not find path.


      if (this.m_openList.isEmpty()) {
        this.m_query.status = _Status["default"].PARTIAL_RESULT;
      }

      return new _UpdateSlicedPathResult["default"](this.m_query.status, iter);
    } /// Finalizes and returns the results of a sliced path query.
    ///  @param[out]	path		An ordered list of polygon references representing the path. (Start to end.) 
    ///  							[(polyRef) * @p pathCount]
    /// @returns The status flags for the query.

  }, {
    key: "finalizeSlicedFindPath",
    value: function finalizeSlicedFindPath() {
      var path = [];

      if (this.m_query.status == _Status["default"].FAILURE) {
        // Reset query.
        this.m_query = new _QueryData["default"]();
        return new _FindPathResult["default"](_Status["default"].FAILURE, path);
      }

      if (this.m_query.startRef == this.m_query.endRef) {
        // Special case: the search starts and ends at same poly.
        path.push(this.m_query.startRef);
      } else {
        // Reverse the path.
        if (this.m_query.lastBestNode.id != this.m_query.endRef) this.m_query.status = _Status["default"].PARTIAL_RESULT;
        var prev = null;
        var node = this.m_query.lastBestNode;
        var prevRay = 0;

        do {
          var next = this.m_nodePool.getNodeAtIdx(node.pidx);
          node.pidx = this.m_nodePool.getNodeIdx(prev);
          prev = node;
          var nextRay = node.flags & _Node["default"].DT_NODE_PARENT_DETACHED; // keep track of whether parent is not adjacent (i.e. due to raycast shortcut)

          node.flags = this.or(this.and(node.flags, ~_Node["default"].DT_NODE_PARENT_DETACHED), prevRay); // and store it in the reversed path's node

          prevRay = nextRay;
          node = next;
        } while (node != null); // Store path


        node = prev;

        do {
          var _next = this.m_nodePool.getNodeAtIdx(node.pidx);

          if ((node.flags & _Node["default"].DT_NODE_PARENT_DETACHED) != 0) {
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

      this.m_query = new _QueryData["default"]();
      return new _FindPathResult["default"](status, path);
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
        return new _FindPathResult["default"](_Status["default"].FAILURE, path);
      }

      if (this.m_query.status == _Status["default"].FAILURE) {
        // Reset query.
        this.m_query = new _QueryData["default"]();
        return new _FindPathResult["default"](_Status["default"].FAILURE, path);
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
          this.m_query.status = _Status["default"].PARTIAL_RESULT;
          node = this.m_query.lastBestNode;
        } // Reverse the path.


        var prevRay = 0;

        do {
          var next = this.m_nodePool.getNodeAtIdx(node.pidx);
          node.pidx = this.m_nodePool.getNodeIdx(prev);
          prev = node;
          var nextRay = this.and(node.flags, _Node["default"].DT_NODE_PARENT_DETACHED); // keep track of whether parent is not adjacent (i.e. due to raycast shortcut)

          node.flags = this.or(this.and(node.flags & ~_Node["default"].DT_NODE_PARENT_DETACHED), prevRay); // and store it in the reversed path's node

          prevRay = nextRay;
          node = next;
        } while (node != null); // Store path


        node = prev;

        do {
          var _next2 = this.m_nodePool.getNodeAtIdx(node.pidx);

          if ((node.flags & _Node["default"].DT_NODE_PARENT_DETACHED) != 0) {
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

      this.m_query = new _QueryData["default"]();
      return new _FindPathResult["default"](status, path);
    }
  }, {
    key: "appendVertex",
    value: function appendVertex(pos, flags, ref, straightPath, maxStraightPath) {
      if (straightPath.length > 0 && _DetourCommon["default"].vEqual(straightPath[straightPath.length - 1].pos, pos)) {
        // The vertices are equal, update flags and poly.
        straightPath[straightPath.length - 1].flags = flags;
        straightPath[straightPath.length - 1].ref = ref;
      } else {
        if (straightPath.length < maxStraightPath) {
          // Append new vertex.
          straightPath.push(new _StraightPathItem["default"](pos, flags, ref));
        } // If reached end of path or there is no space to append more vertices, return.


        if (flags == NavMeshQuery.DT_STRAIGHTPATH_END || straightPath.length >= maxStraightPath) {
          return _Status["default"].SUCCSESS;
        }
      }

      return _Status["default"].IN_PROGRESS;
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


        var interect = _DetourCommon["default"].intersectSegSeg2D(startPos, endPos, left, right);

        if (interect[0]) {
          t = interect.third;

          var pt = _DetourCommon["default"].vLerp3(left, right, t);

          stat = this.appendVertex(pt, 0, path[i + 1], straightPath, maxStraightPath);
          if (!stat == _Status["default"].IN_PROGRESS) return stat;
        }
      }

      return _Status["default"].IN_PROGRESS;
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
      if (!stat == _Status["default"].IN_PROGRESS) return straightPath;

      if (path.length > 1) {
        var portalApex = _DetourCommon["default"].vCopy_return(closestStartPos);

        var portalLeft = _DetourCommon["default"].vCopy_return(portalApex);

        var portalRight = _DetourCommon["default"].vCopy_return(portalApex);

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
                if (!stat == _Status["default"].IN_PROGRESS) return straightPath;
              }

              this.appendVertex(closestEndPos, 0, path[i], straightPath, maxStraightPath);
              return straightPath;
            } // If starting really close the portal, advance.


            if (i == 0) {
              var dt = _DetourCommon["default"].distancePtSegSqr2D3(portalApex, left, right);

              if (dt[0] < _DetourCommon["default"].sqr(0.001)) continue;
            }
          } else {
            // End of the path.
            left = _DetourCommon["default"].vCopy_return(closestEndPos);
            right = _DetourCommon["default"].vCopy_return(closestEndPos);
            toType = _Poly["default"].DT_POLYTYPE_GROUND;
          } // Right vertex.


          if (_DetourCommon["default"].triArea2D3(portalApex, portalRight, right) <= 0.0) {
            if (_DetourCommon["default"].vEqual(portalApex, portalRight) || _DetourCommon["default"].triArea2D3(portalApex, portalLeft, right) > 0.0) {
              portalRight = _DetourCommon["default"].vCopy_return(right);
              rightPolyRef = i + 1 < path.length ? path[i + 1] : 0;
              rightPolyType = toType;
              rightIndex = i;
            } else {
              // Append portals aPoly the current straight path segment.
              if (this.and(options, this.or(NavMeshQuery.DT_STRAIGHTPATH_AREA_CROSSINGS, NavMeshQuery.DT_STRAIGHTPATH_ALL_CROSSINGS)) != 0) {
                stat = appendPortals(apexIndex, leftIndex, portalLeft, path, straightPath, options, maxStraightPath);
                if (!stat == _Status["default"].IN_PROGRESS) return straightPath;
              }

              portalApex = _DetourCommon["default"].vCopy_return(portalLeft);
              apexIndex = leftIndex;
              var flags = 0;
              if (leftPolyRef == 0) flags = NavMeshQuery.DT_STRAIGHTPATH_END;else if (leftPolyType == _Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION) flags = NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION;
              var ref = leftPolyRef; // Append or update vertex

              stat = this.appendVertex(portalApex, flags, ref, straightPath, maxStraightPath);
              if (!stat == _Status["default"].IN_PROGRESS) return straightPath;
              portalLeft = _DetourCommon["default"].vCopy_return(portalApex);
              portalRight = _DetourCommon["default"].vCopy_return(portalApex);
              leftIndex = apexIndex;
              rightIndex = apexIndex; // Restart

              i = apexIndex;
              continue;
            }
          } // Left vertex.


          if (_DetourCommon["default"].triArea2D3(portalApex, portalLeft, left) >= 0.0) {
            if (_DetourCommon["default"].vEqual(portalApex, portalLeft) || _DetourCommon["default"].triArea2D3(portalApex, portalRight, left) < 0.0) {
              portalLeft = _DetourCommon["default"].vCopy_return(left);
              leftPolyRef = i + 1 < path.length ? path[i + 1] : 0;
              leftPolyType = toType;
              leftIndex = i;
            } else {
              // Append portals aPoly the current straight path segment.
              if (this.and(options & this.or(NavMeshQuery.DT_STRAIGHTPATH_AREA_CROSSINGS, NavMeshQuery.DT_STRAIGHTPATH_ALL_CROSSINGS)) != 0) {
                stat = appendPortals(apexIndex, rightIndex, portalRight, path, straightPath, options, maxStraightPath);
                if (!stat == _Status["default"].IN_PROGRESS) return straightPath;
              }

              portalApex = _DetourCommon["default"].vCopy_return(portalRight);
              apexIndex = rightIndex;
              var _flags = 0;
              if (rightPolyRef == 0) _flags = NavMeshQuery.DT_STRAIGHTPATH_END;else if (rightPolyType == _Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION) _flags = NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION;
              var _ref2 = rightPolyRef; // Append or update vertex

              stat = this.appendVertex(portalApex, _flags, _ref2, straightPath, maxStraightPath);
              if (!stat == _Status["default"].IN_PROGRESS) return straightPath;
              portalLeft = _DetourCommon["default"].vCopy_return(portalApex);
              portalRight = _DetourCommon["default"].vCopy_return(portalApex);
              leftIndex = apexIndex;
              rightIndex = apexIndex; // Restart

              i = apexIndex;
              continue;
            }
          }
        } // Append portals aPoly the current straight path segment.


        if (this.and(options & this.or(NavMeshQuery.DT_STRAIGHTPATH_AREA_CROSSINGS, NavMeshQuery.DT_STRAIGHTPATH_ALL_CROSSINGS)) != 0) {
          stat = appendPortals(apexIndex, path.length - 1, closestEndPos, path, straightPath, options, maxStraightPath);
          if (!stat == _Status["default"].IN_PROGRESS) return straightPath;
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
      this.m_tinyNodePool = new _NodePool["default"]();
      var startNode = this.m_tinyNodePool.getNode(startRef);
      startNode.pidx = 0;
      startNode.cost = 0;
      startNode.total = 0;
      startNode.id = startRef;
      startNode.flags = _Node["default"].DT_NODE_CLOSED;
      var stack = [];
      stack.push(startNode);
      var bestPos = new Array(3);
      var bestDist = Number.MAX_VALUE;
      var bestNode = null;

      _DetourCommon["default"].vCopy(bestPos, startPos); // Search constraints


      var searchPos = _DetourCommon["default"].vLerp3(startPos, endPos, 0.5);

      var searchRadSqr = _DetourCommon["default"].sqr(_DetourCommon["default"].vDist2(startPos, endPos) / 2.0 + 0.001);

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
          arraycopy(curTile.data.verts, curPoly.verts[i] * 3, verts, i * 3, 3);
        } // If target is inside the poly, stop search.


        if (_DetourCommon["default"].pointInPolygon(endPos, verts, nverts)) {
          bestNode = curNode;

          _DetourCommon["default"].vCopy(bestPos, endPos);

          break;
        } // Find wall edges and find nearest poPoly inside the walls.


        for (var _i4 = 0, j = curPoly.vertCount - 1; _i4 < curPoly.vertCount; j = _i4++) {
          // Find links to neighbours.
          var MAX_NEIS = 8;
          var nneis = 0;
          var neis = new Array(MAX_NEIS);

          if ((curPoly.neis[j] & _NavMesh["default"].DT_EXT_LINK) != 0) {
            // Tile border.
            for (var k = curPoly.firstLink; k != _NavMesh["default"].DT_NULL_LINK; k = curTile.links[k].next) {
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

            var distSeg = _DetourCommon["default"].distancePtSegSqr2D4(endPos, verts, vj, vi);

            var _distSqr = distSeg[0];
            var _tseg = distSeg[1];

            if (_distSqr < bestDist) {
              // Update nearest distance.
              bestPos = _DetourCommon["default"].vLerp4(verts, vj, vi, _tseg);
              bestDist = _distSqr;
              bestNode = curNode;
            }
          } else {
            for (var _k = 0; _k < nneis; ++_k) {
              // Skip if no node can be allocated.
              var neighbourNode = this.m_tinyNodePool.getNode(neis[_k]);
              if (neighbourNode == null) continue; // Skip if already visited.

              if ((neighbourNode.flags & _Node["default"].DT_NODE_CLOSED) != 0) continue; // Skip the link if it is too far from search constraint.
              // TODO: Maybe should use getPortalPoints(), but this one is way faster.

              var _vj = j * 3;

              var _vi = _i4 * 3;

              var distseg = _DetourCommon["default"].distancePtSegSqr2D4(searchPos, verts, _vj, _vi);

              var _distSqr2 = distseg[0];
              if (_distSqr2 > searchRadSqr) continue; // Mark as the node as visited and push to queue.

              neighbourNode.pidx = this.m_tinyNodePool.getNodeIdx(curNode);
              neighbourNode.flags |= _Node["default"].DT_NODE_CLOSED;
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

      return new _MoveAlongSurfaceResult["default"](bestPos, visited);
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

      for (var i = fromPoly.firstLink; i != _NavMesh["default"].DT_NULL_LINK; i = fromTile.links[i].next) {
        if (fromTile.links[i].ref == to) {
          link = fromTile.links[i];
          break;
        }
      }

      if (link == null) throw new IllegalArgumentException("Null link"); // Handle off-mesh connections.

      if (fromPoly.getType() == _Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION) {
        // Find link that points to first vertex.
        for (var _i5 = fromPoly.firstLink; _i5 != _NavMesh["default"].DT_NULL_LINK; _i5 = fromTile.links[_i5].next) {
          if (fromTile.links[_i5].ref == to) {
            var v = fromTile.links[_i5].edge;
            arraycopy(fromTile.data.verts, fromPoly.verts[v] * 3, left, 0, 3);
            arraycopy(fromTile.data.verts, fromPoly.verts[v] * 3, right, 0, 3);
            return new PortalResult(left, right, fromType, toType);
          }
        }

        throw new IllegalArgumentException("Invalid offmesh from connection");
      }

      if (toPoly.getType() == _Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION) {
        for (var _i6 = toPoly.firstLink; _i6 != _NavMesh["default"].DT_NULL_LINK; _i6 = toTile.links[_i6].next) {
          if (toTile.links[_i6].ref == from) {
            var _v3 = toTile.links[_i6].edge;
            arraycopy(toTile.data.verts, toPoly.verts[_v3] * 3, left, 0, 3);
            arraycopy(toTile.data.verts, toPoly.verts[_v3] * 3, right, 0, 3);
            return new PortalResult(left, right, fromType, toType);
          }
        }

        throw new IllegalArgumentException("Invalid offmesh to connection");
      } // Find portal vertices.


      var v0 = fromPoly.verts[link.edge];
      var v1 = fromPoly.verts[(link.edge + 1) % fromPoly.vertCount];
      arraycopy(fromTile.data.verts, v0 * 3, left, 0, 3);
      arraycopy(fromTile.data.verts, v1 * 3, right, 0, 3); // If the link is at tile boundary, dtClamp the vertices to
      // the link width.

      if (link.side != 0xff) {
        // Unpack portal limits.
        if (link.bmin != 0 || link.bmax != 255) {
          s = 1.0 / 255.0;
          tmin = link.bmin * s;
          tmax = link.bmax * s;
          left = _DetourCommon["default"].vLerp3(fromTile.data.verts, v0 * 3, v1 * 3, tmin);
          right = _DetourCommon["default"].vLerp3(fromTile.data.verts, v0 * 3, v1 * 3, tmax);
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
      var hit = new _RaycastHit["default"]();
      var verts = new Array(this.m_nav.getMaxVertsPerPoly() * 3 + 3);
      var curPos = new Array(3),
          lastPos = new Array(3);

      _DetourCommon["default"].vCopy(curPos, startPos);

      var dir = _DetourCommon["default"].vSub(endPos, startPos);

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
          arraycopy(tile.data.verts, poly.verts[i] * 3, verts, nv * 3, 3);
          nv++;
        }

        var iresult = _DetourCommon["default"].intersectSegmentPoly2D(startPos, endPos, verts, nv);

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

        for (var _i7 = poly.firstLink; _i7 != _NavMesh["default"].DT_NULL_LINK; _i7 = tile.links[_i7].next) {
          var link = tile.links[_i7]; // Find link which contains this edge.

          if (link.edge != iresult.segMax) continue; // Get pointer to the next polygon.

          tileAndPolyUns = this.m_nav.getTileAndPolyByRefUnsafe(link.ref);
          nextTile = tileAndPolyUns[0];
          nextPoly = tileAndPolyUns[1]; // Skip off-mesh connections.

          if (nextPoly.getType() == _Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION) continue; // Skip links based on filter.

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
          _DetourCommon["default"].vCopy(lastPos, curPos);

          curPos = _DetourCommon["default"].vMad(startPos, dir, hit.t);
          var e1 = new VectorPtr(verts, iresult.segMax * 3);
          var e2 = new VectorPtr(verts, (iresult.segMax + 1) % nv * 3);

          var eDir = _DetourCommon["default"].vSub(e2, e1);

          var diff = _DetourCommon["default"].vSub(new VectorPtr(curPos), e1);

          s = _DetourCommon["default"].sqr(eDir[0]) > _DetourCommon["default"].sqr(eDir[2]) ? diff[0] / eDir[0] : diff[2] / eDir[2];
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

      _DetourCommon["default"].vCopy(startNode.pos, centerPos);

      startNode.pidx = 0;
      startNode.cost = 0;
      startNode.total = 0;
      startNode.id = startRef;
      startNode.flags = _Node["default"].DT_NODE_OPEN;
      this.m_openList.push(startNode);
      radiusSqr = _DetourCommon["default"].sqr(radius);

      while (!this.m_openList.length == 0) {
        var bestNode = m_openList.pop();
        bestNode.flags &= ~_Node["default"].DT_NODE_OPEN;
        bestNode.flags |= _Node["default"].DT_NODE_CLOSED; // Get let and tile.
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

        for (var i = bestPoly.firstLink; i != _NavMesh["default"].DT_NULL_LINK; i = bestTile.links[i].next) {
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

          var distseg = _DetourCommon["default"].distancePtSegSqr2D3(centerPos, va, vb);

          distSqr = distseg[0];
          if (distSqr > radiusSqr) continue;
          var neighbourNode = this.m_nodePool.getNode(neighbourRef);
          if ((neighbourNode.flags & _Node["default"].DT_NODE_CLOSED) != 0) continue; // Cost

          if (neighbourNode.flags == 0) neighbourNode.pos = vLerp3(va, vb, 0.5);
          cost = filter.getCost(bestNode.pos, neighbourNode.pos, parentRef, parentTile, parentPoly, bestRef, bestTile, bestPoly, neighbourRef, neighbourTile, neighbourPoly);
          total = bestNode.total + cost; // The node is already in open list and the new result is worse, skip.

          if ((neighbourNode.flags & _Node["default"].DT_NODE_OPEN) != 0 && total >= neighbourNode.total) continue;
          neighbourNode.id = neighbourRef;
          neighbourNode.pidx = this.m_nodePool.getNodeIdx(bestNode);
          neighbourNode.total = total;

          if ((neighbourNode.flags & _Node["default"].DT_NODE_OPEN) != 0) {
            this.m_openList.modify(neighbourNode);
          } else {
            neighbourNode.flags = _Node["default"].DT_NODE_OPEN;
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

      _DetourCommon["default"].vCopy(startNode.pos, centerPos);

      startNode.pidx = 0;
      startNode.cost = 0;
      startNode.total = 0;
      startNode.id = startRef;
      startNode.flags = _Node["default"].DT_NODE_OPEN;
      this.m_openList.push(startNode);

      while (!this.m_openList.length == 0) {
        var bestNode = this.m_openList.pop();
        bestNode.flags &= ~_Node["default"].DT_NODE_OPEN;
        bestNode.flags |= _Node["default"].DT_NODE_CLOSED; // Get let and tile.
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

        for (var _i8 = bestPoly.firstLink; _i8 != _NavMesh["default"].DT_NULL_LINK; _i8 = bestTile.links[_i8].next) {
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

          var ir = _DetourCommon["default"].intersectSegmentPoly2D(va, vb, verts, nverts);

          if (!ir.intersects) continue;
          if (ir.tmin > 1.0 || ir.tmax < 0.0) continue;
          var neighbourNode = this.m_nodePool.getNode(neighbourRef);
          if ((neighbourNode.flags & _Node["default"].DT_NODE_CLOSED) != 0) continue; // Cost

          if (neighbourNode.flags == 0) neighbourNode.pos = _DetourCommon["default"].vLerp3(va, vb, 0.5);
          cost = filter.getCost(bestNode.pos, neighbourNode.pos, parentRef, parentTile, parentPoly, bestRef, bestTile, bestPoly, neighbourRef, neighbourTile, neighbourPoly);
          total = bestNode.total + cost; // The node is already in open list and the new result is worse, skip.

          if ((neighbourNode.flags & _Node["default"].DT_NODE_OPEN) != 0 && total >= neighbourNode.total) continue;
          neighbourNode.id = neighbourRef;
          neighbourNode.pidx = this.m_nodePool.getNodeIdx(bestNode);
          neighbourNode.total = total;

          if ((neighbourNode.flags & _Node["default"].DT_NODE_OPEN) != 0) {
            this.m_openList.modify(neighbourNode);
          } else {
            neighbourNode.flags = _Node["default"].DT_NODE_OPEN;
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
      this.m_tinyNodePool = new _NodePool["default"]();
      var startNode = this.m_tinyNodePool.getNode(startRef);
      startNode.pidx = 0;
      startNode.id = startRef;
      startNode.flags = _Node["default"].DT_NODE_CLOSED;
      var stack = [];
      stack.push(startNode);
      resultRef.push(startNode.id);
      resultParent.push(0);

      var radiusSqr = _DetourCommon["default"].sqr(radius);

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

        for (var i = curPoly.firstLink; i != _NavMesh["default"].DT_NULL_LINK; i = curTile.links[i].next) {
          var link = curTile.links[i];
          var neighbourRef = link.ref; // Skip invalid neighbours.

          if (neighbourRef == 0) continue; // Skip if cannot alloca more nodes.

          var neighbourNode = this.m_tinyNodePool.getNode(neighbourRef);
          if (neighbourNode == null) continue; // Skip visited.

          if ((neighbourNode.flags & _Node["default"].DT_NODE_CLOSED) != 0) continue; // Expand to neighbour

          tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
          var neighbourTile = tileAndPoly[0];
          var neighbourPoly = tileAndPoly[1]; // Skip off-mesh connections.

          if (neighbourPoly.getType() == _Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION) continue; // Do not advance if the polygon is excluded by the filter.

          if (!filter.passFilter(neighbourRef, neighbourTile, neighbourPoly)) continue; // Find edge and calc distance to the edge.

          var pp = this.getPortalPoints7(curRef, curPoly, curTile, neighbourRef, neighbourPoly, neighbourTile, 0, 0);
          var va = pp.left;
          var vb = pp.right; // If the circle is not touching the next polygon, skip it.

          var distseg = _DetourCommon["default"].distancePtSegSqr2D3(centerPos, va, vb);

          var _distSqr3 = distseg[0];
          if (_distSqr3 > radiusSqr) continue; // Mark node visited, this is done before the overlap test so that
          // we will not visit the let again if the test fails.

          neighbourNode.flags |= _Node["default"].DT_NODE_CLOSED;
          neighbourNode.pidx = this.m_tinyNodePool.getNodeIdx(curNode); // Check that the polygon does not collide with existing polygons.
          // Collect vertices of the neighbour poly.

          var npa = neighbourPoly.vertCount;

          for (var k = 0; k < npa; ++k) {
            arraycopy(neighbourTile.data.verts, neighbourPoly.verts[k] * 3, pa, k * 3, 3);
          }

          var overlap = false;

          for (var j = 0; j < resultRef.length; ++j) {
            var pastRef = resultRef[j]; // Connected polys do not overlap.

            var connected = false;

            for (var _k2 = curPoly.firstLink; _k2 != _NavMesh["default"].DT_NULL_LINK; _k2 = curTile.links[_k2].next) {
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
              arraycopy(pastTile.data.verts, pastPoly.verts[_k3] * 3, pb, _k3 * 3, 3);
            }

            if (_DetourCommon["default"].overlapPolyPoly2D(pa, npa, pb, npb)) {
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

      return new _FindLocalNeighbourhoodResult["default"](resultRef, resultParent);
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

        if ((poly.neis[j] & _NavMesh["default"].DT_EXT_LINK) != 0) {
          // Tile border.
          for (var k = poly.firstLink; k != _NavMesh["default"].DT_NULL_LINK; k = tile.links[k].next) {
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
          arraycopy(tile.data.verts, _vj2, seg, 0, 3);
          arraycopy(tile.data.verts, _vi2, seg, 3, 3);
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

            arraycopy(vLerp4(tile.data.verts, vj, vi, tmin), 0, _seg, 0, 3);
            arraycopy(vLerp4(tile.data.verts, vj, vi, tmax), 0, _seg, 3, 3);
            segmentVerts.push(_seg);
            segmentRefs.push(ints[_k4].ref);
          } // Wall segment.


          var imin = ints[_k4 - 1].tmax;
          var imax = ints[_k4].tmin;

          if (imin != imax) {
            tmin = imin / 255.0;
            tmax = imax / 255.0;

            var _seg2 = new Array(6);

            arraycopy(vLerp4(tile.data.verts, vj, vi, tmin), 0, _seg2, 0, 3);
            arraycopy(vLerp4(tile.data.verts, vj, vi, tmax), 0, _seg2, 3, 3);
            segmentVerts.push(_seg2);
            segmentRefs.push(0);
          }
        }
      }

      return new _GetPolyWallSegmentsResult["default"](segmentVerts, segmentRefs);
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

      _DetourCommon["default"].vCopy(startNode.pos, centerPos);

      startNode.pidx = 0;
      startNode.cost = 0;
      startNode.total = 0;
      startNode.id = startRef;
      startNode.flags = _Node["default"].DT_NODE_OPEN;
      this.m_openList.push(startNode);
      radiusSqr = _DetourCommon["default"].sqr(maxRadius);
      var hitPos = new Array(3);
      var bestvj = null;
      var bestvi = null;

      while (!this.m_openList.length == 0) {
        var bestNode = this.m_openList.pop();
        bestNode.flags &= ~_Node["default"].DT_NODE_OPEN;
        bestNode.flags |= _Node["default"].DT_NODE_CLOSED; // Get let and tile.
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
          if ((bestPoly.neis[j] & _NavMesh["default"].DT_EXT_LINK) != 0) {
            // Tile border.
            var solid = true;

            for (var k = bestPoly.firstLink; k != _NavMesh["default"].DT_NULL_LINK; k = bestTile.links[k].next) {
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

          var distseg = _DetourCommon["default"].distancePtSegSqr2D4(centerPos, bestTile.data.verts, vj, vi);

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

        for (var _i9 = bestPoly.firstLink; _i9 != _NavMesh["default"].DT_NULL_LINK; _i9 = bestTile.links[_i9].next) {
          var _link = bestTile.links[_i9];
          var neighbourRef = _link.ref; // Skip invalid neighbours and do not follow back to parent.

          if (neighbourRef == 0 || neighbourRef == parentRef) continue; // Expand to neighbour.

          tileAndPoly = this.m_nav.getTileAndPolyByRefUnsafe(neighbourRef);
          var neighbourTile = tileAndPoly[0];
          var neighbourPoly = tileAndPoly[1]; // Skip off-mesh connections.

          if (neighbourPoly.getType() == _Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION) continue; // Calc distance to the edge.

          var va = bestPoly.verts[_link.edge] * 3;
          var vb = bestPoly.verts[(_link.edge + 1) % bestPoly.vertCount] * 3;

          var _distseg = _DetourCommon["default"].distancePtSegSqr2D4(centerPos, bestTile.data.verts, va, vb);

          distSqr = _distseg[0]; // If the circle is not touching the next polygon, skip it.

          if (distSqr > radiusSqr) continue;
          if (!filter.passFilter(neighbourRef, neighbourTile, neighbourPoly)) continue;
          var neighbourNode = this.m_nodePool.getNode(neighbourRef);
          if ((neighbourNode.flags & _Node["default"].DT_NODE_CLOSED) != 0) continue; // Cost

          if (neighbourNode.flags == 0) {
            neighbourNode.pos = this.getEdgeMidPoint6(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile);
          }

          total = bestNode.total + _DetourCommon["default"].vDist2(bestNode.pos, neighbourNode.pos); // The node is already in open list and the new result is worse, skip.

          if ((neighbourNode.flags & _Node["default"].DT_NODE_OPEN) != 0 && total >= neighbourNode.total) continue;
          neighbourNode.id = neighbourRef;
          neighbourNode.flags = neighbourNode.flags & ~_Node["default"].DT_NODE_CLOSED;
          neighbourNode.pidx = this.m_nodePool.getNodeIdx(bestNode);
          neighbourNode.total = total;

          if ((neighbourNode.flags & _Node["default"].DT_NODE_OPEN) != 0) {
            this.m_openList.modify(neighbourNode);
          } else {
            neighbourNode.flags |= _Node["default"].DT_NODE_OPEN;
            this.m_openList.push(neighbourNode);
          }
        }
      } // Calc hit normal.


      var hitNormal = new Array(3);

      if (bestvi != null && bestvj != null) {
        var tangent = _DetourCommon["default"].vSub(bestvi, bestvj);

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

_defineProperty(NavMeshQuery, "SegInternval", (_temp = function SegInterval(ref, tmin, tmax) {
  _classCallCheck(this, SegInterval);

  _defineProperty(this, "ref", void 0);

  _defineProperty(this, "tmin", void 0);

  _defineProperty(this, "tmax", void 0);

  this.ref = ref;
  this.tmin = tmin;
  this.tmax = tmax;
}, _temp));

var _default = NavMeshQuery;
exports["default"] = _default;