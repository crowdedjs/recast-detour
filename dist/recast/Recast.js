"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _RecastVectors = _interopRequireDefault(require("./RecastVectors.js"));

var _CompactHeightfield = _interopRequireDefault(require("./CompactHeightfield.js"));

var _RecastConstants = _interopRequireDefault(require("./RecastConstants.js"));

var _CompactCell = _interopRequireDefault(require("./CompactCell.js"));

var _CompactSpan = _interopRequireDefault(require("./CompactSpan.js"));

var _RecastCommon = _interopRequireDefault(require("./RecastCommon.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

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

        if (norm[1] <= walkableThr) areas[i] = _RecastConstants["default"].RC_NULL_AREA;
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

      _RecastVectors["default"].subA(e0, verts, v1 * 3, v0 * 3);

      _RecastVectors["default"].subA(e1, verts, v2 * 3, v0 * 3);

      _RecastVectors["default"].cross(norm, e0, e1);

      _RecastVectors["default"].normalize(norm);
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
            if (s.area != _RecastConstants["default"].RC_NULL_AREA) spanCount++;
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
      var chf = new _CompactHeightfield["default"]();
      var w = hf.width;
      var h = hf.height;
      var spanCount = this.getHeightFieldSpanCount(ctx, hf); // Fill in header.

      chf.width = w;
      chf.height = h;
      chf.spanCount = spanCount;
      chf.walkableHeight = walkableHeight;
      chf.walkableClimb = walkableClimb;
      chf.maxRegions = 0;

      _RecastVectors["default"].copy2(chf.bmin, hf.bmin);

      _RecastVectors["default"].copy2(chf.bmax, hf.bmax);

      chf.bmax[1] += walkableHeight * hf.ch;
      chf.cs = hf.cs;
      chf.ch = hf.ch; // let bigSize = w*h;

      chf.cells = new Array(Math.floor(w) * Math.floor(h));
      chf.spans = new Array(spanCount);
      chf.areas = new Array(spanCount);
      var MAX_HEIGHT = 0xffff;

      for (var i = 0; i < chf.cells.length; i++) {
        chf.cells[i] = new _CompactCell["default"]();
      }

      for (var _i2 = 0; _i2 < chf.spans.length; _i2++) {
        chf.spans[_i2] = new _CompactSpan["default"]();
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
            if (s.area != _RecastConstants["default"].RC_NULL_AREA) {
              var bot = s.smax;
              var top = s.next != null ? Math.floor(s.next.smin) : MAX_HEIGHT;
              chf.spans[idx].y = _RecastCommon["default"].clamp(bot, 0, 0xffff);
              chf.spans[idx].h = _RecastCommon["default"].clamp(top - bot, 0, 0xff);
              chf.areas[idx] = s.area;
              idx++; // if(idx == 450240)
              // 	console.log("here")

              c.count++;
            }

            s = s.next;
          }
        }
      } // Find neighbour connections.


      var MAX_LAYERS = _RecastConstants["default"].RC_NOT_CONNECTED - 1;
      var tooHighNeighbour = 0;

      for (var y = 0; y < h; ++y) {
        for (var x = 0; x < w; ++x) {
          var _c = chf.cells[x + y * w];

          for (var _i3 = _c.index, ni = _c.index + _c.count; _i3 < ni; ++_i3) {
            // if(i == 450240)
            // 	console.log("stop")
            var _s = chf.spans[_i3];

            for (var dir = 0; dir < 4; ++dir) {
              _RecastCommon["default"].SetCon(_s, dir, _RecastConstants["default"].RC_NOT_CONNECTED);

              var nx = x + _RecastCommon["default"].GetDirOffsetX(dir);

              var ny = y + _RecastCommon["default"].GetDirOffsetY(dir); // First check that the neighbour cell is in bounds.


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

                  _RecastCommon["default"].SetCon(_s, dir, lidx);

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

var _default = Recast;
exports["default"] = _default;