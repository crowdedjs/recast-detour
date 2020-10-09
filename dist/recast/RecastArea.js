"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _RecastConstants = _interopRequireDefault(require("./RecastConstants.js"));

var _RecastCommon = _interopRequireDefault(require("./RecastCommon.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

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

            if (chf.areas[i] == _RecastConstants["default"].RC_NULL_AREA) {
              areas[i] = chf.areas[i];
              continue;
            }

            nei = new Array(9);

            for (var j = 0; j < 9; ++j) {
              nei[j] = chf.areas[i];
            }

            for (var dir = 0; dir < 4; ++dir) {
              if (_RecastCommon["default"].GetCon(s, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                var ax = x + _RecastCommon["default"].GetDirOffsetX(dir);

                var ay = y + _RecastCommon["default"].GetDirOffsetY(dir);

                var ai = chf.cells[ax + ay * w].index + _RecastCommon["default"].GetCon(s, dir);

                if (chf.areas[ai] != _RecastConstants["default"].RC_NULL_AREA) nei[dir * 2 + 0] = chf.areas[ai];
                var as = chf.spans[ai];
                var dir2 = dir + 1 & 0x3;

                if (_RecastCommon["default"].GetCon(as, dir2) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                  var ax2 = ax + _RecastCommon["default"].GetDirOffsetX(dir2);

                  var ay2 = ay + _RecastCommon["default"].GetDirOffsetY(dir2);

                  var ai2 = chf.cells[ax2 + ay2 * w].index + _RecastCommon["default"].GetCon(as, dir2);

                  if (chf.areas[ai2] != _RecastConstants["default"].RC_NULL_AREA) nei[dir * 2 + 1] = chf.areas[ai2];
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
              if (chf.areas[i] != _RecastConstants["default"].RC_NULL_AREA) chf.areas[i] = areaMod.apply(chf.areas[i]);
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
            if (chf.areas[i] == _RecastConstants["default"].RC_NULL_AREA) continue;

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
            if (chf.areas[i] == _RecastConstants["default"].RC_NULL_AREA) {
              dist[i] = 0;
            } else {
              var _s3 = chf.spans[i];
              var nc = 0;

              for (var dir = 0; dir < 4; ++dir) {
                if (_RecastCommon["default"].GetCon(_s3, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                  var nx = x + _RecastCommon["default"].GetDirOffsetX(dir);

                  var ny = y + _RecastCommon["default"].GetDirOffsetY(dir);

                  var nidx = chf.cells[nx + ny * w].index + _RecastCommon["default"].GetCon(_s3, dir);

                  if (chf.areas[nidx] != _RecastConstants["default"].RC_NULL_AREA) {
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

            if (_RecastCommon["default"].GetCon(_s4, 0) != _RecastConstants["default"].RC_NOT_CONNECTED) {
              // (-1,0)
              var ax = _x + _RecastCommon["default"].GetDirOffsetX(0);

              var ay = _y + _RecastCommon["default"].GetDirOffsetY(0);

              var ai = chf.cells[ax + ay * w].index + _RecastCommon["default"].GetCon(_s4, 0);

              var as = chf.spans[ai];
              nd = Math.min(dist[ai] + 2, 255);
              if (nd < dist[_i]) dist[_i] = nd; // (-1,-1)

              if (_RecastCommon["default"].GetCon(as, 3) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                var aax = ax + _RecastCommon["default"].GetDirOffsetX(3);

                var aay = ay + _RecastCommon["default"].GetDirOffsetY(3);

                var aai = chf.cells[aax + aay * w].index + _RecastCommon["default"].GetCon(as, 3);

                nd = Math.min(dist[aai] + 3, 255);
                if (nd < dist[_i]) dist[_i] = nd;
              }
            }

            if (_RecastCommon["default"].GetCon(_s4, 3) != _RecastConstants["default"].RC_NOT_CONNECTED) {
              // (0,-1)
              var _ax = _x + _RecastCommon["default"].GetDirOffsetX(3);

              var _ay = _y + _RecastCommon["default"].GetDirOffsetY(3);

              var _ai = chf.cells[_ax + _ay * w].index + _RecastCommon["default"].GetCon(_s4, 3);

              var _as = chf.spans[_ai];
              nd = Math.min(dist[_ai] + 2, 255);
              if (nd < dist[_i]) dist[_i] = nd; // (1,-1)

              if (_RecastCommon["default"].GetCon(_as, 2) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                var _aax = _ax + _RecastCommon["default"].GetDirOffsetX(2);

                var _aay = _ay + _RecastCommon["default"].GetDirOffsetY(2);

                var _aai = chf.cells[_aax + _aay * w].index + _RecastCommon["default"].GetCon(_as, 2);

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

            if (_RecastCommon["default"].GetCon(_s5, 2) != _RecastConstants["default"].RC_NOT_CONNECTED) {
              // (1,0)
              var _ax2 = _x2 + _RecastCommon["default"].GetDirOffsetX(2);

              var _ay2 = _y2 + _RecastCommon["default"].GetDirOffsetY(2);

              var _ai2 = chf.cells[_ax2 + _ay2 * w].index + _RecastCommon["default"].GetCon(_s5, 2);

              if (_ai2 == 450241) console.log("Here");
              var _as2 = chf.spans[_ai2];
              nd = Math.min(dist[_ai2] + 2, 255);
              if (nd < dist[_i2]) dist[_i2] = nd; // (1,1)

              if (_RecastCommon["default"].GetCon(_as2, 1) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                var _aax2 = _ax2 + _RecastCommon["default"].GetDirOffsetX(1);

                var _aay2 = _ay2 + _RecastCommon["default"].GetDirOffsetY(1);

                var _aai2 = chf.cells[_aax2 + _aay2 * w].index + _RecastCommon["default"].GetCon(_as2, 1);

                nd = Math.min(dist[_aai2] + 3, 255);
                if (nd < dist[_i2]) dist[_i2] = nd;
              }
            }

            if (_RecastCommon["default"].GetCon(_s5, 1) != _RecastConstants["default"].RC_NOT_CONNECTED) {
              // (0,1)
              var _ax3 = _x2 + _RecastCommon["default"].GetDirOffsetX(1);

              var _ay3 = _y2 + _RecastCommon["default"].GetDirOffsetY(1);

              var _ai3 = chf.cells[_ax3 + _ay3 * w].index + _RecastCommon["default"].GetCon(_s5, 1);

              var _as3 = chf.spans[_ai3];
              nd = Math.min(dist[_ai3] + 2, 255);
              if (nd < dist[_i2]) dist[_i2] = nd; // (-1,1)

              if (_RecastCommon["default"].GetCon(_as3, 0) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                var _aax3 = _ax3 + _RecastCommon["default"].GetDirOffsetX(0);

                var _aay3 = _ay3 + _RecastCommon["default"].GetDirOffsetY(0);

                var _aai3 = chf.cells[_aax3 + _aay3 * w].index + _RecastCommon["default"].GetCon(_as3, 0);

                nd = Math.min(dist[_aai3] + 3, 255);
                if (nd < dist[_i2]) dist[_i2] = nd;
              }
            }
          }
        }
      }

      var thr = radius * 2;

      for (var _i3 = 0; _i3 < chf.spanCount; ++_i3) {
        if (dist[_i3] < thr) chf.areas[_i3] = _RecastConstants["default"].RC_NULL_AREA;
      }

      ctx.stopTimer("ERODE_AREA");
    }
  }, {
    key: "pointInPoly",
    value: function pointInPoly(verts, p) {
      var c = false;
      var i, j;

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
            if (chf.areas[_i5] == _RecastConstants["default"].RC_NULL_AREA) continue;

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

var _default = RecastArea;
exports["default"] = _default;