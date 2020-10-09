"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _RecastVectors = _interopRequireDefault(require("./RecastVectors.js"));

var _RecastCommon = _interopRequireDefault(require("./RecastCommon.js"));

var _RecastConstants = _interopRequireDefault(require("./RecastConstants.js"));

var _Span = _interopRequireDefault(require("./Span.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

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
      var s = new _Span["default"]();
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

          _RecastVectors["default"].copy4(buf, out2 + n * 3, buf, out1 + m * 3);

          m++;
          n++; // add the i'th poPoly to the right polygon. Do NOT add points that are on the dividing line
          // since these were already added above

          if (d[_i] > 0) {
            _RecastVectors["default"].copy4(buf, out1 + m * 3, buf, _in + _i * 3);

            m++;
          } else if (d[_i] < 0) {
            _RecastVectors["default"].copy4(buf, out2 + n * 3, buf, _in + _i * 3);

            n++;
          }
        } else // same side
          {
            // add the i'th poPoly to the right polygon. Addition is done even for points on the dividing line
            if (d[_i] >= 0) {
              _RecastVectors["default"].copy4(buf, out1 + m * 3, buf, _in + _i * 3);

              m++;
              if (d[_i] != 0) continue;
            }

            _RecastVectors["default"].copy4(buf, out2 + n * 3, buf, _in + _i * 3);

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

      _RecastVectors["default"].copy3(tmin, verts, v0 * 3);

      _RecastVectors["default"].copy3(tmax, verts, v0 * 3);

      _RecastVectors["default"].min(tmin, verts, v1 * 3);

      _RecastVectors["default"].min(tmin, verts, v2 * 3);

      _RecastVectors["default"].max(tmax, verts, v1 * 3);

      _RecastVectors["default"].max(tmax, verts, v2 * 3); // If the triangle does not touch the bbox of the heightfield, skip the triagle.


      if (!RecastRasterization.overlapBounds(bmin, bmax, tmin, tmax)) return; // Calculate the footprPoly of the triangle on the grid's y-axis

      var y0 = Math.floor((tmin[2] - bmin[2]) * ics);
      var y1 = Math.floor((tmax[2] - bmin[2]) * ics);
      y0 = _RecastCommon["default"].clamp(y0, 0, h - 1);
      y1 = _RecastCommon["default"].clamp(y1, 0, h - 1); // Clip the triangle into all grid cells it touches.

      var buf = new Array(7 * 3 * 4).fill(0);
      var _in = 0;
      var inrow = 7 * 3;
      var p1 = inrow + 7 * 3;
      var p2 = p1 + 7 * 3;

      _RecastVectors["default"].copy4(buf, 0, verts, v0 * 3);

      _RecastVectors["default"].copy4(buf, 3, verts, v1 * 3);

      _RecastVectors["default"].copy4(buf, 6, verts, v2 * 3);

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
        x0 = _RecastCommon["default"].clamp(x0, 0, w - 1);
        x1 = _RecastCommon["default"].clamp(x1, 0, w - 1);
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

          var ismin = _RecastCommon["default"].clamp(Math.floor(smin * ich), 0, _RecastConstants["default"].RC_SPAN_MAX_HEIGHT);

          var ismax = _RecastCommon["default"].clamp(Math.ceil(smax * ich), ismin + 1, _RecastConstants["default"].RC_SPAN_MAX_HEIGHT);

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

var _default = RecastRasterization;
exports["default"] = _default;