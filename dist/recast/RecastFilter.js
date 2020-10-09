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
          var previousArea = _RecastConstants["default"].RC_NULL_AREA;

          for (var s = solid.spans[x + y * w]; s != null; ps = s, s = s.next) {
            var walkable = s.area != _RecastConstants["default"].RC_NULL_AREA; // If current span is not walkable, but there is walkable
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
            if (s.area == _RecastConstants["default"].RC_NULL_AREA) continue;
            var bot = s.smax;
            var top = s.next != null ? s.next.smin : MAX_HEIGHT; // Find neighbours minimum height.

            var minh = MAX_HEIGHT; // Min and max height of accessible neighbours.

            var asmin = s.smax;
            var asmax = s.smax;

            for (var dir = 0; dir < 4; ++dir) {
              var dx = x + _RecastCommon["default"].GetDirOffsetX(dir);

              var dy = y + _RecastCommon["default"].GetDirOffsetY(dir); // Skip neighbours which are out of bounds.


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


            if (minh < -walkableClimb) s.area = _RecastConstants["default"].RC_NULL_AREA; // If the difference between all neighbours is too large,
            // we are at steep slope, mark the span as ledge.

            if (asmax - asmin > walkableClimb) {
              s.area = _RecastConstants["default"].RC_NULL_AREA;
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
            if (top - bot <= walkableHeight) s.area = _RecastConstants["default"].RC_NULL_AREA;
          }
        }
      }

      ctx.stopTimer("FILTER_WALKABLE");
    }
  }]);

  return RecastFilter;
}();

var _default = RecastFilter;
exports["default"] = _default;