"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _DetourCommon = _interopRequireDefault(require("../DetourCommon.js"));

var _temp;

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _createForOfIteratorHelper(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

function arraycopy(one, oneStart, two, twoStart, len) {
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
      arraycopy(s, 0, seg.s, 0, 6);
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

      _DetourCommon["default"].vCopy(this.m_center, pos); // First query non-overlapping polygons.


      var res = navquery.findLocalNeighbourhood(ref, pos, collisionQueryRange, filter);
      this.m_polys = res.getRefs();
      this.m_segs = []; // Secondly, store all polygon edges.

      for (var j = 0; j < this.m_polys.length; ++j) {
        var gpws = navquery.getPolyWallSegments(this.m_polys[j], false, filter);

        for (var k = 0; k < gpws.getSegmentRefs().length; ++k) {
          var s = gpws.getSegmentVerts()[k]; // Skip too distant segments.

          var distseg = _DetourCommon["default"].distancePtSegSqr2D4(pos, s, 0, 3);

          if (distseg[0] > _DetourCommon["default"].sqr(collisionQueryRange)) continue;
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

_defineProperty(LocalBoundary, "Segment", (_temp = function Segment() {
  _classCallCheck(this, Segment);

  _defineProperty(this, "s", new Array(6).fill(0));

  _defineProperty(this, "d", void 0);
}, _temp));

var _default = LocalBoundary;
exports["default"] = _default;