"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _NavMesh = _interopRequireDefault(require("./NavMesh.js"));

var _DetourCommon = _interopRequireDefault(require("./DetourCommon.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

// import QueryFilter from "./QueryFilter.js"
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
    var areaCost = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : new Array(_NavMesh["default"].DT_MAX_AREAS).fill(1, 0, _NavMesh["default"].DT_MAX_AREAS);

    _classCallCheck(this, DefaultQueryFilter);

    _defineProperty(this, "m_excludeFlags", 0);

    _defineProperty(this, "m_includeFlags", 0);

    _defineProperty(this, "m_areaCost", new Array(_NavMesh["default"].DT_MAX_AREAS));

    this.m_includeFlags = includeFlags;
    this.m_excludeFlags = excludeFlags;

    for (var i = 0; i < Math.min(_NavMesh["default"].DT_MAX_AREAS, areaCost.length); ++i) {
      this.m_areaCost[i] = areaCost[i];
    }

    for (var _i = areaCost.length; _i < _NavMesh["default"].DT_MAX_AREAS; ++_i) {
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
      return _DetourCommon["default"].vDist2(pa, pb) * this.m_areaCost[curPoly.getArea()];
    }
  }]);

  return DefaultQueryFilter;
}();

var _default = DefaultQueryFilter;
exports["default"] = _default;