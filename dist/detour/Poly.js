"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _NavMesh = _interopRequireDefault(require("./NavMesh.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

/** Defines a polyogn within a dtPoly object. */
var Poly = /*#__PURE__*/function () {
  /** The polygon is a standard convex polygon that is part of the surface of the mesh. */

  /** The polygon is an off-mesh connection consisting of two vertices. */

  /** Index to first link in linked list. (Or #DT_NULL_LINK if there is no link.) */

  /** The indices of the polygon's vertices. The actual vertices are located in MeshTile::verts. */

  /** Packed data representing neighbor polygons references and flags for each edge. */

  /** The user defined polygon flags. */

  /** The number of vertices in the polygon. */

  /**
   * The bit packed area id and polygon type.
   * 
   * @note Use the structure's set and get methods to access this value.
   */
  function Poly(index, maxVertsPerPoly) {
    _classCallCheck(this, Poly);

    _defineProperty(this, "index", void 0);

    _defineProperty(this, "firstLink", 0);

    _defineProperty(this, "verts", []);

    _defineProperty(this, "neis", []);

    _defineProperty(this, "flags", 0);

    _defineProperty(this, "vertCount", 0);

    _defineProperty(this, "areaAndtype", void 0);

    this.index = index;
    this.firstLink = _NavMesh["default"].DT_NULL_LINK;
    this.verts = new Array(maxVertsPerPoly);
    this.neis = new Array(maxVertsPerPoly);

    for (var i = 0; i < this.verts.length; i++) {
      this.verts[i] = 0;
    }

    for (var _i = 0; _i < this.neis.length; _i++) {
      this.neis[_i] = 0;
    }
  }

  _createClass(Poly, [{
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
    /** Sets the user defined area id. [Limit: < #DT_MAX_AREAS] */

  }, {
    key: "setArea",
    value: function setArea(a) {
      this.areaAndtype = this.or(this.and(this.areaAndtype, 0xc0), this.and(a, 0x3f));
    }
    /** Sets the polygon type. (See: #dtPolyTypes.) */

  }, {
    key: "setType",
    value: function setType(t) {
      this.areaAndtype = this.or(this.and(this.areaAndtype, 0x3f) | t << 6);
    }
    /** Gets the user defined area id. */

  }, {
    key: "getArea",
    value: function getArea() {
      return this.areaAndtype & 0x3;
    }
    /** Gets the polygon type. (See: #dtPolyTypes) */

  }, {
    key: "getType",
    value: function getType() {
      return this.areaAndtype >> 6;
    }
  }]);

  return Poly;
}();

_defineProperty(Poly, "DT_POLYTYPE_GROUND", 0);

_defineProperty(Poly, "DT_POLYTYPE_OFFMESH_CONNECTION", 1);

var _default = Poly;
exports["default"] = _default;