"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = void 0;

var _temp;

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

/// Represents a set of heightfield layers.
/// @ingroup recast
/// @see rcAllocHeightfieldLayerSet, rcFreeHeightfieldLayerSet 
var HeightfieldLayerSet = function HeightfieldLayerSet() {
  _classCallCheck(this, HeightfieldLayerSet);

  _defineProperty(this, "layers", void 0);
} ///< The layers in the set. [Size: #nlayers]
;

_defineProperty(HeightfieldLayerSet, "HeightfieldLayer", (_temp = function HeightfieldLayer() {
  _classCallCheck(this, HeightfieldLayer);

  _defineProperty(this, "bmin", new Array(3));

  _defineProperty(this, "bmax", new Array(3));

  _defineProperty(this, "cs", void 0);

  _defineProperty(this, "ch", void 0);

  _defineProperty(this, "width", void 0);

  _defineProperty(this, "height", void 0);

  _defineProperty(this, "minx", void 0);

  _defineProperty(this, "maxx", void 0);

  _defineProperty(this, "miny", void 0);

  _defineProperty(this, "maxy", void 0);

  _defineProperty(this, "hmin", void 0);

  _defineProperty(this, "hmax", void 0);

  _defineProperty(this, "heights", void 0);

  _defineProperty(this, "areas", void 0);

  _defineProperty(this, "cons", void 0);
} ///< Packed neighbor connection information. [Size: Same as #heights]
, _temp));

var _default = HeightfieldLayerSet;
exports["default"] = _default;