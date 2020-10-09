"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var ConvexVolume = function ConvexVolume() {
  _classCallCheck(this, ConvexVolume);

  _defineProperty(this, "verts", void 0);

  _defineProperty(this, "hmin", void 0);

  _defineProperty(this, "hmax", void 0);

  _defineProperty(this, "areaMod", void 0);
};

var _default = ConvexVolume;
exports["default"] = _default;