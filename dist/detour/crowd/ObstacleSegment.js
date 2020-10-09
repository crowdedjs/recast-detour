"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var ObstacleSegment = function ObstacleSegment() {
  _classCallCheck(this, ObstacleSegment);

  _defineProperty(this, "p", new Array(3));

  _defineProperty(this, "q", new Array(3));

  _defineProperty(this, "touch", false);
};

var _default = ObstacleSegment;
exports["default"] = _default;