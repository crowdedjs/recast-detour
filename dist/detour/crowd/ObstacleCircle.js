"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var ObstacleCircle = function ObstacleCircle() {
  _classCallCheck(this, ObstacleCircle);

  _defineProperty(this, "p", new Array(3));

  _defineProperty(this, "vel", new Array(3));

  _defineProperty(this, "dvel", new Array(3));

  _defineProperty(this, "rad", void 0);

  _defineProperty(this, "dp", new Array(3));

  _defineProperty(this, "np", new Array(3));
};

var _default = ObstacleCircle;
exports["default"] = _default;