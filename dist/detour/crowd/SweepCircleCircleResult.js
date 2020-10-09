"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var SweepCircleCircleResult = function SweepCircleCircleResult(intersection, htmin, htmax) {
  _classCallCheck(this, SweepCircleCircleResult);

  _defineProperty(this, "intersection", false);

  _defineProperty(this, "htmin", 0);

  _defineProperty(this, "htmax", 0);

  this.intersection = intersection;
  this.htmin = htmin;
  this.htmax = htmax;
};

var _default = SweepCircleCircleResult;
exports["default"] = _default;