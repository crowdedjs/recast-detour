"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var IntersectResult = function IntersectResult() {
  _classCallCheck(this, IntersectResult);

  _defineProperty(this, "intersects", false);

  _defineProperty(this, "tmin", 0);

  _defineProperty(this, "tmax", 1);

  _defineProperty(this, "segMin", -1);

  _defineProperty(this, "segMax", -1);
};

var _default = IntersectResult;
exports["default"] = _default;