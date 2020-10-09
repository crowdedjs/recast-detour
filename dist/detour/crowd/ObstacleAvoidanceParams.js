"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var ObstacleAvoidanceParams = ///< grid
///< adaptive
///< adaptive
///< adaptive
function ObstacleAvoidanceParams() {
  _classCallCheck(this, ObstacleAvoidanceParams);

  _defineProperty(this, "velBias", 0);

  _defineProperty(this, "weightDesVel", 0);

  _defineProperty(this, "weightCurVel", 0);

  _defineProperty(this, "weightSide", 0);

  _defineProperty(this, "weightToi", 0);

  _defineProperty(this, "horizTime", 0);

  _defineProperty(this, "gridSize", 0);

  _defineProperty(this, "adaptiveDivs", 0);

  _defineProperty(this, "adaptiveRings", 0);

  _defineProperty(this, "adaptiveDepth", 0);

  this.velBias = 0.4;
  this.weightDesVel = 2.0;
  this.weightCurVel = 0.75;
  this.weightSide = 0.75;
  this.weightToi = 2.5;
  this.horizTime = 2.5;
  this.gridSize = 33;
  this.adaptiveDivs = 7;
  this.adaptiveRings = 2;
  this.adaptiveDepth = 5;
};

;
var _default = ObstacleAvoidanceParams;
exports["default"] = _default;