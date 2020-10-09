"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

/*
Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
Recast4J Copyright (c) 2015 Piotr Piastucki piotr@jtilia.org

This software is provided 'as-is', without any express or implied
warranty.  In no event will the authors be held liable for any damages
arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:
1. The origin of this software must not be misrepresented; you must not
 claim that you wrote the original software. If you use this software
 in a product, an acknowledgment in the product documentation would be
 appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
 misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
/// Configuration parameters for a crowd agent.
/// @ingroup crowd
var CrowdAgentParams = function CrowdAgentParams() {
  _classCallCheck(this, CrowdAgentParams);

  _defineProperty(this, "radius", 0);

  _defineProperty(this, "height", 0);

  _defineProperty(this, "maxAcceleration", 0);

  _defineProperty(this, "maxSpeed", 0);

  _defineProperty(this, "collisionQueryRange", 0);

  _defineProperty(this, "pathOptimizationRange", 0);

  _defineProperty(this, "separationWeight", 0);

  _defineProperty(this, "updateFlags", 0);

  _defineProperty(this, "obstacleAvoidanceType", 0);

  _defineProperty(this, "queryFilterType", 0);

  _defineProperty(this, "userData", null);
};

_defineProperty(CrowdAgentParams, "DT_CROWD_ANTICIPATE_TURNS", 1);

_defineProperty(CrowdAgentParams, "DT_CROWD_OBSTACLE_AVOIDANCE", 2);

_defineProperty(CrowdAgentParams, "DT_CROWD_SEPARATION", 4);

_defineProperty(CrowdAgentParams, "DT_CROWD_OPTIMIZE_VIS", 8);

_defineProperty(CrowdAgentParams, "DT_CROWD_OPTIMIZE_TOPO", 16);

var _default = CrowdAgentParams;
exports["default"] = _default;