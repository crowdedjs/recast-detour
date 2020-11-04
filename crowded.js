(function (factory) {
  typeof define === 'function' && define.amd ? define(factory) :
  factory();
}((function () { 'use strict';

  Object.defineProperty(exports, "__esModule", {
    value: true
  });
  Object.defineProperty(exports, "CrowdAgentParams", {
    enumerable: true,
    get: function get() {
      return _CrowdAgentParams["default"];
    }
  });
  Object.defineProperty(exports, "RecastTestMeshBuilder", {
    enumerable: true,
    get: function get() {
      return _RecastTestMeshBuilder["default"];
    }
  });
  Object.defineProperty(exports, "NavMesh", {
    enumerable: true,
    get: function get() {
      return _NavMesh["default"];
    }
  });
  Object.defineProperty(exports, "NavMeshQuery", {
    enumerable: true,
    get: function get() {
      return _NavMeshQuery["default"];
    }
  });
  Object.defineProperty(exports, "Crowd", {
    enumerable: true,
    get: function get() {
      return _Crowd["default"];
    }
  });
  Object.defineProperty(exports, "ObstacleAvoidanceParams", {
    enumerable: true,
    get: function get() {
      return _ObstacleAvoidanceParams["default"];
    }
  });

  var _CrowdAgentParams = _interopRequireDefault(require("./detour/crowd/CrowdAgentParams.js"));

  var _RecastTestMeshBuilder = _interopRequireDefault(require("./detour/RecastTestMeshBuilder.js"));

  var _NavMesh = _interopRequireDefault(require("./detour/NavMesh.js"));

  var _NavMeshQuery = _interopRequireDefault(require("./detour/NavMeshQuery.js"));

  var _Crowd = _interopRequireDefault(require("./detour/crowd/Crowd.js"));

  var _ObstacleAvoidanceParams = _interopRequireDefault(require("./detour/crowd/ObstacleAvoidanceParams.js"));

  function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

})));
