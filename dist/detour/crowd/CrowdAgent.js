"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _DetourCommon = _interopRequireDefault(require("../DetourCommon.js"));

var _NavMeshQuery = _interopRequireDefault(require("../NavMeshQuery.js"));

var _PathCorridor = _interopRequireDefault(require("./PathCorridor.js"));

var _LocalBoundary = _interopRequireDefault(require("./LocalBoundary.js"));

var _CrowdAgentAnimation = _interopRequireDefault(require("./CrowdAgentAnimation.js"));

var _PathQueue = _interopRequireDefault(require("./PathQueue.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

/// Represents an agent managed by a #dt object.
/// @ingroup crowd
var CrowdAgent = /*#__PURE__*/function () {
  /// The type of navigation mesh polygon the agent is currently traversing.
  /// @ingroup crowd
  /// True if the agent is active, false if the agent is in an unused slot in the agent pool.
  /// The type of mesh polygon the agent is traversing. (See: #CrowdAgent)
  /// True if the agent has valid path (targetState == DT_CROWDAGENT_TARGET_VALID) and the path does not lead to the requested position, else false.
  /// The path corridor the agent is using.
  /// The local boundary data for the agent.
  /// Time since the agent's path corridor was optimized.
  /// The known neighbors of the agent.
  /// The desired speed.
  ///< The current agent position. [(x, y, z)]
  ///< A temporary value used to accumulate agent displacement during iterative collision resolution. [(x, y, z)]
  ///< The desired velocity of the agent. Based on the current path, calculated from scratch each frame. [(x, y, z)]
  ///< The desired velocity adjusted by obstacle avoidance, calculated from scratch each frame. [(x, y, z)]
  ///< The actual velocity of the agent. The change from nvel => vel is constrained by max acceleration. [(x, y, z)]
  /// The agent's configuration parameters.
  /// The local path corridor corners for the agent.
  ///< State of the movement request.
  ///< Target polyref of the movement request.
  ///< Target position of the movement request (or velocity in case of DT_CROWDAGENT_TARGET_VELOCITY).
  ///< Path finder ref.
  ///< Flag indicating that the current path is being replanned.
  /// <Time since the agent's target was replanned.
  function CrowdAgent(idx) {
    _classCallCheck(this, CrowdAgent);

    _defineProperty(this, "idx", void 0);

    _defineProperty(this, "active", void 0);

    _defineProperty(this, "state", void 0);

    _defineProperty(this, "partial", void 0);

    _defineProperty(this, "corridor", void 0);

    _defineProperty(this, "boundary", void 0);

    _defineProperty(this, "topologyOptTime", void 0);

    _defineProperty(this, "neis", []);

    _defineProperty(this, "desiredSpeed", void 0);

    _defineProperty(this, "npos", new Array(3));

    _defineProperty(this, "disp", new Array(3));

    _defineProperty(this, "dvel", new Array(3));

    _defineProperty(this, "nvel", new Array(3));

    _defineProperty(this, "vel", new Array(3));

    _defineProperty(this, "params", void 0);

    _defineProperty(this, "corners", []);

    _defineProperty(this, "targetState", void 0);

    _defineProperty(this, "targetRef", void 0);

    _defineProperty(this, "targetPos", new Array(3));

    _defineProperty(this, "targetPathqReq", void 0);

    _defineProperty(this, "targetReplan", void 0);

    _defineProperty(this, "targetReplanTime", void 0);

    _defineProperty(this, "animation", void 0);

    this.idx = idx;
    this.corridor = new _PathCorridor["default"]();
    this.boundary = new _LocalBoundary["default"]();
    this.animation = new _CrowdAgentAnimation["default"]();
  }

  _createClass(CrowdAgent, [{
    key: "integrate",
    value: function integrate(dt) {
      // Fake dynamic constraint.
      var maxDelta = this.params.maxAcceleration * dt;

      var dv = _DetourCommon["default"].vSub(this.nvel, this.vel);

      var ds = _DetourCommon["default"].vLen(dv);

      if (ds > maxDelta) dv = _DetourCommon["default"].vScale(dv, maxDelta / ds);
      this.vel = _DetourCommon["default"].vAdd(this.vel, dv); // Integrate

      if (_DetourCommon["default"].vLen(this.vel) > 0.0001) this.npos = _DetourCommon["default"].vMad(this.npos, this.vel, dt);else _DetourCommon["default"].vSet(this.vel, 0, 0, 0);
    }
  }, {
    key: "overOffmeshConnection",
    value: function overOffmeshConnection(radius) {
      if (this.corners.length == 0) return false;
      var offMeshConnection = (this.corners[this.corners.length - 1].getFlags() & _NavMeshQuery["default"].DT_STRAIGHTPATH_OFFMESH_CONNECTION) != 0 ? true : false;

      if (offMeshConnection) {
        distSq = _DetourCommon["default"].vDist2D(this.npos, this.corners[this.corners.length - 1].getPos());
        if (distSq < radius * radius) return true;
      }

      return false;
    }
  }, {
    key: "getDistanceToGoal",
    value: function getDistanceToGoal(range) {
      if (this.corners.length == 0) return range;
      var endOfPath = (this.corners[this.corners.length - 1].getFlags() & _NavMeshQuery["default"].DT_STRAIGHTPATH_END) != 0 ? true : false;
      if (endOfPath) return Math.min(_DetourCommon["default"].vDist2D(this.npos, this.corners[this.corners.length - 1].getPos()), range);
      return range;
    }
  }, {
    key: "calcSmoothSteerDirection",
    value: function calcSmoothSteerDirection() {
      var dir = new Array(3);

      if (!this.corners.length == 0) {
        var ip0 = 0;
        var ip1 = Math.min(1, this.corners.length - 1);
        var p0 = this.corners[ip0].getPos();
        var p1 = this.corners[ip1].getPos();

        var dir0 = _DetourCommon["default"].vSub(p0, this.npos);

        var dir1 = _DetourCommon["default"].vSub(p1, this.npos);

        dir0[1] = 0;
        dir1[1] = 0;

        var len0 = _DetourCommon["default"].vLen(dir0);

        var len1 = _DetourCommon["default"].vLen(dir1);

        if (len1 > 0.001) dir1 = _DetourCommon["default"].vScale(dir1, 1.0 / len1);
        dir[0] = dir0[0] - dir1[0] * len0 * 0.5;
        dir[1] = 0;
        dir[2] = dir0[2] - dir1[2] * len0 * 0.5;

        _DetourCommon["default"].vNormalize(dir);
      }

      return dir;
    }
  }, {
    key: "calcStraightSteerDirection",
    value: function calcStraightSteerDirection() {
      var dir = new Array(3);

      if (!this.corners.length == 0) {
        dir = _DetourCommon["default"].vSub(this.corners[0].getPos(), this.npos);
        dir[1] = 0;

        _DetourCommon["default"].vNormalize(dir);
      }

      return dir;
    }
  }, {
    key: "setTarget",
    value: function setTarget(ref, pos) {
      this.targetRef = ref;

      _DetourCommon["default"].vCopy(this.targetPos, pos);

      this.targetPathqRef = _PathQueue["default"].DT_PATHQ_INVALID;
      if (this.targetRef != 0) this.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_REQUESTING;else this.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_FAILED;
    }
  }, {
    key: "getAgentIndex",
    value: function getAgentIndex() {
      return this.idx;
    }
  }, {
    key: "isActive",
    value: function isActive() {
      return this.active;
    }
  }]);

  return CrowdAgent;
}();

_defineProperty(CrowdAgent, "DT_CROWDAGENT_STATE_INVALID", 0);

_defineProperty(CrowdAgent, "DT_CROWDAGENT_STATE_WALKING", 1);

_defineProperty(CrowdAgent, "DT_CROWDAGENT_STATE_OFFMESH", 2);

_defineProperty(CrowdAgent, "DT_CROWDAGENT_TARGET_NONE", 0);

_defineProperty(CrowdAgent, "DT_CROWDAGENT_TARGET_FAILED", 1);

_defineProperty(CrowdAgent, "DT_CROWDAGENT_TARGET_VALID", 2);

_defineProperty(CrowdAgent, "DT_CROWDAGENT_TARGET_REQUESTING", 3);

_defineProperty(CrowdAgent, "DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE", 4);

_defineProperty(CrowdAgent, "DT_CROWDAGENT_TARGET_WAITING_FOR_PATH", 5);

_defineProperty(CrowdAgent, "DT_CROWDAGENT_TARGET_VELOCITY", 6);

var _default = CrowdAgent;
exports["default"] = _default;