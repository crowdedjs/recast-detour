"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _ObstacleAvoidanceQuery = _interopRequireDefault(require("./ObstacleAvoidanceQuery.js"));

var _DetourCommon = _interopRequireDefault(require("../DetourCommon.js"));

var _ProximityGrid = _interopRequireDefault(require("./ProximityGrid.js"));

var _PathQueue = _interopRequireDefault(require("./PathQueue.js"));

var _CrowdAgent = _interopRequireDefault(require("./CrowdAgent.js"));

var _NavMeshQuery = _interopRequireDefault(require("../NavMeshQuery.js"));

var _DefaultQueryFilter = _interopRequireDefault(require("../DefaultQueryFilter.js"));

var _PriorityQueue = _interopRequireDefault(require("./PriorityQueue.js"));

var _Status = _interopRequireDefault(require("../Status.js"));

var _ObstacleAvoidanceParams = _interopRequireDefault(require("./ObstacleAvoidanceParams.js"));

var _CrowdAgentParams = _interopRequireDefault(require("./CrowdAgentParams.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _toConsumableArray(arr) { return _arrayWithoutHoles(arr) || _iterableToArray(arr) || _unsupportedIterableToArray(arr) || _nonIterableSpread(); }

function _nonIterableSpread() { throw new TypeError("Invalid attempt to spread non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); }

function _iterableToArray(iter) { if (typeof Symbol !== "undefined" && Symbol.iterator in Object(iter)) return Array.from(iter); }

function _arrayWithoutHoles(arr) { if (Array.isArray(arr)) return _arrayLikeToArray(arr); }

function _createForOfIteratorHelper(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var Crowd = /*#__PURE__*/function () {
  _createClass(Crowd, [{
    key: "tween",
    /// The maximum number of neighbors that a crowd agent can take into account
    /// for steering decisions.
    /// @ingroup crowd
    /// The maximum number of corners a crowd agent will look ahead in the path.
    /// This value is used for sizing the crowd agent corner buffers.
    /// Due to the behavior of the crowd manager, the actual number of useful
    /// corners will be one less than this number.
    /// @ingroup crowd
    /// The maximum number of crowd avoidance configurations supported by the
    /// crowd manager.
    /// @ingroup crowd
    /// @see dtObstacleAvoidanceParams, dtCrowd::setObstacleAvoidanceParams(), dtCrowd::getObstacleAvoidanceParams(),
    /// dtCrowdAgentParams::obstacleAvoidanceType
    /// The maximum number of query filter types supported by the crowd manager.
    /// @ingroup crowd
    /// @see dtQueryFilter, dtCrowd::getFilter() dtCrowd::getEditableFilter(),
    /// dtCrowdAgentParams::queryFilterType
    /// Provides neighbor data for agents managed by the crowd.
    /// @ingroup crowd
    /// @see dtCrowdAgent::neis, dtCrowd
    value: function tween(t, t0, t1) {
      return _DetourCommon["default"].clamp((t - t0) / (t1 - t0), 0.0, 1.0);
    }
  }, {
    key: "getNeighbours",
    value: function getNeighbours(pos, height, range, skip, agents, grid) {
      var result = [];
      var ids = grid.queryItems(pos[0] - range, pos[2] - range, pos[0] + range, pos[2] + range);

      var _iterator = _createForOfIteratorHelper(ids),
          _step;

      try {
        for (_iterator.s(); !(_step = _iterator.n()).done;) {
          var id = _step.value;
          var _ag = agents[id];

          if (_ag == skip || !_ag.active) {
            // console.log("Testing");
            continue;
          } // Check for overlap.


          var diff = _DetourCommon["default"].vSub(pos, _ag.npos);

          if (Math.abs(diff[1]) >= (height + _ag.params.height) / 2.0) {
            continue;
          }

          diff[1] = 0;

          var _distSqr = _DetourCommon["default"].vLenSqr(diff);

          if (_distSqr > _DetourCommon["default"].sqr(range)) {
            continue;
          }

          this.addNeighbour(id, _distSqr, result);
        }
      } catch (err) {
        _iterator.e(err);
      } finally {
        _iterator.f();
      }

      return result;
    }
  }, {
    key: "addNeighbour",
    value: function addNeighbour(idx, dist, neis) {
      // Insert neighbour based on the distance.
      var nei = new this.CrowdNeighbor(idx, dist);
      neis.push(nei);
      neis.sort(function (o1, o2) {
        return o1.dist - o2.dist;
      });
    }
  }, {
    key: "addToOptQueue",
    value: function addToOptQueue(newag, agents) {
      // Insert neighbour based on greatest time.
      agents.push(newag);
    } // Insert neighbour based on greatest time.

  }, {
    key: "addToPathQueue",
    value: function addToPathQueue(newag, agents) {
      agents.push(newag);
    } ///
    /// Initializes the crowd.
    /// May be called more than once to purge and re-initialize the crowd.
    /// @param[in] maxAgents The maximum number of agents the crowd can manage. [Limit: >= 1]
    /// @param[in] maxAgentRadius The maximum radius of any agent that will be added to the crowd. [Limit: > 0]
    /// @param[in] nav The navigation mesh to use for planning.
    /// @return True if the initialization succeeded.

  }]);

  function Crowd(maxAgents, maxAgentRadius, nav) {
    var _temp;

    var queryFilterFactory = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : function (i) {
      return new _DefaultQueryFilter["default"]();
    };

    _classCallCheck(this, Crowd);

    _defineProperty(this, "CrowdNeighbor", (_temp = /// < The index of the neighbor in the crowd.
    /// < The distance between the current agent and the neighbor.
    function CrowdNeighbour(idx, dist) {
      _classCallCheck(this, CrowdNeighbour);

      _defineProperty(this, "idx", void 0);

      _defineProperty(this, "dist", void 0);

      this.idx = idx;
      this.dist = dist;
    }, _temp));

    _defineProperty(this, "m_maxAgents", void 0);

    _defineProperty(this, "m_agents", []);

    _defineProperty(this, "m_activeAgents", []);

    _defineProperty(this, "m_pathq", void 0);

    _defineProperty(this, "m_obstacleQueryParams", new Array(Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS));

    _defineProperty(this, "m_obstacleQuery", void 0);

    _defineProperty(this, "m_grid", void 0);

    _defineProperty(this, "m_ext", new Array(3));

    _defineProperty(this, "m_filters", new Array(Crowd.DT_CROWD_MAX_QUERY_FILTER_TYPE));

    _defineProperty(this, "m_maxAgentRadius", void 0);

    _defineProperty(this, "m_velocitySampleCount", void 0);

    _defineProperty(this, "m_navquery", void 0);

    this.m_maxAgents = maxAgents;
    this.m_maxAgentRadius = maxAgentRadius;

    _DetourCommon["default"].vSet(this.m_ext, this.m_maxAgentRadius * 2.0, this.m_maxAgentRadius * 1.5, this.m_maxAgentRadius * 2.0);

    this.m_grid = new _ProximityGrid["default"](this.m_maxAgents * 4, this.m_maxAgentRadius * 3);
    this.m_obstacleQuery = new _ObstacleAvoidanceQuery["default"](6, 8);

    for (var i = 0; i < Crowd.DT_CROWD_MAX_QUERY_FILTER_TYPE; i++) {
      this.m_filters[i] = queryFilterFactory.apply(i);
    } // Init obstacle query params.


    for (var _i = 0; _i < Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS; ++_i) {
      this.m_obstacleQueryParams[_i] = new _ObstacleAvoidanceParams["default"]();
    } // Allocate temp buffer for merging paths.


    this.m_pathq = new _PathQueue["default"](Crowd.MAX_PATHQUEUE_NODES, nav);
    this.m_agents = new Array(this.m_maxAgents);
    this.m_activeAgents = [];

    for (var _i2 = 0; _i2 < this.m_maxAgents; ++_i2) {
      this.m_agents[_i2] = new _CrowdAgent["default"](_i2);
      this.m_agents[_i2].active = false;
    } // The navquery is mostly used for local searches, no need for large
    // node pool.


    this.m_navquery = new _NavMeshQuery["default"](nav);
  } /// Sets the shared avoidance configuration for the specified index.
  /// @param[in] idx The index. [Limits: 0 <= value <
  /// #DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS]
  /// @param[in] params The new configuration.


  _createClass(Crowd, [{
    key: "setObstacleAvoidanceParams",
    value: function setObstacleAvoidanceParams(idx, params) {
      if (idx >= 0 && idx < Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS) {
        this.m_obstacleQueryParams[idx] = params;
      }
    } /// Gets the shared avoidance configuration for the specified index.
    /// @param[in] idx The index of the configuration to retreive.
    /// [Limits: 0 <= value < #DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS]
    /// @return The requested configuration.

  }, {
    key: "getObstacleAvoidanceParams",
    value: function getObstacleAvoidanceParams(idx) {
      if (idx >= 0 && idx < Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS) {
        return this.m_obstacleQueryParams[idx];
      }

      return null;
    } /// The maximum number of agents that can be managed by the object.
    /// @return The maximum number of agents.

  }, {
    key: "getAgentCount",
    value: function getAgentCount() {
      return this.m_maxAgents;
    } /// Gets the specified agent from the pool.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
    /// @return The requested agent.
    /// Agents in the pool may not be in use. Check #dtCrowdAgent.active before using the returned object.

  }, {
    key: "getAgent",
    value: function getAgent(idx) {
      return idx < 0 || idx >= this.m_agents.length ? null : this.m_agents[idx];
    } ///
    /// Gets the specified agent from the pool.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
    /// @return The requested agent.
    /// Agents in the pool may not be in use. Check #dtCrowdAgent.active before using the returned object.

  }, {
    key: "getEditableAgent",
    value: function getEditableAgent(idx) {
      return idx < 0 || idx >= this.m_agents.length ? null : this.m_agents[idx];
    } /// Updates the specified agent's configuration.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
    /// @param[in] params The new agent configuration.

  }, {
    key: "updateAgentParameters",
    value: function updateAgentParameters(idx, params) {
      if (idx < 0 || idx >= this.m_maxAgents) {
        return;
      }

      this.m_agents[idx].params = params;
    } /// Adds a new agent to the crowd.
    /// @param[in] pos The requested position of the agent. [(x, y, z)]
    /// @param[in] params The configutation of the agent.
    /// @return The index of the agent in the agent pool. Or -1 if the agent
    /// could not be added.

  }, {
    key: "addAgent",
    value: function addAgent(pos, params) {
      // Find empty slot.
      var idx = -1;

      for (var i = 0; i < this.m_maxAgents; ++i) {
        if (!this.m_agents[i].active) {
          idx = i;
          break;
        }
      }

      if (idx == -1) {
        return -1;
      }

      var ag = this.m_agents[idx];
      this.updateAgentParameters(idx, params); // Find nearest position on navmesh and place the agent there.

      var nearest = this.m_navquery.findNearestPoly(pos, this.m_ext, this.m_filters[ag.params.queryFilterType]);
      ag.corridor.reset(nearest.getNearestRef(), nearest.getNearestPos());
      ag.boundary.reset();
      ag.partial = false;
      ag.topologyOptTime = 0;
      ag.targetReplanTime = 0;

      _DetourCommon["default"].vSet(ag.dvel, 0, 0, 0);

      _DetourCommon["default"].vSet(ag.nvel, 0, 0, 0);

      _DetourCommon["default"].vSet(ag.vel, 0, 0, 0);

      _DetourCommon["default"].vCopy(ag.npos, nearest.getNearestPos());

      ag.desiredSpeed = 0;

      if (nearest.getNearestRef() != 0) {
        ag.state = _CrowdAgent["default"].DT_CROWDAGENT_STATE_WALKING;
      } else {
        ag.state = _CrowdAgent["default"].DT_CROWDAGENT_STATE_INVALID;
      }

      ag.targetState = _CrowdAgent["default"].DT_CROWDAGENT_TARGET_NONE;
      ag.active = true;
      return idx;
    } /// Removes the agent from the crowd.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
    ///
    /// The agent is deactivated and will no longer be processed. Its
    /// #dt object
    /// is not removed from the pool. It is marked as inactive so that it is
    /// available for reuse.
    /// Removes the agent from the crowd.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]

  }, {
    key: "removeAgent",
    value: function removeAgent(idx) {
      if (idx >= 0 && idx < this.m_maxAgents) {
        this.m_agents[idx].active = false;
      }
    }
  }, {
    key: "requestMoveTargetReplan",
    value: function requestMoveTargetReplan(ag, ref, pos) {
      ag.setTarget(ref, pos);
      ag.targetReplan = true;
      return true;
    } /// Submits a new move request for the specified agent.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
    /// @param[in] ref The position's polygon reference.
    /// @param[in] pos The position within the polygon. [(x, y, z)]
    /// @return True if the request was successfully submitted.
    ///
    /// This method is used when a new target is set.
    ///
    /// The position will be constrained to the surface of the navigation mesh.
    ///
    /// The request will be processed during the next #update().

  }, {
    key: "requestMoveTarget",
    value: function requestMoveTarget(idx, ref, pos) {
      if (idx < 0 || idx >= this.m_maxAgents) {
        return false;
      }

      if (ref == 0) {
        return false;
      }

      var ag = this.m_agents[idx]; // Initialize request.

      ag.setTarget(ref, pos);
      ag.targetReplan = false;
      return true;
    } /// Submits a new move request for the specified agent.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
    /// @param[in] vel The movement velocity. [(x, y, z)]
    /// @return True if the request was successfully submitted.

  }, {
    key: "requestMoveVelocity",
    value: function requestMoveVelocity(idx, vel) {
      if (idx < 0 || idx >= this.m_maxAgents) {
        return false;
      }

      ag = this.m_agents[idx]; // Initialize request.

      ag.targetRef = 0;

      _DetourCommon["default"].vCopy(ag.targetPos, vel);

      ag.targetPathqRef = _PathQueue["default"].DT_PATHQ_INVALID;
      ag.targetReplan = false;
      ag.targetState = _CrowdAgent["default"].DT_CROWDAGENT_TARGET_VELOCITY;
      return true;
    } /// Resets any request for the specified agent.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
    /// @return True if the request was successfully reseted.

  }, {
    key: "resetMoveTarget",
    value: function resetMoveTarget(idx) {
      if (idx < 0 || idx >= this.m_maxAgents) {
        return false;
      }

      ag = this.m_agents[idx]; // Initialize request.

      ag.targetRef = 0;

      _DetourCommon["default"].vSet(ag.targetPos, 0, 0, 0);

      _DetourCommon["default"].vSet(ag.dvel, 0, 0, 0);

      ag.targetPathqRef = _PathQueue["default"].DT_PATHQ_INVALID;
      ag.targetReplan = false;
      ag.targetState = _CrowdAgent["default"].DT_CROWDAGENT_TARGET_NONE;
      return true;
    } /// Gets the active agents let the agent pool.
    /// @param[out] agents An array of agent pointers. [(#dt *) * maxAgents]
    /// @param[in] maxAgents The size of the crowd agent array.
    /// @return The number of agents returned in @p agents.

  }, {
    key: "getActiveAgents",
    value: function getActiveAgents() {
      var agents = []; //new Array(this.m_maxAgents);

      for (var i = 0; i < this.m_maxAgents; ++i) {
        if (this.m_agents[i].active) {
          agents.push(this.m_agents[i]);
        }
      }

      return agents;
    }
  }, {
    key: "updateMoveRequest",
    value: function updateMoveRequest() {
      var queue = new _PriorityQueue["default"](function (a1, a2) {
        return a2.targetReplanTime - a1.targetReplanTime;
      }); // Fire off new requests.

      for (var i = 0; i < this.m_maxAgents; ++i) {
        // if(i==12)
        //     console.log("Bad agent.")
        var _ag2 = this.m_agents[i];

        if (!_ag2.active) {
          continue;
        }

        if (_ag2.state == _CrowdAgent["default"].DT_CROWDAGENT_STATE_INVALID) {
          continue;
        }

        if (_ag2.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_NONE || _ag2.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_VELOCITY) {
          continue;
        }

        if (_ag2.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_REQUESTING) {
          var path = _ag2.corridor.getPath();

          if (path.length == 0) {
            throw new IllegalArgumentException("Empty path");
          } // Quick search towards the goal.


          this.m_navquery.initSlicedFindPath(path[0], _ag2.targetRef, _ag2.npos, _ag2.targetPos, this.m_filters[_ag2.params.queryFilterType], 0);
          this.m_navquery.updateSlicedFindPath(Crowd.MAX_ITER);
          var pathFound = void 0;

          if (_ag2.targetReplan) // && npath > 10)
            {
              // Try to use existing steady path during replan if
              // possible.
              pathFound = this.m_navquery.finalizeSlicedFindPathPartial(path);
            } else {
            // Try to move towards target when goal changes.
            pathFound = this.m_navquery.finalizeSlicedFindPath();
          }

          var reqPath = pathFound.getRefs();
          var reqPos = new Array(3);

          if (!pathFound.getStatus() == _Status["default"].FAILURE && reqPath.length > 0) {
            // In progress or succeed.
            if (reqPath[reqPath.length - 1] != _ag2.targetRef) {
              // Partial path, constrain target position inside the
              // last polygon.
              var cr = this.m_navquery.closestPointOnPoly(reqPath[reqPath.length - 1], _ag2.targetPos);
              reqPos = cr.getClosest();
            } else {
              _DetourCommon["default"].vCopy(reqPos, _ag2.targetPos);
            }
          } else {
            // Could not find path, start the request from current
            // location.
            _DetourCommon["default"].vCopy(reqPos, _ag2.npos);

            reqPath = [];
            reqPath.push(path[0]);
          }

          _ag2.corridor.setCorridor(reqPos, reqPath);

          _ag2.boundary.reset();

          _ag2.partial = false;

          if (reqPath[reqPath.length - 1] == _ag2.targetRef) {
            _ag2.targetState = _CrowdAgent["default"].DT_CROWDAGENT_TARGET_VALID;
            _ag2.targetReplanTime = 0.0;
          } else {
            // The path is longer or potentially unreachable, full plan.
            _ag2.targetState = _CrowdAgent["default"].DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE;
          }
        }

        if (_ag2.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE) {
          this.addToPathQueue(_ag2, queue);
        }
      }

      while (!queue.isEmpty()) {
        var _ag3 = queue.poll();

        _ag3.targetPathqRef = this.m_pathq.request(_ag3.corridor.getLastPoly(), _ag3.targetRef, _ag3.corridor.getTarget(), _ag3.targetPos, this.m_filters[_ag3.params.queryFilterType]);

        if (_ag3.targetPathqRef != _PathQueue["default"].DT_PATHQ_INVALID) {
          _ag3.targetState = _CrowdAgent["default"].DT_CROWDAGENT_TARGET_WAITING_FOR_PATH;
        }
      } // Update requests.


      this.m_pathq.update(Crowd.MAX_ITERS_PER_UPDATE); // Process path results.

      for (var _i3 = 0; _i3 < this.m_maxAgents; ++_i3) {
        var _ag4 = this.m_agents[_i3];

        if (!_ag4.active) {
          continue;
        }

        if (_ag4.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_NONE || _ag4.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_VELOCITY) {
          continue;
        }

        if (_ag4.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_WAITING_FOR_PATH) {
          // Poll path queue.
          var status = this.m_pathq.getRequestStatus(_ag4.targetPathqRef);

          if (status != null && status == _Status["default"].FAILURE) {
            // Path find failed, retry if the target location is still
            // valid.
            _ag4.targetPathqRef = _PathQueue["default"].DT_PATHQ_INVALID;

            if (_ag4.targetRef != 0) {
              _ag4.targetState = _CrowdAgent["default"].DT_CROWDAGENT_TARGET_REQUESTING;
            } else {
              _ag4.targetState = _CrowdAgent["default"].DT_CROWDAGENT_TARGET_FAILED;
            }

            _ag4.targetReplanTime = 0.0;
          } else if (status != null && (status == _Status["default"].SUCCSESS || status == _Status["default"].PARTIAL_RESULT)) {
            var _path = _ag4.corridor.getPath();

            if (_path.length == 0) {
              throw new IllegalArgumentException("Empty path");
            } // Apply results.


            var targetPos = _ag4.targetPos;
            var valid = true;

            var _pathFound = this.m_pathq.getPathResult(_ag4.targetPathqRef);

            var res = _pathFound.getRefs();

            status = _pathFound.getStatus();

            if (status == _Status["default"].FAILURE || res.length == 0) {
              valid = false;
            }

            if (status != null & status == _Status["default"].PARTIAL_RESULT) {
              _ag4.partial = true;
            } else {
              _ag4.partial = false;
            } // Merge result and existing path.
            // The agent might have moved whilst the request is
            // being processed, so the path may have changed.
            // We assume that the end of the path is at the same
            // location
            // where the request was issued.
            // The last ref in the old path should be the same as
            // the location where the request was issued..


            if (valid && _path[_path.length - 1] != res[0]) {
              // if (valid && path[path.length - 1].longValue() != res[0].longValue()) {
              valid = false;
            }

            if (valid) {
              // Put the old path infront of the old path.
              if (_path.length > 1) {
                // path.remove(path.length - 1);
                _path.splice(_path.length - 1, 1);

                _path.push.apply(_path, _toConsumableArray(res));

                res = _path; // Remove trackbacks

                for (var j = 1; j < res.length - 1; ++j) {
                  if (j - 1 >= 0 && j + 1 < res.length) {
                    if (res[j - 1] == res[j + 1]) {
                      res.splice(j + 1, 1);
                      res.splice(j, 1);
                      j -= 2;
                    }
                  }
                }
              } // Check for partial path.


              if (res[res.length - 1] != _ag4.targetRef) {
                // Partial path, constrain target position inside
                // the last polygon.
                var _cr = this.m_navquery.closestPointOnPoly(res[res.length - 1], targetPos);

                targetPos = _cr.getClosest();
              }
            }

            if (valid) {
              // Set current corridor.
              _ag4.corridor.setCorridor(targetPos, res); // Force to update boundary.


              _ag4.boundary.reset();

              _ag4.targetState = _CrowdAgent["default"].DT_CROWDAGENT_TARGET_VALID;
            } else {
              // Something went wrong.
              _ag4.targetState = _CrowdAgent["default"].DT_CROWDAGENT_TARGET_FAILED;
            }

            _ag4.targetReplanTime = 0.0;
          }
        }
      }
    }
  }, {
    key: "updateTopologyOptimization",
    // seconds
    value: function updateTopologyOptimization(agents, dt) {
      if (!agents.length == 0) {
        return;
      }

      var queue = new _PriorityQueue["default"](function (a1, a2) {
        return Float.compare(a2.topologyOptTime, a1.topologyOptTime);
      });

      for (var i = 0; i < agents.length; ++i) {
        ag = agents[i];

        if (ag.state != _CrowdAgent["default"].DT_CROWDAGENT_STATE_WALKING) {
          continue;
        }

        if (ag.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_NONE || ag.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_VELOCITY) {
          continue;
        }

        if ((ag.params.updateFlags & _CrowdAgentParams["default"].DT_CROWD_OPTIMIZE_TOPO) == 0) {
          continue;
        }

        ag.topologyOptTime += dt;

        if (ag.topologyOptTime >= OPT_TIME_THR) {
          addToOptQueue(ag, queue);
        }
      }

      while (!queue.length == 0) {
        ag = queue.poll();
        ag.corridor.optimizePathTopology(this.m_navquery, this.m_filters[ag.params.queryFilterType]);
        ag.topologyOptTime = 0;
      }
    }
  }, {
    key: "checkPathValidity",
    // seconds
    value: function checkPathValidity(agents, dt) {
      for (var i = 0; i < agents.length; ++i) {
        var _ag5 = agents[i];

        if (_ag5.state != _CrowdAgent["default"].DT_CROWDAGENT_STATE_WALKING) {
          continue;
        }

        _ag5.targetReplanTime += dt;
        var replan = false; // First check that the current location is valid.

        var agentPos = new Array(3);

        var agentRef = _ag5.corridor.getFirstPoly();

        _DetourCommon["default"].vCopy(agentPos, _ag5.npos);

        if (!this.m_navquery.isValidPolyRef(agentRef, this.m_filters[_ag5.params.queryFilterType])) {
          // Current location is not valid, try to reposition.
          // TODO: this can snap agents, how to handle that?
          var fnp = this.m_navquery.findNearestPoly(_ag5.npos, this.m_ext, this.m_filters[_ag5.params.queryFilterType]);
          agentRef = fnp.getNearestRef();

          if (fnp.getNearestPos() != null) {
            _DetourCommon["default"].vCopy(agentPos, fnp.getNearestPos());
          }

          if (agentRef == 0) {
            // Could not find location in navmesh, set state to invalid.
            _ag5.corridor.reset(0, agentPos);

            _ag5.partial = false;

            _ag5.boundary.reset();

            _ag5.state = _CrowdAgent["default"].DT_CROWDAGENT_STATE_INVALID;
            continue;
          } // Make sure the first polygon is valid, but leave other valid
          // polygons in the path so that replanner can adjust the path
          // better.


          _ag5.corridor.fixPathStart(agentRef, agentPos); // ag.corridor.trimInvalidPath(agentRef, agentPos, m_navquery,
          // &m_filter);


          _ag5.boundary.reset();

          _DetourCommon["default"].vCopy(_ag5.npos, agentPos);

          replan = true;
        } // If the agent does not have move target or is controlled by
        // velocity, no need to recover the target nor replan.


        if (_ag5.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_NONE || _ag5.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_VELOCITY) {
          continue;
        } // Try to recover move request position.


        if (_ag5.targetState != _CrowdAgent["default"].DT_CROWDAGENT_TARGET_NONE && _ag5.targetState != _CrowdAgent["default"].DT_CROWDAGENT_TARGET_FAILED) {
          if (!this.m_navquery.isValidPolyRef(_ag5.targetRef, this.m_filters[_ag5.params.queryFilterType])) {
            // Current target is not valid, try to reposition.
            var _fnp = this.m_navquery.findNearestPoly(_ag5.targetPos, this.m_ext, this.m_filters[_ag5.params.queryFilterType]);

            _ag5.targetRef = _fnp.getNearestRef();

            if (_fnp.getNearestPos() != null) {
              _DetourCommon["default"].vCopy(_ag5.targetPos, _fnp.getNearestPos());
            }

            replan = true;
          }

          if (_ag5.targetRef == 0) {
            // Failed to reposition target, fail moverequest.
            _ag5.corridor.reset(agentRef, agentPos);

            _ag5.partial = false;
            _ag5.targetState = _CrowdAgent["default"].DT_CROWDAGENT_TARGET_NONE;
          }
        } // If nearby corridor is not valid, replan.


        if (!_ag5.corridor.isValid(Crowd.CHECK_LOOKAHEAD, this.m_navquery, this.m_filters[_ag5.params.queryFilterType])) {
          // Fix current path.
          // ag.corridor.trimInvalidPath(agentRef, agentPos, m_navquery,
          // &m_filter);
          // ag.boundary.reset();
          replan = true;
        } // If the end of the path is near and it is not the requested
        // location, replan.


        if (_ag5.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_VALID) {
          if (_ag5.targetReplanTime > Crowd.TARGET_REPLAN_DELAY && _ag5.corridor.getPathCount() < Crowd.CHECK_LOOKAHEAD && _ag5.corridor.getLastPoly() != _ag5.targetRef) {
            replan = true;
          }
        } // Try to replan path to goal.


        if (replan) {
          if (_ag5.targetState != _CrowdAgent["default"].DT_CROWDAGENT_TARGET_NONE) {
            this.requestMoveTargetReplan(_ag5, _ag5.targetRef, _ag5.targetPos);
          }
        }
      }
    }
  }, {
    key: "update",
    value: function update(dt, debug, frame) {
      if (frame == 114) console.log("Bug");
      this.m_velocitySampleCount = 0;
      var debugIdx = debug != null ? debug.idx : -1;
      var agents = this.getActiveAgents(); // Check that all agents still have valid paths.

      this.checkPathValidity(agents, dt); // Update async move request and path finder.

      this.updateMoveRequest(); // Optimize path topology.

      this.updateTopologyOptimization(agents, dt); // Register agents to proximity grid.

      this.m_grid.clear();

      for (var i = 0; i < agents.length; ++i) {
        var _ag6 = agents[i];
        var p = _ag6.npos;
        var r = _ag6.params.radius;
        this.m_grid.addItem(i, p[0] - r, p[2] - r, p[0] + r, p[2] + r);
      } // Get nearby navmesh segments and agents to collide with.


      var _iterator2 = _createForOfIteratorHelper(agents),
          _step2;

      try {
        for (_iterator2.s(); !(_step2 = _iterator2.n()).done;) {
          var _ag13 = _step2.value;

          if (_ag13.state != _CrowdAgent["default"].DT_CROWDAGENT_STATE_WALKING) {
            continue;
          } // Update the collision boundary after certain distance has been passed or
          // if it has become invalid.


          var updateThr = _ag13.params.collisionQueryRange * 0.25;

          if (_DetourCommon["default"].vDist2DSqr(_ag13.npos, _ag13.boundary.getCenter()) > _DetourCommon["default"].sqr(updateThr) || !_ag13.boundary.isValid(this.m_navquery, this.m_filters[_ag13.params.queryFilterType])) {
            _ag13.boundary.update(_ag13.corridor.getFirstPoly(), _ag13.npos, _ag13.params.collisionQueryRange, this.m_navquery, this.m_filters[_ag13.params.queryFilterType]);
          } // Query neighbour agents


          _ag13.neis = this.getNeighbours(_ag13.npos, _ag13.params.height, _ag13.params.collisionQueryRange, _ag13, agents, this.m_grid);
        } // Find next corner to steer to.

      } catch (err) {
        _iterator2.e(err);
      } finally {
        _iterator2.f();
      }

      for (var _i4 = 0; _i4 < agents.length; ++_i4) {
        var _ag7 = agents[_i4];

        if (_ag7.state != _CrowdAgent["default"].DT_CROWDAGENT_STATE_WALKING) {
          continue;
        }

        if (_ag7.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_NONE || _ag7.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_VELOCITY) {
          continue;
        } // Find corners for steering


        _ag7.corners = _ag7.corridor.findCorners(Crowd.DT_CROWDAGENT_MAX_CORNERS, this.m_navquery, this.m_filters[_ag7.params.queryFilterType]); // Check to see if the corner after the next corner is directly visible,
        // and short cut to there.

        if ((_ag7.params.updateFlags & _CrowdAgentParams["default"].DT_CROWD_OPTIMIZE_VIS) != 0 && _ag7.corners.length > 0) {
          var target = _ag7.corners[Math.min(1, _ag7.corners.length - 1)].getPos();

          _ag7.corridor.optimizePathVisibility(target, _ag7.params.pathOptimizationRange, this.m_navquery, this.m_filters[_ag7.params.queryFilterType]); // Copy data for debug purposes.


          if (debugIdx == _i4) {
            _DetourCommon["default"].vCopy(debug.optStart, _ag7.corridor.getPos());

            _DetourCommon["default"].vCopy(debug.optEnd, target);
          }
        } else {
          // Copy data for debug purposes.
          if (debugIdx == _i4) {
            _DetourCommon["default"].vSet(debug.optStart, 0, 0, 0);

            _DetourCommon["default"].vSet(debug.optEnd, 0, 0, 0);
          }
        }
      } // Trigger off-mesh connections (depends on corners).


      var _iterator3 = _createForOfIteratorHelper(agents),
          _step3;

      try {
        for (_iterator3.s(); !(_step3 = _iterator3.n()).done;) {
          var _ag14 = _step3.value;

          if (_ag14.state != _CrowdAgent["default"].DT_CROWDAGENT_STATE_WALKING) {
            continue;
          }

          if (_ag14.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_NONE || _ag14.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_VELOCITY) {
            continue;
          } // Check


          var triggerRadius = _ag14.params.radius * 2.25;

          if (_ag14.overOffmeshConnection(triggerRadius)) {
            // Prepare to off-mesh connection.
            anim = _ag14.animation; // Adjust the path over the off-mesh connection.

            var refs = new long[2]();

            if (_ag14.corridor.moveOverOffmeshConnection(_ag14.corners[_ag14.corners.length - 1].getRef(), refs, anim.startPos, anim.endPos, this.m_navquery)) {
              _DetourCommon["default"].vCopy(anim.initPos, _ag14.npos);

              anim.polyRef = refs[1];
              anim.active = true;
              anim.t = 0.0;
              anim.tmax = _DetourCommon["default"].vDist2D(anim.startPos, anim.endPos) / _ag14.params.maxSpeed * 0.5;
              _ag14.state = _CrowdAgent["default"].DT_CROWDAGENT_STATE_OFFMESH;
              _ag14.corners = [];
              _ag14.neis = [];
              continue;
            } else {// Path validity check will ensure that bad/blocked connections will be replanned.
            }
          }
        } // Calculate steering.

      } catch (err) {
        _iterator3.e(err);
      } finally {
        _iterator3.f();
      }

      var _iterator4 = _createForOfIteratorHelper(agents),
          _step4;

      try {
        for (_iterator4.s(); !(_step4 = _iterator4.n()).done;) {
          var _ag15 = _step4.value;

          if (_ag15.state != _CrowdAgent["default"].DT_CROWDAGENT_STATE_WALKING) {
            continue;
          }

          if (_ag15.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_NONE) {
            continue;
          }

          var dvel = new Array(3);

          if (_ag15.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_VELOCITY) {
            _DetourCommon["default"].vCopy(dvel, _ag15.targetPos);

            _ag15.desiredSpeed = _DetourCommon["default"].vLen(_ag15.targetPos);
          } else {
            // Calculate steering direction.
            if ((_ag15.params.updateFlags & _CrowdAgentParams["default"].DT_CROWD_ANTICIPATE_TURNS) != 0) {
              dvel = _ag15.calcSmoothSteerDirection();
            } else {
              dvel = _ag15.calcStraightSteerDirection();
            } // Calculate speed scale, which tells the agent to slowdown at the end of the path.


            var slowDownRadius = _ag15.params.radius * 2; // TODO: make less hacky.

            var speedScale = _ag15.getDistanceToGoal(slowDownRadius) / slowDownRadius;
            _ag15.desiredSpeed = _ag15.params.maxSpeed;
            dvel = _DetourCommon["default"].vScale(dvel, _ag15.desiredSpeed * speedScale);
          } // Separation


          if ((_ag15.params.updateFlags & _CrowdAgentParams["default"].DT_CROWD_SEPARATION) != 0) {
            separationDist = _ag15.params.collisionQueryRange;
            invSeparationDist = 1.0 / separationDist;
            separationWeight = _ag15.params.separationWeight;
            w = 0;
            var disp = new Array(3);

            for (var _j3 = 0; _j3 < _ag15.neis.length; ++_j3) {
              nei = agents.get(_ag15.neis[_j3].idx);

              var _diff = _DetourCommon["default"].vSub(_ag15.npos, nei.npos);

              _diff[1] = 0;
              distSqr = _DetourCommon["default"].vLenSqr(_diff);

              if (distSqr < 0.00001) {
                continue;
              }

              if (distSqr > _DetourCommon["default"].sqr(separationDist)) {
                continue;
              }

              dist = Math.sqrt(distSqr);
              weight = separationWeight * (1.0 - _DetourCommon["default"].sqr(dist * invSeparationDist));
              disp = _DetourCommon["default"].vMad(disp, _diff, weight / dist);
              w += 1.0;
            }

            if (w > 0.0001) {
              // Adjust desired velocity.
              dvel = _DetourCommon["default"].vMad(dvel, disp, 1.0 / w); // Clamp desired velocity to desired speed.

              speedSqr = _DetourCommon["default"].vLenSqr(dvel);
              desiredSqr = _DetourCommon["default"].sqr(_ag15.desiredSpeed);

              if (speedSqr > desiredSqr) {
                dvel = _DetourCommon["default"].vScale(dvel, desiredSqr / speedSqr);
              }
            }
          } // Set the desired velocity.


          _DetourCommon["default"].vCopy(_ag15.dvel, dvel);
        } // Velocity planning.

      } catch (err) {
        _iterator4.e(err);
      } finally {
        _iterator4.f();
      }

      for (var _i5 = 0; _i5 < agents.length; ++_i5) {
        var _ag8 = agents[_i5];

        if (_ag8.state != _CrowdAgent["default"].DT_CROWDAGENT_STATE_WALKING) {
          continue;
        }

        if ((_ag8.params.updateFlags & _CrowdAgentParams["default"].DT_CROWD_OBSTACLE_AVOIDANCE) != 0) {
          this.m_obstacleQuery.reset(); // Add neighbours as obstacles.

          for (var j = 0; j < _ag8.neis.length; ++j) {
            var _nei = agents[_ag8.neis[j].idx];
            this.m_obstacleQuery.addCircle(_nei.npos, _nei.params.radius, _nei.vel, _nei.dvel);
          }

          if (_ag8.neis.length > 0 && _i5 == 0) {
            console.log("Frame " + frame);
          } // Append neighbour segments as obstacles.


          for (var _j = 0; _j < _ag8.boundary.getSegmentCount(); ++_j) {
            var s = _ag8.boundary.getSegment(_j); //let s3 = Arrays.copyOfRange(s, 3, 6);
            //let s3 = Arrays.copyOfRange(s, 3, 6);


            var s3 = s.slice(3, 6);

            if (_DetourCommon["default"].triArea2D3(_ag8.npos, s, s3) < 0.0) {
              continue;
            }

            this.m_obstacleQuery.addSegment(s, s3);
          }

          var vod = null;

          if (debugIdx == _i5) {
            vod = debug.vod;
          } // Sample new safe velocity.


          var adaptive = true;
          var ns = 0;
          var params = this.m_obstacleQueryParams[_ag8.params.obstacleAvoidanceType];

          if (adaptive) {
            var nsnvel = this.m_obstacleQuery.sampleVelocityAdaptive(_ag8.npos, _ag8.params.radius, _ag8.desiredSpeed, _ag8.vel, _ag8.dvel, params, vod);
            ns = nsnvel[0];
            _ag8.nvel = nsnvel[1];
          } else {
            var _nsnvel = this.m_obstacleQuery.sampleVelocityGrid(_ag8.npos, _ag8.params.radius, _ag8.desiredSpeed, _ag8.vel, _ag8.dvel, params, vod);

            ns = _nsnvel[0];
            _ag8.nvel = _nsnvel[1];
          }

          this.m_velocitySampleCount += ns;
        } else {
          // If not using velocity planning, new velocity is directly the desired velocity.
          _DetourCommon["default"].vCopy(_ag8.nvel, _ag8.dvel);
        }
      } // Integrate.


      for (var _i6 = 0; _i6 < agents.length; ++_i6) {
        var _ag9 = agents[_i6];

        if (_ag9.state != _CrowdAgent["default"].DT_CROWDAGENT_STATE_WALKING) {
          continue;
        }

        _ag9.integrate(dt);
      } // Handle collisions.


      for (var iter = 0; iter < 4; ++iter) {
        for (var _i7 = 0; _i7 < agents.length; ++_i7) {
          var _ag10 = agents[_i7];

          var idx0 = _ag10.getAgentIndex();

          if (_ag10.state != _CrowdAgent["default"].DT_CROWDAGENT_STATE_WALKING) {
            continue;
          }

          _DetourCommon["default"].vSet(_ag10.disp, 0, 0, 0);

          var _w = 0;

          for (var _j2 = 0; _j2 < _ag10.neis.length; ++_j2) {
            var _nei2 = agents[_ag10.neis[_j2].idx];

            var idx1 = _nei2.getAgentIndex();

            var diff = _DetourCommon["default"].vSub(_ag10.npos, _nei2.npos);

            diff[1] = 0;

            var _dist = _DetourCommon["default"].vLenSqr(diff);

            if (_dist > _DetourCommon["default"].sqr(_ag10.params.radius + _nei2.params.radius)) {
              continue;
            }

            _dist = Math.sqrt(_dist);
            var pen = _ag10.params.radius + _nei2.params.radius - _dist;

            if (_dist < 0.0001) {
              // Agents on top of each other, try to choose diverging separation directions.
              if (idx0 > idx1) {
                _DetourCommon["default"].vSet(diff, -_ag10.dvel[2], 0, _ag10.dvel[0]);
              } else {
                _DetourCommon["default"].vSet(diff, _ag10.dvel[2], 0, -_ag10.dvel[0]);
              }

              pen = 0.01;
            } else {
              pen = 1.0 / _dist * (pen * 0.5) * Crowd.COLLISION_RESOLVE_FACTOR;
            }

            _ag10.disp = _DetourCommon["default"].vMad(_ag10.disp, diff, pen);
            _w += 1.0;
          }

          if (_w > 0.0001) {
            var iw = 1.0 / _w;
            _ag10.disp = _DetourCommon["default"].vScale(_ag10.disp, iw);
          }
        }

        for (var _i8 = 0; _i8 < agents.length; ++_i8) {
          var _ag11 = agents[_i8];

          if (_ag11.state != _CrowdAgent["default"].DT_CROWDAGENT_STATE_WALKING) {
            continue;
          }

          _ag11.npos = _DetourCommon["default"].vAdd(_ag11.npos, _ag11.disp);
        }
      }

      for (var _i9 = 0; _i9 < agents.length; ++_i9) {
        if (frame == 492) console.log("Bad agent");
        var _ag12 = agents[_i9];

        if (_ag12.state != _CrowdAgent["default"].DT_CROWDAGENT_STATE_WALKING) {
          continue;
        } // Move along navmesh.


        _ag12.corridor.movePosition(_ag12.npos, this.m_navquery, this.m_filters[_ag12.params.queryFilterType]); // Get valid constrained position back.


        _DetourCommon["default"].vCopy(_ag12.npos, _ag12.corridor.getPos()); // If not using path, truncate the corridor to just one poly.


        if (_ag12.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_NONE || _ag12.targetState == _CrowdAgent["default"].DT_CROWDAGENT_TARGET_VELOCITY) {
          _ag12.corridor.reset(_ag12.corridor.getFirstPoly(), _ag12.npos);

          _ag12.partial = false;
        }
      } // Update agents using off-mesh connection.


      for (var _i10 = 0; _i10 < this.m_maxAgents; ++_i10) {
        var _anim = this.m_agents[_i10].animation;

        if (!_anim.active) {
          continue;
        }

        ag = this.m_agents[_i10];
        _anim.t += dt;

        if (_anim.t > _anim.tmax) {
          // Reset animation
          _anim.active = false; // Prepare agent for walking.

          ag.state = _CrowdAgent["default"].DT_CROWDAGENT_STATE_WALKING;
          continue;
        } // Update position


        ta = _anim.tmax * 0.15;
        tb = _anim.tmax;

        if (_anim.t < ta) {
          u = tween(_anim.t, 0.0, ta);
          ag.npos = _DetourCommon["default"].vLerp3(_anim.initPos, _anim.startPos, u);
        } else {
          u = tween(_anim.t, ta, tb);
          ag.npos = _DetourCommon["default"].vLerp3(_anim.startPos, _anim.endPos, u);
        } // Update velocity.


        _DetourCommon["default"].vSet(ag.vel, 0, 0, 0);

        _DetourCommon["default"].vSet(ag.dvel, 0, 0, 0);
      }
    }
  }, {
    key: "getQueryExtents",
    value: function getQueryExtents() {
      return this.m_ext;
    }
  }, {
    key: "getFilter",
    value: function getFilter(i) {
      return i >= 0 && i < Crowd.DT_CROWD_MAX_QUERY_FILTER_TYPE ? this.m_filters[i] : null;
    }
  }]);

  return Crowd;
}();

_defineProperty(Crowd, "MAX_ITERS_PER_UPDATE", 100);

_defineProperty(Crowd, "MAX_PATHQUEUE_NODES", 4096);

_defineProperty(Crowd, "MAX_COMMON_NODES", 512);

_defineProperty(Crowd, "DT_CROWDAGENT_MAX_NEIGHBOURS", 6);

_defineProperty(Crowd, "DT_CROWDAGENT_MAX_CORNERS", 4);

_defineProperty(Crowd, "DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS", 8);

_defineProperty(Crowd, "DT_CROWD_MAX_QUERY_FILTER_TYPE", 16);

_defineProperty(Crowd, "MAX_ITER", 20);

_defineProperty(Crowd, "OPT_TIME_THR", 0.5);

_defineProperty(Crowd, "CHECK_LOOKAHEAD", 10);

_defineProperty(Crowd, "TARGET_REPLAN_DELAY", 1.0);

_defineProperty(Crowd, "COLLISION_RESOLVE_FACTOR", 0.7);

var _default = Crowd;
exports["default"] = _default;