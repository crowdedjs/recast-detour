"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _PathQuery = _interopRequireDefault(require("./PathQuery.js"));

var _NavMeshQuery = _interopRequireDefault(require("../NavMeshQuery.js"));

var _DetourCommon = _interopRequireDefault(require("../DetourCommon.js"));

var _Status = _interopRequireDefault(require("../Status.js"));

var _FindPathResult = _interopRequireDefault(require("../FindPathResult.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var PathQueue = /*#__PURE__*/function () {
  // in update ticks.
  function PathQueue(maxSearchNodeCount, nav) {
    _classCallCheck(this, PathQueue);

    _defineProperty(this, "m_queue", new Array(PathQueue.MAX_QUEUE));

    _defineProperty(this, "m_nextHandle", 1);

    _defineProperty(this, "m_queueHead", void 0);

    _defineProperty(this, "m_navquery", void 0);

    this.m_navquery = new _NavMeshQuery["default"](nav);

    for (var i = 0; i < PathQueue.MAX_QUEUE; ++i) {
      this.m_queue[i] = new _PathQuery["default"]();
      this.m_queue[i].ref = PathQueue.DT_PATHQ_INVALID;
      this.m_queue[i].path = new Array(256);
    }

    this.m_queueHead = 0;
  }

  _createClass(PathQueue, [{
    key: "update",
    value: function update(maxIters) {
      // Update path request until there is nothing to update
      // or upto maxIters pathfinder iterations has been consumed.
      var iterCount = maxIters;

      for (var i = 0; i < PathQueue.MAX_QUEUE; ++i) {
        var q = this.m_queue[this.m_queueHead % PathQueue.MAX_QUEUE]; // Skip inactive requests.

        if (q.ref == PathQueue.DT_PATHQ_INVALID) {
          this.m_queueHead++;
          continue;
        } // Handle compPolyed request.


        if (q.status != null && (q.status == _Status["default"].SUCCSESS || q.status == _Status["default"].PARTIAL_RESULT || q.status == _Status["default"].FAILURE)) {
          // If the path result has not been read in few frames, free the slot.
          q.keepAlive++;

          if (q.keepAlive > PathQueue.MAX_KEEP_ALIVE) {
            q.ref = PathQueue.DT_PATHQ_INVALID;
            q.status = null;
          }

          this.m_queueHead++;
          continue;
        } // Handle query start.


        if (q.status == null) {
          q.status = this.m_navquery.initSlicedFindPath(q.startRef, q.endRef, q.startPos, q.endPos, q.filter, 0);
        } // Handle query in progress.


        if (q.status == _Status["default"].IN_PROGRESS) {
          var iters = 0;
          var res = this.m_navquery.updateSlicedFindPath(iterCount);
          iters = res.getIterations();
          q.status = res.getStatus();
          iterCount -= iters;
        }

        if (q.status == _Status["default"].SUCCSESS || q.status == _Status["default"].PARTIAL_RESULT) {
          var path = this.m_navquery.finalizeSlicedFindPath();
          q.status = path.getStatus();
          q.path = path.getRefs();
        }

        if (iterCount <= 0) break;
        this.m_queueHead++;
      }
    }
  }, {
    key: "request",
    value: function request(startRef, endRef, startPos, endPos, filter) {
      // Find empty slot
      var slot = -1;

      for (var i = 0; i < PathQueue.MAX_QUEUE; ++i) {
        if (this.m_queue[i].ref == PathQueue.DT_PATHQ_INVALID) {
          slot = i;
          break;
        }
      } // Could not find slot.


      if (slot == -1) return PathQueue.DT_PATHQ_INVALID;
      var ref = this.m_nextHandle++;
      if (this.m_nextHandle == PathQueue.DT_PATHQ_INVALID) this.m_nextHandle++;
      var q = this.m_queue[slot];
      q.ref = ref;

      _DetourCommon["default"].vCopy(q.startPos, startPos);

      q.startRef = startRef;

      _DetourCommon["default"].vCopy(q.endPos, endPos);

      q.endRef = endRef;
      q.status = null;
      q.filter = filter;
      q.keepAlive = 0;
      return ref;
    }
  }, {
    key: "getRequestStatus",
    value: function getRequestStatus(ref) {
      for (var i = 0; i < PathQueue.MAX_QUEUE; ++i) {
        if (this.m_queue[i].ref == ref) return this.m_queue[i].status;
      }

      return _Status["default"].FAILURE;
    }
  }, {
    key: "getPathResult",
    value: function getPathResult(ref) {
      for (var i = 0; i < PathQueue.MAX_QUEUE; ++i) {
        if (this.m_queue[i].ref == ref) {
          var q = this.m_queue[i]; // Free request for reuse.

          q.ref = PathQueue.DT_PATHQ_INVALID;
          q.status = null;
          return new _FindPathResult["default"](_Status["default"].SUCCSESS, q.path);
        }
      }

      return new _FindPathResult["default"](_Status["default"].FAILURE, null);
    }
  }]);

  return PathQueue;
}();

_defineProperty(PathQueue, "MAX_QUEUE", 8);

_defineProperty(PathQueue, "DT_PATHQ_INVALID", 0);

_defineProperty(PathQueue, "MAX_KEEP_ALIVE", 2);

var _default = PathQueue;
exports["default"] = _default;