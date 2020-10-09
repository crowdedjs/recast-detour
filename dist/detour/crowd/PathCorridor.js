"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _DetourCommon = _interopRequireDefault(require("../DetourCommon.js"));

var _NavMeshQuery = _interopRequireDefault(require("../NavMeshQuery.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _toConsumableArray(arr) { return _arrayWithoutHoles(arr) || _iterableToArray(arr) || _unsupportedIterableToArray(arr) || _nonIterableSpread(); }

function _nonIterableSpread() { throw new TypeError("Invalid attempt to spread non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _iterableToArray(iter) { if (typeof Symbol !== "undefined" && Symbol.iterator in Object(iter)) return Array.from(iter); }

function _arrayWithoutHoles(arr) { if (Array.isArray(arr)) return _arrayLikeToArray(arr); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

/**
 * Represents a dynamic polygon corridor used to plan agent movement.
 * 
 * The corridor is loaded with a path, usually obtained from a
 * #NavMeshQuery::findPath() query. The corridor is then used to plan local
 * movement, with the corridor automatically updating as needed to deal with
 * inaccurate agent locomotion.
 * 
 * Example of a common use case:
 * 
 * -# Construct the corridor object and call 
 * -# Obtain a path from a #dtNavMeshQuery object. 
 * -# Use #reset() to set the agent's current position. (At the beginning of the path.) -# Use
 * #setCorridor() to load the path and target. -# Use #findCorners() to plan
 * movement. (This handles dynamic path straightening.) -# Use #movePosition()
 * to feed agent movement back into the corridor. (The corridor will
 * automatically adjust as needed.) -# If the target is moving, use
 * #moveTargetPosition() to update the end of the corridor. (The corridor will
 * automatically adjust as needed.) -# Repeat the previous 3 steps to continue
 * to move the agent.
 * 
 * The corridor position and target are always constrained to the navigation
 * mesh.
 * 
 * One of the difficulties in maintaining a path is that ing poPoly errors,
 * locomotion inaccuracies, and/or local steering can result in the agent
 * crossing the boundary of the path corridor, temporarily invalidating the
 * path. This class uses local mesh queries to detect and update the corridor as
 * needed to handle these types of issues.
 * 
 * The fact that local mesh queries are used to move the position and target
 * locations results in two beahviors that need to be considered:
 * 
 * Every time a move function is used there is a chance that the path will
 * become non-optimial. Basically, the further the target is moved from its
 * original location, and the further the position is moved outside the original
 * corridor, the more likely the path will become non-optimal. This issue can be
 * addressed by periodically running the #optimizePathTopology() and
 * #optimizePathVisibility() methods.
 * 
 * All local mesh queries have distance limitations. (Review the #dtNavMeshQuery
 * methods for details.) So the most accurate use case is to move the position
 * and target in small increments. If a large increment is used, then the
 * corridor may not be able to accurately find the new location. Because of this
 * limiation, if a position is moved in a large increment, then compare the
 * desired and resulting polygon references. If the two do not match, then path
 * replanning may be needed. E.g. If you move the target, check #getLastPoly()
 * to see if it is the expected polygon.
 *
 */
var PathCorridor = /*#__PURE__*/function () {
  _createClass(PathCorridor, [{
    key: "mergeCorridorStartMoved",
    value: function mergeCorridorStartMoved(path, visited) {
      var furthestPath = -1;
      var furthestVisited = -1; // Find furthest common polygon.

      for (var i = path.length - 1; i >= 0; --i) {
        var found = false;

        for (var j = visited.length - 1; j >= 0; --j) {
          if (path[i] == visited[j]) {
            furthestPath = i;
            furthestVisited = j;
            found = true;
          }
        }

        if (found) break;
      } // If no intersection found just return current path.


      if (furthestPath == -1 || furthestVisited == -1) return path; // Concatenate paths.
      // Adjust beginning of the buffer to include the visited.

      var result = []; // Store visited

      for (var _i = visited.length - 1; _i > furthestVisited; --_i) {
        result.push(visited[_i]);
      }

      result.push.apply(result, _toConsumableArray(path.slice(furthestPath, path.length)));
      return result;
    }
  }, {
    key: "mergeCorridorEndMoved",
    value: function mergeCorridorEndMoved(path, visited) {
      var furthestPath = -1;
      var furthestVisited = -1; // Find furthest common polygon.

      for (var i = 0; i < path.length; ++i) {
        var found = false;

        for (var j = visited.length - 1; j >= 0; --j) {
          if (path[i] == visited[j]) {
            furthestPath = i;
            furthestVisited = j;
            found = true;
          }
        }

        if (found) break;
      } // If no intersection found just return current path.


      if (furthestPath == -1 || furthestVisited == -1) return path; // Concatenate paths.

      var result = path.subList(0, furthestPath);
      result.addAll(visited.subList(furthestVisited, visited.length));
      return result;
    }
  }, {
    key: "mergeCorridorStartShortcut",
    value: function mergeCorridorStartShortcut(path, visited) {
      var furthestPath = -1;
      var furthestVisited = -1; // Find furthest common polygon.

      for (var i = path.length - 1; i >= 0; --i) {
        var found = false;

        for (var j = visited.length - 1; j >= 0; --j) {
          if (path[i] == visited[j]) {
            furthestPath = i;
            furthestVisited = j;
            found = true;
          }
        }

        if (found) break;
      } // If no intersection found just return current path.


      if (furthestPath == -1 || furthestVisited <= 0) return path; // Concatenate paths.
      // Adjust beginning of the buffer to include the visited.

      var result = visited.subList(0, furthestVisited);
      result.addAll(path.subList(furthestPath, path.length));
      return result;
    }
    /**
     * Allocates the corridor's path buffer.
     */

  }]);

  function PathCorridor() {
    _classCallCheck(this, PathCorridor);

    _defineProperty(this, "m_pos", new Array(3));

    _defineProperty(this, "m_target", new Array(3));

    _defineProperty(this, "m_path", []);

    this.m_path = [];
  }
  /**
   * Resets the path corridor to the specified position.
   * @param ref The polygon reference containing the position.
   * @param pos The new position in the corridor. [(x, y, z)]
   */


  _createClass(PathCorridor, [{
    key: "reset",
    value: function reset(ref, pos) {
      this.m_path = [];
      this.m_path.push(ref);

      _DetourCommon["default"].vCopy(this.m_pos, pos);

      _DetourCommon["default"].vCopy(this.m_target, pos);
    }
    /**
     * Finds the corners in the corridor from the position toward the target.
     * (The straightened path.)
     * 
     * This is the function used to plan local movement within the corridor. One
     * or more corners can be detected in order to plan movement. It performs
     * essentially the same function as #dtNavMeshQuery::findStraightPath.
     * 
     * Due to internal optimizations, the maximum number of corners returned
     * will be (@p maxCorners - 1) For example: If the buffers are sized to hold
     * 10 corners, the function will never return more than 9 corners. So if 10
     * corners are needed, the buffers should be sized for 11 corners.
     * 
     * If the target is within range, it will be the last corner and have a
     * polygon reference id of zero.
     * @param filter 
     * 
     * @param[in] navquery The query object used to build the corridor.
     * @return Corners
     */

  }, {
    key: "findCorners",
    value: function findCorners(maxCorners, navquery, filter) {
      var MIN_TARGET_DIST = _DetourCommon["default"].sqr(0.01);

      var path = navquery.findStraightPath(this.m_pos, this.m_target, this.m_path, maxCorners, 0); // Prune points in the beginning of the path which are too close.
      // for (let iter = path.iterator(); iter.hasNext();) {
      // 	let spi = iter.next();
      // 	if ((spi.getFlags() & NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION) != 0
      // 		|| DetourCommon.vDist2D(spi.getPos(), this.m_pos) > MIN_TARGET_DIST) {
      // 		break;
      // 	}
      // 	iter.remove();
      // }
      // TODO

      var newPath = [];
      var i = 0;
      var l = path.length;

      for (; i < path.length; i++) {
        var spi = path[i];

        if ((spi.getFlags() & _NavMeshQuery["default"].DT_STRAIGHTPATH_OFFMESH_CONNECTION) != 0 || _DetourCommon["default"].vDist2D(spi.getPos(), this.m_pos) > MIN_TARGET_DIST) {
          break;
        }
      }

      var j = i; // for(; j < l; j++);
      // {
      // 	newPath.push(path[j])
      // }

      while (j < l) {
        newPath.push(path[j]);
        j++;
      }

      path = newPath;
      return path;
    }
    /**
     * Attempts to optimize the path if the specified poPoly is visible from the
     * current position.
     * 
     * Inaccurate locomotion or dynamic obstacle avoidance can force the agent
     * position significantly outside the original corridor. Over time this can
     * result in the formation of a non-optimal corridor. Non-optimal paths can
     * also form near the corners of tiles.
     * 
     * This function uses an efficient local visibility search to try to
     * optimize the corridor between the current position and @p next.
     * 
     * The corridor will change only if @p next is visible from the current
     * position and moving directly toward the poPoly is better than following
     * the existing path.
     * 
     * The more inaccurate the agent movement, the more beneficial this function
     * becomes. Simply adjust the frequency of the call to match the needs to
     * the agent.
     * 
     * This function is not suitable for let distance searches.
     * 
     * @param next
     *            The poPoly to search toward. [(x, y, z])
     * @param pathOptimizationRange
     *            The maximum range to search. [Limit: > 0]
     * @param navquery
     *            The query object used to build the corridor.
     * @param filter
     *            The filter to apply to the operation.
     */

  }, {
    key: "optimizePathVisibility",
    value: function optimizePathVisibility(next, pathOptimizationRange, navquery, filter) {
      // Clamp the ray to max distance.
      var dist = _DetourCommon["default"].vDist2D(this.m_pos, next); // If too close to the goal, do not try to optimize.


      if (dist < 0.01) return; // Overshoot a little. This helps to optimize open fields in tiled
      // meshes.

      dist = Math.min(dist + 0.01, pathOptimizationRange); // Adjust ray length.

      var delta = _DetourCommon["default"].vSub(next, this.m_pos);

      var goal = _DetourCommon["default"].vMad(this.m_pos, delta, pathOptimizationRange / dist);

      var rc = navquery.raycast(this.m_path[0], this.m_pos, goal, filter, 0, 0);

      if (rc.path.length > 1 && rc.t > 0.99) {
        this.m_path = this.mergeCorridorStartShortcut(this.m_path, rc.path);
      }
    }
    /**
     * Attempts to optimize the path using a local area search. (Partial
     * replanning.)
     * 
     * Inaccurate locomotion or dynamic obstacle avoidance can force the agent
     * position significantly outside the original corridor. Over time this can
     * result in the formation of a non-optimal corridor. This function will use
     * a local area path search to try to re-optimize the corridor.
     * 
     * The more inaccurate the agent movement, the more beneficial this function
     * becomes. Simply adjust the frequency of the call to match the needs to
     * the agent.
     * 
     * @param navquery The query object used to build the corridor.
     * @param filter The filter to apply to the operation.
     * 
     */

  }, {
    key: "optimizePathTopology",
    value: function optimizePathTopology(navquery, filter) {
      if (this.m_path.length < 3) return false;
      var MAX_ITER = 32;
      navquery.initSlicedFindPath(this.m_path[0], this.m_path[this.m_path.length - 1], this.m_pos, this.m_target, filter, 0);
      navquery.updateSlicedFindPath(MAX_ITER);
      var fpr = navquery.finalizeSlicedFindPathPartial(this.m_path);

      if (fpr.getStatus().isSuccess() && fpr.getRefs().length > 0) {
        this.m_path = mergeCorridorStartShortcut(this.m_path, fpr.getRefs());
        return true;
      }

      return false;
    }
  }, {
    key: "moveOverOffmeshConnection",
    value: function moveOverOffmeshConnection(offMeshConRef, refs, start, end, navquery) {
      // Advance the path up to and over the off-mesh connection.
      var prevRef = 0;
      var polyRef = this.m_path[0];
      var npos = 0;

      while (npos < this.m_path.length && polyRef != offMeshConRef) {
        prevRef = polyRef;
        polyRef = this.m_path[npos];
        npos++;
      }

      if (npos == this.m_path.length) {
        // Could not find offMeshConRef
        return false;
      } // Prune path


      this.m_path = this.m_path.subList(npos, this.m_path.length);
      refs[0] = prevRef;
      refs[1] = polyRef;
      var nav = navquery.getAttachedNavMesh();
      var startEnd = nav.getOffMeshConnectionPolyEndPoints(refs[0], refs[1]);

      _DetourCommon["default"].vCopy(this.m_pos, startEnd.second);

      _DetourCommon["default"].vCopy(start, startEnd.first);

      _DetourCommon["default"].vCopy(end, startEnd.second);

      return true;
    }
    /**
     * Moves the position from the current location to the desired location,
     * adjusting the corridor as needed to reflect the change.
     * 
     * Behavior:
     * 
     * - The movement is constrained to the surface of the navigation mesh. -
     * The corridor is automatically adjusted (shorted or lengthened) in order
     * to remain valid. - The new position will be located in the adjusted
     * corridor's first polygon.
     * 
     * The expected use case is that the desired position will be 'near' the
     * current corridor. What is considered 'near' depends on local polygon
     * density, query search extents, etc.
     * 
     * The resulting position will differ from the desired position if the
     * desired position is not on the navigation mesh, or it can't be reached
     * using a local search.
     * 
     * @param npos
     *            The desired new position. [(x, y, z)]
     * @param navquery
     *            The query object used to build the corridor.
     * @param filter
     *            The filter to apply to the operation.
     */

  }, {
    key: "movePosition",
    value: function movePosition(npos, navquery, filter) {
      // Move aPoly navmesh and update new position.
      var masResult = navquery.moveAlongSurface(this.m_path[0], this.m_pos, npos, filter);
      this.m_path = this.mergeCorridorStartMoved(this.m_path, masResult.getVisited()); // Adjust the position to stay on top of the navmesh.

      _DetourCommon["default"].vCopy(this.m_pos, masResult.getResultPos());

      try {
        this.m_pos[1] = navquery.getPolyHeight(this.m_path[0], masResult.getResultPos());
      } catch (e) {// Silently disregard the returned status of DT_FAILURE | DT_INVALID_PARAM to stay
        // consistent with what the library in C++ does, see DetourPathCorridor.cpp.
      }
    }
    /**
     * Moves the target from the curent location to the desired location,
     * adjusting the corridor as needed to reflect the change. Behavior: - The
     * movement is constrained to the surface of the navigation mesh. - The
     * corridor is automatically adjusted (shorted or lengthened) in order to
     * remain valid. - The new target will be located in the adjusted corridor's
     * last polygon.
     * 
     * The expected use case is that the desired target will be 'near' the
     * current corridor. What is considered 'near' depends on local polygon
     * density, query search extents, etc. The resulting target will differ from
     * the desired target if the desired target is not on the navigation mesh,
     * or it can't be reached using a local search.
     *
     * @param npos
     *            The desired new target position. [(x, y, z)]
     * @param navquery
     *            The query object used to build the corridor.
     * @param filter
     *            The filter to apply to the operation.
     */

  }, {
    key: "moveTargetPosition",
    value: function moveTargetPosition(npos, navquery, filter) {
      // Move aPoly navmesh and update new position.
      var masResult = navquery.moveAlongSurface(this.m_path[this.m_path.length - 1], this.m_target, npos, filter);
      this.m_path = mergeCorridorEndMoved(this.m_path, masResult.getVisited()); // TODO: should we do that?
      // Adjust the position to stay on top of the navmesh.

      /*
       *  h = m_target[1]; navquery=>getPolyHeight(this.m_path[m_npath-1],
       * result, &h); result[1] = h;
       */

      _DetourCommon["default"].vCopy(this.m_target, masResult.getResultPos());
    }
    /**
     * Loads a new path and target into the corridor. The current corridor
     * position is expected to be within the first polygon in the path. The
     * target is expected to be in the last polygon.
     * 
     * @warning The size of the path must not exceed the size of corridor's path
     *          buffer set during #init().
     * @param target
     *            The target location within the last polygon of the path. [(x,
     *            y, z)]
     * @param path
     *            The path corridor.
     */

  }, {
    key: "setCorridor",
    value: function setCorridor(target, path) {
      _DetourCommon["default"].vCopy(this.m_target, target);

      this.m_path = _toConsumableArray(path);
    }
  }, {
    key: "fixPathStart",
    value: function fixPathStart(safeRef, safePos) {
      _DetourCommon["default"].vCopy(this.m_pos, safePos);

      if (this.m_path.length < 3 && this.m_path.length > 0) {
        var p = this.m_path[this.m_path.length - 1];
        this.m_path = [];
        this.m_path.push(safeRef);
        this.m_path.push(0);
        this.m_path.push(p);
      } else {
        this.m_path = [];
        this.m_path.push(safeRef);
        this.m_path.push(0);
      }
    }
  }, {
    key: "trimInvalidPath",
    value: function trimInvalidPath(safeRef, safePos, navquery, filter) {
      // Keep valid path as far as possible.
      var n = 0;

      while (n < this.m_path.length && navquery.isValidPolyRef(this.m_path[n], filter)) {
        n++;
      }

      if (n == 0) {
        // The first polyref is bad, use current safe values.
        _DetourCommon["default"].vCopy(this.m_pos, safePos);

        this.m_path = [];
        this.m_path.push(safeRef);
      } else if (n < this.m_path.length) {
        this.m_path = this.m_path.subList(0, n); // The path is partially usable.
      } // Clamp target pos to last poly


      _DetourCommon["default"].vCopy(this.m_target, navquery.closestPointOnPolyBoundary(this.m_path[this.m_path.length - 1], this.m_target));
    }
    /**
     * Checks the current corridor path to see if its polygon references remain
     * valid. The path can be invalidated if there are structural changes to the
     * underlying navigation mesh, or the state of a polygon within the path
     * changes resulting in it being filtered out. (E.g. An exclusion or
     * inclusion flag changes.)
     * 
     * @param maxLookAhead
     *            The number of polygons from the beginning of the corridor to
     *            search.
     * @param navquery
     *            The query object used to build the corridor.
     * @param filter
     *            The filter to apply to the operation.
     * @return
     */

  }, {
    key: "isValid",
    value: function isValid(maxLookAhead, navquery, filter) {
      // Check that all polygons still pass query filter.
      var n = Math.min(this.m_path.length, maxLookAhead);

      for (var i = 0; i < n; ++i) {
        if (!navquery.isValidPolyRef(this.m_path[i], filter)) return false;
      }

      return true;
    }
    /**
     * Gets the current position within the corridor. (In the first polygon.)
     * @return The current position within the corridor.
     */

  }, {
    key: "getPos",
    value: function getPos() {
      return this.m_pos;
    }
    /**
     * Gets the current target within the corridor. (In the last polygon.)
     * @return The current target within the corridor.
     */

  }, {
    key: "getTarget",
    value: function getTarget() {
      return this.m_target;
    }
    /**
     * The polygon reference id of the first polygon in the corridor, the polygon containing the position.
     * @return The polygon reference id of the first polygon in the corridor. (Or zero if there is no path.)
     */

  }, {
    key: "getFirstPoly",
    value: function getFirstPoly() {
      return this.m_path.length == 0 ? 0 : this.m_path[0];
    }
    /**
     * The polygon reference id of the last polygon in the corridor, the polygon containing the target.
     * @return The polygon reference id of the last polygon in the corridor. (Or zero if there is no path.)
     */

  }, {
    key: "getLastPoly",
    value: function getLastPoly() {
      return this.m_path.length == 0 ? 0 : this.m_path[this.m_path.length - 1];
    }
    /**
     * The corridor's path. 
     */

  }, {
    key: "getPath",
    value: function getPath() {
      return this.m_path;
    }
    /**
     * The number of polygons in the current corridor path.
     * @return The number of polygons in the current corridor path.
     */

  }, {
    key: "getPathCount",
    value: function getPathCount() {
      return this.m_path.length;
    }
  }]);

  return PathCorridor;
}();

var _default = PathCorridor;
exports["default"] = _default;