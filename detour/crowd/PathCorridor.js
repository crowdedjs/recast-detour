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

import DetourCommon from "../DetourCommon.js"
import NavMeshQuery from "../NavMeshQuery.js"


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
class PathCorridor {

	m_pos = new Array(3);
	m_target = new Array(3);
	m_path = [];

	mergeCorridorStartMoved(path, visited) {
		let furthestPath = -1;
		let furthestVisited = -1;

		// Find furthest common polygon.
		for (let i = path.length - 1; i >= 0; --i) {
			let found = false;
			for (let j = visited.length - 1; j >= 0; --j) {
				if (path[i] == visited[j]) {
					furthestPath = i;
					furthestVisited = j;
					found = true;
				}
			}
			if (found)
				break;
		}

		// If no intersection found just return current path.
		if (furthestPath == -1 || furthestVisited == -1)
			return path;

		// Concatenate paths.

		// Adjust beginning of the buffer to include the visited.
		let result = [];
		// Store visited
		for (let i = visited.length - 1; i > furthestVisited; --i)
			result.push(visited[i]);
		result.push(...path.slice(furthestPath, path.length));
		return result;
	}

	mergeCorridorEndMoved(path, visited) {
		let furthestPath = -1;
		let furthestVisited = -1;

		// Find furthest common polygon.
		for (let i = 0; i < path.length; ++i) {
			let found = false;
			for (let j = visited.length - 1; j >= 0; --j) {
				if (path[i] == visited[j]) {
					furthestPath = i;
					furthestVisited = j;
					found = true;
				}
			}
			if (found)
				break;
		}

		// If no intersection found just return current path.
		if (furthestPath == -1 || furthestVisited == -1)
			return path;

		// Concatenate paths.
		let result = path.subList(0, furthestPath);
		result.addAll(visited.subList(furthestVisited, visited.length));
		return result;
	}

	mergeCorridorStartShortcut(path, visited) {

		let furthestPath = -1;
		let furthestVisited = -1;

		// Find furthest common polygon.
		for (let i = path.length - 1; i >= 0; --i) {
			let found = false;
			for (let j = visited.length - 1; j >= 0; --j) {
				if (path[i] == visited[j]) {
					furthestPath = i;
					furthestVisited = j;
					found = true;
				}
			}
			if (found)
				break;
		}

		// If no intersection found just return current path.
		if (furthestPath == -1 || furthestVisited <= 0)
			return path;

		// Concatenate paths.

		// Adjust beginning of the buffer to include the visited.
		let result = visited.subList(0, furthestVisited);
		result.addAll(path.subList(furthestPath, path.length));
		return result;
	}

	/**
	 * Allocates the corridor's path buffer.
	 */
	constructor() {
		this.m_path = [];
	}

	/**
	 * Resets the path corridor to the specified position.
	 * @param ref The polygon reference containing the position.
	 * @param pos The new position in the corridor. [(x, y, z)]
	 */
	reset(ref, pos) {
		this.m_path = [];
		this.m_path.push(ref);
		DetourCommon.vCopy(this.m_pos, pos);
		DetourCommon.vCopy(this.m_target, pos);
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
	findCorners(maxCorners, navquery, filter) {
		let MIN_TARGET_DIST = DetourCommon.sqr(0.01);

		let path = navquery.findStraightPath(this.m_pos, this.m_target, this.m_path, maxCorners, 0);
		// Prune points in the beginning of the path which are too close.

		// for (let iter = path.iterator(); iter.hasNext();) {
		// 	let spi = iter.next();
		// 	if ((spi.getFlags() & NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION) != 0
		// 		|| DetourCommon.vDist2D(spi.getPos(), this.m_pos) > MIN_TARGET_DIST) {
		// 		break;
		// 	}
		// 	iter.remove();
		// }
		// TODO
		let newPath = [];
		let i = 0;
		let l = path.length;
		for (; i < path.length; i++) {
			let spi = path[i];
			if ((spi.getFlags() & NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION) != 0
				|| DetourCommon.vDist2D(spi.getPos(), this.m_pos) > MIN_TARGET_DIST) {
				break;

			}
		}
		let j = i;
		// for(; j < l; j++);
		// {
		// 	newPath.push(path[j])
		// }
		while(j < l){
			newPath.push(path[j])
			j++
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

		optimizePathVisibility(next, pathOptimizationRange, navquery,
			filter) {
			// Clamp the ray to max distance.
			let dist = DetourCommon.vDist2D(this.m_pos, next);

			// If too close to the goal, do not try to optimize.
			if (dist < 0.01)
				return;

			// Overshoot a little. This helps to optimize open fields in tiled
			// meshes.
			dist = Math.min(dist + 0.01, pathOptimizationRange);

			// Adjust ray length.
			let delta = DetourCommon.vSub(next, this.m_pos);
			let goal = DetourCommon.vMad(this.m_pos, delta, pathOptimizationRange / dist);

			let rc = navquery.raycast(this.m_path[0], this.m_pos, goal, filter, 0, 0);
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
		optimizePathTopology(navquery, filter) {
			if (this.m_path.length < 3)
				return false;

			let MAX_ITER = 32;

			navquery.initSlicedFindPath(this.m_path[0], this.m_path[this.m_path.length - 1], this.m_pos, this.m_target, filter, 0);
			navquery.updateSlicedFindPath(MAX_ITER);
			let fpr = navquery.finalizeSlicedFindPathPartial(this.m_path);

			if (fpr.getStatus().isSuccess() && fpr.getRefs().length > 0) {
				this.m_path = mergeCorridorStartShortcut(this.m_path, fpr.getRefs());
				return true;
			}

			return false;
		}

		moveOverOffmeshConnection(offMeshConRef, refs, start, end, navquery) {
			// Advance the path up to and over the off-mesh connection.
			let prevRef = 0;
			let polyRef = this.m_path[0];
			let npos = 0;
			while (npos < this.m_path.length && polyRef != offMeshConRef) {
				prevRef = polyRef;
				polyRef = this.m_path[npos];
				npos++;
			}
			if (npos == this.m_path.length) {
				// Could not find offMeshConRef
				return false;
			}

			// Prune path
			this.m_path = this.m_path.subList(npos, this.m_path.length);
			refs[0] = prevRef;
			refs[1] = polyRef;

			let nav = navquery.getAttachedNavMesh();
			let startEnd = nav.getOffMeshConnectionPolyEndPoints(refs[0], refs[1]);
			DetourCommon.vCopy(this.m_pos, startEnd.second);
			DetourCommon.vCopy(start, startEnd.first);
			DetourCommon.vCopy(end, startEnd.second);
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
		movePosition(npos, navquery, filter) {
			// Move aPoly navmesh and update new position.
			let masResult = navquery.moveAlongSurface(this.m_path[0], this.m_pos, npos, filter);
			this.m_path = this.mergeCorridorStartMoved(this.m_path, masResult.getVisited());
			// Adjust the position to stay on top of the navmesh.
			DetourCommon.vCopy(this.m_pos, masResult.getResultPos());
			try {
				this.m_pos[1] = navquery.getPolyHeight(this.m_path[0], masResult.getResultPos());
			} catch (e) {
				// Silently disregard the returned status of DT_FAILURE | DT_INVALID_PARAM to stay
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
		moveTargetPosition(npos, navquery, filter) {
			// Move aPoly navmesh and update new position.
			let masResult = navquery.moveAlongSurface(this.m_path[this.m_path.length - 1], this.m_target, npos,
				filter);
			this.m_path = mergeCorridorEndMoved(this.m_path, masResult.getVisited());
			// TODO: should we do that?
			// Adjust the position to stay on top of the navmesh.
			/*
			 *  h = m_target[1]; navquery=>getPolyHeight(this.m_path[m_npath-1],
			 * result, &h); result[1] = h;
			 */
			DetourCommon.vCopy(this.m_target, masResult.getResultPos());
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

		setCorridor(target, path) {
			DetourCommon.vCopy(this.m_target, target);
			this.m_path = [...path];
		}

		fixPathStart(safeRef, safePos) {
			DetourCommon.vCopy(this.m_pos, safePos);
			if (this.m_path.length < 3 && this.m_path.length > 0) {
				let p = this.m_path[this.m_path.length - 1];
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

		trimInvalidPath(safeRef, safePos, navquery, filter) {
			// Keep valid path as far as possible.
			let n = 0;
			while (n < this.m_path.length && navquery.isValidPolyRef(this.m_path[n], filter)) {
				n++;
			}

			if (n == 0) {
				// The first polyref is bad, use current safe values.
				DetourCommon.vCopy(this.m_pos, safePos);
				this.m_path = [];
				this.m_path.push(safeRef);
			} else if (n < this.m_path.length) {
				this.m_path = this.m_path.subList(0, n);
				// The path is partially usable.
			}
			// Clamp target pos to last poly
			DetourCommon.vCopy(this.m_target, navquery.closestPointOnPolyBoundary(this.m_path[this.m_path.length - 1], this.m_target));
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
		isValid(maxLookAhead, navquery, filter) {
			// Check that all polygons still pass query filter.
			let n = Math.min(this.m_path.length, maxLookAhead);
			for (let i = 0; i < n; ++i) {
				if (!navquery.isValidPolyRef(this.m_path[i], filter))
					return false;
			}

			return true;
		}

		/**
		 * Gets the current position within the corridor. (In the first polygon.)
		 * @return The current position within the corridor.
		 */
		getPos() { return this.m_pos; }

		/**
		 * Gets the current target within the corridor. (In the last polygon.)
		 * @return The current target within the corridor.
		 */
		getTarget() { return this.m_target; }

		/**
		 * The polygon reference id of the first polygon in the corridor, the polygon containing the position.
		 * @return The polygon reference id of the first polygon in the corridor. (Or zero if there is no path.)
		 */
		getFirstPoly() {
			return this.m_path.length == 0 ? 0 : this.m_path[0];
		}

		/**
		 * The polygon reference id of the last polygon in the corridor, the polygon containing the target.
		 * @return The polygon reference id of the last polygon in the corridor. (Or zero if there is no path.)
		 */
		getLastPoly() { return this.m_path.length == 0 ? 0 : this.m_path[this.m_path.length - 1]; }

		/**
		 * The corridor's path. 
		 */
		getPath() { return this.m_path; }

		/**
		 * The number of polygons in the current corridor path.
		 * @return The number of polygons in the current corridor path.
		 */
		getPathCount() { return this.m_path.length; }
	}

	export default PathCorridor;
