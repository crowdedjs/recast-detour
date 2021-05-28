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






import ObstacleAvoidanceQuery from "./ObstacleAvoidanceQuery.js"
import DetourCommon from "../DetourCommon.js"
import ProximityGrid from "./ProximityGrid.js"
import PathQueue from "./PathQueue.js"
import CrowdAgent from "./CrowdAgent.js"
import NavMeshQuery from "../NavMeshQuery.js"
import DefaultQueryFilter from "../DefaultQueryFilter.js"
import PriorityQueue from "./PriorityQueue.js"
import Status from "../Status.js"
import ObstacleAvoidanceParams from "./ObstacleAvoidanceParams.js"
import CrowdAgentParams from "./CrowdAgentParams.js"


class Crowd {

    static MAX_ITERS_PER_UPDATE = 100;

    static MAX_PATHQUEUE_NODES = 4096;
    static MAX_COMMON_NODES = 512;

    /// The maximum number of neighbors that a crowd agent can take into account
    /// for steering decisions.
    /// @ingroup crowd
    static DT_CROWDAGENT_MAX_NEIGHBOURS = 6;

    /// The maximum number of corners a crowd agent will look ahead in the path.
    /// This value is used for sizing the crowd agent corner buffers.
    /// Due to the behavior of the crowd manager, the actual number of useful
    /// corners will be one less than this number.
    /// @ingroup crowd
    static DT_CROWDAGENT_MAX_CORNERS = 4;

    /// The maximum number of crowd avoidance configurations supported by the
    /// crowd manager.
    /// @ingroup crowd
    /// @see dtObstacleAvoidanceParams, dtCrowd::setObstacleAvoidanceParams(), dtCrowd::getObstacleAvoidanceParams(),
    /// dtCrowdAgentParams::obstacleAvoidanceType
    static DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS = 8;

    /// The maximum number of query filter types supported by the crowd manager.
    /// @ingroup crowd
    /// @see dtQueryFilter, dtCrowd::getFilter() dtCrowd::getEditableFilter(),
    /// dtCrowdAgentParams::queryFilterType
    static DT_CROWD_MAX_QUERY_FILTER_TYPE = 16;

    /// Provides neighbor data for agents managed by the crowd.
    /// @ingroup crowd
    /// @see dtCrowdAgent::neis, dtCrowd
    CrowdNeighbor = class CrowdNeighbour {
        idx; /// < The index of the neighbor in the crowd.
        dist; /// < The distance between the current agent and the neighbor.

        constructor(idx, dist) {
            this.idx = idx;
            this.dist = dist;
        }
    };

    m_maxAgents;
    m_agents = [];
    m_activeAgents = [];
    m_pathq;

    m_obstacleQueryParams = new Array(Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS);
    m_obstacleQuery;

    m_grid;

    m_ext = new Array(3);

    m_filters = new Array(Crowd.DT_CROWD_MAX_QUERY_FILTER_TYPE);

    m_maxAgentRadius;

    m_velocitySampleCount;

    m_navquery;

    tween(t, t0, t1) {
        return DetourCommon.clamp((t - t0) / (t1 - t0), 0.0, 1.0);
    }

    getNeighbours(pos, height, range, skip, agents, grid) {

        let result = [];
        let ids = grid.queryItems(pos[0] - range, pos[2] - range, pos[0] + range, pos[2] + range);

        for (let id of ids) {
            let ag = agents[id];

            if (ag == skip || !ag.active) {
                // console.log("Testing");
                continue;
            }

            // Check for overlap.
            let diff = DetourCommon.vSub(pos, ag.npos);
            if (Math.abs(diff[1]) >= (height + ag.params.height) / 2.0) {
                continue;
            }
            diff[1] = 0;
            let distSqr = DetourCommon.vLenSqr(diff);
            if (distSqr > DetourCommon.sqr(range)) {
                continue;
            }

            this.addNeighbour(id, distSqr, result);
        }
        return result;

    }

    addNeighbour(idx, dist, neis) {
        // Insert neighbour based on the distance.
        let nei = new this.CrowdNeighbor(idx, dist);
        neis.push(nei);
        neis.sort((o1, o2) => o1.dist - o2.dist);
    }

    addToOptQueue(newag, agents) {
        // Insert neighbour based on greatest time.
        agents.push(newag);
    }

    // Insert neighbour based on greatest time.
    addToPathQueue(newag, agents) {
        agents.push(newag);
    }

    ///
    /// Initializes the crowd.
    /// May be called more than once to purge and re-initialize the crowd.
    /// @param[in] maxAgents The maximum number of agents the crowd can manage. [Limit: >= 1]
    /// @param[in] maxAgentRadius The maximum radius of any agent that will be added to the crowd. [Limit: > 0]
    /// @param[in] nav The navigation mesh to use for planning.
    /// @return True if the initialization succeeded.


    constructor(maxAgents, maxAgentRadius, nav, queryFilterFactory = (i => new DefaultQueryFilter())) {

        this.m_maxAgents = maxAgents;
        this.m_maxAgentRadius = maxAgentRadius;
        DetourCommon.vSet(this.m_ext, this.m_maxAgentRadius * 2.0, this.m_maxAgentRadius * 1.5, this.m_maxAgentRadius * 2.0);

        this.m_grid = new ProximityGrid(this.m_maxAgents * 4, this.m_maxAgentRadius * 3);
        this.m_obstacleQuery = new ObstacleAvoidanceQuery(6, 8);

        for (let i = 0; i < Crowd.DT_CROWD_MAX_QUERY_FILTER_TYPE; i++) {
            this.m_filters[i] = queryFilterFactory.apply(i);
        }
        // Init obstacle query params.
        for (let i = 0; i < Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS; ++i) {
            this.m_obstacleQueryParams[i] = new ObstacleAvoidanceParams();
        }

        // Allocate temp buffer for merging paths.
        this.m_pathq = new PathQueue(Crowd.MAX_PATHQUEUE_NODES, nav);
        this.m_agents = new Array(this.m_maxAgents);
        this.m_activeAgents = [];
        for (let i = 0; i < this.m_maxAgents; ++i) {
            this.m_agents[i] = new CrowdAgent(i);
            this.m_agents[i].active = false;
        }

        // The navquery is mostly used for local searches, no need for large
        // node pool.
        this.m_navquery = new NavMeshQuery(nav);
    }

    /// Sets the shared avoidance configuration for the specified index.
    /// @param[in] idx The index. [Limits: 0 <= value <
    /// #DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS]
    /// @param[in] params The new configuration.
    setObstacleAvoidanceParams(idx, params) {
        if (idx >= 0 && idx < Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS) {
            this.m_obstacleQueryParams[idx] = params;
        }
    }

    /// Gets the shared avoidance configuration for the specified index.
    /// @param[in] idx The index of the configuration to retreive.
    /// [Limits: 0 <= value < #DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS]
    /// @return The requested configuration.
    getObstacleAvoidanceParams(idx) {
        if (idx >= 0 && idx < Crowd.DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS) {
            return this.m_obstacleQueryParams[idx];
        }
        return null;
    }

    /// The maximum number of agents that can be managed by the object.
    /// @return The maximum number of agents.
    getAgentCount() {
        return this.m_maxAgents;
    }

    /// Gets the specified agent from the pool.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
    /// @return The requested agent.
    /// Agents in the pool may not be in use. Check #dtCrowdAgent.active before using the returned object.
    getAgent(idx) {
        return idx < 0 || idx >= this.m_agents.length ? null : this.m_agents[idx];
    }

    ///
    /// Gets the specified agent from the pool.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
    /// @return The requested agent.
    /// Agents in the pool may not be in use. Check #dtCrowdAgent.active before using the returned object.
    getEditableAgent(idx) {
        return idx < 0 || idx >= this.m_agents.length ? null : this.m_agents[idx];
    }

    /// Updates the specified agent's configuration.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
    /// @param[in] params The new agent configuration.
    updateAgentParameters(idx, params) {
        if (idx < 0 || idx >= this.m_maxAgents) {
            return;
        }
        this.m_agents[idx].params = params;
    }

    /// Adds a new agent to the crowd.
    /// @param[in] pos The requested position of the agent. [(x, y, z)]
    /// @param[in] params The configutation of the agent.
    /// @return The index of the agent in the agent pool. Or -1 if the agent
    /// could not be added.
    addAgent(pos, params) {
        // Find empty slot.
        let idx = -1;
        for (let i = 0; i < this.m_maxAgents; ++i) {
            if (!this.m_agents[i].active) {
                idx = i;
                break;
            }
        }
        if (idx == -1) {
            return -1;
        }

        let ag = this.m_agents[idx];

        this.updateAgentParameters(idx, params);

        // Find nearest position on navmesh and place the agent there.
        let nearest = this.m_navquery.findNearestPoly(pos, this.m_ext, this.m_filters[ag.params.queryFilterType]);

        ag.corridor.reset(nearest.getNearestRef(), nearest.getNearestPos());
        ag.boundary.reset();
        ag.partial = false;

        ag.topologyOptTime = 0;
        ag.targetReplanTime = 0;

        DetourCommon.vSet(ag.dvel, 0, 0, 0);
        DetourCommon.vSet(ag.nvel, 0, 0, 0);
        DetourCommon.vSet(ag.vel, 0, 0, 0);
        DetourCommon.vCopy(ag.npos, nearest.getNearestPos());

        ag.desiredSpeed = 0;

        if (nearest.getNearestRef() != 0) {
            ag.state = CrowdAgent.DT_CROWDAGENT_STATE_WALKING;
        } else {
            ag.state = CrowdAgent.DT_CROWDAGENT_STATE_INVALID;
        }

        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_NONE;

        ag.active = true;

        return idx;
    }

    /// Removes the agent from the crowd.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
    ///
    /// The agent is deactivated and will no longer be processed. Its
    /// #dt object
    /// is not removed from the pool. It is marked as inactive so that it is
    /// available for reuse.
    /// Removes the agent from the crowd.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
    removeAgent(idx) {
        if (idx >= 0 && idx < this.m_maxAgents) {
            this.m_agents[idx].active = false;
        }
    }

    requestMoveTargetReplan(ag, ref, pos) {
        ag.setTarget(ref, pos);
        ag.targetReplan = true;
        return true;
    }

    /// Submits a new move request for the specified agent.
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
    requestMoveTarget(idx, ref, pos) {
        if (idx < 0 || idx >= this.m_maxAgents) {
            return false;
        }
        if (ref == 0) {
            return false;
        }

        let ag = this.m_agents[idx];

        // Initialize request.
        ag.setTarget(ref, pos);
        ag.targetReplan = false;

        return true;
    }

    /// Submits a new move request for the specified agent.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
    /// @param[in] vel The movement velocity. [(x, y, z)]
    /// @return True if the request was successfully submitted.
    requestMoveVelocity(idx, vel) {
        if (idx < 0 || idx >= this.m_maxAgents) {
            return false;
        }

        ag = this.m_agents[idx];

        // Initialize request.
        ag.targetRef = 0;
        DetourCommon.vCopy(ag.targetPos, vel);
        ag.targetPathqRef = PathQueue.DT_PATHQ_INVALID;
        ag.targetReplan = false;
        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY;

        return true;
    }

    /// Resets any request for the specified agent.
    /// @param[in] idx The agent index. [Limits: 0 <= value < #getAgentCount()]
    /// @return True if the request was successfully reseted.
    resetMoveTarget(idx) {
        if (idx < 0 || idx >= this.m_maxAgents) {
            return false;
        }

        ag = this.m_agents[idx];

        // Initialize request.
        ag.targetRef = 0;
        DetourCommon.vSet(ag.targetPos, 0, 0, 0);
        DetourCommon.vSet(ag.dvel, 0, 0, 0);
        ag.targetPathqRef = PathQueue.DT_PATHQ_INVALID;
        ag.targetReplan = false;
        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_NONE;
        return true;
    }

    /// Gets the active agents let the agent pool.
    /// @param[out] agents An array of agent pointers. [(#dt *) * maxAgents]
    /// @param[in] maxAgents The size of the crowd agent array.
    /// @return The number of agents returned in @p agents.
    getActiveAgents() {
        let agents = []//new Array(this.m_maxAgents);
        for (let i = 0; i < this.m_maxAgents; ++i) {
            if (this.m_agents[i].active) {
                agents.push(this.m_agents[i]);
            }
        }
        return agents;
    }

    static MAX_ITER = 20;

    updateMoveRequest() {
        let queue = new PriorityQueue((a1, a2) => a2.targetReplanTime - a1.targetReplanTime);

        // Fire off new requests.
        for (let i = 0; i < this.m_maxAgents; ++i) {
            // if(i==12)
            //     console.log("Bad agent.")
            let ag = this.m_agents[i];
            if (!ag.active) {
                continue;
            }
            if (ag.state == CrowdAgent.DT_CROWDAGENT_STATE_INVALID) {
                continue;
            }
            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE
                || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
                continue;
            }

            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_REQUESTING) {
                let path = ag.corridor.getPath();
                if (path.length == 0) {
                    throw new IllegalArgumentException("Empty path");
                }
                // Quick search towards the goal.
                this.m_navquery.initSlicedFindPath(path[0], ag.targetRef, ag.npos, ag.targetPos, this.m_filters[ag.params.queryFilterType], 0);
                this.m_navquery.updateSlicedFindPath(Crowd.MAX_ITER);
                let pathFound;
                if (ag.targetReplan) // && npath > 10)
                {
                    // Try to use existing steady path during replan if
                    // possible.
                    pathFound = this.m_navquery.finalizeSlicedFindPathPartial(path);
                } else {
                    // Try to move towards target when goal changes.
                    pathFound = this.m_navquery.finalizeSlicedFindPath();
                }
                let reqPath = pathFound.getRefs();
                let reqPos = new Array(3);
                if (!pathFound.getStatus() == Status.FAILURE && reqPath.length > 0) {
                    // In progress or succeed.
                    if (reqPath[reqPath.length - 1] != ag.targetRef) {
                        // Partial path, constrain target position inside the
                        // last polygon.
                        let cr = this.m_navquery.closestPointOnPoly(reqPath[reqPath.length - 1], ag.targetPos);
                        reqPos = cr.getClosest();
                    } else {
                        DetourCommon.vCopy(reqPos, ag.targetPos);
                    }
                } else {
                    // Could not find path, start the request from current
                    // location.
                    DetourCommon.vCopy(reqPos, ag.npos);
                    reqPath = [];
                    reqPath.push(path[0]);
                }

                ag.corridor.setCorridor(reqPos, reqPath);
                ag.boundary.reset();
                ag.partial = false;

                if (reqPath[reqPath.length - 1] == ag.targetRef) {
                    ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_VALID;
                    ag.targetReplanTime = 0.0;
                } else {
                    // The path is longer or potentially unreachable, full plan.
                    ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE;
                }
            }

            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE) {
                this.addToPathQueue(ag, queue);
            }
        }

        while (!queue.isEmpty()) {
            let ag = queue.poll();
            ag.targetPathqRef = this.m_pathq.request(ag.corridor.getLastPoly(), ag.targetRef, ag.corridor.getTarget(), ag.targetPos,
                this.m_filters[ag.params.queryFilterType]);
            if (ag.targetPathqRef != PathQueue.DT_PATHQ_INVALID) {
                ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_WAITING_FOR_PATH;
            }
        }

        // Update requests.
        this.m_pathq.update(Crowd.MAX_ITERS_PER_UPDATE);

        // Process path results.
        for (let i = 0; i < this.m_maxAgents; ++i) {
            let ag = this.m_agents[i];
            if (!ag.active) {
                continue;
            }
            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE
                || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
                continue;
            }

            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_WAITING_FOR_PATH) {
                // Poll path queue.
                let status = this.m_pathq.getRequestStatus(ag.targetPathqRef);
                if (status != null && status == Status.FAILURE) {
                    // Path find failed, retry if the target location is still
                    // valid.
                    ag.targetPathqRef = PathQueue.DT_PATHQ_INVALID;
                    if (ag.targetRef != 0) {
                        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_REQUESTING;
                    } else {
                        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_FAILED;
                    }
                    ag.targetReplanTime = 0.0;
                } else if (status != null && (status == Status.SUCCSESS || status == Status.PARTIAL_RESULT)) {
                    let path = ag.corridor.getPath();
                    if (path.length == 0) {
                        throw new IllegalArgumentException("Empty path");
                    }

                    // Apply results.
                    let targetPos = ag.targetPos;

                    let valid = true;
                    let pathFound = this.m_pathq.getPathResult(ag.targetPathqRef);
                    let res = pathFound.getRefs();
                    status = pathFound.getStatus();
                    if (status == Status.FAILURE || res.length == 0) {
                        valid = false;
                    }

                    if (status != null & status == Status.PARTIAL_RESULT) {
                        ag.partial = true;
                    } else {
                        ag.partial = false;
                    }

                    // Merge result and existing path.
                    // The agent might have moved whilst the request is
                    // being processed, so the path may have changed.
                    // We assume that the end of the path is at the same
                    // location
                    // where the request was issued.

                    // The last ref in the old path should be the same as
                    // the location where the request was issued..
                    if (valid && path[path.length - 1] != res[0]) {
                    // if (valid && path[path.length - 1].longValue() != res[0].longValue()) {
                        valid = false;
                    }

                    if (valid) {
                        // Put the old path infront of the old path.
                        if (path.length > 1) {
                            // path.remove(path.length - 1);
                            path.splice(path.length - 1, 1);
                            path.push(...res);
                            res = path;
                            // Remove trackbacks
                            for (let j = 1; j < res.length - 1; ++j) {
                                if (j - 1 >= 0 && j + 1 < res.length) {
                                    if (res[j - 1] == res[j + 1]) {
                                        res.splice(j + 1,1);
                                        res.splice(j,1);
                                        j -= 2;
                                    }
                                }
                            }
                        }

                        // Check for partial path.
                        if (res[res.length - 1] != ag.targetRef) {
                            // Partial path, constrain target position inside
                            // the last polygon.
                            let cr = this.m_navquery.closestPointOnPoly(res[res.length - 1], targetPos);
                            targetPos = cr.getClosest();
                        }
                    }

                    if (valid) {
                        // Set current corridor.
                        ag.corridor.setCorridor(targetPos, res);
                        // Force to update boundary.
                        ag.boundary.reset();
                        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_VALID;
                    } else {
                        // Something went wrong.
                        ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_FAILED;
                    }

                    ag.targetReplanTime = 0.0;
                }
            }
        }
    }

    static OPT_TIME_THR = 0.5; // seconds

    updateTopologyOptimization(agents, dt) {
        if (!agents.length == 0) {
            return;
        }

        let queue = new PriorityQueue((a1, a2) => Float.compare(a2.topologyOptTime, a1.topologyOptTime));

        for (let i = 0; i < agents.length; ++i) {
            ag = agents[i];
            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
                continue;
            }
            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE
                || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
                continue;
            }
            if ((ag.params.updateFlags & CrowdAgentParams.DT_CROWD_OPTIMIZE_TOPO) == 0) {
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

    static CHECK_LOOKAHEAD = 10;
    static TARGET_REPLAN_DELAY = 1.0; // seconds

    checkPathValidity(agents, dt) {

        for (let i = 0; i < agents.length; ++i) {
            let ag = agents[i];

            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
                continue;
            }

            ag.targetReplanTime += dt;

            let replan = false;

            // First check that the current location is valid.
            let agentPos = new Array(3);
            let agentRef = ag.corridor.getFirstPoly();
            DetourCommon.vCopy(agentPos, ag.npos);
            if (!this.m_navquery.isValidPolyRef(agentRef, this.m_filters[ag.params.queryFilterType])) {
                // Current location is not valid, try to reposition.
                // TODO: this can snap agents, how to handle that?
                let fnp = this.m_navquery.findNearestPoly(ag.npos, this.m_ext, this.m_filters[ag.params.queryFilterType]);
                agentRef = fnp.getNearestRef();
                if (fnp.getNearestPos() != null) {
                    DetourCommon.vCopy(agentPos, fnp.getNearestPos());
                }

                if (agentRef == 0) {
                    // Could not find location in navmesh, set state to invalid.
                    ag.corridor.reset(0, agentPos);
                    ag.partial = false;
                    ag.boundary.reset();
                    ag.state = CrowdAgent.DT_CROWDAGENT_STATE_INVALID;
                    continue;
                }

                // Make sure the first polygon is valid, but leave other valid
                // polygons in the path so that replanner can adjust the path
                // better.
                ag.corridor.fixPathStart(agentRef, agentPos);
                // ag.corridor.trimInvalidPath(agentRef, agentPos, m_navquery,
                // &m_filter);
                ag.boundary.reset();
                DetourCommon.vCopy(ag.npos, agentPos);

                replan = true;
            }

            // If the agent does not have move target or is controlled by
            // velocity, no need to recover the target nor replan.
            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE
                || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
                continue;
            }

            // Try to recover move request position.
            if (ag.targetState != CrowdAgent.DT_CROWDAGENT_TARGET_NONE
                && ag.targetState != CrowdAgent.DT_CROWDAGENT_TARGET_FAILED) {
                if (!this.m_navquery.isValidPolyRef(ag.targetRef, this.m_filters[ag.params.queryFilterType])) {
                    // Current target is not valid, try to reposition.
                    let fnp = this.m_navquery.findNearestPoly(ag.targetPos, this.m_ext, this.m_filters[ag.params.queryFilterType]);
                    ag.targetRef = fnp.getNearestRef();
                    if (fnp.getNearestPos() != null) {
                        DetourCommon.vCopy(ag.targetPos, fnp.getNearestPos());
                    }
                    replan = true;
                }
                if (ag.targetRef == 0) {
                    // Failed to reposition target, fail moverequest.
                    ag.corridor.reset(agentRef, agentPos);
                    ag.partial = false;
                    ag.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_NONE;
                }
            }

            // If nearby corridor is not valid, replan.
            if (!ag.corridor.isValid(Crowd.CHECK_LOOKAHEAD, this.m_navquery, this.m_filters[ag.params.queryFilterType])) {
                // Fix current path.
                // ag.corridor.trimInvalidPath(agentRef, agentPos, m_navquery,
                // &m_filter);
                // ag.boundary.reset();
                replan = true;
            }

            // If the end of the path is near and it is not the requested
            // location, replan.
            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VALID) {
                if (ag.targetReplanTime > Crowd.TARGET_REPLAN_DELAY && ag.corridor.getPathCount() < Crowd.CHECK_LOOKAHEAD
                    && ag.corridor.getLastPoly() != ag.targetRef) {
                    replan = true;
                }
            }

            // Try to replan path to goal.
            if (replan) {
                if (ag.targetState != CrowdAgent.DT_CROWDAGENT_TARGET_NONE) {
                    this.requestMoveTargetReplan(ag, ag.targetRef, ag.targetPos);
                }
            }
        }
    }

    static COLLISION_RESOLVE_FACTOR = 0.7;

    update(dt, debug, frame) {
        // if (frame == 114)
        //     console.log("Bug");
        this.m_velocitySampleCount = 0;

        let debugIdx = debug != null ? debug.idx : -1;

        let agents = this.getActiveAgents();

        // Check that all agents still have valid paths.
        this.checkPathValidity(agents, dt);

        // Update async move request and path finder.
        this.updateMoveRequest();

        // Optimize path topology.
        this.updateTopologyOptimization(agents, dt);

        // Register agents to proximity grid.
        this.m_grid.clear();
        for (let i = 0; i < agents.length; ++i) {
            let ag = agents[i];
            let p = ag.npos;
            let r = ag.params.radius;
            this.m_grid.addItem(i, p[0] - r, p[2] - r, p[0] + r, p[2] + r);
        }

        // Get nearby navmesh segments and agents to collide with.
        for (let ag of agents) {
            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
                continue;
            }

            // Update the collision boundary after certain distance has been passed or
            // if it has become invalid.
            let updateThr = ag.params.collisionQueryRange * 0.25;
            if (DetourCommon.vDist2DSqr(ag.npos, ag.boundary.getCenter()) > DetourCommon.sqr(updateThr)
                || !ag.boundary.isValid(this.m_navquery, this.m_filters[ag.params.queryFilterType])) {
                ag.boundary.update(ag.corridor.getFirstPoly(), ag.npos, ag.params.collisionQueryRange, this.m_navquery,
                    this.m_filters[ag.params.queryFilterType]);
            }
            // Query neighbour agents
            ag.neis = this.getNeighbours(ag.npos, ag.params.height, ag.params.collisionQueryRange, ag, agents, this.m_grid);
        }

        // Find next corner to steer to.
        for (let i = 0; i < agents.length; ++i) {
            let ag = agents[i];

            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
                continue;
            }
            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE
                || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
                continue;
            }

            // Find corners for steering
            ag.corners = ag.corridor.findCorners(Crowd.DT_CROWDAGENT_MAX_CORNERS, this.m_navquery, this.m_filters[ag.params.queryFilterType]);

            // Check to see if the corner after the next corner is directly visible,
            // and short cut to there.
            if ((ag.params.updateFlags & CrowdAgentParams.DT_CROWD_OPTIMIZE_VIS) != 0 && ag.corners.length > 0) {
                let target = ag.corners[Math.min(1, ag.corners.length - 1)].getPos();
                ag.corridor.optimizePathVisibility(target, ag.params.pathOptimizationRange, this.m_navquery,
                    this.m_filters[ag.params.queryFilterType]);

                // Copy data for debug purposes.
                if (debugIdx == i) {
                    DetourCommon.vCopy(debug.optStart, ag.corridor.getPos());
                    DetourCommon.vCopy(debug.optEnd, target);
                }
            } else {
                // Copy data for debug purposes.
                if (debugIdx == i) {
                    DetourCommon.vSet(debug.optStart, 0, 0, 0);
                    DetourCommon.vSet(debug.optEnd, 0, 0, 0);
                }
            }
        }

        // Trigger off-mesh connections (depends on corners).
        for (let ag of agents) {

            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
                continue;
            }
            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE
                || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
                continue;
            }

            // Check
            let triggerRadius = ag.params.radius * 2.25;
            if (ag.overOffmeshConnection(triggerRadius)) {
                // Prepare to off-mesh connection.
                anim = ag.animation;

                // Adjust the path over the off-mesh connection.
                let refs = new long[2];
                if (ag.corridor.moveOverOffmeshConnection(ag.corners[ag.corners.length - 1].getRef(), refs, anim.startPos, anim.endPos,
                    this.m_navquery)) {
                    DetourCommon.vCopy(anim.initPos, ag.npos);
                    anim.polyRef = refs[1];
                    anim.active = true;
                    anim.t = 0.0;
                    anim.tmax = (DetourCommon.vDist2D(anim.startPos, anim.endPos) / ag.params.maxSpeed) * 0.5;

                    ag.state = CrowdAgent.DT_CROWDAGENT_STATE_OFFMESH;
                    ag.corners = [];
                    ag.neis = [];
                    continue;
                } else {
                    // Path validity check will ensure that bad/blocked connections will be replanned.
                }
            }
        }

        // Calculate steering.
        for (let ag of agents) {

            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
                continue;
            }
            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE) {
                continue;
            }

            let dvel = new Array(3);

            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
                DetourCommon.vCopy(dvel, ag.targetPos);
                ag.desiredSpeed = DetourCommon.vLen(ag.targetPos);
            } else {
                // Calculate steering direction.
                if ((ag.params.updateFlags & CrowdAgentParams.DT_CROWD_ANTICIPATE_TURNS) != 0) {
                    dvel = ag.calcSmoothSteerDirection();
                } else {
                    dvel = ag.calcStraightSteerDirection();
                }
                // Calculate speed scale, which tells the agent to slowdown at the end of the path.
                let slowDownRadius = ag.params.radius * 2; // TODO: make less hacky.
                let speedScale = ag.getDistanceToGoal(slowDownRadius) / slowDownRadius;

                ag.desiredSpeed = ag.params.maxSpeed;
                dvel = DetourCommon.vScale(dvel, ag.desiredSpeed * speedScale);
            }

            // Separation
            if ((ag.params.updateFlags & CrowdAgentParams.DT_CROWD_SEPARATION) != 0) {
                separationDist = ag.params.collisionQueryRange;
                invSeparationDist = 1.0 / separationDist;
                separationWeight = ag.params.separationWeight;

                w = 0;
                let disp = new Array(3);

                for (let j = 0; j < ag.neis.length; ++j) {
                    nei = agents.get(ag.neis[j].idx);

                    let diff = DetourCommon.vSub(ag.npos, nei.npos);
                    diff[1] = 0;

                    distSqr = DetourCommon.vLenSqr(diff);
                    if (distSqr < 0.00001) {
                        continue;
                    }
                    if (distSqr > DetourCommon.sqr(separationDist)) {
                        continue;
                    }
                    dist = Math.sqrt(distSqr);
                    weight = separationWeight * (1.0 - DetourCommon.sqr(dist * invSeparationDist));

                    disp = DetourCommon.vMad(disp, diff, weight / dist);
                    w += 1.0;
                }

                if (w > 0.0001) {
                    // Adjust desired velocity.
                    dvel = DetourCommon.vMad(dvel, disp, 1.0 / w);
                    // Clamp desired velocity to desired speed.
                    speedSqr = DetourCommon.vLenSqr(dvel);
                    desiredSqr = DetourCommon.sqr(ag.desiredSpeed);
                    if (speedSqr > desiredSqr) {
                        dvel = DetourCommon.vScale(dvel, desiredSqr / speedSqr);
                    }
                }
            }

            // Set the desired velocity.
            DetourCommon.vCopy(ag.dvel, dvel);
        }

        // Velocity planning.
        for (let i = 0; i < agents.length; ++i) {
            let ag = agents[i];

            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
                continue;
            }

            if ((ag.params.updateFlags & CrowdAgentParams.DT_CROWD_OBSTACLE_AVOIDANCE) != 0) {
                this.m_obstacleQuery.reset();

                // Add neighbours as obstacles.
                for (let j = 0; j < ag.neis.length; ++j) {
                    let nei = agents[ag.neis[j].idx];
                    this.m_obstacleQuery.addCircle(nei.npos, nei.params.radius, nei.vel, nei.dvel);
                }

                // if (ag.neis.length > 0 && i == 0) {
                //     console.log("Frame " + frame)
                // }

                // Append neighbour segments as obstacles.
                for (let j = 0; j < ag.boundary.getSegmentCount(); ++j) {
                    let s = ag.boundary.getSegment(j);
                    //let s3 = Arrays.copyOfRange(s, 3, 6);
                    //let s3 = Arrays.copyOfRange(s, 3, 6);
                    let s3 = s.slice(3,6);
                    if (DetourCommon.triArea2D3(ag.npos, s, s3) < 0.0) {
                        continue;
                    }
                    this.m_obstacleQuery.addSegment(s, s3);
                }

                let vod = null;
                if (debugIdx == i) {
                    vod = debug.vod;
                }

                // Sample new safe velocity.
                let adaptive = true;
                let ns = 0;

                let params = this.m_obstacleQueryParams[ag.params.obstacleAvoidanceType];

                if (adaptive) {
                    let nsnvel = this.m_obstacleQuery.sampleVelocityAdaptive(ag.npos, ag.params.radius, ag.desiredSpeed,
                        ag.vel, ag.dvel, params, vod);
                    ns = nsnvel[0];
                    ag.nvel = nsnvel[1];
                } else {
                    let nsnvel = this.m_obstacleQuery.sampleVelocityGrid(ag.npos, ag.params.radius, ag.desiredSpeed,
                        ag.vel, ag.dvel, params, vod);
                    ns = nsnvel[0];
                    ag.nvel = nsnvel[1];
                }
                this.m_velocitySampleCount += ns;
            } else {
                // If not using velocity planning, new velocity is directly the desired velocity.
                DetourCommon.vCopy(ag.nvel, ag.dvel);
            }
        }

        // Integrate.
        for (let i = 0; i < agents.length; ++i) {
            let ag = agents[i];
            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
                continue;
            }
            ag.integrate(dt);
        }

        // Handle collisions.

        for (let iter = 0; iter < 4; ++iter) {
            for (let i = 0; i < agents.length; ++i) {
                let ag = agents[i];
                let idx0 = ag.getAgentIndex();
                if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
                    continue;
                }

                DetourCommon.vSet(ag.disp, 0, 0, 0);

                let w = 0;

                for (let j = 0; j < ag.neis.length; ++j) {
                    let nei = agents[ag.neis[j].idx];
                    let idx1 = nei.getAgentIndex();
                    let diff = DetourCommon.vSub(ag.npos, nei.npos);
                    diff[1] = 0;

                    let dist = DetourCommon.vLenSqr(diff);
                    if (dist > DetourCommon.sqr(ag.params.radius + nei.params.radius)) {
                        continue;
                    }
                    dist = Math.sqrt(dist);
                    let pen = (ag.params.radius + nei.params.radius) - dist;
                    if (dist < 0.0001) {
                        // Agents on top of each other, try to choose diverging separation directions.
                        if (idx0 > idx1) {
                            DetourCommon.vSet(diff, -ag.dvel[2], 0, ag.dvel[0]);
                        } else {
                            DetourCommon.vSet(diff, ag.dvel[2], 0, -ag.dvel[0]);
                        }
                        pen = 0.01;
                    } else {
                        pen = (1.0 / dist) * (pen * 0.5) * Crowd.COLLISION_RESOLVE_FACTOR;
                    }

                    ag.disp = DetourCommon.vMad(ag.disp, diff, pen);

                    w += 1.0;
                }

                if (w > 0.0001) {
                    let iw = 1.0 / w;
                    ag.disp = DetourCommon.vScale(ag.disp, iw);
                }
            }

            for (let i = 0; i < agents.length; ++i) {
                let ag = agents[i];
                if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
                    continue;
                }

                ag.npos = DetourCommon.vAdd(ag.npos, ag.disp);
            }
        }

        for (let i = 0; i < agents.length; ++i) {
            // if (frame == 492)
            //     console.log("Bad agent")
            let ag = agents[i];
            if (ag.state != CrowdAgent.DT_CROWDAGENT_STATE_WALKING) {
                continue;
            }

            // Move along navmesh.
            ag.corridor.movePosition(ag.npos, this.m_navquery, this.m_filters[ag.params.queryFilterType]);
            // Get valid constrained position back.
            DetourCommon.vCopy(ag.npos, ag.corridor.getPos());

            // If not using path, truncate the corridor to just one poly.
            if (ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_NONE
                || ag.targetState == CrowdAgent.DT_CROWDAGENT_TARGET_VELOCITY) {
                ag.corridor.reset(ag.corridor.getFirstPoly(), ag.npos);
                ag.partial = false;
            }

        }

        // Update agents using off-mesh connection.
        for (let i = 0; i < this.m_maxAgents; ++i) {
            let anim = this.m_agents[i].animation;
            if (!anim.active) {
                continue;
            }
            ag = this.m_agents[i];

            anim.t += dt;
            if (anim.t > anim.tmax) {
                // Reset animation
                anim.active = false;
                // Prepare agent for walking.
                ag.state = CrowdAgent.DT_CROWDAGENT_STATE_WALKING;
                continue;
            }

            // Update position
            ta = anim.tmax * 0.15;
            tb = anim.tmax;
            if (anim.t < ta) {
                u = tween(anim.t, 0.0, ta);
                ag.npos = DetourCommon.vLerp3(anim.initPos, anim.startPos, u);
            } else {
                u = tween(anim.t, ta, tb);
                ag.npos = DetourCommon.vLerp3(anim.startPos, anim.endPos, u);
            }

            // Update velocity.
            DetourCommon.vSet(ag.vel, 0, 0, 0);
            DetourCommon.vSet(ag.dvel, 0, 0, 0);
        }
    }

    getQueryExtents() {
        return this.m_ext;
    }

    getFilter(i) {
        return i >= 0 && i < Crowd.DT_CROWD_MAX_QUERY_FILTER_TYPE ? this.m_filters[i] : null;
    }

}

export default Crowd;
