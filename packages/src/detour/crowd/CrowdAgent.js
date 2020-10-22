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
import PathCorridor from "./PathCorridor.js"
import LocalBoundary from "./LocalBoundary.js"
import CrowdAgentAnimation from "./CrowdAgentAnimation.js"
import PathQueue from "./PathQueue.js"


/// Represents an agent managed by a #dt object.
/// @ingroup crowd
class CrowdAgent {

	/// The type of navigation mesh polygon the agent is currently traversing.
	/// @ingroup crowd
	static DT_CROWDAGENT_STATE_INVALID = 0;
	static DT_CROWDAGENT_STATE_WALKING = 1;
	static DT_CROWDAGENT_STATE_OFFMESH = 2;



	static DT_CROWDAGENT_TARGET_NONE = 0;
	static DT_CROWDAGENT_TARGET_FAILED = 1;
	static DT_CROWDAGENT_TARGET_VALID = 2;
	static DT_CROWDAGENT_TARGET_REQUESTING = 3;
	static DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE = 4;
	static DT_CROWDAGENT_TARGET_WAITING_FOR_PATH = 5;
	static DT_CROWDAGENT_TARGET_VELOCITY = 6;


	idx;

	/// True if the agent is active, false if the agent is in an unused slot in the agent pool.
	active;

	/// The type of mesh polygon the agent is traversing. (See: #CrowdAgent)
	state;

	/// True if the agent has valid path (targetState == DT_CROWDAGENT_TARGET_VALID) and the path does not lead to the requested position, else false.
	partial;

	/// The path corridor the agent is using.
	corridor;

	/// The local boundary data for the agent.
	boundary;

	/// Time since the agent's path corridor was optimized.
	topologyOptTime;

	/// The known neighbors of the agent.
	neis = [];

	/// The desired speed.
	desiredSpeed;

	npos = new Array(3); ///< The current agent position. [(x, y, z)]
	disp = new Array(3); ///< A temporary value used to accumulate agent displacement during iterative collision resolution. [(x, y, z)]
	dvel = new Array(3); ///< The desired velocity of the agent. Based on the current path, calculated from scratch each frame. [(x, y, z)]
	nvel = new Array(3); ///< The desired velocity adjusted by obstacle avoidance, calculated from scratch each frame. [(x, y, z)]
	vel = new Array(3); ///< The actual velocity of the agent. The change from nvel => vel is constrained by max acceleration. [(x, y, z)]

	/// The agent's configuration parameters.
	params;
	/// The local path corridor corners for the agent.
	corners = [];

	targetState; ///< State of the movement request.
	targetRef; ///< Target polyref of the movement request.
	targetPos = new Array(3); ///< Target position of the movement request (or velocity in case of DT_CROWDAGENT_TARGET_VELOCITY).
	targetPathqReq; ///< Path finder ref.
	targetReplan; ///< Flag indicating that the current path is being replanned.
	targetReplanTime; /// <Time since the agent's target was replanned.

	animation;

	constructor(idx) {
		this.idx = idx;
		this.corridor = new PathCorridor();
		this.boundary = new LocalBoundary();
		this.animation = new CrowdAgentAnimation();
	}

	integrate(dt) {
		// Fake dynamic constraint.
		let maxDelta = this.params.maxAcceleration * dt;
		let dv = DetourCommon.vSub(this.nvel, this.vel);
		let ds = DetourCommon.vLen(dv);
		if (ds > maxDelta)
			dv = DetourCommon.vScale(dv, maxDelta / ds);
		this.vel = DetourCommon.vAdd(this.vel, dv);

		// Integrate
		if (DetourCommon.vLen(this.vel) > 0.0001)
			this.npos = DetourCommon.vMad(this.npos, this.vel, dt);
		else
			DetourCommon.vSet(this.vel, 0, 0, 0);
	}

	overOffmeshConnection(radius) {
		if (this.corners.length == 0)
			return false;

		let offMeshConnection = ((this.corners[this.corners.length - 1].getFlags() & NavMeshQuery.DT_STRAIGHTPATH_OFFMESH_CONNECTION) != 0)
			? true : false;
		if (offMeshConnection) {
			distSq = DetourCommon.vDist2D(this.npos, this.corners[this.corners.length - 1].getPos());
			if (distSq < radius * radius)
				return true;
		}

		return false;
	}

	getDistanceToGoal(range) {
		if (this.corners.length == 0)
			return range;

		let endOfPath = ((this.corners[this.corners.length - 1].getFlags() & NavMeshQuery.DT_STRAIGHTPATH_END) != 0) ? true : false;
		if (endOfPath)
			return Math.min(DetourCommon.vDist2D(this.npos, this.corners[this.corners.length - 1].getPos()), range);

		return range;
	}

	calcSmoothSteerDirection() {
		let dir = new Array(3);
		if (!this.corners.length == 0) {

			let ip0 = 0;
			let ip1 = Math.min(1, this.corners.length - 1);
			let p0 = this.corners[ip0].getPos();
			let p1 = this.corners[ip1].getPos();

			let dir0 = DetourCommon.vSub(p0, this.npos);
			let dir1 = DetourCommon.vSub(p1, this.npos);
			dir0[1] = 0;
			dir1[1] = 0;

			let len0 = DetourCommon.vLen(dir0);
			let len1 = DetourCommon.vLen(dir1);
			if (len1 > 0.001)
				dir1 = DetourCommon.vScale(dir1, 1.0 / len1);

			dir[0] = dir0[0] - dir1[0] * len0 * 0.5;
			dir[1] = 0;
			dir[2] = dir0[2] - dir1[2] * len0 * 0.5;

			DetourCommon.vNormalize(dir);
		}
		return dir;
	}

	calcStraightSteerDirection() {
		let dir = new Array(3);
		if (!this.corners.length == 0) {
			dir = DetourCommon.vSub(this.corners[0].getPos(), this.npos);
			dir[1] = 0;
			DetourCommon.vNormalize(dir);
		}
		return dir;
	}


	setTarget(ref, pos) {
		this.targetRef = ref;
		DetourCommon.vCopy(this.targetPos, pos);
		this.targetPathqRef = PathQueue.DT_PATHQ_INVALID;
		if (this.targetRef != 0)
			this.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_REQUESTING;
		else
			this.targetState = CrowdAgent.DT_CROWDAGENT_TARGET_FAILED;
	}

	getAgentIndex() {
		return this.idx;
	}

	isActive() {
		return this.active;
	}
}

export default CrowdAgent;