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
class CrowdAgentParams {
	radius = 0;						///<  radius. [Limit: >= 0]
	height = 0;						///<  height. [Limit: > 0]
	maxAcceleration = 0;			///< Maximum allowed acceleration. [Limit: >= 0]
	maxSpeed = 0;						///< Maximum allowed speed. [Limit: >= 0]

	/// Defines how close a collision element must be before it is considered for steering behaviors. [Limits: > 0]
	collisionQueryRange = 0;

	pathOptimizationRange = 0;		///< The path visibility optimization range. [Limit: > 0]

	/// How aggresive the agent manager should be at avoiding collisions with this agent. [Limit: >= 0]
	separationWeight = 0;

	///  agent update flags.
	static DT_CROWD_ANTICIPATE_TURNS = 1;
	static DT_CROWD_OBSTACLE_AVOIDANCE = 2;
	static DT_CROWD_SEPARATION = 4;
	static DT_CROWD_OPTIMIZE_VIS = 8;			///< Use #dtPathCorridor::optimizePathVisibility() to optimize the agent path.
	static DT_CROWD_OPTIMIZE_TOPO = 16;		///< Use dtPathCorridor::optimizePathTopology() to optimize the agent path.

	/// Flags that impact steering behavior. (See: #UpdateFlags)
	updateFlags = 0;

	/// The index of the avoidance configuration to use for the agent. 
	/// [Limits: 0 <= value <= #DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS]
	obstacleAvoidanceType = 0;

	/// The index of the query filter used by this agent.
	queryFilterType = 0;

	/// User defined data attached to the agent.
	userData = null;
}

export default CrowdAgentParams;