import { CrowdAgentParams, RecastTestMeshBuilder, NavMesh, NavMeshQuery, Crowd, ObstacleAvoidanceParams } from "../src/Main.js"

class NodeApp {
  updateFlags = CrowdAgentParams.DT_CROWD_ANTICIPATE_TURNS | CrowdAgentParams.DT_CROWD_OPTIMIZE_VIS
    | CrowdAgentParams.DT_CROWD_OPTIMIZE_TOPO | CrowdAgentParams.DT_CROWD_OBSTACLE_AVOIDANCE;
  query;
  crowd;
  agents = [];
  ext;
  filter;
  ap;
  navmesh;

  constructor(objString, agentsString, ticks) {
    this.agents = [];
    this.objString = objString;
    this.agentsString = agentsString;
    this.ticks = ticks;
    this.result = "";
  }
  go() {
    this.bootMesh(this.objString);

    //Path is the path to the file where we will store our results
    
    let stream = this.agentsString.split('\n');
    stream.forEach((l,index) => l.trim().length > 0 ? this.agents.push(this.newAgent(l,index)) : 0 == 0);

    let currentMillisecond = 0; //The current time
    let millisecondsBetweenFrames = 40; //40ms between frames, or 25fps
    let secondsOfSimulation = this.ticks; //How long should the simulation run? Change this to whatever you want.
    for (let i = 0; i < secondsOfSimulation; i++) {
      
      for (let j = 0; j < this.agents.length; j++) {
        let agent = this.agents[j]; //Grab each agent in the list

        //Ignore agents that have come into the simulation and exited
        if (agent.hasEntered && !agent.inSimulation) continue;

        //See if we need to add the agent to the simulation
        if (!agent.hasEntered && agent.startMSec <= currentMillisecond) {
          let start = this.getStart(agent);//Get the agent's starting point as a  array
          let idx = this.crowd.addAgent(start, this.getAgentParams(this.updateFlags)); //Assign that poPoly to the agent
          agent.idx = idx;

          //Now find the nearest valid location to the agent's desired destination
          //and assign that nearest point.
          let nearest = this.query.findNearestPoly(this.getEnd(agent), this.ext, this.filter);
          this.crowd.requestMoveTarget(agent.idx, nearest.getNearestRef(), nearest.getNearestPos());
          agent.hasEntered = true;
          agent.inSimulation = true;
        }
        if (agent.hasEntered) {
          let agentCurPos = [this.crowd.getAgent(j).npos[0], this.crowd.getAgent(j).npos[1], this.crowd.getAgent(j).npos[2]];
          let agentDes = this.getEnd(agent);

          let _x = agentCurPos.x - agent.destX;
          let _z = agentCurPos.z - agent.destZ;

          let distanceToDestination = Math.sqrt(_x * _x + _z * _z);
          if (distanceToDestination < 2) {
            this.crowd.removeAgent(agent.idx);
          }
          if (this.comparePos(agentCurPos, agentDes)) {
            agent.inSimulation = false;
          }
        }
      }
      this.crowd.update(1 / 25.0, null, i);

      //Update the current simulation time
      this.writeAgentPosition(currentMillisecond, this.objFilename, this.agentStartsFilename, this.ticks);
      currentMillisecond += millisecondsBetweenFrames;
    }
    return this.result;
  }
  
  truncate(num) {
    if (num > 0)
      return Math.floor(num)
    else
      return Math.ceil(num);
  }
  comparePos(pos1, pos2) {
    let res = false;

    let endX1 = pos1[0];
    let endY1 = pos1[2];

    let endX2 = pos2[0];
    let endY2 = pos2[2];

    let absEndX = Math.abs(this.truncate(endX1 - endX2));
    let absEndY = Math.abs(this.truncate(endY1 - endY2));

    if (absEndX <= 2 && absEndY <= 2) {
      res = true;
    }
    return res;
  }
  bootMesh(objFileContents) {
    this.nmd = RecastTestMeshBuilder.fromFile(objFileContents).getMeshData();
    this.navmesh = new NavMesh(this.nmd, 6, 0);
    this.query = new NavMeshQuery(this.navmesh);
    this.crowd = new Crowd(500, 0.6, this.navmesh);
    let params = new ObstacleAvoidanceParams();
    params.velBias = 0.5;
    params.adaptiveDivs = 5;
    params.adaptiveRings = 2;
    params.adaptiveDepth = 1;
    this.crowd.setObstacleAvoidanceParams(0, params);


    this.ap = this.getAgentParams(this.updateFlags);
    this.ext = this.crowd.getQueryExtents();
    this.filter = this.crowd.getFilter(0);
  }

  getAgentParams(updateFlags) {
    let ap = new CrowdAgentParams();
    ap.radius = 0.6;
    ap.height = 2;
    ap.maxAcceleration = 8.0;
    ap.maxSpeed = 2.5; //Originally 3.5f
    ap.collisionQueryRange = ap.radius * 12;
    ap.pathOptimizationRange = ap.radius * 30;
    ap.updateFlags = updateFlags;
    ap.obstacleAvoidanceType = 0;
    ap.separationWeight = 1; //Originally 2f
    return ap;
  }

  writeAgentPosition(currentMillisecond, objFilename, agentStartsFilename, ticks) {
   
    for (let j = 0; j < this.agents.length; j++) {
      let agent = this.agents[j];
      if (agent.idx != 0 && !agent.idx) continue;
      let idx = agent.idx;
      let x = this.crowd.getAgent(idx).npos[0];
      let y = this.crowd.getAgent(idx).npos[1];
      let z = this.crowd.getAgent(idx).npos[2];

      if (this.agents[j].inSimulation)
        this.result += ("" + j + "," + currentMillisecond + "," + x + "," + y + "," + z + "\n");
    }
  }
  newAgent(l,index) {
    let splits = l.split(",");
    if (splits.length == 0 || splits.length < 8)
      console.log("Error")

    let toReturn = {id:index};
    toReturn.startX = parseFloat(splits[2]);
    toReturn.startZ = parseFloat(splits[3]);
    toReturn.startY = parseFloat(splits[4]);

    toReturn.destX = parseFloat(splits[5]);
    toReturn.destZ = parseFloat(splits[6]);
    toReturn.destY = parseFloat(splits[7]);

    toReturn.startMSec = Math.floor(parseFloat(splits[1]));
    return toReturn;
  }
  getStart(agent) {return [agent.startX, agent.startZ, agent.startY];}
  getEnd(agent) { return [agent.destX, agent.destZ, agent.destY]; }
}

//module.exports = NodeApp;
export default NodeApp;