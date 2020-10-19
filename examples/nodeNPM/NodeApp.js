const Agent = require("./sim/Agent.js")
const fs = require("fs");
const path = require("path");
const util = require("util");
const stream = require("stream");
const { CrowdAgentParams, RecastTestMeshBuilder, NavMesh, NavMeshQuery, Crowd, ObstacleAvoidanceParams } = require("crowded")




class NodeApp {
  static updateFlags = CrowdAgentParams.DT_CROWD_ANTICIPATE_TURNS | CrowdAgentParams.DT_CROWD_OPTIMIZE_VIS
    | CrowdAgentParams.DT_CROWD_OPTIMIZE_TOPO | CrowdAgentParams.DT_CROWD_OBSTACLE_AVOIDANCE;
  static query;
  crowd;
  static agents = [];
  static ext;
  static filter;
  ap;
  /**
  * Fields required to run the open source backend. There is no need to touch these.
  */
  //md;
  navmesh;

  outStream;
  constructor(objString, agentStartsString, ticks) {
    NodeApp.agents = [];
    this.objFilename = objString;
    this.agentStartsFilename = agentStartsString;
    this.ticks = ticks;


  }
  async go() {
    let obj = fs.readFileSync(path.join(__dirname, "../objs/" + this.objFilename), "utf-8");
    //Boot simulation tells Recast to load the scene

    this.bootMesh(obj);

    //Path is the path to the file where we will store our results
    let result = fs.readFileSync(path.join(process.cwd(), "examples/agentStarts/" + this.agentStartsFilename), "utf-8");

    let stream = result.split('\n');
    Agent.index = 0;
    stream.forEach(l => l.trim().length > 0 ? NodeApp.agents.push(new Agent(l)) : 0 == 0);

    let currentMillisecond = 0; //The current time
    let millisecondsBetweenFrames = 40; //40ms between frames, or 25fps
    let secondsOfSimulation = this.ticks; //How long should the simulation run? Change this to whatever you want.
    for (let i = 0; i < secondsOfSimulation; i++) {
      if (i < 1) {
        // initialize all agent's id
        for (let id = 0; id < NodeApp.agents.length; id++) {
          let agent = NodeApp.agents[id];
          agent.setId(id);
        }
      }

      for (let j = 0; j < NodeApp.agents.length; j++) {
        let agent = NodeApp.agents[j]; //Grab each agent in the list

        //Ignore agents that have come into the simulation and exited
        if (agent.hasEntered && !agent.inSimulation) continue;

        //See if we need to add the agent to the simulation
        if (!agent.hasEntered && agent.startMSec <= currentMillisecond) {
          let start = agent.getStart();//Get the agent's starting point as a  array
          let idx = this.crowd.addAgent(start, this.getAgentParams(NodeApp.updateFlags)); //Assign that poPoly to the agent
          agent.idx = idx;

          //Now find the nearest valid location to the agent's desired destination
          //and assign that nearest point.
          let nearest = this.query.findNearestPoly(agent.getEnd(), this.ext, this.filter);
          this.crowd.requestMoveTarget(agent.idx, nearest.getNearestRef(), nearest.getNearestPos());
          agent.hasEntered = true;
          agent.inSimulation = true;
        }
        if (agent.hasEntered) {

          agent.setActive(true);
          agent.setActive(false);
          let agentCurPos = [this.crowd.getAgent(j).npos[0], this.crowd.getAgent(j).npos[1], this.crowd.getAgent(j).npos[2]];
          let agentDes = agent.getEnd();

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
    await this.finish();
  }
  async finish() {
    //= require( https://nodesource.com/blog/understanding-streams-in-nodejs/
    const finished = util.promisify(stream.finished); // (A)

    this.outStream.end();
    await finished(this.outStream);
    console.log("Ended outStream")
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
    if (!this.outStream) {
      this.outStream = fs.createWriteStream(`${objFilename}-${agentStartsFilename}-${ticks}-out.csv`, { emitClose: true })
    }
    for (let j = 0; j < NodeApp.agents.length; j++) {
      let agent = NodeApp.agents[j];
      if (agent.idx != 0 && !agent.idx) continue;
      let idx = agent.idx;
      let x = this.crowd.getAgent(idx).npos[0];
      let y = this.crowd.getAgent(idx).npos[1];
      let z = this.crowd.getAgent(idx).npos[2];

      if (NodeApp.agents[j].inSimulation)
        this.outStream.write("" + j + "," + currentMillisecond + "," + x + "," + y + "," + z + "\n");
    }
  }
}

module.exports = NodeApp;