import CrowdSimApp from "./sim/CrowdSimApp.js"
import Agent from "./sim/Agent.js"
import Vector3 from "./behavior/Vector3.js"
import fs from "fs";


class NodeApp extends CrowdSimApp {
  constructor() {
    super();
    let obj = fs.readFileSync("./exampleNode/hospital.obj", "utf-8");
    //Boot simulation tells Recast to load the scene
       
    this.bootMesh(obj);

    //Path is the path to the file where we will store our results
    //Currently, this is out.csv.
    let result = fs.readFileSync("./exampleNode/in.csv", "utf-8");

    let stream = result.split('\n');
    stream.forEach(l => l.trim().length > 0 ? CrowdSimApp.agents.push(new Agent(l)) : 0 == 0);
            


    let currentMillisecond = 0; //The current time
    let millisecondsBetweenFrames = 40; //40ms between frames, or 25fps
    let secondsOfSimulation = 500; //How long should the simulation run? Change this to whatever you want.
    for (let i = 0; i < secondsOfSimulation; i++) {
      if (i == 29)
        console.log("bog")
      // if(i%10 == 0)
      //   console.log("Pause")
      // console.log(i);
      // if (i == 491)
      //   console.log("Bug")
      if (i < 1) {
        // initialize all agent's id
        for (let id = 0; id < CrowdSimApp.agents.length; id++) {
          let agent = CrowdSimApp.agents[id];
          agent.setId(id);
          // agent.setChecked(false);
        }
      }

      // for(let  a : agents) System.out.println("who is checked? "+a.getChecked() + " who is waiting? " + a.isWaiting()  + a.getId() + " " + i);
      for (let j = 0; j < CrowdSimApp.agents.length; j++) {
        if (j == 12 && currentMillisecond == 42520)
          console.log("Bug")

        let agent = CrowdSimApp.agents[j]; //Grab each agent in the list

        //Ignore agents that have come into the simulation and exited
        if (agent.hasEntered && !agent.inSimulation) continue;

        //See if we need to add the agent to the simulation
        if (!agent.hasEntered && agent.startMSec <= currentMillisecond) {
          let start = agent.getStart();//Get the agent's starting point as a  array
          let idx = this.crowd.addAgent(start, this.getAgentParams(CrowdSimApp.updateFlags)); //Assign that poPoly to the agent
          agent.idx = idx;

          //Now find the nearest valid location to the agent's desired destination
          //and assign that nearest point.
          let nearest = this.query.findNearestPoly(agent.getEnd(), this.ext, this.filter);
          this.crowd.requestMoveTarget(agent.idx, nearest.getNearestRef(), nearest.getNearestPos());
          agent.hasEntered = true;
          agent.inSimulation = true;
        }

        if (agent.hasEntered) {
          let newDestination = new Vector3(0, 0.318020731, 12);

          let limit = currentMillisecond % 450;
          if (/*limit > 400*/true) {
            agent.setActive(true);
            let newDestination = agent.behavior.update(CrowdSimApp.agents, this.crowd, currentMillisecond);
            agent.setActive(false);
            let agentCurPos = [this.crowd.getAgent(j).npos[0], this.crowd.getAgent(j).npos[1], this.crowd.getAgent(j).npos[2]];
            let agentDes = agent.getEnd();

            // newDestination = new  (f[0], f[1], f[2]);
            if (newDestination != null) {
              let checkTempName = newDestination.asArray();
              let nearest = query.findNearestPoly(newDestination.asArray(), ext, filter);
              this.crowd.requestMoveTarget(agent.idx, nearest.getNearestRef(), nearest.getNearestPos());
            } else if (agent.inSimulation == false) {
              this.crowd.removeAgent(agent.idx);
              console.log("removing " + j);
            } else {

            }
            if (this.comparePos(agentCurPos, agentDes)) {
              agent.inSimulation = false;
              //                            System.out.println("==========removing " + j);
            }
          }
        }
        if (j == 0) {
          if(CrowdSimApp.agents[j].idx != 0 && !CrowdSimApp.agents[j].idx) continue
          if(!CrowdSimApp.agents[j].inSimulation) continue;
          let idx = CrowdSimApp.agents[j].idx;
          let x = this.crowd.getAgent(idx).npos[0];
          if(Number.isNaN(x)) continue;
          let y = this.crowd.getAgent(idx).npos[1];
          let z = this.crowd.getAgent(idx).npos[2];
          // console.log(`${i} ${this.crowd.m_agents[0].boundary.m_center} ${this.crowd.m_agents[0].boundary.m_segs.length} ${x} ${y} ${z}`);
          console.log(`${i} ${this.crowd.m_agents[0].boundary.m_center[0]} ${this.crowd.m_agents[0].boundary.m_center[2]} ${x} ${z}`);
        }
      }
      this.crowd.update(1 / 25.0, null, i);

      //Write the agents' positions in a file
      //this.writeAgentPosition(currentMillisecond);
      // console.log("Good")
      //Update the current simulation time
      currentMillisecond += millisecondsBetweenFrames;


    }


    //this.writeJSFile();

  }
  truncate(num){
    if(num > 0)
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
}

export default NodeApp;