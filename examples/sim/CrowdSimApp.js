import CrowdAgentParams from "../../src/detour/crowd/CrowdAgentParams.js"
import RecastTestMeshBuilder from "../../src/detour/RecastTestMeshBuilder.js"
import NavMesh from "../../src/detour/NavMesh.js"
import NavMeshQuery from "../../src/detour/NavMeshQuery.js"
import Crowd from "../../src/detour/crowd/Crowd.js"
import ObstacleAvoidanceParams from "../../src/detour/crowd/ObstacleAvoidanceParams.js"

class CrowdSimApp {

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
    md;
    navmesh;

    outStream;



    bootMesh(objFileContents) {


        this.nmd = RecastTestMeshBuilder.fromFile(objFileContents).getMeshData();
        //.log(JSON.stringify(this.nmd, null, 2));
        this.navmesh = new NavMesh(this.nmd, 6, 0);
        this.query = new NavMeshQuery(this.navmesh);
        this.crowd = new Crowd(500, 0.6, this.navmesh);
        let params = new ObstacleAvoidanceParams();
        params.velBias = 0.5;
        params.adaptiveDivs = 5;
        params.adaptiveRings = 2;
        params.adaptiveDepth = 1;
        this.crowd.setObstacleAvoidanceParams(0, params);
        params = new ObstacleAvoidanceParams();
        params.velBias = 0.5;
        params.adaptiveDivs = 5;
        params.adaptiveRings = 2;
        params.adaptiveDepth = 2;
        this.crowd.setObstacleAvoidanceParams(1, params);
        params = new ObstacleAvoidanceParams();
        params.velBias = 0.5;
        params.adaptiveDivs = 7;
        params.adaptiveRings = 2;
        params.adaptiveDepth = 3;
        this.crowd.setObstacleAvoidanceParams(2, params);
        params = new ObstacleAvoidanceParams();
        params.velBias = 0.5;
        params.adaptiveDivs = 7;
        params.adaptiveRings = 3;
        params.adaptiveDepth = 3;
        this.crowd.setObstacleAvoidanceParams(3, params);

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

    writeAgentPosition(currentMillisecond) {
        for(let j = 0; j < CrowdSimApp.agents.length; j++) {
            let agent = CrowdSimApp.agents[j];
            if(agent.idx != 0 && !agent.idx) continue;
            let idx = agent.idx;
            let x = this.crowd.getAgent(idx).npos[0];
            let y = this.crowd.getAgent(idx).npos[1];
            let z = this.crowd.getAgent(idx).npos[2];

            if (CrowdSimApp.agents[j].inSimulation)
                this.outStream.write("" + j + "," + currentMillisecond + "," + x + "," + y + "," + z + "\n");
        }
    }

}



export default CrowdSimApp;