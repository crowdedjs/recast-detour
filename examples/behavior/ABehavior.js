


/**
 * Created by bricks on 6/19/2018.
 */

 import Vector3 from "./Vector3.js"

class ABehavior {

    index;
    location;
    agents;
    crowd;
    msec;


    constructor(myIndex) {
        this.index = myIndex;
    }

    update(agents, crowd, msec) {
        //Do nothing since this is the default "none" behavior.

        let idx = agents[this.index].idx;
        let x = crowd.getAgent(idx).npos[0];
        let y = crowd.getAgent(idx).npos[1];
        let z = crowd.getAgent(idx).npos[2];

        this.location = new Vector3(x, y, z);
        this.agents = agents;
        this.crowd = crowd;
        this.msec = msec;

        return this.checkEndOfSimulation();

    }


}
export default ABehavior;