const ABehavior = require( "./ABehavior.js")


/**
 * Created by bricks on 6/19/2018.
 */
class None extends ABehavior{


    constructor(myIndex) {
        super(myIndex);
    }

    checkEndOfSimulation() {
        let _x = this.location.x - this.agents[this.index].destX;
        let _z = this.location.z - this.agents[this.index].destZ;
        let _y = this.location.y - this.agents[this.index].destY;

        let distanceToDestination = Math.sqrt(_x * _x + _z * _z);

        if (distanceToDestination < 2) {
            this.agents[this.index].inSimulation = false;
            return null;
        }
        else {
            return this.mainBehavior();
        }
    }

    mainBehavior() {
        //Nothing, because we just head to the destination.
        //Meant to be overridden by child classes.

        return null;
    }
}

module.exports = None;
