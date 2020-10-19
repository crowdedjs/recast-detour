

class Agent {
  startX;
  startY;
  startZ;
  destX;
  destY;
  destZ;
  startMSec;
  inSimulation = false;
  hasEntered = false;
  id;
  active;

  static index = 0;
  idx; //Corresponds to the internal idx number used by recast

  constructor(l) {
    let splits = l.split(",");
    if (splits.length == 0 || splits.length < 8)
      console.log("Error")

    this.startX = parseFloat(splits[2]);
    this.startZ = parseFloat(splits[3]);
    this.startY = parseFloat(splits[4]);

    this.destX = parseFloat(splits[5]);
    this.destZ = parseFloat(splits[6]);
    this.destY = parseFloat(splits[7]);

    this.startMSec = Math.floor(parseFloat(splits[1]));
  }

  getStart() { return [this.startX, this.startZ, this.startY]; }
  getEnd() { return [this.destX, this.destZ, this.destY]; }

  setId(i) { this.id = i; }
  getId() { return this.id; }

  isActive() { return active; }
  setActive(active) { this.active = active; }
}

module.exports = Agent;