
import None from "../behavior/None.js"

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
  behavior;

  ////////////zhicheng////////////
  headToDes;
  paired;
  friendList = [];
  id;
  checked;
  waiting;
  holdTrash;
  active;
  lead;
  follow;
  prevPos = [];
  followTo;
  followers = {};
  gatePick;
  gateIndex;


  ////////////zhicheng////////////
  static  index = 0;
   idx; //Corresponds to the internal idx number used by recast

  constructor(l) {
    
    let splits = l.split(",");
    if(splits.length == 0 || splits.length < 8)
      console.log("Error")

    this.startX = parseFloat(splits[2]);
    this.startZ = parseFloat(splits[3]);
    this.startY = parseFloat(splits[4]);

    this.destX = parseFloat(splits[5]);
    this.destZ = parseFloat(splits[6]);
    this.destY = parseFloat(splits[7]);

    this.startMSec = Math.floor(parseFloat(splits[1]));

    ///Now check to see what behavior this agent should receive
    if (splits.length == 8 || splits[8].trim() == "" || splits[8].trim().toLowerCase() == "none") //i.e. there is no behavior specification
      this.behavior = new None(Agent.index++);
    else {
      //In this case we have to figure it out
      //The behavior is determined by the first word after the comma. Everything else can be a argument in the constructor
      let behave = splits[8].trim().toLowerCase();
      if (behave == "queue")
        this.behavior = new Queue(CrowdSimApp.query, CrowdSimApp.crowd, CrowdSimApp.agents, CrowdSimApp.ext, CrowdSimApp.filter);
      else if (behave =="flee") {
        this.behavior = new PerfectFlee(Agent.index++, new Vector3(0, 0, 0));
      }
    }

    ////////////zhicheng////////////
    this.checked = false;
    this.waiting = false;
    this.paired = false;
    this.friendList = [];
    this.holdTrash = false;
    this.active = false;
    this.follow = true;
    this.lead = false;
    this.followTo = -1;
    this.prevPos = [];
    this.followers = {};
    this.gatePick = false;
    ////////////zhicheng////////////
  }
  getFollowers() { return this.followers; }
	
	 addFollower(a) { this.followers.push(a.getId()); }

getFollowToId() { return this.followTo; }

setFollowTo(leadId) { this.followTo = leadId; }

getPrevPos() { return this.prevPos; }

setPrevPos( prevPos) {
  this.prevPos[0] = prevPos[0];
  this.prevPos[1] = prevPos[1];
  this.prevPos[2] = prevPos[2];
}

isFollow() { return follow; }

setFollow(follow) { this.follow = follow; }

isLead() { return this.lead; }

 setLead(lead) { this.lead = lead; }

getStart() { return [ this.startX, this.startZ, this.startY ]; }

getEnd() { return [ this.destX, this.destZ, this.destY ]; }

////////////////////////zhicheng////////////////////////
 setChecked() { this.checked = true; }

getChecked() { return this.checked; }

 setId(i) { this.id = i; }

getId() { return this.id; }

 setWaiting() { this.waiting = true; }

 setWaitingFalse() { this.waiting = false; }

isWaiting() { return this.waiting; }

//      addFriend(Agent agentFriend) { this.friendList.push(agentFriend); }
//    
//     int getNumOfFriend() { return this.friendList.length; }
//    
//     ArrayList<Agent> getAgentFriendList() { return this.friendList; }

getStartMSec() { return this.startMSec; }

 setPairState(state) { this.paired = state; }

getPairState() { return this.paired; }

 setHasTrash() { this.holdTrash = true; }

 setHasNoTrash() { this.holdTrash = false; }
/**
 * return treu if agent holding trash
 * return false if agent not holding trash.
 * @return boolean
 */
trashCondition() { return this.holdTrash; }

/**
 * If current agent is active, update its newDestination to App() class's update method
 */
isActive() { return active; }
 setActive(active) { this.active = active; }

 getHeadToDes() { return headToDes; }

 setHeadToDes(headToDes) { this.headToDes = headToDes; }


    ////////////////////////zhicheng////////////////////////
    
}

export default Agent;