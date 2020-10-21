const NodeApp = require("./NodeApp")
const fs = require("fs");
const path = require("path");


function crowded(objFilename, agentStartsFilename, ticks){
  let obj = fs.readFileSync(path.join(__dirname, "../../objs/" + objFilename), "utf-8");
  let agents = fs.readFileSync(path.join(__dirname, "../../agentStarts/" + agentStartsFilename), "utf-8");

  let app = new NodeApp(obj, agents, ticks);
  let result = app.go();
  fs.writeFileSync(`${objFilename}-${agentStartsFilename}-${ticks}-out.csv`, result);
}

module.exports = crowded;