import NodeApp from "./NodeApp.js"
import fs from "fs";
import path from "path";

let backtrack = ["../../../examples/objs/", "../../../examples/agentStarts/"]
crowded("hospital.obj", "nine.csv", 500)

function crowded(objFilename, agentStartsFilename, ticks){
  let obj = fs.readFileSync(path.join(process.cwd(), "/examples/objs/" + backtrack[0] + objFilename), "utf-8");
  let agents = fs.readFileSync(path.join(process.cwd(), "/examples/agentStarts/" + backtrack[1] + agentStartsFilename), "utf-8");

  let app = new NodeApp(obj, agents, ticks);
  let result = app.go();
  fs.writeFileSync(`${objFilename}-${agentStartsFilename}-${ticks}-out.csv`, result);
  return result;
}

//module.exports = crowded;
export default crowded;