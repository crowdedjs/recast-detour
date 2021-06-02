import NodeApp from "./NodeApp.js"
import fs from "fs";
import path from "path";


function crowded(objFilename, agentStartsFilename, ticks){
  let obj = fs.readFileSync(path.join(process.cwd(), "/examples/objs/" + objFilename), "utf-8");
  let agents = fs.readFileSync(path.join(process.cwd(), "/examples/agentStarts/" + agentStartsFilename), "utf-8");

  let app = new NodeApp(obj, agents, ticks);
  let result = app.go();
  //fs.writeFileSync(`${objFilename}-${agentStartsFilename}-${ticks}-out.csv`, result);
  return result;
}

//module.exports = crowded;
export default crowded;