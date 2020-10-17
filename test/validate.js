import chai from "chai";
const assert = chai.assert;
import NodeApp from "../examples/node/NodeApp.js"
import fs from "fs";
import path from "path";


describe("Crowd Simulation", function(){
  let objs = ["hospital.obj"];
  let csvs = ["nine.csv", "one.csv"];
  let ticks = [500];
  for(let obj of objs){
    for(let csv of csvs){
      for(let tick of ticks){
        it(`Runs ${obj} ${csv} ${tick}`, async function(){
          this.timeout(10000);
      
          let objFilename = obj;
          let agentStartsFilename = csv;
          let _ticks = tick;
      
          let app = new NodeApp(objFilename, agentStartsFilename, _ticks);
          await app.go();
      
          let baseFilename = `${objFilename}-${agentStartsFilename}-${_ticks}`;
      
          let outFilename = `${baseFilename}-out.csv`;
          let fullFilename = path.join(process.cwd(),  "/" + outFilename);
      
          let out = fs.readFileSync(fullFilename, "utf-8");
          let outExpected = fs.readFileSync(path.join(process.cwd(), `/test/expected/${baseFilename}-expected.csv`), "utf-8");
      
          assert.equal(out,outExpected.replace(/\r\n/g, "\n"));
          fs.unlinkSync(fullFilename);
        })
      }
    }
  }
  
})