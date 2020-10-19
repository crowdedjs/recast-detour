import chai from "chai";
const assert = chai.assert;
// import NodeApp from "../examples/nodeNPM/NodeApp.js"
import crowded from "../examples/nodeNPM/crowded.js"
import fs from "fs";
import path from "path";


describe("Crowd Simulation", function () {
  let runs = [
    { obj: "hospital.obj", csv: "nine.csv", ticks: 500 },
    { obj: "hospital.obj", csv: "one.csv", ticks: 500 },

  ]
  for (let run of runs) {
    it(`Runs ${run.obj} ${run.csv} ${run.ticks}`, async function () {
      this.timeout(10000);

      let objFilename = run.obj;
      let agentStartsFilename = run.csv;
      let _ticks = run.ticks;

      // let app = new NodeApp(objFilename, agentStartsFilename, _ticks);
      // await app.go();
      crowded(objFilename, agentStartsFilename, _ticks);

      let baseFilename = `${objFilename}-${agentStartsFilename}-${_ticks}`;

      let outFilename = `${baseFilename}-out.csv`;
      let fullFilename = path.join(process.cwd(), "/" + outFilename);

      let out = fs.readFileSync(fullFilename, "utf-8");
      let outExpected = fs.readFileSync(path.join(process.cwd(), `/test/expected/${baseFilename}-expected.csv`), "utf-8");

      assert.equal(out, outExpected.replace(/\r\n/g, "\n"));
      fs.unlinkSync(fullFilename);
    })

  }

})