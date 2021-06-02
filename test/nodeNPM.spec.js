import chai from "chai";
const assert = chai.assert;
const expect = chai.expect;
// import NodeApp from "../examples/nodeNPM/NodeApp.js"
//import crowded from "../examples/node/nodeNPM/crowded.js"
import crowded from "./crowded.js"
import fs from "fs";
import path from "path";
import { performance } from "perf_hooks";


describe("Crowd Simulation", function () {
  let runs = [
    { obj: "hospital.obj", csv: "nine.csv", ticks: 1 },
    { obj: "hospital.obj", csv: "one.csv", ticks: 1 },

  ]
  for (let run of runs) {

    it(`Runs ${run.obj} ${run.csv} ${run.ticks}`, async function () {
      this.timeout(10000);

      let objFilename = run.obj;
      let agentStartsFilename = run.csv;
      let _ticks = run.ticks;

      let start = performance.now();
      let result = crowded(objFilename, agentStartsFilename, _ticks);
      let end = performance.now();
      console.log(end - start);

      let baseFilename = `${objFilename}-${agentStartsFilename}-${_ticks}`;

      let outFilename = `${baseFilename}-out.csv`;
      let fullFilename = path.join(process.cwd(), "/" + outFilename);

      //let out = fs.readFileSync(fullFilename, "utf-8");
      let out = result;
      let outExpected = fs.readFileSync(path.join(process.cwd(), `/test/expected/${baseFilename}-expected.csv`), "utf-8");

      //Make sure the out file matches the expected one
      //assert.equal(out, outExpected.replace(/\r\n/g, "\n"));
      //expect('1').to.equal('2')
      try {
        expect(out).to.equal(outExpected.replace(/\r\n/g, "\n"))
      } catch (e) {
        assert.fail();
        console.error("Didn't match")
        //console.error(e);
      }
      //Delete the out file.
      //fs.unlinkSync(fullFilename);
    })

  }

})