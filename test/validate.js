import chai from "chai";
const assert = chai.assert;
import NodeApp from "../examples/node/NodeApp.js"
import fs from "fs";
import path from "path";


describe("Crowd Simulation", function(){
  it("Runs", function(){
    this.timeout(10000);
    new NodeApp();
    let out = fs.readFileSync(path.join(process.cwd(), "/examples/node/out.csv"), "utf-8");
    let outExpected = fs.readFileSync(path.join(process.cwd(), "/test/outExpected.csv"), "utf-8");

    assert.equal(out,outExpected)
  })
})