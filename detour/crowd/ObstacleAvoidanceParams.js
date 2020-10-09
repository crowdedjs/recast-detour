class ObstacleAvoidanceParams {
  velBias=0;
  weightDesVel=0;
  weightCurVel=0;
  weightSide=0;
  weightToi=0;
  horizTime=0;
  gridSize=0; ///< grid
  adaptiveDivs=0; ///< adaptive
  adaptiveRings=0; ///< adaptive
  adaptiveDepth=0; ///< adaptive

  constructor() {
    this.velBias = 0.4;
    this.weightDesVel = 2.0;
    this.weightCurVel = 0.75;
    this.weightSide = 0.75;
    this.weightToi = 2.5;
    this.horizTime = 2.5;
    this.gridSize = 33;
    this.adaptiveDivs = 7;
    this.adaptiveRings = 2;
    this.adaptiveDepth = 5;
  }
};

export default ObstacleAvoidanceParams;