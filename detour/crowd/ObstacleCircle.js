class ObstacleCircle {
  /** Position of the obstacle */
  p = new Array(3);
  /** Velocity of the obstacle */
  vel = new Array(3);
  /** Velocity of the obstacle */
  dvel = new Array(3);
  /** Radius of the obstacle */
  rad;
  /** Use for side selection during sampling. */
  dp = new Array(3);
  /** Use for side selection during sampling. */
  np = new Array(3);
}

export default ObstacleCircle;