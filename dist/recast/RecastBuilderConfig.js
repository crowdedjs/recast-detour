"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _RecastVectors = _interopRequireDefault(require("./RecastVectors.js"));

var _Recast = _interopRequireDefault(require("./Recast.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var RecastBuilderConfig =
/** The width of the field aPoly the x-axis. [Limit: >= 0] [Units: vx] **/

/** The height of the field aPoly the z-axis. [Limit: >= 0] [Units: vx] **/

/** The minimum bounds of the field's AABB. [(x, y, z)] [Units: wu] **/

/** The maximum bounds of the field's AABB. [(x, y, z)] [Units: wu] **/

/** The size of the non-navigable border around the heightfield. [Limit: >=0] [Units: vx] **/

/** Set to true for tiled build **/

/** Set to false to disable building detailed mesh **/
// RecastBuilderConfig(RecastConfig cfg, bmin, bmax) {
//         this(cfg, bmin, bmax, 0, 0, false);
//     }
// RecastBuilderConfig(RecastConfig cfg, bmin, bmax,  tx,  ty,  tiled) {
//         this(cfg, bmin, bmax, tx, ty, tiled, true);
//     }
function RecastBuilderConfig(cfg, bmin, bmax) {
  var tx = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : 0;
  var ty = arguments.length > 4 && arguments[4] !== undefined ? arguments[4] : 0;
  var tiled = arguments.length > 5 && arguments[5] !== undefined ? arguments[5] : false;
  var buildMeshDetail = arguments.length > 6 && arguments[6] !== undefined ? arguments[6] : true;

  _classCallCheck(this, RecastBuilderConfig);

  _defineProperty(this, "cfg", null);

  _defineProperty(this, "width", 0);

  _defineProperty(this, "height", 0);

  _defineProperty(this, "bmin", new Array(3));

  _defineProperty(this, "bmax", new Array(3));

  _defineProperty(this, "borderSize", 0);

  _defineProperty(this, "tiled", 0);

  _defineProperty(this, "buildMeshDetail", false);

  this.cfg = cfg;
  this.tiled = tiled;
  this.buildMeshDetail = buildMeshDetail;

  _RecastVectors["default"].copy2(this.bmin, bmin);

  _RecastVectors["default"].copy2(this.bmax, bmax);

  if (tiled) {
    ts = cfg.tileSize * cfg.cs;
    this.bmin[0] += tx * ts;
    this.bmin[2] += ty * ts;
    this.bmax[0] = this.bmin[0] + ts;
    this.bmax[2] = this.bmin[2] + ts; // Expand the heighfield bounding box by border size to find the extents of geometry we need to build this
    // tile.
    //
    // This is done in order to make sure that the navmesh tiles connect correctly at the borders,
    // and the obstacles close to the border work correctly with the dilation process.
    // No polygons (or contours) will be created on the border area.
    //
    // IMPORTANT!
    //
    // :''''''''':
    // : +-----+ :
    // : | | :
    // : | |<--- tile to build
    // : | | :
    // : +-----+ :<-- geometry needed
    // :.........:
    //
    // You should use this bounding box to query your input geometry.
    //
    // For example if you build a navmesh for terrain, and want the navmesh tiles to match the terrain tile size
    // you will need to pass in data from neighbour terrain tiles too! In a simple case, just pass in all the 8
    // neighbours,
    // or use the bounding box below to only pass in a sliver of each of the 8 neighbours.

    this.borderSize = cfg.walkableRadius + 3; // Reserve enough padding.

    this.bmin[0] -= this.borderSize * cfg.cs;
    this.bmin[2] -= this.borderSize * cfg.cs;
    this.bmax[0] += this.borderSize * cfg.cs;
    this.bmax[2] += this.borderSize * cfg.cs;
    this.width = cfg.tileSize + this.borderSize * 2;
    this.height = cfg.tileSize + this.borderSize * 2;
  } else {
    var wh = _Recast["default"].calcGridSize(this.bmin, this.bmax, cfg.cs);

    this.width = wh[0];
    this.height = wh[1];
    this.borderSize = 0;
  } //ADDED:


  this.width = Math.floor(this.width);
  this.height = Math.floor(this.height);
};

var _default = RecastBuilderConfig;
exports["default"] = _default;