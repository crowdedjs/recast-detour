"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _NavMesh = _interopRequireDefault(require("./NavMesh.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

/**
 * Defines a navigation mesh tile.
 */
var MeshTile =
/** Counter describing modifications to the tile. */

/** The tile data. */

/** The tile links. */

/** Index to the next free link. */

/** Tile flags. (See: #dtTileFlags) */

/** The next free tile, or the next tile in the spatial grid. */
function MeshTile(index) {
  _classCallCheck(this, MeshTile);

  _defineProperty(this, "index", 0);

  _defineProperty(this, "salt", 0);

  _defineProperty(this, "data", null);

  _defineProperty(this, "links", []);

  _defineProperty(this, "linksFreeList", _NavMesh["default"].DT_NULL_LINK);

  _defineProperty(this, "flags", 0);

  _defineProperty(this, "next", null);

  this.index = index;
};

var _default = MeshTile;
exports["default"] = _default;