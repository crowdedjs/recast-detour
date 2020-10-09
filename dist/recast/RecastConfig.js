"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

/*
Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
Recast4J Copyright (c) 2015 Piotr Piastucki piotr@jtilia.org

This software is provided 'as-is', without any express or implied
warranty.  In no event will the authors be held liable for any damages
arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:
1. The origin of this software must not be misrepresented; you must not
 claim that you wrote the original software. If you use this software
 in a product, an acknowledgment in the product documentation would be
 appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
 misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
var RecastConfig =
/** The width/height size of tile's on the xz-plane. [Limit: >= 0] [Units: vx] **/

/** The xz-plane cell size to use for fields. [Limit: > 0] [Units: wu] **/

/** The y-axis cell size to use for fields. [Limit: > 0] [Units: wu] **/

/** The maximum slope that is considered walkable. [Limits: 0 <= value < 90] [Units: Degrees] **/

/**
 * Minimum floor to 'ceiling' height that will still allow the floor area to be considered walkable. [Limit: >= 3]
 * [Units: vx]
 **/

/** Maximum ledge height that is considered to still be traversable. [Limit: >=0] [Units: vx] **/

/**
 * The distance to erode/shrink the walkable area of the heightfield away from obstructions. [Limit: >=0] [Units:
 * vx]
 **/

/** The maximum allowed length for contour edges aPoly the border of the mesh. [Limit: >=0] [Units: vx] **/

/**
 * The maximum distance a simplfied contour's border edges should deviate the original raw contour. [Limit: >=0]
 * [Units: vx]
 **/

/** The minimum number of cells allowed to form isolated island areas. [Limit: >=0] [Units: vx] **/

/**
 * Any regions with a span count smaller than this value will, if possible, be merged with larger regions. [Limit:
 * >=0] [Units: vx]
 **/

/**
 * The maximum number of vertices allowed for polygons generated during the contour to polygon conversion process.
 * [Limit: >= 3]
 **/

/**
 * Sets the sampling distance to use when generating the detail mesh. (For height detail only.) [Limits: 0 or >=
 * 0.9] [Units: wu]
 **/

/**
 * The maximum distance the detail mesh surface should deviate from heightfield data. (For height detail only.)
 * [Limit: >=0] [Units: wu]
 **/
function RecastConfig(partitionType, cellSize, cellHeight, agentHeight, agentRadius, agentMaxClimb, agentMaxSlope, regionMinSize, regionMergeSize, edgeMaxLen, edgeMaxError, vertsPerPoly, detailSampleDist, detailSampleMaxError, tileSize, walkableAreaMod) {
  _classCallCheck(this, RecastConfig);

  _defineProperty(this, "partitionType", null);

  _defineProperty(this, "tileSize", void 0);

  _defineProperty(this, "cs", void 0);

  _defineProperty(this, "ch", void 0);

  _defineProperty(this, "walkableSlopeAngle", void 0);

  _defineProperty(this, "walkableHeight", void 0);

  _defineProperty(this, "walkableClimb", void 0);

  _defineProperty(this, "walkableRadius", void 0);

  _defineProperty(this, "maxEdgeLen", void 0);

  _defineProperty(this, "maxSimplificationError", void 0);

  _defineProperty(this, "minRegionArea", void 0);

  _defineProperty(this, "mergeRegionArea", void 0);

  _defineProperty(this, "maxVertsPerPoly", 0);

  _defineProperty(this, "detailSampleDist", void 0);

  _defineProperty(this, "detailSampleMaxError", void 0);

  _defineProperty(this, "walkableAreaMod", void 0);

  this.partitionType = partitionType;
  this.cs = cellSize;
  this.ch = cellHeight;
  this.walkableSlopeAngle = agentMaxSlope;
  this.walkableHeight = Math.ceil(agentHeight / this.ch);
  this.walkableClimb = Math.floor(agentMaxClimb / this.ch);
  this.walkableRadius = Math.ceil(agentRadius / this.cs);
  this.maxEdgeLen = Math.floor(edgeMaxLen / cellSize);
  this.maxSimplificationError = edgeMaxError;
  this.minRegionArea = regionMinSize * regionMinSize; // Note: area = size*size

  this.mergeRegionArea = regionMergeSize * regionMergeSize; // Note: area = size*size

  this.maxVertsPerPoly = vertsPerPoly;
  this.detailSampleDist = detailSampleDist < 0.9 ? 0 : cellSize * detailSampleDist;
  this.detailSampleMaxError = cellHeight * detailSampleMaxError;
  this.tileSize = tileSize;
  this.walkableAreaMod = walkableAreaMod;
};

var _default = RecastConfig;
exports["default"] = _default;