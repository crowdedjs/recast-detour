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

/** Provides high level information related to a dt object.*/
var MeshHeader = function MeshHeader() {
  _classCallCheck(this, MeshHeader);

  _defineProperty(this, "magic", void 0);

  _defineProperty(this, "version", void 0);

  _defineProperty(this, "x", void 0);

  _defineProperty(this, "y", void 0);

  _defineProperty(this, "layer", void 0);

  _defineProperty(this, "userId", void 0);

  _defineProperty(this, "polyCount", void 0);

  _defineProperty(this, "vertCount", void 0);

  _defineProperty(this, "maxLinkCount", void 0);

  _defineProperty(this, "detailMeshCount", void 0);

  _defineProperty(this, "detailVertCount", void 0);

  _defineProperty(this, "detailTriCount", void 0);

  _defineProperty(this, "bvNodeCount", void 0);

  _defineProperty(this, "offMeshConCount", void 0);

  _defineProperty(this, "offMeshBase", void 0);

  _defineProperty(this, "walkableHeight", void 0);

  _defineProperty(this, "walkableRadius", void 0);

  _defineProperty(this, "walkableClimb", void 0);

  _defineProperty(this, "bmin", new Array(3));

  _defineProperty(this, "bmax", new Array(3));

  _defineProperty(this, "bvQuantFactor", void 0);
};

_defineProperty(MeshHeader, "DT_NAVMESH_MAGIC", 'D' << 24 | 'N' << 16 | 'A' << 8 | 'V');

_defineProperty(MeshHeader, "DT_NAVMESH_VERSION", 7);

_defineProperty(MeshHeader, "DT_NAVMESH_VERSION_RECAST4J", 0x8807);

_defineProperty(MeshHeader, "DT_NAVMESH_STATE_MAGIC", 'D' << 24 | 'N' << 16 | 'M' << 8 | 'S');

_defineProperty(MeshHeader, "DT_NAVMESH_STATE_VERSION", 1);

var _default = MeshHeader;
exports["default"] = _default;