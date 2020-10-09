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
/// Represents the source data used to build an navigation mesh tile.
var NavMeshDataCreateParams = function NavMeshDataCreateParams() {
  _classCallCheck(this, NavMeshDataCreateParams);

  _defineProperty(this, "verts", []);

  _defineProperty(this, "vertCount", 0);

  _defineProperty(this, "polys", []);

  _defineProperty(this, "polyFlags", []);

  _defineProperty(this, "polyAreas", []);

  _defineProperty(this, "polyCount", 0);

  _defineProperty(this, "nvp", 0);

  _defineProperty(this, "detailMeshes", []);

  _defineProperty(this, "detailVerts", []);

  _defineProperty(this, "detailVertsCount", 0);

  _defineProperty(this, "detailTris", []);

  _defineProperty(this, "detailTriCount", 0);

  _defineProperty(this, "offMeshConVerts", []);

  _defineProperty(this, "offMeshConRad", []);

  _defineProperty(this, "offMeshConFlags", []);

  _defineProperty(this, "offMeshConAreas", []);

  _defineProperty(this, "offMeshConDir", []);

  _defineProperty(this, "offMeshConUserID", []);

  _defineProperty(this, "offMeshConCount", 0);

  _defineProperty(this, "userId", 0);

  _defineProperty(this, "tileX", 0);

  _defineProperty(this, "tileY", 0);

  _defineProperty(this, "tileLayer", 0);

  _defineProperty(this, "bmin", []);

  _defineProperty(this, "bmax", []);

  _defineProperty(this, "walkableHeight", 0);

  _defineProperty(this, "walkableRadius", 0);

  _defineProperty(this, "walkableClimb", 0);

  _defineProperty(this, "cs", 0);

  _defineProperty(this, "ch", 0);

  _defineProperty(this, "buildBvTree", false);
} /// @}
;

var _default = NavMeshDataCreateParams;
exports["default"] = _default;