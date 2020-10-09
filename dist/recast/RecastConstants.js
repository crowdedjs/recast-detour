"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

/*
Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
Recast4J Copyright (c) 2015-2018 Piotr Piastucki piotr@jtilia.org

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
var RecastConstants = function RecastConstants() {
  _classCallCheck(this, RecastConstants);
};

_defineProperty(RecastConstants, "RC_NULL_AREA", 0);

_defineProperty(RecastConstants, "RC_NOT_CONNECTED", 0x3f);

_defineProperty(RecastConstants, "RC_WALKABLE_AREA", 63);

_defineProperty(RecastConstants, "RC_SPAN_HEIGHT_BITS", 13);

_defineProperty(RecastConstants, "RC_SPAN_MAX_HEIGHT", (1 << RecastConstants.RC_SPAN_HEIGHT_BITS) - 1);

_defineProperty(RecastConstants, "RC_BORDER_REG", 0x8000);

_defineProperty(RecastConstants, "RC_MULTIPLE_REGS", 0);

_defineProperty(RecastConstants, "RC_BORDER_VERTEX", 0x10000);

_defineProperty(RecastConstants, "RC_AREA_BORDER", 0x20000);

_defineProperty(RecastConstants, "RC_CONTOUR_REG_MASK", 0xffff);

_defineProperty(RecastConstants, "RC_MESH_NULL_IDX", 0xffff);

_defineProperty(RecastConstants, "RC_CONTOUR_TESS_WALL_EDGES", 0x01);

_defineProperty(RecastConstants, "RC_CONTOUR_TESS_AREA_EDGES", 0x02);

_defineProperty(RecastConstants, "WATERSHED", 10);

_defineProperty(RecastConstants, "MONOTONE", 20);

_defineProperty(RecastConstants, "LAYERS", 30);

_defineProperty(RecastConstants, "RC_LOG_WARNING", 1);

var _default = RecastConstants;
exports["default"] = _default;