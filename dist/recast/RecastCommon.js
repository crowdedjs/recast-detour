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
var RecastCommon = function RecastCommon() {
  _classCallCheck(this, RecastCommon);
};

_defineProperty(RecastCommon, "GetCon", function (s, dir) {
  var shift = dir * 6;
  return s.con >> shift & 0x3f;
});

_defineProperty(RecastCommon, "GetDirOffsetX", function (dir) {
  var offset = [-1, 0, 1, 0];
  return offset[dir & 0x03];
});

_defineProperty(RecastCommon, "GetDirOffsetY", function (dir) {
  var offset = [0, 1, 0, -1];
  return offset[dir & 0x03];
});

_defineProperty(RecastCommon, "rcGetDirForOffset", function (x, y) {
  var dirs = [3, 0, -1, 2, 1];
  return dirs[(y + 1 << 1) + x];
});

_defineProperty(RecastCommon, "SetCon", function (s, dir, i) {
  var shift = dir * 6;
  var con = s.con;
  s.con = con & ~(0x3f << shift) | (i & 0x3f) << shift;
});

_defineProperty(RecastCommon, "clamp", function (v, min, max) {
  return Math.max(Math.min(max, v), min);
});

var _default = RecastCommon;
exports["default"] = _default;