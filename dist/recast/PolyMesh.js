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

/** Represents a polygon mesh suitable for use in building a navigation mesh. */
var PolyMesh = function PolyMesh() {
  _classCallCheck(this, PolyMesh);

  _defineProperty(this, "verts", []);

  _defineProperty(this, "polys", []);

  _defineProperty(this, "regs", []);

  _defineProperty(this, "areas", []);

  _defineProperty(this, "nverts", 0);

  _defineProperty(this, "npolys", 0);

  _defineProperty(this, "nvp", 0);

  _defineProperty(this, "maxpolys", 0);

  _defineProperty(this, "flags", []);

  _defineProperty(this, "bmin", new Array(3));

  _defineProperty(this, "bmax", new Array(3));

  _defineProperty(this, "cs", 0);

  _defineProperty(this, "ch", 0);

  _defineProperty(this, "borderSize", 0);

  _defineProperty(this, "maxEdgeError", 0);
};

var _default = PolyMesh;
exports["default"] = _default;