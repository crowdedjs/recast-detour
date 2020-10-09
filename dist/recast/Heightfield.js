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

/** Represents a heightfield layer within a layer set. */
var Heightfield =
/** The width of the heightfield. (APoly the x-axis in cell units.) */

/** The height of the heightfield. (APoly the z-axis in cell units.) */

/** The minimum bounds in world space. [(x, y, z)] */

/** The maximum bounds in world space. [(x, y, z)] */

/** The size of each cell. (On the xz-plane.) */

/** The height of each cell. (The minimum increment aPoly the y-axis.) */

/** Heightfield of spans (width*height). */
function Heightfield(width, height, bmin, bmax, cs, ch) {
  _classCallCheck(this, Heightfield);

  _defineProperty(this, "width", 0);

  _defineProperty(this, "height", 0);

  _defineProperty(this, "bmin", void 0);

  _defineProperty(this, "bmax", void 0);

  _defineProperty(this, "cs", void 0);

  _defineProperty(this, "ch", void 0);

  _defineProperty(this, "spans", []);

  this.width = width;
  this.height = height;
  this.bmin = bmin;
  this.bmax = bmax;
  this.cs = cs;
  this.ch = ch;
  this.spans = new Array(width * height);
};

var _default = Heightfield;
exports["default"] = _default;