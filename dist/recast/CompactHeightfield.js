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

/** A compact, static heightfield representing unobstructed space. */
var CompactHeightfield = function CompactHeightfield() {
  _classCallCheck(this, CompactHeightfield);

  _defineProperty(this, "width", 0);

  _defineProperty(this, "height", 0);

  _defineProperty(this, "spanCount", 0);

  _defineProperty(this, "walkableHeight", 0);

  _defineProperty(this, "walkableClimb", 0);

  _defineProperty(this, "borderSize", 0);

  _defineProperty(this, "maxDistance", 0);

  _defineProperty(this, "maxRegions", void 0);

  _defineProperty(this, "bmin", new Array(3));

  _defineProperty(this, "bmax", new Array(3));

  _defineProperty(this, "cs", 0);

  _defineProperty(this, "ch", 0);

  _defineProperty(this, "cells", []);

  _defineProperty(this, "spans", []);

  _defineProperty(this, "dist", []);

  _defineProperty(this, "areas", []);
};

var _default = CompactHeightfield;
exports["default"] = _default;