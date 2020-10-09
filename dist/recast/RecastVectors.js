"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

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
 _in a product, an acknowledgment _in the product documentation would be
 appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
 misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
var RecastVectors = /*#__PURE__*/function () {
  function RecastVectors() {
    _classCallCheck(this, RecastVectors);
  }

  _createClass(RecastVectors, null, [{
    key: "min",
    value: function min(a, b, i) {
      a[0] = Math.min(a[0], b[i + 0]);
      a[1] = Math.min(a[1], b[i + 1]);
      a[2] = Math.min(a[2], b[i + 2]);
    }
  }, {
    key: "max",
    value: function max(a, b, i) {
      a[0] = Math.max(a[0], b[i + 0]);
      a[1] = Math.max(a[1], b[i + 1]);
      a[2] = Math.max(a[2], b[i + 2]);
    }
  }, {
    key: "copy3",
    value: function copy3(out, _in, i) {
      RecastVectors.copy4(out, 0, _in, i);
    }
  }, {
    key: "copy2",
    value: function copy2(out, _in) {
      RecastVectors.copy4(out, 0, _in, 0);
    }
  }, {
    key: "copy4",
    value: function copy4(out, n, _in, m) {
      out[n] = _in[m];
      out[n + 1] = _in[m + 1];
      out[n + 2] = _in[m + 2];
    }
  }, {
    key: "add",
    value: function add(e0, a, verts, i) {
      e0[0] = a[0] + verts[i];
      e0[1] = a[1] + verts[i + 1];
      e0[2] = a[2] + verts[i + 2];
    }
  }, {
    key: "subA",
    value: function subA(e0, verts, i, j) {
      e0[0] = verts[i] - verts[j];
      e0[1] = verts[i + 1] - verts[j + 1];
      e0[2] = verts[i + 2] - verts[j + 2];
    }
  }, {
    key: "subB",
    value: function subB(e0, i, verts, j) {
      e0[0] = i[0] - verts[j];
      e0[1] = i[1] - verts[j + 1];
      e0[2] = i[2] - verts[j + 2];
    }
  }, {
    key: "cross",
    value: function cross(dest, v1, v2) {
      dest[0] = v1[1] * v2[2] - v1[2] * v2[1];
      dest[1] = v1[2] * v2[0] - v1[0] * v2[2];
      dest[2] = v1[0] * v2[1] - v1[1] * v2[0];
    }
  }, {
    key: "normalize",
    value: function normalize(v) {
      var d = 1.0 / Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
      v[0] *= d;
      v[1] *= d;
      v[2] *= d;
    }
  }]);

  return RecastVectors;
}();

var _default = RecastVectors;
exports["default"] = _default;