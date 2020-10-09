"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = void 0;

var _temp, _temp2;

function _createForOfIteratorHelper(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

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
var ChunkyTriMesh = /*#__PURE__*/function () {
  _createClass(ChunkyTriMesh, [{
    key: "calcExtends",
    value: function calcExtends(items, imin, imax, bmin, bmax) {
      bmin[0] = items[imin].bmin[0];
      bmin[1] = items[imin].bmin[1];
      bmax[0] = items[imin].bmax[0];
      bmax[1] = items[imin].bmax[1];

      for (var i = imin + 1; i < imax; ++i) {
        var it = items[i];

        if (it.bmin[0] < bmin[0]) {
          bmin[0] = it.bmin[0];
        }

        if (it.bmin[1] < bmin[1]) {
          bmin[1] = it.bmin[1];
        }

        if (it.bmax[0] > bmax[0]) {
          bmax[0] = it.bmax[0];
        }

        if (it.bmax[1] > bmax[1]) {
          bmax[1] = it.bmax[1];
        }
      }
    }
  }, {
    key: "longestAxis",
    value: function longestAxis(x, y) {
      return y > x ? 1 : 0;
    } // https://stackoverflow.com/a/45245772/10047920

  }, {
    key: "partialSort",
    value: function partialSort(arr, start, end, sortFx) {
      var preSorted = arr.slice(0, start),
          postSorted = arr.slice(end);
      var sorted = arr.slice(start, end).sort(sortFx);
      arr.length = 0;
      arr.push.apply(arr, preSorted.concat(sorted).concat(postSorted));
      return arr;
    }
  }, {
    key: "subdivide",
    value: function subdivide(items, imin, imax, trisPerChunk, nodes, inTris) {
      var _this = this;

      var inum = imax - imin;
      var node = new ChunkyTriMesh.ChunkyTriMeshNode();
      nodes.push(node);

      if (inum <= trisPerChunk) {
        // Leaf
        this.calcExtends(items, imin, imax, node.bmin, node.bmax); // Copy triangles.

        node.i = nodes.length;
        node.tris = new Array(inum * 3);
        var dst = 0;

        for (var i = imin; i < imax; ++i) {
          var src = items[i].i * 3;
          node.tris[dst++] = inTris[src];
          node.tris[dst++] = inTris[src + 1];
          node.tris[dst++] = inTris[src + 2];
        }
      } else {
        // Split
        this.calcExtends(items, imin, imax, node.bmin, node.bmax);
        var axis = this.longestAxis(node.bmax[0] - node.bmin[0], node.bmax[1] - node.bmin[1]);

        if (axis == 0) {
          //Arrays.sort(items, imin, imax, new CompareItemX());
          this.partialSort(items, imin, imax, function (a, b) {
            return _this.CompareItemX(a, b);
          }); // let subarray = items.slice(imin, imax);
          // subarray.sort((a, b) => (this.CompareItemX(a, b)))
          // // let newitems = items.slice(0, imin);
          // // newitems.push(...subarray)
          // // newitems.push(...items.slice(imax, items.length));
          // for (let i = imin; i < imax; i++) {
          //     items[i].bmin[0] = subarray[i].bmin[0];
          //     items[i].bmin[1] = subarray[i].bmin[1];
          //     items[i].bmax[0] = subarray[i].bmax[0];
          //     items[i].bmax[1] = subarray[i].bmax[1];
          // }
          //items = newitems;
          // Sort aPoly x-axis
        } else if (axis == 1) {
          // Arrays.sort(items, imin, imax, new CompareItemY());
          // Sort aPoly y-axis
          this.partialSort(items, imin, imax, function (a, b) {
            return _this.CompareItemY(a, b);
          }); // let subarray = items.slice(imin, imax);
          // subarray.sort((a, b) => (this.CompareItemY(a, b)))
          // // let newitems = items.slice(0, imin);
          // // newitems.push(...subarray)
          // // newitems.push(...items.slice(imax, items.length));
          // for (let i = imin; i < imax; i++) {
          //     items[i].bmin[0] = subarray[i].bmin[0];
          //     items[i].bmin[1] = subarray[i].bmin[1];
          //     items[i].bmax[0] = subarray[i].bmax[0];
          //     items[i].bmax[1] = subarray[i].bmax[1];
          // }
          //items = newitems;
        }

        var isplit = Math.floor(imin + Math.floor(inum / 2));
        var s = "";

        for (var _i = 0; _i < items.length; _i++) {
          var item = items[_i];
          s += "" + _i + " " + item.bmin[0] + " " + item.bmin[1] + "\n";
          s += "" + _i + " " + item.bmax[0] + " " + item.bmax[1] + "\n";
        } // fs.writeFileSync("items_" + imin + "_" + imax + ".txt", s);


        console.log("done"); // Left

        console.log("Before left " + imin + " " + isplit);
        console.log(items[279].bmin[1]);
        this.subdivide(items, imin, isplit, trisPerChunk, nodes, inTris);
        console.log("Done left " + imin + " " + isplit);
        console.log(items[279].bmin[1]); // Right

        console.log("Before right " + isplit + " " + imax);
        console.log(items[279].bmin[1]);
        this.subdivide(items, isplit, imax, trisPerChunk, nodes, inTris);
        console.log("Done right " + isplit + " " + imax);
        console.log(items[279].bmin[1]); // Negative index means escape.

        node.i = -nodes.length;
      }
    }
  }]);

  function ChunkyTriMesh(verts, tris, ntris, trisPerChunk) {
    _classCallCheck(this, ChunkyTriMesh);

    _defineProperty(this, "CompareItemX", function compare(a, b) {
      if (a.bmin[0] < b.bmin[0]) {
        return -1;
      }

      if (a.bmin[0] > b.bmin[0]) {
        return 1;
      }

      return 0;
    });

    _defineProperty(this, "CompareItemY", function compare(a, b) {
      if (a.bmin[1] < b.bmin[1]) {
        return -1;
      }

      if (a.bmin[1] > b.bmin[1]) {
        return 1;
      }

      return 0;
    });

    _defineProperty(this, "nodes", []);

    _defineProperty(this, "ntris", void 0);

    _defineProperty(this, "maxTrisPerChunk", void 0);

    var nchunks = Math.floor((ntris + trisPerChunk - 1) / trisPerChunk);
    this.nodes = [];
    this.ntris = ntris; // Build tree

    var items = []; //new Array(ntris);

    for (var i = 0; i < ntris; i++) {
      var t = i * 3;
      var it = items[i] = new ChunkyTriMesh.BoundsItem();
      it.i = i; // Calc triangle XZ bounds.

      it.bmin[0] = it.bmax[0] = verts[tris[t] * 3 + 0];
      it.bmin[1] = it.bmax[1] = verts[tris[t] * 3 + 2];

      for (var j = 1; j < 3; ++j) {
        var v = tris[t + j] * 3;

        if (verts[v] < it.bmin[0]) {
          it.bmin[0] = verts[v];
        }

        if (verts[v + 2] < it.bmin[1]) {
          it.bmin[1] = verts[v + 2];
        }

        if (verts[v] > it.bmax[0]) {
          it.bmax[0] = verts[v];
        }

        if (verts[v + 2] > it.bmax[1]) {
          it.bmax[1] = verts[v + 2];
        }
      }
    }

    this.subdivide(items, 0, ntris, trisPerChunk, this.nodes, tris); // Calc max tris per node.

    this.maxTrisPerChunk = 0;

    var _iterator = _createForOfIteratorHelper(this.nodes),
        _step;

    try {
      for (_iterator.s(); !(_step = _iterator.n()).done;) {
        var _node = _step.value;
        var isLeaf = _node.i >= 0;

        if (!isLeaf) {
          continue;
        }

        if (_node.tris.length / 3 > this.maxTrisPerChunk) {
          this.maxTrisPerChunk = _node.tris.length / 3;
        }
      }
    } catch (err) {
      _iterator.e(err);
    } finally {
      _iterator.f();
    }
  }

  _createClass(ChunkyTriMesh, [{
    key: "checkOverlapRect",
    value: function checkOverlapRect(amin, amax, bmin, bmax) {
      var overlap = true;
      overlap = amin[0] > bmax[0] || amax[0] < bmin[0] ? false : overlap;
      overlap = amin[1] > bmax[1] || amax[1] < bmin[1] ? false : overlap;
      return overlap;
    }
  }, {
    key: "getChunksOverlappingRect",
    value: function getChunksOverlappingRect(bmin, bmax) {
      // Traverse tree
      ids = [];
      var i = 0;

      while (i < nodes.length) {
        node = nodes[i];
        var overlap = checkOverlapRect(bmin, bmax, node.bmin, node.bmax);
        var isLeafNode = node.i >= 0;

        if (isLeafNode && overlap) {
          ids.push(node);
        }

        if (overlap || isLeafNode) {
          i++;
        } else {
          i = -node.i;
        }
      }

      return ids;
    }
  }]);

  return ChunkyTriMesh;
}();

_defineProperty(ChunkyTriMesh, "BoundsItem", (_temp = function BoundsItem() {
  _classCallCheck(this, BoundsItem);

  _defineProperty(this, "bmin", new Array(2));

  _defineProperty(this, "bmax", new Array(2));

  _defineProperty(this, "i", void 0);
}, _temp));

_defineProperty(ChunkyTriMesh, "ChunkyTriMeshNode", (_temp2 = function ChunkyTriMeshNode() {
  _classCallCheck(this, ChunkyTriMeshNode);

  _defineProperty(this, "bmin", new Array(2));

  _defineProperty(this, "bmax", new Array(2));

  _defineProperty(this, "i", void 0);

  _defineProperty(this, "tris", []);
}, _temp2));

var _default = ChunkyTriMesh;
exports["default"] = _default;