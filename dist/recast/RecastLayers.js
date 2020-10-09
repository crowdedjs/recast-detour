"use strict";

var _temp;

function _createForOfIteratorHelper(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var RecastLayers = /*#__PURE__*/function () {
  function RecastLayers() {
    _classCallCheck(this, RecastLayers);
  }

  _createClass(RecastLayers, null, [{
    key: "addUnique",
    value: function addUnique(a, v) {
      if (!a.contains(v)) {
        a.push(v);
      }
    }
  }, {
    key: "contains",
    value: function contains(a, v) {
      return a.contains(v);
    }
  }, {
    key: "overlapRange",
    value: function overlapRange(amin, amax, bmin, bmax) {
      return amin > bmax || amax < bmin ? false : true;
    }
  }, {
    key: "buildHeightfieldLayers",
    value: function buildHeightfieldLayers(ctx, chf, borderSize, walkableHeight) {
      ctx.startTimer("RC_TIMER_BUILD_LAYERS");
      var w = chf.width;
      var h = chf.height;
      var srcReg = new Array(chf.spanCount);
      Arrays.fill(srcReg, 0xFF);
      var nsweeps = chf.width; // Math.max(chf.width, chf.height);

      var sweeps = new SweepSpan[nsweeps]();

      for (var i = 0; i < sweeps.length; i++) {
        sweeps[i] = new SweepSpan();
      } // Partition walkable area into monotone regions.


      var prevCount = new Array(256);
      var regId = 0; // Sweep one line at a time.

      for (var y = borderSize; y < h - borderSize; ++y) {
        // Collect spans from this row.
        Arrays.fill(prevCount, 0, regId, 0);
        var sweepId = 0;

        for (var x = borderSize; x < w - borderSize; ++x) {
          var c = chf.cells[x + y * w];

          for (var _i = c.index, ni = c.index + c.count; _i < ni; ++_i) {
            var s = chf.spans[_i];
            if (chf.areas[_i] == RecastConstants.RC_NULL_AREA) continue;
            var sid = 0xFF; // -x

            if (GetCon(s, 0) != RecastConstants.RC_NOT_CONNECTED) {
              var ax = x + RecastCommon.GetDirOffsetX(0);
              var ay = y + RecastCommon.GetDirOffsetY(0);
              var ai = chf.cells[ax + ay * w].index + GetCon(s, 0);
              if (chf.areas[ai] != RecastConstants.RC_NULL_AREA && srcReg[ai] != 0xff) sid = srcReg[ai];
            }

            if (sid == 0xff) {
              sid = sweepId++;
              sweeps[sid].nei = 0xff;
              sweeps[sid].ns = 0;
            } // -y


            if (GetCon(s, 3) != RecastConstants.RC_NOT_CONNECTED) {
              var _ax = x + RecastCommon.GetDirOffsetX(3);

              var _ay = y + RecastCommon.GetDirOffsetY(3);

              var _ai = chf.cells[_ax + _ay * w].index + GetCon(s, 3);

              var nr = srcReg[_ai];

              if (nr != 0xff) {
                // Set neighbour when first valid neighbour is
                // encoutered.
                if (sweeps[sid].ns == 0) sweeps[sid].nei = nr;

                if (sweeps[sid].nei == nr) {
                  // Update existing neighbour
                  sweeps[sid].ns++;
                  prevCount[nr]++;
                } else {
                  // This is hit if there is nore than one
                  // neighbour.
                  // Invalidate the neighbour.
                  sweeps[sid].nei = 0xff;
                }
              }
            }

            srcReg[_i] = sid;
          }
        } // Create unique ID.


        for (var _i2 = 0; _i2 < sweepId; ++_i2) {
          // If the neighbour is set and there is only one continuous
          // connection to it,
          // the sweep will be merged with the previous one, else new
          // region is created.
          if (sweeps[_i2].nei != 0xff && prevCount[sweeps[_i2].nei] == sweeps[_i2].ns) {
            sweeps[_i2].id = sweeps[_i2].nei;
          } else {
            if (regId == 255) {
              throw new RuntimeException("rcBuildHeightfieldLayers: Region ID overflow.");
            }

            sweeps[_i2].id = regId++;
          }
        } // Remap local sweep ids to region ids.


        for (var _x = borderSize; _x < w - borderSize; ++_x) {
          var _c = chf.cells[_x + y * w];

          for (var _i3 = _c.index, _ni = _c.index + _c.count; _i3 < _ni; ++_i3) {
            if (srcReg[_i3] != 0xff) srcReg[_i3] = sweeps[srcReg[_i3]].id;
          }
        }
      }

      var nregs = regId;
      var regs = new LayerRegion[nregs](); // Construct regions

      for (var _i4 = 0; _i4 < nregs; ++_i4) {
        regs[_i4] = new LayerRegion(_i4);
      } // Find region neighbours and overlapping regions.


      var lregs = [];

      for (var _y = 0; _y < h; ++_y) {
        for (var _x2 = 0; _x2 < w; ++_x2) {
          var _c2 = chf.cells[_x2 + _y * w];
          lregs = [];

          for (var _i5 = _c2.index, _ni2 = _c2.index + _c2.count; _i5 < _ni2; ++_i5) {
            var _s = chf.spans[_i5];
            var ri = srcReg[_i5];
            if (ri == 0xff) continue;
            regs[ri].ymin = Math.min(regs[ri].ymin, _s.y);
            regs[ri].ymax = Math.max(regs[ri].ymax, _s.y); // Collect all region layers.

            lregs.push(ri); // Update neighbours

            for (var dir = 0; dir < 4; ++dir) {
              if (GetCon(_s, dir) != RecastConstants.RC_NOT_CONNECTED) {
                var _ax2 = _x2 + RecastCommon.GetDirOffsetX(dir);

                var _ay2 = _y + RecastCommon.GetDirOffsetY(dir);

                var _ai2 = chf.cells[_ax2 + _ay2 * w].index + GetCon(_s, dir);

                var rai = srcReg[_ai2];
                if (rai != 0xff && rai != ri) addUnique(regs[ri].neis, rai);
              }
            }
          } // Update overlapping regions.


          for (var _i6 = 0; _i6 < lregs.length - 1; ++_i6) {
            for (var j = _i6 + 1; j < lregs.length; ++j) {
              if (lregs[_i6].intValue() != lregs[j].intValue()) {
                var _ri = regs[lregs[_i6]];
                var rj = regs[lregs[j]];
                addUnique(_ri.layers, lregs[j]);
                addUnique(rj.layers, lregs[_i6]);
              }
            }
          }
        }
      } // Create 2D layers from regions.


      var layerId = 0;
      var stack = [];

      for (var _i7 = 0; _i7 < nregs; ++_i7) {
        var root = regs[_i7]; // Skip already visited.

        if (root.layerId != 0xff) continue; // Start search.

        root.layerId = layerId;
        root.base = true;
        stack.push(_i7);

        while (!stack.length == 0) {
          // Pop front
          var reg = regs[stack.remove(0)];

          var _iterator = _createForOfIteratorHelper(reg.neis),
              _step;

          try {
            for (_iterator.s(); !(_step = _iterator.n()).done;) {
              var nei = _step.value;
              var regn = regs[nei]; // Skip already visited.

              if (regn.layerId != 0xff) continue; // Skip if the neighbour is overlapping root region.

              if (contains(root.layers, nei)) continue; // Skip if the height range would become too large.

              var ymin = Math.min(root.ymin, regn.ymin);
              var ymax = Math.max(root.ymax, regn.ymax);
              if (ymax - ymin >= 255) continue; // Deepen

              stack.push(nei); // Mark layer id

              regn.layerId = layerId; // Merge current layers to root.

              var _iterator2 = _createForOfIteratorHelper(regn.layers),
                  _step2;

              try {
                for (_iterator2.s(); !(_step2 = _iterator2.n()).done;) {
                  var layer = _step2.value;
                  addUnique(root.layers, layer);
                }
              } catch (err) {
                _iterator2.e(err);
              } finally {
                _iterator2.f();
              }

              root.ymin = Math.min(root.ymin, regn.ymin);
              root.ymax = Math.max(root.ymax, regn.ymax);
            }
          } catch (err) {
            _iterator.e(err);
          } finally {
            _iterator.f();
          }
        }

        layerId++;
      } // Merge non-overlapping regions that are close in height.


      var mergeHeight = walkableHeight * 4;

      for (var _i8 = 0; _i8 < nregs; ++_i8) {
        var _ri2 = regs[_i8];
        if (!_ri2.base) continue;
        var newId = _ri2.layerId;

        for (;;) {
          var oldId = 0xff;

          for (var _j = 0; _j < nregs; ++_j) {
            if (_i8 == _j) continue;
            var _rj = regs[_j];
            if (!_rj.base) continue; // Skip if the regions are not close to each other.

            if (!overlapRange(_ri2.ymin, _ri2.ymax + mergeHeight, _rj.ymin, _rj.ymax + mergeHeight)) continue; // Skip if the height range would become too large.

            var _ymin = Math.min(_ri2.ymin, _rj.ymin);

            var _ymax = Math.max(_ri2.ymax, _rj.ymax);

            if (_ymax - _ymin >= 255) continue; // Make sure that there is no overlap when merging 'ri' and
            // 'rj'.

            var overlap = false; // Iterate over all regions which have the same layerId as
            // 'rj'

            for (var k = 0; k < nregs; ++k) {
              if (regs[k].layerId != _rj.layerId) continue; // Check if region 'k' is overlapping region 'ri'
              // Index to 'regs' is the same as region id.

              if (contains(_ri2.layers, k)) {
                overlap = true;
                break;
              }
            } // Cannot merge of regions overlap.


            if (overlap) continue; // Can merge i and j.

            oldId = _rj.layerId;
            break;
          } // Could not find anything to merge with, stop.


          if (oldId == 0xff) break; // Merge

          for (var _j2 = 0; _j2 < nregs; ++_j2) {
            var _rj2 = regs[_j2];

            if (_rj2.layerId == oldId) {
              _rj2.base = false; // Remap layerIds.

              _rj2.layerId = newId; // Add overlaid layers from 'rj' to 'ri'.

              var _iterator3 = _createForOfIteratorHelper(_rj2.layers),
                  _step3;

              try {
                for (_iterator3.s(); !(_step3 = _iterator3.n()).done;) {
                  var _layer = _step3.value;
                  addUnique(_ri2.layers, _layer);
                } // Update height bounds.

              } catch (err) {
                _iterator3.e(err);
              } finally {
                _iterator3.f();
              }

              _ri2.ymin = Math.min(_ri2.ymin, _rj2.ymin);
              _ri2.ymax = Math.max(_ri2.ymax, _rj2.ymax);
            }
          }
        }
      } // Compact layerIds


      var remap = new Array(256); // Find number of unique layers.

      layerId = 0;

      for (var _i9 = 0; _i9 < nregs; ++_i9) {
        remap[regs[_i9].layerId] = 1;
      }

      for (var _i10 = 0; _i10 < 256; ++_i10) {
        if (remap[_i10] != 0) remap[_i10] = layerId++;else remap[_i10] = 0xff;
      } // Remap ids.


      for (var _i11 = 0; _i11 < nregs; ++_i11) {
        regs[_i11].layerId = remap[regs[_i11].layerId];
      } // No layers, return empty.


      if (layerId == 0) {
        // ctx.stopTimer(RC_TIMER_BUILD_LAYERS);
        return null;
      } // Create layers.
      // rcAssert(lset.layers == 0);


      var lw = w - borderSize * 2;
      var lh = h - borderSize * 2; // Build contracted bbox for layers.

      var bmin = new Array(3);
      var bmax = new Array(3);
      copy(bmin, chf.bmin);
      copy(bmax, chf.bmax);
      bmin[0] += borderSize * chf.cs;
      bmin[2] += borderSize * chf.cs;
      bmax[0] -= borderSize * chf.cs;
      bmax[2] -= borderSize * chf.cs;
      var lset = new HeightfieldLayerSet();
      lset.layers = new HeightfieldLayer[layerId]();

      for (var _i12 = 0; _i12 < lset.layers.length; _i12++) {
        lset.layers[_i12] = new HeightfieldLayer();
      } // Store layers.


      for (var _i13 = 0; _i13 < lset.layers.length; ++_i13) {
        var curId = _i13;
        var _layer2 = lset.layers[_i13];
        var gridSize = lw * lh;
        _layer2.heights = new Array(gridSize);
        Arrays.fill(_layer2.heights, 0xFF);
        _layer2.areas = new Array(gridSize);
        _layer2.cons = new Array(gridSize); // Find layer height bounds.

        var hmin = 0,
            hmax = 0;

        for (var _j3 = 0; _j3 < nregs; ++_j3) {
          if (regs[_j3].base && regs[_j3].layerId == curId) {
            hmin = regs[_j3].ymin;
            hmax = regs[_j3].ymax;
          }
        }

        _layer2.width = lw;
        _layer2.height = lh;
        _layer2.cs = chf.cs;
        _layer2.ch = chf.ch; // Adjust the bbox to fit the heightfield.

        copy(_layer2.bmin, bmin);
        copy(_layer2.bmax, bmax);
        _layer2.bmin[1] = bmin[1] + hmin * chf.ch;
        _layer2.bmax[1] = bmin[1] + hmax * chf.ch;
        _layer2.hmin = hmin;
        _layer2.hmax = hmax; // Update usable data region.

        _layer2.minx = _layer2.width;
        _layer2.maxx = 0;
        _layer2.miny = _layer2.height;
        _layer2.maxy = 0; // Copy height and area from compact heightfield.

        for (var _y2 = 0; _y2 < lh; ++_y2) {
          for (var _x3 = 0; _x3 < lw; ++_x3) {
            var cx = borderSize + _x3;
            var cy = borderSize + _y2;
            var _c3 = chf.cells[cx + cy * w];

            for (var _j4 = _c3.index, nj = _c3.index + _c3.count; _j4 < nj; ++_j4) {
              var _s2 = chf.spans[_j4]; // Skip unassigned regions.

              if (srcReg[_j4] == 0xff) continue; // Skip of does nto bePoly to current layer.

              var lid = regs[srcReg[_j4]].layerId;
              if (lid != curId) continue; // Update data bounds.

              _layer2.minx = Math.min(_layer2.minx, _x3);
              _layer2.maxx = Math.max(_layer2.maxx, _x3);
              _layer2.miny = Math.min(_layer2.miny, _y2);
              _layer2.maxy = Math.max(_layer2.maxy, _y2); // Store height and area type.

              var idx = _x3 + _y2 * lw;
              _layer2.heights[idx] = char(_s2.y - hmin);
              _layer2.areas[idx] = chf.areas[_j4]; // Check connection.

              var portal = 0;
              var con = 0;

              for (var _dir = 0; _dir < 4; ++_dir) {
                if (GetCon(_s2, _dir) != RecastConstants.RC_NOT_CONNECTED) {
                  var _ax3 = cx + RecastCommon.GetDirOffsetX(_dir);

                  var _ay3 = cy + RecastCommon.GetDirOffsetY(_dir);

                  var _ai3 = chf.cells[_ax3 + _ay3 * w].index + GetCon(_s2, _dir);

                  var alid = srcReg[_ai3] != 0xff ? regs[srcReg[_ai3]].layerId : 0xff; // Portal mask

                  if (chf.areas[_ai3] != RecastConstants.RC_NULL_AREA && lid != alid) {
                    portal |= 1 << _dir; // Update height so that it matches on both
                    // sides of the portal.

                    var as = chf.spans[_ai3];
                    if (as.y > hmin) _layer2.heights[idx] = Math.max(_layer2.heights[idx], char(as.y - hmin));
                  } // Valid connection mask


                  if (chf.areas[_ai3] != RecastConstants.RC_NULL_AREA && lid == alid) {
                    var nx = _ax3 - borderSize;
                    var ny = _ay3 - borderSize;
                    if (nx >= 0 && ny >= 0 && nx < lw && ny < lh) con |= 1 << _dir;
                  }
                }
              }

              _layer2.cons[idx] = portal << 4 | con;
            }
          }
        }

        if (_layer2.minx > _layer2.maxx) _layer2.minx = _layer2.maxx = 0;
        if (_layer2.miny > _layer2.maxy) _layer2.miny = _layer2.maxy = 0;
      } // ctx=>stopTimer(RC_TIMER_BUILD_LAYERS);


      return lset;
    }
  }]);

  return RecastLayers;
}();

_defineProperty(RecastLayers, "RC_MAX_LAYERS", RecastConstants.RC_NOT_CONNECTED);

_defineProperty(RecastLayers, "RC_MAX_NEIS", 16);

_defineProperty(RecastLayers, "LayerRegion", (_temp = function LayerRegion(i) {
  _classCallCheck(this, LayerRegion);

  _defineProperty(this, "id", void 0);

  _defineProperty(this, "layerId", void 0);

  _defineProperty(this, "base", void 0);

  _defineProperty(this, "ymin", void 0);

  _defineProperty(this, "ymax", void 0);

  _defineProperty(this, "layers", void 0);

  _defineProperty(this, "neis", void 0);

  this.id = i;
  this.ymin = 0xFFFF;
  this.layerId = 0xff;
  this.layers = [];
  this.neis = [];
}, _temp));