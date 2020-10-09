"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _RecastConstants = _interopRequireDefault(require("./RecastConstants.js"));

var _RecastCommon = _interopRequireDefault(require("./RecastCommon.js"));

var _temp, _temp2;

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var RecastRegion = /*#__PURE__*/function () {
  function RecastRegion() {
    _classCallCheck(this, RecastRegion);
  }

  _createClass(RecastRegion, null, [{
    key: "calculateDistanceField",
    value: function calculateDistanceField(chf, src) {
      var maxDist;
      var w = chf.width;
      var h = chf.height; // Init distance and points.

      for (var i = 0; i < chf.spanCount; ++i) {
        src[i] = 0xffff;
      } // Mark boundary cells.


      for (var y = 0; y < h; ++y) {
        for (var _x = 0; _x < w; ++_x) {
          var c = chf.cells[_x + y * w];

          for (var _i = c.index, ni = c.index + c.count; _i < ni; ++_i) {
            var s = chf.spans[_i];
            var area = chf.areas[_i];
            var nc = 0;

            for (var dir = 0; dir < 4; ++dir) {
              if (_RecastCommon["default"].GetCon(s, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                var ax = _x + _RecastCommon["default"].GetDirOffsetX(dir);

                var ay = y + _RecastCommon["default"].GetDirOffsetY(dir);

                var ai = chf.cells[ax + ay * w].index + _RecastCommon["default"].GetCon(s, dir);

                if (area == chf.areas[ai]) {
                  nc++;
                }
              }
            }

            if (nc != 4) {
              src[_i] = 0;
            }
          }
        }
      } // Pass 1


      for (var _y = 0; _y < h; ++_y) {
        for (var _x2 = 0; _x2 < w; ++_x2) {
          var _c = chf.cells[_x2 + _y * w];

          for (var _i2 = _c.index, _ni = _c.index + _c.count; _i2 < _ni; ++_i2) {
            var _s = chf.spans[_i2];

            if (_RecastCommon["default"].GetCon(_s, 0) != _RecastConstants["default"].RC_NOT_CONNECTED) {
              // (-1,0)
              var _ax = _x2 + _RecastCommon["default"].GetDirOffsetX(0);

              var _ay = _y + _RecastCommon["default"].GetDirOffsetY(0);

              var _ai = chf.cells[_ax + _ay * w].index + _RecastCommon["default"].GetCon(_s, 0);

              var as = chf.spans[_ai];

              if (src[_ai] + 2 < src[_i2]) {
                src[_i2] = src[_ai] + 2;
              } // (-1,-1)


              if (_RecastCommon["default"].GetCon(as, 3) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                var aax = _ax + _RecastCommon["default"].GetDirOffsetX(3);

                var aay = _ay + _RecastCommon["default"].GetDirOffsetY(3);

                var aai = chf.cells[aax + aay * w].index + _RecastCommon["default"].GetCon(as, 3);

                if (src[aai] + 3 < src[_i2]) {
                  src[_i2] = src[aai] + 3;
                }
              }
            }

            if (_RecastCommon["default"].GetCon(_s, 3) != _RecastConstants["default"].RC_NOT_CONNECTED) {
              // (0,-1)
              var _ax2 = _x2 + _RecastCommon["default"].GetDirOffsetX(3);

              var _ay2 = _y + _RecastCommon["default"].GetDirOffsetY(3);

              var _ai2 = chf.cells[_ax2 + _ay2 * w].index + _RecastCommon["default"].GetCon(_s, 3);

              var _as = chf.spans[_ai2];

              if (src[_ai2] + 2 < src[_i2]) {
                src[_i2] = src[_ai2] + 2;
              } // (1,-1)


              if (_RecastCommon["default"].GetCon(_as, 2) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                var _aax = _ax2 + _RecastCommon["default"].GetDirOffsetX(2);

                var _aay = _ay2 + _RecastCommon["default"].GetDirOffsetY(2);

                var _aai = chf.cells[_aax + _aay * w].index + _RecastCommon["default"].GetCon(_as, 2);

                if (src[_aai] + 3 < src[_i2]) {
                  src[_i2] = src[_aai] + 3;
                }
              }
            }
          }
        }
      } // Pass 2


      for (var _y2 = h - 1; _y2 >= 0; --_y2) {
        for (var _x3 = w - 1; _x3 >= 0; --_x3) {
          var _c2 = chf.cells[_x3 + _y2 * w];

          for (var _i3 = _c2.index, _ni2 = _c2.index + _c2.count; _i3 < _ni2; ++_i3) {
            var _s2 = chf.spans[_i3];

            if (_RecastCommon["default"].GetCon(_s2, 2) != _RecastConstants["default"].RC_NOT_CONNECTED) {
              // (1,0)
              var _ax3 = _x3 + _RecastCommon["default"].GetDirOffsetX(2);

              var _ay3 = _y2 + _RecastCommon["default"].GetDirOffsetY(2);

              var _ai3 = chf.cells[_ax3 + _ay3 * w].index + _RecastCommon["default"].GetCon(_s2, 2);

              var _as2 = chf.spans[_ai3];

              if (src[_ai3] + 2 < src[_i3]) {
                src[_i3] = src[_ai3] + 2;
              } // (1,1)


              if (_RecastCommon["default"].GetCon(_as2, 1) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                var _aax2 = _ax3 + _RecastCommon["default"].GetDirOffsetX(1);

                var _aay2 = _ay3 + _RecastCommon["default"].GetDirOffsetY(1);

                var _aai2 = chf.cells[_aax2 + _aay2 * w].index + _RecastCommon["default"].GetCon(_as2, 1);

                if (src[_aai2] + 3 < src[_i3]) {
                  src[_i3] = src[_aai2] + 3;
                }
              }
            }

            if (_RecastCommon["default"].GetCon(_s2, 1) != _RecastConstants["default"].RC_NOT_CONNECTED) {
              // (0,1)
              var _ax4 = _x3 + _RecastCommon["default"].GetDirOffsetX(1);

              var _ay4 = _y2 + _RecastCommon["default"].GetDirOffsetY(1);

              var _ai4 = chf.cells[_ax4 + _ay4 * w].index + _RecastCommon["default"].GetCon(_s2, 1);

              var _as3 = chf.spans[_ai4];

              if (src[_ai4] + 2 < src[_i3]) {
                src[_i3] = src[_ai4] + 2;
              } // (-1,1)


              if (_RecastCommon["default"].GetCon(_as3, 0) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                var _aax3 = _ax4 + _RecastCommon["default"].GetDirOffsetX(0);

                var _aay3 = _ay4 + _RecastCommon["default"].GetDirOffsetY(0);

                var _aai3 = chf.cells[_aax3 + _aay3 * w].index + _RecastCommon["default"].GetCon(_as3, 0);

                if (src[_aai3] + 3 < src[_i3]) {
                  src[_i3] = src[_aai3] + 3;
                }
              }
            }
          }
        }
      }

      maxDist = 0;

      for (var _i4 = 0; _i4 < chf.spanCount; ++_i4) {
        maxDist = Math.max(src[_i4], maxDist);
      }

      return maxDist;
    }
  }, {
    key: "boxBlur",
    value: function boxBlur(chf, thr, src) {
      var w = chf.width;
      var h = chf.height;
      var dst = new Array(chf.spanCount);
      thr *= 2;

      for (var y = 0; y < h; ++y) {
        for (var _x4 = 0; _x4 < w; ++_x4) {
          var c = chf.cells[_x4 + y * w];

          for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
            var s = chf.spans[i];
            var cd = src[i];

            if (cd <= thr) {
              dst[i] = cd;
              continue;
            }

            var d = cd;

            for (var dir = 0; dir < 4; ++dir) {
              if (_RecastCommon["default"].GetCon(s, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                var ax = _x4 + _RecastCommon["default"].GetDirOffsetX(dir);

                var ay = y + _RecastCommon["default"].GetDirOffsetY(dir);

                var ai = chf.cells[ax + ay * w].index + _RecastCommon["default"].GetCon(s, dir);

                d += src[ai];
                var as = chf.spans[ai];
                var dir2 = dir + 1 & 0x3;

                if (_RecastCommon["default"].GetCon(as, dir2) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                  var ax2 = ax + _RecastCommon["default"].GetDirOffsetX(dir2);

                  var ay2 = ay + _RecastCommon["default"].GetDirOffsetY(dir2);

                  var ai2 = chf.cells[ax2 + ay2 * w].index + _RecastCommon["default"].GetCon(as, dir2);

                  d += src[ai2];
                } else {
                  d += cd;
                }
              } else {
                d += cd * 2;
              }
            }

            dst[i] = Math.floor((d + 5) / 9);
          }
        }
      }

      return dst;
    }
  }, {
    key: "floodRegion",
    value: function floodRegion(x, y, i, level, r, chf, srcReg, srcDist, stack) {
      var w = chf.width;
      var area = chf.areas[i]; // Flood fill mark region.

      stack = [];
      stack.push(x);
      stack.push(y);
      stack.push(i);
      srcReg[i] = r;
      srcDist[i] = 0;
      var lev = level >= 2 ? level - 2 : 0;
      var count = 0;

      while (stack.length > 0) {
        var ci = stack.splice(stack.length - 1, 1)[0];
        var cy = stack.splice(stack.length - 1, 1)[0];
        var cx = stack.splice(stack.length - 1, 1)[0];
        var cs = chf.spans[ci]; // Check if any of the neighbours already have a valid region set.

        var ar = 0;

        for (var dir = 0; dir < 4; ++dir) {
          // 8 connected
          if (_RecastCommon["default"].GetCon(cs, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
            var ax = cx + _RecastCommon["default"].GetDirOffsetX(dir);

            var ay = cy + _RecastCommon["default"].GetDirOffsetY(dir);

            var ai = chf.cells[ax + ay * w].index + _RecastCommon["default"].GetCon(cs, dir);

            if (chf.areas[ai] != area) {
              continue;
            }

            var nr = srcReg[ai];

            if ((nr & _RecastConstants["default"].RC_BORDER_REG) != 0) {
              continue;
            }

            if (nr != 0 && nr != r) {
              ar = nr;
              break;
            }

            var as = chf.spans[ai];
            var dir2 = dir + 1 & 0x3;

            if (_RecastCommon["default"].GetCon(as, dir2) != _RecastConstants["default"].RC_NOT_CONNECTED) {
              var ax2 = ax + _RecastCommon["default"].GetDirOffsetX(dir2);

              var ay2 = ay + _RecastCommon["default"].GetDirOffsetY(dir2);

              var ai2 = chf.cells[ax2 + ay2 * w].index + _RecastCommon["default"].GetCon(as, dir2);

              if (chf.areas[ai2] != area) {
                continue;
              }

              var nr2 = srcReg[ai2];

              if (nr2 != 0 && nr2 != r) {
                ar = nr2;
                break;
              }
            }
          }
        }

        if (ar != 0) {
          srcReg[ci] = 0;
          continue;
        }

        count++; // Expand neighbours.

        for (var _dir = 0; _dir < 4; ++_dir) {
          if (_RecastCommon["default"].GetCon(cs, _dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
            var _ax5 = cx + _RecastCommon["default"].GetDirOffsetX(_dir);

            var _ay5 = cy + _RecastCommon["default"].GetDirOffsetY(_dir);

            var _ai5 = chf.cells[_ax5 + _ay5 * w].index + _RecastCommon["default"].GetCon(cs, _dir);

            if (chf.areas[_ai5] != area) {
              continue;
            }

            if (chf.dist[_ai5] >= lev && srcReg[_ai5] == 0) {
              srcReg[_ai5] = r;
              srcDist[_ai5] = 0;
              stack.push(_ax5);
              stack.push(_ay5);
              stack.push(_ai5);
            }
          }
        }
      }

      return count > 0;
    }
  }, {
    key: "expandRegions",
    value: function expandRegions(maxIter, level, chf, srcReg, srcDist, stack, fillStack) {
      var w = chf.width;
      var h = chf.height;

      if (fillStack) {
        // Find cells revealed by the raised level.
        //stack = [];
        stack = [];

        for (var y = 0; y < h; ++y) {
          for (var _x5 = 0; _x5 < w; ++_x5) {
            var c = chf.cells[_x5 + y * w];

            for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
              if (chf.dist[i] >= level && srcReg[i] == 0 && chf.areas[i] != _RecastConstants["default"].RC_NULL_AREA) {
                stack.push(_x5);
                stack.push(y);
                stack.push(i);
              }
            }
          }
        }
      } else // use cells in the input stack
        {
          // mark all cells which already have a region
          for (var j = 0; j < stack.length; j += 3) {
            var _i5 = stack[j + 2];

            if (srcReg[_i5] != 0) {
              stack[j + 2] = -1;
            }
          }
        }

      var dirtyEntries = [];
      var iter = 0;

      while (stack.length > 0) {
        var failed = 0; // dirtyEntries = [];

        dirtyEntries = [];

        for (var _j = 0; _j < stack.length; _j += 3) {
          var _x6 = stack[_j + 0];
          var _y3 = stack[_j + 1];
          var _i6 = stack[_j + 2];

          if (_i6 < 0) {
            failed++;
            continue;
          }

          var r = srcReg[_i6];
          var d2 = 0xfff;
          var area = chf.areas[_i6];
          var s = chf.spans[_i6];

          for (var dir = 0; dir < 4; ++dir) {
            if (_RecastCommon["default"].GetCon(s, dir) == _RecastConstants["default"].RC_NOT_CONNECTED) {
              continue;
            }

            var ax = _x6 + _RecastCommon["default"].GetDirOffsetX(dir);

            var ay = _y3 + _RecastCommon["default"].GetDirOffsetY(dir);

            var ai = chf.cells[ax + ay * w].index + _RecastCommon["default"].GetCon(s, dir);

            if (chf.areas[ai] != area) {
              continue;
            }

            if (srcReg[ai] > 0 && (srcReg[ai] & _RecastConstants["default"].RC_BORDER_REG) == 0) {
              if (srcDist[ai] + 2 < d2) {
                r = srcReg[ai];
                d2 = srcDist[ai] + 2;
              }
            }
          }

          if (r != 0) {
            stack[_j + 2] = -1; // mark as used

            dirtyEntries.push(_i6);
            dirtyEntries.push(r);
            dirtyEntries.push(d2);
          } else {
            failed++;
          }
        } // Copy entries that differ between src and dst to keep them in sync.


        for (var _i7 = 0; _i7 < dirtyEntries.length; _i7 += 3) {
          var idx = dirtyEntries[_i7]; // if (idx == 1344)
          //     console.log("dirty")

          srcReg[idx] = dirtyEntries[_i7 + 1];
          srcDist[idx] = dirtyEntries[_i7 + 2];
        }

        if (failed * 3 == stack.length) {
          break;
        }

        if (level > 0) {
          ++iter;

          if (iter >= maxIter) {
            break;
          }
        }
      }

      return srcReg;
    }
  }, {
    key: "sortCellsByLevel",
    value: function sortCellsByLevel(startLevel, chf, srcReg, nbStacks, stacks, loglevelsPerStack) // the levels per stack (2 in our case) as a bit shift
    {
      var w = chf.width;
      var h = chf.height;
      startLevel = startLevel >> loglevelsPerStack;

      for (var j = 0; j < nbStacks; ++j) {
        // stacks[j] = new Array(1024);
        // stacks[j] = [];
        stacks[j] = [];
      }

      ; // put all cells in the level range into the appropriate stacks

      for (var y = 0; y < h; ++y) {
        for (var _x7 = 0; _x7 < w; ++_x7) {
          var c = chf.cells[_x7 + y * w];

          for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
            if (chf.areas[i] == _RecastConstants["default"].RC_NULL_AREA || srcReg[i] != 0) {
              continue;
            }

            var level = chf.dist[i] >> loglevelsPerStack;
            var sId = startLevel - level;

            if (sId >= nbStacks) {
              continue;
            }

            if (sId < 0) {
              sId = 0;
            }

            stacks[sId].push(_x7);
            stacks[sId].push(y);
            stacks[sId].push(i);
          }
        }
      }
    }
  }, {
    key: "appendStacks",
    value: function appendStacks(srcStack, dstStack, srcReg) {
      for (var j = 0; j < srcStack.length; j += 3) {
        var i = srcStack[j + 2];

        if (i < 0 || srcReg[i] != 0) {
          continue;
        }

        dstStack.push(srcStack[j]);
        dstStack.push(srcStack[j + 1]);
        dstStack.push(srcStack[j + 2]);
      }
    }
  }, {
    key: "removeAdjacentNeighbours",
    value: function removeAdjacentNeighbours(reg) {
      // Remove adjacent duplicates.
      for (var i = 0; i < reg.connections.length && reg.connections.length > 1;) {
        var ni = (i + 1) % reg.connections.length;

        if (reg.connections[i] == reg.connections[ni] && reg.connections[i] >= -128 && reg.connections[i] <= 127) {
          reg.connections.splice(i, 1);
        } else {
          ++i;
        }
      }
    }
  }, {
    key: "replaceNeighbour",
    value: function replaceNeighbour(reg, oldId, newId) {
      var neiChanged = false;

      for (var i = 0; i < reg.connections.length; ++i) {
        if (reg.connections[i] == oldId) {
          reg.connections[i] = newId;
          neiChanged = true;
        }
      }

      for (var _i8 = 0; _i8 < reg.floors.length; ++_i8) {
        if (reg.floors[_i8] == oldId) {
          reg.floors.set(_i8, newId);
        }
      }

      if (neiChanged) {
        RecastRegion.removeAdjacentNeighbours(reg);
      }
    }
  }, {
    key: "canMergeWithRegion",
    value: function canMergeWithRegion(rega, regb) {
      if (rega.areaType != regb.areaType) {
        return false;
      }

      var n = 0;

      for (var i = 0; i < rega.connections.length; ++i) {
        if (rega.connections[i] == regb.id) {
          n++;
        }
      }

      if (n > 1) {
        return false;
      }

      for (var _i9 = 0; _i9 < rega.floors.length; ++_i9) {
        if (rega.floors[_i9] == regb.id) {
          return false;
        }
      }

      return true;
    }
  }, {
    key: "addUniqueFloorRegion",
    value: function addUniqueFloorRegion(reg, n) {
      if (!reg.floors.includes(n)) {
        reg.floors.push(n);
      }
    }
  }, {
    key: "mergeRegions",
    value: function mergeRegions(rega, regb) {
      var aid = rega.id;
      var bid = regb.id; // Duplicate current neighbourhood.

      var acon = rega.connections;
      var bcon = regb.connections; // Find insertion poPoly on A.

      var insa = -1;

      for (var i = 0; i < acon.length; ++i) {
        if (acon[i] == bid) {
          insa = i;
          break;
        }
      }

      if (insa == -1) {
        return false;
      } // Find insertion poPoly on B.


      var insb = -1;

      for (var _i10 = 0; _i10 < bcon.length; ++_i10) {
        if (bcon[_i10] == aid) {
          insb = _i10;
          break;
        }
      }

      if (insb == -1) {
        return false;
      } // Merge neighbours.


      rega.connections = [];

      for (var _i11 = 0, ni = acon.length; _i11 < ni - 1; ++_i11) {
        rega.connections.push(acon[(insa + 1 + _i11) % ni]);
      }

      for (var _i12 = 0, _ni3 = bcon.length; _i12 < _ni3 - 1; ++_i12) {
        rega.connections.push(bcon[(insb + 1 + _i12) % _ni3]);
      }

      RecastRegion.removeAdjacentNeighbours(rega);

      for (var j = 0; j < regb.floors.length; ++j) {
        RecastRegion.addUniqueFloorRegion(rega, regb.floors[j]);
      }

      rega.spanCount += regb.spanCount;
      regb.spanCount = 0;
      regb.connections = [];
      return true;
    }
  }, {
    key: "isRegionConnectedToBorder",
    value: function isRegionConnectedToBorder(reg) {
      // Region is connected to border if
      // one of the neighbours is null id.
      return reg.connections.includes(0);
    }
  }, {
    key: "isSolidEdge",
    value: function isSolidEdge(chf, srcReg, x, y, i, dir) {
      var s = chf.spans[i];
      var r = 0;

      if (_RecastCommon["default"].GetCon(s, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
        var ax = x + _RecastCommon["default"].GetDirOffsetX(dir);

        var ay = y + _RecastCommon["default"].GetDirOffsetY(dir);

        var ai = chf.cells[ax + ay * chf.width].index + _RecastCommon["default"].GetCon(s, dir);

        r = srcReg[ai];
      }

      if (r == srcReg[i]) {
        return false;
      }

      return true;
    }
  }, {
    key: "walkContour",
    value: function walkContour(x, y, i, dir, chf, srcReg, cont) {
      var startDir = dir;
      var starti = i;
      var ss = chf.spans[i];
      var curReg = 0;

      if (_RecastCommon["default"].GetCon(ss, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
        var ax = x + _RecastCommon["default"].GetDirOffsetX(dir);

        var ay = y + _RecastCommon["default"].GetDirOffsetY(dir);

        var ai = chf.cells[ax + ay * chf.width].index + _RecastCommon["default"].GetCon(ss, dir);

        curReg = srcReg[ai];
      }

      cont.push(curReg);
      var iter = 0;

      while (++iter < 40000) {
        var s = chf.spans[i];

        if (RecastRegion.isSolidEdge(chf, srcReg, x, y, i, dir)) {
          // Choose the edge corner
          var r = 0;

          if (_RecastCommon["default"].GetCon(s, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
            var _ax6 = x + _RecastCommon["default"].GetDirOffsetX(dir);

            var _ay6 = y + _RecastCommon["default"].GetDirOffsetY(dir);

            var _ai6 = chf.cells[_ax6 + _ay6 * chf.width].index + _RecastCommon["default"].GetCon(s, dir);

            r = srcReg[_ai6];
          }

          if (r != curReg) {
            curReg = r;
            cont.push(curReg);
          }

          dir = dir + 1 & 0x3; // Rotate CW
        } else {
          var ni = -1;

          var nx = x + _RecastCommon["default"].GetDirOffsetX(dir);

          var ny = y + _RecastCommon["default"].GetDirOffsetY(dir);

          if (_RecastCommon["default"].GetCon(s, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
            var nc = chf.cells[nx + ny * chf.width];
            ni = nc.index + _RecastCommon["default"].GetCon(s, dir);
          }

          if (ni == -1) {
            // Should not happen.
            return;
          }

          x = nx;
          y = ny;
          i = ni;
          dir = dir + 3 & 0x3; // Rotate CCW
        }

        if (starti == i && startDir == dir) {
          break;
        }
      } // Remove adjacent duplicates.


      if (cont.length > 1) {
        for (var j = 0; j < cont.length;) {
          var nj = (j + 1) % cont.length;

          if (cont[j] == cont[nj] && cont[j] >= -128 && cont[j] <= 127) {
            cont.splice(j, 1);
          } else {
            ++j;
          }
        }
      }
    }
  }, {
    key: "mergeAndFilterRegions",
    value: function mergeAndFilterRegions(ctx, minRegionArea, mergeRegionSize, maxRegionId, chf, srcReg, overlaps) {
      var w = chf.width;
      var h = chf.height;
      var nreg = maxRegionId + 1;
      var regions = new Array(nreg); // Construct regions

      for (var i = 0; i < nreg; ++i) {
        regions[i] = new RecastRegion.Region(i);
      } // Find edge of a region and find connections around the contour.


      for (var y = 0; y < h; ++y) {
        for (var _x8 = 0; _x8 < w; ++_x8) {
          var c = chf.cells[_x8 + y * w];

          for (var _i13 = c.index, ni = c.index + c.count; _i13 < ni; ++_i13) {
            var r = srcReg[_i13];

            if (r == 0 || r >= nreg) {
              continue;
            }

            var reg = regions[r];
            reg.spanCount++; // Update floors.

            for (var j = c.index; j < ni; ++j) {
              if (_i13 == j) {
                continue;
              }

              var floorId = srcReg[j];

              if (floorId == 0 || floorId >= nreg) {
                continue;
              }

              if (floorId == r) {
                reg.overlap = true;
              }

              RecastRegion.addUniqueFloorRegion(reg, floorId);
            } // Have found contour


            if (reg.connections.length > 0) {
              continue;
            }

            reg.areaType = chf.areas[_i13]; // Check if this cell is next to a border.

            var ndir = -1;

            for (var dir = 0; dir < 4; ++dir) {
              if (RecastRegion.isSolidEdge(chf, srcReg, _x8, y, _i13, dir)) {
                ndir = dir;
                break;
              }
            }

            if (ndir != -1) {
              // The cell is at border.
              // Walk around the contour to find all the neighbours.
              RecastRegion.walkContour(_x8, y, _i13, ndir, chf, srcReg, reg.connections);
            }
          }
        }
      } // Remove too small regions.


      var stack = new Array(32);
      var trace = new Array(32);

      for (var _i14 = 0; _i14 < nreg; ++_i14) {
        var _reg = regions[_i14];

        if (_reg.id == 0 || (_reg.id & _RecastConstants["default"].RC_BORDER_REG) != 0) {
          continue;
        }

        if (_reg.spanCount == 0) {
          continue;
        }

        if (_reg.visited) {
          continue;
        } // Count the total size of all the connected regions.
        // Also keep track of the regions connects to a tile border.


        var connectsToBorder = false;
        var spanCount = 0; //stack = [];
        //trace = [];

        stack = [];
        trace = [];
        _reg.visited = true;
        stack.push(_i14);

        while (stack.length > 0) {
          // Pop
          var ri = stack.splice(stack.length - 1, 1);
          var creg = regions[ri];
          spanCount += creg.spanCount;
          trace.push(ri);

          for (var _j2 = 0; _j2 < creg.connections.length; ++_j2) {
            if ((creg.connections[_j2] & _RecastConstants["default"].RC_BORDER_REG) != 0) {
              connectsToBorder = true;
              continue;
            }

            var neireg = regions[creg.connections[_j2]];

            if (neireg.visited) {
              continue;
            }

            if (neireg.id == 0 || (neireg.id & _RecastConstants["default"].RC_BORDER_REG) != 0) {
              continue;
            } // Visit


            stack.push(neireg.id);
            neireg.visited = true;
          }
        } // If the accumulated regions size is too small, remove it.
        // Do not remove areas which connect to tile borders
        // as their size cannot be estimated correctly and removing them
        // can potentially remove necessary areas.


        if (spanCount < minRegionArea && !connectsToBorder) {
          // Kill all visited regions.
          for (var _j3 = 0; _j3 < trace.length; ++_j3) {
            regions[trace[_j3]].spanCount = 0;
            regions[trace[_j3]].id = 0;
          }
        }
      } // Merge too small regions to neighbour regions.


      var mergeCount = 0;

      do {
        mergeCount = 0;

        for (var _i15 = 0; _i15 < nreg; ++_i15) {
          var _reg2 = regions[_i15];

          if (_reg2.id == 0 || (_reg2.id & _RecastConstants["default"].RC_BORDER_REG) != 0) {
            continue;
          }

          if (_reg2.overlap) {
            continue;
          }

          if (_reg2.spanCount == 0) {
            continue;
          } // Check to see if the region should be merged.


          if (_reg2.spanCount > mergeRegionSize && RecastRegion.isRegionConnectedToBorder(_reg2)) {
            continue;
          } // Small region with more than 1 connection.
          // Or region which is not connected to a border at all.
          // Find smallest neighbour region that connects to this one.


          var smallest = 0xfffffff;
          var mergeId = _reg2.id;

          for (var _j4 = 0; _j4 < _reg2.connections.length; ++_j4) {
            if ((_reg2.connections[_j4] & _RecastConstants["default"].RC_BORDER_REG) != 0) {
              continue;
            }

            var mreg = regions[_reg2.connections[_j4]];

            if (mreg.id == 0 || (mreg.id & _RecastConstants["default"].RC_BORDER_REG) != 0 || mreg.overlap) {
              continue;
            }

            if (mreg.spanCount < smallest && RecastRegion.canMergeWithRegion(_reg2, mreg) && RecastRegion.canMergeWithRegion(mreg, _reg2)) {
              smallest = mreg.spanCount;
              mergeId = mreg.id;
            }
          } // Found new id.


          if (mergeId != _reg2.id) {
            var oldId = _reg2.id;
            var target = regions[mergeId]; // Merge neighbours.

            if (RecastRegion.mergeRegions(target, _reg2)) {
              // Fixup regions pointing to current region.
              for (var _j5 = 0; _j5 < nreg; ++_j5) {
                if (regions[_j5].id == 0 || (regions[_j5].id & _RecastConstants["default"].RC_BORDER_REG) != 0) {
                  continue;
                } // If another region was already merged into current region
                // change the nid of the previous region too.


                if (regions[_j5].id == oldId) {
                  regions[_j5].id = mergeId;
                } // Replace the current region with the new one if the
                // current regions is neighbour.


                RecastRegion.replaceNeighbour(regions[_j5], oldId, mergeId);
              }

              mergeCount++;
            }
          }
        }
      } while (mergeCount > 0); // Compress region Ids.


      for (var _i16 = 0; _i16 < nreg; ++_i16) {
        regions[_i16].remap = false;

        if (regions[_i16].id == 0) {
          continue; // Skip nil regions.
        }

        if ((regions[_i16].id & _RecastConstants["default"].RC_BORDER_REG) != 0) {
          continue; // Skip external regions.
        }

        regions[_i16].remap = true;
      }

      var regIdGen = 0;

      for (var _i17 = 0; _i17 < nreg; ++_i17) {
        if (!regions[_i17].remap) {
          continue;
        }

        var _oldId = regions[_i17].id;
        var newId = ++regIdGen;

        for (var _j6 = _i17; _j6 < nreg; ++_j6) {
          if (regions[_j6].id == _oldId) {
            regions[_j6].id = newId;
            regions[_j6].remap = false;
          }
        }
      }

      maxRegionId = regIdGen; // Remap regions.

      for (var _i18 = 0; _i18 < chf.spanCount; ++_i18) {
        if ((srcReg[_i18] & _RecastConstants["default"].RC_BORDER_REG) == 0) {
          srcReg[_i18] = regions[srcReg[_i18]].id;
        }
      } // Return regions that we found to be overlapping.


      for (var _i19 = 0; _i19 < nreg; ++_i19) {
        if (regions[_i19].overlap) {
          overlaps.push(regions[_i19].id);
        }
      }

      return maxRegionId;
    }
  }, {
    key: "addUniqueConnection",
    value: function addUniqueConnection(reg, n) {
      if (!reg.connections.contains(n)) {
        reg.connections.push(n);
      }
    }
  }, {
    key: "mergeAndFilterLayerRegions",
    value: function mergeAndFilterLayerRegions(ctx, minRegionArea, maxRegionId, chf, srcReg, overlaps) {
      var w = chf.width;
      var h = chf.height;
      var nreg = maxRegionId + 1;
      var regions = new Array(nreg); // Construct regions

      for (var i = 0; i < nreg; ++i) {
        regions[i] = new Region(i);
      } // Find region neighbours and overlapping regions.


      var lregs = new Array(32);

      for (var y = 0; y < h; ++y) {
        for (var _x9 = 0; _x9 < w; ++_x9) {
          var c = chf.cells[_x9 + y * w];
          lregs = [];

          for (var _i20 = c.index, ni = c.index + c.count; _i20 < ni; ++_i20) {
            var s = chf.spans[_i20];
            var ri = srcReg[_i20];

            if (ri == 0 || ri >= nreg) {
              continue;
            }

            var reg = regions[ri];
            reg.spanCount++;
            reg.ymin = Math.min(reg.ymin, s.y);
            reg.ymax = Math.max(reg.ymax, s.y); // Collect all region layers.

            lregs.push(ri); // Update neighbours

            for (var dir = 0; dir < 4; ++dir) {
              if (_RecastCommon["default"].GetCon(s, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                var ax = _x9 + _RecastCommon["default"].GetDirOffsetX(dir);

                var ay = y + _RecastCommon["default"].GetDirOffsetY(dir);

                var ai = chf.cells[ax + ay * w].index + _RecastCommon["default"].GetCon(s, dir);

                var rai = srcReg[ai];

                if (rai > 0 && rai < nreg && rai != ri) {
                  addUniqueConnection(reg, rai);
                }

                if ((rai & _RecastConstants["default"].RC_BORDER_REG) != 0) {
                  reg.connectsToBorder = true;
                }
              }
            }
          } // Update overlapping regions.


          for (var _i21 = 0; _i21 < lregs.length - 1; ++_i21) {
            for (var j = _i21 + 1; j < lregs.length; ++j) {
              if (lregs[_i21] != lregs[j]) {
                var _ri = regions[lregs[_i21]];
                var rj = regions[lregs[j]];
                RecastRegion.addUniqueFloorRegion(_ri, lregs[j]);
                RecastRegion.addUniqueFloorRegion(rj, lregs[_i21]);
              }
            }
          }
        }
      } // Create 2D layers from regions.


      var layerId = 1;

      for (var _i22 = 0; _i22 < nreg; ++_i22) {
        regions[_i22].id = 0;
      } // Merge montone regions to create non-overlapping areas.


      var stack = new Array(32);

      for (var _i23 = 1; _i23 < nreg; ++_i23) {
        var root = regions[_i23]; // Skip already visited.

        if (root.id != 0) {
          continue;
        } // Start search.


        root.id = layerId;
        stack = [];
        stack.push(_i23);

        while (stack.length > 0) {
          // Pop front
          var _reg3 = regions[stack.remove(0)];
          var ncons = _reg3.connections.length;

          for (var _j7 = 0; _j7 < ncons; ++_j7) {
            var nei = _reg3.connections[_j7];
            var regn = regions[nei]; // Skip already visited.

            if (regn.id != 0) {
              continue;
            } // Skip if the neighbour is overlapping root region.


            var overlap = false;

            for (var k = 0; k < root.floors.length; k++) {
              if (root.floors[k] == nei) {
                overlap = true;
                break;
              }
            }

            if (overlap) {
              continue;
            } // Deepen


            stack.push(nei); // Mark layer id

            regn.id = layerId; // Merge current layers to root.

            for (var _k = 0; _k < regn.floors.length; ++_k) {
              RecastRegion.addUniqueFloorRegion(root, regn.floors[_k]);
            }

            root.ymin = Math.min(root.ymin, regn.ymin);
            root.ymax = Math.max(root.ymax, regn.ymax);
            root.spanCount += regn.spanCount;
            regn.spanCount = 0;
            root.connectsToBorder = root.connectsToBorder || regn.connectsToBorder;
          }
        }

        layerId++;
      } // Remove small regions


      for (var _i24 = 0; _i24 < nreg; ++_i24) {
        if (regions[_i24].spanCount > 0 && regions[_i24].spanCount < minRegionArea && !regions[_i24].connectsToBorder) {
          var _reg4 = regions[_i24].id;

          for (var _j8 = 0; _j8 < nreg; ++_j8) {
            if (regions[_j8].id == _reg4) {
              regions[_j8].id = 0;
            }
          }
        }
      } // Compress region Ids.


      for (var _i25 = 0; _i25 < nreg; ++_i25) {
        regions[_i25].remap = false;

        if (regions[_i25].id == 0) {
          continue; // Skip nil regions.
        }

        if ((regions[_i25].id & _RecastConstants["default"].RC_BORDER_REG) != 0) {
          continue; // Skip external regions.
        }

        regions[_i25].remap = true;
      }

      var regIdGen = 0;

      for (var _i26 = 0; _i26 < nreg; ++_i26) {
        if (!regions[_i26].remap) {
          continue;
        }

        var oldId = regions[_i26].id;
        var newId = ++regIdGen;

        for (var _j9 = _i26; _j9 < nreg; ++_j9) {
          if (regions[_j9].id == oldId) {
            regions[_j9].id = newId;
            regions[_j9].remap = false;
          }
        }
      }

      maxRegionId = regIdGen; // Remap regions.

      for (var _i27 = 0; _i27 < chf.spanCount; ++_i27) {
        if ((srcReg[_i27] & _RecastConstants["default"].RC_BORDER_REG) == 0) {
          srcReg[_i27] = regions[srcReg[_i27]].id;
        }
      }

      return maxRegionId;
    } /// @par
    ///
    /// This is usually the second to the last step in creating a fully built
    /// compact heightfield. This step is required before regions are built
    /// using #rcBuildRegions or #rcBuildRegionsMonotone.
    ///
    /// After this step, the distance data is available via the rcCompactHeightfield::maxDistance
    /// and rcCompactHeightfield::dist fields.
    ///
    /// @see rcCompactHeightfield, rcBuildRegions, rcBuildRegionsMonotone

  }, {
    key: "buildDistanceField",
    value: function buildDistanceField(ctx, chf) {
      ctx.startTimer("BUILD_DISTANCEFIELD");
      var src = new Array(chf.spanCount);
      ctx.startTimer("DISTANCEFIELD_DIST");
      var maxDist = this.calculateDistanceField(chf, src);
      chf.maxDistance = maxDist;
      ctx.stopTimer("DISTANCEFIELD_DIST");
      ctx.startTimer("DISTANCEFIELD_BLUR"); // Blur

      src = this.boxBlur(chf, 1, src); // Store distance.

      chf.dist = src;
      ctx.stopTimer("DISTANCEFIELD_BLUR");
      ctx.stopTimer("BUILD_DISTANCEFIELD");
    }
  }, {
    key: "paintRectRegion",
    value: function paintRectRegion(minx, maxx, miny, maxy, regId, chf, srcReg) {
      var w = chf.width;

      for (var y = miny; y < maxy; ++y) {
        for (var _y4 = minx; x < maxx; ++x) {
          var c = chf.cells[x + _y4 * w];

          for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
            if (chf.areas[i] != _RecastConstants["default"].RC_NULL_AREA) {
              srcReg[i] = regId;
            }
          }
        }
      }
    } /// @par
    ///
    /// Non-null regions will consist of connected, non-overlapping walkable spans that form a single contour.
    /// Contours will form simple polygons.
    ///
    /// If multiple regions form an area that is smaller than @p minRegionArea, then all spans will be
    /// re-assigned to the zero (null) region.
    ///
    /// Partitioning can result in smaller than necessary regions. @p mergeRegionArea helps
    /// reduce unecessarily small regions.
    ///
    /// See the #rcConfig documentation for more information on the configuration parameters.
    ///
    /// The region data will be available via the rcCompactHeightfield::maxRegions
    /// and rcCompactSpan::reg fields.
    ///
    /// @warning The distance field must be created using #rcBuildDistanceField before attempting to build regions.
    ///
    /// @see rcCompactHeightfield, rcCompactSpan, rcBuildDistanceField, rcBuildRegionsMonotone, rcConfig

  }, {
    key: "buildRegionsMonotone",
    value: function buildRegionsMonotone(ctx, chf, borderSize, minRegionArea, mergeRegionArea) {
      ctx.startTimer("BUILD_REGIONS");
      var w = chf.width;
      var h = chf.height;
      var id = 1;
      var srcReg = new Array(chf.spanCount);
      var nsweeps = Math.max(chf.width, chf.height);
      sweeps = new Array(nsweeps);

      for (var i = 0; i < sweeps.length; i++) {
        sweeps[i] = new SweepSpan();
      } // Mark border regions.


      if (borderSize > 0) {
        // Make sure border will not overflow.
        var bw = Math.min(w, borderSize);
        var bh = Math.min(h, borderSize); // PaPoly regions

        paintRectRegion(0, bw, 0, h, id | _RecastConstants["default"].RC_BORDER_REG, chf, srcReg);
        id++;
        paintRectRegion(w - bw, w, 0, h, id | _RecastConstants["default"].RC_BORDER_REG, chf, srcReg);
        id++;
        paintRectRegion(0, w, 0, bh, id | _RecastConstants["default"].RC_BORDER_REG, chf, srcReg);
        id++;
        paintRectRegion(0, w, h - bh, h, id | _RecastConstants["default"].RC_BORDER_REG, chf, srcReg);
        id++;
      }

      chf.borderSize = borderSize;
      var prev = new Array(256); // Sweep one line at a time.

      for (var y = borderSize; y < h - borderSize; ++y) {
        // Collect spans from this row.
        prev.fill(0, 0, id);
        var rid = 1;

        for (var _y5 = borderSize; x < w - borderSize; ++x) {
          var c = chf.cells[x + _y5 * w];

          for (var _i28 = c.index, ni = c.index + c.count; _i28 < ni; ++_i28) {
            var s = chf.spans[_i28];

            if (chf.areas[_i28] == _RecastConstants["default"].RC_NULL_AREA) {
              continue;
            } // -x


            var previd = 0;

            if (_RecastCommon["default"].GetCon(s, 0) != _RecastConstants["default"].RC_NOT_CONNECTED) {
              var ax = x + _RecastCommon["default"].GetDirOffsetX(0);

              var ay = _y5 + _RecastCommon["default"].GetDirOffsetY(0);

              var ai = chf.cells[ax + ay * w].index + _RecastCommon["default"].GetCon(s, 0);

              if ((srcReg[ai] & _RecastConstants["default"].RC_BORDER_REG) == 0 && chf.areas[_i28] == chf.areas[ai]) {
                previd = srcReg[ai];
              }
            }

            if (previd == 0) {
              previd = rid++;
              sweeps[previd].rid = previd;
              sweeps[previd].ns = 0;
              sweeps[previd].nei = 0;
            } // -y


            if (_RecastCommon["default"].GetCon(s, 3) != _RecastConstants["default"].RC_NOT_CONNECTED) {
              var _ax7 = x + _RecastCommon["default"].GetDirOffsetX(3);

              var _ay7 = _y5 + _RecastCommon["default"].GetDirOffsetY(3);

              var _ai7 = chf.cells[_ax7 + _ay7 * w].index + _RecastCommon["default"].GetCon(s, 3);

              if (srcReg[_ai7] != 0 && (srcReg[_ai7] & _RecastConstants["default"].RC_BORDER_REG) == 0 && chf.areas[_i28] == chf.areas[_ai7]) {
                var nr = srcReg[_ai7];

                if (sweeps[previd].nei == 0 || sweeps[previd].nei == nr) {
                  sweeps[previd].nei = nr;
                  sweeps[previd].ns++;
                  prev[nr]++;
                } else {
                  sweeps[previd].nei = RC_NULL_NEI;
                }
              }
            }

            srcReg[_i28] = previd;
          }
        } // Create unique ID.


        for (var _i29 = 1; _i29 < rid; ++_i29) {
          if (sweeps[_i29].nei != RC_NULL_NEI && sweeps[_i29].nei != 0 && prev[sweeps[_i29].nei] == sweeps[_i29].ns) {
            sweeps[_i29].id = sweeps[_i29].nei;
          } else {
            sweeps[_i29].id = id++;
          }
        } // Remap IDs


        for (var _y6 = borderSize; x < w - borderSize; ++x) {
          var _c3 = chf.cells[x + _y6 * w];

          for (var _i30 = _c3.index, _ni4 = _c3.index + _c3.count; _i30 < _ni4; ++_i30) {
            if (srcReg[_i30] > 0 && srcReg[_i30] < rid) {
              srcReg[_i30] = sweeps[srcReg[_i30]].id;
            }
          }
        }
      }

      ctx.startTimer("BUILD_REGIONS_FILTER"); // Merge regions and filter out small regions.

      var overlaps = [];
      chf.maxRegions = mergeAndFilterRegions(ctx, minRegionArea, mergeRegionArea, id, chf, srcReg, overlaps); // Monotone partitioning does not generate overlapping regions.

      ctx.stopTimer("BUILD_REGIONS_FILTER"); // Store the result out.

      for (var _i31 = 0; _i31 < chf.spanCount; ++_i31) {
        // if (i == 1344)
        //     console.log("4431")
        chf.spans[_i31].reg = srcReg[_i31];
      }

      ctx.stopTimer("BUILD_REGIONS");
    } /// @par
    ///
    /// Non-null regions will consist of connected, non-overlapping walkable spans that form a single contour.
    /// Contours will form simple polygons.
    ///
    /// If multiple regions form an area that is smaller than @p minRegionArea, then all spans will be
    /// re-assigned to the zero (null) region.
    ///
    /// Watershed partitioning can result in smaller than necessary regions, especially in diagonal corridors.
    /// @p mergeRegionArea helps reduce unecessarily small regions.
    ///
    /// See the #rcConfig documentation for more information on the configuration parameters.
    ///
    /// The region data will be available via the rcCompactHeightfield::maxRegions
    /// and rcCompactSpan::reg fields.
    ///
    /// @warning The distance field must be created using #rcBuildDistanceField before attempting to build regions.
    ///
    /// @see rcCompactHeightfield, rcCompactSpan, rcBuildDistanceField, rcBuildRegionsMonotone, rcConfig

  }, {
    key: "buildRegions",
    value: function buildRegions(ctx, chf, borderSize, minRegionArea, mergeRegionArea) {
      ctx.startTimer("BUILD_REGIONS");
      var w = chf.width;
      var h = chf.height;
      ctx.startTimer("REGIONS_WATERSHED");
      var LOG_NB_STACKS = 3;
      var NB_STACKS = 1 << LOG_NB_STACKS;
      var lvlStacks = [];

      for (var i = 0; i < NB_STACKS; ++i) {
        lvlStacks.push([]);
      }

      var stack = new Array(1024);
      var srcReg = new Array(chf.spanCount).fill(0);
      var srcDist = new Array(chf.spanCount).fill(0);
      var regionId = 1;
      var level = chf.maxDistance + 1 & ~1; // TODO: Figure better formula, expandIters defines how much the
      // watershed "overflows" and simplifies the regions. Tying it to
      // agent radius was usually good indication how greedy it could be.
      // const let expandIters = 4 + walkableRadius * 2;

      var expandIters = 8;

      if (borderSize > 0) {
        // Make sure border will not overflow.
        var bw = Math.min(w, borderSize);
        var bh = Math.min(h, borderSize); // PaPoly regions

        paintRectRegion(0, bw, 0, h, regionId | _RecastConstants["default"].RC_BORDER_REG, chf, srcReg);
        regionId++;
        paintRectRegion(w - bw, w, 0, h, regionId | _RecastConstants["default"].RC_BORDER_REG, chf, srcReg);
        regionId++;
        paintRectRegion(0, w, 0, bh, regionId | _RecastConstants["default"].RC_BORDER_REG, chf, srcReg);
        regionId++;
        paintRectRegion(0, w, h - bh, h, regionId | _RecastConstants["default"].RC_BORDER_REG, chf, srcReg);
        regionId++;
      }

      chf.borderSize = borderSize;
      var sId = -1;

      while (level > 0) {
        level = level >= 2 ? level - 2 : 0;
        sId = sId + 1 & NB_STACKS - 1; // ctx=>startTimer(RC_TIMER_DIVIDE_TO_LEVELS);

        if (sId == 0) {
          RecastRegion.sortCellsByLevel(level, chf, srcReg, NB_STACKS, lvlStacks, 1);
        } else {
          RecastRegion.appendStacks(lvlStacks[sId - 1], lvlStacks[sId], srcReg); // copy left overs from last level
        } // ctx=>stopTimer(RC_TIMER_DIVIDE_TO_LEVELS);


        ctx.startTimer("BUILD_REGIONS_EXPAND"); // Expand current regions until no empty connected cells found.

        RecastRegion.expandRegions(expandIters, level, chf, srcReg, srcDist, lvlStacks[sId], false);
        ctx.stopTimer("BUILD_REGIONS_EXPAND");
        ctx.startTimer("BUILD_REGIONS_FLOOD"); // Mark new regions with IDs.

        for (var j = 0; j < lvlStacks[sId].length; j += 3) {
          var _x10 = lvlStacks[sId][j];
          var y = lvlStacks[sId][j + 1];
          var _i32 = lvlStacks[sId][j + 2];

          if (_i32 >= 0 && srcReg[_i32] == 0) {
            if (RecastRegion.floodRegion(_x10, y, _i32, level, regionId, chf, srcReg, srcDist, stack)) {
              regionId++;
            }
          }
        }

        ctx.stopTimer("BUILD_REGIONS_FLOOD");
      } // Expand current regions until no empty connected cells found.


      RecastRegion.expandRegions(expandIters * 8, 0, chf, srcReg, srcDist, stack, true);
      ctx.stopTimer("BUILD_REGIONS_WATERSHED");
      ctx.startTimer("BUILD_REGIONS_FILTER"); // Merge regions and filter out smalle regions.

      var overlaps = [];
      chf.maxRegions = this.mergeAndFilterRegions(ctx, minRegionArea, mergeRegionArea, regionId, chf, srcReg, overlaps); // If overlapping regions were found during merging, split those regions.

      if (overlaps.length > 0) {
        throw new RuntimeException("rcBuildRegions: " + overlaps.length + " overlapping regions.");
      }

      ctx.stopTimer("BUILD_REGIONS_FILTER"); // Write the result out.

      for (var _i33 = 0; _i33 < chf.spanCount; ++_i33) {
        // if (i == 1344)
        //     console.log("4431")
        chf.spans[_i33].reg = srcReg[_i33];
      }

      ctx.stopTimer("BUILD_REGIONS");
    }
  }, {
    key: "buildLayerRegions",
    value: function buildLayerRegions(ctx, chf, borderSize, minRegionArea) {
      ctx.startTimer("BUILD_REGIONS");
      var w = chf.width;
      var h = chf.height;
      var id = 1;
      var srcReg = new Array(chf.spanCount);
      var nsweeps = Math.max(chf.width, chf.height);
      var sweeps = new ArrayList(nsweeps);

      for (var i = 0; i < sweeps.length; i++) {
        sweeps[i] = new SweepSpan();
      } // Mark border regions.


      if (borderSize > 0) {
        // Make sure border will not overflow.
        var bw = Math.min(w, borderSize);
        var bh = Math.min(h, borderSize); // PaPoly regions

        paintRectRegion(0, bw, 0, h, id | _RecastConstants["default"].RC_BORDER_REG, chf, srcReg);
        id++;
        paintRectRegion(w - bw, w, 0, h, id | _RecastConstants["default"].RC_BORDER_REG, chf, srcReg);
        id++;
        paintRectRegion(0, w, 0, bh, id | _RecastConstants["default"].RC_BORDER_REG, chf, srcReg);
        id++;
        paintRectRegion(0, w, h - bh, h, id | _RecastConstants["default"].RC_BORDER_REG, chf, srcReg);
        id++;
      }

      chf.borderSize = borderSize;
      var prev = new Array(256); // Sweep one line at a time.

      for (var y = borderSize; y < h - borderSize; ++y) {
        // Collect spans from this row.
        prev.fill(0, 0, id);
        var rid = 1;

        for (var _y7 = borderSize; x < w - borderSize; ++x) {
          var c = chf.cells[x + _y7 * w];

          for (var _i34 = c.index, ni = c.index + c.count; _i34 < ni; ++_i34) {
            var s = chf.spans[_i34];

            if (chf.areas[_i34] == _RecastConstants["default"].RC_NULL_AREA) {
              continue;
            } // -x


            var previd = 0;

            if (_RecastCommon["default"].GetCon(s, 0) != _RecastConstants["default"].RC_NOT_CONNECTED) {
              var ax = x + _RecastCommon["default"].GetDirOffsetX(0);

              var ay = _y7 + _RecastCommon["default"].GetDirOffsetY(0);

              var ai = chf.cells[ax + ay * w].index + _RecastCommon["default"].GetCon(s, 0);

              if ((srcReg[ai] & _RecastConstants["default"].RC_BORDER_REG) == 0 && chf.areas[_i34] == chf.areas[ai]) {
                previd = srcReg[ai];
              }
            }

            if (previd == 0) {
              previd = rid++;
              sweeps[previd].rid = previd;
              sweeps[previd].ns = 0;
              sweeps[previd].nei = 0;
            } // -y


            if (_RecastCommon["default"].GetCon(s, 3) != _RecastConstants["default"].RC_NOT_CONNECTED) {
              var _ax8 = x + _RecastCommon["default"].GetDirOffsetX(3);

              var _ay8 = _y7 + _RecastCommon["default"].GetDirOffsetY(3);

              var _ai8 = chf.cells[_ax8 + _ay8 * w].index + _RecastCommon["default"].GetCon(s, 3);

              if (srcReg[_ai8] != 0 && (srcReg[_ai8] & _RecastConstants["default"].RC_BORDER_REG) == 0 && chf.areas[_i34] == chf.areas[_ai8]) {
                var nr = srcReg[_ai8];

                if (sweeps[previd].nei == 0 || sweeps[previd].nei == nr) {
                  sweeps[previd].nei = nr;
                  sweeps[previd].ns++;
                  prev[nr]++;
                } else {
                  sweeps[previd].nei = RC_NULL_NEI;
                }
              }
            }

            srcReg[_i34] = previd;
          }
        } // Create unique ID.


        for (var _i35 = 1; _i35 < rid; ++_i35) {
          if (sweeps[_i35].nei != RC_NULL_NEI && sweeps[_i35].nei != 0 && prev[sweeps[_i35].nei] == sweeps[_i35].ns) {
            sweeps[_i35].id = sweeps[_i35].nei;
          } else {
            sweeps[_i35].id = id++;
          }
        } // Remap IDs


        for (var _y8 = borderSize; x < w - borderSize; ++x) {
          var _c4 = chf.cells[x + _y8 * w];

          for (var _i36 = _c4.index, _ni5 = _c4.index + _c4.count; _i36 < _ni5; ++_i36) {
            if (srcReg[_i36] > 0 && srcReg[_i36] < rid) {
              srcReg[_i36] = sweeps[srcReg[_i36]].id;
            }
          }
        }
      }

      ctx.startTimer("BUILD_REGIONS_FILTER"); // Merge monotone regions to layers and remove small regions.

      var overlaps = [];
      chf.maxRegions = mergeAndFilterLayerRegions(ctx, minRegionArea, id, chf, srcReg, overlaps);
      ctx.stopTimer("BUILD_REGIONS_FILTER"); // Store the result out.

      for (var _i37 = 0; _i37 < chf.spanCount; ++_i37) {
        chf.spans[_i37].reg = srcReg[_i37];
      }

      ctx.stopTimer("BUILD_REGIONS");
    }
  }]);

  return RecastRegion;
}();

_defineProperty(RecastRegion, "RC_NULL_NEI", 0xfff);

_defineProperty(RecastRegion, "SweepSpan", (_temp = function SweepSpan() {
  _classCallCheck(this, SweepSpan);

  _defineProperty(this, "rid", void 0);

  _defineProperty(this, "id", void 0);

  _defineProperty(this, "ns", void 0);

  _defineProperty(this, "nei", void 0);
} // neighbour id
, _temp));

_defineProperty(RecastRegion, "Region", (_temp2 = // Number of spans belonging to this region
// ID of the region
// Are type.
function Region(i) {
  _classCallCheck(this, Region);

  _defineProperty(this, "spanCount", 0);

  _defineProperty(this, "id", 0);

  _defineProperty(this, "areaType", 0);

  _defineProperty(this, "remap", false);

  _defineProperty(this, "visited", false);

  _defineProperty(this, "overlap", false);

  _defineProperty(this, "connectsToBorder", false);

  _defineProperty(this, "ymin", 0);

  _defineProperty(this, "ymax", 0);

  _defineProperty(this, "connections", []);

  _defineProperty(this, "floors", []);

  this.id = i;
  this.ymin = 0xFFFF;
  this.connections = [];
  this.floors = [];
}, _temp2));

var _default = RecastRegion;
exports["default"] = _default;