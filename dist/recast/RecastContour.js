"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _RecastVectors = _interopRequireDefault(require("./RecastVectors.js"));

var _RecastConstants = _interopRequireDefault(require("./RecastConstants.js"));

var _RecastCommon = _interopRequireDefault(require("./RecastCommon.js"));

var _ContourSet = _interopRequireDefault(require("./ContourSet.js"));

var _RecastMesh = _interopRequireDefault(require("./RecastMesh.js"));

var _Contour = _interopRequireDefault(require("./Contour.js"));

var _temp, _temp2, _temp3;

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var RecastContour = /*#__PURE__*/function () {
  function RecastContour() {
    _classCallCheck(this, RecastContour);

    _defineProperty(this, "CompareHoles", function (a, b) {
      if (a.minx == b.minx) {
        if (a.minz < b.minz) return -1;
        if (a.minz > b.minz) return 1;
      } else {
        if (a.minx < b.minx) return -1;
        if (a.minx > b.minx) return 1;
      }

      return 0;
    });

    _defineProperty(this, "CompareDiagDist", function (va, vb) {
      var a = va;
      var b = vb;
      if (a.dist < b.dist) return -1;
      if (a.dist > b.dist) return 1;
      return 0;
    });
  }

  _createClass(RecastContour, [{
    key: "mergeRegionHoles",
    value: function mergeRegionHoles(ctx, region) {
      // Sort holes from left to right.
      for (var i = 0; i < region.nholes; i++) {
        var minleft = findLeftMostVertex(region.holes[i].contour);
        region.holes[i].minx = minleft[0];
        region.holes[i].minz = minleft[1];
        region.holes[i].leftmost = minleft[2];
      }

      Arrays.sort(region.holes, new CompareHoles());
      var maxVerts = region.outline.nverts;

      for (var _i = 0; _i < region.nholes; _i++) {
        maxVerts += region.holes[_i].contour.nverts;
      }

      diags = new Array(maxVerts);

      for (var pd = 0; pd < maxVerts; pd++) {
        diags[pd] = new PotentialDiagonal();
      }

      var outline = region.outline; // Merge holes into the outline one by one.

      for (var _i2 = 0; _i2 < region.nholes; _i2++) {
        var hole = region.holes[_i2].contour;
        var index = -1;
        var bestVertex = region.holes[_i2].leftmost;

        for (var iter = 0; iter < hole.nverts; iter++) {
          // Find potential diagonals.
          // The 'best' vertex must be in the cone described by 3 cosequtive vertices of the outline.
          // ..o j-1
          //   |
          //   |   * best
          //   |
          // j o-----o j+1
          //         :
          var ndiags = 0;
          var corner = bestVertex * 4;

          for (var j = 0; j < outline.nverts; j++) {
            if (inCone(j, outline.nverts, outline.verts, corner, hole.verts)) {
              var dx = outline.verts[j * 4 + 0] - hole.verts[corner + 0];
              var dz = outline.verts[j * 4 + 2] - hole.verts[corner + 2];
              diags[ndiags].vert = j;
              diags[ndiags].dist = dx * dx + dz * dz;
              ndiags++;
            }
          } // Sort potential diagonals by distance, we want to make the connection as short as possible.


          Arrays.sort(diags, 0, ndiags, new CompareDiagDist()); // Find a diagonal that is not intersecting the outline not the remaining holes.

          index = -1;

          for (var _j = 0; _j < ndiags; _j++) {
            var pt = diags[_j].vert * 4;
            var intersect = intersectSegCountour(pt, corner, diags[_i2].vert, outline.nverts, outline.verts, outline.verts, hole.verts);

            for (var k = _i2; k < region.nholes && !intersect; k++) {
              intersect |= intersectSegCountour(pt, corner, -1, region.holes[k].contour.nverts, region.holes[k].contour.verts, outline.verts, hole.verts);
            }

            if (!intersect) {
              index = diags[_j].vert;
              break;
            }
          } // If found non-intersecting diagonal, stop looking.


          if (index != -1) break; // All the potential diagonals for the current vertex were intersecting, try next vertex.

          bestVertex = (bestVertex + 1) % hole.nverts;
        }

        if (index == -1) {
          ctx.warn("mergeHoles: Failed to find merge points for");
          continue;
        }

        mergeContours(region.outline, hole, index, bestVertex);
      }
    } /// @par
    ///
    /// The raw contours will match the region outlines exactly. The @p maxError and @p maxEdgeLen
    /// parameters control how closely the simplified contours will match the raw contours.
    ///
    /// Simplified contours are generated such that the vertices for portals between areas match up.
    /// (They are considered mandatory vertices.)
    ///
    /// Setting @p maxEdgeLength to zero will disabled the edge length feature.
    ///
    /// See the #rcConfig documentation for more information on the configuration parameters.
    ///
    /// @see rcAllocContourSet, rcCompactHeightfield, rcContourSet, rcConfig

  }], [{
    key: "getCornerHeight",
    value: function getCornerHeight(x, y, i, dir, chf, isBorderVertex) {
      var s = chf.spans[i];
      var ch = s.y;
      var dirp = dir + 1 & 0x3;
      var regs = [0, 0, 0, 0]; // Combine region and area codes in order to prevent
      // border vertices which are in between two areas to be removed.

      regs[0] = chf.spans[i].reg | chf.areas[i] << 16;

      if (_RecastCommon["default"].GetCon(s, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
        var ax = x + _RecastCommon["default"].GetDirOffsetX(dir);

        var ay = y + _RecastCommon["default"].GetDirOffsetY(dir);

        var ai = chf.cells[ax + ay * chf.width].index + _RecastCommon["default"].GetCon(s, dir);

        var as = chf.spans[ai];
        ch = Math.max(ch, as.y);
        regs[1] = chf.spans[ai].reg | chf.areas[ai] << 16;

        if (_RecastCommon["default"].GetCon(as, dirp) != _RecastConstants["default"].RC_NOT_CONNECTED) {
          var ax2 = ax + _RecastCommon["default"].GetDirOffsetX(dirp);

          var ay2 = ay + _RecastCommon["default"].GetDirOffsetY(dirp);

          var ai2 = chf.cells[ax2 + ay2 * chf.width].index + _RecastCommon["default"].GetCon(as, dirp);

          var as2 = chf.spans[ai2];
          ch = Math.max(ch, as2.y);
          regs[2] = chf.spans[ai2].reg | chf.areas[ai2] << 16;
        }
      }

      if (_RecastCommon["default"].GetCon(s, dirp) != _RecastConstants["default"].RC_NOT_CONNECTED) {
        var _ax = x + _RecastCommon["default"].GetDirOffsetX(dirp);

        var _ay = y + _RecastCommon["default"].GetDirOffsetY(dirp);

        var _ai = chf.cells[_ax + _ay * chf.width].index + _RecastCommon["default"].GetCon(s, dirp);

        var _as = chf.spans[_ai];
        ch = Math.max(ch, _as.y);
        regs[3] = chf.spans[_ai].reg | chf.areas[_ai] << 16;

        if (_RecastCommon["default"].GetCon(_as, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
          var _ax2 = _ax + _RecastCommon["default"].GetDirOffsetX(dir);

          var _ay2 = _ay + _RecastCommon["default"].GetDirOffsetY(dir);

          var _ai2 = chf.cells[_ax2 + _ay2 * chf.width].index + _RecastCommon["default"].GetCon(_as, dir);

          var _as2 = chf.spans[_ai2];
          ch = Math.max(ch, _as2.y);
          regs[2] = chf.spans[_ai2].reg | chf.areas[_ai2] << 16;
        }
      } // Check if the vertex is special edge vertex, these vertices will be removed later.


      for (var j = 0; j < 4; ++j) {
        var a = j;
        var b = j + 1 & 0x3;
        var c = j + 2 & 0x3;
        var d = j + 3 & 0x3; // The vertex is a border vertex there are two same exterior cells in a row,
        // followed by two interior cells and none of the regions are out of bounds.

        var twoSameExts = (regs[a] & regs[b] & _RecastConstants["default"].RC_BORDER_REG) != 0 && regs[a] == regs[b];
        var twoInts = ((regs[c] | regs[d]) & _RecastConstants["default"].RC_BORDER_REG) == 0;
        var intsSameArea = regs[c] >> 16 == regs[d] >> 16;
        var noZeros = regs[a] != 0 && regs[b] != 0 && regs[c] != 0 && regs[d] != 0;

        if (twoSameExts && twoInts && intsSameArea && noZeros) {
          isBorderVertex = true;
          break;
        }
      }

      return ch;
    }
  }, {
    key: "walkContour",
    value: function walkContour(x, y, i, chf, flags, points) {
      // Choose the first non-connected edge
      var dir = 0;

      while ((flags[i] & 1 << dir) == 0) {
        dir++;
      }

      var startDir = dir;
      var starti = i;
      var area = chf.areas[i];
      var iter = 0;

      while (++iter < 40000) {
        if ((flags[i] & 1 << dir) != 0) {
          // Choose the edge corner
          var isBorderVertex = false;
          var isAreaBorder = false;
          var px = x;
          var py = RecastContour.getCornerHeight(x, y, i, dir, chf, isBorderVertex);
          var pz = y;

          switch (dir) {
            case 0:
              pz++;
              break;

            case 1:
              px++;
              pz++;
              break;

            case 2:
              px++;
              break;
          }

          var r = 0;
          var s = chf.spans[i];

          if (_RecastCommon["default"].GetCon(s, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
            var ax = x + _RecastCommon["default"].GetDirOffsetX(dir);

            var ay = y + _RecastCommon["default"].GetDirOffsetY(dir);

            var ai = chf.cells[ax + ay * chf.width].index + _RecastCommon["default"].GetCon(s, dir);

            r = chf.spans[ai].reg;
            if (area != chf.areas[ai]) isAreaBorder = true;
          }

          if (isBorderVertex) r |= _RecastConstants["default"].RC_BORDER_VERTEX;
          if (isAreaBorder) r |= _RecastConstants["default"].RC_AREA_BORDER;
          points.push(px);
          points.push(py);
          points.push(pz);
          points.push(r);
          flags[i] &= ~(1 << dir); // Remove visited edges

          dir = dir + 1 & 0x3; // Rotate CW
        } else {
          var ni = -1;

          var nx = x + _RecastCommon["default"].GetDirOffsetX(dir);

          var ny = y + _RecastCommon["default"].GetDirOffsetY(dir);

          var _s = chf.spans[i];

          if (_RecastCommon["default"].GetCon(_s, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
            var nc = chf.cells[nx + ny * chf.width];
            ni = nc.index + _RecastCommon["default"].GetCon(_s, dir);
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
      }
    }
  }, {
    key: "distancePtSeg",
    value: function distancePtSeg(x, z, px, pz, qx, qz) {
      var pqx = qx - px;
      var pqz = qz - pz;
      var dx = x - px;
      var dz = z - pz;
      var d = pqx * pqx + pqz * pqz;
      var t = pqx * dx + pqz * dz;
      if (d > 0) t /= d;
      if (t < 0) t = 0;else if (t > 1) t = 1;
      dx = px + t * pqx - x;
      dz = pz + t * pqz - z;
      return dx * dx + dz * dz;
    }
  }, {
    key: "simplifyContour",
    value: function simplifyContour(points, simplified, maxError, maxEdgeLen, buildFlags) {
      // Add initial points.
      var hasConnections = false;

      for (var i = 0; i < points.length; i += 4) {
        if ((points[i + 3] & _RecastConstants["default"].RC_CONTOUR_REG_MASK) != 0) {
          hasConnections = true;
          break;
        }
      }

      if (hasConnections) {
        // The contour has some portals to other regions.
        // Add a new poPoly to every location where the region changes.
        for (var _i3 = 0, ni = points.length / 4; _i3 < ni; ++_i3) {
          var ii = (_i3 + 1) % ni;
          var differentRegs = (points[_i3 * 4 + 3] & _RecastConstants["default"].RC_CONTOUR_REG_MASK) != (points[ii * 4 + 3] & _RecastConstants["default"].RC_CONTOUR_REG_MASK);
          var areaBorders = (points[_i3 * 4 + 3] & _RecastConstants["default"].RC_AREA_BORDER) != (points[ii * 4 + 3] & _RecastConstants["default"].RC_AREA_BORDER);

          if (differentRegs || areaBorders) {
            simplified.push(points[_i3 * 4 + 0]);
            simplified.push(points[_i3 * 4 + 1]);
            simplified.push(points[_i3 * 4 + 2]);
            simplified.push(_i3);
          }
        }
      }

      if (simplified.length == 0) {
        // If there is no connections at all,
        // create some initial points for the simplification process.
        // Find lower-left and upper-right vertices of the contour.
        var llx = points[0];
        var lly = points[1];
        var llz = points[2];
        var lli = 0;
        var urx = points[0];
        var ury = points[1];
        var urz = points[2];
        var uri = 0;

        for (var _i4 = 0; _i4 < points.length; _i4 += 4) {
          var x = points[_i4 + 0];
          var y = points[_i4 + 1];
          var z = points[_i4 + 2];

          if (x < llx || x == llx && z < llz) {
            llx = x;
            lly = y;
            llz = z;
            lli = _i4 / 4;
          }

          if (x > urx || x == urx && z > urz) {
            urx = x;
            ury = y;
            urz = z;
            uri = _i4 / 4;
          }
        }

        simplified.push(llx);
        simplified.push(lly);
        simplified.push(llz);
        simplified.push(lli);
        simplified.push(urx);
        simplified.push(ury);
        simplified.push(urz);
        simplified.push(uri);
      } // Add points until all raw points are within
      // error tolerance to the simplified shape.


      var pn = points.length / 4;

      for (var _i5 = 0; _i5 < simplified.length / 4;) {
        // console.log(simplified.length)
        var _ii = (_i5 + 1) % (simplified.length / 4);

        var ax = simplified[_i5 * 4 + 0];
        var az = simplified[_i5 * 4 + 2];
        var ai = simplified[_i5 * 4 + 3];
        var bx = simplified[_ii * 4 + 0];
        var bz = simplified[_ii * 4 + 2];
        var bi = simplified[_ii * 4 + 3]; // Find maximum deviation from the segment.

        var maxd = 0;
        var maxi = -1;
        var ci = void 0,
            cinc = void 0,
            endi = void 0; // Traverse the segment in lexilogical order so that the
        // max deviation is calculated similarly when traversing
        // opposite segments.

        if (bx > ax || bx == ax && bz > az) {
          cinc = 1;
          ci = (ai + cinc) % pn;
          endi = bi;
        } else {
          cinc = pn - 1;
          ci = (bi + cinc) % pn;
          endi = ai;
          var temp = ax;
          ax = bx;
          bx = temp;
          temp = az;
          az = bz;
          bz = temp;
        } // Tessellate only outer edges or edges between areas.
        // console.log("start")


        if ((points[ci * 4 + 3] & _RecastConstants["default"].RC_CONTOUR_REG_MASK) == 0 || (points[ci * 4 + 3] & _RecastConstants["default"].RC_AREA_BORDER) != 0) {
          while (ci != endi) {
            // if(Math.random() < .01) console.log(`${ci} ${endi}`);
            var d = RecastContour.distancePtSeg(points[ci * 4 + 0], points[ci * 4 + 2], ax, az, bx, bz);

            if (d > maxd) {
              maxd = d;
              maxi = ci;
            }

            ci = (ci + cinc) % pn;
          }
        } // console.log("stop")
        // If the max deviation is larger than accepted error,
        // add new point, else continue to next segment.


        if (maxi != -1 && maxd > maxError * maxError) {
          // Add the point.
          simplified.splice((_i5 + 1) * 4 + 0, 0, points[maxi * 4 + 0]);
          simplified.splice((_i5 + 1) * 4 + 1, 0, points[maxi * 4 + 1]);
          simplified.splice((_i5 + 1) * 4 + 2, 0, points[maxi * 4 + 2]);
          simplified.splice((_i5 + 1) * 4 + 3, 0, maxi);
        } else {
          ++_i5;
        }
      } // Split too let edges.


      if (maxEdgeLen > 0 && (buildFlags & (_RecastConstants["default"].RC_CONTOUR_TESS_WALL_EDGES | _RecastConstants["default"].RC_CONTOUR_TESS_AREA_EDGES)) != 0) {
        for (var _i6 = 0; _i6 < simplified.length / 4;) {
          var _ii2 = (_i6 + 1) % (simplified.length / 4);

          var _ax3 = simplified[_i6 * 4 + 0];
          var _az = simplified[_i6 * 4 + 2];
          var _ai3 = simplified[_i6 * 4 + 3];
          var _bx = simplified[_ii2 * 4 + 0];
          var _bz = simplified[_ii2 * 4 + 2];
          var _bi = simplified[_ii2 * 4 + 3]; // Find maximum deviation from the segment.

          var _maxi = -1;

          var _ci = (_ai3 + 1) % pn; // Tessellate only outer edges or edges between areas.


          var tess = false; // Wall edges.

          if ((buildFlags & _RecastConstants["default"].RC_CONTOUR_TESS_WALL_EDGES) != 0 && (points[_ci * 4 + 3] & _RecastConstants["default"].RC_CONTOUR_REG_MASK) == 0) tess = true; // Edges between areas.

          if ((buildFlags & _RecastConstants["default"].RC_CONTOUR_TESS_AREA_EDGES) != 0 && (points[_ci * 4 + 3] & _RecastConstants["default"].RC_AREA_BORDER) != 0) tess = true;

          if (tess) {
            var dx = _bx - _ax3;
            var dz = _bz - _az;

            if (dx * dx + dz * dz > maxEdgeLen * maxEdgeLen) {
              // Round based on the segments in lexilogical order so that the
              // max tesselation is consistent regardles in which direction
              // segments are traversed.
              var n = _bi < _ai3 ? _bi + pn - _ai3 : _bi - _ai3;

              if (n > 1) {
                if (_bx > _ax3 || _bx == _ax3 && _bz > _az) _maxi = Math.floor((_ai3 + n / 2) % pn);else _maxi = Math.floor((_ai3 + (n + 1) / 2) % pn);
              }
            }
          } // If the max deviation is larger than accepted error,
          // add new point, else continue to next segment.


          if (_maxi != -1) {
            // Add the point.
            simplified.splice((_i6 + 1) * 4 + 0, 0, points[_maxi * 4 + 0]);
            simplified.splice((_i6 + 1) * 4 + 1, 0, points[_maxi * 4 + 1]);
            simplified.splice((_i6 + 1) * 4 + 2, 0, points[_maxi * 4 + 2]);
            simplified.splice((_i6 + 1) * 4 + 3, 0, _maxi);
          } else {
            ++_i6;
          }
        }
      }

      for (var _i7 = 0; _i7 < simplified.length / 4; ++_i7) {
        // The edge vertex flag is take from the current raw point,
        // and the neighbour region is take from the next raw point.
        var _ai4 = (simplified[_i7 * 4 + 3] + 1) % pn;

        var _bi2 = simplified[_i7 * 4 + 3];
        simplified[_i7 * 4 + 3] = points[_ai4 * 4 + 3] & (_RecastConstants["default"].RC_CONTOUR_REG_MASK | _RecastConstants["default"].RC_AREA_BORDER) | points[_bi2 * 4 + 3] & _RecastConstants["default"].RC_BORDER_VERTEX;
      }
    }
  }, {
    key: "calcAreaOfPolygon2D",
    value: function calcAreaOfPolygon2D(verts, nverts) {
      var area = 0;

      for (var i = 0, j = nverts - 1; i < nverts; j = i++) {
        var vi = i * 4;
        var vj = j * 4;
        area += verts[vi + 0] * verts[vj + 2] - verts[vj + 0] * verts[vi + 2];
      }

      return (area + 1) / 2;
    }
  }, {
    key: "intersectSegCountour",
    value: function intersectSegCountour(d0, d1, i, n, verts, d0verts, d1verts) {
      // For each edge (k,k+1) of P
      var pverts = new Array(4 * 4);

      for (var g = 0; g < 4; g++) {
        pverts[g] = d0verts[d0 + g];
        pverts[4 + g] = d1verts[d1 + g];
      }

      d0 = 0;
      d1 = 4;

      for (var k = 0; k < n; k++) {
        var k1 = _RecastMesh["default"].next(k, n); // Skip edges incident to i.


        if (i == k || i == k1) continue;
        var p0 = k * 4;
        var p1 = k1 * 4;

        for (var _g = 0; _g < 4; _g++) {
          pverts[8 + _g] = verts[p0 + _g];
          pverts[12 + _g] = verts[p1 + _g];
        }

        p0 = 8;
        p1 = 12;
        if (_RecastMesh["default"].vequal(pverts, d0, p0) || _RecastMesh["default"].vequal(pverts, d1, p0) || _RecastMesh["default"].vequal(pverts, d0, p1) || _RecastMesh["default"].vequal(pverts, d1, p1)) continue;
        if (_RecastMesh["default"].intersect(pverts, d0, d1, p0, p1)) return true;
      }

      return false;
    }
  }, {
    key: "inCone",
    value: function inCone(i, n, verts, pj, vertpj) {
      pi = i * 4;
      pi1 = _RecastMesh["default"].next(i, n) * 4;
      pin1 = _RecastMesh["default"].prev(i, n) * 4;
      pverts = new Array(4 * 4);

      for (var g = 0; g < 4; g++) {
        pverts[g] = verts[pi + g];
        pverts[4 + g] = verts[pi1 + g];
        pverts[8 + g] = verts[pin1 + g];
        pverts[12 + g] = vertpj[pj + g];
      }

      pi = 0;
      pi1 = 4;
      pin1 = 8;
      pj = 12; // If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].

      if (_RecastMesh["default"].leftOn(pverts, pin1, pi, pi1)) return _RecastMesh["default"].left(pverts, pi, pj, pin1) && _RecastMesh["default"].left(pverts, pj, pi, pi1); // Assume (i-1,i,i+1) not collinear.
      // else P[i] is reflex.

      return !(_RecastMesh["default"].leftOn(pverts, pi, pj, pi1) && _RecastMesh["default"].leftOn(pverts, pj, pi, pin1));
    }
  }, {
    key: "removeDegenerateSegments",
    value: function removeDegenerateSegments(simplified) {
      // Remove adjacent vertices which are equal on xz-plane,
      // or else the triangulator will get confused.
      var npts = simplified.length / 4;

      for (var i = 0; i < npts; ++i) {
        var ni = _RecastMesh["default"].next(i, npts); //			if (vequal(&simplified[i*4], &simplified[ni*4]))


        if (simplified[i * 4] == simplified[ni * 4] && simplified[i * 4 + 2] == simplified[ni * 4 + 2]) {
          // Degenerate segment, remove.
          simplified.splice(i * 4, 1);
          simplified.splice(i * 4, 1);
          simplified.splice(i * 4, 1);
          simplified.splice(i * 4, 1);
          npts--;
        }
      }
    }
  }, {
    key: "mergeContours",
    value: function mergeContours(ca, cb, ia, ib) {
      var maxVerts = ca.nverts + cb.nverts + 2;
      var verts = new Array(maxVerts * 4);
      var nv = 0; // Copy contour A.

      for (var i = 0; i <= ca.nverts; ++i) {
        var dst = nv * 4;
        var src = (ia + i) % ca.nverts * 4;
        verts[dst + 0] = ca.verts[src + 0];
        verts[dst + 1] = ca.verts[src + 1];
        verts[dst + 2] = ca.verts[src + 2];
        verts[dst + 3] = ca.verts[src + 3];
        nv++;
      } // Copy contour B


      for (var _i8 = 0; _i8 <= cb.nverts; ++_i8) {
        var _dst = nv * 4;

        var _src = (ib + _i8) % cb.nverts * 4;

        verts[_dst + 0] = cb.verts[_src + 0];
        verts[_dst + 1] = cb.verts[_src + 1];
        verts[_dst + 2] = cb.verts[_src + 2];
        verts[_dst + 3] = cb.verts[_src + 3];
        nv++;
      }

      ca.verts = verts;
      ca.nverts = nv;
      cb.verts = null;
      cb.nverts = 0;
    } // Finds the lowest leftmost vertex of a contour.

  }, {
    key: "findLeftMostVertex",
    value: function findLeftMostVertex(contour) {
      var minx = contour.verts[0];
      var minz = contour.verts[2];
      var leftmost = 0;

      for (var i = 1; i < contour.nverts; i++) {
        var x = contour.verts[i * 4 + 0];
        var z = contour.verts[i * 4 + 2];

        if (x < minx || x == minx && z < minz) {
          minx = x;
          minz = z;
          leftmost = i;
        }
      }

      return [minx, minz, leftmost];
    }
  }, {
    key: "buildContours",
    value: function buildContours(ctx, chf, maxError, maxEdgeLen, buildFlags) {
      var w = chf.width;
      var h = chf.height;
      var borderSize = chf.borderSize;
      var cset = new _ContourSet["default"]();
      ctx.startTimer("BUILD_CONTOURS");

      _RecastVectors["default"].copy3(cset.bmin, chf.bmin, 0);

      _RecastVectors["default"].copy3(cset.bmax, chf.bmax, 0);

      if (borderSize > 0) {
        // If the heightfield was build with bordersize, remove the offset.
        pad = borderSize * chf.cs;
        cset.bmin[0] += pad;
        cset.bmin[2] += pad;
        cset.bmax[0] -= pad;
        cset.bmax[2] -= pad;
      }

      cset.cs = chf.cs;
      cset.ch = chf.ch;
      cset.width = chf.width - chf.borderSize * 2;
      cset.height = chf.height - chf.borderSize * 2;
      cset.borderSize = chf.borderSize;
      cset.maxError = maxError;
      var flags = new Array(chf.spanCount).fill(0);
      ctx.startTimer("BUILD_CONTOURS_TRACE"); // Mark boundaries.

      for (var y = 0; y < h; ++y) {
        for (var x = 0; x < w; ++x) {
          var c = chf.cells[x + y * w];

          for (var i = c.index, ni = c.index + c.count; i < ni; ++i) {
            // if(y == 3 && x == 3 && i == 1344)
            // 	console.log("earlier");
            var res = 0;
            var s = chf.spans[i];

            if (chf.spans[i].reg == 0 || (chf.spans[i].reg & _RecastConstants["default"].RC_BORDER_REG) != 0) {
              flags[i] = 0;
              continue;
            }

            for (var dir = 0; dir < 4; ++dir) {
              var r = 0;

              if (_RecastCommon["default"].GetCon(s, dir) != _RecastConstants["default"].RC_NOT_CONNECTED) {
                var ax = x + _RecastCommon["default"].GetDirOffsetX(dir);

                var ay = y + _RecastCommon["default"].GetDirOffsetY(dir);

                var ai = chf.cells[ax + ay * w].index + _RecastCommon["default"].GetCon(s, dir);

                r = chf.spans[ai].reg;
              }

              if (r == chf.spans[i].reg) res |= 1 << dir;
            }

            flags[i] = res ^ 0xf; // Inverse, mark non connected edges.
          }
        }
      }

      ctx.stopTimer("BUILD_CONTOURS_TRACE");
      var verts = new Array(256);
      var simplified = new Array(64);
      ;

      for (var _y = 0; _y < h; ++_y) {
        for (var _x = 0; _x < w; ++_x) {
          var _c = chf.cells[_x + _y * w];

          for (var _i9 = _c.index, _ni = _c.index + _c.count; _i9 < _ni; ++_i9) {
            // if(y==3 && x == 3 && i == 1344)
            // 	console.log("dak");
            if (flags[_i9] == 0 || flags[_i9] == 0xf) {
              flags[_i9] = 0;
              continue;
            }

            var reg = chf.spans[_i9].reg;
            if (reg == 0 || (reg & _RecastConstants["default"].RC_BORDER_REG) != 0) continue;
            var area = chf.areas[_i9];
            verts = [];
            simplified = [];
            ctx.startTimer("BUILD_CONTOURS_TRACE");
            RecastContour.walkContour(_x, _y, _i9, chf, flags, verts);
            ctx.stopTimer("BUILD_CONTOURS_TRACE");
            ctx.startTimer("BUILD_CONTOURS_SIMPLIFY");
            RecastContour.simplifyContour(verts, simplified, maxError, maxEdgeLen, buildFlags);
            RecastContour.removeDegenerateSegments(simplified);
            ctx.stopTimer("BUILD_CONTOURS_SIMPLIFY"); // Store region=>contour remap info.
            // Create contour.

            if (simplified.length / 4 >= 3) {
              var cont = new _Contour["default"]();
              cset.conts.push(cont);
              cont.nverts = simplified.length / 4;
              cont.verts = new Array(simplified.length);

              for (var l = 0; l < cont.verts.length; l++) {
                cont.verts[l] = simplified[l];
              }

              if (borderSize > 0) {
                // If the heightfield was build with bordersize, remove the offset.
                for (var j = 0; j < cont.nverts; ++j) {
                  cont.verts[j * 4] -= borderSize;
                  cont.verts[j * 4 + 2] -= borderSize;
                }
              }

              cont.nrverts = verts.length / 4;
              cont.rverts = new Array(verts.length);

              for (var _l = 0; _l < cont.rverts.length; _l++) {
                cont.rverts[_l] = verts[_l];
              }

              if (borderSize > 0) {
                // If the heightfield was build with bordersize, remove the offset.
                for (var _j2 = 0; _j2 < cont.nrverts; ++_j2) {
                  cont.rverts[_j2 * 4] -= borderSize;
                  cont.rverts[_j2 * 4 + 2] -= borderSize;
                }
              }

              cont.reg = reg;
              cont.area = area;
            }
          }
        }
      } // Merge holes if needed.


      if (cset.conts.length > 0) {
        // Calculate winding of all polygons.
        var winding = new Array(cset.conts.length);
        var nholes = 0;

        for (var _i10 = 0; _i10 < cset.conts.length; ++_i10) {
          var _cont = cset.conts[_i10]; // If the contour is wound backwards, it is a hole.

          winding[_i10] = RecastContour.calcAreaOfPolygon2D(_cont.verts, _cont.nverts) < 0 ? -1 : 1;
          if (winding[_i10] < 0) nholes++;
        }

        if (nholes > 0) {
          // Collect outline contour and holes contours per region.
          // We assume that there is one outline and multiple holes.
          var nregions = chf.maxRegions + 1;
          var regions = new Array(nregions);

          for (var _i11 = 0; _i11 < nregions; _i11++) {
            regions[_i11] = new ContourRegion();
          }

          for (var _i12 = 0; _i12 < cset.conts.length; ++_i12) {
            var _cont2 = cset.conts[_i12]; // Positively would contours are outlines, negative holes.

            if (winding[_i12] > 0) {
              if (regions[_cont2.reg].outline != null) {
                throw new RuntimeException("rcBuildContours: Multiple outlines for region " + _cont2.reg + ".");
              }

              regions[_cont2.reg].outline = _cont2;
            } else {
              regions[_cont2.reg].nholes++;
            }
          }

          for (var _i13 = 0; _i13 < nregions; _i13++) {
            if (regions[_i13].nholes > 0) {
              regions[_i13].holes = new ContourHole[regions[_i13].nholes]();

              for (var nh = 0; nh < regions[_i13].nholes; nh++) {
                regions[_i13].holes[nh] = new ContourHole();
              }

              regions[_i13].nholes = 0;
            }
          }

          for (var _i14 = 0; _i14 < cset.conts.length; ++_i14) {
            var _cont3 = cset.conts[_i14];
            var _reg = regions[_cont3.reg];
            if (winding[_i14] < 0) _reg.holes[_reg.nholes++].contour = _cont3;
          } // Finally merge each regions holes into the outline.


          for (var _i15 = 0; _i15 < nregions; _i15++) {
            var _reg2 = regions[_i15];
            if (_reg2.nholes == 0) continue;

            if (_reg2.outline != null) {
              mergeRegionHoles(ctx, _reg2);
            } else {
              // The region does not have an outline.
              // This can happen if the contour becaomes selfoverlapping because of
              // too aggressive simplification settings.
              throw new RuntimeException("rcBuildContours: Bad outline for region " + _i15 + ", contour simplification is likely too aggressive.");
            }
          }
        }
      }

      ctx.stopTimer("BUILD_CONTOURS");
      return cset;
    }
  }]);

  return RecastContour;
}();

_defineProperty(RecastContour, "ContourRegion", (_temp = function ContourRegion() {
  _classCallCheck(this, ContourRegion);

  _defineProperty(this, "outline", void 0);

  _defineProperty(this, "holes", []);

  _defineProperty(this, "nholes", void 0);
}, _temp));

_defineProperty(RecastContour, "ContourHole", (_temp2 = function ContourHole() {
  _classCallCheck(this, ContourHole);

  _defineProperty(this, "leftmost", void 0);

  _defineProperty(this, "minx", void 0);

  _defineProperty(this, "minz", void 0);

  _defineProperty(this, "contour", void 0);
}, _temp2));

_defineProperty(RecastContour, "PotentialDiagonal", (_temp3 = function PotentialDiagonal() {
  _classCallCheck(this, PotentialDiagonal);

  _defineProperty(this, "dist", void 0);

  _defineProperty(this, "vert", void 0);
}, _temp3));

var _default = RecastContour;
exports["default"] = _default;