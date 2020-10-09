"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _NavMeshParams = _interopRequireDefault(require("./NavMeshParams.js"));

var _DetourCommon = _interopRequireDefault(require("./DetourCommon.js"));

var _MeshTile = _interopRequireDefault(require("./MeshTile.js"));

var _Poly = _interopRequireDefault(require("./Poly.js"));

var _Link = _interopRequireDefault(require("./Link.js"));

var _ClosestPointOnPolyResult = _interopRequireDefault(require("./ClosestPointOnPolyResult.js"));

var _FindNearestPolyResult = _interopRequireDefault(require("./FindNearestPolyResult.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _createForOfIteratorHelper(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

function arraycopy(one, oneStart, two, twoStart, len) {
  for (var i = 0; i < len; i++) {
    two[twoStart + i] = one[oneStart + i];
  }
}

var NavMesh = /*#__PURE__*/function () {
  _createClass(NavMesh, [{
    key: "getMaxTiles",
    /// A flag that indicates that an entity links to an external entity.
    /// (E.g. A polygon edge is a portal that links to another polygon.)
    /// A value that indicates the entity does not link to anything.
    // static DT_NULL_LINK = 0xffffffff;
    /// A flag that indicates that an off-mesh connection can be traversed in
    /// both directions. (Is bidirectional.)
    /// The maximum number of user defined area ids.
    /// Limit raycasting during any angle pahfinding
    /// The limit is given as a multiple of the character radius
    /// < Current initialization params. TODO: do not store this info twice.
    /// < Origin of the tile (0,0)
    //  m_orig[3]; ///< Origin of the tile (0,0)
    /// < Dimensions of each tile.
    /// < Max number of tiles.
    /// < Tile hash lookup size (must be pot).
    /// < Tile hash lookup mask.
    /// < Tile hash lookup.
    /// < Freelist of tiles.
    /// < List of tiles.

    /** The maximum number of vertices per navigation polygon. */

    /**
     * The maximum number of tiles supported by the navigation mesh.
     * 
     * @return The maximum number of tiles supported by the navigation mesh.
     */
    value: function getMaxTiles() {
      return this.m_maxTiles;
    }
    /**
     * Returns tile in the tile array.
     */

  }, {
    key: "getTile",
    value: function getTile(i) {
      return this.m_tiles[i];
    }
    /**
     * Gets the polygon reference for the tile's base polygon.
     * 
     * @param tile
     *            The tile.
     * @return The polygon reference for the base polygon in the specified tile.
     */

  }, {
    key: "getPolyRefBase",
    value: function getPolyRefBase(tile) {
      if (tile == null) return 0;
      var it = tile.index;
      return NavMesh.encodePolyId(tile.salt, it, 0);
    }
    /**
     * Derives a standard polygon reference.
     * 
     * @note This function is generally meant for internal use only.
     * @param salt
     *            The tile's salt value.
     * @param it
     *            The index of the tile.
     * @param ip
     *            The index of the polygon within the tile.
     * @return encoded polygon reference
     */
    //https://stackoverflow.com/a/337572/10047920

  }, {
    key: "allocLink",
    value: function allocLink(tile) {
      if (tile.linksFreeList == NavMesh.DT_NULL_LINK) {
        var _link = new _Link["default"]();

        _link.next = NavMesh.DT_NULL_LINK;
        tile.links.push(_link);
        return tile.links.length - 1;
      }

      var link = tile.linksFreeList;
      tile.linksFreeList = tile.links[link].next;
      return link;
    }
  }, {
    key: "freeLink",
    value: function freeLink(tile, link) {
      tile.links[link].next = tile.linksFreeList;
      tile.linksFreeList = link;
    }
    /**
     * Calculates the tile grid location for the specified world position.
     * 
     * @param pos
     *            The world position for the query. [(x, y, z)]
     * @return 2-element let array with (tx,ty) tile location
     */

  }, {
    key: "calcTileLoc",
    value: function calcTileLoc(pos) {
      var tx = Math.floor((pos[0] - this.m_orig[0]) / this.m_tileWidth);
      var ty = Math.floor((pos[2] - this.m_orig[2]) / this.m_tileHeight);
      return [tx, ty];
    }
  }, {
    key: "getTileAndPolyByRef",
    value: function getTileAndPolyByRef(ref) {
      if (ref == 0) {
        throw new IllegalArgumentException("ref = 0");
      }

      var saltitip = NavMesh.decodePolyId(ref);
      var salt = saltitip[0];
      var it = saltitip[1];
      var ip = saltitip[2];
      if (it >= this.m_maxTiles) throw new IllegalArgumentException("tile > m_maxTiles");
      if (this.m_tiles[it].salt != salt || this.m_tiles[it].data.header == null) throw new IllegalArgumentException("Invalid salt or header");
      if (ip >= this.m_tiles[it].data.header.polyCount) throw new IllegalArgumentException("poly > polyCount");
      return [this.m_tiles[it], this.m_tiles[it].data.polys[ip]];
    } /// @par
    ///
    /// @warning Only use this function if it is known that the provided polygon
    /// reference is valid. This function is faster than #getTileAndPolyByRef,
    /// but
    /// it does not validate the reference.

  }, {
    key: "getTileAndPolyByRefUnsafe",
    value: function getTileAndPolyByRefUnsafe(ref) {
      var saltitip = NavMesh.decodePolyId(ref);
      var it = saltitip[1];
      var ip = saltitip[2];
      return [this.m_tiles[it], this.m_tiles[it].data.polys[ip]];
    }
  }, {
    key: "isValidPolyRef",
    value: function isValidPolyRef(ref) {
      if (ref == 0) return false;
      var saltitip = NavMesh.decodePolyId(ref);
      var salt = saltitip[0];
      var it = saltitip[1];
      var ip = saltitip[2];
      if (it >= this.m_maxTiles) return false;
      if (this.m_tiles[it].salt != salt || this.m_tiles[it].data == null) return false;
      if (ip >= this.m_tiles[it].data.header.polyCount) return false;
      return true;
    }
  }, {
    key: "getParams",
    value: function getParams() {
      return this.m_params;
    }
  }], [{
    key: "lshift",
    value: function lshift(num, bits) {
      return num * Math.pow(2, bits);
    }
  }, {
    key: "rshift",
    value: function rshift(num, bits) {
      return Math.floor(num / Math.pow(2, bits));
    } //https://stackoverflow.com/a/43666199/10047920

  }, {
    key: "and",
    value: function and(v1, v2) {
      var hi = 0x80000000;
      var low = 0x7fffffff;
      var hi1 = ~~(v1 / hi);
      var hi2 = ~~(v2 / hi);
      var low1 = v1 & low;
      var low2 = v2 & low;
      var h = hi1 & hi2;
      var l = low1 & low2;
      return h * hi + l;
    }
  }, {
    key: "or",
    value: function or(v1, v2) {
      var hi = 0x80000000;
      var low = 0x7fffffff;
      var hi1 = ~~(v1 / hi);
      var hi2 = ~~(v2 / hi);
      var low1 = v1 & low;
      var low2 = v2 & low;
      var h = hi1 | hi2;
      var l = low1 | low2;
      return h * hi + l;
    }
  }, {
    key: "encodePolyId",
    value: function encodePolyId(salt, it, ip) {
      var a = NavMesh.lshift(salt, NavMesh.DT_POLY_BITS + NavMesh.DT_TILE_BITS);
      var b = NavMesh.lshift(it, NavMesh.DT_POLY_BITS);
      var c = ip;
      return NavMesh.or(NavMesh.or(a, b), ip);
    } /// Decodes a standard polygon reference.
    /// @note This function is generally meant for internal use only.
    /// @param[in] ref The polygon reference to decode.
    /// @param[out] salt The tile's salt value.
    /// @param[out] it The index of the tile.
    /// @param[out] ip The index of the polygon within the tile.
    /// @see #encodePolyId

  }, {
    key: "decodePolyId",
    value: function decodePolyId(ref) {
      var salt;
      var it;
      var ip;
      var saltMask = NavMesh.lshift(1, NavMesh.DT_SALT_BITS) - 1;
      var tileMask = NavMesh.lshift(1, NavMesh.DT_TILE_BITS) - 1;
      var polyMask = NavMesh.lshift(1, NavMesh.DT_POLY_BITS) - 1;
      salt = Math.floor(NavMesh.and(NavMesh.rshift(ref, NavMesh.DT_POLY_BITS + NavMesh.DT_TILE_BITS), saltMask));
      it = Math.floor(NavMesh.and(NavMesh.rshift(ref, NavMesh.DT_POLY_BITS), tileMask));
      ip = Math.floor(NavMesh.and(ref, polyMask));
      return [salt, it, ip];
    } /// Extracts a tile's salt value from the specified polygon reference.
    /// @note This function is generally meant for internal use only.
    /// @param[in] ref The polygon reference.
    /// @see #encodePolyId

  }, {
    key: "decodePolyIdSalt",
    value: function decodePolyIdSalt(ref) {
      var saltMask = (1 << NavMesh.DT_SALT_BITS) - 1;
      return Math.floor(ref >> NavMesh.DT_POLY_BITS + NavMesh.DT_TILE_BITS & saltMask);
    } /// Extracts the tile's index from the specified polygon reference.
    /// @note This function is generally meant for internal use only.
    /// @param[in] ref The polygon reference.
    /// @see #encodePolyId

  }, {
    key: "decodePolyIdTile",
    value: function decodePolyIdTile(ref) {
      var tileMask = (1 << NavMesh.DT_TILE_BITS) - 1;
      return Math.floor(ref >> NavMesh.DT_POLY_BITS & tileMask);
    } /// Extracts the polygon's index (within its tile) from the specified
    /// polygon reference.
    /// @note This function is generally meant for internal use only.
    /// @param[in] ref The polygon reference.
    /// @see #encodePolyId

  }, {
    key: "decodePolyIdPoly",
    value: function decodePolyIdPoly(ref) {
      var polyMask = (1 << NavMesh.DT_POLY_BITS) - 1;
      return Math.floor(ref & polyMask);
    }
  }]);

  function NavMesh(one, maxVertsPerPoly, flags) {
    _classCallCheck(this, NavMesh);

    _defineProperty(this, "m_params", null);

    _defineProperty(this, "m_orig", []);

    _defineProperty(this, "m_tileWidth", 0);

    _defineProperty(this, "m_tileHeight", 0);

    _defineProperty(this, "m_maxTiles", 0);

    _defineProperty(this, "m_tileLutSize", 0);

    _defineProperty(this, "m_tileLutMask", 0);

    _defineProperty(this, "m_posLookup", []);

    _defineProperty(this, "m_nextFree", null);

    _defineProperty(this, "m_tiles", []);

    _defineProperty(this, "m_maxVertPerPoly", 0);

    _defineProperty(this, "m_tileCount", 0);

    if (flags || flags == 0) {
      this._constructor(NavMesh.getNavMeshParams(one), maxVertsPerPoly);

      this.addTile(one, flags, 0);
    } else {
      this._constructor(one, maxVertsPerPoly);
    }
  }

  _createClass(NavMesh, [{
    key: "_constructor",
    value: function _constructor(params, maxVertsPerPoly) {
      this.m_params = params;
      this.m_orig = params.orig;
      this.m_tileWidth = params.tileWidth;
      this.m_tileHeight = params.tileHeight; // Init tiles

      this.m_maxTiles = params.maxTiles;
      this.m_maxVertPerPoly = maxVertsPerPoly;

      var lutsize = _DetourCommon["default"].nextPow2(params.maxTiles / 4);

      if (lutsize == 0) lutsize = 1;
      this.m_tileLutSize = lutsize;
      this.m_tileLutMask = this.m_tileLutSize - 1;
      this.m_tiles = new Array(this.m_maxTiles);
      this.m_posLookup = new Array(this.m_tileLutSize);
      this.m_nextFree = null;

      for (var i = this.m_maxTiles - 1; i >= 0; --i) {
        this.m_tiles[i] = new _MeshTile["default"](i);
        this.m_tiles[i].salt = 1;
        this.m_tiles[i].next = this.m_nextFree;
        this.m_nextFree = this.m_tiles[i];
      }
    }
  }, {
    key: "queryPolygonsInTile",
    // TODO: These methods are duplicates from dtNavMeshQuery, but are needed
    // for off-mesh connection finding.
    value: function queryPolygonsInTile(tile, qmin, qmax) {
      var polys = [];

      if (tile.data.bvTree != null) {
        var nodeIndex = 0;
        var tbmin = tile.data.header.bmin;
        var tbmax = tile.data.header.bmax;
        var qfac = tile.data.header.bvQuantFactor; // Calculate quantized box

        var _bmin = new Array(3);

        var _bmax = new Array(3); // dtClamp query box to world box.


        var _minx = _DetourCommon["default"].clamp(qmin[0], tbmin[0], tbmax[0]) - tbmin[0];

        var miny = _DetourCommon["default"].clamp(qmin[1], tbmin[1], tbmax[1]) - tbmin[1];
        var minz = _DetourCommon["default"].clamp(qmin[2], tbmin[2], tbmax[2]) - tbmin[2];

        var _maxx = _DetourCommon["default"].clamp(qmax[0], tbmin[0], tbmax[0]) - tbmin[0];

        var maxy = _DetourCommon["default"].clamp(qmax[1], tbmin[1], tbmax[1]) - tbmin[1];
        var maxz = _DetourCommon["default"].clamp(qmax[2], tbmin[2], tbmax[2]) - tbmin[2]; // Quantize

        _bmin[0] = NavMesh.and(Math.floor(qfac * _minx), 0xfffe);
        _bmin[1] = NavMesh.and(Math.floor(qfac * miny), 0xfffe);
        _bmin[2] = NavMesh.and(Math.floor(qfac * minz), 0xfffe);
        _bmax[0] = NavMesh.or(Math.floor(qfac * _maxx + 1), 1);
        _bmax[1] = NavMesh.or(Math.floor(qfac * maxy + 1), 1);
        _bmax[2] = NavMesh.or(Math.floor(qfac * maxz + 1), 1); // Traverse tree

        var base = this.getPolyRefBase(tile);
        var end = tile.data.header.bvNodeCount;

        while (nodeIndex < end) {
          var node = tile.data.bvTree[nodeIndex];

          var overlap = _DetourCommon["default"].overlapQuantBounds(_bmin, _bmax, node.bmin, node.bmax);

          var isLeafNode = node.i >= 0;

          if (isLeafNode && overlap) {
            polys.push(NavMesh.or(base, node.i));
          }

          if (overlap || isLeafNode) nodeIndex++;else {
            var escapeIndex = -node.i;
            nodeIndex += escapeIndex;
          }
        }

        return polys;
      } else {
        bmin = [null, null, null];
        bmax = [null, null, null];

        var _base = this.getPolyRefBase(tile);

        for (var i = 0; i < tile.data.header.polyCount; ++i) {
          var p = tile.data.polys[i]; // Do not return off-mesh connection polygons.

          if (p.getType() == _Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION) continue; // Calc polygon bounds.

          var v = p.verts[0] * 3;

          _DetourCommon["default"].vCopy(bmin, tile.data.verts, v);

          _DetourCommon["default"].vCopy(bmax, tile.data.verts, v);

          for (var j = 1; j < p.vertCount; ++j) {
            v = p.verts[j] * 3;

            _DetourCommon["default"].vMin(bmin, tile.data.verts, v);

            _DetourCommon["default"].vMax(bmax, tile.data.verts, v);
          }

          if (overlapBounds(qmin, qmax, bmin, bmax)) {
            polys.push(NavMesh.or(_base, i));
          }
        }

        return polys;
      }
    } /// Adds a tile to the navigation mesh.
    /// @param[in] data Data for the new tile mesh. (See: #dtCreateNavMeshData)
    /// @param[in] dataSize Data size of the new tile mesh.
    /// @param[in] flags Tile flags. (See: #dtTileFlags)
    /// @param[in] lastRef The desired reference for the tile. (When reloading a
    /// tile.) [opt] [Default: 0]
    /// @param[out] result The tile reference. (If the tile was succesfully
    /// added.) [opt]
    /// @return The status flags for the operation.
    /// @par
    ///
    /// The add operation will fail if the data is in the wrong format, the
    /// allocated tile
    /// space is full, or there is a tile already at the specified reference.
    ///
    /// The lastRef parameter is used to restore a tile with the same tile
    /// reference it had previously used. In this case the #dtPolyRef's for the
    /// tile will be restored to the same values they were before the tile was
    /// removed.
    ///
    /// The nav mesh assumes exclusive access to the data passed and will make
    /// changes to the dynamic portion of the data. For that reason the data
    /// should not be reused in other nav meshes until the tile has been successfully
    /// removed from this nav mesh.
    ///
    /// @see dtCreateNavMeshData, #removeTile

  }, {
    key: "addTile",
    value: function addTile(data, flags, lastRef) {
      // Make sure the data is in right format.
      var header = data.header; // Make sure the location is free.

      if (this.getTileAt(header.x, header.y, header.layer) != null) throw new RuntimeException("Tile already exists"); // Allocate a tile.

      var tile = null;

      if (lastRef == 0) {
        if (this.m_nextFree != null) {
          tile = this.m_nextFree;
          this.m_nextFree = tile.next;
          tile.next = null;
          this.m_tileCount++;
        }
      } else {
        // Try to relocate the tile to specific index with same salt.
        var tileIndex = decodePolyIdTile(lastRef);
        if (tileIndex >= this.m_maxTiles) throw new RuntimeException("Tile index too high"); // Try to find the specific tile id from the free list.

        var target = this.m_tiles[tileIndex];
        var prev = null;
        tile = this.m_nextFree;

        while (tile != null && tile != target) {
          prev = tile;
          tile = tile.next;
        } // Could not find the correct location.


        if (tile != target) throw new RuntimeException("Could not find tile"); // Remove from freelist

        if (prev == null) this.m_nextFree = tile.next;else prev.next = tile.next; // Restore salt.

        tile.salt = decodePolyIdSalt(lastRef);
      } // Make sure we could allocate a tile.


      if (tile == null) throw new RuntimeException("Could not allocate a tile");
      tile.data = data;
      tile.flags = flags;
      tile.links = []; // Insert tile into the position lut.

      var h = NavMesh.computeTileHash(header.x, header.y, this.m_tileLutMask);
      tile.next = this.m_posLookup[h];
      this.m_posLookup[h] = tile; // Patch header pointers.
      // If there are no items in the bvtree, reset the tree pointer.

      if (tile.data.bvTree != null && tile.data.bvTree.length == 0) tile.data.bvTree = null; // Init tile.

      this.connectIntLinks(tile); // Base off-mesh connections to their starting polygons and connect connections inside the tile.

      this.baseOffMeshLinks(tile);
      this.connectExtOffMeshLinks(tile, tile, -1); // Connect with layers in current tile.

      var neis = this.getTilesAt(header.x, header.y);

      for (var j = 0; j < neis.length; ++j) {
        if (neis[j] == tile) {
          continue;
        }

        connectExtLinks(tile, neis[j], -1);
        connectExtLinks(neis[j], tile, -1);
        this.connectExtOffMeshLinks(tile, neis[j], -1);
        this.connectExtOffMeshLinks(neis[j], tile, -1);
      } // Connect with neighbour tiles.


      for (var i = 0; i < 8; ++i) {
        neis = this.getNeighbourTilesAt(header.x, header.y, i);

        for (var _j = 0; _j < neis.length; ++_j) {
          connectExtLinks(tile, neis[_j], i);
          connectExtLinks(neis[_j], tile, _DetourCommon["default"].oppositeTile(i));
          this.connectExtOffMeshLinks(tile, neis[_j], i);
          this.connectExtOffMeshLinks(neis[_j], tile, _DetourCommon["default"].oppositeTile(i));
        }
      }

      return this.getTileRef(tile);
    } /// Removes the specified tile from the navigation mesh.
    /// @param[in] ref The reference of the tile to remove.
    /// @param[out] data Data associated with deleted tile.
    /// @param[out] dataSize Size of the data associated with deleted tile.
    /// @return The status flags for the operation.
    // dtStatus removeTile(dtTileRef ref, char** data, int* dataSize);
    /// @par
    ///
    /// This function returns the data for the tile so that, if desired,
    /// it can be added back to the navigation mesh at a later point.
    ///
    /// @see #addTile

  }, {
    key: "removeTile",
    value: function removeTile(ref) {
      if (ref == 0) {
        return null;
      }

      var tileIndex = decodePolyIdTile(ref);
      var tileSalt = decodePolyIdSalt(ref);
      if (tileIndex >= this.m_maxTiles) throw new RuntimeException("Invalid tile index");
      var tile = this.m_tiles[tileIndex];
      if (tile.salt != tileSalt) throw new RuntimeException("Invalid tile salt"); // Remove tile from hash lookup.

      var h = NavMesh.computeTileHash(tile.data.header.x, tile.data.header.y, this.m_tileLutMask);
      var prev = null;
      var cur = this.m_posLookup[h];

      while (cur != null) {
        if (cur == tile) {
          if (prev != null) prev.next = cur.next;else this.m_posLookup[h] = cur.next;
          break;
        }

        prev = cur;
        cur = cur.next;
      } // Remove connections to neighbour tiles.
      // Create connections with neighbour tiles.
      // Disconnect from other layers in current tile.


      nneis = this.getTilesAt(tile.data.header.x, tile.data.header.y);

      var _iterator = _createForOfIteratorHelper(nneis),
          _step;

      try {
        for (_iterator.s(); !(_step = _iterator.n()).done;) {
          var _j2 = _step.value;
          if (_j2 == tile) continue;
          unconnectLinks(_j2, tile);
        } // Disconnect from neighbour tiles.

      } catch (err) {
        _iterator.e(err);
      } finally {
        _iterator.f();
      }

      for (var i = 0; i < 8; ++i) {
        nneis = this.getNeighbourTilesAt(tile.data.header.x, tile.data.header.y, i);

        var _iterator2 = _createForOfIteratorHelper(nneis),
            _step2;

        try {
          for (_iterator2.s(); !(_step2 = _iterator2.n()).done;) {
            var j = _step2.value;
            unconnectLinks(j, tile);
          }
        } catch (err) {
          _iterator2.e(err);
        } finally {
          _iterator2.f();
        }
      }

      var data = tile.data; // Reset tile.

      tile.data = null;
      tile.flags = 0;
      tile.links = []; // Update salt, salt should never be zero.

      tile.salt = tile.salt + 1 & (1 << NavMesh.DT_SALT_BITS) - 1;
      if (tile.salt == 0) tile.salt++; // Add to free list.

      tile.next = this.m_nextFree;
      this.m_nextFree = tile;
      this.m_tileCount--;
      return data;
    } /// Builds internal polygons links for a tile.

  }, {
    key: "connectIntLinks",
    value: function connectIntLinks(tile) {
      if (tile == null) return;
      var base = this.getPolyRefBase(tile);

      for (var i = 0; i < tile.data.header.polyCount; ++i) {
        var poly = tile.data.polys[i];
        poly.firstLink = NavMesh.DT_NULL_LINK;
        if (poly.getType() == _Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION) continue; // Build edge links backwards so that the links will be
        // in the linked list from lowest index to highest.

        for (var j = poly.vertCount - 1; j >= 0; --j) {
          // Skip hard and non-internal edges.
          if (poly.neis[j] == 0 || (poly.neis[j] & NavMesh.DT_EXT_LINK) != 0) continue;
          var idx = this.allocLink(tile);
          var link = tile.links[idx];
          link.ref = NavMesh.or(base, poly.neis[j] - 1);
          link.edge = j;
          link.side = 0xff;
          link.bmin = link.bmax = 0; // Add to linked list.

          link.next = poly.firstLink;
          poly.firstLink = idx;
        }
      }
    }
  }, {
    key: "unconnectLinks",
    value: function unconnectLinks(tile, target) {
      if (tile == null || target == null) return;
      var targetNum = decodePolyIdTile(this.getTileRef(target));

      for (var i = 0; i < tile.data.header.polyCount; ++i) {
        var poly = tile.data.polys[i];
        var j = poly.firstLink;
        var pj = NavMesh.DT_NULL_LINK;

        while (j != NavMesh.DT_NULL_LINK) {
          if (decodePolyIdTile(tile.links[j].ref) == targetNum) {
            // Remove link.
            var nj = tile.links[j].next;
            if (pj == NavMesh.DT_NULL_LINK) poly.firstLink = nj;else tile.links[pj].next = nj;
            freeLink(tile, j);
            j = nj;
          } else {
            // Advance
            pj = j;
            j = tile.links[j].next;
          }
        }
      }
    }
  }, {
    key: "connectExtLinks",
    value: function connectExtLinks(tile, target, side) {
      if (tile == null) return; // Connect border links.

      for (var i = 0; i < tile.data.header.polyCount; ++i) {
        var poly = tile.data.polys[i]; // Create new links.
        // short m = NavMesh.DT_EXT_LINK | (short)side;

        var nv = poly.vertCount;

        for (var j = 0; j < nv; ++j) {
          // Skip non-portal edges.
          if ((poly.neis[j] & NavMesh.DT_EXT_LINK) == 0) continue;
          var dir = poly.neis[j] & 0xff;
          if (side != -1 && dir != side) continue; // Create new links

          var va = poly.verts[j] * 3;
          var vb = poly.verts[(j + 1) % nv] * 3;
          connectedPolys = findConnectingPolys(tile.data.verts, va, vb, target, _DetourCommon["default"].oppositeTile(dir), 4);
          nei = connectedPolys[0];
          neia = connectedPolys[1];
          var nnei = connectedPolys.third;

          for (var k = 0; k < nnei; ++k) {
            var idx = this.allocLink(tile);
            var link = tile.links[idx];
            link.ref = nei[k];
            link.edge = j;
            link.side = dir;
            link.next = poly.firstLink;
            poly.firstLink = idx; // Compress portal limits to a byte value.

            if (dir == 0 || dir == 4) {
              tmin = (neia[k * 2 + 0] - tile.data.verts[va + 2]) / (tile.data.verts[vb + 2] - tile.data.verts[va + 2]);
              tmax = (neia[k * 2 + 1] - tile.data.verts[va + 2]) / (tile.data.verts[vb + 2] - tile.data.verts[va + 2]);

              if (tmin > tmax) {
                temp = tmin;
                tmin = tmax;
                tmax = temp;
              }

              link.bmin = Math.floor(_DetourCommon["default"].clamp(tmin, 0.0, 1.0) * 255.0);
              link.bmax = Math.floor(_DetourCommon["default"].clamp(tmax, 0.0, 1.0) * 255.0);
            } else if (dir == 2 || dir == 6) {
              tmin = (neia[k * 2 + 0] - tile.data.verts[va]) / (tile.data.verts[vb] - tile.data.verts[va]);
              tmax = (neia[k * 2 + 1] - tile.data.verts[va]) / (tile.data.verts[vb] - tile.data.verts[va]);

              if (tmin > tmax) {
                temp = tmin;
                tmin = tmax;
                tmax = temp;
              }

              link.bmin = Math.floor(_DetourCommon["default"].clamp(tmin, 0.0, 1.0) * 255.0);
              link.bmax = Math.floor(_DetourCommon["default"].clamp(tmax, 0.0, 1.0) * 255.0);
            }
          }
        }
      }
    }
  }, {
    key: "connectExtOffMeshLinks",
    value: function connectExtOffMeshLinks(tile, target, side) {
      if (tile == null) return; // Connect off-mesh links.
      // We are interested on links which land from target tile to this tile.

      var oppositeSide = side == -1 ? 0xff : _DetourCommon["default"].oppositeTile(side);

      for (var i = 0; i < target.data.header.offMeshConCount; ++i) {
        var targetCon = target.data.offMeshCons[i];
        if (targetCon.side != oppositeSide) continue;
        var targetPoly = target.data.polys[targetCon.poly]; // Skip off-mesh connections which start location could not be
        // connected at all.

        if (targetPoly.firstLink == NavMesh.DT_NULL_LINK) continue;
        var ext = [targetCon.rad, target.data.header.walkableClimb, targetCon.rad]; // Find polygon to connect to.

        var p = new Array(3);
        p[0] = targetCon.pos[3];
        p[1] = targetCon.pos[4];
        p[2] = targetCon.pos[5];
        var nearest = this.findNearestPolyInTile(tile, p, ext);
        var ref = nearest.getNearestRef();
        if (ref == 0) continue;
        var nearestPt = nearest.getNearestPos(); // findNearestPoly may return too optimistic results, further check
        // to make sure.

        if (_DetourCommon["default"].sqr(nearestPt[0] - p[0]) + _DetourCommon["default"].sqr(nearestPt[2] - p[2]) > _DetourCommon["default"].sqr(targetCon.rad)) continue; // Make sure the location is on curren mesh.

        target.data.verts[targetPoly.verts[1] * 3] = nearestPt[0];
        target.data.verts[targetPoly.verts[1] * 3 + 1] = nearestPt[1];
        target.data.verts[targetPoly.verts[1] * 3 + 2] = nearestPt[2]; // let off-mesh connection to target poly.

        var idx = this.allocLink(target);
        var link = target.links[idx];
        link.ref = ref;
        link.edge = 1;
        link.side = oppositeSide;
        link.bmin = link.bmax = 0; // Add to linked list.

        link.next = targetPoly.firstLink;
        targetPoly.firstLink = idx; // let target poly to off-mesh connection.

        if ((targetCon.flags & NavMesh.DT_OFFMESH_CON_BIDIR) != 0) {
          var tidx = this.allocLink(tile);
          var landPolyIdx = NavMesh.decodePolyIdPoly(ref);
          var landPoly = tile.data.polys[landPolyIdx];
          link = tile.links[tidx];
          link.ref = NavMesh.or(this.getPolyRefBase(target), targetCon.poly);
          link.edge = 0xff;
          link.side = side == -1 ? 0xff : side;
          link.bmin = link.bmax = 0; // Add to linked list.

          link.next = landPoly.firstLink;
          landPoly.firstLink = tidx;
        }
      }
    }
  }, {
    key: "findConnectingPolys",
    value: function findConnectingPolys(verts, va, vb, tile, side, maxcon) {
      if (tile == null) return [null, null, 0];
      var con = new Array(maxcon);
      var conarea = new Array(maxcon * 2);
      var amin = new Array(2);
      var amax = new Array(2);
      calcSlabEndPoints(verts, va, vb, amin, amax, side);
      apos = getSlabCoord(verts, va, side); // Remove links pointing to 'side' and compact the links array.

      var bmin = new Array(2);
      var bmax = new Array(2);
      var m = NavMesh.or(NavMesh.DT_EXT_LINK, side);
      var n = 0;
      var base = this.getPolyRefBase(tile);

      for (var i = 0; i < tile.data.header.polyCount; ++i) {
        var poly = tile.data.polys[i];
        var nv = poly.vertCount;

        for (var j = 0; j < nv; ++j) {
          // Skip edges which do not poPoly to the right side.
          if (poly.neis[j] != m) continue;
          var vc = poly.verts[j] * 3;
          var vd = poly.verts[(j + 1) % nv] * 3;
          bpos = getSlabCoord(tile.data.verts, vc, side); // Segments are not close enough.

          if (Math.abs(apos - bpos) > 0.01) continue; // Check if the segments touch.

          calcSlabEndPoints(tile.data.verts, vc, vd, bmin, bmax, side);
          if (!overlapSlabs(amin, amax, bmin, bmax, 0.01, tile.data.header.walkableClimb)) continue; // Add return value.

          if (n < maxcon) {
            conarea[n * 2 + 0] = Math.max(amin[0], bmin[0]);
            conarea[n * 2 + 1] = Math.min(amax[0], bmax[0]);
            con[n] = NavMesh.or(base, i);
            n++;
          }

          break;
        }
      }

      return [con, conarea, n];
    }
  }, {
    key: "overlapSlabs",
    value: function overlapSlabs(amin, amax, bmin, bmax, px, py) {
      // Check for horizontal overlap.
      // The segment is shrunken a little so that slabs which touch
      // at end points are not connected.
      minx = Math.max(amin[0] + px, bmin[0] + px);
      maxx = Math.min(amax[0] - px, bmax[0] - px);
      if (minx > maxx) return false; // Check vertical overlap.

      ad = (amax[1] - amin[1]) / (amax[0] - amin[0]);
      ak = amin[1] - ad * amin[0];
      bd = (bmax[1] - bmin[1]) / (bmax[0] - bmin[0]);
      bk = bmin[1] - bd * bmin[0];
      aminy = ad * minx + ak;
      amaxy = ad * maxx + ak;
      bminy = bd * minx + bk;
      bmaxy = bd * maxx + bk;
      dmin = bminy - aminy;
      dmax = bmaxy - amaxy; // Crossing segments always overlap.

      if (dmin * dmax < 0) return true; // Check for overlap at endpoints.

      thr = py * 2 * (py * 2);
      if (dmin * dmin <= thr || dmax * dmax <= thr) return true;
      return false;
    }
    /**
     * Builds internal polygons links for a tile.
     * 
     * @param tile
     */

  }, {
    key: "baseOffMeshLinks",
    value: function baseOffMeshLinks(tile) {
      if (tile == null) return;
      var base = this.getPolyRefBase(tile); // Base off-mesh connection start points.

      for (var i = 0; i < tile.data.header.offMeshConCount; ++i) {
        var con = tile.data.offMeshCons[i];
        var poly = tile.data.polys[con.poly];
        var ext = [con.rad, tile.data.header.walkableClimb, con.rad]; // Find polygon to connect to.

        var nearestPoly = this.findNearestPolyInTile(tile, con.pos, ext);
        var ref = nearestPoly.getNearestRef();
        if (ref == 0) continue;
        var p = con.pos; // First vertex

        var nearestPt = nearestPoly.getNearestPos(); // findNearestPoly may return too optimistic results, further check
        // to make sure.

        if (_DetourCommon["default"].sqr(nearestPt[0] - p[0]) + _DetourCommon["default"].sqr(nearestPt[2] - p[2]) > _DetourCommon["default"].sqr(con.rad)) continue; // Make sure the location is on current mesh.

        tile.data.verts[poly.verts[0] * 3] = nearestPt[0];
        tile.data.verts[poly.verts[0] * 3 + 1] = nearestPt[1];
        tile.data.verts[poly.verts[0] * 3 + 2] = nearestPt[2]; // let off-mesh connection to target poly.

        var idx = this.allocLink(tile);
        var link = tile.links[idx];
        link.ref = ref;
        link.edge = 0;
        link.side = 0xff;
        link.bmin = link.bmax = 0; // Add to linked list.

        link.next = poly.firstLink;
        poly.firstLink = idx; // Start end-poPoly is always connect back to off-mesh connection.

        var tidx = this.allocLink(tile);
        var landPolyIdx = NavMesh.decodePolyIdPoly(ref);
        var landPoly = tile.data.polys[landPolyIdx];
        link = tile.links[tidx];
        link.ref = NavMesh.or(base, con.poly);
        link.edge = 0xff;
        link.side = 0xff;
        link.bmin = link.bmax = 0; // Add to linked list.

        link.next = landPoly.firstLink;
        landPoly.firstLink = tidx;
      }
    }
    /**
     * Returns closest poPoly on polygon.
     * 
     * @param ref
     * @param pos
     * @return
     */

  }, {
    key: "closestPointOnPoly",
    value: function closestPointOnPoly(ref, pos) {
      var tileAndPoly = this.getTileAndPolyByRefUnsafe(ref);
      var tile = tileAndPoly[0];
      var poly = tileAndPoly[1]; // Off-mesh connections don't have detail polygons.

      if (poly.getType() == _Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION) {
        var v0 = poly.verts[0] * 3;
        var v1 = poly.verts[1] * 3;

        var d0 = _DetourCommon["default"].vDist3(pos, tile.data.verts, v0);

        var d1 = _DetourCommon["default"].vDist3(pos, tile.data.verts, v1);

        var u = d0 / (d0 + d1);

        var _closest = _DetourCommon["default"].vLerp4(tile.data.verts, v0, v1, u);

        return new _ClosestPointOnPolyResult["default"](false, _closest);
      } // Clamp poPoly to be inside the polygon.


      var verts = new Array(this.m_maxVertPerPoly * 3);
      var edged = new Array(this.m_maxVertPerPoly);
      var edget = new Array(this.m_maxVertPerPoly);
      var nv = poly.vertCount;

      for (var i = 0; i < nv; ++i) {
        arraycopy(tile.data.verts, poly.verts[i] * 3, verts, i * 3, 3);
      }

      var posOverPoly = false;
      var closest = new Array(3);

      _DetourCommon["default"].vCopy(closest, pos);

      if (!_DetourCommon["default"].distancePtPolyEdgesSqr(pos, verts, nv, edged, edget)) {
        // PoPoly is outside the polygon, dtClamp to nearest edge.
        var _dmin = edged[0];
        var imin = 0;

        for (var _i = 1; _i < nv; ++_i) {
          if (edged[_i] < _dmin) {
            _dmin = edged[_i];
            imin = _i;
          }
        }

        var va = imin * 3;
        var vb = (imin + 1) % nv * 3;
        closest = _DetourCommon["default"].vLerp4(verts, va, vb, edget[imin]);
        posOverPoly = false;
      } else {
        posOverPoly = true;
      } // Find height at the location.


      var ip = poly.index;

      if (tile.data.detailMeshes != null && tile.data.detailMeshes.length > ip) {
        var pd = tile.data.detailMeshes[ip];

        for (var j = 0; j < pd.triCount; ++j) {
          var t = (pd.triBase + j) * 4;
          var v = []; //Was [3][]

          for (var k = 0; k < 3; ++k) {
            if (tile.data.detailTris[t + k] < poly.vertCount) {
              var index = poly.verts[tile.data.detailTris[t + k]] * 3;
              v[k] = [tile.data.verts[index], tile.data.verts[index + 1], tile.data.verts[index + 2]];
            } else {
              var _index = (pd.vertBase + (tile.data.detailTris[t + k] - poly.vertCount)) * 3;

              v[k] = [tile.data.detailVerts[_index], tile.data.detailVerts[_index + 1], tile.data.detailVerts[_index + 2]];
            }
          }

          var heightResult = _DetourCommon["default"].closestHeightPointTriangle(closest, v[0], v[1], v[2]);

          if (heightResult[0]) {
            closest[1] = heightResult[1];
            break;
          }
        }
      }

      return new _ClosestPointOnPolyResult["default"](posOverPoly, closest);
    }
  }, {
    key: "findNearestPolyInTile",
    value: function findNearestPolyInTile(tile, center, extents) {
      var nearestPt = null;

      var bmin = _DetourCommon["default"].vSub(center, extents);

      var bmax = _DetourCommon["default"].vAdd(center, extents); // Get nearby polygons from proximity grid.


      var polys = this.queryPolygonsInTile(tile, bmin, bmax); // Find nearest polygon amongst the nearby polygons.

      var nearest = 0;
      var nearestDistanceSqr = Number.MAX_VALUE;

      for (var i = 0; i < polys.length; ++i) {
        var ref = polys[i];
        var d = 0;
        var cpp = this.closestPointOnPoly(ref, center);
        var posOverPoly = cpp.isPosOverPoly();
        var closestPtPoly = cpp.getClosest(); // If a poPoly is directly over a polygon and closer than
        // climb height, favor that instead of straight line nearest point.

        var diff = _DetourCommon["default"].vSub(center, closestPtPoly);

        if (posOverPoly) {
          d = Math.abs(diff[1]) - tile.data.header.walkableClimb;
          d = d > 0 ? d * d : 0;
        } else {
          d = _DetourCommon["default"].vLenSqr(diff);
        }

        if (d < nearestDistanceSqr) {
          nearestPt = closestPtPoly;
          nearestDistanceSqr = d;
          nearest = ref;
        }
      }

      return new _FindNearestPolyResult["default"](nearest, nearestPt);
    }
  }, {
    key: "getTileAt",
    value: function getTileAt(x, y, layer) {
      // Find tile based on hash.
      var h = NavMesh.computeTileHash(x, y, this.m_tileLutMask);
      var tile = this.m_posLookup[h];

      while (tile != null) {
        if (tile.data.header != null && tile.data.header.x == x && tile.data.header.y == y && tile.data.header.layer == layer) {
          return tile;
        }

        tile = tile.next;
      }

      return null;
    }
  }, {
    key: "getNeighbourTilesAt",
    value: function getNeighbourTilesAt(x, y, side) {
      var nx = x,
          ny = y;

      switch (side) {
        case 0:
          nx++;
          break;

        case 1:
          nx++;
          ny++;
          break;

        case 2:
          ny++;
          break;

        case 3:
          nx--;
          ny++;
          break;

        case 4:
          nx--;
          break;

        case 5:
          nx--;
          ny--;
          break;

        case 6:
          ny--;
          break;

        case 7:
          nx++;
          ny--;
          break;
      }

      return this.getTilesAt(nx, ny);
    }
  }, {
    key: "getTilesAt",
    value: function getTilesAt(x, y) {
      var tiles = []; // Find tile based on hash.

      var h = NavMesh.computeTileHash(x, y, this.m_tileLutMask);
      var tile = this.m_posLookup[h];

      while (tile != null) {
        if (tile.data.header != null && tile.data.header.x == x && tile.data.header.y == y) {
          tiles.push(tile);
        }

        tile = tile.next;
      }

      return tiles;
    }
  }, {
    key: "getTileRefAt",
    value: function getTileRefAt(x, y, layer) {
      // Find tile based on hash.
      var h = NavMesh.computeTileHash(x, y, this.m_tileLutMask);
      var tile = this.m_posLookup[h];

      while (tile != null) {
        if (tile.data.header != null && tile.data.header.x == x && tile.data.header.y == y && tile.data.header.layer == layer) {
          return this.getTileRef(tile);
        }

        tile = tile.next;
      }

      return 0;
    }
  }, {
    key: "getTileByRef",
    value: function getTileByRef(ref) {
      if (ref == 0) return null;
      var tileIndex = decodePolyIdTile(ref);
      var tileSalt = decodePolyIdSalt(ref);
      if (tileIndex >= this.m_maxTiles) return null;
      var tile = this.m_tiles[tileIndex];
      if (tile.salt != tileSalt) return null;
      return tile;
    }
  }, {
    key: "getTileRef",
    value: function getTileRef(tile) {
      if (tile == null) return 0;
      var it = tile.index;
      return NavMesh.encodePolyId(tile.salt, it, 0);
    }
  }, {
    key: "getOffMeshConnectionPolyEndPoints",
    /// @par
    ///
    /// Off-mesh connections are stored in the navigation mesh as special
    /// 2-vertex
    /// polygons with a single edge. At least one of the vertices is expected to
    /// be
    /// inside a normal polygon. So an off-mesh connection is "entered" from a
    /// normal polygon at one of its endpoints. This is the polygon identified
    /// by
    /// the prevRef parameter.
    value: function getOffMeshConnectionPolyEndPoints(prevRef, polyRef) {
      if (polyRef == 0) throw new IllegalArgumentException("polyRef = 0"); // Get current polygon

      var saltitip = NavMesh.decodePolyId(polyRef);
      var salt = saltitip[0];
      var it = saltitip[1];
      var ip = saltitip[2];

      if (it >= this.m_maxTiles) {
        throw new IllegalArgumentException("Invalid tile ID > max tiles");
      }

      if (this.m_tiles[it].salt != salt || this.m_tiles[it].data.header == null) {
        throw new IllegalArgumentException("Invalid salt or missing tile header");
      }

      var tile = this.m_tiles[it];

      if (ip >= tile.data.header.polyCount) {
        throw new IllegalArgumentException("Invalid poly ID > poly count");
      }

      var poly = tile.data.polys[ip]; // Make sure that the current poly is indeed off-mesh link.

      if (poly.getType() != _Poly["default"].DT_POLYTYPE_OFFMESH_CONNECTION) throw new IllegalArgumentException("Invalid poly type"); // Figure out which way to hand out the vertices.

      var idx0 = 0,
          idx1 = 1; // Find link that points to first vertex.

      for (var i = poly.firstLink; i != NavMesh.DT_NULL_LINK; i = tile.links[i].next) {
        if (tile.links[i].edge == 0) {
          if (tile.links[i].ref != prevRef) {
            idx0 = 1;
            idx1 = 0;
          }

          break;
        }
      }

      var startPos = new Array(3);
      var endPos = new Array(3);

      _DetourCommon["default"].vCopy(startPos, tile.data.verts, poly.verts[idx0] * 3);

      _DetourCommon["default"].vCopy(endPos, tile.data.verts, poly.verts[idx1] * 3);

      return [startPos, endPos];
    }
  }, {
    key: "getMaxVertsPerPoly",
    value: function getMaxVertsPerPoly() {
      return this.m_maxVertPerPoly;
    }
  }, {
    key: "getTileCount",
    value: function getTileCount() {
      return this.m_tileCount;
    }
  }], [{
    key: "getNavMeshParams",
    value: function getNavMeshParams(data) {
      var params = new _NavMeshParams["default"]();

      _DetourCommon["default"].vCopy(params.orig, data.header.bmin);

      params.tileWidth = data.header.bmax[0] - data.header.bmin[0];
      params.tileHeight = data.header.bmax[2] - data.header.bmin[2];
      params.maxTiles = 1;
      params.maxPolys = data.header.polyCount;
      return params;
    }
  }, {
    key: "getSlabCoord",
    value: function getSlabCoord(verts, va, side) {
      if (side == 0 || side == 4) return verts[va];else if (side == 2 || side == 6) return verts[va + 2];
      return 0;
    }
  }, {
    key: "calcSlabEndPoints",
    value: function calcSlabEndPoints(verts, va, vb, bmin, bmax, side) {
      if (side == 0 || side == 4) {
        if (verts[va + 2] < verts[vb + 2]) {
          bmin[0] = verts[va + 2];
          bmin[1] = verts[va + 1];
          bmax[0] = verts[vb + 2];
          bmax[1] = verts[vb + 1];
        } else {
          bmin[0] = verts[vb + 2];
          bmin[1] = verts[vb + 1];
          bmax[0] = verts[va + 2];
          bmax[1] = verts[va + 1];
        }
      } else if (side == 2 || side == 6) {
        if (verts[va + 0] < verts[vb + 0]) {
          bmin[0] = verts[va + 0];
          bmin[1] = verts[va + 1];
          bmax[0] = verts[vb + 0];
          bmax[1] = verts[vb + 1];
        } else {
          bmin[0] = verts[vb + 0];
          bmin[1] = verts[vb + 1];
          bmax[0] = verts[va + 0];
          bmax[1] = verts[va + 1];
        }
      }
    }
  }, {
    key: "computeTileHash",
    value: function computeTileHash(x, y, mask) {
      var h1 = 0x8da6b343; // Large multiplicative constants;

      var h2 = 0xd8163841; // here arbitrarily chosen primes

      var n = h1 * x + h2 * y;
      return n & mask;
    }
  }]);

  return NavMesh;
}();

_defineProperty(NavMesh, "DT_SALT_BITS", 16);

_defineProperty(NavMesh, "DT_TILE_BITS", 28);

_defineProperty(NavMesh, "DT_POLY_BITS", 20);

_defineProperty(NavMesh, "DT_EXT_LINK", 0x8000);

_defineProperty(NavMesh, "DT_NULL_LINK", -1);

_defineProperty(NavMesh, "DT_OFFMESH_CON_BIDIR", 1);

_defineProperty(NavMesh, "DT_MAX_AREAS", 64);

_defineProperty(NavMesh, "DT_RAY_CAST_LIMIT_PROPORTIONS", 50.0);

var _default = NavMesh;
exports["default"] = _default;