"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _Context = _interopRequireDefault(require("./Context.js"));

var _Heightfield = _interopRequireDefault(require("./Heightfield.js"));

var _Recast = _interopRequireDefault(require("../recast/Recast.js"));

var _RecastRasterization = _interopRequireDefault(require("./RecastRasterization.js"));

var _RecastFilter = _interopRequireDefault(require("./RecastFilter.js"));

var _RecastArea = _interopRequireDefault(require("./RecastArea.js"));

var _RecastConstants = _interopRequireDefault(require("./RecastConstants.js"));

var _RecastRegion = _interopRequireDefault(require("./RecastRegion.js"));

var _RecastContour = _interopRequireDefault(require("./RecastContour.js"));

var _RecastMesh = _interopRequireDefault(require("./RecastMesh.js"));

var _RecastMeshDetail = _interopRequireDefault(require("./RecastMeshDetail.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _createForOfIteratorHelper(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var RecastBuilderResult = /*#__PURE__*/function () {
  function RecastBuilderResult(pmesh, dmesh) {
    _classCallCheck(this, RecastBuilderResult);

    _defineProperty(this, "pmesh", void 0);

    _defineProperty(this, "dmesh", void 0);

    this.pmesh = pmesh;
    this.dmesh = dmesh;
  }

  _createClass(RecastBuilderResult, [{
    key: "getMesh",
    value: function getMesh() {
      return this.pmesh;
    }
  }, {
    key: "getMeshDetail",
    value: function getMeshDetail() {
      return this.dmesh;
    }
  }]);

  return RecastBuilderResult;
}();

var RecastBuilder = /*#__PURE__*/function () {
  // class RecastBuilderProgressListener {
  //     onProgress(completed, total);
  // }
  // RecastBuilder() {
  //     progressListener = null;
  // }
  function RecastBuilder() {
    var progressListener = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : null;

    _classCallCheck(this, RecastBuilder);

    _defineProperty(this, "progressListener", void 0);

    this.progressListener = progressListener;
  }

  _createClass(RecastBuilder, [{
    key: "buildTiles",
    value: function buildTiles(geom, cfg, threads) {
      var bmin = geom.getMeshBoundsMin();
      var bmax = geom.getMeshBoundsMax();

      var twh = _Recast["default"].calcTileCount(bmin, bmax, cfg.cs, cfg.tileSize);

      var tw = twh[0];
      var th = twh[1];
      result = null;

      if (threads == 1) {
        result = buildSingleThread(geom, cfg, bmin, bmax, tw, th);
      } else {
        result = buildMultiThread(geom, cfg, bmin, bmax, tw, th, threads);
      }

      return result;
    }
  }, {
    key: "buildSingleThread",
    value: function buildSingleThread(geom, cfg, bmin, bmax, tw, th) {
      result = new RecastBuilderResult[tw][th]();
      counter = new AtomicInteger();

      for (var x = 0; x < tw; ++x) {
        for (var y = 0; y < th; ++y) {
          result[x][y] = buildTile(geom, cfg, bmin, bmax, x, y, counter, tw * th);
        }
      }

      return result;
    }
  }, {
    key: "buildMultiThread",
    value: function buildMultiThread(geom, cfg, bmin, bmax, tw, th, threads) {
      ec = Executors.newFixedThreadPool(threads);
      result = new RecastBuilderResult[tw][th]();
      counter = new AtomicInteger();

      for (var x = 0; x < tw; ++x) {
        var _loop = function _loop(y) {
          var tx = x;
          var ty = y;
          ec.submit(function () {
            result[tx][ty] = buildTile(geom, cfg, bmin, bmax, tx, ty, counter, tw * th);
          });
        };

        for (var y = 0; y < th; ++y) {
          _loop(y);
        }
      }

      ec.shutdown();

      try {
        ec.awaitTermination(1000, TimeUnit.HOURS);
      } catch (e) {}

      return result;
    }
  }, {
    key: "buildTile",
    value: function buildTile(geom, cfg, bmin, bmax, tx, ty, counter, total) {
      result = build(geom, new RecastBuilderConfig(cfg, bmin, bmax, tx, ty, true));

      if (progressListener != null) {
        progressListener.onProgress(counter.incrementAndGet(), total);
      }

      return result;
    }
  }, {
    key: "build",
    value: function build(geom, builderCfg) {
      var cfg = builderCfg.cfg;
      var ctx = new _Context["default"]();
      var chf = this.buildCompactHeightfield(geom, builderCfg, ctx); // Partition the heightfield so that we can use simple algorithm later
      // to triangulate the walkable areas.
      // There are 3 martitioning methods, each with some pros and cons:
      // 1) Watershed partitioning
      // - the classic Recast partitioning
      // - creates the nicest tessellation
      // - usually slowest
      // - partitions the heightfield into nice regions without holes or
      // overlaps
      // - the are some corner cases where this method creates produces holes
      // and overlaps
      // - holes may appear when a small obstacles is close to large open area
      // (triangulation can handle this)
      // - overlaps may occur if you have narrow spiral corridors (i.e
      // stairs), this make triangulation to fail
      // * generally the best choice if you precompute the nacmesh, use this
      // if you have large open areas
      // 2) Monotone partioning
      // - fastest
      // - partitions the heightfield into regions without holes and overlaps
      // (guaranteed)
      // - creates let thin polygons, which sometimes causes paths with
      // detours
      // * use this if you want fast navmesh generation
      // 3) Layer partitoining
      // - quite fast
      // - partitions the heighfield into non-overlapping regions
      // - relies on the triangulation code to cope with holes (thus slower
      // than monotone partitioning)
      // - produces better triangles than monotone partitioning
      // - does not have the corner cases of watershed partitioning
      // - can be slow and create a bit ugly tessellation (still better than
      // monotone)
      // if you have large open areas with small obstacles (not a problem if
      // you use tiles)
      // * good choice to use for tiled navmesh with medium and small sized
      // tiles

      if (cfg.partitionType == _RecastConstants["default"].WATERSHED) {
        // Prepare for region partitioning, by calculating distance field
        // aPoly the walkable surface.
        _RecastRegion["default"].buildDistanceField(ctx, chf); // Partition the walkable surface into simple regions without holes.


        _RecastRegion["default"].buildRegions(ctx, chf, builderCfg.borderSize, cfg.minRegionArea, cfg.mergeRegionArea);
      } else if (cfg.partitionType == PartitionType.MONOTONE) {
        // Partition the walkable surface into simple regions without holes.
        // Monotone partitioning does not need distancefield.
        _RecastRegion["default"].buildRegionsMonotone(ctx, chf, builderCfg.borderSize, cfg.minRegionArea, cfg.mergeRegionArea);
      } else {
        // Partition the walkable surface into simple regions without holes.
        _RecastRegion["default"].buildLayerRegions(ctx, chf, builderCfg.borderSize, cfg.minRegionArea);
      } //
      // Step 5. Trace and simplify region contours.
      //
      // Create contours.


      var cset = _RecastContour["default"].buildContours(ctx, chf, cfg.maxSimplificationError, cfg.maxEdgeLen, _RecastConstants["default"].RC_CONTOUR_TESS_WALL_EDGES); //
      // Step 6. Build polygons mesh from contours.
      //


      var pmesh = _RecastMesh["default"].buildPolyMesh(ctx, cset, cfg.maxVertsPerPoly); //
      // Step 7. Create detail mesh which allows to access approximate height
      // on each polygon.
      //


      var dmesh = builderCfg.buildMeshDetail ? _RecastMeshDetail["default"].buildPolyMeshDetail(ctx, pmesh, chf, cfg.detailSampleDist, cfg.detailSampleMaxError) : null;
      return new RecastBuilderResult(pmesh, dmesh);
    }
  }, {
    key: "buildCompactHeightfield",
    value: function buildCompactHeightfield(geomProvider, builderCfg, ctx) {
      var cfg = builderCfg.cfg; //
      // Step 2. Rasterize input polygon soup.
      //
      // Allocate voxel heightfield where we rasterize our input data to.

      var solid = new _Heightfield["default"](builderCfg.width, builderCfg.height, builderCfg.bmin, builderCfg.bmax, cfg.cs, cfg.ch); // Allocate array that can hold triangle area types.
      // If you have multiple meshes you need to process, allocate
      // and array which can hold the max number of triangles you need to
      // process.
      // Find triangles which are walkable based on their slope and rasterize
      // them.
      // If your input data is multiple meshes, you can transform them here,
      // calculate
      // the are type for each of the meshes and rasterize them.

      var meshes = geomProvider.meshes();

      var _iterator = _createForOfIteratorHelper(meshes),
          _step;

      try {
        for (_iterator.s(); !(_step = _iterator.n()).done;) {
          var geom = _step.value;
          var verts = geom.getVerts();
          var tiled = cfg.tileSize > 0;
          var totaltris = 0;

          if (tiled) {
            var tbmin = new Array(2);
            var tbmax = new Array(2);
            tbmin[0] = builderCfg.bmin[0];
            tbmin[1] = builderCfg.bmin[2];
            tbmax[0] = builderCfg.bmax[0];
            tbmax[1] = builderCfg.bmax[2];
            var nodes = geom.getChunksOverlappingRect(tbmin, tbmax);

            var _iterator3 = _createForOfIteratorHelper(nodes),
                _step3;

            try {
              for (_iterator3.s(); !(_step3 = _iterator3.n()).done;) {
                var node = _step3.value;
                var tris = node.tris;
                var ntris = tris.length / 3;
                totaltris += ntris;

                var m_triareas = _Recast["default"].markWalkableTriangles(ctx, cfg.walkableSlopeAngle, verts, tris, ntris, cfg.walkableAreaMod);

                _RecastRasterization["default"].rasterizeTriangles(ctx, verts, tris, m_triareas, ntris, solid, cfg.walkableClimb);
              }
            } catch (err) {
              _iterator3.e(err);
            } finally {
              _iterator3.f();
            }
          } else {
            var _tris = geom.getTris();

            var _ntris = _tris.length / 3;

            var _m_triareas = _Recast["default"].markWalkableTriangles(ctx, cfg.walkableSlopeAngle, verts, _tris, _ntris, cfg.walkableAreaMod);

            totaltris = _ntris;

            _RecastRasterization["default"].rasterizeTrianglesA(ctx, verts, _tris, _m_triareas, _ntris, solid, cfg.walkableClimb);
          }
        } // console.log(solid.spans[0])
        //
        // Step 3. Filter walkables surfaces.
        //
        // Once all geometry is rasterized, we do initial pass of filtering to
        // remove unwanted overhangs caused by the conservative rasterization
        // as well as filter spans where the character cannot possibly stand.

      } catch (err) {
        _iterator.e(err);
      } finally {
        _iterator.f();
      }

      _RecastFilter["default"].filterLowHangingWalkableObstacles(ctx, cfg.walkableClimb, solid);

      _RecastFilter["default"].filterLedgeSpans(ctx, cfg.walkableHeight, cfg.walkableClimb, solid);

      _RecastFilter["default"].filterWalkableLowHeightSpans(ctx, cfg.walkableHeight, solid); //
      // Step 4. Partition walkable surface to simple regions.
      //
      // Compact the heightfield so that it is faster to handle from now on.
      // This will result more cache coherent data as well as the neighbours
      // between walkable cells will be calculated.


      var chf = _Recast["default"].buildCompactHeightfield(ctx, cfg.walkableHeight, cfg.walkableClimb, solid); // console.log(chf.spans[450240])
      // Erode the walkable area by agent radius.


      _RecastArea["default"].erodeWalkableArea(ctx, cfg.walkableRadius, chf); // (Optional) Mark areas.


      var _iterator2 = _createForOfIteratorHelper(geomProvider.getConvexVolumes()),
          _step2;

      try {
        for (_iterator2.s(); !(_step2 = _iterator2.n()).done;) {
          var vol = _step2.value;

          _RecastArea["default"].markConvexPolyArea(ctx, vol.verts, vol.hmin, vol.hmax, vol.areaMod, chf);
        }
      } catch (err) {
        _iterator2.e(err);
      } finally {
        _iterator2.f();
      }

      return chf;
    }
  }, {
    key: "buildLayers",
    value: function buildLayers(geom, cfg) {
      var ctx = new _Context["default"]();
      var chf = buildCompactHeightfield(geom, cfg, ctx);
      return RecastLayers.buildHeightfieldLayers(ctx, chf, cfg.borderSize, cfg.cfg.walkableHeight);
    }
  }]);

  return RecastBuilder;
}();

var _default = RecastBuilder;
exports["default"] = _default;