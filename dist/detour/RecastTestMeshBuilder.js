"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _NavMesh = _interopRequireDefault(require("./NavMesh.js"));

var _ObjImporter = _interopRequireDefault(require("../recast/ObjImporter.js"));

var _RecastConstants = _interopRequireDefault(require("../recast/RecastConstants.js"));

var _RecastConfig = _interopRequireDefault(require("../recast/RecastConfig.js"));

var _SampleAreaModifications = _interopRequireDefault(require("../detour/SampleAreaModifications.js"));

var _RecastBuilderConfig = _interopRequireDefault(require("../recast/RecastBuilderConfig.js"));

var _RecastBuilder = _interopRequireDefault(require("../recast/RecastBuilder.js"));

var _NavMeshDataCreateParams = _interopRequireDefault(require("./NavMeshDataCreateParams.js"));

var _NavMeshBuilder = _interopRequireDefault(require("./NavMeshBuilder.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var RecastTestMeshBuilder = /*#__PURE__*/function () {
  function RecastTestMeshBuilder(m_geom, m_partitionType, m_cellSize, m_cellHeight, m_agentHeight, m_agentRadius, m_agentMaxClimb, m_agentMaxSlope, m_regionMinSize, m_regionMergeSize, m_edgeMaxLen, m_edgeMaxError, m_vertsPerPoly, m_detailSampleDist, m_detailSampleMaxError) {
    _classCallCheck(this, RecastTestMeshBuilder);

    _defineProperty(this, "meshData", void 0);

    var cfg = new _RecastConfig["default"](m_partitionType, m_cellSize, m_cellHeight, m_agentHeight, m_agentRadius, m_agentMaxClimb, m_agentMaxSlope, m_regionMinSize, m_regionMergeSize, m_edgeMaxLen, m_edgeMaxError, m_vertsPerPoly, m_detailSampleDist, m_detailSampleMaxError, 0, _SampleAreaModifications["default"].SAMPLE_AREAMOD_GROUND);
    var bcfg = new _RecastBuilderConfig["default"](cfg, m_geom.getMeshBoundsMin(), m_geom.getMeshBoundsMax());
    var rcBuilder = new _RecastBuilder["default"]();
    var rcResult = rcBuilder.build(m_geom, bcfg);
    var m_pmesh = rcResult.getMesh();

    for (var i = 0; i < m_pmesh.npolys; ++i) {
      m_pmesh.flags[i] = 1;
    }

    var m_dmesh = rcResult.getMeshDetail();
    var params = new _NavMeshDataCreateParams["default"]();
    params.verts = m_pmesh.verts;
    params.vertCount = m_pmesh.nverts;
    params.polys = m_pmesh.polys;
    params.polyAreas = m_pmesh.areas;
    params.polyFlags = m_pmesh.flags;
    params.polyCount = m_pmesh.npolys;
    params.nvp = m_pmesh.nvp;
    params.detailMeshes = m_dmesh.meshes;
    params.detailVerts = m_dmesh.verts;
    params.detailVertsCount = m_dmesh.nverts;
    params.detailTris = m_dmesh.tris;
    params.detailTriCount = m_dmesh.ntris;
    params.walkableHeight = m_agentHeight;
    params.walkableRadius = m_agentRadius;
    params.walkableClimb = m_agentMaxClimb;
    params.bmin = m_pmesh.bmin;
    params.bmax = m_pmesh.bmax;
    params.cs = m_cellSize;
    params.ch = m_cellHeight;
    params.buildBvTree = true;
    params.offMeshConVerts = new Array(6);
    params.offMeshConVerts[0] = 0.1;
    params.offMeshConVerts[1] = 0.2;
    params.offMeshConVerts[2] = 0.3;
    params.offMeshConVerts[3] = 0.4;
    params.offMeshConVerts[4] = 0.5;
    params.offMeshConVerts[5] = 0.6;
    params.offMeshConRad = new Array(1);
    params.offMeshConRad[0] = 0.1;
    params.offMeshConDir = new Array(1);
    params.offMeshConDir[0] = 1;
    params.offMeshConAreas = new Array(1);
    params.offMeshConAreas[0] = 2;
    params.offMeshConFlags = new Array(1);
    params.offMeshConFlags[0] = 12;
    params.offMeshConUserID = new Array(1);
    params.offMeshConUserID[0] = 0x4567;
    params.offMeshConCount = 1;
    this.meshData = _NavMeshBuilder["default"].createNavMeshData(params);
  }

  _createClass(RecastTestMeshBuilder, [{
    key: "getMeshData",
    value: function getMeshData() {
      return this.meshData;
    }
  }]);

  return RecastTestMeshBuilder;
}();

_defineProperty(RecastTestMeshBuilder, "m_cellSize", 0.3);

_defineProperty(RecastTestMeshBuilder, "m_cellHeight", 0.2);

_defineProperty(RecastTestMeshBuilder, "m_agentHeight", 2.0);

_defineProperty(RecastTestMeshBuilder, "m_agentRadius", 0.6);

_defineProperty(RecastTestMeshBuilder, "m_agentMaxClimb", 0.9);

_defineProperty(RecastTestMeshBuilder, "m_agentMaxSlope", 45.0);

_defineProperty(RecastTestMeshBuilder, "m_regionMinSize", 8);

_defineProperty(RecastTestMeshBuilder, "m_regionMergeSize", 20);

_defineProperty(RecastTestMeshBuilder, "m_edgeMaxLen", 12.0);

_defineProperty(RecastTestMeshBuilder, "m_edgeMaxError", 1.3);

_defineProperty(RecastTestMeshBuilder, "m_vertsPerPoly", 6);

_defineProperty(RecastTestMeshBuilder, "m_detailSampleDist", 6.0);

_defineProperty(RecastTestMeshBuilder, "m_detailSampleMaxError", 1.0);

_defineProperty(RecastTestMeshBuilder, "fromFile", function (objFileContents) {
  return new RecastTestMeshBuilder(new _ObjImporter["default"]().load(objFileContents), _RecastConstants["default"].WATERSHED, this.m_cellSize, this.m_cellHeight, this.m_agentHeight, this.m_agentRadius, this.m_agentMaxClimb, this.m_agentMaxSlope, this.m_regionMinSize, this.m_regionMergeSize, this.m_edgeMaxLen, this.m_edgeMaxError, this.m_vertsPerPoly, this.m_detailSampleDist, this.m_detailSampleMaxError);
});

var _default = RecastTestMeshBuilder;
exports["default"] = _default;