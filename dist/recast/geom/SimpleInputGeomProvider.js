"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _RecastVectors = _interopRequireDefault(require("../RecastVectors.js"));

var _TriMesh = _interopRequireDefault(require("./TriMesh.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var SimpleInputGeomProvider
/*extends InputGeomProvider*/
= /*#__PURE__*/function () {
  _createClass(SimpleInputGeomProvider, null, [{
    key: "mapFaces",
    value: function mapFaces(meshFaces) {
      var faces = new Array(meshFaces.length);

      for (var i = 0; i < faces.length; i++) {
        faces[i] = meshFaces[i];
      }

      return faces;
    }
  }, {
    key: "mapVertices",
    value: function mapVertices(vertexPositions) {
      var vertices = new Array(vertexPositions.length);

      for (var i = 0; i < vertices.length; i++) {
        vertices[i] = vertexPositions[i];
      }

      return vertices;
    }
  }]);

  function SimpleInputGeomProvider(vertices, faces) {
    _classCallCheck(this, SimpleInputGeomProvider);

    _defineProperty(this, "vertices", void 0);

    _defineProperty(this, "faces", void 0);

    _defineProperty(this, "bmin", void 0);

    _defineProperty(this, "bmax", void 0);

    _defineProperty(this, "volumes", []);

    this.vertices = vertices;
    this.faces = faces;
    this.bmin = new Array(3);
    this.bmax = new Array(3);

    _RecastVectors["default"].copy3(this.bmin, vertices, 0);

    _RecastVectors["default"].copy3(this.bmax, vertices, 0);

    for (var i = 1; i < vertices.length / 3; i++) {
      _RecastVectors["default"].min(this.bmin, vertices, i * 3);

      _RecastVectors["default"].max(this.bmax, vertices, i * 3);
    }
  }

  _createClass(SimpleInputGeomProvider, [{
    key: "getMeshBoundsMin",
    value: function getMeshBoundsMin() {
      return this.bmin;
    }
  }, {
    key: "getMeshBoundsMax",
    value: function getMeshBoundsMax() {
      return this.bmax;
    }
  }, {
    key: "getConvexVolumes",
    value: function getConvexVolumes() {
      return [];
    }
  }, {
    key: "addConvexVolume",
    value: function addConvexVolume(verts, minh, maxh, areaMod) {
      var vol = new ConvexVolume();
      vol.hmin = minh;
      vol.hmax = maxh;
      vol.verts = verts;
      vol.areaMod = areaMod;
      this.volumes.push(vol);
    }
  }, {
    key: "meshes",
    value: function meshes() {
      // return Collections.singPolyonList(new TriMesh(vertices, faces));
      return [new _TriMesh["default"](this.vertices, this.faces)];
    }
  }]);

  return SimpleInputGeomProvider;
}();

_defineProperty(SimpleInputGeomProvider, "fromIndeces", function (vertexPositions, meshFaces) {
  return new SimpleInputGeomProvider(this.mapVertices(vertexPositions), this.mapFaces(meshFaces));
});

var _default = SimpleInputGeomProvider;
exports["default"] = _default;