"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _ChunkyTriMesh = _interopRequireDefault(require("./ChunkyTriMesh.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var TriMesh = /*#__PURE__*/function () {
  function TriMesh(vertices, faces) {
    _classCallCheck(this, TriMesh);

    _defineProperty(this, "vertices", void 0);

    _defineProperty(this, "faces", void 0);

    _defineProperty(this, "chunkyTriMesh", void 0);

    this.vertices = vertices;
    this.faces = faces;
    this.chunkyTriMesh = new _ChunkyTriMesh["default"](vertices, faces, faces.length / 3, 256);
  }

  _createClass(TriMesh, [{
    key: "getTris",
    value: function getTris() {
      return this.faces;
    }
  }, {
    key: "getVerts",
    value: function getVerts() {
      return this.vertices;
    }
  }, {
    key: "getChunksOverlappingRect",
    value: function getChunksOverlappingRect(bmin, bmax) {
      return chunkyTriMesh.getChunksOverlappingRect(bmin, bmax);
    }
  }]);

  return TriMesh;
}();

var _default = TriMesh;
exports["default"] = _default;