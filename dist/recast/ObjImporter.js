"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _SimpleInputGeomProvider = _interopRequireDefault(require("./geom/SimpleInputGeomProvider.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _createForOfIteratorHelper(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var ObjImporterContext = function ObjImporterContext() {
  _classCallCheck(this, ObjImporterContext);

  _defineProperty(this, "vertexPositions", []);

  _defineProperty(this, "meshFaces", []);
};

var ObjImporter = /*#__PURE__*/function () {
  function ObjImporter() {
    _classCallCheck(this, ObjImporter);
  }

  _createClass(ObjImporter, [{
    key: "load",
    // OBJImporterContext = 
    value: function load(is) {
      var context = new ObjImporterContext();
      var reader = null;

      try {
        var slurp = is;
        var lines = slurp.split(/\r?\n/);

        for (var i = 0; i < lines.length; i++) {
          var line = lines[i];
          this.readLine(line, context);
        } // reader = new BufferedReader(new InputStreamReader(is));
        // let line;
        // while ((line = reader.readLine()) != null) {
        //     line = line.trim();
        //     readLine(line, context);
        // }

      } catch (e) {
        throw e;
      } finally {
        if (reader != null) {
          try {
            reader.close();
          } catch (e) {
            throw new RuntimeException(e.getMessage(), e);
          }
        }
      }

      return _SimpleInputGeomProvider["default"].fromIndeces(context.vertexPositions, context.meshFaces);
    }
  }, {
    key: "readLine",
    value: function readLine(line, context) {
      if (line.startsWith("v")) {
        this.readVertex(line, context);
      } else if (line.startsWith("f")) {
        this.readFace(line, context);
      }
    }
  }, {
    key: "readVertex",
    value: function readVertex(line, context) {
      if (line.startsWith("v ")) {
        var vert = this.readVector3f(line);

        var _iterator = _createForOfIteratorHelper(vert),
            _step;

        try {
          for (_iterator.s(); !(_step = _iterator.n()).done;) {
            var vp = _step.value;
            context.vertexPositions.push(vp);
          }
        } catch (err) {
          _iterator.e(err);
        } finally {
          _iterator.f();
        }
      }
    }
  }, {
    key: "readVector3f",
    value: function readVector3f(line) {
      var v = line.split(/\s+/);

      if (v.length < 4) {
        throw new RuntimeException("Invalid vector, expected 3 coordinates, found " + (v.length - 1));
      }

      return [parseFloat(v[1]), parseFloat(v[2]), parseFloat(v[3])];
    }
  }, {
    key: "readFace",
    value: function readFace(line, context) {
      var v = line.split(/\s+/);

      if (v.length < 4) {
        throw new RuntimeException("Invalid number of face vertices: 3 coordinates expected, found " + v.length);
      }

      for (var j = 0; j < v.length - 3; j++) {
        context.meshFaces.push(this.readFaceVertex(v[1], context));

        for (var i = 0; i < 2; i++) {
          context.meshFaces.push(this.readFaceVertex(v[2 + j + i], context));
        }
      }
    }
  }, {
    key: "readFaceVertex",
    value: function readFaceVertex(face, context) {
      var v = face.split(/\//);
      return this.getIndex(parseInt(v[0]), context.vertexPositions.length);
    }
  }, {
    key: "getIndex",
    value: function getIndex(posi, size) {
      if (posi > 0) {
        posi--;
      } else if (posi < 0) {
        posi = size + posi;
      } else {
        throw new RuntimeException("0 vertex index");
      }

      return posi;
    }
  }]);

  return ObjImporter;
}();

var _default = ObjImporter;
exports["default"] = _default;