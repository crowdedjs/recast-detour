"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

var DetourBuilder = /*#__PURE__*/function () {
  function DetourBuilder() {
    _classCallCheck(this, DetourBuilder);
  }

  _createClass(DetourBuilder, [{
    key: "build",
    value: function build(params, tileX, tileY) {
      console.log("NOT BUILD");
      var data = NavMeshBuilder.createNavMeshData(params);

      if (data != null) {
        data.header.x = tileX;
        data.header.y = tileY;
      }

      return data;
    }
  }]);

  return DetourBuilder;
}();

var _default = DetourBuilder;
exports["default"] = _default;