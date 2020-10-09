"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _DetourCommon = _interopRequireDefault(require("./DetourCommon.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

//TODO: (PP) Add comments
var StraightPathItem = /*#__PURE__*/function () {
  function StraightPathItem(pos, flags, ref) {
    _classCallCheck(this, StraightPathItem);

    _defineProperty(this, "pos", []);

    _defineProperty(this, "flags", 0);

    _defineProperty(this, "ref", 0);

    this.pos = _DetourCommon["default"].vCopy_return(pos);
    this.flags = flags;
    this.ref = ref;
  }

  _createClass(StraightPathItem, [{
    key: "getPos",
    value: function getPos() {
      return this.pos;
    }
  }, {
    key: "getFlags",
    value: function getFlags() {
      return this.flags;
    }
  }, {
    key: "getRef",
    value: function getRef() {
      return this.ref;
    }
  }]);

  return StraightPathItem;
}();

var _default = StraightPathItem;
exports["default"] = _default;