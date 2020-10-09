"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var AreaModification = /*#__PURE__*/function () {
  /**
   * Mask is set to all available bits, which means value is fully applied
   * 
   * @param value
   *            The area id to apply. [Limit: <= #RC_AREA_FLAGS_MASK]
   */
  // constructor(value) {
  // 		this.value = value;
  // 		this.mask = AreaModification.RC_AREA_FLAGS_MASK;
  // 	}

  /**
   * 
   * @param value
   *            The area id to apply. [Limit: <= #RC_AREA_FLAGS_MASK]
   * @param mask
   *            Bitwise mask used when applying value. [Limit: <= #RC_AREA_FLAGS_MASK]
   */
  function AreaModification(value) {
    var mask = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : AreaModification.RC_AREA_FLAGS_MASK;

    _classCallCheck(this, AreaModification);

    _defineProperty(this, "value", void 0);

    _defineProperty(this, "mask", void 0);

    this.value = value;
    this.mask = mask;
  }

  _createClass(AreaModification, [{
    key: "getMaskedValue",
    value: function getMaskedValue() {
      return this.value & this.mask;
    }
  }, {
    key: "apply",
    value: function apply(area) {
      return this.value & this.mask | area & ~this.mask;
    }
  }]);

  return AreaModification;
}();

_defineProperty(AreaModification, "RC_AREA_FLAGS_MASK", 0x3F);

_defineProperty(AreaModification, "clone", function (other) {
  return new AreaModification(other.value, other.mask); // this.value = other.value;
  // this.mask = other.mask;
});

var _default = AreaModification;
exports["default"] = _default;