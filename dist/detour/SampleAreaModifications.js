"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _AreaModification = _interopRequireDefault(require("../recast/AreaModification.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var SampleAreaModifications = function SampleAreaModifications() {
  _classCallCheck(this, SampleAreaModifications);
};

_defineProperty(SampleAreaModifications, "SAMPLE_POLYAREA_TYPE_MASK", 0x07);

_defineProperty(SampleAreaModifications, "SAMPLE_POLYAREA_TYPE_GROUND", 0x1);

_defineProperty(SampleAreaModifications, "SAMPLE_POLYAREA_TYPE_WATER", 0x2);

_defineProperty(SampleAreaModifications, "SAMPLE_POLYAREA_TYPE_ROAD", 0x3);

_defineProperty(SampleAreaModifications, "SAMPLE_POLYAREA_TYPE_DOOR", 0x4);

_defineProperty(SampleAreaModifications, "SAMPLE_POLYAREA_TYPE_GRASS", 0x5);

_defineProperty(SampleAreaModifications, "SAMPLE_POLYAREA_TYPE_JUMP", 0x6);

_defineProperty(SampleAreaModifications, "SAMPLE_AREAMOD_GROUND", new _AreaModification["default"](SampleAreaModifications.SAMPLE_POLYAREA_TYPE_GROUND, SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK));

_defineProperty(SampleAreaModifications, "SAMPLE_AREAMOD_WATER", new _AreaModification["default"](SampleAreaModifications.SAMPLE_POLYAREA_TYPE_WATER, SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK));

_defineProperty(SampleAreaModifications, "SAMPLE_AREAMOD_ROAD", new _AreaModification["default"](SampleAreaModifications.SAMPLE_POLYAREA_TYPE_ROAD, SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK));

_defineProperty(SampleAreaModifications, "SAMPLE_AREAMOD_GRASS", new _AreaModification["default"](SampleAreaModifications.SAMPLE_POLYAREA_TYPE_GRASS, SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK));

_defineProperty(SampleAreaModifications, "SAMPLE_AREAMOD_DOOR", new _AreaModification["default"](SampleAreaModifications.SAMPLE_POLYAREA_TYPE_DOOR, SampleAreaModifications.SAMPLE_POLYAREA_TYPE_DOOR));

_defineProperty(SampleAreaModifications, "SAMPLE_AREAMOD_JUMP", new _AreaModification["default"](SampleAreaModifications.SAMPLE_POLYAREA_TYPE_JUMP, SampleAreaModifications.SAMPLE_POLYAREA_TYPE_JUMP));

_defineProperty(SampleAreaModifications, "SAMPLE_POLYFLAGS_WALK", 0x01);

_defineProperty(SampleAreaModifications, "SAMPLE_POLYFLAGS_SWIM", 0x02);

_defineProperty(SampleAreaModifications, "SAMPLE_POLYFLAGS_DOOR", 0x04);

_defineProperty(SampleAreaModifications, "SAMPLE_POLYFLAGS_JUMP", 0x08);

_defineProperty(SampleAreaModifications, "SAMPLE_POLYFLAGS_DISABLED", 0x10);

_defineProperty(SampleAreaModifications, "SAMPLE_POLYFLAGS_ALL", 0xfff);

var _default = SampleAreaModifications;
exports["default"] = _default;