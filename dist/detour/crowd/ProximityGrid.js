"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = void 0;

var _temp;

function _toConsumableArray(arr) { return _arrayWithoutHoles(arr) || _iterableToArray(arr) || _unsupportedIterableToArray(arr) || _nonIterableSpread(); }

function _nonIterableSpread() { throw new TypeError("Invalid attempt to spread non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _iterableToArray(iter) { if (typeof Symbol !== "undefined" && Symbol.iterator in Object(iter)) return Array.from(iter); }

function _arrayWithoutHoles(arr) { if (Array.isArray(arr)) return _arrayLikeToArray(arr); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

/*
Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
Recast4J Copyright (c) 2015 Piotr Piastucki piotr@jtilia.org

This software is provided 'as-is', without any express or implied
warranty.  In no event will the authors be held liable for any damages
arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:
1. The origin of this software must not be misrepresented; you must not
 claim that you wrote the original software. If you use this software
 in a product, an acknowledgment in the product documentation would be
 appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
 misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
var ProximityGrid = /*#__PURE__*/function () {
  function ProximityGrid(m_cellSize, m_invCellSize) {
    _classCallCheck(this, ProximityGrid);

    _defineProperty(this, "m_cellSize", void 0);

    _defineProperty(this, "m_invCellSize", void 0);

    _defineProperty(this, "items", {});

    _defineProperty(this, "m_bounds", new Array(4));

    this.m_cellSize = m_cellSize;
    this.m_invCellSize = m_invCellSize;
    this.items = {};
  }

  _createClass(ProximityGrid, [{
    key: "clear",
    value: function clear() {
      this.items = {};
      this.m_bounds[0] = 0xfff;
      this.m_bounds[1] = 0xfff;
      this.m_bounds[2] = -0xfff;
      this.m_bounds[3] = -0xfff;
    }
  }, {
    key: "addItem",
    value: function addItem(id, minx, miny, maxx, maxy) {
      var iminx = Math.floor(minx * this.m_invCellSize);
      var iminy = Math.floor(miny * this.m_invCellSize);
      var imaxx = Math.floor(maxx * this.m_invCellSize);
      var imaxy = Math.floor(maxy * this.m_invCellSize);
      this.m_bounds[0] = Math.min(this.m_bounds[0], iminx);
      this.m_bounds[1] = Math.min(this.m_bounds[1], iminy);
      this.m_bounds[2] = Math.min(this.m_bounds[2], imaxx);
      this.m_bounds[3] = Math.min(this.m_bounds[3], imaxy);

      for (var _y = iminy; _y <= imaxy; ++_y) {
        for (var _x = iminx; _x <= imaxx; ++_x) {
          var _key = new ProximityGrid.ItemKey(_x, _y);

          var _ids = this.items[JSON.stringify(_key)];

          if (_ids == null) {
            _ids = [];
            this.items[JSON.stringify(_key)] = _ids;
          }

          _ids.push(id);
        }
      }
    }
  }, {
    key: "queryItems",
    value: function queryItems(minx, miny, maxx, maxy) {
      var iminx = Math.floor(minx * this.m_invCellSize);
      var iminy = Math.floor(miny * this.m_invCellSize);
      var imaxx = Math.floor(maxx * this.m_invCellSize);
      var imaxy = Math.floor(maxy * this.m_invCellSize);
      var result = new Set();

      for (var _y2 = iminy; _y2 <= imaxy; ++_y2) {
        for (var _x2 = iminx; _x2 <= imaxx; ++_x2) {
          var _key2 = new ProximityGrid.ItemKey(_x2, _y2);

          var _ids2 = this.items[JSON.stringify(_key2)];

          if (_ids2 != null) {
            result.add.apply(result, _toConsumableArray(_ids2));
          }
        }
      }

      return result;
    }
  }, {
    key: "getItemCountAt",
    value: function getItemCountAt(x, y) {
      key = new ItemKey(x, y);
      ids = this.items[key];
      return ids != null ? ids.length : 0;
    }
  }]);

  return ProximityGrid;
}();

_defineProperty(ProximityGrid, "ItemKey", (_temp = /*#__PURE__*/function () {
  function ItemKey(x, y) {
    _classCallCheck(this, ItemKey);

    _defineProperty(this, "x", void 0);

    _defineProperty(this, "y", void 0);

    this.x = x;
    this.y = y;
  }

  _createClass(ItemKey, [{
    key: "hashCode",
    value: function hashCode() {
      prime = 31;
      result = 1;
      result = prime * result + x;
      result = prime * result + y;
      return result;
    }
  }, {
    key: "equals",
    value: function equals(obj) {
      if (this == obj) return true;
      if (obj == null) return false;
      if (getClass() != obj.getClass()) return false;
      other = obj;
      if (x != other.x) return false;
      if (y != other.y) return false;
      return true;
    }
  }]);

  return ItemKey;
}(), _temp));

var _default = ProximityGrid;
exports["default"] = _default;