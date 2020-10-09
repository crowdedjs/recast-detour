"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _Node = _interopRequireDefault(require("./Node.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _createForOfIteratorHelper(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

var NodePool = /*#__PURE__*/function () {
  function NodePool() {
    _classCallCheck(this, NodePool);

    _defineProperty(this, "m_map", []);

    _defineProperty(this, "m_nodes", []);
  }

  _createClass(NodePool, [{
    key: "clear",
    value: function clear() {
      this.m_nodes = [];
      this.m_map = [];
    }
  }, {
    key: "findNodes",
    value: function findNodes(id) {
      var nodes = this.m_map[id];

      if (nodes == null) {
        nodes = [];
      }

      return nodes;
    }
  }, {
    key: "findNode",
    value: function findNode(id) {
      var nodes = this.m_map[id];

      if (nodes != null && !nodes.length == 0) {
        return nodes[0];
      }

      return null;
    }
  }, {
    key: "findNode",
    value: function findNode(id, state) {
      var nodes = this.m_map[id];

      if (nodes != null) {
        var _iterator = _createForOfIteratorHelper(nodes),
            _step;

        try {
          for (_iterator.s(); !(_step = _iterator.n()).done;) {
            var node = _step.value;

            if (node.state == state) {
              return node;
            }
          }
        } catch (err) {
          _iterator.e(err);
        } finally {
          _iterator.f();
        }
      }

      return null;
    }
  }, {
    key: "getNode",
    value: function getNode(id) {
      var state = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
      var nodes = this.m_map[id];

      if (nodes != null) {
        var _iterator2 = _createForOfIteratorHelper(nodes),
            _step2;

        try {
          for (_iterator2.s(); !(_step2 = _iterator2.n()).done;) {
            var node = _step2.value;

            if (node.state == state) {
              return node;
            }
          }
        } catch (err) {
          _iterator2.e(err);
        } finally {
          _iterator2.f();
        }
      }

      return this.create(id, state);
    }
  }, {
    key: "create",
    value: function create(id, state) {
      var node = new _Node["default"](this.m_nodes.length + 1);
      node.id = id;
      node.state = state;
      this.m_nodes.push(node);
      var nodes = this.m_map[id];

      if (nodes == null) {
        nodes = [];
        this.m_map[id] = nodes;
      }

      nodes.push(node);
      return node;
    }
  }, {
    key: "getNodeIdx",
    value: function getNodeIdx(node) {
      return node != null ? node.index : 0;
    }
  }, {
    key: "getNodeAtIdx",
    value: function getNodeAtIdx(idx) {
      return idx != 0 ? this.m_nodes[idx - 1] : null;
    }
  }, {
    key: "getNodeCount",
    value: function getNodeCount() {
      return this.m_nodes.length;
    } // getNode(ref) {
    // return this.getNode(ref, 0);
    // }

    /*
    
    inline let getMaxNodes() const { return m_maxNodes; }
    inline dtNodeIndex getFirst(bucket) const { return m_first[bucket]; }
    inline dtNodeIndex getNext(i) const { return m_next[i]; }
    */

  }]);

  return NodePool;
}();

var _default = NodePool;
exports["default"] = _default;