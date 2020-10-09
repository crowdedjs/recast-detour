"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

var _PriorityQueue = _interopRequireDefault(require("./crowd/PriorityQueue.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { "default": obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function compare(one, two) {
  if (one == two) return 0;
  if (one < two) return -1;else return 1;
}

var NodeQueue = /*#__PURE__*/function () {
  function NodeQueue() {
    _classCallCheck(this, NodeQueue);

    this.m_heap = new _PriorityQueue["default"](function (n1, n2) {
      return compare(n1.total, n2.total);
    });
  }

  _createClass(NodeQueue, [{
    key: "clear",
    value: function clear() {
      this.m_heap = new _PriorityQueue["default"](function (n1, n2) {
        return compare(n1.total, n2.total);
      });
    }
  }, {
    key: "top",
    value: function top() {
      return this.m_heap.peek();
    }
  }, {
    key: "pop",
    value: function pop() {
      return this.m_heap.poll();
    }
  }, {
    key: "push",
    value: function push(node) {
      this.m_heap.insert(node);
    }
  }, {
    key: "modify",
    value: function modify(node) {
      this.m_heap.remove(node);
      this.m_heap.offer(node);
    }
  }, {
    key: "isEmpty",
    value: function isEmpty() {
      return this.m_heap.isEmpty();
    }
  }]);

  return NodeQueue;
}();

var _default = NodeQueue;
exports["default"] = _default;