"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

//https://truetocode.com/binary-treemax-heap-priority-queue-and-implementation-using-javascript/427/
var PriorityQueue = /*#__PURE__*/function () {
  function PriorityQueue(comparator) {
    _classCallCheck(this, PriorityQueue);

    this.comparator = comparator;
    this.array = [];
  }

  _createClass(PriorityQueue, [{
    key: "isEmpty",
    value: function isEmpty() {
      return this.array.length == 0;
    }
  }, {
    key: "poll",
    value: function poll() {
      return this.array.splice(0, 1)[0];
    }
  }, {
    key: "insert",
    value: function insert(element) {
      this.push(element);
    }
  }, {
    key: "push",
    value: function push(element) {
      this.array.push(element);
      this.array.sort(this.comparator);
    }
  }, {
    key: "offer",
    value: function offer(element) {
      this.push(element);
    }
  }, {
    key: "remove",
    value: function remove(element) {
      var index = this.array.indexOf(element);
      if (index >= 0) this.array.splice(index, 1);
    }
  }]);

  return PriorityQueue;
}();

var _default = PriorityQueue;
exports["default"] = _default;