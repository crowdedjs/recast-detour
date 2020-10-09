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
    this.array = [null];
    this.size = 0; //this.array = [null, ...array];
    //this.size = array.length
  }

  _createClass(PriorityQueue, [{
    key: "isEmpty",
    value: function isEmpty() {
      return this.size == 0;
    }
  }, {
    key: "arrange",
    value: function arrange(idx) {
      while (idx > 1 && this.compare(Math.floor(idx / 2), idx)) {
        this.swap(idx, Math.floor(idx / 2));
        idx = Math.floor(idx / 2);
      }
    }
  }, {
    key: "heaper",
    value: function heaper(idx1) {
      while (2 * idx1 <= this.size) {
        var idx2 = 2 * idx1;
        if (idx2 < this.size && this.compare(idx2, idx2 + 1)) idx2++;
        if (!this.compare(idx1, idx2)) break;
        this.swap(idx1, idx2);
        idx1 = idx2;
      }
    }
  }, {
    key: "poll",
    value: function poll() {
      return this.rootdelete();
    }
  }, {
    key: "rootdelete",
    value: function rootdelete() {
      var max = this.array[1];
      this.swap(1, this.size--);
      this.heaper(1);
      this.array[this.size + 1] = null;
      return max;
    }
  }, {
    key: "push",
    value: function push(element) {
      this.insert(element);
    }
  }, {
    key: "insert",
    value: function insert(element) {
      this.array[++this.size] = element;
      this.arrange(this.size);
    }
  }, {
    key: "compare",
    value: function compare(idx1, idx2) {
      return this.comparator(this.array[idx1], this.array[idx2]); // return this.array[idx1].priority < this.array[idx2].priority;
    }
  }, {
    key: "swap",
    value: function swap(idx1, idx2) {
      var temp = this.array[idx1];
      this.array[idx1] = this.array[idx2];
      this.array[idx2] = temp;
    }
  }, {
    key: "returnall",
    value: function returnall() {
      return this.array;
    }
  }]);

  return PriorityQueue;
}();

var _default = PriorityQueue;
exports["default"] = _default;