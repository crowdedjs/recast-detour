"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports["default"] = void 0;

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
var Node = /*#__PURE__*/function () {
  /** parent of the node is not adjacent. Found using raycast. */

  /** Position of the node. */

  /** Cost from previous node to current node. */

  /** Cost up to the node. */

  /** Index to parent node. */

  /** extra state information. A polyRef can have multiple nodes with different extra info. see DT_MAX_STATES_PER_NODE */

  /** Node flags. A combination of dtNodeFlags. */

  /** Polygon ref the node corresponds to. */
  function Node(index) {
    _classCallCheck(this, Node);

    _defineProperty(this, "index", 0);

    _defineProperty(this, "pos", new Array(3));

    _defineProperty(this, "cost", 0);

    _defineProperty(this, "total", 0);

    _defineProperty(this, "pidx", 0);

    _defineProperty(this, "state", 0);

    _defineProperty(this, "flags", 0);

    _defineProperty(this, "id", 0);

    this.index = index;
  }

  _createClass(Node, [{
    key: "toString",
    value: function toString() {
      return "Node [id=" + id + "]";
    }
  }]);

  return Node;
}();

_defineProperty(Node, "DT_NODE_OPEN", 0x01);

_defineProperty(Node, "DT_NODE_CLOSED", 0x02);

_defineProperty(Node, "DT_NODE_PARENT_DETACHED", 0x04);

var _default = Node;
exports["default"] = _default;