# Recast Detour ported to JavaScript

This is a JavaScript implementation of the Recast Navigation [library](https://github.com/recastnavigation/recastnavigation), also called Recast Detour. This implementation is not a direct implementation, but a port of the Java [port](https://github.com/ppiastucki/recast4j).

This project is unrelated to the [recast-detour](https://www.npmjs.com/package/recast-detour) project found on npm. That project uses emscripten to create a wrapper around the original Recast Navigation code. This project is a full rewrite in JavaScript.

# License information

This has the same license ([zlib](https://opensource.org/licenses/Zlib)) as both the original recast navigation [library](https://github.com/recastnavigation/recastnavigation) and the Java [port](https://github.com/ppiastucki/recast4j) on which this is based.

In addition to the core library, this repository contains an examples folder. The examples folder has copies of the [axios](https://github.com/axios/axios) library and parts of the [three.js](https://github.com/mrdoob/three.js/) library in the examples folder. They are both MIT licensed. The code in the examples folder that is unique to the project is also MIT licensed.