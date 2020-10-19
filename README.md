# Recast Detour ported to JavaScript

This is a JavaScript implementation of the Recast Navigation [library](https://github.com/recastnavigation/recastnavigation), also called Recast Detour. This implementation is not a direct implementation, but a port of the Java [port](https://github.com/ppiastucki/recast4j).

This project is unrelated to the [recast-detour](https://www.npmjs.com/package/recast-detour) project found on npm. That project uses emscripten to create a wrapper around the original Recast Navigation code. This project is a full rewrite in JavaScript.

# Getting Started

## In the browser

The easiest way to use this project in a browser is through a CDN, 
```https://cdn.jsdelivr.net/npm/@ricksteam/recastdetourjs/recastdetourjs.js```

or 
```https://cdn.jsdelivr.net/npm/@ricksteam/recastdetourjs/recastdetourjs.min.js```

for the minified version.

## For node projects

To access this library in a node project, first install it as a dependency,

```npm i @ricksteam/recastdetourjs```

and then include it in your code,

```const recastdetourjs = require('recastdetourjs')```

## Building the project
The project is built use babel and rollup. The package.json file has helper scripts for this. To build the project from scratch, run
```npm run build```

followed by
```npm run rollup```

This will create a new version of recastdetourjs.js in the main directory.

# License information

This has the same license ([zlib](https://opensource.org/licenses/Zlib)) as both the original recast navigation [library](https://github.com/recastnavigation/recastnavigation) and the Java [port](https://github.com/ppiastucki/recast4j) on which this is based.

In addition to the core library, this repository contains an examples folder. The examples folder has copies of the [axios](https://github.com/axios/axios) library and parts of the [three.js](https://github.com/mrdoob/three.js/) library in the examples folder. They are both MIT licensed by their authors. The code in the examples folder that is unique to this project is also MIT licensed.