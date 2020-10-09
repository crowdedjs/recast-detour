/*
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

import SimpleInputGeomProvider from "./geom/SimpleInputGeomProvider.js"

class ObjImporterContext {
    vertexPositions = [];
    meshFaces = [];
}


class ObjImporter {

    // OBJImporterContext = 

    load(is) {
        let context = new ObjImporterContext();
        let reader = null;
        try {
            let slurp = is
            let lines = slurp.split(/\r?\n/);
            for (let i = 0; i < lines.length; i++) {
                let line = lines[i];
                this.readLine(line, context);
            }
            // reader = new BufferedReader(new InputStreamReader(is));
            // let line;
            // while ((line = reader.readLine()) != null) {
            //     line = line.trim();
            //     readLine(line, context);
            // }
        } catch (e) {
            throw e;
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (e) {
                    throw new RuntimeException(e.getMessage(), e);
                }
            }
        }
        return SimpleInputGeomProvider.fromIndeces(context.vertexPositions, context.meshFaces);

    }

    readLine(line, context) {
        if (line.startsWith("v")) {
            this.readVertex(line, context);
        } else if (line.startsWith("f")) {
            this.readFace(line, context);
        }
    }

    readVertex(line, context) {
        if (line.startsWith("v ")) {
            let vert = this.readVector3f(line);
            for (let vp of vert) {
                context.vertexPositions.push(vp);
            }
        }
    }

    readVector3f(line) {
        let v = line.split(/\s+/);
        if (v.length < 4) {
            throw new RuntimeException("Invalid vector, expected 3 coordinates, found " + (v.length - 1));
        }
        return [parseFloat(v[1]), parseFloat(v[2]), parseFloat(v[3])];
    }

    readFace(line, context) {
        let v = line.split(/\s+/);
        if (v.length < 4) {
            throw new RuntimeException("Invalid number of face vertices: 3 coordinates expected, found "
                + v.length);
        }
        for (let j = 0; j < v.length - 3; j++) {
            context.meshFaces.push(this.readFaceVertex(v[1], context));
            for (let i = 0; i < 2; i++) {
                context.meshFaces.push(this.readFaceVertex(v[2 + j + i], context));
            }
        }
    }

    readFaceVertex(face, context) {
        let v = face.split(/\//);
        return this.getIndex(parseInt(v[0]), context.vertexPositions.length);
    }

    getIndex(posi, size) {
        if (posi > 0) {
            posi--;
        } else if (posi < 0) {
            posi = size + posi;
        } else {
            throw new RuntimeException("0 vertex index");
        }
        return posi;
    }

}

export default ObjImporter;