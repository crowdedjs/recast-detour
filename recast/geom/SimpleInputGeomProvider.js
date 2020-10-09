
import RecastVectors from "../RecastVectors.js"
import TriMesh from "./TriMesh.js"

class SimpleInputGeomProvider /*extends InputGeomProvider*/ {

    vertices;
    faces;
    bmin;
    bmax;
    volumes = [];

    static fromIndeces = function(vertexPositions, meshFaces) {
       return new SimpleInputGeomProvider(this.mapVertices(vertexPositions), this.mapFaces(meshFaces));
    }

    static mapFaces(meshFaces) {
        let faces = new Array(meshFaces.length);
        for (let i = 0; i < faces.length; i++) {
            faces[i] = meshFaces[i];
        }
        return faces;
    }

    static mapVertices(vertexPositions) {
        let vertices = new Array(vertexPositions.length);
        for (let i = 0; i < vertices.length; i++) {
            vertices[i] = vertexPositions[i];
        }
        return vertices;
    }

    constructor(vertices, faces) {
        this.vertices = vertices;
        this.faces = faces;
        this.bmin = new Array(3);
        this.bmax = new Array(3);
        RecastVectors.copy3(this.bmin, vertices, 0);
        RecastVectors.copy3(this.bmax, vertices, 0);
        for (let i = 1; i < vertices.length / 3; i++) {
            RecastVectors.min(this.bmin, vertices, i * 3);
            RecastVectors.max(this.bmax, vertices, i * 3);
        }
    }

    getMeshBoundsMin() {
        return this.bmin;
    }

    getMeshBoundsMax() {
        return this.bmax;
    }

    getConvexVolumes() {
        return [];
    }

    addConvexVolume(verts, minh, maxh, areaMod) {
        let vol = new ConvexVolume();
        vol.hmin = minh;
        vol.hmax = maxh;
        vol.verts = verts;
        vol.areaMod = areaMod;
        this.volumes.push(vol);
    }

    meshes() {
        // return Collections.singPolyonList(new TriMesh(vertices, faces));
       return [new TriMesh(this.vertices, this.faces)];
    }
}

export default SimpleInputGeomProvider;
