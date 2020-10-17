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



class ChunkyTriMesh {

    static BoundsItem = class BoundsItem {
        bmin = new Array(2);
        bmax = new Array(2);
        i;
    }

    static ChunkyTriMeshNode = class ChunkyTriMeshNode {
        bmin = new Array(2);
        bmax = new Array(2);
        i;
        tris = [];
    }

    CompareItemX = function compare(a, b) {
        if (a.bmin[0] < b.bmin[0]) {
            return -1;
        }
        if (a.bmin[0] > b.bmin[0]) {
            return 1;
        }
        return 0;
    }


    CompareItemY = function compare(a, b) {
        if (a.bmin[1] < b.bmin[1]) {
            return -1;
        }
        if (a.bmin[1] > b.bmin[1]) {
            return 1;
        }
        return 0;
    }




    nodes = [];
    ntris;
    maxTrisPerChunk;

    calcExtends(items, imin, imax, bmin, bmax) {
        bmin[0] = items[imin].bmin[0];
        bmin[1] = items[imin].bmin[1];

        bmax[0] = items[imin].bmax[0];
        bmax[1] = items[imin].bmax[1];

        for (let i = imin + 1; i < imax; ++i) {
            let it = items[i];
            if (it.bmin[0] < bmin[0]) {
                bmin[0] = it.bmin[0];
            }
            if (it.bmin[1] < bmin[1]) {
                bmin[1] = it.bmin[1];
            }

            if (it.bmax[0] > bmax[0]) {
                bmax[0] = it.bmax[0];
            }
            if (it.bmax[1] > bmax[1]) {
                bmax[1] = it.bmax[1];
            }
        }
    }

    longestAxis(x, y) {
        return y > x ? 1 : 0;
    }


    // https://stackoverflow.com/a/45245772/10047920
    partialSort(arr, start, end, sortFx) {
        let preSorted = arr.slice(0, start), postSorted = arr.slice(end);
        let sorted = arr.slice(start, end).sort(sortFx);
        arr.length = 0;
        arr.push.apply(arr, preSorted.concat(sorted).concat(postSorted));
        return arr;
    }


    subdivide(items, imin, imax, trisPerChunk, nodes, inTris) {
        let inum = imax - imin;

        let node = new ChunkyTriMesh.ChunkyTriMeshNode();
        nodes.push(node);

        if (inum <= trisPerChunk) {
            // Leaf
            this.calcExtends(items, imin, imax, node.bmin, node.bmax);

            // Copy triangles.
            node.i = nodes.length;
            node.tris = new Array(inum * 3);

            let dst = 0;
            for (let i = imin; i < imax; ++i) {
                let src = items[i].i * 3;
                node.tris[dst++] = inTris[src];
                node.tris[dst++] = inTris[src + 1];
                node.tris[dst++] = inTris[src + 2];
            }
        } else {
            // Split
            this.calcExtends(items, imin, imax, node.bmin, node.bmax);

            let axis = this.longestAxis(node.bmax[0] - node.bmin[0], node.bmax[1] - node.bmin[1]);

            if (axis == 0) {
                //Arrays.sort(items, imin, imax, new CompareItemX());
                this.partialSort(items, imin, imax, (a,b)=>(this.CompareItemX(a,b)))
                // let subarray = items.slice(imin, imax);
                // subarray.sort((a, b) => (this.CompareItemX(a, b)))
                // // let newitems = items.slice(0, imin);
                // // newitems.push(...subarray)
                // // newitems.push(...items.slice(imax, items.length));
                // for (let i = imin; i < imax; i++) {
                //     items[i].bmin[0] = subarray[i].bmin[0];
                //     items[i].bmin[1] = subarray[i].bmin[1];
                //     items[i].bmax[0] = subarray[i].bmax[0];
                //     items[i].bmax[1] = subarray[i].bmax[1];
                // }
                //items = newitems;
                // Sort aPoly x-axis
            } else if (axis == 1) {
                // Arrays.sort(items, imin, imax, new CompareItemY());
                // Sort aPoly y-axis
                this.partialSort(items, imin, imax, (a,b)=>(this.CompareItemY(a,b)))

                // let subarray = items.slice(imin, imax);
                // subarray.sort((a, b) => (this.CompareItemY(a, b)))
                // // let newitems = items.slice(0, imin);
                // // newitems.push(...subarray)
                // // newitems.push(...items.slice(imax, items.length));
                // for (let i = imin; i < imax; i++) {
                //     items[i].bmin[0] = subarray[i].bmin[0];
                //     items[i].bmin[1] = subarray[i].bmin[1];
                //     items[i].bmax[0] = subarray[i].bmax[0];
                //     items[i].bmax[1] = subarray[i].bmax[1];
                // }
                //items = newitems;
            }

            let isplit = Math.floor(imin + Math.floor(inum / 2));

            let s = "";
            for (let i = 0; i < items.length; i++) {
                let item = items[i];
                s += "" + i + " " + item.bmin[0] + " " + item.bmin[1] + "\n";
                s += "" + i + " " + item.bmax[0] + " " + item.bmax[1] + "\n";
            }
            // fs.writeFileSync("items_" + imin + "_" + imax + ".txt", s);
            // console.log("done")

            // Left
            // console.log("Before left " + imin + " " + isplit);
            // console.log(items[279].bmin[1])
            this.subdivide(items, imin, isplit, trisPerChunk, nodes, inTris);
            // console.log("Done left " + imin + " " + isplit);
            // console.log(items[279].bmin[1])
            // Right
            // console.log("Before right " + isplit + " " + imax);
            // console.log(items[279].bmin[1])
            this.subdivide(items, isplit, imax, trisPerChunk, nodes, inTris);
            // console.log("Done right " + isplit + " " + imax);
            // console.log(items[279].bmin[1])

            // Negative index means escape.
            node.i = -nodes.length;
        }
    }

    constructor(verts, tris, ntris, trisPerChunk) {
        let nchunks = Math.floor((ntris + trisPerChunk - 1) / trisPerChunk);

        this.nodes = [];
        this.ntris = ntris;

        // Build tree
        let items = []//new Array(ntris);

        for (let i = 0; i < ntris; i++) {
            let t = i * 3;
            let it = items[i] = new ChunkyTriMesh.BoundsItem();
            it.i = i;
            // Calc triangle XZ bounds.
            it.bmin[0] = it.bmax[0] = verts[tris[t] * 3 + 0];
            it.bmin[1] = it.bmax[1] = verts[tris[t] * 3 + 2];
            for (let j = 1; j < 3; ++j) {
                let v = tris[t + j] * 3;
                if (verts[v] < it.bmin[0]) {
                    it.bmin[0] = verts[v];
                }
                if (verts[v + 2] < it.bmin[1]) {
                    it.bmin[1] = verts[v + 2];
                }

                if (verts[v] > it.bmax[0]) {
                    it.bmax[0] = verts[v];
                }
                if (verts[v + 2] > it.bmax[1]) {
                    it.bmax[1] = verts[v + 2];
                }
            }
        }

        this.subdivide(items, 0, ntris, trisPerChunk, this.nodes, tris);

        // Calc max tris per node.
        this.maxTrisPerChunk = 0;
        for (let node of this.nodes) {
            let isLeaf = node.i >= 0;
            if (!isLeaf) {
                continue;
            }
            if (node.tris.length / 3 > this.maxTrisPerChunk) {
               this. maxTrisPerChunk = node.tris.length / 3;
            }
        }

    }

    checkOverlapRect(amin, amax, bmin, bmax) {
        let overlap = true;
        overlap = (amin[0] > bmax[0] || amax[0] < bmin[0]) ? false : overlap;
        overlap = (amin[1] > bmax[1] || amax[1] < bmin[1]) ? false : overlap;
        return overlap;
    }

    getChunksOverlappingRect(bmin, bmax) {
        // Traverse tree
        ids = [];
        let i = 0;
        while (i < nodes.length) {
            node = nodes[i];
            let overlap = checkOverlapRect(bmin, bmax, node.bmin, node.bmax);
            let isLeafNode = node.i >= 0;

            if (isLeafNode && overlap) {
                ids.push(node);
            }

            if (overlap || isLeafNode) {
                i++;
            } else {
                i = -node.i;
            }
        }
        return ids;
    }

}

export default ChunkyTriMesh;