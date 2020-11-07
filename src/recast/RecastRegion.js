/*
Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
Recast4J Copyright (c) 2015-2018 Piotr Piastucki piotr@jtilia.org

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

import RecastConstants from "./RecastConstants.js"
import RecastCommon from "./RecastCommon.js"

class RecastRegion {

    static RC_NULL_NEI = 0xfff;

    static SweepSpan = class SweepSpan {
        rid; // row id
        id; // region id
        ns; // number samples
        nei; // neighbour id
    }

    static calculateDistanceField(chf, src) {
        let maxDist;
        let w = chf.width;
        let h = chf.height;

        // Init distance and points.
        for (let i = 0; i < chf.spanCount; ++i) {
            src[i] = 0xffff;
        }

        // Mark boundary cells.
        for (let y = 0; y < h; ++y) {
            for (let x = 0; x < w; ++x) {
                let c = chf.cells[x + y * w];
                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
                    let s = chf.spans[i];
                    let area = chf.areas[i];

                    let nc = 0;
                    for (let dir = 0; dir < 4; ++dir) {
                        if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
                            let ax = x + RecastCommon.GetDirOffsetX(dir);
                            let ay = y + RecastCommon.GetDirOffsetY(dir);
                            let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);
                            if (area == chf.areas[ai]) {
                                nc++;
                            }
                        }
                    }
                    if (nc != 4) {
                        src[i] = 0;
                    }
                }
            }
        }

        // Pass 1
        for (let y = 0; y < h; ++y) {
            for (let x = 0; x < w; ++x) {
                let c = chf.cells[x + y * w];
                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
                    let s = chf.spans[i];

                    if (RecastCommon.GetCon(s, 0) != RecastConstants.RC_NOT_CONNECTED) {
                        // (-1,0)
                        let ax = x + RecastCommon.GetDirOffsetX(0);
                        let ay = y + RecastCommon.GetDirOffsetY(0);
                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 0);
                        let as = chf.spans[ai];
                        if (src[ai] + 2 < src[i]) {
                            src[i] = src[ai] + 2;
                        }

                        // (-1,-1)
                        if (RecastCommon.GetCon(as, 3) != RecastConstants.RC_NOT_CONNECTED) {
                            let aax = ax + RecastCommon.GetDirOffsetX(3);
                            let aay = ay + RecastCommon.GetDirOffsetY(3);
                            let aai = chf.cells[aax + aay * w].index + RecastCommon.GetCon(as, 3);
                            if (src[aai] + 3 < src[i]) {
                                src[i] = src[aai] + 3;
                            }
                        }
                    }
                    if (RecastCommon.GetCon(s, 3) != RecastConstants.RC_NOT_CONNECTED) {
                        // (0,-1)
                        let ax = x + RecastCommon.GetDirOffsetX(3);
                        let ay = y + RecastCommon.GetDirOffsetY(3);
                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 3);
                        let as = chf.spans[ai];
                        if (src[ai] + 2 < src[i]) {
                            src[i] = src[ai] + 2;
                        }

                        // (1,-1)
                        if (RecastCommon.GetCon(as, 2) != RecastConstants.RC_NOT_CONNECTED) {
                            let aax = ax + RecastCommon.GetDirOffsetX(2);
                            let aay = ay + RecastCommon.GetDirOffsetY(2);
                            let aai = chf.cells[aax + aay * w].index + RecastCommon.GetCon(as, 2);
                            if (src[aai] + 3 < src[i]) {
                                src[i] = src[aai] + 3;
                            }
                        }
                    }
                }
            }
        }

        // Pass 2
        for (let y = h - 1; y >= 0; --y) {
            for (let x = w - 1; x >= 0; --x) {
                let c = chf.cells[x + y * w];
                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
                    let s = chf.spans[i];

                    if (RecastCommon.GetCon(s, 2) != RecastConstants.RC_NOT_CONNECTED) {
                        // (1,0)
                        let ax = x + RecastCommon.GetDirOffsetX(2);
                        let ay = y + RecastCommon.GetDirOffsetY(2);
                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 2);
                        let as = chf.spans[ai];
                        if (src[ai] + 2 < src[i]) {
                            src[i] = src[ai] + 2;
                        }

                        // (1,1)
                        if (RecastCommon.GetCon(as, 1) != RecastConstants.RC_NOT_CONNECTED) {
                            let aax = ax + RecastCommon.GetDirOffsetX(1);
                            let aay = ay + RecastCommon.GetDirOffsetY(1);
                            let aai = chf.cells[aax + aay * w].index + RecastCommon.GetCon(as, 1);
                            if (src[aai] + 3 < src[i]) {
                                src[i] = src[aai] + 3;
                            }
                        }
                    }
                    if (RecastCommon.GetCon(s, 1) != RecastConstants.RC_NOT_CONNECTED) {
                        // (0,1)
                        let ax = x + RecastCommon.GetDirOffsetX(1);
                        let ay = y + RecastCommon.GetDirOffsetY(1);
                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 1);
                        let as = chf.spans[ai];
                        if (src[ai] + 2 < src[i]) {
                            src[i] = src[ai] + 2;
                        }

                        // (-1,1)
                        if (RecastCommon.GetCon(as, 0) != RecastConstants.RC_NOT_CONNECTED) {
                            let aax = ax + RecastCommon.GetDirOffsetX(0);
                            let aay = ay + RecastCommon.GetDirOffsetY(0);
                            let aai = chf.cells[aax + aay * w].index + RecastCommon.GetCon(as, 0);
                            if (src[aai] + 3 < src[i]) {
                                src[i] = src[aai] + 3;
                            }
                        }
                    }
                }
            }
        }

        maxDist = 0;
        for (let i = 0; i < chf.spanCount; ++i) {
            maxDist = Math.max(src[i], maxDist);
        }

        return maxDist;
    }

    static boxBlur(chf, thr, src) {
        let w = chf.width;
        let h = chf.height;
        let dst = new Array(chf.spanCount);

        thr *= 2;

        for (let y = 0; y < h; ++y) {
            for (let x = 0; x < w; ++x) {
                let c = chf.cells[x + y * w];
                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
                    let s = chf.spans[i];
                    let cd = src[i];
                    if (cd <= thr) {
                        dst[i] = cd;
                        continue;
                    }

                    let d = cd;
                    for (let dir = 0; dir < 4; ++dir) {
                        if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
                            let ax = x + RecastCommon.GetDirOffsetX(dir);
                            let ay = y + RecastCommon.GetDirOffsetY(dir);
                            let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);
                            d += src[ai];

                            let as = chf.spans[ai];
                            let dir2 = (dir + 1) & 0x3;
                            if (RecastCommon.GetCon(as, dir2) != RecastConstants.RC_NOT_CONNECTED) {
                                let ax2 = ax + RecastCommon.GetDirOffsetX(dir2);
                                let ay2 = ay + RecastCommon.GetDirOffsetY(dir2);
                                let ai2 = chf.cells[ax2 + ay2 * w].index + RecastCommon.GetCon(as, dir2);
                                d += src[ai2];
                            } else {
                                d += cd;
                            }
                        } else {
                            d += cd * 2;
                        }
                    }
                    dst[i] = Math.floor((d + 5) / 9);
                }
            }
        }
        return dst;
    }

    static floodRegion(x, y, i, level, r, chf, srcReg,
        srcDist, stack) {
        let w = chf.width;

        let area = chf.areas[i];

        // Flood fill mark region.
        stack = [];
        stack.push(x);
        stack.push(y);
        stack.push(i);
        srcReg[i] = r;
        srcDist[i] = 0;

        let lev = level >= 2 ? level - 2 : 0;
        let count = 0;

        while (stack.length > 0) {
            let ci = stack.splice(stack.length - 1,1)[0];
            let cy = stack.splice(stack.length - 1,1)[0];
            let cx = stack.splice(stack.length - 1,1)[0];

            let cs = chf.spans[ci];

            // Check if any of the neighbours already have a valid region set.
            let ar = 0;
            for (let dir = 0; dir < 4; ++dir) {
                // 8 connected
                if (RecastCommon.GetCon(cs, dir) != RecastConstants.RC_NOT_CONNECTED) {
                    let ax = cx + RecastCommon.GetDirOffsetX(dir);
                    let ay = cy + RecastCommon.GetDirOffsetY(dir);
                    let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(cs, dir);
                    if (chf.areas[ai] != area) {
                        continue;
                    }
                    let nr = srcReg[ai];
                    if ((nr & RecastConstants.RC_BORDER_REG) != 0) {
                        continue;
                    }
                    if (nr != 0 && nr != r) {
                        ar = nr;
                        break;
                    }

                    let as = chf.spans[ai];

                    let dir2 = (dir + 1) & 0x3;
                    if (RecastCommon.GetCon(as, dir2) != RecastConstants.RC_NOT_CONNECTED) {
                        let ax2 = ax + RecastCommon.GetDirOffsetX(dir2);
                        let ay2 = ay + RecastCommon.GetDirOffsetY(dir2);
                        let ai2 = chf.cells[ax2 + ay2 * w].index + RecastCommon.GetCon(as, dir2);
                        if (chf.areas[ai2] != area) {
                            continue;
                        }
                        let nr2 = srcReg[ai2];
                        if (nr2 != 0 && nr2 != r) {
                            ar = nr2;
                            break;
                        }
                    }
                }
            }
            if (ar != 0) {
                srcReg[ci] = 0;
                continue;
            }

            count++;

            // Expand neighbours.
            for (let dir = 0; dir < 4; ++dir) {
                if (RecastCommon.GetCon(cs, dir) != RecastConstants.RC_NOT_CONNECTED) {
                    let ax = cx + RecastCommon.GetDirOffsetX(dir);
                    let ay = cy + RecastCommon.GetDirOffsetY(dir);
                    let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(cs, dir);
                    if (chf.areas[ai] != area) {
                        continue;
                    }
                    if (chf.dist[ai] >= lev && srcReg[ai] == 0) {
                        srcReg[ai] = r;
                        srcDist[ai] = 0;
                        stack.push(ax);
                        stack.push(ay);
                        stack.push(ai);
                    }
                }
            }
        }

        return count > 0;
    }

    static expandRegions(maxIter, level, chf, srcReg, srcDist,
        stack, fillStack) {
        let w = chf.width;
        let h = chf.height;

        if (fillStack) {
            // Find cells revealed by the raised level.
            //stack = [];
            stack = []
            for (let y = 0; y < h; ++y) {
                for (let x = 0; x < w; ++x) {
                    let c = chf.cells[x + y * w];
                    for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
                        if (chf.dist[i] >= level && srcReg[i] == 0 && chf.areas[i] != RecastConstants.RC_NULL_AREA) {
                            stack.push(x);
                            stack.push(y);
                            stack.push(i);
                        }
                    }
                }
            }
        } else // use cells in the input stack
        {
            // mark all cells which already have a region
            for (let j = 0; j < stack.length; j += 3) {
                let i = stack[j + 2];
                if (srcReg[i] != 0) {
                    stack[j + 2] = -1;
                }
            }
        }

        let dirtyEntries = [];
        let iter = 0;
        while (stack.length > 0) {
            let failed = 0;
            // dirtyEntries = [];
            dirtyEntries = [];

            for (let j = 0; j < stack.length; j += 3) {
                let x = stack[j + 0];
                let y = stack[j + 1];
                let i = stack[j + 2];
                if (i < 0) {
                    failed++;
                    continue;
                }

                let r = srcReg[i];
                let d2 = 0xfff;
                let area = chf.areas[i];
                let s = chf.spans[i];
                for (let dir = 0; dir < 4; ++dir) {
                    if (RecastCommon.GetCon(s, dir) == RecastConstants.RC_NOT_CONNECTED) {
                        continue;
                    }
                    let ax = x + RecastCommon.GetDirOffsetX(dir);
                    let ay = y + RecastCommon.GetDirOffsetY(dir);
                    let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);
                    if (chf.areas[ai] != area) {
                        continue;
                    }
                    if (srcReg[ai] > 0 && (srcReg[ai] & RecastConstants.RC_BORDER_REG) == 0) {
                        if (srcDist[ai] + 2 < d2) {
                            r = srcReg[ai];
                            d2 = srcDist[ai] + 2;
                        }
                    }
                }
                if (r != 0) {
                    stack[j + 2]= -1; // mark as used
                    dirtyEntries.push(i);
                    dirtyEntries.push(r);
                    dirtyEntries.push(d2);
                } else {
                    failed++;
                }
            }

            // Copy entries that differ between src and dst to keep them in sync.
            for (let i = 0; i < dirtyEntries.length; i += 3) {
                let idx = dirtyEntries[i];
                // if (idx == 1344)
                //     console.log("dirty")
                srcReg[idx] = dirtyEntries[i + 1];
                srcDist[idx] = dirtyEntries[i + 2];
            }

            if (failed * 3 == stack.length) {
                break;
            }

            if (level > 0) {
                ++iter;
                if (iter >= maxIter) {
                    break;
                }
            }
        }

        return srcReg;
    }

    static sortCellsByLevel(startLevel, chf, srcReg, nbStacks,
        stacks, loglevelsPerStack) // the levels per stack (2 in our case) as a bit shift
    {
        let w = chf.width;
        let h = chf.height;
        startLevel = startLevel >> loglevelsPerStack;

        for (let j = 0; j < nbStacks; ++j) {
            // stacks[j] = new Array(1024);
            // stacks[j] = [];
            stacks[j] = []
        }
        ;

        // put all cells in the level range into the appropriate stacks
        for (let y = 0; y < h; ++y) {
            for (let x = 0; x < w; ++x) {
                let c = chf.cells[x + y * w];
                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
                    if (chf.areas[i] == RecastConstants.RC_NULL_AREA || srcReg[i] != 0) {
                        continue;
                    }

                    let level = chf.dist[i] >> loglevelsPerStack;
                    let sId = startLevel - level;
                    if (sId >= nbStacks) {
                        continue;
                    }
                    if (sId < 0) {
                        sId = 0;
                    }

                    stacks[sId].push(x);
                    stacks[sId].push(y);
                    stacks[sId].push(i);
                }
            }
        }
    }

    static appendStacks(srcStack, dstStack, srcReg) {
        for (let j = 0; j < srcStack.length; j += 3) {
            let i = srcStack[j + 2];
            if ((i < 0) || (srcReg[i] != 0)) {
                continue;
            }
            dstStack.push(srcStack[j]);
            dstStack.push(srcStack[j + 1]);
            dstStack.push(srcStack[j + 2]);
        }
    }

    static Region = class Region {
        spanCount = 0; // Number of spans belonging to this region
        id = 0; // ID of the region
        areaType = 0; // Are type.
        remap = false;
        visited = false
        overlap = false;
        connectsToBorder = false;
        ymin = 0;
        ymax = 0;
        connections = [];
        floors = [];

        constructor(i) {
            this.id = i;
            this.ymin = 0xFFFF;
            this.connections = [];
            this.floors = [];
        }

    };

    static removeAdjacentNeighbours(reg) {
        // Remove adjacent duplicates.
        for (let i = 0; i < reg.connections.length && reg.connections.length > 1;) {
            let ni = (i + 1) % reg.connections.length;
            if (reg.connections[i] == reg.connections[ni] && reg.connections[i] >=-128 && reg.connections[i] <= 127) {
                reg.connections.splice(i,1);
            } else {
                ++i;
            }
        }
    }

    static replaceNeighbour(reg, oldId, newId) {
        let neiChanged = false;
        for (let i = 0; i < reg.connections.length; ++i) {
            if (reg.connections[i] == oldId) {
                reg.connections[i] = newId;
                neiChanged = true;
            }
        }
        for (let i = 0; i < reg.floors.length; ++i) {
            if (reg.floors[i] == oldId) {
                reg.floors.set(i, newId);
            }
        }
        if (neiChanged) {
            RecastRegion.removeAdjacentNeighbours(reg);
        }
    }

    static canMergeWithRegion(rega, regb) {
        if (rega.areaType != regb.areaType) {
            return false;
        }
        let n = 0;
        for (let i = 0; i < rega.connections.length; ++i) {
            if (rega.connections[i] == regb.id) {
                n++;
            }
        }
        if (n > 1) {
            return false;
        }
        for (let i = 0; i < rega.floors.length; ++i) {
            if (rega.floors[i] == regb.id) {
                return false;
            }
        }
        return true;
    }

    static addUniqueFloorRegion(reg, n) {
        if (!reg.floors.includes(n)) {
            reg.floors.push(n);
        }
    }

    static mergeRegions(rega, regb) {
        let aid = rega.id;
        let bid = regb.id;

        // Duplicate current neighbourhood.
        let acon = rega.connections;
        let bcon = regb.connections;

        // Find insertion poPoly on A.
        let insa = -1;
        for (let i = 0; i < acon.length; ++i) {
            if (acon[i] == bid) {
                insa = i;
                break;
            }
        }
        if (insa == -1) {
            return false;
        }

        // Find insertion poPoly on B.
        let insb = -1;
        for (let i = 0; i < bcon.length; ++i) {
            if (bcon[i] == aid) {
                insb = i;
                break;
            }
        }
        if (insb == -1) {
            return false;
        }

        // Merge neighbours.
        rega.connections = [];
        for (let i = 0, ni = acon.length; i < ni - 1; ++i) {
            rega.connections.push(acon[(insa + 1 + i) % ni]);
        }

        for (let i = 0, ni = bcon.length; i < ni - 1; ++i) {
            rega.connections.push(bcon[(insb + 1 + i) % ni]);
        }

        RecastRegion.removeAdjacentNeighbours(rega);

        for (let j = 0; j < regb.floors.length; ++j) {
            RecastRegion.addUniqueFloorRegion(rega, regb.floors[j]);
        }
        rega.spanCount += regb.spanCount;
        regb.spanCount = 0;
        regb.connections = [];

        return true;
    }

    static isRegionConnectedToBorder(reg) {
        // Region is connected to border if
        // one of the neighbours is null id.
        return reg.connections.includes(0);
    }

    static isSolidEdge(chf, srcReg, x, y, i, dir) {
        let s = chf.spans[i];
        let r = 0;
        if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
            let ax = x + RecastCommon.GetDirOffsetX(dir);
            let ay = y + RecastCommon.GetDirOffsetY(dir);
            let ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(s, dir);
            r = srcReg[ai];
        }
        if (r == srcReg[i]) {
            return false;
        }
        return true;
    }

    static walkContour(x, y, i, dir, chf, srcReg,
        cont) {
        let startDir = dir;
        let starti = i;

        let ss = chf.spans[i];
        let curReg = 0;
        if (RecastCommon.GetCon(ss, dir) != RecastConstants.RC_NOT_CONNECTED) {
            let ax = x + RecastCommon.GetDirOffsetX(dir);
            let ay = y + RecastCommon.GetDirOffsetY(dir);
            let ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(ss, dir);
            curReg = srcReg[ai];
        }
        cont.push(curReg);

        let iter = 0;
        while (++iter < 40000) {
            let s = chf.spans[i];

            if (RecastRegion.isSolidEdge(chf, srcReg, x, y, i, dir)) {
                // Choose the edge corner
                let r = 0;
                if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
                    let ax = x + RecastCommon.GetDirOffsetX(dir);
                    let ay = y + RecastCommon.GetDirOffsetY(dir);
                    let ai = chf.cells[ax + ay * chf.width].index + RecastCommon.GetCon(s, dir);
                    r = srcReg[ai];
                }
                if (r != curReg) {
                    curReg = r;
                    cont.push(curReg);
                }

                dir = (dir + 1) & 0x3; // Rotate CW
            } else {
                let ni = -1;
                let nx = x + RecastCommon.GetDirOffsetX(dir);
                let ny = y + RecastCommon.GetDirOffsetY(dir);
                if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
                    let nc = chf.cells[nx + ny * chf.width];
                    ni = nc.index + RecastCommon.GetCon(s, dir);
                }
                if (ni == -1) {
                    // Should not happen.
                    return;
                }
                x = nx;
                y = ny;
                i = ni;
                dir = (dir + 3) & 0x3; // Rotate CCW
            }

            if (starti == i && startDir == dir) {
                break;
            }
        }

        // Remove adjacent duplicates.
        if (cont.length > 1) {
            for (let j = 0; j < cont.length;) {
                let nj = (j + 1) % cont.length;
                if (cont[j] == cont[nj] && cont[j]>=-128 && cont[j] <=127) {
                    cont.splice(j,1);
                } else {
                    ++j;
                }
            }
        }
    }

    static mergeAndFilterRegions(ctx, minRegionArea, mergeRegionSize, maxRegionId,
        chf, srcReg, overlaps) {
        let w = chf.width;
        let h = chf.height;

        let nreg = maxRegionId + 1;
        let regions = new Array(nreg);

        // Construct regions
        for (let i = 0; i < nreg; ++i) {
            regions[i] = new RecastRegion.Region(i);
        }

        // Find edge of a region and find connections around the contour.
        for (let y = 0; y < h; ++y) {
            for (let x = 0; x < w; ++x) {
                let c = chf.cells[x + y * w];
                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
                    let r = srcReg[i];
                    if (r == 0 || r >= nreg) {
                        continue;
                    }

                    let reg = regions[r];
                    reg.spanCount++;

                    // Update floors.
                    for (let j = c.index; j < ni; ++j) {
                        if (i == j) {
                            continue;
                        }
                        let floorId = srcReg[j];
                        if (floorId == 0 || floorId >= nreg) {
                            continue;
                        }
                        if (floorId == r) {
                            reg.overlap = true;
                        }
                        RecastRegion.addUniqueFloorRegion(reg, floorId);
                    }

                    // Have found contour
                    if (reg.connections.length > 0) {
                        continue;
                    }

                    reg.areaType = chf.areas[i];

                    // Check if this cell is next to a border.
                    let ndir = -1;
                    for (let dir = 0; dir < 4; ++dir) {
                        if (RecastRegion.isSolidEdge(chf, srcReg, x, y, i, dir)) {
                            ndir = dir;
                            break;
                        }
                    }

                    if (ndir != -1) {
                        // The cell is at border.
                        // Walk around the contour to find all the neighbours.
                        RecastRegion.walkContour(x, y, i, ndir, chf, srcReg, reg.connections);
                    }
                }
            }
        }

        // Remove too small regions.
        let stack = new Array(32);
        let trace = new Array(32);
        for (let i = 0; i < nreg; ++i) {
            let reg = regions[i];
            if (reg.id == 0 || (reg.id & RecastConstants.RC_BORDER_REG) != 0) {
                continue;
            }
            if (reg.spanCount == 0) {
                continue;
            }
            if (reg.visited) {
                continue;
            }

            // Count the total size of all the connected regions.
            // Also keep track of the regions connects to a tile border.
            let connectsToBorder = false;
            let spanCount = 0;
            //stack = [];
            //trace = [];
            stack = [];
            trace = [];

            reg.visited = true;
            stack.push(i);

            while (stack.length > 0) {
                // Pop
                let ri = stack.splice(stack.length - 1, 1);

                let creg = regions[ri];

                spanCount += creg.spanCount;
                trace.push(ri);

                for (let j = 0; j < creg.connections.length; ++j) {
                    if ((creg.connections[j] & RecastConstants.RC_BORDER_REG) != 0) {
                        connectsToBorder = true;
                        continue;
                    }
                    let neireg = regions[creg.connections[j]];
                    if (neireg.visited) {
                        continue;
                    }
                    if (neireg.id == 0 || (neireg.id & RecastConstants.RC_BORDER_REG) != 0) {
                        continue;
                    }
                    // Visit
                    stack.push(neireg.id);
                    neireg.visited = true;
                }
            }

            // If the accumulated regions size is too small, remove it.
            // Do not remove areas which connect to tile borders
            // as their size cannot be estimated correctly and removing them
            // can potentially remove necessary areas.
            if (spanCount < minRegionArea && !connectsToBorder) {
                // Kill all visited regions.
                for (let j = 0; j < trace.length; ++j) {
                    regions[trace[j]].spanCount = 0;
                    regions[trace[j]].id = 0;
                }
            }
        }

        // Merge too small regions to neighbour regions.
        let mergeCount = 0;
        do {
            mergeCount = 0;
            for (let i = 0; i < nreg; ++i) {
                let reg = regions[i];
                if (reg.id == 0 || (reg.id & RecastConstants.RC_BORDER_REG) != 0) {
                    continue;
                }
                if (reg.overlap) {
                    continue;
                }
                if (reg.spanCount == 0) {
                    continue;
                }

                // Check to see if the region should be merged.
                if (reg.spanCount > mergeRegionSize && RecastRegion.isRegionConnectedToBorder(reg)) {
                    continue;
                }

                // Small region with more than 1 connection.
                // Or region which is not connected to a border at all.
                // Find smallest neighbour region that connects to this one.
                let smallest = 0xfffffff;
                let mergeId = reg.id;
                for (let j = 0; j < reg.connections.length; ++j) {
                    if ((reg.connections[j] & RecastConstants.RC_BORDER_REG) != 0) {
                        continue;
                    }
                    let mreg = regions[reg.connections[j]];
                    if (mreg.id == 0 || (mreg.id & RecastConstants.RC_BORDER_REG) != 0 || mreg.overlap) {
                        continue;
                    }
                    if (mreg.spanCount < smallest && RecastRegion.canMergeWithRegion(reg, mreg) && RecastRegion.canMergeWithRegion(mreg, reg)) {
                        smallest = mreg.spanCount;
                        mergeId = mreg.id;
                    }
                }
                // Found new id.
                if (mergeId != reg.id) {
                    let oldId = reg.id;
                    let target = regions[mergeId];

                    // Merge neighbours.
                    if (RecastRegion.mergeRegions(target, reg)) {
                        // Fixup regions pointing to current region.
                        for (let j = 0; j < nreg; ++j) {
                            if (regions[j].id == 0 || (regions[j].id & RecastConstants.RC_BORDER_REG) != 0) {
                                continue;
                            }
                            // If another region was already merged into current region
                            // change the nid of the previous region too.
                            if (regions[j].id == oldId) {
                                regions[j].id = mergeId;
                            }
                            // Replace the current region with the new one if the
                            // current regions is neighbour.
                            RecastRegion.replaceNeighbour(regions[j], oldId, mergeId);
                        }
                        mergeCount++;
                    }
                }
            }
        } while (mergeCount > 0);

        // Compress region Ids.
        for (let i = 0; i < nreg; ++i) {
            regions[i].remap = false;
            if (regions[i].id == 0) {
                continue; // Skip nil regions.
            }
            if ((regions[i].id & RecastConstants.RC_BORDER_REG) != 0) {
                continue; // Skip external regions.
            }
            regions[i].remap = true;
        }

        let regIdGen = 0;
        for (let i = 0; i < nreg; ++i) {
            if (!regions[i].remap) {
                continue;
            }
            let oldId = regions[i].id;
            let newId = ++regIdGen;
            for (let j = i; j < nreg; ++j) {
                if (regions[j].id == oldId) {
                    regions[j].id = newId;
                    regions[j].remap = false;
                }
            }
        }
        maxRegionId = regIdGen;

        // Remap regions.
        for (let i = 0; i < chf.spanCount; ++i) {
            if ((srcReg[i] & RecastConstants.RC_BORDER_REG) == 0) {
                srcReg[i] = regions[srcReg[i]].id;
            }
        }

        // Return regions that we found to be overlapping.
        for (let i = 0; i < nreg; ++i) {
            if (regions[i].overlap) {
                overlaps.push(regions[i].id);
            }
        }

        return maxRegionId;
    }

    static addUniqueConnection(reg, n) {
        if (!reg.connections.contains(n)) {
            reg.connections.push(n);
        }
    }

    static mergeAndFilterLayerRegions(ctx, minRegionArea, maxRegionId,
        chf, srcReg, overlaps) {
        let w = chf.width;
        let h = chf.height;

        let nreg = maxRegionId + 1;
        let regions = new Array(nreg);

        // Construct regions
        for (let i = 0; i < nreg; ++i) {
            regions[i] = new Region(i);
        }

        // Find region neighbours and overlapping regions.
        let lregs = new Array(32);
        for (let y = 0; y < h; ++y) {
            for (let x = 0; x < w; ++x) {
                let c = chf.cells[x + y * w];

                lregs = [];

                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
                    let s = chf.spans[i];
                    let ri = srcReg[i];
                    if (ri == 0 || ri >= nreg) {
                        continue;
                    }
                    let reg = regions[ri];

                    reg.spanCount++;

                    reg.ymin = Math.min(reg.ymin, s.y);
                    reg.ymax = Math.max(reg.ymax, s.y);
                    // Collect all region layers.
                    lregs.push(ri);

                    // Update neighbours
                    for (let dir = 0; dir < 4; ++dir) {
                        if (RecastCommon.GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
                            let ax = x + RecastCommon.GetDirOffsetX(dir);
                            let ay = y + RecastCommon.GetDirOffsetY(dir);
                            let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, dir);
                            let rai = srcReg[ai];
                            if (rai > 0 && rai < nreg && rai != ri) {
                                addUniqueConnection(reg, rai);
                            }
                            if ((rai & RecastConstants.RC_BORDER_REG) != 0) {
                                reg.connectsToBorder = true;
                            }
                        }
                    }

                }

                // Update overlapping regions.
                for (let i = 0; i < lregs.length - 1; ++i) {
                    for (let j = i + 1; j < lregs.length; ++j) {
                        if (lregs[i] != lregs[j]) {
                            let ri = regions[lregs[i]];
                            let rj = regions[lregs[j]];
                            RecastRegion.addUniqueFloorRegion(ri, lregs[j]);
                            RecastRegion.addUniqueFloorRegion(rj, lregs[i]);
                        }
                    }
                }

            }
        }

        // Create 2D layers from regions.
        let layerId = 1;

        for (let i = 0; i < nreg; ++i) {
            regions[i].id = 0;
        }

        // Merge montone regions to create non-overlapping areas.
        let stack = new Array(32);
        for (let i = 1; i < nreg; ++i) {
            let root = regions[i];
            // Skip already visited.
            if (root.id != 0) {
                continue;
            }

            // Start search.
            root.id = layerId;

            stack = [];
            stack.push(i);

            while (stack.length > 0) {
                // Pop front
                let reg = regions[stack.remove(0)];

                let ncons = reg.connections.length;
                for (let j = 0; j < ncons; ++j) {
                    let nei = reg.connections[j];
                    let regn = regions[nei];
                    // Skip already visited.
                    if (regn.id != 0) {
                        continue;
                    }
                    // Skip if the neighbour is overlapping root region.
                    let overlap = false;
                    for(let k = 0; k < root.floors.length; k++) {
                        if (root.floors[k] == nei) {
                            overlap = true;
                            break;
                        }
                    }
                    if (overlap) {
                        continue;
                    }

                    // Deepen
                    stack.push(nei);

                    // Mark layer id
                    regn.id = layerId;
                    // Merge current layers to root.
                    for(let k = 0; k < regn.floors.length; ++k) {
                        RecastRegion.addUniqueFloorRegion(root, regn.floors[k]);
                    }
                    root.ymin = Math.min(root.ymin, regn.ymin);
                    root.ymax = Math.max(root.ymax, regn.ymax);
                    root.spanCount += regn.spanCount;
                    regn.spanCount = 0;
                    root.connectsToBorder = root.connectsToBorder || regn.connectsToBorder;
                }
            }

            layerId++;
        }

        // Remove small regions
        for (let i = 0; i < nreg; ++i) {
            if (regions[i].spanCount > 0 && regions[i].spanCount < minRegionArea && !regions[i].connectsToBorder) {
                let reg = regions[i].id;
                for (let j = 0; j < nreg; ++j) {
                    if (regions[j].id == reg) {
                        regions[j].id = 0;
                    }
                }
            }
        }

        // Compress region Ids.
        for (let i = 0; i < nreg; ++i) {
            regions[i].remap = false;
            if (regions[i].id == 0) {
                continue; // Skip nil regions.
            }
            if ((regions[i].id & RecastConstants.RC_BORDER_REG) != 0) {
                continue; // Skip external regions.
            }
            regions[i].remap = true;
        }

        let regIdGen = 0;
        for (let i = 0; i < nreg; ++i) {
            if (!regions[i].remap) {
                continue;
            }
            let oldId = regions[i].id;
            let newId = ++regIdGen;
            for (let j = i; j < nreg; ++j) {
                if (regions[j].id == oldId) {
                    regions[j].id = newId;
                    regions[j].remap = false;
                }
            }
        }
        maxRegionId = regIdGen;

        // Remap regions.
        for (let i = 0; i < chf.spanCount; ++i) {
            if ((srcReg[i] & RecastConstants.RC_BORDER_REG) == 0) {
                srcReg[i] = regions[srcReg[i]].id;
            }
        }

        return maxRegionId;
    }

    /// @par
    ///
    /// This is usually the second to the last step in creating a fully built
    /// compact heightfield. This step is required before regions are built
    /// using #rcBuildRegions or #rcBuildRegionsMonotone.
    ///
    /// After this step, the distance data is available via the rcCompactHeightfield::maxDistance
    /// and rcCompactHeightfield::dist fields.
    ///
    /// @see rcCompactHeightfield, rcBuildRegions, rcBuildRegionsMonotone
    static buildDistanceField(ctx, chf) {

        ctx.startTimer("BUILD_DISTANCEFIELD");
        let src = new Array(chf.spanCount);
        ctx.startTimer("DISTANCEFIELD_DIST");

        let maxDist = this.calculateDistanceField(chf, src);
        chf.maxDistance = maxDist;

        ctx.stopTimer("DISTANCEFIELD_DIST");

        ctx.startTimer("DISTANCEFIELD_BLUR");

        // Blur
        src = this.boxBlur(chf, 1, src);

        // Store distance.
        chf.dist = src;

        ctx.stopTimer("DISTANCEFIELD_BLUR");

        ctx.stopTimer("BUILD_DISTANCEFIELD");

    }

    static paintRectRegion(minx, maxx, miny, maxy, regId, chf,
        srcReg) {
        let w = chf.width;
        for (let y = miny; y < maxy; ++y) {
            for (let y = minx; x < maxx; ++x) {
                let c = chf.cells[x + y * w];
                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
                    if (chf.areas[i] != RecastConstants.RC_NULL_AREA) {
                        srcReg[i] = regId;
                    }
                }
            }
        }
    }

    /// @par
    ///
    /// Non-null regions will consist of connected, non-overlapping walkable spans that form a single contour.
    /// Contours will form simple polygons.
    ///
    /// If multiple regions form an area that is smaller than @p minRegionArea, then all spans will be
    /// re-assigned to the zero (null) region.
    ///
    /// Partitioning can result in smaller than necessary regions. @p mergeRegionArea helps
    /// reduce unecessarily small regions.
    ///
    /// See the #rcConfig documentation for more information on the configuration parameters.
    ///
    /// The region data will be available via the rcCompactHeightfield::maxRegions
    /// and rcCompactSpan::reg fields.
    ///
    /// @warning The distance field must be created using #rcBuildDistanceField before attempting to build regions.
    ///
    /// @see rcCompactHeightfield, rcCompactSpan, rcBuildDistanceField, rcBuildRegionsMonotone, rcConfig
    static buildRegionsMonotone(ctx, chf, borderSize, minRegionArea,
        mergeRegionArea) {
        ctx.startTimer("BUILD_REGIONS");

        let w = chf.width;
        let h = chf.height;
        let id = 1;

        let srcReg = new Array(chf.spanCount);

        let nsweeps = Math.max(chf.width, chf.height);
        sweeps = new Array(nsweeps);
        for (let i = 0; i < sweeps.length; i++) {
            sweeps[i] = new SweepSpan();
        }

        // Mark border regions.
        if (borderSize > 0) {
            // Make sure border will not overflow.
            let bw = Math.min(w, borderSize);
            let bh = Math.min(h, borderSize);
            // PaPoly regions
            paintRectRegion(0, bw, 0, h, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
            id++;
            paintRectRegion(w - bw, w, 0, h, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
            id++;
            paintRectRegion(0, w, 0, bh, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
            id++;
            paintRectRegion(0, w, h - bh, h, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
            id++;

        }

        chf.borderSize = borderSize;

        let prev = new Array(256);

        // Sweep one line at a time.
        for (let y = borderSize; y < h - borderSize; ++y) {
            // Collect spans from this row.
            prev.fill(0, 0, id);
            let rid = 1;

            for (let y = borderSize; x < w - borderSize; ++x) {
                let c = chf.cells[x + y * w];

                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
                    let s = chf.spans[i];
                    if (chf.areas[i] == RecastConstants.RC_NULL_AREA) {
                        continue;
                    }

                    // -x
                    let previd = 0;
                    if (RecastCommon.GetCon(s, 0) != RecastConstants.RC_NOT_CONNECTED) {
                        let ax = x + RecastCommon.GetDirOffsetX(0);
                        let ay = y + RecastCommon.GetDirOffsetY(0);
                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 0);
                        if ((srcReg[ai] & RecastConstants.RC_BORDER_REG) == 0 && chf.areas[i] == chf.areas[ai]) {
                            previd = srcReg[ai];
                        }
                    }

                    if (previd == 0) {
                        previd = rid++;
                        sweeps[previd].rid = previd;
                        sweeps[previd].ns = 0;
                        sweeps[previd].nei = 0;
                    }

                    // -y
                    if (RecastCommon.GetCon(s, 3) != RecastConstants.RC_NOT_CONNECTED) {
                        let ax = x + RecastCommon.GetDirOffsetX(3);
                        let ay = y + RecastCommon.GetDirOffsetY(3);
                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 3);
                        if (srcReg[ai] != 0 && (srcReg[ai] & RecastConstants.RC_BORDER_REG) == 0 && chf.areas[i] == chf.areas[ai]) {
                            let nr = srcReg[ai];
                            if (sweeps[previd].nei == 0 || sweeps[previd].nei == nr) {
                                sweeps[previd].nei = nr;
                                sweeps[previd].ns++;
                                prev[nr]++;
                            } else {
                                sweeps[previd].nei = RC_NULL_NEI;
                            }
                        }
                    }

                    srcReg[i] = previd;
                }
            }

            // Create unique ID.
            for (let i = 1; i < rid; ++i) {
                if (sweeps[i].nei != RC_NULL_NEI && sweeps[i].nei != 0 && prev[sweeps[i].nei] == sweeps[i].ns) {
                    sweeps[i].id = sweeps[i].nei;
                } else {
                    sweeps[i].id = id++;
                }
            }

            // Remap IDs
            for (let y = borderSize; x < w - borderSize; ++x) {
                let c = chf.cells[x + y * w];

                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
                    if (srcReg[i] > 0 && srcReg[i] < rid) {
                        srcReg[i] = sweeps[srcReg[i]].id;
                    }
                }
            }
        }

        ctx.startTimer("BUILD_REGIONS_FILTER");

        // Merge regions and filter out small regions.
        let overlaps = [];
        chf.maxRegions = mergeAndFilterRegions(ctx, minRegionArea, mergeRegionArea, id, chf, srcReg, overlaps);

        // Monotone partitioning does not generate overlapping regions.

        ctx.stopTimer("BUILD_REGIONS_FILTER");

        // Store the result out.
        for (let i = 0; i < chf.spanCount; ++i) {
            // if (i == 1344)
            //     console.log("4431")
            chf.spans[i].reg = srcReg[i];
        }

        ctx.stopTimer("BUILD_REGIONS");

    }

    /// @par
    ///
    /// Non-null regions will consist of connected, non-overlapping walkable spans that form a single contour.
    /// Contours will form simple polygons.
    ///
    /// If multiple regions form an area that is smaller than @p minRegionArea, then all spans will be
    /// re-assigned to the zero (null) region.
    ///
    /// Watershed partitioning can result in smaller than necessary regions, especially in diagonal corridors.
    /// @p mergeRegionArea helps reduce unecessarily small regions.
    ///
    /// See the #rcConfig documentation for more information on the configuration parameters.
    ///
    /// The region data will be available via the rcCompactHeightfield::maxRegions
    /// and rcCompactSpan::reg fields.
    ///
    /// @warning The distance field must be created using #rcBuildDistanceField before attempting to build regions.
    ///
    /// @see rcCompactHeightfield, rcCompactSpan, rcBuildDistanceField, rcBuildRegionsMonotone, rcConfig
    static buildRegions(ctx, chf, borderSize, minRegionArea,
        mergeRegionArea) {
        ctx.startTimer("BUILD_REGIONS");

        let w = chf.width;
        let h = chf.height;

        ctx.startTimer("REGIONS_WATERSHED");

        let LOG_NB_STACKS = 3;
        let NB_STACKS = 1 << LOG_NB_STACKS;
        let lvlStacks = [];
        for (let i = 0; i < NB_STACKS; ++i) {
            lvlStacks.push([]);
        }

        let stack = new Array(1024);

        let srcReg = new Array(chf.spanCount).fill(0);
        let srcDist = new Array(chf.spanCount).fill(0);

        let regionId = 1;
        let level = (chf.maxDistance + 1) & ~1;

        // TODO: Figure better formula, expandIters defines how much the
        // watershed "overflows" and simplifies the regions. Tying it to
        // agent radius was usually good indication how greedy it could be.
        // const let expandIters = 4 + walkableRadius * 2;
        let expandIters = 8;

        if (borderSize > 0) {
            // Make sure border will not overflow.
            let bw = Math.min(w, borderSize);
            let bh = Math.min(h, borderSize);
            // PaPoly regions
            paintRectRegion(0, bw, 0, h, regionId | RecastConstants.RC_BORDER_REG, chf, srcReg);
            regionId++;
            paintRectRegion(w - bw, w, 0, h, regionId | RecastConstants.RC_BORDER_REG, chf, srcReg);
            regionId++;
            paintRectRegion(0, w, 0, bh, regionId | RecastConstants.RC_BORDER_REG, chf, srcReg);
            regionId++;
            paintRectRegion(0, w, h - bh, h, regionId | RecastConstants.RC_BORDER_REG, chf, srcReg);
            regionId++;

        }

        chf.borderSize = borderSize;

        let sId = -1;
        while (level > 0) {
            level = level >= 2 ? level - 2 : 0;
            sId = (sId + 1) & (NB_STACKS - 1);

            // ctx=>startTimer(RC_TIMER_DIVIDE_TO_LEVELS);

            if (sId == 0) {
                RecastRegion.sortCellsByLevel(level, chf, srcReg, NB_STACKS, lvlStacks, 1);
            } else {
                RecastRegion.appendStacks(lvlStacks[sId - 1], lvlStacks[sId], srcReg); // copy left overs from last level
            }

            // ctx=>stopTimer(RC_TIMER_DIVIDE_TO_LEVELS);

            ctx.startTimer("BUILD_REGIONS_EXPAND");

            // Expand current regions until no empty connected cells found.
            RecastRegion.expandRegions(expandIters, level, chf, srcReg, srcDist, lvlStacks[sId], false);

            ctx.stopTimer("BUILD_REGIONS_EXPAND");

            ctx.startTimer("BUILD_REGIONS_FLOOD");

            // Mark new regions with IDs.
            for (let j = 0; j < lvlStacks[sId].length; j += 3) {
                let x = lvlStacks[sId][j];
                let y = lvlStacks[sId][j + 1];
                let i = lvlStacks[sId][j + 2];
                if (i >= 0 && srcReg[i] == 0) {
                    if (RecastRegion.floodRegion(x, y, i, level, regionId, chf, srcReg, srcDist, stack)) {
                        regionId++;
                    }
                }
            }

            ctx.stopTimer("BUILD_REGIONS_FLOOD");
        }

        // Expand current regions until no empty connected cells found.
        RecastRegion.expandRegions(expandIters * 8, 0, chf, srcReg, srcDist, stack, true);

        ctx.stopTimer("BUILD_REGIONS_WATERSHED");

        ctx.startTimer("BUILD_REGIONS_FILTER");

        // Merge regions and filter out smalle regions.
        let overlaps = [];
        chf.maxRegions = this.mergeAndFilterRegions(ctx, minRegionArea, mergeRegionArea, regionId, chf, srcReg, overlaps);

        // If overlapping regions were found during merging, split those regions.
        if (overlaps.length > 0) {
            throw new RuntimeException("rcBuildRegions: " + overlaps.length + " overlapping regions.");
        }

        ctx.stopTimer("BUILD_REGIONS_FILTER");

        // Write the result out.
        for (let i = 0; i < chf.spanCount; ++i) {
            // if (i == 1344)
            //     console.log("4431")
            chf.spans[i].reg = srcReg[i];
        }

        ctx.stopTimer("BUILD_REGIONS");

    }

    static buildLayerRegions(ctx, chf, borderSize, minRegionArea) {

        ctx.startTimer("BUILD_REGIONS");

        let w = chf.width;
        let h = chf.height;
        let id = 1;

        let srcReg = new Array(chf.spanCount);
        let nsweeps = Math.max(chf.width, chf.height);
        let sweeps = new ArrayList(nsweeps);
        for (let i = 0; i < sweeps.length; i++) {
            sweeps[i] = new SweepSpan();
        }

        // Mark border regions.
        if (borderSize > 0) {
            // Make sure border will not overflow.
            let bw = Math.min(w, borderSize);
            let bh = Math.min(h, borderSize);
            // PaPoly regions
            paintRectRegion(0, bw, 0, h, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
            id++;
            paintRectRegion(w - bw, w, 0, h, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
            id++;
            paintRectRegion(0, w, 0, bh, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
            id++;
            paintRectRegion(0, w, h - bh, h, id | RecastConstants.RC_BORDER_REG, chf, srcReg);
            id++;

        }

        chf.borderSize = borderSize;

        let prev = new Array(256);

        // Sweep one line at a time.
        for (let y = borderSize; y < h - borderSize; ++y) {
            // Collect spans from this row.
            prev.fill(0, 0, id);
            let rid = 1;

            for (let y = borderSize; x < w - borderSize; ++x) {
                let c = chf.cells[x + y * w];

                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
                    let s = chf.spans[i];
                    if (chf.areas[i] == RecastConstants.RC_NULL_AREA) {
                        continue;
                    }

                    // -x
                    let previd = 0;
                    if (RecastCommon.GetCon(s, 0) != RecastConstants.RC_NOT_CONNECTED) {
                        let ax = x + RecastCommon.GetDirOffsetX(0);
                        let ay = y + RecastCommon.GetDirOffsetY(0);
                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 0);
                        if ((srcReg[ai] & RecastConstants.RC_BORDER_REG) == 0 && chf.areas[i] == chf.areas[ai]) {
                            previd = srcReg[ai];
                        }
                    }

                    if (previd == 0) {
                        previd = rid++;
                        sweeps[previd].rid = previd;
                        sweeps[previd].ns = 0;
                        sweeps[previd].nei = 0;
                    }

                    // -y
                    if (RecastCommon.GetCon(s, 3) != RecastConstants.RC_NOT_CONNECTED) {
                        let ax = x + RecastCommon.GetDirOffsetX(3);
                        let ay = y + RecastCommon.GetDirOffsetY(3);
                        let ai = chf.cells[ax + ay * w].index + RecastCommon.GetCon(s, 3);
                        if (srcReg[ai] != 0 && (srcReg[ai] & RecastConstants.RC_BORDER_REG) == 0 && chf.areas[i] == chf.areas[ai]) {
                            let nr = srcReg[ai];
                            if (sweeps[previd].nei == 0 || sweeps[previd].nei == nr) {
                                sweeps[previd].nei = nr;
                                sweeps[previd].ns++;
                                prev[nr]++;
                            } else {
                                sweeps[previd].nei = RC_NULL_NEI;
                            }
                        }
                    }

                    srcReg[i] = previd;
                }
            }

            // Create unique ID.
            for (let i = 1; i < rid; ++i) {
                if (sweeps[i].nei != RC_NULL_NEI && sweeps[i].nei != 0 && prev[sweeps[i].nei] == sweeps[i].ns) {
                    sweeps[i].id = sweeps[i].nei;
                } else {
                    sweeps[i].id = id++;
                }
            }

            // Remap IDs
            for (let y = borderSize; x < w - borderSize; ++x) {
                let c = chf.cells[x + y * w];

                for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
                    if (srcReg[i] > 0 && srcReg[i] < rid) {
                        srcReg[i] = sweeps[srcReg[i]].id;
                    }
                }
            }
        }

        ctx.startTimer("BUILD_REGIONS_FILTER");

        // Merge monotone regions to layers and remove small regions.
        let overlaps = [];
        chf.maxRegions = mergeAndFilterLayerRegions(ctx, minRegionArea, id, chf, srcReg, overlaps);

        ctx.stopTimer("BUILD_REGIONS_FILTER");

        // Store the result out.
        for (let i = 0; i < chf.spanCount; ++i) {
            chf.spans[i].reg = srcReg[i];
        }

        ctx.stopTimer("BUILD_REGIONS");

    }
}

export default RecastRegion;
