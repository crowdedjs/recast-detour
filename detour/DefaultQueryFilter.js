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


/**
 * <b>The Default Implementation</b>
 * 
 * At construction: All area costs default to 1.0. All flags are included and none are excluded.
 * 
 * If a polygon has both an include and an exclude flag, it will be excluded.
 * 
 * The way filtering works, a navigation mesh polygon must have at least one flag set to ever be considered by a query.
 * So a polygon with no flags will never be considered.
 * 
 * Setting the include flags to 0 will result in all polygons being excluded.
 * 
 * <b>Custom Implementations</b>
 * 
 * Implement a custom query filter by overriding the virtual passFilter() and getCost() functions. If this is done, both
 * functions should be as fast as possible. Use cached local copies of data rather than accessing your own objects where
 * possible.
 * 
 * Custom implementations do not need to adhere to the flags or cost logic used by the default implementation.
 * 
 * In order for A* searches to work properly, the cost should be proportional to the travel distance. Implementing a
 * cost modifier less than 1.0 is likely to lead to problems during pathfinding.
 * 
 * @see NavMeshQuery
 */

import NavMesh from "./NavMesh.js"
import DetourCommon from "./DetourCommon.js"
// import QueryFilter from "./QueryFilter.js"

 class DefaultQueryFilter /*extends QueryFilter*/ {

m_excludeFlags = 0;
m_includeFlags = 0;
m_areaCost = new Array(NavMesh.DT_MAX_AREAS);

// DefaultQueryFilter() {
// 		this.m_includeFlags = 0xfff;
// 		this.m_excludeFlags = 0;
// 		for (let i = 0; i < NavMesh.DT_MAX_AREAS; ++i) {
// 			m_areaCost[i] = 1.0;
// 		}
// 	}

constructor(includeFlags = 0xffff,  excludeFlags = 0, areaCost = new Array(NavMesh.DT_MAX_AREAS).fill(1,0,NavMesh.DT_MAX_AREAS)) {
		this.m_includeFlags = includeFlags;
		this.m_excludeFlags = excludeFlags;
		for (let i = 0; i < Math.min(NavMesh.DT_MAX_AREAS, areaCost.length); ++i) {
			this.m_areaCost[i] = areaCost[i];
		}
		for (let i = areaCost.length; i < NavMesh.DT_MAX_AREAS; ++i) {
			this.m_areaCost[i] = 1.0;
		}
	}

	and(v1, v2) {
    var hi = 0x80000000;
    var low = 0x7fffffff;
    var hi1 = ~~(v1 / hi);
    var hi2 = ~~(v2 / hi);
    var low1 = v1 & low;
    var low2 = v2 & low;
    var h = hi1 & hi2;
    var l = low1 & low2;
    return h*hi + l;
}
	or(v1, v2) {
    var hi = 0x80000000;
    var low = 0x7fffffff;
    var hi1 = ~~(v1 / hi);
    var hi2 = ~~(v2 / hi);
    var low1 = v1 & low;
    var low2 = v2 & low;
    var h = hi1 | hi2;
    var l = low1 | low2;
    return h*hi + l;
}
	

passFilter(ref, tile, poly) {
		return this.and(poly.flags , this.m_includeFlags) != 0 && this.and(poly.flags , this.m_excludeFlags) == 0;
	}

 getCost( pa, pb,  prevRef, prevTile, prevPoly,  curRef,
			curTile, curPoly,  nextRef, nextTile, nextPoly) {
		return DetourCommon.vDist2(pa, pb) * this.m_areaCost[curPoly.getArea()];
	}

}

export default DefaultQueryFilter;
