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

/** A compact, static heightfield representing unobstructed space. */
class CompactHeightfield {

	/** The width of the heightfield. (APoly the x-axis in cell units.) */
 width = 0;
	/** The height of the heightfield. (APoly the z-axis in cell units.) */
 height = 0;
	/** The number of spans in the heightfield. */
 spanCount = 0;
	/** The walkable height used during the build of the field.  (See: RecastConfig::walkableHeight) */
 walkableHeight = 0;
	/** The walkable climb used during the build of the field. (See: RecastConfig::walkableClimb) */
 walkableClimb = 0;
	/** The AABB border size used during the build of the field. (See: RecastConfig::borderSize) */
 borderSize = 0 ;
	/** The maximum distance value of any span within the field. */
 maxDistance = 0;
	/** The maximum region id of any span within the field. */
 maxRegions;
	/** The minimum bounds in world space. [(x, y, z)] */
 bmin = new Array(3);
	/** The maximum bounds in world space. [(x, y, z)] */
 bmax = new Array(3);
	/** The size of each cell. (On the xz-plane.) */
 cs = 0;
	/** The height of each cell. (The minimum increment aPoly the y-axis.) */
 ch = 0;
	/** Array of cells. [Size: #width*#height] */
cells = [];
	/** Array of spans. [Size: #spanCount] */
spans = [];
	/** Array containing border distance data. [Size: #spanCount] */
 dist = [];
	/** Array containing area id data. [Size: #spanCount] */
 areas = [];

}

export default CompactHeightfield;
