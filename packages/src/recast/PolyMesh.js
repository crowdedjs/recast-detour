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

/** Represents a polygon mesh suitable for use in building a navigation mesh. */
class PolyMesh {

	/** The mesh vertices. [Form: (x, y, z) coordinates * #nverts] */
	verts = [];
	/** Polygon and neighbor data. [Length: #maxpolys * 2 * #nvp] */
	polys = [];
	/** The region id assigned to each polygon. [Length: #maxpolys] */
	regs = [];
	/** The area id assigned to each polygon. [Length: #maxpolys] */
	areas = [];
	/** The number of vertices. */
	nverts = 0;
	/** The number of polygons. */
	npolys = 0;
	/** The maximum number of vertices per polygon. */
	nvp = 0;
	/** The number of allocated polygons. */
	maxpolys = 0;
	/** The user defined flags for each polygon. [Length: #maxpolys] */
	flags = [];
	/** The minimum bounds in world space. [(x, y, z)] */
	bmin = new Array(3);
	/** The maximum bounds in world space. [(x, y, z)] */
	bmax = new Array(3);
	/** The size of each cell. (On the xz-plane.) */
	cs = 0;
	/** The height of each cell. (The minimum increment aPoly the y-axis.) */
	ch = 0;
	/** The AABB border size used to generate the source data from which the mesh was derived. */
	borderSize = 0;
	/** The max error of the polygon edges in the mesh. */
	maxEdgeError = 0;
}

export default PolyMesh;
