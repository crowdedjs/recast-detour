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

/** Provides high level information related to a dt object.*/
class MeshHeader {
	/** A magic number used to detect compatibility of navigation tile data. */
static  DT_NAVMESH_MAGIC = 'D' << 24 | 'N' << 16 | 'A' << 8 | 'V';
	/** A version number used to detect compatibility of navigation tile data.*/
static  DT_NAVMESH_VERSION = 7;
static  DT_NAVMESH_VERSION_RECAST4J = 0x8807;
	/** A magic number used to detect the compatibility of navigation tile states.*/
static  DT_NAVMESH_STATE_MAGIC = 'D' << 24 | 'N' << 16 | 'M' << 8 | 'S';
	/** A version number used to detect compatibility of navigation tile states.*/
static  DT_NAVMESH_STATE_VERSION = 1;

	/** < Tile magic number. (Used to identify the data format.)*/
magic;
	/** < Tile data format version number.*/
version;
	/** < The x-position of the tile within the dtNavMesh tile grid. (x, y, layer)*/
x;
	/** < The y-position of the tile within the dtNavMesh tile grid. (x, y, layer)*/
y;
	/** < The layer of the tile within the dtNavMesh tile grid. (x, y, layer)*/
layer;
	/** < The user defined id of the tile.*/
userId;
	/** < The number of polygons in the tile.*/
 polyCount;
	/** < The number of vertices in the tile.*/
 vertCount;
	/** < The number of allocated links.*/
 maxLinkCount;
	/** < The number of sub-meshes in the detail mesh.*/
 detailMeshCount;
	/** The number of unique vertices in the detail mesh. (In addition to the polygon vertices.)*/
 detailVertCount;
	/** < The number of triangles in the detail mesh.*/
 detailTriCount;
	/** < The number of bounding volume nodes. (Zero if bounding volumes are disabled.)*/
 bvNodeCount;
	/** < The number of off-mesh connections.*/
 offMeshConCount;
	/** < The index of the first polygon which is an off-mesh connection.*/
 offMeshBase;
	/** < The height of the agents using the tile.*/
 walkableHeight;
	/** < The radius of the agents using the tile.*/
 walkableRadius;
	/** < The maximum climb height of the agents using the tile.*/
 walkableClimb;
	/** < The minimum bounds of the tile's AABB. [(x, y, z)]*/
 bmin = new Array(3);
	/** < The maximum bounds of the tile's AABB. [(x, y, z)]*/
 bmax = new Array(3);
	/** The bounding volume quantization factor.*/
 bvQuantFactor;
}

export default MeshHeader;
