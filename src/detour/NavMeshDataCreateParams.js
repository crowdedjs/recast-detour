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

/// Represents the source data used to build an navigation mesh tile.
class NavMeshDataCreateParams {

	/// @name Polygon Mesh Attributes
	/// Used to create the base navigation graph.
	/// See #rcPolyMesh for details related to these attributes.
	/// @{

 verts = [];			///< The polygon mesh vertices. [(x, y, z) * #vertCount] [Unit: vx]
 vertCount = 0;		///< The number vertices in the polygon mesh. [Limit: >= 3]
 polys = [];			///< The polygon data. [Size: #polyCount * 2 * #nvp]
 polyFlags = [];		///< The user defined flags assigned to each polygon. [Size: #polyCount]
 polyAreas = [];		///< The user defined area ids assigned to each polygon. [Size: #polyCount]
 polyCount = 0;		///< Number of polygons in the mesh. [Limit: >= 1]
 nvp = 0;				///< Number maximum number of vertices per polygon. [Limit: >= 3]

	/// @}
	/// @name Height Detail Attributes (Optional)
	/// See #rcPolyMeshDetail for details related to these attributes.
	/// @{

 detailMeshes = [];			///< The height detail sub-mesh data. [Size: 4 * #polyCount]
 detailVerts = [];		///< The detail mesh vertices. [Size: 3 * #detailVertsCount] [Unit: wu]
 detailVertsCount = 0;		///< The number of vertices in the detail mesh.
 detailTris = [];			///< The detail mesh triangles. [Size: 4 * #detailTriCount]
 detailTriCount = 0;			///< The number of triangles in the detail mesh.

	/// @}
	/// @name Off-Mesh Connections Attributes (Optional)
	/// Used to define a custom point-to-po edge within the navigation graph, an 
	/// off-mesh connection is a user defined traversable connection made up to two vertices, 
	/// at least one of which resides within a navigation mesh polygon.
	/// @{

	/// Off-mesh connection vertices. [(ax, ay, az, bx, by, bz) * #offMeshConCount] [Unit: wu]
 offMeshConVerts = [];
	/// Off-mesh connection radii. [Size: #offMeshConCount] [Unit: wu]
 offMeshConRad = [];
	/// User defined flags assigned to the off-mesh connections. [Size: #offMeshConCount]
 offMeshConFlags = [];
	/// User defined area ids assigned to the off-mesh connections. [Size: #offMeshConCount]
 offMeshConAreas = [];
	/// The permitted travel direction of the off-mesh connections. [Size: #offMeshConCount]
	///
	/// 0 = Travel only from endpo A to endpo B.<br/>
	/// #DT_OFFMESH_CON_BIDIR = Bidirectional travel.
 offMeshConDir = [];	
	/// The user defined ids of the off-mesh connection. [Size: #offMeshConCount]
 offMeshConUserID = [];
	/// The number of off-mesh connections. [Limit: >= 0]
 offMeshConCount = 0;

	/// @}
	/// @name Tile Attributes
	/// @note The tile grid/layer data can be left at zero if the destination is a single tile mesh.
	/// @{

 userId = 0;	///< The user defined id of the tile.
 tileX = 0;				///< The tile's x-grid location within the multi-tile destination mesh. (A the x-axis.)
 tileY = 0;				///< The tile's y-grid location within the multi-tile desitation mesh. (A the z-axis.)
 tileLayer = 0;			///< The tile's layer within the layered destination mesh. [Limit: >= 0] (A the y-axis.)
 bmin = [];			///< The minimum bounds of the tile. [(x, y, z)] [Unit: wu]
 bmax = [];			///< The maximum bounds of the tile. [(x, y, z)] [Unit: wu]

	/// @}
	/// @name General Configuration Attributes
	/// @{

 walkableHeight = 0;	///< The agent height. [Unit: wu]
 walkableRadius= 0;	///< The agent radius. [Unit: wu]
 walkableClimb = 0;	///< The agent maximum traversable ledge. (Up/Down) [Unit: wu]
 cs= 0 ;				///< The xz-plane cell size of the polygon mesh. [Limit: > 0] [Unit: wu]
 ch = 0;				///< The y-axis cell height of the polygon mesh. [Limit: > 0] [Unit: wu]

	/// True if a bounding volume tree should be built for the tile.
	/// @note The BVTree is not normally needed for layered navigation meshes.
 buildBvTree = false;

	/// @}

}

export default NavMeshDataCreateParams;