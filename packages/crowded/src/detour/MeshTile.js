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

import NavMesh from "./NavMesh.js"

/**
 * Defines a navigation mesh tile.
 */
class MeshTile {
	index = 0;
	/** Counter describing modifications to the tile. */
	salt = 0;
	/** The tile data. */
	data = null;
	/** The tile links. */
	links = [];
	/** Index to the next free link. */
	linksFreeList = NavMesh.DT_NULL_LINK;
	/** Tile flags. (See: #dtTileFlags) */
	flags = 0;
	/** The next free tile, or the next tile in the spatial grid. */
	next = null;
	constructor(index) {
		this.index = index;
	}

}
export default MeshTile;
