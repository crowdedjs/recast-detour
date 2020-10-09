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

/**
 * Defines a link between polygons.
 * 
 * @note This structure is rarely if ever used by the end user.
 * @see MeshTile
 */
class Link {
	/** Neighbour reference. (The neighbor that is linked to.) */
	 ref = 0;
	/** Index of the next link. */
	next = 0;
	/** Index of the polygon edge that owns this link. */
	 edge = 0;
	/** If a boundary link, defines on which side the link is. */
	 side = 0;
	/** If a boundary link, defines the minimum sub-edge area. */
	bmin = 0;
	/** If a boundary link, defines the maximum sub-edge area. */
	bmax = 0;

}

export default Link;