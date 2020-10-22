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


/** Represents a group of related contours. */
class ContourSet {

	/** A list of the contours in the set. */
  conts = [];
	/** The minimum bounds in world space. [(x, y, z)] */
 bmin = new Array(3);
	/** The maximum bounds in world space. [(x, y, z)] */
 bmax = new Array(3);
	/** The size of each cell. (On the xz-plane.) */
 cs = 0;
	/** The height of each cell. (The minimum increment aPoly the y-axis.) */
 ch = 0;
	/** The width of the set. (APoly the x-axis in cell units.) */
 width = 0;
	/** The height of the set. (APoly the z-axis in cell units.) */
 height = 0;
	/** The AABB border size used to generate the source data from which the contours were derived. */
 borderSize = 0;
	/** The max edge error that this contour set was simplified with. */
 maxError = 0;
}

export default ContourSet;
