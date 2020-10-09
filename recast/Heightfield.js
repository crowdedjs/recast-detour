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

/** Represents a heightfield layer within a layer set. */
class Heightfield {

	/** The width of the heightfield. (APoly the x-axis in cell units.) */
  width = 0;
	/** The height of the heightfield. (APoly the z-axis in cell units.) */
  height = 0;
	/** The minimum bounds in world space. [(x, y, z)] */
bmin;
	/** The maximum bounds in world space. [(x, y, z)] */
 bmax;
	/** The size of each cell. (On the xz-plane.) */
  cs;
	/** The height of each cell. (The minimum increment aPoly the y-axis.) */
  ch;
	/** Heightfield of spans (width*height). */
 spans = [];

constructor(width,  height, bmin, bmax,  cs,  ch) {
		this.width = width;
		this.height = height;
		this.bmin = bmin;
		this.bmax = bmax;
		this.cs = cs;
		this.ch = ch;
		this.spans = new Array(width * height);

	}
}

export default Heightfield;
