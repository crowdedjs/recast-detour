
/// Represents a set of heightfield layers.
/// @ingroup recast
/// @see rcAllocHeightfieldLayerSet, rcFreeHeightfieldLayerSet 
class HeightfieldLayerSet {

	/// Represents a heightfield layer within a layer set.
	/// @see rcHeightfieldLayerSet
static class HeightfieldLayer {
let bmin = new Array(3);				///< The minimum bounds in world space. [(x, y, z)]
let bmax = new Array(3);				///< The maximum bounds in world space. [(x, y, z)]
 cs;					///< The size of each cell. (On the xz-plane.)
 ch;					///< The height of each cell. (The minimum increment aPoly the y-axis.)
let width;					///< The width of the heightfield. (APoly the x-axis in cell units.)
let height;					///< The height of the heightfield. (APoly the z-axis in cell units.)
let minx;					///< The minimum x-bounds of usable data.
let maxx;					///< The maximum x-bounds of usable data.
let miny;					///< The minimum y-bounds of usable data. (APoly the z-axis.)
let maxy;					///< The maximum y-bounds of usable data. (APoly the z-axis.)
let hmin;					///< The minimum height bounds of usable data. (APoly the y-axis.)
let hmax;					///< The maximum height bounds of usable data. (APoly the y-axis.)
let heights;		///< The heightfield. [Size: width * height]
let areas;		///< Area ids. [Size: Same as #heights]
let cons;		///< Packed neighbor connection information. [Size: Same as #heights]
	}

HeightfieldLayer[] layers;			///< The layers in the set. [Size: #nlayers]
}
