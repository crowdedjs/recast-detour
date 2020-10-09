
/// Represents a set of heightfield layers.
/// @ingroup recast
/// @see rcAllocHeightfieldLayerSet, rcFreeHeightfieldLayerSet 
class HeightfieldLayerSet {

	/// Represents a heightfield layer within a layer set.
	/// @see rcHeightfieldLayerSet
	static HeightfieldLayer = class HeightfieldLayer {
		bmin = new Array(3);				///< The minimum bounds in world space. [(x, y, z)]
		bmax = new Array(3);				///< The maximum bounds in world space. [(x, y, z)]
		cs;					///< The size of each cell. (On the xz-plane.)
		ch;					///< The height of each cell. (The minimum increment aPoly the y-axis.)
		width;					///< The width of the heightfield. (APoly the x-axis in cell units.)
		height;					///< The height of the heightfield. (APoly the z-axis in cell units.)
		minx;					///< The minimum x-bounds of usable data.
		maxx;					///< The maximum x-bounds of usable data.
		miny;					///< The minimum y-bounds of usable data. (APoly the z-axis.)
		maxy;					///< The maximum y-bounds of usable data. (APoly the z-axis.)
		hmin;					///< The minimum height bounds of usable data. (APoly the y-axis.)
		hmax;					///< The maximum height bounds of usable data. (APoly the y-axis.)
		heights;		///< The heightfield. [Size: width * height]
		areas;		///< Area ids. [Size: Same as #heights]
		cons;		///< Packed neighbor connection information. [Size: Same as #heights]
	}

	layers;			///< The layers in the set. [Size: #nlayers]
}

export default HeightfieldLayerSet;