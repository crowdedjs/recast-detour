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

import NavMeshParams from "./NavMeshParams.js"
import DetourCommon from "./DetourCommon.js"
import MeshTile from "./MeshTile.js"
import Poly from "./Poly.js"
import Link from "./Link.js"
import ClosestPointOnPolyResult from "./ClosestPointOnPolyResult.js"
import FindNearestPolyResult from "./FindNearestPolyResult.js"

function arraycopy(one, oneStart, two, twoStart, len) {
	for (let i = 0; i < len; i++) {
		two[twoStart + i] = one[oneStart + i];
	}
}

class NavMesh {

	static DT_SALT_BITS = 16;
	static DT_TILE_BITS = 28;
	static DT_POLY_BITS = 20;

	/// A flag that indicates that an entity links to an external entity.
	/// (E.g. A polygon edge is a portal that links to another polygon.)
	static DT_EXT_LINK = 0x8000;

	/// A value that indicates the entity does not link to anything.
	// static DT_NULL_LINK = 0xffffffff;
	static DT_NULL_LINK = -1;

	/// A flag that indicates that an off-mesh connection can be traversed in
	/// both directions. (Is bidirectional.)
	static DT_OFFMESH_CON_BIDIR = 1;

	/// The maximum number of user defined area ids.
	static DT_MAX_AREAS = 64;

	/// Limit raycasting during any angle pahfinding
	/// The limit is given as a multiple of the character radius
	static DT_RAY_CAST_LIMIT_PROPORTIONS = 50.0;

	m_params = null/// < Current initialization params. TODO: do not store this info twice.
	m_orig = []; /// < Origin of the tile (0,0)
	//  m_orig[3]; ///< Origin of the tile (0,0)
	m_tileWidth = 0;
	m_tileHeight = 0; /// < Dimensions of each tile.
	m_maxTiles = 0; /// < Max number of tiles.
	m_tileLutSize = 0; /// < Tile hash lookup size (must be pot).
	m_tileLutMask = 0; /// < Tile hash lookup mask.
	m_posLookup = []; /// < Tile hash lookup.
	m_nextFree = null; /// < Freelist of tiles.
	m_tiles = []; /// < List of tiles.
	/** The maximum number of vertices per navigation polygon. */
	m_maxVertPerPoly = 0;
	m_tileCount = 0;

	/**
	 * The maximum number of tiles supported by the navigation mesh.
	 * 
	 * @return The maximum number of tiles supported by the navigation mesh.
	 */
	getMaxTiles() {
		return this.m_maxTiles;
	}

	/**
	 * Returns tile in the tile array.
	 */
	getTile(i) {
		return this.m_tiles[i];
	}

	/**
	 * Gets the polygon reference for the tile's base polygon.
	 * 
	 * @param tile
	 *            The tile.
	 * @return The polygon reference for the base polygon in the specified tile.
	 */
	getPolyRefBase(tile) {
		if (tile == null)
			return 0;
		let it = tile.index;
		return NavMesh.encodePolyId(tile.salt, it, 0);
	}

	/**
	 * Derives a standard polygon reference.
	 * 
	 * @note This function is generally meant for internal use only.
	 * @param salt
	 *            The tile's salt value.
	 * @param it
	 *            The index of the tile.
	 * @param ip
	 *            The index of the polygon within the tile.
	 * @return encoded polygon reference
	 */
	//https://stackoverflow.com/a/337572/10047920
	static lshift(num, bits) {
		return num * Math.pow(2, bits);
	}
	static rshift(num, bits) {
		return Math.floor(num / Math.pow(2, bits));
	}
	//https://stackoverflow.com/a/43666199/10047920
	static and(v1, v2) {
		var hi = 0x80000000;
		var low = 0x7fffffff;
		var hi1 = ~~(v1 / hi);
		var hi2 = ~~(v2 / hi);
		var low1 = v1 & low;
		var low2 = v2 & low;
		var h = hi1 & hi2;
		var l = low1 & low2;
		return h * hi + l;
	}
	static or(v1, v2) {
		var hi = 0x80000000;
		var low = 0x7fffffff;
		var hi1 = ~~(v1 / hi);
		var hi2 = ~~(v2 / hi);
		var low1 = v1 & low;
		var low2 = v2 & low;
		var h = hi1 | hi2;
		var l = low1 | low2;
		return h * hi + l;
	}
	static encodePolyId(salt, it, ip) {
		let a = (NavMesh.lshift(salt, NavMesh.DT_POLY_BITS + NavMesh.DT_TILE_BITS));
		let b = (NavMesh.lshift(it, NavMesh.DT_POLY_BITS));
		let c = ip;
		return NavMesh.or(NavMesh.or(a, b), ip);
	}

	/// Decodes a standard polygon reference.
	/// @note This function is generally meant for internal use only.
	/// @param[in] ref The polygon reference to decode.
	/// @param[out] salt The tile's salt value.
	/// @param[out] it The index of the tile.
	/// @param[out] ip The index of the polygon within the tile.
	/// @see #encodePolyId
	static decodePolyId(ref) {
		let salt;
		let it;
		let ip;
		let saltMask = NavMesh.lshift(1, NavMesh.DT_SALT_BITS) - 1;
		let tileMask = NavMesh.lshift(1, NavMesh.DT_TILE_BITS) - 1;
		let polyMask = NavMesh.lshift(1, NavMesh.DT_POLY_BITS) - 1;
		salt = Math.floor(NavMesh.and(NavMesh.rshift(ref, (NavMesh.DT_POLY_BITS + NavMesh.DT_TILE_BITS)), saltMask));
		it = Math.floor(NavMesh.and(NavMesh.rshift(ref, NavMesh.DT_POLY_BITS), tileMask));
		ip = Math.floor(NavMesh.and(ref, polyMask));
		return [salt, it, ip];
	}

	/// Extracts a tile's salt value from the specified polygon reference.
	/// @note This function is generally meant for internal use only.
	/// @param[in] ref The polygon reference.
	/// @see #encodePolyId
	static decodePolyIdSalt(ref) {
		let saltMask = (1 << NavMesh.DT_SALT_BITS) - 1;
		return Math.floor((ref >> (NavMesh.DT_POLY_BITS + NavMesh.DT_TILE_BITS)) & saltMask);
	}

	/// Extracts the tile's index from the specified polygon reference.
	/// @note This function is generally meant for internal use only.
	/// @param[in] ref The polygon reference.
	/// @see #encodePolyId
	static decodePolyIdTile(ref) {
		let tileMask = (1 << NavMesh.DT_TILE_BITS) - 1;
		return Math.floor((ref >> NavMesh.DT_POLY_BITS) & tileMask);
	}

	/// Extracts the polygon's index (within its tile) from the specified
	/// polygon reference.
	/// @note This function is generally meant for internal use only.
	/// @param[in] ref The polygon reference.
	/// @see #encodePolyId
	static decodePolyIdPoly(ref) {
		let polyMask = (1 << NavMesh.DT_POLY_BITS) - 1;
		return Math.floor(ref & polyMask);
	}

	allocLink(tile) {
		if (tile.linksFreeList == NavMesh.DT_NULL_LINK) {
			let link = new Link();
			link.next = NavMesh.DT_NULL_LINK;
			tile.links.push(link);
			return tile.links.length - 1;
		}
		let link = tile.linksFreeList;
		tile.linksFreeList = tile.links[link].next;
		return link;
	}

	freeLink(tile, link) {
		tile.links[link].next = tile.linksFreeList;
		tile.linksFreeList = link;
	}

	/**
	 * Calculates the tile grid location for the specified world position.
	 * 
	 * @param pos
	 *            The world position for the query. [(x, y, z)]
	 * @return 2-element let array with (tx,ty) tile location
	 */
	calcTileLoc(pos) {
		let tx = Math.floor((pos[0] - this.m_orig[0]) / this.m_tileWidth);
		let ty = Math.floor((pos[2] - this.m_orig[2]) / this.m_tileHeight);
		return [tx, ty];
	}

	getTileAndPolyByRef(ref) {
		if (ref == 0) {
			throw new IllegalArgumentException("ref = 0");
		}
		let saltitip = NavMesh.decodePolyId(ref);
		let salt = saltitip[0];
		let it = saltitip[1];
		let ip = saltitip[2];
		if (it >= this.m_maxTiles)
			throw new IllegalArgumentException("tile > m_maxTiles");
		if (this.m_tiles[it].salt != salt || this.m_tiles[it].data.header == null)
			throw new IllegalArgumentException("Invalid salt or header");
		if (ip >= this.m_tiles[it].data.header.polyCount)
			throw new IllegalArgumentException("poly > polyCount");
		return [this.m_tiles[it], this.m_tiles[it].data.polys[ip]];
	}

	/// @par
	///
	/// @warning Only use this function if it is known that the provided polygon
	/// reference is valid. This function is faster than #getTileAndPolyByRef,
	/// but
	/// it does not validate the reference.
	getTileAndPolyByRefUnsafe(ref) {
		let saltitip = NavMesh.decodePolyId(ref);
		let it = saltitip[1];
		let ip = saltitip[2];
		return [this.m_tiles[it], this.m_tiles[it].data.polys[ip]];
	}

	isValidPolyRef(ref) {
		if (ref == 0)
			return false;
		let saltitip = NavMesh.decodePolyId(ref);
		let salt = saltitip[0];
		let it = saltitip[1];
		let ip = saltitip[2];
		if (it >= this.m_maxTiles)
			return false;
		if (this.m_tiles[it].salt != salt || this.m_tiles[it].data == null)
			return false;
		if (ip >= this.m_tiles[it].data.header.polyCount)
			return false;
		return true;
	}

	getParams() {
		return this.m_params;
	}

	constructor(one, maxVertsPerPoly, flags) {

		if (flags || flags == 0) {
			this._constructor(NavMesh.getNavMeshParams(one), maxVertsPerPoly);
			this.addTile(one, flags, 0);
		}
		else {
			this._constructor(one, maxVertsPerPoly);
		}
	}

	_constructor(params, maxVertsPerPoly) {
		this.m_params = params;
		this.m_orig = params.orig;
		this.m_tileWidth = params.tileWidth;
		this.m_tileHeight = params.tileHeight;
		// Init tiles
		this.m_maxTiles = params.maxTiles;
		this.m_maxVertPerPoly = maxVertsPerPoly;
		let lutsize = DetourCommon.nextPow2(params.maxTiles / 4);
		if (lutsize == 0)
			lutsize = 1;
		this.m_tileLutSize = lutsize;
		this.m_tileLutMask = this.m_tileLutSize - 1;
		this.m_tiles = new Array(this.m_maxTiles);
		this.m_posLookup = new Array(this.m_tileLutSize);
		this.m_nextFree = null;
		for (let i = this.m_maxTiles - 1; i >= 0; --i) {
			this.m_tiles[i] = new MeshTile(i);
			this.m_tiles[i].salt = 1;
			this.m_tiles[i].next = this.m_nextFree;
			this.m_nextFree = this.m_tiles[i];
		}

	}

	static getNavMeshParams(data) {
		let params = new NavMeshParams();
		DetourCommon.vCopy(params.orig, data.header.bmin);
		params.tileWidth = data.header.bmax[0] - data.header.bmin[0];
		params.tileHeight = data.header.bmax[2] - data.header.bmin[2];
		params.maxTiles = 1;
		params.maxPolys = data.header.polyCount;
		return params;
	}

	// TODO: These methods are duplicates from dtNavMeshQuery, but are needed
	// for off-mesh connection finding.

	queryPolygonsInTile(tile, qmin, qmax) {
		let polys = [];
		if (tile.data.bvTree != null) {
			let nodeIndex = 0;
			let tbmin = tile.data.header.bmin;
			let tbmax = tile.data.header.bmax;
			let qfac = tile.data.header.bvQuantFactor;
			// Calculate quantized box
			let bmin = new Array(3);
			let bmax = new Array(3);
			// dtClamp query box to world box.
			let minx = DetourCommon.clamp(qmin[0], tbmin[0], tbmax[0]) - tbmin[0];
			let miny = DetourCommon.clamp(qmin[1], tbmin[1], tbmax[1]) - tbmin[1];
			let minz = DetourCommon.clamp(qmin[2], tbmin[2], tbmax[2]) - tbmin[2];
			let maxx = DetourCommon.clamp(qmax[0], tbmin[0], tbmax[0]) - tbmin[0];
			let maxy = DetourCommon.clamp(qmax[1], tbmin[1], tbmax[1]) - tbmin[1];
			let maxz = DetourCommon.clamp(qmax[2], tbmin[2], tbmax[2]) - tbmin[2];
			// Quantize
			bmin[0] = NavMesh.and(Math.floor((qfac * minx)), 0xfffe);
			bmin[1] = NavMesh.and(Math.floor((qfac * miny)), 0xfffe);
			bmin[2] = NavMesh.and(Math.floor((qfac * minz)), 0xfffe);
			bmax[0] = NavMesh.or(Math.floor((qfac * maxx + 1)), 1);
			bmax[1] = NavMesh.or(Math.floor((qfac * maxy + 1)), 1);
			bmax[2] = NavMesh.or(Math.floor((qfac * maxz + 1)), 1);

			// Traverse tree
			let base = this.getPolyRefBase(tile);
			let end = tile.data.header.bvNodeCount;
			while (nodeIndex < end) {
				let node = tile.data.bvTree[nodeIndex];
				let overlap = DetourCommon.overlapQuantBounds(bmin, bmax, node.bmin, node.bmax);
				let isLeafNode = node.i >= 0;

				if (isLeafNode && overlap) {
					polys.push(NavMesh.or(base, node.i));
				}

				if (overlap || isLeafNode)
					nodeIndex++;
				else {
					let escapeIndex = -node.i;
					nodeIndex += escapeIndex;
				}
			}

			return polys;
		} else {
			bmin = [null, null, null];
			bmax = [null, null, null];
			let base = this.getPolyRefBase(tile);
			for (let i = 0; i < tile.data.header.polyCount; ++i) {
				let p = tile.data.polys[i];
				// Do not return off-mesh connection polygons.
				if (p.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION)
					continue;
				// Calc polygon bounds.
				let v = p.verts[0] * 3;
				DetourCommon.vCopy(bmin, tile.data.verts, v);
				DetourCommon.vCopy(bmax, tile.data.verts, v);
				for (let j = 1; j < p.vertCount; ++j) {
					v = p.verts[j] * 3;
					DetourCommon.vMin(bmin, tile.data.verts, v);
					DetourCommon.vMax(bmax, tile.data.verts, v);
				}
				if (overlapBounds(qmin, qmax, bmin, bmax)) {
					polys.push(NavMesh.or(base, i));
				}
			}
			return polys;
		}
	}

	/// Adds a tile to the navigation mesh.
	/// @param[in] data Data for the new tile mesh. (See: #dtCreateNavMeshData)
	/// @param[in] dataSize Data size of the new tile mesh.
	/// @param[in] flags Tile flags. (See: #dtTileFlags)
	/// @param[in] lastRef The desired reference for the tile. (When reloading a
	/// tile.) [opt] [Default: 0]
	/// @param[out] result The tile reference. (If the tile was succesfully
	/// added.) [opt]
	/// @return The status flags for the operation.
	/// @par
	///
	/// The add operation will fail if the data is in the wrong format, the
	/// allocated tile
	/// space is full, or there is a tile already at the specified reference.
	///
	/// The lastRef parameter is used to restore a tile with the same tile
	/// reference it had previously used. In this case the #dtPolyRef's for the
	/// tile will be restored to the same values they were before the tile was
	/// removed.
	///
	/// The nav mesh assumes exclusive access to the data passed and will make
	/// changes to the dynamic portion of the data. For that reason the data
	/// should not be reused in other nav meshes until the tile has been successfully
	/// removed from this nav mesh.
	///
	/// @see dtCreateNavMeshData, #removeTile
	addTile(data, flags, lastRef) {
		// Make sure the data is in right format.
		let header = data.header;

		// Make sure the location is free.
		if (this.getTileAt(header.x, header.y, header.layer) != null)
			throw new RuntimeException("Tile already exists");

		// Allocate a tile.
		let tile = null;
		if (lastRef == 0) {
			if (this.m_nextFree != null) {
				tile = this.m_nextFree;
				this.m_nextFree = tile.next;
				tile.next = null;
				this.m_tileCount++;
			}
		} else {
			// Try to relocate the tile to specific index with same salt.
			let tileIndex = decodePolyIdTile(lastRef);
			if (tileIndex >= this.m_maxTiles)
				throw new RuntimeException("Tile index too high");
			// Try to find the specific tile id from the free list.
			let target = this.m_tiles[tileIndex];
			let prev = null;
			tile = this.m_nextFree;
			while (tile != null && tile != target) {
				prev = tile;
				tile = tile.next;
			}
			// Could not find the correct location.
			if (tile != target)
				throw new RuntimeException("Could not find tile");
			// Remove from freelist
			if (prev == null)
				this.m_nextFree = tile.next;
			else
				prev.next = tile.next;

			// Restore salt.
			tile.salt = decodePolyIdSalt(lastRef);
		}

		// Make sure we could allocate a tile.
		if (tile == null)
			throw new RuntimeException("Could not allocate a tile");

		tile.data = data;
		tile.flags = flags;
		tile.links = [];

		// Insert tile into the position lut.
		let h = NavMesh.computeTileHash(header.x, header.y, this.m_tileLutMask);
		tile.next = this.m_posLookup[h];
		this.m_posLookup[h] = tile;

		// Patch header pointers.

		// If there are no items in the bvtree, reset the tree pointer.
		if (tile.data.bvTree != null && tile.data.bvTree.length == 0)
			tile.data.bvTree = null;

		// Init tile.

		this.connectIntLinks(tile);
		// Base off-mesh connections to their starting polygons and connect connections inside the tile.
		this.baseOffMeshLinks(tile);
		this.connectExtOffMeshLinks(tile, tile, -1);

		// Connect with layers in current tile.
		let neis = this.getTilesAt(header.x, header.y);
		for (let j = 0; j < neis.length; ++j) {
			if (neis[j] == tile) {
				continue;
			}
			connectExtLinks(tile, neis[j], -1);
			connectExtLinks(neis[j], tile, -1);
			this.connectExtOffMeshLinks(tile, neis[j], -1);
			this.connectExtOffMeshLinks(neis[j], tile, -1);
		}

		// Connect with neighbour tiles.
		for (let i = 0; i < 8; ++i) {
			neis = this.getNeighbourTilesAt(header.x, header.y, i);
			for (let j = 0; j < neis.length; ++j) {
				connectExtLinks(tile, neis[j], i);
				connectExtLinks(neis[j], tile, DetourCommon.oppositeTile(i));
				this.connectExtOffMeshLinks(tile, neis[j], i);
				this.connectExtOffMeshLinks(neis[j], tile, DetourCommon.oppositeTile(i));
			}
		}

		return this.getTileRef(tile);
	}

	/// Removes the specified tile from the navigation mesh.
	/// @param[in] ref The reference of the tile to remove.
	/// @param[out] data Data associated with deleted tile.
	/// @param[out] dataSize Size of the data associated with deleted tile.
	/// @return The status flags for the operation.
	// dtStatus removeTile(dtTileRef ref, char** data, int* dataSize);
	/// @par
	///
	/// This function returns the data for the tile so that, if desired,
	/// it can be added back to the navigation mesh at a later point.
	///
	/// @see #addTile
	removeTile(ref) {
		if (ref == 0) {
			return null;
		}
		let tileIndex = decodePolyIdTile(ref);
		let tileSalt = decodePolyIdSalt(ref);
		if (tileIndex >= this.m_maxTiles)
			throw new RuntimeException("Invalid tile index");
		let tile = this.m_tiles[tileIndex];
		if (tile.salt != tileSalt)
			throw new RuntimeException("Invalid tile salt");

		// Remove tile from hash lookup.
		let h = NavMesh.computeTileHash(tile.data.header.x, tile.data.header.y, this.m_tileLutMask);
		let prev = null;
		let cur = this.m_posLookup[h];
		while (cur != null) {
			if (cur == tile) {
				if (prev != null)
					prev.next = cur.next;
				else
					this.m_posLookup[h] = cur.next;
				break;
			}
			prev = cur;
			cur = cur.next;
		}

		// Remove connections to neighbour tiles.
		// Create connections with neighbour tiles.

		// Disconnect from other layers in current tile.
		nneis = this.getTilesAt(tile.data.header.x, tile.data.header.y);
		for (let j of nneis) {
			if (j == tile)
				continue;
			unconnectLinks(j, tile);
		}

		// Disconnect from neighbour tiles.
		for (let i = 0; i < 8; ++i) {
			nneis = this.getNeighbourTilesAt(tile.data.header.x, tile.data.header.y, i);
			for (let j of nneis)
				unconnectLinks(j, tile);
		}
		let data = tile.data;
		// Reset tile.
		tile.data = null;

		tile.flags = 0;
		tile.links = [];

		// Update salt, salt should never be zero.
		tile.salt = (tile.salt + 1) & ((1 << NavMesh.DT_SALT_BITS) - 1);
		if (tile.salt == 0)
			tile.salt++;

		// Add to free list.
		tile.next = this.m_nextFree;
		this.m_nextFree = tile;
		this.m_tileCount--;
		return data;
	}

	/// Builds internal polygons links for a tile.
	connectIntLinks(tile) {
		if (tile == null)
			return;

		let base = this.getPolyRefBase(tile);

		for (let i = 0; i < tile.data.header.polyCount; ++i) {
			let poly = tile.data.polys[i];
			poly.firstLink = NavMesh.DT_NULL_LINK;

			if (poly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION)
				continue;

			// Build edge links backwards so that the links will be
			// in the linked list from lowest index to highest.
			for (let j = poly.vertCount - 1; j >= 0; --j) {
				// Skip hard and non-internal edges.
				if (poly.neis[j] == 0 || (poly.neis[j] & NavMesh.DT_EXT_LINK) != 0)
					continue;

				let idx = this.allocLink(tile);
				let link = tile.links[idx];
				link.ref = NavMesh.or(base, (poly.neis[j] - 1));
				link.edge = j;
				link.side = 0xff;
				link.bmin = link.bmax = 0;
				// Add to linked list.
				link.next = poly.firstLink;
				poly.firstLink = idx;
			}
		}
	}
	unconnectLinks(tile, target) {
		if (tile == null || target == null)
			return;

		let targetNum = decodePolyIdTile(this.getTileRef(target));

		for (let i = 0; i < tile.data.header.polyCount; ++i) {
			let poly = tile.data.polys[i];
			let j = poly.firstLink;
			let pj = NavMesh.DT_NULL_LINK;
			while (j != NavMesh.DT_NULL_LINK) {
				if (decodePolyIdTile(tile.links[j].ref) == targetNum) {
					// Remove link.
					let nj = tile.links[j].next;
					if (pj == NavMesh.DT_NULL_LINK)
						poly.firstLink = nj;
					else
						tile.links[pj].next = nj;
					freeLink(tile, j);
					j = nj;
				} else {
					// Advance
					pj = j;
					j = tile.links[j].next;
				}
			}
		}
	}

	connectExtLinks(tile, target, side) {
		if (tile == null)
			return;

		// Connect border links.
		for (let i = 0; i < tile.data.header.polyCount; ++i) {
			let poly = tile.data.polys[i];

			// Create new links.
			// short m = NavMesh.DT_EXT_LINK | (short)side;

			let nv = poly.vertCount;
			for (let j = 0; j < nv; ++j) {
				// Skip non-portal edges.
				if ((poly.neis[j] & NavMesh.DT_EXT_LINK) == 0)
					continue;

				let dir = poly.neis[j] & 0xff;
				if (side != -1 && dir != side)
					continue;

				// Create new links
				let va = poly.verts[j] * 3;
				let vb = poly.verts[(j + 1) % nv] * 3;
				connectedPolys = findConnectingPolys(tile.data.verts, va, vb, target,
					DetourCommon.oppositeTile(dir), 4);
				nei = connectedPolys[0];
				neia = connectedPolys[1];
				let nnei = connectedPolys.third;
				for (let k = 0; k < nnei; ++k) {
					let idx = this.allocLink(tile);
					let link = tile.links[idx];
					link.ref = nei[k];
					link.edge = j;
					link.side = dir;

					link.next = poly.firstLink;
					poly.firstLink = idx;

					// Compress portal limits to a byte value.
					if (dir == 0 || dir == 4) {
						tmin = (neia[k * 2 + 0] - tile.data.verts[va + 2])
							/ (tile.data.verts[vb + 2] - tile.data.verts[va + 2]);
						tmax = (neia[k * 2 + 1] - tile.data.verts[va + 2])
							/ (tile.data.verts[vb + 2] - tile.data.verts[va + 2]);
						if (tmin > tmax) {
							temp = tmin;
							tmin = tmax;
							tmax = temp;
						}
						link.bmin = Math.floor(DetourCommon.clamp(tmin, 0.0, 1.0) * 255.0);
						link.bmax = Math.floor(DetourCommon.clamp(tmax, 0.0, 1.0) * 255.0);
					} else if (dir == 2 || dir == 6) {
						tmin = (neia[k * 2 + 0] - tile.data.verts[va])
							/ (tile.data.verts[vb] - tile.data.verts[va]);
						tmax = (neia[k * 2 + 1] - tile.data.verts[va])
							/ (tile.data.verts[vb] - tile.data.verts[va]);
						if (tmin > tmax) {
							temp = tmin;
							tmin = tmax;
							tmax = temp;
						}
						link.bmin = Math.floor(DetourCommon.clamp(tmin, 0.0, 1.0) * 255.0);
						link.bmax = Math.floor(DetourCommon.clamp(tmax, 0.0, 1.0) * 255.0);
					}
				}
			}
		}
	}

	connectExtOffMeshLinks(tile, target, side) {
		if (tile == null)
			return;

		// Connect off-mesh links.
		// We are interested on links which land from target tile to this tile.
		let oppositeSide = (side == -1) ? 0xff : DetourCommon.oppositeTile(side);

		for (let i = 0; i < target.data.header.offMeshConCount; ++i) {
			let targetCon = target.data.offMeshCons[i];
			if (targetCon.side != oppositeSide)
				continue;

			let targetPoly = target.data.polys[targetCon.poly];
			// Skip off-mesh connections which start location could not be
			// connected at all.
			if (targetPoly.firstLink == NavMesh.DT_NULL_LINK)
				continue;

			let ext = [targetCon.rad, target.data.header.walkableClimb, targetCon.rad];

			// Find polygon to connect to.
			let p = new Array(3);
			p[0] = targetCon.pos[3];
			p[1] = targetCon.pos[4];
			p[2] = targetCon.pos[5];
			let nearest = this.findNearestPolyInTile(tile, p, ext);
			let ref = nearest.getNearestRef();
			if (ref == 0)
				continue;
			let nearestPt = nearest.getNearestPos();
			// findNearestPoly may return too optimistic results, further check
			// to make sure.

			if (DetourCommon.sqr(nearestPt[0] - p[0]) + DetourCommon.sqr(nearestPt[2] - p[2]) > DetourCommon.sqr(targetCon.rad))
				continue;
			// Make sure the location is on curren mesh.
			target.data.verts[targetPoly.verts[1] * 3] = nearestPt[0];
			target.data.verts[targetPoly.verts[1] * 3 + 1] = nearestPt[1];
			target.data.verts[targetPoly.verts[1] * 3 + 2] = nearestPt[2];

			// let off-mesh connection to target poly.
			let idx = this.allocLink(target);
			let link = target.links[idx];
			link.ref = ref;
			link.edge = 1;
			link.side = oppositeSide;
			link.bmin = link.bmax = 0;
			// Add to linked list.
			link.next = targetPoly.firstLink;
			targetPoly.firstLink = idx;

			// let target poly to off-mesh connection.
			if ((targetCon.flags & NavMesh.DT_OFFMESH_CON_BIDIR) != 0) {
				let tidx = this.allocLink(tile);
				let landPolyIdx = NavMesh.decodePolyIdPoly(ref);
				let landPoly = tile.data.polys[landPolyIdx];
				link = tile.links[tidx];
				link.ref = NavMesh.or(this.getPolyRefBase(target), (targetCon.poly));
				link.edge = 0xff;
				link.side = (side == -1 ? 0xff : side);
				link.bmin = link.bmax = 0;
				// Add to linked list.
				link.next = landPoly.firstLink;
				landPoly.firstLink = tidx;
			}
		}
	}

	findConnectingPolys(verts, va, vb, tile, side, maxcon) {
		if (tile == null)
			return [null, null, 0];
		let con = new Array(maxcon);
		let conarea = new Array(maxcon * 2);
		let amin = new Array(2);
		let amax = new Array(2);
		calcSlabEndPoints(verts, va, vb, amin, amax, side);
		apos = getSlabCoord(verts, va, side);

		// Remove links pointing to 'side' and compact the links array.
		let bmin = new Array(2);
		let bmax = new Array(2);
		let m = NavMesh.or(NavMesh.DT_EXT_LINK, side);
		let n = 0;
		let base = this.getPolyRefBase(tile);

		for (let i = 0; i < tile.data.header.polyCount; ++i) {
			let poly = tile.data.polys[i];
			let nv = poly.vertCount;
			for (let j = 0; j < nv; ++j) {
				// Skip edges which do not poPoly to the right side.
				if (poly.neis[j] != m)
					continue;
				let vc = poly.verts[j] * 3;
				let vd = poly.verts[(j + 1) % nv] * 3;
				bpos = getSlabCoord(tile.data.verts, vc, side);
				// Segments are not close enough.
				if (Math.abs(apos - bpos) > 0.01)
					continue;

				// Check if the segments touch.
				calcSlabEndPoints(tile.data.verts, vc, vd, bmin, bmax, side);

				if (!overlapSlabs(amin, amax, bmin, bmax, 0.01, tile.data.header.walkableClimb))
					continue;

				// Add return value.
				if (n < maxcon) {
					conarea[n * 2 + 0] = Math.max(amin[0], bmin[0]);
					conarea[n * 2 + 1] = Math.min(amax[0], bmax[0]);
					con[n] = NavMesh.or(base, i);
					n++;
				}
				break;
			}
		}
		return [con, conarea, n];
	}

	static getSlabCoord(verts, va, side) {
		if (side == 0 || side == 4)
			return verts[va];
		else if (side == 2 || side == 6)
			return verts[va + 2];
		return 0;
	}

	static calcSlabEndPoints(verts, va, vb, bmin, bmax, side) {
		if (side == 0 || side == 4) {
			if (verts[va + 2] < verts[vb + 2]) {
				bmin[0] = verts[va + 2];
				bmin[1] = verts[va + 1];
				bmax[0] = verts[vb + 2];
				bmax[1] = verts[vb + 1];
			} else {
				bmin[0] = verts[vb + 2];
				bmin[1] = verts[vb + 1];
				bmax[0] = verts[va + 2];
				bmax[1] = verts[va + 1];
			}
		} else if (side == 2 || side == 6) {
			if (verts[va + 0] < verts[vb + 0]) {
				bmin[0] = verts[va + 0];
				bmin[1] = verts[va + 1];
				bmax[0] = verts[vb + 0];
				bmax[1] = verts[vb + 1];
			} else {
				bmin[0] = verts[vb + 0];
				bmin[1] = verts[vb + 1];
				bmax[0] = verts[va + 0];
				bmax[1] = verts[va + 1];
			}
		}
	}

	overlapSlabs(amin, amax, bmin, bmax, px, py) {
		// Check for horizontal overlap.
		// The segment is shrunken a little so that slabs which touch
		// at end points are not connected.
		minx = Math.max(amin[0] + px, bmin[0] + px);
		maxx = Math.min(amax[0] - px, bmax[0] - px);
		if (minx > maxx)
			return false;

		// Check vertical overlap.
		ad = (amax[1] - amin[1]) / (amax[0] - amin[0]);
		ak = amin[1] - ad * amin[0];
		bd = (bmax[1] - bmin[1]) / (bmax[0] - bmin[0]);
		bk = bmin[1] - bd * bmin[0];
		aminy = ad * minx + ak;
		amaxy = ad * maxx + ak;
		bminy = bd * minx + bk;
		bmaxy = bd * maxx + bk;
		dmin = bminy - aminy;
		dmax = bmaxy - amaxy;

		// Crossing segments always overlap.
		if (dmin * dmax < 0)
			return true;

		// Check for overlap at endpoints.
		thr = (py * 2) * (py * 2);
		if (dmin * dmin <= thr || dmax * dmax <= thr)
			return true;

		return false;
	}

	/**
	 * Builds internal polygons links for a tile.
	 * 
	 * @param tile
	 */
	baseOffMeshLinks(tile) {
		if (tile == null)
			return;

		let base = this.getPolyRefBase(tile);

		// Base off-mesh connection start points.
		for (let i = 0; i < tile.data.header.offMeshConCount; ++i) {
			let con = tile.data.offMeshCons[i];
			let poly = tile.data.polys[con.poly];

			let ext = [con.rad, tile.data.header.walkableClimb, con.rad];

			// Find polygon to connect to.
			let nearestPoly = this.findNearestPolyInTile(tile, con.pos, ext);
			let ref = nearestPoly.getNearestRef();
			if (ref == 0)
				continue;
			let p = con.pos; // First vertex
			let nearestPt = nearestPoly.getNearestPos();
			// findNearestPoly may return too optimistic results, further check
			// to make sure.
			if (DetourCommon.sqr(nearestPt[0] - p[0]) + DetourCommon.sqr(nearestPt[2] - p[2]) > DetourCommon.sqr(con.rad))
				continue;
			// Make sure the location is on current mesh.
			tile.data.verts[poly.verts[0] * 3] = nearestPt[0];
			tile.data.verts[poly.verts[0] * 3 + 1] = nearestPt[1];
			tile.data.verts[poly.verts[0] * 3 + 2] = nearestPt[2];

			// let off-mesh connection to target poly.
			let idx = this.allocLink(tile);
			let link = tile.links[idx];
			link.ref = ref;
			link.edge = 0;
			link.side = 0xff;
			link.bmin = link.bmax = 0;
			// Add to linked list.
			link.next = poly.firstLink;
			poly.firstLink = idx;

			// Start end-poPoly is always connect back to off-mesh connection.
			let tidx = this.allocLink(tile);
			let landPolyIdx = NavMesh.decodePolyIdPoly(ref);
			let landPoly = tile.data.polys[landPolyIdx];
			link = tile.links[tidx];
			link.ref = NavMesh.or(base, (con.poly));
			link.edge = 0xff;
			link.side = 0xff;
			link.bmin = link.bmax = 0;
			// Add to linked list.
			link.next = landPoly.firstLink;
			landPoly.firstLink = tidx;
		}
	}

	/**
	 * Returns closest poPoly on polygon.
	 * 
	 * @param ref
	 * @param pos
	 * @return
	 */
	closestPointOnPoly(ref, pos) {
		let tileAndPoly = this.getTileAndPolyByRefUnsafe(ref);
		let tile = tileAndPoly[0];
		let poly = tileAndPoly[1];
		// Off-mesh connections don't have detail polygons.
		if (poly.getType() == Poly.DT_POLYTYPE_OFFMESH_CONNECTION) {
			let v0 = poly.verts[0] * 3;
			let v1 = poly.verts[1] * 3;
			let d0 = DetourCommon.vDist3(pos, tile.data.verts, v0);
			let d1 = DetourCommon.vDist3(pos, tile.data.verts, v1);
			let u = d0 / (d0 + d1);
			let closest = DetourCommon.vLerp4(tile.data.verts, v0, v1, u);
			return new ClosestPointOnPolyResult(false, closest);
		}

		// Clamp poPoly to be inside the polygon.
		let verts = new Array(this.m_maxVertPerPoly * 3);
		let edged = new Array(this.m_maxVertPerPoly);
		let edget = new Array(this.m_maxVertPerPoly);
		let nv = poly.vertCount;
		for (let i = 0; i < nv; ++i)
			arraycopy(tile.data.verts, poly.verts[i] * 3, verts, i * 3, 3);

		let posOverPoly = false;
		let closest = new Array(3);
		DetourCommon.vCopy(closest, pos);
		if (!DetourCommon.distancePtPolyEdgesSqr(pos, verts, nv, edged, edget)) {
			// PoPoly is outside the polygon, dtClamp to nearest edge.
			let dmin = edged[0];
			let imin = 0;
			for (let i = 1; i < nv; ++i) {
				if (edged[i] < dmin) {
					dmin = edged[i];
					imin = i;
				}
			}
			let va = imin * 3;
			let vb = ((imin + 1) % nv) * 3;
			closest = DetourCommon.vLerp4(verts, va, vb, edget[imin]);
			posOverPoly = false;
		} else {
			posOverPoly = true;
		}

		// Find height at the location.
		let ip = poly.index;
		if (tile.data.detailMeshes != null && tile.data.detailMeshes.length > ip) {
			let pd = tile.data.detailMeshes[ip];
			for (let j = 0; j < pd.triCount; ++j) {
				let t = (pd.triBase + j) * 4;
				let v = [];//Was [3][]
				for (let k = 0; k < 3; ++k) {
					if (tile.data.detailTris[t + k] < poly.vertCount) {
						let index = poly.verts[tile.data.detailTris[t + k]] * 3;
						v[k] = [
							tile.data.verts[index], tile.data.verts[index + 1],
							tile.data.verts[index + 2]
						];
					} else {
						let index = (pd.vertBase + (tile.data.detailTris[t + k] - poly.vertCount)) * 3;
						v[k] = [
							tile.data.detailVerts[index], tile.data.detailVerts[index + 1],
							tile.data.detailVerts[index + 2]
						];
					}
				}
				let heightResult = DetourCommon.closestHeightPointTriangle(closest, v[0], v[1], v[2]);
				if (heightResult[0]) {
					closest[1] = heightResult[1];
					break;
				}
			}
		}
		return new ClosestPointOnPolyResult(posOverPoly, closest);
	}

	findNearestPolyInTile(tile, center, extents) {
		let nearestPt = null;
		let bmin = DetourCommon.vSub(center, extents);
		let bmax = DetourCommon.vAdd(center, extents);

		// Get nearby polygons from proximity grid.
		let polys = this.queryPolygonsInTile(tile, bmin, bmax);

		// Find nearest polygon amongst the nearby polygons.
		let nearest = 0;
		let nearestDistanceSqr = Number.MAX_VALUE;
		for (let i = 0; i < polys.length; ++i) {
			let ref = polys[i];
			let d = 0;
			let cpp = this.closestPointOnPoly(ref, center);
			let posOverPoly = cpp.isPosOverPoly();
			let closestPtPoly = cpp.getClosest();

			// If a poPoly is directly over a polygon and closer than
			// climb height, favor that instead of straight line nearest point.
			let diff = DetourCommon.vSub(center, closestPtPoly);
			if (posOverPoly) {
				d = Math.abs(diff[1]) - tile.data.header.walkableClimb;
				d = d > 0 ? d * d : 0;
			} else {
				d = DetourCommon.vLenSqr(diff);
			}
			if (d < nearestDistanceSqr) {
				nearestPt = closestPtPoly;
				nearestDistanceSqr = d;
				nearest = ref;
			}
		}
		return new FindNearestPolyResult(nearest, nearestPt);
	}

	getTileAt(x, y, layer) {
		// Find tile based on hash.
		let h = NavMesh.computeTileHash(x, y, this.m_tileLutMask);
		let tile = this.m_posLookup[h];
		while (tile != null) {
			if (tile.data.header != null && tile.data.header.x == x && tile.data.header.y == y
				&& tile.data.header.layer == layer) {
				return tile;
			}
			tile = tile.next;
		}
		return null;
	}

	getNeighbourTilesAt(x, y, side) {
		let nx = x, ny = y;
		switch (side) {
			case 0:
				nx++;
				break;
			case 1:
				nx++;
				ny++;
				break;
			case 2:
				ny++;
				break;
			case 3:
				nx--;
				ny++;
				break;
			case 4:
				nx--;
				break;
			case 5:
				nx--;
				ny--;
				break;
			case 6:
				ny--;
				break;
			case 7:
				nx++;
				ny--;
				break;
		}
		return this.getTilesAt(nx, ny);
	}

	getTilesAt(x, y) {
		let tiles = [];
		// Find tile based on hash.
		let h = NavMesh.computeTileHash(x, y, this.m_tileLutMask);
		let tile = this.m_posLookup[h];
		while (tile != null) {
			if (tile.data.header != null && tile.data.header.x == x && tile.data.header.y == y) {
				tiles.push(tile);
			}
			tile = tile.next;
		}
		return tiles;
	}

	getTileRefAt(x, y, layer) {
		// Find tile based on hash.
		let h = NavMesh.computeTileHash(x, y, this.m_tileLutMask);
		let tile = this.m_posLookup[h];
		while (tile != null) {
			if (tile.data.header != null && tile.data.header.x == x && tile.data.header.y == y
				&& tile.data.header.layer == layer) {
				return this.getTileRef(tile);
			}
			tile = tile.next;
		}
		return 0;
	}

	getTileByRef(ref) {
		if (ref == 0)
			return null;
		let tileIndex = decodePolyIdTile(ref);
		let tileSalt = decodePolyIdSalt(ref);
		if (tileIndex >= this.m_maxTiles)
			return null;
		let tile = this.m_tiles[tileIndex];
		if (tile.salt != tileSalt)
			return null;
		return tile;
	}

	getTileRef(tile) {
		if (tile == null)
			return 0;
		let it = tile.index;
		return NavMesh.encodePolyId(tile.salt, it, 0);
	}

	static computeTileHash(x, y, mask) {
		let h1 = 0x8da6b343; // Large multiplicative constants;
		let h2 = 0xd8163841; // here arbitrarily chosen primes
		let n = h1 * x + h2 * y;
		return n & mask;
	}

	/// @par
	///
	/// Off-mesh connections are stored in the navigation mesh as special
	/// 2-vertex
	/// polygons with a single edge. At least one of the vertices is expected to
	/// be
	/// inside a normal polygon. So an off-mesh connection is "entered" from a
	/// normal polygon at one of its endpoints. This is the polygon identified
	/// by
	/// the prevRef parameter.
	getOffMeshConnectionPolyEndPoints(prevRef, polyRef) {
		if (polyRef == 0)
			throw new IllegalArgumentException("polyRef = 0");

		// Get current polygon
		let saltitip = NavMesh.decodePolyId(polyRef);
		let salt = saltitip[0];
		let it = saltitip[1];
		let ip = saltitip[2];
		if (it >= this.m_maxTiles) {
			throw new IllegalArgumentException("Invalid tile ID > max tiles");
		}
		if (this.m_tiles[it].salt != salt || this.m_tiles[it].data.header == null) {
			throw new IllegalArgumentException("Invalid salt or missing tile header");
		}
		let tile = this.m_tiles[it];
		if (ip >= tile.data.header.polyCount) {
			throw new IllegalArgumentException("Invalid poly ID > poly count");
		}
		let poly = tile.data.polys[ip];

		// Make sure that the current poly is indeed off-mesh link.
		if (poly.getType() != Poly.DT_POLYTYPE_OFFMESH_CONNECTION)
			throw new IllegalArgumentException("Invalid poly type");

		// Figure out which way to hand out the vertices.
		let idx0 = 0, idx1 = 1;

		// Find link that points to first vertex.
		for (let i = poly.firstLink; i != NavMesh.DT_NULL_LINK; i = tile.links[i].next) {
			if (tile.links[i].edge == 0) {
				if (tile.links[i].ref != prevRef) {
					idx0 = 1;
					idx1 = 0;
				}
				break;
			}
		}
		let startPos = new Array(3);
		let endPos = new Array(3);
		DetourCommon.vCopy(startPos, tile.data.verts, poly.verts[idx0] * 3);
		DetourCommon.vCopy(endPos, tile.data.verts, poly.verts[idx1] * 3);
		return [startPos, endPos];

	}

	getMaxVertsPerPoly() {
		return this.m_maxVertPerPoly;
	}

	getTileCount() {
		return this.m_tileCount;
	}
}

export default NavMesh;