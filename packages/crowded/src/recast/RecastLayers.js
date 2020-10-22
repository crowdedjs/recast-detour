



class RecastLayers {

	static RC_MAX_LAYERS = RecastConstants.RC_NOT_CONNECTED;
	static RC_MAX_NEIS = 16;

	static LayerRegion = class LayerRegion {
		id;
		layerId;
		base;
		ymin;
		ymax;
		layers;
		neis;

		constructor(i) {
			this.id = i;
			this.ymin = 0xFFFF;
			this.layerId = 0xff;
			this.layers = [];
			this.neis = [];
		}

	};

	static addUnique(a, v) {
		if (!a.contains(v)) {
			a.push(v);
		}
	}

	static contains(a, v) {
		return a.contains(v);
	}

	static overlapRange(amin, amax, bmin, bmax) {
		return (amin > bmax || amax < bmin) ? false : true;
	}

	static buildHeightfieldLayers(ctx, chf, borderSize, walkableHeight) {

		ctx.startTimer("RC_TIMER_BUILD_LAYERS");
		let w = chf.width;
		let h = chf.height;
		let srcReg = new Array(chf.spanCount);
		Arrays.fill(srcReg, 0xFF);
		let nsweeps = chf.width;// Math.max(chf.width, chf.height);
		let sweeps = new SweepSpan[nsweeps];
		for (let i = 0; i < sweeps.length; i++) {
			sweeps[i] = new SweepSpan();
		}
		// Partition walkable area into monotone regions.
		let prevCount = new Array(256);
		let regId = 0;
		// Sweep one line at a time.
		for (let y = borderSize; y < h - borderSize; ++y) {
			// Collect spans from this row.
			Arrays.fill(prevCount, 0, regId, 0);
			let sweepId = 0;

			for (let x = borderSize; x < w - borderSize; ++x) {
				let c = chf.cells[x + y * w];

				for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
					let s = chf.spans[i];
					if (chf.areas[i] == RecastConstants.RC_NULL_AREA)
						continue;
					let sid = 0xFF;
					// -x

					if (GetCon(s, 0) != RecastConstants.RC_NOT_CONNECTED) {
						let ax = x + RecastCommon.GetDirOffsetX(0);
						let ay = y + RecastCommon.GetDirOffsetY(0);
						let ai = chf.cells[ax + ay * w].index + GetCon(s, 0);
						if (chf.areas[ai] != RecastConstants.RC_NULL_AREA && srcReg[ai] != 0xff)
							sid = srcReg[ai];
					}

					if (sid == 0xff) {
						sid = sweepId++;
						sweeps[sid].nei = 0xff;
						sweeps[sid].ns = 0;
					}

					// -y
					if (GetCon(s, 3) != RecastConstants.RC_NOT_CONNECTED) {
						let ax = x + RecastCommon.GetDirOffsetX(3);
						let ay = y + RecastCommon.GetDirOffsetY(3);
						let ai = chf.cells[ax + ay * w].index + GetCon(s, 3);
						let nr = srcReg[ai];
						if (nr != 0xff) {
							// Set neighbour when first valid neighbour is
							// encoutered.
							if (sweeps[sid].ns == 0)
								sweeps[sid].nei = nr;

							if (sweeps[sid].nei == nr) {
								// Update existing neighbour
								sweeps[sid].ns++;
								prevCount[nr]++;
							} else {
								// This is hit if there is nore than one
								// neighbour.
								// Invalidate the neighbour.
								sweeps[sid].nei = 0xff;
							}
						}
					}

					srcReg[i] = sid;
				}
			}

			// Create unique ID.
			for (let i = 0; i < sweepId; ++i) {
				// If the neighbour is set and there is only one continuous
				// connection to it,
				// the sweep will be merged with the previous one, else new
				// region is created.
				if (sweeps[i].nei != 0xff && prevCount[sweeps[i].nei] == sweeps[i].ns) {
					sweeps[i].id = sweeps[i].nei;
				} else {
					if (regId == 255) {
						throw new RuntimeException("rcBuildHeightfieldLayers: Region ID overflow.");
					}
					sweeps[i].id = regId++;
				}
			}

			// Remap local sweep ids to region ids.
			for (let x = borderSize; x < w - borderSize; ++x) {
				let c = chf.cells[x + y * w];
				for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
					if (srcReg[i] != 0xff)
						srcReg[i] = sweeps[srcReg[i]].id;
				}
			}
		}
		let nregs = regId;
		let regs = new LayerRegion[nregs];

		// Construct regions
		for (let i = 0; i < nregs; ++i) {
			regs[i] = new LayerRegion(i);
		}

		// Find region neighbours and overlapping regions.
		let lregs = [];
		for (let y = 0; y < h; ++y) {
			for (let x = 0; x < w; ++x) {
				let c = chf.cells[x + y * w];

				lregs = [];

				for (let i = c.index, ni = c.index + c.count; i < ni; ++i) {
					let s = chf.spans[i];
					let ri = srcReg[i];
					if (ri == 0xff)
						continue;

					regs[ri].ymin = Math.min(regs[ri].ymin, s.y);
					regs[ri].ymax = Math.max(regs[ri].ymax, s.y);

					// Collect all region layers.
					lregs.push(ri);

					// Update neighbours
					for (let dir = 0; dir < 4; ++dir) {
						if (GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
							let ax = x + RecastCommon.GetDirOffsetX(dir);
							let ay = y + RecastCommon.GetDirOffsetY(dir);
							let ai = chf.cells[ax + ay * w].index + GetCon(s, dir);
							let rai = srcReg[ai];
							if (rai != 0xff && rai != ri)
								addUnique(regs[ri].neis, rai);
						}
					}

				}

				// Update overlapping regions.
				for (let i = 0; i < lregs.length - 1; ++i) {
					for (let j = i + 1; j < lregs.length; ++j) {
						if (lregs[i].intValue() != lregs[j].intValue()) {
							let ri = regs[lregs[i]];
							let rj = regs[lregs[j]];
							addUnique(ri.layers, lregs[j]);
							addUnique(rj.layers, lregs[i]);
						}
					}
				}

			}
		}

		// Create 2D layers from regions.
		let layerId = 0;

		let stack = [];

		for (let i = 0; i < nregs; ++i) {
			let root = regs[i];
			// Skip already visited.
			if (root.layerId != 0xff)
				continue;

			// Start search.
			root.layerId = layerId;
			root.base = true;

			stack.push(i);

			while (!stack.length == 0) {
				// Pop front
				let reg = regs[stack.remove(0)];

				for (let nei of reg.neis) {
					let regn = regs[nei];
					// Skip already visited.
					if (regn.layerId != 0xff)
						continue;
					// Skip if the neighbour is overlapping root region.
					if (contains(root.layers, nei))
						continue;
					// Skip if the height range would become too large.
					let ymin = Math.min(root.ymin, regn.ymin);
					let ymax = Math.max(root.ymax, regn.ymax);
					if ((ymax - ymin) >= 255)
						continue;

					// Deepen
					stack.push(nei);

					// Mark layer id
					regn.layerId = layerId;
					// Merge current layers to root.
					for (let layer of regn.layers)
						addUnique(root.layers, layer);
					root.ymin = Math.min(root.ymin, regn.ymin);
					root.ymax = Math.max(root.ymax, regn.ymax);
				}
			}

			layerId++;
		}

		// Merge non-overlapping regions that are close in height.
		let mergeHeight = walkableHeight * 4;

		for (let i = 0; i < nregs; ++i) {
			let ri = regs[i];
			if (!ri.base)
				continue;

			let newId = ri.layerId;

			for ( ; ;) {
				let oldId = 0xff;

				for (let j = 0; j < nregs; ++j) {
					if (i == j)
						continue;
					let rj = regs[j];
					if (!rj.base)
						continue;

					// Skip if the regions are not close to each other.
					if (!overlapRange(ri.ymin, ri.ymax + mergeHeight, rj.ymin, rj.ymax + mergeHeight))
						continue;
					// Skip if the height range would become too large.
					let ymin = Math.min(ri.ymin, rj.ymin);
					let ymax = Math.max(ri.ymax, rj.ymax);
					if ((ymax - ymin) >= 255)
						continue;

					// Make sure that there is no overlap when merging 'ri' and
					// 'rj'.
					let overlap = false;
					// Iterate over all regions which have the same layerId as
					// 'rj'
					for (let k = 0; k < nregs; ++k) {
						if (regs[k].layerId != rj.layerId)
							continue;
						// Check if region 'k' is overlapping region 'ri'
						// Index to 'regs' is the same as region id.
						if (contains(ri.layers, k)) {
							overlap = true;
							break;
						}
					}
					// Cannot merge of regions overlap.
					if (overlap)
						continue;

					// Can merge i and j.
					oldId = rj.layerId;
					break;
				}

				// Could not find anything to merge with, stop.
				if (oldId == 0xff)
					break;

				// Merge
				for (let j = 0; j < nregs; ++j) {
					let rj = regs[j];
					if (rj.layerId == oldId) {
						rj.base = false;
						// Remap layerIds.
						rj.layerId = newId;
						// Add overlaid layers from 'rj' to 'ri'.
						for (let layer of rj.layers)
						addUnique(ri.layers, layer);
						// Update height bounds.
						ri.ymin = Math.min(ri.ymin, rj.ymin);
						ri.ymax = Math.max(ri.ymax, rj.ymax);
					}
				}
			}
		}

		// Compact layerIds
		let remap = new Array(256);

		// Find number of unique layers.
		layerId = 0;
		for (let i = 0; i < nregs; ++i)
			remap[regs[i].layerId] = 1;
		for (let i = 0; i < 256; ++i) {
			if (remap[i] != 0)
				remap[i] = layerId++;
			else
				remap[i] = 0xff;
		}
		// Remap ids.
		for (let i = 0; i < nregs; ++i)
			regs[i].layerId = remap[regs[i].layerId];

		// No layers, return empty.
		if (layerId == 0) {
			// ctx.stopTimer(RC_TIMER_BUILD_LAYERS);
			return null;
		}

		// Create layers.
		// rcAssert(lset.layers == 0);

		let lw = w - borderSize * 2;
		let lh = h - borderSize * 2;

		// Build contracted bbox for layers.
		let bmin = new Array(3);
		let bmax = new Array(3);
		copy(bmin, chf.bmin);
		copy(bmax, chf.bmax);
		bmin[0] += borderSize * chf.cs;
		bmin[2] += borderSize * chf.cs;
		bmax[0] -= borderSize * chf.cs;
		bmax[2] -= borderSize * chf.cs;

		let lset = new HeightfieldLayerSet();
		lset.layers = new HeightfieldLayer[layerId];
		for (let i = 0; i < lset.layers.length; i++) {
			lset.layers[i] = new HeightfieldLayer();
		}

		// Store layers.
		for (let i = 0; i < lset.layers.length; ++i) {
			let curId = i;

			let layer = lset.layers[i];

			let gridSize = lw * lh;

			layer.heights = new Array(gridSize);
			Arrays.fill(layer.heights, 0xFF);
			layer.areas = new Array(gridSize);
			layer.cons = new Array(gridSize);

			// Find layer height bounds.
			let hmin = 0, hmax = 0;
			for (let j = 0; j < nregs; ++j) {
				if (regs[j].base && regs[j].layerId == curId) {
					hmin = regs[j].ymin;
					hmax = regs[j].ymax;
				}
			}

			layer.width = lw;
			layer.height = lh;
			layer.cs = chf.cs;
			layer.ch = chf.ch;

			// Adjust the bbox to fit the heightfield.
			copy(layer.bmin, bmin);
			copy(layer.bmax, bmax);
			layer.bmin[1] = bmin[1] + hmin * chf.ch;
			layer.bmax[1] = bmin[1] + hmax * chf.ch;
			layer.hmin = hmin;
			layer.hmax = hmax;

			// Update usable data region.
			layer.minx = layer.width;
			layer.maxx = 0;
			layer.miny = layer.height;
			layer.maxy = 0;

			// Copy height and area from compact heightfield.
			for (let y = 0; y < lh; ++y) {
				for (let x = 0; x < lw; ++x) {
					let cx = borderSize + x;
					let cy = borderSize + y;
					let c = chf.cells[cx + cy * w];
					for (let j = c.index, nj = c.index + c.count; j < nj; ++j) {
						let s = chf.spans[j];
						// Skip unassigned regions.
						if (srcReg[j] == 0xff)
							continue;
						// Skip of does nto bePoly to current layer.
						let lid = regs[srcReg[j]].layerId;
						if (lid != curId)
							continue;

						// Update data bounds.
						layer.minx = Math.min(layer.minx, x);
						layer.maxx = Math.max(layer.maxx, x);
						layer.miny = Math.min(layer.miny, y);
						layer.maxy = Math.max(layer.maxy, y);

						// Store height and area type.
						let idx = x + y * lw;
						layer.heights[idx] = (char)(s.y - hmin);
						layer.areas[idx] = chf.areas[j];

						// Check connection.
						let portal = 0;
						let con = 0;
						for (let dir = 0; dir < 4; ++dir) {
							if (GetCon(s, dir) != RecastConstants.RC_NOT_CONNECTED) {
								let ax = cx + RecastCommon.GetDirOffsetX(dir);
								let ay = cy + RecastCommon.GetDirOffsetY(dir);
								let ai = chf.cells[ax + ay * w].index + GetCon(s, dir);
								let alid = srcReg[ai] != 0xff ? regs[srcReg[ai]].layerId : 0xff;
								// Portal mask
								if (chf.areas[ai] != RecastConstants.RC_NULL_AREA && lid != alid) {
									portal |= (1 << dir);
									// Update height so that it matches on both
									// sides of the portal.
									let as = chf.spans[ai];
									if (as.y > hmin)
										layer.heights[idx] = Math.max(layer.heights[idx], (char)(as.y - hmin));
								}
								// Valid connection mask
								if (chf.areas[ai] != RecastConstants.RC_NULL_AREA && lid == alid) {
									let nx = ax - borderSize;
									let ny = ay - borderSize;
									if (nx >= 0 && ny >= 0 && nx < lw && ny < lh)
										con |= (1 << dir);
								}
							}
						}
						layer.cons[idx] = (portal << 4) | con;
					}
				}
			}

			if (layer.minx > layer.maxx)
				layer.minx = layer.maxx = 0;
			if (layer.miny > layer.maxy)
				layer.miny = layer.maxy = 0;
		}

		// ctx=>stopTimer(RC_TIMER_BUILD_LAYERS);
		return lset;
	}
}
