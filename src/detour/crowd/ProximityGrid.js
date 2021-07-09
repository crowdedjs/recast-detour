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




class ProximityGrid {

	static ItemKey = class ItemKey {

		x;
		y;

		constructor(x, y) {
			this.x = x;
			this.y = y;
		}

		hashCode() {
			prime = 31;
			result = 1;
			result = prime * result + x;
			result = prime * result + y;
			return result;
		}

		equals(obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			other = obj;
			if (x != other.x)
				return false;
			if (y != other.y)
				return false;
			return true;
		}

	};

	m_cellSize;
	m_invCellSize;
	items = {};
	m_bounds = new Array(4);

	constructor(m_cellSize, m_invCellSize) {
		this.m_cellSize = m_cellSize;
		this.m_invCellSize = m_invCellSize;
		this.items = {};
	}

	clear() {
		this.items = {};
		this.m_bounds[0] = 0xfff;
		this.m_bounds[1] = 0xfff;
		this.m_bounds[2] = -0xfff;
		this.m_bounds[3] = -0xfff;
	}

	addItem(id, minx, miny, maxx, maxy) {
		let iminx = Math.floor(minx * this.m_invCellSize);
		let iminy = Math.floor(miny * this.m_invCellSize);
		let imaxx = Math.floor(maxx * this.m_invCellSize);
		let imaxy = Math.floor(maxy * this.m_invCellSize);

		this.m_bounds[0] = Math.min(this.m_bounds[0], iminx);
		this.m_bounds[1] = Math.min(this.m_bounds[1], iminy);
		this.m_bounds[2] = Math.min(this.m_bounds[2], imaxx);
		this.m_bounds[3] = Math.min(this.m_bounds[3], imaxy);

		for (let y = iminy; y <= imaxy; ++y) {
			for (let x = iminx; x <= imaxx; ++x) {
				//let key = new ProximityGrid.ItemKey(x, y);
				let _string = this.stringIt(x,y);
				let ids = this.items[_string];
				if (ids == null) {
					ids = [];
					this.items[_string]= ids;
				}
				ids.push(id);
			}
		}
	}

	stringIt(x,y){
		return x + "," + y;
	}

	queryItems(minx, miny, maxx, maxy) {
		let iminx = Math.floor(minx * this.m_invCellSize);
		let iminy = Math.floor(miny * this.m_invCellSize);
		let imaxx = Math.floor(maxx * this.m_invCellSize);
		let imaxy = Math.floor(maxy * this.m_invCellSize);

		let result = new Set();
		for (let y = iminy; y <= imaxy; ++y) {
			for (let x = iminx; x <= imaxx; ++x) {
				let _string = this.stringIt(x,y);
				
				//let key = new ProximityGrid.ItemKey(x, y);
				let ids = this.items[_string];
				if (ids != null) {
					result.add(...ids);
				}
			}
		}

		return result;
	}

	getItemCountAt(x, y) {
		key = new ItemKey(x, y);
		ids = this.items[key];
		return ids != null ? ids.length : 0;
	}
}

export default ProximityGrid;