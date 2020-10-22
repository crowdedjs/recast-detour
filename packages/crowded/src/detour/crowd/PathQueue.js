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

import PathQuery from "./PathQuery.js"
import NavMeshQuery from "../NavMeshQuery.js"
import DetourCommon from "../DetourCommon.js"
import Status from "../Status.js";
import FindPathResult from "../FindPathResult.js"


class PathQueue {

	static MAX_QUEUE = 8;
	static DT_PATHQ_INVALID = 0;
	static MAX_KEEP_ALIVE = 2; // in update ticks.

	m_queue = new Array(PathQueue.MAX_QUEUE);
	m_nextHandle = 1;
	m_queueHead;
	m_navquery;

	constructor(maxSearchNodeCount, nav) {
		this.m_navquery = new NavMeshQuery(nav);
		for (let i = 0; i < PathQueue.MAX_QUEUE; ++i) {
			this.m_queue[i] = new PathQuery();
			this.m_queue[i].ref = PathQueue.DT_PATHQ_INVALID;
			this.m_queue[i].path = new Array(256);
		}
		this.m_queueHead = 0;
	}

	update(maxIters) {
		// Update path request until there is nothing to update
		// or upto maxIters pathfinder iterations has been consumed.
		let iterCount = maxIters;

		for (let i = 0; i < PathQueue.MAX_QUEUE; ++i) {
			let q = this.m_queue[this.m_queueHead % PathQueue.MAX_QUEUE];

			// Skip inactive requests.
			if (q.ref == PathQueue.DT_PATHQ_INVALID) {
				this.m_queueHead++;
				continue;
			}

			// Handle compPolyed request.
			if (q.status != null && (q.status == Status.SUCCSESS || q.status == Status.PARTIAL_RESULT || q.status == Status.FAILURE)) {
				// If the path result has not been read in few frames, free the slot.
				q.keepAlive++;
				if (q.keepAlive > PathQueue.MAX_KEEP_ALIVE) {
					q.ref = PathQueue.DT_PATHQ_INVALID;
					q.status = null;
				}

				this.m_queueHead++;
				continue;
			}

			// Handle query start.
			if (q.status == null) {
				q.status = this.m_navquery.initSlicedFindPath(q.startRef, q.endRef, q.startPos, q.endPos, q.filter, 0);
			}
			// Handle query in progress.
			if (q.status == Status.IN_PROGRESS) {
				let iters = 0;
				let res = this.m_navquery.updateSlicedFindPath(iterCount);
				iters = res.getIterations();
				q.status = res.getStatus();
				iterCount -= iters;
			}
			if (q.status == Status.SUCCSESS || q.status == Status.PARTIAL_RESULT) {
				let path = this.m_navquery.finalizeSlicedFindPath();
				q.status = path.getStatus();
				q.path = path.getRefs();
			}

			if (iterCount <= 0)
				break;

			this.m_queueHead++;
		}

	}

	request(startRef, endRef, startPos, endPos, filter) {
		// Find empty slot
		let slot = -1;
		for (let i = 0; i < PathQueue.MAX_QUEUE; ++i) {
			if (this.m_queue[i].ref == PathQueue.DT_PATHQ_INVALID) {
				slot = i;
				break;
			}
		}
		// Could not find slot.
		if (slot == -1)
			return PathQueue.DT_PATHQ_INVALID;

		let ref = this.m_nextHandle++;
		if (this.m_nextHandle == PathQueue.DT_PATHQ_INVALID)
			this.m_nextHandle++;

		let q = this.m_queue[slot];
		q.ref = ref;
		DetourCommon.vCopy(q.startPos, startPos);
		q.startRef = startRef;
		DetourCommon.vCopy(q.endPos, endPos);
		q.endRef = endRef;
		q.status = null;
		q.filter = filter;
		q.keepAlive = 0;
		return ref;

	}

	getRequestStatus(ref) {
		for (let i = 0; i < PathQueue.MAX_QUEUE; ++i) {
			if (this.m_queue[i].ref == ref)
				return this.m_queue[i].status;
		}
		return Status.FAILURE;

	}

	getPathResult(ref) {
		for (let i = 0; i < PathQueue.MAX_QUEUE; ++i) {
			if (this.m_queue[i].ref == ref) {
				let q = this.m_queue[i];
				// Free request for reuse.
				q.ref = PathQueue.DT_PATHQ_INVALID;
				q.status = null;
				return new FindPathResult(Status.SUCCSESS, q.path);
			}
		}
		return new FindPathResult(Status.FAILURE, null);
	}
}

export default PathQueue;
