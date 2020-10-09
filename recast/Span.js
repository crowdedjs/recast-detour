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

/** Represents a span in a heightfield. */
class Span {

	/** The lower limit of the span. [Limit: < #smax] */
smin = 0;
	/** The upper limit of the span. [Limit: <= #RC_SPAN_MAX_HEIGHT] */
smax = 0;
	/** The area id assigned to the span. */
area = 0;
	/** The next span higher up in column. */
next = null;

}

export default Span;