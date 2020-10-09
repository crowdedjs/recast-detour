
class AreaModification {

static  RC_AREA_FLAGS_MASK = 0x3F;
value;
mask;

	/**
	 * Mask is set to all available bits, which means value is fully applied
	 * 
	 * @param value
	 *            The area id to apply. [Limit: <= #RC_AREA_FLAGS_MASK]
	 */
// constructor(value) {
// 		this.value = value;
// 		this.mask = AreaModification.RC_AREA_FLAGS_MASK;
// 	}

	/**
	 * 
	 * @param value
	 *            The area id to apply. [Limit: <= #RC_AREA_FLAGS_MASK]
	 * @param mask
	 *            Bitwise mask used when applying value. [Limit: <= #RC_AREA_FLAGS_MASK]
	 */
constructor(value,  mask = AreaModification.RC_AREA_FLAGS_MASK) {
		this.value = value;
		this.mask = mask;
	}

static clone = function(other) {
	return new AreaModification(other.value, other.mask);
		// this.value = other.value;
		// this.mask = other.mask;
	}

getMaskedValue() {
		return this.value & this.mask;
	}

apply(area) {
		return ((this.value & this.mask) | (area & ~this.mask));
	}
}

export default AreaModification
