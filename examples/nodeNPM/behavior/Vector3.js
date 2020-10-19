
/**
 * Created by bricksphd on 6/19/2018.
 */
class Vector3 {

    static Zero = new Vector3(0, 0, 0);

    x = 0;
    y = 0;
    z = 0;

    constructor(x, y, z) {
        if (!y || !z) {
            this.x = x[0];
            this.y = x[1];
            this.z = x[2];
        }
        else {
            this.x = x;
            this.y = y;
            this.z = z;
        }
    }

    asArray() {
        return new [this.x, this.y, this.z];
    }
    toString() {
        return "" + this.x + "," + this.y + "," + this.z;
    }

    static subtract(end, start) {
        return new Vector3(end.x - start.x, end.y - start.y, end.z - start.z);
    }

    normalize() {
        let length = length();

        this.x /= length;
        this.y /= length;
        this.z /= length;

        return this;
    }

    length() {
        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    }

    scale(f) {
        this.x *= f;
        this.y *= f;
        this.z *= f;

        return this;
    }

    scaleXZ(f) {
        this.x *= f;
        this.z *= f;
        return this;
    }

    add(other) {
        this.x += other.x;
        this.y += other.y;
        this.z += other.z;

        return this;
    }
}

module.exports = Vector3