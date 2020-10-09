//https://truetocode.com/binary-treemax-heap-priority-queue-and-implementation-using-javascript/427/
class PriorityQueue {
  constructor(comparator) {
    this.comparator = comparator;
    this.array = [null];
    this.size = 0;
    //this.array = [null, ...array];
    //this.size = array.length
  }
  isEmpty(){
    return this.size == 0;
  }
  arrange(idx) {
    while (idx > 1 && this.compare(Math.floor(idx / 2), idx)) {

      this.swap(idx, Math.floor(idx / 2));
      idx = Math.floor(idx / 2);
    }
  }

  heaper(idx1) {
    while (2 * idx1 <= this.size) {

      let idx2 = 2 * idx1;

      if (idx2 < this.size && this.compare(idx2, idx2 + 1)) idx2++;

      if (!this.compare(idx1, idx2)) break;


      this.swap(idx1, idx2);

      idx1 = idx2;
    }
  }
  poll(){
    return this.rootdelete();
  }
  rootdelete() {

    let max = this.array[1]

    this.swap(1, this.size--);
    this.heaper(1);

    this.array[this.size + 1] = null;

    return max;
  }
  push(element){
    this.insert(element);
  }

  insert(element) {

    this.array[++this.size] = element;

    this.arrange(this.size);
  }

  compare(idx1, idx2) {
    return this.comparator(this.array[idx1], this.array[idx2]);
    // return this.array[idx1].priority < this.array[idx2].priority;
  }

  swap(idx1, idx2) {
    const temp = this.array[idx1];
    this.array[idx1] = this.array[idx2];
    this.array[idx2] = temp;
  }
  returnall() {
    return this.array;
  }
}

export default PriorityQueue;