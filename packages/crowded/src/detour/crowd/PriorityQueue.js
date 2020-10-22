//https://truetocode.com/binary-treemax-heap-priority-queue-and-implementation-using-javascript/427/
class PriorityQueue {
  constructor(comparator) {
    this.comparator = comparator;
    this.array = [];
    
  }
  isEmpty(){
    return this.array.length == 0;
  }
  
  poll(){
    return this.array.splice(0,1)[0];
  }
  insert(element){
    this.push(element);
  }
  push(element){
    this.array.push(element);
    this.array.sort(this.comparator);
  }
  offer(element){
    this.push(element);
  }
  remove(element){
    let index =  this.array.indexOf(element)
    if(index >= 0)
      this.array.splice(index,1);
  }

  
}

export default PriorityQueue;