// this list contains all the objects that will collide and the time of collision
// it also keeps track of the statistics of the collisions
export class collisionClass {

    constructor() {

        this.list = []
        this.count = 0
        this.period = 0.0
    }

    resetCount() {
        this.count = 0
        this.period = 0.0
    }

    resetList() {
        this.list = []
    }

    add( a, b, t) {

        this.list.push({a,b,t})
    }

    collide(interval) { 

        let list = this.list
        
        // check
        if (list.length < 1) return

        // sort the list 
        list.sort( (r1, r2) => r1.t - r2.t )

        // do all collisions
        let len = list.length
        let col = null
        for (let i=0; i < len ; i++) {

            // get the collision record
            col = list[i]

            // parameters are the other object, the collision time (relative to the start if the interval) and the interval duration
            col.a.collide(col.b, col.t, interval)
        }

        this.period += interval
        this.count += len
    }

    get() {
        return {count: this.count, period: this.period}
    }
}