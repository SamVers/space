import * as THREE from 'three/build/three.module.js'

// colors used for bounding boxes
const boxFreeClr = new THREE.Color("#00FF00")
const boxHitClr =  new THREE.Color("#FF0000")

export class objectClass {

    constructor(mass, radius) {

        this.pos = new THREE.Vector3()
        this.speed = new THREE.Vector3()
        this.radius = radius
        this.mass = mass

        // the bounding box when moving from one position to another
        this.bBox = new THREE.Box3

        // for debugging we sometimes need to see the bounding boxes
        this.boxHelper = null

        // objects can be added to a linked list (eg in the collisionlist)
        this.next = null
    }

    // calculate the size of the bounding box based on the new position q and the size of tha ball
    setBBox(q) {

        // notation
        let p = this.pos, r = this.radius, min = this.bBox.min, max = this.bBox.max

        // we subtract/add r to make sure the ball lies completely in the box
        min.x = p.x < q.x ? p.x-r : q.x-r
        min.y = p.y < q.y ? p.y-r : q.y-r
        min.z = p.z < q.z ? p.z-r : q.z-r

        max.x = p.x > q.x ? p.x+r : q.x+r
        max.y = p.y > q.y ? p.y+r : q.y+r
        max.z = p.z > q.z ? p.z+r : q.z+r

        // a box that has collided changes gradually to normal
        if (this.boxHelper) this.boxHelper.material.color.lerp(boxFreeClr,0.1)
    }

    newPosition(t) {

        // calculate the new position
        let q = this.pos.clone().addScaledVector(this.speed,t)

        // adjust the bounding box
        this.setBBox(q)

        // copy the new position
        this.pos.copy(q)
    }

    // a function to test for intersection of trajectories
    intersects( other, interval ) {

        if (this.bBox.intersectsBox(other.bBox)) {

            // we have a candidate - calculate relative position and speed
            // p = p2-p1, v= v2-v1 - note that p is already the new position t sec ago
            let p = other.pos.clone().sub(this.pos)
            let v = other.speed.clone().sub(this.speed)

            // calculate a,b and c of the quadratic distance equation
            let a = v.lengthSq()
            let b = 2*p.dot(v)
            let c = p.lengthSq() - (other.radius+this.radius)**2

            // the determinant
            let det = b*b - 4*a*c
            
            // if the determinant is positive there are intersections
            if ( det > 0) {

                // show that the objects will collide
                if (this.boxHelper) this.boxHelper.material.color.copy(boxHitClr)
                if (other.boxHelper) other.boxHelper.material.color.copy(boxHitClr)

                // calculate the roots - nota that t1 is the smallest root
                let sqrt = Math.sqrt(det)
                let t1 = (-b - sqrt)/(2*a)

                /*** 
                there is a second 'collision' at t2 when the balls intersect when moving apart
                t2 = (-b + sqrt)/(2*a)
                t2 = -tCol - b/a
                t2 > 0 becomes then t1 < -b/a
                so if we have missed the first intersection (at t1) but t2 is still in the future, we still consider t1
                ***/

                // the collision should have happened in the interval except when the second collision is still in the future
                if ( ( t1 < 0 ) && ((t1 > -interval) || (t1 < -b/a)) ) return t1
                else  return Infinity
            }
            else return Infinity
        }
        return Infinity
    }

    collide( other, tCol, interval ) {

        // bring the two objects back to the moment of the collision - tCol is always negative !
        this.pos.addScaledVector(this.speed, tCol)
        other.pos.addScaledVector(other.speed, tCol)

        // formula to be found in wikipedia https://en.wikipedia.org/wiki/Elastic_collision#Two-dimensional
        let p = this.pos.clone().sub(other.pos)
        let v = this.speed.clone().sub(other.speed)

        // the factors for the speed calculation
        let f1 =  ( 2 * other.mass / (this.mass + other.mass) ) * ( v.dot(p) / p.lengthSq() )
        let f2 = this.mass/other.mass

        // calculate the new speeds - keep the order !
        this.speed.sub( p.multiplyScalar(f1) )
        other.speed.add( p.multiplyScalar(f2) )

        // the two objects have still tCol seconds to travel with the new speed (-tCol is positive !)
        this.pos.addScaledVector(this.speed, -tCol)
        other.pos.addScaledVector(other.speed, -tCol)      
    }
}

