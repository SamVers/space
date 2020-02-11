import * as THREE from 'three/build/three.module.js'
import {singleLinkedList} from '../util/util.js'

export class octreeClass {

    // creates the octant-tree recursively
    constructor(p1, p2, level) {

        // the octant is a box defined by two diagonal points
        this.space = new THREE.Box3(p1, p2)

        // octant is an array that will contain the sub-octants if any
        this.octant = null //[null, null, null, null, null, null, null, null]

        // the linked list that will contain all objects in this octant that do not fit in a sub-octant
        this.objectList = new singleLinkedList

        // if not at the lowest level (note 1 = 8, 2 = 64, 3 = 512 suboctants)
        if (level > 0) {

            // it contains 8 new octants
            this.octant = [null, null, null, null, null, null, null, null]

            // dp is a vector with the size of the octant at this level
            let dp = new THREE.Vector3( (p2.x-p1.x)/2,(p2.y-p1.y)/2,(p2.z-p1.z)/2 ) 
            
            // get a new vector
            let p = new THREE.Vector3()

            // make a total of eight new octants
            let n = 0
            for (let i=0; i<2; i++)
                for (let j=0; j<2; j++)
                    for (let k=0; k<2; k++) {

                        // p defines one point of the octant..
                        p.set( p1.x + i*dp.x, p1.y + j*dp.y , p1.z + k*dp.z)

                        // ..and p+dp is the diagonally opposed point 
                        this.octant[n++] = new octreeClass(p, p.clone().add( dp ), level-1)
                    }
        }
    }

    // when the size of the container changes, the sizes of the octants have to be adapted as well
    sizeChange(p1,p2) {

        // change the size of the box
        this.space.min = {...p1}
        this.space.max = {...p2}

        // if there are there still sub-octants..
        if (this.octant) {

            // dp is a vector with the size of the octant at this level
            let dp = new THREE.Vector3( (p2.x-p1.x)/2,(p2.y-p1.y)/2,(p2.z-p1.z)/2 ) 
            
            // a copy of p1
            let p = new THREE.Vector3()

            // adjust the size for for all octants
            let n = 0
            for (let i=0; i<2; i++)
                for (let j=0; j<2; j++)
                    for (let k=0; k<2; k++) {

                        // p defines one point of the octant..
                        p.set( p1.x + i*dp.x, p1.y + j*dp.y , p1.z + k*dp.z)

                        // ..and p+dp is the diagonally opposed point 
                        this.octant[n++].sizeChange(p, p.clone().add( dp ))
                    }
        }
    }

    // all the object lists are reset
    reset() {
        // reset the object list
        this.objectList.reset()

        // check the sub-octants
        if (this.octant) 
            for (let i=0; i < this.octant.length; i++) 
                if ( this.octant[i]) this.octant[i].reset()
    }

    // object is checked against the other objects already in the octree
    // secTime is the duration of the interval in seconds 
    // collisions is an array that will contain all collision events in the interval
    collisionCheck( object, secTime, collisions) {

        // take the first object on the objects in this octant
        let other = this.objectList.first
        let when = 0.0

        // check against collisions with objects stored at this level
        while (other != null ) {

            // calculate if the objects intersect and if so at what time in the interval - inifinity means no collision
            if ( ( when = object.intersects(other, secTime)) != Infinity ) collisions.add(object,other,when)
            
            // take the next object
            other = other.next
        }

        // check the sub-octants
        if (this.octant) {

            // check for each octant..
            for (let i=0; i<this.octant.length; i++) {

                // ..if the object bbox has an intersection with the octant space 
                if (object.bBox.intersectsBox(this.octant[i].space)) 

                    // ..then check for collisions in the octant - 'true' is returned if the object fits completely in a (sub) octant
                    if ( this.octant[i].collisionCheck( object, secTime, collisions ) ) return true
            }
        }
        // The object did not fit in the sub-octants - check if the object fits completely in this octant 
        if ( ! this.space.containsBox( object.bBox )) return false
        
        // ...it fits - put it on the object list for this octant
        this.objectList.push(object)

        // we're done
        return true            
    }

    // count the nr of objects in the octree - just for debugging
    count() {
        let n = this.objectList.length()

        if (this.octant) 
            for (let i=0; i<this.octant.length; i++)  
                n += this.octant[i].count()
        return n
    }

    // count the total nr of sub-octants - just for debugging
    size() {
        let n = 1

        if (this.octant) {
            for (let i=0; i<this.octant.length; i++)  
                n += this.octant[i].size()
        }
        return n
    }
}// end of octree class

