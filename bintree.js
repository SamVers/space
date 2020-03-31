// import * as THREE from 'three/build/three.min.js'
import * as THREE from 'three'
import {singleLinkedList} from '../util/util.js'

export class bintreeClass {

    // creates the bin-tree recursively
    constructor(p1, p2, level) {

        // the space is a box defined by two diagonal points
        this.space = new THREE.Box3(p1, p2)

        // octant is an array that will contain the sub-octants if any
        this.bintree = null 

        // the linked list that will contain all objects in this octant that do not fit in a sub-octant
        this.objectList = new singleLinkedList

        // if not at the lowest level (note 1 = 8, 2 = 64, 3 = 512 suboctants)
        if (level > 0) {

            // it contains 2 new octants
            this.bintree = [null, null]

            // dp is a vector with the size of the octant at this level
            let delta = new THREE.Vector3( p2.x-p1.x,p2.y-p1.y,p2.z-p1.z ) 
            
            let q1 = p1.clone()
            let q2 = p2.clone()
            q2.z = q1.z + delta.z/2
            this.bintree[0] = new bintreeClass(q1,q2, level-1 )
            q1.z += delta.z/2
            q2.z += delta.z/2
            this.bintree[1] = new bintreeClass(q1,q2, level-1 )
        }
    }

    // when the size of the container changes, the sizes of the octants have to be adapted as well
    sizeChange(p1,p2) {

        // change the size of the box
        this.space.min = {...p1}
        this.space.max = {...p2}

        // if there are there still sub-octants..
        if (this.bintree) {

            let q1 = p1.clone()
            let q2 = p2.clone()
            q2.z = q1.z + delta.z/2
            this.bintree[0].sizeChange(q1,q2)
            q1.z += delta.z/2
            q2.z += delta.z/2
            this.bintree[1].sizeChange(q1,q2)
        }
    }

    // all the object lists are reset
    reset() {
        // reset the object list
        this.objectList.reset()

        // check the sub-octants
        if (this.bintree) {

            this.bintree[0].reset()
            this.bintree[1].reset()
        }
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
        if (this.bintree) {

            // check for each octant..
            for (let i=0; i<this.bintree.length; i++) {

                // ..if the object bbox has an intersection with the octant space 
                if (object.bBox.intersectsBox(this.bintree[i].space)) 

                    // ..then check for collisions in the octant - 'true' is returned if the object fits completely in a (sub) octant
                    if ( this.bintree[i].collisionCheck( object, secTime, collisions ) ) return true
            }
        }
        // The object did not fit in the subtrees - check if the object fits completely in this bin
        if ( ! this.space.containsBox( object.bBox )) return false
        
        // ...it fits - put it on the object list for this tree
        this.objectList.push(object)

        // we're done
        return true            
    }

    // count the nr of objects in the octree - just for debugging
    count() {
        let n = this.objectList.length()

        if (this.bintree) 
            for (let i=0; i<this.bintree.length; i++)  
                n += this.octant[i].count()
        return n
    }

    // count the total nr of sub-octants - just for debugging
    size() {
        let n = 1

        if (this.bintree) {
            for (let i=0; i<this.bintree.length; i++)  
                n += this.octant[i].size()
        }
        return n
    }
}// end of bintree class

