import * as THREE from 'three'
import {singleLinkedList} from '../util/util.js'

export class longColliderClass {

    constructor(p1, p2, level) {

        let q1=p1.clone(), q2=p1.clone
        let n = 2**level
        this.zStep = p2.z - p1.z / n

        // Make a simple list of boxes - each box has a list of the particles contained in the box
        for (let i = 0; i < n; i++) {       
            q1.z = q2.z
            q2.z += zStep
            this.slices[i] = { box:new THREE.Box3(q1,q2), particleList: new singleLinkedList }
        }
    }

    // when the size of the container changes, the sizes of the octants have to be adapted as well
    sizeChange(p1,p2) {

        let q1=p1.clone(), q2=p1.clone
        let n = this.slices.length
        this.zStep = p2.z-p1.z/n
        for (let i = 0; i < n; i++) {       
            q1.z = q2.z
            q2.z += zStep
            this.slices[i].box.min = {...q1}
            this.slices[i].box.max = {...q2}
        }
    }

    // all the particle lists are reset
    reset() {
        // reset the particle lis
        this.slices.forEach( (slice) => slice.particleList.reset())
    }

    // particle is checked against the other particles already in the octree
    // secTime is the duration of the interval in seconds 
    // collisions is an array that will contain all collision events in the interval
    collisionCheck( particle, secTime, collisions) {

        let n = this.slices.length
        let other = null
        zMin = particle.bBox.min.z
        zMax = particle.bBox.max.z

        // check all slices
        for (let i = 0; i < n; i++) {     
            slice = this.slices[i]

            // check if particle passes through this slice
            if ((zMin < slice.box.max.z ) && (zMax > slice.box.min.z)) {
                other = slice.particleList.first
                while (other) {
                    let when = 0.0
                    if ( ( when = particle.intersects(other, secTime)) != Infinity ) {
                        collisions.add(particle,other,when)
                        slice.particleList.attach(particle)
                        return true;
                    }
                    other = other.next   
                }      
            }
            // there are no more slices where the particle passes through
            else return false
        }
        return false
    }
}// end of class

