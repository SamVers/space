//import * as THREE from 'three/build/three.min.js'
import * as THREE from 'three'
import * as SI from '../../lib/util/SI-units.js'

// The box container
export class cuboidContainerClass {

constructor(scene, spec) {

    // spec = {vector: (x,y,z) and unit: "m"}

    // size of the bounding box - everything in meter 
    let factor = SI.factor("length", spec.unit, "m") 
    let s = {x: factor*spec.vector.x, y:factor*spec.vector.y, z:factor*spec.vector.z}

    // positioned around 0
    let p = {x:0, y:0, z:0}

    // geometry, material, mesh
    this.geometry = this.material = this.mesh = null

    // make the bounding box
    this.bBox = new THREE.Box3( new THREE.Vector3(p.x-s.x/2, p.y-s.y/2, p.z-s.z/2), 
                                new THREE.Vector3(p.x+s.x/2, p.y+s.y/2, p.z+s.z/2))

    // the material settings
    this.materialSettings = { color: "#aabbdd",wireframe: true,transparent: true,opacity: 0.5}

    // the force exerted on the walls of the container 
    this.force = {  x1:0.0, x2:0.0,
                    y1:0.0, y2:0.0,
                    z1:0.0, z2:0.0 }

    // the nr of collisions on each wall
    this.collisions = { nX1:0, nX2:0,
                        nY1:0, nY2:0,
                        nZ1:0, nZ2:0}

    // we only need the rest if we have a scene to draw into
    if (!scene) return

    // make the material
    this.makeMaterial()

    // make the mesh of the box
    this.makeGeometry(s)

    // make the mesh
    this.addMesh(scene)
}

destructor(scene) {
    if (this.mesh) scene.remove(this.mesh), this.mesh=null
    if (this.geometry) this.geometry.dispose(), this.geometry = null
    if (this.material) this.material.dispose(), this.material = null
}

changeSize(scene, spec) {

    // size of the bounding box - everything in meter 
    let factor = SI.factor("length", spec.unit, "m") 
    let s = {x: factor*spec.vector.x, y:factor*spec.vector.y, z:factor*spec.vector.z}

    // positioned around 0
    let p = {x:0, y:0, z:0}

    // make the bounding box
    this.bBox = new THREE.Box3( new THREE.Vector3(p.x-s.x/2, p.y-s.y/2, p.z-s.z/2), 
                                new THREE.Vector3(p.x+s.x/2, p.y+s.y/2, p.z+s.z/2))

    // we only need the rest if we have a scene to draw into
    if (!scene) return   

    // make the mesh of the box
    this.makeGeometry(s)

    // make the mesh
    this.addMesh(scene)
}

makeGeometry(s) {

    // we make some subdivisions on the cuboid
    let min = s.x < s.y ? ( s.x < s.z ? s.x : s.z) : ( s.y < s.z ? s.y : s.z)

    // calculate the nr of subdivisions 
    let [nx,ny,nz] = [4*Math.round(s.x/min), 4*Math.round(s.y/min), 4*Math.round(s.z/min)]

    // max nr of divisions
    const maxDiv = 16
    nx = nx > maxDiv ? maxDiv : nx,  ny= ny > maxDiv ? maxDiv : ny,  nz = nz > maxDiv ? maxDiv : nz

    // get rid of old geometry if any and create the new one
    if (this.geometry) this.geometry.dispose()
    this.geometry = new THREE.CubeGeometry(s.x,s.y,s.z, nx, ny, nz)
}

makeMaterial() {
    if (this.material) this.material.dispose()
    this.material = new THREE.MeshLambertMaterial(this.materialSettings)
}

addMesh(scene) {
    if (this.mesh) scene.remove(this.mesh)
    this.mesh = new THREE.Mesh(this.geometry, this.material)
    scene.add(this.mesh)
}

magnitude() {
    return this.bBox.max.clone().sub(this.bBox.min).length()
}

// we do not care about the collision time - reflection is simple to calculate
// we set collision time at 100 - after sorting collisions with the container will be at the end
collisionCheck(object, interval, collisions) {

    // notation
    let p = object.pos, box = this.bBox, r = object.radius
    
    // check if the object is beyond the box - in the x, y or z direction
    if ((p.x < box.min.x + r) || (p.x > box.max.x - r)||
        (p.y < box.min.y + r) || (p.y > box.max.y - r)|| 
        (p.z < box.min.z + r) || (p.z > box.max.z - r)  ) {

        // if it is push it on the collision-array 
        collisions.add(this, object, 100.0)
    }
}

// Calculates the new position of the object after collision with the cuboid
// Note that the position of the object is already the new position but without the collision
// Here we simply calculate the reflected position and the reflected speed
collide( object, tCol, interval ) {

    let p = object.pos,  v = object.speed,  r = object.radius, m= object.mass
    let box = this.bBox
    
    if (p.x < box.min.x + r) {
        p.x = -p.x + 2*(box.min.x + r)
        v.x = -v.x
        this.force.x1 += 2*m*v.x
        this.collisions.nX1++
    }
    else if ( p.x > box.max.x - r ) {
        p.x = -p.x + 2*(box.max.x - r)
        v.x = -v.x
        this.force.x2 -= 2*m*v.x
        this.collisions.nX2++
    }
    
    if (p.y < box.min.y + r) {
        p.y = -p.y + 2*(box.min.y + r)
        v.y = -v.y
        this.force.y1 += 2*m*v.y
        this.collisions.nY1++
    }
    else if (p.y > box.max.y - r) {
        p.y = -p.y + 2*(box.max.y - r)
        v.y = -v.y
        this.force.y2 -= 2*m*v.y
        this.collisions.nY2++
    }

    if ( p.z < box.min.z + r) {
        p.z = -p.z + 2*(box.min.z + r)
        v.z = -v.z
        this.force.z1 += 2*m*v.z
        this.collisions.nZ1++
    }
    else if ( p.z > box.max.z - r ) {
        p.z = -p.z + 2*(box.max.z - r)
        v.z = -v.z
        this.force.z2 -= 2*m*v.z
        this.collisions.nZ2++
    }
}

getForceAndArea() {

    // notation
    let f = this.force
    let s = this.bBox.max.clone().sub(this.bBox.min)

    // calculate
    let force = f.x1 + f.x2 + f.y1 + f.y2 + f.z1 + f.z2
    let area = 2*(s.x*s.y + s.x*s.z + s.y*s.z)
    
    // reset the force
    f.x1 = f.x2 = f.y1 = f.y2 = f.z1 = f.z2 = 0.0

    // done
    return [force, area]
}

getCollisions() {

    let c = this.collisions
    let count = c.nX1 + c.nX2 + c.nY1 + c.nY2 + c.nZ1 + c.nZ2
    c.nX1 = c.nX2 = c.nY1 = c.nY2 = c.nZ1 = c.nZ2 = 0
    return count
}

reset() {
    // reset the nr of collisions
    c.nX1 = c.nX2 = c.nY1 = c.nY2 = c.nZ1 = c.nZ2 = 0

    // reset the force
    f.x1 = f.x2 = f.y1 = f.y2 = f.z1 = f.z2 = 0.0
}

} // end of container class