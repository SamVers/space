import * as THREE from 'three/build/three.min.js'
import * as SI  from '../util/util.js'
import {objectClass} from './object.js'

// colors used for bounding boxes
const boxFreeClr = new THREE.Color("#00FF00")
const boxHitClr =  new THREE.Color("#FF0000")

const cmTom = 0.01


// This function generates random numbers along a Normal or Log-normal distribution 
// using the Marsaglia polar method. the mean is 0 and the deviation is 1
// the function generates two values at each invocation - that is why we have the spare
let spareRandom = null
function normalRandom() {

    let val, u, v, s, mul;

    if(spareRandom !== null) {
        val = spareRandom;
        spareRandom = null;
    }
    else {
        do {
            u = Math.random()*2-1;
            v = Math.random()*2-1;
            s = u*u+v*v;
        } while(s === 0 || s >= 1);
        mul = Math.sqrt(-2 * Math.log(s) / s);
        val = u * mul;
        spareRandom = v * mul;
    }
    return val;
}

export class objectGroupClass {

    constructor(init) {

        // save the name
        this.name = String(init.name)

        // we also save the color
        this.color = init.color

        // the objects positions etc
        this.objects = []

        // the helper group
        this.boxHelperGroup = null

        // we keep the mesh, geometry and material
        this.mesh = this.geometry = this.material = null

        // the color for the objects of this group
        let color = new THREE.Color(this.color)

        // the material
        this.material = new THREE.MeshLambertMaterial({color: color})

        // convert to meter
        let radius = init.objectRadius ?
                     +init.objectRadius.value * SI.factor("length",init.objectRadius.unit, "m")
                     : 0.0

        // object geometry
        this.geometry = new THREE.SphereGeometry(radius , 12, 12)
    }

    newObjectCount(scene, count, mass, radius) {

        // remove the box helpers if any
        if (this.boxHelperGroup != null) {

            scene.remove( this.boxHelperGroup )
            this.boxHelperGroup = null
        }

        // get a new array for the objects
        this._setObjectArrays(count, mass, radius)

        // if we still have a mesh - release it
        if (this.mesh) {
            scene.remove(this.mesh)
            this.mesh = null
        }
    }

    _setObjectArrays(count, mass, radius) {

        // if the array is too big -
        if (this.objects.length >= count) {

            //...trim to count
            this.objects.length = count
        }
        else {
            // if too small - free
            this.objects.length = 0
        
            // allocate a new array
            this.objects = new Array(count)

            // initialize
            for (let i=0; i<count; i++) this.objects[i] = new objectClass(mass, radius)
        }
    }

    checkMassAndRadius(container, objectCount, objectMass, objectRadius) {

        // the number of objects
        let count = +objectCount.current

		// convert the mass of the object to g
		let mass = +objectMass.value * SI.factor("mass",objectMass.unit,"g")

		// convert the radius to meter
        let radius = +objectRadius.value * SI.factor("length",objectRadius.unit, "m")
        
        return [count, mass, radius]
    }

    // set the balls in the cube
    initPosition(container, ballBox, location) {

        // check
        if (this.objects.length < 1) return

        // the box where the balls have to fit
        let cuboid = ballBox.vector

        // the displacement of the box
        let slide = (container.bBox.max.x - container.bBox.min.x - cuboid.x) * location.current/location.max

        // the center of the box is
        let center = new THREE.Vector3()
        center.x = container.bBox.min.x + cuboid.x/2 + slide
        center.y = (container.bBox.max.y + container.bBox.min.y)/2
        center.z = (container.bBox.max.z + container.bBox.min.z)/2

        this.randomInCuboid(cuboid, center)
    }

    // initialise the bounding boxes of the objects
    initBBox() {
        let object = null
        for (let i=0; i< this.objects.length; i++) {
            object = this.objects[i]
            object.bBox.min.copy(object.pos).subScalar(object.radius)
            object.bBox.max.copy(object.pos).addScalar(object.radius)
        }
    }

    randomInCuboid(cuboid, center){

        let object =null
        let r = this.objects[0].radius

        let x = cuboid.x - 2*r, dx = - cuboid.x/2 + center.x
        let y = cuboid.y - 2*r, dy = - cuboid.y/2 + center.y
        let z = cuboid.z - 2*r, dz = - cuboid.z/2 + center.z

        // everything is centered around 0 - so we shift half
        for (let i=0; i< this.objects.length; i++) {

            object = this.objects[i]
            object.pos.set( x * Math.random() + dx, 
                            y * Math.random() + dy, 
                            z * Math.random() + dz)
        }
    }

    slideObjects(container, ballBox, location) {

        // check
        if (this.objects.length < 1) return

        let pos = null
        let r = this.objects[0].radius
        let max = container.bBox.max, min = container.bBox.min

        // the displacement of the box
        let slide = (max.x - min.x - ballBox.vector.x) * location.delta/location.max
        for (let i=0; i< this.objects.length; i++) {
            pos = this.objects[i].pos
            pos.x += slide
            if (pos.x < min.x + r ) pos.x = min.x + r
            if (pos.x > max.x - r ) pos.x = max.x - r
        }
    }

    orderedInCuboid(cuboid){

        // the nr of vectors we need
        let nVectors = this.objects.length
        let r = this.objects[0].radius

        // the raster we will use - note that nx*ny*nz > nVectors !
        let nx = Math.ceil( (nVectors*(cuboid.x/cuboid.y)*(cuboid.x/cuboid.z))**(1/3))
        let ny = Math.ceil( (nVectors*(cuboid.y/cuboid.x)*(cuboid.y/cuboid.z))**(1/3))
        let nz = Math.ceil( (nVectors*(cuboid.z/cuboid.x)*(cuboid.z/cuboid.y))**(1/3))
    
        // the spacing between objects
        // let dx=3*r,dy=3*r,dz=3*r
        // if (dx*nx > cuboid.x) dx = cuboid.x/nx
        // if (dy*ny > cuboid.y) dy = cuboid.y/ny
        // if (dz*nz > cuboid.z) dz = cuboid.z/nz

        let dx = cuboid.x/nx
        let dy = cuboid.y/ny
        let dz = cuboid.z/nz
    
        // the indices
        let ix = 0, iy = 0, iz = 0, count = 0

        // we center around 0
        let offset = {x: -dx*(nx-1)/2, y: -dy*(ny-1)/2, z: -dz*(nz-1)/2}
 
        // caculate the positions of all objects
        for (iz=0; iz<nz; iz++) {
            for (iy=0; iy<ny; iy++) {
                for (ix=0; ix<nx; ix++) {              
                    this.objects[count].pos.set(ix*dx + offset.x, 
                                                iy*dy + offset.y, 
                                                iz*dz + offset.z)
                    if (++count == nVectors) return
                }
            }
        }
    }

    randomInSphere(radius) {

        let phi=0.0, theta=0.0, R=0.0
        for (let i=0; i< this.objects.length; i++) {
            phi   = Math.random()*Math.PI*2
            theta = Math.random()*Math.PI
            R     = Math.random()*radius
            this.objects[i].pos.set(R*Math.sin(theta)*Math.cos(phi), 
                                    R*Math.sin(theta)*Math.sin(phi),
                                    R*Math.cos(theta))
        }
    }

    orderedInSphere(radius) {

        let R=0.0
        let dphi=2*Math.PI
        let dtheta=Math.PI
        let r = 2*this.objects[0].radius
        let count = 0

        for (let shell=1; shell < 10; shell++) {

            R = shell * r
            dphi = 2*Math.PI/(shell*6)
            dtheta = Math.PI/(shell*6)
            for (let i = 0; i <= shell*6; i++) {
                for (let j=0; j<= shell*6;j++) {

                    this.objects[count].pos.set(
                        R*Math.sin(dtheta*j)*Math.cos(dphi*i), 
                        R*Math.sin(dtheta*j)*Math.sin(dphi*i),
                        R*Math.cos(dtheta*j))

                    if (++count == this.objects.length) return
                }
            }
        }
    }

    placeObjectsInScene(scene) {

        // the matrix to position the objects
        let matrix = new THREE.Matrix4()

        // check if we have to place a new mesh
        if (this.mesh == null) {

            // create an instanced mesh with the geometry - it will be instanced nObjects time
            this.mesh = new THREE.InstancedMesh(this.geometry, this.material, this.objects.length)

            // add to the scene
            scene.add(this.mesh)
        }

        // set the update bit
		this.mesh.instanceMatrix.needsUpdate = true

        // place every instance...
        for (let i=0; i<this.objects.length; i++) {

            // make the matrix to position the object
            matrix.setPosition(this.objects[i].pos)

            // position the object
            this.mesh.setMatrixAt(i,matrix)
        }
    }

    // {   text:"fixed value",         value:{type:"fixed",fixedSpeed}},
    // {   text:"normal distributed",  value:{type:"rnd-normal", mu, sigma}},
    // {   text:"uniform distributed", value:{type:"rnd-uniform",a, b}},

    initSpeed(direction, amplitude) {

        let object = null
        let A=0, fixedSpeed=0, sigma=0, mu=0, a=0, b=0
        let phi = 0.0, theta = 0.0

        // to speed up the calculations
        switch(amplitude.type) {

            case "fixed":  
                fixedSpeed = +amplitude.fixedSpeed.value
                break

            case "rnd-normal": 
                sigma = +amplitude.sigma.value 
                mu = +amplitude.mu.value
                break

            case "rnd-uniform": 
                a = +amplitude.a.value
                b = +amplitude.b.value
                break
        }
        
        // calculate for each object
        for(let i=0;i<this.objects.length;i++) {

            // notation
            object = this.objects[i]

            // first get the amplitue
            switch(amplitude.type) {
                case "fixed": A = fixedSpeed
                    break

                case "rnd-normal": A = normalRandom()*sigma + mu
                    break

                case "rnd-uniform": A = (Math.random() * (b-a)) + a
                    break
            }

            // then the vector
            switch(direction.type) {
        
                case "random":
                    // get the random vectors using polar coordinates
                    phi = Math.random()*Math.PI*2
                    theta = Math.random()*Math.PI
                    object.speed.set(   A*Math.sin(theta)*Math.cos(phi), 
                                        A*Math.sin(theta)*Math.sin(phi),
                                        A*Math.cos(theta))
                    break
        
                case "radial":
                    object.speed.copy(object.pos).sub(direction.center).setLength(A)
                    break
        
                case "parallel":
                    object.speed.copy(direction.parallel).setLength(A)
                    break
        
                default: 
                    console.log("[UNKNOWN] unknown speed vector type in function speedVectors", direction.type)
                    return
            }
        }
        //console.log("objects",this.objects )
    }

    // put the helpers in a THREE.Group !
    showAABB(scene) {

        let object = null

        // check if we have to create the box helper group
        if (this.boxHelperGroup == null ) {

            this.boxHelperGroup = new THREE.Group()

            for (let i=0; i < this.objects.length; i++) {

                // notation
                object = this.objects[i]

                // create a boxhelper object
                object.boxHelper = new THREE.Box3Helper( object.bBox, boxFreeClr )

                // add the helpers to the box helper group
                this.boxHelperGroup.add(object.boxHelper)
            }
        }

        // add the helper group to the scene
        scene.add(this.boxHelperGroup)
    }

    hideAABB(scene) {

        if (this.boxHelperGroup) scene.remove( this.boxHelperGroup )
    }

    changeMass(value, unit) {

        // convert the mass of the object to g
        let mass = (+value) * SI.factor("mass",unit,"g")

        // change all masses
        for (let i=0; i < this.objects.length; i++) 
            this.objects[i].mass = mass
    }

    changeRadius(value, unit) {

        // recalculate the radius
        let radius = (+value) * SI.factor("length",unit,"m")

        // we have to change the geometry
        this.geometry.dispose()

        // get a new geometry
        this.geometry = new THREE.SphereGeometry(radius , 12, 12) 
        
        // and put it in the mesh
        this.mesh.geometry = this.geometry

        // change all objects
        for (let i=0; i < this.objects.length; i++) 
            this.objects[i].radius = radius
    }
}