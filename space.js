define(['exports', 'three/build/three.module.js', 'three/examples/jsm/controls/OrbitControls.js'], function (exports, THREE, OrbitControls_js) { 'use strict';

    // this list contains all the objects that will collide and the time of collision
    // it also keeps track of the statistics of the collisions
    class collisionClass {

        constructor() {

            this.list = [];
            this.count = 0;
            this.period = 0.0;
        }

        resetCount() {
            this.count = 0;
            this.period = 0.0;
        }

        resetList() {
            this.list = [];
        }

        add( a, b, t) {

            this.list.push({a,b,t});
        }

        collide(interval) { 

            let list = this.list;
            
            // check
            if (list.length < 1) return

            // sort the list 
            list.sort( (r1, r2) => r1.t - r2.t );

            // do all collisions
            let len = list.length;
            let col = null;
            for (let i=0; i < len ; i++) {

                // get the collision record
                col = list[i];

                // parameters are the other object, the collision time (relative to the start if the interval) and the interval duration
                col.a.collide(col.b, col.t, interval);
            }

            this.period += interval;
            this.count += len;
        }

        get() {
            return {count: this.count, period: this.period}
        }
    }

    // the units that are 
    const allSeries = [
        {type:"unknown", units: [ {str:"?", exp:0}] },
        {type:"length",  units: [ {str:"km", exp:3}, {str:"m", exp:0}, {str:"cm", exp:-2}, {str:"mm", exp:-3}, {str:"µm", exp:-6}, {str:"nm", exp:-9}, {str:"pm", exp:-12}] },
        {type:"mass",    units: [ {str:"Gt", exp:12}, {str:"Mt",exp:9}, {str:"kt",exp:6}, {str:"t",exp:3}, {str:"kg",exp:0}, {str:"g",exp:-3}, {str:"µg",exp:-6}, {str:"ng",exp:-9}, {str:"pg",exp:-12},{str:"Da",exp:-27},]},
        {type:"speed",   units: [ {str:"km/s", exp:3},{str:"km/h", exp:0}, {str:"m/s", exp:0}, {str:"cm/s", exp:-2}, {str:"mm/s", exp:-3}, {str:"µm/s", exp:-6}, {str:"nm/s", exp:-9}] },
    ];

    // conversion factor from one unit to the other
    function exponent(type, unitFrom, unitTo) {

        // find the units
        let series = allSeries.find( (s) => s.type === type);

        // check 
        if (series === undefined) return null

        // find the two units in the series
        let F = series.units.find( (u) => u.str === unitFrom );
        let T = series.units.find( (u) => u.str === unitTo );

        // check
        if (F === undefined || T === undefined) return null

        // subtract the exponents
        return (F.exp - T.exp)
    }

    function factor(type, unitFrom, unitTo) {
        
        let exp = exponent(type, unitFrom, unitTo);

        if (exp == null) return null
        else return 10**exp
    }

    // The box container
    class cuboidContainerClass {

    constructor(scene, spec) {

        // spec = {vector: (x,y,z) and unit: "m"}

        // size of the bounding box - everything in meter 
        let factor$1 = factor("length", spec.unit, "m"); 
        let s = {x: factor$1*spec.vector.x, y:factor$1*spec.vector.y, z:factor$1*spec.vector.z};

        // positioned around 0
        let p = {x:0, y:0, z:0};

        // geometry, material, mesh
        this.geometry = this.material = this.mesh = null;

        // make the bounding box
        this.bBox = new THREE.Box3( new THREE.Vector3(p.x-s.x/2, p.y-s.y/2, p.z-s.z/2), 
                                    new THREE.Vector3(p.x+s.x/2, p.y+s.y/2, p.z+s.z/2));

        // the material settings
        this.materialSettings = { color: "#aabbdd",wireframe: true,transparent: true,opacity: 0.5};

        // the force exerted on the walls of the container 
        this.force = {  x1:0.0, x2:0.0,
                        y1:0.0, y2:0.0,
                        z1:0.0, z2:0.0 };

        // the nr of collisions on each wall
        this.collisions = { nX1:0, nX2:0,
                            nY1:0, nY2:0,
                            nZ1:0, nZ2:0};

        // we only need the rest if we have a scene to draw into
        if (!scene) return

        // make the material
        this.makeMaterial();

        // make the mesh of the box
        this.makeGeometry(s);

        // make the mesh
        this.addMesh(scene);
    }

    changeSize(scene, spec) {

        // size of the bounding box - everything in meter 
        let factor$1 = factor("length", spec.unit, "m"); 
        let s = {x: factor$1*spec.vector.x, y:factor$1*spec.vector.y, z:factor$1*spec.vector.z};

        // positioned around 0
        let p = {x:0, y:0, z:0};

        // make the bounding box
        this.bBox = new THREE.Box3( new THREE.Vector3(p.x-s.x/2, p.y-s.y/2, p.z-s.z/2), 
                                    new THREE.Vector3(p.x+s.x/2, p.y+s.y/2, p.z+s.z/2));

        // we only need the rest if we have a scene to draw into
        if (!scene) return   

        // make the mesh of the box
        this.makeGeometry(s);

        // make the mesh
        this.addMesh(scene);
    }

    makeGeometry(s) {

        // we make some subdivisions on the cuboid
        let min = s.x < s.y ? ( s.x < s.z ? s.x : s.z) : ( s.y < s.z ? s.y : s.z);

        // calculate the nr of subdivisions and limit to 32
        let [nx,ny,nz] = [4*Math.round(s.x/min), 4*Math.round(s.y/min), 4*Math.round(s.z/min)];

        // max nr of divisions
        const maxDiv = 16;
        nx = nx > maxDiv ? maxDiv : nx,  ny= ny > maxDiv ? maxDiv : ny,  nz = nz > maxDiv ? maxDiv : nz;

        // get rid of old geometry if any and create the new one
        if (this.geometry) this.geometry.dispose();
        this.geometry = new THREE.CubeGeometry(s.x,s.y,s.z, nx, ny, nz);
    }

    makeMaterial() {
        if (this.material) this.material.dispose();
        this.material = new THREE.MeshLambertMaterial(this.materialSettings);
    }

    addMesh(scene) {
        if (this.mesh) scene.remove(this.mesh);
        this.mesh = new THREE.Mesh(this.geometry, this.material);
        scene.add(this.mesh);
    }

    magnitude() {
        return this.bBox.max.clone().sub(this.bBox.min).length()
    }

    // we do not care about the collision time - reflection is simple to calculate
    // we set collision time at 100 - after sorting collisions with the container will be at the end
    collisionCheck(object, interval, collisions) {

        // notation
        let p = object.pos, box = this.bBox, r = object.radius;
        
        // check if the object is beyond the box - in the x, y or z direction
        if ((p.x < box.min.x + r) || (p.x > box.max.x - r)||
            (p.y < box.min.y + r) || (p.y > box.max.y - r)|| 
            (p.z < box.min.z + r) || (p.z > box.max.z - r)  ) {

            // if it is push it on the collision-array 
            collisions.add(this, object, 100.0);
        }
    }

    // Calculates the new position of the object after collision with the cuboid
    // Note that the position of the object is already the new position but without the collision
    // Here we simply calculate the reflected position and the reflected speed
    collide( object, tCol, interval ) {

        let p = object.pos,  v = object.speed,  r = object.radius, m= object.mass;
        let box = this.bBox;
        
        if (p.x < box.min.x + r) {
            p.x = -p.x + 2*(box.min.x + r);
            v.x = -v.x;
            this.force.x1 += 2*m*v.x;
            this.collisions.nX1++;
        }
        else if ( p.x > box.max.x - r ) {
            p.x = -p.x + 2*(box.max.x - r);
            v.x = -v.x;
            this.force.x2 -= 2*m*v.x;
            this.collisions.nX2++;
        }
        
        if (p.y < box.min.y + r) {
            p.y = -p.y + 2*(box.min.y + r);
            v.y = -v.y;
            this.force.y1 += 2*m*v.y;
            this.collisions.nY1++;
        }
        else if (p.y > box.max.y - r) {
            p.y = -p.y + 2*(box.max.y - r);
            v.y = -v.y;
            this.force.y2 -= 2*m*v.y;
            this.collisions.nY2++;
        }

        if ( p.z < box.min.z + r) {
            p.z = -p.z + 2*(box.min.z + r);
            v.z = -v.z;
            this.force.z1 += 2*m*v.z;
            this.collisions.nZ1++;
        }
        else if ( p.z > box.max.z - r ) {
            p.z = -p.z + 2*(box.max.z - r);
            v.z = -v.z;
            this.force.z2 -= 2*m*v.z;
            this.collisions.nZ2++;
        }
    }

    getForceAndArea() {

        // notation
        let f = this.force;
        let s = this.bBox.max.clone().sub(this.bBox.min);

        // calculate
        let force = f.x1 + f.x2 + f.y1 + f.y2 + f.z1 + f.z2;
        let area = 2*(s.x*s.y + s.x*s.z + s.y*s.z);
        
        // reset the force
        f.x1 = f.x2 = f.y1 = f.y2 = f.z1 = f.z2 = 0.0;

        // done
        return [force, area]
    }

    getCollisions() {

        let c = this.collisions;
        let count = c.nX1 + c.nX2 + c.nY1 + c.nY2 + c.nZ1 + c.nZ2;
        c.nX1 = c.nX2 = c.nY1 = c.nY2 = c.nZ1 = c.nZ2 = 0;
        return count
    }

    } // end of container class

    // a doubly linked list of elements
    class singleLinkedList {

        constructor() {

            this.n = 0;
            this.first = null;
            this.last = null;
        }

        reset() {
            this.n = 0;
            this.first = null;
            this.last = null;        
        }

        length() { return this.n }

        // attach detach at the end of the list
        attach(p) {
            p.next = null;
            if (this.last == null)  {
                this.first = this.last = p;
            }
            else {
                this.last.next = p;
                this.last = p; 
            }
            this.n++;
        }

        detach() {
            if (this.last==null) return null
            let p = this.last;
            if (this.last==this.first)
                this.first = this.last = null;
            else {
                this.last = this.first;
                while (this.last.next) this.last = this.last.next;
            }
            this.n--;
            return p
        }

        // push and pop at the front of the list
        push(p) {
            p.next = this.first;
            if (this.first == null) {
                this.first = this.last = p;
            }
            else {
                this.first = p;
            }
            this.n++;
        }

        pop() {
            if (this.first == null) return null
            let p = this.first;
            this.first = p.next;
            p.next = null;
            if (this.first == null) this.last = null;
            this.n--;
            return p
        }
    }

    // colors used for bounding boxes
    const boxFreeClr = new THREE.Color("#00FF00");
    const boxHitClr =  new THREE.Color("#FF0000");

    class objectClass {

        constructor(mass, radius) {

            this.pos = new THREE.Vector3();
            this.speed = new THREE.Vector3();
            this.radius = radius;
            this.mass = mass;

            // the bounding box when moving from one position to another
            this.bBox = new THREE.Box3;

            // for debugging we sometimes need to see the bounding boxes
            this.boxHelper = null;

            // objects can be added to a linked list (eg in the collisionlist)
            this.next = null;
        }

        // calculate the size of the bounding box based on the new position q and the size of tha ball
        setBBox(q) {

            // notation
            let p = this.pos, r = this.radius, min = this.bBox.min, max = this.bBox.max;

            // we subtract/add r to make sure the ball lies completely in the box
            min.x = p.x < q.x ? p.x-r : q.x-r;
            min.y = p.y < q.y ? p.y-r : q.y-r;
            min.z = p.z < q.z ? p.z-r : q.z-r;

            max.x = p.x > q.x ? p.x+r : q.x+r;
            max.y = p.y > q.y ? p.y+r : q.y+r;
            max.z = p.z > q.z ? p.z+r : q.z+r;

            // a box that has collided changes gradually to normal
            if (this.boxHelper) this.boxHelper.material.color.lerp(boxFreeClr,0.1);
        }

        newPosition(t) {

            // calculate the new position
            let q = this.pos.clone().addScaledVector(this.speed,t);

            // adjust the bounding box
            this.setBBox(q);

            // copy the new position
            this.pos.copy(q);
        }

        // a function to test for intersection of trajectories
        intersects( other, interval ) {

            if (this.bBox.intersectsBox(other.bBox)) {

                // we have a candidate - calculate relative position and speed
                // p = p2-p1, v= v2-v1 - note that p is already the new position t sec ago
                let p = other.pos.clone().sub(this.pos);
                let v = other.speed.clone().sub(this.speed);

                // calculate a,b and c of the quadratic distance equation
                let a = v.lengthSq();
                let b = 2*p.dot(v);
                let c = p.lengthSq() - (other.radius+this.radius)**2;

                // the determinant
                let det = b*b - 4*a*c;
                
                // if the determinant is positive there are intersections
                if ( det > 0) {

                    // show that the objects will collide
                    if (this.boxHelper) this.boxHelper.material.color.copy(boxHitClr);
                    if (other.boxHelper) other.boxHelper.material.color.copy(boxHitClr);

                    // calculate the roots - nota that t1 is the smallest root
                    let sqrt = Math.sqrt(det);
                    let t1 = (-b - sqrt)/(2*a);

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
            this.pos.addScaledVector(this.speed, tCol);
            other.pos.addScaledVector(other.speed, tCol);

            // formula to be found in wikipedia https://en.wikipedia.org/wiki/Elastic_collision#Two-dimensional
            let p = this.pos.clone().sub(other.pos);
            let v = this.speed.clone().sub(other.speed);

            // the factors for the speed calculation
            let f1 =  ( 2 * other.mass / (this.mass + other.mass) ) * ( v.dot(p) / p.lengthSq() );
            let f2 = this.mass/other.mass;

            // calculate the new speeds - keep the order !
            this.speed.sub( p.multiplyScalar(f1) );
            other.speed.add( p.multiplyScalar(f2) );

            // the two objects have still tCol seconds to travel with the new speed (-tCol is positive !)
            this.pos.addScaledVector(this.speed, -tCol);
            other.pos.addScaledVector(other.speed, -tCol);      
        }
    }

    // colors used for bounding boxes
    const boxFreeClr$1 = new THREE.Color("#00FF00");
    const boxHitClr$1 =  new THREE.Color("#FF0000");


    // This function generates random numbers along a Normal or Log-normal distribution 
    // using the Marsaglia polar method. the mean is 0 and the deviation is 1
    // the function generates two values at each invocation - that is why we have the spare
    let spareRandom = null;
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

    class objectGroupClass {

        constructor(init) {

            // save the name
            this.name = String(init.name);

            // we also save the color
            this.color = init.color;

            // the objects positions etc
            this.objects = [];

            // the helper group
            this.boxHelperGroup = null;

            // we keep the mesh, geometry and material
            this.mesh = this.geometry = this.material = null;

            // the color for the objects of this group
            let color = new THREE.Color(this.color);

            // the material
            this.material = new THREE.MeshLambertMaterial({color: color});

            // convert to meter
            let radius = init.objectRadius ?
                         +init.objectRadius.value * factor("length",init.objectRadius.unit, "m")
                         : 0.0;

            // object geometry
            this.geometry = new THREE.SphereGeometry(radius , 12, 12);
        }

        newObjectCount(scene, count, mass, radius) {

            // remove the box helpers if any
            if (this.boxHelperGroup != null) {

                scene.remove( this.boxHelperGroup );
                this.boxHelperGroup = null;
            }

            // get a new array for the objects
            this._setObjectArrays(count, mass, radius);

            // if we still have a mesh - release it
            if (this.mesh) {
                scene.remove(this.mesh);
                this.mesh = null;
            }
        }

        _setObjectArrays(count, mass, radius) {

            // if the array is too big -
            if (this.objects.length >= count) {

                //...trim to count
                this.objects.length = count;
            }
            else {
                // if too small - free
                this.objects.length = 0;
            
                // allocate a new array
                this.objects = new Array(count);

                // initialize
                for (let i=0; i<count; i++) this.objects[i] = new objectClass(mass, radius);
            }
        }

        checkMassAndRadius(container, objectCount, objectMass, objectRadius) {

            // the number of objects
            let count = +objectCount.current;

    		// convert the mass of the object to g
    		let mass = +objectMass.value * factor("mass",objectMass.unit,"g");

    		// convert the radius to meter
            let radius = +objectRadius.value * factor("length",objectRadius.unit, "m");
            
            return [count, mass, radius]
        }

        // set the balls in the cube
        initPosition(container, ballBox, location) {

            // check
            if (this.objects.length < 1) return

            // the box where the balls have to fit
            let cuboid = ballBox.vector;

            // the displacement of the box
            let slide = (container.bBox.max.x - container.bBox.min.x - cuboid.x) * location.current/location.max;

            // the center of the box is
            let center = new THREE.Vector3();
            center.x = container.bBox.min.x + cuboid.x/2 + slide;
            center.y = (container.bBox.max.y + container.bBox.min.y)/2;
            center.z = (container.bBox.max.z + container.bBox.min.z)/2;

            this.randomInCuboid(cuboid, center);
        }

        // initialise the bounding boxes of the objects
        initBBox() {
            let object = null;
            for (let i=0; i< this.objects.length; i++) {
                object = this.objects[i];
                object.bBox.min.copy(object.pos).subScalar(object.radius);
                object.bBox.max.copy(object.pos).addScalar(object.radius);
            }
        }

        randomInCuboid(cuboid, center){

            let object =null;
            let r = this.objects[0].radius;

            let x = cuboid.x - 2*r, dx = - cuboid.x/2 + center.x;
            let y = cuboid.y - 2*r, dy = - cuboid.y/2 + center.y;
            let z = cuboid.z - 2*r, dz = - cuboid.z/2 + center.z;

            // everything is centered around 0 - so we shift half
            for (let i=0; i< this.objects.length; i++) {

                object = this.objects[i];
                object.pos.set( x * Math.random() + dx, 
                                y * Math.random() + dy, 
                                z * Math.random() + dz);
            }
        }

        slideObjects(container, ballBox, location) {

            // check
            if (this.objects.length < 1) return

            let pos = null;
            let r = this.objects[0].radius;
            let max = container.bBox.max, min = container.bBox.min;

            // the displacement of the box
            let slide = (max.x - min.x - ballBox.vector.x) * location.delta/location.max;
            for (let i=0; i< this.objects.length; i++) {
                pos = this.objects[i].pos;
                pos.x += slide;
                if (pos.x < min.x + r ) pos.x = min.x + r;
                if (pos.x > max.x - r ) pos.x = max.x - r;
            }
        }

        orderedInCuboid(cuboid){

            // the nr of vectors we need
            let nVectors = this.objects.length;
            let r = this.objects[0].radius;

            // the raster we will use - note that nx*ny*nz > nVectors !
            let nx = Math.ceil( (nVectors*(cuboid.x/cuboid.y)*(cuboid.x/cuboid.z))**(1/3));
            let ny = Math.ceil( (nVectors*(cuboid.y/cuboid.x)*(cuboid.y/cuboid.z))**(1/3));
            let nz = Math.ceil( (nVectors*(cuboid.z/cuboid.x)*(cuboid.z/cuboid.y))**(1/3));
        
            // the spacing between objects
            // let dx=3*r,dy=3*r,dz=3*r
            // if (dx*nx > cuboid.x) dx = cuboid.x/nx
            // if (dy*ny > cuboid.y) dy = cuboid.y/ny
            // if (dz*nz > cuboid.z) dz = cuboid.z/nz

            let dx = cuboid.x/nx;
            let dy = cuboid.y/ny;
            let dz = cuboid.z/nz;
        
            // the indices
            let ix = 0, iy = 0, iz = 0, count = 0;

            // we center around 0
            let offset = {x: -dx*(nx-1)/2, y: -dy*(ny-1)/2, z: -dz*(nz-1)/2};
     
            // caculate the positions of all objects
            for (iz=0; iz<nz; iz++) {
                for (iy=0; iy<ny; iy++) {
                    for (ix=0; ix<nx; ix++) {              
                        this.objects[count].pos.set(ix*dx + offset.x, 
                                                    iy*dy + offset.y, 
                                                    iz*dz + offset.z);
                        if (++count == nVectors) return
                    }
                }
            }
        }

        randomInSphere(radius) {

            let phi=0.0, theta=0.0, R=0.0;
            for (let i=0; i< this.objects.length; i++) {
                phi   = Math.random()*Math.PI*2;
                theta = Math.random()*Math.PI;
                R     = Math.random()*radius;
                this.objects[i].pos.set(R*Math.sin(theta)*Math.cos(phi), 
                                        R*Math.sin(theta)*Math.sin(phi),
                                        R*Math.cos(theta));
            }
        }

        orderedInSphere(radius) {

            let R=0.0;
            let dphi=2*Math.PI;
            let dtheta=Math.PI;
            let r = 2*this.objects[0].radius;
            let count = 0;

            for (let shell=1; shell < 10; shell++) {

                R = shell * r;
                dphi = 2*Math.PI/(shell*6);
                dtheta = Math.PI/(shell*6);
                for (let i = 0; i <= shell*6; i++) {
                    for (let j=0; j<= shell*6;j++) {

                        this.objects[count].pos.set(
                            R*Math.sin(dtheta*j)*Math.cos(dphi*i), 
                            R*Math.sin(dtheta*j)*Math.sin(dphi*i),
                            R*Math.cos(dtheta*j));

                        if (++count == this.objects.length) return
                    }
                }
            }
        }

        placeObjectsInScene(scene) {

            // the matrix to position the objects
            let matrix = new THREE.Matrix4();

            // check if we have to place a new mesh
            if (this.mesh == null) {

                // create an instanced mesh with the geometry - it will be instanced nObjects time
                this.mesh = new THREE.InstancedMesh(this.geometry, this.material, this.objects.length);

                // add to the scene
                scene.add(this.mesh);
            }

            // set the update bit
    		this.mesh.instanceMatrix.needsUpdate = true;

            // place every instance...
            for (let i=0; i<this.objects.length; i++) {

                // make the matrix to position the object
                matrix.setPosition(this.objects[i].pos);

                // position the object
                this.mesh.setMatrixAt(i,matrix);
            }
        }

        // {   text:"fixed value",         value:{type:"fixed",fixedSpeed}},
        // {   text:"normal distributed",  value:{type:"rnd-normal", mu, sigma}},
        // {   text:"uniform distributed", value:{type:"rnd-uniform",a, b}},

        initSpeed(direction, amplitude) {

            let object = null;
            let A=0, fixedSpeed=0, sigma=0, mu=0, a=0, b=0;
            let phi = 0.0, theta = 0.0;

            // to speed up the calculations
            switch(amplitude.type) {

                case "fixed":  
                    fixedSpeed = +amplitude.fixedSpeed.value;
                    break

                case "rnd-normal": 
                    sigma = +amplitude.sigma.value; 
                    mu = +amplitude.mu.value;
                    break

                case "rnd-uniform": 
                    a = +amplitude.a.value;
                    b = +amplitude.b.value;
                    break
            }
            
            // calculate for each object
            for(let i=0;i<this.objects.length;i++) {

                // notation
                object = this.objects[i];

                // first get the amplitue
                switch(amplitude.type) {
                    case "fixed": A = fixedSpeed;
                        break

                    case "rnd-normal": A = normalRandom()*sigma + mu;
                        break

                    case "rnd-uniform": A = (Math.random() * (b-a)) + a;
                        break
                }

                // then the vector
                switch(direction.type) {
            
                    case "random":
                        // get the random vectors using polar coordinates
                        phi = Math.random()*Math.PI*2;
                        theta = Math.random()*Math.PI;
                        object.speed.set(   A*Math.sin(theta)*Math.cos(phi), 
                                            A*Math.sin(theta)*Math.sin(phi),
                                            A*Math.cos(theta));
                        break
            
                    case "radial":
                        object.speed.copy(object.pos).sub(direction.center).setLength(A);
                        break
            
                    case "parallel":
                        object.speed.copy(direction.parallel).setLength(A);
                        break
            
                    default: 
                        console.log("[UNKNOWN] unknown speed vector type in function speedVectors", direction.type);
                        return
                }
            }
            //console.log("objects",this.objects )
        }

        // put the helpers in a THREE.Group !
        showAABB(scene) {

            let object = null;

            // check if we have to create the box helper group
            if (this.boxHelperGroup == null ) {

                this.boxHelperGroup = new THREE.Group();

                for (let i=0; i < this.objects.length; i++) {

                    // notation
                    object = this.objects[i];

                    // create a boxhelper object
                    object.boxHelper = new THREE.Box3Helper( object.bBox, boxFreeClr$1 );

                    // add the helpers to the box helper group
                    this.boxHelperGroup.add(object.boxHelper);
                }
            }

            // add the helper group to the scene
            scene.add(this.boxHelperGroup);
        }

        hideAABB(scene) {

            if (this.boxHelperGroup) scene.remove( this.boxHelperGroup );
        }

        changeMass(value, unit) {

            // convert the mass of the object to g
            let mass = (+value) * factor("mass",unit,"g");

            // change all masses
            for (let i=0; i < this.objects.length; i++) 
                this.objects[i].mass = mass;
        }

        changeRadius(value, unit) {

            // recalculate the radius
            let radius = (+value) * factor("length",unit,"m");

            // we have to change the geometry
            this.geometry.dispose();

            // get a new geometry
            this.geometry = new THREE.SphereGeometry(radius , 12, 12); 
            
            // and put it in the mesh
            this.mesh.geometry = this.geometry;

            // change all objects
            for (let i=0; i < this.objects.length; i++) 
                this.objects[i].radius = radius;
        }
    }

    class octreeClass {

        // creates the octant-tree recursively
        constructor(p1, p2, level) {

            // the octant is a box defined by two diagonal points
            this.space = new THREE.Box3(p1, p2);

            // octant is an array that will contain the sub-octants if any
            this.octant = null; //[null, null, null, null, null, null, null, null]

            // the linked list that will contain all objects in this octant that do not fit in a sub-octant
            this.objectList = new singleLinkedList;

            // if not at the lowest level (note 1 = 8, 2 = 64, 3 = 512 suboctants)
            if (level > 0) {

                // it contains 8 new octants
                this.octant = [null, null, null, null, null, null, null, null];

                // dp is a vector with the size of the octant at this level
                let dp = new THREE.Vector3( (p2.x-p1.x)/2,(p2.y-p1.y)/2,(p2.z-p1.z)/2 ); 
                
                // get a new vector
                let p = new THREE.Vector3();

                // make a total of eight new octants
                let n = 0;
                for (let i=0; i<2; i++)
                    for (let j=0; j<2; j++)
                        for (let k=0; k<2; k++) {

                            // p defines one point of the octant..
                            p.set( p1.x + i*dp.x, p1.y + j*dp.y , p1.z + k*dp.z);

                            // ..and p+dp is the diagonally opposed point 
                            this.octant[n++] = new octreeClass(p, p.clone().add( dp ), level-1);
                        }
            }
        }

        // when the size of the container changes, the sizes of the octants have to be adapted as well
        sizeChange(p1,p2) {

            // change the size of the box
            this.space.min = {...p1};
            this.space.max = {...p2};

            // if there are there still sub-octants..
            if (this.octant) {

                // dp is a vector with the size of the octant at this level
                let dp = new THREE.Vector3( (p2.x-p1.x)/2,(p2.y-p1.y)/2,(p2.z-p1.z)/2 ); 
                
                // a copy of p1
                let p = new THREE.Vector3();

                // adjust the size for for all octants
                let n = 0;
                for (let i=0; i<2; i++)
                    for (let j=0; j<2; j++)
                        for (let k=0; k<2; k++) {

                            // p defines one point of the octant..
                            p.set( p1.x + i*dp.x, p1.y + j*dp.y , p1.z + k*dp.z);

                            // ..and p+dp is the diagonally opposed point 
                            this.octant[n++].sizeChange(p, p.clone().add( dp ));
                        }
            }
        }

        // all the object lists are reset
        reset() {
            // reset the object list
            this.objectList.reset();

            // check the sub-octants
            if (this.octant) 
                for (let i=0; i < this.octant.length; i++) 
                    if ( this.octant[i]) this.octant[i].reset();
        }

        // object is checked against the other objects already in the octree
        // secTime is the duration of the interval in seconds 
        // collisions is an array that will contain all collision events in the interval
        collisionCheck( object, secTime, collisions) {

            // take the first object on the objects in this octant
            let other = this.objectList.first;
            let when = 0.0;

            // check against collisions with objects stored at this level
            while (other != null ) {

                // calculate if the objects intersect and if so at what time in the interval - inifinity means no collision
                if ( ( when = object.intersects(other, secTime)) != Infinity ) collisions.add(object,other,when);
                
                // take the next object
                other = other.next;
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
            this.objectList.push(object);

            // we're done
            return true            
        }

        // count the nr of objects in the octree - just for debugging
        count() {
            let n = this.objectList.length();

            if (this.octant) 
                for (let i=0; i<this.octant.length; i++)  
                    n += this.octant[i].count();
            return n
        }

        // count the total nr of sub-octants - just for debugging
        size() {
            let n = 1;

            if (this.octant) {
                for (let i=0; i<this.octant.length; i++)  
                    n += this.octant[i].size();
            }
            return n
        }
    }// end of octree class

    class studioClass {

    constructor(canvas) {

    	// the width and height of the canvas
    	let W = canvas.clientWidth;
    	let H = canvas.clientHeight;

    	let nearPlane=1, farPlane=10000;

    	// the camera
    	this.camera = new THREE.PerspectiveCamera(45,W / H, nearPlane, farPlane);
    	this.camera.position.set(1,1,1);
    	this.camera.lookAt(0,0,0);

    	//the scene
    	this.scene = new THREE.Scene();

    	// we will create an ambient and a point light
    	this.lights = {ambient:null, point:null};

    	// add an ambient light
    	this.lights.ambient = new THREE.AmbientLight(0xffffff, 0.2);
    	this.scene.add(this.lights.ambient);

    	// add a point light
    	this.lights.point = new THREE.PointLight(0xffffff, 0.9, 0, 2);
    	this.lights.point.position.set(0,0,0);
    	this.scene.add(this.lights.point);
    	
    	// set up the renderer
    	this.renderer = new THREE.WebGLRenderer({canvas: canvas, antialias: true});
    	this.renderer.setClearColor(0x000000);
    	this.renderer.setPixelRatio(window.devicePixelRatio);
    	this.renderer.setSize(W, H);
    	
    	// set up the orbit controls
    	// NOTE BAD REACTION FROM ORBITCONTROLS WHEN HEIGHT +  VH
    	this.orbitControls = new OrbitControls_js.OrbitControls( this.camera, canvas );
    	this.orbitControls.minDistance = nearPlane;
    	this.orbitControls.maxDistance = farPlane;
    	this.orbitControls.maxPolarAngle = Math.PI;

        // add an axes helper to the scene
    	this.scene.add( new THREE.AxesHelper( 1 ) );
    }

    adjustToScene(d) {
    	if (this.camera) {
    		this.camera.position.set(d/2,d/2,d);
    		this.camera.lookAt(0,0,0);
    	}
    	if (this.lights.point) this.lights.point.position.set( 2*d, 2*d, 2*d);
    }

    render() {

    	// just render the scene
    	this.renderer.render(this.scene, this.camera);
    }

    }//end of studio class

    class _3DTimerClass {

    	/**
    	 * fTime is the start of the frame
    	 * uTime is the start of the update
    	 * rTime is the start of the render
    	 * 
    	 * flextime means update will use real time difference between two frames (could 17ms but also 20ms etc.)
    	 * if flextime = false then each frame is considered to have the same fixed length, eg 16.666 ms
    	 * 
    	 * The timefactor is how much one framesecond is in real time - eg timefactor = 0.1 slows down the action by 10
    	 */
         constructor(fixedLength = 0.0, timeFactor = 1.0, reportCount = 600) {

    		//if a fixed frame time is given (in seconds)
    		this.flexTime = fixedLength > 0.0 ? false : true;

    		// save the fixed time
    		this.fixedLength = fixedLength;

    		// the time factor
    		this.timeFactor = timeFactor;

    		// the integral, update and render time
    		this.fTime  = this.uTime  = this.rTime = 0.0;
    		
    		// the time of one loop
    		this.fDelta = this.uDelta = this.rDelta = 0.0;
    		
    		// the sum 
    		this.fCumul = this.uCumul = this.rCumul = 0.0;
    		
    		// the nr of times render/update were done
            this.fCount = this.uCount = this.rCount = 0;

            // report frame count
            this.reportCount = reportCount;
        }
        
    reset() {
    	// set the counts to zero
    	this.fCount = this.uCount = this.rCount = 0;
    	
    	// reset the cumulative time
    	this.fCumul = this.uCumul = this.rCumul = 0;

    	// set the starting time
    	this.fTime = window.performance.now();
    }

    frameStart() {
    	// get the time
    	let now = window.performance.now();

    	// Get the duration of the previous loop
    	this.fCumul += (this.fDelta = now - this.fTime);

    	// save the start time of the integral loop / also the start time of the update loop
    	this.fTime = now;

    	// counter
    	this.fCount++;

    	// the time lapsed since the last update - used by update
    	if (this.flexTime) 
    		return this.fDelta * this.timeFactor
    	else 
    		return this.fixedLength * this.timeFactor
    }

    beforeUpdate() {
    	// get the update start time
    	this.uTime = window.performance.now();
    }

    afterUpdate(){
    	// get the duration of the update loop
    	this.uCumul += (this.uDelta = window.performance.now() - this.uTime);

    	// counter
    	this.uCount++;
    }

    beforeRender() {
    	// get the render start time
    	this.rTime =  window.performance.now();
    }

    afterRender(){
    	// get the duraion of the render loop
    	this.rCumul += (this.rDelta = window.performance.now() - this.rTime);

    	// counter
        this.rCount++;

    	// report if required
    	if (this.fCount > this.reportCount) this.reportTiming();
    }

    reportTiming() {

    	// calculate the mean values
    	let f = this.fCount > 0 ? this.fCumul/this.fCount : 0.001;
    	let u = this.uCount > 0 ? this.uCumul/this.uCount : 0.0;
    	let r = this.rCount > 0 ? this.rCumul/this.rCount : 0.0;

        // report the timings
        console.log(`[TIMING] ${(1000/f).toFixed(0)} fps - frame ${f.toFixed(3)} update ${u.toFixed(3)} render ${r.toFixed(3)}`);

        // set the counts and timings to zero
        this.fCumul = this.uCumul = this.rCumul = 0.0;
        this.fCount = this.uCount = this.rCount = 0;
    }

    } // end of timer class

    exports._3DTimerClass = _3DTimerClass;
    exports.collisionClass = collisionClass;
    exports.cuboidContainerClass = cuboidContainerClass;
    exports.objectClass = objectClass;
    exports.objectGroupClass = objectGroupClass;
    exports.octreeClass = octreeClass;
    exports.studioClass = studioClass;

    Object.defineProperty(exports, '__esModule', { value: true });

});
//# sourceMappingURL=space.js.map
