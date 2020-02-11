import * as THREE from 'three/build/three.module.js'
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js'

export class studioClass {

constructor(canvas) {

	// the width and height of the canvas
	let W = canvas.clientWidth
	let H = canvas.clientHeight

	let nearPlane=1, farPlane=10000

	// the camera
	this.camera = new THREE.PerspectiveCamera(45,W / H, nearPlane, farPlane)
	this.camera.position.set(1,1,1)
	this.camera.lookAt(0,0,0)

	//the scene
	this.scene = new THREE.Scene();

	// we will create an ambient and a point light
	this.lights = {ambient:null, point:null}

	// add an ambient light
	this.lights.ambient = new THREE.AmbientLight(0xffffff, 0.2)
	this.scene.add(this.lights.ambient)

	// add a point light
	this.lights.point = new THREE.PointLight(0xffffff, 0.9, 0, 2);
	this.lights.point.position.set(0,0,0)
	this.scene.add(this.lights.point)
	
	// set up the renderer
	this.renderer = new THREE.WebGLRenderer({canvas: canvas, antialias: true})
	this.renderer.setClearColor(0x000000)
	this.renderer.setPixelRatio(window.devicePixelRatio);
	this.renderer.setSize(W, H)
	
	// set up the orbit controls
	// NOTE BAD REACTION FROM ORBITCONTROLS WHEN HEIGHT +  VH
	this.orbitControls = new OrbitControls( this.camera, canvas )
	this.orbitControls.minDistance = nearPlane
	this.orbitControls.maxDistance = farPlane
	this.orbitControls.maxPolarAngle = Math.PI;

    // add an axes helper to the scene
	this.scene.add( new THREE.AxesHelper( 1 ) )
}

adjustToScene(d) {
	if (this.camera) {
		this.camera.position.set(d/2,d/2,d)
		this.camera.lookAt(0,0,0)
	}
	if (this.lights.point) this.lights.point.position.set( 2*d, 2*d, 2*d)
}

render() {

	// just render the scene
	this.renderer.render(this.scene, this.camera);
}

}//end of studio class