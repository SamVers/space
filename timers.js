export class _3DTimerClass {

	//  fTime is the start of the frame
	//  uTime is the start of the update
	//  rTime is the start of the render
	 
	//  flextime means update will use real time difference between two frames (could 17ms but also 20ms etc.)
	//  if flextime = false then each frame is considered to have the same fixed length, eg 16.666 ms
	 
	//  The timefactor is how much one framesecond is in real time - eg timefactor = 0.1 slows down the action by 10
	 
    constructor(fixedLength = 0.0, timeFactor = 1.0, reportCount = 600) {

		//if a fixed frame time is given (in seconds)
		this.flexTime = fixedLength > 0.0 ? false : true

		// save the fixed time
		this.fixedLength = fixedLength

		// the time factor
		this.timeFactor = timeFactor

		// the integral, update and render time
		this.fTime  = this.uTime  = this.rTime = 0.0
		
		// the time of one loop
		this.fDelta = this.uDelta = this.rDelta = 0.0
		
		// the sum 
		this.fCumul = this.uCumul = this.rCumul = 0.0
		
		// the nr of times render/update were done
        this.fCount = this.uCount = this.rCount = 0

        // report frame count
        this.reportCount = reportCount
    }
    
reset() {
	// set the counts to zero
	this.fCount = this.uCount = this.rCount = 0
	
	// reset the cumulative time
	this.fCumul = this.uCumul = this.rCumul = 0

	// set the starting time
	this.fTime = window.performance.now()
}

frameStart() {
	// get the time
	let now = window.performance.now()

	// Get the duration of the previous loop
	this.fCumul += (this.fDelta = now - this.fTime)

	// save the start time of the integral loop / also the start time of the update loop
	this.fTime = now

	// counter
	this.fCount++

	// the time lapsed since the last update - used by update
	if (this.flexTime) 
		return this.fDelta * this.timeFactor
	else 
		return this.fixedLength * this.timeFactor
}

beforeUpdate() {
	// get the update start time
	this.uTime = window.performance.now()
}

afterUpdate(){
	// get the duration of the update loop
	this.uCumul += (this.uDelta = window.performance.now() - this.uTime)

	// counter
	this.uCount++
}

beforeRender() {
	// get the render start time
	this.rTime =  window.performance.now()
}

afterRender(){
	// get the duraion of the render loop
	this.rCumul += (this.rDelta = window.performance.now() - this.rTime)

	// counter
    this.rCount++

	// report if required
	if (this.fCount > this.reportCount) this.reportTiming()
}

reportTiming() {

	// calculate the mean values
	let f = this.fCount > 0 ? this.fCumul/this.fCount : 0.001
	let u = this.uCount > 0 ? this.uCumul/this.uCount : 0.0
	let r = this.rCount > 0 ? this.rCumul/this.rCount : 0.0

    // report the timings
    console.log(`[TIMING] ${(1000/f).toFixed(0)} fps - frame ${f.toFixed(3)} update ${u.toFixed(3)} render ${r.toFixed(3)}`)

    // set the counts and timings to zero
    this.fCumul = this.uCumul = this.rCumul = 0.0
    this.fCount = this.uCount = this.rCount = 0
}

} // end of timer class