(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){
'use strict'

module.exports = IntegratorFactory

function defaultErrorScaleFunction( i, dt, y, dydt ) {
  return Math.abs(y) + Math.abs(dt * dydt) + 1e-32
}

function defaultErrorReduceFunction( i, accumulatedError, errorEstimate ) {
  return Math.max( accumulatedError, Math.abs(errorEstimate))
}

function defaultErrorPostFunction( accumulatedError ) {
  return accumulatedError
}

function minMag (a, b) {
  return (a > 0 ? Math.min : Math.max)(a, b);
}

function maxMag (a, b) {
  return (a > 0 ? Math.max : Math.min)(a, b);
}

var Integrator = function Integrator( y0, deriv, t0, dt0, options ) {
  var opts = options || {}
  this.tol = opts.tol===undefined ? 1e-8 : opts.tol
  this.maxIncreaseFactor = opts.maxIncreaseFactor===undefined ? 10 : opts.maxIncreaseFactor
  this.maxDecreaseFactor = opts.maxDecreaseFactor===undefined ? 10 : opts.maxDecreaseFactor
  this.dtMinMag = opts.dtMinMag===undefined ? 0 : Math.abs(opts.dtMinMag)
  this.dtMaxMag = opts.dtMaxMag===undefined ? undefined : Math.abs(opts.dtMaxMag)
  this.verbose = opts.verbose===undefined ? true : !!opts.verbose;

  var logCnt = 0
  var maxLogs = 10
  var maxLogWarningIssued = false
  this.__log = function (method, msg) {
    if (!this.verbose) return;
    if (logCnt < maxLogs) {
      console.log('ode45-cash-karp::' + method + '(): ' + msg)
      logCnt++
    } else {
      if (!maxLogWarningIssued) {
        console.log('ode45-cash-karp: too many warnings. Silencing further output')
        maxLogWarningIssued = true
      }
    }
  }.bind(this)

  this.errorScaleFunction = opts.errorScaleFunction === undefined ? defaultErrorScaleFunction : opts.errorScaleFunction
  this.errorReduceFunction = opts.errorReduceFunction === undefined ? defaultErrorReduceFunction : opts.errorReduceFunction
  this.errorPostFunction = opts.errorPostFunction === undefined ? defaultErrorPostFunction : opts.errorPostFunction

  // This is technically a parameter, but I think the value of leaving this undocumented exceeds the
  // value of documenting this and only adding confusion. I can't imagine this will even need to be
  // modified.
  this.safetyFactor = opts.safetyFactor===undefined ? 0.9 : opts.safetyFactor

  // Bind variables to this:
  this.deriv = deriv
  this.y = y0
  this.n = this.y.length
  this.dt = dt0
  this.t = t0

  // Create a scratch array into which we compute the derivative:
  this._ctor = this.y.constructor

  this._errorScale = new this._ctor( this.n )
  this._w = new this._ctor( this.n )
  this._k1 = new this._ctor( this.n )
  this._k2 = new this._ctor( this.n )
  this._k3 = new this._ctor( this.n )
  this._k4 = new this._ctor( this.n )
  this._k5 = new this._ctor( this.n )
  this._k6 = new this._ctor( this.n )
}

Integrator.prototype._calculateK1 = function() {
  this.deriv( this._k1, this.y, this.t )

  return this
}

Integrator.prototype._calculateKs = function(dt) {
  var i

  //var a21 =  0.200000000000000000 // 1/5
  //var a31 =  0.075000000000000000 // 3/40
  //var a32 =  0.225000000000000000 // 9/40
  //var a41 =  0.300000000000000000 // 3/10
  //var a42 = -0.900000000000000000 // -9/10
  //var a43 =  1.200000000000000000 // 6/5
  //var a51 = -0.203703703703703703 // -11/54
  //var a52 =  2.500000000000000000 // 5/2
  //var a53 = -2.592592592592592592 // -70/27
  //var a54 =  1.296296296296296296 // 35/27
  //var a61 =  0.029495804398148148 // 1631/55296
  //var a62 =  0.341796875000000000 // 175/512
  //var a63 =  0.041594328703703703 // 575/13824
  //var a64 =  0.400345413773148148 // 44275/110592
  //var a65 =  0.061767578125000000 // 253/4096

  //var b1  =  0.000000000000000000 // 0
  //var b2  =  0.200000000000000000 // 1/5
  //var b3  =  0.300000000000000000 // 3/10
  //var b4  =  0.600000000000000000 // 3/5
  //var b5  =  1.000000000000000000 // 1
  //var b6  =  0.875000000000000000 // 7/8

  // Same for every step, so don't repeat:
  //this.deriv( this._k1, this.y, this.t )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
      0.2 * this._k1[i] )
  }

  this.deriv( this._k2, this._w, this.t + dt * 0.2 )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
      0.075 * this._k1[i] +
      0.225 * this._k2[i] )
  }

  this.deriv( this._k3, this._w, this.t + dt * 0.3 )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
       0.3 * this._k1[i] +
      -0.9 * this._k2[i] +
       1.2 * this._k3[i] )
  }

  this.deriv( this._k4, this._w, this.t + dt * 0.6 )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
      -0.203703703703703703 * this._k1[i] +
       2.5                  * this._k2[i] +
      -2.592592592592592592 * this._k3[i] +
       1.296296296296296296 * this._k4[i] )
  }

  this.deriv( this._k5, this._w, this.t + dt /* * b5 */ )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
      0.029495804398148148 * this._k1[i] +
      0.341796875          * this._k2[i] +
      0.041594328703703703 * this._k3[i] +
      0.400345413773148148 * this._k4[i] +
      0.061767578125       * this._k5[i] )
  }

  this.deriv( this._k6, this._w, this.t + dt * 0.875 )

  return this
}

Integrator.prototype._calculateError = function(dt) {
  //var cs1 =  0.102177372685185185 // 2825/27648
  //var cs2 =  0.000000000000000000 // 0
  //var cs3 =  0.383907903439153439 // 18575/48384
  //var cs4 =  0.244592737268518518 // 13525/55296
  //var cs5 =  0.019321986607142857 // 277/14336
  //var cs6 =  0.250000000000000000 // 1/4

  //var dc1 =  0.004293774801587301 // cs1 - c1
  //var dc2 =  0.000000000000000000 // cs2 - c2
  //var dc3 = -0.018668586093857832 // cs3 - c3
  //var dc4 =  0.034155026830808080 // cs4 - c4
  //var dc5 =  0.019321986607142857 // cs5 - c5
  //var dc6 = -0.039102202145680406 // cs6 - c6

  var error = 0
  for(var i=0; i<this.n; i++) {
    error =  this.errorReduceFunction( i, error,
      dt * (
         0.004293774801587301 * this._k1[i] +
        -0.018668586093857832 * this._k3[i] +
         0.034155026830808080 * this._k4[i] +
         0.019321986607142857 * this._k5[i] +
        -0.039102202145680406 * this._k6[i]
      ) / this._errorScale[i]
    )
  }

  return this.errorPostFunction(error)
}

Integrator.prototype._update = function(dt) {
  //var c1  =  0.097883597883597883 // 37/378
  //var c2  =  0.000000000000000000 // 0
  //var c3  =  0.402576489533011272 // 250/621
  //var c4  =  0.210437710437710437 // 125/594
  //var c5  =  0.000000000000000000 // 0
  //var c6  =  0.289102202145680406 // 512/1771

  for(var i=0; i<this.n; i++) {
    this.y[i] += dt * (
      0.097883597883597883 * this._k1[i] +
      0.402576489533011272 * this._k3[i] +
      0.210437710437710437 * this._k4[i] +
      0.289102202145680406 * this._k6[i]
    )
  }
  this.t += dt
  return this
}

Integrator.prototype._calculateErrorScale = function(dt) {
  for(var i=0; i<this.n; i++) {
    this._errorScale[i] = this.errorScaleFunction(i, dt, this.y[i], this._k1[i])
  }
  return this
}

Integrator.prototype.step = function( tLimit ) {
  // Bail out early if we're *at* the limit:
  if (Math.abs(this.t - tLimit) < this.dt * 1e-10) {
    return false;
  }

  var thisDt = this.dt;

  // Don't integrate past a tLimit, if provided:
  if( tLimit !== undefined ) {
    thisDt = thisDt > 0 ? Math.min( tLimit - this.t, thisDt ) : Math.max( tLimit - this.t, thisDt )
  }

  // Limit the magnitude of dt to dtMaxMag
  if( this.dtMaxMag !== undefined && Math.abs( thisDt ) > this.dtMaxMag ) {
    this.__log('step', 'step greater than maximum stepsize requested. dt magnitude has been limited.')
    thisDt = thisDt > 0 ? this.dtMaxMag : -this.dtMaxMag
  }

  // Limit the magnitude of dt to dtMinMag
  if( this.dtMinMag !== undefined && Math.abs( thisDt ) < this.dtMinMag ) {
    this.__log('step', 'step smaller than minimum stepsize requested. dt magnitude has been limited.')
    thisDt = thisDt > 0 ? this.dtMinMag : -this.dtMinMag
  }

  // The first derivative doesn't change even if dt does, so only calculate this once:
  this._calculateK1()

  // The scale factor per-dimension probably doesn't need to change either across a single adaptive step:
  this._calculateErrorScale(thisDt)

  var error = Infinity
  var maxError = 0
  var nextDt
  var lowerDtLimitReached = false

  while(true) {

    // Calculate intermediate k's for the proposed step:
    this._calculateKs(thisDt)

    // Calculate the max error of the proposed step:
    error = this._calculateError(thisDt)

    if( error < this.tol || lowerDtLimitReached ) {
      // Success! Exit:
      break
    }

    if( ! Number.isFinite(error) ) {
      throw new Error('ode45-cash-karp::step() NaN encountered while integrating.')
    }

    // Failure. Adapt the timestep:
    nextDt = this.safetyFactor * thisDt * Math.pow( this.tol / error, 0.2 )

    // Cut the timestep, but not by more than maxDecreaseFactor
    thisDt = maxMag( thisDt / this.maxDecreaseFactor, nextDt )

    // If stepsize too small, finish off by taking the currently proposed step and logging a warning:
    if( this.dtMinMag !== undefined && Math.abs(thisDt) < this.dtMinMag ) {
      thisDt = this.dtMinMag * (thisDt > 0 ? 1 : -1);
      this.__log('step', 'minimum stepsize reached.')
      lowerDtLimitReached = true
    }
  }

  // Apply this update:
  this._update(thisDt)

  // Calculate the next timestep size:
  nextDt = this.safetyFactor * thisDt * Math.pow( this.tol / error, 0.25 )

  // Increase the timestep for the next time around, but not by more than the maxIncreaseFactor:
  this.dt = maxMag(this.dt / this.maxDecreaseFactor, minMag( this.dt * this.maxIncreaseFactor, nextDt ));

  if( tLimit !== undefined ) {
    return Math.abs(this.t - tLimit) > this.dt * 1e-8;
  } else {
    return true
  }
}

Integrator.prototype.steps = function( n, tLimit ) {
  for(var step=0; step<n; step++) {
    if( ! this.step(tLimit) ) return false;
  }
}

function IntegratorFactory( y0, deriv, t, dt, options ) {
  return new Integrator( y0, deriv, t, dt, options )
}

},{}],2:[function(require,module,exports){
const solver = require('ode45-cash-karp');
window.require = require;

},{"ode45-cash-karp":1}]},{},[2])
//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIkM6L1VzZXJzL2JlbmphL0FwcERhdGEvTG9jYWwvWWFybi9EYXRhL2dsb2JhbC9ub2RlX21vZHVsZXMvYnJvd3Nlci1wYWNrL19wcmVsdWRlLmpzIiwibm9kZV9tb2R1bGVzL29kZTQ1LWNhc2gta2FycC9saWIvaW5kZXguanMiLCJvZGU0NS5qcyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiQUFBQTtBQ0FBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN2VEE7QUFDQSIsImZpbGUiOiJnZW5lcmF0ZWQuanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlc0NvbnRlbnQiOlsiKGZ1bmN0aW9uKCl7ZnVuY3Rpb24gcihlLG4sdCl7ZnVuY3Rpb24gbyhpLGYpe2lmKCFuW2ldKXtpZighZVtpXSl7dmFyIGM9XCJmdW5jdGlvblwiPT10eXBlb2YgcmVxdWlyZSYmcmVxdWlyZTtpZighZiYmYylyZXR1cm4gYyhpLCEwKTtpZih1KXJldHVybiB1KGksITApO3ZhciBhPW5ldyBFcnJvcihcIkNhbm5vdCBmaW5kIG1vZHVsZSAnXCIraStcIidcIik7dGhyb3cgYS5jb2RlPVwiTU9EVUxFX05PVF9GT1VORFwiLGF9dmFyIHA9bltpXT17ZXhwb3J0czp7fX07ZVtpXVswXS5jYWxsKHAuZXhwb3J0cyxmdW5jdGlvbihyKXt2YXIgbj1lW2ldWzFdW3JdO3JldHVybiBvKG58fHIpfSxwLHAuZXhwb3J0cyxyLGUsbix0KX1yZXR1cm4gbltpXS5leHBvcnRzfWZvcih2YXIgdT1cImZ1bmN0aW9uXCI9PXR5cGVvZiByZXF1aXJlJiZyZXF1aXJlLGk9MDtpPHQubGVuZ3RoO2krKylvKHRbaV0pO3JldHVybiBvfXJldHVybiByfSkoKSIsIid1c2Ugc3RyaWN0J1xuXG5tb2R1bGUuZXhwb3J0cyA9IEludGVncmF0b3JGYWN0b3J5XG5cbmZ1bmN0aW9uIGRlZmF1bHRFcnJvclNjYWxlRnVuY3Rpb24oIGksIGR0LCB5LCBkeWR0ICkge1xuICByZXR1cm4gTWF0aC5hYnMoeSkgKyBNYXRoLmFicyhkdCAqIGR5ZHQpICsgMWUtMzJcbn1cblxuZnVuY3Rpb24gZGVmYXVsdEVycm9yUmVkdWNlRnVuY3Rpb24oIGksIGFjY3VtdWxhdGVkRXJyb3IsIGVycm9yRXN0aW1hdGUgKSB7XG4gIHJldHVybiBNYXRoLm1heCggYWNjdW11bGF0ZWRFcnJvciwgTWF0aC5hYnMoZXJyb3JFc3RpbWF0ZSkpXG59XG5cbmZ1bmN0aW9uIGRlZmF1bHRFcnJvclBvc3RGdW5jdGlvbiggYWNjdW11bGF0ZWRFcnJvciApIHtcbiAgcmV0dXJuIGFjY3VtdWxhdGVkRXJyb3Jcbn1cblxuZnVuY3Rpb24gbWluTWFnIChhLCBiKSB7XG4gIHJldHVybiAoYSA+IDAgPyBNYXRoLm1pbiA6IE1hdGgubWF4KShhLCBiKTtcbn1cblxuZnVuY3Rpb24gbWF4TWFnIChhLCBiKSB7XG4gIHJldHVybiAoYSA+IDAgPyBNYXRoLm1heCA6IE1hdGgubWluKShhLCBiKTtcbn1cblxudmFyIEludGVncmF0b3IgPSBmdW5jdGlvbiBJbnRlZ3JhdG9yKCB5MCwgZGVyaXYsIHQwLCBkdDAsIG9wdGlvbnMgKSB7XG4gIHZhciBvcHRzID0gb3B0aW9ucyB8fCB7fVxuICB0aGlzLnRvbCA9IG9wdHMudG9sPT09dW5kZWZpbmVkID8gMWUtOCA6IG9wdHMudG9sXG4gIHRoaXMubWF4SW5jcmVhc2VGYWN0b3IgPSBvcHRzLm1heEluY3JlYXNlRmFjdG9yPT09dW5kZWZpbmVkID8gMTAgOiBvcHRzLm1heEluY3JlYXNlRmFjdG9yXG4gIHRoaXMubWF4RGVjcmVhc2VGYWN0b3IgPSBvcHRzLm1heERlY3JlYXNlRmFjdG9yPT09dW5kZWZpbmVkID8gMTAgOiBvcHRzLm1heERlY3JlYXNlRmFjdG9yXG4gIHRoaXMuZHRNaW5NYWcgPSBvcHRzLmR0TWluTWFnPT09dW5kZWZpbmVkID8gMCA6IE1hdGguYWJzKG9wdHMuZHRNaW5NYWcpXG4gIHRoaXMuZHRNYXhNYWcgPSBvcHRzLmR0TWF4TWFnPT09dW5kZWZpbmVkID8gdW5kZWZpbmVkIDogTWF0aC5hYnMob3B0cy5kdE1heE1hZylcbiAgdGhpcy52ZXJib3NlID0gb3B0cy52ZXJib3NlPT09dW5kZWZpbmVkID8gdHJ1ZSA6ICEhb3B0cy52ZXJib3NlO1xuXG4gIHZhciBsb2dDbnQgPSAwXG4gIHZhciBtYXhMb2dzID0gMTBcbiAgdmFyIG1heExvZ1dhcm5pbmdJc3N1ZWQgPSBmYWxzZVxuICB0aGlzLl9fbG9nID0gZnVuY3Rpb24gKG1ldGhvZCwgbXNnKSB7XG4gICAgaWYgKCF0aGlzLnZlcmJvc2UpIHJldHVybjtcbiAgICBpZiAobG9nQ250IDwgbWF4TG9ncykge1xuICAgICAgY29uc29sZS5sb2coJ29kZTQ1LWNhc2gta2FycDo6JyArIG1ldGhvZCArICcoKTogJyArIG1zZylcbiAgICAgIGxvZ0NudCsrXG4gICAgfSBlbHNlIHtcbiAgICAgIGlmICghbWF4TG9nV2FybmluZ0lzc3VlZCkge1xuICAgICAgICBjb25zb2xlLmxvZygnb2RlNDUtY2FzaC1rYXJwOiB0b28gbWFueSB3YXJuaW5ncy4gU2lsZW5jaW5nIGZ1cnRoZXIgb3V0cHV0JylcbiAgICAgICAgbWF4TG9nV2FybmluZ0lzc3VlZCA9IHRydWVcbiAgICAgIH1cbiAgICB9XG4gIH0uYmluZCh0aGlzKVxuXG4gIHRoaXMuZXJyb3JTY2FsZUZ1bmN0aW9uID0gb3B0cy5lcnJvclNjYWxlRnVuY3Rpb24gPT09IHVuZGVmaW5lZCA/IGRlZmF1bHRFcnJvclNjYWxlRnVuY3Rpb24gOiBvcHRzLmVycm9yU2NhbGVGdW5jdGlvblxuICB0aGlzLmVycm9yUmVkdWNlRnVuY3Rpb24gPSBvcHRzLmVycm9yUmVkdWNlRnVuY3Rpb24gPT09IHVuZGVmaW5lZCA/IGRlZmF1bHRFcnJvclJlZHVjZUZ1bmN0aW9uIDogb3B0cy5lcnJvclJlZHVjZUZ1bmN0aW9uXG4gIHRoaXMuZXJyb3JQb3N0RnVuY3Rpb24gPSBvcHRzLmVycm9yUG9zdEZ1bmN0aW9uID09PSB1bmRlZmluZWQgPyBkZWZhdWx0RXJyb3JQb3N0RnVuY3Rpb24gOiBvcHRzLmVycm9yUG9zdEZ1bmN0aW9uXG5cbiAgLy8gVGhpcyBpcyB0ZWNobmljYWxseSBhIHBhcmFtZXRlciwgYnV0IEkgdGhpbmsgdGhlIHZhbHVlIG9mIGxlYXZpbmcgdGhpcyB1bmRvY3VtZW50ZWQgZXhjZWVkcyB0aGVcbiAgLy8gdmFsdWUgb2YgZG9jdW1lbnRpbmcgdGhpcyBhbmQgb25seSBhZGRpbmcgY29uZnVzaW9uLiBJIGNhbid0IGltYWdpbmUgdGhpcyB3aWxsIGV2ZW4gbmVlZCB0byBiZVxuICAvLyBtb2RpZmllZC5cbiAgdGhpcy5zYWZldHlGYWN0b3IgPSBvcHRzLnNhZmV0eUZhY3Rvcj09PXVuZGVmaW5lZCA/IDAuOSA6IG9wdHMuc2FmZXR5RmFjdG9yXG5cbiAgLy8gQmluZCB2YXJpYWJsZXMgdG8gdGhpczpcbiAgdGhpcy5kZXJpdiA9IGRlcml2XG4gIHRoaXMueSA9IHkwXG4gIHRoaXMubiA9IHRoaXMueS5sZW5ndGhcbiAgdGhpcy5kdCA9IGR0MFxuICB0aGlzLnQgPSB0MFxuXG4gIC8vIENyZWF0ZSBhIHNjcmF0Y2ggYXJyYXkgaW50byB3aGljaCB3ZSBjb21wdXRlIHRoZSBkZXJpdmF0aXZlOlxuICB0aGlzLl9jdG9yID0gdGhpcy55LmNvbnN0cnVjdG9yXG5cbiAgdGhpcy5fZXJyb3JTY2FsZSA9IG5ldyB0aGlzLl9jdG9yKCB0aGlzLm4gKVxuICB0aGlzLl93ID0gbmV3IHRoaXMuX2N0b3IoIHRoaXMubiApXG4gIHRoaXMuX2sxID0gbmV3IHRoaXMuX2N0b3IoIHRoaXMubiApXG4gIHRoaXMuX2syID0gbmV3IHRoaXMuX2N0b3IoIHRoaXMubiApXG4gIHRoaXMuX2szID0gbmV3IHRoaXMuX2N0b3IoIHRoaXMubiApXG4gIHRoaXMuX2s0ID0gbmV3IHRoaXMuX2N0b3IoIHRoaXMubiApXG4gIHRoaXMuX2s1ID0gbmV3IHRoaXMuX2N0b3IoIHRoaXMubiApXG4gIHRoaXMuX2s2ID0gbmV3IHRoaXMuX2N0b3IoIHRoaXMubiApXG59XG5cbkludGVncmF0b3IucHJvdG90eXBlLl9jYWxjdWxhdGVLMSA9IGZ1bmN0aW9uKCkge1xuICB0aGlzLmRlcml2KCB0aGlzLl9rMSwgdGhpcy55LCB0aGlzLnQgKVxuXG4gIHJldHVybiB0aGlzXG59XG5cbkludGVncmF0b3IucHJvdG90eXBlLl9jYWxjdWxhdGVLcyA9IGZ1bmN0aW9uKGR0KSB7XG4gIHZhciBpXG5cbiAgLy92YXIgYTIxID0gIDAuMjAwMDAwMDAwMDAwMDAwMDAwIC8vIDEvNVxuICAvL3ZhciBhMzEgPSAgMC4wNzUwMDAwMDAwMDAwMDAwMDAgLy8gMy80MFxuICAvL3ZhciBhMzIgPSAgMC4yMjUwMDAwMDAwMDAwMDAwMDAgLy8gOS80MFxuICAvL3ZhciBhNDEgPSAgMC4zMDAwMDAwMDAwMDAwMDAwMDAgLy8gMy8xMFxuICAvL3ZhciBhNDIgPSAtMC45MDAwMDAwMDAwMDAwMDAwMDAgLy8gLTkvMTBcbiAgLy92YXIgYTQzID0gIDEuMjAwMDAwMDAwMDAwMDAwMDAwIC8vIDYvNVxuICAvL3ZhciBhNTEgPSAtMC4yMDM3MDM3MDM3MDM3MDM3MDMgLy8gLTExLzU0XG4gIC8vdmFyIGE1MiA9ICAyLjUwMDAwMDAwMDAwMDAwMDAwMCAvLyA1LzJcbiAgLy92YXIgYTUzID0gLTIuNTkyNTkyNTkyNTkyNTkyNTkyIC8vIC03MC8yN1xuICAvL3ZhciBhNTQgPSAgMS4yOTYyOTYyOTYyOTYyOTYyOTYgLy8gMzUvMjdcbiAgLy92YXIgYTYxID0gIDAuMDI5NDk1ODA0Mzk4MTQ4MTQ4IC8vIDE2MzEvNTUyOTZcbiAgLy92YXIgYTYyID0gIDAuMzQxNzk2ODc1MDAwMDAwMDAwIC8vIDE3NS81MTJcbiAgLy92YXIgYTYzID0gIDAuMDQxNTk0MzI4NzAzNzAzNzAzIC8vIDU3NS8xMzgyNFxuICAvL3ZhciBhNjQgPSAgMC40MDAzNDU0MTM3NzMxNDgxNDggLy8gNDQyNzUvMTEwNTkyXG4gIC8vdmFyIGE2NSA9ICAwLjA2MTc2NzU3ODEyNTAwMDAwMCAvLyAyNTMvNDA5NlxuXG4gIC8vdmFyIGIxICA9ICAwLjAwMDAwMDAwMDAwMDAwMDAwMCAvLyAwXG4gIC8vdmFyIGIyICA9ICAwLjIwMDAwMDAwMDAwMDAwMDAwMCAvLyAxLzVcbiAgLy92YXIgYjMgID0gIDAuMzAwMDAwMDAwMDAwMDAwMDAwIC8vIDMvMTBcbiAgLy92YXIgYjQgID0gIDAuNjAwMDAwMDAwMDAwMDAwMDAwIC8vIDMvNVxuICAvL3ZhciBiNSAgPSAgMS4wMDAwMDAwMDAwMDAwMDAwMDAgLy8gMVxuICAvL3ZhciBiNiAgPSAgMC44NzUwMDAwMDAwMDAwMDAwMDAgLy8gNy84XG5cbiAgLy8gU2FtZSBmb3IgZXZlcnkgc3RlcCwgc28gZG9uJ3QgcmVwZWF0OlxuICAvL3RoaXMuZGVyaXYoIHRoaXMuX2sxLCB0aGlzLnksIHRoaXMudCApXG5cbiAgZm9yKGk9MDsgaTx0aGlzLm47IGkrKykge1xuICAgIHRoaXMuX3dbaV0gPSB0aGlzLnlbaV0gKyBkdCAqIChcbiAgICAgIDAuMiAqIHRoaXMuX2sxW2ldIClcbiAgfVxuXG4gIHRoaXMuZGVyaXYoIHRoaXMuX2syLCB0aGlzLl93LCB0aGlzLnQgKyBkdCAqIDAuMiApXG5cbiAgZm9yKGk9MDsgaTx0aGlzLm47IGkrKykge1xuICAgIHRoaXMuX3dbaV0gPSB0aGlzLnlbaV0gKyBkdCAqIChcbiAgICAgIDAuMDc1ICogdGhpcy5fazFbaV0gK1xuICAgICAgMC4yMjUgKiB0aGlzLl9rMltpXSApXG4gIH1cblxuICB0aGlzLmRlcml2KCB0aGlzLl9rMywgdGhpcy5fdywgdGhpcy50ICsgZHQgKiAwLjMgKVxuXG4gIGZvcihpPTA7IGk8dGhpcy5uOyBpKyspIHtcbiAgICB0aGlzLl93W2ldID0gdGhpcy55W2ldICsgZHQgKiAoXG4gICAgICAgMC4zICogdGhpcy5fazFbaV0gK1xuICAgICAgLTAuOSAqIHRoaXMuX2syW2ldICtcbiAgICAgICAxLjIgKiB0aGlzLl9rM1tpXSApXG4gIH1cblxuICB0aGlzLmRlcml2KCB0aGlzLl9rNCwgdGhpcy5fdywgdGhpcy50ICsgZHQgKiAwLjYgKVxuXG4gIGZvcihpPTA7IGk8dGhpcy5uOyBpKyspIHtcbiAgICB0aGlzLl93W2ldID0gdGhpcy55W2ldICsgZHQgKiAoXG4gICAgICAtMC4yMDM3MDM3MDM3MDM3MDM3MDMgKiB0aGlzLl9rMVtpXSArXG4gICAgICAgMi41ICAgICAgICAgICAgICAgICAgKiB0aGlzLl9rMltpXSArXG4gICAgICAtMi41OTI1OTI1OTI1OTI1OTI1OTIgKiB0aGlzLl9rM1tpXSArXG4gICAgICAgMS4yOTYyOTYyOTYyOTYyOTYyOTYgKiB0aGlzLl9rNFtpXSApXG4gIH1cblxuICB0aGlzLmRlcml2KCB0aGlzLl9rNSwgdGhpcy5fdywgdGhpcy50ICsgZHQgLyogKiBiNSAqLyApXG5cbiAgZm9yKGk9MDsgaTx0aGlzLm47IGkrKykge1xuICAgIHRoaXMuX3dbaV0gPSB0aGlzLnlbaV0gKyBkdCAqIChcbiAgICAgIDAuMDI5NDk1ODA0Mzk4MTQ4MTQ4ICogdGhpcy5fazFbaV0gK1xuICAgICAgMC4zNDE3OTY4NzUgICAgICAgICAgKiB0aGlzLl9rMltpXSArXG4gICAgICAwLjA0MTU5NDMyODcwMzcwMzcwMyAqIHRoaXMuX2szW2ldICtcbiAgICAgIDAuNDAwMzQ1NDEzNzczMTQ4MTQ4ICogdGhpcy5fazRbaV0gK1xuICAgICAgMC4wNjE3Njc1NzgxMjUgICAgICAgKiB0aGlzLl9rNVtpXSApXG4gIH1cblxuICB0aGlzLmRlcml2KCB0aGlzLl9rNiwgdGhpcy5fdywgdGhpcy50ICsgZHQgKiAwLjg3NSApXG5cbiAgcmV0dXJuIHRoaXNcbn1cblxuSW50ZWdyYXRvci5wcm90b3R5cGUuX2NhbGN1bGF0ZUVycm9yID0gZnVuY3Rpb24oZHQpIHtcbiAgLy92YXIgY3MxID0gIDAuMTAyMTc3MzcyNjg1MTg1MTg1IC8vIDI4MjUvMjc2NDhcbiAgLy92YXIgY3MyID0gIDAuMDAwMDAwMDAwMDAwMDAwMDAwIC8vIDBcbiAgLy92YXIgY3MzID0gIDAuMzgzOTA3OTAzNDM5MTUzNDM5IC8vIDE4NTc1LzQ4Mzg0XG4gIC8vdmFyIGNzNCA9ICAwLjI0NDU5MjczNzI2ODUxODUxOCAvLyAxMzUyNS81NTI5NlxuICAvL3ZhciBjczUgPSAgMC4wMTkzMjE5ODY2MDcxNDI4NTcgLy8gMjc3LzE0MzM2XG4gIC8vdmFyIGNzNiA9ICAwLjI1MDAwMDAwMDAwMDAwMDAwMCAvLyAxLzRcblxuICAvL3ZhciBkYzEgPSAgMC4wMDQyOTM3NzQ4MDE1ODczMDEgLy8gY3MxIC0gYzFcbiAgLy92YXIgZGMyID0gIDAuMDAwMDAwMDAwMDAwMDAwMDAwIC8vIGNzMiAtIGMyXG4gIC8vdmFyIGRjMyA9IC0wLjAxODY2ODU4NjA5Mzg1NzgzMiAvLyBjczMgLSBjM1xuICAvL3ZhciBkYzQgPSAgMC4wMzQxNTUwMjY4MzA4MDgwODAgLy8gY3M0IC0gYzRcbiAgLy92YXIgZGM1ID0gIDAuMDE5MzIxOTg2NjA3MTQyODU3IC8vIGNzNSAtIGM1XG4gIC8vdmFyIGRjNiA9IC0wLjAzOTEwMjIwMjE0NTY4MDQwNiAvLyBjczYgLSBjNlxuXG4gIHZhciBlcnJvciA9IDBcbiAgZm9yKHZhciBpPTA7IGk8dGhpcy5uOyBpKyspIHtcbiAgICBlcnJvciA9ICB0aGlzLmVycm9yUmVkdWNlRnVuY3Rpb24oIGksIGVycm9yLFxuICAgICAgZHQgKiAoXG4gICAgICAgICAwLjAwNDI5Mzc3NDgwMTU4NzMwMSAqIHRoaXMuX2sxW2ldICtcbiAgICAgICAgLTAuMDE4NjY4NTg2MDkzODU3ODMyICogdGhpcy5fazNbaV0gK1xuICAgICAgICAgMC4wMzQxNTUwMjY4MzA4MDgwODAgKiB0aGlzLl9rNFtpXSArXG4gICAgICAgICAwLjAxOTMyMTk4NjYwNzE0Mjg1NyAqIHRoaXMuX2s1W2ldICtcbiAgICAgICAgLTAuMDM5MTAyMjAyMTQ1NjgwNDA2ICogdGhpcy5fazZbaV1cbiAgICAgICkgLyB0aGlzLl9lcnJvclNjYWxlW2ldXG4gICAgKVxuICB9XG5cbiAgcmV0dXJuIHRoaXMuZXJyb3JQb3N0RnVuY3Rpb24oZXJyb3IpXG59XG5cbkludGVncmF0b3IucHJvdG90eXBlLl91cGRhdGUgPSBmdW5jdGlvbihkdCkge1xuICAvL3ZhciBjMSAgPSAgMC4wOTc4ODM1OTc4ODM1OTc4ODMgLy8gMzcvMzc4XG4gIC8vdmFyIGMyICA9ICAwLjAwMDAwMDAwMDAwMDAwMDAwMCAvLyAwXG4gIC8vdmFyIGMzICA9ICAwLjQwMjU3NjQ4OTUzMzAxMTI3MiAvLyAyNTAvNjIxXG4gIC8vdmFyIGM0ICA9ICAwLjIxMDQzNzcxMDQzNzcxMDQzNyAvLyAxMjUvNTk0XG4gIC8vdmFyIGM1ICA9ICAwLjAwMDAwMDAwMDAwMDAwMDAwMCAvLyAwXG4gIC8vdmFyIGM2ICA9ICAwLjI4OTEwMjIwMjE0NTY4MDQwNiAvLyA1MTIvMTc3MVxuXG4gIGZvcih2YXIgaT0wOyBpPHRoaXMubjsgaSsrKSB7XG4gICAgdGhpcy55W2ldICs9IGR0ICogKFxuICAgICAgMC4wOTc4ODM1OTc4ODM1OTc4ODMgKiB0aGlzLl9rMVtpXSArXG4gICAgICAwLjQwMjU3NjQ4OTUzMzAxMTI3MiAqIHRoaXMuX2szW2ldICtcbiAgICAgIDAuMjEwNDM3NzEwNDM3NzEwNDM3ICogdGhpcy5fazRbaV0gK1xuICAgICAgMC4yODkxMDIyMDIxNDU2ODA0MDYgKiB0aGlzLl9rNltpXVxuICAgIClcbiAgfVxuICB0aGlzLnQgKz0gZHRcbiAgcmV0dXJuIHRoaXNcbn1cblxuSW50ZWdyYXRvci5wcm90b3R5cGUuX2NhbGN1bGF0ZUVycm9yU2NhbGUgPSBmdW5jdGlvbihkdCkge1xuICBmb3IodmFyIGk9MDsgaTx0aGlzLm47IGkrKykge1xuICAgIHRoaXMuX2Vycm9yU2NhbGVbaV0gPSB0aGlzLmVycm9yU2NhbGVGdW5jdGlvbihpLCBkdCwgdGhpcy55W2ldLCB0aGlzLl9rMVtpXSlcbiAgfVxuICByZXR1cm4gdGhpc1xufVxuXG5JbnRlZ3JhdG9yLnByb3RvdHlwZS5zdGVwID0gZnVuY3Rpb24oIHRMaW1pdCApIHtcbiAgLy8gQmFpbCBvdXQgZWFybHkgaWYgd2UncmUgKmF0KiB0aGUgbGltaXQ6XG4gIGlmIChNYXRoLmFicyh0aGlzLnQgLSB0TGltaXQpIDwgdGhpcy5kdCAqIDFlLTEwKSB7XG4gICAgcmV0dXJuIGZhbHNlO1xuICB9XG5cbiAgdmFyIHRoaXNEdCA9IHRoaXMuZHQ7XG5cbiAgLy8gRG9uJ3QgaW50ZWdyYXRlIHBhc3QgYSB0TGltaXQsIGlmIHByb3ZpZGVkOlxuICBpZiggdExpbWl0ICE9PSB1bmRlZmluZWQgKSB7XG4gICAgdGhpc0R0ID0gdGhpc0R0ID4gMCA/IE1hdGgubWluKCB0TGltaXQgLSB0aGlzLnQsIHRoaXNEdCApIDogTWF0aC5tYXgoIHRMaW1pdCAtIHRoaXMudCwgdGhpc0R0IClcbiAgfVxuXG4gIC8vIExpbWl0IHRoZSBtYWduaXR1ZGUgb2YgZHQgdG8gZHRNYXhNYWdcbiAgaWYoIHRoaXMuZHRNYXhNYWcgIT09IHVuZGVmaW5lZCAmJiBNYXRoLmFicyggdGhpc0R0ICkgPiB0aGlzLmR0TWF4TWFnICkge1xuICAgIHRoaXMuX19sb2coJ3N0ZXAnLCAnc3RlcCBncmVhdGVyIHRoYW4gbWF4aW11bSBzdGVwc2l6ZSByZXF1ZXN0ZWQuIGR0IG1hZ25pdHVkZSBoYXMgYmVlbiBsaW1pdGVkLicpXG4gICAgdGhpc0R0ID0gdGhpc0R0ID4gMCA/IHRoaXMuZHRNYXhNYWcgOiAtdGhpcy5kdE1heE1hZ1xuICB9XG5cbiAgLy8gTGltaXQgdGhlIG1hZ25pdHVkZSBvZiBkdCB0byBkdE1pbk1hZ1xuICBpZiggdGhpcy5kdE1pbk1hZyAhPT0gdW5kZWZpbmVkICYmIE1hdGguYWJzKCB0aGlzRHQgKSA8IHRoaXMuZHRNaW5NYWcgKSB7XG4gICAgdGhpcy5fX2xvZygnc3RlcCcsICdzdGVwIHNtYWxsZXIgdGhhbiBtaW5pbXVtIHN0ZXBzaXplIHJlcXVlc3RlZC4gZHQgbWFnbml0dWRlIGhhcyBiZWVuIGxpbWl0ZWQuJylcbiAgICB0aGlzRHQgPSB0aGlzRHQgPiAwID8gdGhpcy5kdE1pbk1hZyA6IC10aGlzLmR0TWluTWFnXG4gIH1cblxuICAvLyBUaGUgZmlyc3QgZGVyaXZhdGl2ZSBkb2Vzbid0IGNoYW5nZSBldmVuIGlmIGR0IGRvZXMsIHNvIG9ubHkgY2FsY3VsYXRlIHRoaXMgb25jZTpcbiAgdGhpcy5fY2FsY3VsYXRlSzEoKVxuXG4gIC8vIFRoZSBzY2FsZSBmYWN0b3IgcGVyLWRpbWVuc2lvbiBwcm9iYWJseSBkb2Vzbid0IG5lZWQgdG8gY2hhbmdlIGVpdGhlciBhY3Jvc3MgYSBzaW5nbGUgYWRhcHRpdmUgc3RlcDpcbiAgdGhpcy5fY2FsY3VsYXRlRXJyb3JTY2FsZSh0aGlzRHQpXG5cbiAgdmFyIGVycm9yID0gSW5maW5pdHlcbiAgdmFyIG1heEVycm9yID0gMFxuICB2YXIgbmV4dER0XG4gIHZhciBsb3dlckR0TGltaXRSZWFjaGVkID0gZmFsc2VcblxuICB3aGlsZSh0cnVlKSB7XG5cbiAgICAvLyBDYWxjdWxhdGUgaW50ZXJtZWRpYXRlIGsncyBmb3IgdGhlIHByb3Bvc2VkIHN0ZXA6XG4gICAgdGhpcy5fY2FsY3VsYXRlS3ModGhpc0R0KVxuXG4gICAgLy8gQ2FsY3VsYXRlIHRoZSBtYXggZXJyb3Igb2YgdGhlIHByb3Bvc2VkIHN0ZXA6XG4gICAgZXJyb3IgPSB0aGlzLl9jYWxjdWxhdGVFcnJvcih0aGlzRHQpXG5cbiAgICBpZiggZXJyb3IgPCB0aGlzLnRvbCB8fCBsb3dlckR0TGltaXRSZWFjaGVkICkge1xuICAgICAgLy8gU3VjY2VzcyEgRXhpdDpcbiAgICAgIGJyZWFrXG4gICAgfVxuXG4gICAgaWYoICEgTnVtYmVyLmlzRmluaXRlKGVycm9yKSApIHtcbiAgICAgIHRocm93IG5ldyBFcnJvcignb2RlNDUtY2FzaC1rYXJwOjpzdGVwKCkgTmFOIGVuY291bnRlcmVkIHdoaWxlIGludGVncmF0aW5nLicpXG4gICAgfVxuXG4gICAgLy8gRmFpbHVyZS4gQWRhcHQgdGhlIHRpbWVzdGVwOlxuICAgIG5leHREdCA9IHRoaXMuc2FmZXR5RmFjdG9yICogdGhpc0R0ICogTWF0aC5wb3coIHRoaXMudG9sIC8gZXJyb3IsIDAuMiApXG5cbiAgICAvLyBDdXQgdGhlIHRpbWVzdGVwLCBidXQgbm90IGJ5IG1vcmUgdGhhbiBtYXhEZWNyZWFzZUZhY3RvclxuICAgIHRoaXNEdCA9IG1heE1hZyggdGhpc0R0IC8gdGhpcy5tYXhEZWNyZWFzZUZhY3RvciwgbmV4dER0IClcblxuICAgIC8vIElmIHN0ZXBzaXplIHRvbyBzbWFsbCwgZmluaXNoIG9mZiBieSB0YWtpbmcgdGhlIGN1cnJlbnRseSBwcm9wb3NlZCBzdGVwIGFuZCBsb2dnaW5nIGEgd2FybmluZzpcbiAgICBpZiggdGhpcy5kdE1pbk1hZyAhPT0gdW5kZWZpbmVkICYmIE1hdGguYWJzKHRoaXNEdCkgPCB0aGlzLmR0TWluTWFnICkge1xuICAgICAgdGhpc0R0ID0gdGhpcy5kdE1pbk1hZyAqICh0aGlzRHQgPiAwID8gMSA6IC0xKTtcbiAgICAgIHRoaXMuX19sb2coJ3N0ZXAnLCAnbWluaW11bSBzdGVwc2l6ZSByZWFjaGVkLicpXG4gICAgICBsb3dlckR0TGltaXRSZWFjaGVkID0gdHJ1ZVxuICAgIH1cbiAgfVxuXG4gIC8vIEFwcGx5IHRoaXMgdXBkYXRlOlxuICB0aGlzLl91cGRhdGUodGhpc0R0KVxuXG4gIC8vIENhbGN1bGF0ZSB0aGUgbmV4dCB0aW1lc3RlcCBzaXplOlxuICBuZXh0RHQgPSB0aGlzLnNhZmV0eUZhY3RvciAqIHRoaXNEdCAqIE1hdGgucG93KCB0aGlzLnRvbCAvIGVycm9yLCAwLjI1IClcblxuICAvLyBJbmNyZWFzZSB0aGUgdGltZXN0ZXAgZm9yIHRoZSBuZXh0IHRpbWUgYXJvdW5kLCBidXQgbm90IGJ5IG1vcmUgdGhhbiB0aGUgbWF4SW5jcmVhc2VGYWN0b3I6XG4gIHRoaXMuZHQgPSBtYXhNYWcodGhpcy5kdCAvIHRoaXMubWF4RGVjcmVhc2VGYWN0b3IsIG1pbk1hZyggdGhpcy5kdCAqIHRoaXMubWF4SW5jcmVhc2VGYWN0b3IsIG5leHREdCApKTtcblxuICBpZiggdExpbWl0ICE9PSB1bmRlZmluZWQgKSB7XG4gICAgcmV0dXJuIE1hdGguYWJzKHRoaXMudCAtIHRMaW1pdCkgPiB0aGlzLmR0ICogMWUtODtcbiAgfSBlbHNlIHtcbiAgICByZXR1cm4gdHJ1ZVxuICB9XG59XG5cbkludGVncmF0b3IucHJvdG90eXBlLnN0ZXBzID0gZnVuY3Rpb24oIG4sIHRMaW1pdCApIHtcbiAgZm9yKHZhciBzdGVwPTA7IHN0ZXA8bjsgc3RlcCsrKSB7XG4gICAgaWYoICEgdGhpcy5zdGVwKHRMaW1pdCkgKSByZXR1cm4gZmFsc2U7XG4gIH1cbn1cblxuZnVuY3Rpb24gSW50ZWdyYXRvckZhY3RvcnkoIHkwLCBkZXJpdiwgdCwgZHQsIG9wdGlvbnMgKSB7XG4gIHJldHVybiBuZXcgSW50ZWdyYXRvciggeTAsIGRlcml2LCB0LCBkdCwgb3B0aW9ucyApXG59XG4iLCJjb25zdCBzb2x2ZXIgPSByZXF1aXJlKCdvZGU0NS1jYXNoLWthcnAnKTtcclxuIl19
