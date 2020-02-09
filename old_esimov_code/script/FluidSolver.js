/**
 * @created by Simo Endre
 * @simo_endre
 */

FS.Solver = (function() {
    // Private variables
    var _NX, _NY, _NX2, _NY2,
        _invNumCells, _dt, _isRGB,
        _solverIterations, _colorDiffusion, _doVorticityConfinement;

    var wrap_x = false;
    var wrap_y = false;

    var _visc, _fadeSpeed, _avgDensity, _uniformity, _avgSpeed;

    var temp = 0;

    // Constants
    var FLUID_DEFAULT_NX = 50,
        FLUID_DEFAULT_NY = 50,
        FLUID_DEFAULT_DT = 0.1,
        FLUID_DEFAULT_VISC = 0.00012,
        FLUID_DEFAULT_COLOR_DIFUSION = 0.1,
        FLUID_DEFAULT_FADESPEED = 0.0014,
        FLUID_DEFAULT_SOLVER_ITERATION = 4,
        FLUID_DEFAULT_VORTICITY_CONFINEMENT = false;
    var self;

    function FSolver() {

        // Public variales
        this.width = 0;
        this.height = 0;
        this.numCells = 0;        

        self = this;

        this.__defineGetter__('wrapX', function() { return wrap_x; });
        this.__defineGetter__('wrapY', function() { return wrap_y; });
		this.__defineGetter__('fadeSpeed', function() { return _fadeSpeed; });
		this.__defineGetter__('solverIterations', function() { return _solverIterations; });
        this.__defineSetter__('rgb', function(value) { _isRGB = value; });
        this.__defineSetter__('fadeSpeed', function(value) { _fadeSpeed = value; });
		this.__defineSetter__('solverIterations', function(value) { _solverIterations = value; });
        this.__defineSetter__('viscosity', function(value) { _visc = value; });
        this.__defineSetter__('deltaT', function(value) { _dt = value; });
        this.__defineSetter__('vorticityConfinement', function(value) { _doVorticityConfinement = value; });
    }

    FSolver.prototype.mainSolver = function(NX, NY) {
        setup(NX, NY);

        this.reset();
    };    

    FSolver.prototype.reset = function() {
        if (typeof Float32Array == null || typeof Float32Array == "undefined")
            Float32Array = Array;        
        this.r = new Float32Array(this.numCells);
        this.g = new Float32Array(this.numCells);
        this.b = new Float32Array(this.numCells);
        this.u = new Float32Array(this.numCells);
        this.v = new Float32Array(this.numCells);                      

        this.rOld = new Float32Array(this.numCells);
        this.gOld = new Float32Array(this.numCells);
        this.bOld = new Float32Array(this.numCells);
        this.uOld = new Float32Array(this.numCells);
        this.vOld = new Float32Array(this.numCells);        

        this.curl_abs = new Float32Array(this.numCells);
        this.curl_orig = new Float32Array(this.numCells);
        this.density = new Float32Array(this.numCells);   
        this.densityOld = new Float32Array(this.numCells);

        this.source = new Float32Array(this.numCells);
        temp = new Float32Array(this.numCells);

        var i = this.numCells;

        while(i-- > -1) {
            this.r[i] = this.rOld[i] = this.g[i] = this.gOld[i] = this.b[i] = this.bOld[i] = 0;
            this.u[i] = this.uOld[i] = this.v[i] = this.vOld[i] = 0;
            this.curl_abs[i] = this.curl_orig[i] = 0;
            this.density[i] = this.densityOld[i] = 0;
            this.source[i] = 0;
        }
    };

    FSolver.prototype.update = function() {
        addDensitySource(this.density, this.densityOld);
        diffuse(0, this.densityOld, this.density, 0);
        advect(0, this.density, this.densityOld, this.u, this.v);

        if (_doVorticityConfinement) {
            calcVorticityConfinement(this.uOld, this.vOld);            
        }

        this.addCellVelocity();
        this.SWAP('uOld', 'u');
        diffuse(1, this.u, this.uOld, _visc);
        this.SWAP('vOld', 'v');
        diffuse(2, this.v, this.vOld, _visc);        
        project(this.u, this.v, this.uOld, this.vOld);


        this.SWAP('uOld', 'u');
        this.SWAP('vOld', 'v');
        advect(1, this.u, this.uOld, this.uOld, this.vOld);
        advect(2, this.v, this.vOld, this.uOld, this.vOld);
        project(this.u, this.v, this.uOld, this.vOld);

        if (_isRGB) {
            addSourceRGB();
            this.SWAP('rOld', 'r');
            
            if (_colorDiffusion != 0 && _dt != 0) {
               diffuseRGB(_colorDiffusion);
               this.SWAP('rOld', 'r');
            }
            
            advectRGB(0, this.u, this.v);
            fadeRGB();
        } else {
            addDensitySource(this.density, this.densityOld);
            swapR();

            diffuse(0, this.r, this.rOld, 0);
            this.SWAP('rOld', 'r');
            advect(0, this.r, this.rOld, this.u, this.v);
            fadeR();
        }
    };


    FSolver.prototype.updateDensity = function() {        
        addDensitySource(this.density, this.densityOld);

        if (_doVorticityConfinement) {            
            calcVorticityConfinement(this.uOld, this.vOld);            
        }
        
        diffuse(0, this.densityOld, this.density, 0);
        advect(0, this.density, this.densityOld, this.u, this.v);
    }
	
    FSolver.prototype.updateVelocity = function() {
        this.addCellVelocity();
        this.SWAP('uOld', 'u');
        diffuse(1, this.u, this.uOld, _visc);
        this.SWAP('vOld', 'v');        
        diffuse(2, this.v, this.vOld, _visc);        
        project(this.u, this.v, this.uOld, this.vOld);


        this.SWAP('uOld', 'u');
        this.SWAP('vOld', 'v');
        advect(1, this.u, this.uOld, this.uOld, this.vOld);
        advect(2, this.v, this.vOld, this.uOld, this.vOld);
        project(this.u, this.v, this.uOld, this.vOld);
    }

    FSolver.prototype.getIndexForCellPosition = function(i, j) {
        i = (i < 1) ? 1 : ((i > _NX) ? _NX : i);
        j = (j < 1) ? 1 : ((j > _NY) ? _NY : i);        
        return FLUID_IX(i, j);
    };

    FSolver.prototype.getIndexForNormalizedPosition = function(x, y) {
        return this.getIndexForCellPosition(parseInt(x * _NX2), parseInt(y * _NY2));
    };

    FSolver.prototype.setWrap = function(x, y) {
        x = x || false;
        y = y || false;
        wrap_x = x;
        wrap_y = y;
    };

    FSolver.prototype.SWAP = function(x0, x) {
        var tmp = this[x0];
        this[x0] = this[x];
        this[x] = tmp;
    }

    FSolver.prototype.addCellVelocity = function() {
        var size = self.numCells;
        while(size-- > -1) {
            if (isNaN(self.u[size])) continue;
            self.u[size] += _dt * self.uOld[size];
            if (isNaN(self.v[size])) continue;
            self.v[size] += _dt * self.vOld[size];
        }
    }
	
	FSolver.prototype.getDensity = function(x, y) {
		return this.density[(x + 1) + (y + 1) * _NY];
	};


    // Private functions

    function setup(NX, NY) {
        _dt = FLUID_DEFAULT_DT;
        _visc = FLUID_DEFAULT_VISC;
        _fadeSpeed = FLUID_DEFAULT_FADESPEED;
        _solverIterations = this.solverIterations;
        _colorDiffusion = FLUID_DEFAULT_COLOR_DIFUSION;
        _doVorticityConfinement = FLUID_DEFAULT_VORTICITY_CONFINEMENT;

        _NX = NX;
        _NY = NY;
        _NX2 = _NX + 2; // cells + extra boundary cells
        _NY2 = _NY + 2;

        self.numCells = _NX2 * _NY2;
        _invNumCells = 1.0 / self.numCells;

        self.width = _NX2;
        self.height = _NY2;

        _isRGB = false;
    }

    function addDensitySource(x, x0) {                
        var size = self.numCells;        
        while(size-- > -1) {
            x[size] += _dt * x0[size];
        }
    }    

    function addSourceRGB() {
        var size = self.numCells;
        while(size-- > -1) {
            if (isNaN(self.r[size])) continue;
            self.r[size] += _dt * self.rOld[size];
            if (isNaN(self.g[size])) continue;
            self.g[size] += _dt * self.gOld[size];
            if (isNaN(self.b[size])) continue;
            self.b[size] += _dt * self.bOld[size];            
        }
    }

    function diffuse(bound, c, c0, _diff) {
        var a = _dt * _diff * _NX * _NY;
        linearSolver(bound, c, c0, a, 1.0 + 4 * a);
    }

    function diffuseUV(bound, _diff) {
        var a = _dt * _diff * _NX * _NY;
        linearSolverUV(bound, a, 1.0 + 4 * a);
    }

    function diffuseRGB(_diff) {
        var a = _dt * _diff * _NX * _NY;
        linearSolverRGB( a, 1.0 + 4 * a);
    }

    function calcVorticityConfinement(x, y) {
        var i, j, index, dx, dy, av, length;

        for (j = _NY; j > 0; --j) {
            index = FLUID_IX(_NX, j);
            for (i = _NX; i > 0; --i) {
                dx = self.u[parseInt(index + _NX2)] - self.u[parseInt(index - _NX2)];
                dy = self.v[parseInt(index + 1)] - self.v[parseInt(index - 1)];

                av = (dy - dx) * 0.5;
                self.curl_orig[parseInt(index)] = av;
                self.curl_abs[parseInt(index)] = (av < 0) ? -av : av;

                --index;
            }
        }

        for (j = _NY-1; j > 1; --j) {
            index = FLUID_IX(_NX-1, j);
            for (i = _NX-1; i > 1; --i) {
                dx = self.curl_abs[parseInt(index + 1)] - self.curl_abs[parseInt(index - 1)];
                dy = self.curl_abs[parseInt(index + _NX2)] - self.curl_abs[parseInt(index - _NX2)];

                length = Math.sqrt(dy * dy + dx * dx) + 0.0001;
                length = 2 / length;                
                dx *= length;
                dy *= length;

                av = self.curl_orig[parseInt(index)];                
                y[parseInt(index)] = dx * av;
                x[parseInt(index)] = dy * -av;

                --index;                
            }
        }
    }

    function fadeR() {        
        var holdAmount = 1 - _fadeSpeed;
        var totalDeviations = 0;
        var currentDeviation = 0;
        var i;

        _avgDensity = 0;
        _avgSpeed = 0;

        for (i = 0; i < self.numCells; i++) {
            self.uOld[i] = 0; self.vOld[i] = 0;
            self.rOld[i] = 0;

            _avgSpeed = self.u[i] * self.u[i] + self.v[i] * self.v[i];
            var density = Math.min(1.0, self.r[i]);
            _avgDensity += density;

            // calc deviation for uniformity
            currentDeviation = density - _avgDensity;
            totalDeviations += currentDeviation * currentDeviation;
            // fade out old
            self.r[i] = density * holdAmount;
        }
        _avgDensity *= _invNumCells;

        _uniformity = 1.0 / (1 + totalDeviations * _invNumCells); // 0: very wide distribution, 1: very uniform
    }

    function fadeRGB() {
        var holdAmount = 1 - _fadeSpeed;
        var totalDeviations = 0;
        var currentDeviation = 0;

        _avgDensity = 0;
        _avgSpeed = 0;

        for (var i = 0; i < self.numCells; i++) {
            self.uOld[i] = 0; self.vOld[i] = 0;
            self.rOld[i] = 0; self.gOld[i] = 0; self.bOld[i] = 0;

            _avgSpeed = self.u[i] * self.u[i] + self.v[i] * self.v[i];
            var dR = Math.min(1.0, self.r[i]);
            var dG = Math.min(1.0, self.g[i]);
            var dB = Math.min(1.0, self.b[i]);
            var density = Math.max(dR, Math.max(dG, dB));
            _avgDensity += density;
            currentDeviation = density - _avgDensity;

            // calc deviation for uniformity
            totalDeviations += currentDeviation * currentDeviation;
            // fade out old
            self.r[i] = dR * holdAmount;
            self.g[i] = dG * holdAmount;
            self.b[i] = dB * holdAmount;
        }

        _avgDensity *= _invNumCells;
        _avgSpeed *= _invNumCells;
        _uniformity = 1.0 / (1+ totalDeviations * _invNumCells); // 0: very wide distribution, 1: very uniform
    }

    function advect(bound, _d, d0, du, dv) {
        var i, j, i0, j0, i1, j1;
        var x, y, s0, t0, s1, t1, dt0, dt1;

        dt0 = _dt * _NX;
        dt1 = _dt * _NY;

        for (j = _NY; j > 0; --j) {
            for (i = _NX; i > 0; --i) {
                x = i - dt0 * du[FLUID_IX(i, j)];
                y = j - dt1 * dv[FLUID_IX(i, j)];

                if (x > _NX + 0.5) x = _NX + 0.5;
                if (x < 0.5) x = 0.5;

                i0 = parseInt(~~x);
                i1 = i0 + 1;

                if (y > _NY + 0.5) y = _NY + 0.5;
                if (y < 0.5) y = 0.5;

                j0 = parseInt(~~y);
                j1 = j0 + 1;

                s1 = x - i0;
                s0 = 1 - s1;
                t1 = y - j0;
                t0 = 1 - t1;

                _d[FLUID_IX(i, j)] = s0 * (t0 * d0[FLUID_IX(i0, j0)] + t1 * d0[FLUID_IX(i0, j1)]) +
                    s1 * (t0 * d0[FLUID_IX(i1, j0)] + t1 * d0[FLUID_IX(i1, j1)]);
            }
        }

        setBoundary(bound, _d);
    }

    function advectRGB(du, dv) {
        var i, j, i0, j0, i1, j1;
        var x, y, s0, t0, s1, t1, dt0x, dt0y;
        var index;

        dt0x = _dt * _NX;
        dt0y = _dt * _NY;

        for (j = _NY; j > 0; --j) {
            for (i = _NX; i > 0; --i) {               
                index = FLUID_IX[i, j];
                x = i - dt0x * du[FLUID_IX(i, j)];
                y = j - dt0y * dv[FLUID_IX(i, j)];

                if (x > _NX + 0.5) x = _NX + 0.5;
                if (x < 0.5) x = 0.5;

                i0 = parseInt(x);
                i1 = i0 + 1;

                if (y > _NY + 0.5) y = _NY + 0.5;
                if (y < 0.5) y = 0.5;

                j0 = parseInt(y);
                j1 = j0 + 1;

                s1 = x - i0;
                s0 = 1 - s1;
                t1 = y - j0;
                t0 = 1 - t1;
                j0 = i0 + _NX2;
                
                io = FLUID_IX(i0, j0);
                
                self.r[index] = s0 * (t0 * self.rOld[i0] + t1 * self.rOld[j0]) + s1 * (t0 * self.rOld[i1] + t1 * self.rOld[j0+1]);
                self.g[index] = s0 * (t0 * self.gOld[i0] + t1 * self.gOld[j0]) + s1 * (t0 * self.gOld[i1] + t1 * self.gOld[j0+1]);
                self.b[index] = s0 * (t0 * self.bOld[i0] + t1 * self.bOld[j0]) + s1 * (t0 * self.bOld[i1] + t1 * self.bOld[j0+1]);
            }
        }

        setBoundaryRGB();
    }

    function project(u, v, p, div) {
        var i, j, index;
        var h = - 0.5 / _NX;
        for (j = _NY; j > 0; --j) {
            index = FLUID_IX(_NX, j);
            for (i = _NX; i > 0; --i) {
                div[index] = h * (u[parseInt(index+1)] - u[parseInt(index-1)] + v[parseInt(index + _NX2)] - v[parseInt(index - _NX2)] );
                p[index] = 0;
                --index;
            }
        }

        setBoundary(0, div);
        setBoundary(0, p);

        linearSolver(0, p, div, 1, 4);

        for (j = _NY; j > 0; --j) {
            index = FLUID_IX(_NX, j);
            for (i = _NX; i > 0; --i) {
                u[index] -= 0.5 * _NX * (p[parseInt(index+1)] - p[parseInt(index-1)]);
                v[index] -= 0.5 * _NY * (p[parseInt(index + _NX2)] - p[parseInt(index - _NX2)]);
                --index;
            }
        }
        setBoundary(1, u);
        setBoundary(2, v);
    }

    function linearSolver(bound, x, x0, a, c) {
        var i, j, k, index;

        if (a == 1 && c == 4) {            
            for (k = 0; k < _solverIterations; k++) {
                for (j = _NY; j > 0; --j) {
                    index = FLUID_IX(_NX, j);
                    for (i = _NX; i > 0; --i) {
                        x[index] = ( (x[parseInt(index-1)] + x[parseInt(index+1)] +
                            x[parseInt(index - _NX2)] + x[parseInt(index + _NX2)]) + x0[parseInt(index)]) * .25;
                        --index;
                    }
                }

                setBoundary(bound, x);
            }
        } else {
            c = 1 / c;
            for (k = 0; k < _solverIterations; k++) {
                for (j = _NY; j > 0; --j) {
                    index = FLUID_IX(_NX, j);
                    for (i = _NX; i > 0; --i) {
                        x[index] = (a * (x[parseInt(index-1)] + x[parseInt(index+1)] +
                            x[parseInt(index - _NX2)] + x[parseInt(index + _NX2)]) + x0[parseInt(index)]) * c;
                        --index;
                    }
                }

                setBoundary(bound, x);
            }
        }
    }

    function linearSolverRGB(a, c) {
        var i, j, k, index, index2, index3;
        c = 1 / c;

        for (k = 0; k < _solverIterations; ++k) {
            for (j = _NY; j > 0; --j) {
                index = FLUID_IX(_NX, j);
                index2 = index - _NX2;
                index3 = index + _NX2;
                
                for (i = _NX; i > 0; --i) {
                    self.r[index] = ( (self.r[parseInt(index-1)] + self.r[parseInt(index+1)] +
                        self.r[parseInt(index2)] + self.r[parseInt(index3)]) * a + self.rOld[index] ) * c;

                    self.g[index] = ( (self.g[parseInt(index-1)] + self.g[parseInt(index+1)] +
                        self.g[parseInt(index2)] + self.g[parseInt(index3)]) * a + self.gOld[index] ) * c;

                    self.b[index] = ( (self.b[parseInt(index-1)] + self.b[parseInt(index+1)] +
                        self.b[parseInt(index2)] + self.b[parseInt(index3)]) * a + self.bOld[index] ) * c;

                    --index;
                    --index2;
                    --index3;
                }
            }

            setBoundaryRGB();
        }
    }

    function linearSolverUV(a, c) {
        var i, j, k, index;
        c = 1 / c;

        for (k = 0; k < _solverIterations; ++k) {
            for (j = _NY; j > 0; --j) {
                index = FLUID_IX(_NX, j);                
                for (i = _NX; i > 0; --i) {
                    self.u[index] = ( (self.u[parseInt(index-1)] + self.u[parseInt(index+1)] + 
                              self.u[parseInt(index - _NX2)] + self.u[parseInt(index + _NX2)]) * a + self.uOld[index] ) * c;
                    self.v[index] = ( (self.v[parseInt(index-1)] + self.v[parseInt(index+1)] + 
                              self.v[parseInt(index - _NX2)] + self.v[parseInt(index + _NX2)]) * a + self.vOld[index] ) * c;
                    --index;
                }
            }

            setBoundary(1, self.u);
            setBoundary(2, self.v);
        }
    }

    function setBoundary(bound, x) {
        var dst1, dst2, src1, src2, i;
        var step = FLUID_IX(0, 1) - FLUID_IX(0, 0);
        dst1 = FLUID_IX(0, 1);
        src1 = FLUID_IX(1, 1);
        dst2 = FLUID_IX(_NX+1, 1);
        src2 = FLUID_IX(_NX, 1);

        if (wrap_x) {
            src1 ^= src2;
            src2 ^= src1;
            src1 ^= src2;
        }

        if (bound == 1 && !wrap_x) {
            for (i = _NY; i > 0; --i) {
                x[dst1] = -x[src1]; dst1 += step; src1 += step;
                x[dst2] = -x[src2]; dst2 += step; src2 += step;
            }
        } else {
            for (i = _NY; i > 0; --i) {
                x[dst1] = x[src1]; dst1 += step; src1 += step;
                x[dst2] = x[src2]; dst2 += step; src2 += step;
            }
        }

        dst1 = FLUID_IX(1, 0);
        src1 = FLUID_IX(1, 1);
        dst2 = FLUID_IX(1, _NY+1);
        src2 = FLUID_IX(1, _NY);

        if (wrap_y) {
            src1 ^= src2;
            src2 ^= src1;
            src1 ^= src2;
        }

        if (bound == 2 && !wrap_y) {
            for (i = _NX; i > 0 ; --i) {
                x[dst1++] = -x[src1++];
                x[dst2++] = -x[src2++];
            }
        } else {
            for (i = _NX; i > 0; --i) {
                x[dst1++] = x[src1++];
                x[dst2++] = x[src2++];
            }
        }

        x[FLUID_IX(0, 0)] = .5 * (x[FLUID_IX(1, 0)] + x[FLUID_IX(0, 1)]);
        x[FLUID_IX(0, _NY+1)] = .5 * (x[FLUID_IX(1, _NY+1)] + x[FLUID_IX(0, _NY)]);
        x[FLUID_IX(_NX+1, 0)] = .5 * (x[FLUID_IX(_NX, 0)] + x[FLUID_IX(_NX+1, _NX)]);
        x[FLUID_IX(_NX+1, _NY+1)] = .5 * (x[FLUID_IX(_NX, _NY+1)] + x[FLUID_IX(_NX+1, _NY)]);
    }

    function setBoundaryRGB() {
       
        if (!wrap_x && !wrap_y) return;
       
        var i, src1, src2, dst1, dst2;
        var step = FLUID_IX(0, 1) - FLUID_IX(0, 0);

        if (wrap_x) {
            dst1 = FLUID_IX(0, 1);
            src1 = FLUID_IX(1, 1);
            dst2 = FLUID_IX(_NX+1, 1);
            src2 = FLUID_IX(_NX, 1);
            
            src1 ^= src2;
            src2 ^= src1;
            src1 ^= src2;        
         
            for (i = _NY; i > 0 ; --i) {
                self.r[dst1] = self.r[src1]; self.g[dst1] = self.g[src1]; self.b[dst1] = self.b[src1]; 
                dst1 += step; src1 += step;

                self.r[dst2] = self.r[src2]; self.g[dst2] = self.g[src2]; self.b[dst2] = self.b[src2]; 
                dst2 += step; src2 += step;            
            }
        }

        if (wrap_y) {
            dst1 = FLUID_IX(1, 0);
            src1 = FLUID_IX(1, 1);
            dst2 = FLUID_IX(1, _NY+1);
            src2 = FLUID_IX(1, _NY);
            
            src1 ^= src2;
            src2 ^= src1;
            src1 ^= src2;        
         
            for (i = _NX; i > 0 ; --i) {
                self.r[dst1] = self.r[src1]; self.g[dst1] = self.g[src1]; self.b[dst1] = self.b[src1]; 
                ++dst1; ++src1;

                self.r[dst2] = self.r[src2]; self.g[dst2] = self.g[src2]; self.b[dst2] = self.b[src2]; 
                ++dst2; ++src2;
            }
        }
    }        

    function swapRGB() {
        temp = self.r; self.r = self.rOld; self.rOld = temp;
        temp = self.g; self.g = self.gOld; self.gOld = temp;
        temp = self.b; self.b = self.bOld; self.bOld = temp;
    }

    function FLUID_IX(i, j) {        
        return (parseInt(i + _NX2 * j)); //get the index of a two dimensional array
    }

    return FSolver;
})();