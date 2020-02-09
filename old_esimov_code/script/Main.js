(function() {
    var lastTime = 0;
    var vendors = ['ms', 'moz', 'webkit', 'o'];
    for(var x = 0; x < vendors.length && !window.requestAnimationFrame; ++x) {
        window.requestAnimationFrame = window[vendors[x]+'RequestAnimationFrame'];
        window.cancelAnimationFrame = window[vendors[x]+'CancelAnimationFrame']
            || window[vendors[x]+'CancelRequestAnimationFrame'];
    }

    if (!window.requestAnimationFrame)
        window.requestAnimationFrame = function(callback, element) {
            var currTime = new Date().getTime();
            var timeToCall = Math.max(0, 16 - (currTime - lastTime));
            var id = window.setTimeout(function() { callback(currTime + timeToCall); },
                timeToCall);
            lastTime = currTime + timeToCall;
            return id;
        };

    if (!window.cancelAnimationFrame)
        window.cancelAnimationFrame = function(id) {
            clearTimeout(id);
        };
}());

// Private Constants
var FLUID_WIDTH = 128;
var NUM_PARTICLES = 4;
var VMAX = 0.013;
var VMAX2 = VMAX * VMAX;
var ADD_DENSITY = 60000;
var VEL_MUL = 40;

// Private variables
var rgb, c, a, m, saturation, v2 = 0, vxNorm = 0, vyNorm = 0;
var mainCanvas = document.createElement('canvas');
mainCanvas.setAttribute('width', '900');
mainCanvas.setAttribute('height', '500');
mainCanvas.setAttribute('id', 'canvas');
var mainContext = document.getElementsByClassName('context')[0];
mainContext.appendChild(mainCanvas);  

var canvasWidth = mainCanvas.width;
var canvasHeight = mainCanvas.height;

var sx = mainCanvas.width / mainCanvas.offsetWidth;
var sy = mainCanvas.height / mainCanvas.offsetHeight;
var cx = 0, cy = 0, lastX, lastY;

mainCanvas.width = FLUID_WIDTH;
mainCanvas.height = FLUID_WIDTH;
var ctx = mainCanvas.getContext('2d');
var prevMouse = new Point();

var sw = mainCanvas.width;
var sh = mainCanvas.height;
var isw = 1 / sw;
var ish = 1 / sh;
var mx = 4.0;
var my = 0;
var aspectRatio = sw * ish;
var dx, dy;

// Initialize GUI values
var GUI = new function() {
    this.viscosity = 0.00006;  
    this.fadeSpeed = 0.0012; 
    this.densityMul = 0.98;
    this.density = ADD_DENSITY;
    this.velMul = VEL_MUL;  
    this.fluidWidth = FLUID_WIDTH;
    this.numParticles = NUM_PARTICLES; 
    this.solverIterations = 4;
    this.doVorticity = false;
    
  return this;
}  

var partList = new Array();
var fSolver = new FS.Solver();
    fSolver.mainSolver(GUI.fluidWidth, parseInt(GUI.fluidWidth * sh / sw));

    fSolver.viscosity = GUI.viscosity;
    fSolver.fadeSpeed = GUI.fadeSpeed;
    fSolver.solverIterations = GUI.solverIterations;
    fSolver.vorticityConfinement = GUI.doVorticity;

var fsWidth = fSolver.width;
var fsHeight = fSolver.height;
var image = ctx.createImageData(fsWidth, fsHeight);
var data32 = new Uint32Array(image.data.buffer);

var frameCount = 0;
var drawMode = 0,
    drawParticles = true,
    drawFluid = false,
    drawLines = false,
    mouseDown = false,    
    addVel    = true,
    add = true;

var blurMatrix = [
    0.01, 0.02, 0.04, 0.02, 0.01,
    0.02, 0.04, 0.08, 0.04, 0.02,
    0.04, 0.08, 0.16, 0.08, 0.04,
    0.02, 0.04, 0.08, 0.04, 0.02,
    0.01, 0.02, 0.04, 0.02, 0.01
];

var counter = 0;

var gui = new dat.GUI({autoPlace:false, width: 270, align: "left"});
gui.domElement.style.position = 'relative';
gui.domElement.style.top = '60px';
gui.domElement.style.left = '80px';

var datGUI = document.getElementById('gui');
datGUI.appendChild(gui.domElement);
var gs = gui.addFolder('General settings');
gs.open();
gs.add(GUI,'viscosity').min(0.00001).max(0.001).step(0.00001).name('Viscosity').onChange(function(value) { fSolver.viscosity = value });
gs.add(GUI,'fadeSpeed').min(0.0006).max(0.009).step(0.0001).name('Fade speed').onChange(function(value) { fSolver.fadeSpeed = value });
gs.add(GUI,'density').min(10000).max(90000).step(100).name('Density').onChange(function() { fSolver.updateDensity() });
gs.add(GUI,'velMul').min(10).max(100).step(1).name('Velocity').onChange(function() { fSolver.addCellVelocity() });
gs.add(GUI,'densityMul').min(0.7).max(1.1).step(0.1).name('Density Mult.').onChange(function(value) { GUI.densityMul = value });
gs.add(GUI,'solverIterations').min(3).max(10).step(1).name('Solver Iteration').onChange(function(value) { GUI.solverIterations = fSolver.solverIterations = value;} );
gs.add(GUI,'doVorticity').onFinishChange(function(e) { 
    if (GUI.doVorticity) { 
        fSolver.vorticityConfinement = true;        
    } else {
        fSolver.vorticityConfinement = false;
    }    
});

//gs.add(GUI,'fluidWidth').min(32).max(512).step(64).name('Fluid width').onChange(function(value) { GUI.fluidWidth = value });
gui.add(fSolver, 'reset').name('Clear');

function Particle() {

    // Public constants
    this.MOMENTUM = 0.5;
    this.FLUID_FORCE = 0.5;

    // Public variables
    this.x = 0;
    this.y = 0;
    this.vx = 0;
    this.vy = 0;
    this.radius = 0;
    this.alpha = 0;
    this.mass = 0;
}

Particle.prototype.init = function(x, y) {
    this.x = x;
    this.y = y;

    this.vx = this.vy = 0;
    this.radius = 5;
    this.alpha = Utils.Random.fl(.3, 1);
    this.mass = Utils.Random.fl(.1, 1);
};

Particle.prototype.update = function() {
    if (this.alpha == 0) return;
    var fluidIndex = fSolver.getIndexForNormalizedPosition(this.x * isw, this.y * ish);
    this.vx = fSolver.u[fluidIndex] * sw * this.mass * this.FLUID_FORCE + this.vx * this.MOMENTUM;
    this.vy = fSolver.v[fluidIndex] * sh * this.mass * this.FLUID_FORCE + this.vy * this.MOMENTUM;

    this.x += this.vx;
    this.y += this.vy;

    if (this.x < 0) {
        if (fSolver.wrapX)
            this.x += sw;
        else {
            this.x = 1;
            this.vx *= -1; // inverse velocity
        }
    }

    if (this.y < 0) {
        if (fSolver.wrapY)
            this.y += sh;
        else {
            this.y = 1;
            this.vy *= -1; // inverse velocity
        }
    }

    if (this.x > this.sw) {
        if (fSolver.wrapX)
            this.x -= sw;
        else {
            this.x = 1- sw;
            this.vx *= -1;
        }
    }

    if (this.y > this.sh) {
        if (fSolver.wrapY)
            this.y -= this.y;
        else {
            this.y = 1 - sh;
            this.vy *= -1;
        }
    }

    // makes particles to glitter if they slows down below threshold
    if (this.vx * this.vx + this.vy * this.vy < .5) {
        this.alpha = 0;
        this.vx = Utils.Random.fl (-1, 1);
        this.vy = Utils.Random.fl (-1, 1);
    }

    this.alpha *= 0.0999;
    if (this.alpha < 0.02) { this.alpha = 0; }
};


function addParticles(x, y, n) {
    while (n--) {
        var p = new Particle();
        p.init(x + Utils.Random.fl(-15.0, 15.0), y + Utils.Random.fl(-15.0, 15.0));
        partList.push(p);        
    }
}


function particlesUpdate(bmp) {
    var i = 0;
    var length = partList.length;
    var rgb = null;

    while (i++ < length-1) {       
        var p = partList[i];  
        if (p instanceof Particle) {
            if (p.alpha > 0) {
                p.update();
                a = parseInt(p.alpha * 0xFF + 0.5);

                if (drawFluid) {
                    c = a << 24 | a << 16 | a << 8 | a;                    
                } else {                   
                    vxNorm = p.vx * isw;
                    vyNorm = p.vy * ish;
                    v2 = vxNorm * vxNorm + vyNorm * vyNorm;
                    if (v2 > VMAX) v2 = VMAX;
                    
                    m = p.mass;
                    saturation = (p.mass > 0.5) ? p.mass * p.mass * p.mass : 0;
                    saturation *= saturation;                                        
                    rgb = hsvToRgb(0, Utils.Number.map2(v2, 0, VMAX2, 0, 1) + saturation,
                        Utils.Number.interpolate(m, 0.5, 1) * 1);                        
                        
                    //console.log(v2 + ":" + Utils.Number.map2(v2, 0, VMAX2, 0, 1) + " : " + Utils.Number.interpolate(m, 0.5, 1) * p.alpha);
                    bmp.context.globalAlpha = a;
                    c = (rgb.r & 0xff) << 16 | (rgb.g & 0xff) << 8 | (rgb.b & 0xff);
                }
                line(bmp, parseInt(p.x - p.vx - mx + .5), parseInt(p.y - p.vy - my + .5), parseInt(p.x + .5), parseInt(p.y + .5), c);
                //EFLA(bmp, parseInt(p.x - p.vx - mx + .5), parseInt(p.y - p.vy - my + .5), parseInt(p.x + .5), parseInt(p.y + .5), c);
                //AALine(bmp, parseInt(p.x - p.vx - mx + .5), parseInt(p.y - p.vy - my + .5), parseInt(p.x + .5), parseInt(p.y + .5), c);
            }
        }
    }
}

//init();
render();
function draw() {
    var i=0, j;
    ctx.fillStyle = '#000';
    ctx.fillRect(0, 0, mainCanvas.width, mainCanvas.height);
    ctx.fill();
    
    for (; i < fSolver.numCells; i++) {
        j = fSolver.density[i];                
        if (j > 255) j = 255;
        
        data32[i] =
            (255 << 24) | // alpha, 255 = 100% opacity
            (j << 16) | // blue
            (j << 8) | // green
            j; // red

        image.data[i * 4 + 0] = j;
        image.data[i * 4 + 1] = j;
        image.data[i * 4 + 2] = j;
        image.data[i * 4 + 3] = 255;
    }

    ctx.putImageData(image, 0, 0);
}

function getScrollX() {
    return window.pageXOffset || window.document.documentElement.scrollLeft;
}

function getScrollY() {
    return window.pageYOffset || window.document.documentElement.scrollTop;
}

function getMousePos(event) {    
    event.preventDefault();
    dx = event.pageX - (getScrollX() + mainCanvas.getBoundingClientRect().left)- lastX;
    dy = event.pageY - (getScrollY() + mainCanvas.getBoundingClientRect().top) - lastY;        
    lastX = event.pageX - (getScrollX() + mainCanvas.getBoundingClientRect().left) - cx;
    lastY = event.pageY - (getScrollY() + mainCanvas.getBoundingClientRect().top) - cy;        
    
    return {
        x: lastX,
        y: lastY
    }    
}

function getTouchPos(event) {    
    event.preventDefault();
    dx = event.changedTouches[0].pageX - (getScrollX() + mainCanvas.getBoundingClientRect().left) - lastX;
    dy = event.changedTouches[0].pageY - (getScrollY() + mainCanvas.getBoundingClientRect().top) - lastY;    
    lastX = event.changedTouches[0].pageX - (getScrollX() + mainCanvas.getBoundingClientRect().left);
    lastY = event.changedTouches[0].pageY - (getScrollY() + mainCanvas.getBoundingClientRect().top);        

    lastX = lastX > 0 ? (lastX < canvasWidth ? lastX : canvasWidth) : 0;
    lastY = lastY > 0 ? (lastY < canvasHeight ? lastY : canvasHeight) : 0;        

    return {
        x: lastX,
        y: lastY
    }    
}

function render() {
    var id = requestAnimationFrame(render);
    loop();
}

/** 
 * TODO
 * Render with linked list 
 */

/*function render() {
    var id = requestAnimationFrame(render);
    fSolver.update();

    if (drawLines || drawParticles) {
        renderParticles();
    }

    if (drawFluid) {
        renderFluid();
    }
    frameCount = ++frameCount % 0xffffff;
} */

function uiFeedback() {
    var x = ~~(lastX / GUI.numParticles);
    var y = ~~(lastY / GUI.numParticles);    

    if (x > fSolver.width || y > fSolver.height) return;

    var c = getIndexForCellPosition(x, y);        

    if (!mouseDown) {
        if (dx==undefined) return;
        return;
    } 
        
    if (addVel) {            
        fSolver.uOld[c] += dx / GUI.numParticles * GUI.velMul;      
        fSolver.vOld[c] += dy / GUI.numParticles * GUI.velMul;
    }

    // TODO for more velocity and density, we could pick a radius around the mouse for more fluid injection.

    if (add) {      
        fSolver.source[c] += GUI.density;          
    }    
    
}

function getIndexForCellPosition(i,j) {    
    return i + (fsWidth) * j;
}

mainCanvas.addEventListener('mousemove', onMouseMove);

mainCanvas.addEventListener('mousedown', function(event) {    
    event.preventDefault();
    cx = event.pageX - (getScrollX() + mainCanvas.getBoundingClientRect().left) - lastX;
    cy = event.pageY - (getScrollY() + mainCanvas.getBoundingClientRect().top) - lastY;    
    mouseDown = true;    
});

mainCanvas.addEventListener('mouseup', function(event) {
    event.preventDefault();
    mouseDown = false;
});


// Mobile devices touch events
mainCanvas.addEventListener('touchstart', function(event) {    
    event.preventDefault();    
    lastX = event.changedTouches[0].pageX - (getScrollX() + mainCanvas.getBoundingClientRect().left);
    lastY = event.changedTouches[0].pageY - (getScrollY() + mainCanvas.getBoundingClientRect().top);            
    mouseDown = true;      
});

mainCanvas.addEventListener('touchmove', onTouchMove);

mainCanvas.addEventListener('touchend', function(event) {    
    event.preventDefault();
    mouseDown = false;
})

function loop() {
    for (var i=0; i < fSolver.numCells; i++) {                
        fSolver.vOld[i] = 0;
        fSolver.uOld[i] = 0;        
        fSolver.densityOld[i] = 0;
    }

    uiFeedback();

    for (var i=0; i < fSolver.numCells; i++) {        
        fSolver.densityOld[i] += fSolver.source[i];
        fSolver.source[i] *= GUI.fadeSpeed;
        // additional code to create dissipatation
        fSolver.density[i] *= GUI.densityMul;        
    }

    fSolver.updateVelocity();    
    fSolver.updateDensity();
    draw();
}

function renderParticles() {
    if (drawLines) {
        mainCanvas.colorTransform(fade2alpha);
    } else {
        mainCanvas.colorTransform(fade2black);
    }

    particlesUpdate(mainCanvas);    
    mainCanvas.draw(sparkle.canvas, sparkleMatrix);

    if (frameCount & 1) {
        mainCanvas.draw(fadeImage.canvas, drawMatrix, new ColorTransform(0.1, 0.1, 0.1), BlendMode.ADD);
        fadeImage.applyFilter(fadeImage, new ColorMatrixFilter(blurMatrix));
    }
}

function renderFluid() {
    var fw = fSolver.width;
    var fh = fSolver.height;    
    var y, x;
    for (y = 1; y < fh-1; y++) {
        for (x = 1; x < fw-1; x++) {
            var index = parseInt(y * fw + x);
            var color = ((fSolver.r[index] & 0xff) << 16) |
                        ((fSolver.g[index] & 0xff) << 08) |
                        ( fSolver.b[index] & 0xff);
                
            color = "#" + color.toString(16);                        
            fluidImage.setPixel(x, y, color);                        
        }
    }

    fluidImage.applyFilter(fluidImage, new ColorMatrixFilter(blurMatrix));
}

function onMouseMove(e) {
    if (getMousePos(e).y <= 40) return;    
    //handleForce(getMousePos(e).x, getMousePos(e).y - 40);
}

function onTouchMove(e) {
    if (getTouchPos(e).y <= 40) return;    
}

function handleForce(x, y) {
    var norm_x = x * isw;
    var norm_y = y * ish;
    var vel_x = (x - prevMouse.x) * isw;
    var vel_y = (y - prevMouse.y) * ish;

    addForce(norm_x, norm_y, vel_x, vel_y);

    prevMouse.x = x;
    prevMouse.y = y;
}

function addForce(x, y, dx, dy) {
    var speed = dx * dx + dy * dy * aspectRatio^2;
    if (speed > 0) {
        if (x > 1) x = 1; else if (x < 0) x = 0;
        if (y > 1) y = 1; else if (y < 0) y = 0;
        
        var colorMult = 5;
        var velocityMult = 20.0;
        var index = fSolver.getIndexForNormalizedPosition(x, y);        
        var hue = ( (x + y) * 180 + frameCount) % 360;
        var rgb = hsvToRgb(hue, 1, 1);                
        
        if (fSolver.rgb) {
            fSolver.rOld[index] += rgb.r * colorMult;
            fSolver.gOld[index] += rgb.g * colorMult;
            fSolver.bOld[index] += rgb.b * colorMult;                    
        }

        if (drawFluid || drawLines)
            addParticles(parseInt(x * sw), parseInt(y * sh), NUM_PARTICLES);        
        fSolver.uOld[index] += dx * VEL_MUL;
        fSolver.vOld[index] += dy * VEL_MUL;
    }
}

function restart(e) {   
    if (mainCanvas.mousePos(e).y <= 40) return;
    partList.length = 0;
    fSolver.reset();

    fluidImage.fill("#00ff00");
    mainCanvas.fill('#ffff00');

    drawMode = ++drawMode % 4;    

    if(drawMode == 0) {
        drawFluid = true;
        drawParticles = true;
        drawLines = false;
    } else if (drawMode == 1) {
        drawFluid = true;
        drawParticles = false;
        drawLines = true;
    } else if (drawMode == 2) {
        drawFluid = true;
        drawParticles = false;
        drawLines = false;
    } else if(drawMode == 3){
        drawFluid = false;
        drawParticles = true;
        drawLines = false;
    }
}

function line(bmp, x0, y0, x1, y1, c) {
   var x = x0;
   var y = y0;
   var dx = x1 - x0;
   var dy = y1 - y0;
   var xinc = ( dx > 0 ) ? 1 : -1;
   var yinc = ( dy > 0 ) ? 1 : -1;
   var cumul;
   var i;
   
   dx = (dx ^ (dx >> 31)) - (dx >> 31);
   dy = (dy ^ (dy >> 31)) - (dy >> 31);
   
   bmp.setPixel(x, y, c);
   
   if (dx > dy) {
      cumul = dx >> 1;
      for (i = 1; i <= dx; ++i) {
         x += xinc;
         cumul += dy;
         if (cumul >= dx) {
            cumul -= dx;
            y += yinc;
         }
         bmp.setPixel(x,y, c);
      }
   } else {
      cumul = dy >> 1;
      for (i = 1; i <= dy; ++i) {
         y += yinc;
         cumul += dx;
         if (cumul >= dy) {
            cumul -= dy;
            x += yinc;
         }
         bmp.setPixel(x,y, c);
      }
   }
}

/**
 *    Extremely Fast Line Algorithm
 *    @author     Po-Han Lin (original version: http://www.edepot.com/algorithm.html)
 */
function EFLA(bmp, x1, y1, x2, y2, c) {
    var yLonger = false;
    var shortLen = y2 - y1;
    var longLen = x2 - x1;
    c = "#" + c.toString(16);

    // Tricky and a much faster way to get the absolute value of a number
    if ((shortLen ^ (shortLen >> 31)) - (shortLen >> 31) >  (longLen ^ (longLen >> 31)) - (longLen >> 31)) {
        var swap = shortLen;
        shortLen = longLen;
        longLen = swap;
        yLonger = true;
    } else {
        yLonger = false;
    }

    var inc = (longLen < 0) ? -1 : 1;
    var decInc = (longLen == 0) ? shortLen : shortLen / longLen;
    if (yLonger) {
        for (var i=0; i != longLen; i += inc) {
            bmp.setPixel(x1 + i * decInc, y1+ i, c);
        }
    } else {
        for (i=0; i != longLen; i += inc) {
            bmp.setPixel(x1 + i, y1 + i * decInc, c);
        }
    }
}


/**
 *    Xiaolin Wu's line algorithm with anti aliasing
 *    http://en.wikipedia.org/wiki/Xiaolin_Wu%27s_line_algorithm
 *    http://www.bytearray.org/?p=67
 */

function AALine(bmp, x1, y1, x2, y2, c) {
    var steep = (y2 - y1) < 0 ? -(y2 - y1) : (y2 - y1) > (x2 - x1) < 0 ? -(x2 - x1) : (x2 - x1);
    var dx = x2 - x1;
    var dy = y2 - y1;
    var swap = 0;
    var alpha = 0;

    if (steep) {
        swap = x1; x1 = y1; y1 = swap;
        swap = x2; x2 = y2; y2 = swap;
        swap = dx; dx = dy; dy = swap;
    }
    if (x1 > x2) {
        swap = x1; x1 = x2; x2 = swap;
        swap = y1; y1 = y2; y2 = swap;
    }
    var gradient = dy / dx;
    var xend = x1;
    var yend = y1 + gradient * (xend - x1);
    var xgap = 1.0 - ( (x1 + 0.5) % 1.0);
    var xpx1 = parseInt(xend);
    var ypx1 = yend;

    alpha = (1.0 - (yend % 1.0)) * xgap;
    drawAlpha(bmp, ypx1, xpx1, alpha, c);
    alpha = (yend % 1.0) * xgap;
    drawAlpha(bmp, ypx1, xpx1 + 1, alpha, c);

    var intery = yend + gradient;
    xend = x2;
    yend = y2 + gradient * (xend - x2);
    xgap = (x2 + 0.5) % 1.0;
    var xpx2 = parseInt(xend);
    var ypx2 = yend;

    alpha = (1 - (yend % 1)) * xgap;
    drawAlpha(bmp, ypx2, xpx2, alpha, c);

    alpha = (yend % 1.0) * xgap;
    drawAlpha(bmp, ypx2, xpx2 + 1, alpha, c);

    while (xpx1++ < xpx2) {
        alpha = (1 - (intery % 1.0));
        drawAlpha(bmp, intery, xpx1, alpha, c);
        alpha = (intery % 1.0);
        drawAlpha(bmp, intery + 1, xpx1, alpha, c);
        intery = intery + gradient;
    }
}

function drawAlpha(bmp, x, y, a, c) {
    var hexColor = bmp.getPixel(x, y);
    var rgb = parseInt(hexColor, 16);
    var r0 = (rgb >> 16) & 0xff;
    var g0 = (rgb >> 08) & 0xff;
    var b0 =  rgb & 0xff;
        
    var r1 = (c >> 16) & 0xff;
    var g1 = (c >> 08) & 0xff;
    var b1 = c & 0xff;

    var rc = r1 * a + r0 * (1 - a);
    var gc = g1 * a + g0 * (1 - a);
    var bc = b1 * a + b0 * (1 - a);
    var color = (rc << 16) | (gc << 8) | bc;
    color = "#" + color.toString(16);
    bmp.setPixel(x, y, color);
}

