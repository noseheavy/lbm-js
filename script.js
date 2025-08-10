//
// lattice boltzmann method
// basic reference implementation
//

/* config */

// lattice config
// one pixel per cell
let widthImg  = 200; // dx
let heightImg = 200; // dx
let scaleGrid = 1;   // m

// auto
let scaleCell = scaleGrid / widthImg; // m

// time step
let deltaT = 0.0029; // s, increasing lowers c

// kinematic viscosity
let nu = 0.000299;   // m^2 s^-1

// absolute reference density
// all other densities relative (fractional) values
let rho = 1;         // kg m^-3

let windVel = 0.3;   // mach
let deltaF  = 0.1;   // max deviation from 1
let stabView = false;
let fullView = true;

// relative particle density (relative to ref density)
// 9 directions:
// 6 2 5
// 3 0 1
// 7 4 8
// two sets for swapping
let f = new Float32Array(9 * widthImg * heightImg); // r
let g = new Float32Array(9 * widthImg * heightImg); // r
let swap = false;

// velocity and weights for collision operator
const w = new Float32Array([
     0,  0, 4 /  9,
     1,  0, 1 /  9,
     0,  1, 1 /  9,
    -1,  0, 1 /  9,
     0, -1, 1 /  9,
     1,  1, 1 / 36,
    -1,  1, 1 / 36,
    -1, -1, 1 / 36,
     1, -1, 1 / 36
]);

/* methods */

// getters for solver
function getVelX(i) {
    return w[3 * i    ];
}
function getVelY(i) {
    return w[3 * i + 1];
}
function getW(i) {
    return w[3 * i + 2];
}

// for density
// f for current density
// g for next
function getF(x, y, i) {
    if (swap) {
        return g[9 * widthImg * y + 9 * x + i];
    } else {
        return f[9 * widthImg * y + 9 * x + i];
    }
}
function getG(x, y, i) {
    if (swap) {
        return f[9 * widthImg * y + 9 * x + i];
    } else {
        return g[9 * widthImg * y + 9 * x + i];
    }
}
function setF(x, y, i, val) {
    if (swap) {
        g[9 * widthImg * y + 9 * x + i] = val;
    } else {
        f[9 * widthImg * y + 9 * x + i] = val;
    }
}
function setG(x, y, i, val) {
    if (swap) {
        f[9 * widthImg * y + 9 * x + i] = val;
    } else {
        g[9 * widthImg * y + 9 * x + i] = val;
    }
}
function swapF() {
    swap = !swap;
}

// equil density
function getEq(i, uX, uY) { // u in dx dt^-1, all div by c included in consts
    return getW(i) * (
        1 + 3 * (getVelX(i) * uX + getVelY(i) * uY) +
        4.5 * (getVelX(i) * uX + getVelY(i) * uY) * (getVelX(i) * uX + getVelY(i) * uY) -
        1.5 * (uX * uX + uY * uY)
    );
}

// boundaries
function isWall(i, j) {
    return getF(i, j, 0) < 0;
}
function getOppDir(dir) {
    switch (dir) {
        case 1: return 3;
        case 2: return 4;
        case 3: return 1;
        case 4: return 2;
        case 5: return 7;
        case 6: return 8;
        case 7: return 5;
        case 8: return 6;
    }
}
function applyNoSlip(i, j) {
    // check through all non-wall neighbours and
    // pull particles into center, flip vel
    // particles effectively bounced at boundary
    for (let dir = 1; dir < 9; dir ++) {
        offsetX = getVelX(dir);
        offsetY = getVelY(dir);
        
        if (i + offsetX > -1 && i + offsetX < widthImg &&
            j + offsetY > -1 && j + offsetY < heightImg) {
            if (!isWall(i + offsetX, j + offsetY)) {
                setF(i, j, dir, getF(i + offsetX, j + offsetY, getOppDir(dir)));
            }
        }
    }
}
function applyVel(i, j, uX, uY) { // m s^-1
    // check through all non-wall neighbours and
    // pull particles into center, flip non-eq vel and add
    // particles also bounced at boundary
    let cellX = uX / scaleCell * deltaT;
    let cellY = uY / scaleCell * deltaT;
    
    for (let dir = 1; dir < 9; dir ++) {
        offsetX = getVelX(dir);
        offsetY = getVelY(dir);
        
        if (i + offsetX > -1 && i + offsetX < widthImg &&
            j + offsetY > -1 && j + offsetY < heightImg) {
            if (!isWall(i + offsetX, j + offsetY)) {
                setF(
                    i, j, dir,
                    getF(i + offsetX, j + offsetY, getOppDir(dir)) - getEq(getOppDir(dir), cellX, cellY) +
                    getEq(dir, cellX, cellY)
                );
            }
        }
    }
}

function toRGB(a) {
    a = 1 - Math.min(1, Math.max(0, a));

    let r = 0;
    if (a < 0.25) {
        r = 255;
    } else if (a < 0.5) {
        r = 1 - (a - 0.25) * 4;
        r *= 255;
    }
    
    let gr = 0;
    if (a < 0.25) {
        gr = a * 4;
        gr *= 255;
    } else if (a < 0.75) {
        gr = 255;
    } else {
        gr = 1 - (a - 0.75) * 4;
        gr *= 255;
    }
    
    let b = 0;
    if (a > 0.75) {
        b = 255;
    } else if (a > 0.5) {
        b = (a - 0.5) * 4;
        b *= 255;
    }
    
    return [r, gr, b];
}

/* init */

// calculate speed of sound w lattice config, tau w kinematic viscosity
const c     = 1 / Math.sqrt(3);                    // dx dt^-1, always same for lbm
const realC = scaleCell / deltaT * c;              // m s^-1, notice not just scale divided by time
let tau     = nu / (realC * realC * deltaT) + 0.5; // no units
let invTau  = 1 / tau;

// setup context and image for display
let disp = document.getElementById('disp');
let ctx  = disp.getContext('2d');
let img  = ctx.createImageData(widthImg, heightImg);

// lattice
for (let j = 1; j < heightImg - 1; j ++) {
    for (let i = 1; i < widthImg - 1; i ++) {
        for (let dir = 0; dir < 9; dir ++) {
            setF(i, j, dir, getEq(dir, 0, 0)); // assume all points at equilibrium
        }
    }
}

// mark boundaries
for (let i = 1; i < widthImg - 1; i ++) {
    // mark wall in current
    setF(i,             0, 0, -1);
    setF(i, heightImg - 1, 0, -1);
    
    // mark wall in next
    setG(i,             0, 0, -1);
    setG(i, heightImg - 1, 0, -1);
}
for (let i = 0; i < heightImg; i ++) {
    setF(           0, i, 0, -1);
    setF(widthImg - 1, i, 0, -1);
    
    setG(           0, i, 0, -1);
    setG(widthImg - 1, i, 0, -1);
}

// apply boundaries
for (let i = 1; i < widthImg - 1; i ++) {
    applyVel(i, 0, 0, 0);
    applyNoSlip(i, heightImg - 1);
}
for (let i = 0; i < heightImg; i ++) {
    applyNoSlip(           0, i);
    applyNoSlip(widthImg - 1, i);
}

// sim time
let accumT = 0;             // s
let t0 = Date.now() * 0.001; // real time
let tPrev = t0;

let maxUPrev = 0.001;

// debug
/*for (let dir = 0; dir < 9; dir ++) {
setF(35, 35, dir, getEq(dir, 0.01, 0));
}*/
/* for (let dir = 0; dir < 9; dir ++) {
    // setP(80, 80, dir, getP(10, 10, dir) + 0.2);
    for (let i = 0; i < 30; i ++) {
    for (let j = 0; j < 30; j ++) {
    setP(70 + i, 70 + j, dir, getEq(dir, 0, 0.02));
    }}
    // console.log(getP(0, 1, dir));
    console.log(getP(0, 1, dir));
} */
console.log('width: ' + scaleCell * widthImg + ' m');
console.log('c: ' + Math.round(100 * realC) / 100 + ' m s^-1');
console.log('re: ' + widthImg * (windVel * c) / (nu / scaleCell / scaleCell * deltaT));
// console.log('re: ' + widthImg * (windVel * c) / ((tau - 0.5) * c * c));
console.log('tau: ' + tau);
// console.log('t: ' + t0);
console.log('max lid v: ' + windVel * realC + ' m s^-1');
console.log('init');

/* main */

// sim loop
function render() {
    // timing
    let t = (Date.now() * 0.001) - t0;
    if (t - tPrev > 1 / 30) {
        console.log('t: ' + accumT);
    }
    tPrev = t;

    // render current field
    let maxU = 0.001;
    for (let j = 0; j < heightImg; j ++) {
        for (let i = 0; i < widthImg; i ++) {
            // stored left to right, then top to bottom
            // 4 channels: r, g, b, a
            if (isWall(i, j)) {
                if (stabView) {
                    img.data[4 * (widthImg * j + i)    ] = 255;
                    img.data[4 * (widthImg * j + i) + 1] = 0;
                    img.data[4 * (widthImg * j + i) + 2] = 0;
                } else {
                    img.data[4 * (widthImg * j + i)    ] = 255;
                    img.data[4 * (widthImg * j + i) + 1] = 255;
                    img.data[4 * (widthImg * j + i) + 2] = 255;
                }
                /* if (i == 0) {
                    let col = toRGB(1 - j / (heightImg - 1));
                    img.data[4 * (widthImg * j + i)    ] = col[0];
                    img.data[4 * (widthImg * j + i) + 1] = col[1];
                    img.data[4 * (widthImg * j + i) + 2] = col[2];
                } */
            } else {
                let uX   = 0;
                let uY   = 0;
                let allP = 0;
                
                for (let dir = 0; dir < 9; dir ++) {
                    uX += getF(i, j, dir) * getVelX(dir);
                    uY += getF(i, j, dir) * getVelY(dir);
                    allP += getF(i, j, dir);
                }
                
                if (stabView) {
                    let offsetC = 127.5 * (allP - 1) / deltaF;
                    img.data[4 * (widthImg * j + i)    ] = 127.5 + offsetC;
                    img.data[4 * (widthImg * j + i) + 1] = 127.5 + offsetC;
                    img.data[4 * (widthImg * j + i) + 2] = 127.5 + offsetC;
                } else {
                    let a = Math.sqrt(uY * uY + uX * uX);
                    
                    if (fullView) {
                        if (a > maxU) {
                            maxU = a;
                        }
                        
                        a /= maxUPrev;
                    } else {
                        a /= windVel * c;
                    }
                    
                    let col = toRGB(a);
                    img.data[4 * (widthImg * j + i)    ] = col[0];
                    img.data[4 * (widthImg * j + i) + 1] = col[1];
                    img.data[4 * (widthImg * j + i) + 2] = col[2];
                }
            }
            
            // set alpha
            img.data[4 * (widthImg * j + i) + 3] = 255;
        }
    }
    if (fullView && !stabView) {
        maxUPrev = maxU;
    }
    ctx.putImageData(img, 0, 0);
    
    for (let rep = 0; rep < 1; rep ++) {
        // stream
        for (let j = 1; j < heightImg - 1; j ++) {
            for (let i = 1; i < widthImg - 1; i ++) {
                // take all particles that should be streamed for this cell
                for (let dir = 0; dir < 9; dir ++) {
                    // use micro vel to get offsets
                    setG(
                        i, j, dir,
                        getF(i - getVelX(dir), j - getVelY(dir), dir)
                    );
                }
            }
        }
        
        // bring changes forward
        swapF();
        
        // collide
        for (let j = 1; j < heightImg - 1; j ++) {
            for (let i = 1; i < widthImg - 1; i ++) {
                let uX   = 0;
                let uY   = 0;
                let allP = 0;
                
                for (let dir = 0; dir < 9; dir ++) {
                    uX   += getF(i, j, dir) * getVelX(dir);
                    uY   += getF(i, j, dir) * getVelY(dir);
                    allP += getF(i, j, dir);
                }
                
                // r dx dt^-1 -> dx dt^-1
                uX /= allP;
                uY /= allP;
                
                // temp stability aid
                let scaleP = Math.max(Math.min(allP, 1 + deltaF), 1 - deltaF) / allP;
                
                // bgk approx
                for (let dir = 0; dir < 9; dir ++) {
                    setG(
                        i, j, dir,
                        scaleP * (
                            getF(i, j, dir) + (allP * getEq(dir, uX, uY) - getF(i, j, dir)) * invTau
                        )
                    );
                }
            }
        }
        
        // bring changes forward
        swapF();
        
        // apply boundaries
        for (let i = 1; i < widthImg - 1; i ++) {
            applyVel(i, 0, windVel * realC, 0); // Math.tanh(accumT) * 
            applyNoSlip(i, heightImg - 1);
        }
        for (let i = 0; i < heightImg; i ++) {
            applyNoSlip(           0, i);
            applyNoSlip(widthImg - 1, i);
        }
        
        accumT += deltaT;
    }
    
    // wait for frame
    requestAnimationFrame(render);
    // setTimeout(render, 1000);
}

// begin rendering
render();