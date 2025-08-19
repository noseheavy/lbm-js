//
// lattice boltzmann method
// basic reference implementation
//

/* config */

// lattice dimesions
// one pixel per cell
const widthGrid  = 200; // dx
const heightGrid = 200; // dx
const scaleGrid  = 1;   // m
const scaleCell  = scaleGrid / widthGrid;

// sound speed and kinematic viscosity
const actC  = 1;      // m   s^-1
const actNu = 0.0003; // m^2 s^-1

// absolute reference density
// all other densities fractional values
const actRho = 1; // kg m^-3

/* case setup */

// calculate lbm parameters
const lbmC      = 1 / Math.sqrt(3);                          // dx   dt^-1, always same for lbm (notice not one cell per step)
const lbmDeltaT = lbmC * scaleCell / actC;                   // s    dt^-1
const lbmNu     = actNu * lbmDeltaT / scaleCell / scaleCell; // dx^2 dt^-1

// calculate tau w kinematic viscosity
const tau    = lbmNu / (lbmC * lbmC * 1) + 0.5; // nu / (c * c * dt) + 0.5 -> no units
const invTau = 1 / tau;

let windVel = 0.3; // mach
let deltaF  = 0.1; // max deviation from 1

// standard lbgk 
function applyColl(i, j) {
    let uX   = 0;
    let uY   = 0;
    let allF = 0;

    for (let dir = 0; dir < 9; dir ++) {
        uX   += getF(i, j, dir) * getVelX(dir);
        uY   += getF(i, j, dir) * getVelY(dir);
        allF += getF(i, j, dir);
    }

    // r dx dt^-1 -> dx dt^-1
    uX /= allF;
    uY /= allF;

    // temp stability aid
    let scaleP = Math.max(Math.min(allF, 1 + deltaF), 1 - deltaF) / allF;

    // bgk approx
    for (let dir = 0; dir < 9; dir ++) {
        setG(
            i, j, dir,
            scaleP * (
                getF(i, j, dir) + (allF * getEq(dir, uX, uY) - getF(i, j, dir)) * invTau
            )
        );
    }
}

// setup boundaries
function markBnds() {
    for (let i = 1; i < widthGrid - 1; i ++) {
        markWall(i,              0);
        markWall(i, heightGrid - 1);
    }
    for (let i = 0; i < heightGrid; i ++) {
        markWall(            0, i);
        markWall(widthGrid - 1, i);
    }
}
function applyBnds() {
    for (let i = 1; i < widthGrid - 1; i ++) {
        applyVel(i, 0, windVel * actC, 0);
        
        applyNoSlip(i, heightGrid - 1);
    }
    for (let i = 0; i < heightGrid; i ++) {
        applyNoSlip(            0, i);
        applyNoSlip(widthGrid - 1, i);
    }
}

/* init */

let stabView = false;
let fullView = true;

// relative particle density (to ref density)
// 9 directions (img coord system upside down):
// 6 2 5
// 3 0 1
// 7 4 8
// two sets for swapping
let f = new Float32Array(9 * widthGrid * heightGrid); // r
let g = new Float32Array(9 * widthGrid * heightGrid); // r
let swap = false;

// cell velocity and weights for collision operator
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

// setup context and image for display
let disp = document.getElementById('disp');
let ctx  = disp.getContext('2d');
let img  = ctx.createImageData(widthGrid, heightGrid);

// lattice
for (let j = 1; j < heightGrid - 1; j ++) {
    for (let i = 1; i < widthGrid - 1; i ++) {
        for (let dir = 0; dir < 9; dir ++) {
            setF(i, j, dir, getEq(dir, 0, 0)); // assume all points at equilibrium
        }
    }
}

// init boundaries
markBnds();
applyBnds();

let lbmT     = 0;        // s
let maxUPrev = 0.000001; // dx dt^-1

// debug
console.log('width: ' + scaleGrid + ' m');
console.log('c: ' + actC + ' m s^-1');
console.log('re: ' + scaleGrid * (windVel * actC) / actNu);
console.log('tau: ' + tau);
console.log('max lid v: ' + windVel * actC + ' m s^-1');
// console.log('init');

// timing
const t0    = 0.001 * Date.now();
let   tPrev = t0;

/* methods */

// getters for solver
function getVelX(dir) {
    return w[3 * dir    ];
}
function getVelY(dir) {
    return w[3 * dir + 1];
}
function getW(dir) {
    return w[3 * dir + 2];
}

// for density
// f for current density
// g for next
function getF(i, j, dir) {
    if (swap) {
        return g[9 * widthGrid * j + 9 * i + dir];
    } else {
        return f[9 * widthGrid * j + 9 * i + dir];
    }
}
function getG(i, j, dir) {
    if (swap) {
        return f[9 * widthGrid * j + 9 * i + dir];
    } else {
        return g[9 * widthGrid * j + 9 * i + dir];
    }
}
function setF(i, j, dir, val) {
    if (swap) {
        g[9 * widthGrid * j + 9 * i + dir] = val;
    } else {
        f[9 * widthGrid * j + 9 * i + dir] = val;
    }
}
function setG(i, j, dir, val) {
    if (swap) {
        f[9 * widthGrid * j + 9 * i + dir] = val;
    } else {
        g[9 * widthGrid * j + 9 * i + dir] = val;
    }
}
function swapF() {
    swap = !swap;
}

// equil density
function getEq(dir, uX, uY) { // u in dx dt^-1, all div by c included in consts
    return getW(dir) * (
        1 + 3 * (getVelX(dir) * uX + getVelY(dir) * uY) +
        4.5 * (getVelX(dir) * uX + getVelY(dir) * uY) * (getVelX(dir) * uX + getVelY(dir) * uY) -
        1.5 * (uX * uX + uY * uY)
    );
}

// boundaries
function markWall(i, j) { // mark wall with invalid val
    setF(i, j, 0, -1);
    setG(i, j, 0, -1);
}
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
        
        if (i + offsetX > -1 && i + offsetX < widthGrid &&
            j + offsetY > -1 && j + offsetY < heightGrid) {
            if (!isWall(i + offsetX, j + offsetY)) {
                setF(i, j, dir, getF(i + offsetX, j + offsetY, getOppDir(dir)));
            }
        }
    }
}
function applyVel(i, j, uX, uY) { // u in m s^-1, f implied to be unitary
    // check through all non-wall neighbours and
    // pull particles into center, flip non-eq vel and add
    // particles also bounced at boundary
    let cellX = uX / scaleCell * lbmDeltaT;
    let cellY = uY / scaleCell * lbmDeltaT;
    
    for (let dir = 1; dir < 9; dir ++) {
        offsetX = getVelX(dir);
        offsetY = getVelY(dir);
        
        if (i + offsetX > -1 && i + offsetX < widthGrid &&
            j + offsetY > -1 && j + offsetY < heightGrid) {
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

/* main */

// sim loop
function render() {
    // timing
    let t = 0.001 * Date.now();
    if (t - tPrev > 1 / 30) {
        console.log('t: ' + lbmT);
    }
    tPrev = t;

    // render current field
    let maxU = 0.000001;
    for (let j = 0; j < heightGrid; j ++) {
        for (let i = 0; i < widthGrid; i ++) {
            // stored left to right, then top to bottom
            // 4 channels: r, g, b, a
            if (isWall(i, j)) {
                if (stabView) {
                    img.data[4 * (widthGrid * j + i)    ] = 255;
                    img.data[4 * (widthGrid * j + i) + 1] = 0;
                    img.data[4 * (widthGrid * j + i) + 2] = 0;
                } else {
                    img.data[4 * (widthGrid * j + i)    ] = 255;
                    img.data[4 * (widthGrid * j + i) + 1] = 255;
                    img.data[4 * (widthGrid * j + i) + 2] = 255;
                }
                /* if (i == 0) {
                    let col = toRGB(1 - j / (heightGrid - 1));
                    img.data[4 * (widthGrid * j + i)    ] = col[0];
                    img.data[4 * (widthGrid * j + i) + 1] = col[1];
                    img.data[4 * (widthGrid * j + i) + 2] = col[2];
                } */
            } else {
                let uX   = 0;
                let uY   = 0;
                let allF = 0;
                
                for (let dir = 0; dir < 9; dir ++) {
                    uX   += getF(i, j, dir) * getVelX(dir);
                    uY   += getF(i, j, dir) * getVelY(dir);
                    allF += getF(i, j, dir);
                }
                
                if (stabView) {
                    let offsetC = 127.5 * (allF - 1) / deltaF;
                    img.data[4 * (widthGrid * j + i)    ] = 127.5 + offsetC;
                    img.data[4 * (widthGrid * j + i) + 1] = 127.5 + offsetC;
                    img.data[4 * (widthGrid * j + i) + 2] = 127.5 + offsetC;
                } else {
                    let a = Math.sqrt(uY * uY + uX * uX);
                    // a /= allF; (allF ~= 1, skip)
                    
                    if (fullView) {
                        if (a > maxU) {
                            maxU = a;
                        }
                        
                        a /= maxUPrev;
                    } else {
                        a /= windVel * c;
                    }
                    
                    let col = toRGB(a);
                    img.data[4 * (widthGrid * j + i)    ] = col[0];
                    img.data[4 * (widthGrid * j + i) + 1] = col[1];
                    img.data[4 * (widthGrid * j + i) + 2] = col[2];
                }
            }
            
            // set alpha
            img.data[4 * (widthGrid * j + i) + 3] = 255;
        }
    }
    if (fullView && !stabView) {
        maxUPrev = maxU;
    }
    ctx.putImageData(img, 0, 0);
    
    for (let rep = 0; rep < 1; rep ++) {
        // stream
        for (let j = 1; j < heightGrid - 1; j ++) {
            for (let i = 1; i < widthGrid - 1; i ++) {
                // take all particles that should be streamed into this cell
                for (let dir = 1; dir < 9; dir ++) { // skip stationary parts
                    // use micro vel to get offsets
                    setG(
                        i, j, dir,
                        getF(i - getVelX(dir), j - getVelY(dir), dir)
                    );
                }
                setG(i, j, 0, getF(i, j, 0)); // not to forget
            }
        }
        swapF(); // bring changes forward
        
        // collide
        for (let j = 1; j < heightGrid - 1; j ++) {
            for (let i = 1; i < widthGrid - 1; i ++) {
                applyColl(i, j);
            }
        }
        swapF(); // bring changes forward
        
        // finish step and ready for next
        applyBnds();
        lbmT += lbmDeltaT;
    }
    
    // wait for frame
    requestAnimationFrame(render);
    // setTimeout(render, 1000);
}

// begin rendering
render();