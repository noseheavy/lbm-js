//
// lattice boltzmann method
// basic reference implementation
//

/* config */

// lattice config
// one pixel per cell
let widthImg  = 100; // dx
let heightImg = 100; // dx
let scaleCell = 0.1; // m

// time step
let dt = 0.002; // s

// kinematic viscosity
let nu = 0.0003 // dx^2 dt^-1

// absolute reference density
// all other densities relative (fractional) values
let rho = 1; // kg m^-3

let windVel = 0.0005; // mach

// norm particle density, relative to ref density
// 9 directions:
// 6 2 5
// 3 0 1
// 7 4 8
// two sets for swapping
let p = new Float32Array(9 * widthImg * heightImg);
let q = new Float32Array(9 * widthImg * heightImg);
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
// p for current density
// q for next
function getP(x, y, i) {
    if (swap) {
        return q[9 * widthImg * y + 9 * x + i];
    } else {
        return p[9 * widthImg * y + 9 * x + i];
    }
}
function getQ(x, y, i) {
    if (swap) {
        return p[9 * widthImg * y + 9 * x + i];
    } else {
        return q[9 * widthImg * y + 9 * x + i];
    }
}
function setP(x, y, i, val) {
    if (swap) {
        q[9 * widthImg * y + 9 * x + i] = val;
    } else {
        p[9 * widthImg * y + 9 * x + i] = val;
    }
}
function setQ(x, y, i, val) {
    if (swap) {
        p[9 * widthImg * y + 9 * x + i] = val;
    } else {
        q[9 * widthImg * y + 9 * x + i] = val;
    }
}
function swapP() {
    swap = !swap;
}

// equil density
function getEq(i, uX, uY) { // dx dt^-1, div by c included in consts
    return getW(i) * (
        1 + 3 * (getVelX(i) * uX + getVelY(i) * uY) +
        4.5 * (getVelX(i) * uX + getVelY(i) * uY) * (getVelX(i) * uX + getVelY(i) * uY) -
        1.5 * (uX * uX + uY * uY)
    );
}

// boundaries
function isWall(i, j) {
    return getP(i, j, 0) < 0;
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
                setP(i, j, dir, getP(i + offsetX, j + offsetY, getOppDir(dir)));
            }
        }
    }
}
function applyVel(i, j, uX, uY) { // m s^-1
    // check through all non-wall neighbours and
    // pull particles into center, flip non-eq vel and add
    // particles also bounced at boundary
    let cellX = uX / scaleCell * dt;
    let cellY = uY / scaleCell * dt;
    
    for (let dir = 1; dir < 9; dir ++) {
        offsetX = getVelX(dir);
        offsetY = getVelY(dir);
        
        if (i + offsetX > -1 && i + offsetX < widthImg &&
            j + offsetY > -1 && j + offsetY < heightImg) {
            if (!isWall(i + offsetX, j + offsetY)) {
                setP(
                    i, j, dir,
                    getP(i + offsetX, j + offsetY, getOppDir(dir)) - getEq(getOppDir(dir), cellX, cellY) +
                    getEq(dir, cellX, cellY)
                );
            }
        }
    }
}

/* init */

// calculate speed of sound w lattice config, tau w kinematic viscosity
const c     = 1 / Math.sqrt(3);   // dx dt^-1, always same for lbm
const realC = scaleCell / dt * c; // m s^-1, notice not just scale divided by time
let tau = nu / (c * c * dt) + 0.5;

// setup context and image for display
let disp = document.getElementById('disp');
let ctx  = disp.getContext('2d');
let img  = ctx.createImageData(widthImg, heightImg);

// lattice
for (let j = 1; j < heightImg - 1; j ++) {
    for (let i = 1; i < widthImg - 1; i ++) {
        for (let dir = 0; dir < 9; dir ++) {
            setP(i, j, dir, getEq(dir, 0, 0)); // assume all points at equilibrium
        }
    }
}

// mark boundaries
for (let i = 1; i < widthImg - 1; i ++) {
    // mark wall in current
    setP(i,             0, 0, -1);
    setP(i, heightImg - 1, 0, -1);
    
    // mark wall in next
    setQ(i,             0, 0, -1);
    setQ(i, heightImg - 1, 0, -1);
}
for (let i = 0; i < heightImg; i ++) {
    setP(           0, i, 0, -1);
    setP(widthImg - 1, i, 0, -1);
    
    setQ(           0, i, 0, -1);
    setQ(widthImg - 1, i, 0, -1);
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
let accumT = 0;
let t0 = Date.now() / 1000;

// debug
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
console.log('c: ' + Math.round(10 * realC) / 10 + ' m s^-1');
console.log('re: ' + widthImg * (windVel * c) / nu);
console.log('tau: ' + tau);
// console.log('t: ' + t0);
console.log('max lid v: ' + windVel * realC + ' m s^-1');
console.log('init');

/* main */

// sim loop
function render() {
    // timing
    let t = Date.now() / 1000 - t0;
    // console.log('t: ' + accumT);

    // render current field
    for (let j = 0; j < heightImg; j ++) {
        for (let i = 0; i < widthImg; i ++) {
            // stored left to right, then top to bottom
            // 4 channels: r, g, b, a
            /* img.data[4 * (widthImg * j + i)    ] = 255 * (0.5 + 0.5 * Math.sin(t));
            img.data[4 * (widthImg * j + i) + 1] = 255 * (0.5 + 0.5 * Math.sin(t * 0.9 + 1 * Math.PI / 3));
            img.data[4 * (widthImg * j + i) + 2] = 255 * (0.5 + 0.5 * Math.sin(t * 0.8 + 2 * Math.PI / 3));
            
            img.data[4 * (widthImg * j + i) + 3] = 255; */
            if (isWall(i, j)) {
                img.data[4 * (widthImg * j + i)    ] = 255;
                img.data[4 * (widthImg * j + i) + 1] = 0;
                img.data[4 * (widthImg * j + i) + 2] = 0;
            } else {
                /* img.data[4 * (widthImg * j + i)    ] = 255 * (0.5 + 0.5 * Math.sin(t));
                img.data[4 * (widthImg * j + i) + 1] = 255 * (0.5 + 0.5 * Math.sin(t * 0.9 + 1 * Math.PI / 3));
                img.data[4 * (widthImg * j + i) + 2] = 255 * (0.5 + 0.5 * Math.sin(t * 0.8 + 2 * Math.PI / 3)); */
                let uX   = 0;
                let uY   = 0;
                let allP = 0;
                
                for (let dir = 0; dir < 9; dir ++) {
                    uX += getP(i, j, dir) * getVelX(dir);
                    uY += getP(i, j, dir) * getVelY(dir);
                    allP += getP(i, j, dir);
                }
                
                let offsetC = 127.5 * Math.sqrt(uY * uY + uX * uX) / (windVel * c / scaleCell * dt);
                // let offsetC = 127.5 * uX / (windVel * c / scaleCell * dt);
                // let offsetC = 127.5 * 20000 * (allP - 1);
                img.data[4 * (widthImg * j + i)    ] = 127.5 + offsetC;
                img.data[4 * (widthImg * j + i) + 1] = 127.5 + offsetC;
                img.data[4 * (widthImg * j + i) + 2] = 127.5 + offsetC;
            }
            
            // set alpha
            img.data[4 * (widthImg * j + i) + 3] = 255;
        }
    }
    ctx.putImageData(img, 0, 0);
    
    for (let rep = 0; rep < 5; rep ++) {
    // stream
    for (let j = 1; j < heightImg - 1; j ++) {
        for (let i = 1; i < widthImg - 1; i ++) {
            // take all particles that should be streamed for this cell
            for (let dir = 0; dir < 9; dir ++) {
                // use micro vel to get offsets
                setQ(
                    i, j, dir,
                    getP(i - getVelX(dir), j - getVelY(dir), dir)
                );
            }
        }
    }
    
    // bring changes forward
    swapP();
    
    // collide
    for (let j = 1; j < heightImg - 1; j ++) {
        for (let i = 1; i < widthImg - 1; i ++) {
            let uX   = 0;
            let uY   = 0;
            let allP = 0;
            
            for (let dir = 0; dir < 9; dir ++) {
                uX += getP(i, j, dir) * getVelX(dir);
                uY += getP(i, j, dir) * getVelY(dir);
                allP += getP(i, j, dir);
            }
            
            // temp stability aid
            let scaleP = Math.max(Math.min(allP, 2), 0) / allP;
            allP = scaleP * allP;
            
            // bgk approx
            for (let dir = 0; dir < 9; dir ++) {
                setQ(
                    i, j, dir,
                    scaleP * getP(i, j, dir) + 
                    (allP * getEq(dir, uX, uY) - scaleP * getP(i, j, dir)) / tau
                );
            }
        }
    }
    
    // bring changes forward
    swapP();
    
    // apply boundaries
    for (let i = 1; i < widthImg - 1; i ++) {
        applyVel(i, 0, windVel * Math.tanh(accumT) * c, 0);
        applyNoSlip(i, heightImg - 1);
    }
    for (let i = 0; i < heightImg; i ++) {
        applyNoSlip(           0, i);
        applyNoSlip(widthImg - 1, i);
    }
    
    accumT += dt;
    // console.log('lid v: ' + 0.0003 * Math.tanh(accumT) * c + ' m s^-1');
    } // rep
    
    // wait for frame
    requestAnimationFrame(render);
    // setTimeout(render, 1000);
}

// begin rendering
render();