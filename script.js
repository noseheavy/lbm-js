//
// lattice boltzmann method
// basic reference implementation
//

/* config */

// lattice config
// one pixel per cell
let widthImg  = 150;
let heightImg = 150;
let scaleCell = 0.01;

// time step
let dt = 0.01;

// kinematic viscosity
let nu = 1;

// absolute reference density
// all other densities relative (fractional) values
let rho = 1;

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
let w = new Float32Array([
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

// equil density
function getEq(i, uX, uY) {
    return getW(i) * (
        1 + 3 * (getVelX(i) * uX + getVelY(i) * uY) +
        4.5 * (getVelX(i) * uX + getVelY(i) * uY) * (getVelX(i) * uX + getVelY(i) * uY) -
        1.5 * (uX * uX + uY * uY)
    );
}

/* init */

// calculate speed of sound w lattice config, tau w kinematic viscosity
let c   = scaleCell / dt / Math.sqrt(3); // notice not just scale divided by time
let tau = nu / (c * c * dt) + 0.5;

// setup context and image for display
let disp = document.getElementById('disp');
let ctx  = disp.getContext('2d');
let img  = ctx.createImageData(widthImg, heightImg);

let t0 = Date.now() / 1000;

// lattice
for (let j = 1; j < heightImg - 1; j ++) {
    for (let i = 1; i < widthImg - 1; i ++) {
        for (let dir = 0; dir < 9; dir ++) {
            setP(i, j, dir, getEq(dir, 0, 0)); // assume all points at equilibrium
        }
    }
}

// boundaries
// todo: handle no-slip
for (let i = 0; i < widthImg; i ++) {
    for (let dir = 0; dir < 9; dir ++) {
        setP(i,             0, dir, -1);
        setP(i, heightImg - 1, dir, -1);
    }
}
for (let i = 1; i < heightImg - 1; i ++) {
    for (let dir = 0; dir < 9; dir ++) {
        setP(           0, i, dir, -1);
        setP(widthImg - 1, i, dir, -1);
    }
}

// debug
for (let dir = 0; dir < 9; dir ++) {
    console.log(getP(0, 1, dir));
}
console.log('init');

/* main */

// sim loop
function render() {
    // timing
    let t = Date.now() / 1000 - t0;

    // render current field
    for (let j = 0; j < heightImg; j ++) {
        for (let i = 0; i < widthImg; i ++) {
            // stored left to right, then top to bottom
            // 4 channels: r, g, b, a
            /* img.data[4 * (widthImg * j + i)    ] = 255 * (0.5 + 0.5 * Math.sin(t));
            img.data[4 * (widthImg * j + i) + 1] = 255 * (0.5 + 0.5 * Math.sin(t * 0.9 + 1 * Math.PI / 3));
            img.data[4 * (widthImg * j + i) + 2] = 255 * (0.5 + 0.5 * Math.sin(t * 0.8 + 2 * Math.PI / 3));
            
            img.data[4 * (widthImg * j + i) + 3] = 255; */
            if (getP(i, j, 0) < 0) {
                img.data[4 * (widthImg * j + i)    ] = 155;
                img.data[4 * (widthImg * j + i) + 1] = 155;
                img.data[4 * (widthImg * j + i) + 2] = 155;
            } else {
                img.data[4 * (widthImg * j + i)    ] = 255 * (0.5 + 0.5 * Math.sin(t));
                img.data[4 * (widthImg * j + i) + 1] = 255 * (0.5 + 0.5 * Math.sin(t * 0.9 + 1 * Math.PI / 3));
                img.data[4 * (widthImg * j + i) + 2] = 255 * (0.5 + 0.5 * Math.sin(t * 0.8 + 2 * Math.PI / 3));
            }
            
            img.data[4 * (widthImg * j + i) + 3] = 255;
        }
    }
    ctx.putImageData(img, 0, 0);
    
    // stream
    for (let j = 1; j < heightImg - 1; j ++) {
        for (let i = 1; i < widthImg - 1; i ++) {
        }
    }
    
    // collide
    for (let j = 1; j < heightImg - 1; j ++) {
        for (let i = 1; i < widthImg - 1; i ++) {
        }
    }
    
    // todo: boundaries (before next stream...)
    
    // wait for frame
    requestAnimationFrame(render);
}

// begin rendering
render();