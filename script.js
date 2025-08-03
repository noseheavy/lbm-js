//
// lattice boltzmann method
// basic reference implementation
//

// config

// lattice config
// one pixel per cell
let widthImg  = 500;
let heightImg = 500;
let scaleCell = 1;

// time step
let dt = 0.01;

// kinematic viscosity
let nu = 1;

// reference density
let rho = 1;

// norm particle density, relative to ref density
// 9 directions:
// 6 2 5
// 3 0 1
// 7 4 8
let p = new Float32Array(9 * widthImg * heightImg);

// weights for collision operator
let w = new Float32Array([
     1,  0,
     0,  1,
    -1,  0,
     0, -1,
     1,  1,
    -1,  1
    -1, -1,
     1, -1
]);

// init

// calculate speed of sound with lattice config, tau using kinematic viscosity
let c   = scaleCell / dt / Math.sqrt(3); // notice not simply scale divided by time
let tau = nu / (c * c * dt) + 0.5;

let disp = document.getElementById('disp');
let ctx  = disp.getContext('2d');
let img  = ctx.createImageData(widthImg, heightImg);

console.log('init');

for (let j = 0; j < heightImg; j++) {
    for (let i = 0; i < widthImg; i++) {
        // stored left to right, then top to bottom
        // 4 channels: r, g, b, a
        img.data[4 * (widthImg * j + i)    ] = 255;
        img.data[4 * (widthImg * j + i) + 1] = 0;
        img.data[4 * (widthImg * j + i) + 2] = 0;
        img.data[4 * (widthImg * j + i) + 3] = 255;
    }
}

ctx.putImageData(img, 0, 0);