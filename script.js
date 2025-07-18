let widthImg  = 500;
let heightImg = 500;

let p = new Float32Array(widthImg * heightImg);

let canv = document.getElementById('disp');
let ctx  = canv.getContext('2d');
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