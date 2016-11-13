#!/bin/sh

./build/bin/main pics/00908u.bmp res/align1.bmp --align
./build/bin/main pics/00926u.bmp res/align2.bmp --align
./build/bin/main pics/00932u.bmp res/align3.bmp --align
./build/bin/main pics/00936u.bmp res/align4.bmp --align
./build/bin/main pics/00976u.bmp res/align5.bmp --align

./build/bin/main pics/aston.bmp res/sobel-x.bmp --sobel-x
./build/bin/main pics/aston.bmp res/sobel-y.bmp --sobel-y
./build/bin/main pics/aston.bmp res/unsharp.bmp --unsharp
./build/bin/main pics/aston.bmp res/gray-world.bmp --gray-world
./build/bin/main pics/aston.bmp res/resize-u.bmp --resize 4.3
./build/bin/main pics/aston.bmp res/resize-d.bmp --resize 0.4
./build/bin/main pics/aston.bmp res/autocontrast.bmp --autocontrast 0.1
./build/bin/main pics/aston.bmp res/gaussian.bmp --gaussian 2 5
./build/bin/main pics/aston.bmp res/gaussian-separable.bmp --gaussian-separable 2 5
./build/bin/main pics/aston.bmp res/median.bmp --median 5
./build/bin/main pics/aston.bmp res/median-linear.bmp --median-linear 5
./build/bin/main pics/aston.bmp res/canny.bmp --canny 150 260

./build/bin/main pics/00908u.bmp res/align1-subpixel.bmp --align --subpixel 2
./build/bin/main pics/00908u.bmp res/align1-gray-world.bmp --align --gray-world
./build/bin/main pics/00908u.bmp res/align1-unsharp.bmp --align --unsharp
./build/bin/main pics/00908u.bmp res/align1-autocontrast.bmp --align --autocontrast 0.1
