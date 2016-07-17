cd ../output/animation
pwd
avconv -r 30 -start_number 1 -i theta%05d.tga test.mp4

