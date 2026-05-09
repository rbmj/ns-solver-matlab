#!/bin/sh

ffmpeg -framerate 240 -pattern_type glob -i 'out/*.png' -vf "scale=1600x1200" -c:v libx264 -pix_fmt yuv420p out.mp4
