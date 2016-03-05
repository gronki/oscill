#!/bin/bash

set -e

PREFX=$1

mkdir -p mpeg

LISTF=$(ls png/* | grep "png\/${PREFX}-1[0-9]*\.png")


ffmpeg -i png/${PREFX}-1%06d.png -r 24 -b:v 8M -s 1600x900 -vcodec h264   mpeg/${PREFX}.mpeg

rm $LISTF
