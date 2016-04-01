#!/bin/bash

avconv -framerate 20 -i frame%d.jpg -c:v h264 -crf 1 out.mov
