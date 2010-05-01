#!/bin/bash

PATH=.:$PATH

# Measure error of velocity evaluation
# against analytical Lamb-Oseen velocity
init | rm-comments.sh > PARTICLES
vel-eval < PARTICLES > EVALUATED_VELOCITY
analytical-vel < PARTICLES > ANALYTICAL_VELOCITY
velocity-error.py EVALUATED_VELOCITY ANALYTICAL_VELOCITY

