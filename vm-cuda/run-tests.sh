#!/bin/bash

PATH=.:$PATH

PARTICLES=$(mktemp)
EVALUATED_VELOCITY=$(mktemp)
ANALYTICAL_VELOCITY=$(mktemp)

cat <<EOF
### Measure error of velocity evaluation
### against analytical Lamb-Oseen velocity
EOF
init | rm-comments.sh > $PARTICLES
vel-eval < $PARTICLES > $EVALUATED_VELOCITY
analytical-vel < $PARTICLES > $ANALYTICAL_VELOCITY
velocity-error.py $EVALUATED_VELOCITY $ANALYTICAL_VELOCITY

