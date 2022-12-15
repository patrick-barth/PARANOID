#!/usr/bin/env bash

ACCOUNT="pbarth"
VERSION="1.0"
REPO="hub.docker.com"

DOCKERDIR=$1
OUTPUTDIR=$2

set -e

echo "Directory containing dockerfiles: $DOCKERDIR"
echo "Directory to save images to:      $OUTPUTDIR"

for d in `ls $DOCKERDIR/`; do
    NAME=${d/dockerfile.}
    echo "Pulling image $NAME"
    singularity pull --name $OUTPUTDIR/$ACCOUNT-$NAME-$VERSION.img docker://$ACCOUNT/$NAME:$VERSION > /dev/null
done