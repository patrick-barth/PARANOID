FROM alpine:latest

LABEL tool="STAR"
LABEL tool.version="2.7.10a"

MAINTAINER Patrick Barth <patrick.barth@computational.bio.uni-giessen.de>

ARG TMPDIR=/opt
ARG FILE=2.7.10a.tar.gz

RUN apk update && apk add bash libgcc libstdc++ libgomp g++ wget make zlib-dev musl-dev libc-dev patch && \
	wget -q https://github.com/alexdobin/STAR/archive/$FILE -O $TMPDIR/$FILE && \
	tar xf $TMPDIR/$FILE -C $TMPDIR && \
	rm $TMPDIR/$FILE && \
	cd $TMPDIR/STAR-2.7.10a/source && \
        wget -O - https://patch-diff.githubusercontent.com/raw/alexdobin/STAR/pull/1651.diff | patch -p 2 && \
	make STAR && \
    cp STAR /usr/local/bin/ && \
    cd ../.. && rm -rf STAR-2.7.10a && \
    apk del gcc g++ wget make zlib-dev musl-dev libc-dev patch


ENTRYPOINT ["STAR"]


