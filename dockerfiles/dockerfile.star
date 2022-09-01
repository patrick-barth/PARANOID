#FROM alpine:latest
FROM ubuntu:18.04

LABEL tool="STAR"
LABEL tool.version="2.7.10a"

MAINTAINER Patrick Barth <patrick.barth@computational.bio.uni-giessen.de>

ARG TMPDIR=/opt
ARG FILE=2.7.10a.tar.gz

#RUN apk update && apk add g++ \
RUN apt-get update -y && apt-get install -y g++ \
	wget \
	make \
	zlib1g-dev && \
	wget https://github.com/alexdobin/STAR/archive/$FILE -O $TMPDIR/$FILE && \
	tar -xf $TMPDIR/$FILE -C $TMPDIR && \
	rm $TMPDIR/$FILE && \
	cd $TMPDIR/STAR-2.7.10a/source && \
	make STAR

ENV PATH=$TMPDIR/STAR-2.7.10a/source:$PATH

ENTRYPOINT ["STAR"]


