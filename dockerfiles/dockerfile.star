FROM alpine:3.16

LABEL tool="STAR"
LABEL tool.version="2.7.10b"

MAINTAINER Patrick Barth <patrick.barth@computational.bio.uni-giessen.de>

ARG FILE=2.7.10b.tar.gz

RUN apk update && apk add bash wget && \
	wget -q https://github.com/alexdobin/STAR/archive/$FILE && \
	tar xzf $FILE && \
        cp STAR-2.7.10b/bin/Linux_x86_64_static/STAR /usr/local/bin/ && \
	rm -rf $FILE STAR-2.7.10b && \
        apk del wget

CMD ["STAR"]


