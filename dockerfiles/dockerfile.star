FROM alpine:3.17

LABEL tool="STAR"

LABEL maintainer="patrick.barth@computational.bio.uni-giessen.de"

RUN apk update && apk add bash wget && \
	wget -q https://github.com/alexdobin/STAR/releases/download/2.7.11a/STAR_2.7.11a.zip && \
	unzip STAR_2.7.11a.zip && \
    mv STAR_2.7.11a/Linux_x86_64_static/STAR /usr/local/bin/ && \
	rm STAR_2.7.11a.zip && \
    apk del wget

CMD ["STAR"]


