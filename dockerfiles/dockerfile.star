FROM alpine:3.16

LABEL tool="STAR"

LABEL maintainer="patrick.barth@computational.bio.uni-giessen.de"

RUN apk update && apk add bash wget && \
	wget -q https://github.com/alexdobin/STAR/releases/download/2.7.10b_alpha_220111/STAR_2.7.10b_alpha_230111_Linux_x86_64_static.zip && \
	unzip STAR_2.7.10b_alpha_230111_Linux_x86_64_static.zip && \
        mv STAR /usr/local/bin/ && \
	rm STAR_2.7.10b_alpha_230111_Linux_x86_64_static.zip && \
        apk del wget

CMD ["STAR"]


