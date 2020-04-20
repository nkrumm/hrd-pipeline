FROM "sequenza/sequenza:3.0.0"

WORKDIR /
USER root
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update \
	&& (apt-get install -y procps) \
	&& apt-get clean && apt-get purge \
	&& rm -rf /var/lib/apt/lists/* /tmp/*



