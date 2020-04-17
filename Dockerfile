## Docker file
FROM alpine
MAINTAINER @Toyo-Daichi

# RUN: when build
RUN set -x && \
    apk update && \
    apk add --no-cache gfortran musl-dev

# CMD: when RUN
CMD ["/bin/sh"]

