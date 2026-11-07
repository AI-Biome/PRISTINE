# PRISTINE container: dependencies managed by pixi (see pixi.toml / pixi.lock).
FROM ghcr.io/prefix-dev/pixi:0.81.0-noble AS build

WORKDIR /opt/pristine
COPY pixi.toml pixi.lock pyproject.toml LICENSE README.md ./
COPY src ./src
RUN pixi install --locked \
 && pixi shell-hook -e default > /shell-hook.sh \
 && echo 'exec "$@"' >> /shell-hook.sh

FROM ubuntu:24.04 AS production
COPY --from=build /opt/pristine /opt/pristine
COPY --from=build /shell-hook.sh /shell-hook.sh
WORKDIR /data
# Mount your project directory (containing config.yaml) at /data:
#   docker run --rm -v "$(pwd):/data" pristine --config config.yaml
ENTRYPOINT ["/bin/bash", "/shell-hook.sh", "pristine"]
CMD ["--config", "config.yaml"]
