# EvoNexus runtime image.
#
# The `docker` and `singularity` profiles both point at a single image,
# because RUN_EVOLUTION shells out to MAFFT, IQ-TREE 2 and HyPhy and must find
# all three on PATH alongside the Python package.
#
#   docker build -t ghcr.io/aglucaci/evonexus:0.1.0 .
#   docker push  ghcr.io/aglucaci/evonexus:0.1.0
#
# Until that image is pushed to a registry you can reach, run the container
# profiles against a local build or override the tag:
#
#   nextflow run . -profile full,singularity --container my/evonexus:dev
#
FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="EvoNexus" \
      org.opencontainers.image.description="Evolutionary and ecological networks of viral host sharing in bats" \
      org.opencontainers.image.source="https://github.com/aglucaci/EvoNexus" \
      org.opencontainers.image.licenses="MIT"

USER root
# `ps` is required by Nextflow to collect task runtime metrics; procps is not
# in the micromamba base image.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps ca-certificates \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes

# The pipeline mounts the repository and sets PYTHONPATH itself, so the
# package is not baked in: the image carries the environment, the run carries
# the code. This keeps the image valid across code changes.
ENV PATH=/opt/conda/bin:$PATH \
    PYTHONNOUSERSITE=1 \
    MPLBACKEND=Agg

# A login shell sources /etc/profile, which resets PATH and loses the conda
# prefix -- `bash -lc python` then fails even though ENV PATH is correct.
# Nextflow runs task scripts with a non-login shell, but Singularity and
# interactive `docker run -it` do not, so the PATH is pinned here as well.
USER root
RUN printf 'export PATH=/opt/conda/bin:$PATH\n' > /etc/profile.d/00-conda.sh \
 && chmod 0644 /etc/profile.d/00-conda.sh
USER $MAMBA_USER

SHELL ["/bin/bash", "-c"]
