# Rejuvenation Atlas pipeline container
# Digest-pinned base for reproducibility; update the digest deliberately.
FROM mambaorg/micromamba:1.5.8

COPY environment.yml /tmp/environment.yml
# Cache-bust token: increment when environment.yml dependencies change.
ARG CONDA_CACHE_BUST=1
RUN micromamba create -y -n atlas -f /tmp/environment.yml && \
    micromamba clean --all --yes && \
    aws --version

ENV PATH=/opt/conda/envs/atlas/bin:$PATH \
    MPLBACKEND=Agg

WORKDIR /opt/atlas
COPY src/ ./src/
COPY main.nf nextflow.config ./
COPY modules/ ./modules/
COPY conf/ ./conf/

# Nextflow supplies its own command; do not override ENTRYPOINT.
ENTRYPOINT []
