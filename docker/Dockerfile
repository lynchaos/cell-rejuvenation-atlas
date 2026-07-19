# Rejuvenation Atlas pipeline container
# Digest-pinned base for reproducibility; update the digest deliberately.
FROM mambaorg/micromamba:1.5.8

COPY environment.yml /tmp/environment.yml
RUN micromamba create -y -n atlas -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH=/opt/conda/envs/atlas/bin:$PATH \
    MPLBACKEND=Agg

WORKDIR /opt/atlas
COPY src/ ./src/
COPY main.nf nextflow.config ./
COPY modules/ ./modules/
COPY conf/ ./conf/

# Nextflow supplies its own command; do not override ENTRYPOINT.
ENTRYPOINT []
