# Rejuvenation Atlas pipeline container
# Digest-pinned base for reproducibility; update the digest deliberately.
FROM mambaorg/micromamba:1.5.8

# Install AWS CLI v2 (required by Nextflow for S3/Batch I/O inside tasks)
USER root
RUN apt-get update -qq && apt-get install -y -qq curl unzip && \
    curl -fsSL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscliv2.zip && \
    unzip -q /tmp/awscliv2.zip -d /tmp && \
    /tmp/aws/install --update && \
    rm -rf /tmp/aws /tmp/awscliv2.zip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

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

ENTRYPOINT ["python"]
