# AWS setup (one-time, ~30 minutes)

The container image is built and hosted on GitHub Container Registry (GHCR) via
GitHub Actions (see `.github/workflows/docker.yml`). AWS only provides compute
(S3 + Batch), so you do not need ECR.

1. **S3**: create a bucket for the Nextflow work dir and results.
   ```bash
   aws s3 mb s3://YOUR-BUCKET --region us-east-1
   ```
2. **Batch**: create two compute environments + job queues:
   * `rejuvenation-atlas-queue` — spot, c5/m5 families (CPU jobs)
   * `rejuvenation-atlas-gpu-queue` — g4dn.xlarge (the `gpu`-labeled scVI job)
3. **IAM**:
   * The Batch instance role needs read/write access to the S3 bucket.
   * The instance role (or a secret injected into the job) needs a GitHub token
     with `read:packages` scope so it can pull the GHCR image.
   * Your user needs `batch:SubmitJob`.
4. Edit `conf/awsbatch.config` (replace `YOUR-BUCKET` and `ghcr.io/...` if you
   forked the repo), then run:

   ```bash
   nextflow run main.nf -profile awsbatch --module all
   ```

Costs: a full run is designed to fit in a few USD of spot compute; the scVI job
is the largest (~1 GPU-hour).

> **Note for Altos Labs demo**: start with `--module rejuvenation_clock` to
> reproduce the Gill et al. 2022 epigenetic-age reversal result quickly and
> cheaply on CPU-only spot instances.
