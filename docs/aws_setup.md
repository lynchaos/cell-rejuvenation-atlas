# AWS setup (one-time, ~30 minutes)

1. **S3**: `aws s3 mb s3://YOUR-BUCKET` — used for the Nextflow work dir and results.
2. **ECR + container**:
   ```bash
   aws ecr create-repository --repository-name rejuvenation-atlas
   docker build -t rejuvenation-atlas -f docker/Dockerfile .
   docker tag rejuvenation-atlas YOUR-ACCOUNT.dkr.ecr.us-east-1.amazonaws.com/rejuvenation-atlas:latest
   aws ecr get-login-password | docker login --username AWS --password-stdin YOUR-ACCOUNT.dkr.ecr.us-east-1.amazonaws.com
   docker push YOUR-ACCOUNT.dkr.ecr.us-east-1.amazonaws.com/rejuvenation-atlas:latest
   ```
3. **Batch**: create two compute environments + job queues:
   * `rejuvenation-atlas-queue` — spot, c5/m5 families (CPU jobs)
   * `rejuvenation-atlas-gpu-queue` — g4dn.xlarge (the `gpu`-labeled scVI job)
4. **IAM**: the instance role needs S3 + ECR read/write; your user needs `batch:SubmitJob`.
5. Edit `conf/awsbatch.config` (bucket, queues, container URI), then:

   ```bash
   nextflow run main.nf -profile awsbatch --module all
   ```

Costs: a full run is designed to fit in a few USD of spot compute; the scVI job dominates (~1 GPU-hour).
