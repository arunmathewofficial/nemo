## Instructions:
Add docker images to gitlab container registry:
+ `docker pull --platform linux/amd64 jmackey12345/nemo-ci:latest`
+ `docker tag jmackey12345/nemo-ci:latest registry.git.nrw/jmackey1/nemo/nemo-ci:latest`
+ `docker push registry.git.nrw/jmackey1/nemo/nemo-ci:latest`

This failed on macOS because some layers are not uploaded.  So I tried this, and it worked (thanks to google ai):
+ `docker buildx build  --platform linux/amd64 --provenance=false --sbom=false -f Nemo_Dockerfile -t registry.git.nrw/jmackey1/nemo/nemo-ci:latest --push .`

## Nemo dockerfile:
```
FROM debian:13-slim

# Set non-interactive mode for apt
ENV DEBIAN_FRONTEND=noninteractive

# Update package lists and install packages
RUN apt-get update && apt-get install --no-install-recommends -y \
    clang-format-18 \
    git \
    g++ \
    libsundials-dev \
    openmpi-bin \
    cmake \
    libspdlog-dev \
    libfmt-dev \
    libboost-math-dev \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# Copy application code (if any)
# COPY . .

# Expose port (if needed)
# EXPOSE 8080

# Default command
CMD ["bash"]
```
