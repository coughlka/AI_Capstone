# Stage 1: Build dependencies
FROM python:3.11-slim-bullseye AS builder

# Set working directory
WORKDIR /web-app

# Install build dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir --user -r requirements.txt


# Stage 2: Runtime
FROM python:3.11-slim-bullseye

# Set working directory
WORKDIR /web-app

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy Python dependencies from builder
COPY --from=builder /root/.local /root/.local

# Make sure scripts in .local are usable
ENV PATH=/root/.local/bin:$PATH
ENV PYTHONPATH=/root/.local/lib/python3.11/site-packages

# Copy application code
COPY web-app/ .

# Copy pipeline outputs
COPY outputs/ /outputs/
ENV OUTPUTS_DIR=/outputs

# Expose port (default 8000, can be overridden with BACKEND_PORT env var)
EXPOSE 8000

# Run the application
CMD ["python", "main.py"]
