# Stage 1: Python dependencies and application
FROM --platform=linux/amd64 mcr.microsoft.com/devcontainers/python:3.11-bullseye as builder

# Set working directory
WORKDIR /web-app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir --user -r requirements.txt


# Stage 2: Runtime
FROM --platform=linux/amd64 mcr.microsoft.com/devcontainers/python:3.11-bullseye

# Set working directory
WORKDIR /web-app

# Install runtime dependencies (if any)
RUN apt-get update && apt-get install -y \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy Python dependencies from builder to a location accessible by any user
COPY --from=builder /root/.local /root/.local

# Make sure scripts in .local are usable
ENV PATH=/root/.local/bin:$PATH
ENV PYTHONPATH=/root/.local/lib/python3.11/site-packages:$PYTHONPATH

# Copy application code
COPY web-app/ .

# Expose port (default 8000, can be overridden with BACKEND_PORT env var)
EXPOSE 8000

# Run the application
CMD ["python", "main.py"]

