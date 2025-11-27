
# --- Builder stage: prepare wheels to avoid rebuilding C extensions in final image ---
FROM python:3.11-slim AS builder
ENV PYTHONUNBUFFERED=1
WORKDIR /app

# Install build dependencies (kept in builder only)
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential gcc \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and build wheels into /wheels
COPY requirements.txt /app/requirements.txt
RUN pip install --upgrade pip setuptools wheel \
    && pip wheel --wheel-dir=/wheels -r /app/requirements.txt || true


# --- Runtime stage: smaller image that installs from wheels when possible ---
FROM python:3.11-slim AS runtime
ENV PYTHONUNBUFFERED=1
WORKDIR /app

# Copy prebuilt wheels and requirements
COPY --from=builder /wheels /wheels
COPY requirements.txt /app/requirements.txt

# Install runtime dependencies using local wheels first (falls back to PyPI)
RUN pip install --upgrade pip \
    && pip install --no-cache-dir --no-index --find-links=/wheels -r /app/requirements.txt || pip install --no-cache-dir -r /app/requirements.txt

# Copy project files
COPY . /app

# Clean up wheels and apt lists to reduce image size
RUN rm -rf /wheels /var/lib/apt/lists/* || true

# Use a non-root user for better security
RUN addgroup --system app && adduser --system --ingroup app app || true
USER app

CMD ["python", "-c", "import os,sys; print('ModelHATPC container ready. Files:', os.listdir('.'))"]


# --- Dev target: includes developer requirements ---
FROM runtime AS dev
COPY requirements-dev.txt /app/requirements-dev.txt
RUN pip install --no-cache-dir -r /app/requirements-dev.txt || true


# --- Minimal final image (smaller) ---
FROM gcr.io/distroless/python3 AS minimal

# Copy installed Python libraries from the runtime stage.
# Path is the default for Debian-based official Python images; if using a different
# base change this path accordingly.
COPY --from=runtime /usr/local/lib/python3.11 /usr/local/lib/python3.11
COPY --from=runtime /app /app
WORKDIR /app

# Distroless does not include a shell. Use python as the entrypoint.
CMD ["python", "-c", "import os,sys; print('ModelHATPC minimal container ready. Files:', os.listdir('/app'))"]
