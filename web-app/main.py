"""
CRC Biomarker Evidence Browser - Web Application

Browse and explore colorectal cancer biomarker candidates from the
AI_Capstone pipeline outputs.
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse

from api import biomarkers
from config import config

# Create FastAPI app
app = FastAPI(
    title="CRC Biomarker Evidence Browser",
    description="Browse and explore colorectal cancer biomarker candidates",
    version="1.0.0"
)

# Include biomarker API routes
app.include_router(biomarkers.router)

# CORS middleware - restricted to specific origins for security
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:8000",
        "http://localhost:8080",
        "http://127.0.0.1:8000",
        "http://127.0.0.1:8080",
    ],
    allow_credentials=True,
    allow_methods=["GET", "POST"],  # Read-only API (POST only for reload)
    allow_headers=["Content-Type"],
)

# Mount static files
app.mount("/static", StaticFiles(directory="static"), name="static")


@app.get("/")
async def root():
    """Serve the Evidence Browser UI."""
    return FileResponse("static/index.html")


@app.get("/health")
async def health():
    """Health check endpoint."""
    return {"status": "healthy"}


if __name__ == "__main__":
    import uvicorn

    port = int(config.BACKEND_PORT) if config.BACKEND_PORT else 8000
    uvicorn.run(app, host="0.0.0.0", port=port)
