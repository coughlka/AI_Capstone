"""
Web app for the Biomarker Cancer Project
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from config import config
from api import patient_data

# Create FastAPI app
app = FastAPI(
    title="Biomarker Cancer",
    description="Web app for the Biomarker Cancer app",
    version="0.1.0"
)

app.include_router(patient_data.router)

# Add CORS middleware (allows frontend to connect from different port)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files
app.mount("/static", StaticFiles(directory="static"), name="static")


@app.get("/")
async def root():
    """Serve the UI"""
    return FileResponse("static/index.html")


@app.get("/api")
async def api_info():
    """API info endpoint"""
    return {
        "name": "Biomarker Cancer",
        "version": "0.1.0",
        "status": "running"
    }


@app.get("/health")
async def health():
    """Health check endpoint"""
    return {"status": "healthy"}


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=config.BACKEND_PORT)
