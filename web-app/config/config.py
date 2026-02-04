"""
Configuration file for the Biomarker Cancer Project
"""

import os
from dotenv import load_dotenv

# Load .env file if it exists
load_dotenv()

# Backend port
if os.getenv("BACKEND_PORT"):
    BACKEND_PORT = os.getenv("BACKEND_PORT")
else:
    BACKEND_PORT = 8000