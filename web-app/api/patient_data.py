from fastapi import APIRouter, HTTPException
import random # only for demo purposes

router = APIRouter(prefix="/patient-data", tags=["Patient Data"])

@router.post("/predict", response_model=dict)
async def predict_patient_data(patient_data: dict):
    """
    Predicts the patient's risk of developing cancer based on the patient's biomarker data.
    """

    try:
        result = random.randint(0, 100)
        return {"result": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
