from pydantic import BaseModel

class PatientData(BaseModel):
    age: int
    gender: int
    glucose: int
    cholesterol: int
    hdl: int
    tch: int
    ldl: int
    bmi: int
    smoker: int
    alcohol: int
