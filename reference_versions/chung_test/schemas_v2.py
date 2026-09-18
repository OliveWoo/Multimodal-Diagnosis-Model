# schemas_v2.py
from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field

class KeyLab(BaseModel):
    value: str
    date: Optional[str] = None
    confidence: Optional[float] = None

class Outline(BaseModel):
    case_uid: Optional[str] = None
    version: int = 2
    patient_id: Optional[str] = None
    encounter_date: Optional[str] = None  # YYYY-MM-DD or YYYY/MM
    age: Optional[int] = None
    sex: Optional[str] = None
    chief_complaint: Optional[str] = None
    suspected_dx: List[str] = Field(default_factory=list)
    imaging_flags: List[str] = Field(default_factory=list)
    key_labs: Dict[str, KeyLab] = Field(default_factory=dict)
    pft_trend: Optional[str] = None
    notes: Optional[str] = None
    sources: List[str] = Field(default_factory=list)

class LabItem(BaseModel):
    item: str
    unit: Optional[str] = None
    value: str
    reference: Optional[str] = None
    notes: Optional[str] = None
    source_span: Optional[str] = None
    confidence: Optional[float] = None

class LabBatch(BaseModel):
    date: str
    items: List[LabItem]

class PFT(BaseModel):
    date: str
    FVC_L: Optional[float] = None
    FVC_pct: Optional[float] = None
    FEV1_L: Optional[float] = None
    FEV1_pct: Optional[float] = None
    TLC_L: Optional[float] = None
    TLC_pct: Optional[float] = None
    DLCO_pct: Optional[float] = None
    source_span: Optional[str] = None

class ImagingSection(BaseModel):
    text: str = ""
    flags: List[str] = Field(default_factory=list)

class PathologySection(BaseModel):
    text: str = ""

class Sections(BaseModel):
    laboratory: Dict[str, List[LabBatch]] | Dict[str, Any] | List[LabBatch] | Dict[str, List] = Field(default_factory=dict)
    pft: List[PFT] = Field(default_factory=list)
    imaging: ImagingSection = Field(default_factory=ImagingSection)
    pathology: PathologySection = Field(default_factory=PathologySection)

class CaseCardV2(BaseModel):
    outline: Outline
    sections: Sections
    meta: Dict[str, Any] = Field(default_factory=dict)
