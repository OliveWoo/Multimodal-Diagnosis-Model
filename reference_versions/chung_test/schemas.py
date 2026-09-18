from typing import List, Optional, Union, Dict, Any
from pydantic import BaseModel, Field, field_validator

# ---------- 共用設定 ----------
class _Base(BaseModel):
    # 忽略 LLM 多塞的額外欄位，避免驗證直接噴掉
    model_config = {"extra": "ignore"}

# ---------- Layer 1: text_JSON ----------
class TextOnlySection(_Base):
    text: str = ""

class TextJSONSections(_Base):
    laboratory: TextOnlySection = TextOnlySection()
    imaging:   TextOnlySection = TextOnlySection()
    pft:       TextOnlySection = TextOnlySection()
    pathology: TextOnlySection = TextOnlySection()

class TextJSON(_Base):
    meta: Dict[str, Any] = Field(default_factory=dict)
    sections: TextJSONSections = TextJSONSections()

def empty_textjson(source: str = "") -> "TextJSON":
    return TextJSON(
        meta={"source": source},
        sections=TextJSONSections(),  # 皆有預設值
    )

# ---------- Layer 2: Case Card ----------
class LabItem(_Base):
    item: str
    unit: Optional[str] = None
    value: Union[float, str]
    reference: Optional[str] = None
    notes: Optional[str] = None
    source_span: Optional[str] = None
    confidence: Optional[float] = None

class LabBatch(_Base):
    date: str  # YYYY/MM or YYYY/MM/DD
    items: List[LabItem]

class LaboratorySection(_Base):
    batches: List[LabBatch] = Field(default_factory=list)

class PFTRow(_Base):
    date: str
    FVC_L: Optional[float] = None
    FVC_pct: Optional[float] = None
    FEV1_L: Optional[float] = None
    FEV1_pct: Optional[float] = None
    TLC_L: Optional[float] = None
    TLC_pct: Optional[float] = None
    DLCO_pct: Optional[float] = None
    source_span: Optional[str] = None

class PFTSection(_Base):
    rows: List[PFTRow] = Field(default_factory=list)

class ImagingSection(_Base):
    text: str = ""
    flags: List[str] = Field(default_factory=list)

    @field_validator("flags", mode="before")
    @classmethod
    def _coerce_flags(cls, v):
        if v is None:
            return []
        if isinstance(v, str):
            return [v]
        return list(v)

class PathologySection(_Base):
    text: str = ""

class CaseCardSections(_Base):
    laboratory: LaboratorySection = LaboratorySection()
    pft: PFTSection = PFTSection()
    imaging: ImagingSection = ImagingSection()
    pathology: PathologySection = PathologySection()

class CaseCardSummary(_Base):
    key_labs: Dict[str, Any] = Field(default_factory=dict)
    imaging_flags: List[str] = Field(default_factory=list)
    suspected_dx: List[str] = Field(default_factory=list)

    @field_validator("imaging_flags", "suspected_dx", mode="before")
    @classmethod
    def _coerce_list(cls, v):
        if v is None:
            return []
        if isinstance(v, str):
            return [v]  # 單一字串 → list
        return list(v)

class CaseCard(_Base):
    summary: CaseCardSummary = CaseCardSummary()
    sections: CaseCardSections = CaseCardSections()
    meta: Dict[str, Any] = Field(default_factory=lambda: {"errors": []})

def empty_casecard() -> "CaseCard":
    return CaseCard()  # 全部都有預設值
