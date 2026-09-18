from fastapi import FastAPI, Response
from pydantic import BaseModel
from parsers.text_cleaner import parse_raw_to_textjson
from parsers.casecard_parser import parse_textjson_to_casecard
from schemas import TextJSON, CaseCard
from metrics import METER
from parsers.casecard_parser_v2 import parse_textjson_to_casecard_v2
from fastapi.responses import JSONResponse
import traceback
app = FastAPI(title="Case Card Builder API", version="0.1.1")

class RawIn(BaseModel):
    text: str
    source: str = ""

@app.post("/parse/text", response_model=TextJSON)
def parse_text(payload: RawIn, response: Response):
    out = parse_raw_to_textjson(payload.text, source=payload.source)
    usage = (out.meta or {}).get("usage", {})
    # 把用量寫進回應 Header（方便在前端/工具直接看到）
    response.headers["X-Token-Model"] = str(usage.get("model", ""))
    response.headers["X-Token-Input"] = str(usage.get("input_tokens", ""))
    response.headers["X-Token-Output"] = str(usage.get("output_tokens", ""))
    response.headers["X-Token-Total"] = str(usage.get("total_tokens", ""))
    return out

@app.post("/parse/casecard")
def parse_casecard(payload: TextJSON, response: Response):
    try:
        # 1) 接 payload
        tj = payload.model_dump()   # pydantic v2
        # app.logger.info(f"/parse/casecard payload keys: {list(tj.keys())}")

        # 2) 執行 parser
        out = parse_textjson_to_casecard(tj)

        # 3) 轉成 dict（不管 parser 回傳 model 還是 dict 都統一）
        if isinstance(out, BaseModel):
            out_dict = out.model_dump()
        else:
            out_dict = out

        # 4) 安全取 usage + 寫 header
        meta = (out_dict or {}).get("meta", {}) or {}
        usage = meta.get("usage", {}) or {}
        response.headers["X-Token-Model"]  = str(usage.get("model", ""))
        response.headers["X-Token-Input"]  = str(usage.get("input_tokens", ""))
        response.headers["X-Token-Output"] = str(usage.get("output_tokens", ""))
        response.headers["X-Token-Total"]  = str(usage.get("total_tokens", ""))

        return out_dict

    except Exception as e:
        # 直接把錯誤訊息與 traceback 回給你，Swagger 就看得到，而不是一個空空的 "Internal Server Error"
        return JSONResponse(
            status_code=500,
            content={"error": str(e), "trace": traceback.format_exc()}
        )

@app.post("/parse/casecard.v2")
def parse_casecard_v2(payload: TextJSON, response: Response):
    out = parse_textjson_to_casecard_v2(payload.model_dump())
    usage = (out.meta or {}).get("usage", {})
    response.headers["X-Token-Model"] = str(usage.get("model", ""))
    response.headers["X-Token-Input"] = str(usage.get("input_tokens", ""))
    response.headers["X-Token-Output"] = str(usage.get("output_tokens", ""))
    response.headers["X-Token-Total"] = str(usage.get("total_tokens", ""))
    return out

@app.get("/metrics/tokens")
def token_metrics():
    """回傳服務啟動以來的 token 統計（總量＋依模型）。"""
    return METER.snapshot()
