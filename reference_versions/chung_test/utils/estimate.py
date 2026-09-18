import tiktoken

def estimate_tokens(text: str, model: str = "gpt-4o-mini") -> int:
    """
    簡單估算文字 token 數（方便你送出前做預算）。
    注意：估算非 100% 精準，但相當接近。
    """
    try:
        enc = tiktoken.encoding_for_model(model)
    except Exception:
        enc = tiktoken.get_encoding("cl100k_base")
    return len(enc.encode(text or ""))
