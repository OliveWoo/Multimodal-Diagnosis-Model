from typing import Dict
from threading import Lock
from collections import defaultdict

class TokenMeter:
    """
    簡單的 Token 計數器：可累加 input/output/total，並依模型彙整。
    """
    def __init__(self):
        self.lock = Lock()
        self.totals = defaultdict(int)              # keys: input, output, total
        self.by_model = defaultdict(lambda: defaultdict(int))

    def add(self, model: str, input_tokens: int = 0, output_tokens: int = 0):
        with self.lock:
            total = (input_tokens or 0) + (output_tokens or 0)
            self.totals["input"] += input_tokens or 0
            self.totals["output"] += output_tokens or 0
            self.totals["total"] += total
            self.by_model[model]["input"] += input_tokens or 0
            self.by_model[model]["output"] += output_tokens or 0
            self.by_model[model]["total"] += total

    def snapshot(self) -> Dict:
        with self.lock:
            return {
                "totals": dict(self.totals),
                "by_model": {m: dict(v) for m, v in self.by_model.items()}
            }

METER = TokenMeter()
