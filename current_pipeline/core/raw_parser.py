"""將 Excel 工作簿拆解為各工作表的 Raw JSON 檔案。"""

from __future__ import annotations

import argparse
import importlib
import logging
import pkgutil
from pathlib import Path
from types import ModuleType
from typing import Iterable, Sequence

import parsers
from utils import sanitize_filename

DEFAULT_RAW_EXCEL_DIR = Path("raw_excel")
RAW_SUFFIX = "_json_Raw"

logger = logging.getLogger(__name__)


def _load_parsers() -> dict[str, ModuleType]:
    """載入 parsers 套件底下所有具備 SHEET_NAME 與 parse 的模組。"""
    discovered: dict[str, ModuleType] = {}
    for _loader, module_name, _is_pkg in pkgutil.iter_modules(parsers.__path__):
        module = importlib.import_module(f"parsers.{module_name}")
        if hasattr(module, "SHEET_NAME") and hasattr(module, "parse"):
            sheet_name = getattr(module, "SHEET_NAME")
            discovered[sheet_name] = module
            logger.debug("已載入 parser：%s 對應 %s", module_name, sheet_name)
    return discovered


def _default_output_dir(filepath: Path) -> Path:
    """根據 Excel 檔名產生預設的 `_json_Raw` 資料夾。"""
    sanitized = sanitize_filename(Path(filepath).stem)
    return Path(filepath).with_name(f"{sanitized}{RAW_SUFFIX}")


def convert_workbook(
    filepath: Path | str,
    sheet_names: Iterable[str] | None = None,
    output_dir: Path | None = None,
) -> tuple[list[str], list[tuple[str, Exception | str]], Path]:
    """
    將 Excel 工作簿轉為多個 Raw JSON。

    參數：
        filepath：來源 Excel 路徑。
        sheet_names：要處理的工作表名稱；預設為全部支援的 parser。
        output_dir：Raw JSON 的輸出目錄；預設使用 `_json_Raw` 後綴。

    回傳：
        (成功清單, 錯誤清單, 實際輸出資料夾)
    """
    workbook_path = Path(filepath).resolve()
    if not workbook_path.exists():
        raise FileNotFoundError(f"Excel 檔案不存在：{workbook_path}")

    target_dir = Path(output_dir) if output_dir else _default_output_dir(workbook_path)
    target_dir.mkdir(parents=True, exist_ok=True)

    parsers_map = _load_parsers()
    targets = list(sheet_names) if sheet_names else list(parsers_map.keys())

    succeeded: list[str] = []
    errors: list[tuple[str, Exception | str]] = []

    for sheet in targets:
        parser_module = parsers_map.get(sheet)
        if not parser_module:
            logger.warning("找不到對應 parser：%s", sheet)
            errors.append((sheet, "No parser registered for this sheet"))
            continue
        try:
            parser_module.parse(workbook_path, target_dir)
            succeeded.append(sheet)
            logger.debug("完成 Raw 轉換：%s -> %s", sheet, target_dir)
        except Exception as exc:  # pylint: disable=broad-exception-caught
            logger.exception("轉換工作表 %s 時發生錯誤", sheet, exc_info=exc)
            errors.append((sheet, exc))

    return succeeded, errors, target_dir


def _iter_excel_files(path: Path) -> Iterable[Path]:
    """列出目錄下的 Excel 檔案。"""
    if path.is_file():
        return [path]
    if path.is_dir():
        return sorted(
            file for file in path.iterdir() if file.suffix.lower() in {".xlsx", ".xlsm", ".xls"}
        )
    return []


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """提供簡易 CLI 介面供單獨轉換測試使用。"""
    parser = argparse.ArgumentParser(description="將 Excel 轉為 Raw JSON 檔案。")
    parser.add_argument(
        "target",
        nargs="?",
        help="Excel 檔案或資料夾；未提供時預設掃描 raw_excel/。",
    )
    parser.add_argument(
        "--sheets",
        nargs="*",
        default=None,
        help="指定要處理的工作表名稱（預設處理所有支援的 parser）。",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="指定輸出資料夾；未提供時使用 `_json_Raw` 後綴。",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    """CLI 入口：可用於手動測試 Raw 轉換。"""
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    if args.target:
        targets = _iter_excel_files(Path(args.target))
    else:
        targets = _iter_excel_files(DEFAULT_RAW_EXCEL_DIR)

    if not targets:
        print("找不到可處理的 Excel 檔案。")
        return

    for workbook in targets:
        try:
            succeeded, errors, outdir = convert_workbook(
                workbook,
                sheet_names=args.sheets,
                output_dir=args.output,
            )
        except FileNotFoundError as exc:
            logging.error(str(exc))
            continue

        print(f"來源：{workbook}")
        print(f"輸出目錄：{outdir}")
        print(f"成功工作表：{succeeded}")
        if errors:
            print("失敗工作表：")
            for sheet, err in errors:
                print(f"  - {sheet}: {err}")


if __name__ == "__main__":
    main()
