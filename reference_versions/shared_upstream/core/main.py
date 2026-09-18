from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable, Sequence

from tools.check_excel_duplicates import find_duplicates
from core.llm_parser import parse_with_llm
from core.raw_parser import convert_workbook

PROJECT_ROOT = Path(__file__).resolve().parent.parent
LOG_DIR = PROJECT_ROOT / "outputs" / "logs"
PIPELINE_LOG = LOG_DIR / "pipeline.log"
DUPLICATE_REPORT = LOG_DIR / "duplicate_report.txt"
RAW_SUFFIX = "_json_Raw"
FORMATTED_SUFFIX = "_json"
SUPPORTED_EXTENSIONS = {".xlsx", ".xlsm", ".xls"}


def setup_logging() -> logging.Logger:
    """初始化 logging 設定，確保同時寫入檔案與終端。"""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("pipeline")
    if logger.handlers:
        return logger

    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    file_handler = logging.FileHandler(PIPELINE_LOG, encoding="utf-8")
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """設定命令列參數，支援多個檔案或資料夾以及自動模式。"""
    parser = argparse.ArgumentParser(
        description="批次處理 Excel 檔案，輸出 Raw JSON 並選擇性執行 LLM Parser。",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Excel 檔案或包含 Excel 的資料夾路徑。",
    )
    parser.add_argument(
        "--auto",
        action="store_true",
        help="自動模式：略過互動詢問並直接執行所有階段。",
    )
    return parser.parse_args(argv)


def collect_excel_files(entries: Iterable[str], logger: logging.Logger) -> list[Path]:
    """蒐集輸入路徑中的 Excel 檔案，維持穩定順序並移除重複。"""
    logger.info("開始蒐集 Excel 檔案")
    collected: list[Path] = []
    seen: set[Path] = set()

    for entry in entries:
        path = Path(entry).expanduser()
        if path.is_dir():
            candidates = sorted(
                file for file in path.iterdir() if file.suffix.lower() in SUPPORTED_EXTENSIONS
            )
            for file in candidates:
                resolved = file.resolve()
                if resolved not in seen:
                    collected.append(resolved)
                    seen.add(resolved)
        else:
            if path.exists() and path.suffix.lower() in SUPPORTED_EXTENSIONS:
                resolved = path.resolve()
                if resolved not in seen:
                    collected.append(resolved)
                    seen.add(resolved)
            else:
                logger.warning("忽略無效輸入：%s", path)

    logger.info("結束蒐集 Excel 檔案，共 %d 筆", len(collected))
    return collected


def run_duplicate_check(excel_files: list[Path], logger: logging.Logger) -> bool:
    """執行重複檔案檢查，若發現重複則寫入報告並回傳 True。"""
    logger.info("開始檢查重複檔案")
    duplicates = find_duplicates(excel_files)
    if duplicates:
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        lines: list[str] = []
        for group in duplicates:
            group_paths = [str(Path(item)) for item in group]
            lines.append(", ".join(group_paths))
            logger.warning("發現重複檔案組：%s", group_paths)

        DUPLICATE_REPORT.write_text("\n".join(lines), encoding="utf-8")
        logger.error("偵測到重複檔案，詳細資訊請見 %s", DUPLICATE_REPORT)
        logger.info("結束檢查重複檔案")
        return True

    logger.info("Checked %d files, no duplicates found.", len(excel_files))
    logger.info("結束檢查重複檔案")
    return False


def ensure_formatted_dir(raw_dir: Path) -> Path:
    """建立 `_json` 目錄，儲存 LLM 解析後的整理版 JSON。"""
    if raw_dir.name.endswith(RAW_SUFFIX):
        formatted_name = raw_dir.name[: -len(RAW_SUFFIX)] + FORMATTED_SUFFIX
    else:
        formatted_name = raw_dir.name + FORMATTED_SUFFIX
    formatted_dir = raw_dir.with_name(formatted_name)
    formatted_dir.mkdir(parents=True, exist_ok=True)
    return formatted_dir


def run_raw_conversion(excel_files: Iterable[Path], logger: logging.Logger) -> dict[Path, Path]:
    """呼叫 Raw Parser 將 Excel 轉成原始 JSON，並回傳對應目錄映射。"""
    logger.info("開始 Raw JSON 轉換階段")
    mapping: dict[Path, Path] = {}
    for excel_path in excel_files:
        succeeded, errors, out_dir = convert_workbook(
            excel_path,
            sheet_names=None,
            output_dir=None,
        )
        mapping[excel_path] = out_dir
        logger.info("Raw JSON 轉換完成：%s -> %s", excel_path, out_dir)
        if succeeded:
            logger.debug("成功解析工作表：%s", ", ".join(succeeded))
        if errors:
            for sheet, err in errors:
                logger.error("解析工作表 %s 失敗：%s", sheet, err)
        print(f"[RAW] {excel_path.name} 已輸出至 {out_dir}")
    logger.info("結束 Raw JSON 轉換階段")
    return mapping


def run_llm_pipeline(
    excel_files: Iterable[Path],
    raw_mapping: dict[Path, Path],
    logger: logging.Logger,
    auto_mode: bool,
) -> None:
    """根據互動設定決定是否執行 LLM Parser。"""
    if not raw_mapping:
        logger.warning("沒有可供 LLM Parser 使用的資料，跳過此階段")
        return

    logger.info("開始 LLM Parser 階段")

    if auto_mode:
        for excel_path in excel_files:
            raw_dir = raw_mapping[excel_path]
            formatted_dir = ensure_formatted_dir(raw_dir)  # `_json` 目錄存放整理後結果
            parse_with_llm(raw_dir, formatted_dir)
            logger.info("LLM Parser 轉換完成：%s -> %s", raw_dir, formatted_dir)
            print(f"[LLM] {excel_path.name} 已輸出至 {formatted_dir}")
        logger.info("結束 LLM Parser 階段")
        return

    all_remaining = False
    skip_all = False

    for excel_path in excel_files:
        if skip_all:
            logger.info("使用者選擇跳過後續 LLM Parser，提前結束")
            break

        raw_dir = raw_mapping[excel_path]
        should_convert = all_remaining

        if not all_remaining:
            while True:
                choice = input(
                    f"{excel_path.name} - Convert to formatted JSON using LLM? (y/n/all/skip): "
                ).strip().lower()
                if choice in {"y", "n", "all", "skip"}:
                    break
                print("請輸入 y、n、all 或 skip。")

            if choice == "y":
                should_convert = True
            elif choice == "n":
                should_convert = False
            elif choice == "all":
                should_convert = True
                all_remaining = True
            elif choice == "skip":
                skip_all = True
                logger.info("使用者選擇跳過所有 LLM Parser。")
                break

        if should_convert and not skip_all:
            formatted_dir = ensure_formatted_dir(raw_dir)
            parse_with_llm(raw_dir, formatted_dir)
            logger.info("LLM Parser 轉換完成：%s -> %s", raw_dir, formatted_dir)
            print(f"[LLM] {excel_path.name} 已輸出至 {formatted_dir}")

    logger.info("結束 LLM Parser 階段")


def main(argv: Sequence[str] | None = None) -> None:
    """主程式入口：串接檢查、轉換與互動流程。"""
    args = parse_args(argv)
    logger = setup_logging()
    logger.info("Pipeline 流程啟動")

    excel_files = collect_excel_files(args.inputs, logger)
    if not excel_files:
        logger.error("未找到任何 Excel 檔案，流程終止")
        print("未找到任何 Excel 檔案可以處理。")
        return

    if run_duplicate_check(excel_files, logger):
        print("偵測到重複檔案，流程已終止。請查看 outputs/logs/duplicate_report.txt")
        return

    raw_mapping = run_raw_conversion(excel_files, logger)
    run_llm_pipeline(excel_files, raw_mapping, logger, auto_mode=args.auto)

    logger.info("Pipeline 流程結束")


if __name__ == "__main__":
    main()
