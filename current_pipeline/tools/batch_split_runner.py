from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List

WORK_DIR_TEMPLATE = "work_batch_{idx}"
BACKUP_DIR_NAME = "backup_outputs"
DEFAULT_TASKS = ["analyze", "mngs_big", "summarize"]


def prompt_yes_no(message: str) -> bool:
    answer = input(f"{message} [y/N]: ").strip().lower()
    return answer == "y"


def list_patient_dirs(base: Path) -> list[Path]:
    candidates: list[Path] = []
    for child in base.iterdir():
        if not child.is_dir():
            continue
        # 跳過 Raw 版本的資料夾
        if child.name.lower().endswith("_raw"):
            continue
        candidates.append(child.resolve())

    def _sort_key(path: Path) -> tuple[int, str]:
        # 依檔名中的數字排序，沒有數字就用名稱
        import re

        match = re.search(r"(\d+)", path.name)
        if match:
            return (int(match.group(1)), path.name)
        return (10**9, path.name)

    return sorted(candidates, key=_sort_key)


def split_into_batches(items: List[Path], batch_count: int = 4) -> list[list[Path]]:
    batches: list[list[Path]] = [[] for _ in range(batch_count)]
    for idx, item in enumerate(items):
        batches[idx % batch_count].append(item)
    return batches


def copy_patients_to_workdirs(batches: list[list[Path]], work_root: Path) -> list[Path]:
    work_dirs: list[Path] = []
    work_root.mkdir(parents=True, exist_ok=True)
    for idx, batch in enumerate(batches, start=1):
        work_dir = work_root / WORK_DIR_TEMPLATE.format(idx=idx)
        work_dir.mkdir(parents=True, exist_ok=True)
        for patient_dir in batch:
            dest = work_dir / patient_dir.name
            if dest.exists():
                shutil.rmtree(dest)
            shutil.copytree(patient_dir, dest)
        work_dirs.append(work_dir)
    return work_dirs


def build_task_command(base_cmd: list[str], tasks: list[str], patient_dir: Path) -> list[str]:
    # llm_scheduler 需要位置參數 inputs，這裡直接指定 patient_dir 為 inputs
    # 另外附上 --patients 以便過濾。
    return base_cmd + [str(patient_dir), "--patients", patient_dir.name, "--tasks", *tasks]


def run_command(cmd: list[str]) -> int:
    print(f"Running: {' '.join(cmd)}")
    proc = subprocess.run(cmd, text=True)
    return proc.returncode


def backup_and_merge_outputs(src_patient: Path, dest_patient: Path) -> None:
    backup_root = dest_patient / BACKUP_DIR_NAME
    backup_root.mkdir(parents=True, exist_ok=True)

    for name in ("agent_outputs", "summary_outputs"):
        src = src_patient / name
        if src.exists():
            dst = dest_patient / name
            if dst.exists():
                backup_dst = backup_root / name
                if backup_dst.exists():
                    shutil.rmtree(backup_dst)
                shutil.copytree(dst, backup_dst)
                shutil.rmtree(dst)
            shutil.copytree(src, dst)

    # copy loose JSON files in root (e.g., *_mNGS_max_agent.json)
    for file in src_patient.glob("*.json"):
        target = dest_patient / file.name
        if target.exists():
            shutil.copy2(target, backup_root / file.name)
        shutil.copy2(file, target)


def remove_work_root(work_root: Path) -> None:
    if work_root.exists():
        shutil.rmtree(work_root, ignore_errors=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="將病人資料分成 4 批跑任務，再合併結果（會先備份原輸出）。"
    )
    parser.add_argument(
        "--base",
        required=True,
        type=Path,
        help="病人來源父層目錄（會掃描底下子資料夾）。",
    )
    parser.add_argument(
        "--work-root",
        type=Path,
        default=Path("batch_work"),
        help="暫存工作區根目錄（預設 batch_work）。",
    )
    parser.add_argument(
        "--tasks",
        nargs="+",
        default=DEFAULT_TASKS,
        help="要執行的任務列表，預設 analyze mngs_big summarize。",
    )
    parser.add_argument(
        "--cmd",
        nargs="+",
        default=["python", "-m", "tools.llm_scheduler"],
        help="執行任務的基礎命令（會附加 --tasks/--patients/路徑）。",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=1,
        help="傳給 llm_scheduler 的併發病人數，預設 1（在每批內序列執行）。",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="傳給 llm_scheduler 的跳過已存在輸出旗標。",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    patients = list_patient_dirs(args.base)
    if not patients:
        print("找不到任何病人資料夾，請確認 --base 路徑。")
        sys.exit(1)

    batches = split_into_batches(patients, batch_count=4)
    work_dirs = copy_patients_to_workdirs(batches, args.work_root)

    print("已建立工作目錄：")
    for wd in work_dirs:
        print(f" - {wd}")
    if not prompt_yes_no("是否繼續執行任務？"):
        print("已中止，未執行任何任務，並清除暫存目錄。")
        remove_work_root(args.work_root)
        sys.exit(0)

    for work_dir, batch in zip(work_dirs, batches, strict=True):
        for patient in batch:
            patient_copy = work_dir / patient.name
            cmd = args.cmd + ["--max-workers", str(args.max_workers)]
            if args.skip_existing:
                cmd.append("--skip-existing")
            full_cmd = build_task_command(cmd, args.tasks, patient_copy)
            code = run_command(full_cmd)
            if code != 0:
                print(f"命令失敗，退出碼 {code}：{' '.join(full_cmd)}")
                continue
            backup_and_merge_outputs(patient_copy, patient)

    # 清理暫存工作區
    remove_work_root(args.work_root)


if __name__ == "__main__":
    main()
