#!/usr/bin/env python3
import csv
import re
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATE = "20260825"
VARIANTS = ("no_fusion", "fusion")


def fields(line):
    return dict(re.findall(r"([A-Za-z_][A-Za-z0-9_]*)=([^\s]+)", line))


def number(d, key, default=0.0):
    try:
        return float(d.get(key, default))
    except ValueError:
        return float(default)


def integer(d, key, default=0):
    try:
        return int(d.get(key, default))
    except ValueError:
        return int(default)


def timer(text, label):
    match = re.search(
        rf"^Elapsed time: {re.escape(label)}: total\s+([0-9.]+) sec",
        text,
        re.MULTILINE,
    )
    return float(match.group(1)) if match else None


def wall_seconds(value):
    bits = value.split(":")
    if len(bits) == 2:
        return float(bits[0]) * 60.0 + float(bits[1])
    if len(bits) == 3:
        return float(bits[0]) * 3600.0 + float(bits[1]) * 60.0 + float(bits[2])
    return float(value)


def external_wall(job_text, variant):
    match = re.search(
        rf"PAIR_VARIANT_BEGIN[^\n]*variant={variant}(.*?)"
        rf"PAIR_VARIANT_END[^\n]*variant={variant}",
        job_text,
        re.DOTALL,
    )
    if not match:
        return None
    elapsed = re.search(
        r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([0-9:.]+)",
        match.group(1),
    )
    return wall_seconds(elapsed.group(1)) if elapsed else None


def sum_key(records, key):
    return sum(number(record, key) for record in records)


def precision_records(records, precision):
    return [record for record in records if record.get("precision") == precision]


def fmt(value):
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.12g}"
    return value


run_rows = []
cg_rows = []
dagd_rows = []

for pair_index in range(1, 11):
    pair = f"pair{pair_index:02d}"
    pair_dir = ROOT / pair
    out_path = next(pair_dir.glob("run_pair.pjm.*.out"))
    stats_path = next(pair_dir.glob("run_pair.pjm.*.stats"))
    job_id = re.search(r"\.(\d+)\.out$", out_path.name).group(1)
    job_text = out_path.read_text(errors="replace")
    stats_text = stats_path.read_text(errors="replace")
    exit_match = re.search(r"^\s*EXIT CODE\s*:\s*(\d+)", stats_text, re.MULTILINE)
    exit_code = int(exit_match.group(1)) if exit_match else None
    first = "no_fusion" if pair_index % 2 else "fusion"

    for variant in VARIANTS:
        log_path = pair_dir / variant / "hmc.log"
        text = log_path.read_text(errors="replace")
        cg = [fields(line) for line in text.splitlines() if "CG_PROFILE " in line]
        dagd_fusion = [
            fields(line) for line in text.splitlines() if "DAGD_FUSION_PROFILE " in line
        ]
        dagd_detail = [
            fields(line) for line in text.splitlines() if "DAGD_DETAIL_PROFILE " in line
        ]
        if len(dagd_fusion) != len(dagd_detail):
            raise RuntimeError(f"DdagD record mismatch: {pair} {variant}")

        for record_index, record in enumerate(cg, 1):
            cg_rows.append({"pair": pair, "job_id": job_id, "variant": variant,
                            "record_index": record_index, **record})
        for record_index, (fusion_record, detail_record) in enumerate(
            zip(dagd_fusion, dagd_detail), 1
        ):
            merged = {"pair": pair, "job_id": job_id, "variant": variant,
                      "record_index": record_index}
            merged.update({f"fusion_{key}": value for key, value in fusion_record.items()})
            merged.update({f"detail_{key}": value for key, value in detail_record.items()})
            dagd_rows.append(merged)

        hdiff_match = re.search(r"H\(diff\)\s*=\s*([-+0-9.eE]+)", text)
        row = {
            "pair": pair,
            "job_id": job_id,
            "variant": variant,
            "order_position": 1 if variant == first else 2,
            "exit_code": exit_code,
            "h_diff": float(hdiff_match.group(1)) if hdiff_match else None,
            "trajectory_sec": timer(text, "each trajectory"),
            "update_sec": timer(text, "update_Nf2p1_eo"),
            "main_sec": timer(text, "Main"),
            "external_wall_sec": external_wall(job_text, variant),
            "cg_records": len(cg),
            "cg_iterations_total": sum(integer(record, "iterations") for record in cg),
            "cg_total_sec": sum_key(cg, "total_sec"),
            "cg_steps_sec": sum_key(cg, "steps_sec"),
            "cg_residual_block_sec": sum_key(cg, "step_residual_block_sec"),
            "dagd_records": len(dagd_fusion),
            "dagd_calls_total": sum(integer(record, "calls") for record in dagd_fusion),
            "dagd_total_sec": sum_key(dagd_fusion, "dagd_total_sec"),
            "dagd_d_combine_sec": sum_key(dagd_fusion, "d_combine_total_sec"),
            "dagd_ddag_combine_sec": sum_key(dagd_fusion, "ddag_combine_total_sec"),
            "dagd_combine_sec": sum_key(dagd_fusion, "combine_total_sec"),
        }
        for precision in ("float", "double"):
            pcg = precision_records(cg, precision)
            pdagd = precision_records(dagd_fusion, precision)
            prefix = "sp" if precision == "float" else "dp"
            row[f"cg_{prefix}_records"] = len(pcg)
            row[f"cg_{prefix}_iterations_total"] = sum(
                integer(record, "iterations") for record in pcg
            )
            row[f"cg_{prefix}_total_sec"] = sum_key(pcg, "total_sec")
            row[f"cg_{prefix}_residual_block_sec"] = sum_key(
                pcg, "step_residual_block_sec"
            )
            row[f"dagd_{prefix}_records"] = len(pdagd)
            row[f"dagd_{prefix}_calls_total"] = sum(
                integer(record, "calls") for record in pdagd
            )
            row[f"dagd_{prefix}_total_sec"] = sum_key(pdagd, "dagd_total_sec")
            row[f"dagd_{prefix}_combine_sec"] = sum_key(pdagd, "combine_total_sec")
        run_rows.append(row)


def write_rows(path, rows):
    columns = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(value) for key, value in row.items()})


run_path = ROOT / f"FAIR_COMBINED_RUN_RESULTS_{DATE}.csv"
cg_path = ROOT / f"CG_PROFILE_ALL_RECORDS_{DATE}.csv"
dagd_path = ROOT / f"DAGD_PROFILE_ALL_RECORDS_{DATE}.csv"
write_rows(run_path, run_rows)
write_rows(cg_path, cg_rows)
write_rows(dagd_path, dagd_rows)

metrics = [
    "trajectory_sec", "update_sec", "main_sec", "external_wall_sec",
    "cg_iterations_total", "cg_total_sec", "cg_steps_sec",
    "cg_residual_block_sec", "cg_sp_total_sec", "cg_sp_residual_block_sec",
    "cg_dp_total_sec", "cg_dp_residual_block_sec", "dagd_calls_total",
    "dagd_total_sec", "dagd_combine_sec", "dagd_sp_total_sec",
    "dagd_sp_combine_sec", "dagd_dp_total_sec", "dagd_dp_combine_sec",
    "h_diff",
]
summary_rows = []
for metric in metrics:
    no_values = [float(row[metric]) for row in run_rows if row["variant"] == "no_fusion"]
    yes_values = [float(row[metric]) for row in run_rows if row["variant"] == "fusion"]
    paired_differences = [no_values[i] - yes_values[i] for i in range(10)]
    no_mean = statistics.mean(no_values)
    yes_mean = statistics.mean(yes_values)
    summary_rows.append({
        "metric": metric,
        "no_fusion_n": len(no_values),
        "no_fusion_mean": no_mean,
        "no_fusion_stdev": statistics.stdev(no_values),
        "fusion_n": len(yes_values),
        "fusion_mean": yes_mean,
        "fusion_stdev": statistics.stdev(yes_values),
        "mean_difference_no_minus_fusion": no_mean - yes_mean,
        "reduction_percent": 100.0 * (no_mean - yes_mean) / no_mean if no_mean else 0.0,
        "paired_difference_mean": statistics.mean(paired_differences),
        "paired_difference_stdev": statistics.stdev(paired_differences),
    })

summary_path = ROOT / f"FAIR_COMBINED_SUMMARY_{DATE}.csv"
write_rows(summary_path, summary_rows)

summary = {row["metric"]: row for row in summary_rows}
report_path = ROOT / f"FAIR_COMBINED_REPORT_JA_{DATE}.txt"
with report_path.open("w", encoding="utf-8") as stream:
    stream.write("CG axpy+norm2 / DdagD 融合 公平計測10組 集約結果\n")
    stream.write("Date: 2026-08-25\n\n")
    stream.write("[正常性]\n")
    stream.write(f"実行数: {len(run_rows)}（融合なし10、融合あり10）\n")
    stream.write(f"EXIT CODE=0: {sum(row['exit_code'] == 0 for row in run_rows)}/20\n")
    stream.write(f"CG_PROFILE: {len(cg_rows)} records\n")
    stream.write(f"DAGD_PROFILE: {len(dagd_rows)} records\n")
    stream.write("timing_scheme: single_outer_block\n\n")
    stream.write("[平均値: 融合なし -> 融合あり (削減率)]\n")
    labels = [
        ("CG残差更新ブロック", "cg_residual_block_sec"),
        ("CG全体", "cg_total_sec"),
        ("DdagD combine", "dagd_combine_sec"),
        ("DdagD全体", "dagd_total_sec"),
        ("trajectory", "trajectory_sec"),
        ("update_Nf2p1_eo", "update_sec"),
        ("Main", "main_sec"),
        ("外部wall clock", "external_wall_sec"),
    ]
    for label, metric in labels:
        row = summary[metric]
        stream.write(
            f"{label}: {row['no_fusion_mean']:.9f} s -> "
            f"{row['fusion_mean']:.9f} s "
            f"({row['reduction_percent']:.6f} %), "
            f"paired diff {row['paired_difference_mean']:.9f} "
            f"+/- {row['paired_difference_stdev']:.9f} s\n"
        )
    stream.write("\n[反復回数と妥当性]\n")
    row = summary["cg_iterations_total"]
    stream.write(
        f"CG反復回数合計/実行: {row['no_fusion_mean']:.3f} -> "
        f"{row['fusion_mean']:.3f}\n"
    )
    row = summary["h_diff"]
    stream.write(
        f"H(diff): {row['no_fusion_mean']:.10f} -> {row['fusion_mean']:.10f}\n"
    )
    stream.write("\n[解釈]\n")
    stream.write(
        "CG残差更新ブロックとDdagD combineでは、対応差の標準偏差より平均短縮量が十分大きく、局所的な短縮を確認できる。\n"
    )
    stream.write(
        "一方、trajectory・Main・外部wall clockでは対応差の標準偏差が平均短縮量より大きい。今回の10組だけから全体時間が高速化したとは断定できない。\n"
    )
    stream.write("\n[注意]\n")
    stream.write(
        "fusion条件ではCG融合とDdagD融合を同時に有効化しているため、Mainやtrajectoryの差は両者を合わせた効果である。\n"
    )
    stream.write(
        "局所時間は各PROFILE項目で個別に比較できるが、全体差を片方だけへ帰属するには個別融合ビルドが必要である。\n"
    )

print(report_path.read_text(encoding="utf-8"))
