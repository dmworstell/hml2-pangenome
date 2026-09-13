#!/usr/bin/env python3
"""Render a small, dependency-free PDF from a samtools depth table.

The program is intentionally a static source file.  A run calls it from its
immutable source location and writes to a run-scoped output path; no plotting
program is generated or removed at run time.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable, Iterator


METRIC_LABELS = {
    "raw": "Raw depth",
    "mapq10": "MAPQ >= 10 depth",
}

INTERPRETATION_LINES = (
    "Pooled non-overlapping outer target-flank baseline (descriptive only)",
    "Each assembly's own-flank depth is a one-copy scale on the combined reference.",
    "Body / own-flank normalized depth is descriptive; it is not copy number.",
    "Read-inferred copy number: not estimated",
    "Coverage is diagnostic / assembly support only.",
)


class PlotInputError(ValueError):
    """Raised when an input cannot support an honest coverage plot."""


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def _depth_rows(path: Path) -> Iterator[tuple[str, int, float]]:
    """Yield (contig, one-based position, depth) with strict parsing."""

    with _open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise PlotInputError(
                    f"{path}:{line_number}: expected at least three tab-separated fields"
                )
            contig = fields[0]
            if not contig or any(ch.isspace() for ch in contig):
                raise PlotInputError(f"{path}:{line_number}: invalid contig name")
            try:
                position = int(fields[1])
                depth = float(fields[2])
            except ValueError as exc:
                raise PlotInputError(
                    f"{path}:{line_number}: position and depth must be numeric"
                ) from exc
            if position < 1:
                raise PlotInputError(f"{path}:{line_number}: position must be positive")
            if not math.isfinite(depth) or depth < 0:
                raise PlotInputError(
                    f"{path}:{line_number}: depth must be finite and nonnegative"
                )
            yield contig, position, depth


def _scan_depth(path: Path) -> tuple[list[str], dict[str, dict[str, float]], int]:
    if not path.is_file():
        raise PlotInputError(f"depth table does not exist: {path}")
    order: list[str] = []
    stats: dict[str, dict[str, float]] = {}
    closed: set[str] = set()
    previous_contig: str | None = None
    previous_position: int | None = None
    row_count = 0

    for contig, position, depth in _depth_rows(path):
        row_count += 1
        if contig != previous_contig:
            if previous_contig is not None:
                closed.add(previous_contig)
            if contig in closed:
                raise PlotInputError(
                    f"contig {contig!r} appears in more than one depth-table block"
                )
            order.append(contig)
            stats[contig] = {
                "start": float(position),
                "end": float(position),
                "max_depth": depth,
            }
            previous_contig = contig
            previous_position = position
            continue

        assert previous_position is not None
        if position != previous_position + 1:
            raise PlotInputError(
                f"contig {contig!r} is not contiguous at positions "
                f"{previous_position} and {position}"
            )
        previous_position = position
        stats[contig]["end"] = float(position)
        stats[contig]["max_depth"] = max(stats[contig]["max_depth"], depth)

    if row_count == 0:
        raise PlotInputError(f"depth table has no data rows: {path}")
    return order, stats, row_count


def _binned_depth(
    path: Path,
    order: list[str],
    stats: dict[str, dict[str, float]],
    maximum_bins: int = 620,
) -> dict[str, list[float]]:
    spans = {
        contig: int(stats[contig]["end"] - stats[contig]["start"] + 1)
        for contig in order
    }
    total_span = sum(spans.values())
    bin_counts = {
        contig: max(2, min(spans[contig], round(maximum_bins * spans[contig] / total_span)))
        for contig in order
    }
    sums = {contig: [0.0] * bin_counts[contig] for contig in order}
    counts = {contig: [0] * bin_counts[contig] for contig in order}

    for contig, position, depth in _depth_rows(path):
        start = int(stats[contig]["start"])
        span = spans[contig]
        index = min(bin_counts[contig] - 1, (position - start) * bin_counts[contig] // span)
        sums[contig][index] += depth
        counts[contig][index] += 1

    return {
        contig: [
            sums[contig][index] / counts[contig][index]
            if counts[contig][index]
            else 0.0
            for index in range(bin_counts[contig])
        ]
        for contig in order
    }


def _load_summary(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    try:
        with path.open("r", encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise PlotInputError(f"could not read summary JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise PlotInputError("summary JSON must contain an object")
    estimate = value.get("read_inferred_copy_number")
    if estimate is not None and estimate != "not_estimated":
        raise PlotInputError(
            "read_inferred_copy_number must be the literal string 'not_estimated'"
        )
    warnings = value.get("warnings", [])
    if not isinstance(warnings, list) or not all(isinstance(item, str) for item in warnings):
        raise PlotInputError("summary warnings must be a list of strings")
    return value


def _pdf_string(value: str) -> str:
    safe = value.encode("latin-1", "replace").decode("latin-1")
    return safe.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def _text_command(x: float, y: float, size: float, value: str) -> str:
    return f"BT /F1 {size:.1f} Tf {x:.2f} {y:.2f} Td ({_pdf_string(value)}) Tj ET"


def _short(value: Any, limit: int = 104) -> str:
    text = str(value).replace("\r", " ").replace("\n", " ")
    return text if len(text) <= limit else text[: limit - 3] + "..."


def _summary_lines(summary: dict[str, Any], metric: str) -> list[str]:
    lines: list[str] = []
    pooled_key = f"{metric}_zero_inclusive_median_depth"
    pooled_record = summary.get("pooled_outer_flank_baselines", {})
    pooled = pooled_record.get(pooled_key) if isinstance(pooled_record, dict) else None
    if isinstance(pooled, (int, float)) and math.isfinite(float(pooled)):
        lines.append(
            f"{METRIC_LABELS[metric]} pooled outer target-flank baseline: "
            f"{float(pooled):.4g} depth"
        )

    assemblies = summary.get("assembly_summaries", [])
    if isinstance(assemblies, list):
        for record in assemblies[:2]:
            if not isinstance(record, dict):
                continue
            assembly = _short(record.get("assembly_id", "assembly"), 28)
            metric_record = record.get(metric, {})
            if not isinstance(metric_record, dict):
                metric_record = {}
            own = metric_record.get("own_outer_flank_zero_inclusive_median_depth")
            ratio = metric_record.get("body_own_flank_normalized_depth")
            parts = [assembly]
            if isinstance(own, (int, float)) and math.isfinite(float(own)):
                parts.append(f"own-flank scale={float(own):.4g}")
            if isinstance(ratio, (int, float)) and math.isfinite(float(ratio)):
                parts.append(f"body/own-flank={float(ratio):.4g}")
            lines.append("; ".join(parts))

    for warning in summary.get("warnings", [])[:3]:
        lines.append(f"Warning: {_short(warning, 92)}")
    if not summary.get("warnings"):
        lines.append(
            "Warning context: inspect low-MAPQ and ambiguous-alignment fractions with this plot."
        )
    return lines


def _content_stream(
    order: list[str],
    stats: dict[str, dict[str, float]],
    binned: dict[str, list[float]],
    metric: str,
    metric_label: str,
    title: str,
    summary: dict[str, Any],
) -> bytes:
    commands: list[str] = []
    commands.append(_text_command(48, 572, 16, _short(title, 90)))
    commands.append(_text_command(48, 550, 12, metric_label))

    left, bottom, width, height = 58.0, 214.0, 676.0, 300.0
    commands.extend(
        [
            "0.92 0.92 0.92 RG 0.5 w",
            f"{left:.2f} {bottom:.2f} {width:.2f} {height:.2f} re S",
        ]
    )
    maximum_depth = max(float(stats[contig]["max_depth"]) for contig in order)
    y_scale_max = maximum_depth if maximum_depth > 0 else 1.0
    commands.append(_text_command(4, bottom + height - 4, 8, f"{maximum_depth:.4g}"))
    commands.append(_text_command(30, bottom - 3, 8, "0"))

    spans = {
        contig: int(stats[contig]["end"] - stats[contig]["start"] + 1)
        for contig in order
    }
    total_span = sum(spans.values())
    colors = ((0.08, 0.31, 0.62), (0.76, 0.24, 0.16), (0.16, 0.55, 0.30))
    cursor = left
    for ordinal, contig in enumerate(order):
        segment_width = width * spans[contig] / total_span
        values = binned[contig]
        red, green, blue = colors[ordinal % len(colors)]
        commands.append(f"{red:.2f} {green:.2f} {blue:.2f} RG 0.8 w")
        for index, value in enumerate(values):
            x = cursor + segment_width * index / max(1, len(values) - 1)
            y = bottom + height * value / y_scale_max
            commands.append(f"{x:.2f} {y:.2f} {'m' if index == 0 else 'l'}")
        commands.append("S")
        commands.append(_text_command(cursor + 2, bottom - 14, 7, _short(contig, 42)))
        commands.append(
            _text_command(
                cursor + 2,
                bottom - 25,
                6.5,
                f"{int(stats[contig]['start'])}-{int(stats[contig]['end'])}",
            )
        )
        cursor += segment_width
        if ordinal + 1 < len(order):
            commands.append(f"0.78 0.78 0.78 RG 0.4 w {cursor:.2f} {bottom:.2f} m {cursor:.2f} {bottom + height:.2f} l S")

    y = 165.0
    for line in INTERPRETATION_LINES:
        commands.append(_text_command(48, y, 8.2, line))
        y -= 12.0
    for line in _summary_lines(summary, metric):
        if y < 24:
            break
        commands.append(_text_command(48, y, 7.8, line))
        y -= 11.0
    return ("\n".join(commands) + "\n").encode("latin-1", "replace")


def _pdf_document(content: bytes) -> bytes:
    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        (
            b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 792 612] "
            b"/Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>"
        ),
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
        b"<< /Length " + str(len(content)).encode("ascii") + b" >>\nstream\n" + content + b"endstream",
    ]
    document = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
    offsets = [0]
    for number, body in enumerate(objects, 1):
        offsets.append(len(document))
        document.extend(f"{number} 0 obj\n".encode("ascii"))
        document.extend(body)
        document.extend(b"\nendobj\n")
    xref_offset = len(document)
    document.extend(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
    document.extend(b"0000000000 65535 f \n")
    for offset in offsets[1:]:
        document.extend(f"{offset:010d} 00000 n \n".encode("ascii"))
    document.extend(
        (
            f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\n"
            f"startxref\n{xref_offset}\n%%EOF\n"
        ).encode("ascii")
    )
    return bytes(document)


def validate_pdf(path: Path) -> None:
    try:
        payload = path.read_bytes()
    except OSError as exc:
        raise PlotInputError(f"could not read rendered PDF {path}: {exc}") from exc
    if len(payload) < 256 or not payload.startswith(b"%PDF-"):
        raise PlotInputError(f"rendered file is not a nonempty PDF: {path}")
    if not payload.rstrip().endswith(b"%%EOF") or b"\nxref\n" not in payload:
        raise PlotInputError(f"rendered PDF is incomplete: {path}")


def render_plot(
    depth_path: Path,
    output_path: Path,
    metric: str,
    summary_path: Path | None = None,
    title: str = "CNV target-window coverage diagnostic",
) -> None:
    if metric not in METRIC_LABELS:
        raise PlotInputError(f"unsupported metric: {metric}")
    order, stats, _ = _scan_depth(depth_path)
    binned = _binned_depth(depth_path, order, stats)
    summary = _load_summary(summary_path)
    content = _content_stream(
        order, stats, binned, metric, METRIC_LABELS[metric], title, summary
    )
    payload = _pdf_document(content)

    output_path = output_path.resolve()
    if not output_path.parent.is_dir():
        raise PlotInputError(f"output directory does not exist: {output_path.parent}")
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.", suffix=".tmp", dir=output_path.parent
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        validate_pdf(temporary_path)
        os.replace(temporary_path, output_path)
        validate_pdf(output_path)
    finally:
        try:
            temporary_path.unlink()
        except FileNotFoundError:
            pass


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Render a run-scoped coverage diagnostic as a dependency-free PDF."
    )
    parser.add_argument("--depth-tsv", required=True, type=Path)
    parser.add_argument("--output-pdf", required=True, type=Path)
    parser.add_argument("--metric", required=True, choices=sorted(METRIC_LABELS))
    parser.add_argument("--summary-json", type=Path)
    parser.add_argument("--title", default="CNV target-window coverage diagnostic")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    try:
        render_plot(
            depth_path=arguments.depth_tsv,
            output_path=arguments.output_pdf,
            metric=arguments.metric,
            summary_path=arguments.summary_json,
            title=arguments.title,
        )
    except PlotInputError as exc:
        raise SystemExit(f"plot_cnv_depth_v2.py: {exc}") from exc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
