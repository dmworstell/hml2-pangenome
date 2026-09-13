#!/usr/bin/env python3
"""Minimal, dependency-free extraction of a flat first-sheet XLSX table."""

import csv
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path


NS = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def col_number(cell_ref: str) -> int:
    letters = re.match(r"[A-Z]+", cell_ref).group(0)
    value = 0
    for letter in letters:
        value = value * 26 + ord(letter) - ord("A") + 1
    return value


def extract(source: Path, destination: Path) -> None:
    with zipfile.ZipFile(source) as archive:
        shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        shared = [
            "".join(node.text or "" for node in item.findall(".//x:t", NS))
            for item in shared_root.findall("x:si", NS)
        ]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))

    rows = []
    for row_node in sheet.findall(".//x:sheetData/x:row", NS):
        sparse = {}
        for cell in row_node.findall("x:c", NS):
            value_node = cell.find("x:v", NS)
            value = "" if value_node is None else value_node.text
            if cell.get("t") == "s" and value:
                value = shared[int(value)]
            sparse[col_number(cell.get("r"))] = value
        rows.append([sparse.get(i, "") for i in range(1, max(sparse, default=0) + 1)])

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: extract_xlsx.py SOURCE.xlsx DESTINATION.tsv")
    extract(Path(sys.argv[1]), Path(sys.argv[2]))
