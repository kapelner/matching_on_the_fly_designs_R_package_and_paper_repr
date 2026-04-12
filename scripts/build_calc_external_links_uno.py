#!/usr/bin/env python3
"""Build Calc external-link workbooks via LibreOffice UNO.

Workflow:
1. Import each CSV into a hidden Calc document.
2. Save each imported CSV as an `.ods` source workbook with a named range `DATA`.
3. Create a master `.ods` workbook with one linked sheet per source workbook.
4. Set each area link refresh period to 60 seconds.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import uno


def make_prop(name: str, value):
    prop = uno.createUnoStruct("com.sun.star.beans.PropertyValue")
    prop.Name = name
    prop.Value = value
    return prop


def connect(host: str, port: int):
    local_ctx = uno.getComponentContext()
    resolver = local_ctx.ServiceManager.createInstanceWithContext(
        "com.sun.star.bridge.UnoUrlResolver", local_ctx
    )
    return resolver.resolve(
        f"uno:socket,host={host},port={port};urp;StarOffice.ComponentContext"
    )


def cell_address(sheet_idx: int, col: int, row: int):
    addr = uno.createUnoStruct("com.sun.star.table.CellAddress")
    addr.Sheet = sheet_idx
    addr.Column = col
    addr.Row = row
    return addr


def response_name(csv_path: Path) -> str:
    stem = csv_path.stem
    prefix = "comprehensive_tests_results_nc_1_"
    return stem[len(prefix) :] if stem.startswith(prefix) else stem


def csv_shape(csv_path: Path) -> tuple[int, int]:
    rows = 0
    cols = 0
    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        for row in reader:
            rows += 1
            cols = max(cols, len(row))
    return max(rows, 1), max(cols, 1)


def col_label(index: int) -> str:
    result = []
    n = index + 1
    while n:
        n, rem = divmod(n - 1, 26)
        result.append(chr(ord("A") + rem))
    return "".join(reversed(result))


def absolute_range(sheet_name: str, rows: int, cols: int) -> str:
    return f"${sheet_name}.$A$1:${col_label(cols - 1)}${rows}"


def import_csv_to_source(desktop, csv_path: Path, source_path: Path):
    rows, cols = csv_shape(csv_path)
    csv_url = uno.systemPathToFileUrl(str(csv_path.resolve()))
    source_url = uno.systemPathToFileUrl(str(source_path.resolve()))
    load_props = (
        make_prop("Hidden", True),
        make_prop("FilterName", "Text - txt - csv (StarCalc)"),
        make_prop("FilterOptions", "44,34,76,1"),
    )
    doc = desktop.loadComponentFromURL(csv_url, "_blank", 0, load_props)
    try:
        sheets = doc.getSheets()
        sheet = sheets.getByIndex(0)
        sheet.setName("Data")

        named_ranges = doc.getPropertyValue("NamedRanges")
        base = cell_address(0, 0, 0)
        named_ranges.addNewByName("DATA", absolute_range("Data", rows, cols), base, 0)

        store_props = (
            make_prop("FilterName", "calc8"),
            make_prop("Overwrite", True),
        )
        doc.storeAsURL(source_url, store_props)
    finally:
        doc.close(True)

    return rows, cols, source_url


def fill_manifest(sheet, rows_data):
    data = tuple(tuple(str(v) for v in row) for row in rows_data)
    end_col = len(rows_data[0]) - 1
    end_row = len(rows_data) - 1
    sheet.getCellRangeByPosition(0, 0, end_col, end_row).setDataArray(data)


def freeze_top_row(controller, sheet):
    controller.setActiveSheet(sheet)
    controller.freezeAtPosition(0, 1)


def autosize_columns(sheet, start_col: int, end_col: int):
    cols = sheet.getColumns()
    for col_idx in range(start_col, end_col + 1):
        cols.getByIndex(col_idx).OptimalWidth = True


def build_master(desktop, links, master_path: Path, refresh_seconds: int):
    temp_master_path = master_path.with_name(master_path.stem + ".tmp.ods")
    master_url = uno.systemPathToFileUrl(str(temp_master_path.resolve()))
    doc = desktop.loadComponentFromURL(
        "private:factory/scalc", "_blank", 0, (make_prop("Hidden", True),)
    )
    try:
        sheets = doc.getSheets()
        first = sheets.getByIndex(0)
        first.setName("manifest")

        manifest_rows = [["Sheet", "CSV", "Source ODS"]]
        for idx, link in enumerate(links, start=1):
            sheets.insertNewByName(link["sheet_name"], idx)
            manifest_rows.append([link["sheet_name"], link["csv_path"], link["source_path"]])

        fill_manifest(sheets.getByName("manifest"), manifest_rows)
        controller = doc.getCurrentController()
        freeze_top_row(controller, sheets.getByName("manifest"))
        autosize_columns(sheets.getByName("manifest"), 0, 8)

        area_links = doc.getPropertyValue("AreaLinks")
        for idx, link in enumerate(links, start=1):
            dest = cell_address(idx, 0, 0)
            area_links.insertAtPosition(dest, link["source_url"], "DATA", "", "")

        for i in range(area_links.getCount()):
            area_link = area_links.getByIndex(i)
            try:
                area_link.setPropertyValue("RefreshPeriod", refresh_seconds)
            except Exception:
                pass
            try:
                area_link.setPropertyValue("RefreshDelay", refresh_seconds)
            except Exception:
                pass
            area_link.refresh()

        for idx, link in enumerate(links, start=1):
            sheet = sheets.getByName(link["sheet_name"])
            freeze_top_row(controller, sheet)
            autosize_columns(sheet, 0, 8)

        store_props = (
            make_prop("FilterName", "calc8"),
            make_prop("Overwrite", True),
        )
        doc.storeAsURL(master_url, store_props)
    finally:
        doc.close(True)

    temp_master_path.replace(master_path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "csv_glob",
        nargs="?",
        default="package_tests/comprehensive_tests_results_nc_1_*.csv",
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=2002)
    parser.add_argument("--refresh-seconds", type=int, default=60)
    parser.add_argument("--output-dir", default="package_tests/calc_external_links_uno")
    args = parser.parse_args()

    ctx = connect(args.host, args.port)
    desktop = ctx.ServiceManager.createInstanceWithContext("com.sun.star.frame.Desktop", ctx)

    csv_files = sorted(Path().glob(args.csv_glob))
    if not csv_files:
        raise SystemExit(f"No CSV files matched: {args.csv_glob}")

    output_dir = Path(args.output_dir)
    sources_dir = output_dir / "sources"
    sources_dir.mkdir(parents=True, exist_ok=True)

    links = []
    for csv_path in csv_files:
        source_path = sources_dir / f"{csv_path.stem}.ods"
        rows, cols, source_url = import_csv_to_source(desktop, csv_path, source_path)
        links.append(
            {
                "sheet_name": response_name(csv_path),
                "csv_path": str(csv_path.resolve()),
                "source_path": str(source_path.resolve()),
                "source_url": source_url,
                "rows": rows,
                "cols": cols,
            }
        )

    master_path = output_dir / "comprehensive_tests_results_nc_1_links.ods"
    build_master(desktop, links, master_path, args.refresh_seconds)
    print(master_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
