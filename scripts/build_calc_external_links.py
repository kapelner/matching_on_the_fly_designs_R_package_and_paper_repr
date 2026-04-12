#!/usr/bin/env python3
"""Build Calc workbooks that link comprehensive test CSV outputs.

This creates:
1. One source `.ods` per CSV with a named range `DATA`.
2. One master `.ods` workbook with one sheet per source CSV. Each sheet is
   pre-populated with the current CSV contents and marked as an external link
   to the corresponding source `.ods` named range.

The master workbook can be opened in LibreOffice Calc and refreshed via
Sheet -> Link to External Data / Edit Links.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable
import zipfile
from xml.sax.saxutils import escape


NS = {
    "office": "urn:oasis:names:tc:opendocument:xmlns:office:1.0",
    "table": "urn:oasis:names:tc:opendocument:xmlns:table:1.0",
    "text": "urn:oasis:names:tc:opendocument:xmlns:text:1.0",
    "style": "urn:oasis:names:tc:opendocument:xmlns:style:1.0",
    "fo": "urn:oasis:names:tc:opendocument:xmlns:xsl-fo-compatible:1.0",
    "xlink": "http://www.w3.org/1999/xlink",
    "number": "urn:oasis:names:tc:opendocument:xmlns:datastyle:1.0",
    "svg": "urn:oasis:names:tc:opendocument:xmlns:svg-compatible:1.0",
    "of": "urn:oasis:names:tc:opendocument:xmlns:of:1.2",
}


def xml_text(value: str) -> str:
    return escape(value, {'"': "&quot;"})


def col_label(index: int) -> str:
    result = []
    n = index + 1
    while n:
        n, rem = divmod(n - 1, 26)
        result.append(chr(ord("A") + rem))
    return "".join(reversed(result))


def cell_xml(value: str, link_xml: str = "") -> str:
    return (
        '      <table:table-cell office:value-type="string">'
        f"{link_xml}<text:p>{xml_text(value)}</text:p>"
        "</table:table-cell>"
    )


def empty_cells_xml(count: int) -> str:
    if count <= 0:
        return ""
    return f'      <table:table-cell table:number-columns-repeated="{count}"/>'


def row_xml(values: list[str], first_link_xml: str = "") -> str:
    cells = []
    for idx, value in enumerate(values):
        cells.append(cell_xml(value, first_link_xml if idx == 0 else ""))
    return "    <table:table-row>\n" + "\n".join(cells) + "\n    </table:table-row>"


def named_range_xml(sheet_name: str, n_rows: int, n_cols: int) -> str:
    end_col = col_label(max(n_cols, 1) - 1)
    end_row = max(n_rows, 1)
    return (
        "  <table:named-expressions>\n"
        '    <table:named-range table:name="DATA" '
        f'table:base-cell-address="${xml_text(sheet_name)}.$A$1" '
        f'table:cell-range-address="${xml_text(sheet_name)}.$A$1:.${end_col}${end_row}"/>\n'
        "  </table:named-expressions>"
    )


def source_link_xml(source_uri: str, n_rows: int, n_cols: int, refresh_seconds: int) -> str:
    return (
        "\n        <table:cell-range-source "
        f'table:name="DATA" '
        f'table:last-column-spanned="{max(n_cols, 1)}" '
        f'table:last-row-spanned="{max(n_rows, 1)}" '
        f'table:refresh-delay="PT{refresh_seconds}S" '
        f'xlink:href="{xml_text(source_uri)}" '
        'xlink:type="simple" '
        'xlink:actuate="onLoad"/>'
    )


def table_xml(
    sheet_name: str,
    rows: list[list[str]],
    named_range: bool,
    first_link_xml: str = "",
) -> str:
    normalized_rows = rows or [[""]]
    n_cols = max((len(row) for row in normalized_rows), default=1)
    out = [f'<table:table table:name="{xml_text(sheet_name)}">']
    out.append(
        f'  <table:table-column table:number-columns-repeated="{max(n_cols, 1)}" '
        'table:default-cell-style-name="Default"/>'
    )
    for row_idx, row in enumerate(normalized_rows):
        padded = row + [""] * (n_cols - len(row))
        out.append(row_xml(padded, first_link_xml if row_idx == 0 else ""))
    out.append("</table:table>")
    if named_range:
        out.append(named_range_xml(sheet_name, len(normalized_rows), n_cols))
    return "\n".join(out)


def linked_table_xml(sheet_name: str, n_rows: int, n_cols: int, source_uri: str, refresh_seconds: int) -> str:
    cols = max(n_cols, 1)
    rows = max(n_rows, 1)
    link_xml = source_link_xml(source_uri, rows, cols, refresh_seconds)
    first_row = [""] * cols
    out = [f'<table:table table:name="{xml_text(sheet_name)}">']
    out.append(
        f'  <table:table-column table:number-columns-repeated="{cols}" '
        'table:default-cell-style-name="Default"/>'
    )
    out.append(row_xml(first_row, first_link_xml=link_xml))
    if rows > 1:
        out.append(f'  <table:table-row table:number-rows-repeated="{rows - 1}"/>')
    out.append("</table:table>")
    return "\n".join(out)


def document_xml(tables: Iterable[str]) -> str:
    attrs = " ".join(f'xmlns:{k}="{v}"' for k, v in NS.items())
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<office:document-content {attrs} office:version="1.2">\n'
        "  <office:scripts/>\n"
        "  <office:font-face-decls/>\n"
        "  <office:automatic-styles/>\n"
        "  <office:body>\n"
        "    <office:spreadsheet>\n"
        f"{chr(10).join(tables)}\n"
        "    </office:spreadsheet>\n"
        "  </office:body>\n"
        "</office:document-content>\n"
    )


def manifest_xml() -> str:
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<manifest:manifest '
        'xmlns:manifest="urn:oasis:names:tc:opendocument:xmlns:manifest:1.0" '
        'manifest:version="1.2">\n'
        '  <manifest:file-entry manifest:full-path="/" '
        'manifest:media-type="application/vnd.oasis.opendocument.spreadsheet"/>\n'
        '  <manifest:file-entry manifest:full-path="content.xml" manifest:media-type="text/xml"/>\n'
        '</manifest:manifest>\n'
    )


def write_ods(path: Path, content_xml: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr(
            "mimetype",
            "application/vnd.oasis.opendocument.spreadsheet",
            compress_type=zipfile.ZIP_STORED,
        )
        zf.writestr("content.xml", content_xml, compress_type=zipfile.ZIP_DEFLATED)
        zf.writestr("META-INF/manifest.xml", manifest_xml(), compress_type=zipfile.ZIP_DEFLATED)


def read_csv_rows(path: Path) -> list[list[str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return [list(row) for row in csv.reader(handle)]


def response_name(csv_path: Path) -> str:
    stem = csv_path.stem
    prefix = "comprehensive_tests_results_nc_1_"
    if stem.startswith(prefix):
        return stem[len(prefix):]
    return stem


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "csv_glob",
        nargs="?",
        default="package_tests/comprehensive_tests_results_nc_1_*.csv",
        help="Glob of CSV files to include.",
    )
    parser.add_argument(
        "--output-dir",
        default="package_tests/calc_external_links",
        help="Directory for generated ODS files.",
    )
    parser.add_argument(
        "--refresh-seconds",
        type=int,
        default=60,
        help="External link refresh interval in seconds.",
    )
    args = parser.parse_args()

    csv_files = sorted(Path().glob(args.csv_glob))
    if not csv_files:
        raise SystemExit(f"No CSV files matched: {args.csv_glob}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sources_dir = output_dir / "sources"
    sources_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = [["Sheet", "CSV", "Source ODS"]]
    master_tables = []

    for csv_path in csv_files:
        rows = read_csv_rows(csv_path)
        sheet_name = response_name(csv_path)
        source_path = sources_dir / f"{csv_path.stem}.ods"
        source_uri = source_path.resolve().as_uri()

        source_doc = document_xml([table_xml("Data", rows, named_range=True)])
        write_ods(source_path, source_doc)

        n_rows = max(len(rows), 1)
        n_cols = max((len(row) for row in rows), default=1)
        master_tables.append(
            linked_table_xml(sheet_name, n_rows, n_cols, source_uri, args.refresh_seconds)
        )

        manifest_rows.append([sheet_name, str(csv_path.resolve()), str(source_path.resolve())])

    master_tables.insert(0, table_xml("manifest", manifest_rows, named_range=False))
    master_path = output_dir / "comprehensive_tests_results_nc_1_links.ods"
    master_doc = document_xml(master_tables)
    write_ods(master_path, master_doc)

    print(f"Generated master workbook: {master_path}")
    print(f"Generated {len(csv_files)} source workbooks in: {sources_dir}")
    for csv_path in csv_files:
        print(f" - {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
