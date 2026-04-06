#!/usr/bin/env python3
"""Format generated Rd files to a maximum line width.

This is intended as a post-processing step after roxygen2::roxygenize().
It rewrites common roxygen-generated long-line patterns in man/*.Rd while
preserving parseable Rd syntax.
"""

from __future__ import annotations

import argparse
import re
import textwrap
from pathlib import Path


PRE_BLOCK_RE = re.compile(
    r'\\if\{html\}\{\\out\{<div class="r">\}\}\\preformatted\{([^{}]*)\}\\if\{html\}\{\\out\{</div>\}\}'
)
PKG_LINK_RE = re.compile(
    r"<li><span class=\"pkg-link\"[^>]*><a href='[^']*'><code>([^<]*)</code></a></span></li>"
)
ITEM_HREF_RE = re.compile(r"^(\\item )\\href\{([^}]*)\}\{(\\code\{.*\})\}$")
ITEM_TWOARG_RE = re.compile(r"^(\\item\{[^}]*\})(\{.*\})$")
CODE_LINK_RE = re.compile(r"\\code\{\\link\[([^]:]+):([^\]]+)\]\{([^}]*)\}\}")
LATEX_HYPERTARGET_RE = re.compile(r"^\\if\{latex\}\{\\out\{\\hypertarget\{.*\}\{\}\}\}$")
ITEM_ONLY_RE = re.compile(r"^\\item(?:\s+)?$")


def wrap_plain(line: str, width: int) -> list[str]:
    if len(line) <= width or " " not in line:
        return [line]
    return textwrap.wrap(
        line,
        width=width,
        break_long_words=False,
        break_on_hyphens=False,
        subsequent_indent="  ",
    )


def transform_line(line: str, width: int) -> list[str]:
    line = PRE_BLOCK_RE.sub(lambda m: r"\preformatted{" + m.group(1) + "}", line)
    line = PKG_LINK_RE.sub(lambda m: "<li><code>" + m.group(1) + "</code></li>", line)
    line = line.replace("\\\\preformatted{", r"\preformatted{")

    indent_len = len(line) - len(line.lstrip(" "))
    indent = line[:indent_len]
    stripped = line[indent_len:]

    if len(line) <= width:
        return [line]

    match = ITEM_HREF_RE.match(stripped)
    if match:
        return [
            indent + f"{match.group(1)}\\href{{{match.group(2)}}}",
            indent + match.group(3),
        ]

    match = ITEM_TWOARG_RE.match(stripped)
    if match:
        return [indent + match.group(1), indent + match.group(2)]

    if " -> " in line:
        parts = line.split(" -> ")
        return [parts[0], *["-> " + part for part in parts[1:]]]

    if stripped.startswith("<li><code>") and "</code></li>" in stripped:
        return [indent + "<li>", indent + stripped[4:], indent + "</li>"]

    if stripped.startswith(r"\preformatted{") and stripped.endswith("}") and "$" in stripped:
        left, right = stripped.split("$", 1)
        return [indent + left + "$", indent + right]

    if stripped.startswith(r"\code{") and "$" in stripped:
        left, right = stripped.split("$", 1)
        return [indent + left + "$", indent + right]

    if stripped.startswith("{") and stripped.endswith("}") and " " in stripped:
        wrapped = textwrap.wrap(
            stripped[1:-1],
            width=max(1, width - 2),
            break_long_words=False,
            break_on_hyphens=False,
        )
        if len(wrapped) > 1:
            return [
                indent + "{" + wrapped[0],
                *[indent + piece for piece in wrapped[1:-1]],
                indent + wrapped[-1] + "}",
            ]

    return wrap_plain(line, width)


def transform_file(path: Path, width: int) -> bool:
    lines = path.read_text().splitlines()
    out: list[str] = []
    changed = False
    i = 0

    while i < len(lines):
        line = lines[i]

        if LATEX_HYPERTARGET_RE.match(line):
            changed = True
            i += 1
            continue

        if line == r"\if{html}{\out{<a" and i + 1 < len(lines) and lines[i + 1].lstrip().startswith('id="method-'):
            changed = True
            i += 2
            continue

        if line == r"\if{html}{\out{<div":
            changed = True
            i += 1
            continue

        if line.startswith('  class="r">}}\\preformatted{'):
            changed = True
            line = r"\preformatted{" + line.split("}}\\preformatted{", 1)[1]

        if ITEM_ONLY_RE.match(line) and i + 2 < len(lines) and lines[i + 1].startswith(r"\href{") and lines[i + 2].startswith(r"\code{"):
            pieces = transform_line(r"\item " + lines[i + 2], width)
            out.extend(pieces)
            changed = True
            i += 3
            continue

        if line.startswith(r"\item ") and i + 1 < len(lines) and line[len(r"\item "):].startswith(r"\href{") and lines[i + 1].startswith(r"\code{"):
            pieces = transform_line(r"\item " + lines[i + 1], width)
            out.extend(pieces)
            changed = True
            i += 2
            continue

        if ITEM_ONLY_RE.match(line) and i + 1 < len(lines) and lines[i + 1].startswith(r"\code{"):
            pieces = transform_line(r"\item " + lines[i + 1], width)
            out.extend(pieces)
            changed = True
            i += 2
            continue

        if i + 1 < len(lines) and line.startswith(r"\item{") and line.endswith("}") and lines[i + 1].startswith("{"):
            line = line + lines[i + 1]
            changed = True
            i += 1

        if len(line) > width and line.startswith(r"-> \code{"):
            out.append("->")
            out.append(line[3:])
            changed = True
            continue

        new_line = CODE_LINK_RE.sub(lambda m: r"\code{" + m.group(3) + "}", line)
        if new_line != line:
            changed = True
        line = new_line

        pieces = transform_line(line, width)
        if pieces != [line]:
            changed = True
        out.extend(pieces)
        i += 1

    if changed:
        path.write_text("\n".join(out) + "\n")
    return changed


def count_overlong_lines(man_dir: Path, width: int) -> int:
    count = 0
    for path in sorted(man_dir.glob("*.Rd")):
        for line in path.read_text().splitlines():
            if len(line) > width:
                count += 1
    return count


def main() -> int:
    parser = argparse.ArgumentParser(description="Format generated Rd files to a maximum line width.")
    parser.add_argument("man_dir", nargs="?", default="EDI/man", help="Directory containing .Rd files.")
    parser.add_argument("--width", type=int, default=100, help="Maximum line width.")
    args = parser.parse_args()

    man_dir = Path(args.man_dir)
    if not man_dir.is_dir():
        raise SystemExit(f"{man_dir} is not a directory")

    changed_files = 0
    for path in sorted(man_dir.glob("*.Rd")):
        if transform_file(path, args.width):
            changed_files += 1

    remaining = count_overlong_lines(man_dir, args.width)
    print(f"formatted_files={changed_files}")
    print(f"remaining_overlong_lines={remaining}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
