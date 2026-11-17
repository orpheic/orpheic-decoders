#!/usr/bin/env python3
"""Write the public copy of the decoders (src/) from the Orpheic program's source tree.

    tools/export.py <path to the Orpheic source tree> [--out <directory>]

The program is the single source of truth. This script copies the decoder files as they are, with one
transformation: code the program keeps only for its own comparison tests is compiled there under
`#if ORPHEIC_MP3_PROOF`, and here those blocks are resolved as if the macro were 0 (the `#else` branch
is kept, the directives are dropped) and the guarded default definition of the macro is removed.

It stops with an error, writing nothing, if a directive mentions the macro in a form it does not
understand, if a block is not closed, or if the result still contains a name that belongs to the
removed code.
"""
import os
import re
import sys

FILES = [
    # MP3 (ISO/IEC 11172-3 and 13818-3, Layer III)
    "mp3decoder.cpp", "mp3decoder.h", "mp3exact.h", "mp3huff.cpp", "mp3huff.h", "mp3huff_kernel.inc",
    "mp3hybrid.cpp", "mp3hybrid.h", "mp3hybrid_kernel.inc", "mp3synth.cpp", "mp3synth.h",
    "mp3synth_kernel.inc", "mp3tables.inc", "mp3consts.h", "mp3consts_tables.inc",
    # FLAC
    "flacdecoder.cpp", "flacdecoder.h", "flacdecoder_rice.inc",
    # DST (the lossless compression of DSD on SACD)
    "dstdecoder.cpp", "dstdecoder.h",
]
MACRO = "ORPHEIC_MP3_PROOF"
FORBIDDEN = ["likeAvcodec", "LikeAvcodec", MACRO]

IF_ANY = re.compile(r"^\s*#\s*(if|ifdef|ifndef)\b")
ELSE = re.compile(r"^\s*#\s*else\b")
ELIF = re.compile(r"^\s*#\s*elif\b")
ENDIF = re.compile(r"^\s*#\s*endif\b")
IF_MACRO = re.compile(r"^\s*#\s*if\s+" + MACRO + r"\s*(//.*)?$")
IFNDEF_MACRO = re.compile(r"^\s*#\s*ifndef\s+" + MACRO + r"\s*(//.*)?$")


class ExportError(Exception):
    pass


def strip(name, text):
    """Resolve the marked blocks of one file as if the macro were 0."""
    out = []
    # One entry per open conditional: "keep" (an ordinary one, copied through), "drop" (the macro's
    # #if branch), "else" (the macro's #else branch, kept), "default" (the #ifndef default, dropped).
    stack = []
    lines = text.split("\n")
    for number, line in enumerate(lines, 1):
        where = f"{name}:{number}"
        dropping = any(kind in ("drop", "default") for kind in stack)
        if IF_MACRO.match(line):
            stack.append("drop")
            continue
        if IFNDEF_MACRO.match(line):
            stack.append("default")
            continue
        if MACRO in line and re.match(r"^\s*#", line) and not (stack and stack[-1] == "default"):
            raise ExportError(f"{where}: a directive uses {MACRO} in a form this script does not resolve: {line.strip()}")
        if IF_ANY.match(line):
            stack.append("keep")
            if not dropping:
                out.append(line)
            continue
        if ELIF.match(line):
            if not stack:
                raise ExportError(f"{where}: #elif without #if")
            if stack[-1] != "keep":
                raise ExportError(f"{where}: #elif inside a {MACRO} block is not supported")
            if not dropping:
                out.append(line)
            continue
        if ELSE.match(line):
            if not stack:
                raise ExportError(f"{where}: #else without #if")
            top = stack[-1]
            if top == "drop":
                stack[-1] = "else"
            elif top == "else":
                raise ExportError(f"{where}: second #else in a {MACRO} block")
            elif top == "default":
                raise ExportError(f"{where}: #else in the default definition of {MACRO}")
            elif not dropping:
                out.append(line)
            continue
        if ENDIF.match(line):
            if not stack:
                raise ExportError(f"{where}: #endif without #if")
            top = stack.pop()
            if top == "keep":
                if not dropping:
                    out.append(line)
            elif top == "default":
                # the definition sat between two blank lines: leave one
                following = lines[number] if number < len(lines) else None
                if out and out[-1].strip() == "" and following is not None and following.strip() == "":
                    out.pop()
            continue
        if not dropping:
            out.append(line)
    if stack:
        raise ExportError(f"{name}: {len(stack)} conditional block(s) not closed at the end of the file")
    result = "\n".join(out)
    for word in FORBIDDEN:
        if word in result:
            at = result[: result.index(word)].count("\n") + 1
            raise ExportError(f"{name}: the exported text still contains '{word}' (line {at} of the output)")
    return result


def main(argv):
    if len(argv) < 2 or argv[1] in ("-h", "--help"):
        print(__doc__)
        return 2
    source = argv[1]
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(here, "..", "src")
    if "--out" in argv:
        out_dir = argv[argv.index("--out") + 1]
    out_dir = os.path.abspath(out_dir)
    results = {}
    try:
        for name in FILES:
            path = os.path.join(source, name)
            if not os.path.isfile(path):
                raise ExportError(f"{path}: not found")
            with open(path, encoding="utf-8") as f:
                results[name] = (f.read(), None)
            results[name] = (results[name][0], strip(name, results[name][0]))
    except ExportError as e:
        print(f"export.py: {e}\nNothing was written.", file=sys.stderr)
        return 1
    os.makedirs(out_dir, exist_ok=True)
    removed_total = 0
    for name in FILES:
        before, after = results[name]
        with open(os.path.join(out_dir, name), "w", encoding="utf-8") as f:
            f.write(after)
        removed = before.count("\n") - after.count("\n")
        removed_total += removed
        print(f"  {name}: {after.count(chr(10))} lines" + (f" ({removed} removed)" if removed else ""))
    print(f"export.py: {len(FILES)} files written to {out_dir}, {removed_total} lines removed in all")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
