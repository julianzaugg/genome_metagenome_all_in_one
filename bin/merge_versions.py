#!/usr/bin/env python3
"""
Merge one or more versions.yml files emitted by Nextflow processes into
a single software_versions.tsv with columns: tool, version, process.

Handles both YAML formats used in this pipeline:
  Block:  "PROC_NAME":\n    tool: version
  Flow:   "PROC_NAME": {tool: version}
"""
import re
import sys


def parse(text):
    entries = []
    # Flow-style lines: "PROC": {tool: version[, ...]}
    for m in re.finditer(r'^"([^"]+)":\s*\{([^}]+)\}', text, re.MULTILINE):
        proc = m.group(1)
        for pair in re.finditer(r'([\w][\w-]*):\s*([^,}]+)', m.group(2)):
            entries.append((proc, pair.group(1), pair.group(2).strip()))
    # Block-style: "PROC":\n    tool: version
    for m in re.finditer(r'^"([^"]+)":\s*\n((?:[ \t]+\S[^\n]*\n?)+)', text, re.MULTILINE):
        proc = m.group(1)
        for line in m.group(2).splitlines():
            kv = re.match(r'[ \t]+([\w][\w-]*):\s*(.*)', line)
            if kv:
                entries.append((proc, kv.group(1), kv.group(2).strip()))
    return entries


def main():
    seen = set()
    rows = []
    for path in sys.argv[1:]:
        try:
            text = open(path).read()
        except (FileNotFoundError, OSError):
            continue
        for proc, tool, version in parse(text):
            key = (tool.lower(), version)
            if key not in seen:
                seen.add(key)
                rows.append((tool, version, proc))

    rows.sort(key=lambda r: r[0].lower())
    print("tool\tversion\tprocess")
    for tool, version, proc in rows:
        print(f"{tool}\t{version}\t{proc}")


if __name__ == "__main__":
    main()
