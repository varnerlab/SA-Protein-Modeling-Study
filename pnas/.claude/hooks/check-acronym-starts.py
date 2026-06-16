#!/usr/bin/env python3
"""Flag sentences that begin with an acronym in LaTeX prose.

An "acronym" here is any first-word token containing two or more consecutive
capital letters (SA, DMS, MSA) -- including lowercase-led forms such as
pLDDT or mRNA, which a capital-anchored search misses. Title-case words
(The, Performance) and CamelCase product names (EvoDiff) are NOT flagged,
because they have no run of two consecutive capitals.

Two run modes:

  1. CLI:   python3 check-acronym-starts.py main.tex si-appendix.tex
            Prints a human-readable report. Exit 1 if any violation, else 0.

  2. Hook:  receives PostToolUse JSON on stdin, scans the edited .tex file,
            and emits hook JSON (systemMessage + additionalContext) when it
            finds a violation. Silent and exit 0 otherwise.

To accept a particular sentence-initial token, add it to IGNORE below.
"""

import json
import re
import sys

# Tokens that are acceptable at a sentence start despite the capital run.
IGNORE = {"ORCID"}

# Abbreviations whose trailing period must not be read as a sentence end.
ABBREVS = [
    "vs", "e.g", "i.e", "cf", "al", "etc", "Fig", "Figs", "Eq", "Eqs",
    "Sec", "Secs", "Ref", "Refs", "No", "pp", "Table", "Tab", "Eqn",
    "approx", "resp", "viz", "ca", "Dr", "Prof", "St", "Inc",
]

# Reference-style macros whose braced argument should be dropped wholesale,
# so a label like \ref{SI-fig:...} never masquerades as the word "SI".
_REF_MACROS = r"\\(?:ref|eqref|autoref|cref|Cref|Vref|vref|pageref|label|cite\w*)"

# Environments that are not prose: their bodies are removed before scanning so
# table cells (RRM, SH3, ...) and display-math rows are never read as sentences.
# Captions live OUTSIDE these (in the surrounding table/figure float) and are
# kept, since a sentence-initial acronym in a caption is a real violation.
_NON_PROSE_ENVS = [
    "tabular", "tabular*", "array", "longtable", "longtable*", "tabu",
    "supertabular", "align", "align*", "alignat", "alignat*", "flalign",
    "flalign*", "equation", "equation*", "gather", "gather*", "multline",
    "multline*", "split", "eqnarray", "eqnarray*", "matrix", "pmatrix",
    "bmatrix", "Bmatrix", "vmatrix", "Vmatrix", "smallmatrix", "cases",
    "algorithmic",
]

# Metadata macros whose (possibly nested) braced argument is not prose.
_META_MACROS = [
    "keywords", "title", "author", "affil", "date", "dates", "ORCID",
    "correspondingauthor", "leadauthor", "doi", "equalauthors",
    "authorcontributions",
]


def _strip_balanced(text, macro):
    """Remove every \\macro{...} call, honoring nested braces."""
    out = []
    i = 0
    needle = "\\" + macro
    while i < len(text):
        j = text.find(needle, i)
        if j == -1:
            out.append(text[i:])
            break
        out.append(text[i:j])
        k = j + len(needle)
        while k < len(text) and text[k] in " \t":
            k += 1
        if k < len(text) and text[k] == "{":
            depth = 0
            while k < len(text):
                if text[k] == "{":
                    depth += 1
                elif text[k] == "}":
                    depth -= 1
                    if depth == 0:
                        k += 1
                        break
                k += 1
            i = k
        else:
            out.append(needle)
            i = j + len(needle)
    return "".join(out)


def find_violations(path):
    """Return a list of (line, token, snippet) for sentence-initial acronyms."""
    try:
        with open(path, encoding="utf-8") as fh:
            raw = fh.read()
    except (OSError, UnicodeDecodeError):
        return []

    # Drop whole-line LaTeX comments (a line whose first non-space char is %).
    body_lines = [ln for ln in raw.split("\n") if not ln.lstrip().startswith("%")]
    text = "\n".join(body_lines)

    # Remove non-prose environment bodies (tables, display math, pseudocode).
    for env in _NON_PROSE_ENVS:
        pat = r"\\begin\{" + re.escape(env) + r"\}.*?\\end\{" + re.escape(env) + r"\}"
        text = re.sub(pat, " ", text, flags=re.DOTALL)

    # Remove metadata macro arguments (keyword list, author block, title, ...).
    for macro in _META_MACROS:
        text = _strip_balanced(text, macro)

    # Strip constructs that would inject spurious "words".
    text = re.sub(r"~?" + _REF_MACROS + r"\s*\{[^}]*\}", " ", text)  # \ref{...}
    # Math becomes a lowercase placeholder, not a blank: a sentence that truly
    # begins with a math symbol (e.g. "$d$: PCA dimension") then starts with the
    # placeholder, so the following acronym is not misread as sentence-initial.
    text = re.sub(r"\\\[.*?\\\]", " matheq ", text, flags=re.DOTALL)  # \[ ... \]
    text = re.sub(r"\\\(.*?\\\)", " matheq ", text)                   # \( ... \)
    text = re.sub(r"\$[^$]*\$", " matheq ", text)                     # $ ... $
    text = re.sub(r"\\[a-zA-Z@]+", " ", text)                          # \command
    text = re.sub(r"[{}\[\]]", " ", text)                             # stray braces

    # Protect abbreviation periods from the sentence splitter.
    for ab in ABBREVS:
        text = re.sub(r"\b" + re.escape(ab) + r"\.", ab + "\x00", text)

    out = []
    for sentence in re.split(r"(?<=[.!?])\s+", text):
        sentence = sentence.replace("\x00", ".").strip()
        # First real word: skip leading non-letter noise.
        m = re.search(r"[A-Za-z][A-Za-z0-9]*", sentence)
        if not m:
            continue
        token = m.group(0)
        after = sentence[m.end():m.end() + 1]
        if after == ":":            # field label ("ORCID:"), not prose
            continue
        if token in IGNORE:
            continue
        if not re.search(r"[A-Z]{2,}", token):  # needs a run of >=2 capitals
            continue

        # Best-effort line number: locate the token (plus the following word,
        # to disambiguate common acronyms) in the original file.
        nxt = re.search(r"[A-Za-z][A-Za-z0-9]*", sentence[m.end():])
        if nxt:
            locate = re.compile(re.escape(token) + r"[^A-Za-z]{0,40}" + re.escape(nxt.group(0)))
        else:
            locate = re.compile(re.escape(token))
        hit = locate.search(raw)
        line = raw.count("\n", 0, hit.start()) + 1 if hit else None
        out.append((line, token, sentence[:72]))
    return out


def cli(paths):
    total = 0
    for path in paths:
        for line, token, snippet in find_violations(path):
            total += 1
            loc = f"{path}:{line}" if line else path
            print(f"  {loc}  [{token}]  {snippet}...")
    if total:
        print(f"\n{total} sentence-initial acronym(s). Rewrite so no sentence begins with one.")
        return 1
    print("No sentence-initial acronyms found.")
    return 0


def hook():
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0
    path = (payload.get("tool_input") or {}).get("file_path", "")
    if not path.endswith(".tex"):
        return 0
    violations = find_violations(path)
    if not violations:
        return 0

    name = path.rsplit("/", 1)[-1]
    lines = [f"  L{line} [{tok}] \"{snip}...\"" for line, tok, snip in violations]
    summary = "\n".join(lines)
    msg = (f"Acronym check: {len(violations)} sentence-initial acronym(s) in {name}\n"
           f"{summary}")
    context = (
        f"The style check found {len(violations)} sentence(s) beginning with an "
        f"acronym in {name} (the author forbids starting a sentence with an "
        f"acronym). Rewrite each so the sentence does not begin with the flagged "
        f"token:\n{summary}"
    )
    print(json.dumps({
        "systemMessage": msg,
        "hookSpecificOutput": {
            "hookEventName": "PostToolUse",
            "additionalContext": context,
        },
    }))
    return 0


if __name__ == "__main__":
    sys.exit(cli(sys.argv[1:]) if len(sys.argv) > 1 else hook())
