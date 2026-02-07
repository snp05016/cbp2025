from __future__ import annotations

from pathlib import Path

BEGIN = r"\begin{document}"
END = r"\end{document}"


def extract_document_body(tex: str) -> str:
    """
    Return content strictly between \\begin{document} and \\end{document}.
    If markers aren't found, return the original content unchanged.
    """
    b = tex.find(BEGIN)
    e = tex.rfind(END)
    if b == -1 or e == -1 or e <= b:
        return tex
    return tex[b + len(BEGIN) : e]


def main() -> None:
    here = Path(__file__).resolve().parent
    out = here / "combined.tex"

    tex_files = sorted(
        [p for p in here.glob("*.tex") if p.name != out.name],
        key=lambda p: p.name.lower(),
    )

    parts: list[str] = []
    for p in tex_files:
        src = p.read_text(encoding="utf-8")
        body = extract_document_body(src)
        parts.append(f"% --- BEGIN: {p.name} ---\n{body}\n% --- END: {p.name} ---\n\\clearpage\n")

    combined = (
        "\\documentclass{article}\n"
        "\\usepackage[margin=1in]{geometry}\n"
        "\\usepackage{graphicx}\n"
        "\\usepackage{booktabs}\n"
        "\\usepackage{amsmath}\n"
        "\\usepackage{caption}\n"
        "\n"
        "\\begin{document}\n\n"
        + "".join(parts)
        + "\n\\end{document}\n"
    )

    out.write_text(combined, encoding="utf-8")
    print(f"Wrote {out} (combined {len(tex_files)} files).")


if __name__ == "__main__":
    main()
