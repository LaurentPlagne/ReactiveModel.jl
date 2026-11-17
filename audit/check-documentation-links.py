"""Check local HTML links, anchors and assets after building docs/make.jl.

Run from any directory: python3 audit/check-documentation-links.py
Documenter deployment metadata scripts are optional in a local-only build.
"""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


class Page(HTMLParser):
    def __init__(self, source):
        super().__init__()
        self.refs = []
        self.ids = set()
        self.feed(source)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            self.ids.add(attrs["id"])
        key = "href" if tag in ("a", "link") else "src" if tag in ("img", "script") else None
        if key and attrs.get(key):
            self.refs.append((tag, attrs[key]))


repo = Path(__file__).resolve().parents[1]
root = repo / "docs/build"
pages = {p: Page(p.read_text()) for p in root.rglob("*.html")}
assert pages, "Build the documentation first"
errors = []
checked = 0
optional = 0
for path, page in pages.items():
    for tag, ref in page.refs:
        url = urlsplit(ref)
        if url.scheme or url.netloc:
            if url.scheme == "file":
                errors.append((str(path.relative_to(root)), ref, "machine-specific link"))
            continue
        target = (path.parent / unquote(url.path)).resolve() if url.path else path
        if target.is_dir():
            target /= "index.html"
        # These scripts are inserted by Documenter and populated by deploydocs,
        # which is intentionally not run for a local manual.
        if tag == "script" and target in (root / "siteinfo.js", root.parent / "versions.js"):
            optional += 1
            continue
        checked += 1
        if not target.exists():
            errors.append((str(path.relative_to(root)), ref, "missing file"))
        elif url.fragment and target in pages and unquote(url.fragment) not in pages[target].ids:
            errors.append((str(path.relative_to(root)), ref, "missing anchor"))
blocks = sum(p.read_text().count("```@example ") for p in (repo / "docs/src").rglob("*.md"))
print(f"{len(pages)} HTML pages; {checked} local references/assets checked; {len(errors)} errors.")
print(f"{blocks} executable Documenter example blocks.")
print(f"{optional} optional deployment-script references excluded (siteinfo.js and versions.js).")
for error in errors:
    print(error)
raise SystemExit(bool(errors))
