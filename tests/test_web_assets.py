from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class IdCollector(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()
        self.scripts = []
        self.stylesheets = []

    def handle_starttag(self, tag, attrs):
        values = dict(attrs)
        if "id" in values:
            self.ids.add(values["id"])
        if tag == "script" and values.get("src"):
            self.scripts.append(values["src"])
        if tag == "link" and values.get("rel") == "stylesheet":
            self.stylesheets.append(values.get("href"))


def test_browser_entrypoint_has_required_controls():
    parser = IdCollector()
    parser.feed((ROOT / "index.html").read_text(encoding="utf-8"))

    required = {
        "operation",
        "set-a",
        "set-b",
        "analyze",
        "theme-toggle",
        "result-table",
        "download-results",
    }
    assert required <= parser.ids
    assert "app.js" in parser.scripts
    assert "styles.css" in parser.stylesheets


def test_browser_assets_exist_and_runtime_is_pinned():
    for name in ("app.js", "worker.js", "styles.css", "bed_overlap.py"):
        assert (ROOT / name).is_file(), name

    worker = (ROOT / "worker.js").read_text(encoding="utf-8")
    assert 'PYODIDE_VERSION = "314.0.7"' in worker
    assert 'fetch("bed_overlap.py"' in worker


def test_analyze_button_is_visible_action():
    html = (ROOT / "index.html").read_text(encoding="utf-8")
    assert 'id="analyze"' in html
    assert "Analyze intervals" in html
