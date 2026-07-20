# Rejuvenation Atlas — static Explorer

A publication-styled, **server-free** alternative to the Streamlit companion app
(`app/`). It reads the same `results/` directory and renders the same five
modules, but as a single static HTML page instead of a running Python process.

Why this exists alongside the Streamlit app: Streamlit is convenient while a
module is under active development (hot reload, native widgets, zero build
step), but it requires a live Python server to view — awkward to share with a
collaborator, and its default chrome doesn't read as a finished, citable
artifact. This site is the "ship it" counterpart: build once, then the output
folder is a plain static site you can open locally, email as a zip, or drop
on any static host.

## Build

```bash
python site/build.py --results results/ --out site/public
```

Reads whatever module outputs exist under `results/` (missing ones render an
"run the pipeline" notice, same as the Streamlit app) and writes a
self-contained `site/public/`:

```
site/public/
├── index.html            # everything inlined: CSS, JS, chart data, PNGs (base64)
├── assets/plotly.min.js  # the one non-inlined file, copied from the plotly package
└── data/*.csv            # tables also offered as plain-CSV downloads
```

No Jinja/Node/webpack step — `build.py` only depends on `pandas`, already a
project dependency. Design tokens (palette, Plotly template) live in
`theme_tokens.py`, imported by `build.py`.

## View it

Just open the file — no server needed:

```bash
xdg-open site/public/index.html   # Linux
open site/public/index.html       # macOS
```

Or serve it like any static site (needed only if your browser blocks local
`file://` access, e.g. via `chrome://flags` restrictions):

```bash
python -m http.server -d site/public 8877
```

## Deploy

`site/public/` is a complete static site — push it to GitHub Pages, Netlify,
S3, or any static host. Nothing in it expects a backend.

## Design

Same design tokens as the Streamlit app's `theme.py` (validated categorical /
diverging palette, Source Serif headers over Inter body text, figure/table
captions numbered like a manuscript), implemented independently in
`static/style.css` and `static/app.js` — a plain CSS/JS pair, no framework.
Charts are Plotly.js, themed client-side so the **dark-mode toggle re-themes
the charts themselves**, not just the page chrome (persisted in
`localStorage`, defaults to the OS preference). Tables sort by clicking a
column header.

Edit `static/style.css` / `static/app.js` for design changes, `theme_tokens.py`
for palette/template changes, or `build.py` for what data gets pulled from
`results/` — then re-run the build command.
