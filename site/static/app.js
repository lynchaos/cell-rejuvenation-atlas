(function () {
  "use strict";

  var dataEl = document.getElementById("site-data");
  var DATA = JSON.parse(dataEl.textContent);
  var root = document.documentElement;

  /* ---------------------------------------------------------------
     Theme: explicit toggle wins over prefers-color-scheme, persisted.
  --------------------------------------------------------------- */
  var toggleBtn = document.getElementById("theme-toggle");
  var stored = localStorage.getItem("cra-theme");
  var prefersDark = window.matchMedia("(prefers-color-scheme: dark)").matches;
  var theme = stored || (prefersDark ? "dark" : "light");

  function setToggleLabel() {
    toggleBtn.textContent = theme === "dark" ? "☀ Light" : "● Dark";
  }

  function applyTheme(next, rerender) {
    theme = next;
    root.setAttribute("data-theme", theme);
    localStorage.setItem("cra-theme", theme);
    setToggleLabel();
    if (rerender) renderActiveModuleCharts();
  }

  root.setAttribute("data-theme", theme);
  setToggleLabel();
  toggleBtn.addEventListener("click", function () {
    applyTheme(theme === "dark" ? "light" : "dark", true);
  });

  /* ---------------------------------------------------------------
     Module navigation (single-page, section show/hide)
  --------------------------------------------------------------- */
  var navButtons = Array.prototype.slice.call(document.querySelectorAll(".nav-item"));
  var sections = Array.prototype.slice.call(document.querySelectorAll(".module-section"));
  var renderedOnce = {};

  function showModule(id) {
    sections.forEach(function (s) {
      s.classList.toggle("is-active", s.getAttribute("data-module") === id);
    });
    navButtons.forEach(function (b) {
      b.classList.toggle("is-active", b.getAttribute("data-module") === id);
    });
    renderModuleCharts(id);
    if (history.replaceState) history.replaceState(null, "", "#module-" + id);
  }

  navButtons.forEach(function (btn) {
    btn.addEventListener("click", function () {
      showModule(btn.getAttribute("data-module"));
    });
  });

  /* ---------------------------------------------------------------
     Plotly rendering — theme-aware. Each chart object carries its own
     data/layout plus optional per-theme color / colorscale overrides,
     merged with the shared light/dark template at render time.
  --------------------------------------------------------------- */
  function mergeLayout(base, override) {
    var out = JSON.parse(JSON.stringify(base));
    for (var k in override || {}) {
      if (out[k] && typeof out[k] === "object" && !Array.isArray(out[k]) &&
          override[k] && typeof override[k] === "object" && !Array.isArray(override[k])) {
        out[k] = Object.assign({}, out[k], override[k]);
      } else {
        out[k] = override[k];
      }
    }
    return out;
  }

  function paintChart(containerId, chart) {
    var el = document.getElementById(containerId);
    if (!el || typeof Plotly === "undefined") return;
    var traces = (chart.data || []).map(function (tr) {
      return JSON.parse(JSON.stringify(tr));
    });
    if (chart.color) {
      var c = chart.color[theme] || chart.color.light;
      traces.forEach(function (tr) {
        if (tr.marker) tr.marker = Object.assign({}, tr.marker, { color: c });
        if (tr.line) tr.line = Object.assign({}, tr.line, { color: c });
        if (!tr.marker && !tr.line) tr.marker = { color: c };
      });
    }
    if (chart.colorscale) {
      var cs = chart.colorscale[theme] || chart.colorscale.light;
      traces.forEach(function (tr) { tr.colorscale = cs; });
    }
    var template = DATA.templates[theme] || DATA.templates.light;
    var layout = mergeLayout(template.layout, chart.layout || {});
    Plotly.react(el, traces, layout, { displayModeBar: false, responsive: true });
  }

  function renderModuleCharts(id) {
    var mod = DATA.modules[id];
    if (!mod) return;
    (mod.charts || []).forEach(function (chart) { paintChart(chart.id, chart); });
  }

  function renderActiveModuleCharts() {
    var active = sections.filter(function (s) { return s.classList.contains("is-active"); })[0];
    if (active) renderModuleCharts(active.getAttribute("data-module"));
  }

  /* ---------------------------------------------------------------
     Variant selectors (e.g. module 4 age-group heatmap)
  --------------------------------------------------------------- */
  Array.prototype.slice.call(document.querySelectorAll(".variant-select")).forEach(function (sel) {
    sel.addEventListener("change", function () {
      var modId = sel.getAttribute("data-module");
      var chartId = sel.getAttribute("data-chart");
      var variant = DATA.modules[modId].variants[sel.value];
      if (variant) paintChart(chartId, variant);
    });
  });

  /* ---------------------------------------------------------------
     Click-to-sort tables
  --------------------------------------------------------------- */
  Array.prototype.slice.call(document.querySelectorAll("table.data-table")).forEach(function (table) {
    var headers = Array.prototype.slice.call(table.querySelectorAll("th"));
    headers.forEach(function (th, colIdx) {
      var ascending = true;
      th.addEventListener("click", function () {
        var tbody = table.querySelector("tbody");
        var rows = Array.prototype.slice.call(tbody.querySelectorAll("tr"));
        rows.sort(function (a, b) {
          var av = a.children[colIdx].textContent.trim();
          var bv = b.children[colIdx].textContent.trim();
          var an = parseFloat(av), bn = parseFloat(bv);
          var cmp = (!isNaN(an) && !isNaN(bn)) ? an - bn : av.localeCompare(bv);
          return ascending ? cmp : -cmp;
        });
        rows.forEach(function (r) { tbody.appendChild(r); });
        headers.forEach(function (h) { h.classList.remove("sorted-asc", "sorted-desc"); });
        th.classList.add(ascending ? "sorted-asc" : "sorted-desc");
        ascending = !ascending;
      });
    });
  });

  /* ---------------------------------------------------------------
     Initial render
  --------------------------------------------------------------- */
  var initial = (location.hash.match(/module-(\w+)/) || [])[1];
  if (!initial || !DATA.modules[initial]) initial = Object.keys(DATA.modules)[0];
  showModule(initial);
})();
