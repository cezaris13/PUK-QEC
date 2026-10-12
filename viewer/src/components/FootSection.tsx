import { useState } from "react";
import { runLabel, schemeName, type Row } from "../data";
import { axis, baseLayout, css, distColor } from "../theme";
import { Plot, type Trace } from "./Plot";

// footprint.py's STYLES: line dash and marker per scheme
const STYLE: Record<string, [string, string]> = {
  pairs: ["solid", "circle"], method_a: ["dash", "diamond"], method_b: ["dashdot", "triangle-up"],
  method_c: ["dot", "triangle-down"], pairs_bell: ["longdashdot", "square"],
};
const style = (s: string) => STYLE[s] ?? ["solid", "circle"];

type FootRow = Row & { code?: string; basis?: string };
type Point = { code: string; scheme: string; distance: number; qubits: number; errors: number;
               ler: number; lo: number; hi: number };

// ponytail: Wilson 95% interval, close to footprint.py's Clopper-Pearson without a beta quantile in JS
function wilson(k: number, n: number): [number, number] {
  const z = 1.96, p = k / n, den = 1 + z * z / n;
  const mid = (p + z * z / (2 * n)) / den, half = z * Math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den;
  return [Math.max(0, mid - half), mid + half];
}

/** footprint.py's _combined: one point per (code, scheme, d), LER = 1 - (1 - p_Z)(1 - p_X). Older CSVs with no
 *  basis column are one point per row already. */
function combined(rows: FootRow[]): Point[] {
  const either = (a: number, b: number) => 1 - (1 - a) * (1 - b);
  const by: Record<string, FootRow[]> = {};
  rows.forEach(r => (by[`${r.code ?? "css"}|${r.scheme}|${r.distance}`] ??= []).push(r));
  return Object.values(by).flatMap(g => {
    const z = g.find(r => r.basis === "Z"), x = g.find(r => r.basis === "X");
    const parts = z && x ? [z, x] : g[0].basis ? [] : [g[0]];
    if (!parts.length) return [];
    const ivs = parts.map(r => wilson(r.errors, r.shots));
    return [{ code: g[0].code ?? "css", scheme: g[0].scheme!, distance: g[0].distance, qubits: g[0].qubits!,
              errors: parts.reduce((s, r) => s + r.errors, 0), ler: parts.map(r => r.ler).reduce(either, 0),
              lo: ivs.map(i => i[0]).reduce(either, 0), hi: ivs.map(i => i[1]).reduce(either, 0) }];
  });
}

/** ln LER = a + b d through the points with errors, and where it meets target: footprint.py's _qubits_needed. */
function fit(pts: Point[], target: number) {
  const ok = pts.filter(p => p.errors > 0);
  if (new Set(ok.map(p => p.distance)).size < 2) return null;
  const xs = ok.map(p => p.distance), ys = ok.map(p => Math.log(p.ler)), n = xs.length;
  const mx = xs.reduce((a, b) => a + b) / n, my = ys.reduce((a, b) => a + b) / n;
  const b = xs.reduce((s, x, i) => s + (x - mx) * (ys[i] - my), 0) / xs.reduce((s, x) => s + (x - mx) ** 2, 0);
  if (b >= 0) return null;
  const a = my - b * mx, d = (Math.log(target) - a) / b, nu = pts[0].qubits / pts[0].distance ** 2;
  return { a, b, d, nu, qubits: nu * d * d, lambda: Math.exp(-2 * b) };
}

const NICE = [100, 200, 500, 1000, 1500, 2000, 3000, 4500, 6000, 8000, 10000, 15000, 20000, 30000, 50000];

/** footprint.py's plot: LER vs qubits on a square-root axis, colour = distance, marker = scheme, CSS hollow when
 *  both codes are in the run, grey fits out to the target (star). */
export function FootSection({ foot }: { foot: Record<string, Row[]> }) {
  const paths = Object.keys(foot);
  const [path, setPath] = useState(paths[0] ?? "");
  const [target, setTarget] = useState(1e-8);
  const [byQ, setByQ] = useState(true);
  const pts = combined((foot[path] ?? []) as FootRow[]);
  const codes = [...new Set(pts.map(p => p.code))], schemes = [...new Set(pts.map(p => p.scheme))];
  const both = codes.length > 1;
  const X = (q: number, d: number) => byQ ? Math.sqrt(q) : d;
  const traces: Trace[] = [];
  const legend = (group: string, extra: Trace) =>
    traces.push({ type: "scatter", x: [null], y: [null], legendgroup: group, legendgrouptitle: { text: group }, ...extra });

  const fits = codes.flatMap(code => schemes.map(s => {
    const mine = pts.filter(p => p.code === code && p.scheme === s).sort((a, b) => a.distance - b.distance);
    if (!mine.length) return null;
    const [dash, symbol] = style(s), hollow = both && code === "css", seen = mine.filter(p => p.errors > 0);
    traces.push({ type: "scatter", mode: "markers", showlegend: false,
      x: seen.map(p => X(p.qubits, p.distance)), y: seen.map(p => p.ler),
      error_y: { type: "data", symmetric: false, array: seen.map(p => p.hi - p.ler),
                 arrayminus: seen.map(p => p.ler - p.lo), thickness: 1, width: 3 },
      marker: { size: 9, symbol: hollow ? `${symbol}-open` : symbol, color: seen.map(p => distColor(p.distance)),
                line: { width: hollow ? 1.5 : 0 } },
      text: seen.map(p => `${both ? code + ", " : ""}${schemeName(s)}, d=${p.distance}<br>${p.qubits} qubits`),
      hovertemplate: "%{text}<br>LER %{y:.3g}<extra></extra>" });
    const f = fit(mine, target);
    if (f) {
      const grey = css(hollow ? "--text-muted" : "--text-secondary"), x: number[] = [], y: number[] = [];
      for (let i = 0; i <= 50; i++) {
        const d = mine[0].distance + (f.d - mine[0].distance) * i / 50;
        x.push(X(f.nu * d * d, d)); y.push(Math.exp(f.a + f.b * d));
      }
      const label = `${both ? code + ", " : ""}${schemeName(s)}: ${Math.round(f.qubits).toLocaleString()} qubits, Λ = ${f.lambda.toFixed(1)}`;
      traces.push({ type: "scatter", mode: "lines", x, y, name: label, legendgroup: "Linefits",
                    legendgrouptitle: { text: "Linefits" }, hoverinfo: "skip", line: { color: grey, width: 1, dash } });
      traces.push({ type: "scatter", mode: "markers", x: [X(f.qubits, f.d)], y: [target], showlegend: false,
                    hovertemplate: `${label}<extra></extra>`, marker: { symbol: "star", size: 12, color: grey } });
    }
    return { code, s, f, nu: mine[0].qubits / mine[0].distance ** 2 };
  })).filter(Boolean) as { code: string; s: string; f: ReturnType<typeof fit>; nu: number }[];

  schemes.forEach(s => legend("Schemes", { mode: "markers", name: schemeName(s),
    marker: { symbol: style(s)[1], size: 9, color: css("--text-secondary") } }));
  [...new Set(pts.map(p => p.distance))].sort((a, b) => a - b).forEach(d =>
    legend("Distances", { mode: "markers", name: String(d), marker: { symbol: "square", size: 9, color: distColor(d) } }));

  const allQ = [...pts.map(p => p.qubits), ...fits.map(r => r.f?.qubits ?? 0)];
  const ticks = NICE.filter(q => q <= Math.max(...allQ, 1) * 1.05);
  const r0 = foot[path]?.[0];

  return (
    <>
      <h2>Qubits for a target logical error rate</h2>
      <div className="controls">
        <div className="field">
          <label htmlFor="footSel">Footprint run</label>
          <select id="footSel" value={path} onChange={e => setPath(e.target.value)}>
            {paths.map(p => <option key={p} value={p}>{runLabel(p)}</option>)}
          </select>
        </div>
        <div className="field">
          <label htmlFor="target">Target logical error rate</label>
          <select id="target" value={target} onChange={e => setTarget(+e.target.value)}>
            {([[1e-6, "10⁻⁶"], [1e-8, "10⁻⁸"], [1e-9, "10⁻⁹"], [1e-12, "10⁻¹²"]] as const).map(([v, t]) =>
              <option key={v} value={v}>{t}</option>)}
          </select>
        </div>
        <div className="field">
          <label htmlFor="xaxis">Horizontal axis</label>
          <select id="xaxis" value={byQ ? "q" : "d"} onChange={e => setByQ(e.target.value === "q")}>
            <option value="q">Physical qubits</option><option value="d">Code distance</option>
          </select>
        </div>
      </div>
      <Plot data={traces} tall
        msg={!paths.length && <span>No footprint runs yet — <code>make footprint</code> writes one.</span>}
        layout={{ ...baseLayout(),
          title: { text: r0 ? `${r0.noise} noise, η = ${r0.eta}, p = ${r0.p}; dotted: target ${target}` : "",
                   font: { color: css("--text-primary"), size: 13 }, x: 0, xanchor: "left" },
          legend: { ...baseLayout().legend, groupclick: "toggleitem", font: { size: 11, color: css("--text-primary") } },
          shapes: [{ type: "line", xref: "paper", x0: 0, x1: 1, y0: target, y1: target,
                     line: { color: css("--text-muted"), width: 1, dash: "dot" } }],
          xaxis: byQ
            ? axis("Total physical qubits (√ scale)", { tickvals: ticks.map(Math.sqrt), ticktext: ticks.map(String),
                                                        range: [Math.sqrt(0.8 * Math.min(...allQ)), Math.sqrt(Math.max(...allQ) * 1.05)] })
            : axis("Code distance d", { dtick: 2 }),
          yaxis: axis("Logical error rate, Z or X memory (d rounds)", { type: "log", exponentformat: "power" }) }} />
      <p className="note">
        Z and X memory combined, p_L = 1 − (1 − p_Z)(1 − p_X), with 95% Wilson bars. Points with no logical errors are
        left out. Each fit is ln LER = a + b·d, a straight line on the √ qubit axis because N = ν d²; Λ = LER(d) / LER(d + 2).
        {both && " Hollow markers and light fits: CSS; filled and dark: X3Z3."}
      </p>
      <div className="scroll"><table><tbody>
        <tr><th>Code</th><th>Scheme</th><th>ν = qubits / d²</th><th>Λ</th><th>d for target</th><th>Qubits for target</th></tr>
        {fits.map(({ code, s, f, nu }) => (
          <tr key={code + s}>
            <td>{code}</td><td>{schemeName(s)}</td><td className="num">{nu.toFixed(1)}</td>
            <td className="num">{f ? f.lambda.toFixed(1) : "no fit"}</td>
            <td className="num">{f ? f.d.toFixed(1) : "—"}</td>
            <td className="num">{f ? Math.round(f.qubits).toLocaleString() : "—"}</td>
          </tr>
        ))}
      </tbody></table></div>
    </>
  );
}
