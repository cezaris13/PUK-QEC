import { useState } from "react";
import { runLabel, schemeName, type BiasRow } from "../data";
import { axis, baseLayout, css, isDark } from "../theme";
import { Plot, type Trace } from "./Plot";

// x3z3_bias.py's look: colour = scheme (bright X3Z3, pale CSS), line style and marker = scheme
const COLOUR: Record<string, [string, string]> = {
  pairs: ["#1f77b4", "#8fb4e8"], method_a: ["#d62728", "#f0928f"], method_b: ["#2ca02c", "#8fce8f"],
  method_c: ["#9467bd", "#c5b0d5"], pairs_bell: ["#ff7f0e", "#ffbb78"],
};
const DASH: Record<string, string> = { pairs: "solid", method_a: "dash", method_b: "dashdot", method_c: "dot",
                                       pairs_bell: "longdashdot" };
const MARKER: Record<string, string> = { pairs: "circle", method_a: "diamond", method_b: "triangle-up",
                                         method_c: "triangle-down", pairs_bell: "square" };
const CODES = ["x3z3", "css"], CODE_LABEL: Record<string, string> = { x3z3: "X³Z³", css: "CSS" };
const Z95 = 1.96;

const rate = (r: BiasRow) => r.errors / r.shots;
const either = (a: number, b: number) => 1 - (1 - a) * (1 - b);

/** 95% Wilson interval of errors / shots. */
// ponytail: Wilson, not x3z3_bias.py's Clopper-Pearson; they agree to a few % once a point has tens of errors
function wilson(r: BiasRow): [number, number] {
  const n = r.shots, f = rate(r), z2 = Z95 * Z95;
  const mid = (f + z2 / (2 * n)) / (1 + z2 / n);
  const half = Z95 * Math.sqrt(f * (1 - f) / n + z2 / (4 * n * n)) / (1 + z2 / n);
  return [Math.max(0, mid - half), mid + half];
}

/** x3z3_bias._threshold: the p where log p_L of the largest d first rises above the smallest's, interpolated in
 *  log p, the lower of the Z and X memories', with its 95% interval from each point's binomial sigma(log p_L). */
function threshold(rows: BiasRow[]): [number, number, number] | null {
  const found: [number, number][] = [];
  for (const basis of "ZX") {
    const mine = rows.filter(r => r.basis === basis && r.errors > 0);
    const ds = [...new Set(mine.map(r => r.distance))].sort((a, b) => a - b);
    if (ds.length < 2) continue;
    const at = new Map(mine.map(r => [`${r.distance}|${r.p}`, r]));
    const get = (d: number, p: number) => at.get(`${d}|${p}`)!;
    const ps = [...new Set(mine.map(r => r.p))].filter(p => at.has(`${ds[0]}|${p}`) && at.has(`${ds[ds.length - 1]}|${p}`))
      .sort((a, b) => a - b);
    const hi = ds[ds.length - 1], lo = ds[0];
    const gap = ps.map(p => Math.log(rate(get(hi, p))) - Math.log(rate(get(lo, p))));
    const v = (p: number) => [hi, lo].reduce((s, d) => s + (1 - rate(get(d, p))) / get(d, p).errors, 0);
    for (let i = 0; i + 1 < ps.length; i++) {
      if (gap[i] < 0 && gap[i + 1] >= 0) {
        const g0 = gap[i], g1 = gap[i + 1], span = Math.log(ps[i + 1]) - Math.log(ps[i]);
        const sigma = Z95 * span * Math.hypot(g1 * Math.sqrt(v(ps[i])), g0 * Math.sqrt(v(ps[i + 1]))) / (g0 - g1) ** 2;
        found.push([Math.log(ps[i]) + g0 / (g0 - g1) * span, sigma]);
        break;
      }
    }
  }
  if (!found.length) return null;
  const [m, s] = found.reduce((a, b) => (b[0] < a[0] ? b : a));
  return [Math.exp(m), Math.exp(m - s), Math.exp(m + s)];
}

const rgba = (hex: string, a: number) =>
  `rgba(${[1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16)).join(",")},${a})`;

/** A line with markers and its shaded band, as matplotlib's plot + fill_between. */
function banded(panel: 1 | 2, name: string, xs: number[], ys: (number | null)[], los: (number | null)[],
                his: (number | null)[], code: string, scheme: string, hover: string, inLegend: boolean): Trace[] {
  const colour = COLOUR[scheme]?.[code === "css" ? 1 : 0] ?? "#888888";
  const at = { xaxis: panel === 1 ? "x" : "x2", yaxis: panel === 1 ? "y" : "y2", legendgroup: scheme };
  const edge = { ...at, x: xs, mode: "lines", line: { width: 0 }, hoverinfo: "skip", showlegend: false,
                 type: "scatter", connectgaps: false };
  return [
    { ...edge, y: los },
    { ...edge, y: his, fill: "tonexty", fillcolor: rgba(colour, isDark() ? 0.2 : 0.13) },
    { ...at, type: "scatter", mode: "lines+markers", x: xs, y: ys, connectgaps: false,
      // the legend shows the scheme once (bright); clicking it hides both codes
      name: schemeName(scheme), showlegend: inLegend,
      line: { color: colour, width: 2, dash: DASH[scheme] ?? "solid" },
      marker: { color: code === "css" ? css("--surface-1") : colour, size: 8, symbol: MARKER[scheme] ?? "circle",
                line: { color: colour, width: 1.5 } },
      hovertemplate: hover + `<extra>${name}</extra>` },
  ];
}

/** The Codes legend: bright filled X3Z3, pale hollow CSS, as x3z3_bias.py's second legend block. */
const codeKey = (code: string): Trace => {
  const grey = code === "css" ? "#b3b3b3" : "#404040";
  return { type: "scatter", mode: "lines+markers", x: [null], y: [null], legend: "legend2", hoverinfo: "skip",
           name: `${CODE_LABEL[code]} (${code === "css" ? "pale" : "bright"})`,
           line: { color: grey, width: 2 },
           marker: { size: 8, color: code === "css" ? css("--surface-1") : grey, line: { color: grey, width: 1.5 } } };
};

/** x3z3_bias.py's figure: (a) threshold vs eta and (b) p_L vs eta at one p for the largest distance, every
 *  (code, scheme) of one sweep, with 95% bands. */
export function BiasSection({ bias }: { bias: Record<string, BiasRow[]> }) {
  const paths = Object.keys(bias);
  const [path, setPath] = useState(paths[0] ?? "");
  const rows = bias[path] ?? [];
  const ps = [...new Set(rows.map(r => r.p))].sort((a, b) => a - b);
  const [pPick, setP] = useState<number | null>(null);
  // default: the p closest to 1e-3, x3z3_bias.py's --p-fixed
  const p = pPick !== null && ps.includes(pPick) ? pPick
    : ps.reduce((a, b) => (Math.abs(Math.log(b / 1e-3)) < Math.abs(Math.log(a / 1e-3)) ? b : a), ps[0] ?? 1e-3);
  const etas = [...new Set(rows.map(r => r.eta))].sort((a, b) => a - b);
  const dMax = Math.max(...rows.map(r => r.distance));
  const schemes = [...new Set(rows.map(r => r.scheme))];
  const combos = CODES.flatMap(c => schemes.filter(s => rows.some(r => r.code === c && r.scheme === s)).map(s => [c, s]));

  const traces: Trace[] = [];
  const shown = new Set<string>();
  for (const [code, scheme] of combos) {
    const first = !shown.has(scheme);
    shown.add(scheme);
    const mine = rows.filter(r => r.code === code && r.scheme === scheme);
    const name = `${CODE_LABEL[code] ?? code} · ${schemeName(scheme)}`;
    const th = etas.map(e => threshold(mine.filter(r => r.eta === e)));
    traces.push(...banded(1, name, etas, th.map(t => t?.[0] ?? null), th.map(t => t?.[1] ?? null),
                          th.map(t => t?.[2] ?? null), code, scheme, "η %{x:g}<br>p_th %{y:.3e}", first));
    const pl = etas.map(e => {
      const at = mine.filter(r => r.eta === e && r.p === p && r.distance === dMax);
      const z = at.find(r => r.basis === "Z"), x = at.find(r => r.basis === "X");
      if (!z || !x) return null;
      const [zl, zh] = wilson(z), [xl, xh] = wilson(x);
      return [either(rate(z), rate(x)), either(zl, xl), either(zh, xh)];
    });
    traces.push(...banded(2, name, etas, pl.map(v => v?.[0] || null), pl.map(v => v?.[1] || null),
                          pl.map(v => v?.[2] ?? null), code, scheme, "η %{x:g}<br>p_L %{y:.3e}", false));
  }
  traces.push(...CODES.filter(c => combos.some(([x]) => x === c)).map(codeKey));

  const none = !paths.length && (
    <span>No bias sweeps yet. <code>make x3z3-bias</code> writes one to <code>results/x3z3/</code>.</span>
  );
  const L = baseLayout(), text = css("--text-primary");
  const legendBox = (y: number, title: string) => ({ ...L.legend, x: 1.02, xanchor: "left", y, yanchor: "top",
    title: { text: `<b>${title}</b>`, font: { color: text, size: 12 } }, font: { color: text, size: 12 },
    tracegroupgap: 2 });
  const logX = axis("noise bias η", { type: "log", dtick: 1, exponentformat: "power", showline: true,
                                      linecolor: css("--text-secondary"), ticks: "outside" });
  const panel = (t: string, x: number) => ({ text: t, x, xref: "paper", xanchor: "center", y: 1.04, yref: "paper",
                                            yanchor: "bottom", showarrow: false, font: { size: 15, color: text } });
  return (
    <>
      <h2>X³Z³ vs CSS against noise bias</h2>
      <div className="controls">
        <div className="field">
          <label htmlFor="biasRun">Sweep</label>
          <select id="biasRun" value={path} onChange={e => setPath(e.target.value)}>
            {paths.map(f => <option key={f} value={f}>{runLabel(f)}</option>)}
          </select>
        </div>
        <div className="field">
          <label htmlFor="biasP">p of (b)</label>
          <select id="biasP" value={p} onChange={e => setP(+e.target.value)}>
            {ps.map(v => <option key={v} value={v}>{v.toPrecision(3)}</option>)}
          </select>
        </div>
      </div>
      <Plot data={traces} msg={none}
        layout={{ ...L, font: { ...L.font, size: 13 }, margin: { l: 80, r: 190, t: 50, b: 60 },
          legend: legendBox(1, "Schemes"), legend2: legendBox(0.45, "Codes"),
          annotations: [panel("(a) threshold", 0.21),
                        panel(`(b) d = ${dMax}, ${dMax} rounds, p = ${p.toPrecision(3)}`, 0.79)],
          xaxis: { ...logX, domain: [0, 0.42] }, xaxis2: { ...logX, domain: [0.58, 1], anchor: "y2" },
          yaxis: axis("threshold p<sub>th</sub>", { tickformat: ".4~f", showline: true, ticks: "outside",
                                                   linecolor: css("--text-secondary") }),
          yaxis2: axis("logical failure probability p<sub>L</sub>", { type: "log", anchor: "x2", dtick: 1,
                       exponentformat: "power", showline: true, ticks: "outside", linecolor: css("--text-secondary") }) }} />
      <p className="note">
        As <code>x3z3_bias.py</code>'s figure. Colour, line style and marker are the scheme; X³Z³ is bright and
        filled, CSS pale and hollow. (a) Where the largest distance's p<sub>L</sub> curve crosses the smallest's,
        the lower of the Z and X memories'. Band: 95% from each point's shot noise carried through the crossing.
        (b) p<sub>L</sub> = 1 − (1 − p<sub>Z</sub>)(1 − p<sub>X</sub>) at the largest distance. Band: 95% Wilson
        intervals (close to the script's Clopper-Pearson).
      </p>
    </>
  );
}
