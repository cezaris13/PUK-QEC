import { useState } from "react";
import { binomSE, distanceFor, fitDecay, runLabel, schemeName, type Row } from "../data";
import { axis, baseLayout, css, MARK, SCHEME_COLORS } from "../theme";
import { Plot, type Trace } from "./Plot";

/** footprint.py's LER vs physical qubits at one p, with each scheme's fit run out to a target LER. */
export function FootSection({ foot }: { foot: Record<string, Row[]> }) {
  const paths = Object.keys(foot);
  const [path, setPath] = useState(paths[0] ?? "");
  const [target, setTarget] = useState(1e-6);
  const [byQ, setByQ] = useState(true);
  const rows = foot[path] ?? [];

  const traces: Trace[] = [];
  const fits = [...new Set(rows.map(r => r.scheme!))].map((s, k) => {
    const mine = rows.filter(r => r.scheme === s).sort((a, b) => a.distance - b.distance);
    const color = SCHEME_COLORS[k % SCHEME_COLORS.length], nu = mine[0].qubits! / mine[0].distance ** 2;
    const seen = mine.filter(r => r.ler > 0);
    traces.push({ type: "scatter", mode: "markers", name: schemeName(s),
      x: seen.map(r => byQ ? r.qubits : r.distance), y: seen.map(r => r.ler),
      error_y: { type: "data", array: seen.map(binomSE), visible: true, color, thickness: 1 },
      marker: { size: 10, color, symbol: MARK[k % MARK.length] }, text: seen.map(r => `d=${r.distance}`),
      hovertemplate: `%{text}<br>${byQ ? "%{x} qubits" : "d = %{x}"}<br>LER %{y:.3g}<extra>${schemeName(s)}</extra>` });
    const f = fitDecay(mine), dt = distanceFor(f, target);
    if (f) {
      // linear in d; on the qubit axis it bends only because N = nu d^2
      const x: number[] = [], y: number[] = [];
      for (let d = mine[0].distance; d <= Math.max(dt ?? 0, mine[mine.length - 1].distance + 2); d += 0.25) {
        x.push(byQ ? nu * d * d : d); y.push(10 ** (f.a + f.b * d));
      }
      traces.push({ type: "scatter", mode: "lines", x, y, showlegend: false, hoverinfo: "skip",
                    line: { color, width: 1, dash: "dash" } });
    }
    return { s, nu, f, dt };
  });

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
            <option value={1e-6}>10⁻⁶</option><option value={1e-9}>10⁻⁹</option><option value={1e-12}>10⁻¹²</option>
          </select>
        </div>
        <div className="field">
          <label htmlFor="xaxis">Horizontal axis</label>
          <select id="xaxis" value={byQ ? "q" : "d"} onChange={e => setByQ(e.target.value === "q")}>
            <option value="q">Physical qubits</option><option value="d">Code distance</option>
          </select>
        </div>
      </div>
      <Plot data={traces}
        msg={!paths.length && <span>No footprint runs yet — <code>make footprint</code> writes one.</span>}
        layout={{ ...baseLayout(),
          title: { text: rows.length ? `p = ${rows[0].p}` : "", font: { color: css("--text-primary"), size: 13 },
                   x: 0, xanchor: "left" },
          shapes: [{ type: "line", xref: "paper", x0: 0, x1: 1, y0: target, y1: target,
                     line: { color: css("--text-muted"), width: 1, dash: "dot" } }],
          xaxis: axis(byQ ? "Physical qubits (data + ancillas)" : "Code distance d",
                      { rangemode: "tozero", dtick: byQ ? undefined : 2 }),
          yaxis: axis("Logical error rate (d rounds)", { type: "log" }) }} />
      <p className="note">
        Each point is a direct simulation at one p. The dashed fit is log₁₀ LER = a + b·d per scheme, linear in d;
        on the qubit axis it bends only because N = ν d².
      </p>
      <div className="scroll"><table><tbody>
        <tr><th>Scheme</th><th>ν = qubits / d²</th><th>Fit log₁₀ LER</th><th>d for target</th><th>Qubits for target</th></tr>
        {fits.map(({ s, nu, f, dt }) => (
          <tr key={s}>
            <td>{schemeName(s)}</td><td className="num">{nu.toFixed(1)}</td>
            <td className="num">{f ? `${f.a.toFixed(2)} ${f.b < 0 ? "−" : "+"} ${Math.abs(f.b).toFixed(3)}·d` : "no fit"}</td>
            <td className="num">{dt ?? "—"}</td>
            <td className="num">{dt ? Math.round(nu * dt * dt).toLocaleString() : "—"}</td>
          </tr>
        ))}
      </tbody></table></div>
    </>
  );
}
