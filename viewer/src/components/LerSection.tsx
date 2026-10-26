import { useState } from "react";
import { crossing, dLabel, pct, PCT, roundsLabel, schemeName, type Row, type Series } from "../data";
import { axis, baseLayout, css, distColor, store } from "../theme";
import { combined, style, type FootRow, type Point } from "./FootSection";
import { Plot, type Trace } from "./Plot";

const KEY = "floquet-ler-run", COMBINE = "floquet-ler-combine";

/** physical_to_logical.py's plot for one run: colour = distance, marker and dash = scheme, CSS hollow and paler
 *  when both codes are in the run, Z and X memory combined (1 - (1 - p_Z)(1 - p_X)) or apart (Z solid, X dotted). */
export function LerSection({ ler }: { ler: Series[] }) {
  const [key, setKey] = useState(() => store.get(KEY) || "");
  const [combine, setCombine] = useState(() => store.get(COMBINE) !== "0");
  const s = ler.find(x => x.key === key) ?? ler[0];
  const pts = s ? combined(s.rows as FootRow[], combine) : [];
  const codes = [...new Set(pts.map(p => p.code))], both = codes.length > 1;
  const bases = [...new Set(pts.map(p => p.basis))], apart = bases.length > 1;

  const curves: Record<string, Point[]> = {};
  pts.forEach(p => (curves[`${p.code}|${p.scheme}|${p.basis}`] ??= []).push(p));
  const label = (c: Point) => [both && c.code, schemeName(c.scheme), apart && `${c.basis} memory`].filter(Boolean).join(", ");

  const traces: Trace[] = Object.values(curves).flatMap(curve => {
    const c = curve[0], hollow = both && c.code === "css", [dash, symbol] = style(c.scheme);
    const byD: Record<number, Point[]> = {};
    curve.forEach(p => (byD[p.distance] ??= []).push(p));
    return Object.entries(byD).map(([d, all]) => {
      const seen = all.filter(p => p.errors > 0).sort((a, b) => a.p - b.p);  // dropped: a null y still stretches x
      return { type: "scatter", mode: "lines+markers", showlegend: false, opacity: hollow ? 0.55 : 1,
        x: seen.map(p => p.p), y: seen.map(p => p.ler),
        error_y: { type: "data", symmetric: false, array: seen.map(p => p.hi - p.ler),
                   arrayminus: seen.map(p => p.ler - p.lo), color: distColor(+d), thickness: 1, width: 3 },
        line: { color: distColor(+d), width: 2, dash: apart ? (c.basis === "X" ? "dot" : "solid") : dash },
        marker: { size: 7, color: distColor(+d), symbol: hollow ? `${symbol}-open` : symbol, line: { width: hollow ? 1.5 : 0 } },
        hovertemplate: `p %{x:.2e}<br>LER %{y:.3g}<extra>${label(c)}, d=${dLabel(+d)}</extra>` };
    });
  });
  const legend = (group: string, extra: Trace) =>
    traces.push({ type: "scatter", x: [null], y: [null], legendgroup: group, legendgrouptitle: { text: group }, ...extra });
  const grey = css("--text-secondary");
  [...new Set(pts.map(p => p.scheme))].forEach(sc => legend("Schemes", { mode: apart ? "markers" : "lines+markers",
    name: schemeName(sc), line: { color: grey, dash: style(sc)[0] }, marker: { symbol: style(sc)[1], size: 9, color: grey } }));
  [...new Set(pts.map(p => p.distance))].sort((a, b) => a - b).forEach(d =>
    legend("Distances", { mode: "markers", name: String(dLabel(d)), marker: { symbol: "square", size: 9, color: distColor(d) } }));
  if (both) codes.forEach(code => legend("Codes", { mode: "markers", name: code,
    marker: { symbol: code === "css" ? "circle-open" : "circle", size: 9, color: grey } }));
  if (apart) bases.forEach(b => legend("Memory", { mode: "lines", name: b,
    line: { color: grey, dash: b === "X" ? "dot" : "solid" } }));

  const r0 = s?.rows[0];
  return (
    <>
      <h2>Logical vs physical error rate</h2>
      <div className="controls">
        <div className="field">
          <label htmlFor="lerSel">Run</label>
          <select id="lerSel" value={s?.key ?? ""} onChange={e => { setKey(e.target.value); store.set(KEY, e.target.value); }}>
            {ler.map(x => <option key={x.key} value={x.key}>{x.label}</option>)}
          </select>
        </div>
        <div className="field">
          <span className="lbl">Memory bases</span>
          <label><input type="checkbox" checked={combine}
            onChange={e => { setCombine(e.target.checked); store.set(COMBINE, e.target.checked ? "1" : "0"); }} />{" "}
            combine Z and X into one P<sub>L</sub></label>
        </div>
      </div>
      <Plot data={traces}
        msg={!ler.length && <span>No LER CSVs yet — run <code>make physical-to-logical</code> or <code>make plots</code>.</span>}
        layout={{ ...baseLayout(), title: { text: r0 ? `${r0.noise} noise, η = ${r0.eta}` : "" },
                  legend: { ...baseLayout().legend, groupclick: "toggleitem" },
                  xaxis: axis("Physical error rate p", { type: "log" }),
                  yaxis: axis(`Logical error rate (${roundsLabel()}, ${combine ? "Z or X" : "per"} memory)`, { type: "log" }) }} />
      <p className="note">
        x is p, the two-qubit gate error; the noise model scales the rest off it (spin: readout 5p, reset 2p,
        idling per round 2p, 1-qubit gates p/10). Error bars are 95% Wilson intervals; points with zero logical
        errors are left out. Curves of one scheme cross at its threshold.
      </p>
      <div className="scroll"><table><tbody>
        <tr><th>Curve</th><th>Distances</th><th>Crossing p (smallest vs largest d)</th><th>Bracket</th></tr>
        {Object.entries(curves).map(([k, curve]) => {
          const c = crossing(curve as unknown as Row[]);
          return (
            <tr key={k}>
              <td>{label(curve[0])}</td>
              <td className="num">{c ? c.ds.map(dLabel).join(", ") : "—"}</td>
              <td className="num">{c?.p ? pct(c.p) : "no crossing"}</td>
              <td className="num">{c?.p ? `${(c.lo * PCT).toFixed(3)}–${pct(c.hi)}` : "—"}</td>
            </tr>
          );
        })}
      </tbody></table></div>
    </>
  );
}
