import { useState } from "react";
import { base, binomSE, byDistance, crossing, pct, PCT, schemeName, type Series } from "../data";
import { axis, baseLayout, DASH, distColor, MARK, store } from "../theme";
import { Plot, type Trace } from "./Plot";

const KEY = "floquet-series";

/** Logical vs physical error rate for every scheme and run, with each one's estimated crossing. */
export function LerSection({ ler }: { ler: Series[] }) {
  const [off, setOff] = useState(() => new Set<string>(JSON.parse(store.get(KEY) || "[]")));
  const toggle = (key: string) => {
    const next = new Set(off);
    if (next.has(key)) next.delete(key); else next.add(key);
    store.set(KEY, JSON.stringify([...next]));
    setOff(next);
  };

  const traces: Trace[] = ler.flatMap((s, k) => off.has(s.key) ? [] :
    Object.entries(byDistance(s.rows)).map(([d, rows]) => ({
      type: "scatter", mode: "lines+markers", name: `${s.label} · d=${d}`,
      x: rows.map(r => r.p), y: rows.map(r => r.ler > 0 ? r.ler : null),
      error_y: { type: "data", array: rows.map(binomSE), visible: true,
                 color: distColor(+d), thickness: 1, width: 3 },
      line: { color: distColor(+d), width: 2, dash: DASH[k % DASH.length] },
      marker: { size: 6, color: distColor(+d), symbol: MARK[k % MARK.length] },
      hovertemplate: `p %{x:.2e}<br>LER %{y:.3g}<extra>${s.label} d=${d}</extra>`,
    })));

  return (
    <>
      <h2>Logical vs physical error rate</h2>
      <div className="controls">
        <div className="field" style={{ flex: "1 1 100%" }}>
          <span className="lbl">Series (scheme · run)</span>
          <div className="checks">
            {ler.map((s, i) => (
              <label key={s.key}>
                <input type="checkbox" checked={!off.has(s.key)} onChange={() => toggle(s.key)} />{" "}
                {s.label} <span style={{ color: "var(--text-muted)" }}>({DASH[i % DASH.length]})</span>
              </label>
            ))}
          </div>
        </div>
      </div>
      <Plot data={traces}
        msg={!ler.length && <span>No LER CSVs yet — run <code>make physical-to-logical</code> or <code>make plots</code>.</span>}
        layout={{ ...baseLayout(), showlegend: false,
                  xaxis: axis("Physical error rate p", { type: "log" }),
                  yaxis: axis("Logical error rate (d rounds)", { type: "log" }) }} />
      <p className="note">
        x is p, the two-qubit gate error; the noise model scales the rest off it (spin: readout 5p, reset 2p,
        idling per round 2p, 1-qubit gates p/10). Colour is the code distance, line style and symbol the scheme.
        Error bars are binomial; points with zero logical errors are left out (they sit below 1/shots). Curves
        of one scheme cross at its threshold.
      </p>
      <div className="scroll"><table><tbody>
        <tr><th>Scheme</th><th>Run</th><th>Distances</th><th>Crossing p (smallest vs largest d)</th><th>Bracket</th></tr>
        {ler.map(s => {
          const c = crossing(s.rows);
          return (
            <tr key={s.key}>
              <td>{schemeName(s.scheme)}</td><td>{base(s.run)}</td>
              <td className="num">{c ? c.ds.join(", ") : "—"}</td>
              <td className="num">{c?.p ? pct(c.p) : "no crossing"}</td>
              <td className="num">{c?.p ? `${(c.lo * PCT).toFixed(3)}–${pct(c.hi)}` : "—"}</td>
            </tr>
          );
        })}
      </tbody></table></div>
    </>
  );
}
