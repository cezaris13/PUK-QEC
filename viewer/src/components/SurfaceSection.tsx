import { useState } from "react";
import { pct, PCT, runLabel, schemeName, surfPoints, type ThresholdRun } from "../data";
import { axis, baseLayout, COPPER, css, G_COLOR, T_COLOR } from "../theme";
import { Plot, type Trace } from "./Plot";

type Key = "A" | "B";

/** Run A shaded by p_th, run B a flat overlay, plus the bias-sweep ranges as bars on the p_G and p_T axes. */
function surfaceTraces(path: string, d: ThresholdRun, key: Key): Trace[] {
  const pts = surfPoints(d), shaded = key === "A";
  const c = css(shaded ? "--series-a" : "--series-b");
  const P = (k: "x" | "y" | "z" | "p") => pts.map(q => q[k]);
  const out: Trace[] = [];
  if (pts.length >= 3) out.push({
    type: "mesh3d", x: P("x"), y: P("y"), z: P("z"), delaunayaxis: "z", name: runLabel(path),
    hoverinfo: "skip", flatshading: true, lighting: { ambient: 0.9, diffuse: 0.2, specular: 0 },
    ...(shaded
      ? { intensity: P("p"), colorscale: COPPER, showscale: true, opacity: 0.85,
          colorbar: { title: { text: "p_th (%)", side: "top" }, thickness: 14, len: 0.7,
                      outlinewidth: 0, tickfont: { color: css("--text-secondary") } } }
      : { color: c, opacity: 0.35 }),
  });
  out.push({
    type: "scatter3d", mode: "markers", x: P("x"), y: P("y"), z: P("z"), showlegend: false,
    customdata: pts.map(q => [q.p, q.e]),
    marker: { size: 3, color: shaded ? css("--text-primary") : c },
    hovertemplate: "p_G %{x:.3f} %<br>p_T %{y:.3f} %<br>p_R %{z:.3f} %<br>" +
                   `p_th %{customdata[0]:.3f} ± %{customdata[1]:.3f} %<extra>${runLabel(path)}</extra>`,
  });
  // like the arrows of IBM's panel (b)
  for (const [rows, along, color] of [[d.g_bias, "x", G_COLOR], [d.t_bias, "y", T_COLOR]] as const) {
    const ps = rows.filter(r => r.p_th > 0).map(r => PCT * r.p_th);
    if (ps.length < 2) continue;
    const ends = [Math.min(...ps), Math.max(...ps)];
    out.push({ type: "scatter3d", mode: "lines+markers", showlegend: false, hoverinfo: "skip",
               x: along === "x" ? ends : [0, 0], y: along === "y" ? ends : [0, 0], z: [0, 0],
               line: { color, width: 6 }, marker: { size: 3, color } });
  }
  return out;
}

function biasTraces(path: string, d: ThresholdRun, key: Key): Trace[] {
  const dash = key === "A" ? "solid" : "dash", symbol = key === "A" ? "circle" : "square";
  return ([["g_bias", "eta_g", "p_G^th vs η_G", G_COLOR], ["t_bias", "eta_t", "p_T^th vs η_T", T_COLOR]] as const)
    .map(([k, eta, name, color]) => {
      const rows = d[k].filter(r => r.p_th > 0);
      return { type: "scatter", mode: "lines+markers", name: `${name} · ${key}`,
               x: rows.map(r => r[eta]), y: rows.map(r => PCT * r.p_th),
               error_y: { type: "data", array: rows.map(r => PCT * r.p_th_error), visible: true,
                          color, thickness: 1, width: 3 },
               line: { color, dash, width: 2 }, marker: { color, symbol, size: 7 },
               hovertemplate: `η %{x:.3g}<br>p_th %{y:.3f} %<extra>${runLabel(path)}</extra>` };
    });
}

// an invisible point keeps the empty 3D box on screen until the first run exists
const EMPTY_3D: Trace[] = [{ type: "scatter3d", x: [0], y: [0], z: [0], mode: "markers",
                             marker: { opacity: 0 }, hoverinfo: "skip", showlegend: false }];

/** thresholds.py's (b) surface and (c) bias plot, two runs at a time. */
export function SurfaceSection({ surf }: { surf: Record<string, ThresholdRun> }) {
  const paths = Object.keys(surf);
  const [a, setA] = useState(paths[0] ?? "");
  const [b, setB] = useState("");
  const picked = ([["A", a], ["B", b]] as [Key, string][]).filter(([, p]) => surf[p]);
  const options = paths.map(p => <option key={p} value={p}>{runLabel(p)}</option>);
  const none = !paths.length && (
    <span>No threshold runs yet.<br /><code>make thresholds SCHEME=method_b</code> (or pairs, method_a,
      method_c) writes one to <code>results/…/thresholds/*.json</code>.</span>
  );

  return (
    <>
      <h2>Threshold surface (p<sub>G</sub>, p<sub>T</sub>, p<sub>R</sub>)</h2>
      <div className="controls">
        <div className="field">
          <label htmlFor="selA"><span className="swatch" style={{ background: "var(--series-a)" }} />Run A</label>
          <select id="selA" value={a} onChange={e => setA(e.target.value)}>{options}</select>
        </div>
        <div className="field">
          <label htmlFor="selB"><span className="swatch" style={{ background: "var(--series-b)" }} />Run B (compare)</label>
          <select id="selB" value={b} onChange={e => setB(e.target.value)}>
            <option value="">— none —</option>{options}
          </select>
        </div>
      </div>
      <Plot tall msg={none}
        data={picked.length ? picked.flatMap(([k, p]) => surfaceTraces(p, surf[p], k)) : EMPTY_3D}
        layout={{ ...baseLayout(), margin: { l: 0, r: 0, t: 0, b: 0 },
          legend: { ...baseLayout().legend, x: 0, y: 1 },
          scene: { xaxis: axis("p_G gate (%)", { rangemode: "tozero" }),
                   yaxis: axis("p_T idling (%)", { rangemode: "tozero" }),
                   zaxis: axis("p_R readout (%)", { rangemode: "tozero" }),
                   bgcolor: "rgba(0,0,0,0)", aspectmode: "cube",
                   camera: { eye: { x: 1.6, y: 1.6, z: 1.0 } } } }} />
      <p className="note">
        Each dot is one direction's crossing point p<sub>th</sub>·(cos θ cos φ, cos θ sin φ, sin θ): gate, idling
        and readout error at threshold. Run A is shaded by its length p<sub>th</sub> (copper, as in
        arXiv:2306.17786); run B is a flat overlay. The bars on the p<sub>G</sub> and p<sub>T</sub> axes span the
        thresholds of the bias sweep below.
      </p>
      <Plot msg={none && "Threshold vs bias appears with the first threshold run."}
        data={picked.flatMap(([k, p]) => biasTraces(p, surf[p], k))}
        layout={{ ...baseLayout(),
          shapes: picked.flatMap(([, p]) => ([[surf[p].settings.eta_g, G_COLOR], [surf[p].settings.eta_t, T_COLOR]] as const)
            .map(([x, color]) => ({ type: "line", x0: x, x1: x, yref: "paper", y0: 0, y1: 1,
                                    line: { color, width: 1, dash: "dot" } }))),
          xaxis: axis("Error bias η", { type: "log" }),
          yaxis: axis("Threshold (%)", { rangemode: "tozero" }) }} />
      <p className="note">
        Threshold along pure gate noise against η<sub>G</sub> (2-qubit / 1-qubit gate error) and along pure idling
        against η<sub>T</sub> (dephasing / relaxation). Dotted rules mark the η the surface used.
      </p>
      <div className="scroll"><table><tbody>
        <tr><th>Run</th><th>Scheme</th><th>d</th><th>Decoder</th><th>η<sub>G</sub></th><th>η<sub>T</sub></th>
          <th>p<sub>G</sub><sup>th</sup></th><th>p<sub>T</sub><sup>th</sup></th><th>p<sub>R</sub><sup>th</sup></th>
          <th>Directions with a crossing</th></tr>
        {picked.map(([k, p]) => {
          const d = surf[p], s = d.settings;
          return (
            <tr key={k}>
              <td>{k}</td><td>{schemeName(s.scheme)}</td>
              <td className="num">{s.distances.join(", ")}</td><td>{s.decoder}</td>
              <td className="num">{s.eta_g}</td><td className="num">{s.eta_t}</td>
              {d.axes.map((v, i) => <td key={i} className="num">{pct(v)}</td>)}
              <td className="num">{d.surface.filter(r => r.p_th > 0).length} / {d.surface.length}</td>
            </tr>
          );
        })}
      </tbody></table></div>
    </>
  );
}
