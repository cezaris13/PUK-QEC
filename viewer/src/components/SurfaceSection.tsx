import { useState } from "react";
import { dLabel, pct, PCT, runLabel, schemeName, surfPoints, surfVolume, type ThresholdRun } from "../data";
import { axis, baseLayout, COPPER, css, G_COLOR, T_COLOR } from "../theme";
import { Plot, type Trace } from "./Plot";

type Key = "A" | "B";

// thresholds.py writes .../thresholds/<code>/...; older summaries have no code column
const codeOf = (path: string) => (path.includes("/x3z3/") || /(^|\/)x3z3_/.test(path) ? "x3z3" : "css");

/** The eighth of the sphere of radius `r` in the octant, as points: `alphahull: 0` hulls them into the solid,
 *  with the origin closing its three flat faces on the coordinate planes. */
function octantPoints(r: number, n = 16): [number, number, number][] {
  const pts: [number, number, number][] = [[0, 0, 0]];
  for (let i = 0; i <= n; i++) for (let j = 0; j <= n; j++) {
    const th = (Math.PI / 2 * i) / n, ph = (Math.PI / 2 * j) / n;
    pts.push([r * Math.cos(th) * Math.cos(ph), r * Math.cos(th) * Math.sin(ph), r * Math.sin(th)]);
  }
  return pts;
}

/** Each run's surface and dots shaded by p_th (A copper, B black to its series colour, each with its colourbar),
 * plus the bias-sweep ranges as bars on the p_G and p_T axes. */
function surfaceTraces(path: string, d: ThresholdRun, key: Key, spheres: boolean): Trace[] {
  const pts = surfPoints(d), isA = key === "A";
  const scale = isA ? COPPER : [[0, "rgb(0,0,0)"], [1, css("--series-b")]];
  const P = (k: "x" | "y" | "z" | "p") => pts.map(q => q[k]);
  const range = { cmin: Math.min(...P("p")), cmax: Math.max(...P("p")) };
  const out: Trace[] = [];
  if (pts.length >= 3) out.push({
    type: "mesh3d", x: P("x"), y: P("y"), z: P("z"), delaunayaxis: "z", name: runLabel(path),
    hoverinfo: "skip", flatshading: true, lighting: { ambient: 0.9, diffuse: 0.2, specular: 0 },
    intensity: P("p"), colorscale: scale, ...range, showscale: true, opacity: isA ? 0.5 : 0.35,
    colorbar: { title: { text: `p_th (%) · ${key}`, side: "top" }, thickness: 14, len: 0.7, x: isA ? 1.0 : 1.12,
                outlinewidth: 0, tickfont: { color: css("--text-secondary") } },
  });
  out.push({
    type: "scatter3d", mode: "markers", x: P("x"), y: P("y"), z: P("z"), showlegend: false,
    customdata: pts.map(q => [q.p, q.e]),
    // coloured by p_th on the surface's scale, as thresholds.py's scatter
    marker: { size: 4, color: P("p"), colorscale: scale, ...range, symbol: isA ? "circle" : "diamond" },
    hovertemplate: "p_G %{x:.3f} %<br>p_T %{y:.3f} %<br>p_R %{z:.3f} %<br>" +
                   `p_th %{customdata[0]:.3f} ± %{customdata[1]:.3f} %<extra>${runLabel(path)}</extra>`,
  });
  // the isotropic octant the surface is contained in: the sphere through its largest p_th
  const v = spheres ? surfVolume(d) : null;
  if (v?.radius) {
    const ball = octantPoints(v.radius);
    out.push({
      type: "mesh3d", x: ball.map(q => q[0]), y: ball.map(q => q[1]), z: ball.map(q => q[2]),
      alphahull: 0, color: css(isA ? "--series-a" : "--series-b"), opacity: 0.15, flatshading: true,
      lighting: { ambient: 0.9, diffuse: 0.2, specular: 0 }, hoverinfo: "skip",
      name: `sphere · ${key} (r ${v.radius.toFixed(3)} %, ${v.sphere.toFixed(4)} %³)`,
    });
  }
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
  const [runA, setA] = useState(paths[0] ?? "");
  const [runB, setB] = useState("");
  // By scheme: A is the scheme's CSS run, B its X3Z3 run (the newest of each, as paths are newest first)
  const [byScheme, setByScheme] = useState(false);
  // the isotropic sphere octants, to see the measured volume against them
  const [spheres, setSpheres] = useState(false);
  // equal scale: one % the same length on every axis, so the sphere octant is round; off stretches the box
  // to a cube, which fills the frame but distorts shapes when the axes' ranges differ
  const [equalScale, setEqualScale] = useState(true);
  const schemes = [...new Set(paths.map(p => surf[p].settings.scheme))];
  const [scheme, setScheme] = useState(schemes[0] ?? "");
  const pickedScheme = schemes.includes(scheme) ? scheme : schemes[0];
  const runOf = (code: string) => paths.find(p => surf[p].settings.scheme === pickedScheme && codeOf(p) === code) ?? "";
  const [a, b] = byScheme ? [runOf("css"), runOf("x3z3")] : [runA, runB];
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
          <span className="lbl">Pick</span>
          <label>
            <input type="checkbox" checked={byScheme} onChange={e => setByScheme(e.target.checked)} />{" "}
            by scheme: CSS as A, X³Z³ as B
          </label>
        </div>
        <div className="field">
          <span className="lbl">Show</span>
          <label>
            <input type="checkbox" checked={spheres} onChange={e => setSpheres(e.target.checked)} />{" "}
            sphere octant through the largest p<sub>th</sub>
          </label>
          <label>
            <input type="checkbox" checked={equalScale} onChange={e => setEqualScale(e.target.checked)} />{" "}
            equal scale on all three axes
          </label>
        </div>
        {byScheme ? (
          <div className="field">
            <label htmlFor="selScheme">Scheme</label>
            <select id="selScheme" value={pickedScheme} onChange={e => setScheme(e.target.value)}>
              {schemes.map(s => <option key={s} value={s}>{schemeName(s)}</option>)}
            </select>
          </div>
        ) : (<>
          <div className="field">
            <label htmlFor="selA"><span className="swatch" style={{ background: "var(--series-a)" }} />Run A</label>
            <select id="selA" value={runA} onChange={e => setA(e.target.value)}>{options}</select>
          </div>
          <div className="field">
            <label htmlFor="selB"><span className="swatch" style={{ background: "var(--series-b)" }} />Run B (compare)</label>
            <select id="selB" value={runB} onChange={e => setB(e.target.value)}>
              <option value="">— none —</option>{options}
            </select>
          </div>
        </>)}
      </div>
      <div className="row"><div>
      <Plot tall msg={none}
        data={picked.length ? picked.flatMap(([k, p]) => surfaceTraces(p, surf[p], k, spheres)) : EMPTY_3D}
        layout={{ ...baseLayout(), margin: { l: 0, r: 0, t: 0, b: 0 },
          legend: { ...baseLayout().legend, x: 0, y: 1 },
          scene: { xaxis: axis("p_G gate (%)", { rangemode: "tozero" }),
                   yaxis: axis("p_T idling (%)", { rangemode: "tozero" }),
                   zaxis: axis("p_R readout (%)", { rangemode: "tozero" }),
                   bgcolor: "rgba(0,0,0,0)", aspectmode: equalScale ? "data" : "cube",
                   camera: { eye: { x: 1.6, y: 1.6, z: 1.0 } } } }} />
      <p className="note">
        Each dot is one direction's crossing point p<sub>th</sub>·(cos θ cos φ, cos θ sin φ, sin θ): gate, idling
        and readout error at threshold. Each run, surface and dots, is shaded by its length p<sub>th</sub>:
        run A copper (as in arXiv:2306.17786 and thresholds.py), run B black to its series colour (diamonds). The bars on the p<sub>G</sub> and p<sub>T</sub> axes span the
        thresholds of the bias sweep on the right.
      </p>
      <p className="note">
        Volume under surface: the region between the origin and the surface, in %³, summed as one tetrahedron per
        grid triangle with the origin. Larger means more combinations of gate, idling and readout error stay below
        threshold. Tick "sphere octant" above to draw the isotropic volume it sits inside: the octant of the
        sphere through the largest p<sub>th</sub> (π r³ / 6), what the volume would be if every direction
        tolerated as much as the best one. Measured / flat compares it instead to the
        tetrahedron of the surface's own three axis points (p<sub>G</sub> p<sub>T</sub> p<sub>R</sub> / 6):
        above 1 the surface bulges outward, so mixed noise is tolerated better than the axis thresholds alone
        suggest.
      </p>
      </div><div>
      <Plot tall msg={none && "Threshold vs bias appears with the first threshold run."}
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
      </div></div>
      <div className="scroll"><table><tbody>
        <tr><th>Run</th><th>Scheme</th><th>d</th><th>Decoder</th><th>η<sub>G</sub></th><th>η<sub>T</sub></th>
          <th>p<sub>G</sub><sup>th</sup></th><th>p<sub>T</sub><sup>th</sup></th><th>p<sub>R</sub><sup>th</sup></th>
          <th>Directions with a crossing</th>
          <th>Volume under surface (%³)</th><th>Sphere octant (%³)</th>
          <th>Measured / flat</th></tr>
        {picked.map(([k, p]) => {
          const d = surf[p], s = d.settings, v = surfVolume(d);
          return (
            <tr key={k}>
              <td>{k}</td><td>{schemeName(s.scheme)}{codeOf(p) === "x3z3" ? " · X³Z³" : " · CSS"}</td>
              <td className="num">{s.distances.map(dLabel).join(", ")}</td><td>{s.decoder}</td>
              <td className="num">{s.eta_g}</td><td className="num">{s.eta_t}</td>
              {d.axes.map((v, i) => <td key={i} className="num">{pct(v)}</td>)}
              <td className="num">{d.surface.filter(r => r.p_th > 0).length} / {d.surface.length}</td>
              <td className="num">{v ? v.volume.toFixed(4) + (v.missing ? ` (${v.missing}/${v.total} triangles missing)` : "") : "—"}</td>
              <td className="num">{v?.sphere ? v.sphere.toFixed(4) : "—"}</td>
              <td className="num">{v?.flat ? `× ${v.flat.toFixed(2)}` : "—"}</td>
            </tr>
          );
        })}
      </tbody></table></div>
    </>
  );
}
