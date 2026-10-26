// Loading and reading the result files under results/: LER CSVs (physical_to_logical.py, shared/memory.py),
// footprint CSVs (footprint.py), threshold summaries (thresholds.py's <name>_thresholds.csv) and the X3Z3 vs CSS
// bias sweeps (x3z3_bias.py / x3z3_ler.py).

// The page keeps the results it parsed when it loaded, so hot-swapping this module would pair new code with data
// parsed the old way: reload the page whenever it changes instead.
import.meta.hot?.accept(() => location.reload());

export const PCT = 100;

// Older CSVs record L, the lattice size; newer ones the code distance d = 2L (--distances is d now; a torus of
// 3L x 2L hexagons; floquet.code_distance). App's header switch picks which one labels the plots.
let showCodeDistance = true;
export const setShowCodeDistance = (on: boolean) => { showCodeDistance = on; };
export const codeDistanceShown = () => showCodeDistance;
/** A distance label, the code distance or the lattice size. ponytail: older CSVs record L (always odd: 3, 5, 7),
 *  newer ones the code distance d = 2L (always even), so parity says which a value is. */
export const dLabel = (x: number) => (showCodeDistance ? (x % 2 ? 2 * x : x) : (x % 2 ? x : x / 2));
/** How long a memory experiment runs: `r` periods of the 6-step schedule ("d" for d of them). */
export const roundsLabel = (r: string | number = "d") =>
  !showCodeDistance ? `${r} rounds` : r === "d" ? "d/2 periods" : `${r} periods`;

/** One row of an LER or footprint CSV. `scheme` is missing in shared/memory.py's files, `qubits`
 *  only in footprint.py's. */
export interface Row {
  scheme?: string;
  code?: string;
  basis?: string;
  qubits?: number;
  distance: number;
  rounds: number;
  noise: string;
  eta: number;
  p: number;
  shots: number;
  errors: number;
  ler: number;
}

/** One run's LER rows: every scheme, code, basis and distance in it. */
export interface Series {
  key: string;
  label: string;
  rows: Row[];
}

/** One direction or bias point of thresholds.py, with its crossing p_th. */
export interface ThresholdRow {
  theta: number;
  phi: number;
  eta_g: number;
  eta_t: number;
  p_th: number;
  p_th_error: number;
}

export interface ThresholdRun {
  settings: { scheme: string; code: string; distances: number[]; decoder: string; eta_g: number; eta_t: number };
  axes: [number, number, number];      // p_G, p_T, p_R thresholds along each axis
  surface: ThresholdRow[];
  g_bias: ThresholdRow[];
  t_bias: ThresholdRow[];
}

/** One row of thresholds.py's <name>_thresholds.csv: group is axis, surface, g_bias or t_bias. */
type SummaryRow = ThresholdRow & { group: string; scheme: string; code: string; decoder: string;
                                   distances: string | number };

function thresholdRun(rows: SummaryRow[]): ThresholdRun {
  const axis = rows.filter(r => r.group === "axis"), first = axis[0] ?? rows[0];
  const group = (g: string) => rows.filter(r => r.group === g);
  return {
    settings: { scheme: first.scheme, code: first.code || "css", decoder: first.decoder, eta_g: first.eta_g, eta_t: first.eta_t,
                distances: String(first.distances).split(" ").map(Number) },
    axes: axis.map(r => r.p_th) as [number, number, number],
    surface: group("surface"), g_bias: group("g_bias"), t_bias: group("t_bias"),
  };
}

/** One point of x3z3_ler.py's sweep: one memory basis, errors out of shots, no ler column. */
export interface BiasRow {
  code: string; scheme: string; distance: number; eta: number; p: number; basis: string;
  noise: string; shots: number; errors: number;
}

export interface Results {
  ler: Series[];
  foot: Record<string, Row[]>;
  surf: Record<string, ThresholdRun>;
  bias: Record<string, BiasRow[]>;
}

// ponytail: plain split, fine for the unquoted CSVs this repo writes
export function parseCSV(text: string): Row[] {
  const [head, ...lines] = text.trim().split(/\r?\n/);
  const keys = head.split(",");
  return lines.filter(Boolean).map(l => {
    const v = l.split(","), o: Record<string, string | number> = {};
    keys.forEach((k, i) => o[k] = v[i] === "" || isNaN(+v[i]) ? v[i] : +v[i]);
    return o as unknown as Row;
  });
}

export const base = (path: string) => path.split("/").pop()!;
export const runLabel = (path: string) => base(path).replace(/(_thresholds)?\.csv$/, "");
const SCHEME_NAME: Record<string, string> = {
  pairs: "pairs (2 ancillas)", pairs_bell: "pairs, Bell state", method_a: "method A", method_b: "method B", method_c: "method C",
};
export const schemeName = (s: string) => SCHEME_NAME[s] ?? s;

export async function loadAll(): Promise<Results> {
  const files: string[] = await (await fetch("/api/files")).json();
  const groups: Record<string, Series> = {}, foot: Record<string, Row[]> = {};
  const surf: Record<string, ThresholdRun> = {}, bias: Record<string, BiasRow[]> = {};
  const texts = await Promise.all(files.map(async f => (await fetch("/" + f)).text()));
  files.forEach((f, i) => {                     // in server order, so the selects stay newest first
    const rows = parseCSV(texts[i]);
    if (rows.length && "group" in rows[0]) { surf[f] = thresholdRun(rows as unknown as SummaryRow[]); return; }
    if (rows.length && "basis" in rows[0] && "eta" in rows[0] && !("ler" in rows[0]) && !("theta" in rows[0])) {
      bias[f] = rows as unknown as BiasRow[];
      return;
    }
    // thresholds.py's <name>.csv: its raw points, already summed up in <name>_thresholds.csv
    if (!rows.length || !("ler" in rows[0]) || "theta" in rows[0]) return;
    if ("qubits" in rows[0]) { foot[f] = rows; return; }
    // shared/memory.py writes one file per distance and no scheme column
    const run = f.replace(/_d\d+\.csv$/, "").replace(/\.csv$/, "");
    const memory = rows[0].scheme ? "" : " (memory.py)";
    (groups[run] ??= { key: run, label: base(run) + memory, rows: [] }).rows.push(...rows);
  });
  return { ler: Object.values(groups).sort((a, b) => a.label.localeCompare(b.label)), foot, surf, bias };
}

export const binomSE = (r: Row) => Math.sqrt(r.ler * (1 - r.ler) / r.shots);

export function byDistance(rows: Row[]): Record<number, Row[]> {
  const m: Record<number, Row[]> = {};
  rows.forEach(r => (m[r.distance] ??= []).push(r));
  Object.values(m).forEach(a => a.sort((x, y) => x.p - y.p));
  return m;
}

/** Where the smallest and largest distance cross: last p with the large d better, then log-log
 *  interpolation to the next p. Needs both curves above zero at both points. */
export function crossing(rows: Row[]) {
  const m = byDistance(rows), ds = Object.keys(m).map(Number).sort((a, b) => a - b);
  if (ds.length < 2) return null;
  const lo = m[ds[0]], hi = m[ds[ds.length - 1]];
  const diff = (i: number) => Math.log(hi[i].ler) - Math.log(lo[i].ler);
  for (let i = 0; i + 1 < Math.min(lo.length, hi.length); i++) {
    if (![lo[i], hi[i], lo[i + 1], hi[i + 1]].every(r => r.ler > 0)) continue;
    const a = diff(i), b = diff(i + 1);
    if (a < 0 && b >= 0) {
      const t = a / (a - b);
      return { ds, p: Math.exp(Math.log(lo[i].p) + t * (Math.log(lo[i + 1].p) - Math.log(lo[i].p))),
               lo: lo[i].p, hi: lo[i + 1].p };
    }
  }
  return { ds, p: null, lo: 0, hi: 0 };
}

/** log10 LER = a + b d through the rows with any logical errors. */
export function fitDecay(rows: Row[]) {
  const ok = rows.filter(r => r.ler > 0);
  if (new Set(ok.map(r => r.distance)).size < 2) return null;
  const n = ok.length, xs = ok.map(r => r.distance), ys = ok.map(r => Math.log10(r.ler));
  const mx = xs.reduce((a, b) => a + b) / n, my = ys.reduce((a, b) => a + b) / n;
  const b = xs.reduce((s, x, i) => s + (x - mx) * (ys[i] - my), 0) /
            xs.reduce((s, x) => s + (x - mx) ** 2, 0);
  return { a: my - b * mx, b };
}

/** Smallest odd d whose fitted LER reaches the target. */
export function distanceFor(f: { a: number; b: number } | null, target: number) {
  if (!f || f.b >= 0) return null;
  const d = Math.ceil((Math.log10(target) - f.a) / f.b);
  return d % 2 ? d : d + 1;
}

/** A surface direction's crossing point (p_G, p_T, p_R) in %, with its length p_th and error. */
export const surfPoints = (d: ThresholdRun) => d.surface.filter(r => r.p_th > 0).map(r => ({
  x: PCT * r.p_th * Math.cos(r.theta) * Math.cos(r.phi),
  y: PCT * r.p_th * Math.cos(r.theta) * Math.sin(r.phi),
  z: PCT * r.p_th * Math.sin(r.theta), p: PCT * r.p_th, e: PCT * r.p_th_error,
}));

export const pct = (v: number) => (v * PCT).toFixed(3) + " %";

/** Volume between the origin and the surface, in %³: thresholds.py's grid is a triangle of n(n+1)/2 directions
 *  in (s, t) order, so neighbours form triangles, and each with the origin is a tetrahedron |a·(b×c)| / 6.
 *  `corners` are the surface's own (p_G, p_T, p_R) axis points; `flat` the volume / G T R / 6 of the corner
 *  tetrahedron, 1 when the surface is the flat triangle through them. Triangles with a direction that found no
 *  crossing are left out and counted in `missing`. */
export function surfVolume(d: ThresholdRun) {
  const rows = d.surface, N = rows.length, n = Math.round((Math.sqrt(8 * N + 1) - 1) / 2);
  if (!N || (n * (n + 1)) / 2 !== N) return null;
  const at = (i: number, j: number) => i * n - (i * (i - 1)) / 2 + j;  // row of grid point (i, j), i + j < n
  const vec = (k: number) => {
    const r = rows[k], p = PCT * r.p_th;
    return r.p_th > 0 ? [p * Math.cos(r.theta) * Math.cos(r.phi), p * Math.cos(r.theta) * Math.sin(r.phi),
                         p * Math.sin(r.theta)] : null;
  };
  const tet = (a: number[], b: number[], c: number[]) => Math.abs(
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
  let volume = 0, missing = 0, total = 0;
  for (let i = 0; i < n; i++) for (let j = 0; i + j < n; j++) {
    const tris = [];
    if (i + j + 1 < n) tris.push([at(i, j), at(i + 1, j), at(i, j + 1)]);
    if (i + j + 2 < n) tris.push([at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)]);
    for (const t of tris) {
      total++;
      const [a, b, c] = t.map(vec);
      if (a && b && c) volume += tet(a, b, c); else missing++;
    }
  }
  const corners = [at(0, 0), at(n - 1, 0), at(0, n - 1)].map(k => PCT * rows[k].p_th);  // G, T, R
  const box = corners[0] * corners[1] * corners[2] / 6;
  return { volume, corners, flat: box > 0 ? volume / box : null, missing, total };
}
