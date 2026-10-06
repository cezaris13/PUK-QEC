// Loading and reading the result files under results/: LER CSVs (decode.py, shared/memory.py),
// footprint CSVs (footprint.py) and threshold summaries (thresholds.py's <name>_thresholds.csv).

export const PCT = 100;

/** One row of an LER or footprint CSV. `scheme` is missing in shared/memory.py's files, `qubits`
 *  only in footprint.py's. */
export interface Row {
  scheme?: string;
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

/** One scheme's LER curves from one run, every distance. */
export interface Series {
  key: string;
  label: string;
  scheme: string;
  run: string;
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
  settings: { scheme: string; distances: number[]; decoder: string; eta_g: number; eta_t: number };
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
    settings: { scheme: first.scheme, decoder: first.decoder, eta_g: first.eta_g, eta_t: first.eta_t,
                distances: String(first.distances).split(" ").map(Number) },
    axes: axis.map(r => r.p_th) as [number, number, number],
    surface: group("surface"), g_bias: group("g_bias"), t_bias: group("t_bias"),
  };
}

export interface Results {
  ler: Series[];
  foot: Record<string, Row[]>;
  surf: Record<string, ThresholdRun>;
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
  pairs: "pairs (2 ancillas)", method_a: "method A", method_b: "method B", method_c: "method C",
};
export const schemeName = (s: string) => SCHEME_NAME[s] ?? s;

export async function loadAll(): Promise<Results> {
  const files: string[] = await (await fetch("/api/files")).json();
  const groups: Record<string, Series> = {}, foot: Record<string, Row[]> = {};
  const surf: Record<string, ThresholdRun> = {};
  const texts = await Promise.all(files.map(async f => (await fetch("/" + f)).text()));
  files.forEach((f, i) => {                     // in server order, so the selects stay newest first
    const rows = parseCSV(texts[i]);
    if (rows.length && "group" in rows[0]) { surf[f] = thresholdRun(rows as unknown as SummaryRow[]); return; }
    // thresholds.py's <name>.csv: its raw points, already summed up in <name>_thresholds.csv
    if (!rows.length || !("ler" in rows[0]) || "theta" in rows[0]) return;
    if ("qubits" in rows[0]) { foot[f] = rows; return; }
    // shared/memory.py writes one file per distance and no scheme column
    const run = f.replace(/_d\d+\.csv$/, "").replace(/\.csv$/, "");
    for (const r of rows) {
      const scheme = r.scheme || "pairs";
      const key = run + "|" + scheme + (r.scheme ? "" : "|memory.py");
      (groups[key] ??= { key, run, scheme, rows: [],
        label: `${schemeName(scheme)} · ${base(run)}` + (r.scheme ? "" : " (memory.py)") }).rows.push(r);
    }
  });
  return { ler: Object.values(groups).sort((a, b) => a.label.localeCompare(b.label)), foot, surf };
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
