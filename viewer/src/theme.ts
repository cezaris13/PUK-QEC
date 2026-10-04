// Light/dark theme and the Plotly styling that follows it. Colours come from the CSS tokens in styles.css,
// read at render time, so a re-render after a theme switch picks them up.

const ROOT = document.documentElement;

export const css = (name: string) => getComputedStyle(ROOT).getPropertyValue(name).trim();

export const isDark = () => ROOT.dataset.theme
  ? ROOT.dataset.theme === "dark"
  : matchMedia("(prefers-color-scheme: dark)").matches;

/** localStorage that never throws (private windows, blocked storage). */
export const store = {
  get: (k: string) => { try { return localStorage.getItem(k); } catch { return null; } },
  set: (k: string, v: string) => { try { localStorage.setItem(k, v); } catch { /* per-viewer nicety only */ } },
};

/** "light", "dark", or null to follow the OS. */
export function applyTheme(pref: string | null) {
  if (pref) ROOT.dataset.theme = pref; else delete ROOT.dataset.theme;
}

export const baseLayout = () => ({
  paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
  font: { color: css("--text-secondary"), size: 12 },
  margin: { l: 64, r: 20, t: 30, b: 50 },
  legend: { bgcolor: "rgba(0,0,0,0)", font: { color: css("--text-primary") } },
});

export const axis = (title: string, extra: object = {}) => ({
  title: { text: title }, gridcolor: css("--border"), zerolinecolor: css("--border"),
  color: css("--text-secondary"), ...extra,
});

// Ordered ramp over code distance; dark mode gets its own lighter steps rather than a flip.
const DIST_COLORS = ["#7fc97f", "#1a9850", "#74add1", "#2166ac", "#f4a582", "#d73027"];
const DIST_COLORS_DARK = ["#a3d9a3", "#4fb05a", "#8ec6e8", "#5e9fd6", "#f5b79c", "#e46b63"];
export const distColor = (d: number) => (isDark() ? DIST_COLORS_DARK : DIST_COLORS)[((d - 3) % 6 + 6) % 6];

export const DASH = ["solid", "dash", "dot", "dashdot", "longdash", "longdashdot"];
export const MARK = ["circle", "square", "diamond", "triangle-up", "x", "star"];
export const SCHEME_COLORS = ["#2a78d6", "#eb6834", "#1a9850", "#8e44ad", "#d73027"];
export const G_COLOR = "#2a78d6", T_COLOR = "#1a9850";
export const COPPER = [[0, "rgb(0,0,0)"], [1, "rgb(255,199,127)"]];  // matplotlib "copper", arXiv:2306.17786
