// Dev and preview server for the viewer. Besides the app it serves the repo's results/ folder and lists its
// files at /api/files, so `npm run dev` is all there is to run (no separate Python server).
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import react from "@vitejs/plugin-react";
import { defineConfig, type Connect } from "vite";

const REPO = fileURLToPath(new URL("..", import.meta.url));
const RESULTS = path.join(REPO, "results");

/** Every LER / footprint CSV and every thresholds JSON under results/, newest first. */
function listResults(): string[] {
  if (!fs.existsSync(RESULTS)) return [];
  return (fs.readdirSync(RESULTS, { recursive: true }) as string[])
    .map(f => path.join(RESULTS, f))
    .filter(f => f.endsWith(".csv"))
    .map(f => ({ f, t: fs.statSync(f).mtimeMs }))
    .sort((a, b) => b.t - a.t)
    .map(({ f }) => path.relative(REPO, f).split(path.sep).join("/"));
}

const results: Connect.NextHandleFunction = (req, res, next) => {
  const url = decodeURIComponent((req.url ?? "").split("?")[0]);
  if (url === "/api/files") {
    res.setHeader("Content-Type", "application/json");
    return res.end(JSON.stringify(listResults()));
  }
  if (url.startsWith("/results/")) {
    const file = path.join(REPO, url);           // join normalises "..", so check it stayed inside
    if (!file.startsWith(RESULTS + path.sep) || !fs.existsSync(file)) {
      res.statusCode = 404;
      return res.end();
    }
    return fs.createReadStream(file).pipe(res);
  }
  next();
};

export default defineConfig({
  plugins: [react(), {
    name: "results",
    configureServer: s => { s.middlewares.use(results); },
    configurePreviewServer: s => { s.middlewares.use(results); },
  }],
  server: { port: 8000, open: true },
  preview: { port: 8000, open: true },
  build: { chunkSizeWarningLimit: 6000 },   // Plotly alone is ~4.8 MB
});
