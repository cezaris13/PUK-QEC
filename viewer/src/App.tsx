import { useEffect, useState } from "react";
import { BiasSection } from "./components/BiasSection";
import { FootSection } from "./components/FootSection";
import { LerSection } from "./components/LerSection";
import { SurfaceSection } from "./components/SurfaceSection";
import { codeDistanceShown, loadAll, setShowCodeDistance, type Results } from "./data";
import { applyTheme, isDark, store } from "./theme";

const THEME_KEY = "floquet-theme";
const DIST_KEY = "floquet-distance";  // "lattice" to label with L; anything else, the code distance 2L
setShowCodeDistance(store.get(DIST_KEY) !== "lattice");

export function App() {
  const [data, setData] = useState<Results | null>(null);
  const [error, setError] = useState("");
  const [, setDark] = useState(isDark());     // a re-render reads the CSS tokens afresh
  const [codeD, setCodeD] = useState(codeDistanceShown());
  const pickDistance = (on: boolean) => {
    setShowCodeDistance(on);                  // before the render, so every label reads it
    store.set(DIST_KEY, on ? "code" : "lattice");
    setCodeD(on);
  };

  useEffect(() => { loadAll().then(setData, e => setError(String(e))); }, []);
  useEffect(() => {
    // follow the OS only while no explicit choice is stored
    const mq = matchMedia("(prefers-color-scheme: dark)");
    const on = () => { applyTheme(store.get(THEME_KEY)); setDark(isDark()); };
    mq.addEventListener("change", on);
    return () => mq.removeEventListener("change", on);
  }, []);

  const flip = () => {
    const next = isDark() ? "light" : "dark";
    store.set(THEME_KEY, next);
    applyTheme(next);                         // before the render, so css() sees the new tokens
    setDark(next === "dark");
  };

  return (
    <>
      <div className="hdr">
        <div>
          <h1>Floquet code results</h1>
          <p className="sub">
            Spin-ancilla honeycomb Floquet code: two ancillas per edge (pairs) against one ancilla per edge
            (methods A, B, C). Read straight from the CSV/JSON files under <code>results/</code>.
          </p>
        </div>
        <label className="theme" title="CSVs record L (--distances); the code distance is 2L">
          <input type="checkbox" checked={codeD} onChange={e => pickDistance(e.target.checked)} />{" "}
          {codeD ? "code distance d = 2L" : "lattice size L"}
        </label>
        <button type="button" className="theme" onClick={flip}>{isDark() ? "☀ Light" : "☾ Dark"}</button>
      </div>
      {error ? <p className="empty">Could not load results: {error}</p>
        : !data ? <p className="empty">Loading results…</p>
        : <>
            <LerSection ler={data.ler} />
            <SurfaceSection surf={data.surf} />
            <BiasSection bias={data.bias} />
            <FootSection foot={data.foot} />
          </>}
    </>
  );
}
