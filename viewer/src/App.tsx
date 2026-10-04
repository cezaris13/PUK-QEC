import { useEffect, useState } from "react";
import { FootSection } from "./components/FootSection";
import { LerSection } from "./components/LerSection";
import { SurfaceSection } from "./components/SurfaceSection";
import { loadAll, type Results } from "./data";
import { applyTheme, isDark, store } from "./theme";

const THEME_KEY = "floquet-theme";

export function App() {
  const [data, setData] = useState<Results | null>(null);
  const [error, setError] = useState("");
  const [, setDark] = useState(isDark());     // a re-render reads the CSS tokens afresh

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
        <button type="button" className="theme" onClick={flip}>{isDark() ? "☀ Light" : "☾ Dark"}</button>
      </div>
      {error ? <p className="empty">Could not load results: {error}</p>
        : !data ? <p className="empty">Loading results…</p>
        : <>
            <LerSection ler={data.ler} />
            <SurfaceSection surf={data.surf} />
            <FootSection foot={data.foot} />
          </>}
    </>
  );
}
