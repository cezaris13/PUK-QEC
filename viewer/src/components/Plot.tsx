import Plotly from "plotly.js-dist-min";
import { useEffect, useRef, type ReactNode } from "react";

// ponytail: traces and layouts are loose objects; Plotly's typings miss mesh3d fields (intensity,
// delaunayaxis, lighting), so they are cast once here rather than fought at every call site
export type Trace = Record<string, unknown>;

const CFG = { responsive: true, displaylogo: false };

/** A Plotly chart redrawn on every render, with an optional message laid over it (empty states). */
export function Plot({ data, layout, tall, msg }:
  { data: Trace[]; layout: Record<string, unknown>; tall?: boolean; msg?: ReactNode }) {
  const ref = useRef<HTMLDivElement>(null);
  useEffect(() => {
    Plotly.react(ref.current!, data as Plotly.Data[], layout as Partial<Plotly.Layout>, CFG);
  });
  useEffect(() => {
    const el = ref.current!;
    return () => Plotly.purge(el);
  }, []);
  return (
    <div className="holder">
      <div ref={ref} className={tall ? "plot" : "plot small"} />
      {msg && <div className="msg">{msg}</div>}
    </div>
  );
}
