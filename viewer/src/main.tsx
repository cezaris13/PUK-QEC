import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { App } from "./App";
import { applyTheme, store } from "./theme";
import "./styles.css";

applyTheme(store.get("floquet-theme"));
createRoot(document.getElementById("root")!).render(<StrictMode><App /></StrictMode>);
