import { useState, useMemo, useEffect, useRef, useCallback } from "react";

// ===== KaTeX =====
function useKaTeX() { const [r, setR] = useState(false); useEffect(() => { if (window.katex) { setR(true); return; } const l = document.createElement("link"); l.rel = "stylesheet"; l.href = "https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/katex.min.css"; document.head.appendChild(l); const s = document.createElement("script"); s.src = "https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/katex.min.js"; s.onload = () => setR(true); document.head.appendChild(s); }, []); return r; }
function Tx({ m, d = false }) { const ref = useRef(null), rdy = useKaTeX(); useEffect(() => { if (rdy && ref.current && window.katex) try { window.katex.render(m, ref.current, { displayMode: d, throwOnError: false, trust: true }); } catch (e) { if (ref.current) ref.current.textContent = m; } }, [m, d, rdy]); return <span ref={ref} />; }
function TB({ m }) { return <div style={{ margin: "10px 0", overflowX: "auto" }}><Tx m={m} d={true} /></div>; }
function Def({ show, children }) { if (!show) return null; return <div style={{ background: "#eef2ff", border: "1px solid #a5b4fc", borderRadius: 6, padding: "8px 12px", margin: "6px 0", fontSize: 12, color: "#312e81", lineHeight: 1.7, borderLeft: "4px solid #6366f1" }}>{children}</div>; }

// ===== Matrix =====
function tr(M) { const r = M.length, c = M[0].length, T = []; for (let j = 0; j < c; j++) { T[j] = []; for (let i = 0; i < r; i++)T[j][i] = M[i][j]; } return T; }
function mm(A, B) { const rA = A.length, cA = A[0].length, cB = B[0].length, C = Array.from({ length: rA }, () => Array(cB).fill(0)); for (let i = 0; i < rA; i++)for (let j = 0; j < cB; j++)for (let k = 0; k < cA; k++)C[i][j] += A[i][k] * B[k][j]; return C; }
function inv(M) { const n = M.length, a = M.map((r, i) => [...r, ...Array(n).fill(0).map((_, j) => i === j ? 1 : 0)]); for (let i = 0; i < n; i++) { let mx = i; for (let k = i + 1; k < n; k++)if (Math.abs(a[k][i]) > Math.abs(a[mx][i])) mx = k;[a[i], a[mx]] = [a[mx], a[i]]; const p = a[i][i]; if (Math.abs(p) < 1e-12) return null; for (let j = 0; j < 2 * n; j++)a[i][j] /= p; for (let k = 0; k < n; k++) { if (k === i) continue; const f = a[k][i]; for (let j = 0; j < 2 * n; j++)a[k][j] -= f * a[i][j]; } } return a.map(r => r.slice(n)); }

// ===== Critical value tables =====
const FT = { "2_9": { 0.01: 8.02, 0.05: 4.26, 0.10: 3.01 }, "4_15": { 0.01: 4.89, 0.05: 3.06, 0.10: 2.36 } };
const TT = { 9: { 0.005: 3.250, 0.01: 1.833, 0.025: 2.262, 0.05: 1.833 }, 15: { 0.005: 2.947, 0.01: 2.602, 0.025: 2.131, 0.05: 1.753 } };
function getFc(a, d1, d2) { const k = `${d1}_${d2}`; return FT[k]?.[a] || 4.26; }
function getTc(a, df) { return TT[df]?.[a / 2] || 2.262; }

const ff = (v, d = 4) => Number(v).toFixed(d), f2 = v => Number(v).toFixed(2);
const ALPHAS = [0.01, 0.05, 0.10];

// ===== F-distribution PDF =====
function lnGamma(z) { const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]; if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lnGamma(1 - z); z -= 1; let x = c[0]; for (let i = 1; i < 9; i++) x += c[i] / (z + i); const t = z + 7.5; return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x); }
function fPDF(x, d1, d2) { if (x <= 0) return 0; const lnB = lnGamma(d1 / 2) + lnGamma(d2 / 2) - lnGamma((d1 + d2) / 2); return Math.exp((d1 / 2) * Math.log(d1 / d2) + (d1 / 2 - 1) * Math.log(x) - ((d1 + d2) / 2) * Math.log(1 + d1 * x / d2) - lnB); }

// ===== Default data =====
const DEF1 = [[3.30, 3.42, 3.36, 3.34], [3.25, 3.15, 3.30, 3.20], [3.10, 3.25, 3.18, 3.12]];
const DEF2 = [[9.8, 3.3, 2.8, 3.1, 4.1], [12.6, 4.4, 4.9, 3.5, 3.9], [11.9, 3.9, 5.3, 4.8, 4.7], [13.1, 5.9, 2.6, 3.1, 3.6], [13.3, 4.6, 5.1, 5.0, 4.1], [13.5, 5.2, 3.2, 3.3, 4.3], [10.1, 4.0, 4.0, 3.3, 4.0], [13.1, 4.7, 4.5, 3.5, 3.8], [10.7, 4.5, 4.1, 3.7, 3.6], [11.0, 3.7, 3.6, 3.3, 3.6], [13.0, 4.6, 4.6, 3.6, 3.6], [11.6, 4.7, 3.5, 3.5, 3.7], [12.0, 3.9, 4.6, 3.6, 4.1], [11.4, 4.6, 4.0, 3.4, 3.6], [12.2, 5.1, 3.6, 3.3, 4.0], [12.8, 5.0, 4.4, 3.6, 3.7], [12.4, 4.8, 4.4, 3.4, 3.6], [13.2, 5.3, 3.5, 3.6, 3.7], [10.6, 3.9, 3.8, 3.4, 4.0], [7.9, 3.4, 3.8, 3.4, 3.4]];
const DEFXP = [5.1, 4.7, 4.8, 4.0];

// ===== Calculation engines =====
function calc1(data, sigma, alpha) {
  const nm = ["SwiftVen", "SwiftFast", "SwiftPay"], k = 3, n = 4, N = 12;
  const mn = data.map(g => g.reduce((a, b) => a + b, 0) / n);
  const gm = data.flat().reduce((a, b) => a + b, 0) / N;
  const tau = mn.map(x => x - gm);
  let SST = 0; mn.forEach(x => { SST += n * (x - gm) ** 2; });
  let SSE = 0; data.forEach((g, i) => g.forEach(v => { SSE += (v - mn[i]) ** 2; }));
  const SSTot = SST + SSE, MST = SST / (k - 1), MSE = SSE / (N - k), F0 = MST / MSE;
  const Fc = getFc(alpha, 2, 9);
  const tc = getTc(alpha, 9);
  const LSD = tc * Math.sqrt(2 * MSE / n);
  const cmp = []; for (let i = 0; i < k; i++)for (let j = i + 1; j < k; j++) { const df = Math.abs(mn[i] - mn[j]); cmp.push({ a: nm[i], b: nm[j], i: i + 1, j: j + 1, diff: df, sig: df > LSD }); }
  const vr = data.map((g, i) => g.reduce((a, v) => a + (v - mn[i]) ** 2, 0) / (n - 1));
  const res = []; data.forEach((g, i) => g.forEach(v => res.push(v - mn[i])));
  const resSorted = [...res].sort((a, b) => a - b);
  // Runs test on residuals in observation order
  let runs = 1; for (let i = 1; i < res.length; i++) if ((res[i] >= 0) !== (res[i - 1] >= 0)) runs++;
  const n1 = res.filter(e => e >= 0).length, n2 = res.filter(e => e < 0).length;
  const ER = (2 * n1 * n2) / (n1 + n2) + 1;
  const VR = (2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)) / ((n1 + n2) ** 2 * ((n1 + n2) - 1));
  const sdR = Math.sqrt(VR);
  const Z0runs = sdR > 0 ? (runs - ER) / sdR : 0;
  const Zcrit = alpha === 0.01 ? 2.576 : alpha === 0.05 ? 1.96 : 1.645;
  const runsOk = Math.abs(Z0runs) <= Zcrit;
  return { data, nm, k, n, N, sigma, alpha, mn, gm, tau, SST, SSE, SSTot, MST, MSE, F0, Fc, tc, LSD, cmp, vr, res, resSorted, runs, n1, n2, ER, VR, sdR, Z0runs, Zcrit, runsOk, rejected: F0 > Fc };
}

function calc2(raw, alpha, xp) {
  const no = raw.length, p = 4;
  const Y = raw.map(r => [r[0]]), X = raw.map(r => [1, r[1], r[2], r[3], r[4]]);
  const Xt = tr(X), XtX = mm(Xt, X), Xi = inv(XtX);
  if (!Xi) return null;
  const XtY = mm(Xt, Y), B = mm(Xi, XtY), b = B.map(x => x[0]);
  const Yh = mm(X, B), yv = Y.map(r => r[0]), ym = yv.reduce((a, x) => a + x, 0) / no;
  let SSR = 0, SSE = 0, SST = 0;
  for (let i = 0; i < no; i++) { SSR += (Yh[i][0] - ym) ** 2; SSE += (yv[i] - Yh[i][0]) ** 2; SST += (yv[i] - ym) ** 2; }
  const MSR = SSR / p, MSE = SSE / (no - p - 1), Fr = MSR / MSE, R2 = SSR / SST, R2a = 1 - (1 - R2) * (no - 1) / (no - p - 1);
  const se = b.map((_, i) => Math.sqrt(Math.abs(MSE * Xi[i][i])));
  const tv = b.map((x, i) => se[i] > 0 ? x / se[i] : 0);
  const df = no - p - 1, tc = getTc(alpha, df);
  const Fc = getFc(alpha, 4, df);
  const cl = b.map((x, i) => x - tc * se[i]), ch = b.map((x, i) => x + tc * se[i]);
  const ypr = [1, ...xp].reduce((a, x, i) => a + x * b[i], 0);
  const sigVars = [1, 2, 3, 4].filter(i => Math.abs(tv[i]) > tc);
  const nsVars = [1, 2, 3, 4].filter(i => Math.abs(tv[i]) <= tc);
  return { raw, no, p, b, se, tv, SSR, SSE, SST, MSR, MSE, Fr, R2, R2a, ym, df, tc, Fc, cl, ch, ypr, sigVars, nsVars, rejectedGlobal: Fr > Fc, alpha };
}

// ===== Number input =====
function NI({ value, onChange, w = 58 }) {
  return <input type="number" step="any" value={value} onChange={e => onChange(parseFloat(e.target.value) || 0)}
    style={{ width: w, padding: "3px 4px", border: "1px solid #cbd5e1", borderRadius: 4, fontSize: 12, textAlign: "center", background: "#fffbeb" }} />;
}

// ===== Fisher Chart =====
function FisherChart({ d1, d2, F0, Fc, alpha }) {
  const canvasRef = useRef(null);
  useEffect(() => {
    const canvas = canvasRef.current; if (!canvas) return;
    const ctx = canvas.getContext("2d");
    const W = canvas.width, H = canvas.height;
    const pad = { l: 40, r: 15, t: 15, b: 32 }, pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
    const xMax = Math.max(F0 * 1.4, Fc * 2, 8), steps = 300, dx = xMax / steps;
    let pts = [], yMax = 0;
    for (let i = 0; i <= steps; i++) { const x = i * dx, y = fPDF(x, d1, d2); if (isFinite(y)) { pts.push([x, y]); if (y > yMax) yMax = y; } }
    yMax *= 1.15;
    const toX = x => pad.l + (x / xMax) * pw, toY = y => pad.t + ph - (y / yMax) * ph;
    ctx.clearRect(0, 0, W, H);
    // rejection region fill
    ctx.beginPath(); ctx.moveTo(toX(Fc), toY(0));
    for (const [x, y] of pts) if (x >= Fc) ctx.lineTo(toX(x), toY(y));
    ctx.lineTo(toX(xMax), toY(0)); ctx.closePath();
    ctx.fillStyle = "rgba(233,69,96,0.22)"; ctx.fill();
    // curve
    ctx.beginPath(); pts.forEach(([x, y], i) => i === 0 ? ctx.moveTo(toX(x), toY(y)) : ctx.lineTo(toX(x), toY(y)));
    ctx.strokeStyle = "#0f3460"; ctx.lineWidth = 2.2; ctx.stroke();
    // axes
    ctx.strokeStyle = "#555"; ctx.lineWidth = 1; ctx.beginPath();
    ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, pad.t + ph); ctx.lineTo(pad.l + pw, pad.t + ph); ctx.stroke();
    // Fc line
    ctx.setLineDash([5, 4]); ctx.strokeStyle = "#e94560"; ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(toX(Fc), toY(0)); ctx.lineTo(toX(Fc), toY(fPDF(Fc, d1, d2))); ctx.stroke();
    // F0 line
    ctx.strokeStyle = "#2e7d32";
    ctx.beginPath(); ctx.moveTo(toX(Math.min(F0, xMax * 0.95)), toY(0)); ctx.lineTo(toX(Math.min(F0, xMax * 0.95)), pad.t + 8); ctx.stroke();
    ctx.setLineDash([]);
    // labels
    ctx.font = "11px 'Segoe UI',sans-serif"; ctx.textAlign = "center";
    ctx.fillStyle = "#e94560"; ctx.fillText(`Fc=${Fc}`, toX(Fc), toY(0) + 14);
    ctx.fillStyle = "#2e7d32"; ctx.fillText(`F\u2080=${Number(F0).toFixed(2)}`, toX(Math.min(F0, xMax * 0.95)), pad.t + 2);
    ctx.fillStyle = "#e94560"; ctx.font = "10px 'Segoe UI',sans-serif";
    ctx.fillText(`\u03B1=${alpha}`, Math.min(toX(Fc) + 35, W - 30), toY(fPDF(Fc, d1, d2) * 0.5));
    // title
    ctx.font = "bold 11px 'Segoe UI',sans-serif"; ctx.fillStyle = "#0f3460"; ctx.textAlign = "left";
    ctx.fillText(`Distribución F(${d1}, ${d2})`, pad.l + 5, pad.t + 12);
    ctx.fillStyle = "#333"; ctx.textAlign = "center"; ctx.font = "10px 'Segoe UI',sans-serif";
    ctx.fillText("F", W / 2, H - 3);
  }, [d1, d2, F0, Fc, alpha]);
  return <canvas ref={canvasRef} width={500} height={200} style={{ width: "100%", maxWidth: 500, height: "auto", border: "1px solid #e0e0e0", borderRadius: 6, background: "#fafbfc", marginTop: 8 }} />;
}

// ===== MAIN =====
export default function App() {
  const [tab, setTab] = useState("p1");
  const [s1, setS1] = useState(0);
  const [s2, setS2] = useState(0);
  const [defs, setDefs] = useState(false);
  const [edit, setEdit] = useState(false);
  const kr = useKaTeX();

  // Part 1 state
  const [d1, setD1] = useState(DEF1.map(r => [...r]));
  const [sig, setSig] = useState(0.18);
  const [a1, setA1] = useState(0.05);

  // Part 2 state
  const [d2, setD2] = useState(DEF2.map(r => [...r]));
  const [a2, setA2] = useState(0.05);
  const [xp, setXp] = useState([...DEFXP]);

  const p1 = useMemo(() => calc1(d1, sig, a1), [d1, sig, a1]);
  const p2 = useMemo(() => calc2(d2, a2, xp), [d2, a2, xp]);

  const D = defs;
  const updD1 = (i, j, v) => { const n = d1.map(r => [...r]); n[i][j] = v; setD1(n); };
  const updD2 = (i, j, v) => { const n = d2.map(r => [...r]); n[i][j] = v; setD2(n); };
  const updXp = (i, v) => { const n = [...xp]; n[i] = v; setXp(n); };
  const reset1 = () => { setD1(DEF1.map(r => [...r])); setSig(0.18); setA1(0.05); };
  const reset2 = () => { setD2(DEF2.map(r => [...r])); setA2(0.05); setXp([...DEFXP]); };

  const t1 = ["1. Modelo", "2. ANOVA", "3. Comparaciones", "4. Normalidad", "5. Aleatoriedad", "6. Varianzas", "7. Supuestos"];
  const t2 = ["1. MCO", "2. Predicción", "3. ANOVA", "4. R\u00B2", "5. IC", "6. Signif.", "7. Modelo"];
  const bN = ["\\beta_0", "\\beta_1", "\\beta_2", "\\beta_3", "\\beta_4"];
  const bD = ["Intercepto", "Peticiones/s", "Tam. trama", "Lat. bóveda", "Mem. microserv."];
  const arqN = ["SwiftVen", "SwiftFast", "SwiftPay"];

  const S = {
    bx: { fontFamily: "'Segoe UI',system-ui,sans-serif", maxWidth: 960, margin: "0 auto", padding: 12, color: "#1a1a2e", fontSize: 14, paddingBottom: 65 },
    hd: { textAlign: "center", padding: "16px 0 8px", borderBottom: "3px solid #0f3460" },
    h1: { fontSize: 19, fontWeight: 700, color: "#0f3460", margin: 0 },
    h2: { fontSize: 12, color: "#555", margin: "4px 0 0", fontWeight: 400 },
    tbar: { display: "flex", gap: 4, marginTop: 12, flexWrap: "wrap" },
    mt: a => ({ padding: "7px 16px", border: "none", borderRadius: "6px 6px 0 0", cursor: "pointer", fontWeight: 600, fontSize: 13, background: a ? "#0f3460" : "#ddd", color: a ? "#fff" : "#333" }),
    sb: { display: "flex", gap: 2, flexWrap: "wrap", marginBottom: 8, background: "#f0f0f0", padding: 4, borderRadius: 6 },
    st: a => ({ padding: "4px 8px", border: "none", borderRadius: 4, cursor: "pointer", fontSize: 11, fontWeight: a ? 600 : 400, background: a ? "#16213e" : "transparent", color: a ? "#fff" : "#333", whiteSpace: "nowrap" }),
    cd: { background: "#fff", border: "1px solid #ddd", borderRadius: 8, padding: 14, marginBottom: 10, boxShadow: "0 1px 3px rgba(0,0,0,0.06)" },
    T: { width: "100%", borderCollapse: "collapse", fontSize: 12, marginTop: 6 },
    th: { background: "#0f3460", color: "#fff", padding: "5px 6px", textAlign: "center", fontSize: 11 },
    thF: { background: "#1e3a5f", color: "#fff", padding: "5px 6px", textAlign: "center", fontSize: 11, fontStyle: "italic" },
    td: { padding: "4px 6px", textAlign: "center", borderBottom: "1px solid #eee", fontSize: 12 },
    fm: { background: "#f8f9fa", padding: "8px 12px", borderRadius: 6, margin: "8px 0", borderLeft: "3px solid #0f3460", overflowX: "auto" },
    sp: { background: "#fafbfc", padding: "10px 12px", borderRadius: 6, marginBottom: 7, borderLeft: "4px solid #e94560" },
    spT: { fontWeight: 700, color: "#e94560", fontSize: 12, marginBottom: 3 },
    ok: { background: "#e8f5e9", padding: "10px 12px", borderRadius: 6, fontWeight: 600, color: "#2e7d32", marginTop: 8, fontSize: 13 },
    wn: { background: "#fff3e0", padding: "10px 12px", borderRadius: 6, color: "#e65100", marginTop: 8 },
    q: { background: "#e3f2fd", padding: "8px 12px", borderRadius: 6, fontWeight: 600, color: "#0d47a1", marginBottom: 8, borderLeft: "4px solid #1976d2", fontSize: 13 },
    ed: { background: "#fffbeb", border: "1px solid #f59e0b", borderRadius: 8, padding: 12, marginBottom: 10 },
    edH: { display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 8 },
    btn: (bg) => ({ padding: "5px 12px", border: "none", borderRadius: 5, cursor: "pointer", fontSize: 11, fontWeight: 600, background: bg, color: "#fff" }),
    sel: { padding: "4px 6px", border: "1px solid #cbd5e1", borderRadius: 4, fontSize: 12, background: "#fffbeb" },
    fab: { position: "fixed", bottom: 18, right: 18, display: "flex", gap: 6, zIndex: 999 },
    fabBtn: (bg) => ({ padding: "7px 12px", borderRadius: 18, border: "none", cursor: "pointer", fontSize: 11, fontWeight: 600, boxShadow: "0 2px 8px rgba(0,0,0,0.2)", background: bg, color: "#fff" }),
  };

  if (!kr) return <div style={{ textAlign: "center", padding: 40 }}>Cargando LaTeX...</div>;

  // ===== EDITORS =====
  const Ed1 = () => (<div style={S.ed}>
    <div style={S.edH}><b style={{ color: "#92400e", fontSize: 13 }}>Editar Parámetros — Parte I</b><button style={S.btn("#dc2626")} onClick={reset1}>Restablecer</button></div>
    <div style={{ display: "flex", gap: 16, flexWrap: "wrap", alignItems: "center", marginBottom: 8 }}>
      <label style={{ fontSize: 12 }}><Tx m="\sigma" /> = <NI value={sig} onChange={setSig} w={55} /></label>
      <label style={{ fontSize: 12 }}><Tx m="\alpha" /> = <select style={S.sel} value={a1} onChange={e => setA1(parseFloat(e.target.value))}>{ALPHAS.map(a => <option key={a} value={a}>{a}</option>)}</select></label>
    </div>
    <table style={S.T}><thead><tr><th style={S.th}>Arquitectura</th>{[1, 2, 3, 4].map(j => <th key={j} style={S.th}>Rep {j}</th>)}</tr></thead>
      <tbody>{d1.map((row, i) => <tr key={i}><td style={{ ...S.td, fontWeight: 600 }}>{arqN[i]}</td>{row.map((v, j) => <td key={j} style={S.td}><NI value={v} onChange={x => updD1(i, j, x)} w={52} /></td>)}</tr>)}</tbody></table>
  </div>);

  const Ed2 = () => (<div style={S.ed}>
    <div style={S.edH}><b style={{ color: "#92400e", fontSize: 13 }}>Editar Parámetros — Parte II</b><button style={S.btn("#dc2626")} onClick={reset2}>Restablecer</button></div>
    <div style={{ display: "flex", gap: 16, flexWrap: "wrap", alignItems: "center", marginBottom: 8 }}>
      <label style={{ fontSize: 12 }}><Tx m="\alpha" /> = <select style={S.sel} value={a2} onChange={e => setA2(parseFloat(e.target.value))}>{ALPHAS.map(a => <option key={a} value={a}>{a}</option>)}</select></label>
      <span style={{ fontSize: 12 }}>Predicción:</span>
      {xp.map((v, i) => <label key={i} style={{ fontSize: 12 }}><Tx m={`x_${i + 1}`} />=<NI value={v} onChange={x => updXp(i, x)} w={48} /></label>)}
    </div>
    <div style={{ maxHeight: 220, overflowY: "auto", border: "1px solid #e5e7eb", borderRadius: 4 }}>
      <table style={{ ...S.T, marginTop: 0 }}><thead style={{ position: "sticky", top: 0 }}><tr><th style={S.th}>#</th><th style={S.th}>y</th><th style={S.th}><Tx m="x_1" /></th><th style={S.th}><Tx m="x_2" /></th><th style={S.th}><Tx m="x_3" /></th><th style={S.th}><Tx m="x_4" /></th></tr></thead>
        <tbody>{d2.map((row, i) => <tr key={i}><td style={{ ...S.td, fontWeight: 600, fontSize: 10 }}>{i + 1}</td>{row.map((v, j) => <td key={j} style={S.td}><NI value={v} onChange={x => updD2(i, j, x)} w={46} /></td>)}</tr>)}</tbody></table>
    </div>
  </div>);

  // ===== PART 1 PANELS =====
  const P1 = [
    () => (<div>
      <div style={S.q}>Pregunta 1: Plantee el modelo estadístico y estime <Tx m="\mu" />, <Tx m="\sigma^2" /> y <Tx m="\tau_i" />.</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Modelo Estadístico</h3>
        <div style={S.fm}><TB m={`Y_{ij} = \\mu + \\tau_i + \\varepsilon_{ij}, \\qquad \\varepsilon_{ij} \\sim N(0,\\,\\sigma^2)`} /></div>
        <Def show={D}><b><Tx m="Y_{ij}" />:</b> Observación j del tratamiento i.<br /><b><Tx m="\mu" />:</b> Media global: <Tx m="\mu = \tfrac{1}{k}\sum \mu_i" />.<br /><b><Tx m="\tau_i" />:</b> Efecto del tratamiento: <Tx m="\tau_i = \mu_i - \mu" />, con <Tx m="\sum \tau_i = 0" />.<br /><b><Tx m="\varepsilon_{ij}" />:</b> Error iid <Tx m="\sim N(0,\sigma^2)" />, independiente del tratamiento.<br /><b>Modelo:</b> DCA — 1 factor, k=3 niveles, n=4 réplicas/nivel.</Def>
        <p><Tx m={`i = 1,2,3 \\;\\text{(arquitecturas)},\\quad j = 1,2,3,4 \\;\\text{(réplicas)}`} /></p>
      </div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Tabla de fórmulas</h3>
        <table style={S.T}><thead><tr><th style={S.thF}>Símbolo</th><th style={S.thF}>Fórmula</th><th style={S.thF}>Descripción</th></tr></thead>
          <tbody>
            <tr><td style={S.td}><Tx m="\hat{\mu}" /></td><td style={S.td}><Tx m="\bar{y}_{..} = \tfrac{1}{N}\sum_i\sum_j y_{ij}" /></td><td style={{ ...S.td, textAlign: "left" }}>Media global muestral</td></tr>
            <tr style={{ background: "#f9f9f9" }}><td style={S.td}><Tx m="\bar{y}_i" /></td><td style={S.td}><Tx m="\tfrac{1}{n}\sum_{j=1}^{n}y_{ij}" /></td><td style={{ ...S.td, textAlign: "left" }}>Media del grupo i</td></tr>
            <tr><td style={S.td}><Tx m="\hat{\tau}_i" /></td><td style={S.td}><Tx m="\bar{y}_i - \bar{y}_{..}" /></td><td style={{ ...S.td, textAlign: "left" }}>Efecto estimado del tratamiento</td></tr>
            <tr style={{ background: "#f9f9f9" }}><td style={S.td}><Tx m="\sigma^2" /></td><td style={S.td}><Tx m={`(${ff(sig, 2)})^2`} /></td><td style={{ ...S.td, textAlign: "left" }}>Varianza asumida por diseño</td></tr>
          </tbody></table>
      </div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Valores calculados</h3>
        <table style={S.T}><thead><tr><th style={S.th}>Arquitectura</th>{[1, 2, 3, 4].map(j => <th key={j} style={S.th}>{j}</th>)}<th style={S.th}><Tx m="\bar{y}_i" /></th><th style={S.th}><Tx m="\hat{\tau}_i" /></th></tr></thead>
          <tbody>{p1.nm.map((nm, i) => <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}><td style={{ ...S.td, fontWeight: 600 }}>{nm} (<Tx m={`A_{${i + 1}}`} />)</td>{p1.data[i].map((v, j) => <td key={j} style={S.td}>{ff(v, 2)}</td>)}<td style={{ ...S.td, fontWeight: 700, color: "#0f3460" }}>{ff(p1.mn[i], 4)}</td><td style={{ ...S.td, fontWeight: 600, color: p1.tau[i] > 0 ? "#c62828" : "#2e7d32" }}>{p1.tau[i] > 0 ? "+" : ""}{ff(p1.tau[i], 4)}</td></tr>)}</tbody></table>
        <div style={S.fm}><TB m={`\\hat{\\mu} = ${ff(p1.gm, 4)}\\text{ seg},\\quad \\sigma^2 = ${ff(sig * sig, 4)}\\text{ seg}^2,\\quad \\sigma = ${ff(sig, 2)}\\text{ seg}`} /></div>
        <div style={S.ok}><b>Respuesta:</b> La latencia promedio global del core bancario es <Tx m={`\\hat{\\mu}=${ff(p1.gm, 4)}`} /> seg. SwiftPay presenta el mayor efecto de reducción (<Tx m={`\\hat{\\tau}_3=${ff(p1.tau[2], 4)}`} />), mientras SwiftVen es la de mayor latencia (<Tx m={`\\hat{\\tau}_1=${p1.tau[0] > 0 ? "+" : ""}{${ff(p1.tau[0], 4)}}`} />).</div>
      </div>
    </div>),

    // ANOVA 8 pasos
    () => (<div>
      <div style={S.q}>Pregunta 2: Con ANOVA (8 pasos), ¿hay al menos una arquitectura con desempeño diferente?</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>ANOVA — 8 Pasos</h3>
        <div style={S.sp}><div style={S.spT}>Paso 1: Parámetros</div>
          <p><Tx m={`k=3,\\;n=4,\\;N=12,\\;\\sigma=${ff(sig, 2)}`} /> seg. Parámetros: <Tx m="\mu_1,\\mu_2,\\mu_3" /> (latencias medias de cada arquitectura).</p>
          <Def show={D}><b><Tx m="\mu_i" />:</b> Media poblacional del nivel i: <Tx m="\mu_i = \mu + \tau_i" />. Parámetro desconocido a estimar.</Def>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 2: Hipótesis nula</div>
          <div style={S.fm}><TB m="H_0: \\mu_1 = \\mu_2 = \\mu_3" /></div><p>Las tres arquitecturas tienen igual latencia media.</p>
          <Def show={D}><b><Tx m="H_0" />:</b> Equivale a <Tx m="\tau_1=\tau_2=\cdots=\tau_k=0" />. Se rechaza solo si <Tx m={`F_0 > F_{\\alpha,k-1,N-k}`} />.</Def>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 3: Hipótesis alternativa</div>
          <div style={S.fm}><TB m="H_1: \\text{Al menos un }\\mu_i\\text{ difiere}" /></div>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 4: Nivel de significancia</div>
          <div style={S.fm}><TB m={`\\alpha = ${a1}`} /></div>
          <Def show={D}><b><Tx m="\alpha" />:</b> <Tx m={"P(\\text{rechazar }H_0 \\mid H_0\\text{ verdadera}) \\leq \\alpha"} />. Error Tipo I: concluir diferencia inexistente.</Def>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 5: Estadístico</div>
          <div style={S.fm}><TB m={`F_0 = \\frac{MST}{MSE} \\sim F_{(k-1,\\,N-k)} = F_{(2,\\,9)}`} /></div>
          <Def show={D}><b><Tx m="F_0" />:</b> Cociente de varianzas: <Tx m="MST/MSE" />. Bajo <Tx m="H_0" />, <Tx m={"F_0 \\sim F_{(k-1,\\,N-k)}"} />. Valores grandes → diferencia entre tratamientos.</Def>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 6: Criterio de rechazo</div>
          <div style={S.fm}><TB m={`\\text{Rechazar }H_0\\text{ si }F_0 > F_{${a1},\\,2,\\,9} = ${p1.Fc}`} /></div>
          <Def show={D}><b><Tx m={`F_{\\text{crit}}`} />:</b> Percentil <Tx m="1-\\alpha" /> de <Tx m="F_{(k-1,\\,N-k)}" />. Región de rechazo: cola derecha.</Def>
          <div style={{ marginTop: 8 }}><FisherChart d1={2} d2={9} F0={p1.F0} Fc={p1.Fc} alpha={a1} /></div>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 7: Cálculos</div>
          <h4 style={{ margin: "6px 0 3px", color: "#333", fontSize: 13 }}>Tabla de fórmulas</h4>
          <table style={S.T}><thead><tr><th style={S.thF}>Fuente</th><th style={S.thF}>SS</th><th style={S.thF}>gl</th><th style={S.thF}>MS</th><th style={S.thF}><Tx m="F" /></th></tr></thead>
            <tbody>
              <tr><td style={{ ...S.td, fontWeight: 600 }}>Tratamientos</td><td style={S.td}><Tx m={"n\\sum(\\bar{y}_i-\\bar{y}_{..})^2"} /></td><td style={S.td}><Tx m="k{-}1" /></td><td style={S.td}><Tx m="\\tfrac{SST}{k-1}" /></td><td style={S.td}><Tx m="\\tfrac{MST}{MSE}" /></td></tr>
              <tr style={{ background: "#f9f9f9" }}><td style={{ ...S.td, fontWeight: 600 }}>Error</td><td style={S.td}><Tx m={"\\sum\\!\\sum(y_{ij}-\\bar{y}_i)^2"} /></td><td style={S.td}><Tx m="N{-}k" /></td><td style={S.td}><Tx m="\\tfrac{SSE}{N-k}" /></td><td style={S.td}>—</td></tr>
              <tr><td style={{ ...S.td, fontWeight: 600 }}>Total</td><td style={S.td}><Tx m="SST{+}SSE" /></td><td style={S.td}><Tx m="N{-}1" /></td><td style={S.td}>—</td><td style={S.td}>—</td></tr>
            </tbody></table>
          <Def show={D}><b>SST</b> = <Tx m={"n\\sum(\\bar{y}_i - \\bar{y}_{..})^2"} />: variabilidad entre tratamientos.<br /><b>SSE</b> = <Tx m={"\\sum\\sum(y_{ij} - \\bar{y}_i)^2"} />: variabilidad residual.<br /><b>MS</b> = SS/gl: normalización por grados de libertad.</Def>
          <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>Valores numéricos</h4>
          <table style={S.T}><thead><tr><th style={S.th}>Fuente</th><th style={S.th}>SS</th><th style={S.th}>gl</th><th style={S.th}>MS</th><th style={S.th}><Tx m="F_0" /></th><th style={S.th}><Tx m="F_{\text{crit}}" /></th></tr></thead>
            <tbody>
              <tr><td style={{ ...S.td, fontWeight: 600 }}>Tratamientos</td><td style={S.td}>{ff(p1.SST, 6)}</td><td style={S.td}>2</td><td style={S.td}>{ff(p1.MST, 6)}</td><td style={{ ...S.td, fontWeight: 700, color: "#c62828" }}>{f2(p1.F0)}</td><td style={S.td}>{p1.Fc}</td></tr>
              <tr style={{ background: "#f9f9f9" }}><td style={{ ...S.td, fontWeight: 600 }}>Error</td><td style={S.td}>{ff(p1.SSE, 6)}</td><td style={S.td}>9</td><td style={S.td}>{ff(p1.MSE, 6)}</td><td style={S.td}>—</td><td style={S.td}>—</td></tr>
              <tr style={{ fontWeight: 600 }}><td style={{ ...S.td, fontWeight: 700 }}>Total</td><td style={S.td}>{ff(p1.SSTot, 6)}</td><td style={S.td}>11</td><td style={S.td}>—</td><td style={S.td}>—</td><td style={S.td}>—</td></tr>
            </tbody></table>
          <div style={{ ...S.fm, marginTop: 10 }}><TB m={`F_0 = ${f2(p1.F0)} \\;${p1.rejected ? ">" : "\\leq"}\\; ${p1.Fc} \\;\\Longrightarrow\\; \\boxed{\\text{${p1.rejected ? "Se rechaza" : "No se rechaza"} } H_0}`} /></div>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 8: Conclusión</div>
          <div style={S.ok}>{p1.rejected
            ? <>Con un error del {a1 * 100}%, medias <Tx m={`\\bar{y}_1=${ff(p1.mn[0], 4)}`} />, <Tx m={`\\bar{y}_2=${ff(p1.mn[1], 4)}`} />, <Tx m={`\\bar{y}_3=${ff(p1.mn[2], 4)}`} /> y <Tx m={`\\sigma=${ff(sig, 2)}`} /> seg, se concluye que <b>al menos una arquitectura del core bancario presenta desempeño significativamente diferente</b> en latencia (<Tx m={`F_0=${f2(p1.F0)}>${p1.Fc}`} />).</>
            : <>Con un error del {a1 * 100}%, <b>no hay evidencia suficiente</b> para afirmar que alguna arquitectura difiera en latencia (<Tx m={`F_0=${f2(p1.F0)}\\leq${p1.Fc}`} />).</>
          }</div>
        </div>
      </div>
    </div>),

    // Comparaciones
    () => (<div>
      <div style={S.q}>Pregunta 3: ¿Recomendaría SwiftPay (<Tx m="A_3" />) como estándar para el banco?</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Comparaciones Múltiples — LSD de Fisher</h3>
        <Def show={D}><b>LSD:</b> Prueba post-hoc que compara todos los pares <Tx m="(i,j)" /> tras rechazar <Tx m="H_0" />. Significativa si <Tx m="|\\bar{y}_i - \\bar{y}_j| > LSD" />.</Def>
        <h4 style={{ margin: "6px 0 3px", color: "#333", fontSize: 13 }}>Fórmula y cálculo</h4>
        <div style={S.fm}><TB m={`LSD = t_{\\alpha/2,\\,N-k}\\cdot\\sqrt{\\frac{2\\cdot MSE}{n}} = ${p1.tc}\\times\\sqrt{\\frac{2\\times${ff(p1.MSE, 6)}}{4}} = ${ff(p1.LSD, 4)}`} /></div>
        <Def show={D}><b><Tx m={`t_{${a1 / 2},\\,N-k}`} />:</b> Percentil de Student con <Tx m="N-k" /> gl. <Tx m={"LSD = t_{\\alpha/2}\\sqrt{2 \\cdot MSE / n}"} />: mínima diferencia significativa entre pares.</Def>
        <h4 style={{ margin: "8px 0 3px", color: "#333", fontSize: 13 }}>Resultados</h4>
        <table style={S.T}><thead><tr><th style={S.th}>Comparación</th><th style={S.th}><Tx m={"|\\bar{y}_i-\\bar{y}_j|"} /></th><th style={S.th}>LSD</th><th style={S.th}>Resultado</th></tr></thead>
          <tbody>{p1.cmp.map((x, i) => <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}><td style={{ ...S.td, fontWeight: 600 }}>{x.a} vs {x.b}</td><td style={S.td}>{ff(x.diff, 4)}</td><td style={S.td}>{ff(p1.LSD, 4)}</td><td style={{ ...S.td, fontWeight: 700, color: x.sig ? "#2e7d32" : "#c62828" }}>{x.sig ? "Significativa" : "No significativa"}</td></tr>)}</tbody></table>
        <div style={S.ok}>{p1.cmp[2]?.sig && !p1.cmp[1]?.sig
          ? <><b>Respuesta:</b> SwiftPay (<Tx m="A_3" />) tiene la menor latencia ({ff(p1.mn[2], 4)} seg) y difiere de SwiftVen, pero <b>NO difiere significativamente de SwiftFast</b>. No se puede recomendar como estándar único; se sugiere evaluar costo, estabilidad y escalabilidad entre SwiftFast y SwiftPay.</>
          : p1.cmp.every(x => x.sig)
            ? <><b>Respuesta:</b> Todas las comparaciones son significativas. SwiftPay (<Tx m="A_3" />) con media {ff(p1.mn[2], 4)} seg <b>sí se recomienda</b> como estándar al ser significativamente mejor que ambas alternativas.</>
            : <><b>Respuesta:</b> Con los datos actuales, las diferencias significativas son: {p1.cmp.filter(x => x.sig).map(x => `${x.a} vs ${x.b}`).join(", ") || "ninguna"}. Se requiere evaluar los resultados en contexto para decidir.</>
        }</div>
      </div>
    </div>),

    // Normalidad
    () => (<div>
      <div style={S.q}>Pregunta 4: ¿Apoyaría el supuesto de normalidad?</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Supuesto de Normalidad</h3>
        <Def show={D}><b>Supuesto:</b> <Tx m="\varepsilon_{ij} \sim N(0,\sigma^2)" />. Violación → <Tx m="F_0" /> y p-valor no confiables. Se verifica mediante residuos <Tx m={"e_{ij} = y_{ij} - \\bar{y}_i"} />.</Def>
        <div style={S.fm}><TB m={`e_{ij} = y_{ij} - \\bar{y}_i`} /></div>
        <Def show={D}><b><Tx m="e_{ij}" />:</b> Residuo = estimación de <Tx m="\varepsilon_{ij}" />. Si el modelo es correcto, <Tx m={"e_{ij} \\approx N(0,\\sigma^2)"} /> sin patrones.</Def>
        <p style={{ fontSize: 12 }}>Residuos ordenados:</p>
        <div style={S.fm}><TB m={`\\{${p1.resSorted.map(r => ff(r, 4)).join(',\\;')}\\}`} /></div>
        <p style={{ fontSize: 13 }}>• Simétricos alrededor de 0 • Sin outliers • Rango acotado</p>
        <div style={S.ok}><b>Respuesta:</b> Con <Tx m="n=4" /> por grupo, las pruebas formales tienen bajo poder. Los residuos no muestran asimetrías ni valores extremos, por lo que <b>el supuesto de normalidad se sostiene razonablemente</b> para las mediciones de latencia del core bancario.</div>
      </div>
    </div>),

    // Aleatoriedad — Estadístico Z de rachas
    () => (<div>
      <div style={S.q}>Pregunta 5: ¿Los datos fueron obtenidos al azar?</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Supuesto de Aleatoriedad — Estadístico Z</h3>
        <Def show={D}>
          <b>Racha:</b> Subsecuencia consecutiva de residuos del mismo signo (positivo o negativo) en el orden de observación.<br />
          <b>Residuo:</b> <Tx m={"e_{ij} = y_{ij} - \\bar{y}_i"} /> — desviación de cada observación respecto a la media de su tratamiento.<br />
          <b><Tx m="n_1" />:</b> Número de residuos <Tx m={"\\geq 0"} />. <b><Tx m="n_2" />:</b> Número de residuos <Tx m={"< 0"} />. <b><Tx m="N = n_1 + n_2" />.</b><br />
          <b><Tx m="R" />:</b> Total de rachas observadas.<br />
          <b><Tx m="E(R)" />:</b> Valor esperado: <Tx m={"E(R) = \\frac{2n_1 n_2}{N} + 1"} />.<br />
          <b><Tx m={"\\sigma_R"} />:</b> Desviación estándar: <Tx m={"\\sigma_R = \\sqrt{\\frac{2n_1 n_2(2n_1 n_2 - N)}{N^2(N-1)}}"} />.<br />
          <b><Tx m="Z_0" />:</b> Estadístico normalizado: <Tx m={"Z_0 = \\frac{R - E(R)}{\\sigma_R}"} />. Mide cuántas desviaciones estándar se aleja <Tx m="R" /> de su valor esperado.<br />
          <b>Criterio:</b> Si <Tx m={`|Z_0| \\leq Z_{\\alpha/2}`} />, los residuos son consistentes con independencia.
        </Def>

        <div style={S.cd}>
          <h4 style={{ margin: "0 0 6px", color: "#0f3460", fontSize: 13 }}>1. Residuos en orden de observación</h4>
          <table style={S.T}><thead><tr><th style={S.th}>Obs.</th><th style={S.th}><Tx m="y_{ij}" /></th><th style={S.th}><Tx m={"\\bar{y}_i"} /></th><th style={S.th}><Tx m="e_{ij}" /></th></tr></thead>
            <tbody>{(() => { let idx = 0; return p1.data.flatMap((g, i) => g.map((v, j) => { const e = v - p1.mn[i]; idx++; return <tr key={`${i}-${j}`} style={{ background: idx % 2 ? "#f9f9f9" : "#fff" }}><td style={{ ...S.td, fontSize: 10 }}>{p1.nm[i]}-{j + 1}</td><td style={S.td}>{ff(v, 2)}</td><td style={S.td}>{ff(p1.mn[i], 4)}</td><td style={{ ...S.td, fontWeight: 600, color: e >= 0 ? "#2e7d32" : "#c62828" }}>{ff(e, 4)}</td></tr>; })); })()}</tbody></table>
        </div>

        <div style={S.cd}>
          <h4 style={{ margin: "0 0 6px", color: "#0f3460", fontSize: 13 }}>2. Parámetros de la prueba</h4>
          <table style={S.T}><thead><tr><th style={S.th}>Parámetro</th><th style={S.th}>Fórmula</th><th style={S.th}>Valor</th></tr></thead>
            <tbody>
              <tr><td style={{ ...S.td, fontWeight: 600 }}><Tx m="n_1" /> (residuos <Tx m={"\\geq 0"} />)</td><td style={S.td}>—</td><td style={{ ...S.td, fontWeight: 700 }}>{p1.n1}</td></tr>
              <tr style={{ background: "#f9f9f9" }}><td style={{ ...S.td, fontWeight: 600 }}><Tx m="n_2" /> (residuos <Tx m={"< 0"} />)</td><td style={S.td}>—</td><td style={{ ...S.td, fontWeight: 700 }}>{p1.n2}</td></tr>
              <tr><td style={{ ...S.td, fontWeight: 600 }}><Tx m="N" /></td><td style={S.td}><Tx m="n_1 + n_2" /></td><td style={{ ...S.td, fontWeight: 700 }}>{p1.N}</td></tr>
              <tr style={{ background: "#f9f9f9" }}><td style={{ ...S.td, fontWeight: 600, color: "#0f3460" }}><Tx m="R" /> (rachas)</td><td style={S.td}>—</td><td style={{ ...S.td, fontWeight: 700, color: "#0f3460" }}>{p1.runs}</td></tr>
            </tbody></table>
        </div>

        <div style={S.cd}>
          <h4 style={{ margin: "0 0 6px", color: "#0f3460", fontSize: 13 }}>3. Cálculo del estadístico Z</h4>
          <div style={S.fm}><TB m={`E(R) = \\frac{2(${p1.n1})(${p1.n2})}{${p1.N}} + 1 = ${ff(p1.ER, 4)}`} /></div>
          <div style={S.fm}><TB m={`\\sigma_R = \\sqrt{\\frac{2(${p1.n1})(${p1.n2})\\big[2(${p1.n1})(${p1.n2}) - ${p1.N}\\big]}{${p1.N}^2(${p1.N}-1)}} = ${ff(p1.sdR, 4)}`} /></div>
          <div style={{ ...S.fm, borderLeft: "3px solid #e94560" }}><TB m={`Z_0 = \\frac{R - E(R)}{\\sigma_R} = \\frac{${p1.runs} - ${ff(p1.ER, 4)}}{${ff(p1.sdR, 4)}} = ${ff(p1.Z0runs, 4)}`} /></div>
        </div>

        <div style={S.cd}>
          <h4 style={{ margin: "0 0 6px", color: "#0f3460", fontSize: 13 }}>4. Criterio de rechazo</h4>
          <div style={S.fm}><TB m={`\\text{Independencia si } |Z_0| \\leq Z_{\\alpha/2}`} /></div>
          <div style={S.fm}><TB m={`Z_{\\alpha/2} = Z_{${ff(a1 / 2, 3)}} = ${p1.Zcrit}`} /></div>
          <div style={{ ...S.fm, borderLeft: `3px solid ${p1.runsOk ? "#2e7d32" : "#c62828"}` }}><TB m={`|Z_0| = ${ff(Math.abs(p1.Z0runs), 4)} \\;${p1.runsOk ? "\\leq" : ">"}\\; ${p1.Zcrit} \\;\\Longrightarrow\\; \\boxed{${p1.runsOk ? "\\text{Independencia}" : "\\text{No independencia}"}}`} /></div>
        </div>

        <div style={p1.runsOk ? S.ok : S.wn}><b>Respuesta:</b> Con <Tx m={`R=${p1.runs}`} /> rachas en <Tx m={`N=${p1.N}`} /> observaciones (<Tx m={`n_1=${p1.n1}, n_2=${p1.n2}`} />), se obtiene <Tx m={`Z_0=${ff(p1.Z0runs, 4)}`} />. Dado que <Tx m={`|Z_0|=${ff(Math.abs(p1.Z0runs), 4)}`} /> {p1.runsOk ? "≤" : ">"} <Tx m={`${p1.Zcrit}`} /> (<Tx m={`\\alpha=${a1}`} />), {p1.runsOk ? <><b>los residuos son independientes</b> y el supuesto de aleatoriedad se sostiene.</> : <><b>hay evidencia de dependencia</b> en los residuos; el supuesto de aleatoriedad es cuestionable.</>}</div>
      </div>
    </div>),

    // Varianzas
    () => (<div>
      <div style={S.q}>Pregunta 6: ¿Las tres arquitecturas presentan la misma varianza?</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Homocedasticidad</h3>
        <Def show={D}><b>Homocedasticidad:</b> <Tx m="\sigma_1^2 = \sigma_2^2 = \cdots = \sigma_k^2" />. Violación inflata o deflata <Tx m="F_0" />, invalidando el nivel <Tx m="\alpha" />.</Def>
        <div style={S.fm}><TB m={`s_i^2 = \\frac{1}{n-1}\\sum_{j=1}^{n}(y_{ij}-\\bar{y}_i)^2`} /></div>
        <table style={S.T}><thead><tr><th style={S.th}>Arquitectura</th><th style={S.th}><Tx m="s_i^2" /></th><th style={S.th}><Tx m="s_i" /></th></tr></thead>
          <tbody>{p1.nm.map((nm, i) => <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}><td style={{ ...S.td, fontWeight: 600 }}>{nm}</td><td style={S.td}>{ff(p1.vr[i], 6)}</td><td style={S.td}>{ff(Math.sqrt(p1.vr[i]), 4)}</td></tr>)}</tbody></table>
        {(() => {
          const rat = Math.max(...p1.vr) / Math.min(...p1.vr); return <>
            <div style={S.fm}><TB m={`\\frac{s^2_{\\max}}{s^2_{\\min}} = ${ff(rat, 2)}`} /></div>
            <Def show={D}><b>Regla:</b> Si <Tx m={"s^2_{\\max}/s^2_{\\min} < 3"} /> con <Tx m="n_i" /> iguales, ANOVA es robusto ante heterocedasticidad moderada.</Def>
            <div style={S.ok}><b>Respuesta:</b> La razón de varianzas es {ff(rat, 2)}. {rat < 3 ? "Esto indica una violación leve; " : "Esto sugiere que las varianzas no son exactamente iguales; sin embargo, "}con tamaños iguales (<Tx m="n=4" />), <b>ANOVA es robusto</b> ante esta violación y las conclusiones se mantienen.</div>
          </>;
        })()}
      </div>
    </div>),

    // Supuestos
    () => (<div>
      <div style={S.q}>Pregunta 7: Si no se cumplen los supuestos, ¿se mantienen las conclusiones del inciso 2?</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Robustez de las Conclusiones</h3>
        <Def show={D}><b>Robustez:</b> Con <Tx m="n_i" /> iguales, ANOVA tolera violaciones moderadas de normalidad y homocedasticidad sin alterar significativamente el nivel <Tx m="\alpha" />.</Def>
        <table style={S.T}><thead><tr><th style={S.th}>Supuesto</th><th style={S.th}>Estado</th><th style={S.th}>Impacto</th></tr></thead>
          <tbody>
            <tr><td style={{ ...S.td, fontWeight: 600 }}>Normalidad</td><td style={{ ...S.td, color: "#2e7d32", fontWeight: 600 }}>Razonable</td><td style={{ ...S.td, textAlign: "left" }}>Sin evidencia en contra</td></tr>
            <tr style={{ background: "#f9f9f9" }}><td style={{ ...S.td, fontWeight: 600 }}>Aleatoriedad</td><td style={{ ...S.td, color: "#2e7d32", fontWeight: 600 }}>Razonable</td><td style={{ ...S.td, textAlign: "left" }}>Depende del diseño experimental</td></tr>
            <tr><td style={{ ...S.td, fontWeight: 600 }}>Homocedasticidad</td><td style={{ ...S.td, color: "#e65100", fontWeight: 600 }}>Cuestionable</td><td style={{ ...S.td, textAlign: "left" }}>Razón = {ff(Math.max(...p1.vr) / Math.min(...p1.vr), 2)}</td></tr>
          </tbody></table>
        <p style={{ fontSize: 13, marginTop: 8 }}>1. <b>Diseño balanceado</b> (<Tx m="n_1=n_2=n_3=4" />) → ANOVA robusto.</p>
        <p style={{ fontSize: 13 }}>2. <b><Tx m="F_0" /> holgado:</b> <Tx m={`${f2(p1.F0)}/${p1.Fc} \\approx ${ff(p1.F0 / p1.Fc, 1)}`} /> veces el valor crítico.</p>
        <p style={{ fontSize: 13 }}>3. <b>Normalidad</b> no rechazada por inspección de residuos.</p>
        <Def show={D}><b>Señal fuerte:</b> <Tx m={"F_0/F_c \\gg 1"} /> indica que violaciones moderadas de supuestos no alteran la decisión estadística.</Def>
        <div style={S.ok}><b>Respuesta:</b> Aun si los supuestos no se cumplen estrictamente, <b>las conclusiones se mantienen</b>. El diseño balanceado y <Tx m={`F_0=${f2(p1.F0)}`} /> (ampliamente superior a {p1.Fc}) otorgan robustez. Se recomienda más réplicas para confirmar con mayor rigor.</div>
      </div>
    </div>),
  ];

  // ===== PART 2 PANELS =====
  const P2 = p2 ? [
    () => (<div>
      <div style={S.q}>Pregunta 1: Estimadores del modelo RLM por MCO.</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Mínimos Cuadrados Ordinarios</h3>
        <div style={S.fm}><TB m={"\\hat{\\boldsymbol{\\beta}} = (\\mathbf{X}^\\top\\mathbf{X})^{-1}\\mathbf{X}^\\top\\mathbf{Y}"} /></div>
        <Def show={D}><b>MCO:</b> Minimiza <Tx m={"\\sum e_i^2 = \\sum(y_i - \\hat{y}_i)^2"} />. Solución: <Tx m={"\\hat{\\beta}=(X'X)^{-1}X'Y"} />.<br /><b>SE:</b> <Tx m={"SE(\\hat{\\beta}_i) = \\sqrt{MSE \\cdot c_{ii}}"} />, donde <Tx m={"c_{ii} = [(X'X)^{-1}]_{ii}"} />.</Def>
        <h4 style={{ margin: "8px 0 3px", color: "#333", fontSize: 13 }}>Tabla de símbolos</h4>
        <table style={S.T}><thead><tr><th style={S.thF}>Param.</th><th style={S.thF}>Variable</th><th style={S.thF}>Significado</th></tr></thead>
          <tbody>{bN.map((nm, i) => <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}><td style={{ ...S.td, fontWeight: 600 }}><Tx m={nm} /></td><td style={S.td}>{i === 0 ? "—" : <Tx m={`x_${i}`} />}</td><td style={{ ...S.td, textAlign: "left" }}>{bD[i]}</td></tr>)}</tbody></table>
        <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>Valores estimados</h4>
        <table style={S.T}><thead><tr><th style={S.th}>Param.</th><th style={S.th}><Tx m={"\\hat{\\beta}_i"} /></th><th style={S.th}><Tx m={"SE(\\hat{\\beta}_i)"} /></th></tr></thead>
          <tbody>{bN.map((nm, i) => <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}><td style={{ ...S.td, fontWeight: 600 }}><Tx m={nm} /></td><td style={{ ...S.td, fontWeight: 700 }}>{ff(p2.b[i], 4)}</td><td style={S.td}>{ff(p2.se[i], 4)}</td></tr>)}</tbody></table>
        <div style={S.fm}><TB m={`\\hat{y}=${ff(p2.b[0], 4)}${p2.b[1] >= 0 ? "+" : ""}${ff(p2.b[1], 4)}x_1${p2.b[2] >= 0 ? "+" : ""}${ff(p2.b[2], 4)}x_2${p2.b[3] >= 0 ? "+" : ""}${ff(p2.b[3], 4)}x_3${p2.b[4] >= 0 ? "+" : ""}${ff(p2.b[4], 4)}x_4`} /></div>
        <div style={S.ok}><b>Respuesta:</b> El modelo ajustado por MCO para el uso de CPU de SwiftPay queda definido por la ecuación anterior.</div>
      </div>
    </div>),

    () => (<div>
      <div style={S.q}>Pregunta 2: Prediga CPU para <Tx m={`x_1=${xp[0]},\\,x_2=${xp[1]},\\,x_3=${xp[2]},\\,x_4=${xp[3]}`} />.</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Predicción Puntual</h3>
        <div style={S.fm}><TB m={"\\hat{y}=\\hat{\\beta}_0+\\hat{\\beta}_1x_1+\\hat{\\beta}_2x_2+\\hat{\\beta}_3x_3+\\hat{\\beta}_4x_4"} /></div>
        <div style={S.fm}><TB m={`\\hat{y}=${ff(p2.b[0], 4)}+(${ff(p2.b[1], 4)})(${xp[0]})+(${ff(p2.b[2], 4)})(${xp[1]})+(${ff(p2.b[3], 4)})(${xp[2]})+(${ff(p2.b[4], 4)})(${xp[3]})`} /></div>
        <div style={S.ok}><TB m={`\\boxed{\\hat{y}=${ff(p2.ypr, 4)}\\%}`} /><br /><b>Respuesta:</b> Bajo las condiciones dadas, el modelo predice un uso de CPU del {ff(p2.ypr, 2)}% en los servidores de SwiftPay.</div>
      </div>
    </div>),

    () => (<div>
      <div style={S.q}>Pregunta 3: ANOVA para bondad de ajuste del modelo.</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>ANOVA para RLM — 8 Pasos</h3>
        <Def show={D}><b>Prueba global:</b> <Tx m={"H_0: \\beta_1 = \\cdots = \\beta_p = 0"} />. Si <Tx m={"F_0 > F_c"} />, al menos una variable explica <Tx m="Y" />.</Def>
        <div style={S.sp}><div style={S.spT}>Pasos 1-4</div>
          <p><Tx m={`n=20,\\;p=4,\\;gl_{reg}=4,\\;gl_{error}=15`} />. <Tx m={`\\alpha=${a2}`} /></p>
          <div style={S.fm}><TB m={"H_0:\\beta_1=\\beta_2=\\beta_3=\\beta_4=0\\qquad H_1:\\text{Al menos un }\\beta_i\\neq 0"} /></div>
        </div>
        <div style={S.sp}><div style={S.spT}>Pasos 5-6</div>
          <div style={S.fm}><TB m={`F_0=\\frac{MSR}{MSE}\\sim F_{(4,15)},\\qquad\\text{Rechazar si }F_0>${p2.Fc}`} /></div>
          <Def show={D}><b>MSR</b> = SSR/p: varianza explicada por regresión.<br /><b>MSE</b> = SSE/(n−p−1): varianza residual. <Tx m={"F_0 = MSR/MSE"} />.</Def>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 7</div>
          <h4 style={{ margin: "4px 0 3px", color: "#333", fontSize: 13 }}>Fórmulas</h4>
          <table style={S.T}><thead><tr><th style={S.thF}>Fuente</th><th style={S.thF}>SS</th><th style={S.thF}>gl</th><th style={S.thF}>MS</th></tr></thead>
            <tbody>
              <tr><td style={{ ...S.td, fontWeight: 600 }}>Regresión</td><td style={S.td}><Tx m={"\\sum(\\hat{y}_i-\\bar{y})^2"} /></td><td style={S.td}><Tx m="p" /></td><td style={S.td}><Tx m="SSR/p" /></td></tr>
              <tr style={{ background: "#f9f9f9" }}><td style={{ ...S.td, fontWeight: 600 }}>Error</td><td style={S.td}><Tx m={"\\sum(y_i-\\hat{y}_i)^2"} /></td><td style={S.td}><Tx m="n{-}p{-}1" /></td><td style={S.td}><Tx m="SSE/(n{-}p{-}1)" /></td></tr>
              <tr><td style={{ ...S.td, fontWeight: 600 }}>Total</td><td style={S.td}><Tx m={"\\sum(y_i-\\bar{y})^2"} /></td><td style={S.td}><Tx m="n{-}1" /></td><td style={S.td}>—</td></tr>
            </tbody></table>
          <h4 style={{ margin: "8px 0 3px", color: "#333", fontSize: 13 }}>Valores</h4>
          <table style={S.T}><thead><tr><th style={S.th}>Fuente</th><th style={S.th}>SS</th><th style={S.th}>gl</th><th style={S.th}>MS</th><th style={S.th}><Tx m="F_0" /></th></tr></thead>
            <tbody>
              <tr><td style={{ ...S.td, fontWeight: 600 }}>Regresión</td><td style={S.td}>{ff(p2.SSR, 4)}</td><td style={S.td}>4</td><td style={S.td}>{ff(p2.MSR, 4)}</td><td style={{ ...S.td, fontWeight: 700, color: "#c62828" }}>{f2(p2.Fr)}</td></tr>
              <tr style={{ background: "#f9f9f9" }}><td style={{ ...S.td, fontWeight: 600 }}>Error</td><td style={S.td}>{ff(p2.SSE, 4)}</td><td style={S.td}>15</td><td style={S.td}>{ff(p2.MSE, 4)}</td><td style={S.td}>—</td></tr>
              <tr style={{ fontWeight: 600 }}><td style={{ ...S.td, fontWeight: 700 }}>Total</td><td style={S.td}>{ff(p2.SST, 4)}</td><td style={S.td}>19</td><td style={S.td}>—</td><td style={S.td}>—</td></tr>
            </tbody></table>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 8</div>
          <div style={S.ok}><Tx m={`F_0=${f2(p2.Fr)}\\;${p2.rejectedGlobal ? ">" : "\\leq"}\\;${p2.Fc}`} />.{" "}
            {p2.rejectedGlobal ? "Se rechaza H₀: el modelo es significativo. Al menos un parámetro de tráfico aporta a explicar el uso de CPU de SwiftPay." : "No se rechaza H₀: el modelo no es significativo globalmente."}</div>
        </div>
      </div>
    </div>),

    () => (<div>
      <div style={S.q}>Pregunta 4: Coeficiente de determinación e interpretación.</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Coeficiente de Determinación</h3>
        <h4 style={{ margin: "4px 0 3px", color: "#333", fontSize: 13 }}>Fórmulas</h4>
        <div style={S.fm}><TB m={"R^2=\\frac{SSR}{SST}=1-\\frac{SSE}{SST}"} /><TB m={"R^2_{\\text{adj}}=1-\\frac{(1-R^2)(n-1)}{n-p-1}"} /></div>
        <Def show={D}><b><Tx m="R^2" />:</b> <Tx m={"SSR/SST \\in [0,1]"} /> — fracción de variabilidad explicada por el modelo.<br /><b><Tx m={"R^2_{\\text{adj}}"} />:</b> Penaliza por número de parámetros p; favorece parsimonia al comparar modelos.</Def>
        <h4 style={{ margin: "8px 0 3px", color: "#333", fontSize: 13 }}>Cálculo</h4>
        <div style={S.fm}><TB m={`R^2=\\frac{${ff(p2.SSR, 4)}}{${ff(p2.SST, 4)}}=${ff(p2.R2, 4)}`} /><TB m={`R^2_{\\text{adj}}=${ff(p2.R2a, 4)}`} /></div>
        <div style={S.ok}><b>Respuesta:</b> El modelo explica el <b>{ff(p2.R2 * 100, 2)}%</b> de la variabilidad del uso de CPU. <Tx m={`R^2_{\\text{adj}}=${ff(p2.R2a * 100, 2)}\\%`} />. El {ff((1 - p2.R2) * 100, 2)}% restante corresponde a factores no incluidos o variabilidad aleatoria.</div>
      </div>
    </div>),

    () => (<div>
      <div style={S.q}>Pregunta 5: IC del 95% para cada <Tx m="\beta_i" />.</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Intervalos de Confianza</h3>
        <div style={S.fm}><TB m={`IC_{${(1 - a2) * 100}\\%}:\\;\\hat{\\beta}_i\\pm t_{${a2 / 2},\\,15}\\cdot SE(\\hat{\\beta}_i),\\quad t_{${a2 / 2},15}=${p2.tc}`} /></div>
        <Def show={D}><b>IC:</b> <Tx m={"\\hat{\\beta}_i \\pm t_{\\alpha/2,\\,n-p-1} \\cdot SE(\\hat{\\beta}_i)"} />. Si <Tx m={"0 \\in IC"} /> → <Tx m="\beta_i" /> no es significativamente distinto de 0.</Def>
        <table style={S.T}><thead><tr><th style={S.th}></th><th style={S.th}><Tx m={"\\hat{\\beta}_i"} /></th><th style={S.th}>SE</th><th style={S.th}>Inferior</th><th style={S.th}>Superior</th><th style={S.th}><Tx m={"0\\in IC"} /></th></tr></thead>
          <tbody>{bN.map((nm, i) => { const c0 = p2.cl[i] <= 0 && p2.ch[i] >= 0; return <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}><td style={{ ...S.td, fontWeight: 600 }}><Tx m={nm} /></td><td style={S.td}>{ff(p2.b[i], 4)}</td><td style={S.td}>{ff(p2.se[i], 4)}</td><td style={S.td}>{ff(p2.cl[i], 4)}</td><td style={S.td}>{ff(p2.ch[i], 4)}</td><td style={{ ...S.td, fontWeight: 700, color: c0 ? "#c62828" : "#2e7d32" }}>{c0 ? "Sí" : "No"}</td></tr>; })}</tbody></table>
        <div style={S.ok}><b>Respuesta:</b> Los parámetros cuyo IC no contiene 0 son significativos al {(1 - a2) * 100}% de confianza para el modelo de CPU de SwiftPay.</div>
      </div>
    </div>),

    () => (<div>
      <div style={S.q}>Pregunta 6: ¿El aporte de cada <Tx m="X_i" /> es significativo (5%)?</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Significancia Individual — 8 Pasos</h3>
        <Def show={D}><b>Significancia individual:</b> Prueba <Tx m={"H_0: \\beta_i = 0"} /> para cada variable por separado. Rechazar si <Tx m={"|t_0| > t_c"} />.</Def>
        <div style={S.sp}><div style={S.spT}>Pasos 1-6</div>
          <div style={S.fm}><TB m={`H_0:\\beta_i=0\\quad H_1:\\beta_i\\neq 0\\quad\\alpha=${a2}`} /><TB m={`t_0=\\frac{\\hat{\\beta}_i}{SE(\\hat{\\beta}_i)},\\quad\\text{Rechazar si }|t_0|>t_{${a2 / 2},15}=${p2.tc}`} /></div>
          <Def show={D}><b><Tx m="t_0" />:</b> <Tx m={"t_0 = \\hat{\\beta}_i / SE(\\hat{\\beta}_i)"} /> — distancia estandarizada de <Tx m={"\\hat{\\beta}_i"} /> a 0. Rechazar si <Tx m={`|t_0| > ${p2.tc}`} />.</Def>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 7: Cálculos</div>
          <table style={S.T}><thead><tr><th style={S.th}>Variable</th><th style={S.th}><Tx m={"\\hat{\\beta}_i"} /></th><th style={S.th}>SE</th><th style={S.th}><Tx m="t_0" /></th><th style={S.th}><Tx m={`|t_0|\\text{ vs }${p2.tc}`} /></th><th style={S.th}>Decisión</th></tr></thead>
            <tbody>{bN.slice(1).map((nm, i) => { const ix = i + 1, sg = Math.abs(p2.tv[ix]) > p2.tc; return <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}><td style={{ ...S.td, fontWeight: 600 }}><Tx m={nm} /> ({bD[ix]})</td><td style={S.td}>{ff(p2.b[ix], 4)}</td><td style={S.td}>{ff(p2.se[ix], 4)}</td><td style={{ ...S.td, fontWeight: 700 }}>{ff(p2.tv[ix], 4)}</td><td style={S.td}>{ff(Math.abs(p2.tv[ix]), 4)} {sg ? ">" : "≤"} {p2.tc}</td><td style={{ ...S.td, fontWeight: 700, color: sg ? "#2e7d32" : "#c62828" }}>{sg ? "Signif." : "No signif."}</td></tr>; })}</tbody></table>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 8: Conclusión</div>
          <div style={S.ok}><b>Respuesta:</b> Con <Tx m={`\\alpha=${a2}`} />, las variables significativas son: <b>{p2.sigVars.length > 0 ? p2.sigVars.map(i => bD[i]).join(", ") : "ninguna"}</b>. Las no significativas ({p2.nsVars.map(i => bD[i]).join(", ")}) podrían eliminarse del modelo de telemetría de SwiftPay.</div>
          <Def show={D}><b>Parsimonia:</b> Modelo con menor número de variables que maximice <Tx m={"R^2_{\\text{adj}}"} />.</Def>
        </div>
      </div>
    </div>),

    () => (<div>
      <div style={S.q}>Pregunta 7: ¿Cuál modelo recomendaría? Justifique.</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Modelo Recomendado</h3>
        <table style={S.T}><thead><tr><th style={S.th}>Variable</th><th style={S.th}><Tx m="|t_0|" /></th><th style={S.th}>vs <Tx m={`${p2.tc}`} /></th><th style={S.th}>Decisión</th></tr></thead>
          <tbody>{[1, 2, 3, 4].map(i => { const sg = Math.abs(p2.tv[i]) > p2.tc; return <tr key={i} style={{ background: i % 2 ? "" : "#f9f9f9" }}><td style={{ ...S.td, fontWeight: 600 }}><Tx m={`x_${i}`} /> ({bD[i]})</td><td style={S.td}>{ff(Math.abs(p2.tv[i]), 4)}</td><td style={S.td}>{sg ? ">" : "≤"} {p2.tc}</td><td style={{ ...S.td, fontWeight: 700, color: sg ? "#2e7d32" : "#c62828" }}>{sg ? "Incluir" : "Eliminar"}</td></tr>; })}</tbody></table>
        <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>Modelo reducido propuesto</h4>
        <div style={S.fm}><TB m={`\\hat{y}=\\beta_0${p2.sigVars.map(i => `+\\beta_${i}x_${i}`).join("")}`} /></div>
        <Def show={D}><b>Reducción:</b> Eliminar variables con <Tx m={"|t_0| \\leq t_c"} /> reduce sobreajuste. Verificar con F parcial: <Tx m={"F = \\frac{(SSE_R - SSE_C)/q}{MSE_C}"} />, q = variables eliminadas.</Def>
        <div style={S.ok}><b>Respuesta:</b> Se recomienda un <b>modelo reducido</b> con {p2.sigVars.length > 0 ? <>variables <b>{p2.sigVars.map(i => bD[i]).join(", ")}</b></> : "revisión adicional"}. Las variables {p2.nsVars.map(i => bD[i]).join(", ")} se eliminan. Modelo completo: <Tx m={`R^2=${ff(p2.R2 * 100, 2)}\\%`} />, <Tx m={`R^2_{\\text{adj}}=${ff(p2.R2a * 100, 2)}\\%`} />. El modelo reducido debería mantener un <Tx m={"R^2_{\\text{adj}}"} /> similar, confirmando la simplificación.</div>
      </div>
    </div>),
  ] : [() => <div style={S.wn}>Error en los datos de la Parte II. Verifique los valores ingresados.</div>];

  return (
    <div style={S.bx}>
      <div style={S.hd}>
        <h1 style={S.h1}>Proyecto Final — Estadística Inferencial</h1>
        <p style={S.h2}>Banco de Venezuela — Arquitecturas SwiftPay · Computación FACYT 2026</p>
      </div>
      <div style={S.tbar}>
        <button style={S.mt(tab === "p1")} onClick={() => setTab("p1")}>PARTE I: ANOVA</button>
        <button style={S.mt(tab === "p2")} onClick={() => setTab("p2")}>PARTE II: Regresión</button>
      </div>
      <div style={{ background: "#fff", border: "1px solid #ddd", borderTop: "none", borderRadius: "0 0 8px 8px", padding: 12 }}>
        {edit && (tab === "p1" ? <Ed1 /> : <Ed2 />)}
        {tab === "p1" && (<><div style={S.sb}>{t1.map((t, i) => <button key={i} style={S.st(s1 === i)} onClick={() => setS1(i)}>{t}</button>)}</div>{P1[s1]()}</>)}
        {tab === "p2" && (<><div style={S.sb}>{t2.map((t, i) => <button key={i} style={S.st(s2 === i)} onClick={() => setS2(i)}>{t}</button>)}</div>{P2[s2]()}</>)}
      </div>
      <div style={{ textAlign: "center", padding: "10px 0", color: "#999", fontSize: 10 }}>Estadística Inferencial — FACYT 2026</div>
      <div style={S.fab}>
        <button style={S.fabBtn(edit ? "#f59e0b" : "#92400e")} onClick={() => setEdit(e => !e)}>{edit ? "✕ Cerrar editor" : "⚙️ Editar datos"}</button>
        <button style={S.fabBtn(defs ? "#6366f1" : "#1e1b4b")} onClick={() => setDefs(d => !d)}>{defs ? "✕ Ocultar defs" : "📖 Definiciones"}</button>
      </div>
    </div>
  );
}