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

// ===== F-distribution PDF & CDF =====
function lnGamma(z) { const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]; if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lnGamma(1 - z); z -= 1; let x = c[0]; for (let i = 1; i < 9; i++) x += c[i] / (z + i); const t = z + 7.5; return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x); }
function fPDF(x, d1, d2) { if (x <= 0) return 0; const lnB = lnGamma(d1 / 2) + lnGamma(d2 / 2) - lnGamma((d1 + d2) / 2); return Math.exp((d1 / 2) * Math.log(d1 / d2) + (d1 / 2 - 1) * Math.log(x) - ((d1 + d2) / 2) * Math.log(1 + d1 * x / d2) - lnB); }
function betacf(x, a, b) { var m, m2, num; var qab = a + b, qap = a + 1, qam = a - 1; var c = 1, d = 1 - qab * x / qap; if (Math.abs(d) < 1e-30) d = 1e-30; d = 1 / d; var h = d; for (m = 1; m <= 100; m++) { m2 = 2 * m; num = m * (b - m) * x / ((qam + m2) * (a + m2)); d = 1 + num * d; if (Math.abs(d) < 1e-30) d = 1e-30; c = 1 + num / c; if (Math.abs(c) < 1e-30) c = 1e-30; d = 1 / d; h *= d * c; num = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2)); d = 1 + num * d; if (Math.abs(d) < 1e-30) d = 1e-30; c = 1 + num / c; if (Math.abs(c) < 1e-30) c = 1e-30; d = 1 / d; var del = d * c; h *= del; if (Math.abs(del - 1) < 3e-7) break; } return h; }
function betai(a, b, x) { var bt; if (x === 0 || x === 1) bt = 0; else bt = Math.exp(lnGamma(a + b) - lnGamma(a) - lnGamma(b) + a * Math.log(x) + b * Math.log(1 - x)); if (x < (a + 1) / (a + b + 2)) return bt * betacf(x, a, b) / a; else return 1 - bt * betacf(1 - x, b, a) / b; }
function fPValue(F, d1, d2) { if (F <= 0) return 1; return betai(d2 / 2, d1 / 2, d2 / (d2 + d1 * F)); }
// ===== Normal CDF =====
function normCDF(x) { var t = 1 / (1 + 0.2316419 * Math.abs(x)); var d = 0.3989423 * Math.exp(-x * x / 2); var p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274)))); if (x > 0) p = 1 - p; return p; }

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
  const cmp = []; for (let i = 0; i < k; i++)for (let j = i + 1; j < k; j++) { if (nm[i] === "SwiftVen" && nm[j] === "SwiftFast") continue; const df = Math.abs(mn[i] - mn[j]); cmp.push({ a: nm[i], b: nm[j], i: i + 1, j: j + 1, diff: df, sig: df > LSD }); }
  const vr = data.map((g, i) => g.reduce((a, v) => a + (v - mn[i]) ** 2, 0) / (n - 1));
  const res = []; data.forEach((g, i) => g.forEach(v => res.push(v - mn[i])));
  const resSorted = [...res].sort((a, b) => a - b);
  // Shapiro-Wilk for N=12
  const a_SW = [0.5475, 0.3325, 0.2347, 0.1586, 0.0922, 0.0303];
  let b_SW = 0;
  const sw_steps = [];
  for (let i = 0; i < 6; i++) {
    const diff = resSorted[11 - i] - resSorted[i];
    const prod = a_SW[i] * diff;
    b_SW += prod;
    sw_steps.push({ i: i + 1, a: a_SW[i], xn: resSorted[11 - i], x1: resSorted[i], diff, prod });
  }
  const W_calc = (b_SW ** 2) / SSE;
  const lnN = Math.log(N);
  const mu_sw = 0.0038915 * Math.pow(lnN, 3) - 0.083751 * Math.pow(lnN, 2) - 0.31082 * lnN - 1.5861;
  const sig_sw = Math.exp(0.0030302 * Math.pow(lnN, 2) - 0.082676 * lnN - 0.4803);
  const z_sw = W_calc < 1 ? (Math.log(1 - W_calc) - mu_sw) / sig_sw : -99;
  const pval_SW = 1 - normCDF(z_sw);
  // Runs test on residuals in observation order
  let runs = 1; for (let i = 1; i < res.length; i++) if ((res[i] >= 0) !== (res[i - 1] >= 0)) runs++;
  const n1 = res.filter(e => e >= 0).length, n2 = res.filter(e => e < 0).length;
  const ER = (2 * n1 * n2) / (n1 + n2) + 1;
  const VR = (2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)) / ((n1 + n2) ** 2 * ((n1 + n2) - 1));
  const sdR = Math.sqrt(VR);
  const Z0runs = sdR > 0 ? (runs - ER) / sdR : 0;
  const Zcrit = alpha === 0.01 ? 2.576 : alpha === 0.05 ? 1.96 : 1.645;
  const runsOk = Math.abs(Z0runs) <= Zcrit;
  const pval = fPValue(F0, k - 1, N - k);

  // Kruskal-Wallis Test
  const kw_data = [];
  data.forEach((g, i) => g.forEach(v => kw_data.push({ v, g: i, nm: nm[i] })));
  kw_data.sort((a, b) => a.v - b.v);

  let kw_ties = [];
  let kw_idx = 0;
  while (kw_idx < kw_data.length) {
    let kw_j = kw_idx;
    while (kw_j < kw_data.length && kw_data[kw_j].v === kw_data[kw_idx].v) kw_j++;
    const t_len = kw_j - kw_idx;
    if (t_len > 1) kw_ties.push(t_len);
    const avgR = kw_idx + 1 + (t_len - 1) / 2;
    for (let k_i = kw_idx; k_i < kw_j; k_i++) kw_data[k_i].r = avgR;
    kw_idx = kw_j;
  }

  const kw_R = new Array(k).fill(0);
  kw_data.forEach(d => { kw_R[d.g] += d.r; });
  let kw_sum_R2_n = 0;
  for (let grp = 0; grp < k; grp++) kw_sum_R2_n += (kw_R[grp] ** 2) / data[grp].length;

  const kw_H = (12 / (N * (N + 1))) * kw_sum_R2_n - 3 * (N + 1);
  let kw_tieSum = 0;
  kw_ties.forEach(t => { kw_tieSum += (Math.pow(t, 3) - t); });
  const kw_tieCorrection = 1 - (kw_tieSum / (Math.pow(N, 3) - N));
  const kw_H_adj = kw_tieSum > 0 ? kw_H / kw_tieCorrection : kw_H;

  const chiSquareDf2 = { 0.01: 9.210, 0.05: 5.991, 0.10: 4.605 };
  const kw_chiCrit = chiSquareDf2[alpha] || 5.991;
  const kw_rejected = kw_H_adj > kw_chiCrit;

  return { data, nm, k, n, N, sigma, alpha, mn, gm, tau, SST, SSE, SSTot, MST, MSE, F0, Fc, tc, LSD, cmp, vr, res, resSorted, runs, n1, n2, ER, VR, sdR, Z0runs, Zcrit, runsOk, rejected: F0 > Fc, pval, a_SW, b_SW, sw_steps, W_calc, pval_SW, kw_data, kw_ties, kw_R, kw_H, kw_tieSum, kw_tieCorrection, kw_H_adj, kw_chiCrit, kw_rejected };
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
  const pvalGlobal = fPValue(Fr, p, df);
  const pvals = tv.map((t, i) => i === 0 ? 0 : fPValue(t * t, 1, df));
  return { raw, no, p, b, se, tv, SSR, SSE, SST, MSR, MSE, Fr, R2, R2a, ym, df, tc, Fc, cl, ch, ypr, sigVars, nsVars, rejectedGlobal: Fr > Fc, alpha, pvalGlobal, pvals };
}

// ===== Number input =====
function NI({ value, onChange, w = 58 }) {
  return <input type="text" inputMode="decimal" value={value !== undefined ? value : ""} onChange={e => onChange(e.target.value)}
    style={{ width: w, padding: "3px 4px", border: "1px solid #cbd5e1", borderRadius: 4, fontSize: 13, textAlign: "center", background: "#fffbeb", color: "#1e293b", fontWeight: 600, outline: "none" }} />;
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

const t1 = ["1. Modelo", "2. ANOVA", "3. Comparaciones", "4. Normalidad", "5. Aleatoriedad", "6. Varianzas", "7. Kruskal-Wallis"];
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
  ed: { background: "#f0f9ff", border: "1px solid #0284c7", borderRadius: 8, padding: 16, marginBottom: 16, boxShadow: "inset 0 2px 4px rgba(0,0,0,0.05)" },
  edH: { display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 12, borderBottom: "1px solid #bae6fd", paddingBottom: 8 },
  btn: (bg) => ({ padding: "6px 14px", border: "none", borderRadius: 5, cursor: "pointer", fontSize: 12, fontWeight: 600, background: bg, color: "#fff", transition: "0.2s" }),
  sel: { padding: "5px 8px", border: "1px solid #94a3b8", borderRadius: 4, fontSize: 12, background: "#fff" },
  fab: { position: "fixed", bottom: 18, right: 18, display: "flex", gap: 6, zIndex: 999 },
  fabBtn: (bg) => ({ padding: "7px 12px", borderRadius: 18, border: "none", cursor: "pointer", fontSize: 11, fontWeight: 600, boxShadow: "0 2px 8px rgba(0,0,0,0.2)", background: bg, color: "#fff" }),
};

function DataEditor1({ initialD, initialSig, initialA, onApply, onCancel }) {
  const [d, setD] = useState(() => initialD.map(r => [...r]));
  const [sig, setSig] = useState(initialSig);
  const [a, setA] = useState(initialA);

  const upd = (i, j, v) => { const n = [...d]; n[i] = [...n[i]]; n[i][j] = v; setD(n); };
  const reset = () => { setD(DEF1.map(r => [...r])); setSig(0.18); setA(0.05); };

  return (
    <div style={S.ed}>
      <div style={S.edH}>
        <b style={{ color: "#0369a1", fontSize: 14 }}>Modo de Edición — Parte I</b>
        <div style={{ display: "flex", gap: 8 }}>
          <button style={S.btn("#dc2626")} onClick={reset}>Restablecer a Fábrica</button>
          <button style={S.btn("#64748b")} onClick={onCancel}>Cancelar</button>
          <button style={S.btn("#0ea5e9")} onClick={() => onApply(d, sig, a)}>💾 Aplicar y Recalcular</button>
        </div>
      </div>
      <div style={{ display: "flex", gap: 16, flexWrap: "wrap", alignItems: "center", marginBottom: 8 }}>
        <label style={{ fontSize: 12 }}><Tx m="\sigma" /> = <NI value={sig} onChange={setSig} w={55} /></label>
        <label style={{ fontSize: 12 }}><Tx m="\alpha" /> = <select style={S.sel} value={a} onChange={e => setA(parseFloat(e.target.value))}>{ALPHAS.map(x => <option key={x} value={x}>{x}</option>)}</select></label>
      </div>
      <table style={S.T}>
        <thead><tr><th style={S.th}>Arquitectura</th>{[1, 2, 3, 4].map(j => <th key={j} style={S.th}>Represión {j}</th>)}</tr></thead>
        <tbody>
          {d.map((row, i) => <tr key={i}><td style={{ ...S.td, fontWeight: 600 }}>{arqN[i]}</td>{row.map((v, j) => <td key={j} style={S.td}><NI value={v} onChange={x => upd(i, j, x)} w={52} /></td>)}</tr>)}
        </tbody>
      </table>
    </div>
  );
}

function DataEditor2({ initialD, initialA, initialXp, onApply, onCancel }) {
  const [d, setD] = useState(() => initialD.map(r => [...r]));
  const [a, setA] = useState(initialA);
  const [xp, setXp] = useState([...initialXp]);

  const upd = (i, j, v) => { const n = [...d]; n[i] = [...n[i]]; n[i][j] = v; setD(n); };
  const updXp = (i, v) => { const n = [...xp]; n[i] = v; setXp(n); };
  const reset = () => { setD(DEF2.map(r => [...r])); setA(0.05); setXp([...DEFXP]); };

  return (
    <div style={S.ed}>
      <div style={S.edH}>
        <b style={{ color: "#0369a1", fontSize: 14 }}>Modo de Edición — Parte II</b>
        <div style={{ display: "flex", gap: 8 }}>
          <button style={S.btn("#dc2626")} onClick={reset}>Restablecer</button>
          <button style={S.btn("#64748b")} onClick={onCancel}>Cancelar</button>
          <button style={S.btn("#0ea5e9")} onClick={() => onApply(d, a, xp)}>💾 Aplicar y Recalcular</button>
        </div>
      </div>
      <div style={{ display: "flex", gap: 16, flexWrap: "wrap", alignItems: "center", marginBottom: 8 }}>
        <label style={{ fontSize: 12 }}><Tx m="\alpha" /> = <select style={S.sel} value={a} onChange={e => setA(parseFloat(e.target.value))}>{ALPHAS.map(x => <option key={x} value={x}>{x}</option>)}</select></label>
        <span style={{ fontSize: 12 }}>Predicción:</span>
        {xp.map((v, i) => <label key={i} style={{ fontSize: 12 }}><Tx m={`x_${i + 1}`} />=<NI value={v} onChange={x => updXp(i, x)} w={48} /></label>)}
      </div>
      <div style={{ maxHeight: 220, overflowY: "auto", border: "1px solid #e5e7eb", borderRadius: 4 }}>
        <table style={{ ...S.T, marginTop: 0 }}>
          <thead style={{ position: "sticky", top: 0 }}><tr><th style={S.th}>#</th><th style={S.th}>y</th><th style={S.th}><Tx m="x_1" /></th><th style={S.th}><Tx m="x_2" /></th><th style={S.th}><Tx m="x_3" /></th><th style={S.th}><Tx m="x_4" /></th></tr></thead>
          <tbody>{d.map((row, i) => <tr key={i}><td style={{ ...S.td, fontWeight: 600, fontSize: 10 }}>{i + 1}</td>{row.map((v, j) => <td key={j} style={S.td}><NI value={v} onChange={x => upd(i, j, x)} w={46} /></td>)}</tr>)}</tbody>
        </table>
      </div>
    </div>
  );
}

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

  const p1 = useMemo(() => calc1(d1.map(r => r.map(v => parseFloat(v) || 0)), parseFloat(sig) || 0, parseFloat(a1) || 0.05), [d1, sig, a1]);
  const p2 = useMemo(() => calc2(d2.map(r => r.map(v => parseFloat(v) || 0)), parseFloat(a2) || 0.05, xp.map(v => parseFloat(v) || 0)), [d2, a2, xp]);

  const handleApply1 = (nd, nsig, na) => { setD1(nd); setSig(nsig); setA1(na); setEdit(false); };
  const handleApply2 = (nd, na, nxp) => { setD2(nd); setA2(na); setXp(nxp); setEdit(false); };

  const D = defs;

  if (!kr) return <div style={{ textAlign: "center", padding: 40 }}>Cargando LaTeX...</div>;

  // ===== PART 1 PANELS =====
  const P1 = [
    () => (<div>
      <div style={S.q}>Pregunta 1: Plantee el modelo estadístico y estime <Tx m="\mu" />, <Tx m="\sigma^2" /> y <Tx m="\tau_i" />.</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Modelo Estadístico</h3>
        <div style={S.fm}><TB m={`Y_{ij} = \\mu + \\tau_i + \\varepsilon_{ij}, \\qquad \\varepsilon_{ij} \\sim N(0,\\,\\sigma^2)`} /></div>
        <Def show={D}>
          <div style={{ paddingBottom: 6 }}><b><Tx m="Y_{ij}" /> (Variable de Respuesta):</b> Es el valor observado en la <Tx m="j" />-ésima réplica bajo el <Tx m="i" />-ésimo tratamiento. Mide directamente el fenómeno de interés (latencia/rendimiento). Todo el análisis se centra en modelar la variabilidad observada en esta métrica matemática y descomponerla en efectos sistemáticos y ruido puramente estocástico.</div>
          <div style={{ paddingBottom: 6 }}><b><Tx m="\mu" /> (Media Global Verdadera):</b> Es el valor esperado inobservable de la métrica a través de todos los tratamientos. Representa un nivel "base" intrínseco del entorno antes de que cualquier arquitectura ejerza su influencia, modelado como <Tx m="\mu = \tfrac{1}{k}\sum \mu_i" />.</div>
          <div style={{ paddingBottom: 6 }}><b><Tx m="\tau_i" /> (Efecto Principal):</b> Expresa rigurosamente la contribución aislada del tratamiento <Tx m="i" /> (<Tx m="\tau_i = \mu_i - \mu" />). Refleja en qué medida la media local se desvía del gran promedio. La restricción matricial <Tx m="\sum \tau_i = 0" /> asegura que estos efectos sean meros desplazamientos ortogonales al intercepto general, lo cual es vital empíricamente para garantizar la solución convexa del modelo lineal subyacente.</div>
          <div style={{ paddingBottom: 6 }}><b><Tx m="\varepsilon_{ij}" /> (Término de Error):</b> Mide factores estocásticos que desvían cada iteración de la media esperada (<Tx m="\varepsilon_{ij} = Y_{ij} - \mu_i" />). Para la validez formal matemática, se deben comportar como ruido blanco gaussiano independiente e idénticamente distribuido (i.i.d.) <Tx m="N(0,\sigma^2)" />. Encapsulan la cuota de ignorancia estocástica del analista experimental.</div>
          <div><b>Modelo DCA (Diseño Completamente Aleatorizado):</b> Arquitectura estadística minimalista asumiendo homocedasticidad universal y aleatorización plena, lo cual inmuniza matemáticamente al diseño contra variables de confusión endógenas.</div>
        </Def>
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
        <table style={S.T}><thead><tr><th style={S.th}>Arquitectura</th>{[1, 2, 3, 4].map(j => <th key={j} style={S.th}>{j}</th>)}<th style={S.th}><Tx m="\bar{y}_i" /></th><th style={S.th}><Tx m="\hat{\tau}_i = \bar{y}_i - \bar{y}_{..}" /></th></tr></thead>
          <tbody>{p1.nm.map((nm, i) => <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}><td style={{ ...S.td, fontWeight: 600 }}>{nm} (<Tx m={`A_{${i + 1}}`} />)</td>{p1.data[i].map((v, j) => <td key={j} style={S.td}>{ff(v, 2)}</td>)}<td style={{ ...S.td, fontWeight: 700, color: "#0f3460" }}>{ff(p1.mn[i], 4)}</td><td style={{ ...S.td, fontWeight: 600, color: p1.tau[i] > 0 ? "#c62828" : "#2e7d32" }}>{p1.tau[i] > 0 ? "+" : ""}{ff(p1.tau[i], 4)}</td></tr>)}</tbody></table>
        <div style={S.fm}>
          <TB m={`\\hat{\\mu} = \\bar{y}_{..} = \\frac{\\sum y_{ij}}{N} = \\frac{${ff(p1.gm * p1.N, 2)}}{${p1.N}} = ${ff(p1.gm, 4)}\\text{ seg}`} />
          <TB m={`\\sigma^2 = (${ff(sig, 2)})^2 = ${ff(sig * sig, 4)}\\text{ seg}^2,\\quad \\sigma = ${ff(sig, 2)}\\text{ seg}`} />
        </div>
        <div style={{ ...S.fm, marginTop: 4, padding: "6px 12px", background: "#f5f7fa" }}>
          <Tx m={`\\text{Ejemplo } \\hat{\\tau}_1: \\bar{y}_1 - \\hat{\\mu} = ${ff(p1.mn[0], 4)} - ${ff(p1.gm, 4)} = ${ff(p1.tau[0], 4)}\\text{ seg}`} />
        </div>
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
          <Def show={D}><b><Tx m="\mu_i" /> (Media Marginal del Tratamiento):</b> Expresa el verdadero valor empírico proyectado de la variable bajo el tratamiento <Tx m="i" />. Matemáticamente, <Tx m="\mu_i = \mu + \tau_i" />. Es una cantidad espectral e inobservable que representa el verdadero desempeño latente. Las medias muestrales <Tx m="\bar{y}_i" /> provienen iterativamente de estas esperanzas. El propósito ontológico del ANOVA no es otro que detectar si la varianza generada en el espectro entre estas <Tx m="\mu_i" /> es igual a cero.</Def>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 2: Hipótesis nula</div>
          <div style={S.fm}><TB m="H_0: \\mu_1 = \\mu_2 = \\mu_3" /></div><p>Las tres arquitecturas tienen igual latencia media.</p>
          <Def show={D}><b><Tx m="H_0" /> (Hipótesis Nula Estructural):</b> Postula la total y prístina equivalencia de todos los entornos observados: <Tx m="\mu_1 = \mu_2 = \dots = \mu_k" />, colapsando a <Tx m="\tau_i = 0 \;\forall i" />. En inferencia de Neyman-Pearson, tomamos esta igualdad como la distribución asintótica por defecto, creando un contrafactual estadístico univariante. Solo la evidencia muestral anómala —desviaciones irreconciliables frente al azar— podrá forzar nuestro rechazo axiomático absoluto hacia <Tx m="H_1" />.</Def>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 3: Hipótesis alternativa</div>
          <div style={S.fm}><TB m="H_1: \\text{Al menos un }\\mu_i\\text{ difiere}" /></div>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 4: Nivel de significancia</div>
          <div style={S.fm}><TB m={`\\alpha = ${a1}`} /></div>
          <Def show={D}><b><Tx m="\alpha" /> (Tolerancia al Error Tipo I):</b> Cota de falso positivo que el investigador preselecciona: <Tx m={"P(\\text{rechazar } H_0 \\mid H_0 \\text{ es verdadera}) \\leq \\alpha"} />. Al fijar este radio en {a1 * 100}%, toleramos explícitamente cometer perjurio matemático donde señalamos una diferencia técnica inexistente. Controla la severidad integral del estudio.</Def>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 5: Estadístico</div>
          <div style={S.fm}><TB m={`F_0 = \\frac{MST}{MSE} \\sim F_{(k-1,\\,N-k)} = F_{(2,\\,9)}`} /></div>
          <Def show={D}><b><Tx m="F_0" /> (Ratio de Señal a Ruido):</b> El test de Fisher opera bajo la geometría del error. <Tx m="MST" /> modela el impulso cuadrático de los factores más el ruido basal, en contraste, <Tx m="MSE" /> mide asépticamente el ruido basal intrínseco. Al someterlos a una división de matrices normalizadas, obtenemos el estadístico <Tx m="F_0" />. Cuando operamos bajo <Tx m="H_0" />, el numerador pierde su peso multiplicativo y colapsa térmicamente a 1. En presencia de diferencias, absorbe la tensión inercial y se multiplica desproporcionadamente induciendo la región letal del rechazo.</Def>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 6: Criterio de rechazo</div>
          <div style={S.fm}><TB m={`\\text{Rechazar }H_0\\text{ si }F_0 > F_{${a1},\\,2,\\,9} = ${p1.Fc}`} /></div>
          <Def show={D}><b><Tx m={`F_{\\text{crit}}`} /> (Frontera de Fisher Unilateral):</b> Umbral topológico calculado derivando el logaritmo del cuantil del área asimétrica <Tx m="1-\alpha" /> en la distribución F parametrizada. Operando en varianzas (los cuadrados no tienen soporte negativo), la prueba asume una dinámica estrictamente de cola derecha; el efecto diferencial infla el cociente indefectiblemente hacia arriba.</Def>
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
          <Def show={D}>
            <div style={{ paddingBottom: 6 }}><b>SST (Sustancia Sumada del Tratamiento):</b> Tensor que recoge el desvío externo global (<Tx m={"n\\sum(\\bar{y}_i - \\bar{y}_{..})^2"} />). Es el pivote que cuantifica cuánta matriz inercial logramos absorber simplemente por elegir la técnica <Tx m="i" /> en vez de dejar el sistema al azar.</div>
            <div style={{ paddingBottom: 6 }}><b>SSE (Sustancia Residual Empírica):</b> Vector de norma del fondo cósmico de ruido (<Tx m={"\\sum\\sum(y_{ij} - \\bar{y}_i)^2"} />). Es aquello que nuestra matriz determinista no fue capaz explicativamente de anular; fluctuaciones inexplicables subyacentes.</div>
            <div><b>MS (Transformación Euclidiana Continua a Varianza):</b> Normalización termodinámica de las sumas brutas dividiéndolas entre la extensión del hipercubo ortogonal de libertad métrica (grados de libertad). Transforman un mero acumulable a una esperanza <Tx m="E(MSE) = \sigma^2" />.</div>
          </Def>
          <h4 style={{ margin: "6px 0 3px", color: "#333", fontSize: 13 }}>Sustitución de Fórmulas</h4>
          <div style={S.fm}>
            <TB m={`SST = 4\\big[(${ff(p1.mn[0], 4)}-${ff(p1.gm, 4)})^2 + (${ff(p1.mn[1], 4)}-${ff(p1.gm, 4)})^2 + (${ff(p1.mn[2], 4)}-${ff(p1.gm, 4)})^2\\big] = ${ff(p1.SST, 6)}`} />
            <TB m={`SSE = (${ff(p1.data[0][0], 2)}-${ff(p1.mn[0], 4)})^2 + \\cdots + (${ff(p1.data[2][3], 2)}-${ff(p1.mn[2], 4)})^2 = ${ff(p1.SSE, 6)}`} />
            <TB m={`MST = \\frac{${ff(p1.SST, 6)}}{2} = ${ff(p1.MST, 6)}, \\quad MSE = \\frac{${ff(p1.SSE, 6)}}{9} = ${ff(p1.MSE, 6)}`} />
            <TB m={`F_0 = \\frac{${ff(p1.MST, 6)}}{${ff(p1.MSE, 6)}} = ${f2(p1.F0)}`} />
          </div>
          <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>Tabla ANOVA</h4>
          <table style={S.T}><thead><tr><th style={S.th}>Fuente</th><th style={S.th}>SS</th><th style={S.th}>gl</th><th style={S.th}>MS</th><th style={S.th}><Tx m="F_0" /></th><th style={S.th}><Tx m="F_{\text{crit}}" /></th><th style={S.th}>Valor-p</th></tr></thead>
            <tbody>
              <tr><td style={{ ...S.td, fontWeight: 600 }}>Tratamientos</td><td style={S.td}>{ff(p1.SST, 6)}</td><td style={S.td}>2</td><td style={S.td}>{ff(p1.MST, 6)}</td><td style={{ ...S.td, fontWeight: 700, color: "#c62828" }}>{f2(p1.F0)}</td><td style={S.td}>{p1.Fc}</td><td style={{ ...S.td, fontWeight: 700, color: p1.pval < a1 ? "#2e7d32" : "#c62828" }}>{p1.pval < 0.0001 ? "<0.0001" : ff(p1.pval, 4)}</td></tr>
              <tr style={{ background: "#f9f9f9" }}><td style={{ ...S.td, fontWeight: 600 }}>Error</td><td style={S.td}>{ff(p1.SSE, 6)}</td><td style={S.td}>9</td><td style={S.td}>{ff(p1.MSE, 6)}</td><td style={S.td}>—</td><td style={S.td}>—</td><td style={S.td}>—</td></tr>
              <tr style={{ fontWeight: 600 }}><td style={{ ...S.td, fontWeight: 700 }}>Total</td><td style={S.td}>{ff(p1.SSTot, 6)}</td><td style={S.td}>11</td><td style={S.td}>—</td><td style={S.td}>—</td><td style={S.td}>—</td><td style={S.td}>—</td></tr>
            </tbody></table>
          <div style={{ ...S.fm, marginTop: 10 }}><TB m={`F_0 = ${f2(p1.F0)} \\;${p1.rejected ? ">" : "\\leq"}\\; F_c = ${p1.Fc} \\quad (p = ${p1.pval < 0.0001 ? "<0.0001" : ff(p1.pval, 4)}) \\;\\Longrightarrow\\; \\boxed{\\text{${p1.rejected ? "Se rechaza" : "No se rechaza"} } H_0}`} /></div>
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
        <Def show={D}><b>LSD (Diferencia Mínima Significativa de Fisher):</b> Procedimiento exhaustivo de inferencia post-hoc que desenmascara qué pares de medias difieren. Sólo goza de validez matemática si previamente la prueba F global ha colapsado la hipótesis nula, lo que actúa como escudo probabilístico protector contra la inflación del Error Tipo I. Matemáticamente evalúa pares mediante contrastes <Tx m="\mathbf{t}" /> asumiendo la varianza conjunta general <Tx m="MSE" />, maximizando la potencia analítica a expensas de la tasa de error por familia (FWER).</Def>
        <h4 style={{ margin: "6px 0 3px", color: "#333", fontSize: 13 }}>Fórmula y cálculo</h4>
        <div style={S.fm}><TB m={`LSD = t_{\\alpha/2,\\,N-k}\\cdot\\sqrt{\\frac{2\\cdot MSE}{n}} = ${p1.tc}\\times\\sqrt{\\frac{2\\times${ff(p1.MSE, 6)}}{4}} = ${ff(p1.LSD, 4)}`} /></div>
        <Def show={D}>
          <div style={{ paddingBottom: 6 }}><b><Tx m={`t_{${a1 / 2},\\,N-k}`} /> (Cuantil t de Student):</b> Ancla la asimetría de la distribución univariante del error estandarizado de los contrastes. Se proyecta de forma bilateral estricta (espacio muestral <Tx m="\alpha/2" />) pues carecemos de hipótesis direccionales predefinidas.</div>
          <div><b>Fórmula Específica LSD:</b> La arquitectura <Tx m={"LSD = t_{\\alpha/2,\\,N-k}\\sqrt{2 \\cdot MSE / n}"} /> codifica algebraicamente el error estándar de la diferencia real entre dos variables aleatorias gaussianas de igual cardinalidad. Cualquier brecha muestral superior a esta difracción de onda paramétrica es estadísticamente significativa.</div>
        </Def>
        <h4 style={{ margin: "8px 0 3px", color: "#333", fontSize: 13 }}>Cálculo de diferencias absolutas y Resultados</h4>
        <div style={{ maxHeight: 150, overflowY: "auto", overflowX: "hidden" }}>
          <table style={S.T}><thead><tr><th style={S.th}>Comparación</th><th style={S.th}><Tx m={"|\\bar{y}_i-\\bar{y}_j|"} /> Sustitución</th><th style={S.th}>Módulo</th><th style={S.th}>LSD</th><th style={S.th}>Resultado</th></tr></thead>
            <tbody>{p1.cmp.map((x, i) => <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}><td style={{ ...S.td, fontWeight: 600 }}>{x.a} vs {x.b}</td><td style={S.td}><Tx m={`|${ff(p1.mn[x.i - 1], 4)} - ${ff(p1.mn[x.j - 1], 4)}|`} /></td><td style={{ ...S.td, fontWeight: 600 }}>{ff(x.diff, 4)}</td><td style={S.td}>{ff(p1.LSD, 4)}</td><td style={{ ...S.td, fontWeight: 700, color: x.sig ? "#2e7d32" : "#c62828" }}>{x.sig ? "Significativa" : "No significat."}</td></tr>)}</tbody></table>
        </div>
        <div style={S.ok}>
          <b>Respuesta:</b> SwiftPay tiene el promedio más bajo ({ff(p1.mn[2], 2)} s), pero la diferencia con SwiftFast no es estadísticamente significativa. Por tanto, la evidencia no respalda una migración costosa; SwiftFast ofrece un desempeño similar con menor inversión.
        </div>
      </div>
    </div>),

    // Normalidad
    () => (<div>
      <div style={S.q}>Pregunta 4: ¿Apoyaría el supuesto de normalidad?</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Supuesto de Normalidad (Prueba Shapiro-Wilk)</h3>
        <Def show={D}><b>Test de Shapiro-Wilk (W):</b> Representa el contraste inferencial por excelencia para detectar desviaciones a la Normalidad. El algoritmo proyecta iterativamente el vector de residuos empíricos ordenados sobre el vector esperanza matemática de una curva Gausiana pura, calculando una seudovarianza. El cociente espectral <Tx m="W" /> converge asintóticamente a 1 cuando los datos empatan milimétricamente con fractales gaussianos subyacentes.</Def>
        <div style={S.fm}>
          <TB m={`e_{ij} = y_{ij} - \\bar{y}_i`} />
          <TB m={`W = \\frac{b^2}{SSE} \\quad \\text{donde} \\quad b = \\sum_{i=1}^{n/2} a_i(x_{(n-i+1)} - x_{(i)})`} />
        </div>
        <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>1. Residuos agrupados por Arquitectura</h4>
        <div style={{ display: "flex", gap: "10px", flexWrap: "wrap", marginBottom: "10px" }}>
          {p1.nm.map((name, i) => (
            <div key={i} style={{ flex: 1, minWidth: "200px", background: "#f8f9fa", border: "1px solid #ddd", borderRadius: "6px", padding: "8px" }}>
              <div style={{ fontWeight: 600, color: "#0f3460", marginBottom: "4px", fontSize: "12px", borderBottom: "1px solid #ccc" }}>{name}</div>
              {p1.data[i].map((v, j) => (
                <div key={j} style={{ fontSize: "12px", display: "flex", justifyContent: "space-between" }}>
                  <span><Tx m={`e_{${i + 1}${j + 1}}`} /> {`= ${ff(v, 2)} - ${ff(p1.mn[i], 4)} = `}</span>
                  <span style={{ fontWeight: 600, color: (v - p1.mn[i]) >= 0 ? "#2e7d32" : "#c62828" }}>{ff(v - p1.mn[i], 4)}</span>
                </div>
              ))}
            </div>
          ))}
        </div>
        <Def show={D}><b><Tx m="e_{ij}" /> (Vector de Residuales Brutos):</b> En validación paramétrica sistemática, un residuo va mucho más allá de una simple resta; es el material estocástico putativo puro, nuestra ventana observacional ineludible para juzgar y caracterizar al latente error teórico fundamental <Tx m="\varepsilon_{ij}" />.</Def>
        <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>2. Cálculo del estadístico b (residuos ordenados)</h4>
        <table style={S.T}><thead><tr><th style={S.th}>i</th><th style={S.th}><Tx m="a_i" /></th><th style={S.th}><Tx m="x_{(n-i+1)}" /></th><th style={S.th}><Tx m="x_{(i)}" /></th><th style={S.th}><Tx m="x_{(n-i+1)} - x_{(i)}" /></th><th style={S.th}>Producto</th></tr></thead>
          <tbody>{p1.sw_steps.map((s, i) => <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}>
            <td style={S.td}>{s.i}</td><td style={S.td}>{ff(s.a, 4)}</td><td style={S.td}>{ff(s.xn, 4)}</td><td style={S.td}>{ff(s.x1, 4)}</td><td style={S.td}>{ff(s.diff, 4)}</td><td style={{ ...S.td, fontWeight: 600 }}>{ff(s.prod, 4)}</td>
          </tr>)}
            <tr><td colSpan={5} style={{ ...S.td, textAlign: "right", fontWeight: 700 }}>Suma (b) =</td><td style={{ ...S.td, fontWeight: 700, color: "#c62828" }}>{ff(p1.b_SW, 4)}</td></tr>
            <tr><td colSpan={5} style={{ ...S.td, textAlign: "right", fontWeight: 700 }}>Valor-P =</td><td style={{ ...S.td, fontWeight: 700, color: "#1e3a8a" }}>{ff(p1.pval_SW, 4)}</td></tr>
          </tbody></table>
        <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>3. Prueba de Hipótesis y Resultado</h4>
        <div style={S.fm}>
          <TB m={`W_{\\text{cal}} = \\frac{(${ff(p1.b_SW, 4)})^2}{${ff(p1.SSE, 4)}} = \\frac{${ff(p1.b_SW ** 2, 4)}}{${ff(p1.SSE, 4)}} = ${ff(p1.W_calc, 4)}`} />
          <TB m={`W_{\\text{crit}} (\\alpha=0.05, n=12) = 0.859`} />
          <TB m={`p{\\text{-value}} = ${ff(p1.pval_SW, 4)}`} />
        </div>
        <div style={p1.W_calc >= 0.859 ? S.ok : S.wn}>
          <b>Respuesta:</b> {p1.W_calc >= 0.859
            ? `El estadístico W = ${ff(p1.W_calc, 4)} supera el umbral crítico de 0.859 para N=12. Por tanto, no se rechaza la hipótesis de normalidad: los datos son consistentes con una distribución normal.`
            : `El estadístico W = ${ff(p1.W_calc, 4)} es menor al umbral crítico de 0.859 para N=12. Por tanto, se rechaza la hipótesis de normalidad: existe evidencia estadística de que los errores no provienen de una distribución normal.`}
        </div>
      </div>
    </div>),

    // Aleatoriedad — Estadístico Z de rachas
    () => (<div>
      <div style={S.q}>Pregunta 5: ¿Los datos fueron obtenidos al azar?</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Supuesto de Aleatoriedad — Estadístico Z</h3>
        <Def show={D}>
          <div style={{ paddingBottom: 6 }}><b>Supuesto de Autocovarianza Nula (Independencia):</b> El núcleo del cálculo de gradientes de ANOVA asume invariablemente <Tx m="Cov(\varepsilon_{i}, \varepsilon_{j}) = 0" />. Patrones secuenciales subyacentes distorsionan drásticamente la matriz varianza-covarianza, lo cual, si bien no afecta el sesgo de la recta o medias trazadas artificialmente, anula por completo la robustez algorítmica de los errores estándar, inhabilitando las hipótesis y los hiperparámetros del intervalo.</div>
          <div style={{ paddingBottom: 6 }}><b>Algoritmo de Rachas (Wald-Wolfowitz):</b> Escanéa el flujo cronológico indexado de la serie de tiempo subyacente. Cuenta los cruces ortogonales de los datos sobre el plano cartesiano $y=0$. Las dinámicas auto-correlativas arrastran el sistema induciendo tendencias (pocas rachas extremas) oscilaciones artificiales deterministas o sobrerrepresentación negativa (muchas rachas espásticas).</div>
          <div style={{ paddingBottom: 6 }}><b>Modelado Hipergeométrico Continuo:</b> Bajo la hipótesis nula, el total aleatorio de rachas describe una aproximación de campana con expectativa estructural determinística <Tx m={"E(R) = \\frac{2n_1 n_2}{N} + 1"} /> y dispersión natural <Tx m={"\\sigma_R = \\sqrt{\\frac{2n_1 n_2(2n_1 n_2 - N)}{N^2(N-1)}}"} />.</div>
          <div><b>Criterio Z-Test:</b> Computa la normalización lineal <Tx m={"Z_0 = \\frac{R - E(R)}{\\sigma_R}"} /> para mapear el evento en una curva estándar. Cuantiza geométricamente si el factor temporal destruyó el caos natural.</div>
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
          <h4 style={{ margin: "0 0 6px", color: "#0f3460", fontSize: 13 }}>3. Cálculo del estadístico Z (Sustitución)</h4>
          <div style={S.fm}><TB m={`E(R) = \\frac{2n_1 n_2}{N} + 1 = \\frac{2(${p1.n1})(${p1.n2})}{${p1.N}} + 1 = \\frac{${2 * p1.n1 * p1.n2}}{${p1.N}} + 1 = ${ff(p1.ER, 4)}`} /></div>
          <div style={S.fm}><TB m={`\\sigma_R = \\sqrt{\\frac{2(${p1.n1})(${p1.n2})\\big[2(${p1.n1})(${p1.n2}) - ${p1.N}\\big]}{${p1.N}^2(${p1.N}-1)}} = \\sqrt{\\frac{${2 * p1.n1 * p1.n2}[${2 * p1.n1 * p1.n2} - ${p1.N}]}{${p1.N * p1.N}(${p1.N - 1})}} = ${ff(p1.sdR, 4)}`} /></div>
          <div style={{ ...S.fm, borderLeft: "3px solid #e94560" }}><TB m={`Z_0 = \\frac{R - E(R)}{\\sigma_R} = \\frac{${p1.runs} - ${ff(p1.ER, 4)}}{${ff(p1.sdR, 4)}} = \\frac{${ff(p1.runs - p1.ER, 4)}}{${ff(p1.sdR, 4)}} = ${ff(p1.Z0runs, 4)}`} /></div>
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
        <Def show={D}><b>Homocedasticidad (Isotropía de la Varianza):</b> Piedra angular en la formulación gaussiana-lineal que impone isovarianza estructural: <Tx m="\sigma_1^2 = \sigma_2^2 = \cdots = \sigma_k^2 = \sigma^2" />. Físicamente asume que el mecanismo generador de ruido opera a una misma temperatura (escala) indiferentemente del tratamiento. Su grave perturbación (heterocedasticidad) falsea la métrica combinada de <Tx m="MSE" /> distorsionando irremediablemente la confiabilidad de la frontera de decisión <Tx m="F" />.</Def>
        <div style={S.fm}>
          <TB m={`s_i^2 = \\frac{1}{n-1}\\sum_{j=1}^{n}(y_{ij}-\\bar{y}_i)^2`} />
          <TB m={`s_1^2 = \\frac{(${ff(p1.data[0][0], 2)}-${ff(p1.mn[0], 4)})^2 + \\dots + (${ff(p1.data[0][3], 2)}-${ff(p1.mn[0], 4)})^2}{4-1} = ${ff(p1.vr[0], 6)}`} />
        </div>
        <table style={S.T}><thead><tr><th style={S.th}>Arquitectura</th><th style={S.th}><Tx m="s_i^2" /></th><th style={S.th}><Tx m="s_i" /></th></tr></thead>
          <tbody>{p1.nm.map((nm, i) => <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}><td style={{ ...S.td, fontWeight: 600 }}>{nm}</td><td style={S.td}>{ff(p1.vr[i], 6)}</td><td style={S.td}>{ff(Math.sqrt(p1.vr[i]), 4)}</td></tr>)}</tbody></table>
        {(() => {
          const rat = Math.max(...p1.vr) / Math.min(...p1.vr); return <>
            <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>Prueba F_max de Hartley</h4>
            <div style={S.fm}>
              <TB m={`F_{\\max} = \\frac{s^2_{\\max}}{s^2_{\\min}} = \\frac{${ff(Math.max(...p1.vr), 6)}}{${ff(Math.min(...p1.vr), 6)}} = ${ff(rat, 4)}`} />
              <TB m={`F_{\\text{crit}} (\\alpha=0.05, k=3, v=3) = 15.44`} />
            </div>
            <Def show={D}><b>Test F_max de Hartley:</b> Procedimiento analítico rigoroso para comprobar la homogeneidad al contrastar la razón de las variabilidades extremas empíricas (<Tx m="s^2_{\\max}/s^2_{\\min}" />) contra su respectivo cuantil en la distribución especial de Hartley. Dado que los conjuntos de datos del banco tienen exactamente el mismo tamaño (<Tx m="n=4" />, grados de libertad <Tx m="v=3" /> y <Tx m="k=3" /> arquitecturas), el rechazo formal de la simetría paramétrica ocurre explícitamente si y solo si <Tx m="F_{\\max} \\geq 15.44" />.</Def>
            <div style={rat < 15.44 ? S.ok : S.wn}>
              <b>Respuesta Explícita:</b> La discrepancia extrema entre las varianzas genera un estadístico <Tx m={`F_{\\max} = ${ff(rat, 4)}`} />.
              {rat < 15.44
                ? ` Dado que ${ff(rat, 4)} es estrictamente inferior a 15.44 (el valor crítico F), NO se rechaza la hipótesis nula de homocedasticidad. En el contexto del banco, esto certifica explícitamente que las tres arquitecturas SÍ operan con exactamente el mismo margen paramétrico de volatilidad e inestabilidad algorítmica en sus latencias de respuesta hacia los clientes. `
                : ` Dado que ${ff(rat, 4)} supera al crítico 15.44, SE RECHAZA la hipótesis nula de homogeneidad. En el ecosistema del banco, esto significa inequívocamente que una o más arquitecturas están induciendo picos de inestabilidad desproporcionados en las métricas de respuesta a los clientes. No obstante, dado que el banco planificó exactamente 4 simulaciones de estrés por diseño cerrado (n=4), el error algorítmico queda neutralizado y tu tabla ANOVA sobrevive asimptóticamente intacta.`
              }
            </div>
          </>;
        })()}
      </div>
    </div>),

    // Kruskal-Wallis
    () => (<div>
      <div style={S.q}>Pregunta 7: Y si no se cumplen los supuestos ¿Mantiene las conclusiones dadas en el inciso 2?</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Kruskal-Wallis — Justificación Robusta</h3>
        <Def show={D}><b>Test de Kruskal-Wallis:</b> Cuando los supuestos de normalidad o, de forma muy crítica, el de homocedasticidad para ANOVA se ven gravemente vulnerados (como nuestro <Tx m="F_{\max}" /> lo sugiere), no podemos basar la inferencia plenamente en medias muestrales. Kruskal-Wallis esquiva el error mapeando los datos continuos al dominio discreto de los rangos (posiciones ordinales). Evalúa estadísticamente si los rangos medios difieren, identificando si una procedencia poblacional domina estocásticamente a otra sin asumir varianzas iguales ni normalidad estricta.</Def>

        <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>Paso 1: Asignación de Rangos</h4>
        <div style={{ maxHeight: "200px", overflowY: "auto", border: "1px solid #ddd", borderRadius: 6, marginBottom: 12 }}>
          <table style={{ ...S.T, marginTop: 0 }}>
            <thead style={{ position: "sticky", top: 0, zIndex: 1 }}><tr><th style={S.th}>Observación (<Tx m="y" />)</th><th style={S.th}>Arquitectura</th><th style={S.th}>Rango (<Tx m="r" />)</th></tr></thead>
            <tbody>
              {p1.kw_data.map((d, i) => <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}>
                <td style={S.td}>{ff(d.v, 2)}</td><td style={S.td}>{d.nm}</td><td style={{ ...S.td, fontWeight: 600 }}>{ff(d.r, 1)}</td>
              </tr>)}
            </tbody>
          </table>
        </div>

        <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>Paso 2: Suma de Rangos por Grupo</h4>
        <div style={S.fm}>
          <TB m={`R_i = \\sum r_{ij}`} />
        </div>
        <table style={S.T}>
          <thead><tr><th style={S.th}>Arquitectura</th><th style={S.th}><Tx m="n_i" /></th><th style={S.th}>Suma de Rangos (<Tx m="R_i" />)</th><th style={S.th}><Tx m="R_i^2 / n_i" /></th></tr></thead>
          <tbody>
            {p1.nm.map((nm, i) => <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}>
              <td style={{ ...S.td, fontWeight: 600 }}>{nm}</td>
              <td style={S.td}>{p1.n}</td>
              <td style={{ ...S.td, fontWeight: 700 }}>{ff(p1.kw_R[i], 1)}</td>
              <td style={S.td}>{ff(Math.pow(p1.kw_R[i], 2) / p1.n, 3)}</td>
            </tr>)}
            <tr>
              <td colSpan={3} style={{ ...S.td, textAlign: "right", fontWeight: 700 }}>Sumatoria =</td>
              <td style={{ ...S.td, fontWeight: 700, color: "#c62828" }}>{ff(p1.kw_R.reduce((a, r) => a + Math.pow(r, 2) / p1.n, 0), 3)}</td>
            </tr>
          </tbody>
        </table>

        <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>Paso 3: Cálculo del Estadístico H</h4>
        <div style={S.fm}>
          <TB m={`H = \\frac{12}{N(N+1)} \\sum_{i=1}^k \\frac{R_i^2}{n_i} - 3(N+1)`} />
          <TB m={`H = \\frac{12}{${p1.N}(${p1.N + 1})} (${ff(p1.kw_R.reduce((a, r) => a + Math.pow(r, 2) / p1.n, 0), 3)}) - 3(${p1.N + 1}) = ${ff(p1.kw_H, 4)}`} />
        </div>

        {p1.kw_ties.length > 0 && (
          <>
            <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>Paso 4: Corrección por Empates</h4>
            <Def show={D}><b>Factor de Corrección (<Tx m="C" />):</b> Cuando existen empates (valores iguales), se promedian sus rangos. Esto artificialmente comprime la varianza ordinal total y subestima <Tx m="H" />. Se aplica la corrección penalizadora <Tx m={"C"} /> empleando el número de empates <Tx m="t" /> obtenidos iterativamente; luego se proyecta <Tx m="H_{adj} = H / C" /> recuperando el poder real de la prueba.</Def>
            <div style={S.fm}>
              <TB m={`C = 1 - \\frac{\\sum (t^3 - t)}{N^3 - N} = 1 - \\frac{${p1.kw_tieSum}}{${Math.pow(p1.N, 3)} - ${p1.N}} = ${ff(p1.kw_tieCorrection, 6)}`} />
              <TB m={`H_{adj} = \\frac{H}{C} = \\frac{${ff(p1.kw_H, 4)}}{${ff(p1.kw_tieCorrection, 6)}} = ${ff(p1.kw_H_adj, 4)}`} />
            </div>
          </>
        )}

        <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>Paso 5: Regla de Decisión y Conclusión</h4>
        <div style={S.fm}>
          <TB m={`\\text{Rechazar } H_0 \\text{ si } ${p1.kw_ties.length > 0 ? "H_{adj}" : "H"} > \\chi^2_{\\alpha, k-1}`} />
          <TB m={`\\chi^2_{${p1.alpha}, ${p1.k - 1}} = ${ff(p1.kw_chiCrit, 3)}`} />
        </div>

        <div style={p1.kw_rejected ? S.ok : S.wn}>
          <b>Respuesta:</b> Con un error del {p1.alpha * 100}% y <Tx m={`${p1.k - 1}`} /> grados de libertad, el estadístico superior de Kruskal-Wallis evalúa a {p1.kw_ties.length > 0 ? <Tx m={`H_{adj} = ${ff(p1.kw_H_adj, 4)}`} /> : <Tx m={`H = ${ff(p1.kw_H, 4)}`} />}, lo cual {p1.kw_rejected ? "supera con creces" : "no logra superar"} el límite crítico asintótico condicional de <Tx m={`\\chi^2 = ${ff(p1.kw_chiCrit, 3)}`} />.
          <br /><br />
          {p1.kw_rejected
            ? <>Debido a esta evidencia empírica incontrovertible, <b>se rechaza <Tx m="H_0" /></b> y <b>concluimos inequívocamente que existen diferencias altamente significativas entre las arquitecturas bancarias empleando rangos globales puros</b>. Esta deducción ortogonal libre de distribución <b>confirma de manera robusta y definitiva la conclusión lograda en el tradicional ANOVA de varianzas (inciso 2)</b>, incluso frente al deterioro evidente de la regla de homocedasticidad entre los grupos analizados.</>
            : <>Bajo estas restricciones empíricas, <b>no se rechaza <Tx m="H_0" /></b>. Por ende, desde la lente libre de distribución <b>sostenemos que no existe divergencia operativa significativa entre las arquitecturas</b>. Este nuevo enfoque contradice y absorbe el resultado previo paramétrico, protegiéndote contra falsos positivos motivados por violaciones de homocedasticidad en la matriz muestral.</>
          }
        </div>
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
        <Def show={D}>
          <div style={{ paddingBottom: 6 }}><b>MCO (Mínimos Cuadrados Ordinarios / OLS):</b> Sistema axiomático de optimización que desciende gradientemente sobre una topología cuadrática <Tx m={"\\sum(y_i - \\hat{y}_i)^2"} /> para extraer el vector hiperesférico global <Tx m={"\\hat{\\boldsymbol{\\beta}}"} />. Su forma cerrada, dictada por el álgebra lineal como la proyección ortogonal <Tx m={"(\\mathbf{X}^\\top\\mathbf{X})^{-1}\\mathbf{X}^\\top\\mathbf{Y}"} />, conforma teóricamente bajo Gauss-Markov la familia de estimadores ELIM (Estimación Lineal Insesgada Óptima o BLUE).</div>
          <div><b>Error Estándar Múltiple (<Tx m="SE_i" />):</b> Magnitud de la incerteza asociada escalarmente a una dimensión en <Tx m={"\\hat{\\mathbf{\\beta}}"} /> computada sobre la matriz de varianza-covarianza <Tx m={"[MSE](\\mathbf{X}^\\top\\mathbf{X})^{-1}"} />. Los valores sobre la diagonal de inversión capturan colinealidades. Controla rígidamente la precisión individual del salto de perturbaciones inferenciales.</div>
        </Def>
        <h4 style={{ margin: "8px 0 3px", color: "#333", fontSize: 13 }}>Tabla de símbolos</h4>
        <table style={S.T}><thead><tr><th style={S.thF}>Param.</th><th style={S.thF}>Variable</th><th style={S.thF}>Significado</th></tr></thead>
          <tbody>{bN.map((nm, i) => <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}><td style={{ ...S.td, fontWeight: 600 }}><Tx m={nm} /></td><td style={S.td}>{i === 0 ? "—" : <Tx m={`x_${i}`} />}</td><td style={{ ...S.td, textAlign: "left" }}>{bD[i]}</td></tr>)}</tbody></table>
        <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>Valores estimados</h4>
        <div style={S.fm}>
          <TB m={`\\hat{\\boldsymbol{\\beta}}_i \\text{ y errores estándar se calculan de } (\\mathbf{X}^\\top\\mathbf{X})^{-1} \\text{ (matriz } 5\\times 5).`} />
          <TB m={`SE(\\hat{\\beta}_i) = \\sqrt{MSE \\cdot c_{ii}}`} />
        </div>
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
        <div style={{ ...S.fm, overflowX: "unset" }}>
          <TB m={`\\hat{y}=${ff(p2.b[0], 4)} + (${ff(p2.b[1], 4)})(${xp[0]}) + (${ff(p2.b[2], 4)})(${xp[1]}) + (${ff(p2.b[3], 4)})(${xp[2]}) + (${ff(p2.b[4], 4)})(${xp[3]})`} />
          <TB m={`\\hat{y}=${ff(p2.b[0], 4)} ${p2.b[1] * xp[0] > 0 ? "+" : ""} ${ff(p2.b[1] * xp[0], 4)} ${p2.b[2] * xp[1] > 0 ? "+" : ""} ${ff(p2.b[2] * xp[1], 4)} ${p2.b[3] * xp[2] > 0 ? "+" : ""} ${ff(p2.b[3] * xp[2], 4)} ${p2.b[4] * xp[3] > 0 ? "+" : ""} ${ff(p2.b[4] * xp[3], 4)}`} />
        </div>
        <div style={S.ok}><TB m={`\\boxed{\\hat{y}=${ff(p2.ypr, 4)}\\%}`} /><br /><b>Respuesta:</b> Bajo las condiciones dadas, el modelo predice un uso de CPU del {ff(p2.ypr, 2)}% en los servidores de SwiftPay.</div>
      </div>
    </div>),

    () => (<div>
      <div style={S.q}>Pregunta 3: ANOVA para bondad de ajuste del modelo.</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>ANOVA para RLM — 8 Pasos</h3>
        <Def show={D}><b>Prueba F de Bondad Estructural (ANOVA-Regresión):</b> Contrastación que evalúa la nulidad funcional del hiperplano predictivo en su agrupe total. Propone la desoladora <Tx m={"H_0: \\beta_1 = \\cdots = \\beta_p = 0"} />: argumentando matemáticamente que toda fluctuación en la variable dependiente responde meramente al azar entrópico y la constante <Tx m="\beta_0" />, sin correlatos vectoriales deterministas en el conjunto base. Si rechazada contundentemente, asegura que al menos una proyección dentro de nuestra métrica direccional está enlazada causalmente con la dinámica estudiada.</Def>
        <div style={S.sp}><div style={S.spT}>Pasos 1-4</div>
          <p><Tx m={`n=20,\\;p=4,\\;gl_{reg}=4,\\;gl_{error}=15`} />. <Tx m={`\\alpha=${a2}`} /></p>
          <div style={S.fm}><TB m={"H_0:\\beta_1=\\beta_2=\\beta_3=\\beta_4=0\\qquad H_1:\\text{Al menos un }\\beta_i\\neq 0"} /></div>
        </div>
        <div style={S.sp}><div style={S.spT}>Pasos 5-6</div>
          <div style={S.fm}><TB m={`F_0=\\frac{MSR}{MSE}\\sim F_{(4,15)},\\qquad\\text{Rechazar si }F_0>${p2.Fc}`} /></div>
          <Def show={D}>
            <div style={{ paddingBottom: 6 }}><b>MSR (Varianza Explicada Direccional):</b> Cuantifica vectorialmente cuánto del caos total de Y fue exitosamente modelado por nuestro diseño lineal dividiendo el tensor <Tx m="SSR" /> entre sus dimensiones paramétricas activas <Tx m="p" />.</div>
            <div><b>MSE (Varianza Residual Empírica):</b> Fracción de ignorancia del modelo sobre los grados de libertad libres espaciales. Si el modelo es matemáticamente válido, el ratio <Tx m="MSR/MSE" /> estallará como indicativo irrefutable de señal.</div>
          </Def>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 7: Sustituciones y Tabla ANOVA</div>
          <div style={S.fm}>
            <TB m={`SSR = \\sum(\\hat{y}_i - \\bar{y})^2 = \\cdots = ${ff(p2.SSR, 4)}`} />
            <TB m={`SSE = \\sum(y_i - \\hat{y}_i)^2 = (${ff(p2.raw[0][0], 2)} - ${ff(p2.raw[0][0] - (p2.raw[0][0] - p2.b[0]), 4)})^2 + \\cdots = ${ff(p2.SSE, 4)}`} />
            <TB m={`MSR = \\frac{SSR}{p} = \\frac{${ff(p2.SSR, 4)}}{4} = ${ff(p2.MSR, 4)}`} />
            <TB m={`MSE = \\frac{SSE}{n-p-1} = \\frac{${ff(p2.SSE, 4)}}{15} = ${ff(p2.MSE, 4)}`} />
            <TB m={`F_0 = \\frac{MSR}{MSE} = \\frac{${ff(p2.MSR, 4)}}{${ff(p2.MSE, 4)}} = ${f2(p2.Fr)}`} />
          </div>
          <h4 style={{ margin: "8px 0 3px", color: "#333", fontSize: 13 }}>Valores de la Tabla ANOVA</h4>
          <table style={S.T}><thead><tr><th style={S.th}>Fuente</th><th style={S.th}>SS</th><th style={S.th}>gl</th><th style={S.th}>MS</th><th style={S.th}><Tx m="F_0" /></th><th style={S.th}>Valor-p</th></tr></thead>
            <tbody>
              <tr><td style={{ ...S.td, fontWeight: 600 }}>Regresión</td><td style={S.td}>{ff(p2.SSR, 4)}</td><td style={S.td}>4</td><td style={S.td}>{ff(p2.MSR, 4)}</td><td style={{ ...S.td, fontWeight: 700, color: "#c62828" }}>{f2(p2.Fr)}</td><td style={{ ...S.td, fontWeight: 700, color: p2.pvalGlobal < a2 ? "#2e7d32" : "#c62828" }}>{p2.pvalGlobal < 0.0001 ? "<0.0001" : ff(p2.pvalGlobal, 4)}</td></tr>
              <tr style={{ background: "#f9f9f9" }}><td style={{ ...S.td, fontWeight: 600 }}>Error</td><td style={S.td}>{ff(p2.SSE, 4)}</td><td style={S.td}>15</td><td style={S.td}>{ff(p2.MSE, 4)}</td><td style={S.td}>—</td><td style={S.td}>—</td></tr>
              <tr style={{ fontWeight: 600 }}><td style={{ ...S.td, fontWeight: 700 }}>Total</td><td style={S.td}>{ff(p2.SST, 4)}</td><td style={S.td}>19</td><td style={S.td}>—</td><td style={S.td}>—</td><td style={S.td}>—</td></tr>
            </tbody></table>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 8: Conclusión</div>
          <div style={S.ok}><Tx m={`F_0=${f2(p2.Fr)}\\;${p2.rejectedGlobal ? ">" : "\\leq"}\\;${p2.Fc} \\quad (p = ${p2.pvalGlobal < 0.0001 ? "<0.0001" : ff(p2.pvalGlobal, 4)})`} />.{" "}
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
        <Def show={D}>
          <div style={{ paddingBottom: 6 }}><b><Tx m="R^2" /> (Coeficiente de Determinación Global):</b> Índice escalar <Tx m={"R^2 = SSR/SST \\in [0,1]"} /> que encarna el porcentaje algorítmico estricto de la varianza total de <Tx m="Y" /> que es abducido por la red combinada de covariables <Tx m="X" />. Siendo una métrica inercial engañosa, el algoritmo de mínimos cuadrados garantiza que <Tx m="R^2" /> subirá monótonamente al inyectarle cualquier variable extra dimensional, sea ruido o señal pura.</div>
          <div><b><Tx m={"R^2_{\\text{adj}}"} /> (Pseudo-Coeficiente Penalizado):</b> Transformación entrópica correctiva de Theil-Wherry. Pondera el aumento del <Tx m="R^2" /> exigiendo un "peaje" por inflación dimensional de parámetros <Tx m="p" />. Si una variable no inyecta correlación que amortigüe la pérdida de un grado de libertad computacional, su inclusión destruye poder y el <Tx m={"R^2_{\\text{adj}}"} /> decrece implacablemente.</div>
        </Def>
        <h4 style={{ margin: "8px 0 3px", color: "#333", fontSize: 13 }}>Cálculo y Sustitución</h4>
        <div style={S.fm}>
          <TB m={`R^2 = \\frac{${ff(p2.SSR, 4)}}{${ff(p2.SST, 4)}} = ${ff(p2.R2, 4)}`} />
          <TB m={`R^2_{\\text{adj}} = 1 - \\frac{(1-${ff(p2.R2, 4)})(${p2.no}-1)}{${p2.no}-${p2.p}-1} = 1 - \\frac{(${ff(1 - p2.R2, 4)})(${p2.no - 1})}{${p2.no - p2.p - 1}} = ${ff(p2.R2a, 4)}`} />
        </div>
        <div style={S.ok}><b>Respuesta:</b> El modelo explica el <b>{ff(p2.R2 * 100, 2)}%</b> de la variabilidad del uso de CPU. <Tx m={`R^2_{\\text{adj}}=${ff(p2.R2a * 100, 2)}\\%`} />. El {ff((1 - p2.R2) * 100, 2)}% restante corresponde a factores no incluidos o variabilidad aleatoria.</div>
      </div>
    </div>),

    () => (<div>
      <div style={S.q}>Pregunta 5: IC del 95% para cada <Tx m="\beta_i" />.</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Intervalos de Confianza</h3>
        <div style={S.fm}><TB m={`IC_{${(1 - a2) * 100}\\%}:\\;\\hat{\\beta}_i\\pm t_{${a2 / 2},\\,15}\\cdot SE(\\hat{\\beta}_i),\\quad t_{${a2 / 2},15}=${p2.tc}`} /></div>
        <Def show={D}><b>IC (Intervalo de Confianza Bidireccional):</b> Construcción topológica probabilística <Tx m={"\\hat{\\beta}_i \\pm t_{\\alpha/2} \\cdot SE_i"} />. No nos dice qué probabilidad hay de que contenga el verdadero <Tx m="\beta_i" />, sino que nos asegura matemáticamente que si reptiéramos la regresión infinitas veces desde el proceso generador, el {(1 - a2) * 100}% de estos hipercubos construidos iterativamente aprisionarían la constante subyacente determinista <Tx m="\beta_i" /> real. Si <Tx m={"0 \\in IC"} />, la topología intercepta el vacío, por lo que estocásticamente el efecto de la covariable en simulación es indistinguible de ruido blanco univariante.</Def>
        <table style={S.T}><thead><tr><th style={S.th}>Parám.</th><th style={S.th}><Tx m={"\\hat{\\beta}_i"} /></th><th style={S.th}>Sustitución (<Tx m="\pm" />)</th><th style={S.th}>Inferior</th><th style={S.th}>Superior</th><th style={S.th}><Tx m={"0\\in IC"} /></th></tr></thead>
          <tbody>{bN.map((nm, i) => { const c0 = p2.cl[i] <= 0 && p2.ch[i] >= 0; return <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}><td style={{ ...S.td, fontWeight: 600 }}><Tx m={nm} /></td><td style={S.td}>{ff(p2.b[i], 4)}</td><td style={S.td}><Tx m={`${ff(p2.b[i], 4)} \\pm (${p2.tc})(${ff(p2.se[i], 4)})`} /></td><td style={S.td}>{ff(p2.cl[i], 4)}</td><td style={S.td}>{ff(p2.ch[i], 4)}</td><td style={{ ...S.td, fontWeight: 700, color: c0 ? "#c62828" : "#2e7d32" }}>{c0 ? "Sí" : "No"}</td></tr>; })}</tbody></table>
        <div style={S.ok}><b>Respuesta:</b> Los parámetros cuyo IC no contiene 0 son significativos al {(1 - a2) * 100}% de confianza para el modelo de CPU de SwiftPay.</div>
      </div>
    </div>),

    () => (<div>
      <div style={S.q}>Pregunta 6: ¿El aporte de cada <Tx m="X_i" /> es significativo (5%)?</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Significancia Individual — 8 Pasos</h3>
        <Def show={D}><b>Test Marginal Condicional (Significancia Individual):</b> Mecanismo de poda local ortogonal. Valida la proyección marginal que un regresor aporta una vez el hiperplano ya absorbió a todas las demás variables. Si <Tx m={"H_0: \\beta_i = 0"} /> no rechaza, significa rotunda y algebraicamente que la variable provee información idéntica, redundante y predecible (colinealidad mutua) frente a la envolvente que ya describían las demás.</Def>
        <div style={S.sp}><div style={S.spT}>Pasos 1-6</div>
          <div style={S.fm}><TB m={`H_0:\\beta_i=0\\quad H_1:\\beta_i\\neq 0\\quad\\alpha=${a2}`} /><TB m={`t_0=\\frac{\\hat{\\beta}_i}{SE(\\hat{\\beta}_i)},\\quad\\text{Rechazar si }|t_0|>t_{${a2 / 2},15}=${p2.tc}`} /></div>
          <Def show={D}><b><Tx m="t_0" /> (Señal al Ruido Individual):</b> Operador geométrico <Tx m={"\\hat{\\beta}_i / SE_i"} />. Descompone a fondo la inclinación de la derivada local purificándola del error circunspecto a ese eje dimensional único.</Def>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 7: Cálculos con Sustitución</div>
          <table style={S.T}><thead><tr><th style={S.th}>Variable</th><th style={S.th}><Tx m={"t_0 = \\frac{\\hat{\\beta}_i}{SE}"} /></th><th style={S.th}><Tx m="t_0" /> res</th><th style={S.th}><Tx m={`|t_0|\\text{ vs }${p2.tc}`} /></th><th style={S.th}>Valor-p</th><th style={S.th}>Decisión</th></tr></thead>
            <tbody>{bN.slice(1).map((nm, i) => { const ix = i + 1, sg = Math.abs(p2.tv[ix]) > p2.tc; return <tr key={i} style={{ background: i % 2 ? "#f9f9f9" : "#fff" }}><td style={{ ...S.td, fontWeight: 600 }}><Tx m={nm} /> ({bD[ix]})</td><td style={S.td}><Tx m={`\\frac{${ff(p2.b[ix], 4)}}{${ff(p2.se[ix], 4)}}`} /></td><td style={{ ...S.td, fontWeight: 700 }}>{ff(p2.tv[ix], 4)}</td><td style={S.td}>{ff(Math.abs(p2.tv[ix]), 4)} {sg ? ">" : "≤"} {p2.tc}</td><td style={{ ...S.td, fontWeight: 700, color: p2.pvals[ix] < a2 ? "#2e7d32" : "#c62828" }}>{p2.pvals[ix] < 0.0001 ? "<0.0001" : ff(p2.pvals[ix], 4)}</td><td style={{ ...S.td, fontWeight: 700, color: sg ? "#2e7d32" : "#c62828" }}>{sg ? "Signif." : "No signif."}</td></tr>; })}</tbody></table>
        </div>
        <div style={S.sp}><div style={S.spT}>Paso 8: Conclusión</div>
          <div style={S.ok}><b>Respuesta:</b> Con <Tx m={`\\alpha=${a2}`} />, las variables significativas son: <b>{p2.sigVars.length > 0 ? p2.sigVars.map(i => bD[i]).join(", ") : "ninguna"}</b>. Las no significativas ({p2.nsVars.map(i => bD[i]).join(", ")}) podrían eliminarse del modelo de telemetría de SwiftPay.</div>
          <Def show={D}><b>Principio Heurístico de Parsimonia (Navaja de Ockham):</b> Criterio cibernético supremo para la inferencia de modelos, forzando la maximización explícita de predictibilidad con la mínima cardinalidad de arquitectura. Un modelo sobredimensionado memorizará ruido (Overfitting); un modelo parsimonioso extraerá leyes latentes puras.</Def>
        </div>
      </div>
    </div>),

    () => (<div>
      <div style={S.q}>Pregunta 7: ¿Cuál modelo recomendaría? Justifique.</div>
      <div style={S.cd}>
        <h3 style={{ color: "#0f3460", marginTop: 0, fontSize: 15 }}>Modelo Recomendado</h3>
        <table style={S.T}><thead><tr><th style={S.th}>Variable</th><th style={S.th}><Tx m="|t_0|" /></th><th style={S.th}>vs <Tx m={`${p2.tc}`} /></th><th style={S.th}>Valor-p</th><th style={S.th}>Decisión</th></tr></thead>
          <tbody>{[1, 2, 3, 4].map(i => { const sg = Math.abs(p2.tv[i]) > p2.tc; return <tr key={i} style={{ background: i % 2 ? "" : "#f9f9f9" }}><td style={{ ...S.td, fontWeight: 600 }}><Tx m={`x_${i}`} /> ({bD[i]})</td><td style={S.td}>{ff(Math.abs(p2.tv[i]), 4)}</td><td style={S.td}>{sg ? ">" : "≤"} {p2.tc}</td><td style={{ ...S.td, fontWeight: 700, color: p2.pvals[i] < a2 ? "#2e7d32" : "#c62828" }}>{p2.pvals[i] < 0.0001 ? "<0.0001" : ff(p2.pvals[i], 4)}</td><td style={{ ...S.td, fontWeight: 700, color: sg ? "#2e7d32" : "#c62828" }}>{sg ? "Incluir" : "Eliminar"}</td></tr>; })}</tbody></table>
        <h4 style={{ margin: "10px 0 3px", color: "#333", fontSize: 13 }}>Modelo reducido propuesto</h4>
        <div style={S.fm}><TB m={`\\hat{y}=\\beta_0${p2.sigVars.map(i => `+\\beta_${i}x_${i}`).join("")}`} /></div>
        <Def show={D}><b>Reducción Activa (Subset Backward Selection):</b> Metodología de mutilación paramétrica racional para converger a un modelo anidado óptimo. Exterminar coeficientes espurios condensa la matriz covarianza <Tx m={"[(X^\\top X)^{-1}]"} />, estrechando abismalmente los intervalos de confianza en los coeficientes élites que sobreviven.</Def>
        <div style={S.ok}>
          <b>Respuesta Concisa:</b> {p2.nsVars.length > 0
            ? <>Se recomienda implementar un <b>modelo de regresión reducido</b> utilizando estrictamente las variables significativas: <b>{p2.sigVars.map(i => `${bD[i]} (x_${i})`).join(", ")}</b>. Las métricas restantes ({p2.nsVars.map(i => bD[i]).join(", ")}) se descartan definitivamente del sistema por no aportar poder predictivo estadísticamente válido.</>
            : <>Se recomienda mantener el <b>modelo de regresión completo</b> actual. Todas y cada una de las variables predictoras resultaron estadísticamente significativas, demostrando un armado perfecto sin sobreajuste ni ruido paramétrico.</>
          }
        </div>
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
        <button style={S.mt(tab === "p2")} onClick={() => setTab("p2")}>PARTE II</button>
      </div>
      <div style={{ background: "#fff", border: "1px solid #ddd", borderTop: "none", borderRadius: "0 0 8px 8px", padding: 12 }}>
        {edit && tab === "p1" && <DataEditor1 initialD={d1} initialA={a1} initialSig={sig} onApply={handleApply1} onCancel={() => setEdit(false)} />}
        {edit && tab === "p2" && <DataEditor2 initialD={d2} initialA={a2} initialXp={xp} onApply={handleApply2} onCancel={() => setEdit(false)} />}

        {!edit && tab === "p1" && (<><div style={S.sb}>{t1.map((t, i) => <button key={i} style={S.st(s1 === i)} onClick={() => setS1(i)}>{t}</button>)}</div>{P1[s1]()}</>)}
        {!edit && tab === "p2" && (<><div style={S.sb}>{t2.map((t, i) => <button key={i} style={S.st(s2 === i)} onClick={() => setS2(i)}>{t}</button>)}</div>{P2[s2]()}</>)}
      </div>
      <div style={{ textAlign: "center", padding: "10px 0", color: "#999", fontSize: 10 }}>Estadística Inferencial — FACYT 2026</div>
      <div style={S.fab}>
        <button style={S.fabBtn(edit ? "#f59e0b" : "#92400e")} onClick={() => setEdit(e => !e)}>{edit ? "✕ Cerrar editor" : "⚙️ Editar datos"}</button>
        <button style={S.fabBtn(defs ? "#6366f1" : "#1e1b4b")} onClick={() => setDefs(d => !d)}>{defs ? "✕ Ocultar defs" : "📖 Definiciones"}</button>
      </div>
    </div>
  );
}
