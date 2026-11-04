// splitPathData（`d` のサブパス分割）のテスト。実行: `npm test`
//
// 2 つ目以降のサブパスが相対 `m` で始まると、切り出した単独の `d` では始点が原点基準に
// 化ける（`d` の最初の `m` は絶対として解釈されるため）。splitPathData は現在点を追って
// 始点を絶対 `M x y` に直す。ここでは
//   1. 期待する文字列そのものを確かめる
//   2. splitPathData とは独立に書いた参照の解釈器で「元の `d` の中でのサブパス」と
//      「切り出した単独の `d`」の頂点列が一致することを、ランダムなパスで確かめる
//   3. （任意）ブラウザ自身の解釈と突き合わせる: DERLE_BROWSER_TEST=1 npm test
//      （Playwright と chromium が要る。無ければ飛ばす）
"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const L = require("../logic.js");

// --- 1. 文字列 -------------------------------------------------------------------

test("splitPathData: 絶対 M で始まるサブパスは 1 文字も変えない", () => {
  assert.deepEqual(L.splitPathData("M0 0L1 1M2 2L3 3"), ["M0 0L1 1", "M2 2L3 3"]);
  assert.deepEqual(L.splitPathData("M0,0c1,1,2,2,3,3Z M10,10h5"), ["M0,0c1,1,2,2,3,3Z", "M10,10h5"]);
  assert.deepEqual(L.splitPathData(""), []);
  assert.deepEqual(L.splitPathData(null), []);
  assert.deepEqual(L.splitPathData("  M0 0L1 1  "), ["M0 0L1 1"]);
});

test("splitPathData: 最初のサブパスの m はそのまま（単独でも絶対として解釈されるので同じ）", () => {
  assert.deepEqual(L.splitPathData("m1 2 3 4"), ["m1 2 3 4"]);
});

test("splitPathData: 2 つ目以降の相対 m は、直前の現在点から絶対 M に直す", () => {
  assert.deepEqual(L.splitPathData("M0 0L1 1m2 2l1 0"), ["M0 0L1 1", "M3 3l1 0"]);
  assert.deepEqual(L.splitPathData("m1 2 3 4m1 1"), ["m1 2 3 4", "M5 7"]);
  assert.deepEqual(L.splitPathData("M0,0c1,1,2,2,3,3m-1-1"), ["M0,0c1,1,2,2,3,3", "M2 2"]);
});

test("splitPathData: z の後の現在点はそのサブパスの始点", () => {
  assert.deepEqual(L.splitPathData("M10 10h5v5zm1 1h1"), ["M10 10h5v5z", "M11 11h1"]);
  assert.deepEqual(L.splitPathData("M10 10l5 0 0 5Zm1 1h1"), ["M10 10l5 0 0 5Z", "M11 11h1"]);
});

test("splitPathData: m に続く座標対は暗黙の相対 l として相対のまま残す", () => {
  assert.deepEqual(L.splitPathData("M10 10h5v5zm1 1 2 2 3 3"), ["M10 10h5v5z", "M11 11l2 2 3 3"]);
  assert.deepEqual(L.splitPathData("M0 0m1,1,2,2"), ["M0 0", "M1 1l2,2"]);
  // 暗黙の l の後にコマンドが続く
  assert.deepEqual(L.splitPathData("M0 0m1 1 2 2h3"), ["M0 0", "M1 1l2 2h3"]);
});

test("splitPathData: 相対 m が連続しても現在点を積み上げる", () => {
  assert.deepEqual(L.splitPathData("M0 0m1 1m1 1m1 1"), ["M0 0", "M1 1", "M2 2", "M3 3"]);
  assert.deepEqual(L.splitPathData("M0 0m1 1 1 1m1 1"), ["M0 0", "M1 1l1 1", "M3 3"]);
});

test("splitPathData: 各コマンドの終点（H/V/C/S/Q/T/A、相対・絶対、暗黙の繰り返し）", () => {
  assert.deepEqual(L.splitPathData("M0 0H5m1 1"), ["M0 0H5", "M6 1"]);
  assert.deepEqual(L.splitPathData("M0 0h5 5m1 1"), ["M0 0h5 5", "M11 1"]);
  assert.deepEqual(L.splitPathData("M0 0V5m1 1"), ["M0 0V5", "M1 6"]);
  assert.deepEqual(L.splitPathData("M0 0v5 5m1 1"), ["M0 0v5 5", "M1 11"]);
  assert.deepEqual(L.splitPathData("M0 0C1 1 2 2 30 40m1 1"), ["M0 0C1 1 2 2 30 40", "M31 41"]);
  assert.deepEqual(L.splitPathData("M1 1c1 1 2 2 3 4 1 1 2 2 3 4m1 1"), ["M1 1c1 1 2 2 3 4 1 1 2 2 3 4", "M8 10"]);
  assert.deepEqual(L.splitPathData("M0 0S1 1 7 8m1 1"), ["M0 0S1 1 7 8", "M8 9"]);
  assert.deepEqual(L.splitPathData("M1 1s1 1 7 8m1 1"), ["M1 1s1 1 7 8", "M9 10"]);
  assert.deepEqual(L.splitPathData("M0 0Q1 1 7 8m1 1"), ["M0 0Q1 1 7 8", "M8 9"]);
  assert.deepEqual(L.splitPathData("M1 1q1 1 7 8m1 1"), ["M1 1q1 1 7 8", "M9 10"]);
  assert.deepEqual(L.splitPathData("M0 0T7 8m1 1"), ["M0 0T7 8", "M8 9"]);
  assert.deepEqual(L.splitPathData("M1 1t7 8 1 1m1 1"), ["M1 1t7 8 1 1", "M10 11"]);
  assert.deepEqual(L.splitPathData("M0 0A5 5 0 0 1 10 0m1 1"), ["M0 0A5 5 0 0 1 10 0", "M11 1"]);
  assert.deepEqual(L.splitPathData("M1 1a5 5 0 0 1 10 0m1 1"), ["M1 1a5 5 0 0 1 10 0", "M12 2"]);
});

test("splitPathData: 円弧のフラグは区切りなしで詰めて書ける（\"1010\" は 1, 0, 10）", () => {
  assert.deepEqual(L.splitPathData("M0 0a5 5 0 1010 10m1 1"), ["M0 0a5 5 0 1010 10", "M11 11"]);
  assert.deepEqual(L.splitPathData("M0 0a5,5,30,1,0,10,10m1 1"), ["M0 0a5,5,30,1,0,10,10", "M11 11"]);
  assert.deepEqual(L.splitPathData("M0 0a5 5 0 1 1-10-10m1 1"), ["M0 0a5 5 0 1 1-10-10", "M-9 -9"]);
});

test("splitPathData: 数値の書き方（\"1.5.5\" / \"-1-2\" / 指数 / 改行やタブ）", () => {
  assert.deepEqual(L.splitPathData("M1.5.5L-1-2m.5.5"), ["M1.5.5L-1-2", "M-0.5 -1.5"]);
  assert.deepEqual(L.splitPathData("M1e1 2E-1l+3,-.5m1e1-1e0"), ["M1e1 2E-1l+3,-.5", "M23 -1.3"]);
  assert.deepEqual(L.splitPathData("M0 0 L 10 0 z\n\tm 5,5 l 1 1"), ["M0 0 L 10 0 z", "M5 5l 1 1"]);
});

test("splitPathData: 壊れたパスデータは旧実装と同じ単純な分割に戻す", () => {
  // 引数が足りない
  assert.deepEqual(L.splitPathData("M0 0L1m2 2"), ["M0 0L1", "m2 2"]);
  // M より前に何かある
  assert.deepEqual(L.splitPathData("L1 1m2 2"), ["L1 1", "m2 2"]);
  // 未知の文字
  assert.deepEqual(L.splitPathData("M0 0X1 1m2 2"), ["M0 0X1 1", "m2 2"]);
});

// --- 2. 参照の解釈器との突き合わせ ---------------------------------------------------
//
// splitPathData とは別に書いた、`d` 全体を解釈してサブパスごとの頂点（各コマンドの終点）を
// 返す最小の解釈器。曲線は終点だけを見る（splitPathData が直すのは始点の座標だけなので、
// 頂点列が一致すれば同じ形になる。曲線の制御点は相対のまま残り、始点が合えば一致する）。

function referenceSubpaths(d) {
  const tokens = [];
  const re = /([MmLlHhVvCcSsQqTtAaZz])|([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)|([\s,]+)/g;
  let m;
  let last = 0;
  while ((m = re.exec(d))) {
    if (m.index !== last) throw new Error(`解釈できない文字: ${d.slice(last, m.index)}`);
    last = re.lastIndex;
    if (m[1]) tokens.push(m[1]);
    else if (m[2]) tokens.push(m[2]);
  }
  if (last !== d.length) throw new Error("解釈できない末尾");
  // 円弧のフラグが詰めて書かれていると数値トークンが合体するので、A/a の引数だけ文字単位で読み直す
  const out = [];
  let i = 0;
  let cmd = null;
  let x = 0;
  let y = 0;
  let sx = 0;
  let sy = 0;
  let first = true;
  let flagBuffer = "";
  const nextNumber = () => Number(tokens[i++]);
  const nextFlag = () => {
    if (!flagBuffer) flagBuffer = tokens[i++];
    const c = flagBuffer[0];
    flagBuffer = flagBuffer.slice(1);
    if (c !== "0" && c !== "1") throw new Error("円弧のフラグ");
    return Number(c);
  };
  const nextNumberAfterFlag = () => {
    if (flagBuffer) {
      const v = Number(flagBuffer);
      flagBuffer = "";
      return v;
    }
    return nextNumber();
  };
  const isCmd = (t) => /^[A-Za-z]$/.test(t);
  while (i < tokens.length) {
    if (isCmd(tokens[i])) cmd = tokens[i++];
    const rel = cmd === cmd.toLowerCase();
    switch (cmd.toUpperCase()) {
      case "M": {
        const px = nextNumber();
        const py = nextNumber();
        x = (rel && !first ? x : 0) + px;
        y = (rel && !first ? y : 0) + py;
        first = false;
        sx = x;
        sy = y;
        out.push([[x, y]]);
        cmd = rel ? "l" : "L";
        break;
      }
      case "L":
      case "T": {
        const px = nextNumber();
        const py = nextNumber();
        x = (rel ? x : 0) + px;
        y = (rel ? y : 0) + py;
        out[out.length - 1].push([x, y]);
        break;
      }
      case "H":
        x = (rel ? x : 0) + nextNumber();
        out[out.length - 1].push([x, y]);
        break;
      case "V":
        y = (rel ? y : 0) + nextNumber();
        out[out.length - 1].push([x, y]);
        break;
      case "C":
      case "S":
      case "Q": {
        const n = cmd.toUpperCase() === "C" ? 6 : 4;
        const v = [];
        for (let k = 0; k < n; k += 1) v.push(nextNumber());
        x = (rel ? x : 0) + v[n - 2];
        y = (rel ? y : 0) + v[n - 1];
        out[out.length - 1].push([x, y]);
        break;
      }
      case "A": {
        nextNumber();
        nextNumber();
        nextNumber();
        nextFlag();
        nextFlag();
        const px = nextNumberAfterFlag();
        const py = nextNumber();
        x = (rel ? x : 0) + px;
        y = (rel ? y : 0) + py;
        out[out.length - 1].push([x, y]);
        break;
      }
      case "Z":
        x = sx;
        y = sy;
        out[out.length - 1].push([x, y]);
        // 次が数値なら不正（ここでは生成しない）
        break;
      default:
        throw new Error(`未知のコマンド ${cmd}`);
    }
  }
  return out;
}

/** 決定的な疑似乱数（再現できるように） */
function makeRng(seed) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

/** ランダムなパスを作る。数値の書き方・区切り・相対/絶対・暗黙の繰り返しを混ぜる */
function randomPath(rng) {
  const pick = (arr) => arr[Math.floor(rng() * arr.length)];
  const num = () => {
    const v = (rng() - 0.5) * 200;
    const kind = rng();
    if (kind < 0.4) return String(Math.round(v));
    if (kind < 0.8) return String(Math.round(v * 100) / 100);
    if (kind < 0.9) return (Math.round(v * 10) / 10).toExponential(2);
    return String(Math.round(v * 1000) / 1000).replace(/^(-?)0\./, "$1.");
  };
  let text = "";
  const emitNumbers = (values) => {
    values.forEach((v, k) => {
      if (k === 0) {
        text += pick(["", " "]);
      } else if (v.startsWith("-") || (v.startsWith(".") && /\.\d*$/.test(values[k - 1]) && values[k - 1].includes("."))) {
        text += pick(["", " ", ","]);
      } else {
        text += pick([" ", ","]);
      }
      text += v;
    });
  };
  const subpaths = 2 + Math.floor(rng() * 4);
  for (let sp = 0; sp < subpaths; sp += 1) {
    const mv = sp === 0 ? pick(["M", "m"]) : pick(["m", "m", "M"]);
    text += (sp > 0 ? pick(["", " ", "\n"]) : "") + mv;
    const pairs = 1 + (rng() < 0.3 ? 1 + Math.floor(rng() * 2) : 0);
    const vals = [];
    for (let k = 0; k < pairs * 2; k += 1) vals.push(num());
    emitNumbers(vals);
    const cmds = 1 + Math.floor(rng() * 6);
    for (let c = 0; c < cmds; c += 1) {
      const cmd = pick(["L", "l", "H", "h", "V", "v", "C", "c", "S", "s", "Q", "q", "T", "t", "A", "a"]);
      text += pick(["", " "]) + cmd;
      const reps = rng() < 0.25 ? 2 : 1;
      for (let r = 0; r < reps; r += 1) {
        const upper = cmd.toUpperCase();
        if (upper === "A") {
          const rx = String(1 + Math.round(rng() * 50));
          const ry = String(1 + Math.round(rng() * 50));
          const rot = String(Math.round(rng() * 90));
          const f1 = pick(["0", "1"]);
          const f2 = pick(["0", "1"]);
          const ex = num();
          const ey = num();
          text += (r > 0 ? " " : pick(["", " "])) + `${rx} ${ry} ${rot}`;
          // フラグは詰めて書くこともある
          text += rng() < 0.5 ? ` ${f1}${f2}${ex.startsWith("-") ? "" : " "}${ex}` : ` ${f1} ${f2} ${ex}`;
          text += (ey.startsWith("-") ? "" : " ") + ey;
        } else {
          const n = { L: 2, T: 2, H: 1, V: 1, C: 6, S: 4, Q: 4 }[upper];
          const vals2 = [];
          for (let k = 0; k < n; k += 1) vals2.push(num());
          if (r > 0) text += " ";
          emitNumbers(vals2);
        }
      }
    }
    if (rng() < 0.5) text += pick(["z", "Z", " z"]);
  }
  return text;
}

function assertSameVertices(a, b, d, k) {
  assert.equal(a.length, b.length, `頂点数 (サブパス ${k}): ${d}`);
  for (let i = 0; i < a.length; i += 1) {
    const tol = 1e-9 * Math.max(1, Math.abs(a[i][0]), Math.abs(a[i][1]));
    assert.ok(
      Math.abs(a[i][0] - b[i][0]) <= tol && Math.abs(a[i][1] - b[i][1]) <= tol,
      `サブパス ${k} の頂点 ${i}: ${a[i]} != ${b[i]}\n  d = ${d}`,
    );
  }
}

const FUZZ_PATHS = (() => {
  const rng = makeRng(20261104);
  return Array.from({ length: 2000 }, () => randomPath(rng));
})();

test("参照の解釈器: 自己点検（既知の値）", () => {
  assert.deepEqual(referenceSubpaths("M0 0l1 1m2 2 1 0z"), [[[0, 0], [1, 1]], [[3, 3], [4, 3], [3, 3]]]);
  assert.deepEqual(referenceSubpaths("m1 1a5 5 0 1110 10"), [[[1, 1], [11, 11]]]);
});

test("splitPathData: ランダムなパス 2000 本で、切り出したサブパスの頂点が元の `d` の中と一致する", () => {
  let relativeStarts = 0;
  for (const d of FUZZ_PATHS) {
    const whole = referenceSubpaths(d);
    const parts = L.splitPathData(d);
    assert.equal(parts.length, whole.length, `サブパス数: ${d}`);
    parts.forEach((part, k) => {
      if (k > 0 && /^\s*m/.test(d.split(/(?=[Mm])/)[k] || "")) relativeStarts += 1;
      const alone = referenceSubpaths(part);
      assert.equal(alone.length, 1, `切り出した部分が 1 サブパスでない: ${part}`);
      assertSameVertices(alone[0], whole[k], d, k);
    });
  }
  assert.ok(relativeStarts > 1000, `相対 m の始点を十分に含む (${relativeStarts})`);
});

test("splitPathData → pathDataToPoints: 点列の始点と終点が、元の `d` の中のサブパスと一致する", () => {
  // 円弧は pathDataToPoints が未対応（null）なので、円弧を含まないものだけ比べる
  let compared = 0;
  for (const d of FUZZ_PATHS) {
    if (/[Aa]/.test(d)) continue;
    const parts = L.splitPathData(d);
    // 点列化した始点・終点が、参照の解釈器で元の `d` を読んだときの k 番目のサブパスと一致するか
    const whole = referenceSubpaths(d);
    parts.forEach((part, k) => {
      const pts = L.pathDataToPoints(part, 0.05, 0);
      if (!pts) return;
      assert.ok(Math.abs(pts[0].x - whole[k][0][0]) < 1e-9 && Math.abs(pts[0].y - whole[k][0][1]) < 1e-9, `始点: ${part}`);
      const end = pts[pts.length - 1];
      const refEnd = whole[k][whole[k].length - 1];
      assert.ok(Math.abs(end.x - refEnd[0]) < 1e-9 && Math.abs(end.y - refEnd[1]) < 1e-9, `終点: ${part}`);
      compared += 1;
    });
  }
  assert.ok(compared > 1000, `比べた数 ${compared}`);
});

// --- 3. ブラウザ自身の解釈との突き合わせ（任意） ------------------------------------

let playwright = null;
if (process.env.DERLE_BROWSER_TEST) {
  try {
    playwright = require("playwright");
  } catch (_) {
    playwright = null;
  }
}

test(
  "splitPathData: ブラウザ（chromium）の getPointAtLength で、元の `d` の中のサブパスと切り出したものが同じ形",
  { skip: playwright ? false : "DERLE_BROWSER_TEST=1 と playwright が要る" },
  async () => {
    const browser = await playwright.chromium.launch();
    try {
      const page = await browser.newPage();
      await page.setContent('<svg xmlns="http://www.w3.org/2000/svg" id="s"></svg>');
      const cases = FUZZ_PATHS.slice(0, 400).map((d) => ({
        d,
        // k 番目のサブパスまでを元のまま並べた `d`（旧実装の単純な分割を繋ぎ直すだけ。
        // 元の文字列の接頭辞なので、ブラウザは元の `d` と同じ解釈をする）
        prefixes: (() => {
          const raw = d.trim().split(/(?=[Mm])/);
          return raw.map((_, k) => raw.slice(0, k + 1).join(""));
        })(),
        parts: L.splitPathData(d),
      }));
      const worst = await page.evaluate((cases) => {
        const svg = document.getElementById("s");
        const NS = "http://www.w3.org/2000/svg";
        const mk = (d) => {
          const p = document.createElementNS(NS, "path");
          p.setAttribute("d", d);
          svg.appendChild(p);
          return p;
        };
        let worstDev = 0;
        let worstCase = null;
        let checked = 0;
        for (const c of cases) {
          for (let k = 0; k < c.parts.length; k += 1) {
            const part = mk(c.parts[k]);
            const ctx = mk(c.prefixes[k]);
            const before = k > 0 ? mk(c.prefixes[k - 1]) : null;
            const partLen = part.getTotalLength();
            const offset = before ? before.getTotalLength() : 0;
            const ctxLen = ctx.getTotalLength();
            // 長さが合う（サブパスの長さは位置に依らない）
            const lenDev = Math.abs(ctxLen - offset - partLen);
            if (lenDev > 1e-3 * Math.max(1, ctxLen)) {
              return { worstDev: Infinity, worstCase: { d: c.d, k, lenDev, part: c.parts[k] } };
            }
            if (partLen > 1e-3) {
              for (let i = 1; i < 20; i += 1) {
                const t = (partLen * i) / 20;
                const a = part.getPointAtLength(t);
                const b = ctx.getPointAtLength(offset + t);
                // getPointAtLength は円弧や曲線の長さを近似で求めるので、前のサブパスの長さ（offset）の
                // 誤差がそのまま位置のずれになる。長さに対する比で見る
                const dev = Math.hypot(a.x - b.x, a.y - b.y) / Math.max(1, ctxLen);
                if (dev > worstDev) {
                  worstDev = dev;
                  worstCase = { d: c.d, k, part: c.parts[k], t };
                }
              }
              checked += 1;
            }
            part.remove();
            ctx.remove();
            if (before) before.remove();
          }
        }
        return { worstDev, worstCase, checked };
      }, cases);
      assert.ok(worst.checked > 500, `比べたサブパス数 ${worst.checked}`);
      // 始点を取り違えると、ずれは直前の現在点の大きさ（ここでは数十〜百）になり、比で 1e-2 を超える。
      // 正しく直したときのずれは長さの近似誤差だけ（実測で比 2.4e-4 程度。円弧で長さ 2,000 超のとき 0.58）
      assert.ok(worst.worstDev < 1e-3, `最大ずれ ${worst.worstDev}: ${JSON.stringify(worst.worstCase)}`);
    } finally {
      await browser.close();
    }
  },
);
