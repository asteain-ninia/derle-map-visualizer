// logic.js（DerleLogic）の純粋関数のユニットテスト。実行: `npm test`（= `node --test`）
//
// DOM が要る関数（getCTM で座標変換する点列化など）はここでは扱わない。そちらは
// logic_test.html（ブラウザで開く目視テスト）で確かめる。
"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const L = require("../logic.js");

const P = (x, y) => ({ x, y });

/** 点 p から折れ線 pts までの最短距離 */
function distToPolyline(p, pts) {
  let best = Infinity;
  for (let i = 0; i < pts.length - 1; i += 1) {
    const a = pts[i];
    const b = pts[i + 1];
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const len2 = dx * dx + dy * dy;
    let t = len2 > 0 ? ((p.x - a.x) * dx + (p.y - a.y) * dy) / len2 : 0;
    t = Math.max(0, Math.min(1, t));
    best = Math.min(best, Math.hypot(p.x - (a.x + dx * t), p.y - (a.y + dy * t)));
  }
  return best;
}

function cubicAt(p0, p1, p2, p3, t) {
  const u = 1 - t;
  return P(
    u * u * u * p0.x + 3 * u * u * t * p1.x + 3 * u * t * t * p2.x + t * t * t * p3.x,
    u * u * u * p0.y + 3 * u * u * t * p1.y + 3 * u * t * t * p2.y + t * t * t * p3.y,
  );
}

function assertPointsClose(actual, expected, eps = 1e-9) {
  assert.equal(actual.length, expected.length, `点数 ${actual.length} != ${expected.length}`);
  for (let i = 0; i < expected.length; i += 1) {
    assert.ok(
      Math.abs(actual[i].x - expected[i].x) <= eps && Math.abs(actual[i].y - expected[i].y) <= eps,
      `点 ${i}: (${actual[i].x}, ${actual[i].y}) != (${expected[i].x}, ${expected[i].y})`,
    );
  }
}

test("Node から読み込める（ブラウザ用の <script> と同じファイル）", () => {
  assert.equal(typeof L.pathDataToPoints, "function");
  assert.equal(globalThis.DerleLogic, L);
});

test("nextFrame は requestAnimationFrame の無い環境でも解決する", async () => {
  await L.nextFrame();
});

// --- パスデータの解釈 ------------------------------------------------------------

test("pathDataToPoints: 直線コマンド（絶対・相対・H/V・Z）", () => {
  assertPointsClose(L.pathDataToPoints("M0 0 L10 0 L10 10 Z", 0.05, 0), [P(0, 0), P(10, 0), P(10, 10), P(0, 0)]);
  assertPointsClose(L.pathDataToPoints("m1 1 l9 0 l0 9 z", 0.05, 0), [P(1, 1), P(10, 1), P(10, 10), P(1, 1)]);
  assertPointsClose(L.pathDataToPoints("M0 0H5V5h-5v-5", 0.05, 0), [P(0, 0), P(5, 0), P(5, 5), P(0, 5), P(0, 0)]);
});

test("pathDataToPoints: M/m の後の座標対は暗黙の L/l", () => {
  assertPointsClose(L.pathDataToPoints("M0 0 10 0 10 10", 0.05, 0), [P(0, 0), P(10, 0), P(10, 10)]);
  assertPointsClose(L.pathDataToPoints("m5 5 10 0 0 10", 0.05, 0), [P(5, 5), P(15, 5), P(15, 15)]);
});

test("pathDataToPoints: 数値の書き方（区切りなし・小数点の連続・指数）", () => {
  // "1.5.5" は 1.5 と .5、"-1-2" は -1 と -2
  assertPointsClose(L.pathDataToPoints("M1.5.5L-1-2", 0.05, 0), [P(1.5, 0.5), P(-1, -2)]);
  assertPointsClose(L.pathDataToPoints("M1e1,2E-1l+3,-.5", 0.05, 0), [P(10, 0.2), P(13, -0.3)]);
});

test("pathDataToPoints: 最初のサブパスだけを返す（1 サブパス = 1 点列の契約）", () => {
  assertPointsClose(L.pathDataToPoints("M0 0L1 0M5 5L6 6", 0.05, 0), [P(0, 0), P(1, 0)]);
});

test("pathDataToPoints: 円弧（A/a）は未対応として null（呼び出し側が旧方式へ退避する）", () => {
  assert.equal(L.pathDataToPoints("M0 0A5 5 0 0 1 10 0", 0.05, 0), null);
  assert.equal(L.pathDataToPoints("M0 0a5 5 0 0 1 10 0", 0.05, 0), null);
  assert.equal(L.pathDataToPoints("", 0.05, 0), null);
  assert.equal(L.pathDataToPoints("M0 0", 0.05, 0), null); // 1 点では線にならない
});

test("pathDataToPoints: 三次ベジェは端点が厳密で、曲線からのずれが許容差以内", () => {
  const p0 = P(0, 0);
  const p1 = P(0, 100);
  const p2 = P(100, 100);
  const p3 = P(100, 0);
  const tol = 0.05;
  const pts = L.pathDataToPoints("M0 0C0 100 100 100 100 0", tol, 0);
  assert.deepEqual(pts[0], p0);
  assert.deepEqual(pts[pts.length - 1], p3);
  assert.ok(pts.length > 10, "曲線が分割されている");
  let worst = 0;
  for (let i = 0; i <= 1000; i += 1) {
    worst = Math.max(worst, distToPolyline(cubicAt(p0, p1, p2, p3, i / 1000), pts));
  }
  // 制御点の弦からのずれ ≤ tol は、曲線の弦からのずれ ≤ 0.75·tol を保証する
  assert.ok(worst <= tol, `最大ずれ ${worst}`);
});

test("pathDataToPoints: 相対 c は現在点からの相対", () => {
  const a = L.pathDataToPoints("M10 10C10 110 110 110 110 10", 0.05, 0);
  const b = L.pathDataToPoints("M10 10c0 100 100 100 100 0", 0.05, 0);
  assertPointsClose(b, a, 1e-9);
});

test("pathDataToPoints: S の第 1 制御点は、直前が三次曲線のときだけ反射する", () => {
  // 直前が C → 反射あり。S は C 0 100 ... と同じ形になる
  const viaS = L.pathDataToPoints("M0 0C0 50 50 50 50 0S100 -50 100 0", 0.01, 0);
  const viaC = L.pathDataToPoints("M0 0C0 50 50 50 50 0C50 -50 100 -50 100 0", 0.01, 0);
  assertPointsClose(viaS, viaC, 1e-9);
  // 直前が L → 反射せず現在点を使う（(0,0) を反射して遠方へ飛ばない。AGENT.md「新方式の実装で踏んだ罠」）
  const afterL = L.pathDataToPoints("M0 0L50 0S100 50 100 0", 0.01, 0);
  const expect = L.pathDataToPoints("M0 0L50 0C50 0 100 50 100 0", 0.01, 0);
  assertPointsClose(afterL, expect, 1e-9);
  const maxY = Math.max(...afterL.map((p) => Math.abs(p.y)));
  assert.ok(maxY < 50, `途中が膨らんでいない (max |y| = ${maxY})`);
});

test("pathDataToPoints: Q/T（二次）と T の反射規則", () => {
  const q = L.pathDataToPoints("M0 0Q50 100 100 0", 0.01, 0);
  // 二次 (0,0)-(50,100)-(100,0) の頂点は (50,50)
  const top = Math.max(...q.map((p) => p.y));
  assert.ok(Math.abs(top - 50) < 0.02, `頂点 y=${top}`);
  const viaT = L.pathDataToPoints("M0 0Q50 100 100 0T200 0", 0.01, 0);
  const viaQ = L.pathDataToPoints("M0 0Q50 100 100 0Q150 -100 200 0", 0.01, 0);
  assertPointsClose(viaT, viaQ, 1e-9);
  // 直前が直線の T は制御点 = 現在点 → 直線になる
  const tAfterL = L.pathDataToPoints("M0 0L10 0T20 0", 0.01, 0);
  assert.ok(tAfterL.every((p) => Math.abs(p.y) < 1e-12));
});

test("pathDataToPoints: minChord で間引いても始点と終点は残る", () => {
  const d = "M0 0" + Array.from({ length: 100 }, (_, i) => ` L${(i + 1) * 0.1} 0`).join("");
  const pts = L.pathDataToPoints(d, 0.05, 1);
  assert.deepEqual(pts[0], P(0, 0));
  assert.ok(Math.abs(pts[pts.length - 1].x - 10) < 1e-9);
  for (let i = 1; i < pts.length - 1; i += 1) {
    assert.ok(pts[i].x - pts[i - 1].x >= 1 - 1e-9, "最小間隔以上");
  }
});

test("splitPathData: M/m ごとに分ける", () => {
  assert.deepEqual(L.splitPathData("M0 0L1 1M2 2L3 3"), ["M0 0L1 1", "M2 2L3 3"]);
  assert.deepEqual(L.splitPathData(""), []);
  assert.deepEqual(L.splitPathData(null), []);
  assert.deepEqual(L.splitPathData("  M0 0L1 1  "), ["M0 0L1 1"]);
});

// --- 幾何の部品 ------------------------------------------------------------------

test("distance / distanceSq", () => {
  assert.equal(L.distance(P(0, 0), P(3, 4)), 5);
  assert.equal(L.distanceSq(P(0, 0), P(3, 4)), 25);
  assert.equal(L.distanceSq(P(1, 1), P(1, 1)), 0);
});

test("removeSequentialDuplicates: 連続する近接点だけを落とす", () => {
  const out = L.removeSequentialDuplicates([P(0, 0), P(0, 0), P(1, 0), P(1, 1e-9), P(0, 0)], 1e-6);
  assert.deepEqual(out, [P(0, 0), P(1, 0), P(0, 0)]);
  assert.deepEqual(L.removeSequentialDuplicates([], 1e-6), []);
});

test("segmentIntersections: 交差・平行・同一直線上の重なり", () => {
  assert.deepEqual(L.segmentIntersections(P(0, 0), P(2, 2), P(0, 2), P(2, 0)), [{ t1: 0.5, t2: 0.5 }]);
  assert.equal(L.segmentIntersections(P(0, 0), P(1, 0), P(0, 1), P(1, 1)), null); // 平行
  assert.equal(L.segmentIntersections(P(0, 0), P(1, 1), P(2, 0), P(3, -1)), null); // 延長上でのみ交わる
  // 端点で接する
  const touch = L.segmentIntersections(P(0, 0), P(1, 0), P(1, 0), P(1, 1));
  assert.deepEqual(touch, [{ t1: 1, t2: 0 }]);
  // 同一直線上の重なり → 重なり区間の両端
  const overlap = L.segmentIntersections(P(0, 0), P(4, 0), P(1, 0), P(6, 0));
  assert.ok(overlap && overlap.length >= 1);
  const ts = overlap.map((h) => h.t1).sort((a, b) => a - b);
  assert.ok(Math.abs(ts[0] - 0.25) < 1e-12);
  assert.ok(Math.abs(ts[ts.length - 1] - 1) < 1e-12);
});

test("pointSegmentProjection: 線分の外に落ちる射影は null", () => {
  assert.deepEqual(L.pointSegmentProjection(P(1, 1), P(0, 0), P(2, 0)), { t: 0.5, distance: 1 });
  assert.equal(L.pointSegmentProjection(P(3, 1), P(0, 0), P(2, 0)), null);
  assert.equal(L.pointSegmentProjection(P(3, 1), P(0, 0), P(0, 0)), null);
});

test("polygonArea: 符号付き面積（SVG の y 下向き座標で時計回りが正）", () => {
  const pts = [P(0, 0), P(4, 0), P(4, 3), P(0, 3)];
  assert.equal(L.polygonArea([0, 1, 2, 3], pts), 12);
  assert.equal(L.polygonArea([3, 2, 1, 0], pts), -12);
});

test("pointInRingIndices / ringBounds", () => {
  // L 字形
  const pts = [P(0, 0), P(4, 0), P(4, 1), P(1, 1), P(1, 4), P(0, 4)];
  const ring = [0, 1, 2, 3, 4, 5];
  assert.equal(L.pointInRingIndices(0.5, 0.5, ring, pts), true);
  assert.equal(L.pointInRingIndices(3, 0.5, ring, pts), true);
  assert.equal(L.pointInRingIndices(3, 3, ring, pts), false); // 凹部
  assert.equal(L.pointInRingIndices(-1, 2, ring, pts), false);
  assert.deepEqual(L.ringBounds(ring, pts), { minX: 0, minY: 0, maxX: 4, maxY: 4 });
});

test("properlyCrosses / pointOnOpenSegment", () => {
  assert.equal(L.properlyCrosses(P(0, 0), P(2, 2), P(0, 2), P(2, 0)), true);
  assert.equal(L.properlyCrosses(P(0, 0), P(1, 0), P(1, 0), P(1, 1)), false); // 端点で接するだけ
  assert.equal(L.pointOnOpenSegment(P(1, 0), P(0, 0), P(2, 0), 1e-6), true);
  assert.equal(L.pointOnOpenSegment(P(0, 0), P(0, 0), P(2, 0), 1e-6), false); // 端点は除く
  assert.equal(L.pointOnOpenSegment(P(1, 0.1), P(0, 0), P(2, 0), 1e-6), false);
});

// --- 線分 → グラフ → 面 ---------------------------------------------------------

const seg = (ax, ay, bx, by, source = "border") => ({ a: P(ax, ay), b: P(bx, by), source });

test("splitSegmentsAtIntersections: X 字の交差で 4 本に分ける", async () => {
  const out = await L.splitSegmentsAtIntersections([seg(0, 0, 2, 2), seg(0, 2, 2, 0)], 10, 0);
  assert.equal(out.length, 4);
  for (const s of out) {
    const touchesCenter = [s.a, s.b].some((p) => Math.abs(p.x - 1) < 1e-12 && Math.abs(p.y - 1) < 1e-12);
    assert.ok(touchesCenter, "どの線分も交点 (1,1) で終わる");
  }
});

test("splitSegmentsAtIntersections: 端点が他線分の途中に近接していれば、その位置で分割する", async () => {
  // T 字: 縦線の下端 (1, 0.3) が横線から 0.3 離れている
  const segs = [seg(0, 0, 2, 0), seg(1, 0.3, 1, 2)];
  const noSnap = await L.splitSegmentsAtIntersections(segs, 10, 0);
  assert.equal(noSnap.length, 2);
  const snapped = await L.splitSegmentsAtIntersections(segs, 10, 0.5);
  assert.equal(snapped.length, 3);
});

test("buildGraphFromSegments: 吸着許容値以内の端点を 1 点に併合する", async () => {
  const graph = await L.buildGraphFromSegments([seg(0, 0, 1, 0), seg(1.05, 0, 2, 0)], 0.1);
  assert.equal(graph.points.length, 3);
  assert.equal(graph.edges.length, 2);
  const loose = await L.buildGraphFromSegments([seg(0, 0, 1, 0), seg(1.05, 0, 2, 0)], 0.01);
  assert.equal(loose.points.length, 4);
});

test("extractPolygonsFromGraph: 2 つ並んだ正方形から 2 面を取り、外側の面を除く", async () => {
  const segs = [
    seg(0, 0, 1, 0), seg(1, 0, 2, 0), seg(2, 0, 2, 1), seg(2, 1, 1, 1),
    seg(1, 1, 0, 1), seg(0, 1, 0, 0), seg(1, 0, 1, 1),
  ];
  const graph = await L.buildGraphFromSegments(segs, 1e-6);
  const { polygons, outerFaces } = L.extractPolygonsFromGraph(graph, 0);
  assert.equal(polygons.length, 2);
  assert.equal(outerFaces.length, 1);
  polygons.forEach((p) => assert.ok(Math.abs(Math.abs(p.area) - 1) < 1e-12));
  assert.ok(Math.abs(Math.abs(outerFaces[0].area) - 2) < 1e-12);
});

test("extractPolygonsFromGraph: 連結成分ごとに外側の面を除く（離れた 2 つの島）", async () => {
  const square = (ox) => [seg(ox, 0, ox + 1, 0), seg(ox + 1, 0, ox + 1, 1), seg(ox + 1, 1, ox, 1), seg(ox, 1, ox, 0)];
  const graph = await L.buildGraphFromSegments([...square(0), ...square(5)], 1e-6);
  const { polygons, outerFaces } = L.extractPolygonsFromGraph(graph, 0);
  assert.equal(polygons.length, 2);
  assert.equal(outerFaces.length, 2);
});

test("bridgeOpenBorderGaps: 国境の行き止まりを近くの点へ繋いで面を閉じる", async () => {
  // 正方形の右辺だけ 1.5 だけ届いていない
  const segs = [seg(0, 0, 10, 0), seg(10, 0, 10, 8.5), seg(10, 10, 0, 10), seg(0, 10, 0, 0)];
  const graph = await L.buildGraphFromSegments(segs, 1e-6);
  assert.equal(L.extractPolygonsFromGraph(graph, 0).polygons.length, 0);
  const { bridges } = await L.bridgeOpenBorderGaps(graph, { radius: 2.5, borderRadius: 6, extraBridges: [] });
  assert.equal(bridges.length, 1);
  assert.equal(L.extractPolygonsFromGraph(graph, 0).polygons.length, 1);
});

test("computeContainment / computeOwnership: 内包の親と、領有する側の種別", () => {
  const pts = [
    P(0, 0), P(10, 0), P(10, 10), P(0, 10), // 0-3 外
    P(4, 4), P(6, 4), P(6, 6), P(4, 6), // 4-7 中の島
  ];
  const polys = [
    { points: [0, 1, 2, 3], area: 100 },
    { points: [4, 5, 6, 7], area: 4 },
  ];
  assert.deepEqual(L.computeContainment(polys, pts), [-1, 0]);
  assert.deepEqual(L.computeOwnership(polys, pts, ["海まで国境線", "島国"]).owners, [-1, 0]);
  assert.deepEqual(L.computeOwnership(polys, pts, ["海岸線まで国境界", "島国"]).owners, [-1, -1]);
});

test("computeContainment: 凹部にはまって辺を接しているだけの隣国は親にしない", () => {
  // L 字の国と、その凹部を埋める正方形の国。正方形の外接矩形は L 字の外接矩形に収まるが、
  // 共有していない頂点 (10,10) は L 字の外にある
  const pts = [P(0, 0), P(10, 0), P(10, 5), P(5, 5), P(5, 10), P(0, 10), P(10, 10)];
  const polys = [
    { points: [0, 1, 2, 3, 4, 5], area: 75 },
    { points: [3, 2, 6, 4], area: 25 },
  ];
  assert.deepEqual(L.computeContainment(polys, pts), [-1, -1]);
});

// --- 経度180度の継ぎ目 -----------------------------------------------------------

test("区間演算: merge / intersect / subtract / contain", () => {
  assert.deepEqual(L.mergeIntervals([[5, 6], [0, 2], [1, 3], [3.5, 4]], 0), [[0, 3], [3.5, 4], [5, 6]]);
  assert.deepEqual(L.mergeIntervals([[0, 2], [2.4, 3]], 0.5), [[0, 3]]);
  assert.deepEqual(L.mergeIntervals([[1, 1]], 0), []); // 長さ 0 は捨てる
  assert.deepEqual(L.intersectIntervals([[0, 5], [8, 10]], [[3, 9]]), [[3, 5], [8, 9]]);
  assert.deepEqual(L.subtractIntervals([[0, 10]], [[2, 3], [5, 6]]), [[0, 2], [3, 5], [6, 10]]);
  assert.equal(L.intervalsContain([[0, 1]], 0.5), true);
  assert.equal(L.intervalsContain([[0, 1]], 1), false); // 端は含まない
});

test("seamSideOf / snapPointsToSeam: 外側は内側の 20 倍まで端とみなす", () => {
  const frame = { x: 0, y: 0, width: 360, height: 180 };
  assert.equal(L.seamSideOf(P(0.05, 10), frame, 0.1), "W");
  assert.equal(L.seamSideOf(P(-1.9, 10), frame, 0.1), "W");
  assert.equal(L.seamSideOf(P(-2.1, 10), frame, 0.1), null);
  assert.equal(L.seamSideOf(P(359.95, 10), frame, 0.1), "E");
  assert.equal(L.seamSideOf(P(0.2, 10), frame, 0.1), null);
  const pts = [P(0.05, 1), P(-1, 2), P(180, 3), P(360.5, 4)];
  assert.equal(L.snapPointsToSeam(pts, frame, 0.1), 3);
  assert.deepEqual(pts.map((p) => p.x), [0, 0, 180, 360]);
  assert.equal(L.snapPointsToSeam(pts, frame, 0.1), 0); // 冪等
});

test("seamStrokeParts: 切断区間だけを落とし、本物の海岸線の部分は残す", () => {
  const frame = { x: 0, y: 0, width: 400, height: 200 };
  const seam = { cut: [[20, 70]] };
  // 西端 y=20–80 の辺のうち y=70–80 は本物の海岸線（AGENT.md「球面出力で経度180度の切断線を描かない」）
  const parts = L.seamStrokeParts(P(0, 20), P(0, 80), frame, 0.1, seam);
  assert.equal(parts.length, 1);
  assert.deepEqual(parts[0], [P(0, 70), P(0, 80)]);
  assert.equal(L.isArtificialSeamEdge(P(0, 20), P(0, 80), frame, 0.1, seam), false);
  assert.equal(L.isArtificialSeamEdge(P(0, 30), P(0, 60), frame, 0.1, seam), true);
  // 地図端に乗らない辺はそのまま
  assert.deepEqual(L.seamStrokeParts(P(10, 20), P(10, 80), frame, 0.1, seam), [[P(10, 20), P(10, 80)]]);
});

test("groupSeamPolygons: 西端と東端で同じ切断線に接するピースをまとめる", () => {
  const frame = { x: 0, y: 0, width: 400, height: 200 };
  const pts = [
    // 西端に接する島の半分 (y 50–100)
    P(0, 50), P(20, 50), P(20, 100), P(0, 100),
    // 東端に接するもう半分
    P(380, 50), P(400, 50), P(400, 100), P(380, 100),
    // 東端の別の島 (y 150–160)
    P(390, 150), P(400, 150), P(400, 160), P(390, 160),
  ];
  const polys = [{ points: [0, 1, 2, 3] }, { points: [4, 5, 6, 7] }, { points: [8, 9, 10, 11] }];
  const { groups } = L.groupSeamPolygons(polys, pts, frame, { tol: 0.1 });
  assert.deepEqual(groups, [[0, 1]]);
});

test("isFrameOnlySeaFace: 枠の仮想辺と海岸線だけで閉じた面", () => {
  const src = (pairs) => (a, b) => pairs[`${Math.min(a, b)}-${Math.max(a, b)}`];
  const ring = [0, 1, 2];
  assert.equal(L.isFrameOnlySeaFace(ring, src({ "0-1": new Set(["frame"]), "1-2": new Set(["coast"]), "0-2": new Set(["coast"]) })), true);
  assert.equal(L.isFrameOnlySeaFace(ring, src({ "0-1": new Set(["frame"]), "1-2": new Set(["border"]), "0-2": new Set(["coast"]) })), false);
  assert.equal(L.isFrameOnlySeaFace(ring, src({ "0-1": new Set(["coast"]), "1-2": new Set(["coast"]), "0-2": new Set(["coast"]) })), false);
});
