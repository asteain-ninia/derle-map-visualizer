// 国タイプ判定（classifyPolygons と、その入力を作る部品）のユニットテスト。実行: `npm test`
//
// 判定の順序は index.html の buildGraph / extractPolygons と同じ（下の classifyMap）。
// 地図枠・経度180度の継ぎ目・橋渡しは扱わない（それぞれ logic.test.js で見ている）。
"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const L = require("../logic.js");

const P = (x, y) => ({ x, y });
const seg = (ax, ay, bx, by, source) => ({ a: P(ax, ay), b: P(bx, by), source });
const loop = (pts, source) => pts.map((p, i) => {
  const q = pts[(i + 1) % pts.length];
  return seg(p[0], p[1], q[0], q[1], source);
});
const box = (x0, y0, x1, y1, source) => loop([[x0, y0], [x1, y0], [x1, y1], [x0, y1]], source);

/**
 * 線分（source: "border" / "coast"）から国タイプを求める。index.html の
 * buildGraph（交差分割 → 国境線だけで面 → 共有する海 → 海岸線の除外 → グラフ）と
 * extractPolygons（面 → 辺の索引 → 海岸線の内外 → 国タイプ）の順に呼ぶ。
 * 戻り値は「代表点 → 種別」を引く関数（代表点を含む最小の面の種別）。
 */
async function classifyMap(segments, snap = 0.01) {
  const split = await L.splitSegmentsAtIntersections(segments, 10, snap);
  const rawCoast = split.filter((s) => s.source === "coast");
  const borderCoastPoints = L.collectBorderCoastPoints(split, snap);
  const borderOnly = split.filter((s) => s.source === "border");
  const coastOnly = split.filter((s) => s.source === "coast");
  const borderGraph = await L.buildGraphFromSegments(borderOnly, snap);
  const borderFaces = L.extractPolygonsFromGraph(borderGraph, 0).polygons;
  const shared = L.collectSharedSeaPolygons(borderGraph, borderFaces, snap, coastOnly);
  const keptCoast = L.filterCoastSegments(
    coastOnly, shared.polygonData, borderGraph.points, shared.sharedPolygons, shared.tolerance,
  );
  const graph = await L.buildGraphFromSegments(borderOnly.concat(keptCoast), snap);
  const polygons = L.extractPolygonsFromGraph(graph, 0).polygons.map((f, i) => ({ id: i + 1, points: f.points }));
  const edgeIndex = L.buildPolygonEdgeIndex(polygons, graph);
  const polygonData = L.buildPolygonData(polygons, graph.points);
  const coastInfo = {
    borderCoastPoints: L.mapPointsToGraphIndices(borderCoastPoints, graph.points, snap),
    coastInsidePolygons: L.collectCoastInsidePolygons(
      rawCoast, polygonData, graph.points, Math.max(L.GEOM_EPS, snap * 0.1),
    ),
  };
  const types = L.classifyPolygons(
    polygons, graph, edgeIndex.edgeSources, { ...shared, points: borderGraph.points }, coastInfo,
  );
  const typeAt = (x, y) => {
    let best = null;
    polygons.forEach((poly) => {
      if (!L.pointInPolygon(P(x, y), poly.points, graph.points, 0)) return;
      const area = Math.abs(L.polygonArea(poly.points, graph.points));
      if (!best || area < best.area) best = { area, type: types.get(String(poly.id)) };
    });
    return best ? best.type : null;
  };
  return { typeAt, types, polygons, shared };
}

test("合成した地図で、6 種類の国タイプを判定する", async () => {
  const segments = [
    // 大陸: 海岸線の正方形を国境線 x=50 で東西に分け、西側の中に国境線だけで囲まれた国を置く
    ...box(0, 0, 100, 100, "coast"),
    seg(50, 0, 50, 100, "border"),
    ...box(10, 40, 30, 60, "border"),
    // 海岸線だけの島
    ...box(200, 0, 220, 20, "coast"),
    // 海上の国境線で囲んだ島（島の海岸線は領内として除かれる）
    ...box(300, 0, 320, 20, "coast"),
    ...box(280, -20, 340, 40, "border"),
    // 海上まで引いた国境線で分けた大陸（国境線が海岸線と交わる）
    ...box(400, 0, 500, 100, "coast"),
    ...box(380, -20, 520, 120, "border"),
    seg(450, -20, 450, 120, "border"),
    // 2 国が共有する海: 海岸線の上で終わる国境線が、海上の国境線で囲んだ海の中にある
    ...box(620, 0, 700, 100, "coast"),
    ...box(600, -20, 720, 120, "border"),
    seg(660, 100, 660, 120, "border"),
  ];
  const { typeAt } = await classifyMap(segments);
  assert.equal(typeAt(25, 20), "海岸線まで国境界", "大陸の西側");
  assert.equal(typeAt(75, 50), "海岸線まで国境界", "大陸の東側");
  assert.equal(typeAt(20, 50), "内陸国", "国境線だけで囲まれ、海岸線に触れない");
  assert.equal(typeAt(210, 10), "島国");
  assert.equal(typeAt(310, 10), "完全島国", "海上の国境線の内側に海岸線がある");
  assert.equal(typeAt(285, 30), "完全島国", "島のまわりの海も同じ面");
  assert.equal(typeAt(420, 50), "海まで国境線", "西側");
  assert.equal(typeAt(390, -10), "海まで国境線", "西側の海の部分");
  assert.equal(typeAt(480, 50), "海まで国境線", "東側");
  assert.equal(typeAt(610, 50), "共有する海");
  assert.equal(typeAt(660, 50), "島国", "共有する海の中の陸は海岸線が残るので別の面");
});

test("collectSharedSeaPolygons: 自由端が海岸線の上にある国境線だけが共有する海を作る（「髭」は数えない）", async () => {
  // 回帰: AGENT.md「1 単位の『髭』で国の扱いが反転する」。自由端が海岸線から離れている短い国境線
  // （作図の残り）を「海岸線で終わる国境」と取り違えると、その国が丸ごと共有する海に化ける。
  const coast = box(20, 20, 80, 80, "coast");
  const sea = box(0, 0, 100, 100, "border");
  const onCoast = [...sea, seg(50, 80, 50, 100, "border")];
  // 海上の国境線から内側へ 1 だけ出た線。自由端 (50, 99) は海岸線から 19 離れている。
  const whisker = [...sea, seg(50, 100, 50, 99, "border")];
  for (const [name, borders, expected] of [["海岸線で終わる", onCoast, 1], ["髭（長さ 1）", whisker, 0]]) {
    const graph = await L.buildGraphFromSegments(borders, 0.01);
    const faces = L.extractPolygonsFromGraph(graph, 0).polygons;
    const shared = L.collectSharedSeaPolygons(graph, faces, 0.01, coast);
    assert.equal(shared.sharedPolygons.size, expected, name);
  }
});

test("classifyPolygons: 判定表（入力を直接与える）", () => {
  // 正方形 1 つ。辺の sources と、国境線・海岸線の交点、領内の海岸線を切り替えて判定を見る。
  const graph = { points: [P(0, 0), P(10, 0), P(10, 10), P(0, 10)] };
  const poly = { id: 7, points: [0, 1, 2, 3] };
  const sources = (list) => {
    const m = new Map();
    [[0, 1], [1, 2], [2, 3], [3, 0]].forEach(([a, b], i) => m.set(L.edgeKey(a, b), new Set(list[i])));
    return m;
  };
  const allBorder = sources([["border"], ["border"], ["border"], ["border"]]);
  const mixed = sources([["border"], ["coast"], ["coast"], ["coast"]]);
  const allCoast = sources([["coast"], ["coast"], ["coast"], ["coast"]]);
  const coastInfo = (touch, inside) => ({
    borderCoastPoints: new Set(touch ? [2] : []),
    coastInsidePolygons: new Set(inside ? [0] : []),
  });
  const typeOf = (edgeSources, info, shared = null, isCutEdge) =>
    L.classifyPolygons([poly], graph, edgeSources, shared, info, isCutEdge).get("7");

  assert.equal(typeOf(allCoast, coastInfo(false, false)), "島国");
  assert.equal(typeOf(mixed, coastInfo(false, false)), "海岸線まで国境界");
  assert.equal(typeOf(allBorder, coastInfo(true, false)), "海まで国境線");
  assert.equal(typeOf(allBorder, coastInfo(false, true)), "完全島国");
  assert.equal(typeOf(allBorder, coastInfo(false, false)), "内陸国");
  assert.equal(typeOf(allBorder, null), "海まで国境線", "海岸線レイヤーが無いときは区別しない");
  assert.equal(typeOf(new Map(), coastInfo(false, false)), "不明");

  // 共有する海: 海岸線の辺の中点が、共有する海と判定された面の中にある
  const sharedSea = {
    polygonData: L.buildPolygonData([{ points: [0, 1, 2, 3] }], graph.points),
    points: graph.points,
    sharedPolygons: new Set([0]),
    tolerance: 1e-6,
  };
  assert.equal(typeOf(mixed, coastInfo(false, false), sharedSea), "共有する海");
  assert.equal(typeOf(mixed, coastInfo(false, false), { ...sharedSea, sharedPolygons: new Set() }), "海岸線まで国境界");

  // 地図端の切断線（isCutEdge が真の辺）は境界として数えない
  const cutBorder = (a, b) => L.edgeKey(a, b) === L.edgeKey(0, 1);
  assert.equal(typeOf(mixed, coastInfo(false, false), null, cutBorder), "島国");
});

test("classifyPolygons: 継ぎ目でつないだ地物は全ピースの辺で判定する", () => {
  const graph = { points: [P(0, 0), P(10, 0), P(10, 10), P(0, 10), P(90, 0), P(100, 0), P(100, 10), P(90, 10)] };
  const edgeSources = new Map();
  [[0, 1], [1, 2], [2, 3], [3, 0]].forEach(([a, b]) => edgeSources.set(L.edgeKey(a, b), new Set(["coast"])));
  [[4, 5], [5, 6], [6, 7], [7, 4]].forEach(([a, b]) => edgeSources.set(L.edgeKey(a, b), new Set(["border"])));
  const feature = { id: 1, rings: [[0, 1, 2, 3], [4, 5, 6, 7]] };
  const info = { borderCoastPoints: new Set(), coastInsidePolygons: new Set() };
  assert.equal(L.classifyPolygons([feature], graph, edgeSources, null, info).get("1"), "海岸線まで国境界");
  assert.deepEqual(L.getPolygonRings(feature), feature.rings);
  assert.deepEqual(L.getPolygonRings({ points: [1, 2, 3] }), [[1, 2, 3]]);
});

test("pointInPolygon / findPolygonIndex: 境界から tol 以内は内側、複数リングの地物はどれかに入れば内側", () => {
  const pts = [P(0, 0), P(10, 0), P(10, 10), P(0, 10), P(20, 0), P(30, 0), P(30, 10), P(20, 10)];
  assert.equal(L.pointInPolygon(P(5, 5), [0, 1, 2, 3], pts, 0), true);
  assert.equal(L.pointInPolygon(P(10.05, 5), [0, 1, 2, 3], pts, 0), false);
  assert.equal(L.pointInPolygon(P(10.05, 5), [0, 1, 2, 3], pts, 0.1), true);
  const data = L.buildPolygonData([{ points: [4, 5, 6, 7] }, { rings: [[0, 1, 2, 3], [4, 5, 6, 7]] }], pts);
  assert.equal(L.findPolygonIndex(P(25, 5), data, pts, 0), 0, "先に見つかった面");
  assert.equal(L.findPolygonIndex(P(5, 5), data, pts, 0), 1);
  assert.equal(L.findPolygonIndex(P(15, 5), data, pts, 0), -1);
});

test("collectBorderCoastPoints / mapPointsToGraphIndices: 国境線と海岸線が共有する端点", () => {
  const segments = [seg(0, 0, 10, 0, "border"), seg(10, 0, 10, 10, "coast"), seg(0, 5, 0, 10, "coast")];
  const pts = L.collectBorderCoastPoints(segments, 0.01);
  assert.deepEqual(pts, [P(10, 0)]);
  const graphPoints = [P(0, 0), P(10.005, 0), P(10, 10)];
  assert.deepEqual([...L.mapPointsToGraphIndices(pts, graphPoints, 0.01)], [1]);
  assert.deepEqual([...L.mapPointsToGraphIndices(pts, graphPoints, 0.001)], [], "吸着許容値の外は対応しない");
});
