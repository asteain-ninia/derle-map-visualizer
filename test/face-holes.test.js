// 内包された面の穴（画面表示と .gimoza の書き出しが共有する判定）のユニットテスト。実行: `npm test`
// 現状.md 重要度低 1「包領を画面上で排他にできていない」への対応で、画面の穴は
// DerleLogic.computeFeatureHoles → computeFaceHoles、書き出しの穴は assignContainedHoles で決まる。
"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");

const L = require("../logic.js");
const G = require("../gimoza.js");

const P = (x, y) => ({ x, y });

/** 点 (x, y) のまわりの巻き数（リング群の合計）。SVG の nonzero の塗りは巻き数 ≠ 0 のところ。 */
function windingNumber(x, y, rings, points) {
  let wn = 0;
  for (const ring of rings) {
    for (let i = 0; i < ring.length; i += 1) {
      const a = points[ring[i]];
      const b = points[ring[(i + 1) % ring.length]];
      const cross = (b.x - a.x) * (y - a.y) - (x - a.x) * (b.y - a.y);
      if (a.y <= y) {
        if (b.y > y && cross > 0) wn += 1;
      } else if (b.y <= y && cross < 0) {
        wn -= 1;
      }
    }
  }
  return wn;
}

/** 点 (x, y) を evenodd で数えたときの交差回数の偶奇（選択の暗転は viewBox の矩形＋選択の path を evenodd で塗る）。 */
function evenOddInside(x, y, rings, points) {
  let inside = false;
  for (const ring of rings) if (L.pointInRingIndices(x, y, ring, points)) inside = !inside;
  return inside;
}

/**
 * 合成の地図。SVG 座標（y 下向き）。
 *   A … 0..100 の正方形。中に B（包領）と C を持つ
 *   B … 10..30 の正方形（A の穴）
 *   C … 50..90 の正方形（A の穴）。中に D を持つ
 *   D … 60..70 の正方形（C の穴。A の穴にはしない＝穴の中の飛び地は親の側で扱わない）
 *   F … 200..300 の正方形。橋渡しのスリットで島 G を外したポケットを持つ（自己接触リング）
 *   G … 240..260 の正方形（F の外周には内包されるが、F のポケットが既に抜いているので穴にしない）
 *   H … A と辺 3-2 を共有する隣の面（内包ではない）
 */
function buildScene() {
  const points = [
    P(0, 0), P(100, 0), P(100, 100), P(0, 100), // 0-3: A
    P(10, 10), P(30, 10), P(30, 30), P(10, 30), // 4-7: B
    P(50, 50), P(90, 50), P(90, 90), P(50, 90), // 8-11: C
    P(60, 60), P(70, 60), P(70, 70), P(60, 70), // 12-15: D
    P(200, 0), P(300, 0), P(300, 100), P(200, 100), // 16-19: F
    P(240, 40), P(240, 60), P(260, 60), P(260, 40), // 20-23: G（F のポケットとしては逆回り）
    P(100, 200), P(0, 200), // 24-25: H
  ];
  const ring = (ids) => ({ points: ids, area: L.polygonArea(ids, points) });
  const faces = [
    ring([0, 1, 2, 3]), // 0: A（時計回り＝正）
    ring([4, 5, 6, 7]), // 1: B（A と同じ向き）
    ring([11, 10, 9, 8]), // 2: C（A と逆向き）
    ring([12, 13, 14, 15]), // 3: D
    ring([16, 17, 18, 19, 16, 20, 21, 22, 23, 20]), // 4: F（16 から 20 へスリット、G を逆回りに回って戻る）
    ring([20, 23, 22, 21]), // 5: G
    ring([3, 2, 24, 25]), // 6: H
  ];
  return { points, faces };
}

test("computeFaceHoles: 直接内包された面だけを穴にし、ポケットで抜けている面と隣の面は穴にしない", () => {
  const { points, faces } = buildScene();
  const { holes, parents, holeRings } = L.computeFaceHoles(faces, points);
  assert.deepEqual(parents, [-1, 0, 0, 2, -1, 4, -1]);
  assert.deepEqual(holes, [[1, 2], [], [3], [], [], [], []]);
  // 穴の形は子の外周（ポケットを外して整えたもの）
  assert.deepEqual(holeRings[0], [[4, 5, 6, 7], [11, 10, 9, 8]]);
  assert.deepEqual(holeRings[2], [[12, 13, 14, 15]]);
});

test("computeFaceHoles: 退化した面があっても index は元の並びで返す", () => {
  const { points, faces } = buildScene();
  const withDegenerate = [{ points: [0, 1, 0], area: 0 }, ...faces];
  const { holes, parents } = L.computeFaceHoles(withDegenerate, points);
  assert.deepEqual(parents, [-1, -1, 1, 1, 3, -1, 5, -1]);
  assert.deepEqual(holes, [[], [2, 3], [], [4], [], [], [], []]);
});

test("computeFeatureHoles: 穴は外周と逆回りで、nonzero でも evenodd でも内側の面が抜ける", () => {
  const { points, faces } = buildScene();
  const pieces = faces.map((f, i) => ({ featureId: String(i + 1), points: f.points, area: f.area }));
  const holesOf = L.computeFeatureHoles(pieces, points);
  assert.deepEqual([...holesOf.keys()].sort(), ["1", "3"]);
  assert.deepEqual(holesOf.get("1").map((h) => h.childId), ["2", "3"]);
  assert.deepEqual(holesOf.get("3").map((h) => h.childId), ["4"]);
  assert.equal(holesOf.get("1")[0].area, 400);
  assert.equal(holesOf.get("1")[1].area, 1600);

  const outerSign = Math.sign(L.polygonArea(pieces[0].points, points));
  holesOf.get("1").forEach((h) => {
    assert.equal(Math.sign(L.polygonArea(h.ring, points)), -outerSign, `${h.childId} の穴は外周と逆回り`);
  });

  // 画面の path は「外周＋穴」。A の塗りは B・C（C の中の D も含む）を抜き、A だけの所は塗る
  const pathA = [pieces[0].points, ...holesOf.get("1").map((h) => h.ring)];
  assert.notEqual(windingNumber(40, 20, pathA, points), 0, "A だけの所は A の塗り");
  assert.equal(windingNumber(20, 20, pathA, points), 0, "包領 B の中は A の塗りにならない");
  assert.equal(windingNumber(55, 55, pathA, points), 0, "C の中は A の塗りにならない");
  assert.equal(windingNumber(65, 65, pathA, points), 0, "C の中の D も A の塗りにならない");
  assert.ok(evenOddInside(40, 20, pathA, points));
  assert.ok(!evenOddInside(20, 20, pathA, points));
  assert.ok(!evenOddInside(65, 65, pathA, points));

  // C は D を抜く
  const pathC = [pieces[2].points, ...holesOf.get("3").map((h) => h.ring)];
  assert.notEqual(windingNumber(55, 55, pathC, points), 0);
  assert.equal(windingNumber(65, 65, pathC, points), 0);

  // F は穴を持たない（ポケットが自己接触リングの中で既に逆回りなので、元のリングのまま G が抜ける）
  assert.notEqual(windingNumber(220, 20, [pieces[4].points], points), 0);
  assert.equal(windingNumber(250, 50, [pieces[4].points], points), 0);
});

test("computeFeatureHoles: 選択の暗転（viewBox の矩形＋選んだ path を evenodd）で包領が暗くなり、両方選べば明るい", () => {
  const { points, faces } = buildScene();
  const pieces = faces.map((f, i) => ({ featureId: String(i + 1), points: f.points, area: f.area }));
  const holesOf = L.computeFeatureHoles(pieces, points);
  const frame = [...points, P(-1000, -1000), P(1000, -1000), P(1000, 1000), P(-1000, 1000)];
  const rect = [26, 27, 28, 29];
  const pathA = [pieces[0].points, ...holesOf.get("1").map((h) => h.ring)];
  const pathB = [pieces[1].points];
  const darkened = (x, y, rings) => evenOddInside(x, y, [rect, ...rings], frame);
  // A だけを選ぶ: A は明るく、包領 B は暗い（以前は B まで明るく出ていた）
  assert.ok(!darkened(40, 20, pathA));
  assert.ok(darkened(20, 20, pathA));
  // A と B を選ぶ: どちらも明るい
  assert.ok(!darkened(40, 20, [...pathA, ...pathB]));
  assert.ok(!darkened(20, 20, [...pathA, ...pathB]));
});

test("computeFeatureHoles: 経度180度の継ぎ目で同じ地物になったピースどうしは穴にしない", () => {
  const { points, faces } = buildScene();
  // A と B を同じ地物にする（実データでは起きないが、起きても自分の中に穴を開けない）
  const pieces = faces.map((f, i) => ({ featureId: i === 1 ? "1" : String(i + 1), points: f.points, area: f.area }));
  const holesOf = L.computeFeatureHoles(pieces, points);
  assert.deepEqual(holesOf.get("1").map((h) => h.childId), ["3"]);
});

test("assignContainedHoles: .gimoza の書き出しの穴（ポケットを除く）と画面の穴が同じになる", () => {
  const { points, faces } = buildScene();
  const polygons = faces.map((f, i) => ({ id: i + 1, area: f.area, points: f.points }));
  const edges = [];
  const seen = new Set();
  faces.forEach((f) => {
    const r = L.cleanRing(f.points);
    for (let i = 0; i < r.length; i += 1) {
      const a = r[i];
      const b = r[(i + 1) % r.length];
      const key = L.edgeKey(a, b);
      if (a === b || seen.has(key)) continue;
      seen.add(key);
      edges.push({ a, b, sources: ["border"] });
    }
  });
  // DOM の無い Node では、ラベルと陸地の判定を空にした最小の文書を渡す
  const svgDoc = { getElementById: () => null, querySelector: () => null, getElementsByTagName: () => [] };
  const manifest = G.buildManifest({
    points,
    edges,
    polygons,
    types: faces.map(() => "内陸国"),
    svgDoc,
    // 冗長点の除去を切って、書き出しと画面で同じリングを比べる
    options: { worldRect: { x: -1000, y: -500, width: 2000, height: 1000 }, isLand: () => true, cleanupEpsilon: 0 },
  });
  const exported = new Map();
  manifest.histories.forEach((h) => {
    const feature = h.entries[0].anchor.snapshot.features.find((f) => f.id === h.featureId);
    exported.set(h.featureId, feature.shape.rings.filter((r) => r.kind === "hole").map((r) => r.vertexIds));
  });

  const { holeRings } = L.computeFaceHoles(faces, points);
  const pocketsOf = faces.map((f) => L.normalizeFaceRing(f.points, points).pockets);
  const asIds = (ring) => ring.map((v) => `v${v}`).sort().join(",");
  faces.forEach((f, i) => {
    const fid = `area-${String(i + 1).padStart(3, "0")}`;
    const fromExport = exported.get(fid).map((ids) => [...ids].sort().join(","));
    const expected = [...pocketsOf[i], ...holeRings[i]].map(asIds);
    assert.deepEqual(fromExport, expected, `${fid} の穴`);
  });
  // F（index 4）は G をポケットとして持つので、書き出しでは穴 1 つ・画面の穴は 0
  assert.equal(exported.get("area-005").length, 1);
  assert.equal(holeRings[4].length, 0);
});

test("gimoza.js の cleanRing / decomposeRing は logic.js のものと同じ結果を返す（一本化）", () => {
  const { points, faces } = buildScene();
  faces.forEach((f) => {
    assert.deepEqual(G.cleanRing(f.points), L.cleanRing(f.points));
    const ring = L.cleanRing(f.points);
    assert.deepEqual(G.decomposeRing(ring, points), L.decomposeRing(ring, points));
  });
});
