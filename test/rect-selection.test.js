// 矩形選択（地図上のドラッグで囲んだ面をまとめて選ぶ）の純粋関数のテスト。実行: `npm test`
//
// ポインタ操作・描画は index.html 側（DOM が要る）なので、ここでは判定と選択集合の更新だけを見る。
"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const L = require("../logic.js");

const P = (x, y) => ({ x, y });

/** 点列を足して、その範囲の頂点番号の列（リング）を返す */
function addRing(points, coords) {
  const start = points.length;
  coords.forEach(([x, y]) => points.push(P(x, y)));
  return coords.map((_, i) => start + i);
}

function square(points, x, y, size) {
  return addRing(points, [[x, y], [x + size, y], [x + size, y + size], [x, y + size]]);
}

test("rectFromCorners: どちら向きに引いても同じ矩形になる", () => {
  const expected = { minX: 1, minY: 2, maxX: 5, maxY: 8 };
  assert.deepEqual(L.rectFromCorners(P(1, 2), P(5, 8)), expected);
  assert.deepEqual(L.rectFromCorners(P(5, 8), P(1, 2)), expected);
  assert.deepEqual(L.rectFromCorners(P(5, 2), P(1, 8)), expected);
  assert.deepEqual(L.rectFromCorners(P(1, 8), P(5, 2)), expected);
});

test("ringInsideRect: 完全に収まる面だけが真（境界上は内側、はみ出し・かすめは偽）", () => {
  const points = [];
  const rect = { minX: 0, minY: 0, maxX: 10, maxY: 10 };
  assert.equal(L.ringInsideRect(square(points, 2, 2, 3), points, rect), true, "内側");
  assert.equal(L.ringInsideRect(square(points, 0, 0, 10), points, rect), true, "矩形と一致（境界上）");
  assert.equal(L.ringInsideRect(square(points, 8, 8, 3), points, rect), false, "一部はみ出す");
  assert.equal(L.ringInsideRect(square(points, 20, 20, 3), points, rect), false, "外側");
  assert.equal(L.ringInsideRect(square(points, -5, -5, 20), points, rect), false, "矩形を覆う大きな面");
  assert.equal(L.ringInsideRect([], points, rect), false, "空のリング");
});

test("ringInsideRect: 凹んだ面でも全頂点で判定する（外接矩形と同値）", () => {
  const points = [];
  // コの字形。くぼみの中だけを囲んでも選ばれない
  const ring = addRing(points, [[0, 0], [10, 0], [10, 10], [0, 10], [0, 7], [7, 7], [7, 3], [0, 3]]);
  assert.equal(L.ringInsideRect(ring, points, { minX: 1, minY: 4, maxX: 6, maxY: 6 }), false);
  assert.equal(L.ringInsideRect(ring, points, { minX: -1, minY: -1, maxX: 11, maxY: 11 }), true);
});

test("collectRectSelectionHits: 島の群れを囲むと島だけが選ばれ、周りの海や大陸は選ばれない", () => {
  const points = [];
  const features = [
    { id: 1, points: square(points, 0, 0, 100) }, // 海（島を取り巻く大きな面）
    { id: 2, points: square(points, -50, 40, 52) }, // 大陸の国（矩形を少しかすめる）
    { id: 3, points: square(points, 10, 10, 2) }, // 島
    { id: 4, points: square(points, 15, 12, 1) }, // 島
    { id: 5, points: square(points, 18, 18, 0.5) }, // 島
    { id: 6, points: square(points, 29, 10, 3) }, // 矩形の縁で切れる島
  ];
  const index = L.buildRectSelectionIndex(features, points);
  const rect = L.rectFromCorners(P(30, 25), P(1, 5));
  assert.deepEqual(L.collectRectSelectionHits(index, rect), ["3", "4", "5"]);
});

test("collectRectSelectionHits: 継ぎ目でつないだ地物はピースのどれかが収まれば選ばれる", () => {
  const points = [];
  const features = [
    { id: 7, rings: [square(points, 0, 10, 2), square(points, 998, 10, 2)] }, // 西端と東端のピース
    { id: 8, rings: [square(points, 500, 10, 2)] },
  ];
  const index = L.buildRectSelectionIndex(features, points);
  assert.deepEqual(index.map((e) => e.bounds.length), [2, 1]);
  assert.deepEqual(L.collectRectSelectionHits(index, L.rectFromCorners(P(-1, 9), P(3, 13))), ["7"]);
  assert.deepEqual(L.collectRectSelectionHits(index, L.rectFromCorners(P(990, 0), P(1001, 20))), ["7"]);
  assert.deepEqual(L.collectRectSelectionHits(index, L.rectFromCorners(P(100, 0), P(400, 20))), []);
});

test("rectSelectionMode: 修飾キーなし=置き換え、Shift=追加、Alt/Ctrl/⌘=解除（解除を優先）", () => {
  assert.equal(L.rectSelectionMode({}), "replace");
  assert.equal(L.rectSelectionMode(undefined), "replace");
  assert.equal(L.rectSelectionMode({ shiftKey: true }), "add");
  assert.equal(L.rectSelectionMode({ altKey: true }), "remove");
  assert.equal(L.rectSelectionMode({ ctrlKey: true }), "remove");
  assert.equal(L.rectSelectionMode({ metaKey: true }), "remove");
  assert.equal(L.rectSelectionMode({ shiftKey: true, altKey: true }), "remove");
});

test("applyRectSelection: 置き換え・追加（トグルしない）・解除。元の集合は変えない", () => {
  const current = new Set(["1", "2"]);
  assert.deepEqual([...L.applyRectSelection(current, ["2", "3"], "replace")].sort(), ["2", "3"]);
  assert.deepEqual([...L.applyRectSelection(current, [], "replace")], [], "何も囲まない置き換えは解除");
  assert.deepEqual([...L.applyRectSelection(current, ["2", "3"], "add")].sort(), ["1", "2", "3"], "選択済みの 2 は外れない");
  assert.deepEqual([...L.applyRectSelection(current, ["2", "3"], "remove")], ["1"]);
  assert.deepEqual([...current].sort(), ["1", "2"]);
});

test("exceedsDragThreshold: 閾値以下の小さな動きはクリックのまま", () => {
  assert.equal(L.exceedsDragThreshold(P(0, 0), P(3, 3), 5), false);
  assert.equal(L.exceedsDragThreshold(P(0, 0), P(5, 0), 5), false);
  assert.equal(L.exceedsDragThreshold(P(0, 0), P(4, 4), 5), true);
});

test("collectRectSelectionHits: 800 面・5 万点でも 1 回の判定は軽い（外接矩形の比較だけ）", () => {
  const points = [];
  const features = [];
  for (let i = 0; i < 800; i += 1) {
    const cx = (i % 40) * 50;
    const cy = Math.floor(i / 40) * 50;
    const coords = [];
    for (let k = 0; k < 62; k += 1) {
      const t = (k / 62) * Math.PI * 2;
      coords.push([cx + 20 * Math.cos(t), cy + 20 * Math.sin(t)]);
    }
    features.push({ id: i + 1, points: addRing(points, coords) });
  }
  assert.ok(points.length >= 49000);
  const index = L.buildRectSelectionIndex(features, points);
  const rect = L.rectFromCorners(P(-25, -25), P(225, 125));
  const started = process.hrtime.bigint();
  let hits;
  for (let n = 0; n < 1000; n += 1) hits = L.collectRectSelectionHits(index, rect);
  const perCallMs = Number(process.hrtime.bigint() - started) / 1e6 / 1000;
  assert.equal(hits.length, 5 * 3);
  assert.ok(perCallMs < 2, `1 回 ${perCallMs.toFixed(3)} ms`);
});
