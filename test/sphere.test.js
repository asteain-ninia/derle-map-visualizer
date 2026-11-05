// 球面投影（Wikipedia風出力の正射図法）のユニットテスト。実行: `npm test`（= `node --test`）
//
// 半球クリップは過去に 2 度、塗りが反転している（AGENT.md「球面投影の半球クリッピングで
// 地平線弧が欠落」「地平線弧の向き誤りによる塗りの反転」）。ここでは個々の式の値より、
// 「クリップ結果を円盤上で塗った領域が、球面上の真の領域を正射投影したものと一致する」ことを
// 主に確かめる。真値は円盤上の格子点を球面へ逆投影し、経緯度で領域の内外を判定して求める
// （clipPolygonToHemisphere とは独立した計算）。
"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const L = require("../logic.js");

const DEG = Math.PI / 180;

// 経緯度 1 度 = 地図 1 単位の地図枠（経度 -180〜180 → x 0〜360、緯度 90〜-90 → y 0〜180）。
const FRAME = { x: 0, y: 0, width: 360, height: 180 };
// index.html の getSeamTolerance と同じ式（max(0.1, 枠幅 × 1e-5)）。
const SEAM_TOL = Math.max(0.1, FRAME.width * 1e-5);
const toMap = ([lon, lat]) => ({ x: lon + 180, y: 90 - lat });

/** 経緯度（度）の頂点列を、刻まずにそのまま単位ベクトルへ写す。 */
function rawVectors(lonLats) {
  return lonLats.map((ll) => L.mapPointToVector(toMap(ll), FRAME, SEAM_TOL));
}

/** buildWikipediaCache と同じく、辺を経緯度で SPHERE_EDGE_MAX_DEG 以下に刻んでから写す。 */
function densifiedVectors(lonLats) {
  return L.densifyMapRing(lonLats.map(toMap), FRAME).map((p) => L.mapPointToVector(p, FRAME, SEAM_TOL));
}

/** 経度 -180〜179 を 1 度ずつ回る輪（3D で極を囲む。経度の巻き数 ±1）。lat は経度の関数。 */
function polarRing(latOf, reverse = false) {
  const out = [];
  for (let lon = -180; lon < 180; lon += 1) out.push(L.lonLatDegToVector(lon, latOf(lon)));
  return reverse ? out.reverse() : out;
}

const rotationAt = (lon, lat) => L.buildRotation(L.lonLatDegToVector(lon, lat));

/** 投影座標 (px, py) を、rotation で回した球面上の点へ逆投影し、経緯度（度）を返す。 */
function unprojectToLonLat(px, py, rot) {
  // projectRotated の逆: 回転後 x = py, 回転後 y = px, 見えている側なので z >= 0。
  const x2 = py;
  const y2 = px;
  const z2 = Math.sqrt(Math.max(0, 1 - px * px - py * py));
  // rotateVector の逆回転。
  const x1 = rot.sinLat * x2 + rot.cosLat * z2;
  const z1 = -rot.cosLat * x2 + rot.sinLat * z2;
  const y1 = y2;
  const x = rot.cosLon * x1 - rot.sinLon * y1;
  const y = rot.sinLon * x1 + rot.cosLon * y1;
  return { lon: Math.atan2(y, x) / DEG, lat: Math.asin(Math.max(-1, Math.min(1, z1))) / DEG };
}

/** 多角形の巻き数（SVG の既定の塗り規則 nonzero で、0 以外なら塗られる）。 */
function windingNumber(poly, px, py) {
  let w = 0;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i, i += 1) {
    const a = poly[j];
    const b = poly[i];
    const cross = (b.x - a.x) * (py - a.y) - (px - a.x) * (b.y - a.y);
    if (a.y <= py) {
      if (b.y > py && cross > 0) w += 1;
    } else if (b.y <= py && cross < 0) {
      w -= 1;
    }
  }
  return w;
}

/**
 * 円盤上の格子点で「クリップ結果が塗る」と「真の領域が見えている」を比べ、食い違う割合を返す。
 * truth(lon, lat) は球面上の真の領域（度）。
 */
function mismatchFraction(clipped, truth, rot, n = 120) {
  let total = 0;
  let bad = 0;
  for (let i = 0; i < n; i += 1) {
    for (let j = 0; j < n; j += 1) {
      const px = -1 + (2 * (i + 0.5)) / n;
      const py = -1 + (2 * (j + 0.5)) / n;
      if (px * px + py * py >= 1) continue;
      total += 1;
      const { lon, lat } = unprojectToLonLat(px, py, rot);
      const painted = clipped.length >= 3 && windingNumber(clipped, px, py) !== 0;
      if (painted !== truth(lon, lat)) bad += 1;
    }
  }
  return bad / total;
}

/** 真の領域が円盤に占める割合（何も塗らない結果との食い違い = 真の領域の割合）。 */
const truthFraction = (truth, rot) => mismatchFraction([], truth, rot);

const onUnitCircle = (p) => Math.abs(Math.hypot(p.x, p.y) - 1) < 1e-9;

/** クリップ結果のうち地平線（単位円）に沿って続く区間ごとの、回った角度（符号付き、ラジアン）。 */
function horizonRuns(clipped) {
  const runs = [];
  let current = 0;
  let open = false;
  for (let i = 0; i < clipped.length; i += 1) {
    const a = clipped[i];
    const b = clipped[(i + 1) % clipped.length];
    if (onUnitCircle(a) && onUnitCircle(b)) {
      current += L.normalizeAngle(Math.atan2(b.y, b.x) - Math.atan2(a.y, a.x));
      open = true;
    } else if (open) {
      runs.push(current);
      current = 0;
      open = false;
    }
  }
  if (open) runs.push(current);
  return runs;
}

function backPointInside(vectors, rot) {
  const back = L.getBackLonLat(rot);
  const indices = vectors.map((_, i) => i);
  const lonArr = vectors.map((v) => Math.atan2(v.y, v.x));
  const latArr = vectors.map((v) => Math.atan2(v.z, Math.hypot(v.x, v.y)));
  return L.containsBackPoint(indices, lonArr, latArr, back.lon, back.lat);
}

// 刻んだ多角形での許容差（円盤の格子点のうち食い違ってよい割合）。実測の最大は 0.06%。
// 残る差は地平線弧の 5 度刻み・格子の粗さ・辺を弦で結ぶことによるもの。
const MAX_MISMATCH = 0.003;

// AGENT.md「長い辺を刻む」の例: 経度 -100〜79°×緯度 -6〜66° の矩形を (-165°, 75°) から見る。
const AGENT_RECT = [[-100, -6], [79, -6], [79, 66], [-100, 66]];
const inAgentRect = (lon, lat) => lon >= -100 && lon <= 79 && lat >= -6 && lat <= 66;

// --- 回転と投影の約束 -------------------------------------------------------------

test("buildRotation / rotateVector: 視点は回転後 z = 1、投影の原点に来る", () => {
  for (const [lon, lat] of [[0, 0], [-165, 75], [120, -40], [180, 0], [0, 90], [33, -89.5]]) {
    const center = L.lonLatDegToVector(lon, lat);
    const rot = L.buildRotation(center);
    const r = L.rotateVector(center, rot);
    assert.ok(Math.abs(r.z - 1) < 1e-12, `(${lon}, ${lat}) の z = ${r.z}`);
    const p = L.projectRotated(r);
    assert.ok(Math.hypot(p.x, p.y) < 1e-7, `(${lon}, ${lat}) の投影 = (${p.x}, ${p.y})`);
  }
});

test("projectRotated: 北が上（SVG の y 負）・東が右（AGENT.md「Wikipedia風出力の方位反転」）", () => {
  const rot = rotationAt(20, 10);
  const north = L.projectRotated(L.rotateVector(L.lonLatDegToVector(20, 40), rot));
  assert.ok(north.y < -0.4 && Math.abs(north.x) < 1e-9, `北 = (${north.x}, ${north.y})`);
  const east = L.projectRotated(L.rotateVector(L.lonLatDegToVector(50, 10), rot));
  assert.ok(east.x > 0.4 && east.y < 0, `東 = (${east.x}, ${east.y})`); // 緯線は北へ弧を描くので y は少し負
  const pole = L.projectRotated(L.rotateVector({ x: 0, y: 0, z: 1 }, rotationAt(0, 0)));
  assert.ok(Math.abs(pole.x) < 1e-12 && Math.abs(pole.y + 1) < 1e-12, "赤道から見た北極は円盤の真上");
});

test("rotateVectorArrays は rotateVector / projectRotated と完全に同じ値を返す", () => {
  const vectors = [[0, 0], [10, 20], [-170, -80], [179.5, 45], [90, 0]].map(([a, b]) => L.lonLatDegToVector(a, b));
  const rot = rotationAt(-30, 25);
  const arr = L.rotateVectorArrays(vectors, rot);
  vectors.forEach((v, i) => {
    const r = L.rotateVector(v, rot);
    const p = L.projectRotated(r);
    assert.equal(arr.rotX[i], r.x);
    assert.equal(arr.rotY[i], r.y);
    assert.equal(arr.rotZ[i], r.z);
    assert.equal(arr.projX[i], p.x);
    assert.equal(arr.projY[i], p.y);
    assert.equal(arr.lonArr[i], Math.atan2(v.y, v.x));
    assert.equal(arr.latArr[i], Math.atan2(v.z, Math.hypot(v.x, v.y)));
  });
});

test("lonLatDegToVector / vectorToLonLat / getBackLonLat", () => {
  const v = L.lonLatDegToVector(-120, 35);
  const ll = L.vectorToLonLat(v);
  assert.ok(Math.abs(ll.lon / DEG + 120) < 1e-9 && Math.abs(ll.lat / DEG - 35) < 1e-9);
  const back = L.getBackLonLat(rotationAt(-120, 35));
  assert.ok(Math.abs(back.lon / DEG - 60) < 1e-9 && Math.abs(back.lat / DEG + 35) < 1e-9, "真裏は (60°, -35°)");
  assert.deepEqual(L.normalizeVector({ x: 0, y: 0, z: 0 }), { x: 1, y: 0, z: 0 });
});

test("slerpVectors: 両端と中間（大円上で等分）", () => {
  const a = L.lonLatDegToVector(0, 0);
  const b = L.lonLatDegToVector(90, 0);
  const mid = L.slerpVectors(a, b, 0.5);
  const ll = L.vectorToLonLat(mid);
  assert.ok(Math.abs(ll.lon / DEG - 45) < 1e-9 && Math.abs(ll.lat) < 1e-9);
  const same = L.slerpVectors(a, a, 0.3);
  assert.ok(Math.abs(same.x - 1) < 1e-12);
});

test("normalizeAngle: (-π, π] に収める", () => {
  assert.equal(L.normalizeAngle(Math.PI), Math.PI);
  assert.equal(L.normalizeAngle(-Math.PI), Math.PI);
  assert.ok(Math.abs(L.normalizeAngle(1.5 * Math.PI) + 0.5 * Math.PI) < 1e-12);
  assert.ok(Math.abs(L.normalizeAngle(-1.5 * Math.PI) - 0.5 * Math.PI) < 1e-12);
  assert.equal(L.normalizeAngle(0.25), 0.25);
});

// --- 地図座標 → 球面 -------------------------------------------------------------

test("mapPointToVector: 地図枠の経緯度、地図端（許容値以内）は経度 ±180 度ちょうど", () => {
  const frame = { x: 10, y: 20, width: 720, height: 360 };
  const v = L.mapPointToVector({ x: 10 + 360 + 60, y: 20 + 90 }, frame, 0.1);
  assert.deepEqual(v, L.lonLatDegToVector(30, 45));
  // 西端の 0.05 内側 → -180、東端の 1.5 外側（内側許容値の 20 倍以内）→ +180
  assert.deepEqual(L.mapPointToVector({ x: 10.05, y: 200 }, frame, 0.1), L.lonLatDegToVector(-180, 0));
  assert.deepEqual(L.mapPointToVector({ x: 731.5, y: 200 }, frame, 0.1), L.lonLatDegToVector(180, 0));
  // 緯度は ±90 度に丸める
  assert.deepEqual(L.mapPointToVector({ x: 370, y: 0 }, frame, 0.1), L.lonLatDegToVector(0, 90));
  assert.equal(L.mapPointToVector(null, frame, 0.1), null);
  assert.equal(L.mapPointToVector({ x: 0, y: 0 }, null, 0.1), null);
});

test("densifyMapSegment / densifyMapRing: 経緯度で SPHERE_EDGE_MAX_DEG 以下に刻む（両端は含まない）", () => {
  assert.equal(L.SPHERE_EDGE_MAX_DEG, 1);
  const mids = L.densifyMapSegment(toMap([-100, -6]), toMap([79, -6]), FRAME);
  assert.equal(mids.length, 178); // 179 度 → 179 片
  const chain = [toMap([-100, -6]), ...mids, toMap([79, -6])];
  for (let i = 1; i < chain.length; i += 1) {
    assert.ok(Math.abs(chain[i].x - chain[i - 1].x) <= 1 + 1e-9);
  }
  assert.deepEqual(L.densifyMapSegment(toMap([0, 0]), toMap([0.5, 0.5]), FRAME), []);
  const ring = L.densifyMapRing(AGENT_RECT.map(toMap), FRAME);
  assert.equal(ring.length, 2 * 179 + 2 * 72); // 閉じる辺（最後 → 最初）も刻む
  assert.deepEqual(ring[0], toMap(AGENT_RECT[0]));
});

// --- 地平線弧の部品 --------------------------------------------------------------

test("addHorizonArc: 単位円上を 5 度以下の刻みで回り、両端は含まない", () => {
  const out = [];
  L.addHorizonArc(out, 0, -200 * DEG);
  assert.ok(out.length >= 39);
  out.forEach((p) => assert.ok(onUnitCircle(p)));
  const angles = [0, ...out.map((p) => Math.atan2(p.y, p.x)), -200 * DEG];
  for (let i = 1; i < angles.length; i += 1) {
    const step = Math.abs(L.normalizeAngle(angles[i] - angles[i - 1]));
    assert.ok(step <= 5 * DEG + 1e-12, `刻み ${step / DEG} 度`);
  }
  const none = [];
  L.addHorizonArc(none, 1, 1e-9);
  assert.deepEqual(none, []);
});

test("buildHorizonArc: 不可視区間が基準半直線を横切った回数だけ 1 周ずつ足す（短い弧に限らない）", () => {
  const exit = { x: 1, y: 0 };
  const entry = { x: 0, y: 1 };
  // 横切らない: 短い方の弧（+90 度）
  const short = L.buildHorizonArc(exit, entry, 0);
  assert.ok(Math.abs(L.horizonArcSweep(short) - 90 * DEG) < 1e-12);
  // 1 回横切る: 反対回りの長い弧（-270 度）。旧実装は常に短い弧で結んで塗りが反転した。
  const longArc = L.buildHorizonArc(exit, entry, 1);
  assert.ok(Math.abs(L.horizonArcSweep(longArc) + 270 * DEG) < 1e-12);
  assert.equal(longArc.exitAngle, 0);
});

test("dedupePoints: 連続する重複と、最後 = 最初の閉じ点を落とす", () => {
  const P = (x, y) => ({ x, y });
  assert.deepEqual(
    L.dedupePoints([P(0, 0), P(0, 0), P(1, 0), P(1, 1e-9), P(1, 1), P(0, 0)]),
    [P(0, 0), P(1, 0), P(1, 1)],
  );
});

// --- 真裏の点の内外 --------------------------------------------------------------

test("containsLonLat: 経緯度平面の内外（経度差は短い方で扱う）", () => {
  const vectors = densifiedVectors(AGENT_RECT);
  const indices = vectors.map((_, i) => i);
  const lonArr = vectors.map((v) => Math.atan2(v.y, v.x));
  const latArr = vectors.map((v) => Math.atan2(v.z, Math.hypot(v.x, v.y)));
  assert.equal(L.containsLonLat(indices, lonArr, latArr, 0, 30 * DEG), true);
  assert.equal(L.containsLonLat(indices, lonArr, latArr, 100 * DEG, 30 * DEG), false);
  assert.equal(L.containsLonLat(indices, lonArr, latArr, 0, 70 * DEG), false);
  assert.equal(L.containsLonLat(indices, lonArr, latArr, -150 * DEG, 30 * DEG), false);
});

test("containsBackPoint: 極を囲む輪（巻き数 ±1）は面積の小さい側を内側とする（向きに依らない）", () => {
  const cases = [
    // [輪の緯度, 内側の点, 外側の点]
    [60, [10, 80], [10, 0]], // 北の球冠が小さい
    [-30, [10, -80], [10, 0]], // 南の球冠（緯度 -30 より南）が小さい
    [5, [10, 30], [10, -30]], // 赤道のすぐ北でも北側の方が小さい
  ];
  for (const [ringLat, inPt, outPt] of cases) {
    for (const reverse of [false, true]) {
      const vectors = polarRing(() => ringLat, reverse);
      const indices = vectors.map((_, i) => i);
      const lonArr = vectors.map((v) => Math.atan2(v.y, v.x));
      const latArr = vectors.map((v) => Math.atan2(v.z, Math.hypot(v.x, v.y)));
      const label = `緯度 ${ringLat} の輪${reverse ? "（逆回り）" : ""}`;
      assert.equal(L.containsBackPoint(indices, lonArr, latArr, inPt[0] * DEG, inPt[1] * DEG), true, `${label}: 内側`);
      assert.equal(L.containsBackPoint(indices, lonArr, latArr, outPt[0] * DEG, outPt[1] * DEG), false, `${label}: 外側`);
    }
  }
});

// --- 半球クリップ（clipPolygonToHemisphere / clipPolygonVectors）------------------

test("全体が見えている多角形は、頂点を投影したものがそのまま返る", () => {
  const rot = rotationAt(40, 20);
  const vectors = rawVectors([[30, 10], [50, 10], [50, 30], [30, 30]]);
  const expected = vectors.map((v) => L.projectRotated(L.rotateVector(v, rot)));
  assert.deepEqual(L.clipPolygonVectors(vectors, rot), expected);
});

test("全体が隠れている多角形は空", () => {
  const rot = rotationAt(40, 20);
  // 真裏 (-140°, -20°) のまわり
  assert.deepEqual(L.clipPolygonVectors(densifiedVectors([[-150, -30], [-130, -30], [-130, -10], [-150, -10]]), rot), []);
  assert.deepEqual(L.clipPolygonVectors([], rot), []);
});

test("地平線をまたぐ多角形: 切断点の間を地平線に沿う弧で埋める（弦で結ばない）", () => {
  // 回帰: 「球面投影の半球クリッピングで地平線弧が欠落」。出口と入口を直線で結ぶと弦と弧の間が欠ける。
  const rot = rotationAt(0, 0);
  const vectors = densifiedVectors([[60, -40], [120, -40], [120, 40], [60, 40]]);
  const clipped = L.clipPolygonVectors(vectors, rot);
  clipped.forEach((p) => assert.ok(Math.hypot(p.x, p.y) <= 1 + 1e-12, "円盤の外に出ない"));
  const onCircle = clipped.filter(onUnitCircle);
  assert.ok(onCircle.length >= 10, `地平線上の点 ${onCircle.length} 個`);
  // 地平線上で隣り合う点どうしは 5 度以内（弧の補間点が入っている）
  for (let i = 0; i < clipped.length; i += 1) {
    const a = clipped[i];
    const b = clipped[(i + 1) % clipped.length];
    if (onUnitCircle(a) && onUnitCircle(b)) {
      const step = Math.abs(L.normalizeAngle(Math.atan2(b.y, b.x) - Math.atan2(a.y, a.x)));
      assert.ok(step <= 5 * DEG + 1e-9, `地平線上の隣接点の間隔 ${step / DEG} 度`);
    }
  }
  const truth = (lon, lat) => lon >= 60 && lon <= 120 && lat >= -40 && lat <= 40;
  assert.ok(mismatchFraction(clipped, truth, rot) < MAX_MISMATCH);
});

test("AGENT.md の例（刻んだ場合）: 矩形を (-165°, 75°) から見た塗りが真値と一致する（向き 2 通り）", () => {
  const rot = rotationAt(-165, 75);
  for (const order of [AGENT_RECT, AGENT_RECT.slice().reverse()]) {
    const clipped = L.clipPolygonVectors(densifiedVectors(order), rot);
    const m = mismatchFraction(clipped, inAgentRect, rot);
    assert.ok(m < MAX_MISMATCH, `食い違い ${(m * 100).toFixed(2)}%`);
    // 不可視部分は地平線に沿ってほぼ半周（177 度）回り込む
    const runs = horizonRuns(clipped);
    assert.equal(runs.length, 1);
    assert.ok(Math.abs(Math.abs(runs[0]) - 177 * DEG) < 1 * DEG, `弧 ${runs[0] / DEG} 度`);
  }
  assert.ok(Math.abs(truthFraction(inAgentRect, rot) - 0.305) < 0.01, "真の領域は円盤の約 30%");
});

test("AGENT.md の例（頂点 4 つのまま）: 長い辺を刻まないと真値と大きく食い違う — 既知の制約", () => {
  // これは「正しい結果」ではない。弧の回転角は 3D の弦で、真裏の点の内外は経緯度平面の直線で見るため、
  // 長い辺があると両者が別の図形を見て塗りが反転する（AGENT.md「長い辺を刻む」: 円盤の 90% が食い違う）。
  // だから buildWikipediaCache は球面へ写す辺をすべて SPHERE_EDGE_MAX_DEG 以下に刻んでいる。
  // クリップ側を直して長い辺のままでも正しくなったなら、このテストと AGENT.md を書き換えること。
  const rot = rotationAt(-165, 75);
  for (const order of [AGENT_RECT, AGENT_RECT.slice().reverse()]) {
    const clipped = L.clipPolygonVectors(rawVectors(order), rot);
    const m = mismatchFraction(clipped, inAgentRect, rot);
    assert.ok(m > 0.8, `食い違い ${(m * 100).toFixed(1)}%（記録時 90.3%）`);
  }
});

test("真裏の点が多角形の内側: 弧の 1 本を逆回りにして内外を入れ替える", () => {
  // (170°, -30°) から見ると真裏 (-10°, 30°) が矩形の内側。見えているのは矩形の角 (79°, -6°) 付近だけ。
  const rot = rotationAt(170, -30);
  const vectors = densifiedVectors(AGENT_RECT);
  assert.equal(backPointInside(vectors, rot), true);
  const clipped = L.clipPolygonVectors(vectors, rot);
  assert.ok(clipped.length > 0);
  const m = mismatchFraction(clipped, inAgentRect, rot);
  assert.ok(m < MAX_MISMATCH, `食い違い ${(m * 100).toFixed(2)}%`);
  // 内外を入れ替えないと、見えている小さな角の補集合（円盤のほぼ全面）を塗ってしまう
  assert.ok(truthFraction(inAgentRect, rot) < 0.05);
});

test("不可視部分が地平線の 180 度以上に回り込む: 長い弧で結ぶ", () => {
  // 回帰: 「地平線弧の向き誤りによる塗りの反転」。経度 -170〜170°×緯度 0〜10° の帯を (-105°, -38°) から見ると、
  // 帯の隠れた部分は地平線に沿って 200 度回り込む。短い方の弧（160 度）で結ぶと塗りが反転する。
  const rot = rotationAt(-105, -38);
  const truth = (lon, lat) => lon >= -170 && lon <= 170 && lat >= 0 && lat <= 10;
  const clipped = L.clipPolygonVectors(densifiedVectors([[-170, 0], [170, 0], [170, 10], [-170, 10]]), rot);
  const runs = horizonRuns(clipped);
  assert.ok(runs.some((r) => Math.abs(r) > 190 * DEG), `地平線上の区間 ${runs.map((r) => (r / DEG).toFixed(1)).join(", ")} 度`);
  const m = mismatchFraction(clipped, truth, rot);
  assert.ok(m < MAX_MISMATCH, `食い違い ${(m * 100).toFixed(2)}%`);
});

test("極を囲む輪（巻き数 ±1）: 小さい側の球冠を塗る（向き 2 通り）", () => {
  // 追記「極を囲む輪」の回帰。経緯度平面で閉じないので containsLonLat では判定できない。
  const views = [[0, 0], [30, 50], [-120, -20], [60, 85], [10, -60]];
  for (const [ringLat, truth] of [[60, (lon, lat) => lat > 60], [-30, (lon, lat) => lat < -30]]) {
    for (const reverse of [false, true]) {
      const vectors = polarRing(() => ringLat, reverse);
      for (const [vlon, vlat] of views) {
        const rot = rotationAt(vlon, vlat);
        const m = mismatchFraction(L.clipPolygonVectors(vectors, rot), truth, rot);
        assert.ok(m < MAX_MISMATCH, `緯度 ${ringLat}${reverse ? "（逆回り）" : ""} を (${vlon}, ${vlat}) から: ${(m * 100).toFixed(2)}%`);
      }
    }
  }
});

test("無作為な視点 × いろいろな形で、塗りが真値と一致する", () => {
  let seed = 7;
  const rnd = () => {
    seed = (seed * 1103515245 + 12345) % 2147483648;
    return seed / 2147483648;
  };
  const views = [];
  for (let i = 0; i < 30; i += 1) views.push([rnd() * 360 - 180, Math.asin(rnd() * 2 - 1) / DEG]);

  const wavyLat = (lon) => 50 + 20 * Math.sin(2 * lon * DEG);
  const uShape = [[-150, -40], [150, -40], [150, 40], [-150, 40], [-150, 20], [130, 20], [130, -20], [-150, -20]];
  const shapes = [
    ["AGENT.md の矩形", densifiedVectors(AGENT_RECT), inAgentRect],
    ["AGENT.md の矩形（逆回り）", densifiedVectors(AGENT_RECT.slice().reverse()), inAgentRect],
    ["波打つ極の輪", polarRing(wavyLat), (lon, lat) => lat > wavyLat(lon)],
    ["南の球冠（逆回り）", polarRing(() => -30, true), (lon, lat) => lat < -30],
    ["細い帯", densifiedVectors([[-170, 0], [170, 0], [170, 10], [-170, 10]]),
      (lon, lat) => lon >= -170 && lon <= 170 && lat >= 0 && lat <= 10],
    ["コの字", densifiedVectors(uShape),
      (lon, lat) => lon >= -150 && lon <= 150 && lat >= -40 && lat <= 40 && !(lon < 130 && lat > -20 && lat < 20)],
  ];
  // 検査した場面の内訳（下の assert で、難しい場面を実際に通っていることを確かめる）
  let backInsideAndClipped = 0;
  let longRuns = 0;
  for (const [name, vectors, truth] of shapes) {
    for (const [vlon, vlat] of views) {
      const rot = rotationAt(vlon, vlat);
      const clipped = L.clipPolygonVectors(vectors, rot);
      const m = mismatchFraction(clipped, truth, rot, 100);
      assert.ok(m < MAX_MISMATCH, `${name} を (${vlon.toFixed(1)}, ${vlat.toFixed(1)}) から: ${(m * 100).toFixed(2)}%`);
      const runs = horizonRuns(clipped);
      if (runs.length && backPointInside(vectors, rot)) backInsideAndClipped += 1;
      if (runs.some((r) => Math.abs(r) > Math.PI)) longRuns += 1;
    }
  }
  assert.ok(backInsideAndClipped >= 10, `真裏の点が内側でクリップされた場面 ${backInsideAndClipped}`);
  assert.ok(longRuns >= 5, `地平線の 180 度超を回る場面 ${longRuns}`);
});

// --- 線の描画 --------------------------------------------------------------------

/** パス文字列の座標を (x, y) の組で取り出す。 */
function pathCoords(d) {
  const nums = d.match(/-?\d+\.\d+/g) || [];
  const out = [];
  for (let i = 0; i + 1 < nums.length; i += 2) out.push({ x: Number(nums[i]), y: Number(nums[i + 1]) });
  return out;
}

test("buildGraticulePath: 経緯線はすべて円盤の中、赤道から見た経度 0 度の経線は縦の直線", () => {
  for (const [lon, lat] of [[0, 0], [-165, 75], [100, -45]]) {
    const d = L.buildGraticulePath(rotationAt(lon, lat));
    const coords = pathCoords(d);
    assert.ok(coords.length > 100);
    coords.forEach((p) => assert.ok(Math.hypot(p.x, p.y) <= 1 + 1e-5, `(${p.x}, ${p.y})`));
  }
  const d = L.buildGraticulePath(rotationAt(0, 0));
  const meridian = d.split("M").find((part) => part.startsWith("0.00000 1.00000") || part.startsWith("0.00000 -1.00000"));
  assert.ok(meridian, "経度 0 度の経線が円盤の上端または下端から始まる");
  pathCoords(`M${meridian}`).forEach((p) => assert.equal(Math.abs(p.x), 0));
});

test("buildSegmentPath / buildVisiblePolyline: 隠れた側は地平線で切る", () => {
  const rot = rotationAt(0, 0);
  const a = L.lonLatDegToVector(60, 0);
  const b = L.lonLatDegToVector(120, 0);
  // 見える → 隠れる: 地平線（経度 90 度、投影 (1, 0)）で止まる
  assert.equal(L.buildSegmentPath([{ a, b }], rot), "M0.86603 0.00000 L1.00000 0.00000");
  assert.equal(L.buildSegmentPath([{ a: b, b: L.lonLatDegToVector(150, 0) }], rot), "");
  const line = [0, 30, 60, 120, 150, 210, 240, 300, 330].map((lon) => L.rotateVector(L.lonLatDegToVector(lon, 0), rot));
  const d = L.buildVisiblePolyline(line);
  assert.equal(d.split("M").filter(Boolean).length, 2, d); // 見える区間が 2 つ
  pathCoords(d).forEach((p) => assert.ok(Math.hypot(p.x, p.y) <= 1 + 1e-5));
});

test("buildPathFromPoints / formatCoord: 3 点未満は空、小数 5 桁で閉じたパス", () => {
  assert.equal(L.buildPathFromPoints([{ x: 0, y: 0 }, { x: 1, y: 0 }]), "");
  assert.equal(
    L.buildPathFromPoints([{ x: 0, y: 0 }, { x: 0.5, y: 0 }, { x: 0.5, y: -1 / 3 }]),
    "M0.00000 0.00000 L0.50000 0.00000 L0.50000 -0.33333 Z",
  );
});
