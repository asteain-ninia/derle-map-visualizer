// gimoza.js（.gimoza 書き出し）の純粋関数のユニットテスト。実行: `npm test`
"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const zlib = require("node:zlib");

// gimoza.js は幾何の部品を logic.js（DerleLogic）に頼るので、先に読み込む
require("../logic.js");
const G = require("../gimoza.js");

const P = (x, y) => ({ x, y });

test("crc32: 標準の検査値", () => {
  assert.equal(G.crc32(new TextEncoder().encode("123456789")), 0xcbf43926);
  assert.equal(G.crc32(new Uint8Array(0)), 0);
});

test("buildZip: 各エントリを読み戻すと元のバイト列になる", async () => {
  const enc = new TextEncoder();
  const files = {
    "manifest.json": enc.encode(JSON.stringify({ a: "あいう".repeat(200) })),
    "assets/basemap.img": new Uint8Array([1, 2, 3, 4, 5]),
  };
  const zip = Buffer.from(await G.buildZip(files));
  // 中央ディレクトリ終端の目印
  assert.equal(zip.readUInt32LE(zip.length - 22), 0x06054b50);
  assert.equal(zip.readUInt16LE(zip.length - 22 + 10), 2); // エントリ数

  let off = 0;
  const seen = {};
  while (zip.readUInt32LE(off) === 0x04034b50) {
    const method = zip.readUInt16LE(off + 8);
    const crc = zip.readUInt32LE(off + 14);
    const size = zip.readUInt32LE(off + 18);
    const nameLen = zip.readUInt16LE(off + 26);
    const extraLen = zip.readUInt16LE(off + 28);
    const name = zip.subarray(off + 30, off + 30 + nameLen).toString("utf8");
    const dataStart = off + 30 + nameLen + extraLen;
    const data = zip.subarray(dataStart, dataStart + size);
    const raw = method === 8 ? zlib.inflateRawSync(data) : data;
    assert.equal(G.crc32(raw), crc, `${name} の CRC`);
    seen[name] = Buffer.from(raw);
    off = dataStart + size;
  }
  assert.deepEqual(Object.keys(seen).sort(), Object.keys(files).sort());
  for (const name of Object.keys(files)) {
    assert.ok(seen[name].equals(Buffer.from(files[name])), `${name} の中身`);
  }
});

test("cleanRing: 閉じの重複と連続重複を除き、3 点未満は null", () => {
  assert.deepEqual(G.cleanRing([0, 1, 1, 2, 3, 0]), [0, 1, 2, 3]);
  assert.equal(G.cleanRing([0, 1, 0]), null);
});

test("decomposeRing: 逆回りの部分ループ（ポケット）を穴として分解する", () => {
  // 外周 0-1-2-3（時計回り、SVG 座標）の 0 から内側へスリットを入れて、逆回りの島 4-5-6-7 を回って戻る
  const pts = [
    P(0, 0), P(10, 0), P(10, 10), P(0, 10), // 外周
    P(4, 4), P(4, 6), P(6, 6), P(6, 4), // 島（外周と逆回り）
  ];
  const ring = [0, 1, 2, 3, 0, 4, 5, 6, 7, 4];
  const { outer, pockets } = G.decomposeRing(ring, pts);
  assert.equal(pockets.length, 1);
  assert.deepEqual(new Set(pockets[0]), new Set([4, 5, 6, 7]));
  assert.ok(outer.includes(1) && outer.includes(2) && outer.includes(3));
});
