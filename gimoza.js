/**
 * gimoza 書き出し — 抽出結果を地図アプリ (gimoza) が読める `.gimoza` にする。
 *
 * `.gimoza` は ZIP で、中身は次の 2 つ:
 *   manifest.json        … 地物・設定・ベースマップ参照を持つ本体 (UTF-8 JSON)
 *   assets/basemap.img   … 背景地図の資産。ラスター画像でもベクター JSON でもよい
 *
 * このファイルは「ZIP を書く」ことと「manifest を組み立てる」ことだけを担う。
 * ポリゴンの抽出・分類・領有判定は logic.js / index.html の責務。
 *
 * 依存なし・オフライン完結（要件定義書 §8）。圧縮は CompressionStream("deflate-raw")
 * があれば使い、無ければ無圧縮で書く（どちらも正しい ZIP で、読む側は区別しない）。
 */
(function (global) {
  "use strict";

  // --- CRC32 -----------------------------------------------------------------
  // ZIP の各エントリに必要。テーブルは初回だけ作る。
  let crcTable = null;
  function crc32(bytes) {
    if (crcTable === null) {
      crcTable = new Uint32Array(256);
      for (let n = 0; n < 256; n += 1) {
        let c = n;
        for (let k = 0; k < 8; k += 1) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
        crcTable[n] = c >>> 0;
      }
    }
    let c = 0xffffffff;
    for (let i = 0; i < bytes.length; i += 1) {
      c = crcTable[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
    }
    return (c ^ 0xffffffff) >>> 0;
  }

  /** deflate-raw で圧縮する。使えない環境では null を返す（呼び出し側が無圧縮へ）。 */
  async function deflateRaw(bytes) {
    if (typeof CompressionStream !== "function") return null;
    try {
      const cs = new CompressionStream("deflate-raw");
      const stream = new Blob([bytes]).stream().pipeThrough(cs);
      const buf = await new Response(stream).arrayBuffer();
      return new Uint8Array(buf);
    } catch (err) {
      return null;
    }
  }

  /**
   * ZIP を組み立てる。`files` は { パス: Uint8Array } 。
   *
   * 実装は ZIP の最小仕様どおり: 各エントリの「ローカルヘッダ + データ」を並べ、
   * 最後に「セントラルディレクトリ + 終端レコード」を置く。日時は 1980-01-01 固定
   * （出力を入力だけで決まるものにするため。同じ抽出結果からは常に同じバイト列が出る）。
   */
  async function buildZip(files) {
    const encoder = new TextEncoder();
    const entries = [];
    const chunks = [];
    let offset = 0;

    for (const path of Object.keys(files)) {
      const nameBytes = encoder.encode(path);
      const raw = files[path];
      const crc = crc32(raw);
      const deflated = await deflateRaw(raw);
      // 圧縮して小さくならなければ無圧縮で入れる
      const useDeflate = deflated !== null && deflated.length < raw.length;
      const data = useDeflate ? deflated : raw;
      const method = useDeflate ? 8 : 0;

      const local = new Uint8Array(30 + nameBytes.length);
      const lv = new DataView(local.buffer);
      lv.setUint32(0, 0x04034b50, true); // ローカルヘッダの目印
      lv.setUint16(4, 20, true); // 展開に必要な版
      lv.setUint16(6, 0, true); // フラグ
      lv.setUint16(8, method, true);
      lv.setUint16(10, 0, true); // 時刻
      lv.setUint16(12, 0x0021, true); // 日付 = 1980-01-01
      lv.setUint32(14, crc, true);
      lv.setUint32(18, data.length, true);
      lv.setUint32(22, raw.length, true);
      lv.setUint16(26, nameBytes.length, true);
      lv.setUint16(28, 0, true); // 拡張フィールド長
      local.set(nameBytes, 30);

      chunks.push(local, data);
      entries.push({ nameBytes, crc, method, compSize: data.length, rawSize: raw.length, offset });
      offset += local.length + data.length;
    }

    const centralStart = offset;
    for (const e of entries) {
      const c = new Uint8Array(46 + e.nameBytes.length);
      const cv = new DataView(c.buffer);
      cv.setUint32(0, 0x02014b50, true); // セントラルディレクトリの目印
      cv.setUint16(4, 20, true); // 作成した版
      cv.setUint16(6, 20, true); // 展開に必要な版
      cv.setUint16(8, 0, true);
      cv.setUint16(10, e.method, true);
      cv.setUint16(12, 0, true);
      cv.setUint16(14, 0x0021, true);
      cv.setUint32(16, e.crc, true);
      cv.setUint32(20, e.compSize, true);
      cv.setUint32(24, e.rawSize, true);
      cv.setUint16(28, e.nameBytes.length, true);
      cv.setUint16(30, 0, true);
      cv.setUint16(32, 0, true); // コメント長
      cv.setUint16(34, 0, true); // ディスク番号
      cv.setUint16(36, 0, true); // 内部属性
      cv.setUint32(38, 0, true); // 外部属性
      cv.setUint32(42, e.offset, true);
      c.set(e.nameBytes, 46);
      chunks.push(c);
      offset += c.length;
    }

    const end = new Uint8Array(22);
    const ev = new DataView(end.buffer);
    ev.setUint32(0, 0x06054b50, true); // 終端レコードの目印
    ev.setUint16(4, 0, true);
    ev.setUint16(6, 0, true);
    ev.setUint16(8, entries.length, true);
    ev.setUint16(10, entries.length, true);
    ev.setUint32(12, offset - centralStart, true);
    ev.setUint32(16, centralStart, true);
    ev.setUint16(20, 0, true);
    chunks.push(end);

    let total = 0;
    for (const c of chunks) total += c.length;
    const out = new Uint8Array(total);
    let pos = 0;
    for (const c of chunks) {
      out.set(c, pos);
      pos += c.length;
    }
    return out;
  }

  // ===========================================================================
  // ここから下は derle-import/tools/build-tmap.mjs (Node 用の変換スクリプト) と
  // derle-import/tools/extract-vector-basemap.mjs の移植。
  //
  // 移植にあたって外したもの (visualizer はビルド無し・オフライン完結のため):
  //   fs / path      … 入出力は引数と戻り値でやり取りする
  //   process.argv   … 設定は options で受け取る
  //   fflate         … ZIP は上の buildZip が書く
  //   polygon-clipping … 面の合併が必要だった「隙間の養子化」を、同じ平面分割上の
  //                      面隣接グラフの外側からの到達判定に置き換えた (下の
  //                      findEnclosedFaces のコメントに理由を書いてある)
  //
  // 包含解析は logic.js の DerleLogic.computeContainment を呼ぶ。Node 版は同じ判定を
  // 自前で持っていたが、二重実装は片方だけ直る事故のもとなので一本化した。
  // ===========================================================================

  // デーレ界の測地系: 惑星半径 1.0489 地球半径 (derle-climate PROJECT_CONTEXT.md)
  // → 赤道長 = 40075.017 × 1.0489 ≒ 42035 km。扁平率は未設定のため 0 (球面近似)。
  const EQUATOR_KM = 42035;

  // 面の種別ごとの塗り色。キーは index.html の classifyPolygons が返す種別名。
  const TYPE_STYLE = {
    海まで国境線: "#a8c98a",
    内陸国: "#dbc472",
    完全島国: "#d9a066",
    海岸線まで国境界: "#b79fd0",
    共有する海: "#8fb8dd",
    島国: "#89c9b8",
    不明: "#c0c0c0",
  };
  const NAME_PREFIX = { 共有する海: "海域", 島国: "島" };

  // 同君連合の親名ラベルは集約側で使うため、面への割当から除外する。
  const UNION_PARENT_NAMES = new Set(["ラヴィエンド・ウィル同君王国"]);

  // 点包含で割当できなかった洋上ラベルの手動割当 (ユーザー裁定 2026-08-08)。
  // キー = ラベル文言、値 = 割当先ポリゴンの内部点 (SVG 座標)。
  // ポリゴン id は再抽出のたびに変わるため、座標→点包含で解決する。
  const MANUAL_LABEL_ASSIGN = {
    ピニア: [1713.4, 953.3], // 隣の小島 (島国)
    テレン合議国: [290.3, 395.3], // 北西の大島
    バアカ: [2053.8, 1639.2], // 南方の巨大沿岸国
    "ヤシュー・アレファクト": [2430.8, 763.8], // 最寄りの大島 (洋上ラベル)
    ニィル: [3581.4, 1491.4], // 最寄りの大島 (洋上ラベル)
  };

  // 同君連合 (ユーザー裁定 2026-08-09): 構成国の内部点で国を特定し 1 集約へまとめる。
  const UNION_GROUPS = [
    {
      name: "ラヴィエンド・ウィル同君王国",
      childKindLabel: "構成国",
      members: [
        [509.2, 493.1], // ラヴィエン
        [610.4, 716.5], // ウィル
      ],
    },
  ];

  // 州名ラベルの無い子 (領海など) の手動命名: 面の内部点 (SVG 座標) で指定する。
  // 例: { at: [x, y], name: 'デーレ共和国領海', kindLabel: '領海' }
  const CHILD_NAME_ASSIGN = [];

  // 背景地図に含めるレイヤー (地形のみ。国境・州境・文字はアプリ側の地物なので除く)
  const BASEMAP_LAYER_IDS = [
    "_第一分高面",
    "_第二分高面",
    "_第三分高面",
    "_第四分高面",
    "_第五分高面",
  ];

  const pad = (n, w) => String(n).padStart(w, "0");

  // --- 1. リングの整形 -------------------------------------------------------
  // 閉包重複 (先頭=末尾) と連続重複を除去し、3 頂点未満は捨てる。
  function cleanRing(indices) {
    let r = indices.slice();
    if (r.length >= 2 && r[0] === r[r.length - 1]) r = r.slice(0, -1);
    const out = [];
    for (const i of r) {
      if (out.length === 0 || out[out.length - 1] !== i) out.push(i);
    }
    while (out.length >= 2 && out[0] === out[out.length - 1]) out.pop();
    return new Set(out).size >= 3 ? out : null;
  }

  function ringSignedArea(indices, points) {
    let a = 0;
    for (let i = 0; i < indices.length; i += 1) {
      const p = points[indices[i]];
      const q = points[indices[(i + 1) % indices.length]];
      a += p.x * q.y - q.x * p.y;
    }
    return a / 2;
  }

  /**
   * リングの正規化。
   *
   * 面抽出の一筆書きは「橋渡し辺を往復して島を除外する」自己接触リング (スリット付き)
   * を生むことがある。同一頂点の再訪で閉じる部分ループを検出し、
   *   - 面積 ≈ 0 (往復スパー) → 除去
   *   - 外周と逆回り (ポケット = 島などの除外領域) → 明示的な穴リングに分解
   *   - 外周と同回り (頂点接触のローブ) → そのまま残す
   */
  function decomposeRing(indices, points) {
    const outerSign = Math.sign(ringSignedArea(indices, points));
    let out = [];
    const pockets = [];
    let posOf = new Map();
    for (const v of indices) {
      if (posOf.has(v)) {
        const s = posOf.get(v);
        const loop = out.slice(s);
        const area = ringSignedArea(loop, points);
        if (Math.abs(area) < 1e-9 || new Set(loop).size < 3) {
          // 往復スパー: 捨てる (v の最初の出現だけ残す)
          out = out.slice(0, s + 1);
          posOf = new Map(out.map((u, k) => [u, k]));
          continue;
        }
        if (Math.sign(area) !== outerSign) {
          // ポケット: 穴リングとして分解
          pockets.push(loop);
          out = out.slice(0, s + 1);
          posOf = new Map(out.map((u, k) => [u, k]));
          continue;
        }
        // 同回りのローブ: 残す (頂点重複を許容)
      }
      posOf.set(v, out.length);
      out.push(v);
    }
    return { outer: out, pockets };
  }

  /**
   * 抽出結果の面 (`state.polygons`) を、整形済みの `{id, type, ring, pockets, absArea}`
   * に直す。`typeOf(i)` は面 index → 種別文字列。
   *
   * `rawRing` は整形前 (分解前) のリング。面どうしの隣接 (辺の共有) を数えるのに使う
   * ので、分解で形が変わる前のものを残しておく。
   */
  function normalizePolygons(rawPolygons, points, typeOf, log) {
    const out = [];
    let pocketCount = 0;
    let spurTrimmed = 0;
    for (let i = 0; i < rawPolygons.length; i += 1) {
      const p = rawPolygons[i];
      const ring = cleanRing(p.points);
      if (ring === null) {
        log(`退化した面を除外: id=${p.id} (${p.points.length} 点)`);
        continue;
      }
      const { outer, pockets } = decomposeRing(ring, points);
      const cleaned = cleanRing(outer);
      if (cleaned === null) {
        log(`分解後に退化した面を除外: id=${p.id}`);
        continue;
      }
      if (outer.length !== ring.length) spurTrimmed += 1;
      pocketCount += pockets.length;
      out.push({
        id: p.id,
        type: typeOf(i),
        ring: cleaned,
        rawRing: ring,
        pockets: pockets.map((l) => cleanRing(l)).filter((l) => l !== null),
        absArea: Math.abs(p.area),
      });
    }
    if (pocketCount > 0 || spurTrimmed > 0) {
      log(`リング正規化: ポケット穴 ${pocketCount} 件 / スパー除去 ${spurTrimmed} 面`);
    }
    return out;
  }

  // --- 2. チェーン単位の間引き -----------------------------------------------
  // グラフの次数 2 でない点 = 分岐点。分岐点間の道 (チェーン) を Douglas-Peucker で
  // 間引き、リング列を「残す点集合」でフィルタする。両隣の面が同じチェーンを参照する
  // ため、チェーン単位で間引けば境界の共有が壊れない。
  function simplifyKeepSet(epsilon, points, edges) {
    const NP = points.length;
    const adj = Array.from({ length: NP }, () => []);
    for (const e of edges) {
      adj[e.a].push(e.b);
      adj[e.b].push(e.a);
    }
    const degree = adj.map((l) => l.length);
    const keep = new Uint8Array(NP);

    const perpDist = (p, a, b) => {
      const dx = b.x - a.x;
      const dy = b.y - a.y;
      const len2 = dx * dx + dy * dy;
      if (len2 === 0) return Math.hypot(p.x - a.x, p.y - a.y);
      const t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / len2;
      const cx = a.x + t * dx;
      const cy = a.y + t * dy;
      return Math.hypot(p.x - cx, p.y - cy);
    };

    // chain: 点インデックス列 (両端含む)。両端は必ず keep、内部を DP。
    function dpChain(chain) {
      keep[chain[0]] = 1;
      keep[chain[chain.length - 1]] = 1;
      const stack = [[0, chain.length - 1]];
      while (stack.length > 0) {
        const [s, e] = stack.pop();
        if (e - s < 2) continue;
        const a = points[chain[s]];
        const b = points[chain[e]];
        let maxD = -1;
        let maxI = -1;
        for (let i = s + 1; i < e; i += 1) {
          const d = perpDist(points[chain[i]], a, b);
          if (d > maxD) {
            maxD = d;
            maxI = i;
          }
        }
        if (maxD > epsilon) {
          keep[chain[maxI]] = 1;
          stack.push([s, maxI], [maxI, e]);
        }
      }
    }

    const edgeKey = (a, b) => (a < b ? a * NP + b : b * NP + a);
    const visited = new Set();

    // 分岐点から伸びるチェーンを歩く
    for (let v = 0; v < NP; v += 1) {
      if (degree[v] === 2 || degree[v] === 0) continue;
      for (const first of adj[v]) {
        if (visited.has(edgeKey(v, first))) continue;
        const chain = [v];
        let prev = v;
        let cur = first;
        visited.add(edgeKey(prev, cur));
        chain.push(cur);
        while (degree[cur] === 2) {
          const next = adj[cur][0] === prev ? adj[cur][1] : adj[cur][0];
          visited.add(edgeKey(cur, next));
          prev = cur;
          cur = next;
          chain.push(cur);
        }
        dpChain(chain);
      }
    }
    // 残り = 分岐点を含まない純サイクル (孤島など)。2 分割して DP。
    for (let v = 0; v < NP; v += 1) {
      if (degree[v] !== 2) continue;
      if (visited.has(edgeKey(v, adj[v][0])) && visited.has(edgeKey(v, adj[v][1]))) continue;
      const cycle = [v];
      let prev = v;
      let cur = adj[v][0];
      visited.add(edgeKey(prev, cur));
      while (cur !== v) {
        cycle.push(cur);
        const next = adj[cur][0] === prev ? adj[cur][1] : adj[cur][0];
        visited.add(edgeKey(cur, next));
        prev = cur;
        cur = next;
      }
      // アンカー 2 点 (始点と最遠点) で半分ずつ DP
      let far = 1;
      let farD = -1;
      for (let i = 1; i < cycle.length; i += 1) {
        const d = Math.hypot(points[cycle[i]].x - points[v].x, points[cycle[i]].y - points[v].y);
        if (d > farD) {
          farD = d;
          far = i;
        }
      }
      dpChain(cycle.slice(0, far + 1));
      dpChain([...cycle.slice(far), v]);
      // 極小の島も最低 4 点は残す (退化リング防止)
      const keptInCycle = cycle.filter((i) => keep[i] === 1);
      if (keptInCycle.length < 4) {
        for (let k = 0; k < 4; k += 1) keep[cycle[Math.floor((k * cycle.length) / 4)]] = 1;
      }
    }
    return keep;
  }

  function bboxOfRing(ring, points) {
    let minX = Infinity;
    let minY = Infinity;
    let maxX = -Infinity;
    let maxY = -Infinity;
    for (const i of ring) {
      const q = points[i];
      if (q.x < minX) minX = q.x;
      if (q.y < minY) minY = q.y;
      if (q.x > maxX) maxX = q.x;
      if (q.y > maxY) maxY = q.y;
    }
    return { minX, minY, maxX, maxY };
  }

  function pointInRing(x, y, ring, points) {
    let inside = false;
    for (let i = 0, j = ring.length - 1; i < ring.length; j = i, i += 1) {
      const a = points[ring[i]];
      const b = points[ring[j]];
      if (a.y > y !== b.y > y && x < ((b.x - a.x) * (y - a.y)) / (b.y - a.y) + a.x) {
        inside = !inside;
      }
    }
    return inside;
  }

  /**
   * 各面を直接内包する最小面積の面の index (無ければ -1)。
   * 判定そのものは logic.js の computeContainment に任せる (同じ判定を二重に書かない)。
   */
  function computeParents(polys, points) {
    const logic = global.DerleLogic;
    if (!logic || typeof logic.computeContainment !== "function") {
      throw new Error("gimoza.js は logic.js (DerleLogic.computeContainment) を必要とします");
    }
    return logic.computeContainment(
      polys.map((p) => ({ points: p.ring, area: p.absArea })),
      points,
    );
  }

  /**
   * 冗長点の除去 (既定) / 目標頂点数への間引き (options.targetVertices 指定時)。
   *
   * visualizer は SVG パスを一定間隔でサンプリングするため、元データで 2 点だけの
   * 直線国境も等間隔の点列になる。小さい epsilon の間引きを掛けると直線は元どおり
   * 両端 2 点へ戻る。既定 0.02 地図単位 = 約 190 m で、最大ズーム (50 倍) でも画面上
   * 0.085 px 相当なので見た目は変わらない。
   */
  function applyCleanup(polys, points, edges, worldRect, opts, log) {
    const targetVertices = opts.targetVertices ?? null;
    const cleanupEpsilon = opts.cleanupEpsilon ?? 0.02;
    const keepFullIds = opts.keepFullIds instanceof Set ? opts.keepFullIds : new Set(opts.keepFullIds ?? []);
    const keepNearPairs = opts.keepNearPairs ?? [];
    if (targetVertices === null && !(cleanupEpsilon > 0)) return polys;

    let best = null;
    if (targetVertices !== null) {
      // epsilon を対数スケールで二分探索して目標頂点数へ寄せる
      let lo = 0.01;
      let hi = 20;
      for (let iter = 0; iter < 24; iter += 1) {
        const mid = Math.sqrt(lo * hi);
        const keep = simplifyKeepSet(mid, points, edges);
        const used = new Set();
        for (const p of polys) for (const i of p.ring) if (keep[i] === 1) used.add(i);
        const count = used.size;
        if (best === null || Math.abs(count - targetVertices) < Math.abs(best.count - targetVertices)) {
          best = { epsilon: mid, keep, count };
        }
        if (Math.abs(count - targetVertices) / targetVertices < 0.03) break;
        if (count > targetVertices) lo = mid;
        else hi = mid;
      }
      log(`間引き: epsilon=${best.epsilon.toFixed(4)} → ${best.count} 共有頂点`);
    } else {
      const keep = simplifyKeepSet(cleanupEpsilon, points, edges);
      const used = new Set();
      const before = new Set();
      for (const p of polys) {
        for (const i of p.ring) {
          before.add(i);
          if (keep[i] === 1) used.add(i);
        }
      }
      best = { epsilon: cleanupEpsilon, keep, count: used.size };
      log(`冗長点の除去: epsilon=${cleanupEpsilon} → ${before.size} → ${used.size} 共有頂点`);
    }

    // 保護 (全周リング/ポケット/内包の原型維持) は「積極的な間引き」でだけ必要。
    // 冗長点の除去は形をほぼ動かさないので保護を掛けず、代わりに生成後の整合性
    // チェックで担保する (崩れたら keepFullIds で個別救済)。
    const aggressive = targetVertices !== null;

    // 経度全周に及ぶリング (南極域など bbox が世界幅の 9 割超) はシーム跨ぎ展開の解釈が
    // 間引きで変わり偽の重なりを生むため、軽い除去でも常に間引きを解除する。
    for (const p of polys) {
      const bb = bboxOfRing(p.ring, points);
      const wraps = bb.maxX - bb.minX >= worldRect.width * 0.9;
      if (wraps || keepFullIds.has(p.id)) {
        for (const i of p.ring) best.keep[i] = 1;
        log(`keep-full: 面 id=${p.id} は間引きを解除${wraps ? " (全周リング)" : ""}`);
      }
    }

    if (aggressive) {
      // ポケット穴は外周と独立に間引かれると親からはみ出すため、ポケットの頂点と
      // その近傍の外周頂点は間引かない (ポケットは小さいので費用は僅か)。
      let pocketKept = 0;
      for (const p of polys) {
        for (const loop of p.pockets) {
          const bb = bboxOfRing(loop, points);
          for (const i of loop) {
            if (best.keep[i] !== 1) {
              best.keep[i] = 1;
              pocketKept += 1;
            }
          }
          const margin = best.epsilon * 6;
          for (const i of p.ring) {
            const q = points[i];
            if (
              q.x >= bb.minX - margin &&
              q.x <= bb.maxX + margin &&
              q.y >= bb.minY - margin &&
              q.y <= bb.maxY + margin &&
              best.keep[i] !== 1
            ) {
              best.keep[i] = 1;
              pocketKept += 1;
            }
          }
        }
      }
      if (pocketKept > 0) log(`ポケット穴の保持: ${pocketKept} 点を間引き対象から除外`);

      // 内包関係の穴 (島・包領・内海) も、外周と独立に間引かれると親からはみ出す。
      // 内包ペアを間引き前に求め、内側と、その bbox 近傍の外周頂点を保持する。
      const preParents = computeParents(polys, points);
      // keepFullIds で指定した面は、その内包子 (穴になる面) も併せて原型維持する。
      for (let i = 0; i < polys.length; i += 1) {
        const parent = preParents[i];
        if (parent !== -1 && keepFullIds.has(polys[parent].id)) {
          for (const vi of polys[i].ring) best.keep[vi] = 1;
          for (const loop of polys[i].pockets) for (const vi of loop) best.keep[vi] = 1;
        }
      }
      let containKept = 0;
      for (let i = 0; i < polys.length; i += 1) {
        const parent = preParents[i];
        if (parent === -1) continue;
        const bb = bboxOfRing(polys[i].ring, points);
        // 内側リングも間引かない: 内外の隙間が epsilon 未満だと、内側だけが動いても
        // 外周を跨いでしまう (穴と親の距離は保証されていない)。
        for (const vi of polys[i].ring) {
          if (best.keep[vi] !== 1) {
            best.keep[vi] = 1;
            containKept += 1;
          }
        }
        const margin = best.epsilon * 6;
        for (const vi of polys[parent].ring) {
          const q = points[vi];
          if (
            q.x >= bb.minX - margin &&
            q.x <= bb.maxX + margin &&
            q.y >= bb.minY - margin &&
            q.y <= bb.maxY + margin &&
            best.keep[vi] !== 1
          ) {
            best.keep[vi] = 1;
            containKept += 1;
          }
        }
      }
      if (containKept > 0) log(`内包穴の周辺保持: ${containKept} 点を間引き対象から除外`);
    }

    // keep-near: 面 A の境界点のうち、面 B の bbox 近傍にあるものだけ間引きを解除する
    // (A の簡略化境界が B を横切った場合の局所修復)。
    for (const [aId, bId] of keepNearPairs) {
      const a = polys.find((p) => p.id === aId);
      const b = polys.find((p) => p.id === bId);
      if (!a || !b) continue;
      const bb = bboxOfRing(b.ring, points);
      const margin = best.epsilon * 4;
      let restored = 0;
      for (const i of a.ring) {
        const q = points[i];
        if (
          q.x >= bb.minX - margin &&
          q.x <= bb.maxX + margin &&
          q.y >= bb.minY - margin &&
          q.y <= bb.maxY + margin &&
          best.keep[i] !== 1
        ) {
          best.keep[i] = 1;
          restored += 1;
        }
      }
      log(`keep-near: 面 ${aId} の ${restored} 点を ${bId} 近傍で復元`);
    }

    const out = [];
    for (const p of polys) {
      const ring = cleanRing(p.ring.filter((i) => best.keep[i] === 1));
      if (ring === null) {
        log(`間引きで退化したため原型を維持: 面 id=${p.id}`);
        out.push(p);
      } else {
        const pockets = p.pockets
          .map((l) => cleanRing(l.filter((i) => best.keep[i] === 1)))
          .filter((l) => l !== null);
        out.push({ ...p, ring, pockets });
      }
    }
    return out;
  }

  // --- 4. ラベル -------------------------------------------------------------
  /** `<style>` の中の `.cls-N { ... }` から、クラス名 → 宣言値の対応を作る。 */
  function cssDeclMap(svgDoc, property) {
    const map = new Map();
    const styles = svgDoc.getElementsByTagName("style");
    const ruleRe = /\.([\w-]+(?:\s*,\s*\.[\w-]+)*)\s*\{([^}]*)\}/g;
    const declRe = new RegExp(`(?:^|[;{\\s])${property}:\\s*([^;]+)`);
    for (let s = 0; s < styles.length; s += 1) {
      const css = styles[s].textContent || "";
      for (let m = ruleRe.exec(css); m !== null; m = ruleRe.exec(css)) {
        const d = declRe.exec(m[2]);
        if (d === null) continue;
        // 同じクラスに複数回書かれていたら後勝ち (CSS の詳細度は全て等しいので順序勝負)
        for (const cls of m[1].split(",")) map.set(cls.trim().replace(/^\./, ""), d[1].trim());
      }
    }
    return map;
  }

  /**
   * ラベルレイヤー (`_国名` / `_州名` / `_地名` / `_海名` / `_山岳名`) の `<text>` から
   * 位置・回転角・フォントサイズ・文字列を取り出す。
   *
   * Node 版は生の SVG テキストを正規表現で読んでいたが、こちらは DOMParser で読んだ
   * ドキュメントを直接たどる (同じ属性を見ているので結果は同じ)。フォントサイズだけは
   * 切り離された文書で getComputedStyle が効かないため、`<style>` を正規表現で読む。
   */
  function extractLabels(svgDoc, layerId) {
    if (!svgDoc) return [];
    const layer =
      svgDoc.getElementById(layerId) || svgDoc.querySelector(`[data-name='${layerId.replace(/^_/, "")}']`);
    if (layer === null) return [];
    const fontSizes = cssDeclMap(svgDoc, "font-size");
    const labels = [];
    const texts = layer.getElementsByTagName("text");
    for (let t = 0; t < texts.length; t += 1) {
      const el = texts[t];
      const transform = el.getAttribute("transform") || "";
      const tr = transform.match(/translate\(\s*([-\d.]+)[\s,]+([-\d.]+)\s*\)/);
      if (!tr) continue;
      const rot = transform.match(/rotate\(\s*([-\d.]+)/);
      const cls = el.getAttribute("class");
      const name = (el.textContent || "").replace(/\s+/g, "").trim();
      if (name.length === 0) continue;
      // 1 行目の文字数 (CJK は正方形グリフなので 1 文字 ≒ font-size として長さを見積る)
      const tspan = el.getElementsByTagName("tspan")[0];
      const firstLine = tspan ? (tspan.textContent || "").trim() : "";
      const fsRaw = cls ? fontSizes.get(cls.split(/\s+/)[0]) : undefined;
      const fs = fsRaw ? Number(String(fsRaw).replace(/px\s*$/, "")) : NaN;
      labels.push({
        x: Number(tr[1]),
        y: Number(tr[2]),
        angleDeg: rot ? Number(rot[1]) : 0,
        fontSize: Number.isFinite(fs) && fs > 0 ? fs : 12,
        firstLineLen: firstLine.length || name.length,
        name,
      });
    }
    return labels;
  }

  /** ラベル文字列の中央 (グリフ中央まで持ち上げた位置)。 */
  function labelCenter(label) {
    const rad = (label.angleDeg * Math.PI) / 180;
    const len = label.firstLineLen * label.fontSize;
    const lift = label.fontSize * 0.35;
    return [
      label.x + (Math.cos(rad) * len) / 2 + Math.sin(rad) * lift,
      label.y + (Math.sin(rad) * len) / 2 - Math.cos(rad) * lift,
    ];
  }

  // --- 5. 座標変換 -----------------------------------------------------------
  /**
   * 世界矩形 = SVG の海色レクト (`_海色`)。経度 -180..180 / 緯度 90..-90 に対応する。
   * Node 版はこれを world-rect.json に固定値で持っていたが、SVG から直接読む。
   */
  function worldRectFromSvg(svgDoc) {
    const rect = svgDoc.getElementById("_海色") || svgDoc.querySelector("[data-name='海色']");
    if (rect === null) throw new Error("世界矩形の基準 (_海色 レクト) が SVG にありません");
    const num = (attr, fallback) => {
      const v = Number(rect.getAttribute(attr));
      return Number.isFinite(v) ? v : fallback;
    };
    return { x: num("x", 0), y: num("y", 0), width: num("width", 0), height: num("height", 0) };
  }

  /**
   * SVG 座標 → 経緯度。作図が海色レクトをわずかにはみ出す点があり、経度 +180 を超えると
   * 周回して反対端の地物と極小の重なりを作るため、±180/±90 にクランプする。
   */
  function makeToLonLat(points, worldRect) {
    return (i) => {
      const p = points[i];
      const lon = ((p.x - worldRect.x) / worldRect.width) * 360 - 180;
      const lat = 90 - ((p.y - worldRect.y) / worldRect.height) * 180;
      return {
        lon: Math.min(180, Math.max(-180, Math.round(lon * 1e7) / 1e7)),
        lat: Math.min(90, Math.max(-90, Math.round(lat * 1e7) / 1e7)),
      };
    };
  }

  // --- 6a. 面の隣接と「囲まれた面」の検出 -------------------------------------
  /**
   * 辺 (頂点インデックスの組) → その辺を持つ面 index の一覧。
   * 面はすべて同じ平面分割から取れているので、内部の辺はちょうど 2 面が共有し、
   * 外界に面した辺は 1 面しか持たない。
   */
  function buildFaceAdjacency(polys, pointCount) {
    const edgeToFaces = new Map();
    for (let i = 0; i < polys.length; i += 1) {
      const ring = polys[i].rawRing || polys[i].ring;
      for (let k = 0; k < ring.length; k += 1) {
        const a = ring[k];
        const b = ring[(k + 1) % ring.length];
        if (a === b) continue;
        const key = a < b ? a * pointCount + b : b * pointCount + a;
        const list = edgeToFaces.get(key);
        if (list === undefined) edgeToFaces.set(key, [i]);
        else if (!list.includes(i)) list.push(i);
      }
    }
    // 面 index → 隣接面 index の一覧 / 外界に面しているか
    const neighbors = Array.from({ length: polys.length }, () => new Set());
    const touchesOutside = new Uint8Array(polys.length);
    for (const list of edgeToFaces.values()) {
      if (list.length === 1) {
        touchesOutside[list[0]] = 1;
        continue;
      }
      for (const a of list) for (const b of list) if (a !== b) neighbors[a].add(b);
    }
    return { neighbors, touchesOutside };
  }

  /**
   * `blocked` の面をすべて壁にしたとき、外界から到達できない面の集合を返す
   * (= blocked の面々にぐるりと囲まれた隙間)。
   *
   * Node 版はここで polygon-clipping による面の合併を行い、「合併の穴の中に内部点が
   * 落ちるか」で判定していた。polygon-clipping は使えないので、同じ平面分割上の
   * 面隣接グラフを外側から塗りつぶす方法に置き換えた。面はすべて同一のグラフから
   * 取れた分割なので「合併の穴」と「外側から到達できない面」は一致する。
   *
   * Node 版と意味が変わる点: Node 版は内部点 1 個の位置で判定していたため、判定点の
   * 取り方 (凹形の重心など) 次第で結果がぶれ得たが、こちらは辺の共有だけで決まるので
   * ぶれない。また Node 版は 3 巡まで繰り返して段階的に養子を増やしていたが、
   * 到達判定は 1 回で不動点になる (養子候補は元から外界に到達できないので、
   * それを壁に加えても到達可能性は変わらない)。
   */
  function findEnclosedFaces(adjacency, blocked, faceCount) {
    const visited = new Uint8Array(faceCount);
    const stack = [];
    for (let i = 0; i < faceCount; i += 1) {
      if (blocked.has(i)) continue;
      if (adjacency.touchesOutside[i] === 1) {
        visited[i] = 1;
        stack.push(i);
      }
    }
    while (stack.length > 0) {
      const v = stack.pop();
      for (const w of adjacency.neighbors[v]) {
        if (visited[w] === 1 || blocked.has(w)) continue;
        visited[w] = 1;
        stack.push(w);
      }
    }
    const enclosed = [];
    for (let i = 0; i < faceCount; i += 1) {
      if (!blocked.has(i) && visited[i] === 0) enclosed.push(i);
    }
    return enclosed;
  }

  /** 子地物の塗り分け: 親色を明暗に振る。 */
  function shadeColor(hex, step) {
    const n = parseInt(hex.slice(1), 16);
    const f = (c) => Math.max(0, Math.min(255, Math.round(c + step * 14)));
    const r = f((n >> 16) & 255);
    const g = f((n >> 8) & 255);
    const b = f(n & 255);
    return `#${((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1)}`;
  }

  /** `types` (Map / 配列 / 関数 / 面に埋まった type) を index → 種別名 の関数に均す。 */
  function makeTypeOf(types, polygons) {
    if (typeof types === "function") return (i) => types(i) || "不明";
    if (types instanceof Map) return (i) => types.get(String(polygons[i].id)) || "不明";
    if (Array.isArray(types)) return (i) => types[i] || "不明";
    return (i) => polygons[i].type || "不明";
  }

  // --- 6. manifest の組み立て -------------------------------------------------
  /**
   * 抽出結果から manifest を組み立てる。
   *
   * input:
   *   polygons        … state.polygons ([{id, area, points}])
   *   points          … state.graph.points ([{x,y}])
   *   edges           … state.graph.edges ([{a,b}])。間引きに使う
   *   types           … state.polygonTypes (Map) など。面 index → 種別名 に均される
   *   svgDoc          … DOMParser で読んだ元 SVG (ラベルと世界矩形の取得に使う)
   *   countryPolygons … 階層化用の国レベル抽出 {points, polygons, types}。省略時はフラット
   *   options         … cleanupEpsilon (既定 0.02) / projectName / equatorKm / log ほか
   */
  function buildManifest(input) {
    const opts = input.options || {};
    const log = typeof opts.log === "function" ? opts.log : () => {};
    const points = input.points;
    const svgDoc = input.svgDoc;
    const worldRect = opts.worldRect || worldRectFromSvg(svgDoc);
    const equatorKm = opts.equatorKm ?? EQUATOR_KM;

    // 1. リング整形
    const typeOf = makeTypeOf(input.types, input.polygons);
    const normalized = normalizePolygons(input.polygons, points, typeOf, log);

    // 2. 冗長点の除去 / 間引き
    const activePolygons = applyCleanup(normalized, points, input.edges || [], worldRect, opts, log);

    // 3. 包含解析 → 穴リング
    const bboxes = activePolygons.map((p) => bboxOfRing(p.ring, points));
    const parentOf = computeParents(activePolygons, points);
    log(`包含解析: ${parentOf.filter((v) => v !== -1).length} 面が他の面に内包 → 穴リング化`);

    const toLonLat = makeToLonLat(points, worldRect);

    // --- 世界の継ぎ目をまたぐ長い辺の分割 ---------------------------------------
    //
    // 読み込む側 (アプリ) は経度を「連続量」として復元する。リングを順に辿り、
    // 隣り合う頂点の経度差が 180 度を超えていたら「反対回りの方が近い」と解釈して
    // ±360 度ずらす (src/domain/geo/longitude.ts の unwrapLonSequence)。つまり
    // **辺 1 本で経度 180 度以上を渡ることは表現できない**。
    //
    // 南極域のように世界を一周する帯は、外周が「海岸線を西から東へ → 世界の右端で
    // 南極点まで下る → 南緯 90 度の線を東から西へ戻る → 左端で戻り上がる」という形に
    // なる。この「戻る」1 本が経度 360 度ぶんの辺になる。
    //
    // 旧方式 (getPointAtLength の等長サンプリング) では、この直線も間引き間隔ごとに
    // 刻まれていたので 3,500 点ほどに分かれており、たまたま条件を満たしていた。
    // 新方式 (d 属性の直接解釈) は直線コマンドの端点をそのまま使うため、両端 2 点だけの
    // 1 本の辺になる。すると復元側はこれを「経度 0.006 度の微小な辺」と読み違え、
    // 以降の頂点が丸ごと 1 周ぶん別の枝へ飛ぶ。結果、帯の外周が本来と違う範囲を覆い、
    // 帯の中にある島や領域と重なって整合性チェックが落ちる (実測 10 面 × 2 種 = 20 件)。
    //
    // 点列化の粗さに関係なく成り立たせるべき不変条件なので、書き出しの直前に
    // 「1 辺あたりの経度幅」を上限以下へ割る。地図は正距円筒なので辺は地図平面上の
    // 直線であり、経度・緯度を線形に内挿すれば元の直線そのものになる。
    const MAX_LON_SPAN_DEG = 60;
    /** 分割で増えた頂点の座標。id → { lon, lat }。 */
    const seamExtras = new Map();
    const roundCoord = (lon, lat) => ({
      lon: Math.min(180, Math.max(-180, Math.round(lon * 1e7) / 1e7)),
      lat: Math.min(90, Math.max(-90, Math.round(lat * 1e7) / 1e7)),
    });
    /**
     * リング (点 index の並び) → 頂点 id の並び。経度幅が上限を超える辺には中間点を打つ。
     * 追加頂点の id は辺の両端の点 index から決まり、向きに依らないので、同じ辺を
     * 共有する 2 つの面には同じ id が入る (境界は共有されたまま = 頂点の食い違いを作らない)。
     */
    const ringVertexIds = (ring) => {
      const out = [];
      for (let k = 0; k < ring.length; k += 1) {
        const a = ring[k];
        const b = ring[(k + 1) % ring.length];
        out.push(`v${a}`);
        const ca = toLonLat(a);
        const cb = toLonLat(b);
        const span = Math.abs(cb.lon - ca.lon);
        if (!(span > MAX_LON_SPAN_DEG)) continue;
        const parts = Math.ceil(span / MAX_LON_SPAN_DEG);
        const lo = Math.min(a, b);
        const hi = Math.max(a, b);
        for (let s = 1; s < parts; s += 1) {
          const t = s / parts;
          // id は常に (小さい index → 大きい index) 向きの位置で付ける
          const slot = a === lo ? s : parts - s;
          const id = `v${lo}-${hi}-${parts}-${slot}`;
          if (!seamExtras.has(id)) {
            seamExtras.set(id, roundCoord(ca.lon + (cb.lon - ca.lon) * t, ca.lat + (cb.lat - ca.lat) * t));
          }
          out.push(id);
        }
      }
      return out;
    };
    /** 頂点 id → 経緯度 (元の点はその場で変換、追加頂点は上の表から)。 */
    const coordOfId = (id) => seamExtras.get(id) ?? toLonLat(Number(id.slice(1)));
    /** 頂点 id の並び順 (元の点は index 昇順、追加頂点はその辺の直後)。 */
    const idSortKey = (id) => {
      const m = /^v(\d+)(?:-(\d+)-(\d+)-(\d+))?$/.exec(id);
      return m === null ? [Infinity, 0] : [Number(m[1]), m[4] === undefined ? 0 : Number(m[4])];
    };
    const sortIds = (ids) =>
      [...ids].sort((x, y) => {
        const kx = idSortKey(x);
        const ky = idSortKey(y);
        return kx[0] - ky[0] || kx[1] - ky[1];
      });

    const smallestContaining = (x, y) => {
      let bestI = -1;
      for (let i = 0; i < activePolygons.length; i += 1) {
        const bb = bboxes[i];
        if (x < bb.minX || x > bb.maxX || y < bb.minY || y > bb.maxY) continue;
        if (!pointInRing(x, y, activePolygons[i].ring, points)) continue;
        if (bestI === -1 || activePolygons[i].absArea < activePolygons[bestI].absArea) bestI = i;
      }
      return bestI;
    };

    // 4. 国名ラベルの割り当て
    const countryLabels = extractLabels(svgDoc, "_国名");
    log(`国名ラベル: ${countryLabels.length} 件`);
    const nameOf = new Array(activePolygons.length).fill(null);
    const unmatched = [];
    const hierarchical = input.countryPolygons != null;
    for (const label of countryLabels) {
      if (hierarchical && UNION_PARENT_NAMES.has(label.name)) continue;
      // テキストはアンカーから rotate 方向へ「文字数 × font-size」だけ伸び、グリフは
      // ベースラインの上に乗る。文字列に沿って密にサンプルし、いずれかの点を含む面の
      // うち最小面積のものを選ぶ (極小の飛び国の上を文字列が横切る場合、周囲の大国で
      // なく飛び国に割り当てるため)。
      const rad = (label.angleDeg * Math.PI) / 180;
      const len = label.firstLineLen * label.fontSize;
      const lift = label.fontSize * 0.35; // ベースラインからグリフ中央への持ち上げ
      const atLift = (t, liftScale) => [
        label.x + Math.cos(rad) * len * t + Math.sin(rad) * lift * liftScale,
        label.y + Math.sin(rad) * len * t - Math.cos(rad) * lift * liftScale,
      ];
      const at = (t) => atLift(t, 1);
      const dense = [];
      for (let t = 0; t <= 1.0001; t += 0.1) dense.push(atLift(t, 1), atLift(t, 0));
      // 1) マイクロ国家 (面積 < 100): 文字列がその上を通っていれば最優先で採用
      const MICRO_AREA = 100;
      let bestI = -1;
      for (const [x, y] of dense) {
        const i = smallestContaining(x, y);
        if (
          i !== -1 &&
          activePolygons[i].absArea < MICRO_AREA &&
          (bestI === -1 || activePolygons[i].absArea < activePolygons[bestI].absArea)
        ) {
          bestI = i;
        }
      }
      // 2) 通常国: 文字列の中央 → 1/4・3/4 点 → アンカー の先着順
      if (bestI === -1) {
        for (const [x, y] of [at(0.5), at(0.25), at(0.75), [label.x, label.y], at(1.0)]) {
          bestI = smallestContaining(x, y);
          if (bestI !== -1) break;
        }
      }
      if (bestI === -1) {
        unmatched.push(label.name);
      } else if (nameOf[bestI] === null) {
        nameOf[bestI] = label.name;
      } else {
        nameOf[bestI] = `${nameOf[bestI]}・${label.name}`;
        log(`ラベル衝突: ${nameOf[bestI]} (面 ${activePolygons[bestI].id})`);
      }
    }
    for (const [name, at] of Object.entries(opts.manualLabelAssign || MANUAL_LABEL_ASSIGN)) {
      const idx = unmatched.indexOf(name);
      if (idx === -1) continue;
      const pi = smallestContaining(at[0], at[1]);
      if (pi === -1) {
        log(`手動割当先 (${at[0]},${at[1]}) を含む面がありません (${name})`);
        continue;
      }
      unmatched.splice(idx, 1);
      nameOf[pi] = nameOf[pi] === null ? name : `${nameOf[pi]}・${name}`;
      log(`手動割当: ${name} → 面 ${activePolygons[pi].id}`);
    }
    if (unmatched.length > 0) {
      log(`どの面にも入らなかった国名ラベル (${unmatched.length}): ${unmatched.join(", ")}`);
    }
    log(`名前が付いた面: ${nameOf.filter((n) => n !== null).length}`);

    // 4b. 地域名レイヤーからの仮名格上げ
    //   - 面内に 州名/地名 ラベル → その名前 (島の固有名・州名・地方名)
    //   - 海域の面内に 海名 ラベル → その海名
    //   - 上記なしの島 → 最寄りの海名から「〇〇海の島NN」
    const stateLabels = extractLabels(svgDoc, "_州名");
    const placeLabels = extractLabels(svgDoc, "_地名");
    const mountainLabels = extractLabels(svgDoc, "_山岳名");
    const collectInto = (labels) => {
      const buckets = new Array(activePolygons.length).fill(null).map(() => []);
      for (const label of labels) {
        const [cx, cy] = labelCenter(label);
        const i = smallestContaining(cx, cy);
        if (i !== -1 && !buckets[i].includes(label.name)) buckets[i].push(label.name);
      }
      return buckets;
    };
    const regionNamesOf = collectInto([...stateLabels, ...placeLabels]);
    const stateNamesOf = collectInto(stateLabels);
    const placeNamesOf = collectInto(placeLabels);
    const seaLabels = extractLabels(svgDoc, "_海名").map((l) => {
      const [cx, cy] = labelCenter(l);
      return { name: l.name, cx, cy };
    });
    const seaNamesOf = new Array(activePolygons.length).fill(null).map(() => []);
    for (const sea of seaLabels) {
      const i = smallestContaining(sea.cx, sea.cy);
      if (i !== -1 && !seaNamesOf[i].includes(sea.name)) seaNamesOf[i].push(sea.name);
    }
    const centroids = activePolygons.map((p) => {
      let sx = 0;
      let sy = 0;
      for (const i of p.ring) {
        sx += points[i].x;
        sy += points[i].y;
      }
      return { x: sx / p.ring.length, y: sy / p.ring.length };
    });
    const nearestOf = (i, labels) => {
      const c = centroids[i];
      let best = null;
      let bestD = Infinity;
      for (const l of labels) {
        const d = Math.hypot(l.cx - c.x, l.cy - c.y);
        if (d < bestD) {
          bestD = d;
          best = l.name;
        }
      }
      return best;
    };
    const nearestSeaName = (i) => nearestOf(i, seaLabels);
    // 陸域の位置手掛かり: 州名 + 地名 (山岳名も地形の目印として加える)
    const landmarkLabels = [...stateLabels, ...placeLabels, ...mountainLabels].map((l) => {
      const [cx, cy] = labelCenter(l);
      return { name: l.name, cx, cy };
    });

    // 6a. 国レベルのグループ化 (階層化)
    const groupOfFace = new Array(activePolygons.length).fill(-1);
    const groups = [];
    const adoptedFaces = new Set(); // 養子 (隙間吸収で追加領有された面)
    if (hierarchical) {
      const cIn = input.countryPolygons;
      const cPts = cIn.points;
      const cTypeOf = makeTypeOf(cIn.types, cIn.polygons);
      const cPolys = [];
      for (let i = 0; i < cIn.polygons.length; i += 1) {
        const p = cIn.polygons[i];
        let r = p.points.slice();
        if (r.length >= 2 && r[0] === r[r.length - 1]) r = r.slice(0, -1);
        if (r.length < 3) continue;
        const bb = bboxOfRing(r, cPts);
        cPolys.push({ id: p.id, type: cTypeOf(i), ring: r, absArea: Math.abs(p.area), ...bb });
      }
      const smallestCountryAt = (x, y) => {
        let best = -1;
        for (let i = 0; i < cPolys.length; i += 1) {
          const c = cPolys[i];
          if (x < c.minX || x > c.maxX || y < c.minY || y > c.maxY) continue;
          if (!pointInRing(x, y, c.ring, cPts)) continue;
          if (best === -1 || c.absArea < cPolys[best].absArea) best = i;
        }
        return best;
      };
      // 州レベルの面がどの国に属するかは、**間引き前**の形で決める。
      //
      // どの国の子かは「その面がもともとどこに描かれていたか」という同一性の話で、
      // 間引きの粗さで変わってはいけない。ところが内部点を間引き後の形から取ると、
      // 境界が動いたぶん点が隣の国へはみ出し、所属が変わってしまう。実際、目標 16500 で
      // ラヴィエンド・ウィル同君王国の子が 5 面から 2 面に減り、外れた 2 面が名無しの
      // 「集約04」という別の国に化けていた (フル版では正しく 5 面)。要件で一国と決まって
      // いる箇所なので、階層が崩れるのは通せない。
      //
      // そこで内部点は normalized (間引き前) の形から取り、その対応付けをそのまま
      // 間引き後の面へ持ち越す。normalized と activePolygons は同じ面が同じ並びで
      // 1 対 1 に対応する (applyCleanup は面を増やしも減らしもしない) ので添字は通じる。
      const preParentOf = computeParents(normalized, points);
      const preBboxes = normalized.map((p) => bboxOfRing(p.ring, points));
      const preCentroids = normalized.map((p) => {
        let sx = 0;
        let sy = 0;
        for (const i of p.ring) {
          sx += points[i].x;
          sy += points[i].y;
        }
        return { x: sx / p.ring.length, y: sy / p.ring.length };
      });
      // 面の内部点。重心が外 (凹形) や自面の穴 (内包面) に落ちる場合はグリッド走査。
      const interiorOf = (i) => {
        const holes = [];
        for (let j = 0; j < normalized.length; j += 1) if (preParentOf[j] === i) holes.push(j);
        const ok = (x, y) =>
          pointInRing(x, y, normalized[i].ring, points) &&
          !holes.some((j) => pointInRing(x, y, normalized[j].ring, points));
        const c = preCentroids[i];
        if (ok(c.x, c.y)) return [c.x, c.y];
        const bb = preBboxes[i];
        for (let gy = 1; gy < 32; gy += 1) {
          for (let gx = 1; gx < 32; gx += 1) {
            const x = bb.minX + ((bb.maxX - bb.minX) * gx) / 32;
            const y = bb.minY + ((bb.maxY - bb.minY) * gy) / 32;
            if (ok(x, y)) return [x, y];
          }
        }
        return [c.x, c.y];
      };
      const byCountry = new Map();
      for (let i = 0; i < activePolygons.length; i += 1) {
        const [x, y] = interiorOf(i);
        const ci = smallestCountryAt(x, y);
        if (ci === -1) continue;
        if (!byCountry.has(ci)) byCountry.set(ci, []);
        byCountry.get(ci).push(i);
      }
      for (const u of opts.unionGroups || UNION_GROUPS) {
        const faces = [];
        let type = null;
        for (const [mx, my] of u.members) {
          const ci = smallestCountryAt(mx, my);
          if (ci === -1 || !byCountry.has(ci)) {
            log(`同君連合 ${u.name}: 構成国が見つかりません (${mx},${my})`);
            continue;
          }
          faces.push(...byCountry.get(ci));
          type = type ?? cPolys[ci].type;
          byCountry.delete(ci);
        }
        if (faces.length >= 2) {
          groups.push({ name: u.name, lift: false, childKindLabel: u.childKindLabel, type: type ?? "不明", faces });
        } else {
          log(`同君連合 ${u.name}: 構成面 ${faces.length} — 集約化を見送り`);
        }
      }
      for (const [ci, faces] of byCountry) {
        if (faces.length > 1) {
          groups.push({ name: null, lift: true, childKindLabel: "州", type: cPolys[ci].type, faces });
        }
      }
      groups.forEach((g, gi) => {
        for (const f of g.faces) groupOfFace[f] = gi;
      });
      // 隙間の吸収 (養子縁組): 子の集まりに囲まれるのに、どの子の領域にも属さない面
      // (子同士の境界に挟まれた島・海の断片) を集約の子として追加領有する。
      // visualizer 仕様の「共有海の島々は海領域に領有される」の階層版。
      const adjacency = buildFaceAdjacency(activePolygons, points.length);
      groups.forEach((g, gi) => {
        const blocked = new Set(g.faces);
        for (const f of findEnclosedFaces(adjacency, blocked, activePolygons.length)) {
          // 「他の面に内包されているか」も所属の話なので間引き前の関係で見る
          // (間引き後の関係で見ると、目標 5000 のとき 2 面が養子になり損ねた)。
          if (groupOfFace[f] !== -1 || preParentOf[f] !== -1) continue;
          g.faces.push(f);
          groupOfFace[f] = gi;
          adoptedFaces.add(f);
          log(
            `養子: 面 id=${activePolygons[f].id} (${activePolygons[f].type}) → ${g.name ?? "グループ" + (gi + 1)}`,
          );
        }
      });
      for (const g of groups) {
        log(`階層化: ${g.name ?? "(国名ラベルから命名)"} — 子 ${g.faces.length} 面 [${g.type}]`);
      }
    }

    // 6b. 地物と履歴の構築
    const histories = [];
    const typeCounters = {};
    const groupBuild = groups.map(() => ({ features: [], childIds: [], used: new Set(), childNo: 0 }));
    for (let i = 0; i < activePolygons.length; i += 1) {
      const p = activePolygons[i];
      const fid = `area-${pad(p.id, 3)}`;
      const territoryRingId = `${fid}-r0`;
      const rings = [{ kind: "territory", id: territoryRingId, vertexIds: ringVertexIds(p.ring) }];
      // 正規化で分解したポケット (橋渡し経由で繋がった島などの除外領域) を穴に
      let holeNo = 0;
      for (const loop of p.pockets) {
        holeNo += 1;
        rings.push({
          kind: "hole",
          id: `${fid}-h${holeNo}`,
          vertexIds: ringVertexIds(loop),
          parentRingId: territoryRingId,
        });
      }
      // この面を直接内包親とする面 = 穴。ただし既にポケット穴で覆われている領域
      // (島全体のポケット vs 島を分割した面など) は二重に開けない。
      for (let j = 0; j < activePolygons.length; j += 1) {
        if (parentOf[j] !== i) continue;
        const jRing = activePolygons[j].ring;
        let covered = false;
        for (const loop of p.pockets) {
          const loopSet = new Set(loop);
          const samples = jRing.filter((v) => !loopSet.has(v)).slice(0, 9);
          if (samples.length === 0) {
            covered = true;
            break;
          }
          let inCount = 0;
          for (const v of samples) {
            if (pointInRing(points[v].x, points[v].y, loop, points)) inCount += 1;
          }
          if (inCount * 2 > samples.length) {
            covered = true;
            break;
          }
        }
        if (covered) continue;
        holeNo += 1;
        rings.push({
          kind: "hole",
          id: `${fid}-h${holeNo}`,
          vertexIds: ringVertexIds(jRing),
          parentRingId: territoryRingId,
        });
      }

      // この地物が参照する頂点 (重複除去、座標変換)
      const usedIdx = new Set();
      for (const r of rings) for (const v of r.vertexIds) usedIdx.add(v);
      const vertices = sortIds(usedIdx).map((v) => ({ id: v, coordinate: coordOfId(v) }));

      const gi = groupOfFace[i];
      let kindLabel = p.type;
      let fillColor = TYPE_STYLE[p.type] ?? TYPE_STYLE["不明"];
      let name;
      if (gi !== -1 && adoptedFaces.has(i)) {
        // 養子面: 集約の子だが、命名は単独地物と同じ規則 (島名・海域名など) を使う
        name = nameOf[i];
      } else if (gi !== -1) {
        // 集約の子: lift グループでは国名ラベルを親に使うため、子は州名などから命名。
        const g = groups[gi];
        const gb = groupBuild[gi];
        gb.childNo += 1;
        kindLabel = g.childKindLabel;
        fillColor = shadeColor(TYPE_STYLE[g.type] ?? TYPE_STYLE["不明"], (gb.childNo % 5) - 2);
        name = g.lift ? null : nameOf[i];
        if (name === null && stateNamesOf[i].length > 0) {
          name = stateNamesOf[i].slice(0, 2).join("・");
          kindLabel = "州";
        }
        if (name === null && placeNamesOf[i].length > 0) {
          // 州名は無いが地名が載っている子 (ノルの区分など) はその地名を仮名に使う
          name = placeNamesOf[i].slice(0, 2).join("・");
          kindLabel = "地方";
        }
        if (name === null) {
          const m = (opts.childNameAssign || CHILD_NAME_ASSIGN).find((c) =>
            pointInRing(c.at[0], c.at[1], p.ring, points),
          );
          if (m !== undefined) {
            name = m.name;
            if (m.kindLabel !== undefined) kindLabel = m.kindLabel;
          }
        }
        // 親名はグループ確定後に決まるため、フォールバックは後で置換する
        if (name === null && p.type === "共有する海") {
          name = `__CHILD_SEA__${gb.childNo}`;
          kindLabel = "領海";
        }
        if (name === null) name = `__CHILD_FALLBACK__${gb.childNo}`;
      } else {
        name = nameOf[i];
      }
      if (name === null && (p.type === "島国" || p.type === "完全島国") && regionNamesOf[i].length > 0) {
        // 島は内側の地名 (島の固有名) をそのまま仮名に使う (最大 2 件連結)。島以外の
        // 陸域は直接命名しない (隣の沿岸国に載った島名での誤命名を防ぐ) — 後段の
        // 「〇〇付近の領域NN」形式に落とす。
        name = regionNamesOf[i].slice(0, 2).join("・");
      }
      if (name === null && p.type === "共有する海" && seaNamesOf[i].length > 0) {
        name = seaNamesOf[i].slice(0, 2).join("・");
      }
      if (name === null && (p.type === "島国" || p.type === "完全島国")) {
        const sea = nearestSeaName(i);
        if (sea !== null) {
          const key = `${sea}の島`;
          typeCounters[key] = (typeCounters[key] || 0) + 1;
          name = `${key}${pad(typeCounters[key], 2)}`;
        }
      }
      if (name === null && p.type === "共有する海") {
        const sea = nearestSeaName(i);
        if (sea !== null) {
          const key = `${sea}の海域`;
          typeCounters[key] = (typeCounters[key] || 0) + 1;
          name = `${key}${pad(typeCounters[key], 2)}`;
        }
      }
      if (name === null) {
        // 無名陸域は内包する地名/州名を優先し、無ければ最寄りの地名/州名/山岳名で
        const mark = regionNamesOf[i][0] ?? nearestOf(i, landmarkLabels);
        if (mark !== null) {
          const key = `${mark}付近の領域`;
          typeCounters[key] = (typeCounters[key] || 0) + 1;
          name = `${key}${pad(typeCounters[key], 2)}`;
        }
      }
      if (name === null) {
        const prefix = NAME_PREFIX[p.type] ?? "領域";
        typeCounters[prefix] = (typeCounters[prefix] || 0) + 1;
        name = `${prefix}${pad(typeCounters[prefix], 3)}`;
      }

      const feature = {
        kind: "leaf",
        id: fid,
        attributes: { name, startYear: 0, endYear: 9999 },
        color: {
          autoColorEnabled: false,
          fillColor,
          fillOpacity: 0.45,
          selectedColor: "#00e5ff",
          selectedOpacity: 0.7,
          borderColor: "#4a4a4a",
        },
        kindLabel,
        shape: { rings },
      };

      if (gi === -1) {
        // 単独地物 = 末端地物 1 つ = 単一錨 (year 0 〜 開いた終端) の履歴
        histories.push({
          featureId: fid,
          entries: [
            {
              anchor: { id: `${fid}-anc0`, snapshot: { features: [feature], vertices } },
              start: { granularity: "year", year: 0 },
              end: null,
            },
          ],
        });
      } else {
        const gb = groupBuild[gi];
        feature.parentId = `agg-${pad(gi + 1, 2)}`;
        gb.features.push(feature);
        gb.childIds.push(fid);
        for (const v of usedIdx) gb.used.add(v);
      }
    }

    // 6c. 集約 (親) の履歴
    groups.forEach((g, gi) => {
      const gb = groupBuild[gi];
      const aggId = `agg-${pad(gi + 1, 2)}`;
      let parentName = g.name;
      if (g.lift) {
        const names = [...new Set(g.faces.map((f) => nameOf[f]).filter((n) => n !== null))];
        parentName = names.length > 0 ? names.join("・") : `集約${pad(gi + 1, 2)}`;
      }
      // 子のフォールバック名を親名で確定 (領海/領域を別々に連番)
      let seaNo = 0;
      let restNo = 0;
      for (const f of gb.features) {
        if (/^__CHILD_SEA__\d+$/.test(f.attributes.name)) {
          seaNo += 1;
          f.attributes.name = `${parentName}領海${pad(seaNo, 2)}`;
        } else if (/^__CHILD_FALLBACK__\d+$/.test(f.attributes.name)) {
          restNo += 1;
          f.attributes.name = `${parentName}の領域${pad(restNo, 2)}`;
        }
      }
      const aggregate = {
        kind: "aggregate",
        id: aggId,
        attributes: { name: parentName, startYear: 0, endYear: 9999 },
        color: {
          autoColorEnabled: false,
          fillColor: TYPE_STYLE[g.type] ?? TYPE_STYLE["不明"],
          fillOpacity: 0.45,
          selectedColor: "#00e5ff",
          selectedOpacity: 0.7,
          borderColor: "#4a4a4a",
        },
        kindLabel: g.type,
        childIds: gb.childIds,
      };
      histories.push({
        featureId: aggId,
        entries: [
          {
            anchor: {
              id: `${aggId}-anc0`,
              snapshot: {
                features: [...gb.features, aggregate],
                vertices: sortIds(gb.used).map((v) => ({ id: v, coordinate: coordOfId(v) })),
              },
            },
            start: { granularity: "year", year: 0 },
            end: null,
          },
        ],
      });
      log(`集約: ${parentName} (${aggId}) — 子 ${gb.childIds.length} 地物`);
    });

    const projectName = opts.projectName ?? "デーレ界";
    return {
      version: "tmap/1",
      histories,
      settings: {
        project: {
          calendarId: "gregorian",
          geo: { equatorKm, flattening: 0 },
          gridSpacingDeg: 10,
          snapDistancePx: 50,
        },
        app: { autoColorDefault: true },
      },
      baseMap: {
        source: "assets/basemap.img",
        meta: { lonWest: -180, lonEast: 180, latNorth: 90, latSouth: -90 },
      },
      ui: {
        project: {
          projectName,
          projectDescription:
            opts.projectDescription ??
            `イラレ出界図.svg から derle-map-visualizer のポリゴン抽出でインポート。赤道長 ${equatorKm} km (惑星半径 1.0489 地球半径)。`,
          equatorLengthKm: equatorKm,
          flattening: 0,
          timelineMinYear: 0,
          timelineMaxYear: 10000,
          gridVisible: true,
          gridSpacingDeg: 10,
          gridColor: "#888888",
          gridOpacity: 0.5,
          baseMapId: "assets/basemap.img",
          areaLabelOccupancyThreshold: 0.02,
          newAreaAutoColorDefault: true,
          defaultPaletteTheme: "default",
          customPalette: [],
          zoomMin: 1,
          zoomMax: 50,
        },
        app: {
          snapDistancePx: 50,
          renderFps: 60,
          vertexMarkerCap: 2000,
          showConflictGuide: true,
          undoStackSize: 100,
        },
      },
    };
  }

  // --- 7. 組み上がった manifest の自前整合性チェック ---------------------------
  //
  // なぜ要るか: 間引き (options.targetVertices) は各リングを独立に動かすので、まれに
  //   (a) 穴が親のリングからはみ出す / 親のリングを横切る
  //   (b) 隣り合う面どうしが重なる (境界が相手の内側へ食い込む)
  // が起きる。どの面がそうなるかは幾何次第で事前に読めないため、書き出しの直前に
  // 自分で検査して、落ちた面だけ間引きを解除して組み直す (下の exportGimoza の収束
  // ループ)。以前は Node 側のドライバがアプリの checkConsistency を呼んでいたが、
  // visualizer はブラウザだけで完結する必要があるのでここに同等の判定を持つ。
  //
  // アプリ (time-map-app-4) のどの判定に対応するか (実データで実際に出た 3 種):
  //   ring-hierarchy-invalid / hole-crosses-parent   → 穴リングが親リングと真に交差
  //   ring-hierarchy-invalid / hole-not-contained    → 穴リングの頂点が親リングの外
  //   terminal-exclusion (= anchor-spatial-overlap)  → 末端地物どうしが正の面積で重なる
  // アプリ側はこの 3 種を validateRingHierarchy と polygon-clipping の交差面積で判定
  // するが、ここでは同じ「意味」を素の幾何述語だけで書く (visualizer は依存なし)。
  //
  // 判定を合わせるうえで大事な約束 (アプリ側と同じ):
  //   - 境界の接触 (辺や頂点の共有) は重なりに数えない。だから交差は「真の交差」
  //     (双方の線分が相手の直線を跨ぐ) だけを見る。
  //   - 穴は「親の外に出ていないこと」で判定する。境界上は許す。
  //   - 経度は連続量として扱う (±180 度をまたぐ列を復元してから比較する)。

  /**
   * 経度の列を連続化する。隣どうしの差を「近い方の回り」(±180 度以内) に直して積む。
   * アプリ側 (src/domain/geo/longitude.ts の unwrapLonSequence) と同じ規則。
   */
  function unwrapLons(lons) {
    const out = [];
    let prevRaw = null;
    let prev = 0;
    for (const lon of lons) {
      if (prevRaw === null) {
        out.push(lon);
        prevRaw = lon;
        prev = lon;
        continue;
      }
      let d = (((lon - prevRaw + 180) % 360) + 360) % 360 - 180;
      if (d === -180) d = 180;
      const v = prev + d;
      out.push(v);
      prevRaw = lon;
      prev = v;
    }
    return out;
  }

  /**
   * 頂点 id の並び → 比較用の経緯度列。連続化したうえで、平均経度が [-180,180) に
   * 収まるよう全体を 360 度単位でずらす。こうしておくと、別々に連続化したリングどうしも
   * 同じ「周」に載るので、そのまま比べられる (アプリ側の alignRingTo に相当)。
   *
   * 最後に 1e-7 度の格子へ丸めるのが要 (アプリ側の quantizeCoordinate と同じ刻み)。
   * 連続化は「前の値 + 差分」を積んでいくので、同じ頂点でもリングごとに 13 桁目あたりが
   * ずれる。境界を共有する 2 面はまったく同じ頂点を使うはずなのに、このずれのせいで
   * 「わずかに交差している」と誤判定してしまう (実測: フル版で 913 組の偽の重なり)。
   * 格子へ丸めれば共有頂点はきっちり同じ値に戻り、共有境界は交差ではなく接触になる。
   */
  const QUANTIZE = 1e7;
  const quantize = (v) => Math.round(v * QUANTIZE) / QUANTIZE;
  function ringGeo(vertexIds, coordOf) {
    const raw = vertexIds.map(coordOf);
    const lons = unwrapLons(raw.map((c) => c.lon));
    let sum = 0;
    for (const l of lons) sum += l;
    const mean = lons.length > 0 ? sum / lons.length : 0;
    const shift = -Math.round(mean / 360) * 360;
    return raw.map((c, i) => ({ lon: quantize(lons[i] + shift), lat: quantize(c.lat) }));
  }

  /** 3 点の外積 (符号が向きを表す)。0 なら一直線上。 */
  const crossSign = (a, b, c) => (b.lon - a.lon) * (c.lat - a.lat) - (b.lat - a.lat) * (c.lon - a.lon);

  /**
   * 2 本の線分が「真に交差」するか。端点だけの接触や、重なって並んでいるだけの状態は
   * 交差に数えない (境界の共有は重なりではない、というアプリ側の約束と同じ)。
   *
   * 面はすべて同じ平面分割から取れており、境界を共有する 2 面はまったく同じ頂点列を
   * 使う (間引きの残す/残さないは頂点ごとに決まるので、共有チェーンは両側で同じに
   * なる)。つまり共有辺の外積はきっちり 0 になるので、浮動小数のままで判定できる。
   */
  function properlyCross(p1, p2, p3, p4) {
    const d1 = crossSign(p3, p4, p1);
    const d2 = crossSign(p3, p4, p2);
    const d3 = crossSign(p1, p2, p3);
    const d4 = crossSign(p1, p2, p4);
    return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
  }

  /** 点が線分上に載っているか。 */
  function onSegment(p, a, b) {
    if (crossSign(a, b, p) !== 0) return false;
    return (
      p.lon >= Math.min(a.lon, b.lon) &&
      p.lon <= Math.max(a.lon, b.lon) &&
      p.lat >= Math.min(a.lat, b.lat) &&
      p.lat <= Math.max(a.lat, b.lat)
    );
  }

  /**
   * リングの検索用インデックス。緯度を帯 (バンド) に切り、各辺をその辺が跨ぐ帯すべてに
   * 登録する。内外判定は水平のレイを飛ばすので「その緯度の帯にある辺」だけ見ればよく、
   * 交差判定も相手の辺の緯度範囲の帯だけを見ればよい。全部を総当たりすると
   * 数千頂点の海岸線 × 数百面で現実的な時間に収まらないため。
   */
  function makeRingIndex(geo) {
    const n = geo.length;
    let minLon = Infinity;
    let maxLon = -Infinity;
    let minLat = Infinity;
    let maxLat = -Infinity;
    for (const c of geo) {
      if (c.lon < minLon) minLon = c.lon;
      if (c.lon > maxLon) maxLon = c.lon;
      if (c.lat < minLat) minLat = c.lat;
      if (c.lat > maxLat) maxLat = c.lat;
    }
    const bandCount = Math.max(1, Math.min(512, Math.ceil(Math.sqrt(n))));
    const height = maxLat - minLat;
    const span = height > 0 ? height : 1e-12;
    const bandOf = (lat) =>
      Math.max(0, Math.min(bandCount - 1, Math.floor(((lat - minLat) / span) * bandCount)));
    const bands = Array.from({ length: bandCount }, () => []);
    for (let i = 0; i < n; i += 1) {
      const a = geo[i];
      const b = geo[(i + 1) % n];
      const lo = bandOf(Math.min(a.lat, b.lat));
      const hi = bandOf(Math.max(a.lat, b.lat));
      for (let k = lo; k <= hi; k += 1) bands[k].push(i);
    }
    return { geo, n, minLon, maxLon, minLat, maxLat, bandOf, bands, stamp: new Int32Array(n).fill(-1) };
  }

  const bboxesOverlap = (a, b) =>
    a.minLon <= b.maxLon && b.minLon <= a.maxLon && a.minLat <= b.maxLat && b.minLat <= a.maxLat;

  /** 点とリングの位置関係: "in" (内) / "out" (外) / "on" (境界上)。 */
  function locateInRing(idx, p) {
    if (p.lon < idx.minLon || p.lon > idx.maxLon || p.lat < idx.minLat || p.lat > idx.maxLat) {
      return "out";
    }
    let inside = false;
    for (const i of idx.bands[idx.bandOf(p.lat)]) {
      const a = idx.geo[i];
      const b = idx.geo[(i + 1) % idx.n];
      if (onSegment(p, a, b)) return "on";
      if (
        a.lat > p.lat !== b.lat > p.lat &&
        p.lon < ((b.lon - a.lon) * (p.lat - a.lat)) / (b.lat - a.lat) + a.lon
      ) {
        inside = !inside;
      }
    }
    return inside ? "in" : "out";
  }

  /** リング a の辺が、インデックス b のリングの辺と真に交差するか。 */
  let crossStampSeq = 0;
  function ringsProperlyCross(a, bIdx) {
    if (
      a.minLon > bIdx.maxLon ||
      bIdx.minLon > a.maxLon ||
      a.minLat > bIdx.maxLat ||
      bIdx.minLat > a.maxLat
    ) {
      return false;
    }
    const ga = a.geo;
    for (let i = 0; i < ga.length; i += 1) {
      const p1 = ga[i];
      const p2 = ga[(i + 1) % ga.length];
      const loLat = Math.min(p1.lat, p2.lat);
      const hiLat = Math.max(p1.lat, p2.lat);
      if (hiLat < bIdx.minLat || loLat > bIdx.maxLat) continue;
      const loLon = Math.min(p1.lon, p2.lon);
      const hiLon = Math.max(p1.lon, p2.lon);
      if (hiLon < bIdx.minLon || loLon > bIdx.maxLon) continue;
      // 同じ辺を複数の帯から二重に見ないよう、走査ごとの通し番号で印を付ける
      crossStampSeq += 1;
      const lo = bIdx.bandOf(loLat);
      const hi = bIdx.bandOf(hiLat);
      for (let k = lo; k <= hi; k += 1) {
        for (const j of bIdx.bands[k]) {
          if (bIdx.stamp[j] === crossStampSeq) continue;
          bIdx.stamp[j] = crossStampSeq;
          const q1 = bIdx.geo[j];
          const q2 = bIdx.geo[(j + 1) % bIdx.n];
          if (Math.min(q1.lon, q2.lon) > hiLon || Math.max(q1.lon, q2.lon) < loLon) continue;
          if (Math.min(q1.lat, q2.lat) > hiLat || Math.max(q1.lat, q2.lat) < loLat) continue;
          if (properlyCross(p1, p2, q1, q2)) return true;
        }
      }
    }
    return false;
  }

  /**
   * 点が地物の「実面積」の内側か。外周の内側で、かつどの穴の内側でもないとき "in"。
   * 境界 (外周でも穴でも) に載っていれば "on"。アプリ側の pointInPolygon と同じ扱い。
   */
  function locateInFace(face, p) {
    const outer = locateInRing(face.rings[0], p);
    if (outer !== "in") return outer;
    for (let k = 1; k < face.rings.length; k += 1) {
      const h = locateInRing(face.rings[k], p);
      if (h === "on") return "on";
      if (h === "in") return "out";
    }
    return "in";
  }

  /**
   * リング ring の頂点のうち 1 つでも地物 face の実面積の内側 (境界上は含まない) に
   * 入っているか。
   *
   * 「最初にはっきり決まった 1 点」で打ち切ってはいけない。間引きで、同じ 2 頂点の
   * あいだを面 A は直線 1 本で、面 B は途中の頂点を残した折れ線で結ぶことがある。
   * このとき 2 面はレンズ状に重なるのに、境界は両端の共有頂点で触れているだけなので
   * 「真の交差」にはならない。重なりの証拠は B の途中の頂点が A の内側にあることだけで、
   * それはリングの並びの途中に現れる (実データ: 面 112 と 113 がこの形)。
   */
  function anyVertexStrictlyInside(face, ring) {
    for (const p of ring.geo) {
      if (locateInFace(face, p) === "in") return true;
    }
    return false;
  }

  /**
   * manifest を読んで、アプリの整合性チェックが落とすであろう面を洗い出す。
   * 返り値は `{ kind, faceIds }` の配列。`faceIds` は面 id (area-NNN の NNN) で、
   * そのまま keepFullIds に入れれば「その面の間引きだけ解除」できる。
   */
  function checkManifestConsistency(manifest) {
    // 1. 末端地物 (leaf) を集め、リングを比較用の幾何に直す
    const faces = [];
    const seen = new Set();
    for (const h of manifest.histories) {
      for (const e of h.entries) {
        const snap = e.anchor.snapshot;
        const coordById = new Map(snap.vertices.map((v) => [v.id, v.coordinate]));
        const coordOf = (id) => coordById.get(id);
        for (const f of snap.features) {
          if (f.kind !== "leaf" || seen.has(f.id)) continue;
          seen.add(f.id);
          const m = /^area-(\d+)$/.exec(f.id);
          const territory = f.shape.rings.find((r) => r.kind === "territory");
          if (territory === undefined) continue;
          const holes = f.shape.rings.filter((r) => r.kind === "hole");
          faces.push({
            featureId: f.id,
            faceId: m === null ? null : Number(m[1]),
            // rings[0] が外周、以降が穴 (locateInFace がこの並びを前提にする)
            rings: [territory, ...holes].map((r) => makeRingIndex(ringGeo(r.vertexIds, coordOf))),
          });
        }
      }
    }

    const out = [];
    const add = (kind, ...faceObjs) => {
      const ids = faceObjs.map((f) => f.faceId).filter((v) => v !== null);
      if (ids.length > 0) out.push({ kind, faceIds: ids });
    };

    // 2. 地物ごとのリング階層 (穴が親の中にきちんと入っているか / 穴どうしが重ならないか)
    for (const face of faces) {
      const outer = face.rings[0];
      let bad = false;
      for (let k = 1; k < face.rings.length && !bad; k += 1) {
        const hole = face.rings[k];
        if (ringsProperlyCross(hole, outer)) {
          add("hole-crosses-parent", face);
          bad = true;
          break;
        }
        // 穴の頂点が 1 つでも外周の外に出ていれば違反 (境界に載っているだけの頂点は
        // 接触なので数えない)。アプリ側も交差していない場合はこの見方をしている。
        for (const p of hole.geo) {
          if (locateInRing(outer, p) === "out") {
            add("hole-not-contained", face);
            bad = true;
            break;
          }
        }
      }
      if (bad) continue;
      // 兄弟の穴どうしの重なり (実データでは出ていないが、アプリ側が見ているので揃える)
      for (let a = 1; a < face.rings.length && !bad; a += 1) {
        for (let b = a + 1; b < face.rings.length; b += 1) {
          const ra = face.rings[a];
          const rb = face.rings[b];
          if (!bboxesOverlap(ra, rb)) continue;
          if (
            ringsProperlyCross(ra, rb) ||
            anyVertexStrictlyInside({ rings: [rb] }, ra) ||
            anyVertexStrictlyInside({ rings: [ra] }, rb)
          ) {
            add("siblings-overlap", face);
            bad = true;
            break;
          }
        }
      }
    }

    // 3. 末端地物どうしの重なり (アプリの 末端排他 / 同時刻の空間競合)
    //    重なりの正体は 2 つだけ: 境界が真に交差するか、片方が丸ごと相手の実面積の
    //    内側に入るか。穴として開けてある面は「穴の内側」なので重なりにならない。
    for (let i = 0; i < faces.length; i += 1) {
      const fi = faces[i];
      for (let j = i + 1; j < faces.length; j += 1) {
        const fj = faces[j];
        if (!bboxesOverlap(fi.rings[0], fj.rings[0])) continue;
        let hit = false;
        for (const ra of fi.rings) {
          for (const rb of fj.rings) {
            if (bboxesOverlap(ra, rb) && ringsProperlyCross(ra, rb)) {
              hit = true;
              break;
            }
          }
          if (hit) break;
        }
        if (!hit && (anyVertexStrictlyInside(fi, fj.rings[0]) || anyVertexStrictlyInside(fj, fi.rings[0]))) {
          hit = true;
        }
        if (hit) add("terminal-overlap", fi, fj);
      }
    }
    return out;
  }

  // --- 8. ベクター背景地図 ----------------------------------------------------
  /** rgb(a,b,c) → #rrggbb (JSON を小さく・描画側の扱いを単純に)。 */
  function toHexColor(css) {
    const m = /^rgba?\(\s*([\d.]+)[,\s]+([\d.]+)[,\s]+([\d.]+)/.exec(css);
    if (m === null) return css;
    const h = (v) => Math.round(Number(v)).toString(16).padStart(2, "0");
    return `#${h(m[1])}${h(m[2])}${h(m[3])}`;
  }

  /**
   * 要素の塗り色を決める。style 属性 → fill 属性 → class の CSS → 親から継承 の順。
   *
   * Node 版は本物のブラウザに SVG を載せて getComputedStyle で解決していたが、
   * DOMParser で読んだ文書は生きた DOM に繋がっていないので getComputedStyle が効かない。
   * この SVG の地形レイヤーは全て `class="cls-N"` だけで塗り色が決まっているため、
   * `<style>` を読んで同じ結果になる。どこにも塗りの指定が無い要素は null を返して
   * 描画対象から外す (Node 版では既定の黒として拾われていたが、対象レイヤーには
   * 該当する要素が無いので出力は変わらない)。
   */
  function resolveFill(el, fillMap) {
    let cur = el;
    while (cur && cur.nodeType === 1) {
      const style = cur.getAttribute("style");
      if (style) {
        const m = /(?:^|[;\s])fill:\s*([^;]+)/.exec(style);
        if (m) return m[1].trim();
      }
      const attr = cur.getAttribute("fill");
      if (attr) return attr.trim();
      const cls = cur.getAttribute("class");
      if (cls) {
        for (const c of cls.split(/\s+/)) {
          if (fillMap.has(c)) return fillMap.get(c);
        }
      }
      cur = cur.parentNode;
    }
    return null;
  }

  function isHidden(el) {
    let cur = el;
    while (cur && cur.nodeType === 1) {
      const style = cur.getAttribute("style");
      if (style && /(?:^|[;\s])display:\s*none/.test(style)) return true;
      if (cur.getAttribute("display") === "none") return true;
      cur = cur.parentNode;
    }
    return false;
  }

  /**
   * 出界図 SVG から「ベクター背景地図」を取り出す (assets/basemap.img の中身)。
   *
   * ラスター (PNG) の代わりに .gimoza へ入れるデータで、ズームしてもぼやけない。
   * 形式は SVG のパス定義 (d 属性) をそのまま保持する:
   *   - 曲線がそのまま残るので、点列にサンプリングするより小さく・かつ正確
   *   - 描画側は Path2D でブラウザにパースさせて塗るだけ (実測 世界全体 19.8ms /
   *     10 倍ズーム 5.6ms — 8192×4096 PNG のデコード 162ms より速い)
   */
  function extractVectorBaseMap(svgDoc, worldRect) {
    const rect = worldRect || worldRectFromSvg(svgDoc);
    const fillMap = cssDeclMap(svgDoc, "fill");
    const groups = [];
    for (const id of BASEMAP_LAYER_IDS) {
      const g = svgDoc.getElementById(id) || svgDoc.querySelector(`[data-name='${id.replace(/^_/, "")}']`);
      if (g === null) continue;
      // 同じ塗り色ごとにパスをまとめる (描画時の状態変更を減らす)
      const byFill = new Map();
      const nodes = g.querySelectorAll("path, polygon, rect");
      for (let n = 0; n < nodes.length; n += 1) {
        const p = nodes[n];
        const fill = resolveFill(p, fillMap);
        if (fill === null || fill === "none" || fill === "" || isHidden(p)) continue;
        const tag = p.tagName.toLowerCase();
        let d = null;
        if (tag === "path") {
          d = p.getAttribute("d");
        } else if (tag === "rect") {
          const x = Number(p.getAttribute("x") ?? 0);
          const y = Number(p.getAttribute("y") ?? 0);
          const w = Number(p.getAttribute("width") ?? 0);
          const h = Number(p.getAttribute("height") ?? 0);
          d = `M${x},${y}H${x + w}V${y + h}H${x}Z`;
        } else if (tag === "polygon") {
          const ptsAttr = (p.getAttribute("points") ?? "").trim();
          if (ptsAttr.length > 0) d = `M${ptsAttr}Z`;
        }
        if (d === null || d.length === 0) continue;
        if (!byFill.has(fill)) byFill.set(fill, []);
        byFill.get(fill).push(d);
      }
      for (const [fill, paths] of byFill) groups.push({ fill: toHexColor(fill), paths });
    }
    // 海の色 (背景の下地) — 海色レクトの塗り
    const sea = svgDoc.getElementById("_海色") || svgDoc.querySelector("[data-name='海色']");
    const background = sea === null ? "#4a6fb5" : toHexColor(resolveFill(sea, fillMap) ?? "#4a6fb5");

    return {
      kind: "tmap/basemap-vector@1",
      // SVG ユーザー座標のどの矩形が世界矩形 (経度 -180..180 / 緯度 90..-90) に対応するか。
      // 描画側はこれだけで SVG 座標 → 画面 の変換を決められる。
      viewBox: { x: rect.x, y: rect.y, width: rect.width, height: rect.height },
      background,
      layers: groups,
    };
  }

  // --- 9. 間引き版の収束ループ ------------------------------------------------
  /**
   * 目標頂点数を指定した間引き版の manifest を、整合性が取れるまで組み直して返す。
   *
   * 間引きは各リングを独立に動かすので、まれに穴が親からはみ出したり、隣の面と
   * 重なったりする。どの面がそうなるかは幾何次第で事前に読めないため、
   *   組み立て → 自前チェック → 落ちた面を keepFullIds へ追加 → 組み立て直し
   * を繰り返す。keepFullIds に入った面はその面だけ原型解像度に戻る (局所的な
   * 間引き解除) ので、頂点数はほとんど増えないまま不整合だけが消える。
   *
   * 上限まで回っても消えない場合は、最後に作ったものをそのまま返す (書き出せない
   * よりは書き出して、進捗メッセージで残った件数を伝えるほうがよい)。
   */
  function buildSimplifiedManifest(input, opts) {
    const progress = typeof opts.onProgress === "function" ? opts.onProgress : () => {};
    const log = typeof opts.log === "function" ? opts.log : () => {};
    const keepFullIds = new Set(opts.keepFullIds ?? []);
    const MAX_ROUNDS = 6;
    let manifest = null;
    for (let round = 1; round <= MAX_ROUNDS; round += 1) {
      progress(`間引き ${round} 回目を組み立て中…`);
      manifest = buildManifest({ ...input, options: { ...opts, keepFullIds } });
      progress(`間引き ${round} 回目の整合性を確認中…`);
      const violations = checkManifestConsistency(manifest);
      if (violations.length === 0) {
        progress(
          `整合性チェック: 問題なし (${round} 回目 / 原型維持 ${keepFullIds.size} 面)`,
        );
        return manifest;
      }
      let added = 0;
      for (const v of violations) {
        for (const id of v.faceIds) {
          if (!keepFullIds.has(id)) {
            keepFullIds.add(id);
            added += 1;
          }
        }
      }
      log(`整合性チェック ${round} 回目: ${violations.length} 件 → ${added} 面を原型維持へ`);
      progress(`不整合 ${violations.length} 件 → ${added} 面の間引きを解除して組み直します`);
      if (added === 0) {
        // 新たに救済できる面が無い = これ以上ループしても同じ結果にしかならない
        progress(`整合性チェック: ${violations.length} 件が残りました (救済できる面なし)`);
        return manifest;
      }
    }
    progress(`整合性チェック: ${MAX_ROUNDS} 回で収束しませんでした`);
    return manifest;
  }

  // --- 10. .gimoza の書き出し -------------------------------------------------
  /**
   * 抽出結果から .gimoza (ZIP) のバイト列を作る。
   * `input` は buildManifest と同じ。`input.baseMapAsset` (Uint8Array) を渡すと
   * 背景資産をそれで差し替える (ラスター PNG を入れたい場合など)。
   */
  async function exportGimoza(input) {
    const opts = input.options || {};
    const worldRect = opts.worldRect || worldRectFromSvg(input.svgDoc);
    const manifest =
      opts.targetVertices == null
        ? buildManifest({ ...input, options: { ...opts, worldRect } })
        : buildSimplifiedManifest(input, { ...opts, worldRect });
    const encoder = new TextEncoder();
    const baseMapAsset =
      input.baseMapAsset ?? encoder.encode(JSON.stringify(extractVectorBaseMap(input.svgDoc, worldRect)));
    return buildZip({
      "manifest.json": encoder.encode(JSON.stringify(manifest)),
      "assets/basemap.img": baseMapAsset,
    });
  }

  global.Gimoza = {
    buildZip,
    crc32,
    extractVectorBaseMap,
    worldRectFromSvg,
    buildManifest,
    buildSimplifiedManifest,
    checkManifestConsistency,
    exportGimoza,
    // 下位の部品 (検証・調整用)
    cleanRing,
    decomposeRing,
    simplifyKeepSet,
    extractLabels,
    ringGeo,
    makeRingIndex,
    locateInRing,
    ringsProperlyCross,
  };
})(window);
