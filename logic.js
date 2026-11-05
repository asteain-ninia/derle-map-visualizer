(function (global) {
  "use strict";

  const EPS = 1e-9;
  let measureSvg = null;

  function setMeasureSvg(svgEl) {
    measureSvg = svgEl;
  }

  function nextFrame() {
    // Node（node --test）には requestAnimationFrame が無い。描画の機会を作るための待ちなので、
    // 無い環境ではマクロタスク 1 回ぶん待つだけにする（ブラウザでの挙動は変えない）。
    if (typeof requestAnimationFrame !== "function") {
      return new Promise((resolve) => setTimeout(resolve, 0));
    }
    return new Promise((resolve) => requestAnimationFrame(() => resolve()));
  }

  // --- `d` のサブパス分割 -------------------------------------------------------
  //
  // 点列化はサブパス（M/m で始まる 1 区切り）ごとに行う（AGENT.md「複数サブパスの誤結線」）。
  // 切り出したサブパスは単独の `d` として解釈されるが、SVG では **`d` の最初の `m` は
  // 絶対 `M` として解釈される**。2 つ目以降のサブパスが相対 `m` で始まっていると、
  // 本来「直前のサブパスの現在点からの相対」だった始点が原点からの絶対座標に化け、
  // エラーも警告も出ないまま形だけが移動する（旧実装はこれを数えて警告するだけだった）。
  //
  // そこで分割しながら現在点を追い、相対 `m` で始まるサブパスは始点を絶対 `M x y` に
  // 書き換えて返す。`m` の後ろに続く座標対は暗黙の相対 `l` なので、`M x y l …` として
  // 相対のまま残す。それ以外のコマンドの文字列は元のまま（絶対 `M` で始まるサブパスは
  // 1 文字も変えない）。

  const PATH_PARAM_COUNT = {
    M: 2, L: 2, T: 2, H: 1, V: 1, C: 6, S: 4, Q: 4, A: 7, Z: 0,
  };
  const PATH_NUMBER_RE = /[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?/y;

  /** 旧実装（現在点を追わない単純な分割）。パスデータが壊れているときだけ使う。 */
  function splitPathDataNaive(d) {
    return d.trim().split(/(?=[Mm])/).map((part) => part.trim()).filter(Boolean);
  }

  function isPathSeparator(code) {
    // 空白（SP, TAB, LF, FF, CR）とカンマ
    return code === 32 || code === 9 || code === 10 || code === 12 || code === 13 || code === 44;
  }

  function splitPathData(d) {
    if (!d) return [];
    const src = String(d);
    // 2 文字目以降に相対 m が無ければ、書き換えるサブパスが無いので単純な分割と同じ結果になる。
    // 原図のほとんどの `d` はこれ（全体を読むより 1 桁速い）
    if (src.indexOf("m", src.search(/\S/) + 1) === -1) return splitPathDataNaive(src);
    const n = src.length;
    let pos = 0;

    const skipSeparators = () => {
      while (pos < n && isPathSeparator(src.charCodeAt(pos))) pos += 1;
    };
    // 数値を 1 つ読む。読めなければ null（pos は動かさない）
    const readNumber = () => {
      skipSeparators();
      PATH_NUMBER_RE.lastIndex = pos;
      const m = PATH_NUMBER_RE.exec(src);
      if (!m) return null;
      pos += m[0].length;
      return Number(m[0]);
    };
    // 円弧のフラグは 1 文字の 0/1。"a5 5 0 1010 10" のように区切りなしで詰めて書ける
    // ので、数値として読むと "1010" を 1 つの数にしてしまう
    const readFlag = () => {
      skipSeparators();
      const c = src[pos];
      if (c !== "0" && c !== "1") return null;
      pos += 1;
      return c === "1" ? 1 : 0;
    };
    // 次がコマンド文字か終端か（＝このコマンドの引数が尽きたか）
    const atCommandOrEnd = () => {
      skipSeparators();
      if (pos >= n) return true;
      const c = src.charCodeAt(pos);
      // 'e'/'E' は数値の指数にしか現れないので、ここではコマンドとして扱わない
      return ((c >= 65 && c <= 90) || (c >= 97 && c <= 122)) && c !== 69 && c !== 101;
    };

    const parts = [];
    let x = 0;
    let y = 0;
    let startX = 0;
    let startY = 0;
    let partStart = -1; // 今のサブパスの先頭位置（src 内）
    let prefix = null; // 今のサブパスの先頭を書き換えるとき { text, resumeAt }
    const flush = (end) => {
      if (partStart < 0) return;
      const text = prefix
        ? prefix.text + src.slice(prefix.resumeAt, end)
        : src.slice(partStart, end);
      const trimmed = text.trim();
      if (trimmed) parts.push(trimmed);
    };

    skipSeparators();
    // M/m より前に何かある（壊れた `d`）。ブラウザは描かない。旧実装と同じ扱いにする
    if (pos < n && src[pos] !== "M" && src[pos] !== "m") return splitPathDataNaive(src);

    while (true) {
      skipSeparators();
      if (pos >= n) break;
      const cmd = src[pos];
      const upper = cmd.toUpperCase();
      const count = PATH_PARAM_COUNT[upper];
      if (count === undefined) return splitPathDataNaive(src); // 未知の文字
      const cmdPos = pos;
      pos += 1;
      const rel = cmd !== upper;

      if (upper === "M") {
        flush(cmdPos);
        partStart = cmdPos;
        prefix = null;
        const dx = readNumber();
        const dy = dx === null ? null : readNumber();
        if (dy === null) return splitPathDataNaive(src);
        // `d` の最初の m は絶対として解釈される（SVG 仕様）
        const relStart = rel && parts.length > 0;
        x = (rel ? x : 0) + dx;
        y = (rel ? y : 0) + dy;
        startX = x;
        startY = y;
        if (relStart) {
          // 続く座標対（暗黙の相対 l）は元の文字列のまま、`l` を補って残す
          prefix = { text: `M${x} ${y}`, resumeAt: pos };
          skipSeparators();
          if (!atCommandOrEnd()) prefix.text += "l";
          prefix.resumeAt = pos;
        }
        // 残りの座標対は暗黙の L / l
        while (!atCommandOrEnd()) {
          const px = readNumber();
          const py = px === null ? null : readNumber();
          if (py === null) return splitPathDataNaive(src);
          x = (rel ? x : 0) + px;
          y = (rel ? y : 0) + py;
        }
        continue;
      }

      if (upper === "Z") {
        x = startX;
        y = startY;
        continue;
      }

      // 引数の組を、次のコマンド文字が来るまで繰り返し読む（同じコマンドの暗黙の繰り返し）
      let groups = 0;
      while (!atCommandOrEnd()) {
        const v = [];
        for (let k = 0; k < count; k += 1) {
          const isFlag = upper === "A" && (k === 3 || k === 4);
          const value = isFlag ? readFlag() : readNumber();
          if (value === null) return splitPathDataNaive(src);
          v.push(value);
        }
        if (upper === "H") {
          x = (rel ? x : 0) + v[0];
        } else if (upper === "V") {
          y = (rel ? y : 0) + v[0];
        } else {
          // 終点はどのコマンドでも最後の 2 つ
          x = (rel ? x : 0) + v[count - 2];
          y = (rel ? y : 0) + v[count - 1];
        }
        groups += 1;
      }
      if (groups === 0) return splitPathDataNaive(src); // 引数の無い描画コマンド
    }
    flush(n);
    return parts;
  }

  function applyTransform(point, element) {
    if (!measureSvg) return { x: point.x, y: point.y };
    const ctm = element.getCTM();
    if (!ctm) return { x: point.x, y: point.y };
    const svgPoint = measureSvg.createSVGPoint();
    svgPoint.x = point.x;
    svgPoint.y = point.y;
    const out = svgPoint.matrixTransform(ctm);
    return { x: out.x, y: out.y };
  }

  // --- パスの点列化 -----------------------------------------------------------
  //
  // 方式が 2 つある。既定は "pathdata"（速い方）で、"length" は比較用に残してある。
  //
  //   length   : 旧方式。getTotalLength / getPointAtLength で等間隔にサンプリングする。
  //              **getPointAtLength の 1 回の呼び出しコストがパスの複雑さに比例する**ため、
  //              サンプル数もパス長に比例することと合わせて 1 本のパスで二乗的に増える。
  //              実測（デーレ界の出界図・全 832 パス）:
  //                国境界      5.7 秒 / 45,692 回 =   125 μs/回
  //                第一分高面 171.3 秒 / 77,320 回 = 2,216 μs/回
  //                最悪のパス 1 本（d 属性 64,132 文字）だけで 78.7 秒（7,121 μs/回）
  //              グラフ構築 190 秒のうち 184 秒（97%）がこの段階だった。
  //
  //   pathdata : 新方式。`d` 属性を自前で解釈して点列にする。直線コマンド（L/H/V）は
  //              元の頂点をそのまま使い、曲線（C/S/Q/T）だけ弦からのずれが許容差に入る
  //              まで再帰分割する。ブラウザの API を呼ばないので複雑さに二乗で効かない。
  //              同じ入力で 国境界 16 ミリ秒 / 第一分高面 27 ミリ秒（約 4,000 倍）。
  //              直線が元の 2 点のまま残る利点もある（等間隔サンプリングは直線も刻む）。
  //
  // 切替は setSampleMode() / getSampleMode()。両方式の突き合わせ（形状のずれ・点数・
  // 時間）は logic_test.html から実行できる。
  const SAMPLE_MODE_LENGTH = "length";
  const SAMPLE_MODE_PATHDATA = "pathdata";
  let sampleMode = SAMPLE_MODE_PATHDATA;
  /** 曲線分割の許容差（地図単位）。制御点の弦からのずれがこれ以下なら直線とみなす。 */
  let curveTolerance = 0.05;

  function setSampleMode(mode) {
    if (mode === SAMPLE_MODE_LENGTH || mode === SAMPLE_MODE_PATHDATA) {
      sampleMode = mode;
    }
    return sampleMode;
  }

  function getSampleMode() {
    return sampleMode;
  }

  function setCurveTolerance(value) {
    const v = Number(value);
    if (Number.isFinite(v) && v > 0) curveTolerance = v;
    return curveTolerance;
  }

  function getCurveTolerance() {
    return curveTolerance;
  }

  /** 旧方式: getPointAtLength による等間隔サンプリング。比較用に残す。 */
  function samplePathElementByLength(pathEl, step, pathDataOverride) {
    if (!measureSvg) {
      throw new Error("measureSvg is not set.");
    }
    const clone = pathEl.cloneNode(true);
    if (pathDataOverride) {
      clone.setAttribute("d", pathDataOverride);
    }
    measureSvg.appendChild(clone);
    const length = clone.getTotalLength();
    const count = Math.max(2, Math.ceil(length / step));
    const points = [];
    for (let i = 0; i <= count; i += 1) {
      const raw = clone.getPointAtLength((length * i) / count);
      points.push(applyTransform(raw, clone));
    }
    measureSvg.removeChild(clone);
    return points;
  }

  /**
   * 三次ベジェを折れ線へ。制御点の弦からのずれが tol 以下になるまで半分に割る
   * （適応分割）。曲率の大きいところだけ細かくなるので、等間隔サンプリングより
   * 少ない点で同じ形を表せる。depth は退化入力（弦長 0 付近）での無限再帰よけ。
   *
   * minChord は「これ以上は細かく割らない」下限。**吸着（buildGraphFromSegments）が
   * 許容値以内の点どうしを 1 点に併合する**ため、許容値より密に点を置くと、小さな島の
   * 輪郭がまるごと 1 点に潰れて面が消える。実測（下限なしの版）: 点 99,385 → 50,059、
   * 面 701 → 634（島国が 58 個消失）。旧方式は間引き間隔ちょうどの等間隔だったので
   * この現象が起きなかった。下限を間引き間隔に合わせることで、点の粗さが旧方式を
   * 下回らないこと＝吸着の挙動が変わらないことを保証する。
   */
  function flattenCubic(out, x0, y0, x1, y1, x2, y2, x3, y3, tol, depth, minChord) {
    const dx = x3 - x0;
    const dy = y3 - y0;
    const chord = Math.hypot(dx, dy);
    const dev =
      chord < EPS
        ? Math.max(Math.hypot(x1 - x0, y1 - y0), Math.hypot(x2 - x0, y2 - y0))
        : Math.max(
            Math.abs((x1 - x0) * dy - (y1 - y0) * dx) / chord,
            Math.abs((x2 - x0) * dy - (y2 - y0) * dx) / chord
          );
    if (depth >= 20 || dev <= tol || (minChord > 0 && chord <= minChord)) {
      out.push({ x: x3, y: y3 });
      return;
    }
    const x01 = (x0 + x1) / 2;
    const y01 = (y0 + y1) / 2;
    const x12 = (x1 + x2) / 2;
    const y12 = (y1 + y2) / 2;
    const x23 = (x2 + x3) / 2;
    const y23 = (y2 + y3) / 2;
    const x012 = (x01 + x12) / 2;
    const y012 = (y01 + y12) / 2;
    const x123 = (x12 + x23) / 2;
    const y123 = (y12 + y23) / 2;
    const xm = (x012 + x123) / 2;
    const ym = (y012 + y123) / 2;
    flattenCubic(out, x0, y0, x01, y01, x012, y012, xm, ym, tol, depth + 1, minChord);
    flattenCubic(out, xm, ym, x123, y123, x23, y23, x3, y3, tol, depth + 1, minChord);
  }

  /**
   * `d` 属性を点列にする（変換前のローカル座標）。サブパス（M が複数）は最初の 1 本だけを
   * 返す — 呼び出し側が splitPathData で 1 サブパスずつ渡す契約のため（複数サブパスを
   * 1 本に繋いでグラフを壊した過去の不具合への対処。AGENT.md「複数サブパスの誤結線」）。
   * 未対応コマンド（A 円弧）が現れたら null を返し、呼び出し側が旧方式へ退避できるようにする。
   */
  function pathDataToPoints(d, tol, minChord) {
    const tokens = String(d || "").match(
      /[MmLlHhVvCcSsQqTtAaZz]|[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?/g
    );
    if (!tokens) return null;
    const pts = [];
    let i = 0;
    let cmd = "";
    let x = 0;
    let y = 0;
    let startX = 0;
    let startY = 0;
    // 直前の制御点（S / T の「反射」に使う）。**直前が同種の曲線でなければ反射は行わず
    // 現在点をそのまま制御点とする**のが SVG の規則。ここを取り違えると、反射先が
    // (2×現在点 − 0) のような遠方になり、始点と終点は合っているのに途中だけ地図の外まで
    // 飛ぶ形になる（実測: 最大 1630 地図単位のずれ。新旧比較で検出した）。
    let prevCtrlX = 0;
    let prevCtrlY = 0;
    let prevWasCubic = false; // 直前が C/c/S/s か（S の反射が有効か）
    let prevWasQuad = false; // 直前が Q/q/T/t か（T の反射が有効か）
    let started = false;
    const num = () => {
      const v = Number(tokens[i]);
      i += 1;
      return Number.isFinite(v) ? v : 0;
    };
    const push = (nx, ny) => {
      pts.push({ x: nx, y: ny });
      x = nx;
      y = ny;
    };
    while (i < tokens.length) {
      const token = tokens[i];
      // コマンドは 1 文字。`/[A-Za-z]/` で判定すると指数表記の数値（"1e1"）の e に反応して、
      // 暗黙の繰り返しの先頭に来た指数表記の数値をコマンドと取り違える（読み飛ばして形が欠ける）
      if (token.length === 1 && /[A-Za-z]/.test(token)) {
        // 2 本目以降のサブパスに入ったら打ち切る（1 サブパス = 1 点列の契約）
        if (started && (token === "M" || token === "m")) break;
        cmd = token;
        i += 1;
      }
      if (cmd === "A" || cmd === "a") return null; // 円弧は未対応 → 旧方式へ退避
      switch (cmd) {
        case "M":
        case "m": {
          const nx = (cmd === "m" ? x : 0) + num();
          const ny = (cmd === "m" ? y : 0) + num();
          push(nx, ny);
          startX = nx;
          startY = ny;
          started = true;
          prevWasCubic = false;
          prevWasQuad = false;
          cmd = cmd === "m" ? "l" : "L"; // 以降の座標対は暗黙の LineTo
          break;
        }
        case "L":
        case "l": {
          const nx = (cmd === "l" ? x : 0) + num();
          const ny = (cmd === "l" ? y : 0) + num();
          push(nx, ny);
          prevWasCubic = false;
          prevWasQuad = false;
          break;
        }
        case "H":
        case "h": {
          push((cmd === "h" ? x : 0) + num(), y);
          prevWasCubic = false;
          prevWasQuad = false;
          break;
        }
        case "V":
        case "v": {
          push(x, (cmd === "v" ? y : 0) + num());
          prevWasCubic = false;
          prevWasQuad = false;
          break;
        }
        case "C":
        case "c":
        case "S":
        case "s": {
          const rel = cmd === "c" || cmd === "s";
          const smooth = cmd === "S" || cmd === "s";
          const bx = rel ? x : 0;
          const by = rel ? y : 0;
          // S/s の第 1 制御点は「直前の第 2 制御点を現在点で反射したもの」。ただし
          // 直前が三次曲線でなければ反射せず現在点を使う（SVG 仕様）。
          const x1 = smooth ? (prevWasCubic ? 2 * x - prevCtrlX : x) : bx + num();
          const y1 = smooth ? (prevWasCubic ? 2 * y - prevCtrlY : y) : by + num();
          const x2 = bx + num();
          const y2 = by + num();
          const x3 = bx + num();
          const y3 = by + num();
          if (!started) {
            push(x, y);
            started = true;
          }
          flattenCubic(pts, x, y, x1, y1, x2, y2, x3, y3, tol, 0, minChord);
          prevCtrlX = x2;
          prevCtrlY = y2;
          prevWasCubic = true;
          prevWasQuad = false;
          x = x3;
          y = y3;
          break;
        }
        case "Q":
        case "q":
        case "T":
        case "t": {
          const rel = cmd === "q" || cmd === "t";
          const smooth = cmd === "T" || cmd === "t";
          const bx = rel ? x : 0;
          const by = rel ? y : 0;
          // T/t も同様に、直前が二次曲線でなければ反射せず現在点を制御点にする。
          const cx = smooth ? (prevWasQuad ? 2 * x - prevCtrlX : x) : bx + num();
          const cy = smooth ? (prevWasQuad ? 2 * y - prevCtrlY : y) : by + num();
          const x3 = bx + num();
          const y3 = by + num();
          // 二次 → 三次へ昇格して同じ分割器に載せる（判定を 1 箇所に保つ）
          const x1 = x + (2 / 3) * (cx - x);
          const y1 = y + (2 / 3) * (cy - y);
          const x2 = x3 + (2 / 3) * (cx - x3);
          const y2 = y3 + (2 / 3) * (cy - y3);
          if (!started) {
            push(x, y);
            started = true;
          }
          flattenCubic(pts, x, y, x1, y1, x2, y2, x3, y3, tol, 0, minChord);
          prevCtrlX = cx;
          prevCtrlY = cy;
          prevWasQuad = true;
          prevWasCubic = false;
          x = x3;
          y = y3;
          break;
        }
        case "Z":
        case "z": {
          if (started) push(startX, startY);
          prevWasCubic = false;
          prevWasQuad = false;
          cmd = "";
          break;
        }
        default:
          i += 1;
          break;
      }
    }
    if (pts.length < 2) return null;
    // 最小間隔の間引き。曲線の分割だけでなく、**元データの頂点そのもの**が吸着許容値より
    // 近接していることがある（Illustrator 出力は 1 単位未満の細かい c コマンドを多用する）。
    // それを残すと吸着で併合され、小さな島が面として消える（下限を分割にだけ掛けた版でも
    // 面 701 → 679・島国 17 個消失が残った）。旧方式は間引き間隔ちょうどの等間隔に
    // 打ち直していたので、同じ最小間隔をここで保証する。始点と終点は必ず残す（閉路の維持）。
    if (minChord > 0 && pts.length > 2) {
      const thinned = [pts[0]];
      for (let k = 1; k < pts.length - 1; k += 1) {
        const last = thinned[thinned.length - 1];
        if (Math.hypot(pts[k].x - last.x, pts[k].y - last.y) >= minChord) thinned.push(pts[k]);
      }
      thinned.push(pts[pts.length - 1]);
      return thinned.length >= 2 ? thinned : pts;
    }
    return pts;
  }

  /**
   * 新方式: `d` を直接解釈して点列にする。
   * step（間引き間隔）は曲線分割の下限に使う — 旧方式と同じ粗さを下回らせないため
   * （flattenCubic の minChord 注記を参照。省略時は下限なし＝純粋な精度優先）。
   */
  function samplePathElementFromData(pathEl, pathDataOverride, step) {
    if (!measureSvg) {
      throw new Error("measureSvg is not set.");
    }
    const minChord = Number.isFinite(step) && step > 0 ? step : 0;
    const d = pathDataOverride || pathEl.getAttribute("d") || "";
    const local = pathDataToPoints(d, curveTolerance, minChord);
    if (!local) return null;
    // 変換行列は旧方式と同じ経路で取る（measureSvg 直下のクローンの CTM）ので、
    // 座標系は完全に一致する（AGENT.md「CTM は measureSvg にクローンして取得」）。
    const clone = pathEl.cloneNode(true);
    if (pathDataOverride) clone.setAttribute("d", pathDataOverride);
    measureSvg.appendChild(clone);
    const points = local.map((p) => applyTransform(p, clone));
    measureSvg.removeChild(clone);
    return points;
  }

  function samplePathElement(pathEl, step, pathDataOverride) {
    if (sampleMode === SAMPLE_MODE_PATHDATA) {
      const points = samplePathElementFromData(pathEl, pathDataOverride, step);
      if (points) return points;
      // 未対応コマンドを含むパスだけ旧方式へ退避する（無音で形を落とさない）
    }
    return samplePathElementByLength(pathEl, step, pathDataOverride);
  }

  function samplePolylineElement(polyEl) {
    if (!measureSvg) {
      throw new Error("measureSvg is not set.");
    }
    const clone = polyEl.cloneNode(true);
    measureSvg.appendChild(clone);
    const points = [];
    const raw = clone.getAttribute("points") || "";
    const parts = raw.trim().split(/[\s,]+/).filter(Boolean);
    for (let i = 0; i < parts.length; i += 2) {
      const x = Number(parts[i]);
      const y = Number(parts[i + 1]);
      if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
      points.push(applyTransform({ x, y }, clone));
    }
    measureSvg.removeChild(clone);
    return points;
  }

  function sampleLineElement(lineEl) {
    if (!measureSvg) {
      throw new Error("measureSvg is not set.");
    }
    const clone = lineEl.cloneNode(true);
    measureSvg.appendChild(clone);
    const x1 = Number(clone.getAttribute("x1"));
    const y1 = Number(clone.getAttribute("y1"));
    const x2 = Number(clone.getAttribute("x2"));
    const y2 = Number(clone.getAttribute("y2"));
    const pts = [
      applyTransform({ x: x1, y: y1 }, clone),
      applyTransform({ x: x2, y: y2 }, clone),
    ];
    measureSvg.removeChild(clone);
    return pts;
  }

  /**
   * `<rect>` の 4 隅（変換後の座標）。回転を含む transform でも 4 隅を個別に写すので形は保たれる。
   * index.html（Wikipedia風出力の陸地取得）が使う。点列化の部品はすべてここに置き、
   * index.html に写しを置かない（AGENT.md「重複定義の取り残しで〜」）。
   */
  function sampleRectElement(rectEl) {
    if (!measureSvg) {
      throw new Error("measureSvg is not set.");
    }
    const clone = rectEl.cloneNode(true);
    measureSvg.appendChild(clone);
    const x = Number(clone.getAttribute("x")) || 0;
    const y = Number(clone.getAttribute("y")) || 0;
    const w = Number(clone.getAttribute("width")) || 0;
    const h = Number(clone.getAttribute("height")) || 0;
    const rawPoints = [
      { x, y },
      { x: x + w, y },
      { x: x + w, y: y + h },
      { x, y: y + h },
    ];
    const points = rawPoints.map((pt) => applyTransform(pt, clone));
    measureSvg.removeChild(clone);
    return points;
  }

  function distance(a, b) {
    const dx = a.x - b.x;
    const dy = a.y - b.y;
    return Math.hypot(dx, dy);
  }

  function distanceSq(a, b) {
    const dx = a.x - b.x;
    const dy = a.y - b.y;
    return dx * dx + dy * dy;
  }

  function removeSequentialDuplicates(points, eps) {
    if (points.length === 0) return points;
    const output = [points[0]];
    for (let i = 1; i < points.length; i += 1) {
      if (distance(points[i], output[output.length - 1]) > eps) {
        output.push(points[i]);
      }
    }
    return output;
  }

  async function extractSegmentsFromLayer(layerEl, step, sourceTag, onProgress, options) {
    const segments = [];
    const elements = Array.from(layerEl.querySelectorAll("path, polyline, polygon, line"));
    const total = elements.length;
    const yieldEvery = 20;
    const isDashed = options && typeof options.isDashed === "function" ? options.isDashed : null;
    const dashSnap = options && Number.isFinite(options.dashSnap) ? options.dashSnap : null;
    const groupedPaths = new Set();
    const pathGroupMap = new Map();

    {
      const groupElements = Array.from(layerEl.querySelectorAll("g"));
      groupElements.forEach((groupEl) => {
        if (groupEl.closest && groupEl.closest("defs")) return;
        const children = Array.from(groupEl.children || []);
        if (!children.length) return;
        const pathChildren = children.filter((child) => child.tagName && child.tagName.toLowerCase() === "path");
        if (pathChildren.length < 2 || pathChildren.length !== children.length) return;
        const info = {
          paths: pathChildren,
          processed: false,
        };
        pathChildren.forEach((child) => pathGroupMap.set(child, info));
      });
    }

    const pushSegments = (points, snapOverride) => {
      const cleaned = removeSequentialDuplicates(points, 1e-6);
      for (let i = 0; i < cleaned.length - 1; i += 1) {
        const segment = {
          a: cleaned[i],
          b: cleaned[i + 1],
          source: sourceTag,
        };
        if (Number.isFinite(snapOverride)) {
          segment.snap = snapOverride;
        }
        segments.push(segment);
      }
    };

    const collectSinglePathPoints = (pathEl) => {
      const d = pathEl.getAttribute("d") || "";
      const subpaths = splitPathData(d);
      if (subpaths.length !== 1) return null;
      return samplePathElement(pathEl, step, subpaths[0]);
    };

    const pathLength = (points) => {
      let totalLength = 0;
      for (let i = 0; i < points.length - 1; i += 1) {
        totalLength += distance(points[i], points[i + 1]);
      }
      return totalLength;
    };

    const buildGroupChain = (pathsPoints) => {
      const count = pathsPoints.length;
      if (!count) return null;
      if (count === 1) return pathsPoints[0];

      const oriented = pathsPoints.map((points) => {
        const reversed = points.slice().reverse();
        return {
          forward: points,
          reverse: reversed,
          start: points[0],
          end: points[points.length - 1],
          rStart: reversed[0],
          rEnd: reversed[reversed.length - 1],
          length: pathLength(points),
        };
      });

      const maxExact = 8;
      if (count > maxExact) {
        const used = new Array(count).fill(false);
        let startIndex = 0;
        let maxLength = -Infinity;
        for (let i = 0; i < count; i += 1) {
          if (oriented[i].length > maxLength) {
            maxLength = oriented[i].length;
            startIndex = i;
          }
        }
        used[startIndex] = true;
        const chain = [{ idx: startIndex, reverse: false }];
        let currentEnd = oriented[startIndex].end;
        while (chain.length < count) {
          let best = null;
          for (let i = 0; i < count; i += 1) {
            if (used[i]) continue;
            const forwardCost = distance(currentEnd, oriented[i].start);
            if (!best || forwardCost < best.cost) {
              best = { idx: i, reverse: false, cost: forwardCost };
            }
            const reverseCost = distance(currentEnd, oriented[i].rStart);
            if (reverseCost < best.cost) {
              best = { idx: i, reverse: true, cost: reverseCost };
            }
          }
          if (!best) break;
          used[best.idx] = true;
          chain.push({ idx: best.idx, reverse: best.reverse });
          currentEnd = best.reverse ? oriented[best.idx].rEnd : oriented[best.idx].end;
        }

        const combined = [];
        chain.forEach((entry) => {
          const points = entry.reverse ? oriented[entry.idx].reverse : oriented[entry.idx].forward;
          combined.push(...points);
        });
        return combined;
      }

      const used = new Array(count).fill(false);
      let bestCost = Infinity;
      let bestChain = null;

      const dfs = (lastEnd, chain, cost) => {
        if (chain.length === count) {
          if (cost < bestCost) {
            bestCost = cost;
            bestChain = chain.slice();
          }
          return;
        }
        for (let i = 0; i < count; i += 1) {
          if (used[i]) continue;
          const forwardCost = cost + distance(lastEnd, oriented[i].start);
          if (forwardCost < bestCost) {
            used[i] = true;
            chain.push({ idx: i, reverse: false });
            dfs(oriented[i].end, chain, forwardCost);
            chain.pop();
            used[i] = false;
          }
          const reverseCost = cost + distance(lastEnd, oriented[i].rStart);
          if (reverseCost < bestCost) {
            used[i] = true;
            chain.push({ idx: i, reverse: true });
            dfs(oriented[i].rEnd, chain, reverseCost);
            chain.pop();
            used[i] = false;
          }
        }
      };

      for (let i = 0; i < count; i += 1) {
        used[i] = true;
        dfs(oriented[i].end, [{ idx: i, reverse: false }], 0);
        dfs(oriented[i].rEnd, [{ idx: i, reverse: true }], 0);
        used[i] = false;
      }

      if (!bestChain) return null;
      const combined = [];
      bestChain.forEach((entry) => {
        const points = entry.reverse ? oriented[entry.idx].reverse : oriented[entry.idx].forward;
        combined.push(...points);
      });
      return combined;
    };

    const tryProcessPathGroup = (groupInfo) => {
      if (!groupInfo || groupInfo.paths.length < 2) return false;
      // この一本化は「点線国境が 1 本の線を複数 path に分けて描いてある」のを繋ぎ直すための
      // ものなので、国境レイヤーにだけ効かせる。海岸線に効かせると、独立した島の輪郭どうしが
      // 直線で結ばれて存在しない海岸線ができる（実データの _第一分高面 にも 4 path のグループが
      // あり、繋ぐと 254〜396 地図単位＝世界幅の 6〜9% の人工線ができる形だった。たまたま
      // 1 本がサブパス複数で中止されていただけで、条件が揃えば静かに壊れる）。
      if (sourceTag !== "border") return false;
      const pointsList = [];
      for (let i = 0; i < groupInfo.paths.length; i += 1) {
        const pathEl = groupInfo.paths[i];
        const points = collectSinglePathPoints(pathEl);
        if (!points || points.length < 2) return false;
        pointsList.push(points);
      }

      const combined = buildGroupChain(pointsList);
      if (!combined || combined.length < 2) return false;
      // Dashed class names vary, so isDashed() stays the single source of truth here.
      const snapOverride =
        isDashed && Number.isFinite(dashSnap) && groupInfo.paths.some((pathEl) => isDashed(pathEl)) ? dashSnap : null;
      pushSegments(combined, snapOverride);
      return true;
    };

    for (let index = 0; index < total; index += 1) {
      const item = elements[index];
      const tag = item.tagName.toLowerCase();
      if (groupedPaths.has(item)) {
        if (onProgress && (index % yieldEvery === 0 || index === total - 1)) {
          onProgress(index + 1, total);
          await nextFrame();
        }
        continue;
      }
      if (tag === "path") {
        const groupInfo = pathGroupMap.get(item);
        if (groupInfo && !groupInfo.processed) {
          const grouped = tryProcessPathGroup(groupInfo);
          groupInfo.processed = true;
          if (grouped) {
            groupInfo.paths.forEach((pathEl) => groupedPaths.add(pathEl));
            if (onProgress && (index % yieldEvery === 0 || index === total - 1)) {
              onProgress(index + 1, total);
              await nextFrame();
            }
            continue;
          }
        }
      }
      const dashed = isDashed ? isDashed(item) : false;
      const snapOverride = dashed && Number.isFinite(dashSnap) ? dashSnap : null;
      if (tag === "path") {
        const d = item.getAttribute("d") || "";
        const subpaths = splitPathData(d);
        if (subpaths.length === 1) {
          const points = samplePathElement(item, step, subpaths[0]);
          pushSegments(points, snapOverride);
        } else if (subpaths.length > 1) {
          subpaths.forEach((subpath) => {
            const points = samplePathElement(item, step, subpath);
            pushSegments(points, snapOverride);
          });
        }
      } else if (tag === "polyline" || tag === "polygon") {
        const points = samplePolylineElement(item);
        if (tag === "polygon" && points.length > 2) {
          points.push(points[0]);
        }
        pushSegments(points, snapOverride);
      } else if (tag === "line") {
        const points = sampleLineElement(item);
        pushSegments(points, snapOverride);
      }

      if (onProgress && (index % yieldEvery === 0 || index === total - 1)) {
        onProgress(index + 1, total);
        await nextFrame();
      }
    }
    return segments;
  }

  function lerpPoint(a, b, t) {
    return {
      x: a.x + (b.x - a.x) * t,
      y: a.y + (b.y - a.y) * t,
    };
  }

  function clamp01(value) {
    if (value < 0) return 0;
    if (value > 1) return 1;
    return value;
  }

  function paramAlongSegment(point, a, b) {
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    if (Math.abs(dx) >= Math.abs(dy)) {
      if (Math.abs(dx) < EPS) return 0;
      return (point.x - a.x) / dx;
    }
    if (Math.abs(dy) < EPS) return 0;
    return (point.y - a.y) / dy;
  }

  function pointSegmentProjection(point, a, b) {
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const len2 = dx * dx + dy * dy;
    if (len2 < EPS) return null;
    const t = ((point.x - a.x) * dx + (point.y - a.y) * dy) / len2;
    if (t < 0 || t > 1) return null;
    const proj = {
      x: a.x + dx * t,
      y: a.y + dy * t,
    };
    return { t, distance: distance(point, proj) };
  }

  function segmentIntersections(p1, p2, p3, p4) {
    const dx1 = p2.x - p1.x;
    const dy1 = p2.y - p1.y;
    const dx2 = p4.x - p3.x;
    const dy2 = p4.y - p3.y;
    const denom = dx1 * dy2 - dy1 * dx2;

    if (Math.abs(denom) >= EPS) {
      const dx3 = p1.x - p3.x;
      const dy3 = p1.y - p3.y;
      const t1 = (dx2 * dy3 - dy2 * dx3) / denom;
      const t2 = (dx1 * dy3 - dy1 * dx3) / denom;
      if (t1 < -EPS || t1 > 1 + EPS || t2 < -EPS || t2 > 1 + EPS) return null;
      return [{ t1: clamp01(t1), t2: clamp01(t2) }];
    }

    const cross = (p3.x - p1.x) * dy1 - (p3.y - p1.y) * dx1;
    if (Math.abs(cross) > EPS) return null;

    const len1 = Math.hypot(dx1, dy1);
    const len2 = Math.hypot(dx2, dy2);
    if (len1 < EPS && len2 < EPS) {
      return [{ t1: 0, t2: 0 }];
    }
    if (len1 < EPS) {
      const t2 = paramAlongSegment(p1, p3, p4);
      if (t2 < -EPS || t2 > 1 + EPS) return null;
      return [{ t1: 0, t2: clamp01(t2) }];
    }
    if (len2 < EPS) {
      const t1 = paramAlongSegment(p3, p1, p2);
      if (t1 < -EPS || t1 > 1 + EPS) return null;
      return [{ t1: clamp01(t1), t2: 0 }];
    }

    const tA3 = paramAlongSegment(p3, p1, p2);
    const tA4 = paramAlongSegment(p4, p1, p2);
    const minT = Math.min(tA3, tA4);
    const maxT = Math.max(tA3, tA4);
    const start = Math.max(minT, 0);
    const end = Math.min(maxT, 1);
    if (end < start - EPS) return null;

    const hits = [];
    const addHit = (tA) => {
      const point = lerpPoint(p1, p2, tA);
      const tB = paramAlongSegment(point, p3, p4);
      hits.push({ t1: clamp01(tA), t2: clamp01(tB) });
    };

    if (Math.abs(start - end) <= EPS) {
      addHit(start);
      return hits;
    }

    addHit(start);
    addHit(end);
    return hits;
  }

  async function splitSegmentsAtIntersections(segments, gridSize, endpointSnap, onProgress) {
    let baseSnap = Number(endpointSnap) || 0;
    let progress = onProgress;
    if (typeof endpointSnap === "function") {
      progress = endpointSnap;
      baseSnap = 0;
    }
    const splitEps = 1e-6;
    const segIntersections = new Map();
    const grid = new Map();

    function gridKey(x, y) {
      return `${x},${y}`;
    }

    function segmentSnap(seg) {
      const segSnap = seg && Number.isFinite(seg.snap) ? seg.snap : baseSnap;
      return segSnap > 0 ? segSnap : 0;
    }

    function addToGrid(idx, seg) {
      const pad = segmentSnap(seg);
      const minX = Math.min(seg.a.x, seg.b.x) - pad;
      const maxX = Math.max(seg.a.x, seg.b.x) + pad;
      const minY = Math.min(seg.a.y, seg.b.y) - pad;
      const maxY = Math.max(seg.a.y, seg.b.y) + pad;
      const x0 = Math.floor(minX / gridSize);
      const x1 = Math.floor(maxX / gridSize);
      const y0 = Math.floor(minY / gridSize);
      const y1 = Math.floor(maxY / gridSize);
      for (let gx = x0; gx <= x1; gx += 1) {
        for (let gy = y0; gy <= y1; gy += 1) {
          const key = gridKey(gx, gy);
          if (!grid.has(key)) grid.set(key, []);
          grid.get(key).push(idx);
        }
      }
    }

    function addSplit(idx, t) {
      if (!segIntersections.has(idx)) segIntersections.set(idx, []);
      segIntersections.get(idx).push(t);
    }

    function checkEndpointSnap(point, targetSeg, targetIdx, snapDist) {
      if (!(snapDist > 0)) return;
      const projection = pointSegmentProjection(point, targetSeg.a, targetSeg.b);
      if (!projection) return;
      if (projection.distance <= snapDist) {
        addSplit(targetIdx, projection.t);
      }
    }

    const gridYield = 800;
    for (let i = 0; i < segments.length; i += 1) {
      addToGrid(i, segments[i]);
      if (progress && (i % gridYield === 0 || i === segments.length - 1)) {
        progress("grid", i + 1, segments.length);
        await nextFrame();
      }
    }

    const checked = new Set();
    const buckets = Array.from(grid.values());
    const bucketYield = 40;
    for (let bucketIndex = 0; bucketIndex < buckets.length; bucketIndex += 1) {
      const bucket = buckets[bucketIndex];
      for (let i = 0; i < bucket.length; i += 1) {
        for (let j = i + 1; j < bucket.length; j += 1) {
          const a = bucket[i];
          const b = bucket[j];
          const key = a < b ? `${a}-${b}` : `${b}-${a}`;
          if (checked.has(key)) continue;
          checked.add(key);
          const segA = segments[a];
          const segB = segments[b];
          const hits = segmentIntersections(segA.a, segA.b, segB.a, segB.b);
          if (hits) {
            for (let h = 0; h < hits.length; h += 1) {
              const hit = hits[h];
              addSplit(a, hit.t1);
              addSplit(b, hit.t2);
            }
          }
          const snapA = segmentSnap(segA);
          const snapB = segmentSnap(segB);
          if (snapA > 0) {
            checkEndpointSnap(segA.a, segB, b, snapA);
            checkEndpointSnap(segA.b, segB, b, snapA);
          }
          if (snapB > 0) {
            checkEndpointSnap(segB.a, segA, a, snapB);
            checkEndpointSnap(segB.b, segA, a, snapB);
          }
        }
      }
      if (progress && (bucketIndex % bucketYield === 0 || bucketIndex === buckets.length - 1)) {
        progress("intersections", bucketIndex + 1, buckets.length);
        await nextFrame();
      }
    }

    const output = [];
    const splitYield = 1000;
    for (let idx = 0; idx < segments.length; idx += 1) {
      const seg = segments[idx];
      const splits = segIntersections.get(idx) || [];
      const candidates = splits.filter((t) => t > splitEps && t < 1 - splitEps);
      const sorted = [0, ...candidates, 1].sort((a, b) => a - b);
      const unique = [];
      for (let i = 0; i < sorted.length; i += 1) {
        const t = sorted[i];
        if (unique.length === 0 || Math.abs(t - unique[unique.length - 1]) > splitEps) {
          unique.push(t);
        }
      }
      for (let i = 0; i < unique.length - 1; i += 1) {
        const t0 = unique[i];
        const t1 = unique[i + 1];
        const p0 = lerpPoint(seg.a, seg.b, t0);
        const p1 = lerpPoint(seg.a, seg.b, t1);
        if (distance(p0, p1) > splitEps) {
          const nextSeg = { a: p0, b: p1, source: seg.source };
          if (Number.isFinite(seg.snap)) {
            nextSeg.snap = seg.snap;
          }
          output.push(nextSeg);
        }
      }
      if (progress && (idx % splitYield === 0 || idx === segments.length - 1)) {
        progress("split", idx + 1, segments.length);
        await nextFrame();
      }
    }

    return output;
  }

  async function buildGraphFromSegments(segments, snap, onProgress) {
    const pointGrid = new Map();
    const points = [];
    const pointHits = [];
    const pointSnaps = [];
    const edges = [];
    const edgeMap = new Map();
    const baseSnap = Number(snap) || 0;
    let maxSnap = baseSnap;
    for (let i = 0; i < segments.length; i += 1) {
      const segSnap = Number.isFinite(segments[i].snap) ? segments[i].snap : baseSnap;
      if (segSnap > maxSnap) maxSnap = segSnap;
    }
    const cellSize = maxSnap > 0 ? maxSnap : 1e-6;

    function gridKey(x, y) {
      return `${x}:${y}`;
    }

    function addPoint(pt, snapHint) {
      const snapDist = Number.isFinite(snapHint) ? snapHint : baseSnap;
      const cx = Math.floor(pt.x / cellSize);
      const cy = Math.floor(pt.y / cellSize);
      let bestIdx = null;
      let bestDist = Infinity;

      for (let dx = -1; dx <= 1; dx += 1) {
        for (let dy = -1; dy <= 1; dy += 1) {
          const key = gridKey(cx + dx, cy + dy);
          const candidates = pointGrid.get(key);
          if (!candidates) continue;
          for (let i = 0; i < candidates.length; i += 1) {
            const idx = candidates[i];
            const d = distance(pt, points[idx]);
            const allowed = Math.max(snapDist, pointSnaps[idx] || 0);
            if (d <= allowed && d < bestDist) {
              bestDist = d;
              bestIdx = idx;
            }
          }
        }
      }

      if (bestIdx === null) {
        const idx = points.length;
        points.push({ x: pt.x, y: pt.y });
        pointHits.push(1);
        pointSnaps.push(snapDist);
        const key = gridKey(cx, cy);
        if (!pointGrid.has(key)) pointGrid.set(key, []);
        pointGrid.get(key).push(idx);
        return idx;
      }

      const hits = pointHits[bestIdx] + 1;
      points[bestIdx].x = (points[bestIdx].x * pointHits[bestIdx] + pt.x) / hits;
      points[bestIdx].y = (points[bestIdx].y * pointHits[bestIdx] + pt.y) / hits;
      pointHits[bestIdx] = hits;
      pointSnaps[bestIdx] = Math.max(pointSnaps[bestIdx] || 0, snapDist);
      return bestIdx;
    }

    const buildYield = 2000;
    for (let i = 0; i < segments.length; i += 1) {
      const seg = segments[i];
      const segSnap = Number.isFinite(seg.snap) ? seg.snap : baseSnap;
      const a = addPoint(seg.a, segSnap);
      const b = addPoint(seg.b, segSnap);
      if (a === b) continue;
      const key = a < b ? `${a}-${b}` : `${b}-${a}`;
      let edge = edgeMap.get(key);
      if (!edge) {
        edge = { a, b, sources: new Set() };
        edgeMap.set(key, edge);
        edges.push(edge);
      }
      edge.sources.add(seg.source);

      if (onProgress && (i % buildYield === 0 || i === segments.length - 1)) {
        onProgress(i + 1, segments.length);
        await nextFrame();
      }
    }

    return { points, edges };
  }

  // --- 閉じ損ねた国境の橋渡し -------------------------------------------------
  //
  // 原図の国境線には、作図上の理由で「わずかに届いていない」隙間が実在する。吸着許容値
  // (既定 1) をほんの少し超える程度の隙間でも、面探索から見れば線が途切れているのと同じで、
  // その国は隣国と繋がった 1 つの面になってしまう（中央イノン/東イノン間で 0.95〜、
  // バアカ/コルニア南方の海上国境では 2.6〜5.7 の隙間が多数）。
  //
  // そこで、行き止まり（次数 1）の端点を近傍の点へ繋いで面を閉じる。吸着許容値そのものを
  // 上げて解決しないのは、吸着は全ての点に一律に効くため、小さな島の輪郭が潰れて面ごと
  // 消えるなどの副作用が大きすぎるから（AGENT.md「点の密度と吸着許容値は連動している」）。
  // 橋渡しは「行き止まりの端点」という、本来閉じているはずの箇所だけを狙って直す。
  //
  // 半径が 2 種類あるのは、海の中で自由に終わる共有海式の国境線が、近くを通る海岸線に
  // 誤って吸い付くのを避けるため。大きい半径 (borderRadius) は「自分も相手も国境線の点」
  // のときだけ許す。

  /**
   * 個別指定の橋渡し（このデーレ界の地図に固有の 1 箇所）。
   * ラヴィエンド・ウィル西方の岬: 点線国境の終端が岬の先端にあり、最寄りの海岸線
   * サンプル点まで 4.41 ある（地図の意図は海岸への接続）。座標 ±0.5 に一致する
   * 行き止まりだけ、この半径を許す。他の地図では options.extraBridges で差し替える。
   */
  const DEFAULT_EXTRA_BRIDGES = [{ x: 783.4, y: 739.8, r: 5.0 }];

  function orient2d(a, b, c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  }

  /**
   * 点 p が線分 ab の上（両端を除く）に乗っているか。
   * 橋渡しが既存の頂点を途中で踏むと、そこにノードが無いままグラフが非平面になる。
   * 「真に交差する」判定だけではこれを取りこぼすので別に見る。
   */
  // 「乗っている」と見なす距離。地図の座標は 0〜4200 程度で、吸着で併合されなかった点は
  // 通常 1 以上離れているので、これを超えて緩めると正当な橋渡しまで弾いてしまう。
  const COLLINEAR_EPS = 1e-6;

  function pointOnOpenSegment(p, a, b, eps) {
    const vx = b.x - a.x;
    const vy = b.y - a.y;
    const len2 = vx * vx + vy * vy;
    if (len2 <= 0) return false;
    const len = Math.sqrt(len2);
    // orient2d は |ab| × 点と直線の距離。距離が eps を超えていれば乗っていない
    if (Math.abs(orient2d(a, b, p)) > eps * len) return false;
    const t = ((p.x - a.x) * vx + (p.y - a.y) * vy) / len2;
    return t * len > eps && (1 - t) * len > eps;
  }

  /** 端点を共有しない 2 線分が真に交差するか（端点で触れているだけは交差としない）。 */
  function properlyCrosses(p1, p2, p3, p4) {
    const d1 = orient2d(p3, p4, p1);
    const d2 = orient2d(p3, p4, p2);
    const d3 = orient2d(p1, p2, p3);
    const d4 = orient2d(p1, p2, p4);
    return (
      ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
      ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))
    );
  }

  /**
   * 行き止まりの国境端点を近傍へ繋いで、閉じ損ねた面を閉じる。
   * graph（{ points, edges }）を**その場で**書き換え、追加した橋の一覧を返す。
   *
   * options:
   *   radius       … 通常の橋渡し半径（既定 2.5）
   *   borderRadius … 国境線どうしの橋渡し半径（既定 6.0）。自分も相手も国境線の点のときだけ使う
   *   cellSize     … 交差判定に使う一様グリッドのセル幅（既定 8）
   *   extraBridges … 座標指定で半径を緩める例外（既定 DEFAULT_EXTRA_BRIDGES）
   *   onProgress   … (done, total) を受ける進捗コールバック。指定時は合間に rAF を挟む
   *
   * 戻り値 { bridges: [{ x, y, d, tip, target }], rejected }
   *   rejected … 半径内に候補はあったが、全て交差するため張れなかった行き止まりの数
   */
  async function bridgeOpenBorderGaps(graph, options) {
    const opts = options || {};
    const radiusBase = Number.isFinite(opts.radius) ? opts.radius : 2.5;
    const borderRadiusBase = Number.isFinite(opts.borderRadius) ? opts.borderRadius : 6.0;
    const CELL = Number.isFinite(opts.cellSize) && opts.cellSize > 0 ? opts.cellSize : 8;
    const extraBridges = opts.extraBridges || DEFAULT_EXTRA_BRIDGES;
    const onProgress = typeof opts.onProgress === "function" ? opts.onProgress : null;

    const points = graph.points;
    const edges = graph.edges;
    const deg = new Array(points.length).fill(0);
    for (let i = 0; i < edges.length; i += 1) {
      deg[edges[i].a] += 1;
      deg[edges[i].b] += 1;
    }
    // 次数 1 の点（＝行き止まり）。順序は辺の並び順に揃える — 橋を 1 本張るたびに
    // 次数とグリッドが変わるので、処理順が違うと結果も変わりうる。移植元（Node 版の
    // extract-polygons.mjs）と同じ順序にして、実測値と突き合わせられるようにする。
    const tips = [];
    const tipSeen = new Set();
    for (let i = 0; i < edges.length; i += 1) {
      const e = edges[i];
      if (deg[e.a] === 1 && !tipSeen.has(e.a)) {
        tipSeen.add(e.a);
        tips.push(e.a);
      }
      if (deg[e.b] === 1 && !tipSeen.has(e.b)) {
        tipSeen.add(e.b);
        tips.push(e.b);
      }
    }

    // 隣接リスト（3 ホップ以内の自チェーン点を橋渡し先から除外するため。
    // すぐ隣の点に繋いでも隙間は閉じず、極小の三角形が増えるだけになる）
    const adjAll = new Map();
    for (let i = 0; i < edges.length; i += 1) {
      const e = edges[i];
      if (!adjAll.has(e.a)) adjAll.set(e.a, []);
      if (!adjAll.has(e.b)) adjAll.set(e.b, []);
      adjAll.get(e.a).push(e.b);
      adjAll.get(e.b).push(e.a);
    }

    const isBorderPoint = new Uint8Array(points.length);
    for (let i = 0; i < edges.length; i += 1) {
      const e = edges[i];
      if (e.sources && e.sources.has("border")) {
        isBorderPoint[e.a] = 1;
        isBorderPoint[e.b] = 1;
      }
    }

    // 平面性の維持（必須）: 橋渡し辺が既存の辺と交差すると、交点にノードが無いまま
    // グラフが非平面になり、面探索（左回りの角度順たどり）が破綻する。実測では
    // リングが自己交差して穴が親の外へはみ出す "hole-crosses-parent" が 6 件出た。
    // よって既存辺と真に交差する橋渡しは棄却する。
    // 総当たり（辺 10 万本 × 橋渡し候補）は重すぎるので、一様グリッドで近傍辺に絞る。
    const grid = new Map();
    const cellKey = (cx, cy) => `${cx},${cy}`;
    const addToGrid = (idx, a, b) => {
      const x0 = Math.floor(Math.min(a.x, b.x) / CELL);
      const x1 = Math.floor(Math.max(a.x, b.x) / CELL);
      const y0 = Math.floor(Math.min(a.y, b.y) / CELL);
      const y1 = Math.floor(Math.max(a.y, b.y) / CELL);
      for (let cx = x0; cx <= x1; cx += 1) {
        for (let cy = y0; cy <= y1; cy += 1) {
          const k = cellKey(cx, cy);
          let list = grid.get(k);
          if (list === undefined) {
            list = [];
            grid.set(k, list);
          }
          list.push(idx);
        }
      }
    };
    for (let i = 0; i < edges.length; i += 1) {
      addToGrid(i, points[edges[i].a], points[edges[i].b]);
    }

    // 候補の相手を探すための点のグリッド。セル幅は「ありうる最大の半径」にしておくと、
    // 周囲 1 セルぶんを見るだけで半径内の点を必ず含む。
    let maxRadius = Math.max(radiusBase, borderRadiusBase);
    for (let i = 0; i < extraBridges.length; i += 1) {
      if (extraBridges[i].r > maxRadius) maxRadius = extraBridges[i].r;
    }
    const POINT_CELL = Math.max(maxRadius, 1e-6);
    const pointGrid = new Map();
    for (let i = 0; i < points.length; i += 1) {
      const k = cellKey(Math.floor(points[i].x / POINT_CELL), Math.floor(points[i].y / POINT_CELL));
      let list = pointGrid.get(k);
      if (list === undefined) {
        list = [];
        pointGrid.set(k, list);
      }
      list.push(i);
    }

    /** 橋渡し tip→target が既存辺と交差するなら true（＝採用不可）。 */
    const bridgeCrosses = (ti, tj) => {
      const p1 = points[ti];
      const p2 = points[tj];
      const x0 = Math.floor(Math.min(p1.x, p2.x) / CELL);
      const x1 = Math.floor(Math.max(p1.x, p2.x) / CELL);
      const y0 = Math.floor(Math.min(p1.y, p2.y) / CELL);
      const y1 = Math.floor(Math.max(p1.y, p2.y) / CELL);
      const seen = new Set();
      for (let cx = x0 - 1; cx <= x1 + 1; cx += 1) {
        for (let cy = y0 - 1; cy <= y1 + 1; cy += 1) {
          const list = grid.get(cellKey(cx, cy));
          if (list === undefined) continue;
          for (let n = 0; n < list.length; n += 1) {
            const ei = list[n];
            if (seen.has(ei)) continue;
            seen.add(ei);
            const e = edges[ei];
            const ea = points[e.a];
            const eb = points[e.b];
            const shared = e.a === ti || e.b === ti || e.a === tj || e.b === tj;
            // 端点を共有する辺は「真の交差」にはならない（橋の根元と先で必ず触れるため）。
            // ただし**免除してよいのは交差判定だけ**。共有していても、橋がその辺を
            // なぞって既存の頂点を途中で踏むことがある。
            //   例: 点 (0,0)(1,0)(5,0)、既存辺 0-1、橋 0→(5,0)
            //   → 橋は辺 0-1 に重なり点 (1,0) を踏むが、端点 0 を共有しているので
            //     以前は素通りしていた（外部レビュー指摘）。
            if (!shared && properlyCrosses(p1, p2, ea, eb)) return true;
            // 「真の交差」にはならないが平面性を壊す形。どちらかの端点が相手の
            // 線分の内側（両端を除く）に乗っていれば、そこにノードが無いまま辺が重なる。
            if (pointOnOpenSegment(ea, p1, p2, COLLINEAR_EPS)) return true;
            if (pointOnOpenSegment(eb, p1, p2, COLLINEAR_EPS)) return true;
            if (pointOnOpenSegment(p1, ea, eb, COLLINEAR_EPS)) return true;
            if (pointOnOpenSegment(p2, ea, eb, COLLINEAR_EPS)) return true;
          }
        }
      }
      return false;
    };

    const bridges = [];
    let rejected = 0;
    const yieldEvery = 50;
    for (let t = 0; t < tips.length; t += 1) {
      const tip = tips[t];
      if (onProgress && (t % yieldEvery === 0 || t === tips.length - 1)) {
        onProgress(t + 1, tips.length);
        await nextFrame();
      }
      if (deg[tip] !== 1) continue; // 既に他の橋渡しで次数が増えていたら対象外

      // 3 ホップ以内の自分側の点を候補から外す
      const near = new Set([tip]);
      let frontier = [tip];
      for (let hop = 0; hop < 3; hop += 1) {
        const next = [];
        for (let i = 0; i < frontier.length; i += 1) {
          const list = adjAll.get(frontier[i]);
          if (!list) continue;
          for (let j = 0; j < list.length; j += 1) {
            const w = list[j];
            if (!near.has(w)) {
              near.add(w);
              next.push(w);
            }
          }
        }
        frontier = next;
      }

      const p = points[tip];
      let radius = radiusBase;
      for (let i = 0; i < extraBridges.length; i += 1) {
        const ex = extraBridges[i];
        if (Math.abs(p.x - ex.x) <= 0.5 && Math.abs(p.y - ex.y) <= 0.5) radius = ex.r;
      }
      // 国境線の行き止まりだけ、国境線の点に対して大きい半径を許す
      const borderRadius = isBorderPoint[tip] === 1 ? Math.max(radius, borderRadiusBase) : radius;

      // 候補を距離順に並べ、既存辺と交差しない最も近い相手を採用する。
      // 半径は必ず POINT_CELL 以下なので、周囲 1 セルぶんを見れば取りこぼさない
      // （全点走査だと、行き止まりが大量にある入力で行き止まり数 × 点数まで膨らむ）。
      const candidates = [];
      const pcx = Math.floor(p.x / POINT_CELL);
      const pcy = Math.floor(p.y / POINT_CELL);
      for (let dx = -1; dx <= 1; dx += 1) {
        for (let dy = -1; dy <= 1; dy += 1) {
          const list = pointGrid.get(cellKey(pcx + dx, pcy + dy));
          if (list === undefined) continue;
          for (let n = 0; n < list.length; n += 1) {
            const j = list[n];
            if (near.has(j)) continue;
            const q = points[j];
            const d = Math.hypot(q.x - p.x, q.y - p.y);
            const limit = isBorderPoint[j] === 1 ? borderRadius : radius;
            if (d < limit) candidates.push({ j, d });
          }
        }
      }
      candidates.sort((u, v) => u.d - v.d);
      let best = -1;
      let bestD = 0;
      for (let i = 0; i < candidates.length; i += 1) {
        if (bridgeCrosses(tip, candidates[i].j)) continue;
        best = candidates[i].j;
        bestD = candidates[i].d;
        break;
      }
      if (best === -1) {
        if (candidates.length > 0) rejected += 1;
        continue;
      }

      const newIdx = edges.length;
      edges.push({ a: tip, b: best, sources: new Set(["border"]) });
      addToGrid(newIdx, points[tip], points[best]);
      deg[tip] += 1;
      deg[best] += 1;
      bridges.push({ x: p.x, y: p.y, d: bestD, tip, target: best });
    }

    return { bridges, rejected };
  }

  function polygonArea(indices, points) {
    let sum = 0;
    for (let i = 0; i < indices.length; i += 1) {
      const a = points[indices[i]];
      const b = points[indices[(i + 1) % indices.length]];
      sum += (a.x * b.y - b.x * a.y);
    }
    return sum / 2;
  }

  function extractPolygonsFromGraph(graph, minArea) {
    const { points, edges } = graph;
    const directed = [];
    const outgoing = Array.from({ length: points.length }, () => []);

    edges.forEach((edge, edgeIndex) => {
      const a = edge.a;
      const b = edge.b;
      const angleAB = Math.atan2(points[b].y - points[a].y, points[b].x - points[a].x);
      const angleBA = Math.atan2(points[a].y - points[b].y, points[a].x - points[b].x);
      const idxAB = directed.length;
      directed.push({ from: a, to: b, angle: angleAB, edgeIndex });
      const idxBA = directed.length;
      directed.push({ from: b, to: a, angle: angleBA, edgeIndex });
      outgoing[a].push(idxAB);
      outgoing[b].push(idxBA);
    });

    outgoing.forEach((list) => list.sort((i, j) => directed[i].angle - directed[j].angle));

    const adjacency = Array.from({ length: points.length }, () => []);
    edges.forEach((edge) => {
      adjacency[edge.a].push(edge.b);
      adjacency[edge.b].push(edge.a);
    });

    const components = new Array(points.length).fill(-1);
    let componentId = 0;
    for (let i = 0; i < points.length; i += 1) {
      if (components[i] !== -1) continue;
      if (adjacency[i].length === 0) continue;
      const stack = [i];
      components[i] = componentId;
      while (stack.length) {
        const v = stack.pop();
        const neighbors = adjacency[v];
        for (let j = 0; j < neighbors.length; j += 1) {
          const next = neighbors[j];
          if (components[next] === -1) {
            components[next] = componentId;
            stack.push(next);
          }
        }
      }
      componentId += 1;
    }

    function nextEdgeIndex(dirIndex) {
      const dir = directed[dirIndex];
      const v = dir.to;
      const incomingAngle = Math.atan2(points[dir.from].y - points[v].y, points[dir.from].x - points[v].x);
      const list = outgoing[v];
      if (!list || list.length === 0) return null;
      let chosen = list[list.length - 1];
      for (let i = 0; i < list.length; i += 1) {
        const candidate = list[i];
        if (directed[candidate].angle < incomingAngle - 1e-10) {
          chosen = candidate;
        } else {
          break;
        }
      }
      return chosen;
    }

    const visited = new Array(directed.length).fill(false);
    const faces = [];

    for (let i = 0; i < directed.length; i += 1) {
      if (visited[i]) continue;
      const face = [];
      let current = i;
      let guard = 0;
      while (!visited[current] && guard < directed.length + 10) {
        visited[current] = true;
        const edge = directed[current];
        face.push(edge.from);
        const next = nextEdgeIndex(current);
        if (next === null) break;
        current = next;
        guard += 1;
        if (current === i) break;
      }

      if (face.length < 3) continue;
      const area = polygonArea(face, points);
      if (Math.abs(area) < minArea) continue;
      faces.push({ points: face, area, component: components[face[0]] });
    }

    if (faces.length === 0) {
      return { polygons: [], outerFaces: [] };
    }

    const faceGroups = new Map();
    faces.forEach((face) => {
      const key = face.component;
      if (!faceGroups.has(key)) faceGroups.set(key, []);
      faceGroups.get(key).push(face);
    });

    const outerFaces = [];
    faceGroups.forEach((group) => {
      const largest = group.reduce((max, face) => (Math.abs(face.area) > Math.abs(max.area) ? face : max), group[0]);
      outerFaces.push(largest);
    });

    const outerSet = new Set(outerFaces);
    const polygons = faces.filter((face) => !outerSet.has(face)).map((face) => ({
      area: face.area,
      points: face.points,
    }));

    return { polygons, outerFaces };
  }

  // --- 領有関係（要件定義書 §6.3） ---------------------------------------------
  //
  // 仕様は「どの面が島々を領有するか」を面の種類ごとに定めている:
  //   - 海まで国境線を引く国 … 「国境線の内側にある島々も同時に領有する」
  //   - 共有する海を持つ国々 … 「島々は海領域に領有される」
  //   - 海岸線まで国境界を引く国 … 「島々は領有しない」（領海を主張しないため）
  //   - 内陸国 / 完全島国 … 内側に島を含む状況を想定しない
  //
  // ここは「面 → それを領有する面」の対応だけを返す純粋関数。実際にどう描くか・
  // どう出力するかは呼び出し側の責務（描画は index.html、書き出しは .gimoza 出力）。
  //
  // 未実装だった経緯: 出力が「点・辺・ポリゴン・種別」だけで領有を表せず、国境線で
  // 囲まれた島が独立した面として出ていた（実機で指摘）。

  /** 領有する側になれる面の種類。 */
  const OWNING_TYPES = new Set(["海まで国境線", "共有する海"]);

  function pointInRingIndices(x, y, ringIndices, points) {
    let inside = false;
    for (let i = 0, j = ringIndices.length - 1; i < ringIndices.length; j = i, i += 1) {
      const a = points[ringIndices[i]];
      const b = points[ringIndices[j]];
      if (a.y > y !== b.y > y && x < ((b.x - a.x) * (y - a.y)) / (b.y - a.y) + a.x) {
        inside = !inside;
      }
    }
    return inside;
  }

  function ringBounds(ringIndices, points) {
    let minX = Infinity;
    let minY = Infinity;
    let maxX = -Infinity;
    let maxY = -Infinity;
    for (let i = 0; i < ringIndices.length; i += 1) {
      const q = points[ringIndices[i]];
      if (q.x < minX) minX = q.x;
      if (q.y < minY) minY = q.y;
      if (q.x > maxX) maxX = q.x;
      if (q.y > maxY) maxY = q.y;
    }
    return { minX, minY, maxX, maxY };
  }

  /**
   * 各面を内包する「最も内側の面」を求める（内包の親）。
   *
   * 隣接する面どうしは境界の頂点を共有するので、共有していない頂点だけを標本にして
   * 多数決する（そうしないと「辺を接している」だけの隣国を親と誤判定する）。
   * 戻り値は index → index の配列（親が無ければ -1）。
   */
  function computeContainment(polygons, points) {
    const bounds = polygons.map((p) => ringBounds(p.points, points));
    const areas = polygons.map((p) => Math.abs(p.area));
    const parents = new Array(polygons.length).fill(-1);
    for (let i = 0; i < polygons.length; i += 1) {
      const bi = bounds[i];
      let best = -1;
      for (let j = 0; j < polygons.length; j += 1) {
        if (i === j || areas[j] <= areas[i]) continue;
        const bj = bounds[j];
        if (bj.minX > bi.minX || bj.minY > bi.minY || bj.maxX < bi.maxX || bj.maxY < bi.maxY) {
          continue;
        }
        const outerSet = new Set(polygons[j].points);
        const samples = [];
        for (let k = 0; k < polygons[i].points.length && samples.length < 9; k += 1) {
          const vi = polygons[i].points[k];
          if (!outerSet.has(vi)) samples.push(vi);
        }
        if (samples.length === 0) continue; // 全頂点共有 = 隣接であって内包ではない
        let inCount = 0;
        for (let k = 0; k < samples.length; k += 1) {
          const q = points[samples[k]];
          if (pointInRingIndices(q.x, q.y, polygons[j].points, points)) inCount += 1;
        }
        if (inCount * 2 <= samples.length) continue;
        if (best === -1 || areas[j] < areas[best]) best = j;
      }
      parents[i] = best;
    }
    return parents;
  }

  /**
   * 領有関係を求める。`types` は面 index → 種別文字列。
   * 戻り値 { parents, owners }:
   *   parents … 内包の親（幾何的な包含。種別に関係なく求まる）
   *   owners  … 領有する面（親の種別が領有する側のときだけ設定。無ければ -1）
   */
  function computeOwnership(polygons, points, types) {
    const parents = computeContainment(polygons, points);
    const owners = parents.map((parent, i) => {
      if (parent === -1) return -1;
      return OWNING_TYPES.has(types[parent]) ? parent : -1;
    });
    return { parents, owners };
  }

  // ---- 経度180度の継ぎ目（地図の西端・東端）----
  // frame: 正距円筒図法で経度-180〜+180・緯度+90〜-90に対応する矩形 {x, y, width, height}

  // 地図枠の外側にはみ出した点は、内側の許容値のこの倍数まで端に乗っているとみなす。
  // 作図が枠をわずかにはみ出すことがあり（BUR080809 では南極の切断線が枠の 0.07 外）、
  // 枠の外の点は端に揃えるしかないので、端に斜めに近づく海岸線を巻き込まないよう狭くしている内側より緩くしてよい。
  const SEAM_OUTSIDE_FACTOR = 20;

  // 点が地図枠の西端/東端のどちらに乗っているか（内側tol以内、外側tol×SEAM_OUTSIDE_FACTOR以内）。"W" / "E" / null
  function seamSideOf(point, frame, tol) {
    if (!point || !frame) return null;
    const west = frame.x;
    const east = frame.x + frame.width;
    const outside = tol * SEAM_OUTSIDE_FACTOR;
    if (point.x <= west + tol && point.x >= west - outside) return "W";
    if (point.x >= east - tol && point.x <= east + outside) return "E";
    return null;
  }

  // 西端/東端からtol以内の点のxを枠の端ちょうどに揃える（pointsを直接書き換える）。揃えた点数を返す。
  // indicesを渡すとその点だけを対象にする。
  function snapPointsToSeam(points, frame, tol, indices) {
    if (!points || !frame || !(tol >= 0)) return 0;
    const west = frame.x;
    const east = frame.x + frame.width;
    let count = 0;
    const visit = (p) => {
      if (!p) return;
      const side = seamSideOf(p, frame, tol);
      if (side === "W") {
        if (p.x !== west) count += 1;
        p.x = west;
      } else if (side === "E") {
        if (p.x !== east) count += 1;
        p.x = east;
      }
    };
    if (indices) {
      indices.forEach((i) => visit(points[i]));
    } else {
      points.forEach(visit);
    }
    return count;
  }

  // 両端点が同じ側の地図端に乗っている辺（地図端に沿う辺）か。"W" / "E" / null
  function seamEdgeSide(a, b, frame, tol) {
    const sa = seamSideOf(a, frame, tol);
    if (!sa) return null;
    return seamSideOf(b, frame, tol) === sa ? sa : null;
  }

  // ---- 区間（y範囲）演算: [[y0, y1], ...] を昇順・非重複で扱う ----
  function mergeIntervals(intervals, gap) {
    const list = intervals.filter((iv) => iv[1] > iv[0]).map((iv) => iv.slice()).sort((p, q) => p[0] - q[0]);
    if (!list.length) return [];
    const g = gap || 0;
    const out = [list[0]];
    for (let i = 1; i < list.length; i += 1) {
      const last = out[out.length - 1];
      const cur = list[i];
      if (cur[0] <= last[1] + g) {
        if (cur[1] > last[1]) last[1] = cur[1];
      } else {
        out.push(cur);
      }
    }
    return out;
  }

  function intersectIntervals(a, b) {
    const out = [];
    let i = 0;
    let j = 0;
    while (i < a.length && j < b.length) {
      const lo = Math.max(a[i][0], b[j][0]);
      const hi = Math.min(a[i][1], b[j][1]);
      if (hi > lo) out.push([lo, hi]);
      if (a[i][1] < b[j][1]) i += 1; else j += 1;
    }
    return out;
  }

  function subtractIntervals(a, b) {
    const out = [];
    a.forEach((iv) => {
      let pieces = [iv.slice()];
      b.forEach((cut) => {
        const next = [];
        pieces.forEach((p) => {
          if (cut[1] <= p[0] || cut[0] >= p[1]) {
            next.push(p);
            return;
          }
          if (cut[0] > p[0]) next.push([p[0], cut[0]]);
          if (cut[1] < p[1]) next.push([cut[1], p[1]]);
        });
        pieces = next;
      });
      pieces.forEach((p) => out.push(p));
    });
    return mergeIntervals(out, 0);
  }

  function intervalsContain(intervals, y) {
    for (let i = 0; i < intervals.length; i += 1) {
      if (y > intervals[i][0] && y < intervals[i][1]) return true;
    }
    return false;
  }

  // 地図端に乗る辺の区間を側ごとに集める。filter(edge)がtrueの辺だけ対象。
  function collectSeamEdgeIntervals(edges, points, frame, tol, filter) {
    const W = [];
    const E = [];
    edges.forEach((edge) => {
      if (filter && !filter(edge)) return;
      const a = points[edge.a];
      const b = points[edge.b];
      const side = seamEdgeSide(a, b, frame, tol);
      if (!side) return;
      (side === "W" ? W : E).push([Math.min(a.y, b.y), Math.max(a.y, b.y)]);
    });
    return { W: mergeIntervals(W, tol), E: mergeIntervals(E, tol) };
  }

  // 地図の西端・東端に「枠」の仮想辺（sources: frame）を追加する。
  // 端に乗る点をyで並べ、隣り合う点の間がまだどの辺にも覆われていなければ辺を張る。
  // これにより、地図端で切れた海上国境（海の領域）も面として閉じる。冪等。追加した辺の数を返す。
  function addSeamFrameEdges(graph, frame, tol) {
    const { points, edges } = graph;
    const existing = new Set(edges.map((e) => (e.a < e.b ? `${e.a}-${e.b}` : `${e.b}-${e.a}`)));
    const covered = collectSeamEdgeIntervals(edges, points, frame, tol, null);
    const bySide = { W: [], E: [] };
    points.forEach((p, idx) => {
      const side = seamSideOf(p, frame, tol);
      if (side) bySide[side].push(idx);
    });
    let added = 0;
    ["W", "E"].forEach((side) => {
      const list = bySide[side].sort((i, j) => points[i].y - points[j].y);
      for (let k = 0; k + 1 < list.length; k += 1) {
        const a = list[k];
        const b = list[k + 1];
        const ya = points[a].y;
        const yb = points[b].y;
        if (yb - ya <= EPS) continue;
        if (intervalsContain(covered[side], (ya + yb) / 2)) continue;
        const key = a < b ? `${a}-${b}` : `${b}-${a}`;
        if (existing.has(key)) continue;
        existing.add(key);
        edges.push({ a, b, sources: new Set(["frame"]) });
        added += 1;
      }
    });
    return added;
  }

  // 地図端の区間を分類する。
  // - border: 国境線（sourcesにborder）が端に沿って引かれている区間（西・東どちらか一方にあれば両側に適用）→ 本物の国境
  // - realCoast: 片側だけ海岸線で、反対側の同じyは海 → 経度180度上の本物の海岸線
  // - cut: それ以外（両側とも陸の切断線、または両側とも海）→ 人工的な切断。左右のポリゴンを同一地物として扱う
  // options.mismatchTol: 両側とも陸の区間に隣接する、この長さ以下の「片側だけ陸」の区間は
  //   描画上の位置ずれ（実データでは左右のピースが約3pxずれている）とみなして切断線に含める。既定は枠の高さの0.5%。
  function analyzeSeam(graph, frame, tol, options) {
    const { points, edges } = graph;
    const mismatchTol = options && Number.isFinite(options.mismatchTol) ? options.mismatchTol : frame.height * 0.005;
    const border = collectSeamEdgeIntervals(edges, points, frame, tol, (e) => e.sources.has("border"));
    const coast = collectSeamEdgeIntervals(edges, points, frame, tol,
      (e) => !e.sources.has("border") && e.sources.has("coast"));
    const borderAny = mergeIntervals(border.W.concat(border.E), tol);
    const landBoth = intersectIntervals(coast.W, coast.E);
    const touchesLandBoth = (iv) => landBoth.some((lb) => Math.abs(lb[0] - iv[1]) <= tol || Math.abs(lb[1] - iv[0]) <= tol);
    const dropMismatch = (list) => list.filter((iv) => !(iv[1] - iv[0] <= mismatchTol && touchesLandBoth(iv)));
    const realCoast = {
      W: dropMismatch(subtractIntervals(subtractIntervals(coast.W, coast.E), borderAny)),
      E: dropMismatch(subtractIntervals(subtractIntervals(coast.E, coast.W), borderAny)),
    };
    const blocked = mergeIntervals(borderAny.concat(realCoast.W, realCoast.E), 0);
    const cut = subtractIntervals([[frame.y, frame.y + frame.height]], blocked);
    return { border, coast, borderAny, landBoth, realCoast, cut };
  }

  // 地図端の辺 a→b のうち、人工的な切断区間（seam.cut）に乗らない部分を [[p, q], ...] で返す（向きは a→b のまま）。
  // 1 本の辺が切断区間と本物の海岸線（realCoast）や国境の区間にまたがることがあるので（polygon/polyline の長い辺は
  // 反対側の端点の位置で分割されない）、辺全体を中点で判定せず、区間の境界で切って人工的な部分だけを落とす。
  // 地図端に乗らない辺はそのまま返す。tol 以下の切れ端は捨てる。
  function seamStrokeParts(a, b, frame, tol, seam) {
    if (!seam || !seamEdgeSide(a, b, frame, tol)) return [[a, b]];
    const y0 = Math.min(a.y, b.y);
    const y1 = Math.max(a.y, b.y);
    if (y1 - y0 <= EPS) return [[a, b]];
    const keep = subtractIntervals([[y0, y1]], seam.cut).filter((iv) => iv[1] - iv[0] > tol);
    const at = (y) => {
      if (y === a.y) return a;
      if (y === b.y) return b;
      const t = (y - a.y) / (b.y - a.y);
      return { x: a.x + (b.x - a.x) * t, y };
    };
    const parts = keep.map(([lo, hi]) => (a.y <= b.y ? [at(lo), at(hi)] : [at(hi), at(lo)]));
    if (a.y > b.y) parts.reverse();
    return parts;
  }

  // 辺(a, b)が全体として地図端上の人工的な切断線か（seamはanalyzeSeamの結果）。
  // 一部でも本物の海岸線・国境の区間を含む辺は false（国タイプ判定で境界として数える）。描画は seamStrokeParts で切り分ける。
  function isArtificialSeamEdge(a, b, frame, tol, seam) {
    if (!seam) return false;
    if (!seamEdgeSide(a, b, frame, tol)) return false;
    return seamStrokeParts(a, b, frame, tol, seam).length === 0;
  }

  // 枠の仮想辺だけで海側を閉じた面（国境線を含まない）か。地図端で陸が複数回接するときにできる
  // 「切断線と海岸線に囲まれた海」は実在の領域ではないので除外に使う。getSources(a, b)は辺のsourcesを返す。
  function isFrameOnlySeaFace(ring, getSources) {
    let hasFrameOnly = false;
    for (let i = 0; i < ring.length; i += 1) {
      const sources = getSources(ring[i], ring[(i + 1) % ring.length]);
      if (!sources) continue;
      if (sources.has("border")) return false;
      if (sources.has("frame") && sources.size === 1) hasFrameOnly = true;
    }
    return hasFrameOnly;
  }

  // 各ポリゴン（rings: 点インデックス列の配列）が地図端の人工的な切断区間（cut）に接する区間を求める。
  function collectSeamRuns(polygons, points, frame, tol, cut) {
    return polygons.map((poly) => {
      const west = [];
      const east = [];
      const rings = poly.rings || [poly.points];
      rings.forEach((ring) => {
        for (let i = 0; i < ring.length; i += 1) {
          const a = points[ring[i]];
          const b = points[ring[(i + 1) % ring.length]];
          const side = seamEdgeSide(a, b, frame, tol);
          if (!side) continue;
          (side === "W" ? west : east).push([Math.min(a.y, b.y), Math.max(a.y, b.y)]);
        }
      });
      let W = mergeIntervals(west, tol);
      let E = mergeIntervals(east, tol);
      if (cut) {
        W = intersectIntervals(W, cut);
        E = intersectIntervals(E, cut);
      }
      return { west: W, east: E };
    });
  }

  // 西端の区間と東端の区間が「同じ切断線」とみなせるポリゴン同士をまとめる（union-find、推移的）。
  // 条件: 重なりが短い方の区間長のminOverlapRatio以上、かつ上端・下端の差がどちらもendTol以内。
  // （海上国境が両端でわずかに違うyに当たる場合に、別の国同士をつながないため）
  // options.cut を渡すと人工的な切断区間（analyzeSeam().cut）に限定する。
  // 返り値: groups = 2個以上のポリゴンからなるグループ（ポリゴン添字の配列）、links = 判定根拠。
  function groupSeamPolygons(polygons, points, frame, options) {
    const opts = options || {};
    const tol = Number.isFinite(opts.tol) ? opts.tol : 0;
    const ratio = Number.isFinite(opts.minOverlapRatio) ? opts.minOverlapRatio : 0.5;
    const endTol = Number.isFinite(opts.endTol) ? opts.endTol : frame.height * 0.005;
    const runs = collectSeamRuns(polygons, points, frame, tol, opts.cut || null);
    const parent = polygons.map((_, i) => i);
    const find = (i) => {
      let r = i;
      while (parent[r] !== r) r = parent[r];
      while (parent[i] !== r) {
        const next = parent[i];
        parent[i] = r;
        i = next;
      }
      return r;
    };
    const links = [];
    const eastRuns = [];
    runs.forEach((r, idx) => r.east.forEach((iv) => eastRuns.push({ idx, iv })));
    runs.forEach((r, wIdx) => {
      r.west.forEach((w) => {
        eastRuns.forEach(({ idx: eIdx, iv: e }) => {
          const overlap = Math.min(w[1], e[1]) - Math.max(w[0], e[0]);
          if (overlap <= 0) return;
          const shorter = Math.min(w[1] - w[0], e[1] - e[0]);
          if (overlap < ratio * shorter) return;
          if (Math.abs(w[0] - e[0]) > endTol || Math.abs(w[1] - e[1]) > endTol) return;
          links.push({ west: wIdx, east: eIdx, westRun: w, eastRun: e, overlap });
          if (wIdx !== eIdx) {
            const ra = find(wIdx);
            const rb = find(eIdx);
            if (ra !== rb) parent[Math.max(ra, rb)] = Math.min(ra, rb);
          }
        });
      });
    });
    const byRoot = new Map();
    polygons.forEach((_, i) => {
      const root = find(i);
      if (!byRoot.has(root)) byRoot.set(root, []);
      byRoot.get(root).push(i);
    });
    const groups = Array.from(byRoot.values()).filter((g) => g.length > 1);
    return { groups, links, runs };
  }

  // ---- 国タイプ判定と、その下ごしらえ（面の内外判定・共有する海・海岸線の内外）----
  // index.html の extractPolygons / buildGraph から呼ぶ。DOM と state には触れない
  // （地図端の切断線かどうかだけは地図枠に依存するので、classifyPolygons に関数で渡す）。

  // 面の内外判定・点の照合に使う許容値の下限（地図単位）。
  const GEOM_EPS = 1e-6;

  function pointOnSegment(point, a, b, tol) {
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const len2 = dx * dx + dy * dy;
    const tol2 = tol * tol;
    if (len2 <= tol2) {
      const dist2 = (point.x - a.x) ** 2 + (point.y - a.y) ** 2;
      return dist2 <= tol2;
    }
    const t = ((point.x - a.x) * dx + (point.y - a.y) * dy) / len2;
    if (t < 0 || t > 1) return false;
    const projX = a.x + dx * t;
    const projY = a.y + dy * t;
    const dist2 = (point.x - projX) ** 2 + (point.y - projY) ** 2;
    return dist2 <= tol2;
  }

  // 境界から tol 以内なら内側とする内外判定（クリック選択・海岸線の内外・共有する海の判定用）。
  // pointInRingIndices（境界の扱いを決めない偶奇判定）とは別物なので寄せないこと。
  function pointInPolygon(point, indices, points, tol) {
    let inside = false;
    for (let i = 0, j = indices.length - 1; i < indices.length; j = i++) {
      const a = points[indices[i]];
      const b = points[indices[j]];
      if (!a || !b) continue;
      if (tol > 0 && pointOnSegment(point, a, b, tol)) {
        return true;
      }
      const intersects = (a.y > point.y) !== (b.y > point.y)
        && point.x < ((b.x - a.x) * (point.y - a.y)) / (b.y - a.y) + a.x;
      if (intersects) inside = !inside;
    }
    return inside;
  }

  // ポリゴンの外周リング列。経度180度の継ぎ目でつながった地物は複数のリング（ピース）を持つ。
  function getPolygonRings(poly) {
    return poly.rings || [poly.points];
  }

  function buildPolygonData(polygons, points) {
    return polygons.map((poly) => {
      const rings = getPolygonRings(poly).map((ring) => {
        let minX = Infinity;
        let minY = Infinity;
        let maxX = -Infinity;
        let maxY = -Infinity;
        for (let i = 0; i < ring.length; i += 1) {
          const p = points[ring[i]];
          if (!p) continue;
          if (p.x < minX) minX = p.x;
          if (p.y < minY) minY = p.y;
          if (p.x > maxX) maxX = p.x;
          if (p.y > maxY) maxY = p.y;
        }
        return { points: ring, minX, minY, maxX, maxY };
      });
      return {
        rings,
        minX: Math.min(...rings.map((r) => r.minX)),
        minY: Math.min(...rings.map((r) => r.minY)),
        maxX: Math.max(...rings.map((r) => r.maxX)),
        maxY: Math.max(...rings.map((r) => r.maxY)),
      };
    });
  }

  function findPolygonIndex(point, polygonData, points, tol) {
    for (let i = 0; i < polygonData.length; i += 1) {
      const poly = polygonData[i];
      if (point.x < poly.minX - tol || point.x > poly.maxX + tol
        || point.y < poly.minY - tol || point.y > poly.maxY + tol) {
        continue;
      }
      for (let r = 0; r < poly.rings.length; r += 1) {
        const ring = poly.rings[r];
        if (point.x < ring.minX - tol || point.x > ring.maxX + tol
          || point.y < ring.minY - tol || point.y > ring.maxY + tol) {
          continue;
        }
        if (pointInPolygon(point, ring.points, points, tol)) {
          return i;
        }
      }
    }
    return -1;
  }

  /**
   * 「点が海岸線の上にあるか」を tolerance 以内で答える関数を作る。
   * 線分を一様グリッドに入れて近傍セルだけを見る（海岸線は 4 万本あるので総当たりは無理）。
   */
  function makeCoastProximityTest(coastSegments, tolerance) {
    if (!coastSegments || coastSegments.length === 0) return () => false;
    const cellSize = Math.max(tolerance, GEOM_EPS) * 2;
    const tol2 = tolerance * tolerance;
    const grid = new Map();
    const key = (cx, cy) => `${cx}:${cy}`;
    coastSegments.forEach((seg) => {
      // 線分がまたぐセルすべてに入れる（端点のセルだけだと長い線分を取りこぼす）
      const x0 = Math.floor(Math.min(seg.a.x, seg.b.x) / cellSize);
      const x1 = Math.floor(Math.max(seg.a.x, seg.b.x) / cellSize);
      const y0 = Math.floor(Math.min(seg.a.y, seg.b.y) / cellSize);
      const y1 = Math.floor(Math.max(seg.a.y, seg.b.y) / cellSize);
      for (let cx = x0; cx <= x1; cx += 1) {
        for (let cy = y0; cy <= y1; cy += 1) {
          const k = key(cx, cy);
          if (!grid.has(k)) grid.set(k, []);
          grid.get(k).push(seg);
        }
      }
    });

    const distSqToSegment = (p, a, b) => {
      const vx = b.x - a.x;
      const vy = b.y - a.y;
      const len2 = vx * vx + vy * vy;
      if (len2 <= 0) return distanceSq(p, a);
      let t = ((p.x - a.x) * vx + (p.y - a.y) * vy) / len2;
      t = t < 0 ? 0 : t > 1 ? 1 : t;
      const dx = p.x - (a.x + t * vx);
      const dy = p.y - (a.y + t * vy);
      return dx * dx + dy * dy;
    };

    return (point) => {
      const cx = Math.floor(point.x / cellSize);
      const cy = Math.floor(point.y / cellSize);
      for (let dx = -1; dx <= 1; dx += 1) {
        for (let dy = -1; dy <= 1; dy += 1) {
          const list = grid.get(key(cx + dx, cy + dy));
          if (!list) continue;
          for (let i = 0; i < list.length; i += 1) {
            if (distSqToSegment(point, list[i].a, list[i].b) <= tol2) return true;
          }
        }
      }
      return false;
    };
  }

  function collectSharedSeaPolygons(borderGraph, borderPolygons, snap, coastSegments) {
    if (!borderGraph || borderPolygons.length === 0) {
      return { polygonData: [], sharedPolygons: new Set(), tolerance: 0 };
    }
    const degrees = new Array(borderGraph.points.length).fill(0);
    borderGraph.edges.forEach((edge) => {
      degrees[edge.a] += 1;
      degrees[edge.b] += 1;
    });
    const polygonData = buildPolygonData(borderPolygons, borderGraph.points);
    const tolerance = Math.max(GEOM_EPS, snap * 0.1);
    const sharedPolygons = new Set();
    // 要件定義書 §6.3 の「国境線が海岸線で終了し」を字義どおりに取る。
    // 自由端が海岸線の上にあることまで確かめないと、作図で残った長さ 1 前後の「髭」
    // （2 本の国境パスが V 字に出会う尖端の残り。実データに多数ある）を
    // 「海岸線で終わる国境」と取り違え、その国が丸ごと共有する海に化ける。
    // 実測: この条件を足すと共有する海の面が 64 → 28 に減り、報告された誤りが全部消える。
    const isOnCoast = makeCoastProximityTest(coastSegments, Math.max(snap, GEOM_EPS));

    borderGraph.edges.forEach((edge) => {
      const freeEnds = [];
      if (degrees[edge.a] === 1) freeEnds.push(borderGraph.points[edge.a]);
      if (degrees[edge.b] === 1) freeEnds.push(borderGraph.points[edge.b]);
      if (freeEnds.length === 0) return;
      if (!freeEnds.some((point) => isOnCoast(point))) return;
      const a = borderGraph.points[edge.a];
      const b = borderGraph.points[edge.b];
      const mid = { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 };
      const idx = findPolygonIndex(mid, polygonData, borderGraph.points, tolerance);
      if (idx >= 0) {
        sharedPolygons.add(idx);
      }
    });

    return { polygonData, sharedPolygons, tolerance };
  }

  function filterCoastSegments(segments, polygonData, points, sharedPolygons, tolerance) {
    if (!polygonData.length) return segments;
    const filtered = [];
    for (let i = 0; i < segments.length; i += 1) {
      const seg = segments[i];
      const mid = { x: (seg.a.x + seg.b.x) / 2, y: (seg.a.y + seg.b.y) / 2 };
      const idx = findPolygonIndex(mid, polygonData, points, tolerance);
      if (idx === -1 || sharedPolygons.has(idx)) {
        filtered.push(seg);
      }
    }
    return filtered;
  }

  function edgeKey(a, b) {
    return a < b ? `${a}-${b}` : `${b}-${a}`;
  }

  function buildPointLookup(points, snap) {
    const cellSize = snap > 0 ? snap : GEOM_EPS;
    const grid = new Map();

    function gridKey(x, y) {
      return `${x}:${y}`;
    }

    points.forEach((pt, idx) => {
      const cx = Math.floor(pt.x / cellSize);
      const cy = Math.floor(pt.y / cellSize);
      const key = gridKey(cx, cy);
      if (!grid.has(key)) grid.set(key, []);
      grid.get(key).push(idx);
    });

    function findNearestIndex(pt) {
      const maxDist = snap > 0 ? snap : GEOM_EPS;
      const maxDist2 = maxDist * maxDist;
      const cx = Math.floor(pt.x / cellSize);
      const cy = Math.floor(pt.y / cellSize);
      let bestIdx = null;
      let bestDist2 = Infinity;

      for (let dx = -1; dx <= 1; dx += 1) {
        for (let dy = -1; dy <= 1; dy += 1) {
          const key = gridKey(cx + dx, cy + dy);
          const candidates = grid.get(key);
          if (!candidates) continue;
          for (let i = 0; i < candidates.length; i += 1) {
            const idx = candidates[i];
            const d2 = distanceSq(pt, points[idx]);
            if (d2 <= maxDist2 && d2 < bestDist2) {
              bestDist2 = d2;
              bestIdx = idx;
            }
          }
        }
      }

      return bestIdx;
    }

    return { findNearestIndex };
  }

  function collectBorderCoastPoints(segments, snap) {
    const points = [];
    const flags = [];
    if (!segments || segments.length === 0) return points;

    const cellSize = snap > 0 ? snap : GEOM_EPS;
    const maxDist = snap > 0 ? snap : GEOM_EPS;
    const maxDist2 = maxDist * maxDist;
    const grid = new Map();

    function gridKey(x, y) {
      return `${x}:${y}`;
    }

    function addPoint(pt, flag) {
      const cx = Math.floor(pt.x / cellSize);
      const cy = Math.floor(pt.y / cellSize);
      let bestIdx = null;
      let bestDist2 = Infinity;

      for (let dx = -1; dx <= 1; dx += 1) {
        for (let dy = -1; dy <= 1; dy += 1) {
          const key = gridKey(cx + dx, cy + dy);
          const candidates = grid.get(key);
          if (!candidates) continue;
          for (let i = 0; i < candidates.length; i += 1) {
            const idx = candidates[i];
            const d2 = distanceSq(pt, points[idx]);
            if (d2 <= maxDist2 && d2 < bestDist2) {
              bestDist2 = d2;
              bestIdx = idx;
            }
          }
        }
      }

      if (bestIdx === null) {
        const idx = points.length;
        points.push({ x: pt.x, y: pt.y });
        flags.push(flag);
        const key = gridKey(cx, cy);
        if (!grid.has(key)) grid.set(key, []);
        grid.get(key).push(idx);
        return;
      }

      flags[bestIdx] |= flag;
    }

    for (let i = 0; i < segments.length; i += 1) {
      const seg = segments[i];
      const flag = seg.source === "border" ? 1 : seg.source === "coast" ? 2 : 0;
      if (!flag) continue;
      addPoint(seg.a, flag);
      addPoint(seg.b, flag);
    }

    const borderCoastPoints = [];
    for (let i = 0; i < points.length; i += 1) {
      if ((flags[i] & 1) && (flags[i] & 2)) {
        borderCoastPoints.push(points[i]);
      }
    }

    return borderCoastPoints;
  }

  function mapPointsToGraphIndices(points, graphPoints, snap) {
    const indices = new Set();
    if (!points || points.length === 0) return indices;
    const lookup = buildPointLookup(graphPoints, snap);
    points.forEach((pt) => {
      const idx = lookup.findNearestIndex(pt);
      if (idx !== null && idx !== undefined) {
        indices.add(idx);
      }
    });
    return indices;
  }

  function collectCoastInsidePolygons(coastSegments, polygonData, points, tolerance) {
    const inside = new Set();
    if (!coastSegments || coastSegments.length === 0 || polygonData.length === 0) return inside;

    for (let i = 0; i < coastSegments.length; i += 1) {
      if (inside.size === polygonData.length) break;
      const seg = coastSegments[i];
      const mid = { x: (seg.a.x + seg.b.x) / 2, y: (seg.a.y + seg.b.y) / 2 };
      const idx = findPolygonIndex(mid, polygonData, points, tolerance);
      if (idx >= 0) inside.add(idx);
    }

    return inside;
  }

  function buildPolygonEdgeIndex(polygons, graph) {
    if (!graph) return null;
    const edgeSources = new Map();
    graph.edges.forEach((edge) => {
      edgeSources.set(edgeKey(edge.a, edge.b), edge.sources);
    });

    const edgeToPolys = new Map();
    polygons.forEach((poly) => {
      const polyId = String(poly.id);
      getPolygonRings(poly).forEach((pts) => {
        for (let i = 0; i < pts.length; i += 1) {
          const a = pts[i];
          const b = pts[(i + 1) % pts.length];
          const key = edgeKey(a, b);
          if (!edgeToPolys.has(key)) {
            edgeToPolys.set(key, { a, b, polys: new Set() });
          }
          edgeToPolys.get(key).polys.add(polyId);
        }
      });
    });

    return { edgeSources, edgeToPolys };
  }

  /**
   * 国タイプ判定（要件定義書 §6.3）。戻り値は 地物 id（文字列）→ 種別。
   *
   *   polygons      … 地物の配列（{ id, points } または継ぎ目でつないだ { id, rings }）
   *   graph         … 点の座標を引く（graph.points）
   *   edgeSources   … edgeKey(a, b) → その辺の sources（"border" / "coast" / "frame" の Set）
   *   sharedSeaInfo … collectSharedSeaPolygons の結果に、国境線グラフの points を足したもの（無ければ null）
   *   coastInfo     … { borderCoastPoints: 国境線と海岸線が交わる点の添字の Set,
   *                     coastInsidePolygons: 領内に海岸線がある地物の添字の Set }（海岸線が無ければ null）
   *   isCutEdge     … (a, b, sources) → 境界として数えない辺か（地図端の人工的な切断線・枠の仮想辺）。
   *                   地図枠と継ぎ目の解析結果に依存するので呼び出し側（index.html）が渡す。
   *
   * 国境線と海岸線の両方を持つ: 共有する海に接していれば「共有する海」、そうでなければ「海岸線まで国境界」。
   * 国境線だけ: 国境線と海岸線の交点に触れていれば「海まで国境線」、領内に海岸線があれば「完全島国」、
   * どちらでもなければ「内陸国」（AGENT.md「内陸国の誤分類」「内海・完全島国の誤分類」）。
   * 海岸線だけ: 「島国」。
   */
  function classifyPolygons(polygons, graph, edgeSources, sharedSeaInfo, coastInfo, isCutEdge) {
    const types = new Map();
    const borderCoastPoints = coastInfo && coastInfo.borderCoastPoints;
    const coastInsidePolygons = coastInfo && coastInfo.coastInsidePolygons;
    const coastDataAvailable = !!(borderCoastPoints && coastInsidePolygons);
    polygons.forEach((poly, polyIndex) => {
      const polyId = String(poly.id);
      let hasBorder = false;
      let hasCoast = false;
      let touchesSharedSea = false;
      let touchesCoast = false;

      // 継ぎ目をまたぐ地物は全ピースの辺をまとめて1つの地物として判定する。
      const edgesOfPoly = [];
      getPolygonRings(poly).forEach((ring) => {
        for (let i = 0; i < ring.length; i += 1) edgesOfPoly.push([ring[i], ring[(i + 1) % ring.length]]);
      });
      for (let i = 0; i < edgesOfPoly.length; i += 1) {
        const [aIdx, bIdx] = edgesOfPoly[i];
        const edge = edgeKey(aIdx, bIdx);
        if (!touchesCoast && borderCoastPoints && (borderCoastPoints.has(aIdx) || borderCoastPoints.has(bIdx))) {
          touchesCoast = true;
        }
        const sources = edgeSources.get(edge);
        if (!sources) continue;
        if (isCutEdge && isCutEdge(aIdx, bIdx, sources)) continue;
        if (sources.has("border")) hasBorder = true;
        if (sources.has("coast")) {
          hasCoast = true;
          if (!touchesSharedSea && sharedSeaInfo && sharedSeaInfo.polygonData.length) {
            const a = graph.points[aIdx];
            const b = graph.points[bIdx];
            const mid = { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 };
            const idx = findPolygonIndex(mid, sharedSeaInfo.polygonData, sharedSeaInfo.points, sharedSeaInfo.tolerance);
            if (idx >= 0 && sharedSeaInfo.sharedPolygons.has(idx)) {
              touchesSharedSea = true;
            }
          }
        }
      }

      let type = "不明";
      if (hasBorder && hasCoast) {
        type = touchesSharedSea ? "共有する海" : "海岸線まで国境界";
      } else if (hasBorder) {
        if (!coastDataAvailable) {
          type = "海まで国境線";
        } else if (touchesCoast) {
          type = "海まで国境線";
        } else if (coastInsidePolygons.has(polyIndex)) {
          type = "完全島国";
        } else {
          type = "内陸国";
        }
      } else if (hasCoast) {
        type = "島国";
      }
      types.set(polyId, type);
    });

    return types;
  }
  // ---- 球面投影（Wikipedia風出力の正射図法）----
  // index.html の buildWikipediaCache / buildWikipediaSvgMarkup から呼ぶ。DOM と state には触れない。
  // 座標の約束:
  //   - 球面上の点は単位ベクトル {x, y, z}（lonLatDegToVector）。経緯度はラジアン（vectorToLonLat）。
  //   - buildRotation(視点) で回すと視点が回転後 z = 1 に来る。z >= 0 が見えている半球。
  //   - 投影座標は projectRotated（x = 回転後 y, y = 回転後 x）。単位円盤 [-1, 1]² に収まり、北が上・東が右。
  // 過去に 2 度、塗りが反転している（AGENT.md「球面投影の半球クリッピングで地平線弧が欠落」
  // 「地平線弧の向き誤りによる塗りの反転」）。触るときは test/sphere.test.js を通すこと。

  function degToRad(deg) {
    return (deg * Math.PI) / 180;
  }

  function normalizeVector(vec) {
    const len = Math.hypot(vec.x, vec.y, vec.z);
    if (!len) return { x: 1, y: 0, z: 0 };
    return { x: vec.x / len, y: vec.y / len, z: vec.z / len };
  }

  function vectorToLonLat(vec) {
    const norm = normalizeVector(vec);
    return { lon: Math.atan2(norm.y, norm.x), lat: Math.asin(norm.z) };
  }

  function slerpVectors(a, b, t) {
    const dot = Math.max(-1, Math.min(1, a.x * b.x + a.y * b.y + a.z * b.z));
    const omega = Math.acos(dot);
    if (omega < 1e-6) {
      return normalizeVector({
        x: a.x + (b.x - a.x) * t,
        y: a.y + (b.y - a.y) * t,
        z: a.z + (b.z - a.z) * t,
      });
    }
    const sinOmega = Math.sin(omega);
    const k0 = Math.sin((1 - t) * omega) / sinOmega;
    const k1 = Math.sin(t * omega) / sinOmega;
    return normalizeVector({
      x: a.x * k0 + b.x * k1,
      y: a.y * k0 + b.y * k1,
      z: a.z * k0 + b.z * k1,
    });
  }

  // 経緯度（度）→ 単位ベクトル（x: 経度0度・赤道、y: 東経90度・赤道、z: 北極）。
  function lonLatDegToVector(lon, lat) {
    const lonRad = degToRad(lon);
    const latRad = degToRad(lat);
    const cosLat = Math.cos(latRad);
    return {
      x: cosLat * Math.cos(lonRad),
      y: cosLat * Math.sin(lonRad),
      z: Math.sin(latRad),
    };
  }

  // 地図座標の点 → 単位ベクトル。frame は地図枠（正距円筒図法で経度-180〜+180・緯度+90〜-90）、
  // tol は地図端に乗っているとみなす距離（seamSideOf）。端の点は経度±180度ちょうどにする。
  function mapPointToVector(point, frame, tol) {
    if (!point || !frame) return null;
    let lon = ((point.x - frame.x) / frame.width) * 360 - 180;
    const side = seamSideOf(point, frame, tol);
    if (side === "W") lon = -180;
    else if (side === "E") lon = 180;
    const lat = Math.max(-90, Math.min(90, 90 - ((point.y - frame.y) / frame.height) * 180));
    return lonLatDegToVector(lon, lat);
  }

  // 球面へ写す辺の長さの上限（経度・緯度の度）。地図上の直線は経緯度の直線なので、長い辺を
  // 弦のまま球面へ写すと地図とずれる。地平線弧の回転角と真裏の点の内外判定（経緯度平面）が
  // 食い違うと塗りが反転するため、球面へ写す前にこの長さ以下へ刻む。
  const SPHERE_EDGE_MAX_DEG = 1;

  // 地図座標の辺a→bを経緯度で SPHERE_EDGE_MAX_DEG 以下に刻む中間点（両端は含まない）。
  function densifyMapSegment(a, b, vb) {
    const spanDeg = Math.max(
      (Math.abs(b.x - a.x) / vb.width) * 360,
      (Math.abs(b.y - a.y) / vb.height) * 180,
    );
    const pieces = Math.ceil(spanDeg / SPHERE_EDGE_MAX_DEG);
    const out = [];
    for (let k = 1; k < pieces; k += 1) {
      const t = k / pieces;
      out.push({ x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t });
    }
    return out;
  }

  // 閉じた点列の各辺（最後→最初を含む）を densifyMapSegment で刻む。
  function densifyMapRing(points, vb) {
    const out = [];
    for (let i = 0; i < points.length; i += 1) {
      const a = points[i];
      out.push(a);
      densifyMapSegment(a, points[(i + 1) % points.length], vb).forEach((p) => out.push(p));
    }
    return out;
  }

  function buildRotation(centerVec) {
    const { lon, lat } = vectorToLonLat(centerVec);
    return {
      cosLon: Math.cos(lon),
      sinLon: Math.sin(lon),
      cosLat: Math.cos(lat),
      sinLat: Math.sin(lat),
    };
  }

  function rotateVector(vec, rotation) {
    const x1 = rotation.cosLon * vec.x + rotation.sinLon * vec.y;
    const y1 = -rotation.sinLon * vec.x + rotation.cosLon * vec.y;
    const z1 = vec.z;
    const x2 = rotation.sinLat * x1 - rotation.cosLat * z1;
    const y2 = y1;
    const z2 = rotation.cosLat * x1 + rotation.sinLat * z1;
    return { x: x2, y: y2, z: z2 };
  }

  // 投影座標（正射図法）。回転後の (x, y) を入れ替える: 投影 x = 回転後 y（東が右）、
  // 投影 y = 回転後 x（SVG の y は下向きなので北が上）。
  function projectRotated(vec) {
    return { x: vec.y, y: vec.x };
  }

  // 単位ベクトル列を視点の回転で回し、clipPolygonToHemisphere に渡す配列一式を作る。
  // rotX/rotY/rotZ は rotateVector、projX/projY は projectRotated と同じ値。
  // lonArr/latArr は回す前の経緯度（ラジアン。真裏の点の内外判定に使う）。
  function rotateVectorArrays(vectors, rotation) {
    const count = vectors.length;
    const rotX = new Float64Array(count);
    const rotY = new Float64Array(count);
    const rotZ = new Float64Array(count);
    const projX = new Float64Array(count);
    const projY = new Float64Array(count);
    const lonArr = new Float64Array(count);
    const latArr = new Float64Array(count);
    for (let i = 0; i < count; i += 1) {
      const v = vectors[i];
      lonArr[i] = Math.atan2(v.y, v.x);
      latArr[i] = Math.atan2(v.z, Math.hypot(v.x, v.y));
      // rotateVector と同じ式（点の数が多いのでオブジェクトを作らずに書く）。
      const x1 = rotation.cosLon * v.x + rotation.sinLon * v.y;
      const y1 = -rotation.sinLon * v.x + rotation.cosLon * v.y;
      const z1 = v.z;
      const x2 = rotation.sinLat * x1 - rotation.cosLat * z1;
      const y2 = y1;
      const z2 = rotation.cosLat * x1 + rotation.sinLat * z1;
      rotX[i] = x2;
      rotY[i] = y2;
      rotZ[i] = z2;
      projX[i] = y2;
      projY[i] = x2;
    }
    return { rotX, rotY, rotZ, projX, projY, lonArr, latArr };
  }

  function intersectRotatedIndices(aIdx, bIdx, rotX, rotY, rotZ) {
    const zA = rotZ[aIdx];
    const zB = rotZ[bIdx];
    const denom = zA - zB;
    if (Math.abs(denom) < 1e-9) return null;
    const t = zA / denom;
    const x = rotX[aIdx] + t * (rotX[bIdx] - rotX[aIdx]);
    const y = rotY[aIdx] + t * (rotY[bIdx] - rotY[aIdx]);
    const len = Math.hypot(x, y);
    if (!len) return null;
    const nx = x / len;
    const ny = y / len;
    return { x: ny, y: nx };
  }

  function intersectRotatedVectors(a, b) {
    const denom = a.z - b.z;
    if (Math.abs(denom) < 1e-9) return null;
    const t = a.z / denom;
    const x = a.x + t * (b.x - a.x);
    const y = a.y + t * (b.y - a.y);
    const len = Math.hypot(x, y);
    if (!len) return null;
    const nx = x / len;
    const ny = y / len;
    return { x: ny, y: nx };
  }

  function dedupePoints(points) {
    const out = [];
    const tol = 1e-6;
    for (let i = 0; i < points.length; i += 1) {
      const p = points[i];
      if (!out.length) {
        out.push(p);
        continue;
      }
      const last = out[out.length - 1];
      if (Math.abs(p.x - last.x) > tol || Math.abs(p.y - last.y) > tol) {
        out.push(p);
      }
    }
    if (out.length > 2) {
      const first = out[0];
      const last = out[out.length - 1];
      if (Math.abs(first.x - last.x) < tol && Math.abs(first.y - last.y) < tol) {
        out.pop();
      }
    }
    return out;
  }

  // 地平線上をexitAngleからsweep（符号付き、πを超えてもよい）だけ回る弧の補間点を追加する（両端は含まない）。
  function addHorizonArc(output, exitAngle, sweep) {
    if (Math.abs(sweep) < 1e-6) return;
    const steps = Math.max(2, Math.ceil(Math.abs(sweep) / (Math.PI / 36)));
    for (let s = 1; s < steps; s += 1) {
      const angle = exitAngle + sweep * (s / steps);
      output.push({ x: Math.cos(angle), y: Math.sin(angle) });
    }
  }

  // 不可視区間（出口→入口）を置き換える地平線弧を「短い方の弧の角度差diff + 2π×turns」の形で求める。
  // 回転角は不可視区間の点列が視線軸まわりに回った角度に等しく、出口→入口の角度差に
  // 基準半直線（回転後座標で y=0, x<0）を横切った符号付き回数crossTurns×2πを足して求める。
  function buildHorizonArc(exitPt, entryPt, crossTurns) {
    const exitAngle = Math.atan2(exitPt.y, exitPt.x);
    let diff = Math.atan2(entryPt.y, entryPt.x) - exitAngle;
    if (diff > Math.PI) diff -= 2 * Math.PI;
    if (diff < -Math.PI) diff += 2 * Math.PI;
    // 投影座標(x, y)=(回転後y, 回転後x)。基準半直線が分岐線になる角度 atan2(回転後y, 回転後x) で差を取り、投影座標の向きへ符号を戻す。
    const sweep = -(Math.atan2(entryPt.x, entryPt.y) - Math.atan2(exitPt.x, exitPt.y) + 2 * Math.PI * crossTurns);
    return { exitAngle, diff, turns: Math.round((sweep - diff) / (2 * Math.PI)) };
  }

  function horizonArcSweep(arc) {
    return arc.diff + 2 * Math.PI * arc.turns;
  }

  // 視点中心の真裏の点の経緯度（ラジアン）。
  function getBackLonLat(rotation) {
    return {
      lon: Math.atan2(-rotation.sinLon, -rotation.cosLon),
      lat: -Math.atan2(rotation.sinLat, rotation.cosLat),
    };
  }

  // 経緯度平面での点の内外判定（点から北へ伸ばした経線との交差回数）。
  // 正距円筒図法の2D多角形と同じ内外になる（極や地図の上下端に接する多角形でも判定できる）。経度差は短い方で扱う。
  function containsLonLat(indices, lonArr, latArr, lon, lat) {
    let inside = false;
    let prev = indices[indices.length - 1];
    for (let i = 0; i < indices.length; i += 1) {
      const curr = indices[i];
      const a = normalizeAngle(lonArr[prev] - lon);
      const b = a + normalizeAngle(lonArr[curr] - lonArr[prev]);
      if ((a > 0) !== (b > 0)) {
        const crossLat = latArr[prev] + (latArr[curr] - latArr[prev]) * (a / (a - b));
        if (crossLat > lat) inside = !inside;
      }
      prev = curr;
    }
    return inside;
  }

  // 真裏の点Bの内外判定。2D地図から来た多角形は経緯度平面で閉じているので containsLonLat で判定する。
  // 経度の巻き数が±1周の多角形（3Dで極を囲む輪。経緯度平面では閉じない）は、輪が分ける2領域のうち
  // 面積の小さい方を内側とみなし、その側の極へ向かう経線とBからの交差回数で判定する。
  function containsBackPoint(indices, lonArr, latArr, lon, lat) {
    let winding = 0;
    let capNorth = 0;
    let prev = indices[indices.length - 1];
    for (let i = 0; i < indices.length; i += 1) {
      const curr = indices[i];
      const dLon = normalizeAngle(lonArr[curr] - lonArr[prev]);
      winding += dLon;
      // 北極と辺の間の球面上の面積（符号付き）。
      capNorth += dLon * (1 - Math.sin((latArr[prev] + latArr[curr]) / 2));
      prev = curr;
    }
    if (Math.abs(winding) < Math.PI) {
      return containsLonLat(indices, lonArr, latArr, lon, lat);
    }
    const northArea = Math.abs(capNorth);
    const towardNorth = northArea <= 4 * Math.PI - northArea;
    let crossings = 0;
    prev = indices[indices.length - 1];
    for (let i = 0; i < indices.length; i += 1) {
      const curr = indices[i];
      const a = normalizeAngle(lonArr[prev] - lon);
      const b = a + normalizeAngle(lonArr[curr] - lonArr[prev]);
      if ((a > 0) !== (b > 0)) {
        const crossLat = latArr[prev] + (latArr[curr] - latArr[prev]) * (a / (a - b));
        if (towardNorth ? crossLat > lat : crossLat < lat) crossings += 1;
      }
      prev = curr;
    }
    // 極側の領域の点は、極へ向かう経線が輪を偶数回横切る。
    return crossings % 2 === 0;
  }

  function normalizeAngle(angle) {
    if (angle > Math.PI) return angle - 2 * Math.PI;
    if (angle <= -Math.PI) return angle + 2 * Math.PI;
    return angle;
  }

  // 回転後座標で多角形を可視半球(z>=0)にクリップし、投影座標の点列を返す。
  // 不可視区間は、その区間の点列が視線軸まわりに回った角度だけ地平線に沿って回る弧で置き換える
  // （不可視区間が地平線の180°以上に回り込む場合も正しい向きになる）。
  // さらに視点の真裏の点Bが多角形の内側にある場合は、弧の1本を逆回りにして内外を反転させる。
  // lonArr/latArrは各点の経緯度（ラジアン）、backはBの経緯度。
  function clipPolygonToHemisphere(indices, rotX, rotY, rotZ, projX, projY, lonArr, latArr, back) {
    if (!indices.length) return [];
    // 点と弧（{ arc }）を順に並べ、最後に弧を補間点へ展開する。
    const items = [];
    const arcs = [];
    let lastExit = null;
    let firstEntry = null;
    let crossTurns = 0;
    let leadCrossTurns = 0;
    const pushArc = (exitPt, entryPt, turns) => {
      const arc = buildHorizonArc(exitPt, entryPt, turns);
      arcs.push(arc);
      items.push({ arc });
    };
    let prev = indices[indices.length - 1];
    let prevInside = rotZ[prev] >= 0;
    for (let i = 0; i < indices.length; i += 1) {
      const curr = indices[i];
      const currInside = rotZ[curr] >= 0;
      // 辺が不可視側で基準半直線（y=0, x<0）を横切った向き。
      let crossTurn = 0;
      if (!(prevInside && currInside) && (rotY[prev] > 0) !== (rotY[curr] > 0)) {
        const t = rotY[prev] / (rotY[prev] - rotY[curr]);
        const cx = rotX[prev] + t * (rotX[curr] - rotX[prev]);
        const cz = rotZ[prev] + t * (rotZ[curr] - rotZ[prev]);
        if (cx < 0 && cz < 0) crossTurn = rotY[prev] > 0 ? 1 : -1;
      }
      if (prevInside && currInside) {
        items.push({ x: projX[curr], y: projY[curr] });
      } else if (prevInside && !currInside) {
        const hit = intersectRotatedIndices(prev, curr, rotX, rotY, rotZ);
        if (hit) {
          items.push(hit);
          lastExit = hit;
        }
        crossTurns = crossTurn;
      } else if (!prevInside && !currInside) {
        crossTurns += crossTurn;
      } else {
        crossTurns += crossTurn;
        const hit = intersectRotatedIndices(prev, curr, rotX, rotY, rotZ);
        if (hit) {
          if (lastExit) {
            pushArc(lastExit, hit, crossTurns);
            lastExit = null;
          } else if (!firstEntry) {
            firstEntry = hit;
            leadCrossTurns = crossTurns;
          }
          items.push(hit);
        }
        items.push({ x: projX[curr], y: projY[curr] });
        crossTurns = 0;
      }
      prev = curr;
      prevInside = currInside;
    }
    if (lastExit && firstEntry) {
      pushArc(lastExit, firstEntry, crossTurns + leadCrossTurns);
    }
    if (arcs.length && containsBackPoint(indices, lonArr, latArr, back.lon, back.lat)) {
      // 弧の回転角の総和の符号は、見えている外側領域を囲む向きを表す。
      // その向きの弧のうち最も長いものを逆回り（1周分戻す）にして、内外を入れ替える。
      const sign = arcs.reduce((sum, arc) => sum + horizonArcSweep(arc), 0) >= 0 ? 1 : -1;
      let target = null;
      arcs.forEach((arc) => {
        const sweep = horizonArcSweep(arc) * sign;
        if (sweep > 0 && (!target || sweep > horizonArcSweep(target) * sign)) target = arc;
      });
      if (target) target.turns -= sign;
    }
    const output = [];
    items.forEach((item) => {
      if (item.arc) {
        addHorizonArc(output, item.arc.exitAngle, horizonArcSweep(item.arc));
      } else {
        output.push(item);
      }
    });
    return dedupePoints(output);
  }

  // 単位ベクトルの点列（閉じた多角形）を可視半球にクリップし、投影座標の点列を返す。
  function clipPolygonVectors(vectors, rotation) {
    if (!vectors.length) return [];
    const count = vectors.length;
    const indices = new Array(count);
    for (let i = 0; i < count; i += 1) indices[i] = i;
    const { rotX, rotY, rotZ, projX, projY, lonArr, latArr } = rotateVectorArrays(vectors, rotation);
    return clipPolygonToHemisphere(indices, rotX, rotY, rotZ, projX, projY, lonArr, latArr, getBackLonLat(rotation));
  }

  function formatCoord(value) {
    return value.toFixed(5);
  }

  function buildPathFromPoints(points) {
    if (points.length < 3) return "";
    let d = `M${formatCoord(points[0].x)} ${formatCoord(points[0].y)}`;
    for (let i = 1; i < points.length; i += 1) {
      d += ` L${formatCoord(points[i].x)} ${formatCoord(points[i].y)}`;
    }
    return `${d} Z`;
  }

  function buildSegmentPath(segments, rotation) {
    const parts = [];
    for (let i = 0; i < segments.length; i += 1) {
      const seg = segments[i];
      if (!seg || !seg.a || !seg.b) continue;
      const a = rotateVector(seg.a, rotation);
      const b = rotateVector(seg.b, rotation);
      const zA = a.z;
      const zB = b.z;
      if (zA >= 0 && zB >= 0) {
        const pa = projectRotated(a);
        const pb = projectRotated(b);
        parts.push(`M${formatCoord(pa.x)} ${formatCoord(pa.y)} L${formatCoord(pb.x)} ${formatCoord(pb.y)}`);
        continue;
      }
      if (zA >= 0 || zB >= 0) {
        const hit = intersectRotatedVectors(a, b);
        if (!hit) continue;
        if (zA >= 0) {
          const pa = projectRotated(a);
          parts.push(`M${formatCoord(pa.x)} ${formatCoord(pa.y)} L${formatCoord(hit.x)} ${formatCoord(hit.y)}`);
        } else {
          const pb = projectRotated(b);
          parts.push(`M${formatCoord(hit.x)} ${formatCoord(hit.y)} L${formatCoord(pb.x)} ${formatCoord(pb.y)}`);
        }
      }
    }
    return parts.join(" ");
  }

  function buildGraticulePath(rotation) {
    const parts = [];
    const step = 30;
    const sample = 3;
    for (let lon = -180; lon <= 180; lon += step) {
      const line = [];
      for (let lat = -90; lat <= 90; lat += sample) {
        line.push(rotateVector(lonLatDegToVector(lon, lat), rotation));
      }
      parts.push(buildVisiblePolyline(line));
    }
    for (let lat = -60; lat <= 60; lat += step) {
      const line = [];
      for (let lon = -180; lon <= 180; lon += sample) {
        line.push(rotateVector(lonLatDegToVector(lon, lat), rotation));
      }
      parts.push(buildVisiblePolyline(line));
    }
    return parts.filter(Boolean).join(" ");
  }

  function buildVisiblePolyline(rotatedPoints) {
    if (!rotatedPoints.length) return "";
    let d = "";
    let active = false;
    let prev = rotatedPoints[0];
    for (let i = 1; i < rotatedPoints.length; i += 1) {
      const curr = rotatedPoints[i];
      const prevInside = prev.z >= 0;
      const currInside = curr.z >= 0;
      if (prevInside && currInside) {
        const pp = projectRotated(prev);
        const pc = projectRotated(curr);
        if (!active) {
          d += `M${formatCoord(pp.x)} ${formatCoord(pp.y)}`;
          active = true;
        }
        d += ` L${formatCoord(pc.x)} ${formatCoord(pc.y)}`;
      } else if (prevInside && !currInside) {
        const hit = intersectRotatedVectors(prev, curr);
        if (hit) {
          const pp = projectRotated(prev);
          if (!active) {
            d += `M${formatCoord(pp.x)} ${formatCoord(pp.y)}`;
          }
          d += ` L${formatCoord(hit.x)} ${formatCoord(hit.y)}`;
        }
        active = false;
      } else if (!prevInside && currInside) {
        const hit = intersectRotatedVectors(prev, curr);
        if (hit) {
          const pc = projectRotated(curr);
          d += `M${formatCoord(hit.x)} ${formatCoord(hit.y)} L${formatCoord(pc.x)} ${formatCoord(pc.y)}`;
          active = true;
        }
      } else {
        active = false;
      }
      prev = curr;
    }
    return d;
  }
  // --- 矩形選択（ドラッグで囲んだ面をまとめて選ぶ。要件定義書 §6.7） ---------------
  //
  // 判定は「面が矩形の内側に完全に収まっているか」。矩形は凸なので、リングの全頂点が
  // 矩形内にあることと、リングの外接矩形が矩形内にあることは同値になる。だから面ごとの
  // 外接矩形を 1 度だけ求めておけば、ドラッグ中の判定は外接矩形の比較だけで済む。
  // 「触れたら選ぶ」にしないのは、島を囲むと周りの海の面や、少しかすめた大陸の国まで
  // 拾ってしまうため（細かい島々をまとめて選ぶ、という用途に合わない）。

  /** 2 隅から正規化した矩形 { minX, minY, maxX, maxY } を作る（どちら向きに引いても同じ）。 */
  function rectFromCorners(a, b) {
    return {
      minX: Math.min(a.x, b.x),
      minY: Math.min(a.y, b.y),
      maxX: Math.max(a.x, b.x),
      maxY: Math.max(a.y, b.y),
    };
  }

  /** 外接矩形 bounds が矩形 rect の内側（境界上を含む）に収まるか。 */
  function boundsInsideRect(bounds, rect) {
    return bounds.minX >= rect.minX && bounds.maxX <= rect.maxX
      && bounds.minY >= rect.minY && bounds.maxY <= rect.maxY;
  }

  /** リング（頂点番号の列）が矩形の内側に完全に収まるか。 */
  function ringInsideRect(ringIndices, points, rect) {
    if (!ringIndices || ringIndices.length === 0) return false;
    return boundsInsideRect(ringBounds(ringIndices, points), rect);
  }

  /**
   * 矩形選択の索引を作る。features は [{ id, rings: [頂点番号の列, ...] }]
   * （rings が無ければ points を 1 本のリングとみなす）。戻り値は [{ id, bounds: [各リングの外接矩形] }]。
   * 経度180度の継ぎ目でつないだ地物は複数のリングを持つ。
   */
  function buildRectSelectionIndex(features, points) {
    return features.map((feature) => {
      const rings = feature.rings || [feature.points];
      return {
        id: String(feature.id),
        bounds: rings.filter((ring) => ring && ring.length).map((ring) => ringBounds(ring, points)),
      };
    });
  }

  /**
   * 矩形に収まる地物の id を返す（索引の並び順）。地物のピースのどれか 1 つが収まれば選ぶ。
   * 継ぎ目でつないだ地物はどちらのピースをクリックしても同じ地物が選ばれるので、それに合わせる。
   */
  function collectRectSelectionHits(index, rect) {
    const hits = [];
    for (let i = 0; i < index.length; i += 1) {
      const entry = index[i];
      if (entry.bounds.some((b) => boundsInsideRect(b, rect))) hits.push(entry.id);
    }
    return hits;
  }

  /**
   * 修飾キーから矩形選択のモードを決める。
   * Alt / Ctrl / ⌘ → "remove"（囲んだ面を選択から外す）、Shift → "add"（追加）、なし → "replace"（置き換え）。
   * Alt と Shift を同時に押したときは外す側を優先する（誤って広げるより安全）。
   */
  function rectSelectionMode(modifiers) {
    const m = modifiers || {};
    if (m.altKey || m.ctrlKey || m.metaKey) return "remove";
    if (m.shiftKey) return "add";
    return "replace";
  }

  /**
   * 矩形選択を選択集合に適用した結果を新しい Set で返す（current は変更しない）。
   * replace: 囲んだ面だけにする（何も囲まなければ空。空白クリックで解除するのと同じ）。
   * add: 囲んだ面を足す（Shift+クリックと違ってトグルしない。混ざった範囲を囲んでも全部選ばれる）。
   * remove: 囲んだ面を外す。
   */
  function applyRectSelection(current, hits, mode) {
    if (mode === "replace") return new Set(hits);
    const next = new Set(current);
    if (mode === "add") {
      hits.forEach((id) => next.add(id));
    } else if (mode === "remove") {
      hits.forEach((id) => next.delete(id));
    }
    return next;
  }

  /** ポインタが押した位置から threshold（画面 px）を超えて動いたか。これ以下はクリックとして扱う。 */
  function exceedsDragThreshold(start, current, threshold) {
    return Math.hypot(current.x - start.x, current.y - start.y) > threshold;
  }

  global.DerleLogic = {
    setMeasureSvg,
    nextFrame,
    extractSegmentsFromLayer,
    splitSegmentsAtIntersections,
    buildGraphFromSegments,
    bridgeOpenBorderGaps,
    extractPolygonsFromGraph,
    // 矩形選択（index.html の地図上のドラッグが使う）
    rectFromCorners,
    boundsInsideRect,
    ringInsideRect,
    buildRectSelectionIndex,
    collectRectSelectionHits,
    rectSelectionMode,
    applyRectSelection,
    exceedsDragThreshold,
    // 点列化の方式切替と、新旧の突き合わせに使う下位関数（logic_test.html が使う）
    setSampleMode,
    getSampleMode,
    setCurveTolerance,
    getCurveTolerance,
    samplePathElement,
    samplePathElementByLength,
    samplePathElementFromData,
    pathDataToPoints,
    splitPathData,
    // 領有関係（要件定義書 §6.3）
    computeOwnership,
    computeContainment,
    OWNING_TYPES,
    // 経度180度の継ぎ目（地図の西端・東端）
    seamSideOf,
    seamEdgeSide,
    snapPointsToSeam,
    mergeIntervals,
    intersectIntervals,
    subtractIntervals,
    intervalsContain,
    addSeamFrameEdges,
    analyzeSeam,
    isArtificialSeamEdge,
    seamStrokeParts,
    isFrameOnlySeaFace,
    collectSeamRuns,
    groupSeamPolygons,
    // 要素の点列化と座標変換（index.html もこれを使う。写しを置かないこと）
    applyTransform,
    samplePolylineElement,
    sampleLineElement,
    sampleRectElement,
    // 純粋な幾何の部品（index.html と test/ が使う）
    distance,
    distanceSq,
    removeSequentialDuplicates,
    segmentIntersections,
    pointSegmentProjection,
    polygonArea,
    pointInRingIndices,
    ringBounds,
    pointOnOpenSegment,
    properlyCrosses,
    // 国タイプ判定（index.html と test/ が使う）
    GEOM_EPS,
    pointOnSegment,
    pointInPolygon,
    getPolygonRings,
    buildPolygonData,
    findPolygonIndex,
    makeCoastProximityTest,
    collectSharedSeaPolygons,
    filterCoastSegments,
    edgeKey,
    buildPointLookup,
    collectBorderCoastPoints,
    mapPointsToGraphIndices,
    collectCoastInsidePolygons,
    buildPolygonEdgeIndex,
    classifyPolygons,
    // 球面投影（index.html と test/ が使う。写しを置かないこと）
    SPHERE_EDGE_MAX_DEG,
    degToRad,
    normalizeVector,
    vectorToLonLat,
    slerpVectors,
    lonLatDegToVector,
    mapPointToVector,
    densifyMapSegment,
    densifyMapRing,
    buildRotation,
    rotateVector,
    projectRotated,
    rotateVectorArrays,
    intersectRotatedIndices,
    intersectRotatedVectors,
    dedupePoints,
    addHorizonArc,
    buildHorizonArc,
    horizonArcSweep,
    getBackLonLat,
    normalizeAngle,
    containsLonLat,
    containsBackPoint,
    clipPolygonToHemisphere,
    clipPolygonVectors,
    formatCoord,
    buildPathFromPoints,
    buildSegmentPath,
    buildVisiblePolyline,
    buildGraticulePath,
  };

  // Node（`node --test`）から require できるようにする。ブラウザの <script> では
  // module が無いので何もしない。
  if (typeof module === "object" && module && module.exports) {
    module.exports = global.DerleLogic;
  }
})(typeof window !== "undefined" ? window : globalThis);
