(function (global) {
  "use strict";

  const EPS = 1e-9;
  let measureSvg = null;

  function setMeasureSvg(svgEl) {
    measureSvg = svgEl;
  }

  function nextFrame() {
    return new Promise((resolve) => requestAnimationFrame(() => resolve()));
  }

  function splitPathData(d) {
    if (!d) return [];
    return d.trim().split(/(?=[Mm])/).map((part) => part.trim()).filter(Boolean);
  }

  /**
   * レイヤー内に「2 つ目以降のサブパスが相対 `m` で始まる `d`」がいくつあるかを数える。
   *
   * サブパスは単独の点列として扱うが、SVG では **パスの最初の `m` は絶対 `M` として
   * 解釈される**。そのため相対 `m` で始まる 2 つ目以降のサブパスを切り出すと、本来
   * 「直前のサブパスの終点からの相対」だった始点が原点からの絶対座標に化け、**エラーも
   * 警告も出ないまま図形だけが移動する**。
   *
   * 現在の原図は全サブパスが絶対 `M` なので実害は無い（実測 0 件）。Illustrator の
   * 書き出し設定が変わった日に無症状で発症する種類の問題なので、数えて表に出す。
   */
  function countRelativeSubpathStarts(layerEl) {
    if (!layerEl) return 0;
    let count = 0;
    const paths = layerEl.querySelectorAll("path");
    for (let i = 0; i < paths.length; i += 1) {
      const parts = splitPathData(paths[i].getAttribute("d") || "");
      for (let k = 1; k < parts.length; k += 1) {
        if (parts[k].charCodeAt(0) === 109) count += 1; // 'm'
      }
    }
    return count;
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
      if (/[A-Za-z]/.test(token)) {
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

  function distance(a, b) {
    const dx = a.x - b.x;
    const dy = a.y - b.y;
    return Math.hypot(dx, dy);
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

  global.DerleLogic = {
    setMeasureSvg,
    nextFrame,
    extractSegmentsFromLayer,
    splitSegmentsAtIntersections,
    buildGraphFromSegments,
    bridgeOpenBorderGaps,
    extractPolygonsFromGraph,
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
    countRelativeSubpathStarts,
    // 領有関係（要件定義書 §6.3）
    computeOwnership,
    computeContainment,
    OWNING_TYPES,
  };
})(window);
