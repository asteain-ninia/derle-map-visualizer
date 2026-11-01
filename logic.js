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

  function samplePathElement(pathEl, step, pathDataOverride) {
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

  // ---- 経度180度の継ぎ目（地図の西端・東端）----
  // frame: 正距円筒図法で経度-180〜+180・緯度+90〜-90に対応する矩形 {x, y, width, height}

  // 点が地図枠の西端/東端のどちらに乗っているか（tol以内）。"W" / "E" / null
  function seamSideOf(point, frame, tol) {
    if (!point || !frame) return null;
    if (Math.abs(point.x - frame.x) <= tol) return "W";
    if (Math.abs(point.x - (frame.x + frame.width)) <= tol) return "E";
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
      if (Math.abs(p.x - west) <= tol) {
        if (p.x !== west) count += 1;
        p.x = west;
      } else if (Math.abs(p.x - east) <= tol) {
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

  // 辺(a, b)が地図端上の人工的な切断線か（seamはanalyzeSeamの結果）。
  function isArtificialSeamEdge(a, b, frame, tol, seam) {
    if (!seam) return false;
    if (!seamEdgeSide(a, b, frame, tol)) return false;
    return intervalsContain(seam.cut, (a.y + b.y) / 2);
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

  global.DerleLogic = {
    setMeasureSvg,
    nextFrame,
    extractSegmentsFromLayer,
    splitSegmentsAtIntersections,
    buildGraphFromSegments,
    extractPolygonsFromGraph,
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
    isFrameOnlySeaFace,
    collectSeamRuns,
    groupSeamPolygons,
  };
})(window);
