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
    return d.trim().split(/(?=M)/).map((part) => part.trim()).filter(Boolean);
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

  async function extractSegmentsFromLayer(layerEl, step, sourceTag, onProgress) {
    const segments = [];
    const elements = Array.from(layerEl.querySelectorAll("path, polyline, polygon, line"));
    const total = elements.length;
    const yieldEvery = 20;

    const pushSegments = (points) => {
      const cleaned = removeSequentialDuplicates(points, 1e-6);
      for (let i = 0; i < cleaned.length - 1; i += 1) {
        segments.push({
          a: cleaned[i],
          b: cleaned[i + 1],
          source: sourceTag,
        });
      }
    };

    for (let index = 0; index < total; index += 1) {
      const item = elements[index];
      const tag = item.tagName.toLowerCase();
      if (tag === "path") {
        const d = item.getAttribute("d") || "";
        const subpaths = splitPathData(d);
        if (subpaths.length === 1) {
          const points = samplePathElement(item, step, subpaths[0]);
          pushSegments(points);
        } else if (subpaths.length > 1) {
          subpaths.forEach((subpath) => {
            const points = samplePathElement(item, step, subpath);
            pushSegments(points);
          });
        }
      } else if (tag === "polyline" || tag === "polygon") {
        const points = samplePolylineElement(item);
        if (tag === "polygon" && points.length > 2) {
          points.push(points[0]);
        }
        pushSegments(points);
      } else if (tag === "line") {
        const points = sampleLineElement(item);
        pushSegments(points);
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
    let snap = Number(endpointSnap) || 0;
    let progress = onProgress;
    if (typeof endpointSnap === "function") {
      progress = endpointSnap;
      snap = 0;
    }
    const splitEps = 1e-6;
    const pad = snap > 0 ? snap : 0;
    const segIntersections = new Map();
    const grid = new Map();

    function gridKey(x, y) {
      return `${x},${y}`;
    }

    function addToGrid(idx, seg) {
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

    function checkEndpointSnap(point, targetSeg, targetIdx) {
      const projection = pointSegmentProjection(point, targetSeg.a, targetSeg.b);
      if (!projection) return;
      if (projection.distance <= snap) {
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
          if (snap > 0) {
            checkEndpointSnap(segA.a, segB, b);
            checkEndpointSnap(segA.b, segB, b);
            checkEndpointSnap(segB.a, segA, a);
            checkEndpointSnap(segB.b, segA, a);
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
          output.push({ a: p0, b: p1, source: seg.source });
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
    const edges = [];
    const edgeMap = new Map();
    const cellSize = snap > 0 ? snap : 1e-6;

    function gridKey(x, y) {
      return `${x}:${y}`;
    }

    function addPoint(pt) {
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
            if (d <= snap && d < bestDist) {
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
        const key = gridKey(cx, cy);
        if (!pointGrid.has(key)) pointGrid.set(key, []);
        pointGrid.get(key).push(idx);
        return idx;
      }

      const hits = pointHits[bestIdx] + 1;
      points[bestIdx].x = (points[bestIdx].x * pointHits[bestIdx] + pt.x) / hits;
      points[bestIdx].y = (points[bestIdx].y * pointHits[bestIdx] + pt.y) / hits;
      pointHits[bestIdx] = hits;
      return bestIdx;
    }

    const buildYield = 2000;
    for (let i = 0; i < segments.length; i += 1) {
      const seg = segments[i];
      const a = addPoint(seg.a);
      const b = addPoint(seg.b);
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

  global.DerleLogic = {
    setMeasureSvg,
    nextFrame,
    extractSegmentsFromLayer,
    splitSegmentsAtIntersections,
    buildGraphFromSegments,
    extractPolygonsFromGraph,
  };
})(window);
