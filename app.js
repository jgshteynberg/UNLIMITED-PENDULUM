const defaults = [
  { length: 1.0, angleDeg: 35, omega: 0.0, mass: 1.0, damping: 0.02 },
  { length: 1.0, angleDeg: -20, omega: 0.2, mass: 0.8, damping: 0.02 },
  { length: 0.9, angleDeg: 10, omega: -0.2, mass: 0.6, damping: 0.02 }
];

const gravityInput = document.getElementById('gravity');
const countInput = document.getElementById('count');
const updateCountBtn = document.getElementById('update-count');
const rowsContainer = document.getElementById('pendulum-rows');
const applyBtn = document.getElementById('apply');
const pauseBtn = document.getElementById('pause');
const resetBtn = document.getElementById('reset');
const canvas = document.getElementById('sim-canvas');
const ctx = canvas.getContext('2d');

let gravity = 9.81;
let links = [];
let paused = false;
let lastTime = performance.now();
let bobTrail = [];

const palette = ['#66d9ff', '#ff9f6e', '#9fff9f', '#dba7ff', '#ffe269', '#ff7fb5', '#8ecbff', '#a6f3c6'];

function makeRow(index, values = {}) {
  const row = document.createElement('tr');
  row.innerHTML = `
    <td>${index + 1}</td>
    <td><input type="number" step="0.05" min="0.1" value="${values.length ?? 1.0}" data-field="length"></td>
    <td><input type="number" step="1" value="${values.angleDeg ?? 0}" data-field="angleDeg"></td>
    <td><input type="number" step="0.05" value="${values.omega ?? 0}" data-field="omega"></td>
    <td><input type="number" step="0.1" min="0.1" value="${values.mass ?? 1}" data-field="mass"></td>
    <td><input type="number" step="0.005" min="0" value="${values.damping ?? 0.02}" data-field="damping"></td>
  `;
  return row;
}

function syncRows(count) {
  rowsContainer.innerHTML = '';
  for (let i = 0; i < count; i += 1) {
    rowsContainer.appendChild(makeRow(i, defaults[i] ?? defaults[defaults.length - 1]));
  }
}

function parseInputs() {
  gravity = Number(gravityInput.value);
  if (!Number.isFinite(gravity) || gravity <= 0) {
    gravity = 9.81;
    gravityInput.value = String(gravity);
  }

  links = [...rowsContainer.querySelectorAll('tr')].map((row) => {
    const fields = Object.fromEntries(
      [...row.querySelectorAll('input')].map((input) => [input.dataset.field, Number(input.value)])
    );
    return {
      length: Math.max(0.1, fields.length || 1),
      theta: (fields.angleDeg || 0) * (Math.PI / 180),
      omega: fields.omega || 0,
      mass: Math.max(0.1, fields.mass || 1),
      damping: Math.max(0, fields.damping || 0)
    };
  });

  bobTrail = [];
}

function resizeCanvasToDisplaySize() {
  const scale = window.devicePixelRatio || 1;
  const width = Math.floor(canvas.clientWidth * scale);
  const height = Math.floor(canvas.clientHeight * scale);
  if (canvas.width !== width || canvas.height !== height) {
    canvas.width = width;
    canvas.height = height;
  }
}

function cumulativeMassFrom(i) {
  let total = 0;
  for (let k = i; k < links.length; k += 1) {
    total += links[k].mass;
  }
  return total;
}

function computeMassMatrix() {
  const n = links.length;
  const M = Array.from({ length: n }, () => Array(n).fill(0));

  for (let i = 0; i < n; i += 1) {
    for (let j = 0; j < n; j += 1) {
      const cm = cumulativeMassFrom(Math.max(i, j));
      const li = links[i].length;
      const lj = links[j].length;
      const angle = links[i].theta - links[j].theta;
      M[i][j] = cm * li * lj * Math.cos(angle);
    }
  }

  return M;
}

function dM(i, j, r) {
  if (i === j) {
    return 0;
  }

  const cm = cumulativeMassFrom(Math.max(i, j));
  const li = links[i].length;
  const lj = links[j].length;
  const s = cm * li * lj * Math.sin(links[i].theta - links[j].theta);

  if (r === i) {
    return -s;
  }
  if (r === j) {
    return s;
  }
  return 0;
}

function computeBiasTerms() {
  const n = links.length;
  const h = Array(n).fill(0);
  const g = Array(n).fill(0);
  const damp = Array(n).fill(0);

  for (let i = 0; i < n; i += 1) {
    for (let j = 0; j < n; j += 1) {
      for (let k = 0; k < n; k += 1) {
        const gamma = 0.5 * (dM(i, j, k) + dM(i, k, j) - dM(j, k, i));
        h[i] += gamma * links[j].omega * links[k].omega;
      }
    }

    g[i] = cumulativeMassFrom(i) * gravity * links[i].length * Math.sin(links[i].theta);
    damp[i] = links[i].damping * links[i].omega;
  }

  return { h, g, damp };
}

function solveLinearSystem(A, b) {
  const n = b.length;
  const M = A.map((row, i) => [...row, b[i]]);

  for (let p = 0; p < n; p += 1) {
    let max = p;
    for (let i = p + 1; i < n; i += 1) {
      if (Math.abs(M[i][p]) > Math.abs(M[max][p])) {
        max = i;
      }
    }
    [M[p], M[max]] = [M[max], M[p]];

    const pivot = M[p][p] || 1e-12;
    for (let j = p; j <= n; j += 1) {
      M[p][j] /= pivot;
    }

    for (let i = 0; i < n; i += 1) {
      if (i === p) {
        continue;
      }
      const factor = M[i][p];
      for (let j = p; j <= n; j += 1) {
        M[i][j] -= factor * M[p][j];
      }
    }
  }

  return M.map((row) => row[n]);
}

function stepPhysics(dt) {
  if (links.length === 0) {
    return;
  }

  const M = computeMassMatrix();
  const { h, g, damp } = computeBiasTerms();
  const rhs = links.map((_, i) => -(h[i] + g[i] + damp[i]));
  const qdd = solveLinearSystem(M, rhs);

  for (let i = 0; i < links.length; i += 1) {
    links[i].omega += qdd[i] * dt;
    links[i].theta += links[i].omega * dt;
  }
}

function computePositions() {
  const positions = [];
  let x = 0;
  let y = 0;

  for (const link of links) {
    x += link.length * Math.sin(link.theta);
    y += link.length * Math.cos(link.theta);
    positions.push({ x, y });
  }

  return positions;
}

function render() {
  resizeCanvasToDisplaySize();
  const w = canvas.width;
  const h = canvas.height;

  ctx.clearRect(0, 0, w, h);
  const pivot = { x: w * 0.5, y: h * 0.12 };

  ctx.strokeStyle = '#426384';
  ctx.lineWidth = Math.max(1.2, w * 0.0016);
  ctx.beginPath();
  ctx.moveTo(0, pivot.y);
  ctx.lineTo(w, pivot.y);
  ctx.stroke();

  if (links.length === 0) {
    return;
  }

  const totalLength = links.reduce((sum, link) => sum + link.length, 0);
  const scale = (h * 0.78) / totalLength;
  const positions = computePositions();

  const end = positions[positions.length - 1];
  bobTrail.push({ x: pivot.x + end.x * scale, y: pivot.y + end.y * scale });
  if (bobTrail.length > 220) {
    bobTrail.shift();
  }

  ctx.strokeStyle = 'rgba(102, 217, 255, 0.85)';
  ctx.beginPath();
  bobTrail.forEach((point, i) => {
    if (i === 0) {
      ctx.moveTo(point.x, point.y);
    } else {
      ctx.lineTo(point.x, point.y);
    }
  });
  ctx.stroke();

  let prevX = pivot.x;
  let prevY = pivot.y;
  positions.forEach((pos, i) => {
    const x = pivot.x + pos.x * scale;
    const y = pivot.y + pos.y * scale;
    const color = palette[i % palette.length];

    ctx.strokeStyle = '#f1f7ff';
    ctx.beginPath();
    ctx.moveTo(prevX, prevY);
    ctx.lineTo(x, y);
    ctx.stroke();

    const bobRadius = Math.max(4.5, Math.sqrt(links[i].mass) * 7 * (w / 1200));
    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.arc(x, y, bobRadius, 0, Math.PI * 2);
    ctx.fill();

    prevX = x;
    prevY = y;
  });
}

function animate(time) {
  const dt = Math.min(0.016, (time - lastTime) / 1000);
  lastTime = time;

  if (!paused) {
    for (let i = 0; i < 2; i += 1) {
      stepPhysics(dt / 2);
    }
  }

  render();
  requestAnimationFrame(animate);
}

updateCountBtn.addEventListener('click', () => {
  const count = Math.min(10, Math.max(1, Number(countInput.value) || 1));
  countInput.value = String(count);
  syncRows(count);
});

applyBtn.addEventListener('click', parseInputs);

pauseBtn.addEventListener('click', () => {
  paused = !paused;
  pauseBtn.textContent = paused ? 'Resume' : 'Pause';
});

resetBtn.addEventListener('click', () => {
  parseInputs();
  lastTime = performance.now();
});

window.addEventListener('resize', render);

syncRows(Number(countInput.value));
parseInputs();
requestAnimationFrame((time) => {
  lastTime = time;
  animate(time);
});
