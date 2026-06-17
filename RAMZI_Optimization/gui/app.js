const state = {
  parameters: {},
  bounds: {},
  rings: [],
  frequency: {},
  target: {},
  optimizer: {},
  outputChannel: "cross",
  autoTimer: null,
};

const parameterLabels = {
  fai1: "fai1",
  fai2: "fai2",
  fai3: "fai3",
  fai4: "fai4",
  fait: "fait (faib=-fait)",
  theta13: "theta1 = theta3",
  theta24: "theta2 = theta4",
};

function log(message) {
  const output = document.getElementById("logOutput");
  output.textContent = `${new Date().toLocaleTimeString()}  ${message}\n${output.textContent}`;
}

async function api(path, payload) {
  const response = await fetch(path, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  const data = await response.json();
  if (!response.ok || data.error) throw new Error(data.error || response.statusText);
  return data;
}

function renderRingInfo() {
  const container = document.getElementById("ringInfo");
  container.innerHTML = state.rings.map((ring) => {
    const lengthUm = ring.length_m * 1e6;
    const fsrGhz = ring.fsr_hz / 1e9;
    return `<div><strong>${ring.name}</strong><span>L=${lengthUm.toFixed(0)} um</span><span>ng=${ring.group_index.toFixed(2)}</span><span>FSR=${fsrGhz.toFixed(2)} GHz</span></div>`;
  }).join("");
}

function renderParameters() {
  const panel = document.getElementById("parameterPanel");
  panel.innerHTML = "";
  Object.entries(state.parameters).forEach(([key, value]) => {
    const bounds = state.bounds[key] || [-Math.PI, Math.PI];
    const row = document.createElement("div");
    row.className = "param-row";
    row.innerHTML = `
      <label>${parameterLabels[key] || key}
        <input type="range" min="${bounds[0]}" max="${bounds[1]}" step="0.001" value="${value}" data-slider="${key}">
      </label>
      <div class="param-grid two">
        <label>rad<input type="number" step="0.001" value="${value}" data-param="${key}"></label>
        <label>π 倍数<input type="number" step="0.001" value="${value / Math.PI}" data-pi-param="${key}"></label>
      </div>
    `;
    panel.appendChild(row);
  });

  panel.querySelectorAll("[data-slider]").forEach((input) => {
    input.addEventListener("input", (event) => {
      const key = event.target.dataset.slider;
      setParameter(key, Number(event.target.value));
      syncParameterInputs(key);
      scheduleSimulate();
    });
  });
  panel.querySelectorAll("[data-param]").forEach((input) => {
    input.addEventListener("change", (event) => {
      const key = event.target.dataset.param;
      setParameter(key, Number(event.target.value));
      syncParameterInputs(key);
      scheduleSimulate();
    });
  });
  panel.querySelectorAll("[data-pi-param]").forEach((input) => {
    input.addEventListener("change", (event) => {
      const key = event.target.dataset.piParam;
      setParameter(key, Number(event.target.value) * Math.PI);
      syncParameterInputs(key);
      scheduleSimulate();
    });
  });
}

function setParameter(key, value) {
  const bounds = state.bounds[key] || [-Math.PI, Math.PI];
  state.parameters[key] = Math.min(bounds[1], Math.max(bounds[0], value));
}

function syncParameterInputs(key) {
  const value = state.parameters[key];
  const slider = document.querySelector(`[data-slider="${key}"]`);
  const numeric = document.querySelector(`[data-param="${key}"]`);
  const piNumeric = document.querySelector(`[data-pi-param="${key}"]`);
  if (slider) slider.value = value;
  if (numeric) numeric.value = value.toFixed(6);
  if (piNumeric) piNumeric.value = (value / Math.PI).toFixed(6);
}

function syncGlobalInputs() {
  document.getElementById("outputChannel").value = state.outputChannel;
  document.getElementById("centerThz").value = state.frequency.center_hz / 1e12;
  document.getElementById("spanGhz").value = state.frequency.span_hz / 1e9;
  document.getElementById("passbandGhz").value = state.target.passband_width_hz / 1e9;
  document.getElementById("numPoints").value = state.frequency.num_points;
  document.getElementById("maxIterations").value = state.optimizer.max_iterations;
  document.getElementById("populationSize").value = state.optimizer.population_size;
  document.getElementById("randomSeed").value = state.optimizer.random_seed ?? "";
}

function bindGlobalInputs() {
  document.getElementById("outputChannel").addEventListener("change", (event) => {
    state.outputChannel = event.target.value;
    scheduleSimulate();
  });
  document.getElementById("centerThz").addEventListener("change", (event) => {
    state.frequency.center_hz = Number(event.target.value) * 1e12;
    scheduleSimulate();
  });
  document.getElementById("spanGhz").addEventListener("change", (event) => {
    state.frequency.span_hz = Number(event.target.value) * 1e9;
    scheduleSimulate();
  });
  document.getElementById("passbandGhz").addEventListener("change", (event) => {
    state.target.passband_width_hz = Number(event.target.value) * 1e9;
    scheduleSimulate();
  });
  document.getElementById("numPoints").addEventListener("change", (event) => {
    state.frequency.num_points = Number(event.target.value);
    scheduleSimulate();
  });
  document.getElementById("maxIterations").addEventListener("change", (event) => {
    state.optimizer.max_iterations = Number(event.target.value);
  });
  document.getElementById("populationSize").addEventListener("change", (event) => {
    state.optimizer.population_size = Number(event.target.value);
  });
  document.getElementById("randomSeed").addEventListener("change", (event) => {
    state.optimizer.random_seed = event.target.value === "" ? null : Number(event.target.value);
  });
}

function requestPayload() {
  return {
    parameters: state.parameters,
    frequency: state.frequency,
    target: state.target,
    optimizer: state.optimizer,
    output_channel: state.outputChannel,
  };
}

function scheduleSimulate() {
  clearTimeout(state.autoTimer);
  state.autoTimer = setTimeout(() => simulate({ silent: true }), 180);
}

async function simulate(options = {}) {
  const silent = options.silent === true;
  try {
    if (!silent) log("开始仿真");
    const result = await api("/api/simulate", requestPayload());
    drawChart(result);
    document.getElementById("metrics").textContent = `MSE: ${result.mse.toFixed(6)}`;
    if (!silent) log(`仿真完成，MSE = ${result.mse.toFixed(6)}`);
  } catch (error) {
    log(`仿真失败：${error.message}`);
  }
}

async function optimize() {
  try {
    log("开始优化");
    const result = await api("/api/optimize", requestPayload());
    state.parameters = result.parameters;
    renderParameters();
    drawChart(result);
    document.getElementById("metrics").textContent = `MSE: ${result.best_loss.toFixed(6)}`;
    log(`优化完成，MSE = ${result.best_loss.toFixed(6)}`);
  } catch (error) {
    log(`优化失败：${error.message}`);
  }
}

function drawChart(result) {
  const svg = document.getElementById("chart");
  const width = 900;
  const height = 320;
  const pad = { left: 52, right: 18, top: 18, bottom: 34 };
  const target = downsample(result.target_db, 700);
  const simulated = downsample(result.simulated_db, 700);
  const all = target.concat(simulated);
  const minY = Math.min(...all, -60);
  const maxY = Math.max(...all, 5);
  const x = (index, length) => pad.left + (index / Math.max(1, length - 1)) * (width - pad.left - pad.right);
  const y = (value) => pad.top + ((maxY - value) / Math.max(1, maxY - minY)) * (height - pad.top - pad.bottom);
  const path = (values) => values.map((value, index) => `${index === 0 ? "M" : "L"} ${x(index, values.length).toFixed(2)} ${y(value).toFixed(2)}`).join(" ");
  const frequency = result.frequency_hz || [];
  const center = frequency.length ? 0.5 * (frequency[0] + frequency[frequency.length - 1]) : 0;
  const ticks = [0, 0.25, 0.5, 0.75, 1].map((fraction) => {
    const index = Math.round(fraction * Math.max(0, frequency.length - 1));
    const label = frequency.length ? `${((frequency[index] - center) / 1e9).toFixed(1)} GHz` : "";
    const tickX = pad.left + fraction * (width - pad.left - pad.right);
    return `
      <line class="tick" x1="${tickX.toFixed(2)}" y1="${height - pad.bottom}" x2="${tickX.toFixed(2)}" y2="${height - pad.bottom + 5}"></line>
      <text class="tick-label" x="${tickX.toFixed(2)}" y="${height - 10}" text-anchor="middle">${label}</text>
    `;
  }).join("");
  svg.innerHTML = `
    <line class="axis" x1="${pad.left}" y1="${height - pad.bottom}" x2="${width - pad.right}" y2="${height - pad.bottom}"></line>
    <line class="axis" x1="${pad.left}" y1="${pad.top}" x2="${pad.left}" y2="${height - pad.bottom}"></line>
    <text x="12" y="${pad.top + 8}" fill="#667085">${maxY.toFixed(1)} dB</text>
    <text x="12" y="${height - pad.bottom}" fill="#667085">${minY.toFixed(1)} dB</text>
    ${ticks}
    <text x="${width - pad.right - 170}" y="${pad.top + 12}" fill="#d64545">方波目标</text>
    <text x="${width - pad.right - 76}" y="${pad.top + 12}" fill="#1e66f5">RAMZI</text>
    <path class="target-line" d="${path(target)}"></path>
    <path class="sim-line" d="${path(simulated)}"></path>
  `;
}

function downsample(values, maxPoints) {
  if (!values || values.length <= maxPoints) return values || [];
  const step = (values.length - 1) / (maxPoints - 1);
  return Array.from({ length: maxPoints }, (_, index) => values[Math.round(index * step)]);
}

async function init() {
  const defaults = await fetch("/api/default").then((response) => response.json());
  state.parameters = defaults.parameters;
  state.bounds = defaults.bounds;
  state.rings = defaults.rings;
  state.frequency = defaults.frequency;
  state.target = defaults.target;
  state.optimizer = defaults.optimizer;
  state.outputChannel = defaults.output_channel || "cross";
  renderRingInfo();
  renderParameters();
  syncGlobalInputs();
  bindGlobalInputs();
  document.getElementById("simulateBtn").addEventListener("click", () => simulate());
  document.getElementById("optimizeBtn").addEventListener("click", optimize);
  await simulate();
}

init().catch((error) => log(`初始化失败：${error.message}`));
