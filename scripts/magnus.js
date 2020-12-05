// Ode45 solver
const ode45 = require('ode45-cash-karp');

// Simulation settings
const simulationContainer = d3.select('#simulation');
const simulate = d3.select('#simulationButton')
  .on('click', startSimulation)
const selectBall = d3.select('#ballType')
  .on('change', loadBall)
const maxWidth = parseInt(simulationContainer.style('width'), 10);
const maxHeight = maxWidth * 0.5;
const initialPosition = {'x': 0, 'y': (maxHeight / 2) - 30};
const ballSettings = {
  basketball: {M: 0.61, R: 0.12},
  football: {M: 0.43, R: 0.11},
  volleyball: {M: 0.27, R: 0.105},
  baseball: {M: 0.142, R: 0.036},
  tennis: {M: 0.114, R: 0.033},
  golf: {M: 0.046, R: 0.021},
  tableTennis: {M: 0.0028, R: 0.02}
};

// Speed Scale
const speedScale = d3.scaleLinear()
  .domain([0, 100])
  .range([0, maxWidth * 0.12])

// Simulation variables
let simulating = false;
let ballType = 'football';
let rotation = 0;
let ballHeight = 0;
let ballSpeed = {'x': 0, 'y': 0};
let temperature = 0;
let initialSpeed = 0;
let angle = 0;
let interval = 0;
let integrator = ode45([0, 0], speed, 0, 1);
const position = {'x': initialPosition.x, 'y': initialPosition.y};

// Air density
function getDensity(h, T) {
  const g = 9.81;
  const T0 = 288.15;
  const M = 0.0289654;
  const L = 0.0065;
  const R = 8.31447;
  const p0 = 101325;
  const p = p0 * ((1 - (L * h / T0)) ** ((g * M) / (R * L)));
  const Rs = 287.058;
  return p / (T * Rs);
}

// ODEs for the ball speed
function speed(dydt, y, t) {
  const v = Math.sqrt(y[0] ** 2 + y[1] ** 2);
  const density = getDensity(ballHeight, temperature);
  const ball = ballSettings[ballType];
  const area = Math.PI * (ball.R ** 2);
  const mass = ball.M;
  const Cd = 0.508 + (1 / ((22.503 + 4.196 * ((initialSpeed / Math.abs(rotation)) ** 2.5)) ** 0.4));
  const Cl = 1 / (2.202 + 0.981 * (initialSpeed / Math.abs(rotation)));
  let spin = 1;
  if (rotation < 0) {
    spin = -1;
  }
  dydt[0] = - (1 / (2 * mass)) * density * area * v * (Cd * y[0] + Cl * y[1] * spin);
  dydt[1] = - (1 / (2 * mass)) * density * area * v * (Cd * y[1] - Cl * y[0] * spin);
}

// Load simulation
function loadSimulation() {
  const svg = simulationContainer
    .append('svg')
      .attr('width', maxWidth)
      .attr('height', maxHeight)
  svg.append('rect')
    .attr('width', '100%')
    .attr('height', '100%')
    .attr('fill', 'none')
    .attr('stroke', 'black')
    .attr('stroke-width', 3)
  svg
    .append('image')
      .attr('xlink:href', `data/${ballType}.png`)
      .attr('x', position.x)
      .attr('y', position.y)
      .attr('width', '60px')
      .attr('height', '60px')
      .attr('id', 'ball')
}

// Reset simulation
function resetSimulation() {
  simulating = false;
  clearInterval(interval);
  simulate
    .style('background-color', '#191970')
    .text('Simular')
  position.x = initialPosition.x;
  position.y = initialPosition.y;
  d3.select('#ball')
    .attr('x', position.x)
    .attr('y', position.y)
    .attr('xlink:href', `data/${ballType}.png`)
    .attr('transform-origin', `${position.x + 30} ${position.y + 30}`)
    .attr('transform', `rotate (0)`)
}

// Load ball image
function loadBall() {
  ballType = selectBall.property('value');
  if (!simulating) {
    d3.select('#ball')
      .attr('xlink:href', `data/${ballType}.png`)
  }
}

// Handle start simulation
function startSimulation() {
  if (simulating) {
    return resetSimulation();
  }
  simulate
    .style('background-color', 'red')
    .text('Detener')
  ballSpeed.x = parseInt(d3.select('#ballSpeed').property('value'), 10);
  ballSpeed.y = 0;
  initialSpeed = ballSpeed.x;
  rotation = parseInt(d3.select('#rotationSpeed').property('value'), 10);
  ballHeight = parseInt(d3.select('#height').property('value'), 10);
  temperature = parseInt(d3.select('#temperature').property('value'), 10);
  integrator = ode45([initialSpeed, 0], speed, 0, 1);
  simulating = true;
  interval = setInterval(simulation, 10);
}

// Checks ball collisions
function checkCollision() {
  if (position.x < 0 || position.x > maxWidth) {
    return true;
  } else if (position.y < - 60 || position.y > maxHeight) {
    return true;
  }
  return false;
}

// Calculate next position
function computeNextPosition() {
  try {
    integrator.step();
  } catch {
    resetSimulation();
  }
  ballSpeed.x = integrator.y[0];
  ballSpeed.y = integrator.y[1];
  position.x += speedScale(ballSpeed.x);
  position.y += speedScale(ballSpeed.y);
  angle += (rotation * (180 / Math.PI)) / 100;
}

// Ball movement
function simulation() {
  computeNextPosition();
  d3.select('#ball')
    .attr('x', position.x)
    .attr('y', position.y)
    .attr('transform-origin', `${position.x + 30} ${position.y + 30}`)
    .attr('transform', `rotate (${angle})`)
  if (checkCollision()) {
    resetSimulation();
  }
}

// Initialize:
loadSimulation();
