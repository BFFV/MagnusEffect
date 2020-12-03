// Simulation Settings
const simulationContainer = d3.select('#simulation');
const simulate = d3.select('#simulationButton')
  .on('click', startSimulation)
const selectBall = d3.select('#ballType')
  .on('change', loadBall)
const maxWidth = parseInt(simulationContainer.style('width'), 10);
const maxHeight = maxWidth * 0.5;
const initialPosition = {'x': 0, 'y': (maxHeight / 2) - 30}
const maxSpeed = 50;
const maxSpeedValue = 80;
const ballSettings = {
  football: {M: 1, R: 0.1},
  tennis: {M: 0.5, R: 0.05},
};

// Simulation Variables
let simulating = false;
let angle = 0;
let ballType = d3.select('#ballType').property('value');
let rotation = parseInt(d3.select('#rotationSpeed').property('value'), 10);
let ballHeight = parseInt(d3.select('#height').property('value'), 10);
let ballSpeed = {'x': parseInt(d3.select('#ballSpeed').property('value'), 10), 'y': 0};
let initialSpeed = ballSpeed.x;
let prev = {'x': initialPosition.x, 'y': initialPosition.y};
let interval = 0;
const position = {'x': initialPosition.x, 'y': initialPosition.y}

// Scales
const speedScale = d3.scaleLinear()
  .domain([0, maxSpeedValue])
  .range([0, maxSpeed])
const speedInverse = d3.scaleLinear()
.domain([0, maxSpeed])
.range([0, maxSpeedValue])

// Load Simulation
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

// Reset Simulation
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
}

// Load Ball Image
function loadBall() {
  ballType = selectBall.property('value');
  if (!simulating) {
    d3.select('#ball')
      .attr('xlink:href', `data/${ballType}.png`)
  }
}

// Handle Start Simulation
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
  prev.x = initialPosition.x;
  prev.y = initialPosition.y;
  angle = 0;
  rotation = parseInt(d3.select('#rotationSpeed').property('value'), 10);
  ballHeight = parseInt(d3.select('#height').property('value'), 10);
  simulating = true;
  interval = setInterval(simulation, 20);
}

// Checks Ball Collisions
function checkCollision() {
  if (position.x < 0 || position.x > maxWidth) {
    return true;
  } else if (position.y < - 60 || position.y > maxHeight) {
    return true;
  }
  return false;
}

// Calculate Air Density
function getDensity() {
  return 1;
}

// Calculate Next Position
function computeNextPosition() {
  if (prev.x !== position.x) {
    angle += Math.atan((position.y - prev.y) / (position.x - prev.x));
  }
  prev.y = position.y;
  prev.x = position.x;
  const F = 1 * getDensity() * Math.PHI * rotation * Math.sqrt(ballSpeed.x^2 + ballSpeed.y^2) * ballSettings[ballType].R^3;

  // Calcular velocidad a partir de fuerza de magnus y dinámica
  ballSpeed.x += -F * Math.sin(angle) / ballSettings[ballType].M;
  ballSpeed.y += F * Math.cos(angle) / ballSettings[ballType].M;

  // Buscar forma de realizarlo correctamente

  position.x += speedScale(ballSpeed.x);
  position.y += speedScale(ballSpeed.y);
}

// Ball Movement
function simulation() {
  computeNextPosition();
  d3.select('#ball')
    .attr('x', position.x)
    .attr('y', position.y)
  if (checkCollision()) {
    resetSimulation();
  }
}

// Initialize
loadSimulation();
