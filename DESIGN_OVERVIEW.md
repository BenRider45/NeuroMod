# BCNNeuroModel Design Overview

## System Purpose

BCNNeuroModel is a computational neuroscience framework for simulating biophysically realistic neuron models, specifically targeting HVC-RA (high vocal center robust nucleus of the arcopallium) neurons found in songbird brain circuits. The system implements a two-compartment Hodgkin-Huxley style model capable of simulating the complex electrophysiological dynamics of these neurons.

## Architectural Principles

The system follows a **polymorphic solver architecture** where numerical integration methods (specifically RK4) are decoupled from the specific neuron models being simulated. This design allows:

1. **Flexibility**: Multiple neuron models can use the same solver
2. **Testability**: Simple models (integrate-and-fire) can verify solver correctness before applying to complex models
3. **Extensibility**: New neuron models or integration methods can be added without modifying existing code

## Core Component Architecture

### 1. Abstract ODE System Interface (`OdeAble`)

**Location**: `src/BCNNeuroModel/OdeAble.hpp`

The `OdeAble` base class defines the contract that all dynamical systems must fulfill:

```
┌─────────────────────────────────────────┐
│            OdeAble                      │
├─────────────────────────────────────────┤
│ + ode_num: int                          │
│ + f(t, y) → State                       │
│ + verify_state(y) → bool                │
│ + post_step_process(t, y) → void        │
└─────────────────────────────────────────┘
```

**Key Responsibilities**:
- `f(t, y)`: Defines the system of ODEs (dy/dt = f(t, y))
- `verify_state(y)`: Validates state vector dimensionality
- `post_step_process(t, y)`: Post-integration processing (e.g., spike detection, logging)
- `ode_num`: Number of state variables in the system


### 2. Numerical Integration Engine (`RK4`)

**Location**: `src/BCNNeuroModel/Rk4.hpp`, `Rk4.cpp`

The RK4 class implements the fourth-order Runge-Kutta numerical integration method:

```
┌─────────────────────────────────────────┐
│              RK4                        │
├─────────────────────────────────────────┤
│ - _h: double (timestep)                 │
│ + Simulate(system, y_0, t_0, t_f)      │
│ - Step(system, t, y) → State           │
└─────────────────────────────────────────┘
```

**Integration Algorithm** (src/BCNNeuroModel/Rk4.cpp:29-49):

For each timestep, computes four intermediate slopes:
```
K1 = f(t, y)
K2 = f(t + h/2, y + h/2·K1)
K3 = f(t + h/2, y + h/2·K2)
K4 = f(t + h, y + h·K3)

y_next = y + (h/6)·(K1 + 2·K2 + 2·K3 + K4)
```

**Output Format**: Returns an `Eigen::MatrixXd` where:
- Rows 0 to (ode_num-1): State variables over time
- Row ode_num: Time vector
- Columns: Individual timesteps

### 3. Two-Compartment HVC-RA Neuron Model (`HVCRA`)

**Location**: `src/BCNNeuroModel/HVCRA.hpp`, `HVCRA.cpp`

The HVCRA model implements a detailed biophysical model with somatic and dendritic compartments:

```
┌─────────────────────────────────────────────────────────┐
│                  HVCRA : OdeAble                        │
├─────────────────────────────────────────────────────────┤
│ State Variables (11 total):                             │
│  [0] V_s  - Somatic membrane voltage                    │
│  [1] V_d  - Dendritic membrane voltage                  │
│  [2] h    - Na+ inactivation gate                       │
│  [3] n    - K+ activation gate                          │
│  [4] r    - Ca2+ activation gate                        │
│  [5] c    - Ca-K activation gate                        │
│  [6] Ca   - Intracellular calcium concentration         │
│  [7-10]   - Synaptic conductances (exc/inh, soma/dend) │
├─────────────────────────────────────────────────────────┤
│ Ionic Currents:                                         │
│  Soma: I_L, I_Na, I_Kdr + synaptic + external          │
│  Dendrite: I_L, I_Ca, I_CaK + synaptic + external      │
├─────────────────────────────────────────────────────────┤
│ Key Methods:                                            │
│  + f(t, y) → computes all 11 derivatives                │
│  + post_step_process() → spike detection                │
└─────────────────────────────────────────────────────────┘
```

#### Biophysical Model Details

**Compartment Coupling** (src/BCNNeuroModel/HVCRA.cpp:62-66):

The two compartments are electrically coupled through a resistive connection:
```
dV_s/dt includes: (V_d - V_s)/(R_c · C_m · A_s)
dV_d/dt includes: (V_s - V_d)/(R_c · C_m · A_d)
```

**Somatic Ionic Currents** (HVCRA.cpp:34-45):
- **Leak current**: `I_sL = -G_L·(V_s - E_L)`
- **Fast Na+ current**: `I_sNa = -G_Na·m³·h·(V_s - E_Na)` (action potential upstroke)
- **Delayed rectifier K+**: `I_sKdr = -G_Kdr·n⁴·(V_s - E_K)` (repolarization)
- **Synaptic currents**: Excitatory and inhibitory with exponential decay

**Dendritic Ionic Currents** (HVCRA.cpp:46-55):
- **Leak current**: `I_dL = -G_dL·(V_d - E_L)`
- **High-threshold Ca2+**: `I_dCa = -G_Ca·r²·(V_d - E_Ca)` (dendritic spikes)
- **Ca-activated K+**: `I_dCaK = -(G_CaK·c)/(1 + 6/Ca)·(V_d - E_K)` (afterhyperpolarization)
- **Synaptic currents**: Excitatory and inhibitory

**Gating Variable Dynamics**:

All gating variables follow first-order kinetics: `dg/dt = (g_∞(V) - g)/τ(V)`

Where voltage-dependent functions are defined in HVCRA.hpp:95-112:
- Na+ inactivation: `h_∞(V) = 1/(1 + exp((V + 45)/7))`
- K+ activation: `n_∞(V) = 1/(1 + exp(-(V + 35)/10))`
- Ca2+ activation: `r_∞(V) = 1/(1 + exp(-(V + 5)/10))`
- Ca-K activation: `c_∞(V) = 1/(1 + exp(-(V - 10)/7))`

**External Input Interface**:

The model accepts two function pointers for time-dependent external currents:
- `I_sExt(t)`: Current injection to soma
- `I_dExt(t)`: Current injection to dendrite

This allows flexible stimulus protocols (step currents, ramps, sinusoids, etc.) defined at runtime.

### 4. Integrate-and-Fire Verification Model (`IFNeuron`)

**Location**: `src/BCNNeuroModel/IFNeuron.hpp`

A simple single-compartment integrate-and-fire neuron used to verify the RK4 solver:

```
┌─────────────────────────────────────────┐
│          IFNeuron : OdeAble             │
├─────────────────────────────────────────┤
│ State Variables (1):                    │
│  [0] V - Membrane voltage               │
├─────────────────────────────────────────┤
│ Dynamics:                               │
│  dV/dt = (L - V + I₀·(1 + sin(t)))/τ   │
└─────────────────────────────────────────┘
```

**Implementation** (IFNeuron.hpp:17-24):

The model implements a leaky integrate-and-fire with sinusoidal drive:
```
dV/dt = (L - V + I₀(1 + sin(t))) / τ
```

Where:
- `L`: Leak reversal potential (-70 mV typical)
- `I₀`: Drive current amplitude
- `τ`: Membrane time constant (20 ms typical)

The `main.cpp:62-89` demonstrates testing with multiple timesteps (h = 10⁻⁵, 10⁻⁴, 10⁻³, 10⁻²).

### 5. Data Export Utility (`Exporto`)

**Location**: `src/BCNNeuroModel/Exporto.hpp`, `Exporto.cpp`

Handles exporting simulation data to MATLAB/Octave format:

```
┌─────────────────────────────────────────┐
│            Exporto                      │
├─────────────────────────────────────────┤
│ - _ExportDir: filesystem::path          │
│ + writeMatrixToOctave()                 │
│ + writePlotScript()                     │
│ + writeTimeStepVector()                 │
└─────────────────────────────────────────┘
```

**Output Format**: Creates `.m` files containing MATLAB/Octave arrays for post-simulation analysis and visualization.

## Data Flow Architecture

### Complete Simulation Pipeline

```
┌──────────────────┐
│  User Code       │  - Defines constants, initial conditions
│  (main.cpp)      │  - Creates neuron model & external currents
└────────┬─────────┘
         │
         ├─────────────────────────────────────────────┐
         │                                             │
         v                                             v
┌─────────────────┐                          ┌─────────────────┐
│  HVCRA Model    │                          │  IFNeuron       │
│  (11 ODEs)      │                          │  (1 ODE)        │
└────────┬────────┘                          └────────┬────────┘
         │                                             │
         │ implements OdeAble::f()                     │
         └─────────────────┬───────────────────────────┘
                           │
                           v
                  ┌─────────────────┐
                  │  RK4 Solver     │
                  │  - Timestep: h  │
                  │  - Method: K1-K4│
                  └────────┬────────┘
                           │
                           │ For each timestep:
                           │  1. Compute K1 = f(t, y)
                           │  2. Compute K2 = f(t+h/2, y+h/2·K1)
                           │  3. Compute K3 = f(t+h/2, y+h/2·K2)
                           │  4. Compute K4 = f(t+h, y+h·K3)
                           │  5. Update: y += h/6·(K1+2K2+2K3+K4)
                           │  6. Call post_step_process()
                           │
                           v
                  ┌─────────────────┐
                  │ Eigen::MatrixXd │
                  │  Rows: States   │
                  │  Cols: Time     │
                  └────────┬────────┘
                           │
                           v
                  ┌─────────────────┐
                  │    Exporto      │
                  │ Writes .m files │
                  └────────┬────────┘
                           │
                           v
                  ┌─────────────────┐
                  │ MATLAB/Octave   │
                  │ Visualization   │
                  └─────────────────┘
```

### State Vector Flow in HVCRA

```
Initial State (y_0)
     │
     ├─> V_s, V_d ──┐
     │              │
     ├─> h, n ──────┼──> Voltage-dependent rates
     │              │    (m_∞, h_∞, n_∞, r_∞, c_∞)
     ├─> r, c ──────┤    (τ_h, τ_n, τ_r, τ_c)
     │              │
     ├─> Ca ────────┼──> I_CaK modulation
     │              │
     └─> g_syn ─────┘
                    │
                    v
              ┌──────────────┐
              │  f(t, y)     │
              │  Computes:   │
              │  - Currents  │
              │  - dV/dt     │
              │  - dg/dt     │
              └──────┬───────┘
                     │
                     v
              ┌──────────────┐
              │  RK4 Step    │
              │  y_next      │
              └──────┬───────┘
                     │
                     v
         ┌───────────────────────┐
         │ post_step_process()   │
         │ - Spike detection     │
         │ - V_s > -10 mV?       │
         │ - V_d > -10 mV?       │
         └───────────────────────┘
```

## Dependency Structure

```
┌─────────────────────────────────────────────────────┐
│                    Eigen3                           │
│          (Linear algebra, matrix operations)        │
└──────────────────┬──────────────────────────────────┘
                   │
                   │ provides VectorXd (State)
                   │
         ┌─────────┴─────────┐
         │                   │
         v                   v
┌─────────────────┐   ┌─────────────────┐
│   OdeAble       │   │   Exporto       │
│  (interface)    │   │  (utility)      │
└────────┬────────┘   └─────────────────┘
         │
         │ inherited by
         │
    ┌────┴────┐
    │         │
    v         v
┌─────────┐ ┌─────────┐
│ HVCRA   │ │IFNeuron │
└────┬────┘ └────┬────┘
     │           │
     │           │ both used by
     └─────┬─────┘
           │
           v
    ┌─────────────┐
    │    RK4      │
    │  (solver)   │
    └─────────────┘
```

**External Dependencies**:
- **Eigen3**: All vector/matrix operations (VectorXd = State type)
- **gflags**: Command-line argument parsing (currently unused but linked)
- **GoogleTest**: Unit testing framework

## Key Design Decisions

### 1. Polymorphism Over Templates

The system uses runtime polymorphism with `OdeAble` as a virtual base class. Neuron models inherit from `OdeAble` and the `RK4` solver operates on `OdeAble` references.

### 2. Function Pointers for External Currents

External currents use `std::function<double(double)>` function objects. Example (main.cpp:7-11):
```cpp
double DendriteExternalCurrent(double t) {
  return t >= 20e-3 && t <= 40e-3 ? 1.4e-9 : 0;
}
```

### 3. Post-Step Processing Hook

The `post_step_process()` method is called after each integration step. Current uses:
- Spike detection (HVCRA.cpp:97-104)

### 4. State Validation

Each model implements `verify_state()` which checks that the state vector has the correct dimensionality for that model.

## Verification Strategy: Integrate-and-Fire Testing

The verification workflow in main.cpp:62-89 tests RK4 with `IFNeuron` at multiple timesteps:

```cpp
RK4 solver1(1e-5);   // h = 10 μs
RK4 solver2(1e-4);   // h = 100 μs
RK4 solver3(1e-3);   // h = 1 ms
RK4 solver4(1e-2);   // h = 10 ms
```

The outputs are compared to verify 4th-order convergence (error ∝ h⁴).

## Future Extensions

The architecture supports several planned enhancements:

1. **Network Simulation**: Multiple coupled HVCRA neurons
2. **Additional Models**: Interneurons, other brain regions
3. **Adaptive Timesteps**: Variable h based on solution dynamics
4. **Alternative Solvers**: Higher-order methods, implicit methods for stiff systems
5. **Real-time Visualization**: Integration with plotting libraries

## File Organization Summary

```
src/BCNNeuroModel/
├── OdeAble.hpp          # Abstract ODE system interface
├── Rk4.hpp/cpp          # RK4 numerical solver
├── HVCRA.hpp/cpp        # Two-compartment HH model
├── IFNeuron.hpp         # Verification model
├── Exporto.hpp/cpp      # MATLAB/Octave export
└── main.cpp             # Demo and verification runs

tests/
└── HVCRA_Test.cpp       # Unit tests (GoogleTest)
```

## Typical Usage Pattern

```cpp
// 1. Define physiological constants
HVCRA_CONSTANTS consts(...);

// 2. Define external current functions
auto I_soma = [](double t) { return /* current */; };
auto I_dend = [](double t) { return /* current */; };

// 3. Create neuron model
HVCRA neuron(consts, I_soma, I_dend, G, id);

// 4. Set initial conditions
State y_0(11);
y_0 << -75e-3, -75e-3, 0, 0, 0, 0, 0, 0, 0, 0, 0;

// 5. Create solver with timestep
RK4 solver(1e-5);  // h = 10 μs

// 6. Run simulation
Eigen::MatrixXd data = solver.Simulate(neuron, y_0, 0, 100e-3);

// 7. Export results
Exporto exp("./output");
exp.writeMatrixToOctave(data, "results.m", "sim_data");
```

---

**Document Version**: 1.0
**Last Updated**: 2026-01-13
**Corresponds to**: BCNNeuroModel commit current state
