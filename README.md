# Ultra-Precision Pendulum Simulation

A high-precision numerical simulation of a mathematical pendulum with perfect energy conservation for undamped systems. This C++ implementation achieves machine-precision energy conservation through a combination of advanced numerical techniques.

## Features

- **Perfect Energy Conservation**: Achieves machine-precision energy conservation for undamped pendulum simulations
- **High-Order Integration**: Uses an 8th-order symplectic integrator based on Yoshida coefficients
- **Analytical Energy Enforcement**: Automatically rescales velocities to enforce exact energy conservation
- **Extended Precision**: Uses `long double` for maximum floating-point precision
- **High-Precision Trigonometry**: Custom implementation of trigonometric functions with Taylor series for small angles
- **Adaptive Precision**: Different algorithms for small vs. large angles to maximize precision
- **Damping Support**: Accurate modeling of damped pendulum systems with analytical damping

## Mathematical Background

The pendulum dynamics are governed by the equation:

```
θ'' + (damping/mass) * θ' + (g/L) * sin(θ) = 0
```

Where:
- θ is the angle from the vertical position
- θ' is the angular velocity
- θ'' is the angular acceleration
- g is the gravitational acceleration
- L is the pendulum length
- damping is the damping coefficient

For an undamped pendulum, the total energy should remain constant throughout the simulation:

```
E = KE + PE = (1/2) * m * L² * θ'² + m * g * L * (1 - cos(θ))
```

## Implementation Details

The implementation achieves perfect energy conservation through several key techniques:

1. **8th-Order Symplectic Integration**: Uses high-order symplectic integrators that inherently preserve the Hamiltonian structure of the system.

2. **Analytical Energy Enforcement**: After each integration step, the velocity is rescaled to ensure the total energy matches the initial energy exactly.

3. **Extended Precision Arithmetic**: Uses `long double` (typically 80-bit precision) for all calculations.

4. **Kahan Summation**: Implements compensated summation to minimize round-off errors during energy calculations.

5. **Small-Angle Optimization**: For small angles, uses high-precision Taylor series expansions instead of standard library trigonometric functions.

## Usage

### Basic Usage

```cpp
#include <iostream>

int main() {
    // Simulation parameters
    long double length = 1.0L;         // pendulum length (m)
    long double g = 9.81L;             // gravitational acceleration (m/s²)
    long double mass = 1.0L;           // bob mass (kg)
    long double damping = 0.0L;        // damping coefficient (kg/s)
    long double theta0 = M_PI/4.0L;    // initial angle (rad)
    long double omega0 = 0.0L;         // initial angular velocity (rad/s)
    long double dt = 0.0001L;          // time step (s)
    long double total_time = 10.0L;    // total simulation time (s)
    
    // Create and run the simulation
    UltraPrecisionPendulum pendulum(length, g, mass, damping, theta0, omega0, dt);
    pendulum.runSimulation(total_time, "pendulum_data.csv", 0.01L);
    
    return 0;
}
```

### Simulating a Damped Pendulum

```cpp
// Create a damped pendulum with damping coefficient 0.1
UltraPrecisionPendulum dampedPendulum(1.0L, 9.81L, 1.0L, 0.1L, M_PI/4.0L, 0.0L, 0.0001L);
dampedPendulum.runSimulation(10.0L, "damped_pendulum_data.csv");
```

### Manual Stepping

```cpp
UltraPrecisionPendulum pendulum(1.0L, 9.81L, 1.0L, 0.0L, M_PI/4.0L, 0.0L, 0.0001L);

// Advance the simulation step by step
for (int i = 0; i < 1000; i++) {
    pendulum.step();
    
    // Access the current state
    std::cout << "Time: " << pendulum.getTime() 
              << ", Angle: " << pendulum.getTheta() 
              << ", Energy: " << pendulum.getTotalEnergy() << std::endl;
}
```

## Compiling and Running

Compile with:

```bash
g++ -o pendulum src/UltraPrecisionPendulum.cpp -lm -O3 -std=c++17
```

Run with:

```bash
./pendulum
```

The simulation will generate CSV files with the results.

## Output Format

The output CSV file contains the following columns:

1. **Time**: Simulation time in seconds
2. **Angle**: Pendulum angle in radians
3. **AngularVelocity**: Angular velocity in radians/second
4. **KineticEnergy**: Kinetic energy in Joules
5. **PotentialEnergy**: Potential energy in Joules
6. **TotalEnergy**: Total energy in Joules
7. **RelativeEnergyError**: Relative error in energy conservation

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
