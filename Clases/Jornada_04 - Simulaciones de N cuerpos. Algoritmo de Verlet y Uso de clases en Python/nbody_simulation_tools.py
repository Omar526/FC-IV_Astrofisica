# -*- coding: utf-8 -*-
"""
Created on 2024-09-13 10:32:29

@author: Omar Fernández
@mail: omar.fernandez.o@usach.cl

Description: A module for simulating the motion of celestial bodies under gravitational interaction.

This module contains classes to represent individual celestial bodies and to perform N-body simulations.

It includes:
    - CelestialBody: A class representing a celestial object with properties like mass, position, and velocity.
    - NBodySimulation: A class to handle the numerical simulation of multiple celestial bodies interacting gravitationally.
"""

import time
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation, FFMpegWriter
# import matplotlib.image as mpimg
# plt.style.use('dark_background')



class CelestialBody:
    """A class to represent a celestial body in a 2D space with a given mass, position, and velocity.

    Attributes:
        mass (float): The mass of the celestial body.
        color (str): The color representing the celestial body for visualization purposes.
        name (str): The name of the celestial body. If not provided, a unique default name will be generated.
        position (np.array): The initial position of the celestial body as a 2D vector [x, y].
        velocity (np.array): The initial velocity of the celestial body as a 2D vector [vx, vy].
        positions (np.array): Array to store positions at each simulation step.
        velocities (np.array): Array to store velocities at each simulation step.
        accelerations (np.array): Array to store accelerations at each simulation step.
        current_step (int): Tracks the current step in the simulation.
    
    Methods:
        __init__(mass, x, y, vx, vy, name=None, color=None): Initializes the celestial body with its mass, position, velocity, and optional name and color.
        _initialize_arrays(num_steps): Initializes arrays to store positions, velocities, and accelerations for the simulation.
    """
    
    def __init__(self, mass: float, init_position: tuple, init_velocity: tuple, name: str = None, color: str = None) -> None:
        """Constructor method

        Args:
            mass (float): Celestial body mass.
            x (float): Initial position x.
            y (float): Initial position y.
            vx (float): Initial velocity vx.
            vy (float): Initial velocity vy.
            name (str, optional): Name of the celestial body. Defaults to None, in which case a unique name is generated.
            color (str, optional): Color representing the celestial body. Defaults to None.
        """
        self.mass = mass
        self.color = color
        self.name = name

        # Set initial position and velocity
        self.init_position = np.array(init_position, dtype=float)
        self.init_velocity = np.array(init_velocity, dtype=float)

        # Arrays to store positions, velocities, and accelerations will be initialized later
        self.positions = None
        self.velocities = None
        self.accelerations = None

        # Keep track of current step in the simulation
        self.current_step = 0

    def _initialize_arrays(self, num_steps: int):
        """Initializes arrays to store positions, velocities, and accelerations.

        This method should only be called by the NBodySimulation class to prepare for the simulation.

        Args:
            num_steps (int): Number of simulation steps.
        """
        self.positions = np.zeros((num_steps, 2))  # Stores (x, y) at each step
        self.velocities = np.zeros((num_steps, 2))  # Stores (vx, vy) at each step
        self.accelerations = np.zeros((num_steps, 2))  # Stores (ax, ay) at each step

        # Set initial position and velocity
        self.positions[0] = self.init_position
        self.velocities[0] = self.init_velocity


class NBodySimulation:
    def __init__(self, dt: float, num_steps: int) -> None:
        """Constructor method

        Args:
            bodies (list): List of CelestialBody instances
            dt (float): Time step for the simulation
            num_steps (int): Number of steps for the simulation
        """
        #self.bodies = bodies
        self.bodies = {}
        self.dt = dt
        self.num_steps = num_steps

        self.G = 6.67430e-11  # universal gravitational constant (m^3 kg^-1 s^-2)

        self.current_step = 0


    def add_body(self, body):
        """Add a CelestialBody to the simulation and initialize arrays."""
        self.bodies[body.name] = body
        body._initialize_arrays(self.num_steps)

    def compute_gravitational_forces(self, delta=1e4):
        """Calculates the gravitational forces between all celestial bodies."""
        for i, body1 in enumerate(self.bodies.values()):
            acceleration = np.zeros(2)  # Initialize the total force acting on body1
            for j, body2 in enumerate(self.bodies.values()):
                if i != j:
                    # Calculate the vector distance between body1 and body2
                    r_vec = body1.positions[body2.current_step + 1] - body2.positions[body1.current_step + 1]
                    r_mag = np.linalg.norm(r_vec)  # Magnitude of the distance vector
                
                    if abs(r_mag) > delta:
                        # Compute the gravitational force magnitude and direction
                        acceleration += body2.mass / r_mag**3 * r_vec
            
            acceleration *= - self.G
            
            # Store the resulting acceleration for body1 (F = m * a => a = F / m)
            body1.accelerations[body1.current_step + 1] = acceleration

    def update(self):
        """Updates the positions and velocities of celestial bodies using the Velocity Verlet method."""
        # First, update positions based on the current accelerations
        for body in self.bodies.values():
            current_step = body.current_step
            
            if current_step < self.num_steps - 1:  # Ensure we do not exceed the number of steps
                # 1. Update position: r_{n+1} = r_n + v_n * dt + 0.5 * a_n * dt^2
                body.positions[current_step + 1] = (
                    body.positions[current_step]
                    + body.velocities[current_step] * self.dt
                    + 0.5 * body.accelerations[current_step] * self.dt**2
                )

        # Recompute forces to update accelerations for the next step
        self.compute_gravitational_forces()
        
        # Then, update velocities based on the updated accelerations
        for body in self.bodies.values():
            current_step = body.current_step
            
            if current_step < self.num_steps - 1:  # Ensure we do not exceed the number of steps
                # 2. Calculate the new acceleration for the next step
                new_acceleration = body.accelerations[current_step + 1]

                # 3. Update velocity: v_{n+1} = v_n + 0.5 * (a_n + a_{n+1}) * dt
                body.velocities[current_step + 1] = (
                    body.velocities[current_step]
                    + 0.5 * (body.accelerations[current_step] + new_acceleration) * self.dt
                )

            # 5. Move to the next step
            body.current_step += 1

    def run_simulation(self):
        """Run the simulation for all steps using Velocity Verlet method."""

        print("Start simulation")

        for step in range(self.num_steps - 1):
            # Step 1: Compute gravitational forces to update accelerations
            self.compute_gravitational_forces()

            # Step 2: Update positions and velocities for all bodies
            self.update()

            # Print progress at each step in console
            print(f"Step {step + 2}/{self.num_steps} completed.", end="\r")
        
        print() # this ensures that the next print does not continue to overwrite on the same line
        print("Completed simulation!")


    def compute_center_of_mass(self):
        """Computes the center of mass of the system."""
        total_mass = sum(body.mass for body in self.bodies.values())
        weighted_positions = sum(body.mass * body.positions for body in self.bodies.values())
        center_of_mass = weighted_positions / total_mass

        # print(center_of_mass)
        # print(center_of_mass.shape)
        return center_of_mass
    
    def shift_to_center_of_mass(self):
        """Shifts the reference frame to the center of mass."""
        center_of_mass = self.compute_center_of_mass()
        for body in self.bodies.values():
            body.positions = body.positions - center_of_mass


# ----------------- EXAMPLE OF USE -----------------------

if __name__ == '__main__':

    # Definition of a binary star system
    star1 = CelestialBody(mass=1.989e30, x=0, y=0, vx=0, vy=0, name="Star1")
    star2 = CelestialBody(mass=1.989e30, x=1.5e11, y=1.5e11, vx=0, vy=25_000, name="Star2")

    # Distant asteroid approaching the system.
    asteroid = CelestialBody(mass=1e15, x=5e11, y=-6e11, vx=-5_000, vy=15_000, name="Asteroid")

    # List of celestial bodies
    bodies = [star1, star2, asteroid]

    # Simulation parameters
    dt = 1 * 3600 * 1 * 24 * 2 # timestep de un día (seconds)
    num_steps = 275 # (total time = dt*num_steps)

    # Create an instance to set up the simulation
    simulation = NBodySimulation(dt=dt, num_steps=num_steps)

    # add celestial bodies to simulation
    for body in bodies:
        simulation.add_body(body=body)

    # run simulation
    simulation.run_simulation()

    # Change of reference system to the center of mass
    simulation.shift_to_center_of_mass()
