# -*- coding: utf-8 -*-
"""
Created on 2024-11-07 16:41:31

@author: Omar Fernández
@mail: omar.fernandez.o@usach.cl

Description: To generate galaxies
"""

import numpy as np

def generate_spiral_galaxy(
    Lx,
    Ly,
    Nx,
    Ny,
    rho,
    center_x,
    center_y,
    num_arms,
    arm_spread,
    density_peak,
    decay_rate,
    central_mass_peak,
    background_density,
):
    # Unidades: Lx y Ly en metros, pero x, y serán convertidos a kpc más adelante
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    X, Y = np.meshgrid(x, y)

    # Convertir a coordenadas polares centradas en el centro de la galaxia (en kpc)
    R = np.sqrt((X - center_x) ** 2 + (Y - center_y) ** 2) / 3.086e19  # Convertir a kpc
    theta = np.arctan2(Y - center_y, X - center_x)

    # Añadir densidad de fondo (en kg/m²)
    rho += background_density * np.ones_like(R)

    # Crear la estructura de brazos espirales
    for arm in range(num_arms):
        arm_angle = 2 * np.pi * arm / num_arms
        spiral_theta = arm_angle + np.log(R + 1) * arm_spread

        # Agregar una distribución de masa gaussiana a lo largo de los brazos
        gaussian = (
            density_peak
            * np.exp(-decay_rate * R**2)
            * np.exp(-(((theta - spiral_theta) % (2 * np.pi)) ** 2) / (2 * 0.3**2))
        )
        rho += gaussian

    # Agregar masa al centro
    central_gaussian = central_mass_peak * np.exp(
        -(R**2) / (2 * (0.3**2))
    )  # Distribución gaussiana centrada
    rho += central_gaussian

    return rho, x, y
