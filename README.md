# Mie_Nanoparticle_Concentration_from_Optical_Density_Sphere
Calculate the density of your spherical (nano)particles from their permittivity, size, and optical density using Mie theory and the Beer-Lambert law

This script reads the relative permittivity data of a material and uses Mie theory to compute the extinction cross-section for a spherical particle of that material embedded in a lossless medium. The permittivity file has to be a tab-delimited text file with three columns: energy (in eV), epsilon1, epsilon2.

From the particle size, the script computes the particle concentration for a solution with a specified optical density at resonance.
