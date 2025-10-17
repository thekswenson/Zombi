"""
Global random number generators.  There is one for each type of simulator.
"""
import numpy as np
import random


RNG = random.Random()
NPRNG = np.random.default_rng()

T_RNG = random.Random(11)
T_NPRNG = np.random.default_rng(11)
def seed_T_generators(seed):
    """
    Seed the global random number generators for the tree simulator.
    A seed of 0 means to not reseed.
    """
    global T_RNG, T_NPRNG
    if seed != 0 and seed != "0":
        T_RNG = random.Random(seed)
        T_NPRNG = np.random.default_rng(seed)


G_RNG = random.Random(11)
G_NPRNG = np.random.default_rng(11)
def seed_G_generators(seed):
    """
    Seed the global random number generators for the genome simulator.
    A seed of 0 means to not reseed.
    """
    global G_RNG, G_NPRNG
    if seed != 0 and seed != "0":
        G_RNG = random.Random(seed)
        G_NPRNG = np.random.default_rng(seed)


S_RNG = random.Random(11)
S_NPRNG = np.random.default_rng(11)
def seed_S_generators(seed):
    """
    Seed the global random number generators for the species simulator.
    A seed of 0 means to not reseed.
    """
    global S_RNG, S_NPRNG
    if seed != 0 and seed != "0":
        S_RNG = random.Random(seed)
        S_NPRNG = np.random.default_rng(seed)
