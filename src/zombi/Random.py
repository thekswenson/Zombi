"""
Global random number generators.  There is one for each type of simulator.

In many spots, we sort lists seemingly for no good reason, but this is
done so that the results remain consistant when using the same seed.
(to speed things up, we could sort only when a seed is given)
"""
import numpy as np
import random


RNG = random.Random()
NPRNG = np.random.default_rng()

def T_RNG() -> random.Random:
    """ Global random number generator for the tree simulator. """
    return T_RNGvar
def T_NPRNG() -> np.random.Generator:
    """ Global numpy random number generator for the tree simulator. """
    return T_NPRNGvar

T_RNGvar = random.Random()
T_NPRNGvar = np.random.default_rng()
def seed_T_generators(seed: int):
    """
    Seed the global random number generators for the tree simulator.
    A seed of 0 means to not reseed.
    """
    global T_RNGvar, T_NPRNGvar
    if seed != 0 and seed is not None:
        print("Seeding T generators with seed", seed)
        T_RNGvar = random.Random(seed)
        T_NPRNGvar = np.random.default_rng(seed)


def G_RNG() -> random.Random:
    """ Global random number generator for the genome simulator. """
    return G_RNGvar
def G_NPRNG() -> np.random.Generator:
    """ Global numpy random number generator for the genome simulator. """
    return G_NPRNGvar

G_RNGvar = random.Random()
G_NPRNGvar = np.random.default_rng()
def seed_G_generators(seed: int):
    """
    Seed the global random number generators for the genome simulator.
    A seed of 0 means to not reseed.
    """
    global G_RNGvar, G_NPRNGvar
    if seed != 0 and seed is not None:
        G_RNGvar = random.Random(seed)
        G_NPRNGvar = np.random.default_rng(seed)


def S_RNG() -> random.Random:
    """ Global random number generator for the species simulator. """
    return S_RNGvar
def S_NPRNG() -> np.random.Generator:
    """ Global numpy random number generator for the species simulator. """
    return S_NPRNGvar

S_RNGvar = random.Random()
S_NPRNGvar = np.random.default_rng()
def seed_S_generators(seed: int):
    """
    Seed the global random number generators for the species simulator.
    A seed of 0 means to not reseed.
    """
    global S_RNGvar, S_NPRNGvar
    if seed != 0 and seed is not None:
        S_RNGvar = random.Random(seed)
        S_NPRNGvar = np.random.default_rng(seed)
