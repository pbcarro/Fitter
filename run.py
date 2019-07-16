#!/usr/bin/env python

import main
import line_profiler

molecule = main.AsymmetricMolecule(**{"verbose": False})

results = molecule.simulate_batch(1000, False)

