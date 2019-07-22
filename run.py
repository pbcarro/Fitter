#!/usr/bin/env python

from pyfitter.pyfitter import AsymmetricMolecule

molecule = AsymmetricMolecule(**{"verbose": False})

results = molecule.simulate_batch(1000, False)

