
import os
import numpy as np
import pandas as pd
from pathlib import Path
from ctypes import c_uint, c_int, c_double, create_string_buffer, CDLL, POINTER, byref, Structure


###Structure definition for python
class Level(Structure):
    _fields_ = [
        ("Index", c_uint),
        ("J", c_uint),
        ("Ka", c_uint),
        ("Kc", c_uint),
        ("Energy",c_double)
        ]

class Transition(Structure):
    _fields_ = [
        ("Frequency", c_double),
        ("Upper", c_uint),
        ("Lower", c_uint),
        ("Type", c_uint),
        ("Intensity",c_double)
        ]

class ETauStruct(Structure):
    _fields_ = [
        ("StatePoints", c_uint),
        ("Delta", c_double),
        ("ETVals", POINTER(c_double))
        ]

class Triple(Structure):
    _fields_ = [
        ("TriplesCount", c_uint),
        ("TransitionList", Transition),
        ("TriplesList", POINTER(c_double))
        ]

class Opt_Bundle(Structure):
    #Note this is currently not identical to the C code, TransitionGSL is a pointer rather than a finite block
    _fields_ = [
        ("ETGSL", ETauStruct),
        ("MyDictionary", POINTER(Level)),
        ("TransitionsGSL", POINTER(Transition))
        ]


class AsymmetricMolecule:
    """
    Python interface to Brandon's fast rigid rotor simulator. This provides
    a way to control the program to a certain degree, with the main user
    input being what the rotational constants of the molecule are, and what
    degree of verbosity in the program output.

    The main methods a user will call are `simulate` and `format_result`:
    the former will call the routines to produce a catalog, and the latter
    will format the results of the catalog (which are ctypes) into a more
    user-friendly data type, either as NumPy arrays or Pandas Dataframes.
    """
    def __init__(self, **kwargs):
        self.seed = None
        self.niter = 100
        self.constants = [10000., 5000., 3000.]
        self.verbose = False
        # These are static paths that a user can overwrite
        module_path = Path(__file__).parent
        self.lib_path = module_path.joinpath("Fitter.so")
        self.et_path = module_path.joinpath("etau.dat")
        self.cat_dict_path = module_path.joinpath("base_cat_dict.txt")
        self.cat_path = module_path.joinpath("base_cat.txt")

        # Update the parameters with user defined settings
        self.__dict__.update(**kwargs)

        # Seed the random number generator, unless user provides a value
        np.random.seed(self.seed)

        # Set up a whole bunch of Ctypes stuff
        self._load_library()
        self._init_pointers()
        self._load_tables()

    def _load_tables(self):
        """
        Private method to load the precomputed tables into memory.
        """
        names = ["etau", "catdict", "catalog"]
        self.string_buffers = dict()
        for name, path in zip(names, [self.et_path, self.cat_dict_path, self.cat_path]):
            check_path = Path(path)
            if not check_path.exists():
                raise Exception(f"{path} table not found!")
            self.string_buffers[name] = create_string_buffer(bytes(check_path))
        # Load the base catalog dictionary
        self.FitterLib.Load_Base_Catalog_Dictionary(
            self.string_buffers["catdict"],
            byref(self.levels),
            self._verbose
        )
        # Load the base catalog
        self._statecount = self.FitterLib.Load_Base_Catalog(
            self.string_buffers["catalog"],
            byref(self.catalog),
            self._verbose
        )
        # Load Etau table
        self.FitterLib.Load_ETau_File2(
            self.string_buffers["etau"],
            byref(self.et),
            byref(self._etstatecount)
        )

    def _load_library(self):
        """
        Private method to use ctypes to load in Brandon's program as a static
        library to access its functions. The way that this function is written
        is more verbose, but it is intended to give useful errors when the
        library is not found/missing.
        """
        if not self.lib_path.exists():
            raise Exception("Fitter static library not found; is it compiled?")
        else:
            # Load the statically compiled library
            self.FitterLib = CDLL(self.lib_path)

    def _init_pointers(self):
        """
        Private method to set up the ctypes and pointers prior to spinning up
        the simulations.
        """
        #self._etArray = POINTER(c_double)()
        #self._filedelta = c_double(0)
        self._statepoints = c_int(0)
        self._statecount = c_int(0)
        #self._etstatecount = c_int(0)
        self._verbose = c_int(int(self.verbose))
        self._constantsType = c_double * 3
        self.levels = POINTER(Level)()
        self.catalog = POINTER(Transition)()
        self.et = ETauStruct()

    def simulate(self, constants=None):
        """
        Function to run the program to generate a catalog based on a set of
        rigid rotor constants.

        Parameters
        ----------
        constants : float or None, optional
            3-element iterable consisting of floats, corresponds to the rotational
            constants in MHz.
        """
        if constants is None:
            constants = self.constants
        # Make sure exactly three rotational constants are given
        assert len(constants) == 3
        # Ensure that the values are sorted in the correct order
        constants = list(np.sort(constants))[::-1]
        if self.verbose:
            print(f"Simulating with A, B, C: {constants}")
            print("Debugging information:")
            print(f"StateCount: {self._statecount}")
        c_constants = self._constantsType(*constants)
        # Run the catalog simulation
        self.FitterLib.Get_Catalog(
            self.catalog,
            c_constants,
            self._statecount,
            self._verbose,
            self.et,
            self.levels
        )

    def simulate_batch(self, niter=None, pandas_dataframe=True, **kwargs):
        """
        Function to run a batch of simulations. The constants are randomly generated
        between runs, although the random seed is only changed when this class is
        created so be aware when making comparisons.

        Kwargs are passed into the constant generation, and so a user can provide
        upper and lower boundary values for A, B, C.

        Parameters
        ----------
        niter : int or None, optional
            If a value is provided by the user, sets the number of simulations to run.
        pandas_dataframe : bool, optional
            If True, the results that are returned are organized as Pandas dataframes.
        kwargs
            Kwargs are passed into the random constants generation.

        Returns
        -------
        results : list
            Nested list of two-tuple, corresponding to the constants used
        """
        if not niter:
            niter = self.niter
        results = list()
        for _ in range(niter):
            constants = self.random_constants(**kwargs)
            self.simulate(constants)
            results.append((constants, self.format_results(pandas_dataframe)))
        return results

    def random_constants(self, lower=1000., upper=10000.):
        """
        Function to generate three rotational constants with upper and lower
        boundary values, based on a uniform distribution. The ordering
        of the constants returned is A, B, C

        Parameters
        ----------
        lower, upper : float, optional
            Lower and upper limits to the random number generation, in MHz.

        Returns
        -------
        constants : NumPy 1D array
        """
        constants = np.random.uniform(lower, upper, 3)
        constants = np.sort(constants)[::-1]
        return constants

    def format_results(self, pandas_table=True):
        """
        Function to convert the C-types catalog into Python data types,
        either as NumPy arrays or as a Pandas dataframe. If the NumPy array
        is returned, the headings go as:

        frequency, upper-index, lower-index, type
        J', Ka', Kc', J'', Ka'', Kc''

        Parameters
        ----------
        pandas_table : bool, optional
            If True (default), the output returned is a Pandas dataframe.

        Returns
        -------
        If pandas_table is True, a pandas dataframe, otherwise a NumPy array.
        """
        table = list()
        # This loop is explicitly done non-Pythonically; iterating over a
        # Pointer will continue until it breaks itself in a segfault, and
        # the attribute `_statecount` holds the correct number of transitions.
        for index in range(self._statecount):
            transition = self.catalog[index]
            ustate = self.levels[transition.Upper]
            lstate = self.levels[transition.Lower]
            data = [
                float(transition.Frequency),
                int(transition.Upper),
                int(transition.Lower),
                int(transition.Type),
                int(ustate.J),
                int(ustate.Ka),
                int(ustate.Kc),
                int(lstate.J),
                int(lstate.Ka),
                int(lstate.Kc)
            ]
            table.append(data)
        if pandas_table:
            return pd.DataFrame(
                table, 
                columns=[
                    "frequency", "upper-index", "lower-index", "type", 
                    "J'", "Ka'", "Kc'", "J''", "Ka''", "Kc''"
                    ]
                )
        else:
            return np.array(table)
