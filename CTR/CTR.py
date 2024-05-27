import numpy as np
import json
import traceback
import pandas as pd

################################################################################################################

pi   = np.pi
inf  = np.inf
#
exp  = lambda x: np.exp(x)
#
vec  = lambda *args: np.array(args)
unit = lambda *args: vec(*args)/norm(vec(*args))
span = lambda eval, evec : np.tensordot(eval, evec, axes=0)
#
def norm(v):
    if len(v.shape)>1:
        return np.linalg.norm(v, axis=1)
    else: return np.linalg.norm(v)
def Crit(v):
    n = np.where(v==0, 1, v)
    print(n)
    v = np.where(v==0, 0, 1/n)
    v[(np.where(v==0)[0]+1)%3] *= -1
    return v
def Orthogonal(v):
    z = np.count_nonzero(v)
    if z==1:
        M = vec([0, 1, 0], [0, 0, 1], [1, 0, 0])
        return M@v, M@M@v
    if z==2: o = Crit(v)
    if z==3: o = Crit(vec(*v[0:2], 0))
    return (o, np.cross(v, o))

################################################################################################################

class Xray:
    CuKa1  = 1.5406
    CuKa2  = 1.544390
    def __init__(self, wavelength=CuKa1, Energy=None, TTHETA=np.rad2deg(np.linspace(0, pi, 1024)[1:])):
        (h, c, e) = (6.62607015E-34, 299792458, 1.6021773349E-19)
        if Energy==None:
            self.wavelength = wavelength # Å unit
            self.Energy = h * c / (wavelength * 1E-10) / e
        else:
            self.Energy = Energy
            self.wavelength = h * c / (Energy * 1E-10) / e
        self.k = 2 * pi / self.wavelength
        self.TTHETA = TTHETA
        self.G = 2 * self.k * np.sin(np.deg2rad(self.TTHETA)/2)
    def Q(self, *ref):
        self.q = span(self.G, unit(*ref))
        return self.q
    def pseudoQ(self, *ref):
        ref = vec([1, 1, 0], [-1, 1, 0], [0, 0, 1]) @ vec(*ref) / 2
        return self.Q(*ref)
    PQ = pseudoQ # Alias
    def HKL(self, molecule): return self.q * molecule.abc / 2 / pi
    def I(self, sample): return sample.I(self.q, self.Energy)
    def __call__(self, sample): return self.I(sample)


class Xray2d(Xray):
    def Q(self, *ref):
        ux, uy = Orthogonal(unit(*ref))
        self.X, self.Y = np.meshgrid(self.G,self.G)
        self.q = (span(self.X, ux)+span(self.Y, uy)).reshape(self.X.size, 3)
        return self.q
    def I(self, sample): return sample.I(self.q, self.Energy).reshape(self.X.shape)
    
################################################################################################################

class Atom():
    # File path
    # https://lampz.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
    PATH = "./DATA/sf"
    AFF = pd.read_csv(f"{PATH}/AFF0.csv", delimiter='\s+').applymap(lambda x: x.strip() if isinstance(x, str) else x)
    #
    def __init__(self, Z, def_name=None):
        self.Z = Z
        self.__den__ = 1
        self.f = None
        #
        if def_name == None:
            (_, _, _, text) = traceback.extract_stack()[-2]
            self.def_name = text[:text.find('=')].strip()
    def __truediv__(self, num):
        atom = Atom(self.Z)
        atom.def_name = self.def_name
        atom.__den__ = num
        return atom
    def __call__(self, *args):
        return (self, *args)
    #
    def aff(self, Q, E=Xray().Energy):
        if self.f == None:
            f = pd.read_csv(Atom.PATH+f"/{self.def_name.lower()}.nff", delimiter='\s+')	
            differences = abs(f.iloc[:, 0] - E)
            closest_row = f.loc[[differences.idxmin()]]
            self.f = vec(closest_row.f1 + 1j * closest_row.f2)[0]
            self.COEF = np.float_(Atom.AFF[Atom.AFF['Element'] == self.def_name].iloc[0].values[1:])
        (a1, b1, a2, b2, a3, b3, a4, b4, c) = self.COEF
        f0 = self.f + sum(c + vec(*[a * exp(-1 * b * np.power(norm(Q) / (4 * pi), 2)) for a, b in zip((a1, a2, a3, a4), (b1, b2, b3, b4))]))
        return f0 / self.__den__
#
class Molecule():
    def __init__(self, abc, structure):
        self.abc = vec(*abc)
        self.structure = vec(*structure)
        self.atoms = self.structure[:,0]
        self.RJ = self.structure[:,1:] * self.abc
    def __truediv__(self, substrate):
        substrate, ref = substrate
        return Molecule(
            abc = substrate.const_volume_abc(self.abc, ref),
            structure = self.structure
    )
    def __add__(self, molecule):
        return Molecule(max(self.abc, molecule.abc), [*self.structure, *molecule.structure])
    def __call__(self, *args): return self, args
    #
    def pseudocubic(*abc):
        a, b, c = abc
        ac      = np.sqrt(a**2 + b**2) / 2
        return vec(ac, ac, c/2)
    def const_volume_abc(self, film_abc, ref):
        ab = self.abc * (vec(1,1,1)-unit(*ref))
        return ab + ( np.prod(film_abc) / np.prod(np.where(ab==0, 1, ab)) ) * unit(*ref)
    #
    def Q2HKL(self, Q):
        return Q * self.abc / 2 / pi
    def SF(self, Q, E=Xray().Energy):
        aff = vec(*[atom.aff(Q, E) for atom in self.atoms])
        phase = (1j * self.RJ @ Q.T).astype(dtype=complex)
        return sum(aff * exp(phase))

class SC(Molecule):
    def __init__(self, abc, X):
        structure = [
            X(0, 0, 0)
        ]
        super().__init__(abc, structure)
    #
class BCC(Molecule):
    def __init__(self, abc, AB):
        A, B = AB
        structure = [
            A(0, 0, 0),
            B(0.5, 0.5, 0.5)
        ]
        super().__init__(abc, structure)
    #
class FCC(Molecule):
    def __init__(self, abc, X):
        structure = [
            X(0, 0, 0),
            *[X(*rj) for rj in (np.ones([3, 3]) - np.eye(3))/2]
        ]
        super().__init__(abc, structure)
    #
class Perovskite(Molecule):
    def __init__(self, abc, ABO):
        A, B, O = ABO
        structure = [
            A(0.5, 0.5, 0.5), 	#BCC
            B(0,0,0),
            *[O(*rj) for rj in (np.ones([3, 3]) - np.eye(3))/2],  #FCC
            # *[(B/8)(n >> 2 & 1, n >> 1 & 1, n & 1) for n in range(8)],
            # *[(O/2)(*rj) for rj in (np.ones([6,3])+vec(*np.eye(3), *-np.eye(3)))/2 ]
        ]
        super().__init__(abc, structure)
#
class Film():
	def __init__(self, molecule, N):
		# self.molecule, self.ref = molecule
		self.molecule = molecule
		self.N = vec(*N)
	def __call__(self, *args): return Sample(self, nref=args)
	def __truediv__(self, substrate):
		# Film/Sample
		if 'Sample' in str(substrate.__class__):
			return Sample(self, *substrate.FILMS, nref=substrate.nref)
		else:   # Film/Film
			return Sample(self, substrate, nref=None)
	#
	def SN(self, Q): # X = Q @ molecule.abc
		IX = 1j * Q * self.molecule.abc
		Noinf = np.where(self.N == inf, 0, self.N)
		NUM = np.where(exp(IX)==1, Noinf, 1-exp(IX*Noinf))
		NUM = np.where(self.N == inf, -1, NUM)
		DEN = np.where(exp(IX)==1, 1, 1-exp(IX))
		return np.prod(NUM/DEN, axis=1)
	#
	def F(self, Q, E=Xray().Energy): return self.molecule.SF(Q, E) * self.SN(Q)
	def I(self, Q, E=Xray().Energy): return np.abs(self.F(Q, E)) ** 2

class RoughCut(Film):
	def __init__(self, molecule, N, nref, beta):
		super().__init__(molecule, N)
		self.nref = vec(*nref)
		self.beta = beta
	def SN(self, Q): # X = Q @ molecule.abc
		IX = 1j * Q * self.molecule.abc
		expIX = exp(IX)
		Noinf = np.where(self.N == inf, 0, self.N)
		NUM = np.where(expIX==1, Noinf, 1-exp(IX*Noinf))
		NUM = np.where(self.N == inf, -1, NUM)
		DEN = np.where(expIX==1, 1, 1-expIX)
		# Roughness
		NUMR = exp(Noinf * self.nref) * np.where(self.beta * expIX==1, -1, self.beta * expIX)
		DENR = np.where(self.beta * expIX==1, 1, 1-self.beta*expIX)
		return np.prod(NUM/DEN + NUMR/DENR, axis=1)
#
class Sample():
    def __init__(self, *films, nref):
        self.FILMS = films
        self.substrate = films[-1]
        self.film = films[0:-1]
        self.nref = nref
    def __truediv__(self, substrate):
        # Sample/Sample
        if 'Sample' in str(substrate.__class__):
            return Sample(*self.FILMS, *substrate.FILMS, nref=substrate.nref)
        else:	# Sample/Film
            return Sample(*self.FILMS, substrate, nref=None)
    #
    def F(self, Q, E=Xray().Energy):
        F = self.substrate.F(Q, E)
        PHI = np.zeros_like(Q).astype(dtype=complex)
        for film in self.film:
            IX  = 1j * Q * film.molecule.abc
            expNIX = np.prod(exp(PHI), axis=1)
            F += (expNIX * film.F(Q, E))
            PHI += IX * (np.where(film.N==inf, 0, film.N) * self.nref)
        return vec(*F)
    #
    def I(self, Q, E=Xray().Energy): return np.abs(self.F(Q, E)) ** 2
