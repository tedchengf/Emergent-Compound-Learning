import numpy as np
import random
import copy
import warnings
from itertools import product, combinations, permutations, chain
warnings.filterwarnings("ignore")

class Compound_Model():
	def __init__(self, elements, products, emergent_compound_products = None):
		self._elements = None
		self._compounds = None
		self._E_compounds = None
		self._P_compounds = None
		self._E_compounds_products = None 
		self._compound_to_ind = None
		self._ind_to_compound = None
		self._subset_mat = None
		self._superset_mat = None
		
		self._products = None
		self._products_id = None
		self._prod_to_id = None
		self._id_to_prod = None

		self._initialize(elements, products, emergent_compound_products)

	def __getattr__(self, name):
		if name == "elements": return self._elements.copy()
		if name == "compounds": return self._compounds.copy()
		if name == "E_compounds": return self._E_compounds.copy()
		if name == "P_compounds": return self._P_compounds.copy()
		if name == "E_compounds_products": return self._E_compounds_products.copy()
		if name == "compound_to_ind": return self._compound_to_ind.copy()
		if name == "ind_to_compound": return self._ind_to_compound.copy()
		if name == "subset_mat": return self._subset_mat.copy()
		if name == "superset_mat": return self._superset_mat.copy()
		if name == "compound_products": return self.get_products()

		if name == "products": return self._products.copy()
		if name == "products_id": return self._products_id.copy()
		if name == "prod_to_id": return self._prod_to_id.copy()
		if name == "id_to_prod": return self._prod_to_id.copy()
		raise AttributeError

	# Compounds must be in order (for now)
	def compound_to_ind(self, query):
		return self._query_parser(query, "i")

	def ind_to_compound(self, query):
		return self._query_parser(query, "c")

	# query_A is subset of query_B
	def is_subset(self, query_A, query_B):
		row_ind = self._query_parser(query_A, "i")
		col_ind = self._query_parser(query_B, "i")
		return self._subset_mat[np.ix_(row_ind, col_ind)]

	# query_A is subset of query_B
	def is_superset(self, query_A, query_B):
		row_ind = self._query_parser(query_A, "i")
		col_ind = self._query_parser(query_B, "i")
		return self._superset_mat[np.ix_(row_ind, col_ind)]	

	def find_common_subsets(self, query, target = None, return_type = "c"):
		row_inds = self._query_parser(query, "i")
		if target is None:
			target_compound = self.compounds
			submat = self._superset_mat[np.ix_(row_inds, np.arange(len(self._compounds), dtype = int))]
		else:
			target_compound = self._query_parser(target, "i")
			submat = self._superset_mat[np.ix_(row_inds, target_compound)]

	# find subsets of query (can constrain the search on the target compounds)
	def find_subsets(self, query, target = None, common_subset = False, return_type = "c"):
		row_inds = self._query_parser(query, "i")
		if target is None:
			target_compound = self.compounds
			submat = self._superset_mat[np.ix_(row_inds, np.arange(len(self._compounds), dtype = int))]
		else:
			target_compound = self._query_parser(target, "i")
			submat = self._superset_mat[np.ix_(row_inds, target_compound)]
		if common_subset == False:
			subsets = []
			for r in submat:
				subsets.append(target_compound[r])
			rsp = []
			for subset in subsets:
				rsp.append(self._query_parser(subset, return_type))
			return nparray_convert(rsp)
		else:
			c_sub = np.ones(submat.shape[1], dtype = bool)
			for r in submat: c_sub = np.multiply(c_sub, r)
			return self._query_parser(target_compound[c_sub], return_type)

	# find supersets of query (can constrain the search on the target 
	# compounds)
	def find_supersets(self, query, target = None, common_superset = False, return_type = "c"):
		row_inds = self._query_parser(query, "i")
		if target is None:
			target_compound = self.compounds
			supermat = self._subset_mat[np.ix_(row_inds, np.arange(len(self._compounds), dtype = int))]
		else:
			target_compound = self._query_parser(target, "i")
			supermat = self._subset_mat[np.ix_(row_inds, target_compound)]
		if common_superset == False:
			supersets = []
			for r in supermat:
				supersets.append(target_compound[r])
			rsp = []
			for superset in supersets:
				rsp.append(self._query_parser(superset, return_type))
			return nparray_convert(rsp)
		else:
			c_sup = np.ones(supermat.shape[1], dtype = bool)
			for r in supermat: c_sup = np.multiply(c_sup, r)
			return self._query_parser(target_compound[c_sup], return_type)
	
	def is_subset(self, query_A, query_B):
		row_ind = self._query_parser(query_A, "i")
		col_ind = self._query_parser(query_B, "i")
		return self._subset_mat[np.ix_(row_ind, col_ind)]

	def initialize_products(self, emergent_compound_products):
		self._E_compounds = set(())
		self._P_compounds = set(())
		self._E_compounds_products = dict({})
		self.verify_compounds(emergent_compound_products)
		self.verify_products(emergent_compound_products.values())
		for c in emergent_compound_products:
			p = emergent_compound_products[c]
			self._E_compounds.add(c)
			self._E_compounds_products.update({c: self._prod_to_id[p]})
		for c in self._compounds:
			if c not in self._E_compounds: self._P_compounds.add(c)

	def get_products(self, query, target = "i"):
		if target not in ("p", "i"): raise KeyError("target must be in ('p', 'i')")
		products = []
		for c in self._query_parser(query, 'c'):
			c_partitions = self.find_subsets([c], target = self._E_compounds, return_type = "c")[0]
			if len(c_partitions) == 0: 
				warnings.warn("No valid decomposition found for compound " + str(c) + "; it's reward is assumed to be None.")
				c_prod = [0]
			else:
				c_prod = [self._E_compounds_products[par] for par in c_partitions]
			if target == "i":
				products.append(c_prod)
			else:
				products.append([self._id_to_prod[idx] for idx in c_prod])
		return products

	def verify_compounds(self, compound_list):
		for c in compound_list: 
			if c not in self._compound_to_ind: raise KeyError(
			"Undefined compound: " + str(c))

	def verify_products(self, product_list):
		for p in product_list:
			if p not in self._prod_to_id: 
				print(p)
				raise KeyError("Undefined product: ", + str(p))

	def _initialize(self, elements, products, emergent_compound_products):
		self._elements = np.array(elements, dtype = object)
		self._compounds = nparray_convert(list(powerset(self.elements)))
		self._subset_mat = subset_mat(self.compounds)
		self._superset_mat = superset_mat(self.compounds)
		self._compound_to_ind = dict({})
		self._ind_to_compound = dict({})
		for ind, compound in enumerate(self.compounds):
			self._compound_to_ind.update({compound: ind})
			self._ind_to_compound.update({ind: compound})

		prime_generator = get_primes()
		self._products = [None] + products.copy()
		self._products_id = [0]
		self._prod_to_id = dict({None:0})
		self._id_to_prod = dict({0:None})
		self._id_to_ind = dict({0:0})
		for ind, prod in enumerate(products):
			if prod is None: pass
			prod_id = next(prime_generator)
			self._products_id.append(prod_id)
			self._prod_to_id.update({prod: prod_id})
			self._id_to_prod.update({prod_id: prod})
			self._id_to_ind.update({prod_id: ind + 1})
		self._products = nparray_convert(self._products)
		self._products_id = np.array(self._products_id, dtype = int)
		
		if emergent_compound_products is not None:
			self.initialize_products(emergent_compound_products)
		return 

	def _query_parser(self, query, target):
		if type(query) is tuple: 
			if target == "i": return np.array([self._compound_to_ind[query]], dtype = int)
			else: return nparray_convert([query])
		if type(query) is int or type(query) is np.int64:
			if target == "c": return nparray_convert([self._ind_to_compound[query]])
			else: return np.array([query], dtype = int)
		try:
			iter_test = iter(query)
			if type(query) in (dict, str): raise TypeError("Unacceptable query type: " + str(type(query)))
			if len(query) == 0: return np.empty((0), dtype = object)
			if target == "i":
				return np.array(list(map(lambda x: self._query_parser(x, target)[0], query)), dtype = int)
			else:
				return nparray_convert(list(map(lambda x: self._query_parser(x, target)[0], query)))
		except(TypeError) as te: pass	
		raise TypeError("Unacceptable query type: " + str(type(query)))

from tqdm import tqdm
from scipy.special import logsumexp
class Partition_Space():
	def __init__(self, compound_model, constrain_compounds = None):
		self.hypotheses = None
		self.hypotheses_incidence = None
		self.prior = None
		self._initialize(compound_model, constrain_compounds)

	def set_prior(self, prior_func):
		self.prior = prior_func(self.hypotheses[0].compounds, self.hypotheses_incidence)
		return

	def _initialize(self, compound_model, constrain_compounds):
		compound_ids = np.arange(len(compound_model._compounds), dtype = int)
		compound_len = len(compound_ids)
		if constrain_compounds is not None:
			constrain_compounds_id = compound_model.compound_to_ind	(constrain_compounds.keys())
			compound_ids = list(filter(lambda x: x not in constrain_compounds_id, compound_ids))
			compound_ids = np.array(compound_ids, dtype = int)
		hypotheses = np.zeros((2**(len(compound_ids)), compound_len), dtype = int)
		all_partitions = powerset(compound_ids)
		for i, subset in enumerate(all_partitions): 
			hypotheses[i][list(subset)] = True
		if constrain_compounds is not None:
			hypotheses[:, constrain_compounds_id] = True
		self.hypotheses_incidence = hypotheses
		self.hypotheses = []
		for hypothesis in hypotheses:
			self.hypotheses.append(self._initialize_bayes_model(hypothesis, compound_model, constrain_compounds))
		self.hypotheses = np.array(self.hypotheses, dtype = object)

	def _initialize_bayes_model(self,hypothesis,compound_model,constrain_compounds):
		emergent_compound_products = {}
		for index, compound in enumerate(compound_model._compounds):
			if hypothesis[index] == True:
				emergent_compound_products.update({compound_model._compounds[index]: None})
		if constrain_compounds is not None: emergent_compound_products.update(constrain_compounds)
		return Bayesian_Model(compound_model._elements, list(compound_model._products[1:]), emergent_compound_products)

class Bayesian_Model(Compound_Model):
	def __init__(self, elements, products, emergent_compound_products):
		super().__init__(elements, products, emergent_compound_products)
		self._uninitialized_E_compounds = set({})
		self._product_candidate_table = dict({})
		# update the empty rewards
		for c in emergent_compound_products:
			if emergent_compound_products[c] is None: self._uninitialized_E_compounds.add(c)
		

	def validate(self, compound, observed_products):
		model_pred = self.get_products(compound)[0]
		
		print(model_pred)

		# Depth 1
		prod_encoding = np.product(model_pred)
		print(prod_encoding)
		if prod_encoding == np.product(observed_products): return True
		# Depth 2
		if prod_encoding != 0: return False
		# Depth 3
		if len(model_pred) != len(observed_products): return False
		unknown_mask = model_pred > 0
		print(unknown_mask)
		return True

	def _valid_partition(self, product_1, product_2):
		return

	def add_emergent_compounds(self, emergent_compound_products):
		self.verify_compounds(emergent_compound_products)
		self.verify_products(emergent_compound_products.values())
		for c in emergent_compound_products:
			if c not in self._uninitialized_E_compounds: raise RuntimeError("Unwarrented update; the compound ", c, " is either updated in the past or does does not belong to the emergent compound set")
			self._uninitialized_E_compounds.remove(c)
			self._E_compounds_products.update({c: emergent_compound_products[c]})

# Helper Functions
###############################################################################

def uniform_prior(compounds, incidences):
	return np.ones(incidences.shape[0])*(-np.log(incidences.shape[0]))

# useful when you do not want numpy arrays to collapse iterables
def nparray_convert(arr):
	np_arr = np.empty(len(arr), dtype = object)
	for ind, obj in enumerate(arr):
		np_arr[ind] = obj
	return np_arr

def multi_union(arr):
	uniq_set = set(())
	for compounds in arr: uniq_set = uniq_set.union(compounds)
	return uniq_set

def subset_mat(arr):
	subset = lambda x, y: set(x).issubset(y)
	subset_relmat = np.ones((len(arr), len(arr)), dtype = bool)
	non_symmetrical_matrix_iteration(arr, subset_relmat, subset)
	return subset_relmat

def superset_mat(arr):
	superset = lambda x, y: set(x).issuperset(y)
	superset_relmat = np.ones((len(arr), len(arr)), dtype = bool)
	non_symmetrical_matrix_iteration(arr, superset_relmat, superset)
	return superset_relmat

def incidence_vectors(elements, arr):
	incidence_mat = np.zeros((len(elements), len(arr)), dtype = int)
	for ci,c in enumerate(arr):
		for ei, e in enumerate(elements):
			if e in c: incidence_mat[ei, ci] = 1
	return incidence_mat.T

def powerset(arr):
	s = list(arr)
	return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def non_symmetrical_matrix_iteration(data_array, target_matrix, function):
	mat_dim = target_matrix.shape[0]
	for r in range(mat_dim):
		for c in range(mat_dim):
			target_matrix[r,c] = function(data_array[r], data_array[c])

def symmetrical_matrix_iteration(data_array, target_matrix, function, skip_diagonal=True):
	mat_dim = target_matrix.shape[0]
	if skip_diagonal == True:
		for r in range(mat_dim):
			for c in range(r+1, mat_dim):
				target_matrix[r,c] = function(data_array[r], data_array[c])
				target_matrix[c,r] = target_matrix[r,c]
	else:
		for r in range(mat_dim):
			for c in range(r, mat_dim):
				target_matrix[r,c] = function(data_array[r], data_array[c])
				target_matrix[c,r] = target_matrix[r,c]

# A little generator that yield prime numbers
def get_primes():
	D = {}
	q = 2
	while True:
		if q not in D:
			yield q
			D[q * q] = [q]
		else:
			for p in D[q]:
				D.setdefault(p + q, []).append(p)
			del D[q]
		q += 1