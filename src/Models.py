import numpy as np
import random
import copy
import warnings
from itertools import product, combinations, permutations, chain

class Compound_Model():
	def __init__(self, elements, emergent_compound_rewards = None):
		self.__elements = None
		self.__compounds = None
		self.__E_compounds = None
		self.__P_compounds = None
		#####
		self.__compound_rewards = None
		self.__compound_to_ind = None
		self.__ind_to_compound = None
		self.__compounds_iv = None
		self.__subset_mat = None
		self.__superset_mat = None

		self.__initialize(elements, emergent_compound_rewards)

	# existing emergent compounds' value will be overwritten
	def add_emergent_compounds(self, emergent_compound_rewards):
		self.__verify_compounds(emergent_compound_rewards)
		for c in emergent_compound_rewards:
			self.__E_compounds.add(c)
			self.__compound_rewards.update({c: emergent_compound_rewards[c]})
		affected_compounds = self.find_supersets(emergent_compound_rewards.keys(), return_type = "c")
		# there might be more than one list of compounds; here we convert
		# them into one set
		affected_compounds = multi_union(affected_compounds)
		self.__compound_rewards.update(self.__update_compound_rewards(affected_compounds))
		return

	def del_emergent_compounds(self, emergent_compounds):
		self.__verify_compounds(emergent_compounds)
		for c in emergent_compounds:
			self.__E_compounds.remove(c)
		affected_compounds = self.find_supersets(emergent_compounds, return_type = "c")
		affected_compounds = multi_union(affected_compounds)
		self.__compound_rewards.update(self.__update_compound_rewards(affected_compounds))
		return

	def initialize_rewards(self, emergent_compound_rewards):
		self.__E_compounds = set(())
		self.__P_compounds = set(())
		self.__compound_rewards = {}
		self.__verify_compounds(emergent_compound_rewards)
		for c in emergent_compound_rewards:
			self.__E_compounds.add(c)
			self.__compound_rewards.update({c: emergent_compound_rewards[c]})
		p_compounds = self.__update_compound_rewards(self.__compounds)
		for c in p_compounds: self.__P_compounds.add(c)
		self.__compound_rewards.update(p_compounds)
		return

	def get_rewards(self, query):
		return [self.__compound_rewards[c] for c in self.__query_parser(query, "c")]

	# query_A is subset of query_B
	def is_subset(self, query_A, query_B):
		row_ind = self.__query_parser(query_A, "i")
		col_ind = self.__query_parser(query_B, "i")
		return self.__subset_mat[np.ix_(row_ind, col_ind)]

	# query_A is subset of query_B
	def is_superset(self, query_A, query_B):
		row_ind = self.__query_parser(query_A, "i")
		col_ind = self.__query_parser(query_B, "i")
		return self.__superset_mat[np.ix_(row_ind, col_ind)]	

	# find subsets of query (can constrain the search on the target compounds)
	def find_subsets(self, query, target = None, return_type = "c"):
		row_inds = self.__query_parser(query, "i")
		if target is None:
			target_compound = self.compounds
			submat = self.__superset_mat[np.ix_(row_inds, np.arange(len(self.__compounds), dtype = int))]
		else:
			target_compound = self.__query_parser(target, "i")
			submat = self.__superset_mat[np.ix_(row_inds, target_compound)]
		subsets = []
		for r in submat:
			subsets.append(target_compound[r])
		rsp = []
		for subset in subsets:
			rsp.append(self.__query_parser(subset, return_type))
		return nparray_convert(rsp)

	# find supersets of query (can constrain the search on the target 
	# compounds)
	def find_supersets(self, query, target = None, return_type = "c"):
		row_inds = self.__query_parser(query, "i")
		if target is None:
			target_compound = self.compounds
			submat = self.__subset_mat[np.ix_(row_inds, np.arange(len(self.__compounds), dtype = int))]
		else:
			target_compound = self.__query_parser(target, "i")
			submat = self.__subset_mat[np.ix_(row_inds, target_compound)]
		subsets = []
		for r in submat:
			subsets.append(target_compound[r])
		rsp = []
		for subset in subsets:
			rsp.append(self.__query_parser(subset, return_type))
		return nparray_convert(rsp)

	# Compounds must be in order (for now)
	def compound_to_ind(self, query):
		return self.__query_parser(query, "i")

	def ind_to_compound(self, query):
		return self.__query_parser(query, "c")

	def get_incidence_vectors(self, query):
		return self.__compounds_iv[self.__query_parser(query, "i"), :].copy()

	def __getattr__(self, name):
		if name == "elements": return self.__elements.copy()
		if name == "compounds": return self.__compounds.copy()
		if name == "E_compounds": return self.__E_compounds.copy()
		if name == "P_compounds": return self.__P_compounds.copy()
		if name == "subset_mat": return self.__subset_mat.copy()
		if name == "superset_mat": return self.__superset_mat.copy()
		if name == "compound_rewards": return self.__compound_rewards.copy()
		raise AttributeError

	def __initialize(self, elements, emergent_compound_rewards):
		self.__elements = np.array(elements, dtype = object)
		self.__compounds = nparray_convert(list(powerset(self.elements)))
		self.__compound_to_ind = {}
		self.__ind_to_compound = {}
		for ind, compound in enumerate(self.compounds):
			self.__compound_to_ind.update({compound: ind})
			self.__ind_to_compound.update({ind: compound})
		self.__compounds_iv = incidence_vectors(self.elements, self.compounds)
		self.__subset_mat = subset_mat(self.compounds)
		self.__superset_mat = superset_mat(self.compounds)
		
		if emergent_compound_rewards is not None:
			self.initialize_rewards(emergent_compound_rewards)
		return 

	def __verify_compounds(self, compound_list):
		for c in compound_list: 
			if c not in self.__compound_to_ind: raise ValueError(
			"Unknown compound: " + str(c))

	def __update_compound_rewards(self, compounds):
		new_compound_rewards = {}
		for c in compounds:
			if c not in self.__E_compounds:
				c_partitions = self.find_subsets([c], target = self.__E_compounds, return_type = "c")[0]
				if len(c_partitions) == 0: 
					warnings.warn("No valid decomposition found for compound " + str(c) + "; it's reward is assumed to be 0.")
					c_val = 0
				else:
					c_val = sum([self.__compound_rewards[par] for par in c_partitions])
				new_compound_rewards.update({c: c_val})
		return new_compound_rewards

	def __query_parser(self, query, target):
		if type(query) is tuple: 
			if target == "i": return np.array([self.__compound_to_ind[query]], dtype = int)
			else: return nparray_convert([query])
		if type(query) is int or type(query) is np.int64:
			if target == "c": return nparray_convert([self.__ind_to_compound[query]])
			else: return np.array([query], dtype = int)
		try:
			iter_test = iter(query)
			if type(query) in (dict, str): raise TypeError("Unacceptable query type: " + str(type(query)))
			if len(query) == 0: return np.empty((0), dtype = object)
			if target == "i":
				return np.array(list(map(lambda x: self.__query_parser(x, target)[0], query)), dtype = int)
			else:
				return nparray_convert(list(map(lambda x: self.__query_parser(x, target)[0], query)))
		except(TypeError) as te: pass	
		# if type(query) in (list, set, np.ndarray):
		# 	if len(query) == 0: return np.empty((0), dtype = object)
		# 	if target == "i":
		# 		return np.array(list(map(lambda x: self.__query_parser(x, target)[0], query)), dtype = int)
		# 	else:
		# 		return nparray_convert(list(map(lambda x: self.__query_parser(x, target)[0], query)))
		raise TypeError("Unacceptable query type: " + str(type(query)))

# Helper Functions
###############################################################################

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

def subset_mat(compounds):
	subset = lambda x, y: set(x).issubset(y)
	subset_relmat = np.ones((len(compounds), len(compounds)), dtype = bool)
	non_symmetrical_matrix_iteration(compounds, subset_relmat, subset)
	return subset_relmat

def superset_mat(compounds):
	superset = lambda x, y: set(x).issuperset(y)
	superset_relmat = np.ones((len(compounds), len(compounds)), dtype = bool)
	non_symmetrical_matrix_iteration(compounds, superset_relmat, superset)
	return superset_relmat

def incidence_vectors(elements, compounds):
	incidence_mat = np.zeros((len(elements), len(compounds)), dtype = int)
	for ci,c in enumerate(compounds):
		for ei, e in enumerate(elements):
			if e in c: incidence_mat[ei, ci] = 1
	return incidence_mat.T

def powerset(iterable):
    s = list(iterable)
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