import numpy as np
import random
import copy
import warnings
from itertools import product, combinations, permutations, chain
warnings.filterwarnings("ignore")

class Compound_Model():
	def __init__(self, elements, emergent_compound_rewards = None):
		self._elements = None
		self._compounds = None
		self._E_compounds = None
		self._P_compounds = None
		self._compound_rewards = None
		self._compound_to_ind = None
		self._ind_to_compound = None
		self._compounds_iv = None
		self._subset_mat = None
		self._superset_mat = None

		self._initialize(elements, emergent_compound_rewards)

	# existing emergent compounds' value will be overwritten
	def add_emergent_compounds(self, emergent_compound_rewards):
		self._verify_compounds(emergent_compound_rewards)
		for c in emergent_compound_rewards:
			self._E_compounds.add(c)
			self._compound_rewards.update({c: emergent_compound_rewards[c]})
		affected_compounds = self.find_supersets(emergent_compound_rewards.keys(), return_type = "c")
		# there might be more than one list of compounds; here we convert
		# them into one set
		affected_compounds = multi_union(affected_compounds)
		self._compound_rewards.update(self._update_compound_rewards(affected_compounds))
		return

	def del_emergent_compounds(self, emergent_compounds):
		self._verify_compounds(emergent_compounds)
		for c in emergent_compounds:
			self._E_compounds.remove(c)
		affected_compounds = self.find_supersets(emergent_compounds, return_type = "c")
		affected_compounds = multi_union(affected_compounds)
		self._compound_rewards.update(self._update_compound_rewards(affected_compounds))
		return

	def initialize_rewards(self, emergent_compound_rewards):
		self._E_compounds = set(())
		self._P_compounds = set(())
		self._compound_rewards = {}
		self._verify_compounds(emergent_compound_rewards)
		for c in emergent_compound_rewards:
			self._E_compounds.add(c)
			self._compound_rewards.update({c: emergent_compound_rewards[c]})
		p_compounds = self._update_compound_rewards(self._compounds)
		for c in p_compounds: self._P_compounds.add(c)
		self._compound_rewards.update(p_compounds)
		return

	def get_rewards(self, query):
		return [self._compound_rewards[c] for c in self._query_parser(query, "c")]

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

	# find subsets of query (can constrain the search on the target compounds)
	def find_subsets(self, query, target = None, return_type = "c"):
		row_inds = self._query_parser(query, "i")
		if target is None:
			target_compound = self.compounds
			submat = self._superset_mat[np.ix_(row_inds, np.arange(len(self._compounds), dtype = int))]
		else:
			target_compound = self._query_parser(target, "i")
			submat = self._superset_mat[np.ix_(row_inds, target_compound)]
		subsets = []
		for r in submat:
			subsets.append(target_compound[r])
		rsp = []
		for subset in subsets:
			rsp.append(self._query_parser(subset, return_type))
		return nparray_convert(rsp)

	# find supersets of query (can constrain the search on the target 
	# compounds)
	def find_supersets(self, query, target = None, return_type = "c"):
		row_inds = self._query_parser(query, "i")
		if target is None:
			target_compound = self.compounds
			submat = self._subset_mat[np.ix_(row_inds, np.arange(len(self._compounds), dtype = int))]
		else:
			target_compound = self._query_parser(target, "i")
			submat = self._subset_mat[np.ix_(row_inds, target_compound)]
		subsets = []
		for r in submat:
			subsets.append(target_compound[r])
		rsp = []
		for subset in subsets:
			rsp.append(self._query_parser(subset, return_type))
		return nparray_convert(rsp)

	# Compounds must be in order (for now)
	def compound_to_ind(self, query):
		return self._query_parser(query, "i")

	def ind_to_compound(self, query):
		return self._query_parser(query, "c")

	def get_incidence_vectors(self, query):
		return self._compounds_iv[self._query_parser(query, "i"), :].copy()

	def __getattr__(self, name):
		if name == "elements": return self._elements.copy()
		if name == "compounds": return self._compounds.copy()
		if name == "E_compounds": return self._E_compounds.copy()
		if name == "P_compounds": return self._P_compounds.copy()
		if name == "subset_mat": return self._subset_mat.copy()
		if name == "superset_mat": return self._superset_mat.copy()
		if name == "compound_rewards": return self._compound_rewards.copy()
		raise AttributeError

	def _initialize(self, elements, emergent_compound_rewards):
		self._elements = np.array(elements, dtype = object)
		self._compounds = nparray_convert(list(powerset(self.elements)))
		self._compound_to_ind = {}
		self._ind_to_compound = {}
		for ind, compound in enumerate(self.compounds):
			self._compound_to_ind.update({compound: ind})
			self._ind_to_compound.update({ind: compound})
		self._compounds_iv = incidence_vectors(self.elements, self.compounds)
		self._subset_mat = subset_mat(self.compounds)
		self._superset_mat = superset_mat(self.compounds)
		
		if emergent_compound_rewards is not None:
			self.initialize_rewards(emergent_compound_rewards)
		return 

	def _verify_compounds(self, compound_list):
		for c in compound_list: 
			if c not in self._compound_to_ind: raise ValueError(
			"Unknown compound: " + str(c))

	def _update_compound_rewards(self, compounds):
		new_compound_rewards = {}
		for c in compounds:
			if c not in self._E_compounds:
				c_partitions = self.find_subsets([c], target = self._E_compounds, return_type = "c")[0]
				if len(c_partitions) == 0: 
					warnings.warn("No valid decomposition found for compound " + str(c) + "; it's reward is assumed to be None.")
					c_val = None
				else:
					try:
						c_val = sum([self._compound_rewards[par] for par in c_partitions])
					except TypeError:
						warnings.warn("One or more emergent compound after decomposition of " + str(c) + " are unitialized; it's reward is assumed to be None.")
						c_val = None
				new_compound_rewards.update({c: c_val})
		return new_compound_rewards

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
		# if type(query) in (list, set, np.ndarray):
		# 	if len(query) == 0: return np.empty((0), dtype = object)
		# 	if target == "i":
		# 		return np.array(list(map(lambda x: self.__query_parser(x, target)[0], query)), dtype = int)
		# 	else:
		# 		return nparray_convert(list(map(lambda x: self.__query_parser(x, target)[0], query)))
		raise TypeError("Unacceptable query type: " + str(type(query)))

from tqdm import tqdm
from scipy.special import logsumexp
class Partition_Space():
	def __init__(self, compound_model, constrain_compounds = None):
		self.hypotheses = None
		self.prior = None
		self._initialize(compound_model, constrain_compounds)

	def _initialize(self, compound_model, constrain_compounds):
		compound_ids = np.arange(len(compound_model._compounds), dtype = int)
		compound_len = len(compound_ids)
		if constrain_compounds is not None:
			constrain_compounds_id = compound_model.compound_to_ind		(constrain_compounds.keys())
			compound_ids = list(filter(lambda x: x not in constrain_compounds_id, compound_ids))
			compound_ids = np.array(compound_ids, dtype = int)
		hypotheses = np.zeros((2**(len(compound_ids)), compound_len))
		all_partitions = powerset(compound_ids)
		for i, subset in enumerate(all_partitions): 
			hypotheses[i][list(subset)] = True
		if constrain_compounds is not None:
			hypotheses[:, constrain_compounds_id] = True
		self.hypotheses = hypotheses

	def initialize_bayes_model(self,hypothesis,compound_model,constrain_compounds):
		emergent_compound_rewards = {}
		for index, compound in enumerate(compound_model._compounds):
			if hypothesis[index] == True:
				emergent_compound_rewards.update({compound_model._compounds : None})
		if constrain_compounds is not None: emergent_compound_rewards.update(constrain_compounds)
		return Bayesian_Model(compound_model._elements, emergent_compound_rewards)


class Bayesian_Model(Compound_Model):
	def __init__(self, elements, emergent_compound_rewards):
		super().__init__(elements, emergent_compound_rewards)
		self._uninitialized_E_compounds = set({})
		# update the empty rewards
		for c in emergent_compound_rewards:
			if emergent_compound_rewards[c] is None: self._uninitialized_E_compounds.add(c)

	def add_emergent_compounds(self, emergent_compound_rewards):
		super()._verify_compounds(emergent_compound_rewards)
		for c in emergent_compound_rewards:
			if c not in self._uninitialized_E_compounds: raise RuntimeError("Unwarrented update; the compound ", c, " is either updated in the past or does does not belong to the emergent compound set")
			self._uninitialized_E_compounds.remove(c)
			self._compound_rewards.update({c: emergent_compound_rewards[c]})
		affected_compounds = super().find_supersets(emergent_compound_rewards.keys(), return_type = "c")
		affected_compounds = multi_union(affected_compounds)
		self._compound_rewards.update(self._update_compound_rewards(affected_compounds))

	def del_emergent_compounds(self, emergent_compounds):
		raise RuntimeError("The function del_emergent_compounds cannot be used for an Bayesian_Model instance.")


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