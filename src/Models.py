import numpy as np
import random
import sys, os
import copy
import warnings
from itertools import product, combinations, permutations, chain
from tqdm import tqdm
from scipy.stats import gamma
from scipy.special import logsumexp
from scipy.stats import entropy
from scipy.stats import rankdata
warnings.filterwarnings("ignore")

###############################################################################
# Compound_Model Class
###############################################################################

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
				c_prod = np.array([self._E_compounds_products[par] for par in c_partitions], dtype = int)
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

###############################################################################
# Bayesian_Model Class
###############################################################################

# The bottom level class for the models. The Bayesian_Model is a subclass of
# Compound_Model. It instantiate a bayesian hypothesis and check for consistency
# when updating new observations.
class Bayesian_Model(Compound_Model):
	def __init__(self, elements, products, emergent_compound_products):
		super().__init__(elements, products, emergent_compound_products)
		# A set of emergent compounds whose product have not been determined
		self._uninitialized_E_compounds = set({})
		# The product candidate table that keeps track of possible product
		# assignment to compounds. The product is in keys, and the candidate are
		# in values as lists
		self._product_candidate_table = dict({})
		# A flag that represent whether the model had ran into conflict. 
		self.solve_flag = True
		# update the empty rewards
		for c in emergent_compound_products:
			if emergent_compound_products[c] is None: self._uninitialized_E_compounds.add(c)

	# The main function that update the new observation into the product
	# candidate table. It checks for consistency and return False whenever the
	# prediction is incorrect or the merge operation fails.
	def validate(self, compound, observed_products, update):
		if self.solve_flag == False: return False
		model_pred = self.get_products(compound)[0]
		# Depth 1
		prod_encoding = np.product(model_pred)
		if prod_encoding == np.product(observed_products): 
			return True
		# Depth 2
		if prod_encoding != 0: 
			if update == True:
				self.solve_flag = False
			return False
		# Depth 3
		if len(model_pred) != len(observed_products): 
			if update == True:
				self.solve_flag = False
			return False
		model_pred_set = set(model_pred)
		observed_prod_set = set(observed_products)
		# Note that this holds because we have previously tested for len
		valid_par = model_pred_set - observed_prod_set == {0}
		if valid_par == False: 
			if update == True:
				self.solve_flag = False
			return False
		# Depth 4
		potential_candidates = self.find_subsets([compound], target = self._uninitialized_E_compounds) [0]
		products_dict = ({})
		for prod in observed_prod_set - model_pred_set: 
			products_dict.update({prod: potential_candidates})
		if update == True:
			product_candidate_table = self._product_candidate_table
		else:
			product_candidate_table = copy.deepcopy(self._product_candidate_table)
		res = self._product_merge(products_dict, product_candidate_table, set({}))
		if type(res) is bool:
			if update == True:
				self.solve_flag = False
			return False
		res = self._compound_merge(*res)
		if type(res) is bool: 
			if update == True:
				self.solve_flag = False
			return False
		if update == True:
			self._udpate_states(*res)
		return True

	def _product_merge(self, products_dict, product_candidate_table, changed_rows):
		for prod in products_dict:
			curr_cands = products_dict[prod]
			# simple updates
			if prod not in product_candidate_table:
				product_candidate_table.update({prod: set(curr_cands)})
				changed_rows.add(prod)
			else:
				recorded_cands = product_candidate_table[prod]
				intersection = recorded_cands.intersection(curr_cands)
				# Conflict
				if len(intersection) == 0: return False
				# No update
				elif intersection == recorded_cands: pass
				else:
					product_candidate_table.update({prod: intersection})
					changed_rows.add(prod)
		return product_candidate_table, changed_rows

	def _compound_merge(self, product_candidate_table, changed_rows):
		model_updates = {}
		while len(changed_rows) > 0:
			curr_prod = changed_rows.pop()
			curr_cands = product_candidate_table[curr_prod]
			if len(curr_cands) == 0: return False
			if len(curr_cands) == 1:
				curr_compound = curr_cands.pop()
				# cue the compound for model update
				model_updates.update({curr_compound: curr_prod})
				# delete the prod from the table
				product_candidate_table.pop(curr_prod)
				# delete the compound from all other candidiate lists
				for other_prod in product_candidate_table:
					other_cands = product_candidate_table[other_prod]
					# catch situations where the only candidate have already
					# been deleted
					if len(other_cands) == 0: return False
					if curr_compound in other_cands:
						other_cands.remove(curr_compound)
						product_candidate_table[other_prod] = other_cands
						changed_rows.add(other_prod)
		return product_candidate_table, model_updates

	def _udpate_states(self, product_candidate_table, model_updates):
		# modify the candidate table
		self._product_candidate_table = product_candidate_table
		# pair the compound
		self._E_compounds_products.update(model_updates)
		# delete the compound from unknowns
		for compound in model_updates: self._uninitialized_E_compounds.remove(compound)
		return

###############################################################################
# Partition_Space Class
###############################################################################

# Partition_Space is the middel-level class that manage the entire Bayesian
# hypothesis space. It keep tracks of each hypothesis instantiated as
# instances of Bayesian_Model
class Partition_Space():
	def __init__(self, compound_model, constrain_compounds = None):
		# A list of Bayesian_Models as hypotheses
		self.hypotheses = None
		# A numpy 2D array of shape (hypothesis, compounds). This is an
		# incidence matrix of inclusion of compounds in each model
		self.hypotheses_incidence = None
		# A probability distribution. Each time the model is updated, the
		# posterior will become the new prior
		self.prior = None
		self._initialize(compound_model, constrain_compounds)

	# Set the prior
	def set_prior(self, prior_func):
		self.prior = prior_func(self.hypotheses[0].compounds, self.hypotheses_incidence)
		return

	# Given a compound and the products, calculate the change in posterior
	# distribution by pushing thie product to the hypotheses. Does not change
	# the state of the models 
	# Returns: 
	# 	- the information gain adjusted by P(products) (float)
	def information_gain(self, compound, products, llh_func, ig_func):
		llh = np.empty(len(self.prior))
		for bm_ind, bm in enumerate(self.hypotheses):
			validate_flag = bm.validate(compound, products, update = False)
			llh[bm_ind] = llh_func(validate_flag, bm)
		marginal = logsumexp(np.add(self.prior, llh))
		posterior = llh + self.prior - marginal
		ig = ig_func(np.exp(self.prior)) - ig_func(np.exp(posterior))
		return np.exp(marginal)*ig

	# A function that update the models given the tested compound and the
	# products. Does change the state of the models
	def bayesian_update(self, tested_compound, observed_products, llh_func):
		llh = np.empty(len(self.prior))
		for bm_ind, bm in enumerate(self.hypotheses):
			validate_flag = bm.validate(tested_compound, observed_products, update = True)
			llh[bm_ind] = llh_func(validate_flag, bm)
		marginal = logsumexp(np.add(self.prior, llh))
		posterior = llh + self.prior - marginal
		self.prior = posterior
		return posterior

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

###############################################################################
# Active_Learning Class
###############################################################################

# Active_Learning is the top-level class that manage the active learning part of
# the model. Below is a list of attributes that the user may find useful:
#	- Active_Learning.compound_model
#	- Active_Learning.partition_space
# Below is a list of functions that the user may find useful:
#	- Active_Learning.active_learning
#	- Active_Learning.choose_compound
#	- Active_learning.learning_episode
class Active_Learning():
	def __init__(self, compound_model, constrain_compounds, prior_func, llh_func, ig_func = entropy):
		# The ground turth compound model
		self.compound_model = compound_model
		# The partition space containing bayesian priors
		self.partition_space = Partition_Space(compound_model, constrain_compounds)
		self.partition_space.set_prior(prior_func)
		self.constrain_compounds = constrain_compounds
		self.constrain_compounds_prod = set({})
		for c in self.constrain_compounds: self.constrain_compounds_prod.add(self.compound_model._prod_to_id[self.constrain_compounds[c]])
		# Likelihood Function
		self.llh_func = llh_func
		# Information Gain Function
		self.ig_func = ig_func
		# The following attributes are constantly being updated
		self.novel_compounds = set({})
		self.compound_history = set({})

	# The all-in-one active learning algorithm. Advance trials and update
	# bayesian hypotheses until a model reaches the posterior of 0.99 or all 11
	# compounds have been observed
	# Returns:
	# 	- The best model (an instance of Bayesian_Model)
	#	- The posterior probability of the best model (float)
	#	- The history of tested compounds (list)
	#	- Posterior entrophy after each update (list)
	#	- Expected information gain for each tested compound (list)
	#	- Actual information gain for each tested compound (list)
	def active_learning(self, exit_threshold = 0.99, return_n = 3, verbose = True):
		# Initialize
		counter = 0
		curr_prior = np.exp(self.partition_space.prior)
		top_prior = np.amax(curr_prior)
		choice_history = []
		expected_IG = [0]
		actual_IG = [0]
		posterior_entropy = [self.ig_func(np.exp(self.partition_space.prior))]
		
		# main learning loop
		while top_prior < exit_threshold:
			counter += 1
			if verbose == True:
				print("\n-------------------------------------------------------------------\nCurrent Trial:", counter,"\n-------------------------------------------------------------------\n")
			# Choose a compound
			curr_compound, IG, compound_ranks, compounds, dist = self.choose_compound(verbose)
			choice_history.append(curr_compound)
			expected_IG.append(IG)
			curr_entropy = self.ig_func(np.exp(self.partition_space.prior))
			actual_IG.append(posterior_entropy[-1] - curr_entropy)
			posterior_entropy.append(curr_entropy)
			if verbose == True: print()
			# Observe the compound's actual product and update the bayesian hypotheses
			models, probs, AIG = self.learning_episode(curr_compound, return_n, verbose)
			top_prior = probs[0]

		# end of loop
		if verbose == True:
			print("\n-------------------------------------------------------------------\nResults")
			print(" - Top Model Probability =", top_prior)
			print(" - Learning Accomplished in", counter, "Trials.")
			print("-------------------------------------------------------------------\n")
		
		return models[0], top_prior, choice_history, posterior_entropy, expected_IG, actual_IG

	# The function for choosing the next compound to test basing on the model
	# history. This is the first half of the active_learning function
	# Returns
	#	- The selected compound for testing (tuple)
	#	- The expected information gain for the said compound (float)
	#   - The ranks of each compound basing on their entrophy (smaller means
	#     less information gain) (np.array)
	#	- all compounds (list)
	#	- Each compound's expected information gain (list)
	def choose_compound(self, verbose = False):
		compounds, dist = self.get_test_dist(verbose)
		compound_ranks = rankdata(dist, method = "max")
		best_compound_inds = np.arange(len(compounds), dtype = int)[compound_ranks == len(compound_ranks)]
		comp_ind = np.random.choice(best_compound_inds)
		test_comp = compounds[comp_ind]
		IG = dist[comp_ind]
		return test_comp, IG, compound_ranks, compounds, dist

	# The function fo updating the bayesian hypotheses. This is the second half
	# of the active_learning function
	# Returns:
	#	- The top return_n bayesian models (a list of Bayesian_Model instances)
	#	- The posterior distribution of hypotheses after the update (np.array)
	#	- The actual informaion gain (float)
	def learning_episode(self, curr_compound, return_n = 3, verbose = False):
		curr_prods = self.compound_model.get_products([curr_compound])[0]
		for prod in curr_prods:
			if prod not in self.constrain_compounds_prod: self.novel_compounds.add(prod)
		self.compound_history.add(curr_compound)
		start_entrophy = self.ig_func(np.exp(self.partition_space.prior))
		self.partition_space.bayesian_update(curr_compound, curr_prods, self.llh_func)
		end_entrophy = self.ig_func(np.exp(self.partition_space.prior))
		actual_IG = start_entrophy - end_entrophy
		if verbose == True:
			self.print_compound_prods(curr_compound)
		models, posterior = self.get_best_models(return_n, verbose)
		return models, posterior, actual_IG
	
	def print_compound_prods(self, compound):
		curr_part = self.compound_model.find_subsets(compound, target = self.compound_model._E_compounds)[0]
		print("================================")
		print("Tested Compound:", compound)
		print("Actual Observations:")
		for par in sorted(curr_part):
			print("  - " + str(self.compound_model._E_compounds_products[par]) + " | " + str(par))

	# A function to get the best num Bayesian_Model instances in the current
	# state
	# Returns:
	#	- Top num hyptheses (a list of Bayesian_Model instances)
	#	- The corresponding posterior probability of the hypotheses (np.array)
	def get_best_models(self, num = 1, verbose = False):
		curr_prior = np.exp(self.partition_space.prior)
		top_n_ind = np.argsort(curr_prior)[::-1][:num]
		if verbose == True:
			for rank, model_ind in enumerate(top_n_ind):
				print("================================")
				print("Model rank:", rank + 1)
				print("Current Probability:", np.around(curr_prior[model_ind], decimals = 4))
				print("Emergent Products:")
				for compound in self.partition_space.hypotheses[model_ind]._E_compounds_products:
					if len(compound) > 1:
						print("  - " + str(self.partition_space.hypotheses[model_ind]._E_compounds_products[compound]) + 
								" | " + str(compound))
		return self.partition_space.hypotheses[top_n_ind], curr_prior[top_n_ind]

	# A function to get the expected information gain of each compound that can 
	# be tested
	# Returns:
	#	- The compounds that can be tested (list of tuples)
	#	- The corresponding information gain for each compounds (list)
	def get_test_dist(self, verbose = False):
		if verbose == False:
			sys.stdout = open(os.devnull, 'w')
		valid_compounds = []
		for c in self.compound_model.compounds:
			if len(c) > 1 and c not in self.compound_history:
				valid_compounds.append(c)
		dist = []
		print("Searching for Compounds to Test")
		for i in tqdm(range(len(valid_compounds)), disable = ~verbose):
			c = valid_compounds[i]
			dist.append(self.compound_expected_utility(c))
		if verbose == False:
			sys.stdout = sys.__stdout__
		return valid_compounds, dist

	# A function that compute the expected information gain of a compound
	# Returns:
	#	- The averaged expected information gain (float)
	def compound_expected_utility(self, compound):
		possible_products = self.possible_products(compound)
		utilities = [self.partition_space.information_gain(compound, pp, self.llh_func, self.ig_func) for pp in possible_products]
		return np.nanmean(utilities)

	# A function that generate possible products of a compound given the current
	# state
	# Returns:
	# 	- Products (list of tuples)
	def possible_products(self, compound):
		basic_products = self.compound_model.find_subsets(compound, target = self.constrain_compounds.keys())[0]
		basic_products = [self.compound_model._prod_to_id[self.constrain_compounds[c]] for c in basic_products]
		additional_products = list(powerset(self.novel_compounds))
		possible_products = []
		for comb in additional_products: possible_products.append(basic_products + list(comb))
		return possible_products

# Helper Functions
###############################################################################

def uniform_prior(compounds, incidences):
	return np.ones(incidences.shape[0])*(-np.log(incidences.shape[0]))

def gamma_prior(compounds, incidences, alpha = 1, beta = 1, scaling = 10):
	all_card = len(compounds)
	all_size = 0
	for c in compounds: all_size += len(c)
	all_x = np.empty(incidences.shape[0])
	for ind, inc in enumerate(incidences):
		curr_card = np.sum(inc)
		curr_compounds = compounds[inc.astype(bool)]
		curr_size = 0
		for c in curr_compounds: curr_size += len(c)
		# Here x is in the range of [0,1]
		all_x[ind] = (alpha*(curr_card/all_card) + beta*(curr_size/all_size))/2
	scaled_x = all_x * scaling
	prior = gamma.logpdf(scaled_x, 1)
	prior = prior - logsumexp(prior)
	return prior

def simple_llh(flag, model):
	if flag == False: 
		return -np.inf
	else:
		return -np.log(len(model._uninitialized_E_compounds) + 1)

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