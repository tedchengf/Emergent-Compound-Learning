import numpy as np
from Models import *

# Global Variables
global elements, products, emergent_compound_products
elements = ["A","B","C","D"]
products = ["0", "♁", "♆", "🜺", "🜍", "♄", "☽", "☉"]
emergent_compound_products = {
	(): products[0],
	("A",): products[1],
	("B",): products[2],
	("C",): products[3],
	("D",): products[4],
	("A","B"): products[5],
	("A", "B", "D"): products[6],
	("B", "C", "D"): products[7]
}
constrain_compounds = {
	(): "0",
	("A",): emergent_compound_products[("A",)],
	("B",): emergent_compound_products[("B",)],
	("C",): emergent_compound_products[("C",)],
	("D",): emergent_compound_products[("D",)]
}

# # Global Variables
# global elements, products, emergent_compound_products
# elements = ["A","B","C"]
# products = ["0", "♁", "♆", "🜺", "🜍", "♄"]
# emergent_compound_products = {
# 	(): products[0],
# 	("A",): products[1],
# 	("B",): products[2],
# 	("C",): products[3],
# 	("A","C"): products[4],
# }
# constrain_compounds = {
# 	(): "0",
# 	("A",): emergent_compound_products[("A",)],
# 	("B",): emergent_compound_products[("B",)],
# 	("C",): emergent_compound_products[("C",)],
# }

# Random Seed
# np.random.seed(2023)

def main():

	ground_truth = Compound_Model(elements, products, emergent_compound_products)
	print(ground_truth._id_to_ind)
	ps = Partition_Space(ground_truth, constrain_compounds)
	ps.set_prior(uniform_prior)
	
	print(len(ps.hypotheses))
	# print("----------")
	# for e_c in ground_truth._E_compounds_products:
	# 	print(e_c, "|", ground_truth._E_compounds_products[e_c])
	# print("----------")
	
	hypotheses_masks = np.ones(len(ps.hypotheses), dtype = bool)
	all_compounds = ground_truth.compounds
	np.random.shuffle(all_compounds)
	for idx, c in enumerate(all_compounds):
		if len(c) <= 1: continue
		actual_products = ground_truth.get_products([c])[0]
		actual_partition = ground_truth.find_subsets(c, target = 
					       ground_truth._E_compounds)[0]
		print("-------------------------------------------------------------------")
		print("-------------------------------------------------------------------")
		print("-------------------------------------------------------------------")
		print("Compound num:", idx)
		print("Tested Compound:", c)
		print("Actual Observations:",)
		for par in sorted(actual_partition):
			print("  - " + str(ground_truth._E_compounds_products[par]) + 
	 				" | " + str(par))
		# for bind, bm in enumerate(ps.hypotheses):
		# 	if hypotheses_masks[bind] == False: continue
		# 	hypotheses_masks[bind] = bm.validate(c, actual_products)
		# 	print("=======")
		# 	print("Model num:", bind)
		# 	print("Current Evaluation:", hypotheses_masks[bind])
		# 	print("Emergent Products:")
		# 	for compound in bm._E_compounds_products:
		# 		if len(compound) > 1:
		# 			print("  - " + str(bm._E_compounds_products[compound]) + 
	 	# 					" | " + str(compound))
		ps.prior = ps.posterior(c, actual_products, simple_llh)
		curr_prior = np.exp(ps.prior)
		top_5_ind = np.argsort(curr_prior)[::-1][:5]
		for rank, model_ind in enumerate(top_5_ind):
			print("=======")
			print("Model rank:", rank + 1)
			print("Current Probability:", curr_prior[model_ind])
			print("Emergent Products:")
			for compound in ps.hypotheses[model_ind]._E_compounds_products:
				if len(compound) > 1:
					print("  - " + str(ps.hypotheses[model_ind]._E_compounds_products[compound]) + 
	 						" | " + str(compound))
	# for bm in ps.hypotheses:

	
	
	return


	test_model_1 = Bayesian_Model(elements, products, emergent_compound_products = {
			(): products[0],
			("A",): products[1],
			("B",): products[2],
			("C",): products[3],
			("D",): products[4],
			("A","B"): None,
			("B","D"): None
			})
	test_model_2 = Bayesian_Model(elements, products, emergent_compound_products = {
			(): products[0],
			("A",): products[1],
			("B",): products[2],
			("C",): products[3],
			("D",): products[4],
			("A","B"): None,
			("A","B","D"): None,
			("B","C","D"): None
			})
			
	compound = ("A", "B", "D")
	actual_products = ground_truth.get_products([compound])[0]
	print(actual_products)
	print("-------------")
	print(test_model_1.validate(compound, actual_products))
	print("-------------")
	print(test_model_2.validate(compound, actual_products))

	compound = ("B", "C", "D")
	actual_products = ground_truth.get_products([compound])[0]
	print(actual_products)
	print("-------------")
	print(test_model_1.validate(compound, actual_products))
	print("-------------")
	print(test_model_2.validate(compound, actual_products))

	# print(ps.hypotheses.shape)
	# print(ground_truth._compounds[ps.hypotheses_incidence[50].astype(bool)])
	# print(ps.hypotheses[10]._uninitialized_E_compounds)
	# print(ground_truth._compounds[ps.hypotheses_incidence[200].astype(bool)])
	# print(ps.hypotheses[200]._uninitialized_E_compounds)
	# print(ground_truth._compounds[ps.hypotheses_incidence[980].astype(bool)])
	# print(ps.hypotheses[980]._uninitialized_E_compounds)
	return

if __name__ == "__main__":
	main()