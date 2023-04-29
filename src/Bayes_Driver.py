import numpy as np
from Models_alt import *

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
# Random Seed
np.random.seed(2023)

def main():


	ground_truth = Compound_Model(elements, products, emergent_compound_products)
	print(ground_truth._id_to_ind)
	ps = Partition_Space(ground_truth, constrain_compounds)
	print("----------")
	for e_c in ground_truth._E_compounds_products:
		print(e_c, "|", ground_truth._E_compounds_products[e_c])
	print("----------")

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
			
	compound = ("A", "B", "C")
	actual_products = ground_truth.get_products([compound])[0]
	print(actual_products)
	np.product(actual_products)
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